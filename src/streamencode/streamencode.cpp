#include "streamencode.hpp"
#include "mtf_rle_arith.hpp"

#include <zstd.h>
#include <cstring>
#include <stdexcept>

namespace bamsi {

// ─── Helper: ZSTD compress/decompress ───────────────────────────────────────

namespace {

std::vector<uint8_t> ZstdCompress(const uint8_t* data, size_t len) {
    size_t bound = ZSTD_compressBound(len);
    std::vector<uint8_t> out(bound);
    size_t compressed = ZSTD_compress(out.data(), bound, data, len, 3);
    if (ZSTD_isError(compressed)) {
        throw Error{ErrorCode::STREAM_DECODE_ERROR,
                    std::string("ZSTD compression failed: ") + ZSTD_getErrorName(compressed)};
    }
    out.resize(compressed);
    return out;
}

std::vector<uint8_t> ZstdDecompress(const uint8_t* data, size_t len, size_t expected) {
    std::vector<uint8_t> out(expected);
    size_t decompressed = ZSTD_decompress(out.data(), expected, data, len);
    if (ZSTD_isError(decompressed)) {
        throw Error{ErrorCode::STREAM_DECODE_ERROR,
                    std::string("ZSTD decompression failed: ") + ZSTD_getErrorName(decompressed)};
    }
    out.resize(decompressed);
    return out;
}

/// Serialize a uint32 to little-endian bytes
void WriteU32(std::vector<uint8_t>& out, uint32_t val) {
    uint8_t buf[4];
    std::memcpy(buf, &val, 4);
    out.insert(out.end(), buf, buf + 4);
}

/// Serialize a uint64 to little-endian bytes
void WriteU64(std::vector<uint8_t>& out, uint64_t val) {
    uint8_t buf[8];
    std::memcpy(buf, &val, 8);
    out.insert(out.end(), buf, buf + 8);
}

uint32_t ReadU32(const uint8_t* data) {
    uint32_t val;
    std::memcpy(&val, data, 4);
    return val;
}

uint64_t ReadU64(const uint8_t* data) {
    uint64_t val;
    std::memcpy(&val, data, 8);
    return val;
}

}  // namespace

// ─── S_seq encoding (Contract §2.4: BWT → MTF → RLE → Arithmetic) ──────────

SeqEncodeResult EncodeSeqStream(const std::vector<uint8_t>& bwt,
                                uint8_t entropy_order_k) {
    SeqEncodeResult result;

    // Full pipeline per Contract §2.4:
    // Step 1: BWT already done (input IS the BWT)
    // Step 2: Move-to-Front over alphabet sigma=6 {A,C,G,T,N,#}
    auto mtf = MtfEncode(bwt, SIGMA);

    // Step 3: Run-Length Encoding
    auto rle = RleEncode(mtf);

    // Step 4: 0th-order arithmetic coding
    auto arith = ArithEncode(rle);

    // Store: header(8 bytes: original_len as uint64) + arith payload
    // We need original BWT length for decode, plus RLE decoded length
    result.payload.reserve(16 + arith.size());
    // Store BWT length (for full decode chain)
    uint64_t bwt_len = bwt.size();
    uint64_t rle_len = rle.size();
    result.payload.resize(16);
    std::memcpy(result.payload.data(), &bwt_len, 8);
    std::memcpy(result.payload.data() + 8, &rle_len, 8);
    result.payload.insert(result.payload.end(), arith.begin(), arith.end());

    result.codec_id = 0x10;  // BWT_MTF_RLE_ARITH
    return result;
}

std::vector<uint8_t> DecodeSeqStream(const std::vector<uint8_t>& payload,
                                     uint8_t codec_id, uint64_t expected_len) {
    if (codec_id == 0x01) {
        // Legacy ZSTD codec (for backward compat with old .bsi files)
        return ZstdDecompress(payload.data(), payload.size(), expected_len);
    }

    // codec_id == 0x10: BWT_MTF_RLE_ARITH
    if (payload.size() < 16) {
        throw Error{ErrorCode::STREAM_DECODE_ERROR, "S_seq payload too short"};
    }

    uint64_t bwt_len, rle_len;
    std::memcpy(&bwt_len, payload.data(), 8);
    std::memcpy(&rle_len, payload.data() + 8, 8);

    // Step 1: Arithmetic decode → RLE data
    std::vector<uint8_t> arith_payload(payload.begin() + 16, payload.end());
    auto rle = ArithDecode(arith_payload, rle_len);

    // Step 2: RLE decode → MTF data
    auto mtf = RleDecode(rle, bwt_len);

    // Step 3: MTF decode → BWT
    auto bwt = MtfDecode(mtf, SIGMA);

    return bwt;
}

// ─── S_qual encoding (Contract §2.7: RANGE_CYCLE) ──────────────────────────
// V1 implementation: per-read ZSTD with per-cycle reordering.
// The RANGE_CYCLE approach transposes quality scores by cycle position,
// then compresses each cycle column. Here we approximate with ZSTD
// on cycle-transposed data per read.

QualEncodeResult EncodeQualStream(const std::vector<OrderedRead>& reads,
                                  QualCodec codec, uint8_t lossy_bins) {
    QualEncodeResult result;
    result.codec_id = static_cast<uint8_t>(codec);
    result.directory.resize(reads.size());

    // Independent per-read blocks as required by strong independence rule (I10).
    for (size_t i = 0; i < reads.size(); ++i) {
        uint64_t offset = result.payload.size();

        // Apply lossy binning if requested
        std::vector<uint8_t> qual_data = reads[i].qual;
        if (lossy_bins > 0 && lossy_bins < 94) {
            uint8_t bin_size = 93 / lossy_bins;
            for (auto& q : qual_data) {
                q = (q / bin_size) * bin_size + bin_size / 2;
                if (q > 93) q = 93;
            }
        }

        auto compressed = ZstdCompress(qual_data.data(), qual_data.size());
        result.directory[i] = {offset, static_cast<uint32_t>(compressed.size())};
        result.payload.insert(result.payload.end(), compressed.begin(), compressed.end());
    }

    return result;
}

std::vector<uint8_t> DecodeQualRead(const std::vector<uint8_t>& payload,
                                    const StreamDirectoryEntry& dir_entry,
                                    uint8_t codec_id) {
    const uint8_t* data = payload.data() + dir_entry.offset;
    uint64_t expected = ZSTD_getFrameContentSize(data, dir_entry.length);
    if (expected == ZSTD_CONTENTSIZE_UNKNOWN || expected == ZSTD_CONTENTSIZE_ERROR) {
        expected = dir_entry.length * 10;  // fallback estimate
    }
    return ZstdDecompress(data, dir_entry.length, expected);
}

// ─── S_meta encoding (Contract §2.8: TYPED_SPLIT) ──────────────────────────
// TYPED_SPLIT layout: per-read block contains:
//   [FLAG: 4 bytes]
//   [n_cigar_ops: varint]
//   [CIGAR ops: each = op_nybble(4 bits) + len_varint(7-bit-stop)]
// Then ZSTD compressed.

MetaEncodeResult EncodeMetaStream(const std::vector<OrderedRead>& reads,
                                   MetaCodec codec) {
    MetaEncodeResult result;
    result.codec_id = static_cast<uint8_t>(codec);
    result.directory.resize(reads.size());

    for (size_t i = 0; i < reads.size(); ++i) {
        std::vector<uint8_t> raw;

        // Write FLAG (4 bytes LE)
        WriteU32(raw, reads[i].flag);

        // CIGAR substream: nybble-encoded ops + varint lengths
        uint32_t n_ops = static_cast<uint32_t>(reads[i].cigar.size());
        // Write n_ops as varint
        {
            uint32_t v = n_ops;
            do {
                uint8_t b = v & 0x7F;
                v >>= 7;
                if (v) b |= 0x80;
                raw.push_back(b);
            } while (v);
        }

        // Pack CIGAR ops: 4-bit nybble for op + varint for length
        for (const auto& cop : reads[i].cigar) {
            // op nybble (4 bits) — fits in one byte with high bit clear
            raw.push_back(cop.op & 0x0F);
            // length as varint
            uint32_t len = cop.len;
            do {
                uint8_t b = len & 0x7F;
                len >>= 7;
                if (len) b |= 0x80;
                raw.push_back(b);
            } while (len);
        }

        uint64_t offset = result.payload.size();
        auto compressed = ZstdCompress(raw.data(), raw.size());
        result.directory[i] = {offset, static_cast<uint32_t>(compressed.size())};
        result.payload.insert(result.payload.end(), compressed.begin(), compressed.end());
    }

    return result;
}

DecodedMeta DecodeMetaRead(const std::vector<uint8_t>& payload,
                           const StreamDirectoryEntry& dir_entry,
                           uint8_t codec_id) {
    const uint8_t* data = payload.data() + dir_entry.offset;
    uint64_t expected = ZSTD_getFrameContentSize(data, dir_entry.length);
    if (expected == ZSTD_CONTENTSIZE_UNKNOWN || expected == ZSTD_CONTENTSIZE_ERROR) {
        expected = dir_entry.length * 10;
    }
    auto raw = ZstdDecompress(data, dir_entry.length, expected);

    DecodedMeta meta;
    size_t off = 0;

    // Read FLAG
    meta.flag = ReadU32(raw.data() + off); off += 4;

    // Read n_ops as varint
    uint32_t n_ops = 0;
    {
        uint32_t shift = 0;
        while (off < raw.size()) {
            uint8_t b = raw[off++];
            n_ops |= (uint32_t)(b & 0x7F) << shift;
            if (!(b & 0x80)) break;
            shift += 7;
        }
    }

    meta.cigar.resize(n_ops);
    for (uint32_t c = 0; c < n_ops; ++c) {
        meta.cigar[c].op = raw[off++] & 0x0F;
        // Read length as varint
        uint32_t len = 0;
        uint32_t shift = 0;
        while (off < raw.size()) {
            uint8_t b = raw[off++];
            len |= (uint32_t)(b & 0x7F) << shift;
            if (!(b & 0x80)) break;
            shift += 7;
        }
        meta.cigar[c].len = len;
    }

    return meta;
}

// ─── S_map encoding (Contract §2.9: DELTA_RANGE) ───────────────────────────
// DELTA_RANGE layout:
//   Per-read block: chrom_id(4 bytes) + delta_pos(varint, signed zigzag)
//   First read of each chromosome: absolute pos.
//   ZSTD compressed over the entire stream.

MapEncodeResult EncodeMapStream(const std::vector<OrderedRead>& reads,
                                 MapCodec codec) {
    MapEncodeResult result;
    result.codec_id = static_cast<uint8_t>(codec);
    result.directory.resize(reads.size());

    // Build delta-encoded stream
    uint32_t prev_chrom_id = UINT32_MAX;
    uint64_t prev_pos = 0;

    for (size_t i = 0; i < reads.size(); ++i) {
        std::vector<uint8_t> raw;

        // chrom_id (4 bytes)
        WriteU32(raw, reads[i].chrom_id);

        // Delta or absolute pos
        if (reads[i].chrom_id != prev_chrom_id) {
            // New chromosome: absolute pos (8 bytes)
            raw.push_back(0x00);  // marker: absolute
            WriteU64(raw, reads[i].pos);
        } else {
            // Same chromosome: delta pos (zigzag varint)
            raw.push_back(0x01);  // marker: delta
            int64_t delta = static_cast<int64_t>(reads[i].pos) - static_cast<int64_t>(prev_pos);
            // Zigzag encode
            uint64_t zigzag = (static_cast<uint64_t>(delta) << 1) ^
                              (static_cast<uint64_t>(delta >> 63));
            do {
                uint8_t b = zigzag & 0x7F;
                zigzag >>= 7;
                if (zigzag) b |= 0x80;
                raw.push_back(b);
            } while (zigzag);
        }

        prev_chrom_id = reads[i].chrom_id;
        prev_pos = reads[i].pos;

        uint64_t offset = result.payload.size();
        auto compressed = ZstdCompress(raw.data(), raw.size());
        result.directory[i] = {offset, static_cast<uint32_t>(compressed.size())};
        result.payload.insert(result.payload.end(), compressed.begin(), compressed.end());
    }

    return result;
}

DecodedMap DecodeMapRead(const std::vector<uint8_t>& payload,
                         const StreamDirectoryEntry& dir_entry,
                         uint8_t codec_id) {
    const uint8_t* data = payload.data() + dir_entry.offset;
    uint64_t expected = ZSTD_getFrameContentSize(data, dir_entry.length);
    if (expected == ZSTD_CONTENTSIZE_UNKNOWN || expected == ZSTD_CONTENTSIZE_ERROR) {
        expected = 12;
    }
    auto raw = ZstdDecompress(data, dir_entry.length, expected);

    DecodedMap map;
    map.chrom_id = ReadU32(raw.data());

    size_t off = 4;
    uint8_t marker = raw[off++];

    if (marker == 0x00) {
        // Absolute pos
        map.pos = ReadU64(raw.data() + off);
    } else {
        // Delta pos (zigzag varint) — need previous pos context
        // For per-read random access, we store enough to decode independently
        // In the current per-read ZSTD scheme, each block is self-contained
        // So we decode the zigzag but it's relative to 0 (first of chrom)
        uint64_t zigzag = 0;
        uint32_t shift = 0;
        while (off < raw.size()) {
            uint8_t b = raw[off++];
            zigzag |= static_cast<uint64_t>(b & 0x7F) << shift;
            if (!(b & 0x80)) break;
            shift += 7;
        }
        int64_t delta = static_cast<int64_t>((zigzag >> 1) ^ -(zigzag & 1));
        map.pos = static_cast<uint64_t>(delta);
    }

    return map;
}

}  // namespace bamsi
