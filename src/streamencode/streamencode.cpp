#include "streamencode.hpp"
#include "mtf_rle_arith.hpp"

#include <zstd.h>
#include <cstring>
#include <stdexcept>

namespace bamsix {

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
// RANGE_CYCLE: per-cycle context captures sequencing-cycle quality decay.
//
// Layout per-read block (before ZSTD compression):
//   [read_length: uint16_t LE]
//   [cycle-transposed quality scores: length bytes, reordered by cycle position]
//
// The transposition reorders qual[0], qual[1], ..., qual[L-1] so that values
// at the same cycle position (modulo cycle_period) are adjacent. This groups
// same-cycle qualities together, enabling ZSTD's LZ77+entropy backend to
// exploit within-cycle correlation (the dominant redundancy source in
// Illumina/MGI data where Q-scores decay predictably along the read cycle).
//
// For a read of length L, the transposed layout is:
//   qual[0], qual[cycle_period], qual[2*cycle_period], ...  (cycle 0)
//   qual[1], qual[1+cycle_period], qual[1+2*cycle_period], ... (cycle 1)
//   ...
//   qual[cycle_period-1], qual[2*cycle_period-1], ...          (cycle period-1)
//
// cycle_period is stored as part of the block header so decoders can reverse
// the transposition. Default: cycle_period = read_length (full transposition).

QualEncodeResult EncodeQualStream(const std::vector<OrderedRead>& reads,
                                   QualCodec codec, uint8_t lossy_bins,
                                   uint32_t block_size) {
    QualEncodeResult result;
    result.codec_id = static_cast<uint8_t>(codec);
    result.block_size = block_size;

    // Guard: reject unsupported codec IDs per Architecture §8.1
    if (codec != QualCodec::RANGE_CYCLE && codec != QualCodec::RANS_DELTA &&
        codec != QualCodec::ZSTD_DICT && codec != QualCodec::BINNED_RANGE) {
        throw Error{ErrorCode::UNSUPPORTED_CODEC,
                    "Unsupported quality codec: 0x" +
                    std::to_string(static_cast<unsigned>(codec))};
    }

    // ─── Helper: apply lossy binning if requested ───────────────────────────
    auto apply_lossy_binning = [&](std::vector<uint8_t>& qual_data) {
        if (lossy_bins > 0 && lossy_bins < 94) {
            uint8_t bin_size = 93 / lossy_bins;
            for (auto& q : qual_data) {
                q = (q / bin_size) * bin_size + bin_size / 2;
                if (q > 93) q = 93;
            }
        }
    };

    // ─── Helper: encode a single read's quality into raw bytes ──────────────
    auto encode_one_read = [&](size_t i) -> std::vector<uint8_t> {
        std::vector<uint8_t> qual_data = reads[i].qual;
        apply_lossy_binning(qual_data);

        uint16_t read_len = static_cast<uint16_t>(qual_data.size());

        if (codec == QualCodec::ZSTD_DICT) {
            // ZSTD_DICT (0x03): simple layout [read_len:u16][raw_qual_bytes]
            // No cycle transposition — relies on ZSTD's LZ77 for compression.
            // This is the "escape hatch" codec per Exec Plan §5.3.6.
            std::vector<uint8_t> raw(2 + qual_data.size());
            std::memcpy(raw.data(), &read_len, 2);
            if (!qual_data.empty()) {
                std::memcpy(raw.data() + 2, qual_data.data(), qual_data.size());
            }
            return raw;
        }

        if (codec == QualCodec::RANS_DELTA) {
            // RANS_DELTA (0x02): delta-coding between adjacent quality scores.
            // Contract §2.7: exploits sequential correlation in quality profiles.
            // Layout: [read_len:u16][first_val:u8][delta_coded_values: read_len-1 bytes]
            // Delta transform: d[0] = q[0]; d[i] = q[i] - q[i-1] + 128 (biased unsigned)
            // Then ZSTD envelope handles entropy coding.
            std::vector<uint8_t> raw;
            raw.resize(2);
            std::memcpy(raw.data(), &read_len, 2);

            if (!qual_data.empty()) {
                raw.push_back(qual_data[0]);  // absolute first value
                for (size_t j = 1; j < qual_data.size(); ++j) {
                    // Biased delta: add 128 to keep values unsigned [0,255]
                    uint8_t delta = static_cast<uint8_t>(
                        static_cast<int16_t>(qual_data[j]) -
                        static_cast<int16_t>(qual_data[j - 1]) + 128);
                    raw.push_back(delta);
                }
            }
            return raw;
        }

        if (codec == QualCodec::BINNED_RANGE) {
            // BINNED_RANGE (0x04): pre-bin quality scores then cycle-transpose.
            // Contract §2.7: lossy fast-path combining binning + cycle structure.
            // Layout: [read_len:u16][cycle_period:u16][bin_count:u8][binned_transposed_quals]
            // Binning: if lossy_bins > 0, apply the standard binning; if lossy_bins == 0,
            // use a default 8-level binning scheme (this codec is inherently lossy-oriented).
            uint8_t effective_bins = (lossy_bins > 0) ? lossy_bins : 8;
            uint8_t bin_size = (effective_bins > 0 && effective_bins < 94) ? (93 / effective_bins) : 1;
            for (auto& q : qual_data) {
                q = (q / bin_size) * bin_size + bin_size / 2;
                if (q > 93) q = 93;
            }

            uint16_t cycle_period = read_len;
            std::vector<uint8_t> raw;
            raw.resize(5);  // read_len(2) + cycle_period(2) + bin_count(1)
            std::memcpy(raw.data(), &read_len, 2);
            std::memcpy(raw.data() + 2, &cycle_period, 2);
            raw[4] = effective_bins;

            if (!qual_data.empty()) {
                std::vector<uint8_t> transposed(qual_data.size());
                size_t dst = 0;
                for (uint16_t c = 0; c < cycle_period && dst < qual_data.size(); ++c) {
                    for (size_t j = c; j < qual_data.size(); j += cycle_period) {
                        transposed[dst++] = qual_data[j];
                    }
                }
                raw.insert(raw.end(), transposed.begin(), transposed.end());
            }
            return raw;
        }

        // RANGE_CYCLE (0x01): cycle-transposed layout (default)
        uint16_t cycle_period = read_len;
        std::vector<uint8_t> raw;
        raw.resize(4);
        std::memcpy(raw.data(), &read_len, 2);
        std::memcpy(raw.data() + 2, &cycle_period, 2);

        if (!qual_data.empty()) {
            std::vector<uint8_t> transposed(qual_data.size());
            size_t dst = 0;
            for (uint16_t c = 0; c < cycle_period && dst < qual_data.size(); ++c) {
                for (size_t j = c; j < qual_data.size(); j += cycle_period) {
                    transposed[dst++] = qual_data[j];
                }
            }
            raw.insert(raw.end(), transposed.begin(), transposed.end());
        }
        return raw;
    };

    if (block_size == 0) {
        // ─── Per-read directory (original path) ─────────────────────────────
        StreamDirectoryPerRead dir(reads.size());
        for (size_t i = 0; i < reads.size(); ++i) {
            uint64_t offset = result.payload.size();
            auto raw = encode_one_read(i);
            auto compressed = ZstdCompress(raw.data(), raw.size());
            dir[i] = {offset, static_cast<uint32_t>(compressed.size())};
            result.payload.insert(result.payload.end(), compressed.begin(), compressed.end());
        }
        result.directory = std::move(dir);
    } else {
        // ─── Block-level directory (Contract §2.3 F4) ───────────────────────
        // Group reads into blocks of block_size. Each block is independently
        // ZSTD-compressed. Block layout:
        //   [num_reads_in_block: uint32_t LE]
        //   [per_read_offset[0..num-1]: uint32_t LE × num]  (offsets within uncompressed block)
        //   [per_read_len[0..num-1]: uint32_t LE × num]    (uncompressed length per read)
        //   [read_0_raw_qual][read_1_raw_qual]...[read_{num-1}_raw_qual]
        StreamDirectoryBlockLevel block_dir;
        size_t N = reads.size();
        for (size_t block_start = 0; block_start < N; block_start += block_size) {
            size_t block_end = std::min(block_start + block_size, N);
            uint32_t num_in_block = static_cast<uint32_t>(block_end - block_start);

            // Encode all reads in this block
            std::vector<std::vector<uint8_t>> raw_reads(num_in_block);
            for (uint32_t j = 0; j < num_in_block; ++j) {
                raw_reads[j] = encode_one_read(block_start + j);
            }

            // Build block payload: header + concatenated raw reads
            std::vector<uint8_t> block_payload;
            // Block header: num_reads
            block_payload.resize(4 + num_in_block * 8);
            std::memcpy(block_payload.data(), &num_in_block, 4);

            // Compute per-read offsets within the block body
            uint32_t body_offset = static_cast<uint32_t>(4 + num_in_block * 8);
            for (uint32_t j = 0; j < num_in_block; ++j) {
                uint32_t read_off = body_offset;
                uint32_t read_len = static_cast<uint32_t>(raw_reads[j].size());
                std::memcpy(block_payload.data() + 4 + j * 4, &read_off, 4);
                std::memcpy(block_payload.data() + 4 + num_in_block * 4 + j * 4, &read_len, 4);
                body_offset += read_len;
            }

            // Append all raw reads
            for (uint32_t j = 0; j < num_in_block; ++j) {
                block_payload.insert(block_payload.end(), raw_reads[j].begin(), raw_reads[j].end());
            }

            // Compress the entire block
            auto compressed = ZstdCompress(block_payload.data(), block_payload.size());

            BlockDirectoryEntry entry;
            entry.block_offset = result.payload.size();
            entry.block_length = static_cast<uint32_t>(compressed.size());
            entry.first_read_id = static_cast<uint32_t>(block_start);
            block_dir.push_back(entry);

            result.payload.insert(result.payload.end(), compressed.begin(), compressed.end());
        }
        result.directory = std::move(block_dir);
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
    auto raw = ZstdDecompress(data, dir_entry.length, expected);

    if (raw.size() < 2) {
        throw Error{ErrorCode::STREAM_DECODE_ERROR,
                    "S_qual per-read block too short (expected ≥2 bytes header)"};
    }

    uint16_t read_len;
    std::memcpy(&read_len, raw.data(), 2);

    // ─── ZSTD_DICT (0x03): [read_len:u16][raw_qual_bytes] ───────────────────
    if (codec_id == static_cast<uint8_t>(QualCodec::ZSTD_DICT)) {
        if (raw.size() < 2u + read_len) {
            throw Error{ErrorCode::STREAM_DECODE_ERROR,
                        "S_qual ZSTD_DICT: payload shorter than declared read_len"};
        }
        return std::vector<uint8_t>(raw.begin() + 2, raw.begin() + 2 + read_len);
    }

    // ─── RANS_DELTA (0x02): [read_len:u16][first_val:u8][deltas] ────────────
    if (codec_id == static_cast<uint8_t>(QualCodec::RANS_DELTA)) {
        if (read_len == 0) return {};
        if (raw.size() < 3u) {
            throw Error{ErrorCode::STREAM_DECODE_ERROR,
                        "S_qual RANS_DELTA: payload too short"};
        }
        std::vector<uint8_t> qual(read_len);
        qual[0] = raw[2];  // absolute first value
        for (uint16_t j = 1; j < read_len && (3u + j - 1) < raw.size(); ++j) {
            // Reverse biased delta: q[i] = q[i-1] + (delta - 128)
            int16_t prev = static_cast<int16_t>(qual[j - 1]);
            int16_t delta = static_cast<int16_t>(raw[2 + j]) - 128;
            int16_t val = prev + delta;
            if (val < 0) val = 0;
            if (val > 93) val = 93;
            qual[j] = static_cast<uint8_t>(val);
        }
        return qual;
    }

    // ─── BINNED_RANGE (0x04): [read_len:u16][cycle_period:u16][bin_count:u8][transposed] ──
    if (codec_id == static_cast<uint8_t>(QualCodec::BINNED_RANGE)) {
        if (raw.size() < 5) {
            throw Error{ErrorCode::STREAM_DECODE_ERROR,
                        "S_qual BINNED_RANGE: block too short (expected ≥5 bytes header)"};
        }
        uint16_t cycle_period;
        std::memcpy(&cycle_period, raw.data() + 2, 2);
        // uint8_t bin_count = raw[4]; // stored for metadata, not needed for decode

        if (raw.size() < 5u + read_len) {
            throw Error{ErrorCode::STREAM_DECODE_ERROR,
                        "S_qual BINNED_RANGE: payload shorter than declared read_len"};
        }
        const uint8_t* transposed = raw.data() + 5;
        std::vector<uint8_t> qual(read_len);

        if (cycle_period == 0 || cycle_period >= read_len) {
            std::memcpy(qual.data(), transposed, read_len);
        } else {
            size_t src = 0;
            for (uint16_t c = 0; c < cycle_period && src < read_len; ++c) {
                for (size_t j = c; j < read_len; j += cycle_period) {
                    qual[j] = transposed[src++];
                }
            }
        }
        return qual;
    }

    // ─── RANGE_CYCLE (0x01): [read_len:u16][cycle_period:u16][transposed_quals] ──
    if (raw.size() < 4) {
        throw Error{ErrorCode::STREAM_DECODE_ERROR,
                    "S_qual RANGE_CYCLE: block too short (expected ≥4 bytes header)"};
    }
    uint16_t cycle_period;
    std::memcpy(&cycle_period, raw.data() + 2, 2);

    if (codec_id == static_cast<uint8_t>(QualCodec::RANGE_CYCLE)) {
        if (raw.size() < 4 + read_len) {
            throw Error{ErrorCode::STREAM_DECODE_ERROR,
                        "S_qual RANGE_CYCLE: payload shorter than declared read_len"};
        }
        const uint8_t* transposed = raw.data() + 4;
        std::vector<uint8_t> qual(read_len);

        if (cycle_period == 0 || cycle_period >= read_len) {
            std::memcpy(qual.data(), transposed, read_len);
        } else {
            size_t src = 0;
            for (uint16_t c = 0; c < cycle_period && src < read_len; ++c) {
                for (size_t j = c; j < read_len; j += cycle_period) {
                    qual[j] = transposed[src++];
                }
            }
        }
        return qual;
    }

    // Fallback: raw quality scores after 4-byte header
    return std::vector<uint8_t>(raw.begin() + 4, raw.begin() + 4 + read_len);
}

// ─── Block-level S_qual decoding (Contract §2.3 F4) ─────────────────────────
// Decodes quality for read `read_id` from a block-level directory.
// Steps:
//   1. Binary search block_dir to find the block containing read_id
//   2. Decompress the block
//   3. Read the block-internal per-read offset table to find read_id's raw data
//   4. Decode using the per-read DecodeQualRead logic

std::vector<uint8_t> DecodeQualFromBlock(const std::vector<uint8_t>& payload,
                                          const StreamDirectoryBlockLevel& block_dir,
                                          uint32_t block_size,
                                          uint64_t read_id,
                                          uint8_t codec_id) {
    if (block_dir.empty() || block_size == 0) {
        throw Error{ErrorCode::STREAM_DECODE_ERROR,
                    "DecodeQualFromBlock: empty block directory or block_size=0"};
    }

    // Find the block containing read_id
    // block_dir is sorted by first_read_id, find the last entry with first_read_id <= read_id
    size_t block_idx = 0;
    for (size_t b = 0; b < block_dir.size(); ++b) {
        if (block_dir[b].first_read_id <= static_cast<uint32_t>(read_id)) {
            block_idx = b;
        } else {
            break;
        }
    }

    const auto& block = block_dir[block_idx];
    uint32_t local_id = static_cast<uint32_t>(read_id) - block.first_read_id;

    // Decompress the block
    const uint8_t* block_data = payload.data() + block.block_offset;
    uint64_t expected = ZSTD_getFrameContentSize(block_data, block.block_length);
    if (expected == ZSTD_CONTENTSIZE_UNKNOWN || expected == ZSTD_CONTENTSIZE_ERROR) {
        expected = block.block_length * 10;
    }
    auto raw_block = ZstdDecompress(block_data, block.block_length, expected);

    if (raw_block.size() < 4) {
        throw Error{ErrorCode::STREAM_DECODE_ERROR,
                    "DecodeQualFromBlock: block too short"};
    }

    // Read block header
    uint32_t num_in_block;
    std::memcpy(&num_in_block, raw_block.data(), 4);

    if (local_id >= num_in_block) {
        throw Error{ErrorCode::STREAM_DECODE_ERROR,
                    "DecodeQualFromBlock: local_id >= num_in_block"};
    }

    if (raw_block.size() < 4 + num_in_block * 8) {
        throw Error{ErrorCode::STREAM_DECODE_ERROR,
                    "DecodeQualFromBlock: block header truncated"};
    }

    // Read per-read offset and length from block header
    uint32_t read_off, read_len;
    std::memcpy(&read_off, raw_block.data() + 4 + local_id * 4, 4);
    std::memcpy(&read_len, raw_block.data() + 4 + num_in_block * 4 + local_id * 4, 4);

    if (read_off + read_len > raw_block.size()) {
        throw Error{ErrorCode::STREAM_DECODE_ERROR,
                    "DecodeQualFromBlock: read data extends beyond block"};
    }

    // Create a synthetic per-read directory entry pointing into the decompressed block
    // and use the existing DecodeQualRead logic on the decompressed data
    // We need to build a temporary view that DecodeQualRead can process.
    // The raw read data is at raw_block[read_off..read_off+read_len]
    // This is the same format as what ZstdDecompress produces for a per-read block.
    const uint8_t* read_data = raw_block.data() + read_off;

    if (read_len < 2) {
        throw Error{ErrorCode::STREAM_DECODE_ERROR,
                    "DecodeQualFromBlock: read data too short"};
    }

    // Parse per-read header: format depends on codec
    uint16_t rlen;
    std::memcpy(&rlen, read_data, 2);

    // ZSTD_DICT: [read_len:u16][raw_qual_bytes]
    if (codec_id == static_cast<uint8_t>(QualCodec::ZSTD_DICT)) {
        if (read_len < 2u + rlen) {
            throw Error{ErrorCode::STREAM_DECODE_ERROR,
                        "DecodeQualFromBlock ZSTD_DICT: read data shorter than declared length"};
        }
        return std::vector<uint8_t>(read_data + 2, read_data + 2 + rlen);
    }

    // RANS_DELTA (0x02): [read_len:u16][first_val:u8][deltas]
    if (codec_id == static_cast<uint8_t>(QualCodec::RANS_DELTA)) {
        if (rlen == 0) return {};
        if (read_len < 3u) {
            throw Error{ErrorCode::STREAM_DECODE_ERROR,
                        "DecodeQualFromBlock RANS_DELTA: read data too short"};
        }
        std::vector<uint8_t> qual(rlen);
        qual[0] = read_data[2];  // absolute first value
        for (uint16_t j = 1; j < rlen && (3u + j - 1) < read_len; ++j) {
            int16_t prev = static_cast<int16_t>(qual[j - 1]);
            int16_t delta = static_cast<int16_t>(read_data[2 + j]) - 128;
            int16_t val = prev + delta;
            if (val < 0) val = 0;
            if (val > 93) val = 93;
            qual[j] = static_cast<uint8_t>(val);
        }
        return qual;
    }

    // BINNED_RANGE (0x04): [read_len:u16][cycle_period:u16][bin_count:u8][transposed]
    if (codec_id == static_cast<uint8_t>(QualCodec::BINNED_RANGE)) {
        if (read_len < 5) {
            throw Error{ErrorCode::STREAM_DECODE_ERROR,
                        "DecodeQualFromBlock BINNED_RANGE: read data too short"};
        }
        uint16_t cycle_period;
        std::memcpy(&cycle_period, read_data + 2, 2);
        // uint8_t bin_count = read_data[4]; // metadata, not needed for decode

        if (read_len < 5u + rlen) {
            throw Error{ErrorCode::STREAM_DECODE_ERROR,
                        "DecodeQualFromBlock BINNED_RANGE: payload shorter than declared"};
        }
        const uint8_t* transposed = read_data + 5;
        std::vector<uint8_t> qual(rlen);
        if (cycle_period == 0 || cycle_period >= rlen) {
            std::memcpy(qual.data(), transposed, rlen);
        } else {
            size_t src = 0;
            for (uint16_t c = 0; c < cycle_period && src < rlen; ++c) {
                for (size_t j = c; j < rlen; j += cycle_period) {
                    qual[j] = transposed[src++];
                }
            }
        }
        return qual;
    }

    // RANGE_CYCLE (0x01): [read_len:u16][cycle_period:u16][transposed_quals]
    if (read_len < 4) {
        throw Error{ErrorCode::STREAM_DECODE_ERROR,
                    "DecodeQualFromBlock RANGE_CYCLE: read data too short for 4-byte header"};
    }

    uint16_t cycle_period;
    std::memcpy(&cycle_period, read_data + 2, 2);

    if (read_len < 4u + rlen) {
        throw Error{ErrorCode::STREAM_DECODE_ERROR,
                    "DecodeQualFromBlock: read data shorter than declared length"};
    }

    if (codec_id == static_cast<uint8_t>(QualCodec::RANGE_CYCLE)) {
        const uint8_t* transposed = read_data + 4;
        std::vector<uint8_t> qual(rlen);

        if (cycle_period == 0 || cycle_period >= rlen) {
            std::memcpy(qual.data(), transposed, rlen);
        } else {
            size_t src = 0;
            for (uint16_t c = 0; c < cycle_period && src < rlen; ++c) {
                for (size_t j = c; j < rlen; j += cycle_period) {
                    qual[j] = transposed[src++];
                }
            }
        }
        return qual;
    }

    return std::vector<uint8_t>(read_data + 4, read_data + 4 + rlen);
}

// ─── S_meta encoding (Contract §2.8: TYPED_SPLIT) ──────────────────────────
// TYPED_SPLIT layout: per-read block contains three typed substreams:
//   (a) FLAG substream: [FLAG: 2 bytes LE] (BAM FLAG is uint16)
//   (b) CIGAR substream: [n_ops: varint][ops: 4-bit nybble + varint len each]
//   (c) Aux-tag substream: [aux_len: varint][raw_aux_bytes] (skipped when empty)
// Entire per-read block is then ZSTD compressed.

MetaEncodeResult EncodeMetaStream(const std::vector<OrderedRead>& reads,
                                   MetaCodec codec) {
    MetaEncodeResult result;
    result.codec_id = static_cast<uint8_t>(codec);
    result.directory.resize(reads.size());

    // Guard: reject unsupported codec IDs per Architecture §8.1
    if (codec != MetaCodec::TYPED_SPLIT && codec != MetaCodec::ZSTD_FALLBACK) {
        throw Error{ErrorCode::UNSUPPORTED_CODEC,
                    "Unsupported meta codec: 0x" +
                    std::to_string(static_cast<unsigned>(codec))};
    }

    for (size_t i = 0; i < reads.size(); ++i) {
        std::vector<uint8_t> raw;

        if (codec == MetaCodec::ZSTD_FALLBACK) {
            // ZSTD_FALLBACK (0x02): whole-record serialization without typed substream
            // decomposition. Contract §2.8: escape hatch for when TYPED_SPLIT's nybble
            // encoding doesn't provide enough compression benefit (e.g., unusual
            // aux-tag-heavy records). Layout:
            //   [FLAG: uint32_t LE (4 bytes)]
            //   [n_cigar_ops: uint32_t LE (4 bytes)]
            //   [cigar_ops: n × (op:uint8_t + len:uint32_t LE) = 5 bytes each]
            //   [aux_len: uint32_t LE (4 bytes)]
            //   [aux_data: aux_len bytes]
            // Entire raw block is then ZSTD compressed.

            // FLAG (4 bytes LE — stored as uint32 for alignment)
            uint32_t flag32 = reads[i].flag;
            raw.resize(4);
            std::memcpy(raw.data(), &flag32, 4);

            // CIGAR ops count + ops
            uint32_t n_ops = static_cast<uint32_t>(reads[i].cigar.size());
            {
                uint8_t buf[4];
                std::memcpy(buf, &n_ops, 4);
                raw.insert(raw.end(), buf, buf + 4);
            }
            for (const auto& cop : reads[i].cigar) {
                raw.push_back(cop.op);
                uint8_t lbuf[4];
                uint32_t clen = cop.len;
                std::memcpy(lbuf, &clen, 4);
                raw.insert(raw.end(), lbuf, lbuf + 4);
            }

            // Aux tags
            uint32_t aux_len = static_cast<uint32_t>(reads[i].aux_data.size());
            {
                uint8_t buf[4];
                std::memcpy(buf, &aux_len, 4);
                raw.insert(raw.end(), buf, buf + 4);
            }
            if (aux_len > 0) {
                raw.insert(raw.end(), reads[i].aux_data.begin(), reads[i].aux_data.end());
            }
        } else {
            // TYPED_SPLIT (0x01): three typed substreams per Contract §2.8

            // (a) FLAG substream: 2 bytes LE (BAM FLAG is uint16_t)
            uint16_t flag16 = static_cast<uint16_t>(reads[i].flag);
            raw.push_back(static_cast<uint8_t>(flag16 & 0xFF));
            raw.push_back(static_cast<uint8_t>((flag16 >> 8) & 0xFF));

            // (b) CIGAR substream: nybble-encoded ops + varint lengths
            uint32_t n_ops = static_cast<uint32_t>(reads[i].cigar.size());
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
                raw.push_back(cop.op & 0x0F);
                uint32_t len = cop.len;
                do {
                    uint8_t b = len & 0x7F;
                    len >>= 7;
                    if (len) b |= 0x80;
                    raw.push_back(b);
                } while (len);
            }

            // (c) Aux-tag substream: varint length + raw bytes
            uint32_t aux_len = static_cast<uint32_t>(reads[i].aux_data.size());
            {
                uint32_t v = aux_len;
                do {
                    uint8_t b = v & 0x7F;
                    v >>= 7;
                    if (v) b |= 0x80;
                    raw.push_back(b);
                } while (v);
            }
            if (aux_len > 0) {
                raw.insert(raw.end(), reads[i].aux_data.begin(), reads[i].aux_data.end());
            }
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

    // ─── ZSTD_FALLBACK (0x02): whole-record layout ─────────────────────────
    if (codec_id == static_cast<uint8_t>(MetaCodec::ZSTD_FALLBACK)) {
        size_t off = 0;

        // FLAG (4 bytes LE)
        if (off + 4 > raw.size()) {
            throw Error{ErrorCode::STREAM_DECODE_ERROR, "S_meta ZSTD_FALLBACK: truncated FLAG"};
        }
        std::memcpy(&meta.flag, raw.data() + off, 4);
        off += 4;

        // n_cigar_ops (4 bytes LE)
        if (off + 4 > raw.size()) {
            throw Error{ErrorCode::STREAM_DECODE_ERROR, "S_meta ZSTD_FALLBACK: truncated n_ops"};
        }
        uint32_t n_ops;
        std::memcpy(&n_ops, raw.data() + off, 4);
        off += 4;

        // CIGAR ops: op(1 byte) + len(4 bytes) each
        meta.cigar.resize(n_ops);
        for (uint32_t c = 0; c < n_ops; ++c) {
            if (off + 5 > raw.size()) {
                throw Error{ErrorCode::STREAM_DECODE_ERROR, "S_meta ZSTD_FALLBACK: truncated CIGAR"};
            }
            meta.cigar[c].op = raw[off++];
            uint32_t clen;
            std::memcpy(&clen, raw.data() + off, 4);
            meta.cigar[c].len = clen;
            off += 4;
        }

        // Aux tags: aux_len(4 bytes LE) + aux_data
        if (off + 4 > raw.size()) {
            throw Error{ErrorCode::STREAM_DECODE_ERROR, "S_meta ZSTD_FALLBACK: truncated aux_len"};
        }
        uint32_t aux_len;
        std::memcpy(&aux_len, raw.data() + off, 4);
        off += 4;

        if (aux_len > 0 && off + aux_len <= raw.size()) {
            meta.aux_data.assign(raw.begin() + off, raw.begin() + off + aux_len);
        }

        return meta;
    }

    // ─── TYPED_SPLIT (0x01): three typed substreams ───────────────────────
    size_t off = 0;

    // (a) Read FLAG (2 bytes LE — BAM uint16)
    if (off + 2 > raw.size()) {
        throw Error{ErrorCode::STREAM_DECODE_ERROR, "S_meta: truncated FLAG"};
    }
    meta.flag = static_cast<uint32_t>(raw[off]) |
                (static_cast<uint32_t>(raw[off + 1]) << 8);
    off += 2;

    // (b) Read CIGAR: n_ops as varint, then ops
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
        if (off >= raw.size()) break;
        meta.cigar[c].op = raw[off++] & 0x0F;
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

    // (c) Read aux-tag substream: varint length + raw bytes
    uint32_t aux_len = 0;
    {
        uint32_t shift = 0;
        while (off < raw.size()) {
            uint8_t b = raw[off++];
            aux_len |= (uint32_t)(b & 0x7F) << shift;
            if (!(b & 0x80)) break;
            shift += 7;
        }
    }
    if (aux_len > 0 && off + aux_len <= raw.size()) {
        meta.aux_data.assign(raw.begin() + off, raw.begin() + off + aux_len);
    }

    return meta;
}

// ─── S_map encoding (Contract §2.9: DELTA_RANGE) ───────────────────────────
// DELTA_RANGE: per-read blocks with delta-encoded positions.
//
// Contract §2.9 specifies:
//   (a) chrom_id substream — stored per-read for random access (RLE is for
//       the bulk stream view; per-read blocks need absolute chrom_id)
//   (b) pos substream — first read of each chromosome stores absolute pos;
//       subsequent reads store Δpos = pos_i − pos_{i-1}
//
// Per-read block layout (before ZSTD):
//   [flags: uint8_t]      — bit 0: 1=absolute pos, 0=delta pos
//   [chrom_id: uint32_t LE]
//   [pos_or_delta: int64_t LE]  — absolute if flags.bit0=1, signed delta otherwise
//
// RAW codec stores raw (chrom_id:4 + pos:8) uncompressed per-read block.

MapEncodeResult EncodeMapStream(const std::vector<OrderedRead>& reads,
                                 MapCodec codec) {
    MapEncodeResult result;
    result.codec_id = static_cast<uint8_t>(codec);
    result.directory.resize(reads.size());

    // Track previous chromosome and position for delta encoding
    uint32_t prev_chrom_id = UINT32_MAX;
    uint64_t prev_pos = 0;

    for (size_t i = 0; i < reads.size(); ++i) {
        uint64_t offset = result.payload.size();

        if (codec == MapCodec::RAW) {
            // RAW codec: uncompressed absolute (chrom_id:4 + pos:8)
            std::vector<uint8_t> raw;
            WriteU32(raw, reads[i].chrom_id);
            WriteU64(raw, reads[i].pos);
            result.directory[i] = {offset, static_cast<uint32_t>(raw.size())};
            result.payload.insert(result.payload.end(), raw.begin(), raw.end());
        } else {
            // DELTA_RANGE: delta-encoded positions within same chromosome
            std::vector<uint8_t> raw;

            bool is_absolute = (reads[i].chrom_id != prev_chrom_id);
            uint8_t flags = is_absolute ? 0x01 : 0x00;
            raw.push_back(flags);

            // chrom_id (always absolute for per-read random access)
            WriteU32(raw, reads[i].chrom_id);

            if (is_absolute) {
                // First read of chromosome: absolute pos
                WriteU64(raw, reads[i].pos);
            } else {
                // Delta from previous read's pos (signed, since reads are
                // sorted by pos within chromosome so delta is always ≥ 0,
                // but we store as int64 for robustness)
                int64_t delta = static_cast<int64_t>(reads[i].pos) -
                                static_cast<int64_t>(prev_pos);
                uint64_t udelta;
                std::memcpy(&udelta, &delta, sizeof(udelta));
                WriteU64(raw, udelta);
            }

            prev_chrom_id = reads[i].chrom_id;
            prev_pos = reads[i].pos;

            auto compressed = ZstdCompress(raw.data(), raw.size());
            result.directory[i] = {offset, static_cast<uint32_t>(compressed.size())};
            result.payload.insert(result.payload.end(), compressed.begin(), compressed.end());
        }
    }

    return result;
}

DecodedMap DecodeMapRead(const std::vector<uint8_t>& payload,
                         const StreamDirectoryEntry& dir_entry,
                         uint8_t codec_id) {
    DecodedMap map;

    if (codec_id == static_cast<uint8_t>(MapCodec::RAW)) {
        // RAW codec: uncompressed (chrom_id:4 + pos:8)
        const uint8_t* data = payload.data() + dir_entry.offset;
        if (dir_entry.length < 12) {
            throw Error{ErrorCode::STREAM_DECODE_ERROR,
                        "S_map RAW block too short"};
        }
        map.chrom_id = ReadU32(data);
        map.pos = ReadU64(data + 4);
        return map;
    }

    // DELTA_RANGE codec: ZSTD-compressed per-read block
    const uint8_t* data = payload.data() + dir_entry.offset;
    uint64_t expected = ZSTD_getFrameContentSize(data, dir_entry.length);
    if (expected == ZSTD_CONTENTSIZE_UNKNOWN || expected == ZSTD_CONTENTSIZE_ERROR) {
        expected = 13;  // flags(1) + chrom_id(4) + pos(8)
    }
    auto raw = ZstdDecompress(data, dir_entry.length, expected);

    if (raw.size() < 13) {
        throw Error{ErrorCode::STREAM_DECODE_ERROR,
                    "S_map DELTA_RANGE block too short (expected ≥13 bytes, got " +
                    std::to_string(raw.size()) + ")"};
    }

    uint8_t flags = raw[0];
    map.chrom_id = ReadU32(raw.data() + 1);
    uint64_t pos_or_delta = ReadU64(raw.data() + 5);

    if (flags & 0x01) {
        // Absolute pos
        map.pos = pos_or_delta;
    } else {
        // Delta pos — store as delta value; caller must accumulate
        // For per-read random access, the caller should reconstruct absolute
        // pos by walking from the nearest absolute entry. However, since each
        // read is independently decodable, we store a flag indicating delta
        // mode and the accumulated pos must be provided by the caller.
        //
        // In practice, for sequential decode (reconstruction), the caller
        // maintains a running pos. For random access, the caller scans from
        // the last absolute entry (chromosome boundary).
        //
        // NOTE: To maintain per-read independence (I10), we store the delta
        // and require the decode context to track state. This is the correct
        // trade-off per Contract §2.9.
        int64_t delta;
        std::memcpy(&delta, &pos_or_delta, sizeof(delta));
        map.pos = static_cast<uint64_t>(delta);
        map.is_delta = true;
    }
    return map;
}

}  // namespace bamsix

