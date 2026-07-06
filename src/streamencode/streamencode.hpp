#pragma once

#include "bamsix/types.hpp"
#include <vector>

namespace bamsix {

/// Encode S_seq stream (sequence data).
/// Contract §2.4 mandatory pipeline: BWT → MTF → RLE → 0th-order Arithmetic.
/// Codec ID 0x10 = BWT_MTF_RLE_ARITH (active default).
/// Legacy codec ID 0x01 = ZSTD (decode-only, for backward compat with early builds).
struct SeqEncodeResult {
    std::vector<uint8_t> payload;
    uint8_t              codec_id;
};
SeqEncodeResult EncodeSeqStream(const std::vector<uint8_t>& bwt,
                                uint8_t entropy_order_k);

/// Decode S_seq stream.
std::vector<uint8_t> DecodeSeqStream(const std::vector<uint8_t>& payload,
                                     uint8_t codec_id, uint64_t expected_len);

/// Encode S_qual stream (quality scores).
/// Contract §2.3 F4: dir_qual may be per-read (block_size=0) or block-level.
struct QualEncodeResult {
    std::vector<uint8_t>    payload;
    StreamDirectory         directory;   // per-read OR block-level (F4)
    uint8_t                 codec_id;
    uint32_t                block_size = 0;  // 0 = per-read
};
QualEncodeResult EncodeQualStream(const std::vector<OrderedRead>& reads,
                                  QualCodec codec, uint8_t lossy_bins,
                                  uint32_t block_size = 0);

/// Decode S_qual stream (quality scores for a single read, per-read directory).
std::vector<uint8_t> DecodeQualRead(const std::vector<uint8_t>& payload,
                                    const StreamDirectoryEntry& dir_entry,
                                    uint8_t codec_id);

/// Decode S_qual from a block-level directory (Contract §2.3 F4).
/// Returns quality scores for read `read_id` by decoding the containing block.
std::vector<uint8_t> DecodeQualFromBlock(const std::vector<uint8_t>& payload,
                                          const StreamDirectoryBlockLevel& block_dir,
                                          uint32_t block_size,
                                          uint64_t read_id,
                                          uint8_t codec_id);

/// Encode S_meta stream (CIGAR + FLAG + aux tags).
struct MetaEncodeResult {
    std::vector<uint8_t>    payload;
    StreamDirectoryPerRead  directory;
    uint8_t                 codec_id;
};
MetaEncodeResult EncodeMetaStream(const std::vector<OrderedRead>& reads,
                                   MetaCodec codec);

/// Decode S_meta stream (CIGAR + FLAG + aux tags for a single read).
struct DecodedMeta {
    CigarRecord              cigar;
    uint32_t                 flag;
    std::vector<uint8_t>     aux_data;    // Contract §2.8: optional BAM aux tags
};
DecodedMeta DecodeMetaRead(const std::vector<uint8_t>& payload,
                           const StreamDirectoryEntry& dir_entry,
                           uint8_t codec_id);

/// Encode S_map stream (chrom_id + pos).
struct MapEncodeResult {
    std::vector<uint8_t>    payload;
    StreamDirectoryPerRead  directory;
    uint8_t                 codec_id;
};
MapEncodeResult EncodeMapStream(const std::vector<OrderedRead>& reads,
                                 MapCodec codec);

/// Decode S_map stream (chrom_id + pos for a single read).
struct DecodedMap {
    uint32_t chrom_id;
    uint64_t pos;
    bool     is_delta = false;  // F2: true when pos is a delta (DELTA_RANGE codec)
};
DecodedMap DecodeMapRead(const std::vector<uint8_t>& payload,
                         const StreamDirectoryEntry& dir_entry,
                         uint8_t codec_id);

}  // namespace bamsix
