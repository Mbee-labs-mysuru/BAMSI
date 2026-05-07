#pragma once

#include "bamsi/types.hpp"
#include <vector>

namespace bamsi {

/// Encode S_seq stream (sequence data).
/// V1: ZSTD compression of BWT bytes.
/// Full: BWT → MTF → RLE → Arithmetic with entropy order k.
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
struct QualEncodeResult {
    std::vector<uint8_t>    payload;
    StreamDirectoryPerRead  directory;
    uint8_t                 codec_id;
};
QualEncodeResult EncodeQualStream(const std::vector<OrderedRead>& reads,
                                  QualCodec codec, uint8_t lossy_bins);

/// Decode S_qual stream (quality scores for a single read).
std::vector<uint8_t> DecodeQualRead(const std::vector<uint8_t>& payload,
                                    const StreamDirectoryEntry& dir_entry,
                                    uint8_t codec_id);

/// Encode S_meta stream (CIGAR + FLAG + aux tags).
struct MetaEncodeResult {
    std::vector<uint8_t>    payload;
    StreamDirectoryPerRead  directory;
    uint8_t                 codec_id;
};
MetaEncodeResult EncodeMetaStream(const std::vector<OrderedRead>& reads,
                                   MetaCodec codec);

/// Decode S_meta stream (CIGAR + FLAG for a single read).
struct DecodedMeta {
    CigarRecord cigar;
    uint32_t    flag;
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
};
DecodedMap DecodeMapRead(const std::vector<uint8_t>& payload,
                         const StreamDirectoryEntry& dir_entry,
                         uint8_t codec_id);

}  // namespace bamsi
