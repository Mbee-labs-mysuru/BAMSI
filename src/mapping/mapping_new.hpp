#pragma once

#include "bamsix/types.hpp"
#include "../bitvectors/bitvectors.hpp"
#include "../streamencode/streamencode.hpp"
#include <vector>
#include <cstdint>

namespace bamsix {

/// Direction for CIGAR mapping.
enum class CigarDirection { LEFT, RIGHT };

/// Map a read-offset to a reference position using CIGAR operations.
/// Total and deterministic per Architecture §5.2.
uint64_t CigarRefPos(const CigarRecord& cigar, uint64_t p_anchor,
                     uint64_t offset, CigarDirection direction);

/// Absolute position table entry: one per chromosome boundary (absolute entry).
/// Built once at query load time from the S_map stream via BuildAbsolutePosTable().
struct AbsPosEntry {
    uint64_t read_id;   ///< Index of this absolute-position read in the stream.
    uint64_t abs_pos;   ///< Absolute genomic position of this read.
};

/// Build the sparse absolute-position lookup table from the S_map payload.
/// Called ONCE at index load time. O(N) to scan all read blocks.
/// Stores one AbsPosEntry per chromosome boundary (is_delta==false) read.
/// Result is sorted by read_id for binary search during MapOccurrenceLazy.
std::vector<AbsPosEntry> BuildAbsolutePosTable(
    const std::vector<uint8_t>& map_payload,
    const StreamDirectoryPerRead& dir_map,
    uint8_t map_codec_id);

/// Dense absolute-position table: stores absolute position for EVERY read.
/// Built once at index load time via BuildFullAbsPosTable(). O(N) build,
/// O(1) lookup per occurrence — eliminates the O(K) forward walk that made
/// locate hang for hours on high-occurrence patterns.
/// Memory: N × 8 bytes (16 MB for 2M reads — negligible).
struct FullAbsPosTable {
    std::vector<uint64_t> positions;  ///< positions[read_id] = absolute genomic pos
};

/// Build the dense absolute-position table from the S_map payload.
/// Single O(N) sequential pass: decodes every map block, accumulates deltas,
/// stores the resulting absolute position for each read.
FullAbsPosTable BuildFullAbsPosTable(
    const std::vector<uint8_t>& map_payload,
    const StreamDirectoryPerRead& dir_map,
    uint8_t map_codec_id);

/// Mapping function M_ℓ per Architecture §5.1.
/// Maps an S-position + pattern length to a genomic interval.
/// Uses bitvector operations only for separator detection (I14: no raw S access).
///
/// LEGACY: uses in-memory OrderedRead vector. For small datasets / tests.
MappingResult MapOccurrence(uint64_t pos, uint64_t pattern_len,
                            QueryStrand strand,
                            const SuccinctBitvector& B_read,
                            const std::vector<OrderedRead>& reads);

/// Lazy per-read mapping (Contract §2.3 Table B / §3.6.2 Steps 3-4).
/// Decodes ONLY the needed read's metadata and mapping data from compressed
/// S_meta/S_map streams using per-read directories.
///
/// Delta resolution uses abs_table (pre-built by BuildAbsolutePosTable) for
/// O(log C + K) cost per occurrence, where C = chromosome count and K = reads
/// between consecutive chromosome boundaries (typically a few thousand).
/// Without abs_table, falls back to O(N) scan (legacy behaviour).
///
/// This is the production path for billion-read scalability.
MappingResult MapOccurrenceLazy(uint64_t pos, uint64_t pattern_len,
                                QueryStrand strand,
                                const SuccinctBitvector& B_read,
                                const std::vector<uint8_t>& meta_payload,
                                const std::vector<uint8_t>& map_payload,
                                const StreamDirectoryPerRead& dir_meta,
                                const StreamDirectoryPerRead& dir_map,
                                uint8_t meta_codec_id,
                                uint8_t map_codec_id,
                                const FullAbsPosTable* full_abs_table = nullptr);

/// Check if pos is a separator position.
/// Uses bitvector operations only (I14).
bool IsSeparatorPosition(uint64_t pos, const SuccinctBitvector& B_read);

}  // namespace bamsix
