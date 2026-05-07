#pragma once

#include "bamsi/types.hpp"
#include "../bitvectors/bitvectors.hpp"
#include <vector>

namespace bamsi {

/// Direction for CIGAR mapping.
enum class CigarDirection { LEFT, RIGHT };

/// Map a read-offset to a reference position using CIGAR operations.
/// Total and deterministic per Architecture §5.2.
uint64_t CigarRefPos(const CigarRecord& cigar, uint64_t p_anchor,
                     uint64_t offset, CigarDirection direction);

/// Mapping function M_ℓ per Architecture §5.1.
/// Maps an S-position + pattern length to a genomic interval.
/// Uses bitvector operations only for separator detection (I14: no raw S access).
MappingResult MapOccurrence(uint64_t pos, uint64_t pattern_len,
                            QueryStrand strand,
                            const SuccinctBitvector& B_read,
                            const std::vector<OrderedRead>& reads);

/// Check if pos is a separator position.
/// Uses bitvector operations only (I14).
bool IsSeparatorPosition(uint64_t pos, const SuccinctBitvector& B_read);

}  // namespace bamsi
