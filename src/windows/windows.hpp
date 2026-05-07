#pragma once

#include "bamsi/types.hpp"
#include <vector>

namespace bamsi {

/// Compute the reference span of a CIGAR string.
/// ref_span(read) = Σ len(op) for op ∈ {M, =, X, D, N}
uint64_t CigarRefSpan(const CigarRecord& cigar);

/// Build the WindowTable per Architecture §4.8.
/// Windows partition S into non-overlapping regions of ~T S-characters each,
/// respecting chromosome boundaries and the no-split-read rule.
WindowTable BuildWindows(const std::vector<OrderedRead>& reads,
                         const SequenceBundle& bundle,
                         uint64_t T);

}  // namespace bamsi
