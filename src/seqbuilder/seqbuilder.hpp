#pragma once

#include "bamsi/types.hpp"
#include <vector>

namespace bamsi {

/// Build the concatenated sequence S = r0 # r1 # ... # r_{N-1}
/// and the readStarts array.
/// Input: ordered reads (must already be sorted and have read_id assigned).
/// Output: SequenceBundle with S and readStarts.
SequenceBundle BuildSequence(const std::vector<OrderedRead>& reads);

}  // namespace bamsi
