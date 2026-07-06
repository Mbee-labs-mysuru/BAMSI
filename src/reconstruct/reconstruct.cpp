#include "reconstruct.hpp"

#include <algorithm>
#include <cstring>

namespace bamsix {

// ─── Full text reconstruction via LF-walk (Contract §2.3 Table A) ───────────
//
// The BWT allows us to reconstruct the original text S by reverse LF-walking
// from the first row (which corresponds to the sentinel character '$' at the
// end of S). Starting from the sentinel row:
//   1. The sentinel row's BWT character is S[|S|-1] (last char of S)
//   2. LF(sentinel_row) gives us the row whose suffix starts at position |S|-1
//   3. Continuing LF-walks reconstructs S in reverse order
//
// We reconstruct S[|S|-1], S[|S|-2], ..., S[0] by reading BWT characters
// during the walk, then reverse to get S[0...|S|-1].

std::vector<std::vector<uint8_t>> ExtractAllSequences(
    const FMIndexEngine& fm,
    const SuccinctBitvector& B_read,
    uint64_t N_reads,
    uint64_t S_length) {

    if (S_length == 0 || N_reads == 0) {
        return std::vector<std::vector<uint8_t>>(N_reads);
    }

    // Step 1: Reconstruct the full text S by reverse LF-walk from sentinel row.
    //
    // The FM-index represents the text S$ where $ is the sentinel.
    // Starting from sentinel_row (where SA[sentinel_row] = |S|, conceptual
    // position of $), each LF step walks backward through S:
    //   row_0 = sentinel_row      → SA[row_0] = |S|  (sentinel position)
    //   row_1 = LF(row_0)         → SA[row_1] = |S|-1, BWT[row_0] = S[|S|-1]
    //   row_2 = LF(row_1)         → SA[row_2] = |S|-2, BWT[row_1] = S[|S|-2]
    //   ...
    //   row_k = LF(row_{k-1})     → BWT[row_{k-1}] = S[|S|-k]
    //
    // After |S| steps we have S in reverse order.

    std::vector<uint8_t> S_reversed(S_length);
    uint64_t row = fm.SentinelRow();

    for (uint64_t k = 0; k < S_length; ++k) {
        // BWT[row] = the character at text position |S|-1-k
        uint8_t c = fm.BwtAt(row);
        S_reversed[k] = c;
        row = fm.LF(row);
    }

    // Step 2: Reverse to get S[0...|S|-1]
    std::vector<uint8_t> S(S_length);
    for (uint64_t i = 0; i < S_length; ++i) {
        S[i] = S_reversed[S_length - 1 - i];
    }

    // Step 3: Slice S into per-read sequences using B_read.
    // readStarts[i] = Select1(B_read, i+1)
    // Read i spans S[readStarts[i]..readStarts[i]+len-1]
    // where len = readStarts[i+1] - readStarts[i] - 1 (subtracting the '#' separator)
    // For the last read (i = N-1), len = |S| - readStarts[N-1]

    std::vector<std::vector<uint8_t>> result(N_reads);

    for (uint64_t i = 0; i < N_reads; ++i) {
        uint64_t start = B_read.Select1(i + 1);  // readStarts[i]
        uint64_t end;
        if (i + 1 < N_reads) {
            uint64_t next_start = B_read.Select1(i + 2);
            // The '#' separator is at next_start - 1
            end = next_start - 1;  // exclusive of separator
        } else {
            end = S_length;  // last read goes to end of S
        }

        uint64_t len = end - start;
        result[i].resize(len);
        for (uint64_t j = 0; j < len; ++j) {
            result[i][j] = S[start + j];
        }
    }

    return result;
}

std::vector<uint8_t> ExtractReadSequence(
    const FMIndexEngine& fm,
    const SuccinctBitvector& B_read,
    uint64_t read_id,
    uint64_t read_len) {

    // For single-read extraction, we could do a full reconstruction and slice,
    // but that's wasteful. Instead, use the Locate mechanism in reverse:
    //
    // We know readStarts[read_id] = Select1(B_read, read_id + 1).
    // We need S[readStarts[read_id] ... readStarts[read_id] + read_len - 1].
    //
    // To get S[pos], we need the BWT row r such that SA[r] = pos.
    // This is ISA[pos] = r. If we have ISA samples, we can use them.
    // If not, we must do a full reconstruction.
    //
    // For v1.0, we use the full reconstruction approach for correctness
    // (the per-read extraction via ISA is an optimization for v1.1).
    //
    // However, for a single read, we can still use the LF-walk approach
    // by walking the full text and extracting only the needed range.
    // This is O(|S|) which is expensive but correct.

    uint64_t S_length = fm.SLen();
    uint64_t N_reads_approx = B_read.PopCount();  // number of 1-bits = number of reads

    // Full reconstruction then slice
    auto all = ExtractAllSequences(fm, B_read, N_reads_approx, S_length);
    if (read_id < all.size()) {
        return all[read_id];
    }
    return {};
}

}  // namespace bamsix
