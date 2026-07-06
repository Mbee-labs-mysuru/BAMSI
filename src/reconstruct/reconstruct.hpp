#pragma once

#include "bamsix/types.hpp"
#include "fmindex/fmindex.hpp"
#include "bitvectors/bitvectors.hpp"

#include <vector>
#include <cstdint>

namespace bamsix {

/// C3 fix: Extract a single read's sequence from the BWT via LF-walk.
///
/// Contract §2.3 Table A: "Given only S_seq (which stores the BWT) and B_read,
/// reconstruct read sequences r_i for all i."
///
/// Algorithm:
///   1. Use B_read.Select1(read_id + 1) to find readStarts[read_id] in S
///   2. The read spans [readStarts[read_id], readStarts[read_id] + len - 1]
///   3. To extract S[pos], we need to find the BWT row whose SA value = pos,
///      then the character at that position is determined by the FM-index
///   4. Walk from readStarts[read_id] through consecutive S-positions,
///      using the ISA (inverse suffix array) to map text positions to BWT rows
///
/// Since we don't store the full ISA, we use the LF-walk approach:
///   - Start from row 0 (which corresponds to SA[0], the suffix starting at
///     the lexicographically smallest position)
///   - Walk the entire BWT to reconstruct S in text order
///   - For efficiency, we do a single pass to reconstruct all reads at once
///
/// For per-read extraction, we use a different approach:
///   - Use SA samples to find a nearby sampled row
///   - LF-walk from there to the target position
///   - Read characters at consecutive positions
std::vector<uint8_t> ExtractReadSequence(
    const FMIndexEngine& fm,
    const SuccinctBitvector& B_read,
    uint64_t read_id,
    uint64_t read_len);

/// Batch extraction: reconstruct all read sequences from BWT.
/// Returns a vector of vectors, one per read.
///
/// This performs a single reverse LF-walk over the entire BWT to reconstruct
/// the text S, then slices it into per-read sequences using readStarts from B_read.
std::vector<std::vector<uint8_t>> ExtractAllSequences(
    const FMIndexEngine& fm,
    const SuccinctBitvector& B_read,
    uint64_t N_reads,
    uint64_t S_length);

}  // namespace bamsix
