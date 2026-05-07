#pragma once

#include "bamsi/types.hpp"
#include "../fmindex/fmindex.hpp"
#include "../bitvectors/bitvectors.hpp"
#include "../mapping/mapping.hpp"
#include <functional>
#include <map>
#include <vector>

namespace bamsi {

/// Compute the reverse complement of a pattern (codes 0..4).
std::vector<uint8_t> ReverseComplement(const std::vector<uint8_t>& pattern);

/// Determine Q(P): the set of orientations to search.
/// If P == rc(P), returns {Forward}; else returns {Forward, Reverse}.
struct QueryOrientation {
    std::vector<uint8_t> pattern;
    QueryStrand          strand;
};
std::vector<QueryOrientation> GetOrientations(const std::vector<uint8_t>& P, StrandMode mode);

/// Validate a pattern: must be non-empty and contain only codes 0..4.
ErrorCode ValidatePattern(const std::vector<uint8_t>& pattern);

/// GlobalCount(P) — total occurrences across all orientations.
uint64_t GlobalCount(const std::vector<uint8_t>& P,
                     const FMIndexEngine& fm,
                     StrandMode mode);

/// GlobalExists(P) — true if P occurs at least once.
bool GlobalExists(const std::vector<uint8_t>& P,
                  const FMIndexEngine& fm,
                  StrandMode mode);

/// Locate(P) — return all matches.
/// @param sorted  If true, sort results per §6.2 ordering.
std::vector<Match> Locate(const std::vector<uint8_t>& P,
                          const FMIndexEngine& fm,
                          const SuccinctBitvector& B_read,
                          const std::vector<OrderedRead>& reads,
                          const std::vector<std::string>& chrom_names,
                          StrandMode mode,
                          bool sorted = false);

/// RegionalCount(P, chrom, [a,b]) — BASE tier.
uint64_t RegionalCount(const std::vector<uint8_t>& P,
                       const std::string& chrom,
                       uint64_t a, uint64_t b,
                       const FMIndexEngine& fm,
                       const SuccinctBitvector& B_read,
                       const SuccinctBitvector& B_window,
                       const WindowTable& windows,
                       const std::vector<OrderedRead>& reads,
                       const std::map<std::string, uint32_t>& chrom_to_id,
                       StrandMode mode);

/// RegionalExists(P, T_threshold, chrom, [a,b]) — BASE tier.
bool RegionalExists(const std::vector<uint8_t>& P,
                    uint64_t T_threshold,
                    const std::string& chrom,
                    uint64_t a, uint64_t b,
                    const FMIndexEngine& fm,
                    const SuccinctBitvector& B_read,
                    const SuccinctBitvector& B_window,
                    const WindowTable& windows,
                    const std::vector<OrderedRead>& reads,
                    const std::map<std::string, uint32_t>& chrom_to_id,
                    StrandMode mode);

}  // namespace bamsi
