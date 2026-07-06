#pragma once

#include "bamsix/types.hpp"
#include "../fmindex/fmindex.hpp"
#include "../bitvectors/bitvectors.hpp"
#include "../mapping/mapping.hpp"
#include "../sarange/sarange.hpp"
#include <functional>
#include <map>
#include <vector>

namespace bamsix {

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
/// When overlap == NonOverlapping, resolves SA positions and applies
/// a greedy left-to-right filter (grep -oE semantics). Complexity
/// increases from O(|P|) to O(|P| + occ·s), or O(|P| + occ) with ISA.
uint64_t GlobalCount(const std::vector<uint8_t>& P,
                     const FMIndexEngine& fm,
                     StrandMode mode,
                     OverlapMode overlap = OverlapMode::All);

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

// ─── Lazy query overloads (Contract §2.3 Table B) ──────────────────────────
// These accept compressed stream payloads + directories instead of in-memory
// reads. They decode per-read metadata lazily via MapOccurrenceLazy.
// This is the PRODUCTION path for billion-read scalability.

/// Locate(P) — lazy decode from streams.
/// Pass full_abs_table (built by BuildFullAbsPosTable) to enable O(1)
/// delta resolution; omit or pass nullptr for legacy O(N) fallback.
std::vector<Match> LocateLazy(const std::vector<uint8_t>& P,
                              const FMIndexEngine& fm,
                              const SuccinctBitvector& B_read,
                              const std::vector<uint8_t>& meta_payload,
                              const std::vector<uint8_t>& map_payload,
                              const StreamDirectoryPerRead& dir_meta,
                              const StreamDirectoryPerRead& dir_map,
                              uint8_t meta_codec_id, uint8_t map_codec_id,
                              const std::vector<std::string>& chrom_names,
                              StrandMode mode,
                              bool sorted = false,
                              const FullAbsPosTable* full_abs_table = nullptr);

/// RegionalCount — lazy decode from streams.
uint64_t RegionalCountLazy(const std::vector<uint8_t>& P,
                           const std::string& chrom,
                           uint64_t a, uint64_t b,
                           const FMIndexEngine& fm,
                           const SuccinctBitvector& B_read,
                           const SuccinctBitvector& B_window,
                           const WindowTable& windows,
                           const std::vector<uint8_t>& meta_payload,
                           const std::vector<uint8_t>& map_payload,
                           const StreamDirectoryPerRead& dir_meta,
                           const StreamDirectoryPerRead& dir_map,
                           uint8_t meta_codec_id, uint8_t map_codec_id,
                           const std::map<std::string, uint32_t>& chrom_to_id,
                           StrandMode mode,
                           const FullAbsPosTable* full_abs_table = nullptr);

/// RegionalExists — lazy decode from streams.
bool RegionalExistsLazy(const std::vector<uint8_t>& P,
                        uint64_t T_threshold,
                        const std::string& chrom,
                        uint64_t a, uint64_t b,
                        const FMIndexEngine& fm,
                        const SuccinctBitvector& B_read,
                        const SuccinctBitvector& B_window,
                        const WindowTable& windows,
                        const std::vector<uint8_t>& meta_payload,
                        const std::vector<uint8_t>& map_payload,
                        const StreamDirectoryPerRead& dir_meta,
                        const StreamDirectoryPerRead& dir_map,
                        uint8_t meta_codec_id, uint8_t map_codec_id,
                        const std::map<std::string, uint32_t>& chrom_to_id,
                        StrandMode mode,
                        const FullAbsPosTable* full_abs_table = nullptr);

// ─── ENHANCED tier (SARange) overloads (Architecture §5.3) ─────────────────
// These use the SARange wavelet tree for O(|W_r|·log(|S|/s)) RegionalCount
// instead of O(occ·s) BASE-tier locate-and-filter.

/// RegionalCount — ENHANCED tier using SARange wavelet tree.
uint64_t RegionalCountSARange(const std::vector<uint8_t>& P,
                               const std::string& chrom,
                               uint64_t a, uint64_t b,
                               const FMIndexEngine& fm,
                               const SuccinctBitvector& B_read,
                               const SuccinctBitvector& B_window,
                               const WindowTable& windows,
                               const std::vector<OrderedRead>& reads,
                               const SARange& sarange,
                               const std::map<std::string, uint32_t>& chrom_to_id,
                               StrandMode mode);

/// RegionalExists — ENHANCED tier using SARange wavelet tree.
bool RegionalExistsSARange(const std::vector<uint8_t>& P,
                            uint64_t T_threshold,
                            const std::string& chrom,
                            uint64_t a, uint64_t b,
                            const FMIndexEngine& fm,
                            const SuccinctBitvector& B_read,
                            const SuccinctBitvector& B_window,
                            const WindowTable& windows,
                            const std::vector<OrderedRead>& reads,
                            const SARange& sarange,
                            const std::map<std::string, uint32_t>& chrom_to_id,
                            StrandMode mode);

}  // namespace bamsix

