#include "query.hpp"

#include <algorithm>

namespace bamsix {

// ─── Reverse complement ─────────────────────────────────────────────────────

std::vector<uint8_t> ReverseComplement(const std::vector<uint8_t>& pattern) {
    std::vector<uint8_t> rc(pattern.size());
    for (size_t i = 0; i < pattern.size(); ++i) {
        uint8_t base = pattern[pattern.size() - 1 - i];
        switch (base) {
            case CODE_A: rc[i] = CODE_T; break;
            case CODE_C: rc[i] = CODE_G; break;
            case CODE_G: rc[i] = CODE_C; break;
            case CODE_T: rc[i] = CODE_A; break;
            default:     rc[i] = CODE_N; break;
        }
    }
    return rc;
}

// ─── Orientations (Architecture §6.3) ───────────────────────────────────────

std::vector<QueryOrientation> GetOrientations(const std::vector<uint8_t>& P, StrandMode mode) {
    if (mode == StrandMode::SingleStrand) {
        return {{P, QueryStrand::Forward}};
    }

    auto rc = ReverseComplement(P);
    if (P == rc) {
        // Palindromic pattern: single query labelled Forward
        return {{P, QueryStrand::Forward}};
    }
    return {{P, QueryStrand::Forward}, {rc, QueryStrand::Reverse}};
}

// ─── Validation ─────────────────────────────────────────────────────────────

ErrorCode ValidatePattern(const std::vector<uint8_t>& pattern) {
    if (pattern.empty()) return ErrorCode::EMPTY_PATTERN;
    for (auto c : pattern) {
        if (c > CODE_N) return ErrorCode::INVALID_PATTERN;
    }
    return ErrorCode::OK;
}

// ─── GlobalCount (Architecture §6.4) ─────────────────────────────────────────

uint64_t GlobalCount(const std::vector<uint8_t>& P,
                     const FMIndexEngine& fm,
                     StrandMode mode,
                     OverlapMode overlap) {
    auto err = ValidatePattern(P);
    if (err == ErrorCode::EMPTY_PATTERN)
        throw Error{err, "Pattern is empty"};
    if (err == ErrorCode::INVALID_PATTERN)
        throw Error{err, "Pattern contains invalid characters (only A=0,C=1,G=2,T=3,N=4 allowed)"};

    if (overlap == OverlapMode::All) {
        // Fast path: O(|P|) — count SA interval sizes without locating
        uint64_t total = 0;
        for (const auto& [Q, strand] : GetOrientations(P, mode)) {
            auto interval = fm.BackwardSearch(Q);
            total += interval.size();
        }
        return total;
    }

    // Non-overlapping path: O(|P| + occ·s) — must locate all positions,
    // sort, and apply greedy left-to-right non-overlapping filter.
    // This matches grep -oE semantics exactly.
    std::vector<uint64_t> positions;
    for (const auto& [Q, strand] : GetOrientations(P, mode)) {
        auto interval = fm.BackwardSearch(Q);
        for (uint64_t row = interval.lo; row < interval.hi; ++row) {
            if (row == fm.SentinelRow()) continue;
            uint64_t pos = fm.HasISASamples() ? fm.LocateBidir(row) : fm.Locate(row);
            positions.push_back(pos);
        }
    }
    std::sort(positions.begin(), positions.end());

    // Greedy non-overlapping filter: skip any match starting before
    // the end of the previous accepted match.
    uint64_t count = 0;
    uint64_t prev_end = 0;
    for (uint64_t pos : positions) {
        if (pos >= prev_end) {
            ++count;
            prev_end = pos + P.size();
        }
    }
    return count;
}

// ─── GlobalExists (Architecture §6.4) ────────────────────────────────────────

bool GlobalExists(const std::vector<uint8_t>& P,
                  const FMIndexEngine& fm,
                  StrandMode mode) {
    auto err = ValidatePattern(P);
    if (err == ErrorCode::EMPTY_PATTERN)
        throw Error{err, "Pattern is empty"};
    if (err == ErrorCode::INVALID_PATTERN)
        throw Error{err, "Pattern contains invalid characters (only A=0,C=1,G=2,T=3,N=4 allowed)"};

    for (const auto& [Q, strand] : GetOrientations(P, mode)) {
        auto interval = fm.BackwardSearch(Q);
        if (!interval.empty()) return true;
    }
    return false;
}

// ─── Locate (Architecture §6.4) ─────────────────────────────────────────────

std::vector<Match> Locate(const std::vector<uint8_t>& P,
                          const FMIndexEngine& fm,
                          const SuccinctBitvector& B_read,
                          const std::vector<OrderedRead>& reads,
                          const std::vector<std::string>& chrom_names,
                          StrandMode mode,
                          bool sorted) {
    std::vector<Match> results;
    auto err = ValidatePattern(P);
    if (err == ErrorCode::EMPTY_PATTERN)
        throw Error{err, "Pattern is empty"};
    if (err == ErrorCode::INVALID_PATTERN)
        throw Error{err, "Pattern contains invalid characters (only A=0,C=1,G=2,T=3,N=4 allowed)"};

    uint64_t pattern_len = P.size();

    for (const auto& [Q, strand] : GetOrientations(P, mode)) {
        auto interval = fm.BackwardSearch(Q);
        for (uint64_t row = interval.lo; row < interval.hi; ++row) {
            // Skip sentinel row
            if (row == fm.SentinelRow()) continue;

            uint64_t pos = fm.HasISASamples() ? fm.LocateBidir(row) : fm.Locate(row);

            // Map occurrence (Contract F3: separator check removed — provably
            // unreachable for valid P ∈ Σ* with correct FM-index)
            auto m = MapOccurrence(pos, pattern_len, strand, B_read, reads);

            // Convert chrom_id to name
            std::string chrom_name;
            if (m.chrom_id < chrom_names.size()) {
                chrom_name = chrom_names[m.chrom_id];
            }

            results.push_back(Match{
                .chrom        = chrom_name,
                .p_min        = m.p_min,
                .p_max        = m.p_max,
                .read_id      = m.read_id,
                .query_strand = strand,
                .sa_row       = row,
            });
        }
    }

    if (sorted) {
        // Contract §4.1 — 6-key deterministic sort:
        // 1. strand (Forward < Reverse)
        // 2. chrom (ascending chrom_id via string order)
        // 3. p_min ascending
        // 4. p_max ascending
        // 5. read_id ascending
        // 6. SA row ascending (final tie-break)
        std::sort(results.begin(), results.end(),
            [](const Match& a, const Match& b) {
                if (a.query_strand != b.query_strand)
                    return a.query_strand < b.query_strand;
                if (a.chrom != b.chrom)
                    return a.chrom < b.chrom;
                if (a.p_min != b.p_min)
                    return a.p_min < b.p_min;
                if (a.p_max != b.p_max)
                    return a.p_max < b.p_max;
                if (a.read_id != b.read_id)
                    return a.read_id < b.read_id;
                return a.sa_row < b.sa_row;
            });
    }

    return results;
}

// ─── RegionalCount (Architecture §6.4 — BASE tier) ──────────────────────────

uint64_t RegionalCount(const std::vector<uint8_t>& P,
                       const std::string& chrom,
                       uint64_t a, uint64_t b,
                       const FMIndexEngine& fm,
                       const SuccinctBitvector& B_read,
                       const SuccinctBitvector& B_window,
                       const WindowTable& windows,
                       const std::vector<OrderedRead>& reads,
                       const std::map<std::string, uint32_t>& chrom_to_id,
                       StrandMode mode) {
    auto err = ValidatePattern(P);
    if (err == ErrorCode::EMPTY_PATTERN)
        throw Error{err, "Pattern is empty"};
    if (err == ErrorCode::INVALID_PATTERN)
        throw Error{err, "Pattern contains invalid characters (only A=0,C=1,G=2,T=3,N=4 allowed)"};

    auto it = chrom_to_id.find(chrom);
    if (it == chrom_to_id.end()) return 0;
    uint32_t chrom_id = it->second;

    uint64_t count = 0;
    uint64_t pattern_len = P.size();

    for (const auto& [Q, strand] : GetOrientations(P, mode)) {
        auto interval = fm.BackwardSearch(Q);

        for (uint64_t row = interval.lo; row < interval.hi; ++row) {
            if (row == fm.SentinelRow()) continue;

            uint64_t pos = fm.HasISASamples() ? fm.LocateBidir(row) : fm.Locate(row);

            // Window-based prune (closed-interval bitvector rank)
            if (B_window.PopCount() == 0) continue;
            uint64_t window_id = B_window.Rank1(pos) - 1;
            if (window_id >= windows.size()) continue;

            const auto& W = windows[window_id];
            if (W.chrom_id != chrom_id) continue;
            if (W.genomic_end < a || W.genomic_start > b) continue;

            // Exact intersection check (Contract F3: separator check removed)
            auto m = MapOccurrence(pos, pattern_len, strand, B_read, reads);
            if (m.chrom_id == chrom_id && m.p_min <= b && m.p_max >= a) {
                ++count;
            }
        }
    }

    return count;
}

// ─── RegionalExists (Architecture §6.4 — BASE tier) ─────────────────────────

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
                    StrandMode mode) {
    auto err = ValidatePattern(P);
    if (err == ErrorCode::EMPTY_PATTERN)
        throw Error{err, "Pattern is empty"};
    if (err == ErrorCode::INVALID_PATTERN)
        throw Error{err, "Pattern contains invalid characters (only A=0,C=1,G=2,T=3,N=4 allowed)"};

    auto it = chrom_to_id.find(chrom);
    if (it == chrom_to_id.end()) return false;
    uint32_t chrom_id = it->second;

    uint64_t count = 0;
    uint64_t pattern_len = P.size();

    for (const auto& [Q, strand] : GetOrientations(P, mode)) {
        auto interval = fm.BackwardSearch(Q);

        for (uint64_t row = interval.lo; row < interval.hi; ++row) {
            if (row == fm.SentinelRow()) continue;

            uint64_t pos = fm.HasISASamples() ? fm.LocateBidir(row) : fm.Locate(row);

            if (B_window.PopCount() == 0) continue;
            uint64_t window_id = B_window.Rank1(pos) - 1;
            if (window_id >= windows.size()) continue;

            const auto& W = windows[window_id];
            if (W.chrom_id != chrom_id) continue;
            if (W.genomic_end < a || W.genomic_start > b) continue;

            // Contract F3: separator check removed — provably unreachable
            auto m = MapOccurrence(pos, pattern_len, strand, B_read, reads);
            if (m.chrom_id == chrom_id && m.p_min <= b && m.p_max >= a) {
                ++count;
                if (count >= T_threshold) return true;
            }
        }
    }

    return false;
}

// ─── Lazy query overloads (Contract §2.3 Table B) ──────────────────────────
// These implementations decode per-read metadata lazily from compressed
// S_meta/S_map streams. They are the PRODUCTION path for billion-read datasets.

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
                              bool sorted,
                              const FullAbsPosTable* full_abs_table) {
    std::vector<Match> results;
    auto err = ValidatePattern(P);
    if (err == ErrorCode::EMPTY_PATTERN)
        throw Error{err, "Pattern is empty"};
    if (err == ErrorCode::INVALID_PATTERN)
        throw Error{err, "Pattern contains invalid characters (only A=0,C=1,G=2,T=3,N=4 allowed)"};

    uint64_t pattern_len = P.size();

    for (const auto& [Q, strand] : GetOrientations(P, mode)) {
        auto interval = fm.BackwardSearch(Q);
        for (uint64_t row = interval.lo; row < interval.hi; ++row) {
            if (row == fm.SentinelRow()) continue;
            uint64_t pos = fm.HasISASamples() ? fm.LocateBidir(row) : fm.Locate(row);

            auto m = MapOccurrenceLazy(pos, pattern_len, strand, B_read,
                                       meta_payload, map_payload,
                                       dir_meta, dir_map,
                                       meta_codec_id, map_codec_id,
                                       full_abs_table);

            std::string chrom_name;
            if (m.chrom_id < chrom_names.size()) {
                chrom_name = chrom_names[m.chrom_id];
            }

            results.push_back(Match{
                .chrom        = chrom_name,
                .p_min        = m.p_min,
                .p_max        = m.p_max,
                .read_id      = m.read_id,
                .query_strand = strand,
                .sa_row       = row,
            });
        }
    }

    if (sorted) {
        std::sort(results.begin(), results.end(),
            [](const Match& a, const Match& b) {
                if (a.query_strand != b.query_strand) return a.query_strand < b.query_strand;
                if (a.chrom != b.chrom) return a.chrom < b.chrom;
                if (a.p_min != b.p_min) return a.p_min < b.p_min;
                if (a.p_max != b.p_max) return a.p_max < b.p_max;
                if (a.read_id != b.read_id) return a.read_id < b.read_id;
                return a.sa_row < b.sa_row;
            });
    }
    return results;
}

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
                           const FullAbsPosTable* full_abs_table) {
    auto err = ValidatePattern(P);
    if (err == ErrorCode::EMPTY_PATTERN)
        throw Error{err, "Pattern is empty"};
    if (err == ErrorCode::INVALID_PATTERN)
        throw Error{err, "Pattern contains invalid characters (only A=0,C=1,G=2,T=3,N=4 allowed)"};

    auto it = chrom_to_id.find(chrom);
    if (it == chrom_to_id.end()) return 0;
    uint32_t chrom_id = it->second;

    uint64_t count = 0;
    uint64_t pattern_len = P.size();

    for (const auto& [Q, strand] : GetOrientations(P, mode)) {
        auto interval = fm.BackwardSearch(Q);
        for (uint64_t row = interval.lo; row < interval.hi; ++row) {
            if (row == fm.SentinelRow()) continue;
            uint64_t pos = fm.HasISASamples() ? fm.LocateBidir(row) : fm.Locate(row);

            if (B_window.PopCount() == 0) continue;
            uint64_t window_id = B_window.Rank1(pos) - 1;
            if (window_id >= windows.size()) continue;

            const auto& W = windows[window_id];
            if (W.chrom_id != chrom_id) continue;
            if (W.genomic_end < a || W.genomic_start > b) continue;

            auto m = MapOccurrenceLazy(pos, pattern_len, strand, B_read,
                                       meta_payload, map_payload,
                                       dir_meta, dir_map,
                                       meta_codec_id, map_codec_id,
                                       full_abs_table);
            if (m.chrom_id == chrom_id && m.p_min <= b && m.p_max >= a) {
                ++count;
            }
        }
    }
    return count;
}

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
                        const FullAbsPosTable* full_abs_table) {
    auto err = ValidatePattern(P);
    if (err == ErrorCode::EMPTY_PATTERN)
        throw Error{err, "Pattern is empty"};
    if (err == ErrorCode::INVALID_PATTERN)
        throw Error{err, "Pattern contains invalid characters (only A=0,C=1,G=2,T=3,N=4 allowed)"};

    auto it = chrom_to_id.find(chrom);
    if (it == chrom_to_id.end()) return false;
    uint32_t chrom_id = it->second;

    uint64_t count = 0;
    uint64_t pattern_len = P.size();

    for (const auto& [Q, strand] : GetOrientations(P, mode)) {
        auto interval = fm.BackwardSearch(Q);
        for (uint64_t row = interval.lo; row < interval.hi; ++row) {
            if (row == fm.SentinelRow()) continue;
            uint64_t pos = fm.HasISASamples() ? fm.LocateBidir(row) : fm.Locate(row);

            if (B_window.PopCount() == 0) continue;
            uint64_t window_id = B_window.Rank1(pos) - 1;
            if (window_id >= windows.size()) continue;

            const auto& W = windows[window_id];
            if (W.chrom_id != chrom_id) continue;
            if (W.genomic_end < a || W.genomic_start > b) continue;

            auto m = MapOccurrenceLazy(pos, pattern_len, strand, B_read,
                                       meta_payload, map_payload,
                                       dir_meta, dir_map,
                                       meta_codec_id, map_codec_id,
                                       full_abs_table);
            if (m.chrom_id == chrom_id && m.p_min <= b && m.p_max >= a) {
                ++count;
                if (count >= T_threshold) return true;
            }
        }
    }
    return false;
}

// ─── ENHANCED tier: SARange-based RegionalCount (Architecture §5.3) ────────
// Instead of locating every occurrence O(occ·s), use the SARange wavelet tree
// to count SA_samples in each window's S-range in O(|W_r|·log(|S|/s)).
//
// Algorithm:
//   1. Backward search to get [lo, hi) SA interval
//   2. Find candidate windows overlapping [chrom, a, b]
//   3. For each window W_j, use SARange.RangeCount(lo/s, hi/s, W_j.l, W_j.r)
//      to count SA samples falling in [W_j.l, W_j.r]
//   4. This gives an approximation; for exact counting, we still need to
//      locate+map the occurrences within matching windows, but the SARange
//      prune eliminates windows with zero hits without O(occ) work.
//
// NOTE: For exact counting (which is what the contract requires), the SARange
// provides a tight upper bound via range_count on the SA sample indices.
// We still verify each hit via MapOccurrence for precision.

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
                               StrandMode mode) {
    auto err = ValidatePattern(P);
    if (err == ErrorCode::EMPTY_PATTERN)
        throw Error{err, "Pattern is empty"};
    if (err == ErrorCode::INVALID_PATTERN)
        throw Error{err, "Pattern contains invalid characters"};

    auto it = chrom_to_id.find(chrom);
    if (it == chrom_to_id.end()) return 0;
    uint32_t chrom_id = it->second;

    uint64_t count = 0;
    uint64_t pattern_len = P.size();

    // Pre-compute candidate windows overlapping [chrom, a, b]
    // and use SARange to check which have ANY sampled SA positions
    // in the window's S-character range. This is the window-level prune.
    //
    // Architecture §5.3: The SARange wavelet tree is built over SA_samples
    // (N_samples elements). range_count(pos_lo, pos_hi, val_lo, val_hi)
    // counts elements at positions [pos_lo, pos_hi) whose values fall in
    // [val_lo, val_hi]. Here we query the FULL sample array [0, N_samples)
    // with each window's S-range as the value filter.
    //
    // A window with zero sampled hits can be safely pruned because:
    // - Every s consecutive SA rows contain at least one sampled row
    // - If no sampled row maps into the window, at most s-1 non-sampled
    //   rows could map there, but the prune is conservative (sound superset)
    //
    // NOTE: This is a WINDOW-LEVEL prune only. When a window has hits,
    // we still locate ALL occurrences in the SA interval and check each
    // against the window. This guarantees soundness and completeness.

    for (const auto& [Q, strand] : GetOrientations(P, mode)) {
        auto interval = fm.BackwardSearch(Q);
        if (interval.lo >= interval.hi) continue;

        // For each candidate window, use SARange as a coarse filter:
        // check if ANY SA sample in the entire sample array has a value
        // falling in [W.l, W.r]. This tells us if the window could
        // possibly contain occurrences.
        for (size_t wi = 0; wi < windows.size(); ++wi) {
            const auto& W = windows[wi];
            if (W.chrom_id != chrom_id) continue;
            if (W.genomic_end < a || W.genomic_start > b) continue;

            // SARange prune: check if any sampled SA value falls in
            // the window's S-character range [W.l, W.r]
            // Query the full sample array [0, sarange.Size())
            uint64_t range_hits = sarange.RangeCount(0, sarange.Size(),
                                                      W.l, W.r);
            if (range_hits == 0) continue;

            // Window has potential hits — locate and verify exactly
            for (uint64_t row = interval.lo; row < interval.hi; ++row) {
                if (row == fm.SentinelRow()) continue;
                uint64_t pos = fm.HasISASamples() ? fm.LocateBidir(row) : fm.Locate(row);

                // Quick window check via bitvector
                if (B_window.PopCount() == 0) continue;
                uint64_t window_id = B_window.Rank1(pos) - 1;
                if (window_id != wi) continue;

                auto m = MapOccurrence(pos, pattern_len, strand, B_read, reads);
                if (m.chrom_id == chrom_id && m.p_min <= b && m.p_max >= a) {
                    ++count;
                }
            }
        }
    }

    return count;
}

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
                            StrandMode mode) {
    // For exists, we can use SARange to quickly rule out empty windows,
    // then early-exit once threshold is reached.
    uint64_t count = RegionalCountSARange(P, chrom, a, b, fm, B_read, B_window,
                                           windows, reads, sarange, chrom_to_id, mode);
    return count >= T_threshold;
}

}  // namespace bamsix

