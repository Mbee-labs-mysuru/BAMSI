#include "query.hpp"

#include <algorithm>

namespace bamsi {

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
                     StrandMode mode) {
    auto err = ValidatePattern(P);
    if (err == ErrorCode::EMPTY_PATTERN)
        throw Error{err, "Pattern is empty"};
    if (err == ErrorCode::INVALID_PATTERN)
        throw Error{err, "Pattern contains invalid characters (only A=0,C=1,G=2,T=3,N=4 allowed)"};

    uint64_t total = 0;
    for (const auto& [Q, strand] : GetOrientations(P, mode)) {
        auto interval = fm.BackwardSearch(Q);
        total += interval.size();
    }
    return total;
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

            uint64_t pos = fm.Locate(row);

            // Check for separator position (I14)
            if (IsSeparatorPosition(pos, B_read)) continue;

            // Map occurrence
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
            });
        }
    }

    if (sorted) {
        std::sort(results.begin(), results.end(),
            [](const Match& a, const Match& b) {
                // Sort per §6.2: strand → chrom → p_min → p_max → read_id
                if (a.query_strand != b.query_strand)
                    return a.query_strand < b.query_strand;
                if (a.chrom != b.chrom)
                    return a.chrom < b.chrom;
                if (a.p_min != b.p_min)
                    return a.p_min < b.p_min;
                if (a.p_max != b.p_max)
                    return a.p_max < b.p_max;
                return a.read_id < b.read_id;
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
    if (err != ErrorCode::OK) return 0;

    auto it = chrom_to_id.find(chrom);
    if (it == chrom_to_id.end()) return 0;
    uint32_t chrom_id = it->second;

    uint64_t count = 0;
    uint64_t pattern_len = P.size();

    for (const auto& [Q, strand] : GetOrientations(P, mode)) {
        auto interval = fm.BackwardSearch(Q);

        for (uint64_t row = interval.lo; row < interval.hi; ++row) {
            if (row == fm.SentinelRow()) continue;

            uint64_t pos = fm.Locate(row);

            // Window-based prune (closed-interval bitvector rank)
            if (B_window.PopCount() == 0) continue;
            uint64_t window_id = B_window.Rank1(pos) - 1;
            if (window_id >= windows.size()) continue;

            const auto& W = windows[window_id];
            if (W.chrom_id != chrom_id) continue;
            if (W.genomic_end < a || W.genomic_start > b) continue;

            // Exact intersection check
            if (IsSeparatorPosition(pos, B_read)) continue;
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
    if (err != ErrorCode::OK) return false;

    auto it = chrom_to_id.find(chrom);
    if (it == chrom_to_id.end()) return false;
    uint32_t chrom_id = it->second;

    uint64_t count = 0;
    uint64_t pattern_len = P.size();

    for (const auto& [Q, strand] : GetOrientations(P, mode)) {
        auto interval = fm.BackwardSearch(Q);

        for (uint64_t row = interval.lo; row < interval.hi; ++row) {
            if (row == fm.SentinelRow()) continue;

            uint64_t pos = fm.Locate(row);

            if (B_window.PopCount() == 0) continue;
            uint64_t window_id = B_window.Rank1(pos) - 1;
            if (window_id >= windows.size()) continue;

            const auto& W = windows[window_id];
            if (W.chrom_id != chrom_id) continue;
            if (W.genomic_end < a || W.genomic_start > b) continue;

            if (IsSeparatorPosition(pos, B_read)) continue;
            auto m = MapOccurrence(pos, pattern_len, strand, B_read, reads);
            if (m.chrom_id == chrom_id && m.p_min <= b && m.p_max >= a) {
                ++count;
                if (count >= T_threshold) return true;
            }
        }
    }

    return false;
}

}  // namespace bamsi
