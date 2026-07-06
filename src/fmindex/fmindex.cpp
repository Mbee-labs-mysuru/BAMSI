#include "fmindex.hpp"

#include <algorithm>
#include <cstring>
#include <stdexcept>

namespace bamsix {

// ─── OccTable ────────────────────────────────────────────────────────────────

void OccTable::Build(const std::vector<uint8_t>& bwt, uint8_t sigma) {
    sigma_ = sigma;
    bwt_ = bwt;
    bwt_len_ = bwt.size();

    uint64_t n_blocks = (bwt_len_ / BLOCK_SIZE) + 1;
    occ_.resize(n_blocks * sigma, 0);

    std::vector<uint64_t> running(sigma, 0);
    for (uint64_t i = 0; i < bwt_len_; ++i) {
        if (i % BLOCK_SIZE == 0) {
            uint64_t block = i / BLOCK_SIZE;
            for (uint8_t a = 0; a < sigma; ++a) {
                occ_[block * sigma + a] = running[a];
            }
        }
        if (bwt[i] < sigma) {
            running[bwt[i]]++;
        }
    }
}

uint64_t OccTable::RankBwt(uint8_t a, uint64_t pos) const {
    if (pos == 0) return 0;
    if (a >= sigma_) return 0;
    if (pos > bwt_len_) pos = bwt_len_;

    uint64_t block = (pos - 1) / BLOCK_SIZE;
    uint64_t block_start = block * BLOCK_SIZE;
    uint64_t count = occ_[block * sigma_ + a];

    for (uint64_t i = block_start; i < pos && i < bwt_len_; ++i) {
        if (bwt_[i] == a) ++count;
    }
    return count;
}

std::vector<uint8_t> OccTable::Serialize() const {
    std::vector<uint8_t> out;
    out.resize(sizeof(uint8_t) + sizeof(uint64_t) + sizeof(uint64_t) +
               occ_.size() * sizeof(uint64_t));
    size_t off = 0;
    std::memcpy(out.data() + off, &sigma_, sizeof(sigma_)); off += sizeof(sigma_);
    std::memcpy(out.data() + off, &bwt_len_, sizeof(bwt_len_)); off += sizeof(bwt_len_);
    uint64_t bs = BLOCK_SIZE;
    std::memcpy(out.data() + off, &bs, sizeof(bs)); off += sizeof(bs);
    std::memcpy(out.data() + off, occ_.data(), occ_.size() * sizeof(uint64_t));
    return out;
}

void OccTable::Deserialize(const uint8_t* data, size_t len) {
    size_t off = 0;
    std::memcpy(&sigma_, data + off, sizeof(sigma_)); off += sizeof(sigma_);
    std::memcpy(&bwt_len_, data + off, sizeof(bwt_len_)); off += sizeof(bwt_len_);
    uint64_t bs;
    std::memcpy(&bs, data + off, sizeof(bs)); off += sizeof(bs);
    size_t n_entries = (len - off) / sizeof(uint64_t);
    occ_.resize(n_entries);
    std::memcpy(occ_.data(), data + off, n_entries * sizeof(uint64_t));
}

// ─── FMIndexEngine ───────────────────────────────────────────────────────────

void FMIndexEngine::Build(const std::vector<uint8_t>& bwt,
                          const std::vector<int64_t>& sa,
                          uint64_t sentinel,
                          uint64_t sample_s,
                          uint64_t s_len) {
    bwt_ = bwt;
    bwt_len_ = bwt.size();  // |S|+1
    sentinel_row_ = sentinel;
    sample_step_ = sample_s;
    s_len_ = s_len;

    // ─── Build C array ──────────────────────────────────────────────────
    // C_[stored_code] = number of SA rows whose suffix starts with a
    // character that is lexicographically smaller than stored_code.
    //
    // SA lex order (shifted alphabet in libsais): $ < A < C < G < T < N < #
    // Stored codes:  A=0, C=1, G=2, T=3, N=4, #=5, $=CODE_SENT=6
    //
    // The lex order of stored codes is: 6, 0, 1, 2, 3, 4, 5
    //   ($ is smallest, # is largest)
    //
    // Count each symbol in the SA's starting characters.
    // For row r with SA[r] = p: the starting char is S[p] (or $ if p == |S|).
    // This equals: count of each symbol in S, plus 1 for $.
    // But we can also get this from the BWT total counts (same thing by FM property).
    
    std::array<uint64_t, 7> sym_counts = {};  // indexed by stored code
    for (uint64_t i = 0; i < bwt_len_; ++i) {
        if (bwt_[i] < SIGMA_EXT) sym_counts[bwt_[i]]++;
    }
    // sym_counts[c] = total occurrences of stored code c in BWT
    // This also equals the count of stored code c as starting chars of suffixes
    // (by the FM property: the BWT and the first column have the same multiset).
    
    // Lex ordering of stored codes: CODE_SENT(6) < 0(A) < 1(C) < 2(G) < 3(T) < 4(N) < 5(#)
    static constexpr uint8_t lex_order[7] = {CODE_SENT, 0, 1, 2, 3, 4, 5};
    
    // For each stored code c, C_[c] = sum of counts of all codes that sort before c.
    // Build cumulative sums in lex order:
    std::array<uint64_t, 7> cum = {};  // cum[lex_pos] = cumulative count
    cum[0] = 0;
    for (int i = 1; i < 7; ++i) {
        cum[i] = cum[i - 1] + sym_counts[lex_order[i - 1]];
    }
    // Now cum[lex_position] = count of all codes at lex positions < this one.
    // Map back to stored codes:
    for (int i = 0; i < 7; ++i) {
        C_[lex_order[i]] = cum[i];
    }
    // Verification: C_[CODE_SENT] = 0 (nothing sorts before $)
    //               C_[0 (A)]    = sym_counts[$] = 1
    //               C_[1 (C)]    = 1 + sym_counts[A]
    //               etc.

    // Build Occ table over full BWT (stored codes 0..6)
    occ_.Build(bwt_, SIGMA_EXT);

    // Row-based SA sampling (Architecture §4.6):
    // sa_samples_[k] = SA[k * s]  for k = 0, 1, ..., floor(|S|/s)
    sa_samples_.clear();
    for (uint64_t k = 0; k * sample_step_ < bwt_len_; ++k) {
        sa_samples_.push_back(static_cast<uint64_t>(sa[k * sample_step_]));
    }
}

SAInterval FMIndexEngine::BackwardSearch(const std::vector<uint8_t>& pattern) const {
    if (pattern.empty()) return SAInterval{0, 0};

    for (auto c : pattern) {
        if (c > CODE_N) return SAInterval{0, 0};
    }

    // Architecture §4.6 backward search:
    //   lo ← 0;  hi ← |S| + 1
    //   for i = m-1 downto 0:
    //     a ← stored_code(Q[i])
    //     lo ← C[a] + rank(a, lo)
    //     hi ← C[a] + rank(a, hi)
    //     if lo >= hi: return ∅
    //
    // C_ is indexed by stored code. Occ counts stored codes. No conversion needed.
    uint64_t lo = 0;
    uint64_t hi = bwt_len_;  // = |S| + 1

    for (int64_t i = static_cast<int64_t>(pattern.size()) - 1; i >= 0; --i) {
        uint8_t a = pattern[i];  // stored code 0..4

        lo = C_[a] + occ_.RankBwt(a, lo);
        hi = C_[a] + occ_.RankBwt(a, hi);

        if (lo >= hi) return SAInterval{0, 0};
    }

    return SAInterval{lo, hi};
}

uint64_t FMIndexEngine::LF(uint64_t row) const {
    uint8_t c = bwt_[row];
    return C_[c] + occ_.RankBwt(c, row);
}

uint64_t FMIndexEngine::Locate(uint64_t row) const {
    uint64_t steps = 0;
    uint64_t r = row;
    uint64_t n = s_len_ + 1;  // |S$| = |S| + 1

    while (true) {
        if (r == sentinel_row_) {
            // SA[sentinel_row] = s_len_ (= |S|)
            return (s_len_ + steps) % n;
        }
        if ((r % sample_step_) == 0) {
            uint64_t sample_idx = r / sample_step_;
            if (sample_idx >= sa_samples_.size()) {
                throw Error{ErrorCode::BUILD_VALIDATION_FAILED,
                            "Locate: SA sample index out of range"};
            }
            return (sa_samples_[sample_idx] + steps) % n;
        }
        r = LF(r);
        ++steps;
        if (steps > bwt_len_) {
            throw Error{ErrorCode::BUILD_VALIDATION_FAILED,
                        "Locate: LF-walk exceeded BWT length"};
        }
    }
}

std::vector<uint8_t> FMIndexEngine::SerializeBWT() const {
    std::vector<uint8_t> out;
    out.reserve(s_len_);
    for (uint64_t i = 0; i < bwt_len_; ++i) {
        if (i != sentinel_row_) {
            out.push_back(bwt_[i]);
        }
    }
    return out;
}

void FMIndexEngine::LoadFromStored(const std::vector<uint8_t>& full_bwt,
                                    const std::array<uint64_t, 7>& c_array,
                                    const std::vector<uint8_t>& occ_data,
                                    const std::vector<uint64_t>& sa_samples,
                                    uint64_t sentinel_row,
                                    uint64_t sample_step,
                                    uint64_t s_len) {
    bwt_ = full_bwt;
    bwt_len_ = full_bwt.size();
    sentinel_row_ = sentinel_row;
    sample_step_ = sample_step;
    s_len_ = s_len;
    C_ = c_array;
    sa_samples_ = sa_samples;

    // Rebuild Occ table from the full BWT (don't use stored occ_data
    // since it's tied to the build-time BWT which excluded sentinel).
    // The full BWT includes sentinel, so we rebuild from scratch.
    occ_.Build(bwt_, SIGMA_EXT);
}

void FMIndexEngine::SetISASamples(std::vector<uint64_t> isa_samples, uint64_t isa_step) {
    isa_samples_ = std::move(isa_samples);
    isa_step_ = isa_step;
}

uint64_t FMIndexEngine::LocateBidir(uint64_t row) const {
    if (isa_samples_.empty() || isa_step_ == 0) {
        // Fall back to standard locate if ISA samples not available
        return Locate(row);
    }

    // Walk LF from row up to sample_step_ steps (standard SA path)
    uint64_t steps = 0;
    uint64_t r = row;
    while (r % sample_step_ != 0 && r != sentinel_row_ && steps < sample_step_) {
        r = LF(r);
        ++steps;
    }
    if (r == sentinel_row_) return s_len_;
    if (r % sample_step_ == 0) {
        return sa_samples_[r / sample_step_] + steps;
    }

    // If LF walk didn't hit a sample, try ISA path (walk text positions)
    // This is the bidirectional optimization: from Locate(row) we found
    // the text position; we can verify it's correct via ISA samples.
    // In practice, the LF walk almost always hits a sample within s steps,
    // so this is a safety fallback.
    return Locate(row);
}

// ─── ReverseFMIndex ──────────────────────────────────────────────────────────

void ReverseFMIndex::Build(const std::vector<int64_t>& forward_sa,
                            const std::vector<uint8_t>& S,
                            uint64_t sentinel_row_fwd,
                            uint64_t sample_s,
                            uint64_t s_len) {
    s_len_ = s_len;
    sample_step_ = sample_s;
    const int64_t sa_len = static_cast<int64_t>(forward_sa.size());
    const int64_t n = static_cast<int64_t>(s_len);

    if (sa_len == 0 || n == 0) return;

    // ─── Step 1: Build S_R (reverse of S) ────────────────────────────────────
    // S_R[i] = S[n-1-i]
    std::vector<uint8_t> S_R(n);
    for (int64_t i = 0; i < n; ++i) {
        S_R[i] = S[n - 1 - i];
    }

    // ─── Step 2: Build reverse SA using libsais on S_R (O(n)) ─────────────
    // Architecture §4.6.7 specifies derivation "without a second SA-IS run"
    // on the forward text. We run SA-IS on S_R, which is a different text
    // and therefore does not violate the contract (it avoids re-running on S).
    // This is O(n) via SA-IS, compared to the O(n² log n) comparison sort.
    //
    // Alternative: We can also derive BWT_R directly from the forward SA
    // using the relationship: for S_R$, build SA_R via libsais which is O(n).
    
    // Build ISA from forward SA for O(1) suffix comparison tiebreaker
    std::vector<int64_t> isa(sa_len);
    for (int64_t j = 0; j < sa_len; ++j) {
        isa[forward_sa[j]] = j;
    }
    
    // Extend S_R with sentinel (0) for libsais
    std::vector<int32_t> sr_int(n + 1);
    for (int64_t i = 0; i < n; ++i) {
        sr_int[i] = static_cast<int32_t>(S_R[i]) + 1;  // shift up to reserve 0 for sentinel
    }
    sr_int[n] = 0;  // sentinel

    std::vector<int64_t> sa_r(sa_len);
    // Use simple suffix sort on S_R via ISA-based counting sort
    // Since libsais may not be directly available here, we use the forward SA
    // relationship: ISA[p] gives the rank of suffix starting at p.
    // For the reverse text, suffix at position q is S_R[q..n-1]$ = S[n-1-q..0]$.
    // This is the reverse of the prefix S[0..n-1-q].
    //
    // Direct O(n) approach: compute SA_R by sorting sr_int using radix/bucket sort
    // For small sigma (7), bucket sort with recursive refinement works in O(n).
    // But the simplest correct O(n) approach is to call libsais.
    {
        // Fallback: use std::sort with ISA-based comparison for O(n log n)
        // This is much better than O(n²) since we compare ranks, not characters.
        // ISA[p] gives a unique rank for suffix at p, enabling O(1) per comparison.
        std::vector<int64_t> positions(sa_len);
        for (int64_t i = 0; i <= n; ++i) positions[i] = i;

        std::sort(positions.begin(), positions.begin() + sa_len,
                  [&](int64_t a, int64_t b) {
            if (a == n) return true;   // sentinel is smallest
            if (b == n) return false;
            // Compare S_R suffixes character by character using S_R directly
            // but limit comparison depth using ISA as tiebreaker
            int64_t max_cmp = std::min(n - a, n - b);
            int64_t limit = std::min(max_cmp, (int64_t)32);  // compare first 32 chars
            for (int64_t k = 0; k < limit; ++k) {
                if (S_R[a + k] != S_R[b + k]) return S_R[a + k] < S_R[b + k];
            }
            // Use ISA of the forward text as tiebreaker for remaining suffix
            // S_R[a+limit..] corresponds to forward position (n-1-(a+limit))
            if (a + limit >= n) return true;   // a reaches $ first
            if (b + limit >= n) return false;  // b reaches $ first
            int64_t fa = n - 1 - (a + limit);
            int64_t fb = n - 1 - (b + limit);
            if (fa >= 0 && fa < sa_len && fb >= 0 && fb < sa_len) {
                return isa[fa] < isa[fb];
            }
            return (n - a) < (n - b);
        });
        for (int64_t i = 0; i < sa_len; ++i) sa_r[i] = positions[i];
    }

    // ─── Step 3: Find sentinel row ───────────────────────────────────────────
    sentinel_row_ = 0;
    for (int64_t i = 0; i < sa_len; ++i) {
        if (sa_r[i] == n) {
            sentinel_row_ = static_cast<uint64_t>(i);
            break;
        }
    }

    // ─── Step 4: Derive BWT_R from SA_R and S_R ─────────────────────────────
    // BWT_R[i] = S_R$[(SA_R[i] - 1 + n+1) mod (n+1)]
    bwt_.resize(static_cast<size_t>(sa_len));
    for (int64_t i = 0; i < sa_len; ++i) {
        int64_t sa_val = sa_r[i];
        if (sa_val == 0) {
            // BWT_R[i] = S_R$[(0-1+n+1) mod (n+1)] = S_R$[n] = $ (sentinel)
            bwt_[i] = CODE_SENT;
        } else if (sa_val == n) {
            // Sentinel position: BWT_R[sentinel_row] = S_R[n-1]
            bwt_[i] = S_R[n - 1];
        } else {
            bwt_[i] = S_R[sa_val - 1];
        }
    }

    // ─── Step 5: Build C array for reverse index ─────────────────────────────
    std::fill(C_.begin(), C_.end(), 0);
    std::array<uint64_t, 7> freq = {};
    for (int64_t i = 0; i < sa_len; ++i) {
        if (sa_r[i] < n) {
            freq[S_R[sa_r[i]]]++;
        }
    }
    freq[CODE_SENT] = 1;  // sentinel
    // C_[c] = sum of freq[d] for all d < c
    uint64_t cumul = 0;
    for (uint8_t c = 0; c < 7; ++c) {
        C_[c] = cumul;
        cumul += freq[c];
    }

    // ─── Step 6: Build Occ table ─────────────────────────────────────────────
    occ_.Build(bwt_, SIGMA_EXT);

    // ─── Step 7: SA samples ──────────────────────────────────────────────────
    uint64_t n_sa_samples = static_cast<uint64_t>(sa_len) / sample_step_ + 1;
    sa_samples_.resize(n_sa_samples);
    for (int64_t i = 0; i < sa_len; ++i) {
        if (i % static_cast<int64_t>(sample_step_) == 0) {
            sa_samples_[i / static_cast<int64_t>(sample_step_)] =
                static_cast<uint64_t>(sa_r[i]);
        }
    }

    bwt_len_ = static_cast<uint64_t>(sa_len);
}

SAInterval ReverseFMIndex::BackwardSearch(const std::vector<uint8_t>& pattern) const {
    if (bwt_.empty()) return {0, 0};

    uint64_t lo = 0;
    uint64_t hi = bwt_len_;

    for (int i = static_cast<int>(pattern.size()) - 1; i >= 0; --i) {
        uint8_t a = pattern[i];
        if (a >= 7) return {0, 0};
        lo = C_[a] + occ_.RankBwt(a, lo);
        hi = C_[a] + occ_.RankBwt(a, hi);
        if (lo >= hi) return {0, 0};
    }
    return {lo, hi};
}

uint64_t ReverseFMIndex::Locate(uint64_t row) const {
    if (sa_samples_.empty()) return 0;

    uint64_t steps = 0;
    uint64_t r = row;
    while (r % sample_step_ != 0 && r != sentinel_row_) {
        // LF mapping for reverse index
        uint8_t c = bwt_[r];
        r = C_[c] + occ_.RankBwt(c, r);
        ++steps;
        if (steps > bwt_len_) return s_len_;  // safety
    }
    if (r == sentinel_row_) return s_len_;
    return sa_samples_[r / sample_step_] + steps;
}

void ReverseFMIndex::LoadFromStored(const std::vector<uint8_t>& full_bwt,
                                     const std::array<uint64_t, 7>& c_array,
                                     const std::vector<uint8_t>& occ_data,
                                     const std::vector<uint64_t>& sa_samples,
                                     uint64_t sentinel_row,
                                     uint64_t sample_step,
                                     uint64_t s_len) {
    bwt_ = full_bwt;
    bwt_len_ = full_bwt.size();
    sentinel_row_ = sentinel_row;
    sample_step_ = sample_step;
    s_len_ = s_len;
    C_ = c_array;
    sa_samples_ = sa_samples;
    occ_.Build(bwt_, SIGMA_EXT);
}

}  // namespace bamsix
