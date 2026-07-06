#pragma once

#include "bamsix/types.hpp"
#include <array>
#include <cstdint>
#include <vector>

namespace bamsix {

/// Occ table for FM-index — checkpointed rank over the full BWT.
/// Provides RankBwt(a, pos): half-open count of symbol a in BWT[0..pos)
/// Used ONLY by the FM backward search and LF-mapping (Invariant I15).
class OccTable {
public:
    static constexpr uint64_t BLOCK_SIZE = 64;

    void Build(const std::vector<uint8_t>& bwt, uint8_t sigma);
    uint64_t RankBwt(uint8_t a, uint64_t pos) const;
    std::vector<uint8_t> Serialize() const;
    void Deserialize(const uint8_t* data, size_t len);

private:
    uint8_t sigma_ = 0;
    uint64_t bwt_len_ = 0;
    std::vector<uint8_t> bwt_;
    std::vector<uint64_t> occ_;
};

/// FM-Index engine implementing Architecture §4.6.
///
/// BWT has |S|+1 entries (including sentinel at sentinel_row).
/// The sentinel character is CODE_SENT (6).
///
/// C array is indexed by STORED CODE (0..6).
/// C_[c] = number of SA rows whose suffix starts with a character
///         that is lexicographically smaller than c (in the SA ordering).
///
/// SA ordering (shifted alphabet used by libsais): $ < A < C < G < T < N < #
/// Stored codes:     A=0, C=1, G=2, T=3, N=4, #=5, $=CODE_SENT=6
///
/// The LF formula is: LF(r) = C_[BWT[r]] + Occ(BWT[r], r)
/// This works directly because both C_ and Occ are indexed by stored code.
class FMIndexEngine {
public:
    void Build(const std::vector<uint8_t>& bwt,
               const std::vector<int64_t>& sa,
               uint64_t sentinel,
               uint64_t sample_s,
               uint64_t s_len);

    SAInterval BackwardSearch(const std::vector<uint8_t>& pattern) const;
    uint64_t LF(uint64_t row) const;
    uint64_t Locate(uint64_t row) const;

    /// Locate with ISA samples (bidirectional, Architecture §4.6).
    /// Uses min(SA-walk, ISA-walk) to find position faster.
    /// Only available when ISA samples are loaded.
    uint64_t LocateBidir(uint64_t row) const;
    bool HasISASamples() const { return !isa_samples_.empty(); }

    /// Set ISA samples (called during build or deserialization).
    void SetISASamples(std::vector<uint64_t> isa_samples, uint64_t isa_step);
    const std::vector<uint64_t>& ISASamples() const { return isa_samples_; }
    uint64_t ISAStep() const { return isa_step_; }

    uint8_t BwtAt(uint64_t row) const { return bwt_[row]; }
    uint64_t SentinelRow() const { return sentinel_row_; }
    uint64_t SampleStep() const { return sample_step_; }
    uint64_t SLen() const { return s_len_; }

    /// For serialization: returns BWT with sentinel stripped (Architecture §4.4).
    std::vector<uint8_t> SerializeBWT() const;
    const std::array<uint64_t, 7>& CArray() const { return C_; }
    const OccTable& Occ() const { return occ_; }
    const std::vector<uint64_t>& SASamples() const { return sa_samples_; }

    /// Load FM-index from stored .bsi data (deserialization path).
    void LoadFromStored(const std::vector<uint8_t>& full_bwt,
                        const std::array<uint64_t, 7>& c_array,
                        const std::vector<uint8_t>& occ_data,
                        const std::vector<uint64_t>& sa_samples,
                        uint64_t sentinel_row,
                        uint64_t sample_step,
                        uint64_t s_len);

private:
    std::vector<uint8_t> bwt_;  // Full BWT, length |S|+1
    uint64_t bwt_len_ = 0;
    uint64_t sentinel_row_ = 0;
    uint64_t sample_step_ = 64;
    uint64_t s_len_ = 0;

    // C_[stored_code] = number of SA rows with a lex-smaller starting char.
    // Indexed by stored code (0..6), NOT by lex rank.
    std::array<uint64_t, 7> C_ = {};

    OccTable occ_;
    std::vector<uint64_t> sa_samples_;  // row-based: sa_samples_[k] = SA[k*s]
    std::vector<uint64_t> isa_samples_; // text-based: isa_samples_[k] = ISA[k*s']
    uint64_t isa_step_ = 0;            // ISA sample step s'
};

/// Reverse FM-Index for bidirectional search support (Architecture §4.6.7).
/// Built only when enable_bidirectional = true.
/// The reverse SA is derived from the forward SA without a second SA-IS run:
///   SA_R[i] = |S| - 1 - SA[|S| - 1 - i]
class ReverseFMIndex {
public:
    /// Build from forward SA and original sequence S (Architecture §4.6.7).
    void Build(const std::vector<int64_t>& forward_sa,
               const std::vector<uint8_t>& S,
               uint64_t sentinel_row_fwd,
               uint64_t sample_s,
               uint64_t s_len);

    SAInterval BackwardSearch(const std::vector<uint8_t>& pattern) const;
    uint64_t Locate(uint64_t row) const;

    bool IsBuilt() const { return !bwt_.empty(); }

    // Serialization
    const std::vector<uint8_t>& BWT() const { return bwt_; }
    const std::array<uint64_t, 7>& CArray() const { return C_; }
    const OccTable& Occ() const { return occ_; }
    const std::vector<uint64_t>& SASamples() const { return sa_samples_; }
    uint64_t SentinelRow() const { return sentinel_row_; }

    void LoadFromStored(const std::vector<uint8_t>& full_bwt,
                        const std::array<uint64_t, 7>& c_array,
                        const std::vector<uint8_t>& occ_data,
                        const std::vector<uint64_t>& sa_samples,
                        uint64_t sentinel_row,
                        uint64_t sample_step,
                        uint64_t s_len);

private:
    std::vector<uint8_t> bwt_;
    uint64_t bwt_len_ = 0;
    uint64_t sentinel_row_ = 0;
    uint64_t sample_step_ = 64;
    uint64_t s_len_ = 0;
    std::array<uint64_t, 7> C_ = {};
    OccTable occ_;
    std::vector<uint64_t> sa_samples_;
};

}  // namespace bamsix
