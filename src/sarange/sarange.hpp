#pragma once

#include <cstdint>
#include <vector>

namespace bamsi {

/// SARange wavelet tree — ENHANCED tier (Architecture §5.3).
///
/// A wavelet tree over SA_samples[0..floor(|S|/s)] supporting:
///   range_count(lo, hi, l, r) → count of SA_samples[i] in [l, r] for i in [lo, hi)
///   Complexity: O(log(max_val / min_val))
///
/// Used by ENHANCED-tier RegionalCount to count pattern occurrences
/// within a genomic window without locating all occurrences.
///
/// The wavelet tree is a balanced binary tree where each node stores a
/// bitvector indicating whether each element goes to the left (0) or right (1)
/// child based on the corresponding bit of the value.
class SARange {
public:
    /// Build wavelet tree over the given SA samples.
    /// max_val should be the maximum possible SA value (typically |S|).
    void Build(const std::vector<uint64_t>& sa_samples, uint64_t max_val);

    /// Count elements in sa_samples[lo..hi) that fall within [val_lo, val_hi].
    /// Complexity: O(log(max_val))
    /// Per Architecture §5.3: range_count(lo, hi, l_j, r_j)
    uint64_t RangeCount(uint64_t lo, uint64_t hi,
                        uint64_t val_lo, uint64_t val_hi) const;

    bool IsBuilt() const { return !levels_.empty(); }
    uint64_t Size() const { return n_; }

    /// Serialization support.
    std::vector<uint8_t> Serialize() const;
    void Deserialize(const uint8_t* data, size_t len);

private:
    struct Level {
        std::vector<uint8_t> bits;       // bitvector (packed)
        std::vector<uint64_t> rank0;     // prefix count of 0-bits, sampled
        uint64_t n_zeros = 0;           // total 0-bits at this level
    };

    uint64_t n_ = 0;         // number of elements
    uint64_t max_val_ = 0;   // universe size
    uint32_t n_levels_ = 0;  // number of levels (ceil(log2(max_val)))
    std::vector<Level> levels_;

    // Internal helpers
    bool GetBit(const Level& level, uint64_t pos) const;
    uint64_t Rank0(const Level& level, uint64_t pos) const;
    uint64_t Rank1(const Level& level, uint64_t pos) const;
};

}  // namespace bamsi
