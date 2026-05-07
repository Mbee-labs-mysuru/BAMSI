#pragma once

#include "bamsi/types.hpp"
#include <vector>
#include <cstdint>

namespace bamsi {

/// Succinct bitvector with O(1) rank1 and select1 support.
/// Uses CLOSED-INTERVAL rank convention: rank1(B, pos) = count of 1s in B[0..pos] inclusive.
/// This is DISTINCT from the half-open BWT rank used by OccTable (I15).
class SuccinctBitvector {
public:
    SuccinctBitvector() = default;

    /// Build from a list of positions that should be set to 1.
    /// @param ones  Sorted list of positions where bit = 1.
    /// @param total Total bitvector length.
    void Build(const std::vector<uint64_t>& ones, uint64_t total);

    /// Closed-interval rank: count of 1-bits in B[0..pos] inclusive.
    /// Used for read-id identification and window lookup.
    uint64_t Rank1(uint64_t pos) const;

    /// Select: find the position of the k-th 1-bit (1-indexed).
    /// select1(1) returns the position of the first 1-bit.
    uint64_t Select1(uint64_t k) const;

    /// Access bit at position pos.
    bool Access(uint64_t pos) const;

    /// Total length of the bitvector.
    uint64_t Length() const { return total_; }

    /// Number of 1-bits.
    uint64_t PopCount() const { return ones_count_; }

    /// Serialize to bytes.
    std::vector<uint8_t> Serialize() const;

    /// Deserialize from bytes.
    void Deserialize(const uint8_t* data, size_t len);

private:
    static constexpr uint64_t BLOCK_BITS = 64;
    static constexpr uint64_t SUPERBLOCK_SIZE = 256;  // bits per superblock

    std::vector<uint64_t> blocks_;       // raw bitvector words
    std::vector<uint64_t> superblocks_;  // cumulative rank at each superblock
    uint64_t total_ = 0;
    uint64_t ones_count_ = 0;

    // For O(1) select: store positions of all 1-bits
    std::vector<uint64_t> select_idx_;
};

/// Bitvector pair for the BAMSI build: B_read and B_window.
struct BitVectors {
    SuccinctBitvector B_read;
    SuccinctBitvector B_window;
};

/// Build B_read from readStarts and B_window from window starts.
BitVectors BuildBitvectors(const SequenceBundle& bundle, const WindowTable& windows);

}  // namespace bamsi
