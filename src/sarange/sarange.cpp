#include "sarange.hpp"

#include <algorithm>
#include <cstring>

namespace bamsi {

// ─── Bit helpers ─────────────────────────────────────────────────────────────

bool SARange::GetBit(const Level& level, uint64_t pos) const {
    if (pos >= n_) return false;
    return (level.bits[pos / 8] >> (pos % 8)) & 1;
}

uint64_t SARange::Rank0(const Level& level, uint64_t pos) const {
    // Count 0-bits in [0, pos)
    uint64_t count = 0;
    uint64_t full_bytes = pos / 8;
    for (uint64_t b = 0; b < full_bytes; ++b) {
        count += 8 - __builtin_popcount(level.bits[b]);
    }
    uint64_t rem = pos % 8;
    if (rem > 0 && full_bytes < level.bits.size()) {
        uint8_t byte = level.bits[full_bytes];
        for (uint64_t bit = 0; bit < rem; ++bit) {
            if (!((byte >> bit) & 1)) ++count;
        }
    }
    return count;
}

uint64_t SARange::Rank1(const Level& level, uint64_t pos) const {
    return pos - Rank0(level, pos);
}

// ─── Build ───────────────────────────────────────────────────────────────────

void SARange::Build(const std::vector<uint64_t>& sa_samples, uint64_t max_val) {
    n_ = sa_samples.size();
    max_val_ = max_val;

    if (n_ == 0 || max_val_ == 0) return;

    // Number of levels = ceil(log2(max_val + 1))
    n_levels_ = 0;
    uint64_t mv = max_val_;
    while ((1ULL << n_levels_) <= mv) ++n_levels_;
    if (n_levels_ == 0) n_levels_ = 1;

    levels_.resize(n_levels_);

    // Working arrays: current permutation of values at each level
    std::vector<uint64_t> current = sa_samples;
    std::vector<uint64_t> left_buf, right_buf;
    left_buf.reserve(n_);
    right_buf.reserve(n_);

    for (uint32_t lev = 0; lev < n_levels_; ++lev) {
        uint32_t bit_pos = n_levels_ - 1 - lev;  // MSB first

        // Allocate bitvector
        levels_[lev].bits.resize((n_ + 7) / 8, 0);
        levels_[lev].n_zeros = 0;

        left_buf.clear();
        right_buf.clear();

        for (uint64_t i = 0; i < current.size(); ++i) {
            bool bit = (current[i] >> bit_pos) & 1;
            if (bit) {
                levels_[lev].bits[i / 8] |= (1 << (i % 8));
                right_buf.push_back(current[i]);
            } else {
                levels_[lev].n_zeros++;
                left_buf.push_back(current[i]);
            }
        }

        // Next level: left children first, then right children
        current.clear();
        current.insert(current.end(), left_buf.begin(), left_buf.end());
        current.insert(current.end(), right_buf.begin(), right_buf.end());
    }
}

// ─── RangeCount ──────────────────────────────────────────────────────────────

uint64_t SARange::RangeCount(uint64_t lo, uint64_t hi,
                              uint64_t val_lo, uint64_t val_hi) const {
    if (!IsBuilt() || lo >= hi || val_lo > val_hi) return 0;
    if (hi > n_) hi = n_;

    // Recursive wavelet tree range count using the standard algorithm.
    // We use an iterative stack-based approach for efficiency.

    struct Frame {
        uint64_t lo, hi;       // range in the wavelet tree
        uint64_t node_lo, node_hi;  // value range represented by this node
        uint32_t level;
    };

    uint64_t count = 0;
    std::vector<Frame> stack;
    stack.push_back({lo, hi, 0, (1ULL << n_levels_) - 1, 0});

    while (!stack.empty()) {
        auto f = stack.back();
        stack.pop_back();

        if (f.lo >= f.hi) continue;
        if (f.node_lo > val_hi || f.node_hi < val_lo) continue;

        // If the node range is fully contained in [val_lo, val_hi]
        if (val_lo <= f.node_lo && f.node_hi <= val_hi) {
            count += f.hi - f.lo;
            continue;
        }

        if (f.level >= n_levels_) {
            // Leaf level
            if (f.node_lo >= val_lo && f.node_lo <= val_hi) {
                count += f.hi - f.lo;
            }
            continue;
        }

        const auto& level = levels_[f.level];
        uint64_t mid = (f.node_lo + f.node_hi) / 2;

        // Count 0-bits and 1-bits in [f.lo, f.hi)
        uint64_t r0_lo = Rank0(level, f.lo);
        uint64_t r0_hi = Rank0(level, f.hi);
        uint64_t r1_lo = f.lo - r0_lo;
        uint64_t r1_hi = f.hi - r0_hi;

        // Left child: values in [node_lo, mid] — the 0-bit elements
        if (r0_lo < r0_hi) {
            stack.push_back({r0_lo, r0_hi, f.node_lo, mid, f.level + 1});
        }

        // Right child: values in [mid+1, node_hi] — the 1-bit elements
        if (r1_lo < r1_hi) {
            uint64_t n_zeros = level.n_zeros;
            stack.push_back({n_zeros + r1_lo, n_zeros + r1_hi,
                             mid + 1, f.node_hi, f.level + 1});
        }
    }

    return count;
}

// ─── Serialization ───────────────────────────────────────────────────────────

std::vector<uint8_t> SARange::Serialize() const {
    std::vector<uint8_t> out;
    // Header: n_, max_val_, n_levels_
    out.resize(sizeof(uint64_t) * 2 + sizeof(uint32_t));
    size_t off = 0;
    std::memcpy(out.data() + off, &n_, sizeof(n_)); off += sizeof(n_);
    std::memcpy(out.data() + off, &max_val_, sizeof(max_val_)); off += sizeof(max_val_);
    std::memcpy(out.data() + off, &n_levels_, sizeof(n_levels_)); off += sizeof(n_levels_);

    for (const auto& level : levels_) {
        uint64_t bits_size = level.bits.size();
        uint64_t n_zeros = level.n_zeros;

        size_t old_size = out.size();
        out.resize(old_size + sizeof(uint64_t) * 2 + bits_size);
        std::memcpy(out.data() + old_size, &bits_size, sizeof(bits_size));
        std::memcpy(out.data() + old_size + sizeof(uint64_t), &n_zeros, sizeof(n_zeros));
        std::memcpy(out.data() + old_size + sizeof(uint64_t) * 2,
                    level.bits.data(), bits_size);
    }
    return out;
}

void SARange::Deserialize(const uint8_t* data, size_t len) {
    if (len < sizeof(uint64_t) * 2 + sizeof(uint32_t)) return;

    size_t off = 0;
    std::memcpy(&n_, data + off, sizeof(n_)); off += sizeof(n_);
    std::memcpy(&max_val_, data + off, sizeof(max_val_)); off += sizeof(max_val_);
    std::memcpy(&n_levels_, data + off, sizeof(n_levels_)); off += sizeof(n_levels_);

    levels_.resize(n_levels_);
    for (uint32_t lev = 0; lev < n_levels_; ++lev) {
        if (off + sizeof(uint64_t) * 2 > len) break;

        uint64_t bits_size = 0, n_zeros = 0;
        std::memcpy(&bits_size, data + off, sizeof(bits_size)); off += sizeof(uint64_t);
        std::memcpy(&n_zeros, data + off, sizeof(n_zeros)); off += sizeof(uint64_t);

        levels_[lev].n_zeros = n_zeros;
        levels_[lev].bits.resize(bits_size);
        if (off + bits_size <= len) {
            std::memcpy(levels_[lev].bits.data(), data + off, bits_size);
        }
        off += bits_size;
    }
}

}  // namespace bamsi
