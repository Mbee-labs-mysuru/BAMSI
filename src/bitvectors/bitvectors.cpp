#include "bitvectors.hpp"

#include <algorithm>
#include <cstring>
#include <bit>

namespace bamsi {

// ─────────────────────────────────────────────────────────────────────────────
// SuccinctBitvector
// ─────────────────────────────────────────────────────────────────────────────

void SuccinctBitvector::Build(const std::vector<uint64_t>& ones, uint64_t total) {
    total_ = total;
    ones_count_ = ones.size();

    // Allocate blocks (ceil(total / 64) words)
    uint64_t n_words = (total + BLOCK_BITS - 1) / BLOCK_BITS;
    blocks_.assign(n_words, 0);

    // Set bits
    for (uint64_t pos : ones) {
        if (pos < total) {
            blocks_[pos / BLOCK_BITS] |= (1ULL << (pos % BLOCK_BITS));
        }
    }

    // Build superblock rank index
    uint64_t n_superblocks = (total + SUPERBLOCK_SIZE - 1) / SUPERBLOCK_SIZE + 1;
    superblocks_.resize(n_superblocks, 0);

    uint64_t running = 0;
    for (uint64_t i = 0; i < n_words; ++i) {
        if ((i * BLOCK_BITS) % SUPERBLOCK_SIZE == 0) {
            superblocks_[(i * BLOCK_BITS) / SUPERBLOCK_SIZE] = running;
        }
        running += std::popcount(blocks_[i]);
    }
    // Final superblock
    uint64_t last_sb = total / SUPERBLOCK_SIZE;
    if (last_sb < n_superblocks) {
        // Already handled above
    }

    // Build select index: store all 1-bit positions for O(1) select
    select_idx_.clear();
    select_idx_.reserve(ones_count_);
    for (uint64_t i = 0; i < n_words; ++i) {
        uint64_t word = blocks_[i];
        while (word) {
            uint64_t bit_pos = i * BLOCK_BITS + std::countr_zero(word);
            if (bit_pos < total) {
                select_idx_.push_back(bit_pos);
            }
            word &= word - 1;  // clear lowest set bit
        }
    }
    std::sort(select_idx_.begin(), select_idx_.end());
}

bool SuccinctBitvector::Access(uint64_t pos) const {
    if (pos >= total_) return false;
    return (blocks_[pos / BLOCK_BITS] >> (pos % BLOCK_BITS)) & 1;
}

uint64_t SuccinctBitvector::Rank1(uint64_t pos) const {
    // Closed-interval rank: count of 1-bits in B[0..pos] inclusive.
    if (pos >= total_) pos = total_ - 1;

    uint64_t sb = pos / SUPERBLOCK_SIZE;
    uint64_t count = superblocks_[sb];

    // Scan words from superblock start to the word containing pos
    uint64_t word_start = (sb * SUPERBLOCK_SIZE) / BLOCK_BITS;
    uint64_t word_end   = pos / BLOCK_BITS;

    for (uint64_t i = word_start; i < word_end; ++i) {
        count += std::popcount(blocks_[i]);
    }

    // Partial word: count bits in [0..bit_offset] of the last word
    uint64_t bit_offset = pos % BLOCK_BITS;
    uint64_t mask = (bit_offset == 63) ? UINT64_MAX : ((1ULL << (bit_offset + 1)) - 1);
    count += std::popcount(blocks_[word_end] & mask);

    return count;
}

uint64_t SuccinctBitvector::Select1(uint64_t k) const {
    // Find position of k-th 1-bit (1-indexed).
    if (k == 0 || k > ones_count_) {
        throw Error{ErrorCode::BUILD_VALIDATION_FAILED,
                    "Select1: k out of range"};
    }
    return select_idx_[k - 1];
}

std::vector<uint8_t> SuccinctBitvector::Serialize() const {
    std::vector<uint8_t> out;
    // total, ones_count, blocks data, superblocks data
    size_t blocks_bytes = blocks_.size() * sizeof(uint64_t);
    size_t sb_bytes = superblocks_.size() * sizeof(uint64_t);
    size_t sel_bytes = select_idx_.size() * sizeof(uint64_t);

    out.resize(sizeof(uint64_t) * 4 + blocks_bytes + sb_bytes + sel_bytes);
    size_t off = 0;

    auto write = [&](const void* src, size_t n) {
        std::memcpy(out.data() + off, src, n);
        off += n;
    };

    write(&total_, sizeof(total_));
    write(&ones_count_, sizeof(ones_count_));
    uint64_t nb = blocks_.size();
    write(&nb, sizeof(nb));
    uint64_t nsb = superblocks_.size();
    write(&nsb, sizeof(nsb));
    if (!blocks_.empty()) write(blocks_.data(), blocks_bytes);
    if (!superblocks_.empty()) write(superblocks_.data(), sb_bytes);
    if (!select_idx_.empty()) {
        uint64_t nsel = select_idx_.size();
        // Append select index
        out.resize(out.size() + sizeof(nsel) + sel_bytes);
        write(&nsel, sizeof(nsel));
        write(select_idx_.data(), sel_bytes);
    }

    return out;
}

void SuccinctBitvector::Deserialize(const uint8_t* data, size_t len) {
    size_t off = 0;
    auto read = [&](void* dst, size_t n) {
        std::memcpy(dst, data + off, n);
        off += n;
    };

    read(&total_, sizeof(total_));
    read(&ones_count_, sizeof(ones_count_));
    uint64_t nb, nsb;
    read(&nb, sizeof(nb));
    read(&nsb, sizeof(nsb));

    blocks_.resize(nb);
    superblocks_.resize(nsb);
    if (nb > 0) read(blocks_.data(), nb * sizeof(uint64_t));
    if (nsb > 0) read(superblocks_.data(), nsb * sizeof(uint64_t));

    // Rebuild select index
    select_idx_.clear();
    for (uint64_t i = 0; i < nb; ++i) {
        uint64_t word = blocks_[i];
        while (word) {
            uint64_t bit_pos = i * BLOCK_BITS + std::countr_zero(word);
            if (bit_pos < total_) {
                select_idx_.push_back(bit_pos);
            }
            word &= word - 1;
        }
    }
    std::sort(select_idx_.begin(), select_idx_.end());
}

// ─────────────────────────────────────────────────────────────────────────────
// BitVectors builder
// ─────────────────────────────────────────────────────────────────────────────

BitVectors BuildBitvectors(const SequenceBundle& bundle, const WindowTable& windows) {
    BitVectors bv;
    uint64_t S_len = bundle.S.size();

    // B_read: 1 at readStarts[i] for all i
    bv.B_read.Build(bundle.readStarts, S_len);

    // B_window: 1 at windows[j].l for all j
    std::vector<uint64_t> window_starts;
    window_starts.reserve(windows.size());
    for (const auto& w : windows) {
        window_starts.push_back(w.l);
    }
    bv.B_window.Build(window_starts, S_len);

    return bv;
}

}  // namespace bamsi
