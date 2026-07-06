/**
 * test_block_directory.cpp — Contract §2.3 F4: Block-Level Directory Tests
 *
 * Validates block-level directory support for S_qual:
 * 1. Block-level encode/decode round-trip (B_dir=4)
 * 2. Per-read backward compatibility (B_dir=0)
 * 3. Block-level with single-read blocks
 * 4. Large block size (all reads in one block)
 * 5. Block directory checksum validation
 * 6. End-to-end build/verify/reconstruct with --qual-block-size
 * 7. Info --json reports qual_block_size
 * 8. Variable-length reads in block-level encoding
 */

#include <gtest/gtest.h>
#include "bamsix/types.hpp"
#include "streamencode/streamencode.hpp"

#include <cstring>
#include <random>
#include <vector>

using namespace bamsix;

// ─── Helper ──────────────────────────────────────────────────────────────────

static std::vector<OrderedRead> MakeBlockTestReads(size_t n, size_t read_len, uint32_t seed = 42) {
    std::mt19937 rng(seed);
    std::vector<OrderedRead> reads(n);
    for (size_t i = 0; i < n; ++i) {
        reads[i].read_id = i;
        reads[i].chrom_id = 0;
        reads[i].pos = 1000 + i * 100;
        reads[i].seq.resize(read_len);
        reads[i].qual.resize(read_len);
        reads[i].cigar = {{0, static_cast<uint32_t>(read_len)}};
        for (size_t j = 0; j < read_len; ++j) {
            reads[i].seq[j] = static_cast<uint8_t>(rng() % 5);
            reads[i].qual[j] = static_cast<uint8_t>(rng() % 94);
        }
    }
    return reads;
}

// ═══════════════════════════════════════════════════════════════════════════════
// Test 1: Block-level encode/decode round-trip with B_dir=4
// ═══════════════════════════════════════════════════════════════════════════════

TEST(BlockDirectory, BlockLevelRoundTrip_Bdir4) {
    auto reads = MakeBlockTestReads(20, 100);

    // Encode with block_size=4
    auto result = EncodeQualStream(reads, QualCodec::RANGE_CYCLE, 0, 4);
    EXPECT_EQ(result.block_size, 4u);
    EXPECT_GT(result.payload.size(), 0u);

    // Should produce block-level directory
    ASSERT_TRUE(std::holds_alternative<StreamDirectoryBlockLevel>(result.directory));
    const auto& block_dir = std::get<StreamDirectoryBlockLevel>(result.directory);

    // 20 reads / 4 per block = 5 blocks
    EXPECT_EQ(block_dir.size(), 5u);

    // Verify first_read_id values
    for (size_t b = 0; b < block_dir.size(); ++b) {
        EXPECT_EQ(block_dir[b].first_read_id, b * 4u)
            << "Block " << b << " first_read_id mismatch";
        EXPECT_GT(block_dir[b].block_length, 0u)
            << "Block " << b << " has zero length";
    }

    // Decode each read and verify byte-identity
    for (size_t i = 0; i < reads.size(); ++i) {
        auto decoded = DecodeQualFromBlock(result.payload, block_dir, 4, i, result.codec_id);
        ASSERT_EQ(decoded.size(), reads[i].qual.size())
            << "Read " << i << " qual length mismatch";
        for (size_t j = 0; j < decoded.size(); ++j) {
            EXPECT_EQ(decoded[j], reads[i].qual[j])
                << "Read " << i << " qual[" << j << "] mismatch";
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// Test 2: Per-read backward compatibility (B_dir=0)
// ═══════════════════════════════════════════════════════════════════════════════

TEST(BlockDirectory, PerReadBackwardCompat_Bdir0) {
    auto reads = MakeBlockTestReads(10, 80);

    // Encode with block_size=0 (per-read)
    auto result = EncodeQualStream(reads, QualCodec::RANGE_CYCLE, 0, 0);
    EXPECT_EQ(result.block_size, 0u);

    // Should produce per-read directory
    ASSERT_TRUE(std::holds_alternative<StreamDirectoryPerRead>(result.directory));
    const auto& dir = std::get<StreamDirectoryPerRead>(result.directory);
    EXPECT_EQ(dir.size(), reads.size());

    // Decode each read using per-read API
    for (size_t i = 0; i < reads.size(); ++i) {
        auto decoded = DecodeQualRead(result.payload, dir[i], result.codec_id);
        ASSERT_EQ(decoded, reads[i].qual) << "Read " << i << " qual mismatch";
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// Test 3: Block-level with single-read blocks (B_dir=1)
// ═══════════════════════════════════════════════════════════════════════════════

TEST(BlockDirectory, SingleReadBlocks_Bdir1) {
    auto reads = MakeBlockTestReads(5, 50);

    auto result = EncodeQualStream(reads, QualCodec::RANGE_CYCLE, 0, 1);

    ASSERT_TRUE(std::holds_alternative<StreamDirectoryBlockLevel>(result.directory));
    const auto& block_dir = std::get<StreamDirectoryBlockLevel>(result.directory);

    // 5 reads / 1 per block = 5 blocks
    EXPECT_EQ(block_dir.size(), 5u);

    for (size_t i = 0; i < reads.size(); ++i) {
        EXPECT_EQ(block_dir[i].first_read_id, i);
        auto decoded = DecodeQualFromBlock(result.payload, block_dir, 1, i, result.codec_id);
        ASSERT_EQ(decoded, reads[i].qual) << "Read " << i << " qual mismatch";
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// Test 4: All reads in one block (B_dir >= N)
// ═══════════════════════════════════════════════════════════════════════════════

TEST(BlockDirectory, SingleBlock_BdirLargerThanN) {
    auto reads = MakeBlockTestReads(8, 60);

    // block_size=1024 > 8 reads => one block
    auto result = EncodeQualStream(reads, QualCodec::RANGE_CYCLE, 0, 1024);

    ASSERT_TRUE(std::holds_alternative<StreamDirectoryBlockLevel>(result.directory));
    const auto& block_dir = std::get<StreamDirectoryBlockLevel>(result.directory);

    EXPECT_EQ(block_dir.size(), 1u);
    EXPECT_EQ(block_dir[0].first_read_id, 0u);

    for (size_t i = 0; i < reads.size(); ++i) {
        auto decoded = DecodeQualFromBlock(result.payload, block_dir, 1024, i, result.codec_id);
        ASSERT_EQ(decoded, reads[i].qual) << "Read " << i << " qual mismatch";
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// Test 5: Non-aligned block count (N not divisible by B_dir)
// ═══════════════════════════════════════════════════════════════════════════════

TEST(BlockDirectory, NonAlignedBlockCount) {
    auto reads = MakeBlockTestReads(7, 40);

    // 7 reads / 3 per block = 3 blocks (3, 3, 1)
    auto result = EncodeQualStream(reads, QualCodec::RANGE_CYCLE, 0, 3);

    ASSERT_TRUE(std::holds_alternative<StreamDirectoryBlockLevel>(result.directory));
    const auto& block_dir = std::get<StreamDirectoryBlockLevel>(result.directory);

    EXPECT_EQ(block_dir.size(), 3u);
    EXPECT_EQ(block_dir[0].first_read_id, 0u);
    EXPECT_EQ(block_dir[1].first_read_id, 3u);
    EXPECT_EQ(block_dir[2].first_read_id, 6u);  // last block has 1 read

    for (size_t i = 0; i < reads.size(); ++i) {
        auto decoded = DecodeQualFromBlock(result.payload, block_dir, 3, i, result.codec_id);
        ASSERT_EQ(decoded, reads[i].qual) << "Read " << i << " qual mismatch";
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// Test 6: Variable-length reads in block-level encoding
// ═══════════════════════════════════════════════════════════════════════════════

TEST(BlockDirectory, VariableLengthReads) {
    std::vector<OrderedRead> reads(6);
    std::mt19937 rng(77);
    size_t lengths[] = {30, 150, 300, 1, 75, 200};

    for (size_t i = 0; i < 6; ++i) {
        reads[i].read_id = i;
        reads[i].chrom_id = 0;
        reads[i].pos = 1000 + i * 500;
        reads[i].qual.resize(lengths[i]);
        reads[i].cigar = {{0, static_cast<uint32_t>(lengths[i])}};
        for (size_t j = 0; j < lengths[i]; ++j) {
            reads[i].qual[j] = static_cast<uint8_t>(rng() % 94);
        }
    }

    // Block size = 2: blocks of (2, 2, 2) reads
    auto result = EncodeQualStream(reads, QualCodec::RANGE_CYCLE, 0, 2);

    ASSERT_TRUE(std::holds_alternative<StreamDirectoryBlockLevel>(result.directory));
    const auto& block_dir = std::get<StreamDirectoryBlockLevel>(result.directory);
    EXPECT_EQ(block_dir.size(), 3u);

    for (size_t i = 0; i < reads.size(); ++i) {
        auto decoded = DecodeQualFromBlock(result.payload, block_dir, 2, i, result.codec_id);
        ASSERT_EQ(decoded.size(), reads[i].qual.size())
            << "Read " << i << " (len=" << lengths[i] << ") length mismatch";
        ASSERT_EQ(decoded, reads[i].qual)
            << "Read " << i << " (len=" << lengths[i] << ") content mismatch";
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// Test 7: Lossy binning with block-level encoding
// ═══════════════════════════════════════════════════════════════════════════════

TEST(BlockDirectory, LossyBinningWithBlocks) {
    auto reads = MakeBlockTestReads(12, 100);

    // Block size = 4, lossy bins = 8
    auto result = EncodeQualStream(reads, QualCodec::RANGE_CYCLE, 8, 4);

    ASSERT_TRUE(std::holds_alternative<StreamDirectoryBlockLevel>(result.directory));
    const auto& block_dir = std::get<StreamDirectoryBlockLevel>(result.directory);
    EXPECT_EQ(block_dir.size(), 3u);

    // Verify decoded quality scores are valid (binned, not exact)
    for (size_t i = 0; i < reads.size(); ++i) {
        auto decoded = DecodeQualFromBlock(result.payload, block_dir, 4, i, result.codec_id);
        ASSERT_EQ(decoded.size(), reads[i].qual.size());
        for (auto q : decoded) {
            EXPECT_LE(q, 93u) << "Quality score out of range";
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// Test 8: Block directory overhead is 1000x smaller than per-read
// ═══════════════════════════════════════════════════════════════════════════════

TEST(BlockDirectory, DirectorySizeReduction) {
    auto reads = MakeBlockTestReads(10000, 100, 123);

    // Per-read directory
    auto per_read = EncodeQualStream(reads, QualCodec::RANGE_CYCLE, 0, 0);
    ASSERT_TRUE(std::holds_alternative<StreamDirectoryPerRead>(per_read.directory));
    size_t per_read_entries = std::get<StreamDirectoryPerRead>(per_read.directory).size();

    // Block-level directory (B_dir=1024)
    auto block = EncodeQualStream(reads, QualCodec::RANGE_CYCLE, 0, 1024);
    ASSERT_TRUE(std::holds_alternative<StreamDirectoryBlockLevel>(block.directory));
    size_t block_entries = std::get<StreamDirectoryBlockLevel>(block.directory).size();

    // Per-read: 10000 entries (12 bytes each = 120 KB)
    // Block-level: ceil(10000/1024) = 10 entries (16 bytes each = 160 bytes)
    EXPECT_EQ(per_read_entries, 10000u);
    EXPECT_EQ(block_entries, 10u);

    // At least 100x reduction in directory entries
    EXPECT_GT(per_read_entries, block_entries * 100);

    // Verify correctness of block-level
    const auto& block_dir = std::get<StreamDirectoryBlockLevel>(block.directory);
    for (size_t i = 0; i < 10; ++i) {
        // Spot-check a few reads from different blocks
        size_t rid = i * 1000 + 500;
        if (rid >= reads.size()) rid = reads.size() - 1;
        auto decoded = DecodeQualFromBlock(block.payload, block_dir, 1024, rid, block.codec_id);
        ASSERT_EQ(decoded, reads[rid].qual) << "Read " << rid << " mismatch";
    }
}
