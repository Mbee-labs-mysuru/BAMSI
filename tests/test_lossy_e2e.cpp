/**
 * test_lossy_e2e.cpp — Lossy-mode end-to-end tests
 *
 * Exec Plan §8.3.4: "End-to-end test of every lossy-mode obligation."
 * Contract §4.3.6: Lossy quality binning verification.
 *
 * Tests:
 * 1. Lossy encode → decode round-trip with per-read directory
 * 2. Lossy encode → decode with block-level directory
 * 3. Verify binned values are within valid range
 * 4. Verify lossy output differs from lossless output
 * 5. Both RANGE_CYCLE and ZSTD_DICT codecs under lossy mode
 */

#include <gtest/gtest.h>
#include "bamsix/types.hpp"
#include "streamencode/streamencode.hpp"

#include <random>
#include <algorithm>

using namespace bamsix;

namespace {

std::vector<OrderedRead> MakeLossyTestReads(size_t n, size_t read_len, uint32_t seed = 42) {
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
            // Use diverse quality values 0..93 to exercise binning
            reads[i].qual[j] = static_cast<uint8_t>(rng() % 94);
        }
    }
    return reads;
}

}  // namespace

// ═══════════════════════════════════════════════════════════════════════════════
// Test 1: Lossy per-read RANGE_CYCLE round-trip
// ═══════════════════════════════════════════════════════════════════════════════

TEST(LossyE2E, PerReadRangeCycleLossyRoundTrip) {
    auto reads = MakeLossyTestReads(20, 150);

    // Encode with lossy_bins=8 (8 quality bins), per-read directory
    auto result = EncodeQualStream(reads, QualCodec::RANGE_CYCLE, 8, 0);
    ASSERT_TRUE(std::holds_alternative<StreamDirectoryPerRead>(result.directory));
    const auto& dir = std::get<StreamDirectoryPerRead>(result.directory);

    for (size_t i = 0; i < reads.size(); ++i) {
        auto decoded = DecodeQualRead(result.payload, dir[i], result.codec_id);
        ASSERT_EQ(decoded.size(), reads[i].qual.size())
            << "Read " << i << " length mismatch";

        // All decoded quality values must be valid (0..93)
        for (size_t j = 0; j < decoded.size(); ++j) {
            EXPECT_LE(decoded[j], 93u)
                << "Read " << i << " qual[" << j << "] out of range: " << (int)decoded[j];
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// Test 2: Lossy block-level RANGE_CYCLE round-trip
// ═══════════════════════════════════════════════════════════════════════════════

TEST(LossyE2E, BlockLevelRangeCycleLossyRoundTrip) {
    auto reads = MakeLossyTestReads(20, 100);

    // Block-level directory, lossy bins=4
    auto result = EncodeQualStream(reads, QualCodec::RANGE_CYCLE, 4, 4);
    ASSERT_TRUE(std::holds_alternative<StreamDirectoryBlockLevel>(result.directory));
    const auto& block_dir = std::get<StreamDirectoryBlockLevel>(result.directory);

    for (size_t i = 0; i < reads.size(); ++i) {
        auto decoded = DecodeQualFromBlock(result.payload, block_dir, 4, i, result.codec_id);
        ASSERT_EQ(decoded.size(), reads[i].qual.size());

        for (size_t j = 0; j < decoded.size(); ++j) {
            EXPECT_LE(decoded[j], 93u);
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// Test 3: Lossy output differs from lossless output
// ═══════════════════════════════════════════════════════════════════════════════

TEST(LossyE2E, LossyDiffersFromLossless) {
    auto reads = MakeLossyTestReads(10, 150);

    // Lossless (bins=0)
    auto lossless = EncodeQualStream(reads, QualCodec::RANGE_CYCLE, 0, 0);
    const auto& lossless_dir = std::get<StreamDirectoryPerRead>(lossless.directory);

    // Lossy (bins=8)
    auto lossy = EncodeQualStream(reads, QualCodec::RANGE_CYCLE, 8, 0);
    const auto& lossy_dir = std::get<StreamDirectoryPerRead>(lossy.directory);

    // At least some reads must have different decoded quality values
    int diff_count = 0;
    for (size_t i = 0; i < reads.size(); ++i) {
        auto dec_lossless = DecodeQualRead(lossless.payload, lossless_dir[i], lossless.codec_id);
        auto dec_lossy = DecodeQualRead(lossy.payload, lossy_dir[i], lossy.codec_id);

        ASSERT_EQ(dec_lossless.size(), dec_lossy.size());

        // Lossless must exactly match original
        ASSERT_EQ(dec_lossless, reads[i].qual)
            << "Lossless decode failed at read " << i;

        // Count differences
        for (size_t j = 0; j < dec_lossless.size(); ++j) {
            if (dec_lossy[j] != dec_lossless[j]) {
                diff_count++;
            }
        }
    }

    // With diverse quality values and 8 bins, there MUST be differences
    EXPECT_GT(diff_count, 0)
        << "Lossy output should differ from lossless for diverse quality inputs";
}

// ═══════════════════════════════════════════════════════════════════════════════
// Test 4: Lossy binning produces consistent bins
// ═══════════════════════════════════════════════════════════════════════════════

TEST(LossyE2E, LossyBinningConsistency) {
    auto reads = MakeLossyTestReads(5, 200);

    // Two independent encodes with same lossy_bins must produce identical output
    auto enc1 = EncodeQualStream(reads, QualCodec::RANGE_CYCLE, 8, 0);
    auto enc2 = EncodeQualStream(reads, QualCodec::RANGE_CYCLE, 8, 0);
    const auto& dir1 = std::get<StreamDirectoryPerRead>(enc1.directory);
    const auto& dir2 = std::get<StreamDirectoryPerRead>(enc2.directory);

    for (size_t i = 0; i < reads.size(); ++i) {
        auto dec1 = DecodeQualRead(enc1.payload, dir1[i], enc1.codec_id);
        auto dec2 = DecodeQualRead(enc2.payload, dir2[i], enc2.codec_id);
        ASSERT_EQ(dec1, dec2) << "Inconsistent lossy binning at read " << i;
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// Test 5: ZSTD_DICT codec under lossy mode
// ═══════════════════════════════════════════════════════════════════════════════

TEST(LossyE2E, ZstdDictLossyRoundTrip) {
    auto reads = MakeLossyTestReads(15, 120);

    // ZSTD_DICT codec with lossy_bins=4
    auto result = EncodeQualStream(reads, QualCodec::ZSTD_DICT, 4, 0);
    ASSERT_TRUE(std::holds_alternative<StreamDirectoryPerRead>(result.directory));
    const auto& dir = std::get<StreamDirectoryPerRead>(result.directory);

    for (size_t i = 0; i < reads.size(); ++i) {
        auto decoded = DecodeQualRead(result.payload, dir[i], result.codec_id);
        ASSERT_EQ(decoded.size(), reads[i].qual.size());

        for (size_t j = 0; j < decoded.size(); ++j) {
            EXPECT_LE(decoded[j], 93u);
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// Test 6: ZSTD_DICT lossless round-trip
// ═══════════════════════════════════════════════════════════════════════════════

TEST(LossyE2E, ZstdDictLosslessRoundTrip) {
    auto reads = MakeLossyTestReads(10, 100);

    // ZSTD_DICT codec with lossy_bins=0 (lossless)
    auto result = EncodeQualStream(reads, QualCodec::ZSTD_DICT, 0, 0);
    const auto& dir = std::get<StreamDirectoryPerRead>(result.directory);

    for (size_t i = 0; i < reads.size(); ++i) {
        auto decoded = DecodeQualRead(result.payload, dir[i], result.codec_id);
        ASSERT_EQ(decoded, reads[i].qual)
            << "ZSTD_DICT lossless round-trip failed at read " << i;
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// Test 7: ZSTD_DICT block-level round-trip
// ═══════════════════════════════════════════════════════════════════════════════

TEST(LossyE2E, ZstdDictBlockLevelRoundTrip) {
    auto reads = MakeLossyTestReads(12, 80);

    // ZSTD_DICT + block-level directory
    auto result = EncodeQualStream(reads, QualCodec::ZSTD_DICT, 0, 4);
    ASSERT_TRUE(std::holds_alternative<StreamDirectoryBlockLevel>(result.directory));
    const auto& block_dir = std::get<StreamDirectoryBlockLevel>(result.directory);

    for (size_t i = 0; i < reads.size(); ++i) {
        auto decoded = DecodeQualFromBlock(result.payload, block_dir, 4, i, result.codec_id);
        ASSERT_EQ(decoded, reads[i].qual)
            << "ZSTD_DICT block-level round-trip failed at read " << i;
    }
}
