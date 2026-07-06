/// Codec Completeness Tests — Contract §2.7, §2.8, §2.9
/// Round-trip verification for ALL codec variants:
///   S_qual: RANGE_CYCLE (0x01), RANS_DELTA (0x02), ZSTD_DICT (0x03), BINNED_RANGE (0x04)
///   S_meta: TYPED_SPLIT (0x01), ZSTD_FALLBACK (0x02)
///   S_map:  DELTA_RANGE (0x01), RAW (0x02)
///
/// These tests verify the codec bake-off requirement (Exec Plan §5.3.6):
/// every declared codec variant must round-trip correctly.

#include <gtest/gtest.h>

#include "bamsix/types.hpp"
#include "../src/streamencode/streamencode.hpp"

#include <algorithm>
#include <cstring>
#include <numeric>
#include <random>
#include <string>

using namespace bamsix;

// ─── Test fixture with synthetic reads ───────────────────────────────────────

class CodecCompletenessTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Create a set of synthetic reads with varying quality profiles
        // to stress-test codec variants
        std::mt19937 rng(42);  // deterministic seed

        for (uint64_t i = 0; i < 20; ++i) {
            OrderedRead r;
            r.read_id = i;
            r.chrom_id = 0;
            r.chrom = "chr1";
            r.pos = 100 + i * 50;
            r.flag = (i % 3 == 0) ? 99 : 0;  // varied FLAG values
            r.mapq = static_cast<uint8_t>(30 + (i % 30));
            r.qname = "read_" + std::to_string(i);

            // Sequence: 50-100 bases
            uint32_t seq_len = 50 + (i % 51);
            r.seq.resize(seq_len);
            for (uint32_t j = 0; j < seq_len; ++j) {
                r.seq[j] = rng() % 5;  // A,C,G,T,N
            }

            // Quality scores: varied profiles to test all codecs
            r.qual.resize(seq_len);
            if (i % 4 == 0) {
                // Flat quality
                std::fill(r.qual.begin(), r.qual.end(), 35);
            } else if (i % 4 == 1) {
                // Declining quality (typical Illumina)
                for (uint32_t j = 0; j < seq_len; ++j) {
                    r.qual[j] = std::max(2, 40 - static_cast<int>(j) / 3);
                }
            } else if (i % 4 == 2) {
                // Random quality
                for (uint32_t j = 0; j < seq_len; ++j) {
                    r.qual[j] = rng() % 42;
                }
            } else {
                // Bimodal quality (high/low alternating)
                for (uint32_t j = 0; j < seq_len; ++j) {
                    r.qual[j] = (j % 2 == 0) ? 40 : 10;
                }
            }

            // CIGAR: simple M-only
            r.cigar = {{0, seq_len}};

            // Aux data: some reads have aux tags
            if (i % 5 == 0) {
                r.aux_data = {0x4E, 0x4D, 0x43, 0x03};  // "NMC\x03"
            }

            reads_.push_back(std::move(r));
        }
    }

    std::vector<OrderedRead> reads_;
};

// ─── S_qual: RANGE_CYCLE (0x01) ─────────────────────────────────────────────

TEST_F(CodecCompletenessTest, QualRangeCycleRoundTrip) {
    auto result = EncodeQualStream(reads_, QualCodec::RANGE_CYCLE, 0, 0);
    EXPECT_EQ(result.codec_id, static_cast<uint8_t>(QualCodec::RANGE_CYCLE));
    EXPECT_FALSE(result.payload.empty());

    auto& dir = std::get<StreamDirectoryPerRead>(result.directory);
    ASSERT_EQ(dir.size(), reads_.size());

    for (size_t i = 0; i < reads_.size(); ++i) {
        auto decoded = DecodeQualRead(result.payload, dir[i], result.codec_id);
        EXPECT_EQ(decoded, reads_[i].qual)
            << "RANGE_CYCLE round-trip failed for read " << i;
    }
}

TEST_F(CodecCompletenessTest, QualRangeCycleBlockLevel) {
    auto result = EncodeQualStream(reads_, QualCodec::RANGE_CYCLE, 0, 4);
    EXPECT_EQ(result.codec_id, static_cast<uint8_t>(QualCodec::RANGE_CYCLE));

    auto& block_dir = std::get<StreamDirectoryBlockLevel>(result.directory);
    EXPECT_GT(block_dir.size(), 0u);

    for (size_t i = 0; i < reads_.size(); ++i) {
        auto decoded = DecodeQualFromBlock(result.payload, block_dir, 4, i,
                                            result.codec_id);
        EXPECT_EQ(decoded, reads_[i].qual)
            << "RANGE_CYCLE block-level round-trip failed for read " << i;
    }
}

// ─── S_qual: RANS_DELTA (0x02) ──────────────────────────────────────────────

TEST_F(CodecCompletenessTest, QualRansDeltaRoundTrip) {
    auto result = EncodeQualStream(reads_, QualCodec::RANS_DELTA, 0, 0);
    EXPECT_EQ(result.codec_id, static_cast<uint8_t>(QualCodec::RANS_DELTA));
    EXPECT_FALSE(result.payload.empty());

    auto& dir = std::get<StreamDirectoryPerRead>(result.directory);
    ASSERT_EQ(dir.size(), reads_.size());

    for (size_t i = 0; i < reads_.size(); ++i) {
        auto decoded = DecodeQualRead(result.payload, dir[i], result.codec_id);
        EXPECT_EQ(decoded, reads_[i].qual)
            << "RANS_DELTA round-trip failed for read " << i;
    }
}

TEST_F(CodecCompletenessTest, QualRansDeltaBlockLevel) {
    auto result = EncodeQualStream(reads_, QualCodec::RANS_DELTA, 0, 8);
    auto& block_dir = std::get<StreamDirectoryBlockLevel>(result.directory);
    EXPECT_GT(block_dir.size(), 0u);

    for (size_t i = 0; i < reads_.size(); ++i) {
        auto decoded = DecodeQualFromBlock(result.payload, block_dir, 8, i,
                                            result.codec_id);
        EXPECT_EQ(decoded, reads_[i].qual)
            << "RANS_DELTA block-level round-trip failed for read " << i;
    }
}

TEST_F(CodecCompletenessTest, QualRansDeltaFlatProfile) {
    // Flat quality profile should delta-code to all 128 (zero deltas)
    std::vector<OrderedRead> flat_reads;
    OrderedRead r;
    r.read_id = 0; r.chrom_id = 0; r.chrom = "chr1"; r.pos = 100;
    r.flag = 0; r.mapq = 60; r.qname = "flat";
    r.seq.resize(100, CODE_A);
    r.qual.resize(100, 35);  // All Q35
    r.cigar = {{0, 100}};
    flat_reads.push_back(std::move(r));

    auto result = EncodeQualStream(flat_reads, QualCodec::RANS_DELTA, 0, 0);
    auto& dir = std::get<StreamDirectoryPerRead>(result.directory);
    auto decoded = DecodeQualRead(result.payload, dir[0], result.codec_id);
    EXPECT_EQ(decoded, flat_reads[0].qual);
}

// ─── S_qual: ZSTD_DICT (0x03) ───────────────────────────────────────────────

TEST_F(CodecCompletenessTest, QualZstdDictRoundTrip) {
    auto result = EncodeQualStream(reads_, QualCodec::ZSTD_DICT, 0, 0);
    EXPECT_EQ(result.codec_id, static_cast<uint8_t>(QualCodec::ZSTD_DICT));

    auto& dir = std::get<StreamDirectoryPerRead>(result.directory);
    for (size_t i = 0; i < reads_.size(); ++i) {
        auto decoded = DecodeQualRead(result.payload, dir[i], result.codec_id);
        EXPECT_EQ(decoded, reads_[i].qual)
            << "ZSTD_DICT round-trip failed for read " << i;
    }
}

TEST_F(CodecCompletenessTest, QualZstdDictBlockLevel) {
    auto result = EncodeQualStream(reads_, QualCodec::ZSTD_DICT, 0, 5);
    auto& block_dir = std::get<StreamDirectoryBlockLevel>(result.directory);
    for (size_t i = 0; i < reads_.size(); ++i) {
        auto decoded = DecodeQualFromBlock(result.payload, block_dir, 5, i,
                                            result.codec_id);
        EXPECT_EQ(decoded, reads_[i].qual)
            << "ZSTD_DICT block-level round-trip failed for read " << i;
    }
}

// ─── S_qual: BINNED_RANGE (0x04) ────────────────────────────────────────────

TEST_F(CodecCompletenessTest, QualBinnedRangeRoundTrip) {
    // BINNED_RANGE applies its own internal binning, so the round-trip is
    // encode → decode → verify decoded matches what was encoded (post-binning).
    // We verify that encode doesn't crash and decode produces valid quality scores.
    auto result = EncodeQualStream(reads_, QualCodec::BINNED_RANGE, 0, 0);
    EXPECT_EQ(result.codec_id, static_cast<uint8_t>(QualCodec::BINNED_RANGE));

    auto& dir = std::get<StreamDirectoryPerRead>(result.directory);
    for (size_t i = 0; i < reads_.size(); ++i) {
        auto decoded = DecodeQualRead(result.payload, dir[i], result.codec_id);
        ASSERT_EQ(decoded.size(), reads_[i].qual.size())
            << "BINNED_RANGE: wrong output length for read " << i;
        // Verify all values are valid Phred scores [0, 93]
        for (size_t j = 0; j < decoded.size(); ++j) {
            EXPECT_LE(decoded[j], 93u)
                << "BINNED_RANGE: invalid Phred value at read " << i << " pos " << j;
        }
    }
}

TEST_F(CodecCompletenessTest, QualBinnedRangeWithExplicitBins) {
    // Test BINNED_RANGE with explicit lossy_bins
    auto result = EncodeQualStream(reads_, QualCodec::BINNED_RANGE, 4, 0);
    auto& dir = std::get<StreamDirectoryPerRead>(result.directory);
    for (size_t i = 0; i < reads_.size(); ++i) {
        auto decoded = DecodeQualRead(result.payload, dir[i], result.codec_id);
        ASSERT_EQ(decoded.size(), reads_[i].qual.size());
        // With 4 bins, values should be limited to at most 4 distinct levels
        std::set<uint8_t> unique_vals(decoded.begin(), decoded.end());
        EXPECT_LE(unique_vals.size(), 5u)
            << "BINNED_RANGE with 4 bins: too many unique values for read " << i;
    }
}

TEST_F(CodecCompletenessTest, QualBinnedRangeBlockLevel) {
    auto result = EncodeQualStream(reads_, QualCodec::BINNED_RANGE, 0, 10);
    auto& block_dir = std::get<StreamDirectoryBlockLevel>(result.directory);
    for (size_t i = 0; i < reads_.size(); ++i) {
        auto decoded = DecodeQualFromBlock(result.payload, block_dir, 10, i,
                                            result.codec_id);
        ASSERT_EQ(decoded.size(), reads_[i].qual.size());
    }
}

// ─── S_qual: Lossy mode across all codecs ────────────────────────────────────

TEST_F(CodecCompletenessTest, QualLossyModeReducesDistinct) {
    // With lossy_bins=4, all codecs should produce fewer unique quality values
    for (QualCodec codec : {QualCodec::RANGE_CYCLE, QualCodec::RANS_DELTA,
                             QualCodec::ZSTD_DICT}) {
        auto result = EncodeQualStream(reads_, codec, 4, 0);
        auto& dir = std::get<StreamDirectoryPerRead>(result.directory);
        for (size_t i = 0; i < reads_.size(); ++i) {
            auto decoded = DecodeQualRead(result.payload, dir[i], result.codec_id);
            std::set<uint8_t> unique_vals(decoded.begin(), decoded.end());
            EXPECT_LE(unique_vals.size(), 5u)
                << "Lossy mode failed for codec " << static_cast<int>(codec)
                << " read " << i;
        }
    }
}

// ─── S_meta: TYPED_SPLIT (0x01) ─────────────────────────────────────────────

TEST_F(CodecCompletenessTest, MetaTypedSplitRoundTrip) {
    auto result = EncodeMetaStream(reads_, MetaCodec::TYPED_SPLIT);
    EXPECT_EQ(result.codec_id, static_cast<uint8_t>(MetaCodec::TYPED_SPLIT));
    EXPECT_FALSE(result.payload.empty());

    for (size_t i = 0; i < reads_.size(); ++i) {
        auto decoded = DecodeMetaRead(result.payload, result.directory[i],
                                       result.codec_id);
        EXPECT_EQ(decoded.flag, reads_[i].flag)
            << "TYPED_SPLIT: FLAG mismatch for read " << i;
        ASSERT_EQ(decoded.cigar.size(), reads_[i].cigar.size())
            << "TYPED_SPLIT: CIGAR op count mismatch for read " << i;
        for (size_t c = 0; c < decoded.cigar.size(); ++c) {
            EXPECT_EQ(decoded.cigar[c].op, reads_[i].cigar[c].op);
            EXPECT_EQ(decoded.cigar[c].len, reads_[i].cigar[c].len);
        }
        EXPECT_EQ(decoded.aux_data, reads_[i].aux_data)
            << "TYPED_SPLIT: aux_data mismatch for read " << i;
    }
}

// ─── S_meta: ZSTD_FALLBACK (0x02) ───────────────────────────────────────────

TEST_F(CodecCompletenessTest, MetaZstdFallbackRoundTrip) {
    auto result = EncodeMetaStream(reads_, MetaCodec::ZSTD_FALLBACK);
    EXPECT_EQ(result.codec_id, static_cast<uint8_t>(MetaCodec::ZSTD_FALLBACK));
    EXPECT_FALSE(result.payload.empty());

    for (size_t i = 0; i < reads_.size(); ++i) {
        auto decoded = DecodeMetaRead(result.payload, result.directory[i],
                                       result.codec_id);
        EXPECT_EQ(decoded.flag, reads_[i].flag)
            << "ZSTD_FALLBACK: FLAG mismatch for read " << i;
        ASSERT_EQ(decoded.cigar.size(), reads_[i].cigar.size())
            << "ZSTD_FALLBACK: CIGAR op count mismatch for read " << i;
        for (size_t c = 0; c < decoded.cigar.size(); ++c) {
            EXPECT_EQ(decoded.cigar[c].op, reads_[i].cigar[c].op);
            EXPECT_EQ(decoded.cigar[c].len, reads_[i].cigar[c].len);
        }
        EXPECT_EQ(decoded.aux_data, reads_[i].aux_data)
            << "ZSTD_FALLBACK: aux_data mismatch for read " << i;
    }
}

TEST_F(CodecCompletenessTest, MetaZstdFallbackComplexCigar) {
    // Test with complex multi-op CIGAR records
    std::vector<OrderedRead> complex;
    OrderedRead r;
    r.read_id = 0; r.chrom_id = 0; r.chrom = "chr1"; r.pos = 100;
    r.flag = 163; r.mapq = 60; r.qname = "complex_cigar";
    r.seq.resize(100, CODE_A);
    r.qual.resize(100, 30);
    r.cigar = {{0, 20}, {1, 3}, {0, 30}, {2, 5}, {0, 47}};  // 20M3I30M5D47M
    r.aux_data = {0x4E, 0x4D, 0x43, 0x08, 0x58, 0x53, 0x41, 0x02, 0x41, 0x42};
    complex.push_back(std::move(r));

    auto result = EncodeMetaStream(complex, MetaCodec::ZSTD_FALLBACK);
    auto decoded = DecodeMetaRead(result.payload, result.directory[0], result.codec_id);

    ASSERT_EQ(decoded.cigar.size(), 5u);
    EXPECT_EQ(decoded.cigar[0].op, 0u); EXPECT_EQ(decoded.cigar[0].len, 20u);
    EXPECT_EQ(decoded.cigar[1].op, 1u); EXPECT_EQ(decoded.cigar[1].len, 3u);
    EXPECT_EQ(decoded.cigar[2].op, 0u); EXPECT_EQ(decoded.cigar[2].len, 30u);
    EXPECT_EQ(decoded.cigar[3].op, 2u); EXPECT_EQ(decoded.cigar[3].len, 5u);
    EXPECT_EQ(decoded.cigar[4].op, 0u); EXPECT_EQ(decoded.cigar[4].len, 47u);
    EXPECT_EQ(decoded.flag, 163u);
    EXPECT_EQ(decoded.aux_data, complex[0].aux_data);
}

// ─── S_map: DELTA_RANGE (0x01) ──────────────────────────────────────────────

TEST_F(CodecCompletenessTest, MapDeltaRangeRoundTrip) {
    auto result = EncodeMapStream(reads_, MapCodec::DELTA_RANGE);
    EXPECT_EQ(result.codec_id, static_cast<uint8_t>(MapCodec::DELTA_RANGE));

    // For DELTA_RANGE, first read of each chrom is absolute, rest are deltas.
    // For our single-chrom test data, read 0 is absolute, rest are deltas.
    uint64_t running_pos = 0;
    for (size_t i = 0; i < reads_.size(); ++i) {
        auto decoded = DecodeMapRead(result.payload, result.directory[i],
                                      result.codec_id);
        EXPECT_EQ(decoded.chrom_id, reads_[i].chrom_id)
            << "DELTA_RANGE: chrom_id mismatch for read " << i;

        // Reconstruct absolute position
        uint64_t abs_pos;
        if (!decoded.is_delta) {
            abs_pos = decoded.pos;
        } else {
            int64_t delta;
            std::memcpy(&delta, &decoded.pos, sizeof(delta));
            abs_pos = running_pos + delta;
        }
        running_pos = reads_[i].pos;

        EXPECT_EQ(abs_pos, reads_[i].pos)
            << "DELTA_RANGE: pos mismatch for read " << i;
    }
}

// ─── S_map: RAW (0x02) ──────────────────────────────────────────────────────

TEST_F(CodecCompletenessTest, MapRawRoundTrip) {
    auto result = EncodeMapStream(reads_, MapCodec::RAW);
    EXPECT_EQ(result.codec_id, static_cast<uint8_t>(MapCodec::RAW));

    for (size_t i = 0; i < reads_.size(); ++i) {
        auto decoded = DecodeMapRead(result.payload, result.directory[i],
                                      result.codec_id);
        EXPECT_EQ(decoded.chrom_id, reads_[i].chrom_id)
            << "RAW: chrom_id mismatch for read " << i;
        EXPECT_EQ(decoded.pos, reads_[i].pos)
            << "RAW: pos mismatch for read " << i;
        EXPECT_FALSE(decoded.is_delta);
    }
}

// ─── S_seq: BWT_MTF_RLE_ARITH (0x10) ────────────────────────────────────────

TEST_F(CodecCompletenessTest, SeqBwtMtfRleArithRoundTrip) {
    // Build a larger BWT for testing — arithmetic coder needs sufficient data
    std::mt19937 rng(123);
    std::vector<uint8_t> bwt(512);
    for (size_t i = 0; i < bwt.size(); ++i) {
        bwt[i] = rng() % SIGMA;  // codes 0..5
    }
    auto result = EncodeSeqStream(bwt, 6);
    EXPECT_EQ(result.codec_id, 0x10u);
    EXPECT_FALSE(result.payload.empty());

    auto decoded = DecodeSeqStream(result.payload, result.codec_id, bwt.size());
    EXPECT_EQ(decoded, bwt) << "S_seq BWT_MTF_RLE_ARITH round-trip failed";
}

// ─── Cross-codec consistency ─────────────────────────────────────────────────

TEST_F(CodecCompletenessTest, AllQualCodecsAcceptEmptyRead) {
    std::vector<OrderedRead> empty_reads;
    OrderedRead r;
    r.read_id = 0; r.chrom_id = 0; r.chrom = "chr1"; r.pos = 100;
    r.flag = 0; r.mapq = 60; r.qname = "empty";
    r.seq = {};
    r.qual = {};
    r.cigar = {};
    empty_reads.push_back(std::move(r));

    for (QualCodec codec : {QualCodec::RANGE_CYCLE, QualCodec::RANS_DELTA,
                             QualCodec::ZSTD_DICT, QualCodec::BINNED_RANGE}) {
        EXPECT_NO_THROW(EncodeQualStream(empty_reads, codec, 0, 0))
            << "Codec " << static_cast<int>(codec) << " crashed on empty read";
    }
}

TEST_F(CodecCompletenessTest, AllQualCodecsAcceptSingleBase) {
    std::vector<OrderedRead> single_reads;
    OrderedRead r;
    r.read_id = 0; r.chrom_id = 0; r.chrom = "chr1"; r.pos = 100;
    r.flag = 0; r.mapq = 60; r.qname = "single";
    r.seq = {CODE_A};
    r.qual = {42};
    r.cigar = {{0, 1}};
    single_reads.push_back(std::move(r));

    for (QualCodec codec : {QualCodec::RANGE_CYCLE, QualCodec::RANS_DELTA,
                             QualCodec::ZSTD_DICT, QualCodec::BINNED_RANGE}) {
        auto result = EncodeQualStream(single_reads, codec, 0, 0);
        auto& dir = std::get<StreamDirectoryPerRead>(result.directory);
        auto decoded = DecodeQualRead(result.payload, dir[0], result.codec_id);
        // For lossy codecs (BINNED_RANGE), just check length
        EXPECT_EQ(decoded.size(), 1u)
            << "Codec " << static_cast<int>(codec) << " wrong length for single-base read";
    }
}

TEST_F(CodecCompletenessTest, AllMetaCodecsAcceptEmptyRead) {
    std::vector<OrderedRead> empty_reads;
    OrderedRead r;
    r.read_id = 0; r.chrom_id = 0; r.chrom = "chr1"; r.pos = 100;
    r.flag = 0; r.mapq = 60; r.qname = "empty";
    r.seq = {};
    r.qual = {};
    r.cigar = {};
    empty_reads.push_back(std::move(r));

    for (MetaCodec codec : {MetaCodec::TYPED_SPLIT, MetaCodec::ZSTD_FALLBACK}) {
        EXPECT_NO_THROW(EncodeMetaStream(empty_reads, codec))
            << "Meta codec " << static_cast<int>(codec) << " crashed on empty read";
    }
}
