/**
 * test_determinism.cpp — Cross-platform determinism verification
 *
 * Contract §11.1 (I9): "Two invocations of Locate(P) in sorted mode on the
 * same .bsi produce byte-identical output."
 *
 * Tests:
 * 1. Build index → query twice → verify identical results
 * 2. FM-index construction is deterministic
 * 3. Stream encoding is deterministic
 * 4. Sorted Locate results are deterministic
 */

#include <gtest/gtest.h>
#include "bamsix/types.hpp"
#include "streamencode/streamencode.hpp"
#include "fmindex/fmindex.hpp"
#include "bitvectors/bitvectors.hpp"
#include "sais/sais.hpp"
#include "seqbuilder/seqbuilder.hpp"
#include "query/query.hpp"

#include <cstring>
#include <random>
#include <algorithm>

using namespace bamsix;

namespace {

struct DeterminismDataset {
    std::vector<OrderedRead> reads;
    SequenceBundle bundle;
    SaisResult sais;
    FMIndexEngine fm;
    BitVectors bv;
    std::vector<std::string> chrom_names;
};

DeterminismDataset BuildDeterminismDataset() {
    DeterminismDataset ds;

    ds.reads.resize(8);
    const char* seqs[] = {"ACGTACGT", "GCATAGCG", "TTAACCGG", "ACGTACGT",
                          "GGCCAATT", "TACGTACG", "ATCGATCG", "ACGTACGT"};
    for (int i = 0; i < 8; ++i) {
        auto s = seqs[i];
        size_t len = std::strlen(s);
        ds.reads[i].seq.resize(len);
        ds.reads[i].qual.resize(len);
        for (size_t j = 0; j < len; ++j) {
            switch (s[j]) {
                case 'A': ds.reads[i].seq[j] = 0; break;
                case 'C': ds.reads[i].seq[j] = 1; break;
                case 'G': ds.reads[i].seq[j] = 2; break;
                case 'T': ds.reads[i].seq[j] = 3; break;
                default:  ds.reads[i].seq[j] = 4; break;
            }
            ds.reads[i].qual[j] = 30 + j;
        }
        ds.reads[i].chrom = (i < 4) ? "chr1" : "chr2";
        ds.reads[i].chrom_id = (i < 4) ? 0 : 1;
        ds.reads[i].pos = 100 + i * 50;
        ds.reads[i].read_id = i;
        ds.reads[i].flag = 0;
        ds.reads[i].mapq = 60;
        ds.reads[i].cigar = {{0, static_cast<uint32_t>(len)}};
        ds.reads[i].source_file_id = 0;
        ds.reads[i].bam_offset = i;
    }

    ds.bundle = BuildSequence(ds.reads);
    ds.sais = ComputeSuffixArray(ds.bundle);
    ds.fm.Build(ds.sais.BWT, ds.sais.SA, ds.sais.sentinel_row, 4, ds.bundle.S.size());

    std::vector<uint64_t> starts(ds.bundle.readStarts.begin(), ds.bundle.readStarts.end());
    ds.bv.B_read.Build(starts, ds.bundle.S.size());

    ds.chrom_names = {"chr1", "chr2"};
    return ds;
}

}  // namespace

// ═══════════════════════════════════════════════════════════════════════════════
// Test 1: Locate is deterministic across two invocations (sorted mode)
// ═══════════════════════════════════════════════════════════════════════════════

TEST(Determinism, LocateSortedDeterministic) {
    auto ds = BuildDeterminismDataset();

    std::vector<uint8_t> P = {0, 1, 2, 3}; // ACGT

    auto results1 = Locate(P, ds.fm, ds.bv.B_read, ds.reads,
                           ds.chrom_names, StrandMode::StrandComplete, true);
    auto results2 = Locate(P, ds.fm, ds.bv.B_read, ds.reads,
                           ds.chrom_names, StrandMode::StrandComplete, true);

    ASSERT_EQ(results1.size(), results2.size())
        << "Locate sorted mode produced different result counts";

    for (size_t i = 0; i < results1.size(); ++i) {
        EXPECT_TRUE(results1[i].chrom == results2[i].chrom)
            << "Mismatch at result " << i << " chrom";
        EXPECT_TRUE(results1[i].p_min == results2[i].p_min)
            << "Mismatch at result " << i << " p_min";
        EXPECT_TRUE(results1[i].p_max == results2[i].p_max)
            << "Mismatch at result " << i << " p_max";
        EXPECT_TRUE(results1[i].read_id == results2[i].read_id)
            << "Mismatch at result " << i << " read_id";
        EXPECT_TRUE(static_cast<int>(results1[i].query_strand) ==
                    static_cast<int>(results2[i].query_strand))
            << "Mismatch at result " << i << " strand";
        EXPECT_TRUE(results1[i].sa_row == results2[i].sa_row)
            << "Mismatch at result " << i << " sa_row";
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// Test 2: FM-index construction is deterministic
// ═══════════════════════════════════════════════════════════════════════════════

TEST(Determinism, FMIndexConstructionDeterministic) {
    auto ds1 = BuildDeterminismDataset();
    auto ds2 = BuildDeterminismDataset();

    ASSERT_EQ(ds1.sais.BWT.size(), ds2.sais.BWT.size());
    EXPECT_TRUE(ds1.sais.BWT == ds2.sais.BWT) << "BWT not deterministic";

    ASSERT_EQ(ds1.sais.SA.size(), ds2.sais.SA.size());
    EXPECT_TRUE(ds1.sais.SA == ds2.sais.SA) << "SA not deterministic";

    EXPECT_TRUE(ds1.sais.sentinel_row == ds2.sais.sentinel_row)
        << "Sentinel row not deterministic";
}

// ═══════════════════════════════════════════════════════════════════════════════
// Test 3: Stream encoding is deterministic
// ═══════════════════════════════════════════════════════════════════════════════

TEST(Determinism, StreamEncodingDeterministic) {
    auto ds = BuildDeterminismDataset();

    auto q1 = EncodeQualStream(ds.reads, QualCodec::RANGE_CYCLE, 0, 0);
    auto q2 = EncodeQualStream(ds.reads, QualCodec::RANGE_CYCLE, 0, 0);
    EXPECT_TRUE(q1.payload == q2.payload) << "S_qual payload not deterministic";

    auto m1 = EncodeMetaStream(ds.reads, MetaCodec::TYPED_SPLIT);
    auto m2 = EncodeMetaStream(ds.reads, MetaCodec::TYPED_SPLIT);
    EXPECT_TRUE(m1.payload == m2.payload) << "S_meta payload not deterministic";

    auto mp1 = EncodeMapStream(ds.reads, MapCodec::DELTA_RANGE);
    auto mp2 = EncodeMapStream(ds.reads, MapCodec::DELTA_RANGE);
    EXPECT_TRUE(mp1.payload == mp2.payload) << "S_map payload not deterministic";
}

// ═══════════════════════════════════════════════════════════════════════════════
// Test 4: GlobalCount is deterministic
// ═══════════════════════════════════════════════════════════════════════════════

TEST(Determinism, GlobalCountDeterministic) {
    auto ds = BuildDeterminismDataset();
    std::vector<uint8_t> P = {0, 1, 2, 3}; // ACGT

    auto c1 = GlobalCount(P, ds.fm, StrandMode::StrandComplete);
    auto c2 = GlobalCount(P, ds.fm, StrandMode::StrandComplete);
    EXPECT_TRUE(c1 == c2) << "GlobalCount not deterministic";

    auto e1 = GlobalExists(P, ds.fm, StrandMode::StrandComplete);
    auto e2 = GlobalExists(P, ds.fm, StrandMode::StrandComplete);
    EXPECT_TRUE(e1 == e2) << "GlobalExists not deterministic";
}

// ═══════════════════════════════════════════════════════════════════════════════
// Test 5: Streaming and sorted modes return same multiset
// ═══════════════════════════════════════════════════════════════════════════════

TEST(Determinism, StreamingAndSortedSameMultiset) {
    auto ds = BuildDeterminismDataset();
    std::vector<uint8_t> P = {0, 1, 2, 3}; // ACGT

    auto streaming = Locate(P, ds.fm, ds.bv.B_read, ds.reads,
                            ds.chrom_names, StrandMode::StrandComplete, false);
    auto sorted = Locate(P, ds.fm, ds.bv.B_read, ds.reads,
                         ds.chrom_names, StrandMode::StrandComplete, true);

    ASSERT_EQ(streaming.size(), sorted.size())
        << "Streaming and sorted modes returned different result counts";

    auto streaming_sorted = streaming;
    std::sort(streaming_sorted.begin(), streaming_sorted.end(),
        [](const Match& a, const Match& b) {
            if (a.query_strand != b.query_strand) return a.query_strand < b.query_strand;
            if (a.chrom != b.chrom) return a.chrom < b.chrom;
            if (a.p_min != b.p_min) return a.p_min < b.p_min;
            if (a.p_max != b.p_max) return a.p_max < b.p_max;
            if (a.read_id != b.read_id) return a.read_id < b.read_id;
            return a.sa_row < b.sa_row;
        });

    for (size_t i = 0; i < sorted.size(); ++i) {
        EXPECT_TRUE(streaming_sorted[i].chrom == sorted[i].chrom);
        EXPECT_TRUE(streaming_sorted[i].p_min == sorted[i].p_min);
        EXPECT_TRUE(streaming_sorted[i].p_max == sorted[i].p_max);
        EXPECT_TRUE(streaming_sorted[i].read_id == sorted[i].read_id);
        EXPECT_TRUE(static_cast<int>(streaming_sorted[i].query_strand) ==
                    static_cast<int>(sorted[i].query_strand));
    }
}
