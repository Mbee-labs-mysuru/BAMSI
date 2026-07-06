/**
 * test_locate_sorted.cpp — Locate sorted-mode ordering verification
 *
 * Contract §4.1 (C1): Locate output ordering has two modes — streaming
 * (default, no sort) and sorted (on --sort-output flag); both return the
 * same set.
 *
 * Sorted order (normative):
 *   1. query_strand: Forward before Reverse
 *   2. chrom_id ascending
 *   3. p_min ascending
 *   4. p_max ascending
 *   5. read_id ascending
 *   6. SA row ascending (final tie-break)
 */

#include <gtest/gtest.h>
#include "bamsix/types.hpp"
#include "query/query.hpp"
#include "fmindex/fmindex.hpp"
#include "bitvectors/bitvectors.hpp"
#include "sais/sais.hpp"
#include "seqbuilder/seqbuilder.hpp"

#include <algorithm>
#include <cstring>

using namespace bamsix;

namespace {

struct SortTestDataset {
    std::vector<OrderedRead> reads;
    SequenceBundle bundle;
    SaisResult sais;
    FMIndexEngine fm;
    BitVectors bv;
    std::vector<std::string> chrom_names;
};

SortTestDataset BuildSortTestDataset() {
    SortTestDataset ds;
    ds.reads.resize(6);

    // chr1 reads (chrom_id=0)
    ds.reads[0].seq = {0,1,2,0,1,2};  // ACGACG
    ds.reads[0].chrom = "chr1"; ds.reads[0].chrom_id = 0;
    ds.reads[0].pos = 100; ds.reads[0].read_id = 0;
    ds.reads[0].flag = 0; ds.reads[0].mapq = 60;
    ds.reads[0].cigar = {{0,6}};
    ds.reads[0].qual = {30,31,32,33,34,35};
    ds.reads[0].source_file_id = 0; ds.reads[0].bam_offset = 0;

    ds.reads[1].seq = {3,0,1,2,3,3};  // TACGTT
    ds.reads[1].chrom = "chr1"; ds.reads[1].chrom_id = 0;
    ds.reads[1].pos = 200; ds.reads[1].read_id = 1;
    ds.reads[1].flag = 0; ds.reads[1].mapq = 50;
    ds.reads[1].cigar = {{0,6}};
    ds.reads[1].qual = {25,26,27,28,29,30};
    ds.reads[1].source_file_id = 0; ds.reads[1].bam_offset = 1;

    ds.reads[2].seq = {2,0,1,2,4,4};  // GACGNN
    ds.reads[2].chrom = "chr1"; ds.reads[2].chrom_id = 0;
    ds.reads[2].pos = 150; ds.reads[2].read_id = 2;
    ds.reads[2].flag = 0; ds.reads[2].mapq = 40;
    ds.reads[2].cigar = {{0,6}};
    ds.reads[2].qual = {20,21,22,23,24,25};
    ds.reads[2].source_file_id = 0; ds.reads[2].bam_offset = 2;

    // chr2 reads (chrom_id=1)
    ds.reads[3].seq = {0,0,1,2,3,0};  // AACGTA
    ds.reads[3].chrom = "chr2"; ds.reads[3].chrom_id = 1;
    ds.reads[3].pos = 50; ds.reads[3].read_id = 3;
    ds.reads[3].flag = 0; ds.reads[3].mapq = 60;
    ds.reads[3].cigar = {{0,6}};
    ds.reads[3].qual = {35,36,37,38,39,40};
    ds.reads[3].source_file_id = 0; ds.reads[3].bam_offset = 3;

    ds.reads[4].seq = {0,1,2,0,1,2};  // ACGACG
    ds.reads[4].chrom = "chr2"; ds.reads[4].chrom_id = 1;
    ds.reads[4].pos = 75; ds.reads[4].read_id = 4;
    ds.reads[4].flag = 0; ds.reads[4].mapq = 55;
    ds.reads[4].cigar = {{0,6}};
    ds.reads[4].qual = {40,41,42,43,44,45};
    ds.reads[4].source_file_id = 0; ds.reads[4].bam_offset = 4;

    ds.reads[5].seq = {1,2,3,0,1,2};  // CGTACG
    ds.reads[5].chrom = "chr2"; ds.reads[5].chrom_id = 1;
    ds.reads[5].pos = 300; ds.reads[5].read_id = 5;
    ds.reads[5].flag = 0; ds.reads[5].mapq = 45;
    ds.reads[5].cigar = {{0,6}};
    ds.reads[5].qual = {30,31,32,33,34,35};
    ds.reads[5].source_file_id = 0; ds.reads[5].bam_offset = 5;

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
// Test 1: Sorted results obey 6-key ordering (forward strand only)
// ═══════════════════════════════════════════════════════════════════════════════

TEST(LocateSorted, ForwardStrandSortedOrdering) {
    auto ds = BuildSortTestDataset();
    std::vector<uint8_t> P = {0, 1, 2}; // ACG

    auto results = Locate(P, ds.fm, ds.bv.B_read, ds.reads,
                          ds.chrom_names, StrandMode::SingleStrand, true);

    ASSERT_GT(results.size(), 0u) << "Expected at least one match for ACG";

    // Verify all results are Forward strand (single-strand mode)
    for (const auto& m : results) {
        EXPECT_TRUE(static_cast<int>(m.query_strand) ==
                    static_cast<int>(QueryStrand::Forward));
    }

    // Verify sorted order: chrom < p_min < p_max < read_id < sa_row
    for (size_t i = 1; i < results.size(); ++i) {
        const auto& a = results[i-1];
        const auto& b = results[i];

        bool order_ok = false;
        if (a.chrom < b.chrom) { order_ok = true; }
        else if (a.chrom == b.chrom) {
            if (a.p_min < b.p_min) { order_ok = true; }
            else if (a.p_min == b.p_min) {
                if (a.p_max < b.p_max) { order_ok = true; }
                else if (a.p_max == b.p_max) {
                    if (a.read_id < b.read_id) { order_ok = true; }
                    else if (a.read_id == b.read_id) {
                        order_ok = (a.sa_row <= b.sa_row);
                    }
                }
            }
        }

        EXPECT_TRUE(order_ok)
            << "Sort violation at positions " << (i-1) << " and " << i
            << ": (" << a.chrom << "," << a.p_min << "," << a.p_max << "," << a.read_id << ")"
            << " vs (" << b.chrom << "," << b.p_min << "," << b.p_max << "," << b.read_id << ")";
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// Test 2: Both-strand search sorts Forward before Reverse
// ═══════════════════════════════════════════════════════════════════════════════

TEST(LocateSorted, BothStrandForwardBeforeReverse) {
    auto ds = BuildSortTestDataset();
    std::vector<uint8_t> P = {0, 1, 2}; // ACG — rc is CGT (not self-complementary)

    auto results = Locate(P, ds.fm, ds.bv.B_read, ds.reads,
                          ds.chrom_names, StrandMode::StrandComplete, true);

    bool seen_reverse = false;
    for (size_t i = 0; i < results.size(); ++i) {
        if (results[i].query_strand == QueryStrand::Reverse) {
            seen_reverse = true;
        }
        if (seen_reverse && results[i].query_strand == QueryStrand::Forward) {
            FAIL() << "Forward result at position " << i
                   << " appears after Reverse results — violates Contract §4.1 key 1";
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// Test 3: Multi-chromosome results are in chrom order
// ═══════════════════════════════════════════════════════════════════════════════

TEST(LocateSorted, MultiChromosomeOrdering) {
    auto ds = BuildSortTestDataset();
    std::vector<uint8_t> P = {0, 1, 2}; // ACG

    auto results = Locate(P, ds.fm, ds.bv.B_read, ds.reads,
                          ds.chrom_names, StrandMode::SingleStrand, true);

    std::string last_chrom = "";
    for (const auto& m : results) {
        if (!last_chrom.empty() && m.chrom < last_chrom) {
            FAIL() << "Chromosome order violation: " << m.chrom
                   << " appears after " << last_chrom;
        }
        last_chrom = m.chrom;
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// Test 4: Streaming and sorted modes return same multiset
// ═══════════════════════════════════════════════════════════════════════════════

TEST(LocateSorted, StreamingAndSortedSameMultiset) {
    auto ds = BuildSortTestDataset();
    std::vector<uint8_t> P = {0, 1, 2}; // ACG

    auto streaming = Locate(P, ds.fm, ds.bv.B_read, ds.reads,
                            ds.chrom_names, StrandMode::StrandComplete, false);
    auto sorted = Locate(P, ds.fm, ds.bv.B_read, ds.reads,
                         ds.chrom_names, StrandMode::StrandComplete, true);

    ASSERT_EQ(streaming.size(), sorted.size())
        << "Streaming (" << streaming.size() << ") and sorted ("
        << sorted.size() << ") modes returned different counts";

    auto cmp = [](const Match& a, const Match& b) {
        if (a.query_strand != b.query_strand) return a.query_strand < b.query_strand;
        if (a.chrom != b.chrom) return a.chrom < b.chrom;
        if (a.p_min != b.p_min) return a.p_min < b.p_min;
        if (a.p_max != b.p_max) return a.p_max < b.p_max;
        if (a.read_id != b.read_id) return a.read_id < b.read_id;
        return a.sa_row < b.sa_row;
    };
    std::sort(streaming.begin(), streaming.end(), cmp);

    for (size_t i = 0; i < sorted.size(); ++i) {
        EXPECT_TRUE(streaming[i].chrom == sorted[i].chrom)
            << "Multiset mismatch at " << i << ": chrom";
        EXPECT_TRUE(streaming[i].p_min == sorted[i].p_min)
            << "Multiset mismatch at " << i << ": p_min";
        EXPECT_TRUE(streaming[i].p_max == sorted[i].p_max)
            << "Multiset mismatch at " << i << ": p_max";
        EXPECT_TRUE(streaming[i].read_id == sorted[i].read_id)
            << "Multiset mismatch at " << i << ": read_id";
        EXPECT_TRUE(static_cast<int>(streaming[i].query_strand) ==
                    static_cast<int>(sorted[i].query_strand))
            << "Multiset mismatch at " << i << ": strand";
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// Test 5: GlobalCount matches Locate result count
// ═══════════════════════════════════════════════════════════════════════════════

TEST(LocateSorted, CountMatchesLocate) {
    auto ds = BuildSortTestDataset();
    std::vector<uint8_t> P = {0, 1, 2}; // ACG

    auto count = GlobalCount(P, ds.fm, StrandMode::StrandComplete);
    auto results = Locate(P, ds.fm, ds.bv.B_read, ds.reads,
                          ds.chrom_names, StrandMode::StrandComplete, true);

    EXPECT_TRUE(count == results.size())
        << "GlobalCount (" << count << ") != Locate result count (" << results.size() << ")";
}
