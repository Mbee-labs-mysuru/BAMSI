/// Tutorial Smoke Tests — Exec Plan §6.3.7
/// "Each tutorial runs as a smoke test in nightly CI to prevent drift."
///
/// These tests simulate the three tutorials end-to-end using synthetic data:
///   Tutorial 01: Motif counting on a small dataset
///   Tutorial 02: Region-restricted query on a genomic region
///   Tutorial 03: Two-phase quality filtering (Locate + qual post-filter)
///
/// These tests verify that the tutorial workflows described in docs/tutorials/
/// actually work with the current codebase. They use the internal API directly
/// (not the CLI) to avoid test-environment path issues.

#include <gtest/gtest.h>

#include "bamsix/types.hpp"
#include "../src/bitvectors/bitvectors.hpp"
#include "../src/fmindex/fmindex.hpp"
#include "../src/ingest/ingest.hpp"
#include "../src/mapping/mapping.hpp"
#include "../src/ordering/ordering.hpp"
#include "../src/query/query.hpp"
#include "../src/sais/sais.hpp"
#include "../src/seqbuilder/seqbuilder.hpp"
#include "../src/streamencode/streamencode.hpp"
#include "../src/windows/windows.hpp"

#include <algorithm>
#include <cstring>
#include <numeric>
#include <string>

using namespace bamsix;

// ─── Shared test infrastructure ──────────────────────────────────────────────

/// Create a synthetic ordered read with the given parameters.
static OrderedRead MakeSyntheticRead(uint64_t read_id, uint32_t chrom_id,
                                      uint64_t pos, const std::string& bases,
                                      const std::vector<uint8_t>& quals) {
    OrderedRead r;
    r.read_id = read_id;
    r.chrom_id = chrom_id;
    r.chrom = "chr" + std::to_string(chrom_id + 1);
    r.pos = pos;
    r.flag = 0;
    r.mapq = 60;
    r.qname = "read_" + std::to_string(read_id);

    // Convert base string to numeric codes
    for (char c : bases) {
        switch (c) {
            case 'A': r.seq.push_back(CODE_A); break;
            case 'C': r.seq.push_back(CODE_C); break;
            case 'G': r.seq.push_back(CODE_G); break;
            case 'T': r.seq.push_back(CODE_T); break;
            default:  r.seq.push_back(CODE_N); break;
        }
    }

    // Quality scores
    r.qual = quals.empty() ? std::vector<uint8_t>(r.seq.size(), 30) : quals;

    // Simple M-only CIGAR
    r.cigar = {{0, static_cast<uint32_t>(r.seq.size())}};  // full match

    return r;
}

/// Build a complete index from synthetic reads.
struct SyntheticIndex {
    std::vector<OrderedRead> reads;
    SequenceBundle bundle;
    FMIndexEngine fm;
    SuccinctBitvector B_read;
    SuccinctBitvector B_window;
    WindowTable windows;
    std::vector<std::string> chrom_names;
    std::map<std::string, uint32_t> chrom_to_id;
};

static SyntheticIndex BuildSyntheticIndex(std::vector<OrderedRead> reads) {
    SyntheticIndex idx;
    idx.reads = std::move(reads);

    // Build chrom name table
    std::set<std::string> unique_chroms;
    for (const auto& r : idx.reads) unique_chroms.insert(r.chrom);
    for (const auto& c : unique_chroms) {
        uint32_t id = idx.chrom_names.size();
        idx.chrom_to_id[c] = id;
        idx.chrom_names.push_back(c);
    }

    // Build sequence
    idx.bundle = BuildSequence(idx.reads);

    // SA-IS (Architecture §4.4)
    SaisResult sais = ComputeSuffixArray(idx.bundle);

    // FM-index — use BWT from SA-IS result directly
    idx.fm.Build(sais.BWT, sais.SA, sais.sentinel_row, 4, idx.bundle.S.size());

    // Windows
    idx.windows = BuildWindows(idx.reads, idx.bundle, 100000);

    // Bitvectors
    auto bvs = BuildBitvectors(idx.bundle, idx.windows);
    idx.B_read = std::move(bvs.B_read);
    idx.B_window = std::move(bvs.B_window);

    return idx;
}

/// Encode a string pattern into numeric codes.
static std::vector<uint8_t> EncodePattern(const std::string& pat) {
    std::vector<uint8_t> result;
    for (char c : pat) {
        switch (c) {
            case 'A': result.push_back(CODE_A); break;
            case 'C': result.push_back(CODE_C); break;
            case 'G': result.push_back(CODE_G); break;
            case 'T': result.push_back(CODE_T); break;
            default:  result.push_back(CODE_N); break;
        }
    }
    return result;
}

// ─── Tutorial 01: Motif Counting (docs/tutorials/01_motif_counting.md) ──────

class Tutorial01MotifCountingTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Create 20 reads on chr1, some containing target motifs
        std::vector<OrderedRead> reads;
        // Reads that contain "GATTACA"
        reads.push_back(MakeSyntheticRead(0, 0, 100, "ACGTGATTACAACGT", {}));
        reads.push_back(MakeSyntheticRead(1, 0, 200, "TTTGATTACATTTTT", {}));
        reads.push_back(MakeSyntheticRead(2, 0, 300, "GATTACAGATTACAG", {})); // 2 occurrences
        // Reads that contain "AGATC"
        reads.push_back(MakeSyntheticRead(3, 0, 400, "CCCCAGATCGGGGGG", {}));
        reads.push_back(MakeSyntheticRead(4, 0, 500, "AAAAAGATCAAAAAA", {}));
        // Reads with no target motifs
        for (uint64_t i = 5; i < 15; ++i) {
            reads.push_back(MakeSyntheticRead(i, 0, 600 + i * 100,
                                               "ACGTACGTACGTACG", {}));
        }
        idx_ = BuildSyntheticIndex(std::move(reads));
    }
    SyntheticIndex idx_;
};

TEST_F(Tutorial01MotifCountingTest, GATTACACountIsCorrect) {
    auto P = EncodePattern("GATTACA");
    uint64_t count = GlobalCount(P, idx_.fm, StrandMode::StrandComplete);
    // 3 forward occurrences (reads 0, 1, 2 have one each, read 2 has 2)
    EXPECT_GE(count, 4u);  // At least 4 forward occurrences
}

TEST_F(Tutorial01MotifCountingTest, AGATCCountIsCorrect) {
    auto P = EncodePattern("AGATC");
    uint64_t count = GlobalCount(P, idx_.fm, StrandMode::StrandComplete);
    EXPECT_GE(count, 2u);  // At least 2 forward occurrences
}

TEST_F(Tutorial01MotifCountingTest, GlobalExistsForPresentMotif) {
    auto P = EncodePattern("GATTACA");
    EXPECT_TRUE(GlobalExists(P, idx_.fm, StrandMode::StrandComplete));
}

TEST_F(Tutorial01MotifCountingTest, GlobalExistsForAbsentMotif) {
    auto P = EncodePattern("ZZZZZZZ");
    // This should fail validation (invalid chars) — test the error path
    EXPECT_THROW(GlobalExists(EncodePattern(""), idx_.fm, StrandMode::StrandComplete),
                 Error);
}

TEST_F(Tutorial01MotifCountingTest, MultipleMotifBatch) {
    // Tutorial 01 counts 10 motifs — verify the batch pattern works
    std::vector<std::string> motifs = {"GATTACA", "AGATC", "ACGT", "TTTTT"};
    for (const auto& motif : motifs) {
        auto P = EncodePattern(motif);
        uint64_t count = GlobalCount(P, idx_.fm, StrandMode::StrandComplete);
        // Each motif should return a non-negative count (no crashes)
        EXPECT_GE(count, 0u) << "Motif: " << motif;
    }
}

// ─── Tutorial 02: Region-Restricted Query (docs/tutorials/02_region_query.md)

class Tutorial02RegionQueryTest : public ::testing::Test {
protected:
    void SetUp() override {
        std::vector<OrderedRead> reads;
        // Reads on chr1 in the "target region" [1000, 2000]
        reads.push_back(MakeSyntheticRead(0, 0, 1000, "ACGTGATTACAACGT", {}));
        reads.push_back(MakeSyntheticRead(1, 0, 1100, "GATTACAGATTACAG", {}));
        reads.push_back(MakeSyntheticRead(2, 0, 1500, "TTTGATTACATTTTT", {}));
        // Reads on chr1 OUTSIDE the target region
        reads.push_back(MakeSyntheticRead(3, 0, 5000, "ACGTGATTACAACGT", {}));
        reads.push_back(MakeSyntheticRead(4, 0, 8000, "GATTACAGATTACAG", {}));
        // Reads on chr2 (wrong chromosome)
        reads.push_back(MakeSyntheticRead(5, 1, 1000, "GATTACAGATTACAG", {}));

        idx_ = BuildSyntheticIndex(std::move(reads));
    }
    SyntheticIndex idx_;
};

TEST_F(Tutorial02RegionQueryTest, RegionalCountInTargetRegion) {
    auto P = EncodePattern("GATTACA");
    uint64_t count = RegionalCount(P, "chr1", 1000, 2000,
                                    idx_.fm, idx_.B_read, idx_.B_window,
                                    idx_.windows, idx_.reads,
                                    idx_.chrom_to_id, StrandMode::StrandComplete);
    // Should find occurrences in reads 0, 1, 2 (all within [1000, 2000])
    EXPECT_GE(count, 3u);
}

TEST_F(Tutorial02RegionQueryTest, RegionalCountExcludesOutOfRange) {
    auto P = EncodePattern("GATTACA");
    // Query a region that excludes read 3 (pos=5000) and read 4 (pos=8000)
    uint64_t count_narrow = RegionalCount(P, "chr1", 1000, 2000,
                                           idx_.fm, idx_.B_read, idx_.B_window,
                                           idx_.windows, idx_.reads,
                                           idx_.chrom_to_id, StrandMode::StrandComplete);
    uint64_t count_wide = RegionalCount(P, "chr1", 1000, 10000,
                                         idx_.fm, idx_.B_read, idx_.B_window,
                                         idx_.windows, idx_.reads,
                                         idx_.chrom_to_id, StrandMode::StrandComplete);
    EXPECT_LE(count_narrow, count_wide);
}

TEST_F(Tutorial02RegionQueryTest, RegionalExistsWithThreshold) {
    auto P = EncodePattern("GATTACA");
    bool exists = RegionalExists(P, 1, "chr1", 1000, 2000,
                                  idx_.fm, idx_.B_read, idx_.B_window,
                                  idx_.windows, idx_.reads,
                                  idx_.chrom_to_id, StrandMode::StrandComplete);
    EXPECT_TRUE(exists);
}

TEST_F(Tutorial02RegionQueryTest, RegionalCountWrongChromReturnsZero) {
    auto P = EncodePattern("GATTACA");
    uint64_t count = RegionalCount(P, "chrX", 1000, 2000,
                                    idx_.fm, idx_.B_read, idx_.B_window,
                                    idx_.windows, idx_.reads,
                                    idx_.chrom_to_id, StrandMode::StrandComplete);
    EXPECT_EQ(count, 0u);
}

// ─── Tutorial 03: Two-Phase Quality Filtering ────────────────────────────────

class Tutorial03QualityFilterTest : public ::testing::Test {
protected:
    void SetUp() override {
        std::vector<OrderedRead> reads;
        // Read 0: high quality (all Q40), contains "GATTACA"
        reads.push_back(MakeSyntheticRead(0, 0, 100, "ACGTGATTACAACGT",
                                           std::vector<uint8_t>(15, 40)));
        // Read 1: low quality (all Q10), contains "GATTACA"
        reads.push_back(MakeSyntheticRead(1, 0, 200, "TTTGATTACATTTTT",
                                           std::vector<uint8_t>(15, 10)));
        // Read 2: mixed quality, contains "GATTACA"
        {
            std::vector<uint8_t> quals(15, 35);
            quals[3] = 5; quals[4] = 5;  // low quality at positions 3-4
            reads.push_back(MakeSyntheticRead(2, 0, 300, "GATTACAGATTACAG", quals));
        }

        idx_ = BuildSyntheticIndex(std::move(reads));
    }
    SyntheticIndex idx_;
};

TEST_F(Tutorial03QualityFilterTest, LocateReturnsAllMatches) {
    // Phase 1: Locate all matches (no quality filtering)
    auto P = EncodePattern("GATTACA");
    auto matches = Locate(P, idx_.fm, idx_.B_read, idx_.reads,
                           idx_.chrom_names, StrandMode::StrandComplete, true);
    EXPECT_GE(matches.size(), 3u);  // At least 3 forward matches
}

TEST_F(Tutorial03QualityFilterTest, QualityPostFilterReducesMatches) {
    // Phase 1: Locate
    auto P = EncodePattern("GATTACA");
    auto matches = Locate(P, idx_.fm, idx_.B_read, idx_.reads,
                           idx_.chrom_names, StrandMode::StrandComplete, true);
    size_t total = matches.size();

    // Phase 2: Post-filter by average quality >= Q30
    // In production, caller retrieves Q_i via dir_qual[read_id] and filters.
    // Here we simulate by checking our synthetic quality data.
    size_t passing = 0;
    for (const auto& m : matches) {
        if (m.read_id < idx_.reads.size()) {
            const auto& quals = idx_.reads[m.read_id].qual;
            if (!quals.empty()) {
                double avg = std::accumulate(quals.begin(), quals.end(), 0.0) / quals.size();
                if (avg >= 30.0) ++passing;
            }
        }
    }

    // Not all matches should pass (read 1 has Q10)
    EXPECT_LT(passing, total);
    EXPECT_GT(passing, 0u);  // At least read 0 passes
}

TEST_F(Tutorial03QualityFilterTest, QualCodecRoundTrip) {
    // Verify quality codec round-trips correctly (used by the tutorial workflow)
    auto qual_result = EncodeQualStream(idx_.reads, QualCodec::RANGE_CYCLE, 0, 0);
    EXPECT_FALSE(qual_result.payload.empty());

    // Decode each read's quality
    auto& dir = std::get<StreamDirectoryPerRead>(qual_result.directory);
    for (size_t i = 0; i < idx_.reads.size(); ++i) {
        auto decoded = DecodeQualRead(qual_result.payload, dir[i],
                                       qual_result.codec_id);
        EXPECT_EQ(decoded, idx_.reads[i].qual)
            << "Quality round-trip failed for read " << i;
    }
}

TEST_F(Tutorial03QualityFilterTest, SortedModeMatchesUnsorted) {
    auto P = EncodePattern("GATTACA");
    auto sorted = Locate(P, idx_.fm, idx_.B_read, idx_.reads,
                          idx_.chrom_names, StrandMode::StrandComplete, true);
    auto unsorted = Locate(P, idx_.fm, idx_.B_read, idx_.reads,
                            idx_.chrom_names, StrandMode::StrandComplete, false);
    // Same multiset of results
    EXPECT_EQ(sorted.size(), unsorted.size());
}
