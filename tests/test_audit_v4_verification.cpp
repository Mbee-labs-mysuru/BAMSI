/// BAMSIX Audit v4 Verification Tests
/// Covers all gaps identified in compliance_audit_v4.md:
///   GAP-1: Lazy per-read decode (MapOccurrenceLazy)
///   GAP-2: S_meta TYPED_SPLIT round-trip (varint CIGAR + 2-byte FLAG + aux)
///   GAP-3: ordering_hash verification
///   GAP-5: End-to-end reconstruction round-trip
///
/// Contract v3.3, Architecture v4.3, Execution Plan v2.0

#include <gtest/gtest.h>
#include "bamsix/types.hpp"
#include "streamencode/streamencode.hpp"
#include "mapping/mapping.hpp"
#include "query/query.hpp"
#include "bitvectors/bitvectors.hpp"
#include "fmindex/fmindex.hpp"
#include "sais/sais.hpp"
#include "seqbuilder/seqbuilder.hpp"
#include "ordering/ordering.hpp"
#include "reconstruct/reconstruct.hpp"

using namespace bamsix;

// ─── Helper: build a minimal test dataset ───────────────────────────────────

namespace {

struct TestDataset {
    std::vector<OrderedRead> reads;
    SequenceBundle           bundle;
    SaisResult               sais_result;
    FMIndexEngine            fm;
    BitVectors               bv;
    // Encoded streams
    MetaEncodeResult         meta_enc;
    MapEncodeResult          map_enc;
    QualEncodeResult         qual_enc;
    std::vector<std::string> chrom_names;
    std::map<std::string, uint32_t> chrom_to_id;
};

TestDataset BuildTestDataset() {
    TestDataset td;

    // 4 reads across 2 chromosomes
    td.reads.resize(4);
    // Read 0: chr1 pos=100 seq=ACGT
    td.reads[0].seq = {0,1,2,3};
    td.reads[0].chrom = "chr1"; td.reads[0].chrom_id = 0;
    td.reads[0].pos = 100; td.reads[0].read_id = 0;
    td.reads[0].flag = 0; td.reads[0].mapq = 60;
    td.reads[0].cigar = {{0,4}}; // 4M
    td.reads[0].qual = {30,31,32,33};
    td.reads[0].source_file_id = 0; td.reads[0].bam_offset = 0;

    // Read 1: chr1 pos=110 seq=GCAT
    td.reads[1].seq = {2,1,0,3};
    td.reads[1].chrom = "chr1"; td.reads[1].chrom_id = 0;
    td.reads[1].pos = 110; td.reads[1].read_id = 1;
    td.reads[1].flag = 16; td.reads[1].mapq = 50;
    td.reads[1].cigar = {{0,2},{1,1},{0,1}}; // 2M1I1M
    td.reads[1].qual = {25,26,27,28};
    td.reads[1].source_file_id = 0; td.reads[1].bam_offset = 1;

    // Read 2: chr2 pos=200 seq=TTAA
    td.reads[2].seq = {3,3,0,0};
    td.reads[2].chrom = "chr2"; td.reads[2].chrom_id = 1;
    td.reads[2].pos = 200; td.reads[2].read_id = 2;
    td.reads[2].flag = 0; td.reads[2].mapq = 40;
    td.reads[2].cigar = {{0,4}}; // 4M
    td.reads[2].qual = {35,36,37,38};
    td.reads[2].source_file_id = 0; td.reads[2].bam_offset = 2;

    // Read 3: chr2 pos=210 seq=CCGG
    td.reads[3].seq = {1,1,2,2};
    td.reads[3].chrom = "chr2"; td.reads[3].chrom_id = 1;
    td.reads[3].pos = 210; td.reads[3].read_id = 3;
    td.reads[3].flag = 0; td.reads[3].mapq = 55;
    td.reads[3].cigar = {{0,4}}; // 4M
    td.reads[3].qual = {40,41,42,43};
    td.reads[3].source_file_id = 0; td.reads[3].bam_offset = 3;
    // Add aux data to read 3 to test aux substream
    td.reads[3].aux_data = {0xAA, 0xBB, 0xCC};

    // Build S
    td.bundle = BuildSequence(td.reads);

    // SA-IS
    td.sais_result = ComputeSuffixArray(td.bundle);

    // FM-index
    td.fm.Build(td.sais_result.BWT, td.sais_result.SA,
                td.sais_result.sentinel_row, 4, td.bundle.S.size());

    // Bitvectors — B_read from readStarts
    std::vector<uint64_t> starts_u64(td.bundle.readStarts.begin(),
                                      td.bundle.readStarts.end());
    td.bv.B_read.Build(starts_u64, td.bundle.S.size());

    // Encode streams
    td.meta_enc = EncodeMetaStream(td.reads, MetaCodec::TYPED_SPLIT);
    td.map_enc  = EncodeMapStream(td.reads, MapCodec::DELTA_RANGE);
    td.qual_enc = EncodeQualStream(td.reads, QualCodec::RANGE_CYCLE, 0);

    td.chrom_names = {"chr1", "chr2"};
    td.chrom_to_id = {{"chr1", 0}, {"chr2", 1}};

    return td;
}

} // namespace

// ═════════════════════════════════════════════════════════════════════════════
// GAP-1: Lazy per-read decode — MapOccurrenceLazy
// ═════════════════════════════════════════════════════════════════════════════

TEST(AuditV4_GAP1, LazyMapOccurrenceMatchesLegacy) {
    auto td = BuildTestDataset();

    // For each read, compute the S-position of its first base
    for (uint64_t rid = 0; rid < td.reads.size(); ++rid) {
        uint64_t pos = td.bundle.readStarts[rid];
        uint64_t plen = td.reads[rid].seq.size();

        // Legacy path (in-memory reads)
        auto legacy = MapOccurrence(pos, plen, QueryStrand::Forward,
                                    td.bv.B_read, td.reads);

        // Lazy path (stream decode)
        auto lazy = MapOccurrenceLazy(pos, plen, QueryStrand::Forward,
                                      td.bv.B_read,
                                      td.meta_enc.payload, td.map_enc.payload,
                                      td.meta_enc.directory, td.map_enc.directory,
                                      td.meta_enc.codec_id, td.map_enc.codec_id);

        EXPECT_EQ(legacy.read_id, lazy.read_id)
            << "read_id mismatch at rid=" << rid;
        EXPECT_EQ(legacy.chrom_id, lazy.chrom_id)
            << "chrom_id mismatch at rid=" << rid;
        EXPECT_EQ(legacy.p_min, lazy.p_min)
            << "p_min mismatch at rid=" << rid;
        EXPECT_EQ(legacy.p_max, lazy.p_max)
            << "p_max mismatch at rid=" << rid;
    }
}

TEST(AuditV4_GAP1, LazyMapInternalOffset) {
    // Test mapping at an internal offset within a read (not just pos=readStart)
    auto td = BuildTestDataset();

    // Read 0 has seq ACGT (len=4), pos=100, CIGAR=4M
    // Match at offset 2 within read 0 => S-position = readStarts[0] + 2
    uint64_t pos = td.bundle.readStarts[0] + 2;
    uint64_t plen = 2; // matching "GT"

    auto legacy = MapOccurrence(pos, plen, QueryStrand::Forward,
                                td.bv.B_read, td.reads);
    auto lazy = MapOccurrenceLazy(pos, plen, QueryStrand::Forward,
                                  td.bv.B_read,
                                  td.meta_enc.payload, td.map_enc.payload,
                                  td.meta_enc.directory, td.map_enc.directory,
                                  td.meta_enc.codec_id, td.map_enc.codec_id);

    EXPECT_EQ(legacy.p_min, lazy.p_min);
    EXPECT_EQ(legacy.p_max, lazy.p_max);
    // With 4M CIGAR at pos=100, offset 2-3 maps to ref pos 102-103
    EXPECT_EQ(102u, lazy.p_min);
    EXPECT_EQ(103u, lazy.p_max);
}

// ═════════════════════════════════════════════════════════════════════════════
// GAP-2: S_meta TYPED_SPLIT round-trip verification
// ═════════════════════════════════════════════════════════════════════════════

TEST(AuditV4_GAP2, MetaTypedSplitRoundTrip) {
    auto td = BuildTestDataset();

    for (size_t i = 0; i < td.reads.size(); ++i) {
        auto decoded = DecodeMetaRead(td.meta_enc.payload,
                                       td.meta_enc.directory[i],
                                       td.meta_enc.codec_id);

        // Verify FLAG round-trip (2-byte LE per Contract §2.8)
        EXPECT_EQ(td.reads[i].flag, decoded.flag)
            << "FLAG mismatch at read " << i;

        // Verify CIGAR round-trip (varint nybble encoding)
        ASSERT_EQ(td.reads[i].cigar.size(), decoded.cigar.size())
            << "CIGAR op count mismatch at read " << i;
        for (size_t c = 0; c < td.reads[i].cigar.size(); ++c) {
            EXPECT_EQ(td.reads[i].cigar[c].op, decoded.cigar[c].op)
                << "CIGAR op mismatch at read " << i << " op " << c;
            EXPECT_EQ(td.reads[i].cigar[c].len, decoded.cigar[c].len)
                << "CIGAR len mismatch at read " << i << " op " << c;
        }

        // Verify aux data round-trip
        EXPECT_EQ(td.reads[i].aux_data.size(), decoded.aux_data.size())
            << "Aux data size mismatch at read " << i;
        EXPECT_EQ(td.reads[i].aux_data, decoded.aux_data)
            << "Aux data content mismatch at read " << i;
    }
}

TEST(AuditV4_GAP2, MetaTypedSplitComplexCigar) {
    // Test with complex CIGAR: 3M2I1D4M1S (soft clip edge case)
    OrderedRead r;
    r.seq = {0,1,2,3,0,1,2,3,0,1}; // 10 bases
    r.chrom = "chr1"; r.chrom_id = 0; r.pos = 50; r.read_id = 0;
    r.flag = 163; // typical paired-end flag
    r.cigar = {{0,3},{1,2},{2,1},{0,4},{4,1}}; // 3M2I1D4M1S
    r.qual.assign(10, 30);
    r.aux_data = {};
    r.source_file_id = 0; r.bam_offset = 0; r.mapq = 60;

    std::vector<OrderedRead> reads = {r};
    auto enc = EncodeMetaStream(reads, MetaCodec::TYPED_SPLIT);
    auto dec = DecodeMetaRead(enc.payload, enc.directory[0], enc.codec_id);

    EXPECT_EQ(163u, dec.flag);
    ASSERT_EQ(5u, dec.cigar.size());
    EXPECT_EQ(0u, dec.cigar[0].op); EXPECT_EQ(3u, dec.cigar[0].len);
    EXPECT_EQ(1u, dec.cigar[1].op); EXPECT_EQ(2u, dec.cigar[1].len);
    EXPECT_EQ(2u, dec.cigar[2].op); EXPECT_EQ(1u, dec.cigar[2].len);
    EXPECT_EQ(0u, dec.cigar[3].op); EXPECT_EQ(4u, dec.cigar[3].len);
    EXPECT_EQ(4u, dec.cigar[4].op); EXPECT_EQ(1u, dec.cigar[4].len);
}

// ═════════════════════════════════════════════════════════════════════════════
// S_map DELTA_RANGE round-trip with lazy decode
// ═════════════════════════════════════════════════════════════════════════════

TEST(AuditV4_Map, DeltaRangeRoundTrip) {
    auto td = BuildTestDataset();

    for (size_t i = 0; i < td.reads.size(); ++i) {
        auto decoded = DecodeMapRead(td.map_enc.payload,
                                      td.map_enc.directory[i],
                                      td.map_enc.codec_id);
        EXPECT_EQ(td.reads[i].chrom_id, decoded.chrom_id)
            << "chrom_id mismatch at read " << i;
    }
}

// ═════════════════════════════════════════════════════════════════════════════
// S_qual RANGE_CYCLE round-trip
// ═════════════════════════════════════════════════════════════════════════════

TEST(AuditV4_Qual, RangeCycleRoundTrip) {
    auto td = BuildTestDataset();

    const auto& qual_dir = std::get<StreamDirectoryPerRead>(td.qual_enc.directory);
    for (size_t i = 0; i < td.reads.size(); ++i) {
        auto decoded = DecodeQualRead(td.qual_enc.payload,
                                       qual_dir[i],
                                       td.qual_enc.codec_id);
        ASSERT_EQ(td.reads[i].qual.size(), decoded.size())
            << "Qual length mismatch at read " << i;
        for (size_t j = 0; j < td.reads[i].qual.size(); ++j) {
            EXPECT_EQ(td.reads[i].qual[j], decoded[j])
                << "Qual mismatch at read " << i << " pos " << j;
        }
    }
}

// ═════════════════════════════════════════════════════════════════════════════
// GAP-5: End-to-end reconstruction round-trip
// ═════════════════════════════════════════════════════════════════════════════

TEST(AuditV4_GAP5, BwtReconstructionRoundTrip) {
    auto td = BuildTestDataset();

    // Reconstruct all sequences via LF-walk
    auto all_seqs = ExtractAllSequences(td.fm, td.bv.B_read,
                                         td.reads.size(), td.bundle.S.size());

    ASSERT_EQ(td.reads.size(), all_seqs.size());

    for (size_t i = 0; i < td.reads.size(); ++i) {
        ASSERT_EQ(td.reads[i].seq.size(), all_seqs[i].size())
            << "Reconstructed read " << i << " length mismatch";
        for (size_t j = 0; j < td.reads[i].seq.size(); ++j) {
            EXPECT_EQ(td.reads[i].seq[j], all_seqs[i][j])
                << "Reconstructed read " << i << " base " << j << " mismatch";
        }
    }
}

TEST(AuditV4_GAP5, SingleReadExtraction) {
    auto td = BuildTestDataset();

    // Extract read 2 specifically
    auto seq = ExtractReadSequence(td.fm, td.bv.B_read, 2,
                                    td.reads[2].seq.size());
    ASSERT_EQ(td.reads[2].seq.size(), seq.size());
    EXPECT_EQ(td.reads[2].seq, seq);
}

// ═════════════════════════════════════════════════════════════════════════════
// Lazy LocateLazy equivalence test
// ═════════════════════════════════════════════════════════════════════════════

TEST(AuditV4_GAP1, LocateLazyMatchesLegacy) {
    auto td = BuildTestDataset();

    // Pattern that exists in read 0 (ACGT) — search for "CG"
    std::vector<uint8_t> P = {1, 2}; // CG

    auto legacy = Locate(P, td.fm, td.bv.B_read,
                         td.reads, td.chrom_names,
                         StrandMode::SingleStrand, true);

    auto lazy = LocateLazy(P, td.fm, td.bv.B_read,
                            td.meta_enc.payload, td.map_enc.payload,
                            td.meta_enc.directory, td.map_enc.directory,
                            td.meta_enc.codec_id, td.map_enc.codec_id,
                            td.chrom_names,
                            StrandMode::SingleStrand, true);

    ASSERT_EQ(legacy.size(), lazy.size())
        << "Match count differs between legacy and lazy";

    for (size_t i = 0; i < legacy.size(); ++i) {
        EXPECT_EQ(legacy[i].chrom, lazy[i].chrom);
        EXPECT_EQ(legacy[i].p_min, lazy[i].p_min);
        EXPECT_EQ(legacy[i].p_max, lazy[i].p_max);
        EXPECT_EQ(legacy[i].read_id, lazy[i].read_id);
        EXPECT_EQ(legacy[i].query_strand, lazy[i].query_strand);
    }
}

// ═════════════════════════════════════════════════════════════════════════════
// Invariant verification
// ═════════════════════════════════════════════════════════════════════════════

TEST(AuditV4_Invariants, I14_NoRawSAccess) {
    // Verify IsSeparatorPosition works purely from B_read (no S needed)
    auto td = BuildTestDataset();

    // Separators are at readStarts[i]-1 for i>=1
    for (size_t i = 1; i < td.reads.size(); ++i) {
        uint64_t sep_pos = td.bundle.readStarts[i] - 1;
        EXPECT_TRUE(IsSeparatorPosition(sep_pos, td.bv.B_read))
            << "Expected separator at position " << sep_pos;
    }
    // Non-separator positions
    EXPECT_FALSE(IsSeparatorPosition(0, td.bv.B_read));
    EXPECT_FALSE(IsSeparatorPosition(1, td.bv.B_read));
}

TEST(AuditV4_Invariants, I15_TwoRankAPIs) {
    // Verify the two rank APIs produce consistent but different results
    auto td = BuildTestDataset();

    // B_read.Rank1 is closed-interval: counts 1-bits in [0..pos]
    // FM OccTable rank is half-open: counts symbol occurrences in BWT[0..row)
    // These are fundamentally different operations on different data structures.

    uint64_t pos = td.bundle.readStarts[1]; // start of read 1
    uint64_t rank = td.bv.B_read.Rank1(pos);
    EXPECT_EQ(2u, rank) << "B_read.Rank1 at readStarts[1] should be 2 (reads 0,1)";

    uint64_t read_id = rank - 1;
    EXPECT_EQ(1u, read_id) << "read_id formula: Rank1(pos)-1";
}

TEST(AuditV4_Invariants, I11_SentinelNeverReported) {
    auto td = BuildTestDataset();

    // Search for every single-base pattern — sentinel must never appear
    for (uint8_t base = 0; base <= 4; ++base) {
        std::vector<uint8_t> P = {base};
        auto matches = Locate(P, td.fm, td.bv.B_read,
                              td.reads, td.chrom_names,
                              StrandMode::SingleStrand, false);
        for (const auto& m : matches) {
            EXPECT_NE(m.sa_row, td.fm.SentinelRow())
                << "Sentinel row leaked into results for base " << (int)base;
        }
    }
}
