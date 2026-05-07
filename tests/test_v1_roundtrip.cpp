/// V1 Gate Test: build → load → count → locate → verify
/// Tests the full round-trip pipeline on the synthetic 10-read BAM.

#include "../src/ingest/ingest.hpp"
#include "../src/ordering/ordering.hpp"
#include "../src/seqbuilder/seqbuilder.hpp"
#include "../src/sais/sais.hpp"
#include "../src/fmindex/fmindex.hpp"
#include "../src/streamencode/streamencode.hpp"
#include "../src/windows/windows.hpp"
#include "../src/bitvectors/bitvectors.hpp"
#include "../src/seal/seal.hpp"
#include "../src/format/format.hpp"
#include "../src/query/query.hpp"
#include "bamsi/config.hpp"

#include <cstring>
#include <iostream>
#include <vector>

using namespace bamsi;

static int failures = 0;

void CHECK(bool cond, const char* msg) {
    if (!cond) {
        std::cerr << "FAIL: " << msg << "\n";
        ++failures;
    }
}

int main() {
    const std::string bam_path = "data/test/synthetic_10reads.bam";
    const std::string bsi_path = "data/test/v1_gate_test.bsi";

    std::cerr << "=== V1 Gate Test ===\n\n";

    // ─── Step 1: Build .bsi ─────────────────────────────────────────────────
    std::cerr << "[1] Building .bsi from " << bam_path << " ...\n";
    {
        auto ingest = IngestBams({bam_path});
        auto ordering = OrderReads(ingest);
        auto bundle = BuildSequence(ordering.reads);
        auto sais = ComputeSuffixArray(bundle);

        FMIndexEngine fm;
        fm.Build(sais.BWT, sais.SA, sais.sentinel_row, 64, bundle.S.size());

        auto seq = EncodeSeqStream(sais.BWT, 6);
        auto qual = EncodeQualStream(ordering.reads, QualCodec::RANGE_CYCLE, 0);
        auto meta = EncodeMetaStream(ordering.reads, MetaCodec::TYPED_SPLIT);
        auto map_ = EncodeMapStream(ordering.reads, MapCodec::DELTA_RANGE);
        auto windows = BuildWindows(ordering.reads, bundle, 100000);
        auto bv = BuildBitvectors(bundle, windows);

        BsiHeader header;
        std::memcpy(header.magic, "BMSI", 4);
        header.version = 6;
        std::strncpy(header.bamsi_version, BAMSI_VERSION, 15);
        header.is_lossless = 1;
        header.source_file_count = ingest.source_file_count;
        header.source_manifest_hash = ingest.source_manifest_hash;
        header.ordering_hash = ordering.ordering_hash;
        header.S_length = bundle.S.size();
        header.N_reads = ordering.reads.size();
        header.N_windows = windows.size();
        header.sample_step_s = 64;
        header.sentinel_row = fm.SentinelRow();
        header.chrom_count = ingest.chrom_names.size();
        for (uint32_t i = 0; i < ingest.chrom_names.size(); ++i)
            header.chrom_name_table.push_back({i, ingest.chrom_names[i]});

        WriteBsi(bsi_path, {header, seq, qual, meta, map_, fm, bv, windows});
        std::cerr << "  Built: " << bsi_path << "\n";
    }

    // ─── Step 2: Verify checksum ────────────────────────────────────────────
    std::cerr << "[2] Verifying checksum ...\n";
    CHECK(VerifyBsi(bsi_path), "VerifyBsi should return true");

    // ─── Step 3: Load .bsi ──────────────────────────────────────────────────
    std::cerr << "[3] Loading .bsi ...\n";
    auto idx = ReadBsi(bsi_path);

    CHECK(std::memcmp(idx.header.magic, "BMSI", 4) == 0, "magic == BMSI");
    CHECK(idx.header.version == 6, "version == 6");
    CHECK(idx.header.N_reads == 10, "N_reads == 10");
    CHECK(idx.header.S_length == 149, "S_length == 149");
    CHECK(idx.header.chrom_count == 1, "chrom_count == 1");
    CHECK(idx.chrom_names.size() == 1, "chrom_names.size == 1");
    CHECK(idx.chrom_names[0] == "chr1", "chrom_names[0] == chr1");
    CHECK(idx.header.sentinel_row == 0, "sentinel_row == 0");
    CHECK(idx.fm.SASamples().size() == 3, "SA_samples.size == 3");

    // ─── Step 4: GlobalCount checks ─────────────────────────────────────────
    std::cerr << "[4] Testing GlobalCount ...\n";
    StrandMode mode = StrandMode::StrandComplete;

    // ACGT occurs 13 times (fwd), ACGT is its own RC so no extra
    auto count_acgt = GlobalCount({0,1,2,3}, idx.fm, mode);
    CHECK(count_acgt == 13, "GlobalCount(ACGT) == 13");

    // Single-strand CCGGTTAA = 3
    auto count_ccgg_ss = GlobalCount({1,1,2,2,3,3,0,0}, idx.fm, StrandMode::SingleStrand);
    CHECK(count_ccgg_ss == 3, "GlobalCount(CCGGTTAA, single-strand) == 3");

    // Strand-complete CCGGTTAA = 3 + 2 = 5
    auto count_ccgg_sc = GlobalCount({1,1,2,2,3,3,0,0}, idx.fm, StrandMode::StrandComplete);
    CHECK(count_ccgg_sc == 5, "GlobalCount(CCGGTTAA, strand-complete) == 5");

    // Non-existent pattern
    auto count_zero = GlobalCount({0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, idx.fm, mode);
    // 25 consecutive A's — unlikely in 10 short reads
    // (actually check: longest A run is 4 in "AAAA")
    CHECK(count_zero == 0, "GlobalCount(25xA) == 0");

    // ─── Step 5: GlobalExists checks ────────────────────────────────────────
    std::cerr << "[5] Testing GlobalExists ...\n";
    CHECK(GlobalExists({0,1,2,3}, idx.fm, mode), "GlobalExists(ACGT) == true");
    CHECK(!GlobalExists({4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4}, idx.fm, mode),
          "GlobalExists(16xN) == false");

    // ─── Step 6: Locate S-positions check ───────────────────────────────────
    std::cerr << "[6] Testing Locate (S-positions) ...\n";
    auto interval = idx.fm.BackwardSearch({0,1,2,3}); // ACGT
    CHECK(interval.size() == 13, "ACGT interval size == 13");

    std::vector<uint64_t> positions;
    for (uint64_t row = interval.lo; row < interval.hi; ++row) {
        if (row == idx.fm.SentinelRow()) continue;
        positions.push_back(idx.fm.Locate(row));
    }

    // All positions should be in [0, 148]
    for (auto p : positions) {
        CHECK(p <= 148, "Locate position in valid range");
    }

    // Position 0 should be present (first read starts with ACGT)
    bool has_pos0 = false;
    for (auto p : positions) { if (p == 0) has_pos0 = true; }
    CHECK(has_pos0, "Locate finds ACGT at position 0");

    // ─── Results ────────────────────────────────────────────────────────────
    std::cerr << "\n=== V1 Gate Test: " << (failures == 0 ? "PASSED" : "FAILED")
              << " (" << failures << " failures) ===\n";
    return failures;
}
