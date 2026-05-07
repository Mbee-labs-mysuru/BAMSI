/// Validation Module — TIER 1 Invariant Tests
/// Tests invariants I1–I15 per Architecture §9.1/§9.3.
///
/// Usage: validate <input.bam> [--verbose]
/// Tests the full pipeline and validates all structural invariants.

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
#include "../src/mapping/mapping.hpp"
#include "bamsi/config.hpp"

#include <cstring>
#include <iostream>
#include <set>
#include <algorithm>

using namespace bamsi;

static int pass = 0, fail = 0;

void CHECK(bool cond, const char* invariant, const char* detail) {
    if (cond) {
        ++pass;
    } else {
        std::cerr << "FAIL [" << invariant << "]: " << detail << "\n";
        ++fail;
    }
}

int main(int argc, char** argv) {
    const char* bam = (argc > 1) ? argv[1] : "data/test/synthetic_10reads.bam";
    bool verbose = (argc > 2 && std::string(argv[2]) == "--verbose");

    std::cerr << "=== BAMSI TIER 1 Invariant Validation ===\n";
    std::cerr << "Input: " << bam << "\n\n";

    // ─── Build pipeline ─────────────────────────────────────────────────────
    auto ingest = IngestBams({bam});
    auto ordering = OrderReads(ingest);
    auto bundle = BuildSequence(ordering.reads);
    auto sais_r = ComputeSuffixArray(bundle);
    FMIndexEngine fm;
    fm.Build(sais_r.BWT, sais_r.SA, sais_r.sentinel_row, 64, bundle.S.size());
    auto windows = BuildWindows(ordering.reads, bundle, 100000);
    auto bv = BuildBitvectors(bundle, windows);

    const auto& S = bundle.S;
    const auto& reads = ordering.reads;
    uint64_t N = reads.size();
    uint64_t n = S.size();

    // ─── I1: Read collection ────────────────────────────────────────────────
    CHECK(N > 0, "I1", "Non-empty read collection");
    CHECK(ingest.reads.size() == N, "I1", "All reads pass inclusion rule");

    // ─── I2: Total ordering determinism ─────────────────────────────────────
    {
        bool ordered = true;
        for (uint64_t i = 1; i < N; ++i) {
            const auto& a = reads[i-1];
            const auto& b = reads[i];
            if (a.chrom_id > b.chrom_id) { ordered = false; break; }
            if (a.chrom_id == b.chrom_id && a.pos > b.pos) { ordered = false; break; }
        }
        CHECK(ordered, "I2", "Reads are in total order (chrom_id, pos, ...)");
    }

    // ─── I3: Ordering hash determinism ──────────────────────────────────────
    {
        auto ordering2 = OrderReads(ingest);
        CHECK(ordering.ordering_hash == ordering2.ordering_hash, "I3",
              "Ordering hash is deterministic on same input");
    }

    // ─── I4: S construction ─────────────────────────────────────────────────
    {
        // S = r0 # r1 # ... # r_{N-1}
        // Check separators
        bool seps_ok = true;
        for (uint64_t i = 0; i < N; ++i) {
            uint64_t start = bundle.readStarts[i];
            uint64_t len = reads[i].seq.size();
            // Check read content
            for (uint64_t j = 0; j < len; ++j) {
                if (S[start + j] != reads[i].seq[j]) { seps_ok = false; break; }
            }
            // Check separator (except after last read)
            if (i < N - 1) {
                if (S[start + len] != CODE_SEP) { seps_ok = false; break; }
            }
        }
        CHECK(seps_ok, "I4", "S correctly concatenates reads with # separators");
    }

    // ─── I5: SA correctness ─────────────────────────────────────────────────
    {
        // SA must be a permutation of 0..n
        std::vector<bool> seen(n + 1, false);
        bool perm_ok = true;
        for (auto v : sais_r.SA) {
            if (v < 0 || (uint64_t)v > n || seen[v]) { perm_ok = false; break; }
            seen[v] = true;
        }
        CHECK(perm_ok, "I5", "SA is a valid permutation of 0..|S|");
    }

    // ─── I6: BWT derivation ─────────────────────────────────────────────────
    {
        // BWT[i] = S$[(SA[i] - 1) mod (n+1)]
        // where S$ = S + '$' (conceptual)
        bool bwt_ok = true;
        for (uint64_t i = 0; i < sais_r.BWT.size(); ++i) {
            int64_t sa_val = sais_r.SA[i];
            uint64_t pred_pos = (sa_val == 0) ? n : (sa_val - 1);
            uint8_t expected;
            if (pred_pos == n) {
                // Position n is the sentinel '$' position, but we use CODE_SENT
                expected = CODE_SENT;
            } else {
                expected = S[pred_pos];
            }
            // Special case: sentinel_row BWT is S[n-1]
            if (i == sais_r.sentinel_row) {
                expected = S[n - 1];
            }
            if (sais_r.BWT[i] != expected) {
                bwt_ok = false;
                if (verbose) {
                    std::cerr << "  BWT[" << i << "] = " << (int)sais_r.BWT[i]
                              << " expected " << (int)expected
                              << " (SA=" << sa_val << ")\n";
                }
            }
        }
        CHECK(bwt_ok, "I6", "BWT derivation correct for all rows");
    }

    // ─── I7: LF property ────────────────────────────────────────────────────
    {
        uint64_t lf_errors = 0;
        for (uint64_t r = 0; r < sais_r.BWT.size(); ++r) {
            uint64_t lf_r = fm.LF(r);
            int64_t sa_r = sais_r.SA[r];
            int64_t sa_lf = sais_r.SA[lf_r];
            int64_t expected = ((sa_r - 1) + (int64_t)(n + 1)) % (int64_t)(n + 1);
            if (sa_lf != expected) ++lf_errors;
        }
        CHECK(lf_errors == 0, "I7", "LF property SA[LF(r)] = (SA[r]-1) mod (|S|+1)");
    }

    // ─── I8: Backward search correctness ────────────────────────────────────
    {
        // Test a few patterns against brute-force
        std::vector<std::vector<uint8_t>> test_pats = {
            {0,1,2,3},      // ACGT
            {1,1,2,2},      // CCGG
            {3,3,0,0},      // TTAA
            {0,0,0,0},      // AAAA
            {4,4,4},        // NNN
        };
        bool bs_ok = true;
        for (auto& p : test_pats) {
            auto interval = fm.BackwardSearch(p);
            uint64_t fm_count = interval.size();
            // Brute force
            uint64_t bf = 0;
            for (uint64_t i = 0; i + p.size() <= n; ++i) {
                bool match = true;
                for (uint64_t j = 0; j < p.size(); ++j) {
                    if (S[i + j] != p[j]) { match = false; break; }
                }
                if (match) ++bf;
            }
            if (fm_count != bf) { bs_ok = false; }
        }
        CHECK(bs_ok, "I8", "Backward search matches brute-force for test patterns");
    }

    // ─── I9: B_read popcount ────────────────────────────────────────────────
    CHECK(bv.B_read.PopCount() == N, "I9",
          "B_read has exactly N 1-bits (one per read)");

    // ─── I10: B_read positions match readStarts ─────────────────────────────
    {
        bool starts_ok = true;
        for (uint64_t i = 0; i < N; ++i) {
            uint64_t sel = bv.B_read.Select1(i + 1);
            if (sel != bundle.readStarts[i]) { starts_ok = false; break; }
        }
        CHECK(starts_ok, "I10", "B_read.select1(i+1) == readStarts[i] for all i");
    }

    // ─── I11: Window coverage ───────────────────────────────────────────────
    {
        // Every read must be in at least one window
        std::vector<bool> covered(N, false);
        for (const auto& w : windows) {
            for (uint64_t rid = w.first_read_id; rid <= w.last_read_id; ++rid) {
                if (rid < N) covered[rid] = true;
            }
        }
        bool all_covered = true;
        for (bool c : covered) { if (!c) all_covered = false; }
        CHECK(all_covered, "I11", "All reads covered by at least one window");
    }

    // ─── I12: Window ordering ───────────────────────────────────────────────
    {
        bool win_ordered = true;
        for (size_t i = 1; i < windows.size(); ++i) {
            if (windows[i].chrom_id < windows[i-1].chrom_id) {
                win_ordered = false; break;
            }
            if (windows[i].chrom_id == windows[i-1].chrom_id &&
                windows[i].l <= windows[i-1].l) {
                win_ordered = false; break;
            }
        }
        CHECK(win_ordered, "I12", "Windows are in (chrom_id, l) order");
    }

    // ─── I13: SA sample correctness ─────────────────────────────────────────
    {
        bool samples_ok = true;
        for (uint64_t k = 0; k < fm.SASamples().size(); ++k) {
            uint64_t row = k * 64;  // sample_step = 64
            if (row < sais_r.SA.size()) {
                if (fm.SASamples()[k] != (uint64_t)sais_r.SA[row]) {
                    samples_ok = false; break;
                }
            }
        }
        CHECK(samples_ok, "I13", "SA samples match SA[k*s] for all k");
    }

    // ─── I14: Separator detection without raw S ─────────────────────────────
    {
        bool sep_ok = true;
        for (uint64_t pos = 0; pos < n; ++pos) {
            bool is_sep = IsSeparatorPosition(pos, bv.B_read);
            bool actual_sep = (S[pos] == CODE_SEP);
            if (is_sep != actual_sep) {
                sep_ok = false;
                if (verbose) {
                    std::cerr << "  I14: pos=" << pos << " is_sep=" << is_sep
                              << " actual=" << actual_sep << "\n";
                }
            }
        }
        CHECK(sep_ok, "I14", "IsSeparatorPosition matches S[pos]==# for all pos");
    }

    // ─── I15: Locate correctness ────────────────────────────────────────────
    {
        bool locate_ok = true;
        for (uint64_t r = 0; r < sais_r.BWT.size(); ++r) {
            if (r == fm.SentinelRow()) continue;
            uint64_t pos = fm.Locate(r);
            if (pos != (uint64_t)sais_r.SA[r]) {
                locate_ok = false;
                if (verbose) {
                    std::cerr << "  I15: row " << r << " Locate=" << pos
                              << " SA=" << sais_r.SA[r] << "\n";
                }
            }
        }
        CHECK(locate_ok, "I15", "Locate(r) == SA[r] for all non-sentinel rows");
    }

    // ─── Seal + Verify round-trip ───────────────────────────────────────────
    {
        const std::string bsi_path = "data/test/invariant_test.bsi";

        BsiHeader header;
        std::memcpy(header.magic, "BMSI", 4);
        header.version = 6;
        std::strncpy(header.bamsi_version, BAMSI_VERSION, 15);
        header.is_lossless = 1;
        header.source_file_count = ingest.source_file_count;
        header.source_manifest_hash = ingest.source_manifest_hash;
        header.ordering_hash = ordering.ordering_hash;
        header.S_length = n;
        header.N_reads = N;
        header.N_windows = windows.size();
        header.sample_step_s = 64;
        header.sentinel_row = fm.SentinelRow();
        header.chrom_count = ingest.chrom_names.size();
        for (uint32_t i = 0; i < ingest.chrom_names.size(); ++i)
            header.chrom_name_table.push_back({i, ingest.chrom_names[i]});

        auto seq = EncodeSeqStream(sais_r.BWT, 6);
        auto qual = EncodeQualStream(ordering.reads, QualCodec::RANGE_CYCLE, 0);
        auto meta = EncodeMetaStream(ordering.reads, MetaCodec::TYPED_SPLIT);
        auto map_ = EncodeMapStream(ordering.reads, MapCodec::DELTA_RANGE);

        WriteBsi(bsi_path, {header, seq, qual, meta, map_, fm, bv, windows, ordering.reads});
        CHECK(VerifyBsi(bsi_path), "SEAL", "Global xxHash64 checksum verified");

        auto idx = ReadBsi(bsi_path);
        CHECK(idx.header.N_reads == N, "SEAL", "Round-trip: N_reads preserved");
        CHECK(idx.reads.size() == N, "SEAL", "Round-trip: reads metadata preserved");

        // Verify FM-index round-trip
        uint64_t count_orig = GlobalCount({0,1,2,3}, fm, StrandMode::SingleStrand);
        uint64_t count_loaded = GlobalCount({0,1,2,3}, idx.fm, StrandMode::SingleStrand);
        CHECK(count_orig == count_loaded, "SEAL", "Round-trip: GlobalCount preserved");
    }

    // ─── Summary ────────────────────────────────────────────────────────────
    std::cerr << "\n=== TIER 1 Validation: "
              << pass << " passed, " << fail << " failed ===\n";
    return fail;
}
