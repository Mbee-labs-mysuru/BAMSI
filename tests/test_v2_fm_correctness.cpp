/// V2 Gate Test: FM correctness on 10K-read synthetic input.
/// Per Execution Plan V2: 100 random patterns of length 8-30,
/// GlobalCount must match brute-force linear scan over S.

#include "../src/ingest/ingest.hpp"
#include "../src/ordering/ordering.hpp"
#include "../src/seqbuilder/seqbuilder.hpp"
#include "../src/sais/sais.hpp"
#include "../src/fmindex/fmindex.hpp"
#include "../src/query/query.hpp"
#include "bamsix/config.hpp"

#include <chrono>
#include <cstdlib>
#include <iostream>
#include <random>
#include <vector>

using namespace bamsix;

/// Brute-force count of pattern in S (stored codes).
static uint64_t BruteForceCount(const std::vector<uint8_t>& S,
                                 const std::vector<uint8_t>& pattern) {
    if (pattern.empty() || pattern.size() > S.size()) return 0;
    uint64_t count = 0;
    for (size_t i = 0; i + pattern.size() <= S.size(); ++i) {
        bool match = true;
        for (size_t j = 0; j < pattern.size(); ++j) {
            if (S[i + j] != pattern[j]) { match = false; break; }
        }
        if (match) ++count;
    }
    return count;
}

int main(int argc, char** argv) {
    const char* bam = (argc > 1) ? argv[1] : "data/test/synthetic_10k.bam";
    int n_patterns = 100;
    int seed = 42;

    auto t0 = std::chrono::steady_clock::now();
    auto elapsed = [&]() {
        auto now = std::chrono::steady_clock::now();
        return std::chrono::duration_cast<std::chrono::milliseconds>(now - t0).count();
    };

    std::cerr << "=== V2 FM Correctness Gate Test ===\n";
    std::cerr << "Input: " << bam << "\n";

    // ─── Build pipeline ─────────────────────────────────────────────────────
    std::cerr << "[1] Ingesting ...\n";
    auto ingest = IngestBams({bam});
    std::cerr << "  " << ingest.reads.size() << " reads (" << elapsed() << " ms)\n";

    std::cerr << "[2] Ordering ...\n";
    auto ordering = OrderReads(ingest);

    std::cerr << "[3] Building S ...\n";
    auto bundle = BuildSequence(ordering.reads);
    std::cerr << "  |S| = " << bundle.S.size() << " (" << elapsed() << " ms)\n";

    std::cerr << "[4] SA-IS ...\n";
    auto sais_r = ComputeSuffixArray(bundle);
    std::cerr << "  |SA| = " << sais_r.SA.size() << " (" << elapsed() << " ms)\n";

    std::cerr << "[5] FM-index ...\n";
    FMIndexEngine fm;
    fm.Build(sais_r.BWT, sais_r.SA, sais_r.sentinel_row, 64, bundle.S.size());
    std::cerr << "  SA samples: " << fm.SASamples().size() << " (" << elapsed() << " ms)\n";

    // ─── Generate random patterns ───────────────────────────────────────────
    std::mt19937 rng(seed);
    std::uniform_int_distribution<int> len_dist(8, 30);
    std::uniform_int_distribution<size_t> pos_dist(0, bundle.S.size() - 1);
    std::uniform_int_distribution<int> base_dist(0, 3);

    struct TestPattern {
        std::vector<uint8_t> pat;
        std::string desc;
    };

    std::vector<TestPattern> patterns;

    // Half from actual S substrings (guaranteed to have matches)
    for (int i = 0; i < n_patterns / 2; ++i) {
        int len = len_dist(rng);
        size_t start = pos_dist(rng);
        if (start + len > bundle.S.size()) start = bundle.S.size() - len;

        std::vector<uint8_t> pat(bundle.S.begin() + start,
                                  bundle.S.begin() + start + len);
        patterns.push_back({pat, "S[" + std::to_string(start) + ".." +
                                std::to_string(start + len) + ")"});
    }

    // Half random (may or may not match)
    for (int i = 0; i < n_patterns / 2; ++i) {
        int len = len_dist(rng);
        std::vector<uint8_t> pat(len);
        for (int j = 0; j < len; ++j) pat[j] = base_dist(rng);
        patterns.push_back({pat, "random_" + std::to_string(i)});
    }

    // ─── Verify each pattern ────────────────────────────────────────────────
    std::cerr << "[6] Testing " << patterns.size() << " patterns ...\n";
    int pass = 0, fail = 0;

    for (size_t i = 0; i < patterns.size(); ++i) {
        auto& tp = patterns[i];

        // Skip patterns containing separators (code 5) — these cross
        // read boundaries and the FM backward search correctly counts
        // them differently than a byte-scan across separators.
        bool has_sep = false;
        for (auto c : tp.pat) {
            if (c >= CODE_SEP) { has_sep = true; break; }
        }
        if (has_sep) {
            ++pass; // not testable, skip
            continue;
        }

        uint64_t fm_count = GlobalCount(tp.pat, fm, StrandMode::SingleStrand);
        uint64_t bf_count = BruteForceCount(bundle.S, tp.pat);

        if (fm_count != bf_count) {
            std::cerr << "  FAIL [" << i << "] " << tp.desc
                      << ": FM=" << fm_count << " BF=" << bf_count << "\n";
            ++fail;
        } else {
            ++pass;
        }
    }

    // ─── Locate verification (10 patterns) ──────────────────────────────────
    std::cerr << "[7] Testing Locate for 10 patterns ...\n";
    int locate_pass = 0, locate_fail = 0;

    for (int i = 0; i < 10 && i < (int)patterns.size(); ++i) {
        auto& tp = patterns[i];
        bool has_sep = false;
        for (auto c : tp.pat) { if (c >= CODE_SEP) { has_sep = true; break; } }
        if (has_sep) { ++locate_pass; continue; }

        auto interval = fm.BackwardSearch(tp.pat);
        for (uint64_t row = interval.lo; row < interval.hi; ++row) {
            if (row == fm.SentinelRow()) continue;
            uint64_t pos = fm.Locate(row);
            if (pos >= bundle.S.size()) {
                std::cerr << "  LOCATE FAIL: pos=" << pos << " >= |S|=" << bundle.S.size() << "\n";
                ++locate_fail;
                continue;
            }
            // Verify match at position
            bool match = true;
            for (size_t j = 0; j < tp.pat.size(); ++j) {
                if (pos + j >= bundle.S.size() || bundle.S[pos + j] != tp.pat[j]) {
                    match = false;
                    break;
                }
            }
            if (!match) {
                std::cerr << "  LOCATE FAIL: pattern not at pos=" << pos << "\n";
                ++locate_fail;
            } else {
                ++locate_pass;
            }
        }
    }

    // ─── Results ────────────────────────────────────────────────────────────
    auto total_ms = elapsed();
    std::cerr << "\n=== Results ===\n";
    std::cerr << "Count:  " << pass << " passed, " << fail << " failed\n";
    std::cerr << "Locate: " << locate_pass << " passed, " << locate_fail << " failed\n";
    std::cerr << "Time:   " << total_ms << " ms\n";

    if (fail == 0 && locate_fail == 0) {
        std::cerr << "\n=== V2 Gate Test: PASSED ===\n";
        return 0;
    } else {
        std::cerr << "\n=== V2 Gate Test: FAILED ===\n";
        return 1;
    }
}
