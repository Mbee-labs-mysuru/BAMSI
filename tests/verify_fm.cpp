/*
 * verify_fm.cpp — Quick FM-index correctness check on the synthetic dataset.
 * Ingests the 10-read BAM, builds the pipeline through Stage 5b,
 * and compares FM backward search counts against brute-force
 * linear scan for a set of test patterns.
 *
 * This is the TIER 1 validation check: "FM search on 100 fixed patterns
 * vs. stored expected counts" (Architecture §4.10 TIER 1).
 */
#include "../src/ingest/ingest.hpp"
#include "../src/ordering/ordering.hpp"
#include "../src/seqbuilder/seqbuilder.hpp"
#include "../src/sais/sais.hpp"
#include "../src/fmindex/fmindex.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <cassert>

using namespace bamsix;

/// Brute-force count of pattern in S
static uint64_t BruteForceCount(const std::vector<uint8_t>& S,
                                 const std::vector<uint8_t>& pattern) {
    if (pattern.size() > S.size()) return 0;
    uint64_t count = 0;
    for (size_t i = 0; i <= S.size() - pattern.size(); ++i) {
        bool match = true;
        for (size_t j = 0; j < pattern.size(); ++j) {
            if (S[i + j] != pattern[j]) { match = false; break; }
        }
        if (match) ++count;
    }
    return count;
}

static std::vector<uint8_t> EncodePattern(const std::string& p) {
    std::vector<uint8_t> codes;
    for (char c : p) {
        switch (c) {
            case 'A': codes.push_back(CODE_A); break;
            case 'C': codes.push_back(CODE_C); break;
            case 'G': codes.push_back(CODE_G); break;
            case 'T': codes.push_back(CODE_T); break;
            default:  codes.push_back(CODE_N); break;
        }
    }
    return codes;
}

int main(int argc, char** argv) {
    const char* bam_path = "data/test/synthetic_10reads.bam";
    if (argc > 1) bam_path = argv[1];

    std::cerr << "=== FM-index correctness verification ===\n";

    // Build pipeline through Stage 5b
    auto ingest = IngestBams({bam_path});
    std::cerr << "Ingested " << ingest.reads.size() << " reads\n";

    auto ordering = OrderReads(ingest);
    auto bundle = BuildSequence(ordering.reads);
    std::cerr << "|S| = " << bundle.S.size() << "\n";

    auto sais = ComputeSuffixArray(bundle);
    std::cerr << "|SA| = " << sais.SA.size()
              << ", sentinel_row = " << sais.sentinel_row << "\n";

    FMIndexEngine fm;
    fm.Build(sais.BWT, sais.SA, sais.sentinel_row, 4 /* small step for testing */,
             bundle.S.size());
    std::cerr << "SA samples: " << fm.SASamples().size() << "\n";

    // Test patterns
    std::vector<std::string> patterns = {
        "A", "C", "G", "T",
        "AC", "GT", "AA", "CC", "GG", "TT",
        "ACG", "CGT", "ACGT",
        "AAAA", "CCCC", "GGGG", "TTTT",
        "ACGTACGT",
        "CCGGTTAA",
        "AAAACCCC",
        "GGGGCCCC",
        "AACCGGTT",
        "ACGTACGTACGTACGT",  // full read01/read07
        "NOTFOUND",  // pattern not in data
        "NNNN",  // N codes
    };

    int pass = 0, fail = 0;
    for (const auto& pstr : patterns) {
        auto codes = EncodePattern(pstr);
        uint64_t bf_count = BruteForceCount(bundle.S, codes);
        auto interval = fm.BackwardSearch(codes);
        uint64_t fm_count = interval.size();

        bool ok = (bf_count == fm_count);
        if (ok) {
            pass++;
        } else {
            fail++;
            std::cerr << "FAIL: pattern=\"" << pstr
                      << "\" brute_force=" << bf_count
                      << " fm_count=" << fm_count << "\n";
        }
    }

    // Also test Locate on a known pattern
    {
        auto codes = EncodePattern("ACGT");
        auto interval = fm.BackwardSearch(codes);
        std::cerr << "Locate test: ACGT has " << interval.size() << " occurrences\n";
        for (uint64_t r = interval.lo; r < interval.hi; ++r) {
            if (r == fm.SentinelRow()) continue;
            uint64_t pos = fm.Locate(r);
            std::cerr << "  row " << r << " -> S-pos " << pos << "\n";
            // Verify: S[pos..pos+3] == ACGT
            if (pos + 4 <= bundle.S.size()) {
                bool ok = (bundle.S[pos] == CODE_A && bundle.S[pos+1] == CODE_C &&
                           bundle.S[pos+2] == CODE_G && bundle.S[pos+3] == CODE_T);
                if (!ok) {
                    std::cerr << "    LOCATE FAIL at pos " << pos << ": S[pos..pos+3] = "
                              << (int)bundle.S[pos] << " " << (int)bundle.S[pos+1]
                              << " " << (int)bundle.S[pos+2] << " " << (int)bundle.S[pos+3] << "\n";
                    fail++;
                } else {
                    pass++;
                }
            }
        }
    }

    std::cerr << "\n=== Results: " << pass << " passed, " << fail << " failed ===\n";
    return (fail > 0) ? 1 : 0;
}
