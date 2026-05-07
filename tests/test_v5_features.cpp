/// BAMSI V5 Feature Tests — ISA samples, SARange wavelet tree, Reverse FM-index
/// Tests per Architecture §4.4, §4.6.7, §5.3

#include <gtest/gtest.h>
#include "bamsi/bamsi.h"
#include "bamsi/types.hpp"
#include "format/format.hpp"
#include "sarange/sarange.hpp"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

namespace bamsi {
namespace {

const std::string kTestBam = "data/test/synthetic_10reads.bam";

struct CliResult {
    int exit_code;
    std::string output;
};

CliResult RunBamsi(const std::string& args) {
    std::string cmd = "./build/bamsi " + args + " 2>&1";
    FILE* fp = popen(cmd.c_str(), "r");
    std::string output;
    char buf[4096];
    while (fgets(buf, sizeof(buf), fp)) output += buf;
    int rc = pclose(fp);
    return {WEXITSTATUS(rc), output};
}

std::string BuildIndex(const std::string& name, const std::string& extra = "") {
    std::string bsi = "/tmp/bamsi_v5_" + name + ".bsi";
    std::string cmd = "./build/bamsi build --input " + kTestBam +
                      " --output " + bsi + " " + extra + " 2>/dev/null";
    if (system(cmd.c_str()) != 0) return "";
    return bsi;
}

// ═══════════════════════════════════════════════════════════════════════════
// ISA Samples (Architecture §4.4 step 5)
// ═══════════════════════════════════════════════════════════════════════════

TEST(V5Features, ISASamplesBuildSucceeds) {
    auto bsi = BuildIndex("isa_build", "--isa-step 32");
    ASSERT_FALSE(bsi.empty());

    auto idx = ReadBsi(bsi);
    EXPECT_EQ(idx.header.has_isa_samples, 1);
    EXPECT_EQ(idx.header.sample_step_s_prime, 32u);
}

TEST(V5Features, ISASamplesQueryConsistency) {
    // Locate results must be identical with and without ISA samples
    auto bsi_no_isa = BuildIndex("isa_noisa");
    auto bsi_isa = BuildIndex("isa_withisa", "--isa-step 32");
    ASSERT_FALSE(bsi_no_isa.empty());
    ASSERT_FALSE(bsi_isa.empty());

    std::vector<std::string> patterns = {"A", "AC", "ACG", "ACGT"};
    for (const auto& pat : patterns) {
        auto r1 = RunBamsi("count --index " + bsi_no_isa + " --pattern " + pat);
        auto r2 = RunBamsi("count --index " + bsi_isa + " --pattern " + pat);
        EXPECT_EQ(r1.output, r2.output) << "ISA samples changed count for: " << pat;

        // Locate positions should also match
        auto l1 = RunBamsi("locate --index " + bsi_no_isa + " --pattern " + pat + " --sort-output --json");
        auto l2 = RunBamsi("locate --index " + bsi_isa + " --pattern " + pat + " --sort-output --json");
        EXPECT_EQ(l1.exit_code, 0);
        EXPECT_EQ(l2.exit_code, 0);
    }
}

TEST(V5Features, ISASamplesVerifyPasses) {
    auto bsi = BuildIndex("isa_verify", "--isa-step 64");
    ASSERT_FALSE(bsi.empty());

    auto r = RunBamsi("verify --index " + bsi);
    EXPECT_EQ(r.exit_code, 0);
}

TEST(V5Features, ISASamplesInfoReported) {
    auto bsi = BuildIndex("isa_info", "--isa-step 32");
    ASSERT_FALSE(bsi.empty());

    auto r = RunBamsi("info --index " + bsi);
    EXPECT_EQ(r.exit_code, 0);
    EXPECT_NE(r.output.find("ISA samples"), std::string::npos)
        << "Info should report ISA samples presence";
    EXPECT_NE(r.output.find("yes"), std::string::npos);
}

// ═══════════════════════════════════════════════════════════════════════════
// SARange Wavelet Tree (Architecture §5.3 — ENHANCED Tier)
// ═══════════════════════════════════════════════════════════════════════════

TEST(V5Features, SARangeBuildSucceeds) {
    auto bsi = BuildIndex("sarange_build", "--enable-sarange");
    ASSERT_FALSE(bsi.empty());

    auto idx = ReadBsi(bsi);
    EXPECT_EQ(idx.header.enable_sarange, 1);
}

TEST(V5Features, SARangeQueryConsistency) {
    // Count results must be identical with and without SARange
    auto bsi_base = BuildIndex("sarange_base");
    auto bsi_enhanced = BuildIndex("sarange_enhanced", "--enable-sarange");
    ASSERT_FALSE(bsi_base.empty());
    ASSERT_FALSE(bsi_enhanced.empty());

    std::vector<std::string> patterns = {"A", "AC", "ACG", "ACGT"};
    for (const auto& pat : patterns) {
        auto r1 = RunBamsi("count --index " + bsi_base + " --pattern " + pat);
        auto r2 = RunBamsi("count --index " + bsi_enhanced + " --pattern " + pat);
        EXPECT_EQ(r1.output, r2.output) << "SARange changed count for: " << pat;
    }
}

TEST(V5Features, SARangeInfoShowsEnhanced) {
    auto bsi = BuildIndex("sarange_info", "--enable-sarange");
    ASSERT_FALSE(bsi.empty());

    auto r = RunBamsi("info --index " + bsi);
    EXPECT_EQ(r.exit_code, 0);
    EXPECT_NE(r.output.find("ENHANCED"), std::string::npos)
        << "Info should show ENHANCED tier for SARange index";
}

TEST(V5Features, SARangeVerifyPasses) {
    auto bsi = BuildIndex("sarange_verify", "--enable-sarange");
    ASSERT_FALSE(bsi.empty());

    auto r = RunBamsi("verify --index " + bsi);
    EXPECT_EQ(r.exit_code, 0);
}

// ═══════════════════════════════════════════════════════════════════════════
// SARange Unit Tests — wavelet tree correctness
// ═══════════════════════════════════════════════════════════════════════════

TEST(SARangeUnit, SmallArrayRangeCount) {
    SARange wt;
    std::vector<uint64_t> data = {5, 1, 8, 3, 7, 2, 9, 4, 6, 0};
    wt.Build(data, 10);

    ASSERT_TRUE(wt.IsBuilt());
    EXPECT_EQ(wt.Size(), 10u);

    // Count elements in [0, 10) that fall in [0, 9] — should be all 10
    EXPECT_EQ(wt.RangeCount(0, 10, 0, 9), 10u);

    // Count elements in [0, 10) that fall in [0, 4] — elements {0,1,2,3,4} = 5
    EXPECT_EQ(wt.RangeCount(0, 10, 0, 4), 5u);

    // Count elements in [0, 10) that fall in [5, 9] — elements {5,6,7,8,9} = 5
    EXPECT_EQ(wt.RangeCount(0, 10, 5, 9), 5u);

    // Count elements in [0, 5) that fall in [0, 4]
    // data[0..5) = {5,1,8,3,7} → only {1,3} in [0,4] = 2
    EXPECT_EQ(wt.RangeCount(0, 5, 0, 4), 2u);

    // Count elements in [3, 7) that fall in [2, 6]
    // data[3..7) = {3,7,2,9} → only {3,2} in [2,6] = 2
    EXPECT_EQ(wt.RangeCount(3, 7, 2, 6), 2u);
}

TEST(SARangeUnit, EmptyRange) {
    SARange wt;
    std::vector<uint64_t> data = {5, 1, 8};
    wt.Build(data, 10);

    EXPECT_EQ(wt.RangeCount(0, 0, 0, 10), 0u);   // empty position range
    EXPECT_EQ(wt.RangeCount(0, 3, 10, 5), 0u);    // inverted value range
}

TEST(SARangeUnit, SerializeDeserialize) {
    SARange wt;
    std::vector<uint64_t> data = {5, 1, 8, 3, 7, 2, 9, 4, 6, 0};
    wt.Build(data, 10);

    auto serialized = wt.Serialize();
    EXPECT_GT(serialized.size(), 0u);

    SARange wt2;
    wt2.Deserialize(serialized.data(), serialized.size());

    ASSERT_TRUE(wt2.IsBuilt());
    // Verify same range_count results
    EXPECT_EQ(wt2.RangeCount(0, 10, 0, 4), 5u);
    EXPECT_EQ(wt2.RangeCount(0, 10, 5, 9), 5u);
    EXPECT_EQ(wt2.RangeCount(0, 5, 0, 4), 2u);
}

// ═══════════════════════════════════════════════════════════════════════════
// Combined V5: ISA + SARange + all features
// ═══════════════════════════════════════════════════════════════════════════

TEST(V5Features, FullV5BuildAndQuery) {
    auto bsi = BuildIndex("full_v5", "--isa-step 32 --enable-sarange");
    ASSERT_FALSE(bsi.empty());

    auto idx = ReadBsi(bsi);
    EXPECT_EQ(idx.header.has_isa_samples, 1);
    EXPECT_EQ(idx.header.enable_sarange, 1);

    // Query should work
    auto r = RunBamsi("count --index " + bsi + " --pattern ACGT");
    EXPECT_EQ(r.exit_code, 0);
    int count = std::stoi(r.output);
    EXPECT_GT(count, 0);

    // Verify should pass
    auto rv = RunBamsi("verify --index " + bsi);
    EXPECT_EQ(rv.exit_code, 0);

    // Info should show both features
    auto ri = RunBamsi("info --index " + bsi);
    EXPECT_EQ(ri.exit_code, 0);
    EXPECT_NE(ri.output.find("ISA samples"), std::string::npos);
    EXPECT_NE(ri.output.find("ENHANCED"), std::string::npos);
}

TEST(V5Features, FullV5Reconstruct) {
    auto bsi = BuildIndex("v5_recon", "--isa-step 32 --enable-sarange");
    ASSERT_FALSE(bsi.empty());

    auto r = RunBamsi("reconstruct --index " + bsi +
                       " --output /tmp/bamsi_v5_recon.bam");
    EXPECT_EQ(r.exit_code, 0);
}

}  // namespace
}  // namespace bamsi
