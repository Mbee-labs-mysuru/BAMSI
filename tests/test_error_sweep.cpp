/// BAMSIX Error Handling Sweep — Contract §0.8 compliance
/// Tests every ErrorCode path via C ABI + CLI subprocess.
/// Ensures "no silent failure" for all documented error conditions.

#include <gtest/gtest.h>
#include "bamsix/bamsix.h"

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>

namespace bamsix {
namespace {

const std::string kTestBam = "data/test/synthetic_10reads.bam";

struct CliResult {
    int exit_code;
    std::string output;
};

CliResult RunBamsix(const std::string& args) {
    std::string cmd = "./build/bamsix " + args + " 2>&1";
    FILE* fp = popen(cmd.c_str(), "r");
    std::string output;
    char buf[4096];
    while (fgets(buf, sizeof(buf), fp)) output += buf;
    int rc = pclose(fp);
    return {WEXITSTATUS(rc), output};
}

std::string BuildTestIndex(const std::string& name, const std::string& extra = "") {
    std::string bsi = "/tmp/bamsix_errsweep_" + name + ".bsi";
    std::string cmd = "./build/bamsix build --input " + kTestBam +
                      " --output " + bsi + " " + extra + " 2>/dev/null";
    if (system(cmd.c_str()) != 0) return "";
    return bsi;
}

// ─── CLI-level Error Paths ──────────────────────────────────────────────────

TEST(ErrorSweep, EmptyPatternRejected) {
    auto bsi = BuildTestIndex("empty_pat");
    ASSERT_FALSE(bsi.empty());
    // Count with missing pattern should fail
    auto r = RunBamsix("count --index " + bsi);
    EXPECT_NE(r.exit_code, 0);
    EXPECT_NE(r.output.find("error"), std::string::npos);
}

TEST(ErrorSweep, MissingIndexArgument) {
    auto r = RunBamsix("count --pattern ACGT");
    EXPECT_NE(r.exit_code, 0);
}

TEST(ErrorSweep, NonexistentIndexFile) {
    auto r = RunBamsix("count --index /nonexistent/path.bsi --pattern ACGT");
    EXPECT_NE(r.exit_code, 0);
}

TEST(ErrorSweep, CorruptMagicBytes) {
    std::string corrupt = "/tmp/bamsix_errsweep_corrupt.bsi";
    {
        std::ofstream ofs(corrupt, std::ios::binary);
        ofs << "GARBAGE_DATA_NOT_A_BSI_FILE";
    }
    auto r = RunBamsix("count --index " + corrupt + " --pattern ACGT");
    EXPECT_NE(r.exit_code, 0);
    EXPECT_NE(r.output.find("error"), std::string::npos);
}

TEST(ErrorSweep, VerifyCorruptFile) {
    std::string corrupt = "/tmp/bamsix_errsweep_corrupt_verify.bsi";
    {
        std::ofstream ofs(corrupt, std::ios::binary);
        ofs << "NOT_BSI";
    }
    auto r = RunBamsix("verify --index " + corrupt);
    EXPECT_NE(r.exit_code, 0);
}

TEST(ErrorSweep, LossyReconstructWithoutAllowLossy) {
    auto bsi = BuildTestIndex("lossy_noflag", "--lossy-bins 8");
    ASSERT_FALSE(bsi.empty());
    auto r = RunBamsix("reconstruct --index " + bsi +
                       " --output /tmp/bamsix_errsweep_nope.bam");
    EXPECT_NE(r.exit_code, 0);
    EXPECT_NE(r.output.find("lossy"), std::string::npos);
}

TEST(ErrorSweep, LossyReconstructWithAllowLossy) {
    auto bsi = BuildTestIndex("lossy_flag", "--lossy-bins 8");
    ASSERT_FALSE(bsi.empty());
    auto r = RunBamsix("reconstruct --index " + bsi +
                       " --output /tmp/bamsix_errsweep_ok.bam --allow-lossy");
    EXPECT_EQ(r.exit_code, 0);
}

TEST(ErrorSweep, BuildMissingInput) {
    auto r = RunBamsix("build --output /tmp/bamsix_errsweep_noinput.bsi");
    EXPECT_NE(r.exit_code, 0);
    EXPECT_NE(r.output.find("error"), std::string::npos);
}

TEST(ErrorSweep, BuildInvalidBamPath) {
    auto r = RunBamsix("build --input /nonexistent.bam --output /tmp/bamsix_errsweep_badbam.bsi");
    EXPECT_NE(r.exit_code, 0);
}

TEST(ErrorSweep, InvalidSubcommand) {
    auto r = RunBamsix("foobar");
    EXPECT_NE(r.exit_code, 0);
}

TEST(ErrorSweep, RegionMissingPattern) {
    auto bsi = BuildTestIndex("region_nopat");
    ASSERT_FALSE(bsi.empty());
    auto r = RunBamsix("region-count --index " + bsi + " --region chr1:0-1000");
    EXPECT_NE(r.exit_code, 0);
}

TEST(ErrorSweep, ReconstructMissingOutput) {
    auto bsi = BuildTestIndex("recon_noout");
    ASSERT_FALSE(bsi.empty());
    auto r = RunBamsix("reconstruct --index " + bsi);
    // Should fail gracefully (no --output)
    EXPECT_NE(r.exit_code, 0);
}

// ─── C ABI Error Paths ─────────────────────────────────────────────────────

TEST(ErrorSweep, CAbiNullIndex) {
    uint64_t count = 0;
    EXPECT_EQ(bamsix_global_count(nullptr, nullptr, 0, &count),
              BAMSIX_STATUS_INVALID_ARGUMENT);
}

TEST(ErrorSweep, CAbiNullPattern) {
    auto bsi = BuildTestIndex("cabi_nullpat");
    ASSERT_FALSE(bsi.empty());
    bamsix_index_t* idx = nullptr;
    ASSERT_EQ(bamsix_open(bsi.c_str(), &idx), BAMSIX_STATUS_OK);

    uint64_t count = 0;
    EXPECT_EQ(bamsix_global_count(idx, nullptr, 0, &count),
              BAMSIX_STATUS_INVALID_ARGUMENT);

    bamsix_free(&idx);
}

TEST(ErrorSweep, CAbiZeroLengthPattern) {
    auto bsi = BuildTestIndex("cabi_zerolen");
    ASSERT_FALSE(bsi.empty());
    bamsix_index_t* idx = nullptr;
    ASSERT_EQ(bamsix_open(bsi.c_str(), &idx), BAMSIX_STATUS_OK);

    uint64_t count = 0;
    uint8_t pat[] = {0};
    EXPECT_EQ(bamsix_global_count(idx, pat, 0, &count),
              BAMSIX_STATUS_INVALID_ARGUMENT);

    bamsix_free(&idx);
}

TEST(ErrorSweep, CAbiNullOutputPtr) {
    auto bsi = BuildTestIndex("cabi_nullout");
    ASSERT_FALSE(bsi.empty());
    bamsix_index_t* idx = nullptr;
    ASSERT_EQ(bamsix_open(bsi.c_str(), &idx), BAMSIX_STATUS_OK);

    uint8_t pat[] = {0, 1, 2, 3};
    EXPECT_EQ(bamsix_global_count(idx, pat, 4, nullptr),
              BAMSIX_STATUS_INVALID_ARGUMENT);

    bamsix_free(&idx);
}

TEST(ErrorSweep, CAbiOpenNonexistent) {
    bamsix_index_t* idx = nullptr;
    EXPECT_NE(bamsix_open("/nonexistent.bsi", &idx), BAMSIX_STATUS_OK);
}

TEST(ErrorSweep, CAbiFreeNullSafe) {
    bamsix_index_t* idx = nullptr;
    bamsix_free(&idx);  // Should not crash
    bamsix_free(nullptr);  // Should not crash
    SUCCEED();
}

TEST(ErrorSweep, CAbiApproxStubs) {
    auto bsi = BuildTestIndex("cabi_approx");
    ASSERT_FALSE(bsi.empty());
    bamsix_index_t* idx = nullptr;
    ASSERT_EQ(bamsix_open(bsi.c_str(), &idx), BAMSIX_STATUS_OK);

    uint8_t pat[] = {0, 1, 2, 3};
    bamsix_locate_result_t results[10];
    size_t n = 0;

    EXPECT_EQ(bamsix_approx_locate_hamming(idx, pat, 4, 1, results, 10, &n),
              BAMSIX_STATUS_NOT_IMPLEMENTED_V1);
    EXPECT_EQ(bamsix_approx_locate_edit(idx, pat, 4, 1, results, 10, &n),
              BAMSIX_STATUS_NOT_IMPLEMENTED_V1);

    bamsix_free(&idx);
}

}  // namespace
}  // namespace bamsix
