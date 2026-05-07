/// BAMSI Error Handling Sweep — Contract §0.8 compliance
/// Tests every ErrorCode path via C ABI + CLI subprocess.
/// Ensures "no silent failure" for all documented error conditions.

#include <gtest/gtest.h>
#include "bamsi/bamsi.h"

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>

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

std::string BuildTestIndex(const std::string& name, const std::string& extra = "") {
    std::string bsi = "/tmp/bamsi_errsweep_" + name + ".bsi";
    std::string cmd = "./build/bamsi build --input " + kTestBam +
                      " --output " + bsi + " " + extra + " 2>/dev/null";
    if (system(cmd.c_str()) != 0) return "";
    return bsi;
}

// ─── CLI-level Error Paths ──────────────────────────────────────────────────

TEST(ErrorSweep, EmptyPatternRejected) {
    auto bsi = BuildTestIndex("empty_pat");
    ASSERT_FALSE(bsi.empty());
    // Count with missing pattern should fail
    auto r = RunBamsi("count --index " + bsi);
    EXPECT_NE(r.exit_code, 0);
    EXPECT_NE(r.output.find("error"), std::string::npos);
}

TEST(ErrorSweep, MissingIndexArgument) {
    auto r = RunBamsi("count --pattern ACGT");
    EXPECT_NE(r.exit_code, 0);
}

TEST(ErrorSweep, NonexistentIndexFile) {
    auto r = RunBamsi("count --index /nonexistent/path.bsi --pattern ACGT");
    EXPECT_NE(r.exit_code, 0);
}

TEST(ErrorSweep, CorruptMagicBytes) {
    std::string corrupt = "/tmp/bamsi_errsweep_corrupt.bsi";
    {
        std::ofstream ofs(corrupt, std::ios::binary);
        ofs << "GARBAGE_DATA_NOT_A_BSI_FILE";
    }
    auto r = RunBamsi("count --index " + corrupt + " --pattern ACGT");
    EXPECT_NE(r.exit_code, 0);
    EXPECT_NE(r.output.find("error"), std::string::npos);
}

TEST(ErrorSweep, VerifyCorruptFile) {
    std::string corrupt = "/tmp/bamsi_errsweep_corrupt_verify.bsi";
    {
        std::ofstream ofs(corrupt, std::ios::binary);
        ofs << "NOT_BSI";
    }
    auto r = RunBamsi("verify --index " + corrupt);
    EXPECT_NE(r.exit_code, 0);
}

TEST(ErrorSweep, LossyReconstructWithoutAllowLossy) {
    auto bsi = BuildTestIndex("lossy_noflag", "--lossy-bins 8");
    ASSERT_FALSE(bsi.empty());
    auto r = RunBamsi("reconstruct --index " + bsi +
                       " --output /tmp/bamsi_errsweep_nope.bam");
    EXPECT_NE(r.exit_code, 0);
    EXPECT_NE(r.output.find("lossy"), std::string::npos);
}

TEST(ErrorSweep, LossyReconstructWithAllowLossy) {
    auto bsi = BuildTestIndex("lossy_flag", "--lossy-bins 8");
    ASSERT_FALSE(bsi.empty());
    auto r = RunBamsi("reconstruct --index " + bsi +
                       " --output /tmp/bamsi_errsweep_ok.bam --allow-lossy");
    EXPECT_EQ(r.exit_code, 0);
}

TEST(ErrorSweep, BuildMissingInput) {
    auto r = RunBamsi("build --output /tmp/bamsi_errsweep_noinput.bsi");
    EXPECT_NE(r.exit_code, 0);
    EXPECT_NE(r.output.find("error"), std::string::npos);
}

TEST(ErrorSweep, BuildInvalidBamPath) {
    auto r = RunBamsi("build --input /nonexistent.bam --output /tmp/bamsi_errsweep_badbam.bsi");
    EXPECT_NE(r.exit_code, 0);
}

TEST(ErrorSweep, InvalidSubcommand) {
    auto r = RunBamsi("foobar");
    EXPECT_NE(r.exit_code, 0);
}

TEST(ErrorSweep, RegionMissingPattern) {
    auto bsi = BuildTestIndex("region_nopat");
    ASSERT_FALSE(bsi.empty());
    auto r = RunBamsi("region-count --index " + bsi + " --region chr1:0-1000");
    EXPECT_NE(r.exit_code, 0);
}

TEST(ErrorSweep, ReconstructMissingOutput) {
    auto bsi = BuildTestIndex("recon_noout");
    ASSERT_FALSE(bsi.empty());
    auto r = RunBamsi("reconstruct --index " + bsi);
    // Should fail gracefully (no --output)
    EXPECT_NE(r.exit_code, 0);
}

// ─── C ABI Error Paths ─────────────────────────────────────────────────────

TEST(ErrorSweep, CAbiNullIndex) {
    uint64_t count = 0;
    EXPECT_EQ(bamsi_global_count(nullptr, nullptr, 0, &count),
              BAMSI_STATUS_INVALID_ARGUMENT);
}

TEST(ErrorSweep, CAbiNullPattern) {
    auto bsi = BuildTestIndex("cabi_nullpat");
    ASSERT_FALSE(bsi.empty());
    bamsi_index_t* idx = nullptr;
    ASSERT_EQ(bamsi_open(bsi.c_str(), &idx), BAMSI_STATUS_OK);

    uint64_t count = 0;
    EXPECT_EQ(bamsi_global_count(idx, nullptr, 0, &count),
              BAMSI_STATUS_INVALID_ARGUMENT);

    bamsi_free(&idx);
}

TEST(ErrorSweep, CAbiZeroLengthPattern) {
    auto bsi = BuildTestIndex("cabi_zerolen");
    ASSERT_FALSE(bsi.empty());
    bamsi_index_t* idx = nullptr;
    ASSERT_EQ(bamsi_open(bsi.c_str(), &idx), BAMSI_STATUS_OK);

    uint64_t count = 0;
    uint8_t pat[] = {0};
    EXPECT_EQ(bamsi_global_count(idx, pat, 0, &count),
              BAMSI_STATUS_INVALID_ARGUMENT);

    bamsi_free(&idx);
}

TEST(ErrorSweep, CAbiNullOutputPtr) {
    auto bsi = BuildTestIndex("cabi_nullout");
    ASSERT_FALSE(bsi.empty());
    bamsi_index_t* idx = nullptr;
    ASSERT_EQ(bamsi_open(bsi.c_str(), &idx), BAMSI_STATUS_OK);

    uint8_t pat[] = {0, 1, 2, 3};
    EXPECT_EQ(bamsi_global_count(idx, pat, 4, nullptr),
              BAMSI_STATUS_INVALID_ARGUMENT);

    bamsi_free(&idx);
}

TEST(ErrorSweep, CAbiOpenNonexistent) {
    bamsi_index_t* idx = nullptr;
    EXPECT_NE(bamsi_open("/nonexistent.bsi", &idx), BAMSI_STATUS_OK);
}

TEST(ErrorSweep, CAbiFreeNullSafe) {
    bamsi_index_t* idx = nullptr;
    bamsi_free(&idx);  // Should not crash
    bamsi_free(nullptr);  // Should not crash
    SUCCEED();
}

TEST(ErrorSweep, CAbiApproxStubs) {
    auto bsi = BuildTestIndex("cabi_approx");
    ASSERT_FALSE(bsi.empty());
    bamsi_index_t* idx = nullptr;
    ASSERT_EQ(bamsi_open(bsi.c_str(), &idx), BAMSI_STATUS_OK);

    uint8_t pat[] = {0, 1, 2, 3};
    bamsi_locate_result_t results[10];
    size_t n = 0;

    EXPECT_EQ(bamsi_approx_locate_hamming(idx, pat, 4, 1, results, 10, &n),
              BAMSI_STATUS_NOT_IMPLEMENTED_V1);
    EXPECT_EQ(bamsi_approx_locate_edit(idx, pat, 4, 1, results, 10, &n),
              BAMSI_STATUS_NOT_IMPLEMENTED_V1);

    bamsi_free(&idx);
}

}  // namespace
}  // namespace bamsi
