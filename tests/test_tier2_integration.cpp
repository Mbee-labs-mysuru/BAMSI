/// BAMSI TIER 2 Integration Tests — Contract §S4 / Execution Plan §4.3
/// Tests: Error handling sweep, codec correctness, lossy mode, CLI flags.
/// Build: linked against bamsi-core + gtest

#include <gtest/gtest.h>

#include "bamsi/types.hpp"
#include "bamsi/config.hpp"
#include "bamsi/bamsi.h"
#include "format/format.hpp"
#include "query/query.hpp"
#include "fmindex/fmindex.hpp"

#include <filesystem>
#include <fstream>
#include <cstdlib>

namespace bamsi {
namespace {

const std::string kTestBam = "data/test/synthetic_10reads.bam";
const std::string kTestBam10k = "data/test/synthetic_10k.bam";
const std::string kTmpDir = "/tmp/bamsi_tier2_";

std::string TmpFile(const std::string& name) {
    return kTmpDir + name;
}

// Helper: run bamsi CLI and return exit code + stdout
struct CliResult {
    int exit_code;
    std::string output;
};

CliResult RunBamsi(const std::string& args) {
    std::string cmd = "./build/bamsi " + args + " 2>/dev/null";
    FILE* fp = popen(cmd.c_str(), "r");
    std::string output;
    char buf[4096];
    while (fgets(buf, sizeof(buf), fp)) output += buf;
    int rc = pclose(fp);
    return {WEXITSTATUS(rc), output};
}

// ─── Error Handling Sweep (Contract §0.8) ────────────────────────────────────

class ErrorHandlingTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Build a test index
        bsi_path_ = TmpFile("errors.bsi");
        std::string cmd = "./build/bamsi build --input " + kTestBam +
                          " --output " + bsi_path_;
        ASSERT_EQ(system(cmd.c_str()), 0);
    }
    std::string bsi_path_;
};

TEST_F(ErrorHandlingTest, EmptyPatternCountReturnsError) {
    auto r = RunBamsi("count --index " + bsi_path_ + " --pattern \"\"");
    // Empty pattern should be rejected
    EXPECT_NE(r.exit_code, 0);
}

TEST_F(ErrorHandlingTest, MissingIndexReturnsError) {
    auto r = RunBamsi("count --pattern ACGT");
    EXPECT_NE(r.exit_code, 0);
}

TEST_F(ErrorHandlingTest, NonexistentIndexReturnsError) {
    auto r = RunBamsi("count --index /nonexistent.bsi --pattern ACGT");
    EXPECT_NE(r.exit_code, 0);
}

TEST_F(ErrorHandlingTest, CorruptBsiReturnsError) {
    // Create a corrupt file
    std::string corrupt = TmpFile("corrupt.bsi");
    {
        std::ofstream ofs(corrupt, std::ios::binary);
        ofs << "NOT_A_BSI_FILE";
    }
    auto r = RunBamsi("count --index " + corrupt + " --pattern ACGT");
    EXPECT_NE(r.exit_code, 0);
}

TEST_F(ErrorHandlingTest, VerifyCorruptFileFails) {
    std::string corrupt = TmpFile("corrupt2.bsi");
    {
        std::ofstream ofs(corrupt, std::ios::binary);
        ofs << "NOT_A_BSI";
    }
    auto r = RunBamsi("verify --index " + corrupt);
    EXPECT_NE(r.exit_code, 0);
}

TEST_F(ErrorHandlingTest, LossyReconstructWithoutFlagFails) {
    std::string lossy = TmpFile("lossy_error.bsi");
    std::string cmd = "./build/bamsi build --input " + kTestBam +
                      " --output " + lossy + " --lossy-bins 8";
    ASSERT_EQ(system(cmd.c_str()), 0);
    auto r = RunBamsi("reconstruct --index " + lossy + " --output /tmp/nope.bam");
    EXPECT_NE(r.exit_code, 0);
}

// ─── C ABI Error Codes ───────────────────────────────────────────────────────

TEST(CAbiTest, NullPointerReturnsError) {
    EXPECT_EQ(bamsi_open(nullptr, nullptr), BAMSI_STATUS_INVALID_ARGUMENT);
    EXPECT_EQ(bamsi_verify(nullptr, nullptr), BAMSI_STATUS_INVALID_ARGUMENT);
    EXPECT_EQ(bamsi_global_count(nullptr, nullptr, 0, nullptr),
              BAMSI_STATUS_INVALID_ARGUMENT);
}

TEST(CAbiTest, ApproxQueriesReturnNotImplemented) {
    bamsi_index_t* idx = nullptr;
    std::string bsi = TmpFile("cabi.bsi");
    std::string cmd = "./build/bamsi build --input " + kTestBam +
                      " --output " + bsi + " 2>/dev/null";
    ASSERT_EQ(system(cmd.c_str()), 0);

    ASSERT_EQ(bamsi_open(bsi.c_str(), &idx), BAMSI_STATUS_OK);
    ASSERT_NE(idx, nullptr);

    uint8_t pat[] = {0, 1, 2, 3};
    bamsi_locate_result_t results[10];
    size_t n = 0;

    EXPECT_EQ(bamsi_approx_locate_hamming(idx, pat, 4, 1, results, 10, &n),
              BAMSI_STATUS_NOT_IMPLEMENTED_V1);
    EXPECT_EQ(bamsi_approx_locate_edit(idx, pat, 4, 1, results, 10, &n),
              BAMSI_STATUS_NOT_IMPLEMENTED_V1);

    bamsi_free(&idx);
    EXPECT_EQ(idx, nullptr);
}

TEST(CAbiTest, VersionReturnsNonNull) {
    EXPECT_NE(bamsi_version(), nullptr);
    EXPECT_GT(bamsi_format_version(), 0);
}

TEST(CAbiTest, OpenAndQueryRoundTrip) {
    bamsi_index_t* idx = nullptr;
    std::string bsi = TmpFile("cabi_roundtrip.bsi");
    std::string cmd = "./build/bamsi build --input " + kTestBam +
                      " --output " + bsi + " 2>/dev/null";
    ASSERT_EQ(system(cmd.c_str()), 0);
    ASSERT_EQ(bamsi_open(bsi.c_str(), &idx), BAMSI_STATUS_OK);

    // GlobalCount
    uint8_t pat[] = {0, 1, 2, 3};  // ACGT
    uint64_t count = 0;
    EXPECT_EQ(bamsi_global_count(idx, pat, 4, &count), BAMSI_STATUS_OK);
    EXPECT_GT(count, 0u);

    // GlobalExists
    int exists = 0;
    EXPECT_EQ(bamsi_global_exists(idx, pat, 4, 1, &exists), BAMSI_STATUS_OK);
    EXPECT_EQ(exists, 1);

    // Locate
    bamsi_locate_result_t results[100];
    size_t n_results = 0;
    EXPECT_EQ(bamsi_locate(idx, pat, 4, results, 100, &n_results), BAMSI_STATUS_OK);
    EXPECT_EQ(n_results, count);

    // Info getters
    uint64_t n_reads = 0;
    EXPECT_EQ(bamsi_get_n_reads(idx, &n_reads), BAMSI_STATUS_OK);
    EXPECT_EQ(n_reads, 10u);

    uint64_t s_len = 0;
    EXPECT_EQ(bamsi_get_s_length(idx, &s_len), BAMSI_STATUS_OK);
    EXPECT_EQ(s_len, 149u);

    uint32_t n_windows = 0;
    EXPECT_EQ(bamsi_get_n_windows(idx, &n_windows), BAMSI_STATUS_OK);
    EXPECT_GE(n_windows, 1u);

    int lossless = 0;
    EXPECT_EQ(bamsi_is_lossless(idx, &lossless), BAMSI_STATUS_OK);
    EXPECT_EQ(lossless, 1);

    uint32_t chrom_count = 0;
    EXPECT_EQ(bamsi_get_chrom_count(idx, &chrom_count), BAMSI_STATUS_OK);
    EXPECT_GE(chrom_count, 1u);

    char chrom_buf[256];
    size_t chrom_len = 0;
    EXPECT_EQ(bamsi_get_chrom_name(idx, 0, chrom_buf, 256, &chrom_len), BAMSI_STATUS_OK);
    EXPECT_GT(chrom_len, 0u);

    // Regional query
    uint64_t reg_count = 0;
    EXPECT_EQ(bamsi_regional_count(idx, pat, 4, chrom_buf, 0, UINT64_MAX, &reg_count),
              BAMSI_STATUS_OK);

    int reg_exists = 0;
    EXPECT_EQ(bamsi_regional_exists(idx, pat, 4, chrom_buf, 0, UINT64_MAX, 1, &reg_exists),
              BAMSI_STATUS_OK);

    // Verify
    int valid = 0;
    EXPECT_EQ(bamsi_verify(bsi.c_str(), &valid), BAMSI_STATUS_OK);
    EXPECT_EQ(valid, 1);

    bamsi_free(&idx);
}

// ─── Codec Correctness ──────────────────────────────────────────────────────

TEST(CodecTest, RoundTripLossless) {
    // Build, verify, count — should produce identical results
    std::string bsi1 = TmpFile("codec1.bsi");
    std::string bsi2 = TmpFile("codec2.bsi");

    std::string cmd1 = "./build/bamsi build --input " + kTestBam +
                       " --output " + bsi1 + " 2>/dev/null";
    std::string cmd2 = "./build/bamsi build --input " + kTestBam +
                       " --output " + bsi2 + " 2>/dev/null";
    ASSERT_EQ(system(cmd1.c_str()), 0);
    ASSERT_EQ(system(cmd2.c_str()), 0);

    // Both should verify
    auto r1 = RunBamsi("verify --index " + bsi1);
    auto r2 = RunBamsi("verify --index " + bsi2);
    EXPECT_EQ(r1.exit_code, 0);
    EXPECT_EQ(r2.exit_code, 0);

    // Same count for any pattern
    auto c1 = RunBamsi("count --index " + bsi1 + " --pattern ACGT");
    auto c2 = RunBamsi("count --index " + bsi2 + " --pattern ACGT");
    EXPECT_EQ(c1.output, c2.output);
}

TEST(CodecTest, LossyBuildVerifyReconstruct) {
    std::string bsi = TmpFile("lossy_codec.bsi");
    std::string bam = TmpFile("lossy_recon.bam");

    std::string cmd = "./build/bamsi build --input " + kTestBam +
                      " --output " + bsi + " --lossy-bins 8 2>/dev/null";
    ASSERT_EQ(system(cmd.c_str()), 0);

    // Verify passes
    auto r = RunBamsi("verify --index " + bsi);
    EXPECT_EQ(r.exit_code, 0);

    // Reconstruct with --allow-lossy
    auto r2 = RunBamsi("reconstruct --index " + bsi +
                        " --output " + bam + " --allow-lossy");
    EXPECT_EQ(r2.exit_code, 0);
}

// ─── CLI Flag Tests ─────────────────────────────────────────────────────────

TEST(CliTest, InfoJsonHasAllFields) {
    std::string bsi = TmpFile("info_json.bsi");
    std::string cmd = "./build/bamsi build --input " + kTestBam +
                      " --output " + bsi + " 2>/dev/null";
    ASSERT_EQ(system(cmd.c_str()), 0);

    auto r = RunBamsi("info --index " + bsi + " --json");
    EXPECT_EQ(r.exit_code, 0);

    // Check key fields present
    EXPECT_NE(r.output.find("\"format_version\""), std::string::npos);
    EXPECT_NE(r.output.find("\"N_reads\""), std::string::npos);
    EXPECT_NE(r.output.find("\"source_manifest_hash\""), std::string::npos);
    EXPECT_NE(r.output.find("\"ordering_hash\""), std::string::npos);
    EXPECT_NE(r.output.find("\"host_os_id\""), std::string::npos);
    EXPECT_NE(r.output.find("\"build_timestamp_utc\""), std::string::npos);
    EXPECT_NE(r.output.find("\"chrom_name_table\""), std::string::npos);
}

TEST(CliTest, LocateSortOutput) {
    std::string bsi = TmpFile("sort_locate.bsi");
    std::string cmd = "./build/bamsi build --input " + kTestBam +
                      " --output " + bsi + " 2>/dev/null";
    ASSERT_EQ(system(cmd.c_str()), 0);

    auto r = RunBamsi("locate --index " + bsi +
                       " --pattern ACGT --sort-output");
    EXPECT_EQ(r.exit_code, 0);
    EXPECT_NE(r.output.find("strand"), std::string::npos);  // TSV header
}

TEST(CliTest, RegionFlagParsing) {
    std::string bsi = TmpFile("region_flag.bsi");
    std::string cmd = "./build/bamsi build --input " + kTestBam +
                      " --output " + bsi + " 2>/dev/null";
    ASSERT_EQ(system(cmd.c_str()), 0);

    auto r = RunBamsi("region-count --index " + bsi +
                       " --pattern ACGT --region chr1:0-999999");
    EXPECT_EQ(r.exit_code, 0);
    EXPECT_NE(r.output.find_first_of("0123456789"), std::string::npos);
}

TEST(CliTest, StrandFlagOverride) {
    std::string bsi = TmpFile("strand_flag.bsi");
    std::string cmd = "./build/bamsi build --input " + kTestBam +
                      " --output " + bsi + " 2>/dev/null";
    ASSERT_EQ(system(cmd.c_str()), 0);

    auto r_complete = RunBamsi("count --index " + bsi +
                                " --pattern ACG --strand complete");
    auto r_single = RunBamsi("count --index " + bsi +
                              " --pattern ACG --strand single");
    EXPECT_EQ(r_complete.exit_code, 0);
    EXPECT_EQ(r_single.exit_code, 0);

    // ACG is not palindromic, so complete should be >= single
    int count_complete = std::stoi(r_complete.output);
    int count_single = std::stoi(r_single.output);
    EXPECT_GE(count_complete, count_single);
}

TEST(CliTest, VerifyJson) {
    std::string bsi = TmpFile("verify_json.bsi");
    std::string cmd = "./build/bamsi build --input " + kTestBam +
                      " --output " + bsi + " 2>/dev/null";
    ASSERT_EQ(system(cmd.c_str()), 0);

    auto r = RunBamsi("verify --index " + bsi + " --json");
    EXPECT_EQ(r.exit_code, 0);
    EXPECT_NE(r.output.find("\"valid\""), std::string::npos);
    EXPECT_NE(r.output.find("true"), std::string::npos);
}

// ─── Provenance Tests (Contract §9.2) ───────────────────────────────────────

TEST(ProvenanceTest, BuildTimestampIsNonZero) {
    std::string bsi = TmpFile("provenance.bsi");
    std::string cmd = "./build/bamsi build --input " + kTestBam +
                      " --output " + bsi + " 2>/dev/null";
    ASSERT_EQ(system(cmd.c_str()), 0);

    auto idx = ReadBsi(bsi);
    EXPECT_GT(idx.header.build_timestamp_utc, 0u);
    EXPECT_GT(idx.header.host_os_id, 0u);
    EXPECT_GT(idx.header.cpu_arch_id, 0u);
}

TEST(ProvenanceTest, ManifestHashIsNonZero) {
    std::string bsi = TmpFile("manifest.bsi");
    std::string cmd = "./build/bamsi build --input " + kTestBam +
                      " --output " + bsi + " 2>/dev/null";
    ASSERT_EQ(system(cmd.c_str()), 0);

    auto idx = ReadBsi(bsi);
    // At least one byte in manifest hash should be non-zero
    bool all_zero = true;
    for (auto b : idx.header.source_manifest_hash)
        if (b != 0) { all_zero = false; break; }
    EXPECT_FALSE(all_zero);

    all_zero = true;
    for (auto b : idx.header.ordering_hash)
        if (b != 0) { all_zero = false; break; }
    EXPECT_FALSE(all_zero);
}

}  // namespace
}  // namespace bamsi
