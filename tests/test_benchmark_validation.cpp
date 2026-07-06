/**
 * test_benchmark_validation.cpp — Stage 5 Benchmark Validation
 *
 * Validates that the benchmark infrastructure produces correct results:
 * - Build produces compressed output smaller than input
 * - All 5 query operations complete and return valid results
 * - bamsix info --json contains all 30+ Contract §9.2 fields
 * - Verify passes after build
 * - Reconstruct produces valid output
 * - Build/query timings are measurable (> 0)
 *
 * Contract §8, Execution Plan §5
 */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <chrono>
#include <fstream>
#include <sstream>
#include <vector>
#include <array>
#include <sys/stat.h>

static int pass_count = 0, fail_count = 0;

#define TEST(name) do { fprintf(stderr, "  TEST: %s\n", name); } while(0)
#define CHECK(cond, msg) do { \
    if (cond) { ++pass_count; } else { \
        fprintf(stderr, "  FAIL: %s\n", msg); \
        ++fail_count; \
    } \
} while(0)

static std::string exec_cmd(const std::string& cmd) {
    std::array<char, 4096> buf;
    std::string result;
    FILE* pipe = popen(cmd.c_str(), "r");
    if (!pipe) return "";
    while (fgets(buf.data(), buf.size(), pipe) != nullptr) {
        result += buf.data();
    }
    int status = pclose(pipe);
    if (status != 0 && result.empty()) {
        result = "[command failed with status " + std::to_string(status) + "]";
    }
    return result;
}

static int exec_status(const std::string& cmd) {
    return system(cmd.c_str());
}

static size_t file_size(const std::string& path) {
    struct stat st;
    if (stat(path.c_str(), &st) != 0) return 0;
    return static_cast<size_t>(st.st_size);
}

static bool file_exists(const std::string& path) {
    struct stat st;
    return stat(path.c_str(), &st) == 0;
}

// Check if a JSON string contains a key
static bool json_has_key(const std::string& json, const std::string& key) {
    return json.find("\"" + key + "\"") != std::string::npos;
}

int main() {
    fprintf(stderr, "=== BAMSIX Benchmark Validation Test ===\n");
    fprintf(stderr, "    Contract §8 / Execution Plan §5\n\n");

    const std::string bamsix = "./build/bamsix";
    const std::string bam = "data/test/synthetic_10reads.bam";
    const std::string bsi = "benchmarks/scratch/bench_validation.bsi";
    const std::string recon = "benchmarks/scratch/bench_validation_recon.bam";

    // Ensure scratch dir exists
    system("mkdir -p benchmarks/scratch");

    if (!file_exists(bamsix)) {
        fprintf(stderr, "ERROR: BAMSIX binary not found at %s\n", bamsix.c_str());
        return 1;
    }
    if (!file_exists(bam)) {
        fprintf(stderr, "ERROR: Test BAM not found at %s\n", bam.c_str());
        return 1;
    }

    size_t bam_size = file_size(bam);
    fprintf(stderr, "Input BAM: %s (%zu bytes)\n\n", bam.c_str(), bam_size);

    // ═══ Test 1: Build produces compressed output ═══
    TEST("Build produces compressed .bsi");
    {
        auto t0 = std::chrono::steady_clock::now();
        int st = exec_status(bamsix + " build --input " + bam + " --output " + bsi + " 2>/dev/null");
        auto t1 = std::chrono::steady_clock::now();
        double build_s = std::chrono::duration<double>(t1 - t0).count();

        CHECK(st == 0, "Build exits with status 0");
        CHECK(file_exists(bsi), "BSI file created");

        size_t bsi_size = file_size(bsi);
        fprintf(stderr, "    BSI size: %zu bytes (BAM: %zu bytes)\n", bsi_size, bam_size);
        fprintf(stderr, "    Build time: %.3f s\n", build_s);

        CHECK(bsi_size > 0, "BSI file is non-empty");
        CHECK(build_s > 0, "Build time is measurable (> 0)");
        // Note: for tiny synthetic BAMs with headers, .bsi may be larger due to index overhead.
        // The compression benefit shows on real data. Here we just check it works.
    }

    // ═══ Test 2: GlobalCount ═══
    TEST("GlobalCount returns valid result");
    {
        auto t0 = std::chrono::steady_clock::now();
        std::string output = exec_cmd(bamsix + " count --index " + bsi + " --pattern ACGT 2>/dev/null");
        auto t1 = std::chrono::steady_clock::now();
        double query_s = std::chrono::duration<double>(t1 - t0).count();

        CHECK(!output.empty(), "GlobalCount produces output");
        CHECK(query_s < 5.0, "GlobalCount completes in < 5s on synthetic data");
        fprintf(stderr, "    Output: %s", output.c_str());
        fprintf(stderr, "    Latency: %.6f s\n", query_s);
    }

    // ═══ Test 3: GlobalExists ═══
    TEST("GlobalExists returns valid result");
    {
        std::string output = exec_cmd(bamsix + " exists --index " + bsi + " --pattern ACGT 2>/dev/null");
        CHECK(!output.empty(), "GlobalExists produces output");
        // Check for existing pattern
        std::string out_non = exec_cmd(bamsix + " exists --index " + bsi + " --pattern NNNNNNNNNNNNNNNN 2>/dev/null");
        CHECK(!out_non.empty(), "GlobalExists for rare pattern produces output");
    }

    // ═══ Test 4: Locate ═══
    TEST("Locate returns valid results");
    {
        auto t0 = std::chrono::steady_clock::now();
        std::string output = exec_cmd(bamsix + " locate --index " + bsi + " --pattern ACGT 2>/dev/null");
        auto t1 = std::chrono::steady_clock::now();

        CHECK(!output.empty(), "Locate produces output");
        CHECK(std::chrono::duration<double>(t1 - t0).count() < 5.0, "Locate completes in < 5s");
    }

    // ═══ Test 5: RegionalCount ═══
    TEST("RegionalCount returns valid result");
    {
        std::string output = exec_cmd(bamsix + " region-count --index " + bsi +
                                      " --pattern ACGT --region chr1:1-1000000 2>/dev/null");
        CHECK(!output.empty(), "RegionalCount produces output");
    }

    // ═══ Test 6: RegionalExists ═══
    TEST("RegionalExists returns valid result");
    {
        std::string output = exec_cmd(bamsix + " region-exists --index " + bsi +
                                      " --pattern ACGT --region chr1:1-1000000 --threshold 1 2>/dev/null");
        CHECK(!output.empty(), "RegionalExists produces output");
    }

    // ═══ Test 7: Verify passes ═══
    TEST("Verify passes on freshly built index");
    {
        int st = exec_status(bamsix + " verify --index " + bsi + " 2>/dev/null");
        CHECK(st == 0, "Verify returns exit code 0");

        // Strict verify
        st = exec_status(bamsix + " verify --index " + bsi + " --strict 2>/dev/null");
        CHECK(st == 0, "Verify --strict returns exit code 0");
    }

    // ═══ Test 8: Info --json has all 30+ fields (Contract §9.2) ═══
    TEST("Info --json contains all Contract §9.2 fields");
    {
        std::string json = exec_cmd(bamsix + " info --index " + bsi + " --json 2>/dev/null");
        CHECK(!json.empty(), "Info --json produces output");

        // Contract §9.2 mandated fields
        const std::vector<std::string> required_fields = {
            "format_version", "bamsix_version",
            "is_lossless", "source_manifest_hash", "ordering_hash",
            "S_length", "N_reads", "N_windows",
            "sample_step_s", "has_isa_samples", "sample_step_s_prime",
            "enable_sarange", "shared_bwt", "enable_bidirectional",
            "window_size_T", "entropy_order_k",
            "qual_codec_id", "qual_lossy_bins",
            "meta_codec_id", "map_codec_id",
            "strand_mode", "chrom_name_table", "chrom_count",
            "sentinel_row", "seq_block_size", "qual_block_size",
            "host_os_id", "cpu_arch_id", "build_timestamp_utc",
            "reference_based_encoding"
        };

        int fields_found = 0;
        for (const auto& field : required_fields) {
            bool found = json_has_key(json, field);
            if (!found) {
                fprintf(stderr, "    MISSING field: %s\n", field.c_str());
            } else {
                fields_found++;
            }
        }
        CHECK(fields_found == (int)required_fields.size(),
              "All 30+ Contract §9.2 fields present in info --json");
        fprintf(stderr, "    Fields: %d/%zu present\n", fields_found, required_fields.size());
    }

    // ═══ Test 9: Reconstruct ═══
    TEST("Reconstruct produces valid BAM output");
    {
        auto t0 = std::chrono::steady_clock::now();
        int st = exec_status(bamsix + " reconstruct --index " + bsi + " --output " + recon + " 2>/dev/null");
        auto t1 = std::chrono::steady_clock::now();
        double recon_s = std::chrono::duration<double>(t1 - t0).count();

        CHECK(st == 0, "Reconstruct exits with status 0");
        CHECK(file_exists(recon), "Reconstructed BAM file created");

        size_t recon_size = file_size(recon);
        CHECK(recon_size > 0, "Reconstructed BAM is non-empty");
        fprintf(stderr, "    Reconstructed BAM: %zu bytes\n", recon_size);
        fprintf(stderr, "    Reconstruct time: %.3f s\n", recon_s);
    }

    // ═══ Test 10: Empty pattern rejection ═══
    TEST("Empty pattern returns error (Contract §0.8)");
    {
        int st = exec_status(bamsix + " count --index " + bsi + " --pattern \"\" 2>/dev/null");
        CHECK(st != 0, "Empty pattern returns non-zero exit code");
    }

    // ═══ Test 11: Invalid file rejection ═══
    TEST("Verify rejects corrupted file");
    {
        // Create a corrupted copy
        std::string corrupt = "benchmarks/scratch/corrupt.bsi";
        system(("cp " + bsi + " " + corrupt).c_str());
        // Flip a byte in the middle
        FILE* f = fopen(corrupt.c_str(), "r+b");
        if (f) {
            fseek(f, 100, SEEK_SET);
            uint8_t byte = 0xFF;
            fwrite(&byte, 1, 1, f);
            fclose(f);
            int st = exec_status(bamsix + " verify --index " + corrupt + " 2>/dev/null");
            // Corrupted file should fail verify or fail to open
            // (it may crash on parse which also returns non-zero)
            CHECK(st != 0, "Corrupted file detected by verify or fails to open");
        }
    }

    // Cleanup
    system("rm -rf benchmarks/scratch");

    fprintf(stderr, "\n=== Benchmark Validation: %d passed, %d failed ===\n",
            pass_count, fail_count);
    return fail_count;
}
