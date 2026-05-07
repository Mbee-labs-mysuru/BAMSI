/**
 * C ABI Test Program — Stage 3 §3.3.6
 * Tests every exported C ABI function with valid and invalid inputs.
 * Compile with C compiler to verify pure-C ABI compatibility.
 */

#include "bamsi/bamsi.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

static int pass = 0, fail = 0;

#define CHECK(cond, msg) do { \
    if (cond) { ++pass; } else { \
        fprintf(stderr, "FAIL: %s\n", msg); \
        ++fail; \
    } \
} while(0)

int main(void) {
    fprintf(stderr, "=== BAMSI C ABI Test ===\n\n");

    /* Version */
    const char* ver = bamsi_version();
    CHECK(ver != NULL, "bamsi_version() non-null");
    CHECK(strlen(ver) > 0, "bamsi_version() non-empty");
    fprintf(stderr, "  Version: %s\n", ver);

    uint16_t fmt_ver = bamsi_format_version();
    CHECK(fmt_ver > 0, "bamsi_format_version() > 0");
    fprintf(stderr, "  Format version: %u\n", fmt_ver);

    /* Open with NULL → error */
    bamsi_index_t* idx = NULL;
    bamsi_status_t st = bamsi_open(NULL, &idx);
    CHECK(st != BAMSI_STATUS_OK, "bamsi_open(NULL) returns error");
    CHECK(idx == NULL, "idx remains NULL on error");

    st = bamsi_open("data/test/synthetic_10reads.bsi", NULL);
    CHECK(st != BAMSI_STATUS_OK, "bamsi_open(_, NULL) returns error");

    /* Open valid .bsi */
    st = bamsi_open("data/test/synthetic_10reads.bsi", &idx);
    CHECK(st == BAMSI_STATUS_OK, "bamsi_open() succeeds");
    CHECK(idx != NULL, "idx is non-NULL after open");

    if (st != BAMSI_STATUS_OK || !idx) {
        fprintf(stderr, "Cannot proceed without valid index. Aborting.\n");
        return 1;
    }

    /* Info queries */
    uint64_t n_reads = 0;
    st = bamsi_get_n_reads(idx, &n_reads);
    CHECK(st == BAMSI_STATUS_OK, "bamsi_get_n_reads() succeeds");
    CHECK(n_reads == 10, "N_reads == 10");

    uint64_t s_len = 0;
    st = bamsi_get_s_length(idx, &s_len);
    CHECK(st == BAMSI_STATUS_OK, "bamsi_get_s_length() succeeds");
    CHECK(s_len == 149, "S_length == 149");

    uint32_t n_windows = 0;
    st = bamsi_get_n_windows(idx, &n_windows);
    CHECK(st == BAMSI_STATUS_OK, "bamsi_get_n_windows() succeeds");
    CHECK(n_windows == 1, "N_windows == 1");

    uint32_t chrom_count = 0;
    st = bamsi_get_chrom_count(idx, &chrom_count);
    CHECK(st == BAMSI_STATUS_OK, "bamsi_get_chrom_count() succeeds");
    CHECK(chrom_count == 1, "chrom_count == 1");

    char chrom_buf[64];
    size_t chrom_len = 0;
    st = bamsi_get_chrom_name(idx, 0, chrom_buf, sizeof(chrom_buf), &chrom_len);
    CHECK(st == BAMSI_STATUS_OK, "bamsi_get_chrom_name(0) succeeds");
    CHECK(strcmp(chrom_buf, "chr1") == 0, "chrom_name[0] == 'chr1'");

    st = bamsi_get_chrom_name(idx, 999, chrom_buf, sizeof(chrom_buf), &chrom_len);
    CHECK(st == BAMSI_STATUS_INVALID_ARGUMENT, "bamsi_get_chrom_name(999) invalid");

    int lossless = 0;
    st = bamsi_is_lossless(idx, &lossless);
    CHECK(st == BAMSI_STATUS_OK, "bamsi_is_lossless() succeeds");
    CHECK(lossless == 1, "is_lossless == 1");

    /* GlobalCount */
    uint8_t pat_acgt[] = {0, 1, 2, 3};  /* ACGT */
    uint64_t count = 0;
    st = bamsi_global_count(idx, pat_acgt, 4, &count);
    CHECK(st == BAMSI_STATUS_OK, "bamsi_global_count(ACGT) succeeds");
    CHECK(count == 13, "GlobalCount(ACGT) == 13");

    /* GlobalCount with NULL → error */
    st = bamsi_global_count(idx, NULL, 4, &count);
    CHECK(st == BAMSI_STATUS_INVALID_ARGUMENT, "bamsi_global_count(NULL pattern) invalid");

    st = bamsi_global_count(idx, pat_acgt, 0, &count);
    CHECK(st == BAMSI_STATUS_INVALID_ARGUMENT, "bamsi_global_count(len=0) invalid");

    /* GlobalExists */
    int exists = 0;
    st = bamsi_global_exists(idx, pat_acgt, 4, 1, &exists);
    CHECK(st == BAMSI_STATUS_OK, "bamsi_global_exists(ACGT) succeeds");
    CHECK(exists == 1, "GlobalExists(ACGT) == true");

    uint8_t pat_nonexistent[] = {4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4};
    st = bamsi_global_exists(idx, pat_nonexistent, 16, 1, &exists);
    CHECK(st == BAMSI_STATUS_OK, "bamsi_global_exists(16xN) succeeds");
    CHECK(exists == 0, "GlobalExists(16xN) == false");

    /* Locate */
    bamsi_locate_result_t results[32];
    size_t n_results = 0;
    st = bamsi_locate(idx, pat_acgt, 4, results, 32, &n_results);
    CHECK(st == BAMSI_STATUS_OK, "bamsi_locate(ACGT) succeeds");
    CHECK(n_results > 0, "Locate(ACGT) returns results");
    fprintf(stderr, "  Locate(ACGT): %zu results\n", n_results);

    /* RegionalCount */
    uint64_t rcount = 0;
    st = bamsi_regional_count(idx, pat_acgt, 4, "chr1", 100, 500, &rcount);
    CHECK(st == BAMSI_STATUS_OK, "bamsi_regional_count(chr1:100-500) succeeds");
    fprintf(stderr, "  RegionalCount(ACGT, chr1:100-500): %lu\n", (unsigned long)rcount);

    /* RegionalExists */
    int rexists = 0;
    st = bamsi_regional_exists(idx, pat_acgt, 4, "chr1", 1, 1000000, 1, &rexists);
    CHECK(st == BAMSI_STATUS_OK, "bamsi_regional_exists(chr1) succeeds");
    CHECK(rexists == 1, "RegionalExists(ACGT, chr1) == true");

    /* Approximate stubs → NOT_IMPLEMENTED_V1 */
    size_t approx_n = 0;
    st = bamsi_approx_locate_hamming(idx, pat_acgt, 4, 1, results, 32, &approx_n);
    CHECK(st == BAMSI_STATUS_NOT_IMPLEMENTED_V1, "bamsi_approx_locate_hamming → NOT_IMPLEMENTED_V1");

    st = bamsi_approx_locate_edit(idx, pat_acgt, 4, 1, results, 32, &approx_n);
    CHECK(st == BAMSI_STATUS_NOT_IMPLEMENTED_V1, "bamsi_approx_locate_edit → NOT_IMPLEMENTED_V1");

    /* Verify */
    int valid = 0;
    st = bamsi_verify("data/test/synthetic_10reads.bsi", &valid);
    CHECK(st == BAMSI_STATUS_OK, "bamsi_verify() succeeds");
    CHECK(valid == 1, "bamsi_verify() returns valid");

    st = bamsi_verify(NULL, &valid);
    CHECK(st == BAMSI_STATUS_INVALID_ARGUMENT, "bamsi_verify(NULL) invalid");

    /* Free */
    bamsi_free(&idx);
    CHECK(idx == NULL, "idx is NULL after bamsi_free");
    bamsi_free(&idx);  /* double-free safety */
    CHECK(idx == NULL, "double bamsi_free is safe");

    /* NULL handle queries → error */
    st = bamsi_get_n_reads(NULL, &n_reads);
    CHECK(st == BAMSI_STATUS_INVALID_ARGUMENT, "bamsi_get_n_reads(NULL) invalid");

    st = bamsi_global_count(NULL, pat_acgt, 4, &count);
    CHECK(st == BAMSI_STATUS_INVALID_ARGUMENT, "bamsi_global_count(NULL) invalid");

    fprintf(stderr, "\n=== C ABI Test: %d passed, %d failed ===\n", pass, fail);
    return fail;
}
