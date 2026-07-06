/**
 * C ABI Test Program — Stage 3 §3.3.6
 * Tests every exported C ABI function with valid and invalid inputs.
 * Compile with C compiler to verify pure-C ABI compatibility.
 */

#include "bamsix/bamsix.h"
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
    fprintf(stderr, "=== BAMSIX C ABI Test ===\n\n");

    /* Version */
    const char* ver = bamsix_version();
    CHECK(ver != NULL, "bamsix_version() non-null");
    if (ver != NULL) {
    	CHECK(strlen(ver) > 0, "bamsix_version() non-empty");
    	fprintf(stderr, "  Version: %s\n", ver);
    } 
    else {
    	fprintf(stderr, "  Version: (null)\n");
    }
    fprintf(stderr, "  Version: %s\n", ver);

    uint16_t fmt_ver = bamsix_format_version();
    CHECK(fmt_ver > 0, "bamsix_format_version() > 0");
    fprintf(stderr, "  Format version: %u\n", fmt_ver);

    /* Open with NULL → error */
    bamsix_index_t* idx = NULL;
    bamsix_status_t st = bamsix_open(NULL, &idx);
    CHECK(st != BAMSIX_STATUS_OK, "bamsix_open(NULL) returns error");
    CHECK(idx == NULL, "idx remains NULL on error");

    st = bamsix_open("data/test/synthetic_10reads.bsi", NULL);
    CHECK(st != BAMSIX_STATUS_OK, "bamsix_open(_, NULL) returns error");

    /* Open valid .bsi */
    st = bamsix_open("data/test/synthetic_10reads.bsi", &idx);
    CHECK(st == BAMSIX_STATUS_OK, "bamsix_open() succeeds");
    CHECK(idx != NULL, "idx is non-NULL after open");

    if (st != BAMSIX_STATUS_OK || !idx) {
        fprintf(stderr, "Cannot proceed without valid index. Aborting.\n");
        return 1;
    }

    /* Info queries */
    uint64_t n_reads = 0;
    st = bamsix_get_n_reads(idx, &n_reads);
    CHECK(st == BAMSIX_STATUS_OK, "bamsix_get_n_reads() succeeds");
    CHECK(n_reads == 10, "N_reads == 10");

    uint64_t s_len = 0;
    st = bamsix_get_s_length(idx, &s_len);
    CHECK(st == BAMSIX_STATUS_OK, "bamsix_get_s_length() succeeds");
    CHECK(s_len == 149, "S_length == 149");

    uint32_t n_windows = 0;
    st = bamsix_get_n_windows(idx, &n_windows);
    CHECK(st == BAMSIX_STATUS_OK, "bamsix_get_n_windows() succeeds");
    CHECK(n_windows == 1, "N_windows == 1");

    uint32_t chrom_count = 0;
    st = bamsix_get_chrom_count(idx, &chrom_count);
    CHECK(st == BAMSIX_STATUS_OK, "bamsix_get_chrom_count() succeeds");
    CHECK(chrom_count == 1, "chrom_count == 1");

    char chrom_buf[64];
    size_t chrom_len = 0;
    st = bamsix_get_chrom_name(idx, 0, chrom_buf, sizeof(chrom_buf), &chrom_len);
    CHECK(st == BAMSIX_STATUS_OK, "bamsix_get_chrom_name(0) succeeds");
    CHECK(strcmp(chrom_buf, "chr1") == 0, "chrom_name[0] == 'chr1'");

    st = bamsix_get_chrom_name(idx, 999, chrom_buf, sizeof(chrom_buf), &chrom_len);
    CHECK(st == BAMSIX_STATUS_INVALID_ARGUMENT, "bamsix_get_chrom_name(999) invalid");

    int lossless = 0;
    st = bamsix_is_lossless(idx, &lossless);
    CHECK(st == BAMSIX_STATUS_OK, "bamsix_is_lossless() succeeds");
    CHECK(lossless == 1, "is_lossless == 1");

    /* GlobalCount */
    uint8_t pat_acgt[] = {0, 1, 2, 3};  /* ACGT */
    uint64_t count = 0;
    st = bamsix_global_count(idx, pat_acgt, 4, &count);
    CHECK(st == BAMSIX_STATUS_OK, "bamsix_global_count(ACGT) succeeds");
    CHECK(count == 13, "GlobalCount(ACGT) == 13");

    /* GlobalCount with NULL → error */
    st = bamsix_global_count(idx, NULL, 4, &count);
    CHECK(st == BAMSIX_STATUS_INVALID_ARGUMENT, "bamsix_global_count(NULL pattern) invalid");

    st = bamsix_global_count(idx, pat_acgt, 0, &count);
    CHECK(st == BAMSIX_STATUS_INVALID_ARGUMENT, "bamsix_global_count(len=0) invalid");

    /* GlobalExists */
    int exists = 0;
    st = bamsix_global_exists(idx, pat_acgt, 4, 1, &exists);
    CHECK(st == BAMSIX_STATUS_OK, "bamsix_global_exists(ACGT) succeeds");
    CHECK(exists == 1, "GlobalExists(ACGT) == true");

    uint8_t pat_nonexistent[] = {4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4};
    st = bamsix_global_exists(idx, pat_nonexistent, 16, 1, &exists);
    CHECK(st == BAMSIX_STATUS_OK, "bamsix_global_exists(16xN) succeeds");
    CHECK(exists == 0, "GlobalExists(16xN) == false");

    /* Locate */
    bamsix_locate_result_t results[32];
    size_t n_results = 0;
    st = bamsix_locate(idx, pat_acgt, 4, results, 32, &n_results);
    CHECK(st == BAMSIX_STATUS_OK, "bamsix_locate(ACGT) succeeds");
    CHECK(n_results > 0, "Locate(ACGT) returns results");
    fprintf(stderr, "  Locate(ACGT): %zu results\n", n_results);

    /* RegionalCount */
    uint64_t rcount = 0;
    st = bamsix_regional_count(idx, pat_acgt, 4, "chr1", 100, 500, &rcount);
    CHECK(st == BAMSIX_STATUS_OK, "bamsix_regional_count(chr1:100-500) succeeds");
    fprintf(stderr, "  RegionalCount(ACGT, chr1:100-500): %lu\n", (unsigned long)rcount);

    /* RegionalExists */
    int rexists = 0;
    st = bamsix_regional_exists(idx, pat_acgt, 4, "chr1", 1, 1000000, 1, &rexists);
    CHECK(st == BAMSIX_STATUS_OK, "bamsix_regional_exists(chr1) succeeds");
    CHECK(rexists == 1, "RegionalExists(ACGT, chr1) == true");

    /* Approximate stubs → NOT_IMPLEMENTED_V1 */
    size_t approx_n = 0;
    st = bamsix_approx_locate_hamming(idx, pat_acgt, 4, 1, results, 32, &approx_n);
    CHECK(st == BAMSIX_STATUS_NOT_IMPLEMENTED_V1, "bamsix_approx_locate_hamming → NOT_IMPLEMENTED_V1");

    st = bamsix_approx_locate_edit(idx, pat_acgt, 4, 1, results, 32, &approx_n);
    CHECK(st == BAMSIX_STATUS_NOT_IMPLEMENTED_V1, "bamsix_approx_locate_edit → NOT_IMPLEMENTED_V1");

    /* Verify */
    int valid = 0;
    st = bamsix_verify("data/test/synthetic_10reads.bsi", &valid);
    CHECK(st == BAMSIX_STATUS_OK, "bamsix_verify() succeeds");
    CHECK(valid == 1, "bamsix_verify() returns valid");

    st = bamsix_verify(NULL, &valid);
    CHECK(st == BAMSIX_STATUS_INVALID_ARGUMENT, "bamsix_verify(NULL) invalid");

    /* Free */
    bamsix_free(&idx);
    CHECK(idx == NULL, "idx is NULL after bamsix_free");
    bamsix_free(&idx);  /* double-free safety */
    CHECK(idx == NULL, "double bamsix_free is safe");

    /* NULL handle queries → error */
    st = bamsix_get_n_reads(NULL, &n_reads);
    CHECK(st == BAMSIX_STATUS_INVALID_ARGUMENT, "bamsix_get_n_reads(NULL) invalid");

    st = bamsix_global_count(NULL, pat_acgt, 4, &count);
    CHECK(st == BAMSIX_STATUS_INVALID_ARGUMENT, "bamsix_global_count(NULL) invalid");

    fprintf(stderr, "\n=== C ABI Test: %d passed, %d failed ===\n", pass, fail);
    return fail;
}
