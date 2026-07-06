/**
 * @file bamsix.h
 * @brief BAMSIX Stable C ABI — Contract §10.3
 *
 * This is the ground-truth C ABI for BAMSIX. All functions are
 * extern "C" and use only C-compatible types. Thread-safe for
 * concurrent reads on the same index handle.
 */

#ifndef BAMSIX_BAMSIX_H
#define BAMSIX_BAMSIX_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ─── Opaque types ──────────────────────────────────────────────────── */

typedef struct bamsix_index bamsix_index_t;

/* ─── Status codes ──────────────────────────────────────────────────── */

typedef enum bamsix_status_t {
    BAMSIX_STATUS_OK                    = 0,
    BAMSIX_STATUS_INVALID_ARGUMENT      = 1,
    BAMSIX_STATUS_NOT_IMPLEMENTED_V1    = 2,
    BAMSIX_STATUS_CORRUPT_BSI           = 3,
    BAMSIX_STATUS_FILE_NOT_FOUND        = 4,
    BAMSIX_STATUS_STREAM_DECODE_ERROR   = 5,
    BAMSIX_STATUS_SEPARATOR_POSITION    = 6,
    /* H6 fix: status codes matching all C++ ErrorCode variants */
    BAMSIX_STATUS_CHECKSUM_MISMATCH     = 7,
    BAMSIX_STATUS_ORDERING_HASH_MISMATCH = 8,
    BAMSIX_STATUS_VERSION_MISMATCH      = 9,
    BAMSIX_STATUS_LOSSY_RECONSTRUCTION  = 10,
    BAMSIX_STATUS_UNSUPPORTED_CODEC     = 11,
    BAMSIX_STATUS_BUILD_FAILED          = 12,
    BAMSIX_STATUS_INTERNAL_ERROR        = 99
} bamsix_status_t;

/* ─── Version ───────────────────────────────────────────────────────── */

/**
 * Return the BAMSIX library version string (e.g. "1.0.0").
 * The returned pointer is statically allocated; do NOT free it.
 */
const char* bamsix_version(void);

/**
 * Return the .bsi format version number.
 */
uint16_t bamsix_format_version(void);

/* ─── Index lifecycle ───────────────────────────────────────────────── */

/**
 * Open a .bsi index file and return a handle.
 * @param path    Path to the .bsi file.
 * @param out     On success, receives a non-NULL handle.
 * @return BAMSIX_STATUS_OK on success.
 */
bamsix_status_t bamsix_open(const char* path, bamsix_index_t** out);

/**
 * Free an index handle. After this call, *idx is NULL.
 */
void bamsix_free(bamsix_index_t** idx);

/* ─── Verification ──────────────────────────────────────────────────── */

/**
 * Verify .bsi file integrity (global xxHash64 checksum).
 * @param path      Path to the .bsi file.
 * @param valid_out On success, set to 1 (valid) or 0 (corrupt).
 * @return BAMSIX_STATUS_OK on success.
 */
bamsix_status_t bamsix_verify(const char* path, int* valid_out);

/* ─── Index info ────────────────────────────────────────────────────── */

/**
 * Get the number of reads in the index.
 */
bamsix_status_t bamsix_get_n_reads(const bamsix_index_t* idx, uint64_t* out);

/**
 * Get |S| (concatenated sequence length).
 */
bamsix_status_t bamsix_get_s_length(const bamsix_index_t* idx, uint64_t* out);

/**
 * Get the number of windows.
 */
bamsix_status_t bamsix_get_n_windows(const bamsix_index_t* idx, uint32_t* out);

/**
 * Get the number of chromosomes.
 */
bamsix_status_t bamsix_get_chrom_count(const bamsix_index_t* idx, uint32_t* out);

/**
 * Get a chromosome name by index (0-based).
 * @param buf      Output buffer.
 * @param buf_len  Buffer size.
 * @param out_len  Actual length of the name (may exceed buf_len).
 */
bamsix_status_t bamsix_get_chrom_name(const bamsix_index_t* idx, uint32_t chrom_idx,
                                     char* buf, size_t buf_len, size_t* out_len);

/**
 * Check if the index was built in lossless mode.
 * @param out  Set to 1 (lossless) or 0 (lossy).
 */
bamsix_status_t bamsix_is_lossless(const bamsix_index_t* idx, int* out);

/* ─── Query: GlobalCount ────────────────────────────────────────────── */

/**
 * Count occurrences of a pattern globally.
 * @param idx      Index handle.
 * @param pattern  Pattern as base codes (A=0, C=1, G=2, T=3, N=4).
 * @param pat_len  Length of pattern.
 * @param count    Output count.
 * @return BAMSIX_STATUS_OK on success.
 */
bamsix_status_t bamsix_global_count(const bamsix_index_t* idx,
                                   const uint8_t* pattern, size_t pat_len,
                                   uint64_t* count);

/* ─── Query: GlobalExists ───────────────────────────────────────────── */

/**
 * Check if a pattern exists globally (at least 'threshold' occurrences).
 * @param idx       Index handle.
 * @param pattern   Pattern as base codes.
 * @param pat_len   Length of pattern.
 * @param threshold Minimum count (default 1 for simple existence).
 * @param exists    Output: 1 if exists, 0 otherwise.
 */
bamsix_status_t bamsix_global_exists(const bamsix_index_t* idx,
                                    const uint8_t* pattern, size_t pat_len,
                                    uint64_t threshold, int* exists);

/* ─── Query: Locate ─────────────────────────────────────────────────── */

/**
 * Result of a Locate query.
 */
typedef struct bamsix_locate_result {
    uint32_t chrom_id;
    uint64_t p_min;
    uint64_t p_max;
    uint64_t read_id;
    int      is_reverse;     /* 0 = forward, 1 = reverse complement */
} bamsix_locate_result_t;

/**
 * Locate all occurrences of a pattern.
 * @param idx         Index handle.
 * @param pattern     Pattern as base codes.
 * @param pat_len     Pattern length.
 * @param results     Caller-allocated output array.
 * @param max_results Maximum results to return.
 * @param n_results   Actual number of results returned.
 */
bamsix_status_t bamsix_locate(const bamsix_index_t* idx,
                             const uint8_t* pattern, size_t pat_len,
                             bamsix_locate_result_t* results,
                             size_t max_results, size_t* n_results);

/* ─── Query: RegionalCount ──────────────────────────────────────────── */

/**
 * Count occurrences of a pattern in a genomic region.
 * @param idx      Index handle.
 * @param pattern  Pattern as base codes.
 * @param pat_len  Pattern length.
 * @param chrom    Chromosome name (null-terminated string).
 * @param start    Region start (1-based inclusive).
 * @param end      Region end (1-based inclusive).
 * @param count    Output count.
 */
bamsix_status_t bamsix_regional_count(const bamsix_index_t* idx,
                                     const uint8_t* pattern, size_t pat_len,
                                     const char* chrom,
                                     uint64_t start, uint64_t end,
                                     uint64_t* count);

/* ─── Query: RegionalExists ─────────────────────────────────────────── */

bamsix_status_t bamsix_regional_exists(const bamsix_index_t* idx,
                                      const uint8_t* pattern, size_t pat_len,
                                      const char* chrom,
                                      uint64_t start, uint64_t end,
                                      uint64_t threshold, int* exists);

/* ─── Approximate Query Stubs (V1) ──────────────────────────────────── */

/**
 * Approximate locate with Hamming distance.
 * V1: returns BAMSIX_STATUS_NOT_IMPLEMENTED_V1.
 */
bamsix_status_t bamsix_approx_locate_hamming(const bamsix_index_t* idx,
                                            const uint8_t* pattern, size_t pat_len,
                                            uint32_t max_distance,
                                            bamsix_locate_result_t* results,
                                            size_t max_results, size_t* n_results);

/**
 * Approximate locate with edit distance.
 * V1: returns BAMSIX_STATUS_NOT_IMPLEMENTED_V1.
 */
bamsix_status_t bamsix_approx_locate_edit(const bamsix_index_t* idx,
                                         const uint8_t* pattern, size_t pat_len,
                                         uint32_t max_distance,
                                         bamsix_locate_result_t* results,
                                         size_t max_results, size_t* n_results);
/* ─── Streaming Locate Iterator (Contract §10.3) ───────────────────── */

typedef struct bamsix_locate_iter bamsix_locate_iter_t;

/**
 * Create a streaming locate iterator.
 * @param idx      Index handle.
 * @param pattern  Pattern as base codes.
 * @param pat_len  Pattern length.
 * @param iter     On success, receives a non-NULL iterator handle.
 * @return BAMSIX_STATUS_OK on success.
 */
bamsix_status_t bamsix_locate_iter_create(const bamsix_index_t* idx,
                                         const uint8_t* pattern, size_t pat_len,
                                         bamsix_locate_iter_t** iter);

/**
 * Advance the iterator and return the next result.
 * @param iter     Iterator handle.
 * @param result   Output: the next locate result.
 * @param has_more Output: 1 if more results remain, 0 if exhausted.
 * @return BAMSIX_STATUS_OK on success.
 */
bamsix_status_t bamsix_locate_iter_next(bamsix_locate_iter_t* iter,
                                       bamsix_locate_result_t* result,
                                       int* has_more);

/**
 * Free an iterator handle. After this call, *iter is NULL.
 */
void bamsix_locate_iter_free(bamsix_locate_iter_t** iter);

#ifdef __cplusplus
}
#endif

#endif  /* BAMSIX_BAMSIX_H */
