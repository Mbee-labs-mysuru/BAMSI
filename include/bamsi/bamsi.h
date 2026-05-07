/**
 * @file bamsi.h
 * @brief BAMSI Stable C ABI — Contract §10.3
 *
 * This is the ground-truth C ABI for BAMSI. All functions are
 * extern "C" and use only C-compatible types. Thread-safe for
 * concurrent reads on the same index handle.
 */

#ifndef BAMSI_BAMSI_H
#define BAMSI_BAMSI_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ─── Opaque types ──────────────────────────────────────────────────── */

typedef struct bamsi_index bamsi_index_t;

/* ─── Status codes ──────────────────────────────────────────────────── */

typedef enum bamsi_status_t {
    BAMSI_STATUS_OK                    = 0,
    BAMSI_STATUS_INVALID_ARGUMENT      = 1,
    BAMSI_STATUS_NOT_IMPLEMENTED_V1    = 2,
    BAMSI_STATUS_CORRUPT_BSI           = 3,
    BAMSI_STATUS_FILE_NOT_FOUND        = 4,
    BAMSI_STATUS_STREAM_DECODE_ERROR   = 5,
    BAMSI_STATUS_SEPARATOR_POSITION    = 6,
    BAMSI_STATUS_INTERNAL_ERROR        = 99
} bamsi_status_t;

/* ─── Version ───────────────────────────────────────────────────────── */

/**
 * Return the BAMSI library version string (e.g. "0.1.0").
 * The returned pointer is statically allocated; do NOT free it.
 */
const char* bamsi_version(void);

/**
 * Return the .bsi format version number.
 */
uint16_t bamsi_format_version(void);

/* ─── Index lifecycle ───────────────────────────────────────────────── */

/**
 * Open a .bsi index file and return a handle.
 * @param path    Path to the .bsi file.
 * @param out     On success, receives a non-NULL handle.
 * @return BAMSI_STATUS_OK on success.
 */
bamsi_status_t bamsi_open(const char* path, bamsi_index_t** out);

/**
 * Free an index handle. After this call, *idx is NULL.
 */
void bamsi_free(bamsi_index_t** idx);

/* ─── Verification ──────────────────────────────────────────────────── */

/**
 * Verify .bsi file integrity (global xxHash64 checksum).
 * @param path      Path to the .bsi file.
 * @param valid_out On success, set to 1 (valid) or 0 (corrupt).
 * @return BAMSI_STATUS_OK on success.
 */
bamsi_status_t bamsi_verify(const char* path, int* valid_out);

/* ─── Index info ────────────────────────────────────────────────────── */

/**
 * Get the number of reads in the index.
 */
bamsi_status_t bamsi_get_n_reads(const bamsi_index_t* idx, uint64_t* out);

/**
 * Get |S| (concatenated sequence length).
 */
bamsi_status_t bamsi_get_s_length(const bamsi_index_t* idx, uint64_t* out);

/**
 * Get the number of windows.
 */
bamsi_status_t bamsi_get_n_windows(const bamsi_index_t* idx, uint32_t* out);

/**
 * Get the number of chromosomes.
 */
bamsi_status_t bamsi_get_chrom_count(const bamsi_index_t* idx, uint32_t* out);

/**
 * Get a chromosome name by index (0-based).
 * @param buf      Output buffer.
 * @param buf_len  Buffer size.
 * @param out_len  Actual length of the name (may exceed buf_len).
 */
bamsi_status_t bamsi_get_chrom_name(const bamsi_index_t* idx, uint32_t chrom_idx,
                                     char* buf, size_t buf_len, size_t* out_len);

/**
 * Check if the index was built in lossless mode.
 * @param out  Set to 1 (lossless) or 0 (lossy).
 */
bamsi_status_t bamsi_is_lossless(const bamsi_index_t* idx, int* out);

/* ─── Query: GlobalCount ────────────────────────────────────────────── */

/**
 * Count occurrences of a pattern globally.
 * @param idx      Index handle.
 * @param pattern  Pattern as base codes (A=0, C=1, G=2, T=3, N=4).
 * @param pat_len  Length of pattern.
 * @param count    Output count.
 * @return BAMSI_STATUS_OK on success.
 */
bamsi_status_t bamsi_global_count(const bamsi_index_t* idx,
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
bamsi_status_t bamsi_global_exists(const bamsi_index_t* idx,
                                    const uint8_t* pattern, size_t pat_len,
                                    uint64_t threshold, int* exists);

/* ─── Query: Locate ─────────────────────────────────────────────────── */

/**
 * Result of a Locate query.
 */
typedef struct bamsi_locate_result {
    uint32_t chrom_id;
    uint64_t p_min;
    uint64_t p_max;
    uint64_t read_id;
    int      is_reverse;     /* 0 = forward, 1 = reverse complement */
} bamsi_locate_result_t;

/**
 * Locate all occurrences of a pattern.
 * @param idx         Index handle.
 * @param pattern     Pattern as base codes.
 * @param pat_len     Pattern length.
 * @param results     Caller-allocated output array.
 * @param max_results Maximum results to return.
 * @param n_results   Actual number of results returned.
 */
bamsi_status_t bamsi_locate(const bamsi_index_t* idx,
                             const uint8_t* pattern, size_t pat_len,
                             bamsi_locate_result_t* results,
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
bamsi_status_t bamsi_regional_count(const bamsi_index_t* idx,
                                     const uint8_t* pattern, size_t pat_len,
                                     const char* chrom,
                                     uint64_t start, uint64_t end,
                                     uint64_t* count);

/* ─── Query: RegionalExists ─────────────────────────────────────────── */

bamsi_status_t bamsi_regional_exists(const bamsi_index_t* idx,
                                      const uint8_t* pattern, size_t pat_len,
                                      const char* chrom,
                                      uint64_t start, uint64_t end,
                                      uint64_t threshold, int* exists);

/* ─── Approximate Query Stubs (V1) ──────────────────────────────────── */

/**
 * Approximate locate with Hamming distance.
 * V1: returns BAMSI_STATUS_NOT_IMPLEMENTED_V1.
 */
bamsi_status_t bamsi_approx_locate_hamming(const bamsi_index_t* idx,
                                            const uint8_t* pattern, size_t pat_len,
                                            uint32_t max_distance,
                                            bamsi_locate_result_t* results,
                                            size_t max_results, size_t* n_results);

/**
 * Approximate locate with edit distance.
 * V1: returns BAMSI_STATUS_NOT_IMPLEMENTED_V1.
 */
bamsi_status_t bamsi_approx_locate_edit(const bamsi_index_t* idx,
                                         const uint8_t* pattern, size_t pat_len,
                                         uint32_t max_distance,
                                         bamsi_locate_result_t* results,
                                         size_t max_results, size_t* n_results);

#ifdef __cplusplus
}
#endif

#endif  /* BAMSI_BAMSI_H */
