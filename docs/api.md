# BAMSIX C API Reference

**Reference:** Contract §10.3, Architecture §12

## Overview

The BAMSIX C ABI (`libbamsi.h` / `bamsix.h`) provides a stable interface for FFI integration from any language. All functions are `extern "C"`, use only C-compatible types, and are thread-safe for concurrent reads on the same index handle.

## Header

```c
#include "bamsix/bamsix.h"
```

## Status Codes

```c
typedef enum bamsix_status_t {
    BAMSIX_STATUS_OK                    = 0,
    BAMSIX_STATUS_INVALID_ARGUMENT      = 1,
    BAMSIX_STATUS_NOT_IMPLEMENTED_V1    = 2,
    BAMSIX_STATUS_CORRUPT_BSI           = 3,
    BAMSIX_STATUS_FILE_NOT_FOUND        = 4,
    BAMSIX_STATUS_STREAM_DECODE_ERROR   = 5,
    BAMSIX_STATUS_SEPARATOR_POSITION    = 6,
    BAMSIX_STATUS_INTERNAL_ERROR        = 99
} bamsix_status_t;
```

## Lifecycle

### `bamsix_open`

```c
bamsix_status_t bamsix_open(const char* path, bamsix_index_t** out);
```

Open a `.bsi` file and return an opaque index handle. The caller must free the handle with `bamsix_free`.

### `bamsix_free`

```c
void bamsix_free(bamsix_index_t** idx);
```

Free an index handle. Sets `*idx` to NULL. Safe to call with `*idx == NULL`.

### `bamsix_verify`

```c
bamsix_status_t bamsix_verify(const char* path, int* valid_out);
```

Verify file integrity. Sets `*valid_out` to 1 (valid) or 0 (corrupt).

## Version

```c
const char* bamsix_version(void);
uint16_t bamsix_format_version(void);
```

## Index Information

```c
bamsix_status_t bamsix_get_n_reads(const bamsix_index_t* idx, uint64_t* out);
bamsix_status_t bamsix_get_s_length(const bamsix_index_t* idx, uint64_t* out);
bamsix_status_t bamsix_get_n_windows(const bamsix_index_t* idx, uint32_t* out);
bamsix_status_t bamsix_get_chrom_count(const bamsix_index_t* idx, uint32_t* out);
bamsix_status_t bamsix_get_chrom_name(const bamsix_index_t* idx, uint32_t chrom_idx,
                                     char* buf, size_t buf_len, size_t* out_len);
bamsix_status_t bamsix_is_lossless(const bamsix_index_t* idx, int* out);
```

## Queries

### GlobalCount

```c
bamsix_status_t bamsix_global_count(const bamsix_index_t* idx,
                                   const uint8_t* pattern, size_t pat_len,
                                   uint64_t* count);
```

Pattern encoding: A=0, C=1, G=2, T=3, N=4.

### GlobalExists

```c
bamsix_status_t bamsix_global_exists(const bamsix_index_t* idx,
                                    const uint8_t* pattern, size_t pat_len,
                                    uint64_t threshold, int* exists);
```

### Locate

```c
typedef struct bamsix_locate_result {
    uint32_t chrom_id;
    uint64_t p_min;
    uint64_t p_max;
    uint64_t read_id;
    int      is_reverse;
} bamsix_locate_result_t;

bamsix_status_t bamsix_locate(const bamsix_index_t* idx,
                             const uint8_t* pattern, size_t pat_len,
                             bamsix_locate_result_t* results,
                             size_t max_results, size_t* n_results);
```

### RegionalCount / RegionalExists

```c
bamsix_status_t bamsix_regional_count(const bamsix_index_t* idx,
                                     const uint8_t* pattern, size_t pat_len,
                                     const char* chrom,
                                     uint64_t start, uint64_t end,
                                     uint64_t* count);

bamsix_status_t bamsix_regional_exists(const bamsix_index_t* idx,
                                      const uint8_t* pattern, size_t pat_len,
                                      const char* chrom,
                                      uint64_t start, uint64_t end,
                                      uint64_t threshold, int* exists);
```

### Approximate Queries (V1 Stubs)

```c
// Both return BAMSIX_STATUS_NOT_IMPLEMENTED_V1 in v1.0
bamsix_status_t bamsix_approx_locate_hamming(...);
bamsix_status_t bamsix_approx_locate_edit(...);
```

## Usage Example

```c
#include "bamsix/bamsix.h"
#include <stdio.h>

int main(void) {
    bamsix_index_t* idx = NULL;
    if (bamsix_open("genome.bsi", &idx) != BAMSIX_STATUS_OK) {
        fprintf(stderr, "Failed to open index\n");
        return 1;
    }

    uint8_t pattern[] = {0, 1, 2, 3};  // ACGT
    uint64_t count;
    bamsix_global_count(idx, pattern, 4, &count);
    printf("ACGT count: %lu\n", (unsigned long)count);

    bamsix_free(&idx);
    return 0;
}
```

## Thread Safety

- All `bamsi_*` query functions are safe for concurrent use on the same index handle.
- `bamsix_open` and `bamsix_free` must not be called concurrently on the same handle.
