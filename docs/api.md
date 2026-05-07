# BAMSI C API Reference

**Reference:** Contract §10.3, Architecture §12

## Overview

The BAMSI C ABI (`libbamsi.h` / `bamsi.h`) provides a stable interface for FFI integration from any language. All functions are `extern "C"`, use only C-compatible types, and are thread-safe for concurrent reads on the same index handle.

## Header

```c
#include "bamsi/bamsi.h"
```

## Status Codes

```c
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
```

## Lifecycle

### `bamsi_open`

```c
bamsi_status_t bamsi_open(const char* path, bamsi_index_t** out);
```

Open a `.bsi` file and return an opaque index handle. The caller must free the handle with `bamsi_free`.

### `bamsi_free`

```c
void bamsi_free(bamsi_index_t** idx);
```

Free an index handle. Sets `*idx` to NULL. Safe to call with `*idx == NULL`.

### `bamsi_verify`

```c
bamsi_status_t bamsi_verify(const char* path, int* valid_out);
```

Verify file integrity. Sets `*valid_out` to 1 (valid) or 0 (corrupt).

## Version

```c
const char* bamsi_version(void);
uint16_t bamsi_format_version(void);
```

## Index Information

```c
bamsi_status_t bamsi_get_n_reads(const bamsi_index_t* idx, uint64_t* out);
bamsi_status_t bamsi_get_s_length(const bamsi_index_t* idx, uint64_t* out);
bamsi_status_t bamsi_get_n_windows(const bamsi_index_t* idx, uint32_t* out);
bamsi_status_t bamsi_get_chrom_count(const bamsi_index_t* idx, uint32_t* out);
bamsi_status_t bamsi_get_chrom_name(const bamsi_index_t* idx, uint32_t chrom_idx,
                                     char* buf, size_t buf_len, size_t* out_len);
bamsi_status_t bamsi_is_lossless(const bamsi_index_t* idx, int* out);
```

## Queries

### GlobalCount

```c
bamsi_status_t bamsi_global_count(const bamsi_index_t* idx,
                                   const uint8_t* pattern, size_t pat_len,
                                   uint64_t* count);
```

Pattern encoding: A=0, C=1, G=2, T=3, N=4.

### GlobalExists

```c
bamsi_status_t bamsi_global_exists(const bamsi_index_t* idx,
                                    const uint8_t* pattern, size_t pat_len,
                                    uint64_t threshold, int* exists);
```

### Locate

```c
typedef struct bamsi_locate_result {
    uint32_t chrom_id;
    uint64_t p_min;
    uint64_t p_max;
    uint64_t read_id;
    int      is_reverse;
} bamsi_locate_result_t;

bamsi_status_t bamsi_locate(const bamsi_index_t* idx,
                             const uint8_t* pattern, size_t pat_len,
                             bamsi_locate_result_t* results,
                             size_t max_results, size_t* n_results);
```

### RegionalCount / RegionalExists

```c
bamsi_status_t bamsi_regional_count(const bamsi_index_t* idx,
                                     const uint8_t* pattern, size_t pat_len,
                                     const char* chrom,
                                     uint64_t start, uint64_t end,
                                     uint64_t* count);

bamsi_status_t bamsi_regional_exists(const bamsi_index_t* idx,
                                      const uint8_t* pattern, size_t pat_len,
                                      const char* chrom,
                                      uint64_t start, uint64_t end,
                                      uint64_t threshold, int* exists);
```

### Approximate Queries (V1 Stubs)

```c
// Both return BAMSI_STATUS_NOT_IMPLEMENTED_V1 in v1.0
bamsi_status_t bamsi_approx_locate_hamming(...);
bamsi_status_t bamsi_approx_locate_edit(...);
```

## Usage Example

```c
#include "bamsi/bamsi.h"
#include <stdio.h>

int main(void) {
    bamsi_index_t* idx = NULL;
    if (bamsi_open("genome.bsi", &idx) != BAMSI_STATUS_OK) {
        fprintf(stderr, "Failed to open index\n");
        return 1;
    }

    uint8_t pattern[] = {0, 1, 2, 3};  // ACGT
    uint64_t count;
    bamsi_global_count(idx, pattern, 4, &count);
    printf("ACGT count: %lu\n", (unsigned long)count);

    bamsi_free(&idx);
    return 0;
}
```

## Thread Safety

- All `bamsi_*` query functions are safe for concurrent use on the same index handle.
- `bamsi_open` and `bamsi_free` must not be called concurrently on the same handle.
