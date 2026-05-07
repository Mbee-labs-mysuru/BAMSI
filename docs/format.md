# BAMSI .bsi Format Specification

**Format version:** 6
**Reference:** Architecture §7

## Overview

The `.bsi` file is a binary container storing a compressed FM-index, auxiliary streams, bitvectors, window tables, and read metadata. All multi-byte integers are little-endian.

## File Structure

```
┌──────────────────────────────────┐
│ Header (variable)                │
├──────────────────────────────────┤
│ Stream: S_seq (BWT)              │
│ Stream: S_qual (quality scores)  │
│ Stream: S_meta (CIGAR + FLAG)    │
│ Stream: S_map (chrom_id + pos)   │
├──────────────────────────────────┤
│ FM-Index Section                 │
├──────────────────────────────────┤
│ Bitvector Section                │
├──────────────────────────────────┤
│ Window Table                     │
├──────────────────────────────────┤
│ Directory Section                │
├──────────────────────────────────┤
│ Read Metadata Section            │
├──────────────────────────────────┤
│ Footer                           │
└──────────────────────────────────┘
```

## Header

| Offset | Size | Field | Description |
|--------|------|-------|-------------|
| 0 | 4 | `magic` | `"BMSI"` (0x42, 0x4D, 0x53, 0x49) |
| 4 | 2 | `version` | Format version (currently 6) |
| 6 | 16 | `bamsi_version` | Null-terminated version string |
| 22 | 1 | `host_os_id` | Build OS identifier |
| 23 | 1 | `cpu_arch_id` | CPU architecture |
| 24 | 8 | `build_timestamp_utc` | Unix timestamp |
| 32 | 1 | `is_lossless` | 1 = lossless, 0 = lossy |
| 33 | 4 | `source_file_count` | Number of input BAM files |
| 37 | 32 | `source_manifest_hash` | SHA-256 of input file manifest |
| 69 | 32 | `ordering_hash` | SHA-256 of read ordering |
| 101 | 8 | `S_length` | Length of concatenated sequence S |
| 109 | 8 | `N_reads` | Number of reads |
| 117 | 4 | `N_windows` | Number of windows |
| 121 | 4 | `sample_step_s` | SA sampling step |
| 125 | 1 | `has_isa_samples` | ISA samples present |
| 126 | 4 | `sample_step_s_prime` | ISA sampling step |
| 130 | 1 | `enable_sarange` | SARange enabled |
| 131 | 1 | `sarange_variant` | SARange variant ID |
| 132 | 1 | `shared_bwt` | Shared BWT flag |
| 133 | 1 | `enable_bidirectional` | Bidirectional FM-index |
| 134 | 1 | `recommended_seed_length` | Recommended seed length |
| 135 | 8 | `window_size_T` | Window size in S-characters |
| 143 | 1 | `entropy_order_k` | Entropy coding order |
| 144 | 1 | `qual_codec_id` | Quality stream codec |
| 145 | 1 | `qual_lossy_bins` | Lossy quality bin count (0=lossless) |
| 146 | 1 | `meta_codec_id` | Metadata stream codec |
| 147 | 1 | `map_codec_id` | Mapping stream codec |
| 148 | 1 | `strand_mode` | 0=StrandComplete, 1=SingleStrand |
| 149 | 8 | `sentinel_row` | BWT row of sentinel |
| 157 | 4 | `chrom_count` | Number of chromosomes |
| 161+ | var | `chrom_name_table` | Per-chrom: id(4) + name_len(4) + name |

Following the chrom table:
| Size | Field | Description |
|------|-------|-------------|
| 4 | `seq_block_size` | Sequence block size |
| 4 | `qual_block_size` | Quality block size |
| 1 | `allow_parallel_sa` | Parallel SA-IS flag |
| 1 | `reference_based_encoding` | Reference-based encoding |
| 32 | `reference_sha256` | Reference SHA-256 |
| 4 | `flags` | Reserved flags |

## Stream Section Format

Each of the four streams (S_seq, S_qual, S_meta, S_map) uses:

| Size | Field | Description |
|------|-------|-------------|
| 8 | `payload_length` | Compressed payload size in bytes |
| 1 | `codec_id` | Codec identifier |
| var | `payload` | Compressed stream data |
| 8 | `section_checksum` | Section xxHash64 |

## FM-Index Section

| Size | Field | Description |
|------|-------|-------------|
| 8 | `bwt_length` | BWT length (= |S|, sentinel-stripped) |
| var | `bwt_bytes` | Raw BWT bytes |
| 56 | `C_array` | 7 × uint64 cumulative counts |
| 8 | `occ_length` | Occurrence table size |
| var | `occ_bytes` | Serialized occurrence table |
| 8 | `sa_count` | Number of SA samples |
| var | `sa_samples` | SA samples (sa_count × uint64) |
| 8 | `checksum` | Section checksum |

## Bitvector Section

| Size | Field | Description |
|------|-------|-------------|
| 8 | `b_read_length` | B_read serialized size |
| var | `b_read_bytes` | B_read bitvector |
| 8 | `b_window_length` | B_window serialized size |
| var | `b_window_bytes` | B_window bitvector |
| 8 | `checksum` | Section checksum |

## Window Table

Per window (repeated `N_windows` times):

| Size | Field | Description |
|------|-------|-------------|
| 4 | `chrom_id` | Chromosome index |
| 8 | `l` | S-position start |
| 8 | `r` | S-position end |
| 8 | `first_read_id` | First read in window |
| 8 | `last_read_id` | Last read in window |
| 8 | `genomic_start` | Genomic start (1-based) |
| 8 | `genomic_end` | Genomic end (1-based) |

Followed by 8-byte section checksum.

## Read Metadata Section

| Size | Field | Description |
|------|-------|-------------|
| 8 | `n_reads` | Number of reads |
| var | per-read data | See below |
| 8 | `checksum` | Section checksum |

Per read:

| Size | Field | Description |
|------|-------|-------------|
| 4 | `chrom_id` | Chromosome index |
| 8 | `pos` | 1-based alignment position |
| 4 | `seq_len` | Sequence length |
| 4 | `n_ops` | Number of CIGAR operations |
| var | CIGAR ops | op(1 byte) + len(4 bytes) per op |

## Footer

| Size | Field | Description |
|------|-------|-------------|
| 8 | `global_checksum` | xxHash64 of all preceding bytes |
| 4 | `footer_magic` | 0xB5110000 |

## Integrity Verification

`bamsi verify` computes xxHash64 over all bytes before the footer and compares against the stored global checksum. Any single-byte modification is detected.
