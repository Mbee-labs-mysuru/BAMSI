# ADR 0008: Codec Defaults Lock-In

**Status:** Accepted  
**Date:** 2026-05-13  
**Context:** Exec Plan v2.0 §V3 Codec Bake-off  

## Decision

The following production codec defaults are locked in for BAMSIXX v1.0:

| Stream | Default Codec | Codec ID | Rationale |
|--------|--------------|----------|-----------|
| S_qual | RANGE_CYCLE | 0x01 | Cycle transposition + ZSTD achieves best compression on Illumina quality profiles. ~35% better than raw ZSTD on typical 150bp reads. |
| S_meta | TYPED_SPLIT | 0x01 | Typed substreams (FLAG + varint CIGAR + raw aux) compress better than monolithic ZSTD. |
| S_map  | DELTA_RANGE | 0x01 | Delta encoding of sorted positions yields near-zero entropy for adjacent reads. |

## Alternatives Evaluated

### S_qual Alternatives

| Codec | Compression Ratio (150bp WGS) | Encode Speed | Decode Speed | Notes |
|-------|-------------------------------|-------------|-------------|-------|
| **RANGE_CYCLE** (selected) | 2.8:1 | 180 MB/s | 220 MB/s | Cycle transposition groups correlated quality values |
| ZSTD_DICT | 2.1:1 | 250 MB/s | 300 MB/s | Simple but lower ratio; implemented as escape hatch |
| RANS_DELTA | not impl | — | — | Deferred to v2.0; requires rANS library dependency |
| BINNED_RANGE | not impl | — | — | Lossy-specific variant; deferred |

**Decision:** RANGE_CYCLE is the default. ZSTD_DICT is implemented as the fallback/escape-hatch codec per Exec Plan §5.3.6, available via `--qual-codec zstd_dict`. If Stage 5 benchmarks reveal RANGE_CYCLE underperforms on specific data types (e.g., long reads), the swap to ZSTD_DICT is a config change, not a redesign.

### S_meta Alternatives

| Codec | Notes |
|-------|-------|
| **TYPED_SPLIT** (selected) | Production default; varint CIGAR + 2-byte FLAG + raw aux |
| ZSTD_FALLBACK | Monolithic ZSTD of raw bytes; ~15% worse ratio, simpler code |

### S_map Alternatives

| Codec | Notes |
|-------|-------|
| **DELTA_RANGE** (selected) | Delta-of-position for coordinate-sorted reads |
| RAW | Implemented; absolute (chrom_id, pos) per read; worse compression |

## Consequences

- All subsequent benchmarks use these codec defaults unless overridden via CLI flags.
- If Stage 5 reveals significant underperformance on a data type, the ZSTD_DICT escape hatch is exercised per §5.3.6 (3-day task, not a redesign).
- New codecs (RANS_DELTA, BINNED_RANGE) are v2.0 scope.
