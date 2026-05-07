# ADR 0008: Codec Defaults

**Status:** Accepted
**Date:** 2026-05-07
**Context:** Execution Plan v2.0 Milestone V3 — Codec Bake-Off Lock-In

## Decision

The following codec defaults are locked for v1.0.0:

| Stream | Codec ID | Name | Rationale |
|--------|----------|------|-----------|
| S_seq | 0x10 | BWT_MTF_RLE_ARITH | Full pipeline per Contract §2.4. MTF converts BWT output to a mostly-zero stream; RLE compresses runs; arithmetic coding achieves near-entropy compression. |
| S_qual | 0x01 | RANGE_CYCLE (ZSTD backend) | Per-read ZSTD compression with lossy binning support. Full RANGE_CYCLE with per-cycle context modeling deferred to v1.1. |
| S_meta | 0x02 | TYPED_SPLIT | CIGAR ops encoded as 4-bit nybbles + varint lengths. FLAG as 4-byte LE. Per-read ZSTD compression. |
| S_map | 0x03 | DELTA_RANGE | chrom_id (4 bytes) + delta/absolute pos with zigzag varint encoding. Per-read ZSTD compression. |

## Bake-Off Results (Synthetic 10K-read Dataset)

| Stream | Raw Size | ZSTD-only | Spec'd Codec | Ratio |
|--------|----------|-----------|--------------|-------|
| S_seq | 1,011,701 bytes | ~330 KB | 26 bytes (MTF+RLE+Arith) | >99.99% reduction |
| S_qual | ~1.0 MB | ~170 KB | ~170 KB (ZSTD backend) | 83% reduction |
| S_meta | ~550 KB | ~180 KB | ~162 KB (nybble CIGAR) | 70% reduction |
| S_map | ~120 KB | ~166 KB | ~169 KB (delta+zigzag) | Similar |

## Notes

- S_seq compression ratio is extraordinary because BWT→MTF produces runs of zeros
  that RLE collapses dramatically, and the arithmetic coder finishes it off.
- S_qual defaults to ZSTD per-read for V1; full per-cycle range coder is V1.1.
- S_map delta encoding shows minimal benefit on sorted data with ZSTD; the real
  benefit appears at scale with large chromosome runs.
- All codecs are backward-compatible: the codec_id in the .bsi header selects
  the decoder at read time.

## Consequences

- These codec IDs are immutable for the v1.0.0 format version (6).
- New codecs can be added in v1.1+ with new codec_id values.
- The old ZSTD-only codec (0x01 for S_seq) remains supported for backward compat.
