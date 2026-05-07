# Gate 3 Review — Query Path Hardening Complete

**Date:** 2026-05-07
**Reviewers:** Tech lead
**Status:** PASS

## Checklist (Contract §10.2 line-by-line)

| CLI Command | Status | Notes |
|-------------|--------|-------|
| `bamsi build` | ✅ | All pipeline stages |
| `bamsi count` | ✅ | --json, --strand |
| `bamsi exists` | ✅ | --json, --strand |
| `bamsi locate` | ✅ | --json, --strand, --sort-output |
| `bamsi region-count` | ✅ | --json, --strand, --region chr:a-b |
| `bamsi region-exists` | ✅ | --json, --strand, --threshold, --region chr:a-b |
| `bamsi reconstruct` | ✅ | --output, --streams, --read-ids, --allow-lossy |
| `bamsi info` | ✅ | --json (all 30+ fields per Contract §9.2) |
| `bamsi verify` | ✅ | --json, section checksums, global checksum |

## C ABI (Contract §10.3)

- [x] `bamsi_version`, `bamsi_format_version`
- [x] `bamsi_open`, `bamsi_free`
- [x] `bamsi_verify`
- [x] `bamsi_get_n_reads`, `bamsi_get_s_length`, `bamsi_get_n_windows`
- [x] `bamsi_get_chrom_count`, `bamsi_get_chrom_name`
- [x] `bamsi_is_lossless`
- [x] `bamsi_global_count`, `bamsi_global_exists`
- [x] `bamsi_locate`
- [x] `bamsi_regional_count`, `bamsi_regional_exists`
- [x] `bamsi_approx_locate_hamming` → NOT_IMPLEMENTED_V1
- [x] `bamsi_approx_locate_edit` → NOT_IMPLEMENTED_V1

## Error Handling

Every ErrorCode has at least one code path that generates it.

## Decision

**PASS** — proceed to Stage 4.
