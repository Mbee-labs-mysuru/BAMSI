# BAMSIX Clinical Operations Guide

**Reference:** Contract §6, §10.6

## Provenance Verification

Every `.bsi` file contains verifiable provenance:

| Field | Purpose |
|-------|---------|
| `source_manifest_hash` | SHA-256 of input file manifest |
| `ordering_hash` | SHA-256 of deterministic read ordering |
| `build_timestamp_utc` | When the index was built |
| `source_file_count` | Number of input BAM files |
| `global_checksum` | xxHash64 of entire file |

### Verification Workflow

```bash
# 1. Verify file integrity
bamsix verify --index patient.bsi
# Expected: PASS

# 2. Inspect provenance metadata
bamsix info --index patient.bsi --json | jq '.source_manifest_hash, .ordering_hash'

# 3. Cross-reference with build log
# The source_manifest_hash can be recomputed from the original BAM files
# to confirm the index matches the expected input.
```

## Audit Trail Integration

### Build Audit

For clinical deployments, record:

1. Input BAM file paths and SHA-256 checksums
2. BAMSIX version (`bamsix version`)
3. Build command and parameters
4. Output `.bsi` file SHA-256
5. `bamsix info --json` output (full metadata dump)
6. `bamsix verify` result

### Query Audit

For auditable query results:

```bash
# Record query with full provenance
bamsix count --index patient.bsi --pattern BRCA1MOTIF --json
# Output: {"count":42}
# Log: timestamp, user, index SHA-256, pattern, result
```

## Lossy Mode Caveats

When `is_lossless = 0` (quality scores were binned):

- **Query results are unaffected** — pattern counting uses sequence data only
- **Reconstructed quality scores differ** from originals
- `bamsix reconstruct` warns without `--allow-lossy` flag
- `bamsix info` surfaces the lossy condition prominently

### Checking Lossy Status

```bash
bamsix info --index patient.bsi --json | jq '.is_lossless'
# true = lossless, false = lossy
```

## No Implicit Network Access

BAMSI performs **zero network operations** at any point:

- Index building reads only local files
- Queries operate entirely on the local `.bsi` file
- No telemetry, no update checks, no remote dependencies

This can be verified by running under network isolation:

```bash
# Linux: network namespace isolation
unshare -n bamsix count --index patient.bsi --pattern ACGT
# Must succeed (no network dependency)
```

## Deterministic Builds

The same input BAM(s) produce byte-identical `.bsi` files across:

- Different machines (same architecture)
- Different runs on the same machine
- Docker container vs bare-metal

This is guaranteed by:
- Deterministic read ordering (I2, I3)
- No hash-map iteration in output paths
- No floating-point in the build pipeline
- No timestamp-dependent ordering

### Verification

```bash
# Build twice, compare
bamsix build input.bam -o first.bsi
bamsix build input.bam -o second.bsi
sha256sum first.bsi second.bsi
# Both checksums must match
```

## Clinical Deployment Checklist

- [ ] Verify BAMSIX version matches validated version
- [ ] Verify input BAM integrity (checksums)
- [ ] Build index with explicit parameters
- [ ] Run `bamsix verify` on output
- [ ] Record `bamsix info --json` for audit trail
- [ ] Verify deterministic build (rebuild and compare)
- [ ] Test under network isolation
- [ ] Store build logs with patient identifier linkage
