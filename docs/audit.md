# BAMSIX Audit Document

**Reference:** Execution Plan Stage 4 §4.3.8

Maps every Contract / Architecture clause to its implementing test.

## Invariant → Test Mapping

| Invariant | Description | Test | Status |
|-----------|-------------|------|--------|
| I1 | Read collection non-empty | `test_tier1_invariants.cpp` CHECK I1 | ✅ |
| I2 | Total ordering determinism | `test_tier1_invariants.cpp` CHECK I2 | ✅ |
| I3 | Ordering hash determinism | `test_tier1_invariants.cpp` CHECK I3 | ✅ |
| I4 | S construction (separators) | `test_tier1_invariants.cpp` CHECK I4 | ✅ |
| I5 | SA permutation correctness | `test_tier1_invariants.cpp` CHECK I5 | ✅ |
| I6 | BWT derivation correctness | `test_tier1_invariants.cpp` CHECK I6 | ✅ |
| I7 | LF property | `test_tier1_invariants.cpp` CHECK I7 | ✅ |
| I8 | Backward search vs brute-force | `test_tier1_invariants.cpp` CHECK I8 | ✅ |
| I9 | B_read popcount = N | `test_tier1_invariants.cpp` CHECK I9 | ✅ |
| I10 | B_read select matches readStarts | `test_tier1_invariants.cpp` CHECK I10 | ✅ |
| I11 | Window coverage (all reads) | `test_tier1_invariants.cpp` CHECK I11 | ✅ |
| I12 | Window ordering | `test_tier1_invariants.cpp` CHECK I12 | ✅ |
| I13 | SA sample correctness | `test_tier1_invariants.cpp` CHECK I13 | ✅ |
| I14 | Separator detection via B_read only | `test_tier1_invariants.cpp` CHECK I14 | ✅ |
| I15 | Locate correctness | `test_tier1_invariants.cpp` CHECK I15 | ✅ |

## Gate Test → Milestone Mapping

| Gate Test | Milestone | Test | Status |
|-----------|-----------|------|--------|
| Build → Verify → Load → Count → Locate | V1 | `test_v1_roundtrip.cpp` | ✅ |
| 100 patterns, FM vs brute-force | V2 | `test_v2_fm_correctness.cpp` | ✅ |
| .bsi checksum verification | V1 | `VerifyBsi()` + CLI verify | ✅ |
| C ABI test program | V5 / Stage 3 | `test_c_abi.c` | ✅ |

## Contract Clause Coverage

| Contract § | Description | Implementation | Test |
|------------|-------------|----------------|------|
| §0.1 | System overview | README.md + docs/algorithms.md | — |
| §2.1 | BAM ingestion | `src/ingest/` | TIER 1 I1 |
| §2.2 | Read ordering | `src/ordering/` | TIER 1 I2, I3 |
| §2.3 | Sequence concatenation | `src/seqbuilder/` | TIER 1 I4 |
| §2.4 | SA-IS construction | `src/sais/` | TIER 1 I5, I6, I7 |
| §2.5 | FM-index | `src/fmindex/` | TIER 1 I7, I8 |
| §2.6 | Bitvectors | `src/bitvectors/` | TIER 1 I9, I10 |
| §2.7 | S_seq encoding | `src/streamencode/` | V1 round-trip |
| §2.8 | S_qual encoding | `src/streamencode/` | V1 round-trip |
| §2.9 | S_meta / S_map encoding | `src/streamencode/` | V1 round-trip |
| §3.1 | GlobalCount | `src/query/` | V2 gate (100 patterns) |
| §3.2 | GlobalExists | `src/query/` | V1 round-trip |
| §4.1 | Locate | `src/query/ + mapping/` | V2 gate (15 locates) |
| §4.2 | RegionalCount | `src/query/` | CLI smoke test |
| §4.3 | RegionalExists | `src/query/` | CLI smoke test |
| §4.5 | Approximate search (stub) | C ABI stubs | C ABI test |
| §5.1 | Separator detection | `src/mapping/` | TIER 1 I14 |
| §5.2 | CIGAR mapping | `src/mapping/` | Multi-op reads |
| §7 | .bsi format | `src/seal/ + format/` | V1 round-trip |
| §9.2 | `bamsix info` | CLI dispatch | Smoke test |
| §10.2 | CLI subcommands (9) | `src/cli/dispatch.cpp` | All working |
| §10.3 | C ABI | `src/cabi/cabi.cpp` | `test_c_abi.c` (43/43) |
| §10.6 | Provenance | Checksums in .bsi | `bamsix verify` |

## Architecture Clause Coverage

| Architecture § | Description | Implementation | Test |
|---------------|-------------|----------------|------|
| §3 Data types | Type definitions | `include/bamsix/types.hpp` | Compile-time |
| §4.6 Rank convention | Closed-interval B_read, half-open Occ | Separate functions | I15 |
| §4.8 Windows | Window construction | `src/windows/` | I11, I12 |
| §5.1 MapOccurrence | Position mapping | `src/mapping/` | V2 Locate |
| §5.2 CIGAR mapping | cigar_ref_pos | `src/mapping/` | Multi-CIGAR test |
| §7 File format | .bsi binary layout | `src/seal/ + format/` | Round-trip |
| §9.1 Validation | TIER 1 invariants | `tests/test_tier1_invariants.cpp` | 20/20 |
| §13.1 Source layout | Module structure | `src/` directories | Build system |

## Test Summary

| Test Suite | Tests | Status |
|-----------|-------|--------|
| TIER 1 Invariants (10 reads) | 20 | ✅ All pass |
| TIER 1 Invariants (10K reads) | 20 | ✅ All pass |
| V1 Round-Trip Gate | 15 | ✅ All pass |
| V2 FM Correctness Gate | 115 | ✅ All pass |
| C ABI Test | 43 | ✅ All pass |
| CLI Smoke Tests | 9 commands | ✅ All work |
