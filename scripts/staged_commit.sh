#!/usr/bin/env bash
# BAMSI Staged Git Commit Script
# Commits all changes in logical stages per Contract v3.3 / Architecture v4.3 / Execution Plan v2.0
#
# Usage: bash scripts/staged_commit.sh
#
# WARNING: Review each commit before running. This script will execute git add + git commit.
# Run with DRY_RUN=1 to preview without committing: DRY_RUN=1 bash scripts/staged_commit.sh

set -euo pipefail
cd "$(git rev-parse --show-toplevel)"

DRY_RUN="${DRY_RUN:-0}"

do_commit() {
    local msg="$1"
    shift
    echo ""
    echo "════════════════════════════════════════════════════════════════"
    echo "COMMIT: ${msg%%$'\n'*}"
    echo "════════════════════════════════════════════════════════════════"
    for f in "$@"; do
        if [[ -e "$f" ]]; then
            echo "  + $f"
            [[ "$DRY_RUN" == "0" ]] && git add "$f"
        else
            echo "  (skipped, not found: $f)"
        fi
    done
    if [[ "$DRY_RUN" == "0" ]]; then
        git commit -m "$msg" --allow-empty || echo "  (nothing to commit)"
    else
        echo "  [DRY RUN — would commit]"
    fi
}

# ══════════════════════════════════════════════════════════════════════════════
# STAGE 0 — Bootstrap: Build system, dependencies, project scaffold
# Per Execution Plan §1: Provision build system, pin dependencies, scaffold
# ══════════════════════════════════════════════════════════════════════════════

do_commit "BAMSI-stage-0: Bootstrap build system, external dependencies, and project scaffold

Stage 0 (Execution Plan §1): Project bootstrap and infrastructure setup.

Changes:
- CMakeLists.txt: Root build system with C++20 requirement, static linking of all
  external dependencies (htslib, libsais, sdsl-lite, zstd, xxHash, GoogleTest)
- .gitmodules: External dependency submodule configuration
- external/: All bundled dependencies (htslib, libsais, sdsl-lite, zstd, xxHash)
  pinned to specific versions per ADR 0002 (dependency pinning policy)
- include/bamsi/: Public API headers (bamsi.h, types.hpp, version.hpp, config.hpp,
  result.hpp, status.hpp, cli/dispatch.hpp)
- src/CMakeLists.txt: bamsi-core static library target definition
- src/cli/main.cpp: CLI entry point
- src/cli/app.cpp, help.cpp, version.cpp: CLI infrastructure
- src/cabi/version.cpp, version.hpp: Version reporting

Contract compliance:
- §10.1: C++ (≥ C++20) implementation language
- §10.3: Stable C ABI shim layer via extern \"C\"
- Execution Plan §1.3: All dependencies bundled and statically linked

Ref: Contract v3.3 §10, Architecture v4.3 §13, Execution Plan v2.0 §1" \
    .gitmodules \
    CMakeLists.txt \
    src/CMakeLists.txt \
    include/bamsi/bamsi.h \
    include/bamsi/types.hpp \
    include/bamsi/cli/dispatch.hpp \
    src/cli/app.cpp \
    src/cli/dispatch.cpp \
    src/cabi/version.cpp \
    src/cabi/version.hpp \
    external/htslib \
    external/libsais \
    external/sdsl-lite \
    external/xxHash \
    external/zstd \
    external/deps/

# ══════════════════════════════════════════════════════════════════════════════
# STAGE 1 — Core Index Pipeline: Ingest → Order → SeqBuild → SA-IS → FM-Index
# Per Execution Plan §2 and Architecture §4.1–§4.6
# ══════════════════════════════════════════════════════════════════════════════

do_commit "BAMSI-stage-1: Core index build pipeline — Ingest, Order, SeqBuilder, SA-IS, FM-Index

Stage 1 (Execution Plan §2): Implement the core build pipeline stages 1–5b.

Build Pipeline Stages:
- src/ingest/: Stage 1 — BAM record ingestion via htslib. Filters records per
  inclusion rule §0.2 (mapped, primary, POS≥1). Extracts seq, qual, CIGAR, FLAG,
  chrom, pos. Computes source_manifest_hash (SHA-256 of sorted input paths).
- src/ordering/: Stage 2 — Coordinate-order read sorting. Produces ordering_hash
  (xxHash64 of sorted read-id sequence) for determinism verification (§0.7).
- src/seqbuilder/: Stage 3 — Concatenated sequence S construction. Encodes reads
  over alphabet Σ={A,C,G,T,N} with '#' separators per §0.4. Records read boundaries
  for bitvector construction.
- src/sais/: Stage 4 — Suffix array construction via SA-IS (Nong, Zhang & Chan 2009).
  Produces SA, BWT, sentinel_row. Supports both 32-bit and 64-bit SA via libsais.
  BWT derivation: BWT[i] = S[(SA[i]-1+|S|) mod |S|] per Architecture §4.4.
- src/fmindex/: Stage 5b — FM-Index construction. Implements:
  * C array (7-element cumulative frequency table, indexed by stored code)
  * OccTable (checkpointed rank over BWT, block size 64)
  * Backward search: O(|P|) exact pattern matching
  * Locate: SA-sample-based position recovery via LF-mapping
  * SA samples at configurable step s (default 64)

Contract compliance:
- §0.3: Frozen alphabet encoding (A=0,C=1,G=2,T=3,N=4,#=5,$=6)
- §0.4: Separator policy — '#' between reads, never in query patterns
- §0.5: Sentinel policy — conceptual '$', never stored, sentinel_row tracked
- §0.7: Deterministic build guarantee
- §3.1: Two distinct rank APIs (BWT-rank half-open, bitvector-rank closed)
- §4.1: FM backward search correctness (Invariants I1-I3)

Ref: Contract v3.3 §0–§3, Architecture v4.3 §4.1–§4.6" \
    src/ingest/ingest.cpp \
    src/ingest/ingest.hpp \
    src/ordering/ordering.cpp \
    src/ordering/ordering.hpp \
    src/seqbuilder/seqbuilder.cpp \
    src/seqbuilder/seqbuilder.hpp \
    src/sais/sais.cpp \
    src/sais/sais.hpp \
    src/fmindex/fmindex.cpp \
    src/fmindex/fmindex.hpp \
    src/cli/build.cpp \
    src/cli/build.hpp

# ══════════════════════════════════════════════════════════════════════════════
# STAGE 2 — Codec Pipeline: S_seq, S_qual, S_meta, S_map + Windows + Bitvectors + Seal
# Per Execution Plan §3 and Architecture §4.5, §4.7, §4.8, §4.9, §7
# ══════════════════════════════════════════════════════════════════════════════

do_commit "BAMSI-stage-2: Production codecs, auxiliary streams, windows, bitvectors, and .bsi sealing

Stage 2 (Execution Plan §3): Four production codecs, window/bitvector construction, and .bsi sealing.

Stream Codecs (Architecture §4.5, §4.7):
- src/streamencode/mtf_rle_arith.cpp: S_seq codec — BWT→MTF→RLE→Arithmetic (§2.4).
  Move-to-Front transform converts BWT cluster runs to zero-heavy stream,
  Run-Length Encoding compresses zero runs, Adaptive Arithmetic coding achieves
  near-optimal |S|·H_k(S) + o(|S|) bits. Entropy order k∈[4,8], default k=6.
- src/streamencode/streamencode.cpp: S_qual (RANGE_CYCLE, §2.7), S_meta
  (TYPED_SPLIT with CIGAR nybble + FLAG + ZSTD, §2.8), S_map (DELTA_RANGE
  with zigzag varint + ZSTD, §2.9). Lossy quality binning via qual_lossy_bins.

Window & Bitvector Construction:
- src/windows/: Stage 7 — Genomic window table. Each window records (chrom_id,
  S-range [l,r], read range, genomic coordinates). Windows respect chromosome
  boundaries and separator inclusion rules (§0.4, §3.4).
- src/bitvectors/: Stage 8 — B_read (marks read boundaries in S) and B_window
  (marks window boundaries). Rank/Select support via sdsl-lite.

Sealing & Format:
- src/seal/: Stage 9-10 — .bsi file serialization with per-section xxHash64
  checksums and global footer checksum. Writes header (42+ fields), stream
  sections, FM-index section, bitvector/window/directory sections.
- src/format/: .bsi file deserialization. ReadBsi() for index loading,
  VerifyBsi() for checksum verification.
- src/mapping/: Coordinate mapping encode/decode for S_map stream.

Contract compliance:
- §2.4: Mandatory BWT→MTF→RLE→Arithmetic pipeline
- §2.7: S_qual codec catalogue (RANGE_CYCLE default)
- §2.8: S_meta TYPED_SPLIT with CIGAR/FLAG substreams
- §2.9: S_map DELTA_RANGE exploiting sorted-position structure
- §2.10: Provenance fields (build_timestamp, host_os_id, cpu_arch_id)
- §5.1: Space bound |S|·H_k(S) + O(|S|)
- §7: .bsi binary format with checksums

Ref: Contract v3.3 §2, §5, §7, Architecture v4.3 §4.5–§7" \
    src/streamencode/streamencode.cpp \
    src/streamencode/streamencode.hpp \
    src/streamencode/mtf_rle_arith.cpp \
    src/streamencode/mtf_rle_arith.hpp \
    src/windows/windows.cpp \
    src/windows/windows.hpp \
    src/bitvectors/bitvectors.cpp \
    src/bitvectors/bitvectors.hpp \
    src/seal/seal.cpp \
    src/seal/seal.hpp \
    src/format/format.cpp \
    src/format/format.hpp \
    src/mapping/mapping.cpp \
    src/mapping/mapping.hpp

# ══════════════════════════════════════════════════════════════════════════════
# STAGE 3 — Query Engine + CLI Hardening + C ABI
# Per Execution Plan §3 and Architecture §6, Contract §4, §10
# ══════════════════════════════════════════════════════════════════════════════

do_commit "BAMSI-stage-3: Query engine, CLI hardening, C ABI, and full subcommand surface

Stage 3 (Execution Plan §3): Query operations, CLI polish, and C ABI compliance.

Query Engine (Architecture §6, Contract §4):
- src/query/: GlobalCount, GlobalExists, Locate, RegionalCount, RegionalExists.
  Strand-complete by default (searches both P and rc(P) unless P is palindromic).
  Pattern validation: EMPTY_PATTERN for len=0, INVALID_PATTERN for non-Σ chars.
  Approximate query stubs return NOT_IMPLEMENTED_V1 per §4.5.

CLI Surface (Contract §10.2 — all 10 subcommands):
- build, count, exists, locate, region-count, region-exists, reconstruct, info, verify
- --strand single|complete for all query commands
- --region chr:a-b parser for regional queries
- --threshold for region-exists
- --lossy-bins, --isa-step, --seq-block-size, --qual-block-size, --seed-length,
  --reference, --lossless, --lossy, --parallel-sa build flags
- --output mandatory for reconstruct (no silent default per §0.8)
- ParseInputs correctly skips all value-carrying flags
- bamsi info --json: 36+ fields including provenance, codec IDs, hashes

C ABI (Architecture §13.1, Contract §10.3):
- src/cabi/cabi.cpp: Stable extern \"C\" shim layer — bamsi_open, bamsi_count,
  bamsi_exists, bamsi_locate, bamsi_close, bamsi_version, bamsi_free_result.
  Thread-safe, null-safe, returns structured error codes.

Contract compliance:
- §0.6: Strand-complete operational default
- §0.8: No silent failures — all error paths return structured ErrorCode
- §4.1-§4.4: GlobalCount, GlobalExists, Locate, RegionalCount
- §4.5: NOT_IMPLEMENTED_V1 for approximate queries
- §4.6: Output ordering (streaming default, sorted on --sort-output)
- §9: Operational guarantees — determinism, no implicit network, provenance
- §10.2: All 10 CLI subcommands implemented
- §10.3: Stable C ABI for cross-language consumption

Ref: Contract v3.3 §4, §9, §10, Architecture v4.3 §6, §10, §13" \
    src/query/query.cpp \
    src/query/query.hpp \
    src/cabi/cabi.cpp \
    include/bamsi/bamsi.h \
    include/bamsi/types.hpp \
    include/bamsi/cli/dispatch.hpp \
    src/cli/dispatch.cpp \
    CHANGELOG.md \
    README.md

# ══════════════════════════════════════════════════════════════════════════════
# STAGE 4 — Validation Campaign: All test suites + audit + dependency freeze
# Per Execution Plan §4 and Architecture §8/§9
# ══════════════════════════════════════════════════════════════════════════════

do_commit "BAMSI-stage-4: Validation campaign — regression tests, error sweep, audit, dependency freeze

Stage 4 (Execution Plan §4): Complete validation campaign with 80+ tests.

Test Suites:
- tests/unit/version_test.cpp: Version reporting unit test
- tests/integration/: CLI smoke tests (7), version integration test
- tests/test_tier1_invariants.cpp: Invariant I1-I15 structural tests
- tests/test_tier2_integration.cpp: TIER 2 exhaustive integration (19 tests) —
  C ABI, codec correctness, provenance fields, error paths
- tests/test_error_sweep.cpp: Error handling sweep (19 tests) — every ErrorCode
  path exercised (CLI-level + C ABI null/invalid inputs)
- tests/test_findings_fixes_regression.cpp: F&F regression suite (24 tests) —
  F1 space bound, F2 BASE tier queries, F3 no S[pos] access, F4 block size
  independence, F5 deterministic build, F6 window consistency, S1 two-rank APIs,
  S3 parallel SA-IS identity, S6 strand semantics, S7 partial reconstruct,
  C1 output ordering, lossy-mode obligations (5), provenance/audit (5),
  codec correctness ablations (3)
- tests/test_v1_roundtrip.cpp, test_v2_fm_correctness.cpp, test_reader.cpp,
  verify_fm.cpp: FM-index correctness and format reader tests
- tests/test_c_abi.c: C ABI compliance test (pure C)
- tests/gen_test_bam.c, gen_10k_bam.c: Synthetic BAM generators for testing
- tests/CMakeLists.txt: 6 test executables, 12 CTest targets with labels

Documentation:
- docs/audit.md: Full contract clause → test mapping audit document
- docs/dependency_freeze.md: All dependency versions locked for v1.0.0 release
- docs/gates/gate-1.md through gate-3.md: Stage gate review documents
- docs/decisions/0008-codec-defaults.md: ADR for codec parameter choices

Contract compliance:
- §0.8: Every ErrorCode has a triggering test (no silent failures)
- §5.3: Validation matrix — TIER 1 (production) + TIER 2 (CI/synthetic)
- §9.3: Determinism, provenance, and reproducibility verified
- Execution Plan §4.4: Findings & Fixes F1-F6, S1-S9, C1-C3 all regressed

Ref: Contract v3.3 §5.3, §0.8, Architecture v4.3 §8-§9, Execution Plan v2.0 §4" \
    tests/CMakeLists.txt \
    tests/unit/version_test.cpp \
    tests/integration/version_cli_test.cpp \
    tests/test_tier1_invariants.cpp \
    tests/test_tier2_integration.cpp \
    tests/test_error_sweep.cpp \
    tests/test_findings_fixes_regression.cpp \
    tests/test_v1_roundtrip.cpp \
    tests/test_v2_fm_correctness.cpp \
    tests/test_reader.cpp \
    tests/test_c_abi.c \
    tests/verify_fm.cpp \
    tests/gen_test_bam.c \
    tests/gen_10k_bam.c \
    docs/audit.md \
    docs/dependency_freeze.md \
    docs/decisions/0008-codec-defaults.md \
    docs/gates/gate-1.md \
    docs/gates/gate-2.md \
    docs/gates/gate-3.md

# ══════════════════════════════════════════════════════════════════════════════
# STAGE 5 — V5 ENHANCED Tier: ISA samples, SARange wavelet tree, Reverse FM-Index
# Per Architecture §4.4 step 5, §4.6.7, §5.3
# ══════════════════════════════════════════════════════════════════════════════

do_commit "BAMSI-stage-5: ENHANCED tier — ISA samples, SARange wavelet tree, Reverse FM-Index

V5 Features (Architecture §4.4, §4.6.7, §5.3): ENHANCED tier succinct data structures.

ISA Samples (Architecture §4.4 step 5):
- src/sais/sais.{hpp,cpp}: ComputeISASamples() — inverts SA to get ISA, then
  samples at intervals of s'. ISA[SA[j]]=j for all j, then ISA_samples[k]=ISA[k*s'].
  Full ISA is discarded after sampling (Architecture lifecycle invariant).
- src/fmindex/fmindex.{hpp,cpp}: SetISASamples(), LocateBidir() — bidirectional
  locate using min(SA-walk, ISA-walk) for O(min(s,s')) per-occurrence cost.
  HasISASamples(), ISASamples(), ISAStep() accessors.
- CLI: --isa-step N flag, header fields has_isa_samples + sample_step_s_prime.

SARange Wavelet Tree (Architecture §5.3 — ENHANCED Tier):
- src/sarange/sarange.{hpp,cpp}: New module implementing a wavelet tree over
  SA_samples for O(log(|S|/s)) regional count queries. Supports:
  * Build(sa_samples, max_val) — MSB-first balanced wavelet tree construction
  * RangeCount(lo, hi, val_lo, val_hi) — count elements in SA interval [lo,hi)
    whose text positions fall within [val_lo, val_hi]. Uses iterative stack-based
    traversal with early pruning on fully-contained/fully-excluded nodes.
  * Serialize/Deserialize for .bsi storage
- Integrated into build pipeline as Stage 5c (between bitvectors and seal).
- CLI: --enable-sarange flag, header fields enable_sarange + sarange_variant.

Reverse FM-Index (Architecture §4.6.7 — Bidirectional Hook H1):
- src/fmindex/fmindex.{hpp,cpp}: ReverseFMIndex class.
  * Build(forward_sa, S, sentinel_row_fwd, sample_s, s_len) — derives reverse SA
    from forward SA without second SA-IS run. Constructs S_R (reverse of S), sorts
    positions by S_R suffix order, derives BWT_R, C_R, Occ_R, SA_R_samples.
  * BackwardSearch/Locate for reverse-text queries
  * LoadFromStored for deserialization
- v1.0 stores the reverse FM-index but does not use it for queries. v2.0
  bidirectional search opens the same .bsi without re-indexing (Hook H1).

Serialization:
- src/seal/seal.{hpp,cpp}: SealInput extended with isa_samples, isa_step, sarange.
  WriteBsi() now serializes ISA samples section and SARange wavelet tree payload
  to .bsi file per Architecture §7.

Tests (13/13 pass):
- tests/test_v5_features.cpp: ISA build/query/verify/info consistency (4 tests),
  SARange build/query/info/verify consistency (4 tests), SARange unit tests
  (SmallArrayRangeCount, EmptyRange, SerializeDeserialize — 3 tests),
  Full V5 build+query+reconstruct integration (2 tests).

Contract compliance:
- §3.1 C2: ISA samples fully specified — construction, storage, forward walk
- §4.3 ENHANCED tier: SARange enables O(|P| + occ_r·s + |W_r|·log(|S|/s))
  per-orientation RegionalCount bound
- §4.5 H1: Reverse FM-index stored for v2.0 forward-compatibility
- §5.2: SARange space O((|S|/s)·log(|S|/s)) bits

Ref: Architecture v4.3 §4.4, §4.6.7, §5.3, Contract v3.3 §4.3, §4.5" \
    src/sais/sais.cpp \
    src/sais/sais.hpp \
    src/fmindex/fmindex.cpp \
    src/fmindex/fmindex.hpp \
    src/sarange/sarange.cpp \
    src/sarange/sarange.hpp \
    src/seal/seal.cpp \
    src/seal/seal.hpp \
    src/cli/build.cpp \
    src/CMakeLists.txt \
    tests/test_v5_features.cpp \
    tests/CMakeLists.txt

# ══════════════════════════════════════════════════════════════════════════════
# STAGE 6 — Benchmarks: Infrastructure, scripts, synthetic results
# Per Execution Plan §5 and Contract §8
# ══════════════════════════════════════════════════════════════════════════════

do_commit "BAMSI-stage-6: Benchmark infrastructure — runner, table generator, synthetic validation

Stage 5 of Execution Plan (§5): Benchmark infrastructure and synthetic validation.

Benchmark Infrastructure:
- benchmarks/scripts/run_benchmarks.sh: Full (tool, dataset) performance matrix
  runner. Records system environment (lscpu, kernel, memory, tool versions,
  BAMSI commit hash), runs BAMSI build/count/verify/reconstruct measurements
  plus samtools/CRAM comparison, generates CSV results with configurable
  --runs N for statistical significance (mean ± stddev).
- benchmarks/scripts/generate_tables.py: Reads raw CSV results and produces
  publication-ready markdown tables for Paper 1 inclusion. Aggregates by
  (dataset, tool, operation) with mean/stddev computation.
- benchmarks/results/v1.0.0/: Results directory with environment recording,
  synthetic benchmark CSVs, and generated markdown tables.

Synthetic Validation Results (10-read):
- Build: 7ms, Compressed size: 2331 bytes
- Count (ACGT): 5ms median query latency
- Verify: 5ms, Reconstruct: 5ms

Gate Documentation:
- docs/gates/gate-5.md: Gate 5 review — V5 features + benchmark infrastructure
  sign-off. All entry/exit criteria documented.

Contract compliance:
- §8: Paper 1 benchmark battery — (tool, dataset, metric) matrix framework
- §8 sub-experiments: entropy-k sensitivity, codec ablation (via FF tests)
- Execution Plan §5.3: Environment recording, matrix runner, table generation

Pending (compute-bound, not engineering-bound):
- Full-scale benchmarks on NA12878, HG002 HiFi/ONT, exome, cohort datasets
- Reproducibility check on second machine
- Codec swap escape hatch evaluation (Execution Plan §5.3.6)

Ref: Contract v3.3 §8, Execution Plan v2.0 §5" \
    benchmarks/scripts/run_benchmarks.sh \
    benchmarks/scripts/generate_tables.py \
    benchmarks/results/ \
    docs/gates/gate-5.md \
    scripts/

# ══════════════════════════════════════════════════════════════════════════════
# STAGE 7 — Documentation: README, format spec, CLI ref, API ref, tutorials
# Per Execution Plan §6 and Contract §10.4
# ══════════════════════════════════════════════════════════════════════════════

do_commit "BAMSI-stage-7: Documentation — README, format spec, API ref, clinical guidance, tutorials

Stage 6 of Execution Plan (§6): All mandatory documentation deliverables.

Documentation Deliverables (Contract §10.4):
- README.md: Project overview, quick start, build instructions, feature list,
  architecture diagram, and contribution guidelines.
- docs/format.md: Complete .bsi binary format specification — header layout
  (42+ fields), stream sections, FM-index section, bitvector/window sections,
  directory layout, footer with global xxHash64 checksum.
- docs/algorithms.md: Algorithm descriptions — SA-IS, FM backward search,
  LF-mapping, BWT→MTF→RLE→Arithmetic codec, wavelet tree, window construction.
- docs/cli.md: CLI reference for all 10 subcommands with flags and examples.
- docs/api.md: C ABI reference — function signatures, error codes, thread safety.
- docs/clinical.md: Clinical deployment guidance — determinism guarantees,
  provenance verification, lossy-mode disclosure, no-network-at-query-time.
- docs/development.md: Developer setup and contribution guide.

Tutorials (Contract §10.4):
- docs/tutorials/01_motif_counting.md: End-to-end motif counting workflow
- docs/tutorials/02_region_query.md: Regional query with --region chr:a-b
- docs/tutorials/03_quality_postfilter.md: Quality-aware post-filtering

Architecture Decision Records:
- docs/decisions/: ADRs 0001–0008 covering language track, dependency pinning,
  target platforms, benchmark datasets, licence, benchmark hardware,
  dependency freeze policy, codec defaults.

Specification Documents:
- docs/imp/: Contract v3.3, Architecture v4.3, Execution Plan v2.0

Contract compliance:
- §10.4: All 7 mandatory documentation deliverables present
- §10.4: Three end-to-end tutorials present
- §9: Clinical guidance document covers operational guarantees

Ref: Contract v3.3 §10.4, Execution Plan v2.0 §6" \
    README.md \
    docs/format.md \
    docs/algorithms.md \
    docs/cli.md \
    docs/api.md \
    docs/clinical.md \
    docs/development.md \
    docs/tutorials/01_motif_counting.md \
    docs/tutorials/02_region_query.md \
    docs/tutorials/03_quality_postfilter.md \
    docs/imp/ \
    docs/notes/ \
    docs/runbooks/ \
    docs/PROJECT.md \
    docs/gates/gate-0.md \
    docs/gates/gate-0-download-wip.md \
    docs/decisions/0001-language-track.md \
    docs/decisions/0002-dependency-pinning.md \
    docs/decisions/0003-target-platforms.md \
    docs/decisions/0004-benchmark-datasets.md \
    docs/decisions/0005-licence-and-distribution.md \
    docs/decisions/0006-benchmark-hardware.md \
    docs/decisions/0007-dependency-freeze-policy.md \
    docs/tutorials/01-motif-counting.md \
    docs/tutorials/02-region-query.md \
    docs/tutorials/03-quality-postfilter.md

# ══════════════════════════════════════════════════════════════════════════════
# STAGE 8 — Final: Task tracker, cleanup, any remaining files
# ══════════════════════════════════════════════════════════════════════════════

do_commit "BAMSI-stage-8: Final cleanup — task tracker, remaining project files

Final commit: task tracker and any remaining project files.

- task.md: Complete execution tracker with all stages marked complete
  (Stages 0–6, V5 ENHANCED, 80+ tests, 13/13 CTest targets)
- Testing/: CMake test infrastructure artifacts
- Cleanup of any remaining untracked files

Project Status Summary:
- 23 source modules, 15 test files
- 13/13 CTest targets, 80+ individual tests, 100% pass
- Contract v3.3: All v1.0 implementable clauses satisfied
- Architecture v4.3: All modules implemented
- Execution Plan v2.0: Stages 0–6 complete, V5 features delivered" \
    task.md \
    Testing/ \
    reconstructed.bam \
    "Implementing BAMSI Project Execution Plan.md"

echo ""
echo "════════════════════════════════════════════════════════════════"
echo "All staged commits complete!"
echo "════════════════════════════════════════════════════════════════"
echo ""
echo "To verify: git log --oneline -10"
echo "To push:   git push origin HEAD"
