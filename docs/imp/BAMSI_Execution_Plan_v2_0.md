# BAMSI — Finalized Execution Plan v2.0

**Document version:** 2.0 (supersedes v1.0)
**Reference specifications:**
- BAMSI Contract v3.3 (frozen)
- BAMSI Architecture v4.3 (frozen)
- BAMSI Findings & Fixes Report
- BAMSI Objective Alignment Report

**Why v2.0 supersedes v1.0:** The v1.0 plan was a workable outline but it had gaps that would have caused real schedule slips and risk in execution. This finalized version addresses ten specific weaknesses identified during senior-engineering review. The differences from v1.0 are summarised in Appendix A.

**Audience:** Engineering leads, senior implementers, project sponsors, paper reviewers, and clinical-deployment auditors.

**Scope:** Everything between "Contract and Architecture frozen" and "v1.0.0 release tag with Paper 1 + Paper 2 in submission."

---

## How to read this document

The plan has **ten stages** with **explicit go/no-go gates between every stage** and **five vertical-slice integration milestones inside Stage 2**. Each stage documents:

- **Purpose** — what the stage delivers and why it sits where it does.
- **Entry criteria** — what must be true to start.
- **Work breakdown** — concrete tasks with traceability to Contract / Architecture clauses.
- **Exit criteria** — measurable conditions to declare done.
- **Gate review** — who signs off, what's reviewed, and what triggers a stage rollback.
- **Deliverables** — tangible artefacts.
- **Risks and mitigations** — what derails this stage in practice, mapped to fixes.
- **Time estimate** — wall-clock for a 2-person C++-track team.

The plan is paranoid by design. Every gate review is an opportunity to catch trouble early, not a bureaucratic checkpoint. A stage that fails its gate gets rolled back to the previous gate's exit state — not patched in place.

---

# Stage 0 — Spec Freeze, Decision Records, Hardware Provisioning

## 0.1 Purpose

Lock the spec, record the deliberate choices the spec leaves to the implementing project, and provision benchmark hardware **before any code is written**. The hardware decision is included here because Stage 5 benchmarks at week 14-16 cannot wait for hardware procurement at that point.

## 0.2 Entry criteria

- BAMSI Contract v3.3 finalised.
- BAMSI Architecture v4.3 finalised.
- A team of at least 2 engineers staffed and committed.
- Project sponsor identified.

## 0.3 Work breakdown

**0.3.1 Tag the spec.** Create a `bamsi-spec` repository (or directory in the main repo) and tag both documents as `v3.3` and `v4.3`. Branch protection on the spec files. Any future change requires a revision-history row, a version bump, and explicit reviewer sign-off.

**0.3.2 Write the founding ADRs.** Minimum set, all in `docs/decisions/`:

- `0001-language-track.md` — C++ vs Rust choice, with rationale from the language analysis. Signed off by tech lead.
- `0002-dependency-pinning.md` — Specific patch versions of htslib, libsais, sdsl-lite (or Rust-track succinct-library mix), zstd, xxHash, OpenSSL/libsodium.
- `0003-target-platforms.md` — Linux x86_64 + Linux aarch64 + macOS x86_64 + macOS aarch64 required; Windows best-effort.
- `0004-benchmark-datasets.md` — Accession numbers, FTP URLs, file sizes, SHA-256 checksums for NA12878 30× WGS, HG002 PacBio HiFi, HG002 ONT, exome BAM, 10-sample 1000 Genomes subset.
- `0005-licence-and-distribution.md` — Apache 2.0 confirmed; `LICENSE` and `NOTICE` content; distribution channels per Contract §10.5.
- **`0006-benchmark-hardware.md` (new in v2.0)** — Specifies the benchmark machines: one build/decompression node (high RAM, fast SSD — 256 GB RAM, 4 TB NVMe minimum recommended), one query-latency node (commodity 16-core x86_64, 64 GB RAM). Procurement or cloud-instance reservation completed in Stage 0 — not Stage 5.
- **`0007-dependency-freeze-policy.md` (new in v2.0)** — Pinned versions are immutable from week 10 (start of Stage 4) until v1.0.0 release. After freeze, dependency upgrades require security-advisory CVE justification and a written exception.

**0.3.3 Write the project charter.** `PROJECT.md` summarising scope, non-goals (cite Contract §7 directly), success criteria, the v1.0.0 definition of done, and the gate-review approver list (tech lead + project sponsor minimum).

**0.3.4 Set up the issue tracker.** Pre-populate it with the Contract / Architecture clauses, one issue per clause requiring implementation work. Tag issues with `stage-N` labels mirroring this plan.

**0.3.5 Mirror the benchmark datasets** to internal storage. Verify SHA-256 against the recorded values in ADR `0004`. This must be done in Stage 0, not Stage 5 — public FTP servers go down, accession numbers get retired.

**0.3.6 Provision benchmark hardware.** Either physical machines on order or cloud instances reserved with a clear allocation through Stage 5 (week 14-16). The reservation must include enough storage for input BAMs (~500 GB minimum) plus output `.bsi` files plus comparison-tool outputs.

**0.3.7 Engineering kick-off meeting.** Required attendance. Sign-off that the team has read and understood: Contract v3.3, Architecture v4.3, Findings & Fixes Report, Objective Alignment Report.

## 0.4 Exit criteria

- Spec documents tagged and branch-protected.
- Seven founding ADRs written and merged.
- `PROJECT.md` written and merged with named approvers.
- Issue tracker populated.
- All benchmark datasets mirrored locally and SHA-256-verified.
- Benchmark hardware provisioned (or reservation confirmed) with allocation through week 16.
- Engineering team has read all four reference documents end to end. Confirmed in writing in the kick-off meeting notes.

## 0.5 Gate 0 review — go/no-go to Stage 1

**Approvers:** Tech lead + project sponsor.

**Reviewed:**
- All ADRs are coherent (no contradictions between, e.g., language choice and dependency choice).
- Hardware allocation extends through week 16.
- Datasets verifiably present.
- Team comprehension of specs (each engineer answers three random questions about Contract clauses).

**Rollback trigger:** Any ADR is incomplete or any dataset is missing. Sponsor extends Stage 0 by up to 1 week.

## 0.6 Deliverables

- `bamsi-spec/` directory (or repo) with tagged documents.
- `docs/decisions/0001` through `0007`.
- `PROJECT.md`.
- Populated issue tracker.
- Mirrored benchmark datasets with SHA-256 verification log.
- Hardware provisioning confirmation.
- Engineering kick-off meeting notes with sign-off.

## 0.7 Time estimate

**1.5-2 weeks** for a 2-person team. The hardware-provisioning task is what makes this longer than v1.0's "1-2 weeks" estimate — procurement of physical hardware can take a week alone, and that's not negotiable.

---

# Stage 1 — Repository Skeleton, CI Bootstrap, Determinism Runbook

## 1.1 Purpose

Build the repo scaffold, CI pipeline, and the **determinism-failure runbook** that Stage 2 onwards depends on. The runbook is new in v2.0 — without it, the cross-platform determinism CI test is a tripwire without a defined response.

## 1.2 Entry criteria

- Stage 0 complete and gated through.
- Language track decided per ADR `0001`.

## 1.3 Work breakdown

**1.3.1 Create the source-tree layout** per Architecture §13.1a (C++) or §13.1b (Rust). Empty modules with the correct names.

**1.3.2 Set up the build system** per Architecture §13.2. Verify `cmake --build` (C++) or `cargo build --release --workspace` (Rust) produces expected artefacts on all four target platforms.

**1.3.3 Create the documentation tree** per Architecture §13.1. Stub files for `format.md`, `algorithms.md`, `cli.md`, `api.md`, `clinical.md`, plus three tutorial stubs.

**1.3.4 Create the test layout.** `tests/unit/`, `tests/integration/` (TIER 1), `tests/synthetic/` (TIER 2). Stub a passing "hello world" test in each.

**1.3.5 Create `benchmarks/` with the manifest.** `benchmarks/manifest.json` lists each dataset with name, source URL, expected SHA-256, expected file size. `benchmarks/Makefile` has `make download`, `make verify`, `make benchmarks` (stub).

**1.3.6 Configure CI workflows:**

- `build.yml` — matrix over Linux x86_64, Linux aarch64, macOS x86_64, macOS aarch64.
- `test.yml` — unit + TIER 1 integration tests on every push.
- `lint.yml` — language-appropriate static analysis.
- `sanitisers.yml` — ASan, UBSan (C++) or Miri (Rust) on PRs to `release/*`.
- `determinism.yml` — cross-platform byte-identity test (initially mostly empty; activates as Stage 2 modules land).
- `format-version.yml` — verifies current binary can read previous-major `.bsi` files (scaffolding for v2.0).
- `nightly.yml` — TIER 2 tests + benchmark smoke run.

**1.3.7 Add `LICENSE`, `NOTICE`, `SECURITY.md`, `CHANGELOG.md`, `CONTRIBUTING.md`, `Dockerfile`** at the repo root.

**1.3.8 Set up branch protection.** `main` requires PR review + green CI before merge. `release/*` requires double review.

**1.3.9 Write the determinism-failure runbook** (new in v2.0) at `docs/runbooks/determinism-failure.md`. Defines the response protocol when `determinism.yml` fails:

1. Failed PR cannot merge.
2. The author bisects the failing commit using a synthetic 10K-read input.
3. Common causes (in order of empirical frequency): parallel-iterator collect order in Rayon/parallel STL; hash-map iteration order; floating-point arithmetic; timestamps in headers; uninitialised memory; environment-dependent default thread pool size.
4. Fix: replace non-deterministic operation with a stable-sorted equivalent, or pin the seed/order, or zero the timestamp field.
5. If determinism is broken on `main` (not a PR), revert immediately and reproduce the failure under controlled conditions before re-applying.
6. If three consecutive determinism breaks happen in one week, escalate to tech lead — the underlying design has a non-determinism source that needs structural fix, not patch fixes.

**1.3.10 Set up the gate-review tracking.** A `docs/gates/` directory with one file per gate, recording the review date, approvers, exit-criteria checklist, sign-off signatures, and any conditions for proceeding.

## 1.4 Exit criteria

- All seven CI workflows green on the empty repo.
- Build successful on all four target platforms.
- Documentation tree exists with stubs.
- Benchmark manifest present (datasets mirrored from Stage 0).
- Branch protection enforced.
- Determinism-failure runbook merged.
- Gate-review tracking infrastructure in place.

## 1.5 Gate 1 review — go/no-go to Stage 2

**Approvers:** Tech lead + at least one senior engineer who will work in Stage 2.

**Reviewed:**
- CI is genuinely green on all four platforms (not green on x86_64 with aarch64 timing out — the tracker confirms green status).
- Determinism runbook is comprehensible to someone who didn't write it.
- Module layout actually matches Architecture §13.1 — verify directory listing.

**Rollback trigger:** Any platform fails to build, any CI workflow is missing, branch protection is not enforced.

## 1.6 Deliverables

- The complete repo scaffold matching Architecture §13.1.
- Seven CI workflows.
- `benchmarks/manifest.json`.
- All licence/governance files.
- `docs/runbooks/determinism-failure.md`.
- `docs/gates/gate-1.md` signed.

## 1.7 Time estimate

**1.5-2 weeks** for a 2-person team. CI configuration eats more time than expected; the runbook adds 0.5 days but pays for itself the first time determinism fails.

---

# Stage 2 — Build Pipeline with Five Vertical-Slice Milestones

## 2.1 Purpose

Implement the build path of the indexer end-to-end. **The major v2.0 change is structural:** instead of "implement Stages 1-10 and validate per module," this stage now proceeds in **five vertical-slice milestones** where the entire pipeline (or its current scope) runs end-to-end on a tiny input before any milestone is declared done. This catches inter-module integration bugs at their cheapest fix point — the moment two modules first touch.

## 2.2 Entry criteria

- Stage 1 complete and gated through.
- CI is green and stable on all four platforms.
- Engineering team has read Contract §0-§3 and Architecture §1-§4 end-to-end.

## 2.3 Work breakdown — five vertical-slice milestones

### Milestone V1 — Trivial round-trip (week 1, ~5 days)

**Goal:** A 10-read synthetic BAM goes through the full pipeline and produces a `.bsi` that reconstructs to the same 10 reads.

**Stages implemented (minimum viable):**
- **Stage 1 — ingestion:** htslib BAM open, primary-read filter, RawRead extraction. Stub source manifest hash to a fixed value.
- **Stage 2 — ordering:** sort by `(chrom_id, pos, source_file_id, bam_offset)`, assign `read_id`.
- **Stage 3 — sequence build:** concatenate to S with `#` separators.
- **Stage 4 — SA-IS:** call libsais (or equivalent), produce SA + BWT + sentinel_row.
- **Stage 5a — S_seq:** stub as plain bzip2 of BWT (don't implement MTF+RLE+arithmetic yet).
- **Stage 5b — FM-index:** wavelet tree over BWT, C, Occ, SA samples. Use sdsl-lite for everything.
- **Stage 6 — S_qual / S_meta / S_map:** stub as plain ZSTD per stream.
- **Stage 7 — windows:** single window covering all of S.
- **Stage 8 — bitvectors:** B_read only (B_window trivial).
- **Stage 10 — sealing:** write `.bsi` with checksums.

**Gate test:** synthetic 10-read BAM → build → reconstruct → assert byte-identity vs original BAM (on the in-scope ℛ subset).

**Why this matters:** at the end of week 1 of Stage 2, the team has a working end-to-end pipeline. Every subsequent milestone replaces stubs with the real spec'd implementation. This is the single biggest change from v1.0.

### Milestone V2 — FM correctness on synthetic input (week 2, ~5 days)

**Goal:** GlobalCount and Locate work correctly against brute-force on a synthetic 10K-read input.

**Stages implemented or hardened:**
- Real source manifest hash and ordering hash.
- Real `M_ℓ` mapping function with CIGAR support — start with M-only (no soft clips, no insertions).
- Real B_read with rank/select.
- A working `bamsi count` and `bamsi locate` CLI that calls into the library.

**Gate test:** synthetic 10K-read BAM with random reads, 100 random patterns of length 8-30. For each pattern, BAMSI's GlobalCount must match brute-force linear scan over the reads.

**Risk caught here:** rank-API confusion (Invariant I15). If GlobalCount is correct but Locate returns wrong `read_id`, the closed-vs-half-open rank bug is in play.

### Milestone V3 — Codec bake-off and lock-in (week 3, ~5 days) — NEW in v2.0

**Goal:** Pick the production codec choices for S_qual / S_meta / S_map by measuring compression ratio on a single representative real dataset (not a synthetic one).

**Why this milestone exists in v2.0:** the v1.0 plan deferred codec selection to Stage 5 benchmarks at week 14-16. If defaults underperformed Genozip, the team had no time to swap and re-validate. This milestone forces the choice early so the rest of Stage 2 hardens the chosen codecs, not stubs.

**Procedure:**

1. Pick the smallest representative real BAM — recommend the exome BAM from the dataset manifest, ~10 GB compressed.
2. Implement all alternates from the catalogue: for S_qual implement RANGE_CYCLE, RANS_DELTA, and ZSTD_DICT; for S_meta implement TYPED_SPLIT and ZSTD_FALLBACK; for S_map implement DELTA_RANGE and RAW.
3. Build the `.bsi` with each combination on the test dataset.
4. Measure compression ratio and encoder/decoder time per stream.
5. **Lock in the production defaults** based on the bake-off results, recorded in ADR `0008-codec-defaults.md`. This is a binding decision, not a hint.

**Gate test:** ADR `0008` merged with the chosen defaults and the bake-off table. Subsequent milestones use only the locked-in codec choices for production builds.

### Milestone V4 — Real-world BAM end-to-end (week 4-5, ~10 days)

**Goal:** A real-world 30× WGS BAM (NA12878) builds, queries, and reconstructs correctly. CIGAR edge cases are handled.

**Stages implemented or hardened:**
- Full CIGAR mapping with soft-clip, insertion, deletion, ref-skip handling per Architecture §5.2 (the `10S150M5S` worked example is a test case).
- Real `bwt_rank_open(...)` and `bv_rank_closed(...)` as separately named functions; the C++/Rust type system or named-function-only enforcement prevents aliasing.
- Multi-window construction (Stage 7) with `T = 100,000` S-characters.
- B_window with rank/select.
- `RegionalCount` and `RegionalExists` in BASE tier.
- Full per-stream codec implementations using V3 lock-ins.
- TIER 1 integration tests for every Invariant I1-I15 that's reachable at this scope.

**Gate test:** NA12878 30× WGS BAM. Build completes; cross-platform determinism CI passes; `bamsi reconstruct` produces a BAM byte-identical to the original on ℛ.

**Risk caught here:** memory pressure at scale (the 24 GB directory overhead becomes real), parallel-iteration non-determinism, real-world CIGAR weirdness.

### Milestone V5 — ENHANCED tier and bidirectional FM (week 6-7, ~10 days)

**Goal:** Enable SARange (ENHANCED tier RegionalCount) and `enable_bidirectional` (Hook H1 for v2.0 approximate matching). Verify the build configuration toggles work.

**Stages implemented:**
- SARange — wavelet tree over `SA_samples` with `range_count(lo, hi, l_j, r_j)`.
- Reverse FM-index derived via `SA_R[i] = |S| − 1 − SA[|S| − 1 − i]`.
- ISA samples optional path.
- The complete `bamsi build/count/exists/locate/region-count/region-exists/reconstruct/info/verify` CLI.
- Stable C ABI.

**Gate test:** all five queries pass TIER 1 tests on NA12878 with `enable_sarange = true` and `enable_bidirectional = true`. The ENHANCED-tier RegionalCount returns the same counts as BASE-tier on 50 random regional queries.

## 2.4 Exit criteria for Stage 2 as a whole

- All five vertical-slice milestones complete with their gate tests passing.
- TIER 1 tests cover every reachable Invariant I1-I15.
- TIER 2 synthetic tests exist and pass nightly.
- Cross-platform determinism CI has been green for at least 7 consecutive commits.
- All nine CLI subcommands implemented with `--json` support.
- C ABI test program passes under ASan.
- Codec defaults locked per ADR `0008`.

## 2.5 Gate 2 review — go/no-go to Stage 3

**Approvers:** Tech lead + project sponsor.

**Reviewed:**
- All five milestone gate tests demonstrably passed.
- ADR `0008` codec lock-in present and reasoned.
- Determinism CI green for 7+ consecutive commits.
- TIER 1 coverage of Invariants I1-I15 verified line-by-line against Architecture §9.3.

**Rollback trigger:** any milestone gate test failing, codec defaults not locked, or determinism CI red. Stage 2 extends, Stage 3 does not start.

## 2.6 Deliverables

- All `src/` modules implementing the pipeline.
- Per-module unit tests.
- TIER 1 integration tests.
- TIER 2 synthetic tests.
- Working `bamsi` CLI with all nine subcommands.
- C ABI test program.
- ADR `0008-codec-defaults.md`.
- Five gate-test reports in `docs/gates/`.

## 2.7 Risks and mitigations

- **Risk:** Stage 5b FM-index is the most algorithmically dense module; teams stall here. **Mitigation:** Pair-program through this stage; reference sdsl-lite tutorial code line by line. The V1 milestone uses sdsl-lite throughout, deferring custom code to later.
- **Risk:** Two-rank-API confusion (I15). **Mitigation:** V2 gate test catches read_id mismatches; V4 gate test re-verifies after CIGAR is added.
- **Risk:** Determinism violations from parallel ingestion or hash-map iteration. **Mitigation:** Cross-platform CI test catches it; runbook from Stage 1 defines the response.
- **Risk:** Codec defaults locked at V3 turn out to underperform on long-read data discovered at Stage 5. **Mitigation:** V3 measures on the exome BAM specifically because it's representative without being the worst case; if Stage 5 reveals a significant gap on long-read, the alternate codec is already implemented and a Stage 5.5 codec swap is a 3-day task, not a re-design.
- **Risk:** Implementation drift from spec. **Mitigation:** every PR cites a Contract or Architecture clause; reviewers refuse PRs without citations.

## 2.8 Time estimate

**6-7 weeks** for a 2-person C++ team. Breaks down as: V1 (1 week) + V2 (1 week) + V3 (1 week, codec bake-off) + V4 (2 weeks, real BAM) + V5 (1.5 weeks, ENHANCED + bidirectional). This is slightly longer than v1.0's "6-8 weeks" estimate but the milestone structure means the team finds out at week 1 if things are going wrong, not at week 8.

---

# Stage 3 — Query Path Hardening and CLI Polish

## 3.1 Purpose

Stage 2's V5 milestone delivered a working query layer. Stage 3 polishes it: edge cases, error handling, output formatting, the complete C ABI, and hardening that wasn't necessary at the milestone gate.

This stage is **shorter than v1.0's "3-4 weeks"** because Stage 2 now delivers a more complete query layer at V5 than v1.0's Stage 2 ever did.

## 3.2 Entry criteria

- Stage 2 gated through.
- All five queries pass TIER 1 on NA12878.

## 3.3 Work breakdown

### 3.3.1 Streaming-mode and sorted-mode polish

- Streaming mode is the default for `Locate`. Verify time-to-first-result on a 10-million-occurrence query is under 100 ms.
- Sorted mode is opt-in via `--sort-output`. Verify byte-identity across two runs (deterministic order).
- Document the modes in `docs/cli.md`.

### 3.3.2 `bamsi info` and `bamsi verify` complete

- `bamsi info --json` emits all fields specified in Contract §9.2 — 30+ fields. Verified against a JSON schema in `tests/fixtures/cli_schemas/`.
- `bamsi verify` re-computes `source_manifest_hash`, `ordering_hash`, all section checksums.
- Tampering test: flip a byte, run `bamsi verify`, expect non-zero exit with `CHECKSUM_MISMATCH`.

### 3.3.3 `bamsi reconstruct` for partial reconstruction

- Subset reconstruction: `--streams seq,map` produces a partial output with only sequences and coordinates.
- Per-read reconstruction: `--read-ids 0,1,2,...` produces output for specific reads.
- Lossy-mode handling: refuses without `--allow-lossy` if `is_lossless = 0` in the header.

### 3.3.4 Approximate-query stubs

- `bamsi_approx_locate_hamming` and `bamsi_approx_locate_edit` exposed in `libbamsi.h`.
- Both return `NOT_IMPLEMENTED_V1`.
- The C ABI test program calls them, asserts the error code.

### 3.3.5 Error handling sweep

- Every error code in Architecture §3 ErrorCode enum has at least one test that triggers it.
- Error messages are useful: include the failing field name, the expected value, the actual value.
- No silent fallbacks anywhere in the query path.

### 3.3.6 C ABI hardening

- The C ABI test program in `tests/c_abi/` tests every exported function with valid and invalid inputs.
- All allocations explicitly freed via `bamsi_*_free` to avoid cross-runtime issues.
- The C ABI test program passes under ASan and UBSan.

## 3.4 Exit criteria

- All nine CLI subcommands have `--json` schema validation.
- `bamsi verify` correctly detects every kind of tamper.
- Approximate-query stubs return `NOT_IMPLEMENTED_V1` cleanly.
- C ABI test program passes under ASan.
- Streaming and sorted modes both work with documented latency.
- Every ErrorCode triggered by at least one test.

## 3.5 Gate 3 review — go/no-go to Stage 4

**Approvers:** Tech lead.

**Reviewed:**
- CLI feature-completeness against Contract §10.2 line-by-line.
- C ABI completeness against the function set defined in Contract §10.3.
- Error-handling sweep complete.

**Rollback trigger:** missing CLI subcommand, C ABI gap, silent-fallback found in any code path.

## 3.6 Deliverables

- Hardened query layer.
- Complete CLI with `--json` schemas.
- Complete C ABI.
- C ABI test program.
- `docs/gates/gate-3.md` signed.

## 3.7 Time estimate

**2-3 weeks** for a 2-person team (down from v1.0's 3-4 because Stage 2's milestone structure delivered more in V5).

---

# Stage 4 — Comprehensive Validation Campaign

## 4.1 Purpose

Stage 4 closes every remaining coverage gap before benchmarking begins. The Findings & Fixes Report identified 22 issues; Stage 4 confirms every one is regression-tested. Architecture §9.3's full validation matrix is exercised. **This stage starts the dependency freeze** per ADR `0007`.

## 4.2 Entry criteria

- Stage 3 gated through.
- All five queries hardened.

## 4.3 Work breakdown

### 4.3.1 Dependency freeze

ADR `0007` activates: pinned dependency versions are immutable from this point until v1.0.0 release. Exception requires CVE justification and written approval from tech lead.

### 4.3.2 Map every Invariant I1-I15 to a test

Confirm a test exists, runs in the correct tier, and passes. Architecture §9.3 has the table; Stage 4 walks it line by line and opens issues for any uncovered invariant.

### 4.3.3 Map every Findings & Fixes finding to a regression test

For each F1-F6, S1-S9, C1-C3 finding, confirm a regression test that would catch a regression to the pre-fix behaviour:

- **F1** (space bound): test that builds NA12878 and asserts `|.bsi|` is within the documented `O(|S|)` bound.
- **F2** (SARange tier): tier-comparison test (50 random queries, BASE vs ENHANCED match).
- **F3** (no `S[pos]==#` access): code-grep test in CI confirming no occurrence of `S[` syntax in `src/query/`; RSS-introspection test in TIER 2 verifying S not loaded.
- **F4** (directory split): build with `seq_block_size=1024` and `=0`; reconstruct same content; assert byte-identity.
- **F5** (pipeline topology): assert SA computed exactly once per build via instrumentation.
- **F6** (window unit-consistency): test with different coverage levels validating the `|W_r| = O(L·d/T)` formula.
- **S1** (two-rank APIs): function-uniqueness test.
- **S2** (entropy-k): not a regression test; verify the sensitivity-analysis script exists for Paper 1.
- **S3** (parallel SA-IS): bit-identity test (only if `allow_parallel_sa = true`).
- **S4** (validation tiers): structural test that TIER 1 is ≤5% of build time on the exome BAM.
- **S5** (soft-clip example): the `10S150M5S` worked example as a verbatim test case.
- **S6** (counting semantics): strand-complete test on patterns where P ≠ rc(P).
- **S7** (Tables A and B): partial reconstruction test for each subset.
- **S8** (chrom_id portability): test that two `.bsi` from BAMs with different `@SQ` orderings return same query results via string chrom names.
- **S9** (quality post-filter): tutorial 03 runs end-to-end as a smoke test.
- **C1** (output ordering modes): streaming and sorted return same multiset.
- **C2** (ISA samples): `locate_bidir` returns same positions as `locate`.
- **C3** (motif-scope calibration): test exists; Stage 5 runs it on real data.

### 4.3.4 Cross-platform determinism — sustained green

The CI workflow has been green for 7+ commits at end of Stage 2. Stage 4 monitors it for the entire stage. If it ever turns red:

1. P0 incident.
2. Bisect using the runbook from Stage 1.
3. Fix and reset the consecutive-green counter.

### 4.3.5 Codec correctness ablations

Every codec variant from V3 has a round-trip test on a sampled real BAM. The bake-off table from ADR `0008` is verified to still hold; if results changed, ADR is updated.

### 4.3.6 Lossy-mode obligations

- Build with `qual_lossy_bins = 8` succeeds.
- `is_lossless = 0` in header.
- `bamsi info` surfaces lossy condition.
- `bamsi reconstruct` warns or refuses without `--allow-lossy`.
- Quality reconstruction matches binned form, not byte-identical to original.

### 4.3.7 Provenance and audit

- `bamsi verify` detects single-byte tamper, manifest mismatch, version mismatch.
- No-implicit-network test: run `bamsi count` under network isolation; succeeds.

### 4.3.8 Audit document

`docs/audit.md` maps every Contract / Architecture clause to its implementing test. Useful for paper reviewers and clinical-deployment auditors. **This is a deliverable**, not optional.

## 4.4 Exit criteria

- Every Invariant I1-I15 has a passing test in the correct tier.
- Every Findings & Fixes finding has a regression test.
- `determinism.yml` green throughout the stage.
- Codec ablation table verified.
- Lossy-mode obligations verified.
- `bamsi verify` correctly detects tampering, manifest mismatch, version mismatch.
- `docs/audit.md` complete and merged.
- Dependency freeze in effect; no upgrades since week 10.

## 4.5 Gate 4 review — go/no-go to Stage 5

**Approvers:** Tech lead + project sponsor.

**Reviewed:**
- Audit document complete and accurate (sponsor checks 5 random clauses).
- Dependency freeze respected (no `Cargo.lock` or `vcpkg.json` changes since week 10 except CVE-justified).
- All invariant + finding tests passing.

**Rollback trigger:** any uncovered invariant, any unfreezable dependency change without CVE, any audit-document gap.

## 4.6 Deliverables

- Complete validation matrix in `tests/`.
- `docs/audit.md`.
- `docs/gates/gate-4.md` signed.

## 4.7 Time estimate

**2-3 weeks** for a 2-person team. Most of the work is filling gaps that surface during the audit-document writing. The audit document itself is the forcing function: writing it makes you confront every clause.

---

# Stage 5 — Benchmark Execution

## 5.1 Purpose

Run the full Paper 1 benchmark battery (Contract §8) against the v1.0 implementation. **The major v2.0 change is splitting compute time from engineering time** so the schedule is honest about what's blocking what.

## 5.2 Entry criteria

- Stage 4 gated through.
- All five datasets present and SHA-256-verified.
- Benchmark hardware accessible.

## 5.3 Work breakdown

### 5.3.1 Set up the benchmark environment (engineering, ~3 days)

- Provision dedicated benchmark machines (already procured in Stage 0).
- Pin OS, kernel, CPU governor (performance), drop-caches policy.
- Record `lscpu`, `lsblk`, kernel version, comparison-tool versions.
- Record BAMSI commit hash and Docker image SHA-256.

### 5.3.2 Run the (tool, dataset) matrix (compute, ~7-10 days)

For each of the 5 datasets and each of the 4 comparison tools (CRAM, Genozip, samtools view, NGC), plus BAMSI:

- Compression-ratio measurements: build time + RSS at 1, 8, 16 threads; compressed size in bits per base pair (sequence/quality/total); decompression time + RSS.
- Query-latency measurements (BAMSI vs `samtools view`): 100+ patterns of length 8-32, against five region scales.
- Three independent runs per (tool, dataset) for std-dev.

This is mostly compute time. The team is unblocked on Stage 6 documentation while it runs.

### 5.3.3 Sub-experiments (compute + analysis, ~5 days)

Per Contract §8:
- Entropy-k sensitivity (3 datasets × 5 k values).
- Motif-scope calibration (10 motifs vs `samtools view -c`).
- Soft-clip interval characterisation.
- Directory granularity trade-off.
- Parallel SA-IS speedup (if enabled).
- Window-pruning benefit.
- Stream-codec ablation (refresh from V3, on full datasets now).
- Bidirectional-FM cost (if enabled).

### 5.3.4 Generate publication-ready tables (engineering, ~2 days)

`benchmarks/scripts/generate_tables.py` reads raw CSVs, produces LaTeX or markdown tables for direct Paper 1 inclusion.

### 5.3.5 Reproducibility check (compute + verification, ~2 days)

Independent re-run on a different machine using only the Docker image, dataset manifest, and one `make benchmarks` command. Verify numbers reproduce within recorded std-dev.

### 5.3.6 Codec swap escape hatch (engineering, contingency)

If results show a chosen codec underperforming meaningfully (e.g., S_qual default >10% worse than ZSTD_DICT on long-read), swap the default in ADR `0008`, re-run only the affected sub-experiment (1-2 days), update the audit document. **This is the contingency plan that v1.0 didn't have.**

## 5.4 Exit criteria

- Full (tool, dataset, metric) matrix populated.
- All sub-experiments complete.
- Reproducibility check passes on a different machine.
- LaTeX/markdown tables in `benchmarks/results/v1.0.0/`.
- BAMSI competitive (within 1.10×) compression vs Genozip on NA12878 — informational target.
- Median GlobalCount latency under 50 ms — informational target.

## 5.5 Gate 5 review — go/no-go to Stage 7 (Paper 1)

**Approvers:** Tech lead + project sponsor.

**Reviewed:**
- Numbers genuinely reproduce on a second machine.
- Targets met or explicitly disclosed.
- Codec choices in ADR `0008` still defensible against the data.

**Rollback trigger:** numbers don't reproduce; or significant target miss without disclosed mitigation.

## 5.6 Deliverables

- Populated `benchmarks/results/`.
- Reproducibility-verified Docker image with SHA-256.
- Publication-ready tables.
- `docs/gates/gate-5.md` signed.

## 5.7 Time estimate

**Engineering: ~1.5 weeks. Compute: ~2-3 weeks (overlaps engineering).** The team is unblocked on Stages 6 and 7-Paper-2 during compute time. Total wall-clock: 2-3 weeks.

---

# Stage 6 — Documentation Completion

## 6.1 Purpose

Complete the documentation deliverables so v1.0.0 ships with credible, accurate, maintained docs. **This stage runs in parallel with Stages 4-5.** The stub files were created in Stage 1; they are filled in alongside the modules they describe.

## 6.2 Entry criteria

- Stage 3 complete (CLI and library APIs stable).
- Documentation can be written from current memory, not stale memory.

## 6.3 Work breakdown

### 6.3.1 README.md

One-paragraph summary, quick-start (Docker → build → count, 5 minutes from copy to result), pointers to deeper docs, citation block.

### 6.3.2 docs/format.md

Byte-level `.bsi` format spec mirroring Architecture §7. Section-by-section walkthrough with diagrams. Suitable for an independent re-implementer.

### 6.3.3 docs/algorithms.md

High-level algorithmic overview citing Contract §3-§4, Architecture §4-§6. Pipeline diagram. Glossary.

### 6.3.4 docs/cli.md

Full reference for every subcommand. Recipe section for common workflows.

### 6.3.5 docs/api.md

Stable C ABI reference (the ground truth) plus the chosen-track native-API reference (C++ namespace from `bamsi.hpp` or Rust crate rustdoc).

### 6.3.6 docs/clinical.md

Operational guidance: provenance verification, audit-log integration, lossy-mode caveats, no-implicit-network policy.

### 6.3.7 Three end-to-end tutorials

Each tutorial runs as a smoke test in nightly CI to prevent drift.

- Tutorial 01: motif counting on NA12878.
- Tutorial 02: region query on `chr17:7570000-7580000` (TP53).
- Tutorial 03: two-phase quality filtering.

### 6.3.8 ADRs continued

Continue `docs/decisions/` for any non-trivial implementation choice in Stages 2-5.

## 6.4 Exit criteria

- All seven major docs and three tutorials written and reviewed.
- Tutorials run end-to-end by someone other than the author.
- No "TODO" or "TBD" left.
- Tutorials are smoke tests in nightly CI.

## 6.5 Gate 6 review — go/no-go to Stage 8

**Approvers:** Tech lead.

**Reviewed:**
- All docs are current (last-modified within Stage 6 timeframe).
- Tutorials genuinely run on a fresh machine.

**Rollback trigger:** any major doc with stale information or any tutorial that fails on a fresh machine.

## 6.6 Deliverables

- All `docs/*.md` and `docs/tutorials/*.md`.
- ADR series complete for v1.0.0.
- `docs/gates/gate-6.md` signed.

## 6.7 Time estimate

**1.5-2 weeks** of focused writing, parallel with Stages 4-5. If the docs were filled in alongside their modules (as Stage 1's stubs encourage), this is mostly review and polish. If they were deferred, add 2 more weeks.

---

# Stage 7 — Paper Drafting and Submission

## 7.1 Purpose

Convert the spec, the implementation, and the benchmark numbers into two Q1-quality papers. **The major v2.0 change is splitting Paper 1 and Paper 2 timelines explicitly:** Paper 2 (Theory) starts at week 10 alongside Stage 4 because it doesn't need numbers; Paper 1 (Systems) starts at week 17 after Stage 5 numbers are in.

## 7.2 Entry criteria — Paper 2

- Stage 4 in progress (theorems are implementation-verified).

## 7.3 Entry criteria — Paper 1

- Stage 5 complete with reproducibility-verified numbers.
- Stage 6 complete (Paper 1 methods cite the format spec).

## 7.4 Work breakdown — Paper 2 (Theory) — starts week 10

### 7.4.1 Structure

Following Contract §3-§6 and §8 Paper 2 list:

- Abstract.
- Preliminaries (alphabet, sequence, read collection, reverse-complement, strand-complete counting).
- Index structures (FM-index, bitvectors, windows, sampling, SARange).
- Theorems with proofs (Soundness, Completeness, Rank-convention lemma, BWT lifecycle lemma, Space bound, Tier-aware complexity).
- Approximate-search forward-compatibility discussion (hooks H1-H4).
- Conclusion.

### 7.4.2 Internal review (3 rounds minimum)

Formal proofs benefit from outside eyes. At minimum:
- Round 1: tech lead + one other senior engineer.
- Round 2: external reviewer familiar with succinct data structures.
- Round 3: full team.

### 7.4.3 Submit

Submit to RECOMB, Bioinformatics (Theory & Algorithms track), JCB, or WABI. Preprint on arXiv simultaneously.

## 7.5 Work breakdown — Paper 1 (Systems) — starts week 17

### 7.5.1 Structure

Following Contract §8 Paper 1 list:

- Abstract.
- Introduction with prior art (CRAM, Genozip, NGC, samtools).
- System overview citing Contract §0.1, §2, §4 and Architecture §1.1.
- Stream-specific codecs (Contract §2.7-§2.9 with empirical justification from Stage 5).
- Query architecture (FM-index + bitvectors + windows + SARange; tier-aware complexity).
- Implementation (language track, dependencies, cross-platform determinism).
- Evaluation (the (tool, dataset, metric) matrix from Stage 5; sub-experiments; ablations).
- Discussion (scope, forward compatibility, clinical workflow).
- Related work, conclusion, future work.

### 7.5.2 Reproducibility statement

Docker image SHA-256, `make benchmarks` command, full benchmark scripts in supplementary materials.

### 7.5.3 Internal review (2-3 rounds)

### 7.5.4 Submit

Submit to SIGMOD, VLDB, Bioinformatics (Applications Track), or USENIX ATC. Preprint on bioRxiv simultaneously.

## 7.6 Joint considerations

- Both papers cite each other.
- Both papers cite the Contract and Architecture documents (these can be supplementary or arXiv-hosted technical reports).
- Author list, affiliations, acknowledgments aligned.

## 7.7 Exit criteria

- Paper 2 drafted, reviewed 3+ rounds, submitted, preprinted.
- Paper 1 drafted, reviewed 2-3+ rounds, submitted, preprinted.
- Reproducibility supplementary materials prepared.

## 7.8 Gate 7 review — go/no-go to Stage 8

**Approvers:** Tech lead + project sponsor + paper co-authors.

**Reviewed:**
- Both papers are submission-ready (not "almost done").
- Preprints posted.

**Rollback trigger:** missing review round, claims in paper not backed by the Contract / Architecture / Stage 5 numbers.

## 7.9 Deliverables

- Paper 1 + Paper 2 manuscripts.
- bioRxiv + arXiv preprints.
- `docs/gates/gate-7.md` signed.

## 7.10 Time estimate

- **Paper 2: 4-6 weeks** wall-clock starting week 10. Most of it parallel with Stages 4-5.
- **Paper 1: 4-5 weeks** starting week 17. Bulk of the work happens after Stage 5 numbers are final.

---

# Stage 8 — Pre-Release Hardening

## 8.1 Purpose

Final pass before tagging v1.0.0. Security review, sanitisers, cross-platform CI, lossy-mode end-to-end, audit-trail end-to-end. **The v2.0 change is increasing fuzzing budget from 24 hours to 7 days continuous** — for a clinical-deployment tool, 24 hours is barely a smoke test.

## 8.2 Entry criteria

- Stages 2-7 complete.
- Stage 5 numbers indicate the system is shippable.

## 8.3 Work breakdown

### 8.3.1 Security review

- Threat model: attacker controls a `.bsi` file. What damage can they do?
- Run AFL or libFuzzer (C++) / `cargo-fuzz` (Rust) against the format parser for **7 days continuous** on a dedicated machine.
- Fix every crash; add crashing inputs to regression-test corpus.
- Review every `unsafe` block (Rust) or `reinterpret_cast` (C++).
- Review network-related code paths (there should be none, but verify).

### 8.3.2 Sanitiser sweep

Full TIER 1 + TIER 2 suite under ASan + UBSan (C++) or Miri (Rust). TIER 1 under ThreadSanitizer.

### 8.3.3 Cross-platform CI sustained green

All four target platforms green on full test suite for at least 14 consecutive commits (raised from v1.0's "7"). The release criterion is stricter than the Stage 2 criterion.

### 8.3.4 Lossy-mode end-to-end

End-to-end test of every lossy-mode obligation from Stage 4 §4.3.6, on full datasets now.

### 8.3.5 Audit-trail end-to-end

- Build, record SHA-256.
- `bamsi info --json` parses, all 30+ fields present.
- Re-build on a different machine, verify byte-identical.
- Tamper, verify `bamsi verify` fails non-zero with correct error code.

### 8.3.6 Performance regression sweep

Compare Stage 5 numbers against a baseline machine. If any metric has regressed since the most recent benchmark run, bisect.

## 8.4 Exit criteria

- All security findings fixed; no fuzzer crashes in 7 days continuous.
- Sanitisers clean.
- All four platforms green on full suite for 14+ consecutive commits.
- Lossy-mode tests passing on full datasets.
- Audit-trail test passing.
- No performance regression vs Stage 5 baseline.

## 8.5 Gate 8 review — go/no-go to Stage 9

**Approvers:** Tech lead + project sponsor + security reviewer.

**Reviewed:**
- Fuzzer log shows 7 days clean.
- Audit trail demonstrably complete.
- Performance numbers match Stage 5.

**Rollback trigger:** any new crash, any silent regression, any audit-trail gap.

## 8.6 Deliverables

- Hardened codebase ready for v1.0.0 tag.
- Fuzzer corpus added to regression tests.
- Security review notes in `docs/decisions/`.
- `docs/gates/gate-8.md` signed.

## 8.7 Time estimate

**2-3 weeks** for a 2-person team. The 7-day fuzzing run is wall-clock that the team is unblocked on (final paper polish, release-note writing).

---

# Stage 9 — v1.0.0 Release

## 9.1 Purpose

Tag, build, sign, and publish v1.0.0 across all distribution channels.

## 9.2 Entry criteria

- Stage 8 gated through.
- All papers submitted (or preprinted).

## 9.3 Work breakdown

### 9.3.1 Tag and build

- Bump version to 1.0.0 in build files, CHANGELOG, README.
- Tag `v1.0.0` after the version-bump PR merges.
- Build static-linked binaries for all four platforms.
- Build Docker image; record SHA-256.
- Build source tarball; sign with GPG.

### 9.3.2 Publish

- GitHub: v1.0.0 release with binaries, signed tarball, signed checksums, release notes citing both preprints.
- Docker Hub / GHCR: push `bamsi/bamsi:v1.0.0` and `:latest`. Record SHA-256 in release notes.
- Bioconda: submit recipe to `bioconda-recipes`. **Plan for 1-2 weeks of asynchronous review (see Stage 10).**
- crates.io (Rust track only): `cargo publish` for the workspace crates.

### 9.3.3 Announce

- Post preprint links to community channels.
- Email collaborating labs.
- Update project website.

### 9.3.4 Set up the v1.x maintenance line

- Create `release/1.x` branch.
- Document patch-release process in `CONTRIBUTING.md`.
- Triage incoming issues into v1.0.x patches and v2.0 features.

## 9.4 Exit criteria

- v1.0.0 tag exists.
- Binaries published to GitHub.
- Docker image published.
- Bioconda submission acknowledged (acceptance is asynchronous).
- crates.io publish succeeded (if Rust).
- Preprints publicly visible.
- `release/1.x` branch exists.

## 9.5 Gate 9 review — go/no-go to Stage 10

**Approvers:** Project sponsor.

**Reviewed:**
- Release artefacts genuinely available and downloadable.
- Preprints accessible.
- Maintenance branch ready.

**Rollback trigger:** broken artefact, missing preprint, no maintenance branch.

## 9.6 Deliverables

- v1.0.0 release on GitHub.
- Public Docker image.
- Public crate or Bioconda submission.
- Public preprints.
- `docs/gates/gate-9.md` signed.

## 9.7 Time estimate

**1 week** of engineering. Bioconda review is asynchronous; not on critical path.

---

# Stage 10 — Post-Release Maintenance, Bioconda Iteration, v2.0 Ramp

## 10.1 Purpose

Sustain v1.0 in production, complete the Bioconda review iteration (new in v2.0 — v1.0 plan treated this as fire-and-forget), handle paper review cycles, and ramp v2.0 work.

## 10.2 Work breakdown

### 10.2.1 Bioconda iteration (1-3 weeks asynchronous)

Bioconda reviewers ask for changes. Expected items:
- Recipe metadata adjustments.
- Test-during-build additions.
- Dependency-version constraints.

Each round: respond within 3 days, push fix, request re-review. Typical: 1-3 rounds before acceptance. **Plan for this — don't be surprised by it.**

### 10.2.2 Patch releases (v1.0.x)

Issue triage. Patch cadence: monthly first quarter, slowing to as-needed.

Each patch: changelog entry, SemVer bump, regression tests for the fix, CI green on all platforms.

### 10.2.3 Minor releases (v1.1, v1.2, …)

Backward-compatible additions: new codec choices, new CLI options, performance optimisations.

ABI-compatible only — existing behaviour unchanged.

### 10.2.4 Paper review cycles

Address Paper 1 and Paper 2 reviewer comments. Re-run benchmarks if methodology changes are requested.

### 10.2.5 v2.0 work

Implement approximate-matching overlay using H1-H4 hooks. Multi-sample joint indexing. Quality-aware queries. Secondary/supplementary alignment indexing.

v2.0 has its own Contract version (4.0+) and Architecture version (5.0+). v1.0 `.bsi` files remain readable by v2.0 by design.

### 10.2.6 Community

Public issue tracker. Office-hours or Discourse. Conference talks. Contributor onboarding.

## 10.3 Time estimate

**Ongoing.** Each patch release is 1-3 days. Minor releases: 1-3 weeks. v2.0 effort: multi-month, starting roughly 6 months after v1.0.0 ships.

---

# Total wall-clock to v1.0.0

| Configuration | Wall-clock | Notes |
|---|---|---|
| 2-person C++ team (recommended) | **18-22 weeks (~5 months)** | Honest range with no slack |
| 2-person Rust team | 22-26 weeks (~6 months) | Add ~30% for ecosystem gaps |
| Solo C++ | 28-32 weeks (~7-8 months) | No parallelism on milestones |
| Solo Rust | 32-36 weeks (~8-9 months) | Both penalties combined |

The v1.0 plan estimated 16-20 weeks. The v2.0 estimate is 18-22 weeks because it accounts for:
- Vertical-slice milestone gate testing in Stage 2 (+1 week vs v1.0).
- Codec bake-off in Stage 2 V3 (+1 week vs v1.0).
- Stage-gate review time across all gates (+0.5 weeks total).
- Hardware provisioning in Stage 0 (no schedule impact, but explicitly budgeted).
- 7-day fuzzing in Stage 8 (no engineering wall-clock impact, but explicit).
- Bioconda iteration in Stage 10 (post-release, not on critical path).

The v1.0 estimate assumed ideal execution; the v2.0 estimate assumes realistic execution with the gate-review structure that catches problems early.

---

# Appendix A — Differences from Execution Plan v1.0

| Change | Why |
|---|---|
| Stage 2 restructured into 5 vertical-slice milestones (V1-V5) | v1.0 had 6-8 weeks of parallel module work without integration milestones; integration bugs surfaced at week 8 instead of week 1. |
| Stage 2 V3 = early codec bake-off + ADR `0008` codec lock-in | v1.0 deferred codec selection to Stage 5 benchmarks; if defaults underperformed, the team had no time to swap. V3 forces the choice early. |
| Stage 0 includes ADR `0006` (benchmark hardware) and ADR `0007` (dependency freeze policy) | v1.0 deferred hardware to Stage 5 and didn't address dependency drift; both caused real schedule risk. |
| Stage 1 adds the determinism-failure runbook | v1.0 said "treat as P0" without defining the response protocol. |
| Stage 2 dependency freeze activates at week 10 | v1.0 had no freeze date; dependencies could change mid-project. |
| Stage 4 produces `docs/audit.md` mapping every clause to its test | v1.0 had per-module tests but no central audit document; useful for paper reviewers and clinical auditors. |
| Stage 5 separates compute time from engineering time | v1.0 estimated 2-3 weeks without distinguishing the two. |
| Stage 5 has explicit codec-swap escape hatch | If V3 codec choices underperform on full benchmarks, the swap is a 3-day task, not a re-design. |
| Stage 7 splits Paper 1 and Paper 2 timelines explicitly | Paper 2 starts at week 10 alongside Stage 4; Paper 1 starts at week 17 after Stage 5 numbers. |
| Every stage has a named go/no-go gate review with approvers | v1.0 listed exit criteria but didn't require sign-off; without that, scope creeps. |
| Stage 8 fuzzing budget = 7 days continuous (was 24 hours) | For a clinical-deployment tool, 24 hours of fuzzing is barely a smoke test. |
| Stage 10 has explicit Bioconda iteration sub-stage | v1.0 treated Bioconda submission as fire-and-forget. |
| Total wall-clock estimate revised from 16-20 weeks to 18-22 weeks | The v1.0 estimate assumed ideal execution; v2.0 estimates honestly. |

---

# Appendix B — Discipline that holds the plan together

These rules apply across every stage. They're not stage-specific tasks — they're working norms.

- **Every PR cites a Contract or Architecture clause.** No clause citation = the PR is either v2.0 scope (refuse) or a spec gap (open a doc-revision PR first).
- **Cross-platform determinism CI is P0.** The runbook from Stage 1 defines the response. Any failure on `main` triggers immediate revert and reproduction.
- **Stage gates are blocking.** A stage that fails its gate review rolls back to the previous gate's exit state. Gates are not bureaucracy; they are the only way scope creep is caught early enough to fix.
- **TIER 1 tests live with their modules.** They are not deferred to Stage 4. Stage 4 is the gap-closing audit, not the writing campaign.
- **Documentation lives with its module.** Stage 6 is the closing review, not the writing campaign.
- **Codec defaults locked at V3 are binding.** A late-Stage-5 swap is a contingency, not a default.
- **Determinism trumps performance.** Any optimisation that breaks I9 reverts. Performance optimisations come after determinism is preserved.
- **No silent fallback.** Anywhere. Every error path returns a structured error code with a useful message.
- **Dependency freeze is real.** From week 10, dependency upgrades require CVE justification.
- **Reproducibility is a release criterion, not a nice-to-have.** Stage 5's reproducibility check on a different machine is a gate, not a checkbox.

---

*End of BAMSI Finalized Execution Plan v2.0*
*Reference: Contract v3.3, Architecture v4.3, Findings & Fixes Report, Objective Alignment Report*
*Supersedes: Execution Plan v1.0*
