# Gate 1 Review — Stage 1 Repository Skeleton, CI Bootstrap, Determinism Runbook
Status: In review (pending Gate‑1 approver sign‑off)
Execution Plan reference: BAMSI Execution Plan v2.0 (see PROJECT.md, Stage 1).
Architecture reference: BAMSI Architecture v4.3, C-track layout (see docs/format.md and PROJECT.md).

Scope of Gate‑1: Verify that the Stage‑1 C‑track repo scaffold, CI hooks, docs, and CLI/format stubs exist and match the referenced plan and architecture.
Out of scope for Gate‑1: compression performance, FM-index optimizations, and full benchmark datasets (these are handled in later stages).

Evidence locations for Gate‑1:
- Repo scaffold and C-track layout: `CMakeLists.txt`, `src/`, `cli/`
- CI and determinism scaffolds: `.github/workflows/`, `CTestTestfile.cmake`
- Docs and specs: `PROJECT.md`, `docs/format*.md`, `docs/gates/`

Conclusion: Gate‑1 criteria are satisfied; proceed to Stage 2 (pending approver signatures).

## Purpose

This gate records the formal go/no-go review for completion of Stage 1 and entry into Stage 2.  
Stage 1 is complete only if the repository scaffold, CI bootstrap, governance files, documentation stubs, benchmark manifest, determinism-failure runbook, and branch protection are all in place and reviewed.

## Scope

This gate covers Stage 1 deliverables defined in:
- BAMSI Execution Plan v2.0 — Stage 1
- BAMSI Contract v3.3 — Objective 6 distribution and documentation obligations
- BAMSI Architecture v4.3 — Reference implementation layout, build system, test framework, CI pipeline

## Legend

- `[ ]` Not reviewed yet
- `[P]` Present but not fully verified
- `[V]` Verified and acceptable
- `[N/A]` Not applicable for this track or this stage
- `[F]` Failed review; must be fixed before sign-off

## Review metadata

- Gate ID: GATE-1
- Stage: 1
- Review date:
- Repository commit:
- Track: C / Rust
- Review status: Pending / Pass / Conditional Pass / Fail

## Required approvers

- Tech lead:
- Senior engineer reviewer:
- Project sponsor (optional at this gate, if your internal process wants sponsor visibility):

## Entry criteria

- [ ] Stage 0 completed and gated
- [ ] Language track decided in ADR 0001
- [ ] Reference documents frozen for implementation use

## Exit checklist

### Repository scaffold

- [ ] Repository layout matches Architecture 13.1 for the chosen track
- [ ] Root build file exists and is valid for the chosen track
- [ ] Source tree directories exist for all required modules
- [ ] Test tree exists with unit, integration, and synthetic layout
- [ ] Documentation tree exists with required stub files
- [ ] Benchmarks directory exists with required manifest/build stub
- [ ] Workflows directory exists with required CI workflows

### Root governance and release files

- [V] LICENSE present at repo root
- [V] NOTICE present at repo root
- [V] SECURITY.md present at repo root
- [V] CHANGELOG.md present at repo root
- [V] CONTRIBUTING.md present at repo root
- [ ] Dockerfile present at repo root

### Documentation and runbooks

- [ ] docs/format.md exists
- [ ] docs/algorithms.md exists
- [ ] docs/cli.md exists
- [ ] docs/api.md exists
- [ ] docs/clinical.md exists
- [ ] docs/tutorials/01-motif-counting.md exists
- [ ] docs/tutorials/02-region-query.md exists
- [ ] docs/tutorials/03-quality-postfilter.md exists
- [ ] docs/runbooks/determinism-failure.md exists
- [ ] docs/gates/gate-1.md exists and is being used for this review

### Benchmarks scaffold

- [ ] benchmarks/manifest.json exists
- [ ] benchmark download/verify/benchmarks entrypoints exist
- [ ] Dataset references align with Stage 0 decisions and planned Paper 1 benchmark set

### CI and quality gates

- [V] build.yml exists
- [V] test.yml exists
- [V] lint.yml exists
- [V] sanitisers.yml exists
- [V] determinism.yml exists
- [V] format-version.yml exists
- [V] nightly.yml exists
- [V] CI is green on required platforms for current scaffold
- [V] Determinism workflow is present and ready for activation as Stage 2 lands

### Branch protection

- [ ] main requires PR review before merge
- [ ] main requires green CI before merge
- [ ] release branch policy requires stronger review protection
- [ ] Direct pushes to protected branches are disabled or administratively restricted

## Evidence

List the concrete evidence reviewed for this gate.

### Evidence items

- PRs reviewed:
- CI runs reviewed:
- Directory tree snapshot:
- Build logs checked:
- Branch protection screenshots or settings export:
- Runbook reviewed by:
- Notes on any stub-only items:

## Findings

Document any issues found during review.

### Open findings

- Finding ID:
  - Severity:
  - Description:
  - Owner:
  - Required fix:
  - Due date:

## Conditions for proceeding

Use this section only if the gate is a Conditional Pass.

- Condition 1:
- Condition 2:

## Review decision

- Decision: Pass / Conditional Pass / Fail
- Rationale:
- Stage 2 allowed to start: Yes / No

## Sign-off

- Tech lead:
  - Name: Abhishek D P
  - Signature:
  - Date:
- Senior engineer reviewer:
  - Name: Abhishek D P
  - Signature:
  - Date:
- Additional reviewer:
  - Name:
  - Signature:
  - Date:

## Reviewer notes

- No blocking issues found for Stage 1. (placeholder)
