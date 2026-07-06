#!/usr/bin/env bash
# ═══════════════════════════════════════════════════════════════════════
# BAMSIX v1.0 — Staged GitHub Push Script (v2)
#
# This script creates branches from master and pulls the LATEST files
# from the source branch (story/BAMSI-02) where all development was done.
#
# Usage:
#   bash scripts/github_push_v1.sh <section>
#   Sections: cleanup | build | engine | tests | docs | workflows
# ═══════════════════════════════════════════════════════════════════════

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO_ROOT"

REMOTE="origin"
BASE_BRANCH="master"
SOURCE_BRANCH="story/BAMSI-02"   # Branch with ALL the latest v1.0 code

# Color output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

info()  { echo -e "${BLUE}[INFO]${NC}  $1"; }
ok()    { echo -e "${GREEN}[OK]${NC}    $1"; }
warn()  { echo -e "${YELLOW}[WARN]${NC}  $1"; }
err()   { echo -e "${RED}[ERROR]${NC} $1"; exit 1; }

show_staged_summary() {
    echo ""
    info "═══ Staged files summary ═══"
    git diff --cached --stat
    echo ""
    local count
    count=$(git diff --cached --name-only | wc -l)
    info "Total files staged: $count"
    echo ""
}

# Helper: pull specific files from source branch into current branch
pull_from_source() {
    for f in "$@"; do
        if git show "${SOURCE_BRANCH}:${f}" &>/dev/null; then
            git checkout "$SOURCE_BRANCH" -- "$f"
        else
            warn "File not in $SOURCE_BRANCH: $f (skipping)"
        fi
    done
}

# Helper: pull an entire directory from source branch
pull_dir_from_source() {
    for d in "$@"; do
        git checkout "$SOURCE_BRANCH" -- "$d" 2>/dev/null || warn "Dir not in $SOURCE_BRANCH: $d"
    done
}

# ─────────────────────────────────────────────────────────────────────
# SECTION 0: defect/BAMSIX-01 — Clean up existing repo
# ─────────────────────────────────────────────────────────────────────
do_cleanup() {
    local BRANCH="defect/BAMSIX-01"
    info "═══ $BRANCH: Clean up existing GitHub repo ═══"

    git fetch "$REMOTE"
    git checkout "$BASE_BRANCH"
    git pull "$REMOTE" "$BASE_BRANCH"
    git checkout -b "$BRANCH"

    # Delete files/dirs that should NOT be on GitHub
    info "Removing business-sensitive files..."

    git rm -rf spec/ 2>/dev/null && ok "Removed spec/" || warn "spec/ not tracked"
    git rm -rf images/ 2>/dev/null && ok "Removed images/" || warn "images/ not tracked"
    git rm -f BLANK_README.md 2>/dev/null && ok "Removed BLANK_README.md" || warn "BLANK_README.md not tracked"
    git rm -rf docs/papers/ 2>/dev/null || true
    git rm -rf docs/imp/ 2>/dev/null || true
    git rm -rf docs/gates/ 2>/dev/null || true
    git rm -rf docs/notes/ 2>/dev/null || true
    git rm -rf docs/runbooks/ 2>/dev/null || true
    git rm -f docs/PROJECT.md 2>/dev/null || true
    git rm -f docs/paper2_theory_draft.md 2>/dev/null || true
    git rm -rf benchmarks/ 2>/dev/null || true

    # Pull the updated .gitignore from source branch
    pull_from_source .gitignore
    git add .gitignore
    ok "Updated .gitignore with comprehensive exclusions"

    show_staged_summary

    git commit -m "defect/BAMSIX-01: Remove business-sensitive files from public repo

- Remove spec/ (architecture document, contract)
- Remove images/ (loose screenshots)
- Remove BLANK_README.md (template placeholder)
- Remove docs/papers/, docs/imp/, docs/gates/, docs/notes/, docs/runbooks/
- Remove benchmarks/ (experiment infrastructure)
- Update .gitignore to permanently block sensitive files

These files contain unpublished research, business execution plans,
and proprietary architecture decisions not suitable for public release."

    info "Pushing $BRANCH..."
    git push "$REMOTE" "$BRANCH"

    ok "═══ $BRANCH pushed! Create PR on GitHub, get it reviewed and merged. ═══"
    echo ""
    warn "WAIT: Merge this PR into $BASE_BRANCH before running the next section."
}

# ─────────────────────────────────────────────────────────────────────
# SECTION 1: story/BAMSIX-01 — Build foundation
# ─────────────────────────────────────────────────────────────────────
do_build() {
    local BRANCH="story/BAMSIX-01"
    info "═══ $BRANCH: Build foundation — CMake, headers, deps, Dockerfile ═══"

    git fetch "$REMOTE"
    git checkout "$BASE_BRANCH"
    git pull "$REMOTE" "$BASE_BRANCH"
    git checkout -b "$BRANCH"

    # Pull latest files from source branch
    pull_from_source CMakeLists.txt Dockerfile .gitmodules LICENSE LICENSE.txt NOTICE
    ok "Pulled build system files from $SOURCE_BRANCH"

    pull_dir_from_source include/
    ok "Pulled include/bamsix/ (public headers)"

    pull_dir_from_source external/
    ok "Pulled external/ (dependencies)"

    git add CMakeLists.txt Dockerfile .gitmodules LICENSE LICENSE.txt NOTICE
    git add include/
    git add external/

    show_staged_summary

    git commit -m "story/BAMSIX-01: Build foundation — CMake, headers, dependencies

Build System:
- CMakeLists.txt (C++20, Format v6, project version 1.0.0)
- Dockerfile (containerized build)
- .gitmodules (5 external dependencies)

Public Headers (include/bamsix/):
- bamsix.h: Stable C ABI (18 functions, 13 status codes)
- bamsix.hpp: C++ public API
- types.hpp: Core data structures
- status.hpp, version.hpp, config.hpp.in

External Dependencies:
- htslib 1.21, libsais 2.8.6, sdsl-lite 2.1.1
- zstd 1.5.6, xxHash 0.8.3

License: Apache 2.0"

    info "Pushing $BRANCH..."
    git push "$REMOTE" "$BRANCH"

    ok "═══ $BRANCH pushed! Create PR on GitHub. ═══"
    warn "WAIT: Merge this PR into $BASE_BRANCH before running the next section."
}

# ─────────────────────────────────────────────────────────────────────
# SECTION 2: story/BAMSIX-02 — Core engine
# ─────────────────────────────────────────────────────────────────────
do_engine() {
    local BRANCH="story/BAMSIX-02-engine"
    info "═══ $BRANCH: Core engine — all 19 source modules ═══"

    git fetch "$REMOTE"
    git checkout "$BASE_BRANCH"
    git pull "$REMOTE" "$BASE_BRANCH"
    git checkout -b "$BRANCH"

    pull_dir_from_source src/
    ok "Pulled src/ (19 modules) from $SOURCE_BRANCH"

    git add src/

    show_staged_summary

    git commit -m "story/BAMSIX-02: Core engine — 19 source modules

Build Pipeline (10 stages):
- ingest, ordering, seqbuilder, sais, fmindex
- seqencode, streamencode, windows, bitvectors, validation

Query Engine:
- query, mapping, sarange (ENHANCED tier)

I/O and Integration:
- seal, format, reconstruct, cabi, cli

Reserved for v2.0:
- fmindexreverse (bidirectional FM-index placeholder)"

    info "Pushing $BRANCH..."
    git push "$REMOTE" "$BRANCH"

    ok "═══ $BRANCH pushed! Create PR on GitHub. ═══"
    warn "WAIT: Merge this PR into $BASE_BRANCH before running the next section."
}

# ─────────────────────────────────────────────────────────────────────
# SECTION 3: story/BAMSIX-03 — Test suite + CI
# ─────────────────────────────────────────────────────────────────────
do_tests() {
    local BRANCH="story/BAMSIX-03"
    info "═══ $BRANCH: Test suite, CI, synthetic data, verification scripts ═══"

    git fetch "$REMOTE"
    git checkout "$BASE_BRANCH"
    git pull "$REMOTE" "$BASE_BRANCH"
    git checkout -b "$BRANCH"

    # Tests
    pull_dir_from_source tests/
    git add tests/
    ok "Pulled tests/ (24 test files)"

    # Synthetic test data
    mkdir -p data/test
    pull_from_source data/test/synthetic_10k.bam
    git add data/test/synthetic_10k.bam
    ok "Pulled data/test/synthetic_10k.bam (549 KB)"

    # Safe verification scripts only
    for script in verify_bamsix_operations.sh run_fuzzer.sh verify_bruteforce.sh \
                  verify_datasets.sh install_samtools.sh; do
        pull_from_source "scripts/$script" 2>/dev/null || true
    done
    git add scripts/verify_bamsix_operations.sh scripts/run_fuzzer.sh \
            scripts/verify_bruteforce.sh scripts/verify_datasets.sh \
            scripts/install_samtools.sh 2>/dev/null || true
    ok "Pulled verification scripts"

    # CI workflows
    pull_dir_from_source .github/workflows/
    git add .github/workflows/
    ok "Pulled .github/workflows/ (7 CI pipelines)"

    show_staged_summary

    git commit -m "story/BAMSIX-03: Test suite, CI pipelines, synthetic data

Test Suite (24 tests):
- Tier 1/2 invariants, FM correctness, roundtrip
- Codec completeness, determinism, C ABI, lossy E2E
- Error sweep, audit verification, locate sorted
- Block directory, benchmark validation, tutorial smoke
- V5 features, findings regression, FM verification

Synthetic Data:
- data/test/synthetic_10k.bam (549 KB, 10k reads)

Verification Scripts:
- verify_bamsix_operations.sh (39 end-to-end checks)
- run_fuzzer.sh, verify_bruteforce.sh, verify_datasets.sh

CI (.github/workflows/):
- build, test, lint, determinism, format-version
- sanitisers, nightly"

    info "Pushing $BRANCH..."
    git push "$REMOTE" "$BRANCH"

    ok "═══ $BRANCH pushed! Create PR on GitHub. ═══"
    warn "WAIT: Merge this PR into $BASE_BRANCH before running the next section."
}

# ─────────────────────────────────────────────────────────────────────
# SECTION 4: story/BAMSIX-04 — Documentation
# ─────────────────────────────────────────────────────────────────────
do_docs() {
    local BRANCH="story/BAMSIX-04"
    info "═══ $BRANCH: Documentation ═══"

    git fetch "$REMOTE"
    git checkout "$BASE_BRANCH"
    git pull "$REMOTE" "$BASE_BRANCH"
    git checkout -b "$BRANCH"

    # Root docs
    pull_from_source README.md CHANGELOG.md CONTRIBUTING.md SECURITY.md
    git add README.md CHANGELOG.md CONTRIBUTING.md SECURITY.md
    ok "Pulled root documentation files"

    # Technical docs (public only)
    for doc in format.md algorithms.md cli.md api.md clinical.md audit.md \
               development.md dependency_freeze.md; do
        pull_from_source "docs/$doc" 2>/dev/null || true
    done
    git add docs/format.md docs/algorithms.md docs/cli.md docs/api.md \
            docs/clinical.md docs/audit.md docs/development.md \
            docs/dependency_freeze.md 2>/dev/null || true
    ok "Pulled technical docs"

    # ADRs
    pull_dir_from_source docs/decisions/
    git add docs/decisions/
    ok "Pulled docs/decisions/ (8 ADR files)"

    # Tutorials
    pull_dir_from_source docs/tutorials/
    git add docs/tutorials/
    ok "Pulled docs/tutorials/ (3 tutorials)"

    # README images
    pull_dir_from_source docs/images/
    git add docs/images/
    ok "Pulled docs/images/ (9 figures)"

    show_staged_summary

    git commit -m "story/BAMSIX-04: Documentation — README, docs, tutorials

Root Documentation:
- README.md (2100+ lines, 9 embedded figures)
- CHANGELOG.md, CONTRIBUTING.md, SECURITY.md

Technical Docs:
- format.md, algorithms.md, cli.md, api.md
- clinical.md, audit.md, development.md, dependency_freeze.md

ADRs: 0001-0008
Tutorials: motif counting, region query, quality postfilter
Figures: 9 PNG/SVG from IEEE TCBB paper"

    info "Pushing $BRANCH..."
    git push "$REMOTE" "$BRANCH"

    ok "═══ $BRANCH pushed! Create PR on GitHub. ═══"
    warn "WAIT: Merge this PR into $BASE_BRANCH before running the next section."
}

# ─────────────────────────────────────────────────────────────────────
# SECTION 5: story/BAMSIX-05 — Snakemake workflow
# ─────────────────────────────────────────────────────────────────────
do_workflows() {
    local BRANCH="story/BAMSIX-05"
    info "═══ $BRANCH: Snakemake workflow ═══"

    git fetch "$REMOTE"
    git checkout "$BASE_BRANCH"
    git pull "$REMOTE" "$BASE_BRANCH"
    git checkout -b "$BRANCH"

    pull_from_source workflows/Snakefile workflows/config.yaml workflows/config_synthetic.yaml
    git add workflows/Snakefile workflows/config.yaml workflows/config_synthetic.yaml
    ok "Pulled workflows/ (Snakefile + 2 configs)"

    show_staged_summary

    git commit -m "story/BAMSIX-05: Snakemake reference workflow

- Snakefile: 6-rule pipeline (build, verify, info, motif count,
  regional count, samtools comparison)
- config.yaml: Default config for real datasets
- config_synthetic.yaml: Config for synthetic_10k.bam

Usage:
  pip install snakemake
  cd workflows
  snakemake --cores 4 --configfile config_synthetic.yaml"

    info "Pushing $BRANCH..."
    git push "$REMOTE" "$BRANCH"

    ok "═══ $BRANCH pushed! All done! ═══"
    ok "🎉 ALL 6 BRANCHES PUSHED! Review and merge PRs in order."
}

# ─────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────
case "${1:-}" in
    cleanup)    do_cleanup ;;
    build)      do_build ;;
    engine)     do_engine ;;
    tests)      do_tests ;;
    docs)       do_docs ;;
    workflows)  do_workflows ;;
    *)
        echo ""
        echo "BAMSIX v1.0 — Staged GitHub Push Script"
        echo "========================================"
        echo ""
        echo "PREREQUISITE: Commit all changes on $SOURCE_BRANCH first:"
        echo "  git add -A && git commit -m 'checkpoint: BAMSIX v1.0'"
        echo ""
        echo "Usage: bash scripts/github_push_v1.sh <section>"
        echo ""
        echo "  1. cleanup    → defect/BAMSIX-01"
        echo "  2. build      → story/BAMSIX-01"
        echo "  3. engine     → story/BAMSIX-02-engine"
        echo "  4. tests      → story/BAMSIX-03"
        echo "  5. docs       → story/BAMSIX-04"
        echo "  6. workflows  → story/BAMSIX-05"
        echo ""
        echo "After each: Create PR → Review → Merge → Run next"
        echo ""
        ;;
esac
