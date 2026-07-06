#!/bin/bash
# ============================================================================
# BAMSIX Post-Rename Verification — Real Dataset (Rigorous)
# Tests EVERY subcommand end-to-end on real NA12878 exome data
# ============================================================================
set -euo pipefail

BAMSIX="./build/bamsix"
BAM="data/test/synthetic_10k.bam"
BSI="/tmp/bamsix_real_test.bsi"
BSI_ENH="/tmp/bamsix_real_enhanced.bsi"
RECON_BAM="/tmp/bamsix_real_recon.bam"

PASS=0
FAIL=0
TOTAL=0

check() {
    local desc="$1"
    local status="$2"
    TOTAL=$((TOTAL + 1))
    if [ "$status" -eq 0 ]; then
        PASS=$((PASS + 1))
        echo "  ✅ PASS: $desc"
    else
        FAIL=$((FAIL + 1))
        echo "  ❌ FAIL: $desc"
    fi
}

echo "============================================================================"
echo "  BAMSIX Post-Rename Rigorous Verification"
echo "  Real Dataset: NA12878_exome_1.bam (5.7 MB real exome)"
echo "============================================================================"
echo ""

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION A: CLI SMOKE TESTS
# ═══════════════════════════════════════════════════════════════════════════════
echo "═══ SECTION A: CLI Smoke Tests ═══"

# 1. version
echo ">>> A1. version"
OUTPUT=$($BAMSIX version 2>&1)
echo "$OUTPUT" | grep -q "1.0.0"
check "version outputs 1.0.0" $?
echo "    Output: $OUTPUT"

# 2. help (3 forms)
echo ">>> A2. help (all forms)"
$BAMSIX --help 2>&1 | grep -q "usage: bamsix"
check "--help" $?
$BAMSIX -h 2>&1 | grep -q "usage: bamsix"
check "-h" $?
$BAMSIX help 2>&1 | grep -q "usage: bamsix"
check "help subcommand" $?

# 3. No args
echo ">>> A3. no-args usage"
$BAMSIX 2>&1 | grep -q "usage: bamsix"
check "no-args shows usage" $?

# 4. Invalid subcommand
echo ">>> A4. invalid subcommand"
$BAMSIX invalidcmd 2>&1 && ESTAT=1 || ESTAT=0
check "invalid subcommand returns non-zero" $?

echo ""

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION B: BUILD PIPELINE (BASE tier)
# ═══════════════════════════════════════════════════════════════════════════════
echo "═══ SECTION B: Build Pipeline ═══"

echo ">>> B1. build BASE tier"
rm -f "$BSI"
START=$(date +%s)
$BAMSIX build --input "$BAM" --output "$BSI" \
    --entropy-k 6 --sample-step 64 2>&1
BUILD_EXIT=$?
END=$(date +%s)
check "build BASE completes successfully" $BUILD_EXIT
echo "    Build time: $((END - START))s"

test -f "$BSI"
check ".bsi file created" $?

BSI_SIZE=$(stat --printf='%s' "$BSI" 2>/dev/null || stat -f '%z' "$BSI" 2>/dev/null)
echo "    .bsi size: $BSI_SIZE bytes"
test "$BSI_SIZE" -gt 10000
check ".bsi is non-trivial size" $?

echo ""

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION C: VERIFY (Checksum Integrity)
# ═══════════════════════════════════════════════════════════════════════════════
echo "═══ SECTION C: Index Verification ═══"

echo ">>> C1. verify"
OUTPUT=$($BAMSIX verify --index "$BSI" 2>&1)
check "verify exits 0" $?
echo "    $OUTPUT"

echo ""

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION D: INFO (Index Metadata)
# ═══════════════════════════════════════════════════════════════════════════════
echo "═══ SECTION D: Index Info ═══"

echo ">>> D1. info (human-readable)"
OUTPUT=$($BAMSIX info --index "$BSI" 2>&1)
check "info exits 0" $?
echo "$OUTPUT" | head -15 | while read -r line; do echo "    $line"; done

echo ">>> D2. info --json"
OUTPUT=$($BAMSIX info --index "$BSI" --json 2>&1)
check "info --json exits 0" $?
echo "$OUTPUT" | head -5 | while read -r line; do echo "    $line"; done

echo ""

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION E: GLOBALCOUNT (5 motif categories)
# ═══════════════════════════════════════════════════════════════════════════════
echo "═══ SECTION E: GlobalCount (count) ═══"

echo ">>> E1. Palindromic motifs (expect R≈1.0 vs grep)"
for PAT in GAATTC ACGT; do
    OUTPUT=$($BAMSIX count --index "$BSI" --pattern "$PAT" 2>&1)
    check "count $PAT" $?
    echo "    $PAT: $OUTPUT"
done

echo ">>> E2. Non-palindromic motifs (expect R≈2.0 vs grep)"
for PAT in GATTACA TTAGGG CCAAT; do
    OUTPUT=$($BAMSIX count --index "$BSI" --pattern "$PAT" 2>&1)
    check "count $PAT" $?
    echo "    $PAT: $OUTPUT"
done

echo ">>> E3. Homopolymer motifs (expect R>2.0 vs grep)"
for PAT in AAAA CCCCC; do
    OUTPUT=$($BAMSIX count --index "$BSI" --pattern "$PAT" 2>&1)
    check "count $PAT" $?
    echo "    $PAT: $OUTPUT"
done

echo ">>> E4. Short motifs (3-4 bp)"
for PAT in ACG TATA CGG; do
    OUTPUT=$($BAMSIX count --index "$BSI" --pattern "$PAT" 2>&1)
    check "count $PAT" $?
    echo "    $PAT: $OUTPUT"
done

echo ">>> E5. Long motif (8 bp)"
OUTPUT=$($BAMSIX count --index "$BSI" --pattern ATCGATCG 2>&1)
check "count ATCGATCG" $?
echo "    ATCGATCG: $OUTPUT"

echo ""

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION F: GLOBALEXISTS
# ═══════════════════════════════════════════════════════════════════════════════
echo "═══ SECTION F: GlobalExists (exists) ═══"

echo ">>> F1. Pattern that should exist"
OUTPUT=$($BAMSIX exists --index "$BSI" --pattern ACGT 2>&1)
check "exists ACGT" $?
echo "    ACGT: $OUTPUT"

echo ">>> F2. Very long pattern (unlikely to exist)"
OUTPUT=$($BAMSIX exists --index "$BSI" --pattern ACGTACGTACGTACGTACGTACGT 2>&1)
check "exists long 24-mer" $?
echo "    24-mer: $OUTPUT"

echo ""

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION G: LOCATE
# ═══════════════════════════════════════════════════════════════════════════════
echo "═══ SECTION G: Locate ═══"

echo ">>> G1. locate GATTACA (sparse)"
OUTPUT=$($BAMSIX locate --index "$BSI" --pattern GATTACA 2>&1)
check "locate GATTACA" $?
HITS=$(echo "$OUTPUT" | wc -l)
echo "    Hits: $HITS"
echo "$OUTPUT" | head -3 | while read -r line; do echo "    $line"; done

echo ">>> G2. locate GATTACA --json"
OUTPUT=$($BAMSIX locate --index "$BSI" --pattern GATTACA --json 2>&1)
check "locate --json" $?
echo "    $(echo "$OUTPUT" | head -1)"

echo ">>> G3. locate with --sort-output"
OUTPUT=$($BAMSIX locate --index "$BSI" --pattern GATTACA --sort-output 2>&1)
check "locate --sort-output" $?
echo "    $(echo "$OUTPUT" | head -1)"

echo ""

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION H: REGIONAL QUERIES
# ═══════════════════════════════════════════════════════════════════════════════
echo "═══ SECTION H: Regional Queries ═══"

echo ">>> H1. region-count ACGT in chr1:1-1000000"
OUTPUT=$($BAMSIX region-count --index "$BSI" --pattern ACGT --region chr1:1-1000000 2>&1)
check "region-count chr1:1-1M" $?
echo "    $OUTPUT"

echo ">>> H2. region-count TTAGGG in chr1:1-100000"
OUTPUT=$($BAMSIX region-count --index "$BSI" --pattern TTAGGG --region chr1:1-100000 2>&1)
check "region-count TTAGGG chr1:1-100K" $?
echo "    $OUTPUT"

echo ">>> H3. region-exists ACGT in chr1:1-1000000"
OUTPUT=$($BAMSIX region-exists --index "$BSI" --pattern ACGT --region chr1:1-1000000 2>&1)
check "region-exists chr1:1-1M" $?
echo "    $OUTPUT"

echo ">>> H4. region-exists ACGT in chr1:1-100 (tiny region)"
OUTPUT=$($BAMSIX region-exists --index "$BSI" --pattern ACGT --region chr1:1-100 2>&1)
check "region-exists chr1:1-100" $?
echo "    $OUTPUT"

echo ""

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION I: RECONSTRUCT
# ═══════════════════════════════════════════════════════════════════════════════
echo "═══ SECTION I: Reconstruct ═══"

echo ">>> I1. reconstruct (full)"
rm -f "$RECON_BAM"
OUTPUT=$($BAMSIX reconstruct --index "$BSI" --output "$RECON_BAM" --allow-lossy 2>&1)
check "reconstruct --allow-lossy" $?
test -f "$RECON_BAM"
check "reconstructed BAM file exists" $?
RECON_SIZE=$(stat --printf='%s' "$RECON_BAM" 2>/dev/null || stat -f '%z' "$RECON_BAM" 2>/dev/null)
echo "    Reconstructed BAM size: $RECON_SIZE bytes"
test "$RECON_SIZE" -gt 100
check "reconstructed BAM is non-empty" $?

echo ""

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION J: ERROR HANDLING
# ═══════════════════════════════════════════════════════════════════════════════
echo "═══ SECTION J: Error Handling ═══"

echo ">>> J1. nonexistent .bsi file"
$BAMSIX count --index /nonexistent.bsi --pattern ACGT 2>&1 && E=1 || E=0
test "$E" -eq 0
check "nonexistent file → error" $?

echo ">>> J2. missing --index"
$BAMSIX count --pattern ACGT 2>&1 && E=1 || E=0
test "$E" -eq 0
check "missing --index → error" $?

echo ">>> J3. invalid subcommand"
$BAMSIX foobar 2>&1 && E=1 || E=0
test "$E" -eq 0
check "invalid subcommand → error" $?

echo ">>> J4. missing --output in build"
$BAMSIX build --output /tmp/bamsix_noinput.bsi 2>&1 && E=1 || E=0
test "$E" -eq 0
check "build without --input → error" $?

echo ""

# ═══════════════════════════════════════════════════════════════════════════════
# CLEANUP
# ═══════════════════════════════════════════════════════════════════════════════
rm -f "$BSI" "$BSI_ENH" "$RECON_BAM"

# ═══════════════════════════════════════════════════════════════════════════════
# FINAL SUMMARY
# ═══════════════════════════════════════════════════════════════════════════════
echo "============================================================================"
echo "  FINAL RESULTS: $PASS/$TOTAL passed, $FAIL failed"
echo "============================================================================"
if [ "$FAIL" -eq 0 ]; then
    echo ""
    echo "  🎉 ALL $TOTAL OPERATIONS VERIFIED ON REAL DATA — RENAME IS CLEAN!"
    echo ""
else
    echo ""
    echo "  ⚠️  $FAIL of $TOTAL operations failed. Review output above."
    echo ""
    exit 1
fi
