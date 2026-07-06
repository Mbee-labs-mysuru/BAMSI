#!/usr/bin/env bash
# scripts/run_fuzzer.sh — Run the BAMSIX format parser fuzzer
#
# Exec Plan §8.3.1: 7 days continuous fuzzing on a dedicated machine.
#
# Usage:
#   ./scripts/run_fuzzer.sh [duration_seconds]
#
# Default: 604800 seconds (7 days)
#
# Prerequisites:
#   cmake -DBUILD_FUZZ=ON -DCMAKE_CXX_COMPILER=clang++ ..
#   make fuzz_format_parser

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
BUILD_DIR="${PROJECT_DIR}/build"
CORPUS_DIR="${PROJECT_DIR}/tests/fuzz/corpus"
CRASH_DIR="${PROJECT_DIR}/tests/fuzz/crashes"
LOG_FILE="${PROJECT_DIR}/tests/fuzz/fuzzer.log"

DURATION="${1:-604800}"  # default 7 days = 604800 seconds

FUZZER="${BUILD_DIR}/tests/fuzz/fuzz_format_parser"

if [ ! -f "$FUZZER" ]; then
    echo "ERROR: Fuzzer binary not found at ${FUZZER}"
    echo "Build with: cmake -DBUILD_FUZZ=ON .. && make fuzz_format_parser"
    exit 1
fi

mkdir -p "$CORPUS_DIR" "$CRASH_DIR"

# Seed corpus: create a minimal valid .bsi header
if [ ! -f "${CORPUS_DIR}/seed_header.bsi" ]; then
    echo -n "BMSI" > "${CORPUS_DIR}/seed_header.bsi"
    printf '\x06\x00' >> "${CORPUS_DIR}/seed_header.bsi"  # version = 6
fi

# Also create an empty file seed
if [ ! -f "${CORPUS_DIR}/seed_empty.bsi" ]; then
    touch "${CORPUS_DIR}/seed_empty.bsi"
fi

echo "============================================="
echo "BAMSI Format Parser Fuzzer"
echo "============================================="
echo "Duration:     ${DURATION} seconds ($(echo "scale=1; ${DURATION}/86400" | bc) days)"
echo "Corpus:       ${CORPUS_DIR}"
echo "Crashes:      ${CRASH_DIR}"
echo "Log:          ${LOG_FILE}"
echo "Started at:   $(date -u '+%Y-%m-%dT%H:%M:%SZ')"
echo "============================================="

"$FUZZER" \
    "$CORPUS_DIR" \
    -artifact_prefix="${CRASH_DIR}/" \
    -max_len=65536 \
    -timeout=10 \
    -max_total_time="${DURATION}" \
    -print_final_stats=1 \
    2>&1 | tee "$LOG_FILE"

FUZZ_EXIT=$?

echo ""
echo "============================================="
echo "Fuzzer completed at: $(date -u '+%Y-%m-%dT%H:%M:%SZ')"
echo "Exit code: ${FUZZ_EXIT}"

# Check for crashes
CRASH_COUNT=$(find "$CRASH_DIR" -name 'crash-*' -o -name 'timeout-*' -o -name 'oom-*' 2>/dev/null | wc -l)
echo "Crashes found: ${CRASH_COUNT}"

if [ "$CRASH_COUNT" -gt 0 ]; then
    echo "FAIL: Crashes detected. See ${CRASH_DIR}/ for details."
    echo "Add crash inputs to regression tests after fixing."
    exit 1
else
    echo "PASS: No crashes in ${DURATION} seconds of fuzzing."
fi
echo "============================================="
