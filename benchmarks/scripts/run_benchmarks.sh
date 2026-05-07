#!/usr/bin/env bash
# BAMSI Benchmark Runner — Execution Plan §5.3
# Runs the full (tool, dataset) performance matrix.
#
# Usage: ./benchmarks/scripts/run_benchmarks.sh [--datasets DIR] [--threads N] [--runs N]
#
# Prerequisites:
#   - BAMSI binary built: ./build/bamsi
#   - Comparison tools installed: samtools, genozip (optional)
#   - Datasets in benchmarks/datasets/ (or specified via --datasets)
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
BAMSI="$PROJECT_DIR/build/bamsi"

# Defaults
DATASET_DIR="${PROJECT_DIR}/benchmarks/datasets"
THREADS=1
RUNS=3
OUTPUT_DIR="${PROJECT_DIR}/benchmarks/results/v1.0.0"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)

# Parse args
while [[ $# -gt 0 ]]; do
    case $1 in
        --datasets) DATASET_DIR="$2"; shift 2 ;;
        --threads)  THREADS="$2"; shift 2 ;;
        --runs)     RUNS="$2"; shift 2 ;;
        --output)   OUTPUT_DIR="$2"; shift 2 ;;
        *) echo "Unknown arg: $1"; exit 1 ;;
    esac
done

mkdir -p "$OUTPUT_DIR"

# ─── System Environment Recording (§5.3.1) ───────────────────────────────────

ENV_FILE="$OUTPUT_DIR/environment_${TIMESTAMP}.txt"
{
    echo "=== BAMSI Benchmark Environment ==="
    echo "Date: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "BAMSI version: $($BAMSI version 2>&1 || echo 'unknown')"
    echo "BAMSI commit: $(cd "$PROJECT_DIR" && git rev-parse HEAD 2>/dev/null || echo 'unknown')"
    echo ""
    echo "=== CPU ==="
    lscpu 2>/dev/null || echo "lscpu not available"
    echo ""
    echo "=== Memory ==="
    free -h 2>/dev/null || echo "free not available"
    echo ""
    echo "=== Kernel ==="
    uname -a
    echo ""
    echo "=== Disk ==="
    lsblk 2>/dev/null || echo "lsblk not available"
    echo ""
    echo "=== Comparison Tools ==="
    echo "samtools: $(samtools --version 2>&1 | head -1 || echo 'not installed')"
    echo "genozip:  $(genozip --version 2>&1 | head -1 || echo 'not installed')"
} > "$ENV_FILE"
echo "[benchmark] Environment recorded: $ENV_FILE"

# ─── Utility functions ────────────────────────────────────────────────────────

measure_time_rss() {
    # Returns: wall_time_s, max_rss_kb
    local cmd="$1"
    /usr/bin/time -v bash -c "$cmd" 2>&1 | grep -E "wall clock|Maximum resident" | \
        awk '/wall clock/ {print $NF} /Maximum resident/ {print $NF}'
}

run_bamsi_benchmark() {
    local bam="$1"
    local name="$2"
    local csv="$OUTPUT_DIR/${name}_bamsi_${TIMESTAMP}.csv"

    echo "dataset,tool,operation,run,wall_time_s,file_size_bytes,rss_kb" > "$csv"

    local bsi="/tmp/bamsi_bench_${name}.bsi"
    local recon="/tmp/bamsi_bench_${name}_recon.bam"

    for run in $(seq 1 $RUNS); do
        echo "[benchmark] $name: BAMSI build (run $run/$RUNS)..."

        # Build (compression)
        local t0=$(date +%s%N)
        $BAMSI build --input "$bam" --output "$bsi" 2>/dev/null
        local t1=$(date +%s%N)
        local build_time=$(echo "scale=3; ($t1 - $t0) / 1000000000" | bc)
        local bsi_size=$(stat -c%s "$bsi" 2>/dev/null || stat -f%z "$bsi")

        echo "$name,bamsi,build,$run,$build_time,$bsi_size,0" >> "$csv"

        # Count query (latency)
        local patterns=("ACGT" "ACGTACGT" "TGCATGCA" "AAAA" "CCCC" "GGGG" "TTTT")
        for pat in "${patterns[@]}"; do
            local qt0=$(date +%s%N)
            $BAMSI count --index "$bsi" --pattern "$pat" > /dev/null 2>&1
            local qt1=$(date +%s%N)
            local query_time=$(echo "scale=6; ($qt1 - $qt0) / 1000000000" | bc)
            echo "$name,bamsi,count_${pat},$run,$query_time,0,0" >> "$csv"
        done

        # Verify
        local vt0=$(date +%s%N)
        $BAMSI verify --index "$bsi" > /dev/null 2>&1
        local vt1=$(date +%s%N)
        local verify_time=$(echo "scale=3; ($vt1 - $vt0) / 1000000000" | bc)
        echo "$name,bamsi,verify,$run,$verify_time,0,0" >> "$csv"

        # Reconstruct (decompression)
        local rt0=$(date +%s%N)
        $BAMSI reconstruct --index "$bsi" --output "$recon" 2>/dev/null
        local rt1=$(date +%s%N)
        local recon_time=$(echo "scale=3; ($rt1 - $rt0) / 1000000000" | bc)
        local recon_size=$(stat -c%s "$recon" 2>/dev/null || stat -f%z "$recon")
        echo "$name,bamsi,reconstruct,$run,$recon_time,$recon_size,0" >> "$csv"

        rm -f "$bsi" "$recon"
    done

    echo "[benchmark] $name: Results in $csv"
}

run_samtools_benchmark() {
    local bam="$1"
    local name="$2"
    local csv="$OUTPUT_DIR/${name}_samtools_${TIMESTAMP}.csv"

    if ! command -v samtools &> /dev/null; then
        echo "[benchmark] samtools not found, skipping"
        return
    fi

    echo "dataset,tool,operation,run,wall_time_s,file_size_bytes,rss_kb" > "$csv"

    for run in $(seq 1 $RUNS); do
        echo "[benchmark] $name: samtools (run $run/$RUNS)..."

        # CRAM compression
        local cram="/tmp/bamsi_bench_${name}.cram"
        local t0=$(date +%s%N)
        samtools view -C -o "$cram" "$bam" 2>/dev/null
        local t1=$(date +%s%N)
        local comp_time=$(echo "scale=3; ($t1 - $t0) / 1000000000" | bc)
        local cram_size=$(stat -c%s "$cram" 2>/dev/null || stat -f%z "$cram")
        echo "$name,samtools_cram,compress,$run,$comp_time,$cram_size,0" >> "$csv"

        # Pattern counting via samtools view -c (for comparison)
        local patterns=("ACGT" "ACGTACGT")
        for pat in "${patterns[@]}"; do
            local qt0=$(date +%s%N)
            samtools view "$bam" 2>/dev/null | grep -c "$pat" > /dev/null 2>&1 || true
            local qt1=$(date +%s%N)
            local query_time=$(echo "scale=6; ($qt1 - $qt0) / 1000000000" | bc)
            echo "$name,samtools_grep,count_${pat},$run,$query_time,0,0" >> "$csv"
        done

        rm -f "$cram"
    done

    echo "[benchmark] $name: Results in $csv"
}

# ─── Main ─────────────────────────────────────────────────────────────────────

echo "=== BAMSI Benchmark Suite v1.0.0 ==="
echo "Datasets:  $DATASET_DIR"
echo "Threads:   $THREADS"
echo "Runs:      $RUNS"
echo "Output:    $OUTPUT_DIR"
echo ""

# Check BAMSI binary
if [[ ! -x "$BAMSI" ]]; then
    echo "ERROR: BAMSI binary not found at $BAMSI"
    echo "       Run: cd build && cmake .. && make -j\$(nproc)"
    exit 1
fi

# Run on synthetic test data first (always available)
echo "=== Running on synthetic test data ==="
SYNTHETIC="${PROJECT_DIR}/data/test/synthetic_10reads.bam"
if [[ -f "$SYNTHETIC" ]]; then
    run_bamsi_benchmark "$SYNTHETIC" "synthetic_10"
    run_samtools_benchmark "$SYNTHETIC" "synthetic_10"
fi

# Run on any BAM files in the datasets directory
if [[ -d "$DATASET_DIR" ]]; then
    for bam in "$DATASET_DIR"/*.bam; do
        [[ -f "$bam" ]] || continue
        dataset_name=$(basename "$bam" .bam)
        echo ""
        echo "=== Running on $dataset_name ==="
        run_bamsi_benchmark "$bam" "$dataset_name"
        run_samtools_benchmark "$bam" "$dataset_name"
    done
else
    echo "[benchmark] No datasets directory found at $DATASET_DIR"
    echo "[benchmark] Place BAM files in benchmarks/datasets/ to benchmark real data"
fi

echo ""
echo "=== Benchmark complete ==="
echo "Results in: $OUTPUT_DIR"
echo "Generate tables with: python3 benchmarks/scripts/generate_tables.py"
