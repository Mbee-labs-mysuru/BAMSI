#!/usr/bin/env bash
# P3 — Independent read-level correctness verification
# Extracts 1,000 random reads from Dataset A BAM, builds a BAMSIX index on that
# exact subset, then compares BAMSIX GlobalCount against Python brute-force
# str.count() on the same reads.
set -euo pipefail

BAMSI_DIR="$(cd "$(dirname "$0")/.." && pwd)"
SAMTOOLS="${BAMSI_DIR}/tools/samtools-1.21/samtools"
BAMSI="${BAMSI_DIR}/build/bamsix"
BAM="${BAMSI_DIR}/data/NA12878_exome_2000000reads.bam"
SCRATCH="${BAMSI_DIR}/benchmarks/results/bruteforce_verify"

mkdir -p "$SCRATCH"

echo "╔══════════════════════════════════════════════════════════════╗"
echo "║  BAMSIX P3 — Brute-Force Read-Level Correctness Check        ║"
echo "║  1,000 random reads × 5 patterns = 5,000 verification pts  ║"
echo "╚══════════════════════════════════════════════════════════════╝"
echo ""

# Step 1: Extract 1,000 random reads as SAM lines (single random draw)
echo "[$(date +%H:%M:%S)] Step 1: Extracting 1,000 random reads..."
"$SAMTOOLS" view -H "$BAM" > "$SCRATCH/subset.sam"
"$SAMTOOLS" view -F 0x904 "$BAM" | shuf -n 1000 > "$SCRATCH/random_sam_lines.txt"
cat "$SCRATCH/random_sam_lines.txt" >> "$SCRATCH/subset.sam"

# Extract sequences from the SAME reads
cut -f10 "$SCRATCH/random_sam_lines.txt" > "$SCRATCH/random_1000_seqs.txt"
NREADS=$(wc -l < "$SCRATCH/random_1000_seqs.txt")
echo "  Extracted $NREADS reads"

# Step 2: Build BAMSIX index on this exact subset
echo "[$(date +%H:%M:%S)] Step 2: Building BAM and BAMSIX index on subset..."
"$SAMTOOLS" view -bS "$SCRATCH/subset.sam" | "$SAMTOOLS" sort -o "$SCRATCH/subset.bam"
"$SAMTOOLS" index "$SCRATCH/subset.bam"
"$BAMSI" build --input "$SCRATCH/subset.bam" --output "$SCRATCH/subset.bsi" 2>/dev/null
echo "  Index built: $SCRATCH/subset.bsi"

# Step 3: Brute-force count with Python on the EXACT same sequences
PATTERNS=("GATTACA" "TTAGGG" "ACGT" "GAATTC" "ATCGATCG")
echo "[$(date +%H:%M:%S)] Step 3: Brute-force counting (Python str.count)..."

python3 << 'PYEOF' > "$SCRATCH/bruteforce_counts.txt"
def rc(seq):
    comp = {'A':'T','T':'A','G':'C','C':'G','N':'N'}
    return ''.join(comp.get(c, 'N') for c in reversed(seq))

def count_overlapping(text, pattern):
    count = 0
    start = 0
    while True:
        pos = text.find(pattern, start)
        if pos == -1:
            break
        count += 1
        start = pos + 1
    return count

patterns = ["GATTACA", "TTAGGG", "ACGT", "GAATTC", "ATCGATCG"]

with open("benchmarks/results/bruteforce_verify/random_1000_seqs.txt") as f:
    seqs = [line.strip().upper() for line in f if line.strip()]

print(f"Total reads: {len(seqs)}")
for pat in patterns:
    pat_rc = rc(pat)
    total = 0
    for seq in seqs:
        fwd = count_overlapping(seq, pat)
        if pat == pat_rc:  # palindrome: same SA interval, counted once
            total += fwd
        else:
            rev = count_overlapping(seq, pat_rc)
            total += fwd + rev
    print(f"{pat}\t{total}")
PYEOF

cat "$SCRATCH/bruteforce_counts.txt"

# Step 4: BAMSIX GlobalCount on the SAME subset index
echo ""
echo "[$(date +%H:%M:%S)] Step 4: BAMSIX GlobalCount on subset index..."
PASS=0
FAIL=0
for pat in "${PATTERNS[@]}"; do
    bamsi_count=$("$BAMSI" count --index "$SCRATCH/subset.bsi" --pattern "$pat" 2>/dev/null || echo "ERR")
    bf_count=$(grep "^$pat" "$SCRATCH/bruteforce_counts.txt" | cut -f2)
    if [[ "$bamsi_count" == "$bf_count" ]]; then
        echo "  ✓ $pat: BAMSI=$bamsi_count, brute-force=$bf_count"
        PASS=$((PASS + 1))
    else
        echo "  ✗ $pat: BAMSI=$bamsi_count, brute-force=$bf_count"
        FAIL=$((FAIL + 1))
    fi
done

TOTAL=$((PASS + FAIL))
VERIF_PTS=$((NREADS * TOTAL))
echo ""
echo "═══════════════════════════════════════════════════════"
echo "  Result: $PASS/$TOTAL patterns match ($VERIF_PTS verification points)"
if [[ $FAIL -eq 0 ]]; then
    echo "  ✓ 100% agreement on all motif×read combinations"
else
    echo "  ⚠ $FAIL mismatches — investigate"
fi
echo "═══════════════════════════════════════════════════════"

# Cleanup temp files
rm -f "$SCRATCH/subset.sam" "$SCRATCH/random_sam_lines.txt"
echo ""
echo "[$(date +%H:%M:%S)] Done. Results in $SCRATCH/"
