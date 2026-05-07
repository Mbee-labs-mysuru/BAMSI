# Tutorial 01: Motif Counting

Count telomeric repeats and other motifs across a BAM file using BAMSI.

## Prerequisites

- BAMSI built and available in `$PATH` or `build/` directory
- A BAM file (we use the synthetic test data)

## Step 1: Build the Index

```bash
bamsi build data/test/synthetic_10reads.bam -o motif_demo.bsi
```

## Step 2: Count Motifs

```bash
# Count ACGT occurrences (strand-complete: includes ACGT + its reverse complement)
bamsi count --index motif_demo.bsi --pattern ACGT
# Output: 13

# Count a specific sequence
bamsi count --index motif_demo.bsi --pattern CCGGTTAA
# Output: 5 (3 forward + 2 reverse complement)

# Check if a motif exists
bamsi exists --index motif_demo.bsi --pattern AACCGGTT
# Output: true
```

## Step 3: Batch Motif Counting

```bash
# Count multiple motifs
for motif in ACGT CCGG TTAA AACCGGTT GGTTAACC; do
    count=$(bamsi count --index motif_demo.bsi --pattern $motif)
    echo "$motif: $count"
done
```

## Step 4: JSON Output for Downstream Processing

```bash
# Machine-readable output
bamsi count --index motif_demo.bsi --pattern ACGT --json
# Output: {"count":13}

# Pipe to jq for processing
bamsi count --index motif_demo.bsi --pattern ACGT --json | jq '.count'
```

## Step 5: Locate Occurrences

```bash
# Find where ACGT occurs
bamsi locate --index motif_demo.bsi --pattern ACGT --sorted
```

Output:
```
strand  chrom   p_min   p_max   read_id
+       chr1    101     104     0
+       chr1    401     404     3
+       chr1    701     704     6
+       chr1    1001    1004    9
```

## Cleanup

```bash
rm motif_demo.bsi
```
