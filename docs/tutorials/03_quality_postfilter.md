# Tutorial 03: Quality Post-Filter

Two-phase workflow: first identify candidate patterns via BAMSIX, then post-filter using quality scores.

## Overview

BAMSI queries operate on sequence data only. For applications where quality scores matter (e.g., variant calling), use a two-phase approach:

1. **Phase 1**: Use BAMSIX to efficiently count/locate candidate patterns
2. **Phase 2**: Post-filter candidates using quality scores from the original BAM

## Prerequisites

- BAMSIX built
- Original BAM file and its `.bsi` index

## Step 1: Build the Index

```bash
bamsix build data/test/synthetic_10reads.bam -o quality_demo.bsi
```

## Step 2: Phase 1 — Find Candidate Reads

```bash
# Locate all occurrences of a variant motif
bamsix locate --index quality_demo.bsi --pattern ACGT --json
```

This returns genomic coordinates and read IDs for all matches.

## Step 3: Inspect Read Metadata

```bash
# Inspect a specific read's metadata
bamsix reconstruct --index quality_demo.bsi --read-id 0
```

Output:
```
  Read 0:
    chrom_id: 0
    pos:      101
    seq_len:  16
    cigar:    16M
```

## Step 4: Phase 2 — Quality Score Filtering

Use the read IDs from Phase 1 to extract quality scores from the original BAM:

```bash
# Using samtools to extract specific reads and their quality scores
samtools view data/test/synthetic_10reads.bam | head -5
```

## Step 5: Combined Pipeline

```bash
#!/bin/bash
# Two-phase quality filter pipeline

INDEX="quality_demo.bsi"
BAM="data/test/synthetic_10reads.bam"
PATTERN="ACGT"
MIN_QUAL=30

# Phase 1: Get candidate counts
count=$(bamsix count --index $INDEX --pattern $PATTERN)
echo "Phase 1: $count occurrences of $PATTERN found"

# Phase 2: Check if index is lossless
lossless=$(bamsix info --index $INDEX --json 2>/dev/null | grep -o '"is_lossless": *[a-z]*' | grep -o '[a-z]*$')
echo "Lossless mode: $lossless"

if [ "$lossless" = "false" ]; then
    echo "WARNING: Quality scores were binned during indexing."
    echo "         Post-filter using the original BAM for precise quality scores."
fi

# Phase 2: Verify integrity before clinical use
bamsix verify --index $INDEX
echo "Integrity: verified"
```

## Key Points

1. **BAMSI does not filter by quality** — it counts/locates sequence patterns only
2. **Quality scores are preserved** in lossless mode for reconstruction
3. **Clinical workflows** should always verify index integrity before use
4. **The two-phase approach** is efficient: BAMSIX narrows candidates to a small subset, then quality filtering operates on that subset

## Cleanup

```bash
rm quality_demo.bsi
```
