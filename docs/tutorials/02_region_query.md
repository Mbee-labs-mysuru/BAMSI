# Tutorial 02: Region Query

Query pattern occurrences within specific genomic regions.

## Prerequisites

- BAMSI built
- A BAM file (we use the synthetic 10K-read test data)

## Step 1: Build the Index

```bash
bamsi build data/test/synthetic_10k.bam -o region_demo.bsi
```

## Step 2: Inspect the Index

```bash
bamsi info --index region_demo.bsi
```

Output:
```
BAMSI Index: region_demo.bsi
  Format version:   6
  |S|:              1011701
  Reads:            10000
  Windows:          11
  Chromosomes:      chr1, chr2, chr3, chr4, chr5
```

## Step 3: Global vs Regional Count

```bash
# Global count: how many times does ACGT appear across ALL reads?
bamsi count --index region_demo.bsi --pattern ACGT
# Output: some large number

# Regional count: how many times in chr1:1000-50000?
bamsi region-count --index region_demo.bsi --pattern ACGT \
    --chrom chr1 --start 1000 --end 50000
```

## Step 4: Region Existence Check

```bash
# Does this pattern exist in a specific region?
bamsi region-exists --index region_demo.bsi --pattern ACGTACGT \
    --chrom chr3 --start 1 --end 1000000

# With threshold: at least 5 occurrences?
bamsi region-exists --index region_demo.bsi --pattern ACGT \
    --chrom chr1 --start 1 --end 1000000 --threshold 5
```

## Step 5: Comparing Across Chromosomes

```bash
echo "Pattern: ACGT"
for chrom in chr1 chr2 chr3 chr4 chr5; do
    count=$(bamsi region-count --index region_demo.bsi --pattern ACGT \
        --chrom $chrom --start 1 --end 999999999)
    echo "  $chrom: $count"
done
```

## Step 6: JSON Pipeline

```bash
# Export region counts as JSON for downstream analysis
bamsi region-count --index region_demo.bsi --pattern ACGT \
    --chrom chr1 --start 1000 --end 50000 --json
# Output: {"region_count":42}
```

## Cleanup

```bash
rm region_demo.bsi
```
