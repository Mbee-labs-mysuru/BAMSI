# BAMSIX CLI Reference

**Reference:** Contract §10.2

## Synopsis

```
bamsix <subcommand> [options]
```

## Subcommands

### `bamsix build`

Build a `.bsi` index from one or more BAM files.

```
bamsix build <input.bam> [-o <output.bsi>]
```

| Option | Default | Description |
|--------|---------|-------------|
| `-o, --output` | `<input>.bsi` | Output path |
| `--sample-step` | `64` | SA sampling step |
| `--window-size` | `100000` | Window size T (S-characters) |

**Example:**
```bash
bamsix build NA12878.bam -o NA12878.bsi
```

---

### `bamsix count`

Count global occurrences of a pattern.

```
bamsix count --index <file.bsi> --pattern <ACGT> [--json]
```

| Option | Description |
|--------|-------------|
| `--index` | Path to .bsi file |
| `--pattern` | DNA pattern (ACGTN characters) |
| `--json` | Output as JSON |

**Example:**
```bash
bamsix count --index data.bsi --pattern TTAGGG
# Output: 42
bamsix count --index data.bsi --pattern TTAGGG --json
# Output: {"count":42}
```

---

### `bamsix exists`

Check if a pattern exists (at least one occurrence).

```
bamsix exists --index <file.bsi> --pattern <ACGT> [--json]
```

**Example:**
```bash
bamsix exists --index data.bsi --pattern BRCA1MOTIF
# Output: true
```

---

### `bamsix locate`

Find genomic locations of pattern occurrences.

```
bamsix locate --index <file.bsi> --pattern <ACGT> [--sorted] [--json]
```

| Option | Description |
|--------|-------------|
| `--sorted` | Sort results by (chrom, position) |
| `--json` | Output as JSON |

**Output columns** (TSV):
```
strand  chrom   p_min   p_max   read_id
+       chr17   7577120 7577127 42
```

---

### `bamsix region-count`

Count pattern occurrences in a genomic region.

```
bamsix region-count --index <file.bsi> --pattern <ACGT> --chrom <chr> --start <pos> --end <pos> [--json]
```

| Option | Description |
|--------|-------------|
| `--chrom` | Chromosome name |
| `--start` | Region start (1-based inclusive) |
| `--end` | Region end (1-based inclusive) |

**Example:**
```bash
# Count ACGT occurrences in the TP53 region
bamsix region-count --index data.bsi --pattern ACGT --chrom chr17 --start 7570000 --end 7580000
```

---

### `bamsix region-exists`

Check if a pattern exists in a genomic region.

```
bamsix region-exists --index <file.bsi> --pattern <ACGT> --chrom <chr> --start <pos> --end <pos> [--threshold <N>] [--json]
```

| Option | Default | Description |
|--------|---------|-------------|
| `--threshold` | `1` | Minimum count for existence |

---

### `bamsix reconstruct`

Reconstruct read metadata from a `.bsi` index.

```
bamsix reconstruct --index <file.bsi> [--read-id <N>] [-o <output.bam>]
```

| Option | Description |
|--------|-------------|
| `--read-id` | Inspect a specific read |
| `-o, --output` | Output file path |

**Example:**
```bash
# Show summary of all reads
bamsix reconstruct --index data.bsi

# Inspect read 42
bamsix reconstruct --index data.bsi --read-id 42
```

---

### `bamsix info`

Display index metadata.

```
bamsix info --index <file.bsi> [--json]
```

**Output (plain):**
```
BAMSIX Index: data.bsi
  Format version:   6
  BAMSIX version:    0.1.0
  Lossless:         yes
  |S|:              1011701
  Reads:            10000
  Windows:          11
  Chromosomes:      chr1, chr2, chr3, chr4, chr5
```

**Output (JSON):**
```json
{
  "magic": "BMSI",
  "format_version": 6,
  "bamsix_version": "0.1.0",
  "is_lossless": true,
  "S_length": 1011701,
  "N_reads": 10000,
  "N_windows": 11,
  "sample_step_s": 64,
  "sentinel_row": 0,
  "window_size_T": 100000,
  "chrom_count": 5,
  "chromosomes": ["chr1", "chr2", "chr3", "chr4", "chr5"],
  "strand_mode": 0,
  "source_file_count": 1,
  "enable_sarange": false,
  "enable_bidirectional": false,
  "sa_samples": 15808
}
```

---

### `bamsix verify`

Verify `.bsi` file integrity.

```
bamsix verify --index <file.bsi> [--json]
```

**Exit codes:**
- `0` — File integrity verified
- `1` — Checksum mismatch or corrupt file

---

### `bamsix version`

Print version information.

```
bamsix version
```

## Common Workflows

### Motif Counting

```bash
# Count telomeric repeat TTAGGG across all reads
bamsix count --index genome.bsi --pattern TTAGGG

# Check if a specific variant exists
bamsix exists --index genome.bsi --pattern ATCGATCGATCG
```

### Region-Specific Queries

```bash
# Count BRCA1 region motifs
bamsix region-count --index genome.bsi --pattern ACGT \
    --chrom chr17 --start 43044295 --end 43170245

# Check TP53 region
bamsix region-exists --index genome.bsi --pattern TTAGGG \
    --chrom chr17 --start 7570000 --end 7580000
```

### Quality Verification

```bash
# Verify index integrity after transfer
bamsix verify --index genome.bsi

# Full info dump for audit
bamsix info --index genome.bsi --json > audit.json
```
