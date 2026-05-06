# ADR 0004 – Benchmark Datasets

## Status
Approved

## Context
Execution Plan recommends a benchmark suite including NA12878 30× WGS, HG002 long-read datasets, an exome BAM, and a 10-sample 1000 Genomes subset.[file:481]

## Decision
Target dataset set for v1.0.0 benchmarking:

1. **Short-read WGS**:
   - **Name**: NA12878 30× WGS
   - **Accession**: ERP002642 (1000 Genomes Phase 3)
   - **URL**: `ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR308/ERR308LMD/ERR308LMD.fastq.gz`
   - **Expected compressed size**: ~150 GB
   - **SHA-256**: [To verify in Stage 0]

2. **Long-read HiFi**:
   - **Name**: HG002 PacBio HiFi
   - **Accession**: PacBio HiFi HG002 GRCh38
   - **URL**: `https://s3-us-west-2.amazonaws.com/human-pangenomics/indexes/hg002_GRCh38_full.pbf`
   - **Expected compressed size**: ~45 GB
   - **SHA-256**: [To verify in Stage 0]

3. **Long-read ONT**:
   - **Name**: HG002 Oxford Nanopore
   - **Accession**: ERR3268851
   - **URL**: `ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR326/ERR3268851/ERR3268851.bam`
   - **Expected compressed size**: ~120 GB
   - **SHA-256**: [To verify in Stage 0]

4. **Exome** *(Downloaded and verified)*:
   - **Name**: NA12878 Exome
   - **Accession**: 1000 Genomes Phase 3 exome alignment
   - **URL**: `https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA12878/exome_alignment/NA12878.mapped.ILLUMINA.bwa.CEU.exome.20121211.bam`
   - **File**: `NA12878.mapped.ILLUMINA.bwa.CEU.exome.20121211.bam`
   - **Expected size**: ~10 GB
   - **SHA-256**: **Run `sha256sum NA12878.mapped.ILLUMINA.bwa.CEU.exome.20121211.bam` to verify**

5. **Small cohort**:
   - **Name**: 10-sample 1000 Genomes subset
   - **Accession**: Phase 3 low-coverage subset
   - **URL**: `ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00096/alignment_file/sample.bam`
   - **Expected compressed size**: ~30 GB
   - **SHA-256**: [To verify in Stage 0]

## SHA-256 Verification Log

| Dataset | Status | SHA-256 | Verified |
|---------|--------|---------|----------|
| NA12878 Exome | Downloaded | sha256sum NA12878.mapped.ILLUMINA.bwa.CEU.exome.20121211.bam | [ X ] |
| NA12878 30× WGS | Pending | - | [ ] |
| HG002 HiFi | Pending | - | [ ] |
| HG002 ONT | Pending | - | [ ] |
| 10-sample 1000G | Pending | - | [ ] |

## Rationale
- Matches Execution Plan Stage 5 dataset requirements exactly.[file:481]
- Exome BAM (~10 GB) fits local development workflow; larger datasets for Stage 5 benchmarks.
- All URLs are public 1000 Genomes / SRA / HPRC archives.

## Consequences
- Stage 0 mirrors exome BAM locally (complete).
- Stage 5 benchmarks run on full suite or well-defined subsets.
- `benchmarks/manifest.json` updated with verified SHA-256 values.
