# ADR 0004 – Benchmark Datasets

## Status
Draft

## Context
Execution Plan recommends a benchmark suite including NA12878 30× WGS, HG002 long-read datasets, an exome BAM, and a 10-sample 1000 Genomes subset.[file:2]

## Decision
Target dataset set for v1.0.0 (subject to hardware and storage constraints):

1. Short-read WGS:
   - Name: NA12878 30× WGS
   - Source: <accession / URL>
   - Expected compressed size: ~100–200 GB range (to be verified)
   - SHA-256: <to-fill>

2. Long-read:
   - Name: HG002 PacBio HiFi or HG002 ONT (choose one as primary)
   - Source: <accession / URL>
   - SHA-256: <to-fill>

3. Exome:
   - Name: Representative exome BAM (to be chosen)
   - Source: <accession / URL>
   - SHA-256: <to-fill>

4. Small cohort:
   - Name: “10-sample 1000 Genomes subset”
   - Source: <accession / URL>
   - SHA-256: <to-fill>

If full datasets do not fit comfortably on the available 1 TB disk or exceed practical processing time, we will define deterministic downsampled or region-sliced versions and record them here as updated entries.

## Rationale
- Matches the tool/dataset design in the Execution Plan for Stage 5 benchmarks and Paper 1.[file:2]
- Balances realism (standard public datasets) with the limits of a single 24 GiB / 1 TB laptop.

## Consequences
- Stage 0 mirrors at least the first short-read WGS dataset (or a well-defined subset) locally with SHA-256 verification.
- Stage 5 benchmarks will be defined over this list or any documented downsampled variants.
