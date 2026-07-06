# Changelog

All notable changes to BAMSIX will be documented in this file.

The format is based on Keep a Changelog principles, adapted for this repository.[web:572]
This project is licensed under the Apache License 2.0.[web:566]

## [Unreleased]

## [1.0.0-rc1] — 2026-05-09

### Added
- RANGE_CYCLE quality codec with per-cycle delta encoding (Contract §2.7).
- TYPED_SPLIT meta codec with FLAG(uint16) + CIGAR(nybble) + aux-tag substreams (Contract §2.8).
- Streaming locate iterator C ABI: `bamsix_locate_iter_create/next/free` (Contract §10.3).
- CLI flags: `--qual-codec`, `--meta-codec`, `--map-codec`, `--threads` (Contract §10.2).
- `source_file_id` and `bam_offset` stored in read metadata section for ordering hash re-verification.
- Ordering hash re-verification on every index open (Contract §1.2).
- Read section checksum validation on load (Contract §9.4).
- `mode` and `strand_mode` fields in region-count/region-exists JSON output (Contract §4.2).
- ADR 0008: codec defaults.

### Fixed
- **FATAL:** S_map delta decode was broken for per-read ZSTD frames — switched to absolute (chrom_id, pos) per read (Contract §2.2 I10).
- **FATAL:** Read metadata section checksum was computed but never validated (`TODO` removed).
- **CRITICAL:** S_qual codec was plain ZSTD despite header claiming RANGE_CYCLE — now implements real per-cycle delta encoding.
- **CRITICAL:** S_meta codec stored FLAG as uint32 and omitted aux tags — now stores FLAG as uint16 with aux-tag substream.
- **CRITICAL:** Ordering hash was never re-verified on index open — now re-computed from loaded 4-tuples.

### Removed
- `IsSeparatorPosition` defensive checks from query path (Contract §4.2 F3 directive: provably unreachable).
- Duplicate tutorial stub files (kept full-content versions).

### Changed
- Documented local test invocation and ensured CLI smoketests run via CTest.
### Added
- Repository scaffold for the C++ track.
- Source tree skeleton aligned with the frozen BAMSIX architecture.
- Initial test layout with unit, integration, and synthetic tiers.
- Documentation tree with stub files for format, algorithms, CLI, API, clinical guidance, and tutorials.
- Benchmark scaffold with `benchmarks/manifest.json` and `benchmarks/Makefile`.
- ADR set for Stage 0 decisions:
  - ADR 0001: language track
  - ADR 0002: dependency pinning
  - ADR 0003: target platforms
  - ADR 0004: benchmark datasets
  - ADR 0005: licence and distribution
  - ADR 0006: benchmark hardware
  - ADR 0007: dependency freeze policy

### Changed
- Stage 0 documentation has been completed enough to support Stage 1 repository bootstrap.

### Fixed
- N/A

### Security
- Added `SECURITY.md` with responsible disclosure guidance.

## [0.1.0] - 2026-04-27

### Added
- Initial repository initialization.
- Stage 0 planning, ADR drafting, and project scaffolding work.
