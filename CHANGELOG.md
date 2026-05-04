# Changelog

All notable changes to BAMSI will be documented in this file.

The format is based on Keep a Changelog principles, adapted for this repository.[web:572]
This project is licensed under the Apache License 2.0.[web:566]

## [Unreleased]

### Added
- Repository scaffold for the C++ track.
- Source tree skeleton aligned with the frozen BAMSI architecture.
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
