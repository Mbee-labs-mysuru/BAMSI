# ADR 0002 – Dependency Pinning

## Status
Draft

## Context
BAMSI depends on a small set of core C/C++ libraries (C++ track): htslib, libsais, sdsl-lite, zstd, xxHash, and a cryptographic library such as OpenSSL or libsodium.[file:2]

## Decision
For v1.0.0 we will pin explicit minimum versions for:

- htslib
- libsais
- sdsl-lite
- zstd
- xxHash
- OpenSSL or libsodium

The exact version numbers will be filled in once initial prototypes are built and tested on Fedora 42.

## Rationale
- Reproducible builds and benchmarks.
- Avoid drifting behaviour mid‑project.
- Align with Execution Plan Stage 4 dependency‑freeze expectations.[file:2]

## Consequences
- After the dependency freeze (ADR 0007), upgrades require a documented CVE justification and an explicit decision to change versions.
- CI and benchmark environments will standardise on these pinned versions.
