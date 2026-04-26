# ADR 0003 – Target Platforms

## Status
Accepted

## Context
Execution Plan v2.0 assumes Linux x86_64, Linux aarch64, macOS x86_64, and macOS aarch64 as required, with Windows as best-effort.[file:2]

## Decision
For v1.0.0:

- Required:
  - Linux x86_64 (Fedora Linux 42 on HP Laptop 15-da0xxx – Intel Core i5‑8250U, 24 GiB RAM, 1 TB disk).

- Best-effort (stretch):
  - Linux aarch64
  - macOS x86_64 / aarch64
  - Windows (CLI only, if feasible through cross-compilation or separate build machines).

## Rationale
- Developer has guaranteed access to Linux x86_64 on the described laptop.
- Additional platforms would require extra hardware or cloud resources and are not critical for the first v1.0.0 release.

## Consequences
- CI matrix will start with Linux x86_64 only, with hooks for adding further platforms.
- Paper and documentation will state Linux x86_64 as the primary evaluation platform.
