# ADR 0007 – Dependency Freeze Policy

## Status
Accepted

## Context
Execution Plan: pinned versions become immutable from the start of Stage 4 until v1.0.0; upgrades require CVE justification and written exception.[file:2] In this project, timing is tied to stage progression rather than calendar weeks.

## Decision
- Dependency freeze activates on **entry to Stage 4** (i.e., immediately after Stage 3 gate passes).
- After this point, for core dependencies (htslib, libsais, sdsl-lite, zstd, xxHash, OpenSSL/libsodium):
  - No version changes are allowed unless:
    - There is a documented security advisory (CVE) affecting the frozen version; and
    - A short decision record is added explaining the upgrade.

## Rationale
- Prevents benchmark and correctness drift late in the project.
- Keeps builds reproducible for Stage 5 benchmarks and paper reproduction packages.

## Consequences
- Any dependency changes after Stage 3 must be deliberate and documented.
- Papers and supplementary materials can reference a stable dependency set for v1.0.0.
