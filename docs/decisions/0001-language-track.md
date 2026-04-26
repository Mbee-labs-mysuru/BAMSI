# ADR 0001 – Language Track Choice

## Status
Accepted

## Context
BAMSI reference implementation must be in C++ (≥ C++20) or Rust (≥ 2021) per Contract v3.3 §10.1 and Architecture v4.3 §13.0.[file:1][file:3]

## Decision
We choose the **C++ track** as the reference implementation language.

## Rationale
- Developer experience and existing code in C++.
- sdsl-lite and libsais have mature C/C++ bindings.
- Execution Plan v2.0 gives shorter timelines for a 2-person C++ team than Rust.[file:2]

## Consequences
- Architecture §13.1a (C++ layout) and §13.2a (C++ build system) are normative for this project.[file:3]
- Rust track is considered future work (v1.1 or v2.0), not part of v1.0 reference.
