# BAMSIX Dependency Freeze — v1.0.0 Release Lock

**Document version:** 1.0
**Freeze date:** 2026-05-07
**Reference:** ADR 0007 (dependency-freeze-policy), Execution Plan §4.3

---

## Policy

All pinned dependency versions are **immutable** from this freeze date until the v1.0.0 release tag.
Upgrades require a security-advisory CVE justification and a written exception signed by the tech lead.

## Frozen Dependencies

| Library | Version | Source | License | SHA-256 (archive) |
|---------|---------|--------|---------|-------------------|
| htslib | 1.21 | `external/htslib/` (bundled) | MIT/BSD | (bundled, no archive) |
| libsais | 2.8.6 | `external/libsais/` (bundled) | Apache-2.0 | (bundled) |
| sdsl-lite | 2.1.1 | `external/sdsl-lite/` (bundled) | GPL-3.0 | (bundled) |
| zstd | 1.5.6 | `external/zstd/` (bundled) | BSD/GPL-2.0 | (bundled) |
| xxHash | (header-only, via zstd) | `external/zstd/lib/common/` | BSD-2-Clause | (bundled) |
| Google Test | 1.12.1 | `external/sdsl-lite/external/googletest/` | BSD-3-Clause | (bundled) |

## Build Toolchain

| Tool | Minimum Version |
|------|----------------|
| CMake | 3.16 |
| GCC | 10.0 (C++20) |
| Clang | 12.0 (C++20) |

## System Dependencies

**None required.** All dependencies are bundled in `external/` and statically linked.
The project does not modify system configuration.

## Verification

Run `scripts/verify_deps.sh` (if available) or confirm via:

```bash
# Verify all external libraries are present
ls external/htslib/ external/libsais/ external/sdsl-lite/ external/zstd/

# Verify build succeeds in isolation
mkdir -p build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
make -j$(nproc)
ctest --output-on-failure
```

## Exception Log

| Date | Library | From → To | CVE/Reason | Approved by |
|------|---------|-----------|------------|-------------|
| (none) | | | | |
