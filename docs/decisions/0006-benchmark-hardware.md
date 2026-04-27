# ADR 0006 – Benchmark Hardware

## Status
Draft

## Context
Execution Plan suggests separate benchmark machines: a high-RAM, NVMe build/decompression node and a separate query-latency node.[file:2] In this project, development and benchmarks initially run on a single laptop.

## Decision
For initial v1.0.0 benchmarks, the primary hardware configuration is:

- Build/decompression + query-latency node (same machine):
  - Model: HP Laptop 15-da0xxx
  - CPU: Intel Core i5-8250U (4 cores / 8 threads)
  - RAM: 24.0 GiB
  - Storage: 1.0 TB SSD/HDD (exact type to be documented)
  - OS: Fedora Linux 42 (Workstation, 64-bit)
  - Kernel: Linux 6.19.13-100.fc42.x86_64

If later a cloud or desktop machine with higher RAM/IO is used for full-scale benchmarks, that configuration will be documented as an additional benchmark node in this ADR.

## Rationale
- Reflects the actual hardware available at project start.
- Enables honest reporting of performance and memory constraints.

## Consequences
- Stage 0 ensures enough free disk space for at least one WGS dataset plus `.bsi` outputs and comparison tool outputs.
- Stage 5 benchmark section in the paper will describe this hardware configuration explicitly.
