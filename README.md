# BAMSI — BAM Sequence Indexer

A high-performance genomic sequence index that enables sub-linear pattern queries on BAM/CRAM alignments. BAMSI builds a compressed FM-index over read sequences, supporting global and regional motif counting, existence checks, and genomic location mapping — all without decompressing the underlying data.

## Quick Start

```bash
# Build from source
mkdir build && cd build
cmake .. && make -j$(nproc)

# Build an index
./bamsi build input.bam -o output.bsi

# Query patterns
./bamsi count --index output.bsi --pattern ACGTACGT
./bamsi exists --index output.bsi --pattern TTAGGG
./bamsi locate --index output.bsi --pattern ACGT --sorted
./bamsi region-count --index output.bsi --pattern ACGT --chrom chr17 --start 7570000 --end 7580000

# Verify integrity
./bamsi verify --index output.bsi

# Inspect index
./bamsi info --index output.bsi --json

# Reconstruct read metadata
./bamsi reconstruct --index output.bsi --read-id 42
```

## Features

- **FM-index** over concatenated read sequences with O(m) backward search
- **Global queries**: count, exists, locate across all reads
- **Regional queries**: count/exists within genomic coordinate windows
- **Strand-complete counting**: automatic reverse-complement search
- **CIGAR-aware mapping**: soft clips, insertions, deletions handled per spec
- **Binary .bsi format**: checksummed, versioned, cross-platform deterministic
- **C ABI**: stable `libbamsi.h` for FFI integration
- **9 CLI subcommands**: build, count, exists, locate, region-count, region-exists, reconstruct, info, verify

## Architecture

BAMSI follows a 10-stage build pipeline:

1. **Ingest** — Read BAM files via htslib
2. **Order** — Deterministic total ordering by (chrom_id, pos)
3. **Sequence Build** — Concatenate reads with `#` separators into S
4. **SA-IS** — Suffix array + BWT construction
5. **FM-Index** — C array, occurrence table, SA samples
6. **Stream Encoding** — S_seq, S_qual, S_meta, S_map compression
7. **Windows** — Partition S into genomic coordinate windows
8. **Bitvectors** — B_read (rank/select for read boundaries), B_window
9. **Validation** — TIER 1 invariant checks
10. **Seal** — Write `.bsi` with xxHash64 checksums

## Documentation

- [Format Specification](docs/format.md) — Byte-level `.bsi` file format
- [Algorithms](docs/algorithms.md) — FM-index, backward search, CIGAR mapping
- [CLI Reference](docs/cli.md) — All 9 subcommands with examples
- [C API Reference](docs/api.md) — Stable C ABI documentation
- [Clinical Operations](docs/clinical.md) — Provenance, audit, lossy-mode guidance

## Dependencies

| Library | Version | Purpose |
|---------|---------|---------|
| htslib | 1.21 | BAM/CRAM I/O |
| libsais | 2.8.6 | Suffix array construction |
| sdsl-lite | 2.1.1 | Succinct data structures |
| zstd | 1.5.6 | Stream compression |
| xxHash | 0.8.3 | Checksums |

## Building

```bash
# Prerequisites: cmake >= 3.16, g++ >= 12 (C++20), OpenSSL
git clone <repo-url> && cd BAMSI
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j$(nproc)

# Run tests
../build/test_tier1 ../data/test/synthetic_10reads.bam
../build/test_c_abi
```

## License

Apache 2.0. See [LICENSE](LICENSE) for details.

## Citation

If you use BAMSI in your research, please cite:

```bibtex
@software{bamsi2026,
  title  = {BAMSI: BAM Sequence Indexer},
  year   = {2026},
  url    = {https://github.com/Mbee-labs-mysuru/BAMSI}
}
```
