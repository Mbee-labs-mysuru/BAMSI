# BAMSI Project Charter

## Scope
Implement BAMSI v1.0.0 as specified by:
- BAMSI Contract v3.3
- BAMSI Architecture v4.3

v1.0.0 focuses on:
- Lossless compression of the indexed aligned-read subset ℛ of BAM files (mapped, primary, non-supplementary with valid POS).[file:2]
- Exact substring pattern matching over read sequences via an FM-index, executed in the compressed domain.
- A C++20 reference implementation with a stable C ABI for bindings in other languages.

## Non-goals (v1.0.0)
- Indexing or reconstructing unmapped, secondary, or supplementary alignments (out of scope).
- Implementing approximate matching (k-mismatch / k-edit); only forward-compatible hooks and error codes are provided.
- Implementing quality-aware queries, joint multi-sample indexing, or v2.0 advanced features.

## Success Criteria
- **Correctness**
  - `decode(encode(ℛ)) = ℛ` for in-scope reads.
  - GlobalCount, GlobalExists, Locate, RegionalCount, and RegionalExists match brute-force baselines on validation datasets.
- **Performance**
  - Competitive compression ratio vs CRAM and Genozip on chosen datasets.
  - Acceptable GlobalCount and Locate latency on the benchmark laptop (and any additional benchmark nodes if used).
- **Reproducibility**
  - Bit-identical `.bsi` files across runs and machines for the same BAM inputs and BuildConfig.
- **Research**
  - At least one systems-style paper (Paper 1) and one theory-style paper (Paper 2) submitted to Q1-level venues.

## Gate-review Approvers
- Tech lead: Abhishek D P
- Project sponsor: same
