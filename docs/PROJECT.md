# BAMSI Project Charter

## Scope
Implement BAMSI v1.0.0 as specified by:
- BAMSI Contract v3.3
- BAMSI Architecture v4.3

v1.0.0 focuses on:
- Lossless compression of the indexed aligned-read subset ℛ of BAM files (mapped, primary, non-supplementary reads with valid POS).[file:1][file:3]
- Exact substring pattern matching over read sequences via an FM-index, executed in the compressed domain.[file:1][file:3]
- Stream-aware handling of sequence, quality, metadata, and mapping information in the compressed representation.[file:1][file:3]
- A C++20 reference implementation with a stable C ABI for bindings in other languages.[file:1][file:3]

## Non-goals (v1.0.0)
- Indexing or reconstructing unmapped, secondary, or supplementary alignments.
- Approximate matching (k-mismatch / k-edit); only forward-compatible hooks and error codes are provided.
- Quality-aware queries, joint multi-sample indexing, or other advanced v2.0 features.
- Variant calling, structural variant calling, or general-purpose BAM editing.

## Success Criteria
- **Correctness**
  - `decode(encode(ℛ)) = ℛ` for all in-scope reads.[file:1][file:3]
  - GlobalCount, GlobalExists, Locate, RegionalCount, and RegionalExists match brute-force baselines on validation datasets.[file:1][file:3]
- **Performance**
  - Competitive compression ratio versus CRAM and Genozip on chosen datasets.[file:1][file:2]
  - Acceptable GlobalCount and Locate latency on the benchmark laptop and any additional benchmark nodes used for evaluation.[file:2]
- **Reproducibility**
  - Bit-identical `.bsi` files across runs and machines for the same BAM inputs and BuildConfig.[file:1][file:3]
- **Research**
  - At least one systems-style paper and one theory-style paper submitted to Q1-level venues.[file:2]

## Definition of Done
BAMSI v1.0.0 is done when:
- The frozen Contract and Architecture are implemented faithfully.[file:1][file:3]
- The C++20 reference implementation is complete for the v1.0.0 scope.[file:1][file:3]
- The command-line interface is usable for build, query, and verification tasks.[file:2][file:3]
- The core invariants of correctness, determinism, and losslessness are validated.[file:1][file:2][file:3]
- The benchmark datasets have been downloaded, checksum-verified, and recorded.[file:2]
- The evaluation results are reproducible and documented.[file:2]
- The final codebase, documentation, and results are stable enough for publication and controlled release.[file:2]

## Roles
- **Primary Engineer:** Student developer responsible for implementation, testing, benchmarking, and documentation.
- **Project Sponsor / Approver:** Academic guide or supervisor responsible for reviewing major decisions and approving stage gates.
- **Technical Reviewer:** Any additional expert reviewer who may validate design decisions, experiments, or publication readiness.

## Gate List
- **Gate 0:** Spec freeze and project kickoff.
- **Gate 1:** Repository skeleton, CI, and development environment ready.
- **Gate 2:** Core ingestion and vertical-slice pipeline working.
- **Gate 3:** Compression and reconstruction validated on initial datasets.
- **Gate 4:** Query path implemented and tested.
- **Gate 5:** Full benchmark suite completed.
- **Gate 6:** Paper-ready evaluation and figures prepared.
- **Gate 7:** Release candidate stabilized.
- **Gate 8:** Final submission and publication package prepared.
- **Gate 9:** Post-release maintenance and follow-up improvements.

## Gate-review Approvers
- Tech lead: Abhishek D P
- Project sponsor: same
