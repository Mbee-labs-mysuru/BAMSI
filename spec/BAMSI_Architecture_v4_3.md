# BAMSI Architecture v4.3
## Hostile-Reviewer-Hardened + Objective-Aligned + Language-Constrained • Implementation-Ready
### Aligned to: BAMSI Contract v3.3

---

## Revision History

| Version | Changes |
|---|---|
| 3.3 | Prior version |
| 4.0 | Cross-document audit against Contract v3.0. Read-collection indexing 0..N-1; strand-complete aligned as default; window construction in S-character-space; SA-IS fixed (no "or equivalent"); GlobalExists first-class; source-manifest-hash byte-level; separator-inclusion in windows MUST; S_meta/S_map access separated; entropy pipeline explicit; chrom_id derivation; ISA/SARange normative optional; BWT derivation formula; partial-reconstruction; space-complexity; SEPARATOR_POSITION error; query-output ordering; GlobalCount strand-summing; "may" vs "must" separator rules tightened. |
| 4.1 | Hostile Q1 reviewer simulation fixes (22 findings), in lockstep with Contract v3.1. **Fatal fixes:** (F1) §9.4 space bound corrected to `\|S\|·H_k(S) + O(\|S\|)`. Single-BWT optimisation flag `shared_bwt` added. (F2) §0.11 + §5.3 split into BASE / ENHANCED tiers; SARange formalised. (F3) §6.4 — `S[pos] == '#'` defensive check REMOVED. (F4) §3 + §7.6 directory split: per-read for meta/map, block-level permitted for seq/qual. (F5) §1.1 build pipeline restructured (SA-IS as own stage). (F6) §6.6 window complexity unit-corrected to `\|W_r\| = O(L·d/T)`. **Significant fixes:** (S1) §4.6 + §5.1 two distinct rank APIs explicit. (S2) §3 / §4.5 entropy-k justified. (S3) §8.3 parallel SA-IS variants permitted. (S4) §9.3 validation matrix two-tier. (S5) §5.2 worked CIGAR example. (S6) §6.1 strand counting disambiguated. (S7) §2.7 partial-reconstruction split Table A / Table B. (S8) §2.4 / §7.1 `chrom_name_table` in header. (S9) §11 quality-aware as v2.0. **Clarifications:** (C1) §6.4 / §8.3 streaming + sorted output modes. (C2) §4.4 + §4.6 ISA samples fully specified. (C3) §10 / §11 secondary/supplementary scope + Paper 1 calibration. |
| 4.2 | Objective-alignment pass against the six BAMSI project objectives, in lockstep with Contract v3.2. (O1 — modular lossless compression, reference-free + reference-based, cross-platform/sequencer reproducibility) §1 system-overview expanded; new §3.5 Provenance & Reproducibility module records `bamsi_version`, `host_os_id`, `cpu_arch_id`, `is_lossless` in the `.bsi` header; sequencing-technology-agnostic statement. (O2 — in-compressed-domain + clinical workflows) New §10.1 references the Contract's §9 Operational Guarantees and lists implementation obligations: `bamsi info` machine-parseable output, `bamsi verify`, no implicit network. (O3 — stream-specific compression) §4.7 (S_qual / S_meta / S_map encoding) expanded with normative codec catalogues mirroring Contract §2.7–§2.9: per-cycle range / rANS-delta / ZSTD-dict / binned for S_qual; TYPED_SPLIT (CIGAR/FLAG/aux substreams) for S_meta; DELTA_RANGE for S_map. New error codes added. (O4 — exact + approximate at population scale) §3 `BuildConfig` adds `enable_bidirectional`, `recommended_seed_length`; new §4.6.7 Bidirectional FM Construction (built only when `enable_bidirectional = true`); §6.7 Approximate Query API stub returning `NOT_IMPLEMENTED_V1` per Contract §4.5. (O5 — benchmark vs. industry tools) §11 Paper 1 expanded with named tools (CRAM, Genozip, samtools, NGC), named real datasets, per-(tool,dataset) metrics, sub-experiments. (O6 — open-source CLI / docs / tutorials) New §12 Reference Implementation: directory layout, build system (CMake), test framework, CI pipeline, the normative CLI surface mirroring Contract §10.2, and the documentation deliverables. |
| **4.3** | **This version.** Implementation-language constraint hardened in lockstep with Contract v3.3. §13.1 source-tree layout splits into a C++ track variant and a Rust track variant; the `api/` module name changes from "C and Python bindings" to "stable C ABI shim layer" exposed by either track. §13.2 build system rewritten as a one-of-two choice: CMake ≥ 3.20 + C++20 OR Cargo + Rust 2021 edition; dependency list adjusted per track (htslib via `hts-sys` crate on the Rust track; sdsl-lite has no Rust equivalent so the Rust track uses `vers-vecs`, `fid` or hand-rolled succinct structures noted explicitly). §13.3 test framework: GoogleTest (C++) OR Cargo's built-in test harness + `proptest` (Rust). §13.7 distribution channels: PyPI removed; `crates.io` added for the Rust track; Bioconda binary unchanged (publishes the compiled CLI / library, track-agnostic). The CLI surface (§13.5) and documentation deliverables (§13.6) are unaffected — those are language-implementation-agnostic. |

---

## 0. Frozen Semantic Choices

These choices are **fixed and must not vary** across implementations. Any change requires a version increment of both this document and the Contract.

### 0.1 Scope of Data

BAMSI operates on the **indexed aligned-read subset** of BAM inputs, not on the full raw BAM byte stream.

A BAM record belongs to the indexed set ℛ if and only if:

```
FLAG & 0x4   == 0    (record is mapped)
FLAG & 0x100 == 0    (record is not secondary)
FLAG & 0x800 == 0    (record is not supplementary)
POS >= 1             (SAM 1-based position is valid)
```

ℛ contains only mapped primary alignments.

### 0.2 Query Alphabet

Query patterns are over:

```
Σ = {A, C, G, T, N}
```

- `N` is a **literal** symbol, not a wildcard.
- Codes: A=0, C=1, G=2, T=3, N=4, `#`=5 (separator, not in Σ), `$`=6 (conceptual sentinel, never stored).

### 0.3 Reverse Complement (Frozen)

Given sequence X = x₀x₁…x_{n-1}:

```
rc(X)[i] = complement(X[n-1-i])
complement: A↔T, C↔G, N↔N
```

### 0.4 Separator Policy

- `#` appears only between reads in S.
- Query patterns P never contain `#`.
- `#` is never counted as a match character.
- `#` **MUST** belong to the same window as the preceding read. This is a hard construction rule, not a suggestion.

### 0.5 Sentinel Policy

- `$` is a conceptual sentinel used only during suffix-array construction.
- `$` is lexicographically smaller than all codes 0–5.
- `$` is **not stored** in S or in the BWT payload.
- The FM-index behaves as if `$` were present at position |S| of S.
- The `sentinel_row` (the unique row where SA[row] = |S|) is tracked explicitly and is **never** reported as a match.

### 0.6 Strand Policy — Operational Default

**Strand-complete is the operational deployment default.** All implementations ship with strand-complete mode active unless the user explicitly enables the `SingleStrand` flag.

**Strand-complete definition:**

```
Q(P) = {P}           if P == rc(P)
Q(P) = {P, rc(P)}    otherwise
```

- Hits from P: `query_strand = Forward`
- Hits from rc(P): `query_strand = Reverse`
- If P == rc(P): search once, label `Forward`

`QueryStrand` refers to **query orientation only** — it is independent of the BAM FLAG strand field.

**Formal theorem basis:** All theorems and complexity bounds are stated **per orientation** (i.e., for a single Q ∈ Q(P)). Total cost for strand-complete mode is multiplied by |Q(P)| ≤ 2, which is a constant and does not change asymptotic class.

**Counting semantics (normative — matches Contract v3.1 §0.6):** BAMSI counts occurrences at the **read-sequence level**, not the deduplicated-genomic-locus level. Each occurrence of P in a read contributes 1 to the Forward tally; each occurrence of rc(P) in a read contributes 1 to the Reverse tally. Two matches on the same genomic locus from opposite orientations count as 2. Users requiring deduplicated genomic-position counts post-process `Locate` output externally.

### 0.7 Determinism

Given the same BAM inputs and the same BuildConfig, BAMSI MUST produce a bit-identical `.bsi` file. Required conditions:

- Fixed ordering rule (§2.4)
- Stable tie-breaking via `(source_file_id, bam_offset)` (§2.4)
- Fixed suffix-array algorithm: **SA-IS** (Nong, Zhang & Chan, 2009) is the reference algorithm. A parallel variant MAY be used in production provided it produces **bit-identical SA output** to the reference single-threaded SA-IS on every input in the validation suite. The reference algorithm is always used for TIER 2 validation (§9.3).
- Fixed compression parameters: `k`, coder type, serialization order
- Fixed little-endian byte serialization throughout
- Fixed sampling steps `s` and (optionally) `s'`

### 0.8 Query Semantics

- Queries operate over **occurrences**, not distinct reads.
- Multiple matches in the same read are counted separately.
- `GlobalCount` and `RegionalCount` return totals **summed across all Q ∈ Q(P)**.
- If P == rc(P), the query executes once and is labeled Forward.

### 0.9 Multi-BAM Policy

- BAMSI may ingest one or more BAM files.
- The order of BAMs in the build configuration is frozen as part of the input manifest.
- A read's stable ordering key includes `source_file_id` (zero-based index of BAM file in the input list) and `bam_offset` (stable record index within that file).

### 0.10 Manifest Hash Policy (Byte-Level Specification)

```
source_manifest_hash = SHA-256(
    concatenation over f=0..F-1, in source_file_id order, of:
        uint32_le(len(filename_f)) ||
        utf8_bytes(filename_f)     ||
        SHA-256(BAM_header_bytes_f)   // 32-byte raw digest
)
```

where `BAM_header_bytes_f` = raw bytes of the BAM header block before the first record of file `f`.

### 0.11 Complexity Policy — Two-Tier Guarantee

Worst-case query complexity is stated per orientation (Q ∈ Q(P)) in **two tiers**. Implementations MUST declare their tier in the benchmark header of any paper or report:

**BASE tier — required of every compliant implementation:**

| Query | Per-orientation bound |
|---|---|
| GlobalCount | O(\|P\|) |
| GlobalExists | O(\|P\|) |
| Locate | O(\|P\| + occ · s) |
| RegionalCount | O(\|P\| + occ · s) |
| RegionalExists | O(\|P\| + T_threshold · s) best; O(\|P\| + occ · s) worst |

**ENHANCED tier — only when `enable_sarange = true`:**

| Query | Per-orientation bound |
|---|---|
| GlobalCount | O(\|P\|) (unchanged) |
| GlobalExists | O(\|P\|) (unchanged) |
| Locate | O(\|P\| + occ · s) (unchanged) |
| RegionalCount | O(\|P\| + occ_r · s + \|W_r\| · log(\|S\|/s)) |
| RegionalExists | O(\|P\| + min(T_threshold, occ_r) · s + \|W_r\| · log(\|S\|/s)) |

Total cost for strand-complete: multiply per-orientation bound by |Q(P)| ≤ 2. Paper 2 (Theory) "tight" bound claims require and MUST declare the ENHANCED tier.

---

## 1. System Overview

### 1.1 Build Pipeline (SA-IS as a Distinct Stage — Single BWT Shared)

The pipeline has been restructured (v4.1) so that SA-IS runs as its own stage producing a **single** SA and BWT consumed in parallel by both S_seq encoding and FM-index construction. This enforces the BWT-lifecycle rule (§2.5 of the Contract) at the level of pipeline topology rather than relying on downstream modules to coordinate.

```
BAM inputs
  │
  ▼
Stage 1 — Ingestion (§4.1)
  │    produces RawRead[], chrom_name_table
  ▼
Stage 2 — Ordering (§4.2)
  │    produces OrderedRead[], ordering_hash, source_manifest_hash
  ▼
Stage 3 — Sequence Construction (§4.3)
  │    produces S (transient), readStarts
  ▼
Stage 4 — SA-IS Construction (§4.4)
  │    produces SA, BWT, sentinel_row
  │    (SA and BWT are the ONLY shared intermediates feeding both 5a and 5b)
  │
  ├────────────────┬───────────────────┐
  ▼                ▼                   │
Stage 5a           Stage 5b            │   running in parallel
S_seq Encoding     FM-Index            │   from the same SA/BWT
(§4.5)             Construction (§4.6) │
consumes BWT       consumes SA + BWT   │
produces S_seq     produces FMIndex    │
                                       │
                   After both consume, │
                   SA is discarded     │
  ▼                ▼                   │
  └────────────────┴───────────────────┘
  │
  ▼
Stage 6 — S_qual / S_meta / S_map Encoding (§4.7)
  │    parallel over streams; strong-independence respected;
  │    per-read dir_meta and dir_map (mandatory); block-level dir_seq/dir_qual permitted
  ▼
Stage 7 — Window Management (§4.8)
  │    parallel per chromosome
  ▼
Stage 8 — Bitvector Construction (§4.9)
  │    B_read, B_window from readStarts and window boundaries
  ▼
Stage 9 — Validation (§4.10)
  │    TIER 1 checks always; TIER 2 checks on sampled/synthetic inputs
  ▼
Stage 10 — Sealing (§4.10)
       → sealed '.bsi' file
```

**Why this matters.** In v4.0 the pipeline ran Stream Encoding (including S_seq, which needs BWT) **before** FM-Index Construction (which is where SA-IS was described). The BWT-lifecycle clause mandated a single shared BWT, but the pipeline order implied two independent BWT computations. Moving SA-IS to its own Stage 4 makes the shared-BWT semantics visible in the pipeline topology. Stages 5a and 5b are embarrassingly parallel from the same BWT.

**Section renumbering (v4.1).** Section numbers in §4 have shifted by one for modules downstream of the new Stage 4: what v4.0 called §4.4 (Stream Encoding) is now split into §4.5 (S_seq, consuming BWT from Stage 4) and §4.7 (S_qual/S_meta/S_map). FM-Index Construction is §4.6.

### 1.2 Query Pipeline

```
'.bsi'
    → load index structures (FM-index, bitvectors, windows, directories)
    → FM backward search per Q ∈ Q(P)
    → locate + mapping layer
    → exact query result
```

### 1.3 Module Separation

Each module below has a defined input type, output type, and invariant set. No module may bypass another module's output.

1. Ingestion
2. Ordering
3. Sequence Construction
4. Stream Encoding
5. FM-Index Construction
6. Window Management
7. Bitvector Construction
8. Mapping Layer
9. Query Execution
10. Validation and Sealing
11. Persistence (file format)

---

## 2. Data Model

### 2.1 Aligned Read Collection

$$\mathcal{R} = \{(r_i, c_i, p_i, m_i, Q_i) \mid i = 0, 1, \ldots, N-1\}$$

Where:

- $r_i \in \Sigma^*$, $\Sigma = \{A, C, G, T, N\}$; stored as codes 0–4
- $c_i$: reference identifier (chromosome name string, resolved from BAM header)
- $p_i$: leftmost aligned genomic coordinate; SAM POS, **1-based** (converted from BAM 0-based at ingestion, exactly once)
- $m_i$: metadata — CIGAR ops, FLAG, optional BAM tags
- $Q_i$: Phred quality scores; $|Q_i| = |r_i|$; values in $[0, 93]$

**Index convention:** i runs from 0 to N−1 throughout. 0-based indexing is used for all arrays, bitvectors, SA rows, and read_id values.

**Scope guarantee:** exact reconstruction is guaranteed only for ℛ; excluded BAM records are out of scope.

### 2.2 Inclusion Rule

As in §0.1.

### 2.3 Alphabet Encoding

| Symbol | Code |
|---|---|
| A | 0 |
| C | 1 |
| G | 2 |
| T | 3 |
| N | 4 |
| `#` | 5 |
| `$` | 6 (conceptual only, never stored) |

Only codes {0..5} appear in stored S and stored index structures.

### 2.4 Ordering and Read Identity

**chrom_id derivation (frozen per index):**  
`chrom_id(name)` = the 0-based rank of `name` in the lexicographically sorted list of all distinct reference names across all input BAM headers. This sort is performed once at header-parse time and frozen for the lifetime of the index.

The resulting `(chrom_id → name)` mapping is stored verbatim in the `.bsi` header as `chrom_name_table` (§7.1). **All query APIs accept string chromosome names** (e.g., `"chr10"`); the implementation resolves them to `chrom_id` via the stored table. Lexicographic ordering (under which `"chr10"` sorts before `"chr2"`) is an internal determinism artefact and is transparent to query callers. Two `.bsi` indices built from BAMs with different `@SQ` orderings have different raw `chrom_id` values but identical query behaviour when queries go through string names — which is the only supported interface.

**Total order ≺:**

$$r_i \prec r_j \iff (\text{chrom\_id}(c_i) < \text{chrom\_id}(c_j)) \lor (\text{chrom\_id}(c_i) = \text{chrom\_id}(c_j) \land p_i < p_j)$$

**Tie-breaking (frozen):**

- Single-BAM: ties broken by `bam_offset` (0-based record index within file).
- Multi-BAM: ties broken by `(source_file_id, bam_offset)`.

**Full sort key:** `(chrom_id, pos, source_file_id, bam_offset)` — all ascending.

$$\text{read\_id}(r_i) = \text{0-based rank of } r_i \text{ under } \prec$$

`read_id` is a bijection: ℛ → {0, 1, …, N−1}.

**Ordering Hash (byte-level specification):**

```
ordering_hash = SHA-256(
    for i = 0 to N-1 (in read_id order):
        uint32_le(chrom_id_i) ||
        uint64_le(pos_i)      ||
        uint32_le(source_file_id_i) ||
        uint64_le(bam_offset_i)
    // all fields concatenated, no separators
)
```

### 2.5 Concatenated Sequence S

$$S = r_0 \;\#\; r_1 \;\#\; \ldots \;\#\; r_{N-1}$$

**readStarts:**

$$\text{readStarts}[0] = 0$$
$$\text{readStarts}[i] = \text{readStarts}[i-1] + |r_{i-1}| + 1 \quad \text{for } i \geq 1$$

$$|S| = \sum_{i=0}^{N-1}|r_i| + (N-1)$$

**Properties:**

- `#` appears only at read boundaries
- Patterns P never contain `#`
- No exact match crosses a read boundary
- Each match belongs to exactly one read

**Lifecycle:** S is a transient build-time artifact. It is not stored as an uncompressed contiguous array in the `.bsi` file or materialized at query time.

### 2.6 Sentinel Semantics

For suffix-array construction, conceptually form S$:

- SA is defined over S$
- `$` is lexicographically smallest (less than code 0)
- `$` is **not stored** in S or in the BWT payload
- FM operations behave as if `$` were present at position |S|

**Sentinel row in queries:** `sentinel_row` is the unique logical row where SA[row] = |S|. Backward search is defined over all logical SA rows including `sentinel_row`, but:

- Patterns never contain `$`
- `locate` is never invoked on `sentinel_row`
- `sentinel_row` never produces a reported match

### 2.7 Partial Reconstruction and Query Streams (Disambiguated)

Two distinct access patterns over the sealed `.bsi` file exist and MUST NOT be conflated.

#### Table A — Reconstruction (recovering original ℛ data)

| Streams Available | Data Recoverable |
|---|---|
| $S_{\text{seq}}$ only | Read sequences $r_i$ for all $i$ |
| $S_{\text{map}}$ only | Reference identifiers $c_i$ (via `chrom_name_table`) and positions $p_i$ |
| $S_{\text{seq}} + S_{\text{map}}$ | Read sequences + genomic coordinates |
| $S_{\text{seq}} + S_{\text{qual}} + S_{\text{map}}$ | Above + quality scores |
| $S_{\text{seq}} + S_{\text{meta}} + S_{\text{map}}$ | Reads + CIGAR + FLAG + tags + coordinates (no quality) |
| All four | Full ℛ reconstruction |

Each row achievable from only the listed streams and their embedded headers — no cross-stream dependency (strong independence, I10).

#### Table B — Query Capability (what is loaded at query time)

| Query Operation | Required Structures |
|---|---|
| GlobalCount, GlobalExists | FMIndex (BWT, C, Occ, SA_samples) + chrom_name_table |
| Locate | FMIndex + B_read + dir_meta + dir_map + S_meta (lazy per-read) + S_map (lazy per-read) + chrom_name_table |
| RegionalCount, RegionalExists | Above + B_window + WindowTable + optionally SARange |
| Any query | **Never** requires S_seq, S_qual, raw S, or the uncompressed SA |

**Key distinction:** The FM-index is built from S and **replaces** S at query time. S_seq is a reconstruction-only stream (Table A); it is never loaded to execute any query. Attempting to do so violates I10 and the "raw S never resident at query time" invariant (§8.2). This closes a conflation present in v4.0's single partial-reconstruction table.

---

## 3. Core Data Types

All core types are **immutable** after construction.

```cpp
// ─── Ingestion ───────────────────────────────────────────────────────────────

struct RawRead {
    std::vector<uint8_t> seq;         // encoded bases: codes 0..4
    std::string          chrom;       // chromosome name string
    uint32_t             flag;        // BAM FLAG field
    uint64_t             pos;         // 1-based SAM POS (converted once at ingestion)
    CigarRecord          cigar;       // ordered list of (op, len) pairs
    std::vector<uint8_t> qual;        // Phred scores in [0,93]; |qual| == |seq|
    uint32_t             source_file_id;  // 0-based index in input BAM list
    uint64_t             bam_offset;      // stable record index within source file
};

// ─── Ordering ────────────────────────────────────────────────────────────────

struct OrderedRead : RawRead {
    uint32_t chrom_id;    // 0-based lexicographic rank of chrom name
    uint64_t read_id;     // 0-based rank under total order ≺
};

// ─── Sequence Builder ─────────────────────────────────────────────────────────

struct SequenceBundle {
    std::vector<uint8_t>  S;            // r0#r1#...#r_{N-1}, codes 0..5
    std::vector<uint64_t> readStarts;   // readStarts[i] = start of read i in S
};

// ─── Stream Directories ──────────────────────────────────────────────────────

struct StreamDirectoryEntry {
    uint64_t offset;   // byte offset of this read's block in the stream payload
    uint32_t length;   // byte length of this read's block
};

// Per-read directory: entry[i] corresponds to read_id = i
using StreamDirectoryPerRead = std::vector<StreamDirectoryEntry>;

// Block-level directory: entry[b] corresponds to reads in [b*B_dir, (b+1)*B_dir)
struct BlockDirectoryEntry {
    uint64_t block_offset;    // byte offset of the block's compressed payload
    uint32_t block_length;    // byte length of the block's compressed payload
    uint32_t first_read_id;   // 0-based read_id of the first read in this block
};
using StreamDirectoryBlockLevel = std::vector<BlockDirectoryEntry>;

// Variant type: per-read is mandatory for meta/map; block-level is permitted for seq/qual
using StreamDirectory = std::variant<StreamDirectoryPerRead, StreamDirectoryBlockLevel>;

struct Directories {
    // MANDATORY per-read directories (hot-path for query-time mapping)
    StreamDirectoryPerRead  meta;       // dir_meta[read_id] — required for Locate / RegionalCount
    StreamDirectoryPerRead  map;        // dir_map[read_id]  — required for Locate / RegionalCount

    // OPTIONAL granularity (per-read OR block-level) — seq/qual are reconstruction-only streams
    StreamDirectory         seq;        // dir_seq: per-read or block-level (B_dir)
    StreamDirectory         qual;       // dir_qual: per-read or block-level (B_dir)
    uint32_t                seq_block_size;   // B_dir for seq (0 = per-read); stored in S_seq header
    uint32_t                qual_block_size;  // B_dir for qual (0 = per-read); stored in S_qual header
};

// ─── Encoded Streams ──────────────────────────────────────────────────────────

struct EncodedStreams {
    std::vector<uint8_t> S_seq;    // BWT+MTF+RLE+arithmetic over S
    std::vector<uint8_t> S_qual;   // compressed quality arrays
    std::vector<uint8_t> S_meta;   // compressed CIGAR+FLAG+tags, keyed by read_id
    std::vector<uint8_t> S_map;    // compressed (chrom_id, pos), keyed by read_id
    Directories          dirs;     // per-read directories for all four streams
};

// ─── FM-Index ─────────────────────────────────────────────────────────────────

struct FMIndex {
    std::vector<uint8_t>  BWT;         // stored BWT payload, non-sentinel rows only
    std::array<uint64_t, 7> C;         // C[a] = # SA rows with BWT suffix < a (σ=6 + sentinel)
    OccStructure          occ;         // wavelet tree or RRR bitvectors over BWT
    uint64_t              sample_step_s;           // s ∈ [32,128]
    std::vector<uint64_t> SA_samples;  // SA_samples[k] = SA[k*s]
    uint64_t              sentinel_row;            // unique row where SA[row] == |S|
    // Optional ISA samples
    bool                  has_isa_samples;
    uint64_t              sample_step_s_prime;     // s' (only valid if has_isa_samples)
    std::vector<uint64_t> ISA_samples; // ISA_samples[k] = ISA[k*s'] (only if has_isa_samples)
};

// ─── Windows ──────────────────────────────────────────────────────────────────

struct Window {
    uint32_t chrom_id;
    uint64_t l;                  // inclusive start in S (S-character position)
    uint64_t r;                  // inclusive end in S (S-character position)
    uint64_t first_read_id;      // 0-based
    uint64_t last_read_id;       // 0-based
    uint64_t genomic_start;      // 1-based, for regional query filtering only
    uint64_t genomic_end;        // 1-based, for regional query filtering only
};

// Sorted by (chrom_id, genomic_start) for binary search
using WindowTable = std::vector<Window>;

// ─── Bitvectors ───────────────────────────────────────────────────────────────

class SuccinctBitvector {
public:
    static SuccinctBitvector build(const std::vector<bool>& bits);
    // rank1(pos): count of 1-bits in B[0..pos], CLOSED INTERVAL (both endpoints inclusive)
    uint64_t rank1(uint64_t pos) const;
    // select1(k): position of k-th 1-bit, k >= 1 (1-based rank)
    uint64_t select1(uint64_t k) const;
    void serialize(std::ostream& out) const;
    static SuccinctBitvector deserialize(std::istream& in);
};

struct BitVectors {
    SuccinctBitvector B_read;    // 1 at readStarts[i] for all i; 0 elsewhere
    SuccinctBitvector B_window;  // 1 at windows[j].l for all j; 0 elsewhere
};

// ─── Query Types ──────────────────────────────────────────────────────────────

enum class QueryStrand { Forward, Reverse };

struct MappingResult {
    uint32_t    chrom_id;
    uint64_t    p_min;        // 1-based inclusive
    uint64_t    p_max;        // 1-based inclusive
    uint64_t    read_id;      // 0-based
    QueryStrand query_strand;
};

struct Match {
    std::string chrom;
    uint64_t    p_min;        // 1-based inclusive
    uint64_t    p_max;        // 1-based inclusive
    uint64_t    read_id;      // 0-based
    QueryStrand query_strand;
};

// ─── Build Configuration ──────────────────────────────────────────────────────

enum class StrandMode { StrandComplete, SingleStrand };

// Stream-specific codec IDs (Contract §2.7 / §2.8 / §2.9)
enum class QualCodec : uint8_t {
    RANGE_CYCLE   = 0x01,  // default: range coder with per-cycle context
    RANS_DELTA    = 0x02,  // rANS over delta-encoded Q-scores
    ZSTD_DICT     = 0x03,  // ZSTD with per-dataset dictionary
    BINNED_RANGE  = 0x04,  // range coder over quantised bins (LOSSY when bins < 94)
};

enum class MetaCodec : uint8_t {
    TYPED_SPLIT   = 0x01,  // default: CIGAR + FLAG + aux substreams
    ZSTD_FALLBACK = 0x02,  // baseline ZSTD over raw concat (debug only)
};

enum class MapCodec : uint8_t {
    DELTA_RANGE   = 0x01,  // default: chrom_id RLE + Δpos range coder
    RAW           = 0x02,  // (uint32, uint64) per read (debug only)
};

struct BuildConfig {
    // Index parameters
    uint64_t   window_size_T;            // in S-characters; default 100,000
    uint64_t   sample_step_s;            // SA sampling step; default 64; in [32,128]
    uint64_t   sample_step_s_prime;      // ISA sampling step (0 = disabled)
    uint8_t    entropy_order_k;          // S_seq entropy order; in [4,8]; default 6
    StrandMode strand_mode;              // default: StrandComplete
    bool       enable_sarange;           // ENHANCED tier; default false
    bool       shared_bwt;               // single BWT for S_seq + FM; default true

    // Stream-specific codec choices (Contract §2.7-§2.9)
    QualCodec  qual_codec;               // default: RANGE_CYCLE
    uint8_t    qual_lossy_bins;          // 0 = lossless (94 levels); >0 sets lossy mode
    MetaCodec  meta_codec;               // default: TYPED_SPLIT
    MapCodec   map_codec;                // default: DELTA_RANGE

    // Directory granularity (Contract §2.3)
    uint32_t   seq_block_size;           // B_dir for dir_seq; 0 = per-read; default 1024
    uint32_t   qual_block_size;          // B_dir for dir_qual; 0 = per-read; default 1024

    // Approximate-search forward-compatibility hooks (Contract §4.5)
    bool       enable_bidirectional;     // build FM over S^R also; default false
    uint8_t    recommended_seed_length;  // hint for v2.0 seed-and-extend; default 16

    // Performance / parallelism
    bool       allow_parallel_sa;        // permit deterministic parallel SA-IS; default false
    uint32_t   build_threads;            // 0 = auto; otherwise pinned thread count

    // Reference-based encoding (Contract §2.6)
    std::string reference_fasta_path;    // empty = reference-free (mandatory mode)
    std::array<uint8_t, 32> reference_sha256;  // SHA-256 of reference; zeroed if reference-free
};

// ─── Structured Errors ────────────────────────────────────────────────────────

enum class ErrorCode {
    INVALID_BAM_INPUT,
    CORRUPT_BSI,
    UNSUPPORTED_CIGAR_OP,
    INVALID_PATTERN,        // pattern contains #, $, or characters not in Σ
    EMPTY_PATTERN,          // |P| == 0
    SEPARATOR_POSITION,     // M_ℓ called on a separator position
    CHECKSUM_MISMATCH,
    ORDERING_HASH_MISMATCH,
    MANIFEST_MISMATCH,
    VERSION_MISMATCH,
    BUILD_VALIDATION_FAILED,
    STREAM_DECODE_ERROR,
    NOT_IMPLEMENTED_V1,     // Approximate query attempted on a v1.0 build (Contract §4.5)
    REFERENCE_MISMATCH,     // .bsi was built reference-based; supplied reference SHA-256 differs
    LOSSY_RECONSTRUCTION,   // Reconstruction attempted but is_lossless=0; caller must opt in
    UNSUPPORTED_CODEC,      // qual_codec/meta_codec/map_codec ID is not recognised by this build
};

struct Error {
    ErrorCode   code;
    std::string message;
};
```

---

## 4. Build Architecture

### 4.1 Ingestion Module

**Input:** one or more BAM files in the order specified by BuildConfig  
**Output:** `RawRead[]`, chromosome mapping (name → chrom_id), input manifest

**Responsibilities:**

- Parse BAM records using a SAM/BAM-compliant library (e.g., htslib)
- Apply the inclusion rule (§0.1); discard records not in ℛ — **silently, not as errors**
- Normalize base sequences to codes 0..4; reject any unrecognized base as `INVALID_BAM_INPUT`
- Convert QUAL ASCII to numeric Phred values in [0, 93]
- Convert BAM 0-based POS to SAM 1-based POS exactly once at ingestion time
- Assign `source_file_id` from the frozen BAM file list (0-based index)
- Assign `bam_offset` as the sequential record index within each file (0-based, counting only records that appear in file order — not filtered)
- Preserve and freeze `chrom_id` mapping: sort all distinct reference names lexicographically across all input headers; assign 0-based ranks

**Invariants:**

- Every emitted RawRead satisfies the inclusion rule
- `|qual| == |seq|` for every emitted read
- No excluded record enters ℛ

### 4.2 Ordering Module

**Input:** `RawRead[]`, chromosome mapping, input manifest  
**Output:** `OrderedRead[]`, `ordering_hash`

**Responsibilities:**

- Sort by the 4-tuple key: `(chrom_id, pos, source_file_id, bam_offset)` — all ascending
- Assign `read_id` = 0-based rank in sorted order
- Compute `ordering_hash` per §0.10 byte-level specification
- Compute `source_manifest_hash` per §0.10 byte-level specification

**Invariants:**

- Sort is deterministic and reproducible
- `read_id` is a bijection over {0, 1, …, N−1}
- Both hashes match re-computation on index open; mismatch → `ORDERING_HASH_MISMATCH` or `MANIFEST_MISMATCH`

### 4.3 Sequence Builder Module

**Input:** `OrderedRead[]`  
**Output:** `SequenceBundle`

**Algorithm:**

```
initialize S = []  and  readStarts.resize(N)
for i = 0 to N-1:
    readStarts[i] = S.size()
    append r_i (encoded bases, codes 0..4) to S
    if i < N-1:
        append code 5 (#) to S
```

**Invariants:**

- `S[readStarts[i] .. readStarts[i] + |r_i| - 1] == r_i` for all i
- `#` (code 5) appears only at read boundaries
- `|S| = Σ|r_i| + (N-1)`

**Lifecycle note:** S is a transient build-time array. It must not be written to the `.bsi` payload uncompressed.

### 4.4 SA-IS Construction Module (Stage 4 — Single Shared BWT)

**Input:** `SequenceBundle.S`  
**Output:** `SA`, `BWT`, `sentinel_row`, `ISA` (if ISA samples requested)

This is a **dedicated pipeline stage** whose sole purpose is to compute the suffix array of S once, derive the BWT from it, optionally extract ISA samples, and then hand both the SA and BWT to downstream Stages 5a (S_seq encoding) and 5b (FM-index construction). No other stage may compute SA; no other stage may compute BWT. This enforces the BWT-lifecycle rule (§2.5 of Contract v3.1) at pipeline-topology level.

**Responsibilities:**

1. Conceptually form `S$` (|S|+1 characters; `$` is lexicographically smaller than all codes 0–5).
2. Build SA over `S$` using **SA-IS** (Nong, Zhang & Chan, 2009) — the reference single-threaded algorithm. If `BuildConfig.allow_parallel_sa` is true, a parallel variant MAY be used provided it produces bit-identical SA to the reference algorithm on the validation suite (§8.3).
3. Identify `sentinel_row` = the unique row where SA[row] = |S|.
4. Derive BWT: `BWT[i] = S[(SA[i] − 1 + |S|) mod |S|]` for all i ≠ sentinel_row. (The sentinel row's conceptual BWT character is `$`; it is **not** stored.)
5. If `BuildConfig.sample_step_s_prime > 0` (ISA samples requested): compute ISA by inverting SA (`ISA[SA[j]] = j` for all j), then sample `ISA_samples[k] = ISA[k · s']` for `k = 0..⌊|S|/s'⌋`. **Full ISA is not stored**; only the sampled array.
6. Hand SA (still in memory), BWT, `sentinel_row`, and `ISA_samples` to Stages 5a and 5b. **After both consumers complete**, SA is discarded; ISA (if built) is also discarded (only samples are retained).

**Lifecycle invariant:** Exactly one SA and one BWT are computed per build. Any implementation that recomputes SA or BWT in a downstream stage is non-compliant.

### 4.5 S_seq Encoding Module (Stage 5a — Consumes BWT from §4.4)

**Input:** `BWT` and `sentinel_row` from §4.4  
**Output:** `S_seq` byte payload + per-block/per-read directory

**Mandatory encoding pipeline (normative):**

1. **BWT** is received from §4.4 (not re-computed). The sentinel row's character (`$`) is elided.
2. **MTF (Move-to-Front):** Apply MTF transform to the BWT payload (non-sentinel rows only). Initial MTF alphabet order: [0, 1, 2, 3, 4, 5] (codes in ascending order).
3. **RLE (Run-Length Encoding):** Apply RLE to the MTF output. Zero-run encoding is preferred.
4. **Arithmetic coding:** Encode the RLE output with a 0th-order adaptive arithmetic model.

Entropy order parameter `k ∈ [4, 8]` is stored in the S_seq stream header. Default: `k = 6`. This pipeline achieves `|S_seq| = |S| · H_k(S) + o(|S|)` (the total index bound is `|S|·H_k(S) + O(|S|)` — see §9.4).

**Justification of default k=6 (required of Paper 1):** Paper 1 MUST include a sensitivity analysis showing compression ratio vs. k ∈ {4,5,6,7,8} on at least three representative datasets (short-read WGS, exome, long-read), plus encoder/decoder time. The theoretical bound holds for any fixed k; the empirical section justifies the chosen default.

**Stream header MUST store:** `k`, codec version, BWT `sentinel_row`, `|S|`, `N`, `seq_block_size` (0 = per-read; otherwise the block granularity of `dir_seq`).

**Single-BWT optimisation (optional, `BuildConfig.shared_bwt`).** If `shared_bwt = true`, S_seq does not store its own BWT bytes — it stores only MTF+RLE+arithmetic metadata and references the BWT already stored in the FM-index section. This eliminates the `|S|·H_0(BWT)` term from the total space bound; the bound remains `|S|·H_k(S) + O(|S|)` but with a smaller constant.

### 4.6 FM-Index Builder Module (Stage 5b — Consumes SA + BWT from §4.4)

**Input:** `SA` and `BWT` from §4.4, plus `ISA_samples` (if requested)  
**Output:** `FMIndex`

**Responsibilities:**

1. Receive SA, BWT, `sentinel_row`, and optionally `ISA_samples` from §4.4.
2. Build C array: `C[a] = |{j : BWT[j] < a, j ≠ sentinel_row}|` extended appropriately for the LF mapping.
3. Build Occ structure (wavelet tree or per-symbol RRR bitvectors) over the stored BWT.
4. Sample SA: `SA_samples[k] = SA[k · s]` for `k = 0, 1, …, ⌊|S|/s⌋`.
5. If `has_isa_samples = true`: retain the `ISA_samples` array received from §4.4.
6. If `BuildConfig.enable_sarange = true`: build the SARange wavelet tree over `SA_samples` (see §5.3).

**Two distinct rank APIs (normative — replaces v4.0's misleading "single consistent convention"):**

BAMSI uses two rank APIs with different semantics. They MUST NOT be implemented as a shared function, and implementers MUST NOT interchange them. This is a correctness requirement, not a style preference.

| API | Semantics | Used by | Value at position p |
|---|---|---|---|
| `rank_a(BWT, pos)` | **Half-open** | FM backward search + LF-mapping (this section only) | count of symbol `a` in `BWT[0..pos−1]` |
| `rank1(B, pos)` | **Closed-interval** | B_read, B_window, and all mapping/query code outside the FM inner loop | count of 1-bits in `B[0..pos]` inclusive |

**Equivalence lemma.** At any position $p$, the two conventions differ by at most 1: `rank1_closed(B, p) = rank1_half_open(B, p+1)`. The bitvector-rank choice is closed-interval precisely because this makes `read_id = rank1(B_read, pos) − 1` exact; the `−1` is not a bug, it is an intentional adjustment required by the convention choice. An implementation that applies the half-open convention to `B_read` will produce incorrect `read_id` values.

**BWT formula (normative, as in §4.4):**

```
BWT[i] = S[(SA[i] − 1 + |S|) mod |S|]   for all i where SA[i] ≠ |S|
```

**LF-mapping:**

```
LF(r) = C[BWT[r]] + rank_{BWT[r]}(BWT, r)
```

where `rank_a(BWT, r)` counts occurrences of `a` in `BWT[0..r−1]` (half-open, consistent with the FM backward search recurrence).

**Backward search:**

```
backward_search(Q[0..m-1]):
    lo ← 0;  hi ← |S| + 1     // logical row space for S$ has |S|+1 rows
    for i = m−1 downto 0:
        a ← encode(Q[i])
        lo ← C[a] + rank_a(BWT, lo)    // half-open: counts in [0..lo)
        hi ← C[a] + rank_a(BWT, hi)    // half-open: counts in [0..hi)
        if lo ≥ hi: return ∅
    return [lo, hi)
```

**Locate (forward-only, always available):**

```
locate(row):
    steps ← 0
    r ← row
    while r mod s ≠ 0 and r ≠ sentinel_row:
        r ← LF(r)
        steps ← steps + 1
    if r == sentinel_row: return |S|   // never produces a reported match
    return SA_samples[r / s] + steps
```

**Locate with ISA samples (bidirectional, when `has_isa_samples = true`):**

```
locate_bidir(row):
    // Walk forward via LF up to s steps, or forward via Psi from text position up to s' steps,
    // whichever reaches a sample first. Deterministic: choose LF path on tie.
    candidateA ← walk LF from row; stop at first r with r mod s == 0, or r == sentinel_row, or after s steps
    candidateB ← walk Psi from text position corresponding to row;
                  stop at first text position divisible by s', or after s' steps
    return whichever reached a sample, with accumulated offset added
```

Per-occurrence cost: `O(min(s, s'))`. The equivalence test in §9.3 (TIER 2) verifies that `locate_bidir` returns the same S-position as `locate` for all tested rows.

**Invariants:**

- FM backward search returns the exact SA interval for any query pattern.
- `locate` / `locate_bidir` reconstruct the correct S-position for any non-sentinel row.
- `sentinel_row` is never reported as a match.

#### 4.6.7 Bidirectional FM Construction (Approximate-Search Hook H1)

This sub-module implements forward-compatibility hook **H1** from Contract §4.5: building a second FM-index over the **reverse string** $S^R$ to support v2.0 bidirectional search (k-edit distance, k-mismatch with bidirectional seed extension). It is built **only when** `BuildConfig.enable_bidirectional = true` (default false).

**Construction.** The reverse FM-index does **not** require a second SA-IS run. From the SA already computed in Stage 4 (§4.4), the reverse SA is obtained by:

```
SA_R[i] = |S| − 1 − SA[|S| − 1 − i]    for i ∈ [0, |S|]
```

(See Lam, Sung, Tam, Wong, Yiu 2009, "High-throughput short read alignment via bi-directional BWT".) From `SA_R`, derive `BWT_R`, `C_R`, `Occ_R`, and `SA_R_samples` exactly as in Stage 5b for the forward index.

**Storage.** When present, the reverse FM-index is stored in a separate `.bsi` section (`FM-Index-Reverse`) with the same internal layout as the forward FM-Index section (§7.3). Header flag `enable_bidirectional` (uint8) records its presence.

**Cost.** Build time +O(|S|) for derivation + +O(|S|·H_0(BWT^R)) bits for storage. The reverse SA is derived in linear time and discarded after `BWT_R` and samples are extracted; it is not stored.

**Use in v1.0.** The reverse FM-index is **never used by any v1.0 query operation**. It is built and stored solely so that v2.0 software can open the same `.bsi` and run bidirectional search without rebuilding. v1.0 queries continue to use only the forward FM-index.

**Hooks H2-H4.** Stored in the index header:
- **H2 — `recommended_seed_length`** (uint8, default 16) — hint for v2.0 seed-and-extend k-mismatch search.
- **H3 — SA range API stability** — backward search returning `[lo, hi)` is a stable public API in `libbamsi`.
- **H4 — Locate stability** — `locate(row)` is a stable public API.

This satisfies Contract §4.5 entirely: a v2.0 release implementing `ApproxLocate` opens any v1.0 `.bsi` built with `enable_bidirectional = true` and runs without re-indexing.

### 4.7 S_qual / S_meta / S_map Encoding Module (Stage 6 — Parallel, Independent)

**Input:** `OrderedRead[]`  
**Output:** `EncodedStreams.S_qual`, `S_meta`, `S_map` + respective directories

These three streams are encoded in parallel; none depends on S_seq or FM-index output. The strong-independence rule (§4.7.4 below) is strictly enforced.

#### 4.7.1 S_qual — Quality Stream (Codec Catalogue)

Quality scores have well-known statistical structure (sequencing-cycle decay, smooth transitions, run-length plateaus). Per Contract §2.7, S_qual MUST exploit one of the following codecs, declared in the stream header (`qual_codec_id`):

| Codec ID | Name | Mechanism | Loss |
|---|---|---|---|
| `0x01` | `RANGE_CYCLE` (default) | Range coder with per-cycle context (read position modulo expected length); models cycle decay common to Illumina/MGI | Lossless |
| `0x02` | `RANS_DELTA` | rANS over delta-encoded Q-scores `Q_i[j] − Q_i[j-1]`; exploits smooth transitions | Lossless |
| `0x03` | `ZSTD_DICT` | ZSTD with a per-dataset dictionary trained on a sampled subset of $Q_i$; effective on long-read flat plateaus | Lossless |
| `0x04` | `BINNED_RANGE` | Range coder over quantised bins; **lossy when `qual_lossy_bins ≠ 0`** (e.g., 8-bin Illumina-style); lossless when bins = 94 | LOSSY iff `qual_lossy_bins ≠ 0` |

**Stream-header layout for S_qual:**

```
uint8   qual_codec_id              // {0x01, 0x02, 0x03, 0x04}
uint8   qual_lossy_bins            // 0 = lossless; >0 = lossy bin count
uint32  qual_block_size            // B_dir; 0 = per-read directory
variable codec_metadata            // per-codec params (e.g., dictionary bytes for 0x03)
variable payload                    // compressed quality data
```

**Lossy-mode obligation.** When `qual_lossy_bins ≠ 0`, the top-level `.bsi` header field `is_lossless` MUST be set to 0; `bamsi info` and `bamsi reconstruct` MUST surface the lossy condition prominently.

Compress each $Q_i$ independently per read_id (inside a block if block-level directories are used).

**Directory granularity:** Per-read (`dir_qual[i]`, required if `qual_block_size = 0`) OR block-level (`dir_qual[b]` for each block of `qual_block_size` consecutive reads, default `B_dir = 1024`). Block-level is permitted because S_qual is a reconstruction-only stream (Table B of §2.7).

#### 4.7.2 S_meta — Metadata Stream (Codec Catalogue)

CIGAR strings, FLAG fields, and BAM aux tags have known structure: small CIGAR-op alphabet with run-length, skewed FLAG distribution, sparse repeated-key aux tags. Per Contract §2.8, S_meta MUST exploit this:

| Codec ID | Name | Mechanism |
|---|---|---|
| `0x01` | `TYPED_SPLIT` (default) | Three logical substreams encoded independently within S_meta: (a) **CIGAR substream** — `(op, len)` pairs with op in 4-bit nybble (alphabet of 9 ops), len as 7-bit-stop varint; (b) **FLAG substream** — range coder with adaptive frequency over commonly-set FLAG bits; (c) **aux-tag substream** — key-dictionary + ZSTD over per-tag values, omitted when no aux tags present |
| `0x02` | `ZSTD_FALLBACK` | ZSTD over concatenated raw metadata; permitted as baseline; MUST NOT be used for benchmark headline numbers |

**Stream-header layout for S_meta:**

```
uint8   meta_codec_id              // {0x01, 0x02}
uint8   has_aux_tags               // 0 if no read has aux tags (skips aux substream)
variable codec_metadata
variable payload                    // contains CIGAR + FLAG + (optional) aux substreams
                                    // sub-offsets within the per-read block are stored in dir_meta
```

Stores per-read: CIGAR ops, FLAG field, optional BAM auxiliary tags. Does **not** store `chrom_id` or `pos` — those belong exclusively to S_map.

**Directory granularity: MANDATORY per-read** (`dir_meta[i]`). S_meta is accessed in the hot path by `Locate` / `RegionalCount` and requires O(1) seek per match. Block-level is **not permitted** for S_meta.

#### 4.7.3 S_map — Mapping Stream (Codec Catalogue)

Reads are sorted by `(chrom_id, pos)` (§2.4). The first read of each chromosome carries an absolute `pos`; subsequent reads carry small `Δpos` values clustering around the mean inter-read distance. Per Contract §2.9:

| Codec ID | Name | Mechanism |
|---|---|---|
| `0x01` | `DELTA_RANGE` (default) | (a) `chrom_id` substream — RLE over long runs of identical chrom_id; (b) `pos` substream — first read of each chromosome is absolute; subsequent are `Δpos = pos_i − pos_{i-1}` encoded with a range coder over a small adaptive alphabet of common deltas |
| `0x02` | `RAW` | Raw `(uint32 chrom_id, uint64 pos)` per read; debug / validation baseline only |

**Stream-header layout for S_map:**

```
uint8   map_codec_id               // {0x01, 0x02}
variable codec_metadata
variable payload
```

Stores per-read: `chrom_id` (uint32 conceptually) and `pos` (uint64, 1-based SAM POS). Does **not** store CIGAR, FLAG, or tags.

**Directory granularity: MANDATORY per-read** (`dir_map[i]`). For `DELTA_RANGE`, the per-read directory stores the byte offset within the chromosome chunk; the chunk's prefix-sum of deltas is reconstructed lazily for the queried block, costing O(B_dir) intra-chunk; the per-read mandate ensures O(1) chunk-relative-position lookup, and chunks are sized to keep this practical.

**Directory cost accounting (at 1B reads):**

| Directory | v4.0 (per-read, mandatory) | v4.1+ |
|---|---|---|
| dir_meta | 12 GB | 12 GB (unchanged — required) |
| dir_map | 12 GB | 12 GB (unchanged — required) |
| dir_seq | 12 GB | ~12 MB with B_dir=1024 |
| dir_qual | 12 GB | ~12 MB with B_dir=1024 |
| **Total** | **48 GB** | **~24 GB** |

This closes the hostile-reviewer objection that v4.0 required 48 GB of directory overhead at 1B reads — now down to 24 GB (half), with the reduction applied exactly to the streams that are not needed at query time.

#### 4.7.4 Strong Independence Rule (Frozen)

No entropy model for any one stream may use another stream's bytes or runtime statistics. Any cross-stream value a stream needs must be stored in that stream's own header. This is invariant I10.

### 4.8 Window Manager Module (Stage 7)

**Input:** `OrderedRead[]`, `SequenceBundle.readStarts`, `BuildConfig.window_size_T`  
**Output:** `WindowTable`

#### 4.8.1 Parameter T

`T` is the **target window size in S-characters** — a count of positions in S (each position holds a base code 0..4 or separator code 5). `T` is NOT a genomic coordinate distance.

Default: `T = 100,000` characters. Stored in the index header.

#### 4.8.2 Construction Algorithm (Per Chromosome)

For each chromosome c (in ascending `chrom_id` order):

```
reads_c = reads with chrom_id == c, in ascending read_id order (i.e., sorted order)

idx ← 0
while idx < |reads_c|:
    start_read = reads_c[idx]

    // Window boundaries are determined in S-CHARACTER SPACE
    l = readStarts[start_read.read_id]
    tentative_end_S = l + T − 1        // S-character-space endpoint

    // Extend to include all reads whose S-start falls within [l, tentative_end_S]
    last_idx = idx
    while last_idx + 1 < |reads_c|
      and readStarts[reads_c[last_idx + 1].read_id] ≤ tentative_end_S:
        last_idx ← last_idx + 1

    last_read = reads_c[last_idx]

    // r is the S-position of the last base of last_read
    r = readStarts[last_read.read_id] + |last_read.seq| − 1

    // MUST include the trailing '#' separator if this is not the last read in S
    if last_read.read_id < N − 1:
        r = r + 1

    // Genomic span — derived AFTER S-space construction, for regional query use only
    genomic_start = start_read.pos                          // 1-based
    genomic_end   = max over reads_c[idx..last_idx] of:
                    (read.pos + ref_span(read) − 1)         // 1-based

    emit Window{
        chrom_id      = c,
        l             = l,
        r             = r,
        first_read_id = start_read.read_id,
        last_read_id  = last_read.read_id,
        genomic_start = genomic_start,
        genomic_end   = genomic_end
    }
    idx ← last_idx + 1
```

**ref_span:**

```
ref_span(read) = Σ len(op)  for op ∈ {M, =, X, D, N}
```

**genomic_end fallback:** If `ref_span(read) == 0`, `genomic_end` for that read = `read.pos`.

**Chromosome boundary rule:** Windows MUST NOT cross chromosome boundaries. The last window of chromosome c ends at the last S-character of the last read of c (including its separator, unless that read is the globally last read in S). The first window of chromosome c+1 starts at `readStarts` of the first read of c+1. No gap or overlap at boundaries.

**Edge case — oversized read:** If `T < |r_i.seq|`, that read forms a single-read window. Log a WARNING. No invariant is violated.

#### 4.8.3 Window Table Ordering

Windows are stored sorted by `(chrom_id, genomic_start)` to enable O(log |W|) binary search for regional queries.

#### 4.8.4 Guarantees

- `[0, |S|-1] = ∪_j [l_j, r_j]` — exact partition, no gaps, no overlaps
- `r_j + 1 = l_{j+1}` for all consecutive windows
- No read is split: ∀i, ∃! window j s.t. `[readStarts[i], readStarts[i] + |r_i| - 1] ⊆ [l_j, r_j]`
- `#` always belongs to the preceding read's window (mandatory, not optional)
- Each window belongs to exactly one chromosome

### 4.9 Bitvector Manager Module (Stage 8)

**Input:** `SequenceBundle`, `WindowTable`  
**Output:** `BitVectors`

**Construction:**

```
B_read[pos]   = 1 iff pos == readStarts[i] for some i ∈ {0..N-1}
B_window[pos] = 1 iff pos == windows[j].l  for some j
```

**Rank convention (normative):** `rank1(B, pos)` = count of 1-bits in B[0..pos] inclusive (closed interval).

**Key identities (required invariants):**

```
rank1(B_read,   readStarts[i])  == i + 1        for all i ∈ {0..N-1}
select1(B_read, i + 1)          == readStarts[i] for all i ∈ {0..N-1}
rank1(B_window, windows[j].l)   == j + 1        for all j
select1(B_window, j + 1)        == windows[j].l  for all j
```

**Implementation choices:**

- `B_read`: Jacobson-style static bitvector — O(1) rank/select, O(n) bits + o(n) overhead.
- `B_window`: RRR-compressed bitvector — O(1) rank/select, compressed for sparse bitvectors.

**chrom_id derivation is finalized** in §4.2; `B_read` and `B_window` construction uses the already-determined `readStarts` and window boundaries.

### 4.10 Validation and Sealing Module (Stage 9–10) — Two-Tier Validation

Validation is split into two tiers to resolve the hostile-reviewer objection that v4.0's "MUST" language implied brute-force checks on full production datasets:

#### TIER 1 — Production Build Validation (Fast; runs every build)

Executes in ≤ 5% of total build time. MUST pass before sealing. Failure aborts the build with `BUILD_VALIDATION_FAILED`.

| Check | Validates |
|---|---|
| Stream round-trip on sampled 10,000 reads | encode/decode recovers ℛ on a sample |
| FM search on 100 fixed patterns vs. stored expected counts | backward search correctness on representative queries |
| Bitvector rank/select identities | §4.9 key identities on B_read, B_window |
| Window coverage + non-overlap (sampled check) | I7 on windows |
| No read split (sampled check) | I12 on a sample of reads |
| Separator assignment (sampled) | separator-in-preceding-window rule on a sample |
| Per-section checksum verification | on-disk byte integrity for all sections |
| Ordering hash match | re-computed ordering_hash equals stored value |
| Source manifest hash match | re-computed source_manifest_hash equals stored value |
| Strand-complete query on 10 test patterns | GlobalCount(P) + GlobalCount(rc(P)) equals brute-force on sample |
| Deterministic rebuild (optional on-demand) | two independent builds of same input → bit-identical .bsi |

#### TIER 2 — CI / Research Validation (Exhaustive; runs on synthetic or sampled inputs)

Runs in CI and for research reproducibility; NOT required for every production build. Brute-force and exhaustive tests operate on synthetic small inputs (`|S| ≤ 10^6`) and randomly sampled windows of production datasets — never on full production S.

| Check | Validates |
|---|---|
| FM backward search vs. brute-force linear scan on synthetic S (\|S\| ≤ 10^6) | full I3 correctness |
| Exhaustive occurrence coverage on constructed inputs | Theorem 2 (Completeness) |
| CIGAR edge cases (10S150M5S, 150I, all-soft-clip, …) | I5, I6 on hand-crafted CIGAR patterns |
| Mapping vs. BAM ground-truth on CIGAR corpus | full §5.2 correctness |
| Stream independence: each stream decoded in isolation | I10 with dir_meta/dir_map split |
| Sentinel row identification and skip | I11 end-to-end |
| Partial reconstruction (Tables A and B of §2.7) | each stream subset produces specified output |
| SARange range_count vs. locate-based count | if enable_sarange, ENHANCED tier consistency |
| ISA samples: locate_bidir equals locate | if has_isa_samples, locate equivalence |
| GlobalExists = (GlobalCount > 0) | query equivalence |
| RegionalExists = (RegionalCount > 0) | query equivalence |
| Block-level directory equivalence | dir_seq/dir_qual at per-read and B_dir=1024 return same reconstruction |
| Deterministic parallel SA-IS (if allow_parallel_sa) | bit-identical SA to reference single-threaded SA-IS |

This closes the hostile-reviewer objection about the impracticality of brute-force validation on 30 GB BAMs.

---

## 5. Mapping Architecture

### 5.1 Mapping Function

For an occurrence at S-position `pos` with pattern length `ℓ`, found by searching orientation `Q` (which determines `query_strand`):

$$M_\ell(\text{pos}, \text{query\_strand}) = (\text{chrom\_name}[c],\; [p_{\min}, p_{\max}],\; \text{read\_id},\; \text{query\_strand})$$

**Procedure:**

```
Step 1 — Identify read (closed-interval rank convention):
    read_id    = rank1(B_read, pos) − 1

Step 2 — Read start and match offset:
    read_start   = select1(B_read, read_id + 1)    // == readStarts[read_id]
    offset_start = pos − read_start
    offset_end   = offset_start + ℓ − 1

Step 3 — Retrieve coordinates from S_map (only):
    (chrom_id, p_i) ← decode per-read block from S_map using dir_map[read_id]
    // S_map provides: chrom_id (uint32), pos (uint64, 1-based SAM POS)
    // S_meta is NOT accessed in this step

Step 4 — Retrieve CIGAR from S_meta (only):
    cigar ← decode per-read block from S_meta using dir_meta[read_id]
    // S_meta provides: CIGAR ops, FLAG, optional tags
    // S_map is NOT accessed in this step

Step 5 — Apply CIGAR mapping:
    p_min = cigar_ref_pos(cigar, p_i, offset_start, LEFT)
    p_max = cigar_ref_pos(cigar, p_i, offset_end,   RIGHT)

Step 6 — Return:
    return MappingResult{
        chrom_id    = chrom_id,
        p_min       = p_min,       // 1-based inclusive
        p_max       = p_max,       // 1-based inclusive
        read_id     = read_id,
        query_strand = query_strand
    }
```

**Invariants:**

- If `pos` corresponds to a separator position in S (conceptually `S[pos] == code 5`; detected without materialising S — see below), return `SEPARATOR_POSITION` error; callers MUST check. Separator detection is implemented via the bitvector condition `rank1(B_read, pos) == rank1(B_read, pos-1)` AND `pos > 0` AND `pos ≠ readStarts[rank1(B_read, pos)-1]` — equivalently, `pos` is one position past the end of some read `r_i`. This check uses only bitvector operations; raw S is never materialised.
- `p_min ≥ 1` always
- `p_min ≤ p_max` always
- Deterministic: identical inputs produce identical outputs
- S_map and S_meta are accessed independently via their respective directories

### 5.2 CIGAR Mapping Definition

`cigar_ref_pos(cigar, p_anchor, offset, direction)` is **total** (defined for all valid inputs) and **deterministic** (identical inputs → identical output).

**Traversal:** Maintain `ref_pos = p_anchor`, `read_pos = 0`. Process ops in order:

| CIGAR Op | Read | Ref | Rule |
|---|---|---|---|
| M, =, X | ✓ | ✓ | If `offset ∈ [read_pos, read_pos+len-1]`: return `ref_pos + (offset − read_pos)`. Else: `ref_pos += len; read_pos += len`. |
| I | ✓ | — | If `offset ∈ [read_pos, read_pos+len-1]`: LEFT → last aligned ref base before this I (or `p_anchor` if none exists). RIGHT → first aligned ref base after this I (or last aligned base before if none after; or `p_anchor` if no aligned base at all). Else: `read_pos += len`. |
| D, N | — | ✓ | `ref_pos += len`. `read_pos` unchanged. |
| S | ✓ | — | If `offset ∈ [read_pos, read_pos+len-1]`: map to nearest aligned reference base. Equidistant tie → choose smaller coordinate. If no aligned base exists at all → return `p_anchor`. Else: `read_pos += len`. |
| H, P | — | — | No-op. Neither pointer changes. |

**Aligned base definition:** A position is "aligned" if it is covered by M, =, or X ops.

**No-aligned-base fallback:** If the read has no M/=/X/D/N ops at all (entirely soft/hard-clipped), all offset mappings return `p_anchor`.

**Properties:**

- Total: returns a defined value for every valid `(cigar, p_anchor, offset, direction)` tuple
- Deterministic including tie-breaks (equidistant S → smaller coordinate)
- Returns minimal enclosing 1-based inclusive interval `[p_min, p_max]` when applied to `[offset_start, offset_end]`
- No match is discarded due to insertions, clipping, or any CIGAR structure

**Interpretation note for users of `Locate` output (normative):** `[p_min, p_max]` is the **minimal enclosing reference interval** for the match — it is guaranteed to contain the reference positions of all aligned bases touched by the match. For patterns matched within soft-clipped or inserted bases, this interval may be larger than the pattern length:

- A pattern of length ℓ matched entirely inside a soft-clipped region of length C can map to an interval up to C reference bases wide.
- A pattern straddling an insertion maps so that `p_min` is the last aligned ref base before the insertion and `p_max` is the first aligned ref base after it.

**Worked example — CIGAR `10S150M5S`, p_i = 1000, |r_i| = 165.** Let pattern `P` of length 12 match at read offsets [x=3, y=14]: first 7 bases in the 10S clip, last 5 bases in the 150M region.

| Step | Value |
|---|---|
| `x = 3` falls in S op. Nearest aligned ref base is at read-offset 10 → ref_pos 1000. | `p_min = cigar_ref_pos(·, 1000, 3, LEFT) = 1000` |
| `y = 14` falls in 150M op at (14 − 10) = 4. | `p_max = cigar_ref_pos(·, 1000, 14, RIGHT) = 1000 + 4 = 1004` |
| Reported interval: | `[1000, 1004]`, width 5 for a 12-base pattern |

**Extreme case.** If `P` (length 10) matches entirely inside the 10S clip at read offsets [0, 9]: both LEFT and RIGHT return 1000, yielding `[1000, 1000]` — a 1-base interval for a 10-base pattern.

Users requiring precise sub-read positioning of clipped-region matches should post-process using the returned `read_id`, retrieve the full CIGAR from `S_meta[read_id]` via `dir_meta`, and interpret the clipped region against the pattern offsets independently. BAMSI's interval is guaranteed **sound** and **total** but not biologically tight for patterns dominated by clipped or inserted bases.

### 5.3 SARange Structure (Optional — ENHANCED Tier)

**Purpose.** Enable the ENHANCED tier RegionalCount bound `O(|P| + occ_r · s + |W_r| · log(|S|/s))` per orientation.

**Structure.** A **wavelet tree** over the sampled SA array `SA_samples[0..⌊|S|/s⌋]` (Grossi, Gupta, Vitter 2003). Equivalent 2D range-searchable structures are permitted; the variant used MUST be declared in the header (`sarange_variant` field, §7.1).

**Supported operation:**

```
range_count(lo, hi, l_j, r_j):
    returns: |{ r ∈ [lo, hi) : SA_effective[r] ∈ [l_j, r_j] }|
    complexity: O(log(|S|/s))
```

**Integration with RegionalCount (per orientation):**

```
[lo, hi) ← fm.backward_search(Q)                                        // O(|P|)
for each window j with (chrom_id, genomic_start, genomic_end) overlapping (chrom, [a,b]):
    c_j ← sarange.range_count(lo, hi, windows[j].l, windows[j].r)       // O(log(|S|/s))
    if c_j > 0:
        locate the c_j occurrences in window j                            // O(c_j · s)
        filter by exact [p_min, p_max] ∩ [a,b] intersection
```

Total per orientation: `O(|P| + |W_r| · log(|S|/s) + occ_r · s)` — the ENHANCED tier bound.

**Space.** `O((|S|/s) · log(|S|/s))` bits for the wavelet tree.

**Header flag.** `enable_sarange` (uint8) in `.bsi` header. When false (BASE tier), the structure is absent; RegionalCount uses backward search + locate + window filtering with `O(|P| + occ · s)` worst case.

**Normative status:**

- **BASE tier (enable_sarange=false):** No SARange present. Paper 1 (Systems) may operate in this tier; the ENHANCED bound MUST NOT be claimed.
- **ENHANCED tier (enable_sarange=true):** SARange present. Paper 2 (Theory) claims of tight RegionalCount require and MUST declare this tier.

### 5.4 Metadata Access Invariant

Per-read metadata access MUST always use the dedicated per-stream directory:

```
Coordinate access: dir_map[read_id]  → S_map block only
CIGAR access:      dir_meta[read_id] → S_meta block only
```

No query or mapping operation may access S_map to retrieve CIGAR information, or S_meta to retrieve coordinate information. This separation is required by the strong independence rule and by invariant I10. Note that `dir_meta` and `dir_map` are MANDATORY per-read (§4.7.2, §4.7.3); `dir_seq` and `dir_qual` MAY be block-level since they are not used by this hot path.

---

## 6. Query Architecture

### 6.1 Query Semantics

- Queries operate over **occurrences**, not distinct reads.
- Multiple matches in the same read are counted separately.
- GlobalCount and RegionalCount return totals **summed across all Q ∈ Q(P)**.
- RegionalExists reports true if the summed count ≥ T_threshold.
- If P == rc(P), the query executes once and is labeled Forward.

**Counting semantics (normative — disambiguates strand-complete totals):** BAMSI counts occurrences at the **read-sequence level**, not the deduplicated-genomic-locus level. Each occurrence of P in any `r_i` contributes 1 to the Forward tally; each occurrence of rc(P) in any `r_i` contributes 1 to the Reverse tally. Two matches on the same genomic locus from opposite orientations count as 2 (not deduplicated). Users requiring deduplicated genomic-position counts MUST post-process `Locate` output externally, grouping on `(chrom_id, p_min, p_max)` or by `read_id`. This matches `samtools view -c` on sequence-filtered reads.

### 6.2 Query Definitions

#### GlobalCount(P)

Returns the total number of exact occurrences of P across all reads in ℛ, summed over all orientations Q ∈ Q(P).

$$\text{GlobalCount}(P) = \sum_{Q \in Q(P)} (hi_Q - lo_Q)$$

**Complexity:** O(|P|) per orientation. O(|P|) total. No locate required.

---

#### GlobalExists(P)

Returns `true` if P occurs at least once in any read in ℛ under any orientation in Q(P); `false` otherwise.

$$\text{GlobalExists}(P) = \exists Q \in Q(P) : (hi_Q - lo_Q) > 0$$

**Complexity:** O(|P|) per orientation. O(|P|) total. No locate required.

---

#### Locate(P)

Returns one `Match` per occurrence: `(chrom, p_min, p_max, read_id, query_strand)` across all Q ∈ Q(P).

**Complexity:** O(|P| + occ · s) per orientation (BASE tier, unchanged in ENHANCED tier).

**Output ordering (normative — two modes):** Both modes return the same multiset of results.

- **Streaming mode (default).** Results emitted in SA-row order as computed, without sorting. Minimises time-to-first-result for large occ; deterministic (SA enumeration is fixed) but does not match the formal ordering below.
- **Sorted mode (opt-in via `--sort-output` or equivalent API flag).** Results fully buffered and sorted per this ordering before emission — required for reproducibility tests:
  1. `query_strand`: Forward before Reverse
  2. `chrom_id` ascending
  3. `p_min` ascending
  4. `p_max` ascending
  5. `read_id` ascending
  6. SA row order (final tie-break within read_id)

Implementations MUST support sorted mode for testing and MAY default to streaming mode for production queries. The mode actually used MUST appear in query output metadata so downstream consumers can interpret the order.

---

#### RegionalCount(P, chrom, [a, b])

Counts occurrences of P whose `[p_min, p_max]` intersects `[a, b]` on chromosome `chrom`, summed across Q ∈ Q(P).

**Match condition:** `p_min ≤ b` and `p_max ≥ a`.

**Complexity — two-tier:**

- **BASE tier (enable_sarange=false):** `O(|P| + occ · s)` per orientation worst case.
- **ENHANCED tier (enable_sarange=true):** `O(|P| + occ_r · s + |W_r| · log(|S|/s))` per orientation.

Per orientation; total × |Q(P)| ≤ 2.

---

#### RegionalExists(P, T_threshold, chrom, [a, b])

Returns `true` if `RegionalCount(P, chrom, [a, b]) ≥ T_threshold`; `false` otherwise. Terminates early when the T_threshold-th qualifying hit is found.

**Complexity:** `O(|P| + T_threshold · s)` best case; BASE worst `O(|P| + occ · s)`; ENHANCED worst `O(|P| + min(T_threshold, occ_r) · s + |W_r| · log(|S|/s))`.

---

### 6.3 Strand Execution

```
Q(P):
    if P == rc(P):  Q = {P};               query_strand for all = Forward
    else:           Q = {P, rc(P)};        Forward for P, Reverse for rc(P)
```

### 6.4 Query Execution

#### GlobalCount(P)

```
total ← 0
for each Q ∈ Q(P):
    [lo, hi) ← fm.backward_search(Q)
    total ← total + (hi − lo)
return total
```

#### GlobalExists(P)

```
for each Q ∈ Q(P):
    [lo, hi) ← fm.backward_search(Q)
    if hi − lo > 0: return true
return false
```

#### Locate(P)

```
results ← []                                       // buffered only in sorted mode
for each (Q, strand) ∈ Q(P) with strand labels:
    [lo, hi) ← fm.backward_search(Q)
    for each row r in [lo, hi):
        if r == sentinel_row: continue
        pos ← fm.locate(r)
        // Note: S[pos] cannot be '#' because P ∈ Σ* (enforced by INVALID_PATTERN at
        // query entry); backward search only returns rows whose suffix starts with Q ∈ Σ*.
        // The v4.0 defensive check `if S[pos] == '#': continue` has been REMOVED —
        // it required access to raw S at query time, contradicting §8.2 ("raw S never
        // resident at query time"). The check was provably unreachable for valid queries.
        if M_ℓ(pos) returns SEPARATOR_POSITION: continue  // still checked; M_ℓ uses no raw-S access
        m ← M_ℓ(pos, strand)
        entry ← Match{chrom_name[m.chrom_id], m.p_min, m.p_max, m.read_id, strand}
        if streaming mode: emit entry
        else:             results.append(entry)
if sorted mode:
    sort results per output ordering (§6.2)
    emit results
```

#### RegionalCount(P, chrom, [a, b]) — BASE tier

```
chrom_id ← chrom_name_table.lookup(chrom)     // O(log F) via sorted name table
if chrom_id is not found: return 0
count ← 0

for each (Q, strand) ∈ Q(P) with strand labels:
    [lo, hi) ← fm.backward_search(Q)
    for each row r in [lo, hi):
        if r == sentinel_row: continue
        pos ← fm.locate(r)
        // Window-based prune (closed-interval bitvector rank)
        window_id ← bv.B_window.rank1(pos) − 1
        W ← windows[window_id]
        if W.chrom_id ≠ chrom_id: continue
        if W.genomic_end < a or W.genomic_start > b: continue
        // Exact intersection check (S_meta + S_map accessed via their own dirs)
        m ← M_ℓ(pos, strand)
        if m.chrom_id == chrom_id and m.p_min ≤ b and m.p_max ≥ a:
            count ← count + 1
return count
```

#### RegionalCount(P, chrom, [a, b]) — ENHANCED tier (when enable_sarange=true)

```
chrom_id ← chrom_name_table.lookup(chrom)
if chrom_id is not found: return 0
count ← 0

// Find the set of windows W_r whose genomic span overlaps [a, b] on chrom_id
// This is a binary search in the WindowTable (sorted by (chrom_id, genomic_start))
W_r ← {j : windows[j].chrom_id == chrom_id
           and windows[j].genomic_end ≥ a
           and windows[j].genomic_start ≤ b}

for each (Q, strand) ∈ Q(P) with strand labels:
    [lo, hi) ← fm.backward_search(Q)                          // O(|P|)
    for each j in W_r:                                         // |W_r| iterations
        c_j ← sarange.range_count(lo, hi, windows[j].l, windows[j].r)   // O(log(|S|/s))
        if c_j == 0: continue
        // Locate only the occurrences that actually land in window j
        for each such row r (there are c_j of them):
            pos ← fm.locate(r)
            m ← M_ℓ(pos, strand)
            if m.p_min ≤ b and m.p_max ≥ a:
                count ← count + 1
return count
```

#### RegionalExists(P, T_threshold, chrom, [a, b])

Identical to the active-tier RegionalCount, except the intersection check becomes:

```
if m.chrom_id == chrom_id and m.p_min ≤ b and m.p_max ≥ a:
    count ← count + 1
    if count ≥ T_threshold: return true
return false
```

### 6.5 Complexity Summary — Two-Tier

| Query | BASE tier (per orientation) | ENHANCED tier (per orientation) |
|---|---|---|
| GlobalCount | O(\|P\|) | O(\|P\|) |
| GlobalExists | O(\|P\|) | O(\|P\|) |
| Locate | O(\|P\| + occ · s) | O(\|P\| + occ · s) |
| RegionalCount | O(\|P\| + occ · s) | O(\|P\| + occ_r · s + \|W_r\| · log(\|S\|/s)) |
| RegionalExists | O(\|P\| + T_threshold · s) best; O(\|P\| + occ · s) worst | O(\|P\| + min(T_threshold, occ_r) · s + \|W_r\| · log(\|S\|/s)) |

Total cost in strand-complete mode = per-orientation bound × |Q(P)| ≤ 2 (constant factor).

### 6.6 Practical Regional Behavior — Unit-Consistent Window Formula

Let $d$ be the average read coverage depth (reads per reference base) and $L$ the genomic span of the query region `[a, b]`. A genomic region of length $L$ reference bases contains approximately $L \cdot d$ S-characters (reads × bases, ignoring the O(N) separator contribution which is subsumed). Therefore, with window size $T$ (in S-characters):

$$|W_r| = O\!\left(\frac{L \cdot d}{T}\right) = O\!\left(\frac{L}{T_{\text{genomic}}}\right), \qquad T_{\text{genomic}} = \frac{T}{d}$$

**This is the unit-consistent form.** The v4.0 formula `|W_r| = O(region_length / T)` implicitly assumed unit coverage ($d = 1$) or unit-normalised quantities. For a 30× WGS, 150 bp reads, $T = 100{,}000$: $T_{\text{genomic}} \approx 3{,}333$ bp, so a 1 Mb region hits approximately 300 windows — consistent with the observed locality benefit.

**Key locality property:** `|W_r|` depends on region size $L$ and coverage $d$, not on total dataset size $|S|$ or number of reads $N$. Query cost scales with region size, not dataset size, for fixed $d$ and $T$.

With SARange (ENHANCED tier), the bound tightens to `O(|P| + occ_r · s + |W_r| · log(|S|/s))` per orientation.

### 6.7 Approximate Query API — v1.0 Stub (Contract §4.5)

The Contract specifies that v1.0 implementations MUST refuse approximate queries with a structured `NOT_IMPLEMENTED_V1` error rather than silently returning exact-only results. The reference implementation exposes the following API surface in `libbamsi.h` so that callers can write code today that will compile against v2.0 unchanged:

```c
// v1.0: returns NOT_IMPLEMENTED_V1; v2.0: returns matches with up to k_mismatch Hamming mismatches
bamsi_status_t bamsi_approx_locate_hamming(
    const bamsi_index_t* idx,
    const char* P, size_t P_len,
    int k_mismatch,
    bamsi_match_iterator_t** out);

// v1.0: returns NOT_IMPLEMENTED_V1 (also if enable_bidirectional == false, even in v2.0); v2.0: k-edit
bamsi_status_t bamsi_approx_locate_edit(
    const bamsi_index_t* idx,
    const char* P, size_t P_len,
    int k_edit,
    bamsi_match_iterator_t** out);
```

**v1.0 semantics:**

```
bamsi_approx_locate_hamming(...) {
    return NOT_IMPLEMENTED_V1;
}
bamsi_approx_locate_edit(...) {
    if (k_edit > 0) return NOT_IMPLEMENTED_V1;
    return NOT_IMPLEMENTED_V1;  // v1.0 refuses even k_edit=0 to keep one entry point clean
}
```

This preserves Contract §4.5 — callers get a clear, structured "not in v1.0" signal rather than a misleading exact match — and gives v2.0 a frozen API surface to implement against.

---

## 7. Storage Architecture

`.bsi` is a single sealed binary file. All multi-byte integers are little-endian. The file layout is:

```
[Header]
[S_seq section]
[S_qual section]
[S_meta section]
[S_map section]
[FM-Index section]
[Bitvector section]
[Window section]
[Directory section]
[Footer]
```

### 7.1 Header

| Field | Type | Description |
|---|---|---|
| `magic` | 4 bytes | `"BMSI"` |
| `version` | uint16 | Format version (current: 6 for v4.2 architecture) |
| `bamsi_version` | 16 bytes | Null-padded ASCII string of the BAMSI software version that wrote this file (e.g., `"1.0.0"`) |
| `host_os_id` | uint8 | 0=any, 1=Linux, 2=macOS, 3=Windows (informational; determinism makes it cosmetic) |
| `cpu_arch_id` | uint8 | 0=any, 1=x86_64, 2=aarch64 |
| `build_timestamp_utc` | uint64 | Unix epoch seconds; 0 for time-stripped reproducible builds |
| `is_lossless` | uint8 | 1 if every stream is lossless on v1.0 ℛ scope; 0 otherwise |
| `source_file_count` | uint32 | Number of input BAM files |
| `source_manifest_hash` | 32 bytes | SHA-256 per §0.10 |
| `ordering_hash` | 32 bytes | SHA-256 per §2.4 |
| `S_length` | uint64 | \|S\| |
| `N_reads` | uint64 | N |
| `N_windows` | uint32 | \|W\| |
| `sample_step_s` | uint32 | SA sampling step |
| `has_isa_samples` | uint8 | 0 or 1 |
| `sample_step_s_prime` | uint32 | ISA sampling step (if has_isa_samples) |
| `enable_sarange` | uint8 | 0 or 1; 1 ⇒ ENHANCED tier, SARange section present |
| `sarange_variant` | uint8 | 0 = wavelet tree over SA_samples (reference); reserved for alternatives |
| `shared_bwt` | uint8 | 0 or 1; 1 ⇒ S_seq references FM-index BWT rather than storing its own |
| `enable_bidirectional` | uint8 | 0 or 1; 1 ⇒ FM-Index-Reverse section present (Approximate-Search Hook H1) |
| `recommended_seed_length` | uint8 | Hint for v2.0 seed-and-extend; default 16 |
| `window_size_T` | uint64 | in S-characters |
| `entropy_order_k` | uint8 | in [4,8] |
| `qual_codec_id` | uint8 | per §4.7.1 catalogue |
| `qual_lossy_bins` | uint8 | 0 = lossless; >0 = lossy |
| `meta_codec_id` | uint8 | per §4.7.2 catalogue |
| `map_codec_id` | uint8 | per §4.7.3 catalogue |
| `strand_mode` | uint8 | 0=StrandComplete, 1=SingleStrand |
| `sentinel_row` | uint64 | FM-index sentinel row |
| `chrom_count` | uint32 | Number of distinct chromosomes |
| `chrom_name_table` | variable | `chrom_count` entries, each `(uint32 chrom_id, uint32 name_len, utf8 name_bytes)`, sorted by `chrom_id`. Used by query APIs to resolve string chromosome names. |
| `seq_block_size` | uint32 | B_dir for dir_seq; 0 = per-read |
| `qual_block_size` | uint32 | B_dir for dir_qual; 0 = per-read |
| `allow_parallel_sa` | uint8 | 0 or 1; informational — records whether a parallel SA-IS variant was used |
| `reference_based_encoding` | uint8 | 0 = reference-free, 1 = reference-based; if 1, reference_sha256 is non-zero |
| `reference_sha256` | 32 bytes | SHA-256 of reference FASTA (if reference_based_encoding=1; else zero) |
| `flags` | uint32 | Reserved feature flags |

### 7.2 Stream Sections (S_seq, S_qual, S_meta, S_map)

Each stream section:

```
uint64  payload_length
uint8   codec_id               // identifies the codec used
variable codec_metadata        // codec-specific parameters (e.g., k for S_seq)
variable payload_bytes         // compressed payload
uint64  section_checksum       // xxHash64 over (codec_id || codec_metadata || payload_bytes)
```

### 7.3 FM-Index Section

```
uint64  BWT_length             // number of stored BWT bytes (non-sentinel rows); 0 if shared_bwt=true and payload is absent
variable BWT_payload           // present iff shared_bwt=false OR this IS the owning storage under shared_bwt=true
variable C_array               // 7 × uint64 (for symbols 0..5 and sentinel)
variable Occ_metadata          // wavelet tree or RRR structure description
variable Occ_payload
uint64  SA_samples_count
variable SA_samples            // SA_samples_count × uint64
[if has_isa_samples:]
  uint64  ISA_samples_count
  variable ISA_samples         // ISA_samples_count × uint64
[if enable_sarange:]
  uint64  SARange_length       // bytes of SARange wavelet tree payload
  variable SARange_payload     // wavelet tree over SA_samples
uint64  section_checksum       // xxHash64 over this section's bytes
```

**Shared BWT semantics.** When `shared_bwt = true`, the BWT is stored exactly once in this section and the S_seq section stores only MTF+RLE+arithmetic metadata pointing at the FM-index BWT. When `shared_bwt = false`, S_seq stores its own BWT bytes independently; the FM-index BWT and the S_seq BWT are logically identical (both derived from the same SA-IS run, §4.4) but physically stored twice. Implementations SHOULD set `shared_bwt = true` in production.

### 7.4 Bitvector Section

```
variable B_read_serialized     // length-prefixed succinct bitvector
variable B_window_serialized   // length-prefixed succinct bitvector
uint64   section_checksum
```

### 7.5 Window Section

For each of `N_windows` windows (in sorted order):

```
uint32  chrom_id
uint64  l                      // S-character position, inclusive
uint64  r                      // S-character position, inclusive
uint64  first_read_id
uint64  last_read_id
uint64  genomic_start          // 1-based
uint64  genomic_end            // 1-based
```

Followed by: `uint64 section_checksum`.

### 7.6 Directory Section

Four directories in a fixed order: `dir_seq`, `dir_qual`, `dir_meta`, `dir_map`. Each directory is preceded by a small header indicating its granularity.

**dir_meta and dir_map (mandatory per-read):**

```
uint8   granularity             // must be 0 (PER_READ)
uint64  entry_count             // must equal N
for each i = 0..N-1:
    uint64 offset               // byte offset of read_i's block in the stream payload
    uint32 length               // byte length of read_i's block
uint64  section_checksum        // xxHash64 over the directory bytes
```

Implementations MUST store `dir_meta` and `dir_map` at per-read granularity. Block-level directories for these streams are non-compliant because `Locate` / `RegionalCount` require O(1) seek per match in the hot path.

**dir_seq and dir_qual (per-read OR block-level):**

```
uint8   granularity             // 0 = PER_READ, 1 = BLOCK_LEVEL
if granularity == PER_READ:
    uint64 entry_count          // must equal N
    for each i = 0..N-1:
        uint64 offset
        uint32 length
else:  // BLOCK_LEVEL
    uint32 block_size           // B_dir; must match seq_block_size or qual_block_size in header
    uint64 entry_count          // = ⌈N / B_dir⌉
    for each b = 0..entry_count-1:
        uint64 block_offset     // byte offset of compressed block b in the stream payload
        uint32 block_length     // byte length of compressed block b
        uint32 first_read_id    // 0-based read_id of the first read in this block
uint64  section_checksum
```

**Size accounting at 1B reads.** With `B_dir = 1024` for seq/qual:

- dir_meta + dir_map (required per-read): `2 × 10^9 × 12 B ≈ 24 GB`
- dir_seq + dir_qual (block-level): `2 × (10^9 / 1024) × 16 B ≈ 31 MB`
- Total: `~24 GB` instead of v4.0's `48 GB` — a 50% reduction applied exactly to the streams that don't serve queries.

### 7.7 Footer

```
uint64  global_checksum        // xxHash64 over all preceding bytes in the file
uint32  footer_magic           // 0xBSI1_0000 (BAMSI footer marker)
```

---

## 8. Error, Memory, and Parallelization Model

### 8.1 Error Model

All failures produce structured errors. No silent fallback is allowed.

| Error Code | Trigger |
|---|---|
| `INVALID_BAM_INPUT` | Malformed BAM record, unrecognized base, |qual|≠|seq| |
| `CORRUPT_BSI` | File truncation, structural parse failure |
| `UNSUPPORTED_CIGAR_OP` | CIGAR operation not in {M,I,D,N,S,=,X,H,P} |
| `INVALID_PATTERN` | Pattern contains `#`, `$`, or characters not in Σ |
| `EMPTY_PATTERN` | \|P\| == 0 |
| `SEPARATOR_POSITION` | M_ℓ invoked on a separator position (code 5) |
| `CHECKSUM_MISMATCH` | Section or global checksum fails |
| `ORDERING_HASH_MISMATCH` | Re-computed ordering_hash differs from stored |
| `MANIFEST_MISMATCH` | Re-computed source_manifest_hash differs from stored |
| `VERSION_MISMATCH` | .bsi version ≠ this implementation's supported version |
| `BUILD_VALIDATION_FAILED` | Any validation check in §4.10 fails |
| `STREAM_DECODE_ERROR` | Decompression of any stream fails |

### 8.2 Memory Model

| Structure | Residency at Query Time |
|---|---|
| FMIndex (BWT, C, Occ, SA_samples) | Fully resident or memory-mapped |
| ISA_samples (if present) | Fully resident or memory-mapped |
| SARange (if `enable_sarange = true`) | Fully resident or memory-mapped |
| B_read, B_window | Always resident (required for every query) |
| WindowTable | Always resident (required for regional queries) |
| `dir_meta`, `dir_map` (mandatory per-read) | Always resident (hot path for Locate / RegionalCount mapping) |
| `dir_seq`, `dir_qual` (permitted block-level) | Resident on-demand; only accessed for reconstruction, not for any query |
| `chrom_name_table` | Always resident (required to resolve string chromosome names in query APIs) |
| S_meta, S_map payloads | Lazily decoded per read_id via `dir_meta` / `dir_map` |
| S_seq, S_qual payloads | Lazily decoded only for reconstruction (Table A of §2.7); **never loaded by any query** |
| Raw S (uncompressed) | **Never resident at query time** |

### 8.3 Parallelization

- **Ingestion (Stage 1):** Parallel by BAM block/region with a final ordering barrier.
- **Ordering (Stage 2):** Sort is parallel (e.g., parallel merge-sort); deterministic output.
- **Sequence construction (Stage 3):** Parallel over reads (append to per-chromosome segments, concatenate deterministically at end).
- **SA-IS construction (Stage 4):** Default single-threaded reference SA-IS (Nong, Zhang & Chan 2009). **A deterministic parallel variant is permitted** when `BuildConfig.allow_parallel_sa = true`, provided the variant produces **bit-identical SA** to the reference algorithm on the TIER 2 validation suite. Additional parallelism: SA construction MAY be parallelised **over chromosomes** (each chromosome's suffixes are suffix-sorted independently and the results concatenated in `chrom_id` order) when the chromosome boundaries in S are preserved. This resolves the hostile-reviewer objection that single-threaded SA-IS on a 200 GB BAM (|S| ≈ 30 B characters) takes ~41 hours.
- **S_seq encoding (Stage 5a) and FM-index construction (Stage 5b):** Run in parallel from the same BWT produced by Stage 4.
- **S_qual / S_meta / S_map encoding (Stage 6):** Parallel across streams (strong independence) and within each stream, parallel over read blocks.
- **Windowing (Stage 7):** Parallel per chromosome.
- **Bitvector construction (Stage 8):** Parallel over read blocks for `B_read`; single-pass for `B_window`.
- **Validation (Stage 9):** TIER 1 checks parallel over sampled reads and patterns; TIER 2 checks parallel over test cases.
- **Query execution:** Independent queries trivially parallelizable. Per-query row iteration over `[lo, hi)` is parallelizable; in **streaming mode** results are emitted in SA-row order (parallel workers emit to a local buffer that is drained in row order); in **sorted mode** results are merged into a global sort before emission. All parallelism MUST preserve determinism (identical result multisets).

---

## 9. Invariants and Theorem-Linked Validation

### 9.1 Invariants

| ID | Invariant |
|---|---|
| I1 | Sequence Integrity — S encodes all included reads with `#` at correct boundaries; readStarts is consistent |
| I2 | Read Boundary Correctness — B_read has exactly one 1-bit at readStarts[i] for all i; 0 elsewhere |
| I3 | FM-index Correctness — backward search returns all and only SA rows whose suffixes begin with the query |
| I4 | No Cross-Read Matches — patterns never contain `#`; matches never cross separators (enforced by `P ∈ Σ*` validation + FM-index correctness; no runtime check on S required) |
| I5 | Mapping Soundness — reported intervals match BAM CIGAR semantics |
| I6 | Mapping Completeness — every occurrence position has a defined mapped interval |
| I7 | Window Coverage — windows form a disjoint partition of [0, \|S\|−1] |
| I8 | Bitvector Correctness — rank/select identities hold for B_read and B_window (closed-interval rank convention; distinct from the half-open BWT rank used by the FM recurrence) |
| I9 | Determinism — same input + same config → bit-identical .bsi; permits parallel SA-IS that is bit-identical to reference single-threaded SA-IS |
| I10 | Stream Independence — each stream decodable from its own bytes and embedded header only; `dir_meta` and `dir_map` are distinct per-read directories (mandatory); `dir_seq` and `dir_qual` may be per-read or block-level |
| I11 | Sentinel Consistency — FM-index behaves as if `$` exists; sentinel_row never yields a match |
| I12 | Window Integrity — each read belongs to exactly one window |
| I13 | BWT Lifecycle — SA and BWT are computed exactly once per build (§1.1 Stage 4) and shared by S_seq encoding (Stage 5a) and FM-index construction (Stage 5b); no stage recomputes them |
| I14 | No Raw-S Query Access — no query operation accesses raw S or its uncompressed contents at query time; separator detection uses B_read bitvector operations only |
| I15 | Rank API Separation — the half-open BWT-rank (`rank_a(BWT, pos)` counting in `[0, pos)`) is used only within the FM recurrence and LF-mapping; the closed-interval bitvector-rank (`rank1(B, pos)` counting in `[0, pos]`) is used everywhere else; the two MUST NOT be implemented as a shared function |

### 9.2 Theorem Dependencies

| Theorem | Required Invariants |
|---|---|
| Theorem 1 (Soundness) | I3, I4, I5, I8, I10, I11, I15 |
| Theorem 2 (Completeness) | I3, I6, I7, I11, I12, I15 |
| Space bound | I10, I13 (shared BWT bounds constants) |
| Determinism | I9, I13 (single SA-IS output) |
| Query-time S-free | I14, I15 |

### 9.3 Validation Matrix — Two-Tier Reference

The invariants and properties below are assigned to TIER 1 (production-build) or TIER 2 (CI/research) per §4.10. This cross-references §4.10 and makes the assignment explicit.

| Claim | TIER | Validation Method |
|---|---|---|
| I1 | 1 | Stream round-trip + read-span boundary checks on sampled 10,000 reads |
| I2 | 1 | rank1/select1 identity tests on B_read (sampled positions) |
| I3 | 2 | FM backward search vs. brute-force linear scan on synthetic S (\|S\| ≤ 10^6) + fixed 100-pattern sanity in TIER 1 |
| I4 | 1 | Separator exclusion tests on all `#` positions (cheap O(N)) |
| I5 | 2 | Mapping vs. BAM ground-truth on CIGAR edge-case corpus |
| I6 | 2 | Exhaustive occurrence coverage on small constructed inputs |
| I7 | 1 | Window overlap + coverage check (sampled ranges in TIER 1; full in TIER 2) |
| I8 | 1 | rank1/select1 property tests on B_read and B_window (sampled) |
| I9 | 1 (on-demand) | Rebuild byte-identity test on same input — run explicitly, not every build |
| I10 | 2 | Each stream decoded in isolation; `dir_meta` + S_meta decoded without `dir_map` + S_map and vice versa |
| I11 | 1 | `sentinel_row` identification; verify no query reports sentinel_row |
| I12 | 1 | Read-to-window uniqueness on sampled reads |
| T unit | 1 | Window construction test: verify T is in S-characters, not genomic distance |
| Strand counting | 2 | `GlobalCount(P)` equals brute-force read-sequence-level count of P and rc(P), summed (§6.1 counting semantics) |
| Separator must | 1 | Verify `#` at `readStarts[i]+\|r_i\|` is in same window as read i for all i < N-1 |
| Partial recon (Table A + Table B) | 2 | Each stream subset from §2.7 reconstructs specified data; Table B: S_seq never loaded by any query |
| SARange (if enabled) | 2 | `range_count` matches locate-based count for ≥ 100 regional queries |
| ISA samples (if enabled) | 2 | `locate_bidir` returns same positions as `locate` for all tested rows |
| Block-level dir_seq/dir_qual | 2 | Reconstruction output identical at per-read and B_dir=1024 granularity |
| Deterministic parallel SA-IS (if `allow_parallel_sa`) | 2 | Bit-identical SA to reference single-threaded SA-IS |
| GlobalExists | 1 | `GlobalExists(P) == (GlobalCount(P) > 0)` on test patterns |
| RegionalExists | 1 | `RegionalExists(P, 1, chrom, [a,b]) == (RegionalCount(P, chrom, [a,b]) > 0)` |
| chrom_name_table resolution | 1 | Query APIs accepting string chromosome names return same results across rebuilds with different @SQ orderings (§2.4 portability) |

Any TIER 1 failure MUST abort sealing and produce `BUILD_VALIDATION_FAILED`. TIER 2 failures are diagnostic / CI-blocking and MAY be reported without aborting a production build (unless configured to).

### 9.4 Space Complexity — Corrected Form

$$\text{Total space} = |S| \cdot H_k(S) + O(|S|) \text{ bits}$$

where the $O(|S|)$ term collects the FM-index BWT/Occ structure, bitvectors, auxiliary streams, and directories. This **corrects** the v4.0 statement `|S| H_k(S) + o(|S|) bits total`, which was mathematically unsound: `|S|·H_0(BWT)` for the FM BWT is a separate linear-in-|S| term that cannot be absorbed into the `o(|S|)` of `|S_seq|`. For genomic DNA, `H_0(BWT) ≈ 2.32` bits/symbol and `H_k(S) ≈ 1.7–2.0` bits/symbol for `k ∈ [4,8]`; the FM BWT term is actually larger than the compressed sequence term, making the v4.0 summary especially misleading.

| Component | Space | Status |
|---|---|---|
| S_seq (BWT+MTF+RLE+Arithmetic) | \|S\| H_k(S) + o(\|S\|) bits | dominant storage term |
| BWT in FM-index (Occ structure) | \|S\| H_0(BWT) + o(\|S\|) bits | linear-in-\|S\|; **elided when `shared_bwt=true`** |
| C array | 7 × 8 bytes | negligible |
| SA_samples | (\|S\|/s) × 8 bytes | sublinear |
| ISA_samples (optional) | (\|S\|/s') × 8 bytes | sublinear |
| B_read, B_window | 2\|S\| bits + o(\|S\|) overhead each | linear-in-\|S\| |
| WindowTable | O(\|S\|/T) entries × ~48 bytes/entry | sublinear |
| S_qual, S_meta, S_map | O(N × avg field size) compressed | O(\|S\|) |
| dir_meta, dir_map (required per-read) | 2N × 12 bytes | hot-path query directories |
| dir_seq, dir_qual (permitted block-level) | 2(N/B_dir) × 12 bytes if block; 2N × 12 bytes if per-read | reconstruction only |
| SARange (if enable_sarange=true) | O((\|S\|/s) × log(\|S\|/s)) bits | ENHANCED tier only |

**Dominant term:** `|S| H_k(S)` (the S_seq stream) with a linear-in-|S| tail `O(|S|)` collecting everything else. This is the standard form for a compressed full-text index carrying additional auxiliary data (Ferragina & Manzini 2005; Grossi & Vitter 2005).

---

## 10. Non-Goals (v1.0 — Out of Scope, with Extension Plan)

BAMSI v1.0 deliberately does not address the following. Each is a scope choice with a documented extension path; reviewers and users should understand both the limitation and the upgrade trajectory. (Mirrors Contract §7.)

- **Approximate pattern matching (v1.0 not implemented; v2.0 designed).** v1.0 returns exact substring matches only. The v1.0 build carries forward-compatibility hooks (Contract §4.5; Architecture §4.6.7 for H1, §3 BuildConfig for H2, §6 for H3/H4) so that a v2.0 release adds k-mismatch and k-edit overlays **without re-indexing**. v1.0 implementations refuse approximate queries with `NOT_IMPLEMENTED_V1` (§6.7).
- **Dynamic updates.** Indices are immutable; updates require rebuild. Log-structured overlays: future work (v3.0).
- **Variant calling or genotyping.** BAMSI provides substring counts, not variant inference. Pairs naturally with downstream callers (GATK, DeepVariant, bcftools).
- **Full BAM replacement.** BAMSI is a queryable compressed index, not a BAM viewer.
- **Multi-sample joint indexing.** Each `.bsi` indexes a fixed input manifest. Joint multi-sample indices: future work (v2.0); v1.0 reserves header fields for sample tagging.
- **Reconstruction / indexing of unmapped, secondary, or supplementary records.** ℛ contains only primary mapped reads (§0.1). Consequences:
  - Well-suited to motif finding, coverage QC, k-mer analysis, and any workflow that conventionally filters out secondary/supplementary reads.
  - **NOT** a drop-in replacement for tools analysing structural variants, chimeric reads, or long-read alignments that rely on supplementary alignments.
  - Paper 1 MUST include a motif-scope calibration: `GlobalCount_BAMSI(P)` (primary reads only) vs. `samtools view -c` over all alignments, on 10 representative motifs.
  - Indexing of secondary/supplementary records: future work (v2.0).
- **Quality-aware queries (v1.0 two-phase).** S_qual is compressed with stream-specific codecs (§4.7.1), but quality scores are not used in the v1.0 query model. Two-phase workflow: `Locate(P)` returns matching `read_id`s; the caller retrieves `Q_i` via `dir_qual[read_id]` for quality post-filtering. Quality-threshold extensions to `RegionalCount`: future work (v2.0).
- **Deduplicated genomic-position query semantics.** Strand-complete counts are read-sequence-level, not genomic-position-deduplicated (§6.1 counting semantics). Users requiring deduplication post-process `Locate` output.

---

## 11. Research and System Guarantees

### Paper 1 — Systems Contribution (Aligned to Contract §8 / Objective 5)

**Tools compared.** BAMSI vs. CRAM (samtools 1.21+ with reference), Genozip (latest stable, lossless mode), `samtools view` (query-latency baseline), NGC (Roguski & Deorowicz 2014). Lossy modes of any tool MUST be reported separately and labelled.

**Datasets — public, real-world, mandatory.**

| Dataset | Source | Purpose |
|---|---|---|
| NA12878 30× WGS Illumina HiSeq X BAM | 1000 Genomes / GIAB | short-read WGS benchmark |
| HG002 PacBio HiFi 30× | GIAB | long-read PacBio benchmark |
| HG002 ONT PromethION ~30× | GIAB | long-read ONT benchmark |
| NA12878 exome (Agilent SureSelect or equivalent) | public | capture / target panel benchmark |
| 10-sample 1000 Genomes subset | 1000 Genomes | scalability / cohort benchmark |

**Per-(tool, dataset) metrics:**

- Compression ratio: bits per base pair, separately for sequence, quality, and total.
- Build wall-clock and peak RSS at single-thread and 8/16-thread settings.
- Decompression-to-BAM wall-clock and peak RSS (where applicable).
- Query latency (BAMSI only vs. `samtools view` baseline): GlobalCount, GlobalExists, Locate, RegionalCount, RegionalExists for ≥ 100 patterns spanning length 8–32, region sizes 1 kb / 1 Mb / chromosome / genome.
- Memory footprint at query time (resident set excluding kernel cache).

**Disclosure required per result row:** BAMSI tier (BASE / ENHANCED), `is_lossless`, codec choices (`qual_codec_id`, `meta_codec_id`, `map_codec_id`, `entropy_order_k`), `BuildConfig.shared_bwt`.

**Reproducibility:** All benchmarks reproducible from a public Docker/Singularity container with one command (`make benchmarks`). Container image SHA-256 in the paper.

**Sub-experiments:**

- Entropy-order sensitivity: `k ∈ {4,5,6,7,8}` × three datasets (validates §4.5 default).
- Motif-scope calibration: BAMSI vs. `samtools view -c` (all alignments) on 10 motifs (validates Contract §7 / §10).
- Soft-clip interval characterisation: `[p_min, p_max]` width vs. precise sub-read mapping on CIGAR edge cases (validates §5.2).
- Directory granularity trade-off: per-read vs. block-level `dir_seq`/`dir_qual` (validates §4.7).
- Parallel SA-IS (if `allow_parallel_sa = true`): speedup vs. reference single-threaded; bit-identity test results.
- Window-pruning benefit: RegionalCount latency vs. region size at fixed `T` (validates §6.6).
- **Stream-codec ablation (new in v4.2):** TYPED_SPLIT vs. ZSTD_FALLBACK on S_meta; DELTA_RANGE vs. RAW on S_map; per-cycle range vs. ZSTD-dict on S_qual (validates the codec catalogues of §4.7).
- **Bidirectional FM cost (new in v4.2):** build-time and `.bsi`-size overhead with `enable_bidirectional = true` (informational — quantifies the v2.0 forward-compatibility cost paid in v1.0).

**Performance targets (informational):** On NA12878 30×, total `.bsi` size ≤ 1.10× Genozip output (acknowledging FM-index overhead); median GlobalCount latency < 50 ms on commodity hardware.

### Paper 2 — Theory Contribution

- Formal proofs of Theorems 1 (Soundness) and 2 (Completeness) with invariant citations.
- Entropy bound: `|S_seq| = |S| H_k(S) + o(|S|)` via BWT→MTF→RLE→Arithmetic, parametric in `k`.
- Query complexity per orientation, tier-aware:
  - **BASE tier:** GlobalCount/GlobalExists `O(|P|)`, Locate `O(|P| + occ·s)`, RegionalCount `O(|P| + occ·s)`, RegionalExists as in §6.5.
  - **ENHANCED tier (SARange-enabled):** RegionalCount `O(|P| + occ_r·s + |W_r| log(|S|/s))`; RegionalExists correspondingly tightened. All Paper 2 "tight" claims require and MUST state the ENHANCED tier.
- **Overall space bound:** `|S| H_k(S) + O(|S|)` bits total (§9.4).
- **Rank-convention lemma:** half-open FM rank and closed-interval bitvector rank differ by at most 1 (§4.6).
- **Pipeline-topology lemma:** SA and BWT computed exactly once per build (§1.1 Stage 4) and shared by Stages 5a/5b.
- **Approximate-search forward-compatibility (informational, not Paper 2's core):** the v1.0 index satisfies hooks H1–H4 of Contract §4.5 and is therefore extension-compatible with a v2.0 approximate-matching overlay without re-indexing.

---

## 12. Final Guarantee

Under the above frozen semantics (aligned to Contract v3.3), BAMSI guarantees:

1. **Exact substring search** over indexed reads under strand-complete mode (read-sequence-level counting).
2. **Correct CIGAR-consistent genomic mapping** for every reported occurrence, with sound minimal-enclosing intervals (§5.2 interpretation note for clipped/inserted regions).
3. **Deterministic and reproducible** index construction: same input + same `BuildConfig` → bit-identical `.bsi`.
4. **Lossless reconstruction** of the indexed aligned-read subset ℛ from all four streams (Table A of §2.7), gated by the `is_lossless` header flag.
5. **No undefined behavior** for any valid in-scope input.
6. **Independent stream decodability** — any subset of streams decodes from its own bytes alone (I10).
7. **Invariant-backed theorem soundness and completeness** per Contract v3.3.
8. **Two-tier complexity guarantee** (BASE and ENHANCED) with explicit SARange formalisation (§5.3).
9. **Space bound `|S| H_k(S) + O(|S|)`** — the standard compressed full-text index form with auxiliary data (§9.4).
10. **Two distinct rank APIs** explicitly separated (§4.6): half-open for the FM recurrence, closed-interval for bitvectors.
11. **Stream-specific compression** with declared codec choices for S_qual / S_meta / S_map (§4.7) — Objective 3.
12. **Forward compatibility with v2.0 approximate matching** via hooks H1–H4 (Contract §4.5; Architecture §4.6.7 + §6.7) without re-indexing — Objective 4.
13. **In-compressed-domain query** — no query operation accesses raw S or any of S_seq / S_qual at runtime (I14).
14. **Provenance-complete `.bsi`** suitable for clinical / production audit (`bamsi info`, `bamsi verify`; Contract §9).

---

## 13. Reference Implementation Contract (Objective 6)

The reference implementation is the canonical software realisation of this Architecture. Alternative implementations are encouraged and must conform to the same `.bsi` format, Contract semantics, and observable behaviour, but are not required to mirror the reference's source organisation.

### 13.0 Implementation Language Constraint (Contract §10.1)

The reference implementation MUST be written in **C++ (≥ ISO C++20)** OR **Rust (≥ 2021 edition)**. **No Python implementation of the reference is permitted.** The constraint is binding and is justified in Contract §10.1 (latency targets, clinical workflows, htslib interop, memory safety). Each project picks one of the two languages — the **C++ track** or the **Rust track** — and stays in it; mixing the two within a single reference build is not part of v1.0.

Both tracks expose:
- The CLI surface defined in §13.5 (mirroring Contract §10.2).
- The stable C ABI defined in Contract §10.3 — exposed via `extern "C"` shims (C++ track) or `#[no_mangle] pub extern "C"` functions (Rust track).
- Bit-identical `.bsi` files within the same track across operating systems and CPU architectures (validated in CI per §13.4 step 4). Cross-track byte-identity (a C++ build vs. a Rust build of the same input) is **not** required.

Consumer-language wrappers (Python, Go, Java, R, Julia, …) are user territory and bind against the C ABI. They are not part of the v1.0 reference.

### 13.1 Source Tree Layout

The reference implementation is **single-track**: a project chooses C++ **or** Rust at start-of-work and stays in that track. The two layouts share the same module decomposition (one directory per pipeline stage / runtime concern); only the file naming and build-system metadata differ.

#### 13.1a Layout — C++ Track

```
bamsi/
├── LICENSE                  (Apache 2.0)
├── NOTICE                   (third-party licences)
├── README.md
├── SECURITY.md
├── CHANGELOG.md
├── Dockerfile
├── CMakeLists.txt
├── conanfile.txt | vcpkg.json
├── src/
│   ├── ingest/              (Stage 1)
│   ├── ordering/            (Stage 2)
│   ├── seqbuilder/          (Stage 3)
│   ├── sais/                (Stage 4 — SA-IS + optional parallel variant)
│   ├── seq_encode/          (Stage 5a — S_seq codec)
│   ├── fm_index/            (Stage 5b — forward FM-index)
│   ├── fm_index_reverse/    (Stage 5b' — reverse FM-index when enable_bidirectional)
│   ├── stream_encode/       (Stage 6 — S_qual / S_meta / S_map codecs)
│   ├── windows/             (Stage 7)
│   ├── bitvectors/          (Stage 8)
│   ├── validation/          (Stage 9 — TIER 1 + TIER 2)
│   ├── seal/                (Stage 10)
│   ├── query/               (FM backward search, Locate, Regional* execution)
│   ├── mapping/             (M_ℓ implementation)
│   ├── format/              (.bsi reader/writer)
│   ├── cli/                 (bamsi command-line tool)
│   └── c_abi/               (extern "C" shim — stable C ABI for downstream bindings)
├── include/bamsi/
│   ├── bamsi.hpp            (C++ public header — namespace bamsi::)
│   └── bamsi.h              (C public header — stable C ABI)
├── tests/
│   ├── unit/                (per-module GoogleTest suites)
│   ├── integration/         (TIER 1 validation in CI)
│   └── synthetic/           (TIER 2 validation in CI)
├── benchmarks/              (Paper 1 reproducibility scripts)
├── workflows/               (Snakemake / Nextflow demonstration pipelines)
└── docs/
    ├── format.md
    ├── algorithms.md
    ├── cli.md
    ├── api.md               (C ABI reference + C++ namespace reference)
    ├── clinical.md
    └── tutorials/
        ├── 01_motif_counting.md
        ├── 02_region_query.md
        └── 03_quality_post_filter.md
```

#### 13.1b Layout — Rust Track

```
bamsi/
├── LICENSE                  (Apache 2.0)
├── NOTICE                   (third-party licences)
├── README.md
├── SECURITY.md
├── CHANGELOG.md
├── Dockerfile
├── Cargo.toml               (workspace)
├── crates/
│   ├── bamsi-core/          (library: ingest, ordering, seqbuilder, sais, seq_encode,
│   │                         fm_index, fm_index_reverse, stream_encode, windows,
│   │                         bitvectors, validation, seal, query, mapping, format)
│   ├── bamsi-cli/           (binary: `bamsi` command-line tool, depends on bamsi-core)
│   └── bamsi-c-abi/         (cdylib + staticlib: #[no_mangle] extern "C" shim
│                              over bamsi-core; exports the stable C ABI)
├── include/bamsi/
│   └── bamsi.h              (C public header — stable C ABI; hand-written or
│                              cbindgen-generated from bamsi-c-abi)
├── tests/                   (workspace-level integration + synthetic tests)
│   ├── integration/         (TIER 1)
│   └── synthetic/           (TIER 2)
├── benches/                 (criterion.rs Paper 1 benchmark scripts)
├── workflows/               (Snakemake / Nextflow demonstration pipelines)
└── docs/
    ├── format.md
    ├── algorithms.md
    ├── cli.md
    ├── api.md               (C ABI reference + Rust crate rustdoc reference)
    ├── clinical.md
    └── tutorials/
        ├── 01_motif_counting.md
        ├── 02_region_query.md
        └── 03_quality_post_filter.md
```

The `bamsi-core` crate's module layout mirrors the C++ track's `src/` subdirectory list one-to-one (one Rust module per pipeline stage), so the algorithm-level documentation (`docs/algorithms.md`) is shared between tracks without modification.

### 13.2 Build System

The reference implementation is one-of-two; pick at project inception, do not mix.

#### 13.2a C++ Track Build System

- **Build:** CMake ≥ 3.20, **ISO C++20** (or later). Compilers: GCC ≥ 11, Clang ≥ 14, MSVC ≥ 19.30 (best-effort).
- **Dependencies:** htslib (≥ 1.18) for BAM I/O; libsais (or equivalent SA-IS); sdsl-lite (or equivalent succinct library — wavelet trees, RRR bitvectors); zstd (≥ 1.5); xxHash (≥ 0.8); OpenSSL (≥ 3.0) or libsodium for SHA-256.
- **Dependency management:** vcpkg (manifest mode) or Conan; reference uses vcpkg.
- **Output:** static-linked `bamsi` CLI binary (avoids deployment-environment skew); `libbamsi.a` and `libbamsi.so` (or platform equivalents) exposing both the C++ namespace and the C ABI.

#### 13.2b Rust Track Build System

- **Build:** Cargo, **Rust 2021 edition** (or later), MSRV stated in `Cargo.toml` (recommended ≥ 1.75 for stable async + GATs). Single workspace with three crates as listed in §13.1b.
- **Dependencies (selected, all current at the time of writing — verify versions on each release):**
  - `hts-sys` or `noodles-bam` for BAM I/O (rust-bio or htslib FFI).
  - `libsais-sys` (FFI) or a pure-Rust SA-IS port for Stage 4.
  - **Succinct data structures.** Rust has no single library equivalent to sdsl-lite; the Rust track uses a combination of `sucds`, `vers-vecs`, `fid-rs`, or hand-rolled wavelet-tree / RRR bitvector implementations. Whichever combination is chosen, the validation suite (§9.3) MUST exercise rank/select identities in TIER 1.
  - `zstd` crate, `xxhash-rust`, `sha2` for hashing.
  - `clap` (≥ 4) for the CLI; `serde` + `serde_json` for `--json`; `rayon` for data-parallel sections.
- **Output:** static-linked `bamsi` CLI binary via `cargo build --release`; `bamsi-c-abi` produces both `cdylib` (`libbamsi.so`) and `staticlib` (`libbamsi.a`) with the C ABI; the `bamsi` Rust crate is publishable on `crates.io`.

### 13.3 Test Framework

#### 13.3a C++ Track

- **Unit tests:** GoogleTest, one suite per `src/` module.
- **Property-based tests:** rapidcheck (recommended) for invariant checks (rank/select identities, partial reconstruction).
- **TIER 1 integration:** every push, ≤ 5% of full build time. Round-trip on 10,000 sampled reads; FM on 100 fixed patterns; bitvector identities; checksums; GlobalExists/RegionalExists equivalence.
- **TIER 2 synthetic:** nightly on `|S| ≤ 10^6` synthetic inputs. Brute-force FM cross-check; exhaustive coverage; CIGAR edge cases; stream independence; partial reconstruction (Tables A + B); SARange / ISA equivalence; parallel SA-IS bit-identity; block-level directory equivalence; codec-ablation correctness.

#### 13.3b Rust Track

- **Unit tests:** Cargo's built-in `#[test]` harness, one test module per `bamsi-core` module.
- **Property-based tests:** `proptest` (recommended) for the same invariant checks.
- **TIER 1 / TIER 2 split:** identical to the C++ track in coverage and gating; implemented as `#[cfg(test)]` integration tests in the workspace `tests/` directory and run via `cargo test --release` (TIER 1) and `cargo test --release --features tier2-synthetic` (TIER 2 nightly).

### 13.4 CI Pipeline

GitHub Actions (or equivalent), per-track. For each commit:

1. **Build matrix:** Linux x86_64, Linux aarch64, macOS x86_64, macOS aarch64 (Windows is best-effort).
2. **Static analysis & sanitisers:**
   - **C++ track:** clang-tidy, IWYU, AddressSanitizer, UndefinedBehaviorSanitizer, ThreadSanitizer (TIER 2 only — too slow for every commit).
   - **Rust track:** `cargo clippy -- -D warnings`, `cargo fmt --check`, `cargo deny` (licences + advisories), `cargo miri test` for unsafe-code modules in TIER 2.
3. **Unit + TIER 1 integration tests.**
4. **Determinism test (cross-platform within track):** build the same `.bsi` from the same input on two runners (Linux x86_64 and macOS aarch64) and assert byte-identity. Validates I9 across platforms within the chosen track. Cross-track byte-identity (C++ build vs. Rust build) is **not** required and is **not** tested.
5. **Format-version compatibility:** verify the current binary can read previous-major-version `.bsi` files.

**Nightly:**
- TIER 2 tests.
- Benchmark smoke run (subset of Paper 1's full battery) to detect performance regressions.

### 13.5 CLI Subcommand Surface (Mirrors Contract §10.2)

| Subcommand | Purpose |
|---|---|
| `bamsi build` | Build `.bsi` from BAM input(s) |
| `bamsi count` | GlobalCount(P) |
| `bamsi exists` | GlobalExists(P) |
| `bamsi locate` | Locate(P) |
| `bamsi region-count` | RegionalCount(P, chrom, [a,b]) |
| `bamsi region-exists` | RegionalExists |
| `bamsi reconstruct` | Recover ℛ from streams (full or per-read-id subset) |
| `bamsi info` | Emit provenance and codec disclosure (Contract §9.2) |
| `bamsi verify` | Re-compute hashes / checksums; verify byte integrity |

Every subcommand supports `--json` for machine-parseable output. `--strand`, `--sort-output`, `--bed`, `--threads` available where applicable.

### 13.6 Documentation Deliverables

The `docs/` tree contains, at v1.0 release:

- **format.md** — `.bsi` byte-level format spec (mirrors §7).
- **algorithms.md** — algorithmic overview citing Contract §3, §4 and Architecture §4–§6.
- **cli.md** — full CLI reference with examples.
- **api.md** — Stable C ABI reference (the ground truth for all bindings) plus the chosen-track native-API reference: C++ namespace `bamsi::` from `bamsi.hpp` on the C++ track, Rust crate rustdoc on the Rust track.
- **clinical.md** — operational guidance: provenance verification (`bamsi verify`), audit-log integration, lossy-mode caveats, no-implicit-network policy.
- **tutorials/01_motif_counting.md** — end-to-end: download NA12878, build, count 10 motifs.
- **tutorials/02_region_query.md** — `RegionalCount` over `chr17:7570000-7580000` (TP53).
- **tutorials/03_quality_post_filter.md** — two-phase workflow: `Locate` → fetch qualities → filter at Q ≥ 30.

### 13.7 Distribution Channels

- **GitHub:** tagged SemVer releases (`v1.0.0`, …) with source tarball and signed checksums.
- **Docker Hub / GHCR:** `bamsi/bamsi:v1.0.0` and `bamsi/bamsi:latest`. Image SHA-256 in release notes. The image bundles the compiled CLI binary; track choice (C++ vs. Rust) is not visible to consumers.
- **Bioconda:** `conda install -c bioconda bamsi`. Same compiled-binary contract as the Docker image; track choice is not visible.
- **Track-native package index (Rust track only):** `crates.io` publication of the `bamsi`, `bamsi-cli`, and `bamsi-c-abi` crates. Not applicable to the C++ track (no first-party C++ registry; vcpkg / Conan recipes are user-territory).
- **Cloud-genomics workflow templates:** Snakemake and Nextflow recipes in `workflows/` demonstrating BAMSI alongside `samtools` / `bcftools` in typical pipelines.
- **Not distributed by the reference project:** PyPI, npm, CRAN, Maven Central. These are user-territory; downstream maintainers MAY publish bindings to those registries against the stable C ABI (Contract §10.3) without coupling them to the reference project.

---

*End of BAMSI Architecture v4.3*
*Aligned to: BAMSI Contract v3.3*
