# FINAL FORMAL PROBLEM AND DESIGN SPECIFICATION
## BAMSI — Lossless Compression and Exact Pattern Matching over BAM Data using Succinct Indexing
**Version:** 3.3 — Hostile-Reviewer-Hardened + Objective-Aligned + Language-Constrained  
**Status:** Contract-complete. Cross-aligned with Architecture v4.3. Resolves all 22 hostile-reviewer findings (v3.1), aligns the spec to the six BAMSI project objectives (v3.2), and binds the reference implementation language to **C++ (≥ C++20) or Rust (≥ 2021 edition)** with no Python implementation permitted (v3.3). Suitable for SIGMOD / VLDB / RECOMB / Bioinformatics Q1 submission and clinical / production deployment.

---

## Revision History

| Version | Changes |
|---|---|
| 1.0 | Original specification |
| 1.1 | CIGAR mapping fix, coordinate convention, window determinism |
| 2.0 | rank₁ off-by-one fix; S_meta/S_map lookup separation; T unit consistency; S_seq entropy pipeline; RegionalCount complexity; Exists(P) restored; RegionalExists added; strand policy explicit opt-in; ℛ inclusion rule formal; source_manifest_hash defined; partial reconstruction table; BWT lifecycle; per-section checksums |
| 3.0 | Cross-document audit with Architecture v3.3. Fixes: (1) Strand policy unified — strand-complete operational default; theorems per-orientation. (2) GlobalExists(P) added; (3) GlobalCount strand-summing semantics explicit. (4) T window unit re-affirmed. (5) Source manifest hash byte-level precise. (6) Query output ordering formally specified. (7) ISA samples/SARange normative optional. (8) BWT derivation formula restated. (9) Ordering hash byte encoding specified. (10) chrom_id derivation explicit. (11) Separator-must-include rule tightened. (12) Cross-reference inconsistencies with Architecture v3.3 resolved. |
| 3.1 | Hostile Q1 reviewer simulation fixes (22 findings across Theory, Systems, Bioinformatics panels). **Fatal fixes:** (F1) Space bound corrected to `\|S\|·H_k(S) + O(\|S\|)` — the FM-index BWT and metadata streams form a legitimate additive linear term, not `o(\|S\|)`. (F2) SARange complexity guarantees formally split into BASE tier (all implementations) and ENHANCED tier (SARange-enabled); SARange formally specified in Appendix A. (F3) Removed `S[pos] == '#'` defensive check from Locate/RegionalCount pseudocode — provably impossible for P ∈ Σ* and contradicted the "S never materialized" invariant. (F4) Per-read directories split: dir_meta/dir_map mandatory per-read; dir_seq/dir_qual permitted block-level (default block=1024 reads) to eliminate the 48 GB/1B-read overhead. (F5) Pipeline restructured so SA-IS is its own stage, producing a single BWT consumed by both S_seq encoding and FM-index construction (stages 5a/5b in parallel). (F6) Window complexity formula corrected for unit consistency — `\|W_r\| = O((L·d)/T)` where d is coverage and T is in S-characters. **Significant fixes:** (S1) Two distinct rank APIs explicitly specified (BWT-rank half-open for FM backward search; bitvector-rank closed-interval for B_read/B_window); removed the misleading "single consistent convention" clause. (S2) Entropy order k ∈ [4,8] justified — default k=6 with mandatory sensitivity analysis in Paper 1; k stored in S_seq header. (S3) SA-IS mandate relaxed to "SA-IS or a deterministic parallel variant with bit-identical output on the reference suite." (S4) Validation matrix split into TIER 1 (production builds, fast) and TIER 2 (CI/synthetic, exhaustive). (S5) Soft-clip/insertion interval interpretation note added with worked example. (S6) Strand-complete counting semantics disambiguated (read-sequence-level, not genomic-position deduplicated). (S7) Partial reconstruction table split into Table A (reconstruction) and Table B (query capability) — clarifies S_seq is never loaded at query time. (S8) chrom_id portability clarified — frozen per-index via stored chrom_name_table; resolved by string name in query APIs. (S9) Quality-aware queries positioned as v2.0 future work with two-phase workflow. **Clarifications:** (C1) Locate output ordering has two modes — streaming (default, no sort) and sorted (on `--sort-output` flag); both return the same set. (C2) ISA samples fully specified: construction during SA-IS, s' stored in header, forward walk via LF-mapping. (C3) Secondary/supplementary alignment scope explicitly communicated with Paper 1 empirical comparison. |
| 3.2 | Objective-alignment pass against the six BAMSI project objectives. (O1 — modular lossless compression, reference-free + reference-based, cross-platform/sequencer reproducibility) §0.1 expanded; §2.10 reproducibility/portability/provenance section added; sequencing-technology agnosticism stated explicitly. (O2 — in-compressed-domain pattern matching, research + clinical workflows) §0.1 explicitly states in-compressed-domain semantics; new §9 Operational Guarantees covers determinism, provenance surface, query-time reproducibility, no-implicit-network for clinical workflows. (O3 — stream-specific compression strategies) New §2.7 (S_qual: per-cycle range / rANS-delta / ZSTD-dict / binned codecs), §2.8 (S_meta: TYPED_SPLIT with CIGAR/FLAG/aux substreams), §2.9 (S_map: DELTA_RANGE codec exploiting sorted-pos structure). (O4 — exact AND approximate queries at population scale) New §4.5 Approximate Pattern Matching Extension Point — k-mismatch and k-edit promoted from flat non-goal to formal v2.0 extension with v1.0 forward-compatibility hooks (H1 bidirectional FM, H2 seed-length, H3 SA range API stability, H4 locate stability); v1.0 returns `NOT_IMPLEMENTED_V1` on approximate queries. (O5 — benchmark vs. industry tools on real datasets) §8 Paper 1 fully expanded with named tools (CRAM, Genozip, samtools, NGC), named real datasets (NA12878, HG002 HiFi, HG002 ONT, exome, cohort), per-(tool,dataset) metrics, sub-experiments. (O6 — open-source CLI / docs / tutorials) New §10 Distribution and Software Deliverables: Apache 2.0 licence requirement; **implementation language MUST be C++ (≥ C++20) or Rust (≥ 2021 edition); no Python implementation is permitted for the reference build** — the language constraint is binding because clinical-workflow latency targets (Contract §9, sub-50 ms median GlobalCount) require deterministic memory layout and no GC pauses; normative CLI subcommand surface (`bamsi build/count/exists/locate/region-count/region-exists/reconstruct/info/verify`); stable C ABI for cross-language consumption (other-language bindings — Rust, Go, Java, R — are user-territory but not part of the v1.0 reference); mandatory documentation deliverables (README, format spec, CLI ref, API ref, clinical guidance, three end-to-end tutorials); distribution channels (GitHub/Docker/Bioconda/Snakemake-Nextflow); versioning policy. |
| **3.3** | **This version.** Implementation-language constraint hardened. The v3.2 wording "C and Python bindings" is corrected to **"C++ (≥ C++20) or Rust (≥ 2021 edition) only"** for the reference implementation. Python is removed as a permitted implementation language but **NOT** as a permitted *consumer* language: the reference exposes a stable C ABI (extern "C" from C++ or `#[no_mangle] extern "C"` from Rust) which any FFI-capable language — Python, Go, Java, R, Julia, etc. — MAY bind to as a user-side wrapper outside the v1.0 reference scope. Rationale captured in §10.1 (language constraint with deterministic-layout / no-GC-pause / clinical-latency justification) and §10.3 (library deliverable surface restated as "C++ or Rust core + stable C ABI"). Architecture §13.2 build-system row updated to a one-of-two choice (CMake/C++20 OR Cargo/Rust-2021); §13.7 distribution channels drop PyPI; tutorial deliverables and clinical-workflow `bamsi info --json` output are unaffected because they go through the CLI, which is language-implementation-agnostic. |

---

## 0. Scope and Frozen Semantic Choices

This section establishes choices that are **frozen** — they must not vary across implementations or over the project lifetime. Any change requires a full version increment of this contract.

### 0.1 System Purpose

BAMSI is a system that:

- **Losslessly compresses** the **indexed aligned-read subset** of BAM files into four independent streams, in both **reference-free** (mandatory) and **reference-based** (optional) modes (§2.6).
- Builds a **succinct FM-index** over the concatenated read sequences, enabling **in-compressed-domain pattern matching** — exact search executes directly against compressed structures (FM-index, bitvectors, windows) without decompressing the full BAM (I14, §4.2).
- Supports **exact pattern matching** as the v1.0 query model, with a **formal extension point** for **approximate matching** (k-mismatch / k-edit) defined in §4.5 — the v1.0 index is forward-compatible with a v2.0 approximate-search overlay that does not require re-indexing.
- Provides **stream-specific compression strategies** tuned to the statistical structure of each data class: high-order entropy coding for sequence (§2.4), context-aware Q-score coding for quality (§2.7), CIGAR-run / FLAG-dictionary coding for metadata (§2.8), and delta+chrom-rank coding for mapping coordinates (§2.9).
- Is **scalable to population-genomics dataset sizes** (multi-billion reads per index) via parallel ingestion, parallel SA construction (§0.7), block-level reconstruction directories (§2.3), and ENHANCED-tier output-sensitive query bounds (§4.3).
- Provides formal **soundness, completeness, complexity, and reproducibility** guarantees suitable for Q1 publication.
- **Targets both research and clinical-discovery workflows** under documented provenance and reproducibility guarantees (§9 Operational Guarantees).
- Is **technology-agnostic** with respect to the sequencing platform: short-read (Illumina, MGI/BGI), long-read (PacBio HiFi, ONT) and any platform whose alignments produce SAM/BAM-compliant records satisfying the inclusion rule (§0.2). Compression-ratio behaviour varies by technology (Paper 1 §11 reports per-technology numbers); correctness is invariant.

### 0.2 Scope of ℛ — Inclusion Rule (Formal)

**Definition.** A BAM record belongs to the indexed set ℛ if and only if all of the following hold:

```
FLAG & 0x4   == 0    (record is mapped)
FLAG & 0x100 == 0    (record is not secondary)
FLAG & 0x800 == 0    (record is not supplementary)
POS >= 1             (SAM 1-based position is valid)
```

**Guarantee.** The losslessness guarantee `decode(encode(ℛ)) = ℛ` applies exactly to this set.  
**Non-goal.** Reconstruction of unmapped, secondary, or supplementary records is out of scope.

### 0.3 Alphabet Encoding — Frozen

| Symbol | Code |
|---|---|
| A | 0 |
| C | 1 |
| G | 2 |
| T | 3 |
| N | 4 |
| `#` (separator) | 5 |
| `$` (sentinel) | 6 — conceptual only, never stored |

Only codes `{0..5}` appear in stored S and stored index structures.

### 0.4 Separator Policy

- `#` (code 5) appears only between consecutive reads in S.
- Query patterns P never contain `#`.
- `#` is never counted as a match character.
- `#` belongs to the same window as the **preceding** read. This is a hard rule; no window boundary may fall inside a separator character or between a read's last base and its following separator.

### 0.5 Sentinel Policy

- `$` (code 6) is a conceptual sentinel appended to S to form S$ for suffix-array construction.
- `$` is lexicographically smaller than all other symbols.
- `$` is **not stored** in S or in the BWT payload.
- The FM-index behaves exactly as if `$` were present at position |S| of S.
- The **sentinel row** is the unique SA row where SA[row] = |S|. It is tracked explicitly as `sentinel_row` in the index header.
- No query pattern can contain `$`; the sentinel row never produces a reported match.

### 0.6 Strand Policy — Operational Default and Formal Basis

**Operational default (strand-complete):** All shipped implementations MUST operate in strand-complete mode by default. This is the deployment standard and the mode against which benchmarks, user-facing queries, and system comparisons are made.

**Definition of strand-complete mode:**

- For a query pattern P, define `Q(P) = {P}` if P = rc(P), else `Q(P) = {P, rc(P)}`.
- Each element of Q(P) is searched independently as a distinct pattern over S.
- Each result carries a `QueryStrand` label: `Forward` (from P) or `Reverse` (from rc(P)).
- If P = rc(P), search is performed exactly once and labeled `Forward`.
- `QueryStrand` denotes query orientation only; it is independent of the BAM FLAG strand field.

**Formal theorem basis (per-orientation):** All formal theorems and complexity bounds in this contract are stated **per orientation** — i.e., for a single pattern string Q ∈ Q(P). For strand-complete mode, each bound applies independently to each orientation; total cost is multiplied by `|Q(P)| ≤ 2`.

**Single-strand mode (legacy/optional):** An implementation MAY expose a `SingleStrand` build/query flag. When active, only P itself is searched and no `rc(P)` search is performed. All complexity bounds apply directly without the factor of 2.

**Reverse-complement definition (frozen):** Given sequence X = x₀x₁…x_{n-1} over {A,C,G,T,N}:

```
rc(X)[i] = complement(X[n-1-i])
complement: A↔T, C↔G, N↔N
```

**Counting semantics (normative — disambiguates strand-complete totals):**

BAMSI counts occurrences at the **read-sequence level**, not at the deduplicated-genomic-locus level:

- Each occurrence of `P` inside any `r_i` contributes `1` hit to the `Forward` tally.
- Each occurrence of `rc(P)` inside any `r_i` contributes `1` hit to the `Reverse` tally.
- Two matches on the same genomic locus arising from opposite orientations (e.g., forward match in a read and reverse-complement match in a mate read aligned on the opposite strand) are counted **as 2**, not as 1.
- If `P = rc(P)` (palindrome), the pattern is searched exactly once and only the `Forward` tally is populated; the genomic locus is counted once even if it appears on both strands of the same read.
- This convention matches `samtools view -c` on `SEQ`-filtered reads under sequence-exact matching.

Users requiring **deduplicated genomic-position counts** must post-process `Locate` results by grouping on `(chrom_id, p_min, p_max)` or by read_id as appropriate for their analysis. This is a deliberate scope choice: BAMSI's primary semantic is "how many read-substring occurrences exist," not "how many distinct genomic positions are covered." The chosen semantic preserves exact compositional summability: `GlobalCount(P) = Σ_i occurrences_in_r_i(P) + Σ_i occurrences_in_r_i(rc(P))`.

### 0.7 Determinism Guarantee

Given the same BAM inputs and the same `BuildConfig`, BAMSI MUST produce a bit-identical `.bsi` file on every run. Determinism requires:

- Fixed ordering rule (§1.2)
- Stable tie-breaking via `(source_file_id, bam_offset)` (§1.2)
- Fixed suffix-array algorithm: **SA-IS** (Nong, Zhang & Chan, 2009) is the reference algorithm. A parallel variant MAY be used in implementations provided it produces **bit-identical SA output to the reference single-threaded SA-IS** on every input in the fixed validation suite (§8 Paper 1 / §9.3 of the Architecture). The reference algorithm is used for validation; deployed builds may use a conforming parallel variant for throughput. No non-deterministic SA algorithm is permitted under any circumstance.
- Fixed compression parameters: `k`, coder type, serialization order
- Fixed little-endian byte serialization throughout
- Fixed sampling step `s` (and optionally `s'` for ISA samples)
- Fixed default `T` and `k` values as specified herein

### 0.8 Query Semantics — Occurrences Not Reads

- Queries operate over **occurrences**, not distinct reads.
- Multiple matches of P in the same read are counted separately.
- A pattern P of length 0 is invalid; implementations MUST return a structured error (`EMPTY_PATTERN`).
- A pattern P containing `#`, `$`, or any character outside Σ is invalid; implementations MUST return a structured error (`INVALID_PATTERN`).
- v1.0 query operations are **exact-match only**. Any approximate-match query (k-mismatch, k-edit) MUST return a structured `NOT_IMPLEMENTED_V1` error rather than fall back to exact matching, so that callers cannot silently misinterpret v1.0 results as approximate (§4.5).

### 0.9 Multi-BAM Policy

- BAMSI may ingest one or more BAM files.
- The order of BAM files in the build configuration is frozen as part of the input manifest.
- Each read's stable ordering key includes `source_file_id` (zero-based index of its BAM file in the input list) and `bam_offset` (stable record index within that file).
- `source_file_id` and `bam_offset` together uniquely identify any record across all inputs.

---

## 1. Data Model

### 1.1 Aligned Read Collection

$$\mathcal{R} = \{(r_i, c_i, p_i, m_i, Q_i) \mid i = 0, 1, \ldots, N-1\}$$

Where:

- $r_i \in \Sigma^*$, $\Sigma = \{A, C, G, T, N\}$; codes 0–4 respectively
- $c_i$: reference identifier — chromosome name resolved to a string from the BAM header
- $p_i$: leftmost **aligned** genomic coordinate, SAM POS, **1-based** (converted from BAM's 0-based at ingestion — exactly once at ingestion time; never re-converted)
- $m_i$: metadata record containing CIGAR ops, FLAG, and optional tags
- $Q_i$: Phred score array, $|Q_i| = |r_i|$, values in $[0, 93]$

**Scope:** Only records satisfying the inclusion rule (§0.2) are in ℛ. The index i runs from 0 to N−1 throughout this document.

### 1.2 Ordering, Read Identity, and Hashes

**Total order ≺ over ℛ:**

$$r_i \prec r_j \iff (c_i < c_j) \lor (c_i = c_j \land p_i < p_j)$$

where `<` on chromosome names is defined by `chrom_id` (see below).

**chrom_id derivation (frozen per index):**  
`chrom_id(name)` = the 0-based rank of `name` in the lexicographically sorted list of all distinct reference sequence names appearing in the BAM header(s). This sort is performed once at header-parse time and **frozen** for the lifetime of the index.

The resulting `(chrom_id → chrom_name)` mapping is stored verbatim in the `.bsi` header as the `chrom_name_table`. All query APIs accept **string chromosome names** (e.g., `"chr10"`); the implementation resolves each name to its stored `chrom_id` internally via this table. The lexicographic ordering is used for determinism of the build only and is transparent to query callers — for example, a query `RegionalCount(P, "chr10", …)` is resolved identically regardless of whether, under lexicographic ordering, `"chr10"` sorts before `"chr2"`.

**Portability across multi-BAM builds:** When multiple BAM inputs are combined, the sort is computed over the **union** of all distinct reference names across all input headers. The resulting `chrom_name_table` is stored in the index; any two builds consuming the same `(file set, BuildConfig)` tuple produce identical `chrom_name_table` contents (enforced by I9 determinism and `ordering_hash` / `source_manifest_hash`). Indices built from BAMs with different `@SQ` orderings are **not interchangeable as raw `chrom_id`s**, but their query semantics are identical when queries go through string chromosome names, which is the only supported query interface.

**Tie-breaking (frozen):**

- For single-BAM input: if $(c_i, p_i) = (c_j, p_j)$, order by `bam_offset` (record index in file, 0-based).
- For multi-BAM input: if $(c_i, p_i) = (c_j, p_j)$, order by `(source_file_id, bam_offset)`.

The sort key is the 4-tuple: `(chrom_id, pos, source_file_id, bam_offset)`.

**Definition:**

$$\text{read\_id}(r_i) = \text{0-based rank of } r_i \text{ under } \prec$$

`read_id` is a bijection from ℛ to $\{0, 1, \ldots, N-1\}$.

**Ordering Hash (byte-level specification):**

$$\text{ordering\_hash} = \text{SHA-256}\!\left(\text{little-endian encoding of sorted sequence of } (\text{chrom\_id}_i, p_i, \text{source\_file\_id}_i, \text{bam\_offset}_i)\right)$$

Each 4-tuple is encoded as: `uint32 chrom_id || uint64 pos || uint32 source_file_id || uint64 bam_offset`, all little-endian, concatenated in ascending `read_id` order with no separators. The SHA-256 is computed over these raw bytes.

**Source Manifest Hash (byte-level specification):**

$$\text{source\_manifest\_hash} = \text{SHA-256}\!\left(\text{concatenation of } (\text{len}(f_f) \mathbin\| f_f \mathbin\| \text{SHA-256(BAM\_header\_bytes}_f))_{f=0..F-1}\right)$$

where:
- Files appear in ascending `source_file_id` order.
- $f_f$ = UTF-8 encoded filename bytes of file $f$.
- $\text{len}(f_f)$ = `uint32` little-endian byte count of the filename.
- `BAM_header_bytes`$_f$ = the raw bytes of the BAM header block before the first record of file $f$.
- The SHA-256 of each BAM header is computed as a 32-byte digest.
- The full source_manifest_hash is SHA-256 over the concatenation of all per-file records.

Both hashes are stored in the `.bsi` header and re-verified on every index open. A mismatch MUST cause an immediate, non-silent error (`MANIFEST_MISMATCH` or `ORDERING_MISMATCH`).

### 1.3 Concatenated Sequence S

$$S = r_0 \,\#\, r_1 \,\#\, \ldots \,\#\, r_{N-1}$$

- `#` is code 5; it does not belong to $\Sigma$
- Patterns $P \in \Sigma^*$ never include `#`
- $|S| = \sum_{i=0}^{N-1} |r_i| + (N-1)$

**readStarts:** $\text{readStarts}[i]$ is the 0-based position in S of the first character of $r_i$.

$$\text{readStarts}[0] = 0$$
$$\text{readStarts}[i] = \text{readStarts}[i-1] + |r_{i-1}| + 1 \quad \text{for } i \geq 1$$

**Implications:**

- No exact match crosses a read boundary (enforced by §0.4 and the separator's absence from Σ)
- Each match belongs to exactly one read
- `read_id` uniquely identifies the read containing any given match

### 1.4 CIGAR Mapping Definition (Formal and Total)

#### 1.4.1 Matching Rule

- Pattern matching is performed **only on raw read sequences** $r_i$.
- CIGAR is **ignored during matching**.
- CIGAR is applied **after** a match position is found, to map read offsets to genomic coordinates.

#### 1.4.2 Reference Position Function

Define `cigar_ref_pos(cigar, p_anchor, offset, direction)` where:

- `cigar`: ordered list of `(op, len)` pairs for read $r_i$
- `p_anchor`: $p_i$ (1-based leftmost aligned position)
- `offset`: 0-based index within $r_i$, $0 \leq \text{offset} \leq |r_i| - 1$
- `direction`: LEFT or RIGHT

The function traverses CIGAR ops maintaining `ref_pos` (starting at `p_anchor`) and `read_pos` (starting at 0):

| CIGAR Op | Consumes Read | Consumes Ref | Rule |
|---|---|---|---|
| M, =, X | Yes | Yes | If `offset ∈ [read_pos, read_pos+len-1]`: return `ref_pos + (offset − read_pos)`. Else advance both by `len`. |
| I | Yes | No | If `offset ∈ [read_pos, read_pos+len-1]`: LEFT → last aligned ref base before this insertion (or `p_anchor` if none); RIGHT → first aligned ref base after this insertion (or last aligned base before if no base exists after; or `p_anchor` if no aligned base exists at all). Else advance `read_pos` by `len` only. |
| D, N | No | Yes | Advance `ref_pos` by `len` only. `read_pos` unchanged. |
| S | Yes | No | If `offset ∈ [read_pos, read_pos+len-1]`: map to nearest aligned reference base. If two are equidistant, choose the smaller coordinate. If no aligned base exists at all (no M/=/X/D/N ops in the entire read), return `p_anchor`. Else advance `read_pos` by `len` only. |
| H, P | No | No | Ignore entirely. Neither `ref_pos` nor `read_pos` changes. |

The function is **total**: it returns a defined value for every valid `(cigar, p_anchor, offset, direction)` tuple.

**Equidistance tie-break for S (soft clip):** always choose the smaller coordinate.

**No-aligned-base fallback:** if the read has no M/=/X/D/N operations at all (entire read is soft/hard-clipped or padding), all offset mappings return `p_anchor`.

#### 1.4.3 Interval Mapping

For a match at read interval $[x, y]$ (0-based, inclusive, within $r_i$):

$$p_{\min} = \text{cigar\_ref\_pos}(\text{cigar}_i, p_i, x, \text{LEFT})$$
$$p_{\max} = \text{cigar\_ref\_pos}(\text{cigar}_i, p_i, y, \text{RIGHT})$$

**Properties (invariant):**

- $p_{\min} \geq 1$ always (1-based coordinate)
- $p_{\min} \leq p_{\max}$ always
- The interval $[p_{\min}, p_{\max}]$ is the **minimal enclosing genomic interval**
- Mapping is total and deterministic: identical inputs always produce identical outputs
- No match is ever discarded due to insertions, soft clips, or any CIGAR structure

**Interpretation note (normative for users):** `[p_min, p_max]` is the **minimal enclosing reference interval** for the match — it is guaranteed to contain the reference positions of all aligned bases touched by the match, and for matches that touch no aligned base (e.g., entirely within a soft clip or insertion), it collapses via the rules in §1.4.2 to a well-defined but coarser interval. Specifically:

- A pattern of length $\ell$ matched entirely inside a soft-clipped region of length $C$ can map to an interval whose width is up to $C$ reference bases. The reported interval is the tightest interval that provably encloses the genomic positions of any aligned base adjacent to the clipped region; it is **not** claimed to be the precise sub-base position of the pattern within the clip.
- A pattern straddling an insertion is mapped so that $p_{\min}$ is the last aligned reference base before the insertion (or the insertion anchor) and $p_{\max}$ is the first aligned reference base after it; the reported width reflects the reference span that the insertion occupies (zero) plus any flanking aligned bases.
- A pattern entirely inside a deletion or reference-skip does not arise — these ops do not consume read bases, so there is no read position to match.

**Worked example.** Consider read $r_i$ with CIGAR `10S150M5S`, `p_i = 1000`, `|r_i| = 165`. Let pattern `P` of length 12 match at read offsets $[x, y] = [3, 14]$ — i.e., the first 7 bases are inside the leading 10S clip and the last 5 bases are inside the 150M aligned region:

- `x = 3` falls inside the S op. Per §1.4.2, LEFT maps to the nearest aligned ref base. The nearest aligned M/= /X op starts at read offset 10 with ref_pos 1000. Equidistance tie-break chooses the smaller coordinate, yielding $p_{\min} = 1000$.
- `y = 14` falls inside the 150M op at (offset − read_pos of M) = (14 − 10) = 4. Thus $p_{\max} = 1000 + 4 = 1004$.
- Reported interval: `[1000, 1004]`, width 5 ref bases for a 12-base pattern.

For the extreme case where `P` matches entirely within the 10S clip at read offsets `[0, 9]`: both LEFT and RIGHT return 1000 (nearest aligned base = the first M base), yielding `[1000, 1000]` — a 1-base interval for a 10-base pattern.

**Users requiring precise sub-read positioning** of clipped-region matches should post-process using the returned `read_id`, retrieve the full CIGAR from `S_meta[read_id]`, and interpret the clipped region against the pattern offsets independently. BAMSI's interval is guaranteed **sound** (the true reference positions of aligned bases within the match fall inside `[p_min, p_max]`) and **total** (no match is discarded) but it is not claimed to be biologically tight for patterns dominated by soft-clipped or inserted bases.

#### 1.4.4 Coordinate Convention

$$[p_{\min}, p_{\max}] \text{ is 1-based, inclusive (SAM/BAM convention throughout)}$$

---

## 2. Compression Model

### 2.1 Losslessness

$$\text{decode}(S_{\text{seq}}, S_{\text{qual}}, S_{\text{meta}}, S_{\text{map}}) = \mathcal{R}$$

This guarantee applies to all and only the records in ℛ (primary mapped reads satisfying §0.2). Excluded records are outside scope.

### 2.2 Stream Decomposition (Strong Independence — Frozen)

$$\text{encode}(\mathcal{R}) = (S_{\text{seq}},\; S_{\text{qual}},\; S_{\text{meta}},\; S_{\text{map}})$$

**Stream content (normative):**

| Stream | Content |
|---|---|
| $S_{\text{seq}}$ | Compressed concatenated read sequences S |
| $S_{\text{qual}}$ | Compressed Phred quality score arrays $Q_i$ |
| $S_{\text{meta}}$ | Compressed per-read metadata: CIGAR ops, FLAG, and optional BAM tags |
| $S_{\text{map}}$ | Compressed per-read mapping data: chromosome name (via `chrom_id`) and 1-based POS $p_i$ |

**Strong independence rule (frozen):**

- No entropy model for one stream may depend on another stream's bytes or run-time statistics.
- Any cross-stream statistics that a stream needs MUST be **stored explicitly within that stream's own header** and decoded from there.
- Each stream is **independently decodable** from its own bytes alone.
- Streams may be decoded in any order or individually.
- This rule is invariant I10 and is directly required for Theorem 1 (Soundness).

### 2.3 Partial Reconstruction and Query Streams (Disambiguated)

The following two tables formally separate **reconstruction** (recovering original ℛ data) from **query execution** (answering GlobalCount / Locate / RegionalCount). These are **distinct access patterns** using disjoint subsets of the sealed `.bsi` file; a user must not conflate them.

#### Table A — Reconstruction (recovering original ℛ data)

Given only the listed streams and their embedded headers and directories, the following data is recoverable from the sealed `.bsi`:

| Streams Available | Data Recoverable |
|---|---|
| $S_{\text{seq}}$ only | Read sequences $r_i$ for all $i \in \{0..N-1\}$ |
| $S_{\text{map}}$ only | Reference identifiers $c_i$ (via stored `chrom_name_table`) and positions $p_i$ for all $i$ |
| $S_{\text{seq}} + S_{\text{map}}$ | Read sequences + genomic coordinates — genomically placed reads |
| $S_{\text{seq}} + S_{\text{qual}} + S_{\text{map}}$ | Above + quality scores $Q_i$ |
| $S_{\text{seq}} + S_{\text{meta}} + S_{\text{map}}$ | Read sequences + CIGAR + FLAG + tags + coordinates (without qualities) |
| All four streams | Full ℛ reconstruction — lossless BAM record set |

Each row is achievable using only the listed streams and their embedded headers, with no reference to any other stream. This property is a direct consequence of the strong independence rule (§2.2).

#### Table B — Query Capability (what is loaded at query time)

| Query Operation | Required Structures |
|---|---|
| `GlobalCount(P)`, `GlobalExists(P)` | FM-index (BWT, C, Occ, SA_samples) + `chrom_name_table` |
| `Locate(P)` | FM-index + B_read + Directories{meta, map} + S_meta (per-read lazy) + S_map (per-read lazy) + `chrom_name_table` |
| `RegionalCount(…)`, `RegionalExists(…)` | FM-index + B_read + B_window + WindowTable + Directories{meta, map} + S_meta (per-read lazy) + S_map (per-read lazy) + `chrom_name_table` + (optionally SARange) |
| Any query | **Never** requires $S_{\text{seq}}$, $S_{\text{qual}}$, raw $S$, or the uncompressed SA |

**Key distinction:** The FM-index is built from $S$ and **replaces** $S$ for all query-time operations. $S_{\text{seq}}$ is a reconstruction-oriented stream (Table A), not a query-oriented stream (Table B). An implementation MUST NOT load $S_{\text{seq}}$ to execute any query; attempting to do so violates the "raw $S$ never resident at query time" invariant (§8 of Architecture).

**Per-read vs block-level directory granularity (normative):**

- `dir_meta` and `dir_map` MUST be stored at **per-read granularity** (one entry per `read_id`). They are the hot path for `Locate` / `RegionalCount` mapping and require O(1) seek per match.
- `dir_seq` and `dir_qual` MAY be stored at either per-read granularity or **block-level granularity** (one entry per block of `B_dir` consecutive reads, default `B_dir = 1024`). Block granularity is permitted because $S_{\text{seq}}$ and $S_{\text{qual}}$ are reconstruction-only streams (Table A) and tolerate O(B_dir) intra-block decode cost per read access. Implementations choosing block-level directories MUST store `B_dir` in the corresponding stream header so that decoders can locate and decompress the containing block.
- At 1 billion reads, per-read dirs for meta+map cost 24 GB; with block-level seq/qual (B_dir=1024) the seq+qual dirs drop from 24 GB to approximately 23 MB each. This is deliberately absorbed into the overall `O(|S|)` space overhead.

### 2.4 Entropy Model

**Theoretical guarantee:**

$$|S_{\text{seq}}| = |S| \cdot H_k(S) + o(|S|), \quad k = O(1)$$

where $H_k(S)$ is the $k$-th order empirical entropy of S.

**Mandatory encoding pipeline for $S_{\text{seq}}$ (normative):**

1. **BWT:** Compute the Burrows-Wheeler Transform of S (reusing the SA-IS output from FM-index construction — see §2.5). BWT[i] = S[(SA[i] − 1 + |S|) mod |S|] for non-sentinel rows.
2. **MTF:** Apply Move-to-Front transform to the BWT output.
3. **RLE:** Apply Run-Length Encoding to the MTF output.
4. **Arithmetic coding:** Encode the RLE output with a 0th-order arithmetic model.

This four-step pipeline achieves $H_k(S)$ compression for fixed $k$: the BWT clusters $k$-context runs, enabling the 0th-order arithmetic coder to exploit $k$-th order redundancy. Steps 2–4 are applied to the BWT payload only (the sentinel row's BWT character, which is conceptually `$`, is never stored).

**Entropy order parameter:** $k \in [4, 8]$, fixed at build time and stored in the $S_{\text{seq}}$ stream header (so decoders do not need to guess). Default: $k = 6$.

**Justification for the default $k = 6$ (normative guidance):**

The theoretical space bound $|S|H_k(S) + O(|S|)$ holds for **any fixed $k$**. The default $k=6$ is chosen empirically as a trade-off between compression ratio and encoder complexity on representative genomic datasets; it is **not** claimed that $k=6$ is optimal in any theoretical sense. The contract's space bound is parameterised by $k$ and carries over to any fixed $k \in [4,8]$.

Implementations MUST expose `entropy_order_k` as a `BuildConfig` field, and Paper 1 (Systems) MUST include a sensitivity analysis:

- Compression ratio vs. $k \in \{4, 5, 6, 7, 8\}$ for at least three representative datasets: short-read WGS BAM, exome BAM, and long-read (ONT/PacBio) BAM.
- Encoder/decoder time as a function of $k$.
- A brief justification of the chosen default.

This empirical section exists precisely to pre-empt reviewer objections that $k=6$ is arbitrary. The theory paper (Paper 2) states the bound parametrically and does not commit to any particular $k$.

### 2.5 BWT Lifecycle

The suffix array (SA) over S is computed **exactly once** during the build pipeline using SA-IS. This single SA computation serves two purposes:

1. Derive BWT for $S_{\text{seq}}$ compression (§2.4): `BWT[i] = S[(SA[i] − 1 + |S|) mod |S|]`
2. Construct the FM-index (§3.1)

After both consumers have extracted their outputs, the SA is discarded. S itself is discarded after the FM-index is built. No downstream component requires S or the SA as uncompressed contiguous arrays.

### 2.6 Reference Model

**Mandatory mode:** Reference-free. The FM-index and all query operations are reference-independent.

**Optional mode:** An implementation MAY use a reference genome to improve $S_{\text{seq}}$ encoding (e.g., encoding S as reference-relative differences). If so:

- S itself MUST NOT be modified. The FM-index is always built over the original unreferenced S.
- Query correctness is unaffected.
- The $S_{\text{seq}}$ stream header MUST declare whether reference-based encoding was used and, if so, identify the reference (e.g., by SHA-256 of the reference FASTA).

### 2.7 S_qual — Quality-Score Stream Strategy (Normative)

Quality scores have well-known statistical structure that generic codecs do not fully exploit. To meet Objective 3 (stream-specific compression), `S_qual` MUST use a coder that captures **at least one** of the following redundancy sources, declared in the stream header (`qual_codec_id`):

| Codec ID | Strategy | Captured redundancy |
|---|---|---|
| `0x01 — RANGE_CYCLE` | Range coder with a **per-cycle context** (read position modulo expected read length) | Sequencing-cycle quality decay (universal in Illumina/MGI; less pronounced but present in long-read) |
| `0x02 — RANS_DELTA` | rANS over the delta-encoded Q-score sequence (`Q_i[j] − Q_i[j-1]`) | Smooth quality transitions within reads |
| `0x03 — ZSTD_DICT` | ZSTD with a per-dataset dictionary trained on a sampled subset of $Q_i$ | Long-read flat-quality runs and platform-specific quality histograms |
| `0x04 — BINNED_RANGE` | Range coder with quantised quality bins (Illumina-2-bin / 8-bin / Q-Bin) — **lossy when quantised**; lossless when bin count = 94 | High-redundancy of binned qualities; provides a **declared lossy fast-path** that MUST be opted into via `BuildConfig.qual_lossy_bins` |

**Defaults (lossless mode):** `qual_codec_id = 0x01 (RANGE_CYCLE)`. Implementations MAY support multiple codecs; the choice is per-build, frozen at seal time.

**Lossy quality binning is opt-in only.** When `qual_lossy_bins ≠ 0`, the `.bsi` header field `is_lossless` is set to 0 and the losslessness guarantee of §2.1 is restricted to {sequence, mapping, metadata}; quality scores are reconstructed in their binned form. This is reported prominently on `bamsi info`.

**Independence preserved:** Whichever codec is chosen, S_qual remains independently decodable from its own bytes and embedded header (I10).

### 2.8 S_meta — Metadata Stream Strategy (Normative)

CIGAR strings, FLAG fields, and BAM aux tags have specific structure: CIGAR is a short run of `(op, len)` pairs with a small alphabet of ops (typically {M, I, D, S}); FLAG is a 16-bit bitmask with a heavily skewed distribution (a few values dominate); aux tags are sparse and key-repeated.

S_meta MUST exploit this via a **typed-substream** layout, declared in the stream header (`meta_codec_id`):

| Codec ID | Strategy |
|---|---|
| `0x01 — TYPED_SPLIT` (default) | Split per-read metadata into three logical substreams encoded independently within S_meta: (a) CIGAR substream — variable-length encoding of `(op, len)` pairs with op encoded in a 4-bit nybble (alphabet of 9 ops) and len encoded as a 7-bit-stop varint; (b) FLAG substream — range coder with adaptive frequency over the 16 commonly-set FLAG bits; (c) aux-tag substream — key-dictionary + ZSTD over per-tag values, skipped entirely when no aux tags present. |
| `0x02 — ZSTD_FALLBACK` | ZSTD over the concatenated raw metadata. Permitted as a baseline; MUST NOT be used for benchmark-headline numbers. |

**Default:** `meta_codec_id = 0x01 (TYPED_SPLIT)`. Header records this choice; substream offsets within the S_meta block are stored in the per-read directory (`dir_meta`).

### 2.9 S_map — Mapping Stream Strategy (Normative)

Mapping coordinates are highly compressible because reads are sorted by `(chrom_id, pos)` (§1.2). The first read of each chromosome carries an absolute `pos`; subsequent reads in the same chromosome carry small `Δpos` values clustering around the mean inter-read distance.

S_map MUST exploit this via:

| Codec ID | Strategy |
|---|---|
| `0x01 — DELTA_RANGE` (default) | (a) `chrom_id` substream — run-length encoded (long runs of identical chrom_id); (b) `pos` substream — first read of each chromosome stores absolute `pos`; subsequent reads store `Δpos = pos_i − pos_{i-1}` encoded with a range coder over a small adaptive alphabet of common deltas. |
| `0x02 — RAW` | Raw `(uint32 chrom_id, uint64 pos)` per read. Permitted only as a debug / validation baseline. |

**Default:** `map_codec_id = 0x01 (DELTA_RANGE)`. The header stores the codec ID; chunk boundaries (one chunk per chromosome) are recorded in `dir_map` so per-read O(1) lookup is preserved (the chunk's prefix-sum of deltas is reconstructed lazily for the queried block, costing O(B_dir) for block-level access — the per-read mandate of `dir_map` makes this O(1) per matched read).

### 2.10 Reproducibility, Portability, and Provenance (Normative)

To support reproducibility across platforms (Objective 1) and clinical-discovery workflows (Objective 2), every `.bsi` MUST carry sufficient metadata to be re-verifiable:

- **Build provenance.** The `.bsi` header stores: `bamsi_version` (string), `build_timestamp_utc` (uint64, optional — may be zero for fully-deterministic re-builds), `host_os_id` (uint8: 0=any/unspecified, 1=Linux, 2=macOS, 3=Windows), `cpu_arch_id` (uint8: 0=any, 1=x86_64, 2=aarch64). These are informational; the determinism guarantee (I9) means the same `BuildConfig` produces a bit-identical `.bsi` on any conformant platform. Re-builds verify byte-identity, not platform identity.
- **Source manifest.** `source_manifest_hash` (§1.2) fingerprints all input BAM files by header bytes; mismatches at index-open are non-silent errors (`MANIFEST_MISMATCH`).
- **Reference dependency disclosure.** When `reference_based_encoding = 1` (§2.6), the SHA-256 of the reference FASTA MUST be stored in the `S_seq` header. An implementation opening such a `.bsi` without access to the same reference MAY still answer queries (queries are reference-free) but MUST NOT attempt full reconstruction.
- **Codec disclosure.** All codec IDs (`qual_codec_id`, `meta_codec_id`, `map_codec_id`, `entropy_order_k`, etc.) are stored in stream headers and surfaced by the `bamsi info` command (§9).
- **Lossless / lossy flag.** `is_lossless` (uint8) is `1` iff every stream uses a lossless codec on the v1.0 ℛ scope. Any lossy choice (e.g., `qual_lossy_bins ≠ 0`) sets this to `0` and MUST be reported on every query response and reconstruction operation.
- **Cross-platform endianness.** All multi-byte integers are little-endian throughout (§0.7). On big-endian hosts, decoders byte-swap on read.

---

## 3. Indexing Model

### 3.1 Core Index

- FM-index built over S using alphabet $\sigma = 6$: `{A, C, G, T, N, #}` (codes 0–5)
- Rank/select structures: $O(1)$ amortized per operation
- S is **not materialized** at query time

### 3.2 FM-Index Formalization

- SA is defined over $S\$$ (S with conceptual sentinel appended at position |S|)
- `$` is lexicographically smaller than all codes 0–5
- **BWT derivation formula:** `BWT[i] = S[(SA[i] − 1 + |S|) mod |S|]` for all rows $i$ where $SA[i] \neq |S|$
- The **sentinel row** is the unique row $r^*$ where $SA[r^*] = |S|$. It is stored as `sentinel_row` in the index header.
- The stored BWT payload omits the sentinel row character (`$` is never stored). The sentinel row is tracked by its index `sentinel_row` and is handled as a special case in all algorithms.
- All query algorithms operate over logical SA rows; they skip `sentinel_row` explicitly and never report it as a match.

**C array:** $C[a]$ = number of SA rows whose BWT suffix starts with a character lexicographically less than $a$. Defined over the extended alphabet including `$`, though `$` never appears in queries.

**Backward search (exact, per orientation):**

```
backward_search(Q[0..m-1]):
    lo ← 0;  hi ← |S| + 1   // row space is [0, |S|] for S$ (|S|+1 rows total)
    for i = m−1 downto 0:
        a ← encode(Q[i])
        lo ← C[a] + rank_a(BWT, lo)
        hi ← C[a] + rank_a(BWT, hi)
        if lo ≥ hi: return ∅
    return [lo, hi)
```

Here `rank_a(BWT, pos)` counts occurrences of symbol `a` in BWT[0..pos−1] — **half-open** (classical Ferragina-Manzini convention, used only within the FM recurrence and LF-mapping). This is a different API from bitvector rank (§3.3).

Each step is $O(1)$ via the Occ structure. Total cost per pattern: $O(|Q|)$.

### 3.3 Bitvectors and the Two Rank Conventions (Normative)

BAMSI uses **two distinct rank APIs** with different semantics. They MUST NOT be implemented as a single function with a shared signature, and implementers MUST NOT interchange them.

| API | Semantics | Used by | Return value at position $p$ |
|---|---|---|---|
| `rank_a(BWT, pos)` | **Half-open** | FM-index backward search and LF-mapping (§3.2 only) | number of symbol `a` in `BWT[0..pos−1]` |
| `rank1(B, pos)` | **Closed-interval** | `B_read`, `B_window`, and all mapping/query code outside the FM inner loop | number of 1-bits in `B[0..pos]` inclusive |

**Lemma (equivalence up to boundary shift).** At any position $p$, the two rank conventions differ by at most 1:
$$\text{rank1}_{\text{closed}}(B, p) = \text{rank1}_{\text{half-open}}(B, p+1)$$
and the choice of closed-interval rank for bitvectors is what makes the read-identification formula `read_id = rank1(B_read, pos) − 1` exact. The `−1` is an intentional, provably correct off-by-one adjustment that arises from this choice of convention; it is not an error. An implementation that applies the half-open convention to `B_read` will produce incorrect `read_id` values.

The earlier "single consistent convention" language (v3.0) is retracted; it was a misleading simplification. The two APIs **are** consistent within their own scope, and explicitly differ at the boundary between scopes.

**Bitvector construction:**

- $B_{\text{read}}$: length $|S|$, $B_{\text{read}}[\text{readStarts}[i]] = 1$ for all $i \in \{0..N-1\}$, and 0 elsewhere.
- $B_{\text{window}}$: length $|S|$, $B_{\text{window}}[l_j] = 1$ for all windows $j$ (where $l_j$ is the S-position of the window's first character), and 0 elsewhere.

Both bitvectors support $O(1)$ rank and select.

**Rank convention (normative, closed-interval):**

$$\text{rank}_1(B, \text{pos}) = \text{number of 1-bits in } B[0..\text{pos}] \text{ (both endpoints inclusive)}$$

**Select convention (normative, 1-based rank):**

$$\text{select}_1(B, k) = \text{position of the } k\text{-th 1-bit in } B, \quad k \geq 1$$

**Key identities (required by invariant I8):**

- $\text{rank}_1(B_{\text{read}}, \text{readStarts}[i]) = i + 1$
- $\text{select}_1(B_{\text{read}}, i + 1) = \text{readStarts}[i]$
- $\text{rank}_1(B_{\text{window}}, l_j) = j + 1$
- $\text{select}_1(B_{\text{window}}, j + 1) = l_j$

### 3.4 Locate Sampling

**Sampling step:** $s = \Theta(\log |S|)$ (theoretical); fixed $s \in [32, 128]$ in practice. Default: $s = 64$. Stored in the index header.

**SA samples:** $\text{SA\_samples}[k] = SA[k \cdot s]$ for $k = 0, 1, \ldots, \lfloor|S|/s\rfloor$.

**Locate algorithm (forward-only, always available):**

```
locate(row):
    steps ← 0
    r ← row
    while r mod s ≠ 0 and r ≠ sentinel_row:
        r ← LF(r)      // LF(r) = C[BWT[r]] + rank_{BWT[r]}(BWT, r), half-open FM rank
        steps ← steps + 1
    if r == sentinel_row:
        return |S|      // conceptual position; sentinel row never produces a reported match
    return SA_samples[r / s] + steps
```

Cost per occurrence: $O(s)$. Total locate cost: $O(|P| + \text{occ} \cdot s)$ per orientation.

**Optional ISA samples (normative optional — fully specified):**

- **Construction.** During SA-IS construction, before SA is discarded (§2.5), compute the inverse suffix array `ISA[i] = j` such that `SA[j] = i` for all $i \in [0, |S|]$. Then sample: $\text{ISA\_samples}[k] = \text{ISA}[k \cdot s']$ for $k = 0, 1, \ldots, \lfloor|S|/s'\rfloor$. The full ISA is **not** stored; only the sampled array.
- **Sampling step.** $s'$ is a second sampling step stored in the index header alongside $s$. Default: $s' = s = 64$. $s'$ MAY equal $s$ or differ from it; valid range $s' \in [32, 128]$.
- **Locate with ISA (bidirectional).** From a given BWT row $r$, the implementation may walk forward via LF until hitting a row divisible by $s$ (cost ≤ $s$), OR walk forward from the text position via $\Psi$-mapping until hitting a sampled ISA position (cost ≤ $s'$). The cheaper of the two is chosen per occurrence, giving per-occurrence cost $O(\min(s, s'))$ and total locate cost $O(|P| + \text{occ} \cdot \min(s, s'))$ per orientation.
- **Header flag.** `has_isa_samples` (uint8, 0 or 1) indicates presence; if 1, `sample_step_s_prime` and the `ISA_samples` array are present in the FM-Index section of the `.bsi` file.
- **Equivalence guarantee.** Locate with ISA returns the **same set of S-positions** as forward-only locate; the difference is purely in the number of LF/$\Psi$ steps per occurrence. The validation suite includes a byte-identity test for this equivalence (§9.3 of the Architecture).

### 3.5 Window Partitioning (Deterministic — T in Character-Space)

#### 3.5.1 Definition of T

$T$ is the **target window size in S-characters** — that is, a count of positions in the concatenated string S, where each position holds either a read base (codes 0–4) or a separator `#` (code 5). Default: $T = 100{,}000$ characters. Stored in the index header.

**T is a character-space quantity, not a genomic-coordinate quantity.** The window construction algorithm uses S-positions exclusively; genomic coordinates of windows are derived afterwards for use in regional query filtering only.

For a typical Illumina dataset with 150 bp reads, $T = 100{,}000$ characters yields approximately 650–700 reads per window.

**Relating T (S-characters) to genomic distance (reference bases):** At average coverage depth $d$ (reads per reference base), a genomic region of length $L$ reference bases contains approximately $L \cdot d$ S-characters (ignoring separator overhead, which contributes 1 extra character per read). Therefore, the **effective genomic window size** satisfies:

$$T_{\text{genomic}} \approx \frac{T}{d}$$

and the expected number of windows overlapping a genomic region of length $L$ is:

$$|W_r| = O\!\left(\frac{L \cdot d}{T}\right) = O\!\left(\frac{L}{T_{\text{genomic}}}\right)$$

This is the unit-consistent form of the window complexity claim: for fixed coverage $d$ and window size $T$, query cost scales with **region size, not dataset size**, which is the locality property BAMSI delivers. The complexity tables in §4 state the formula in either form; both are equivalent under the coverage normalisation.

#### 3.5.2 Construction Algorithm (Per Chromosome)

For each chromosome $c$ (processed in ascending `chrom_id` order):

```
Let reads_c = reads with chrom_id = c, in ascending read_id order
idx ← 0
while idx < |reads_c|:
    start_read ← reads_c[idx]
    l ← readStarts[start_read.read_id]          // S-position of window start
    tentative_end_S ← l + T − 1                 // tentative window end in S-character-space

    // Extend to include all reads whose S-start falls within [l, tentative_end_S]
    last_idx ← idx
    while last_idx + 1 < |reads_c|
      and readStarts[reads_c[last_idx + 1].read_id] ≤ tentative_end_S:
        last_idx ← last_idx + 1

    last_read ← reads_c[last_idx]
    r ← readStarts[last_read.read_id] + |last_read.seq| − 1   // S-pos of last read base

    // MUST include the trailing '#' separator in this window if not the last read in S
    if last_read.read_id < N − 1:
        r ← r + 1

    // Compute genomic span (for regional query filtering only — not used in construction)
    genomic_start ← start_read.pos                         // 1-based
    genomic_end   ← max over reads_c[idx..last_idx] of:
                    (read.pos + ref_span(read) − 1)        // 1-based reference end

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

**ref_span:** `ref_span(r_i) = sum of len(op)` for ops in `{M, =, X, D, N}` from the read's CIGAR.  
**genomic_end fallback:** If `ref_span(r_i) = 0`, `genomic_end` for that read = `p_i` (the read has no reference-consuming ops).

**Edge case — oversized read:** If $T < |r_i|$ for some read (read sequence is longer than the window target in S-characters), that read forms a window by itself: `l = readStarts[i]`, `r = readStarts[i] + |r_i| - 1` (plus 1 for the separator if not last). A WARNING is logged to stderr. This does not violate any invariant.

**Chromosome boundary rule:** Windows MUST NOT cross chromosome boundaries. The final window of chromosome $c$ ends at the last S-character of the last read of $c$ (plus its separator if that read is not the last read in S overall). The first window of chromosome $c+1$ starts at `readStarts` of the first read of $c+1$. No gap or overlap may exist at chromosome boundaries.

#### 3.5.3 Guarantees

The window construction algorithm guarantees:

- $[0, |S|-1] = \bigcup_j [l_j, r_j]$ — **exact partition, no gaps, no overlaps**
- $r_j + 1 = l_{j+1}$ for all consecutive windows $j$ and $j+1$
- **No read is split:** $\forall i \in \{0..N-1\}, \exists! j$ such that $[\text{readStarts}[i],\; \text{readStarts}[i] + |r_i| - 1] \subseteq [l_j, r_j]$
- `#` belongs to the preceding read's window (the rule `r ← r + 1` when `last_read.read_id < N−1` is mandatory, not optional)
- Each window belongs to exactly one chromosome
- Windows are sorted by `(chrom_id, genomic_start)` — enabling binary search for regional queries

### 3.6 Mapping Layer (Fully Specified)

#### 3.6.1 Definition

$$M_\ell(pos) = (\text{chrom\_name}[c_j],\; [p_{\min}, p_{\max}],\; \text{read\_id},\; \text{query\_strand})$$

where $pos$ is a 0-based position in $S$, $\ell = |P|$ is the pattern length, and `query_strand` is `Forward` or `Reverse` according to which element of Q(P) produced this match.

#### 3.6.2 Implementation (Step-by-Step)

**Step 1 — Identify read:**

```
read_id = rank1(B_read, pos) − 1
```

Using the closed-interval rank convention: `rank1(B_read, pos)` counts 1-bits in $B_{\text{read}}[0..pos]$ inclusive. Since $B_{\text{read}}[\text{readStarts}[i]] = 1$ and these are the only 1-bits, this gives `rank1 = i+1` at any position within read $i$'s span, yielding `read_id = i`. ✓

**Step 2 — Compute read start and match offset:**

```
read_start   = select1(B_read, read_id + 1)   // = readStarts[read_id]
offset_start = pos − read_start
offset_end   = offset_start + ℓ − 1
```

**Step 3 — Retrieve genomic anchor (from S_map only):**

```
(chrom_id, p_i) ← decode per-read block from S_map using dir_map[read_id]
```

S_map contains and provides: `chrom_id` (uint32) and `p_i` (uint64, 1-based SAM POS).

**Step 4 — Retrieve CIGAR (from S_meta only):**

```
cigar ← decode per-read block from S_meta using dir_meta[read_id]
```

S_meta contains and provides: CIGAR ops, FLAG, and optional tags. The mapping layer uses only CIGAR ops from this step.

**Step 5 — Apply CIGAR mapping:**

```
p_min = cigar_ref_pos(cigar, p_i, offset_start, LEFT)
p_max = cigar_ref_pos(cigar, p_i, offset_end,   RIGHT)
```

**Step 6 — Return:**

```
return (chrom_name[chrom_id], [p_min, p_max], read_id, query_strand)
```

where `chrom_name[chrom_id]` is the string resolved from the frozen BAM header chrom_id mapping.

#### 3.6.3 Invariants

- $M_\ell(pos)$ is defined for all $pos \in [0, |S|-1]$ that do not fall on a separator character.
- If $pos$ corresponds to a separator position in S (conceptually S[pos] = code 5), return a structured error (`SEPARATOR_POSITION`); callers MUST check. Separator detection is implemented via bitvector operations on $B_{\text{read}}$ (a separator position is the one immediately after the last base of some read $r_i$; this is determinable from `rank1(B_read, pos)` and cached read lengths) — raw S is never materialised at query time.
- $p_{\min} \geq 1$ always (1-based)
- $p_{\min} \leq p_{\max}$ always
- Mapping is deterministic and BAM-consistent.
- S_map and S_meta are accessed independently via their respective per-read directories (dir_map and dir_meta). Neither access depends on the other stream.

#### 3.6.4 Space

$$|M| = O(|W|) = O(|S|/T)$$

where $|W|$ is the number of windows.

---

## 4. Query Model

### 4.1 Query Definitions

All queries accept $P \in \Sigma^*$, $|P| \geq 1$. Invalid patterns return a structured error (`EMPTY_PATTERN` or `INVALID_PATTERN`).

In strand-complete mode (operational default), let $Q(P)$ be as defined in §0.6. Each query is stated below in strand-complete terms; for single-strand mode, replace $Q(P)$ with $\{P\}$.

---

#### GlobalCount(P)

Returns the total number of exact occurrences of $P$ across all reads in ℛ, summed over all orientations in $Q(P)$.

$$\text{GlobalCount}(P) = \sum_{Q \in Q(P)} (hi_Q - lo_Q)$$

where $[lo_Q, hi_Q)$ is the SA interval returned by `backward_search(Q)`.

**Complexity:** $O(|P|)$ per orientation; $O(|P|)$ total (|Q(P)| ≤ 2, constant).  
**No locate call is required.**

---

#### GlobalExists(P)

Returns `true` if $P$ occurs at least once in any read in ℛ (under any orientation in $Q(P)$); `false` otherwise.

$$\text{GlobalExists}(P) = \exists Q \in Q(P) : hi_Q - lo_Q > 0$$

**Complexity:** $O(|P|)$ per orientation; $O(|P|)$ total.  
**No locate call is required.**

---

#### Locate(P)

Returns one result $(c, [p_{\min}, p_{\max}], \text{read\_id}, \text{query\_strand})$ per occurrence of $P$, across all orientations in $Q(P)$.

**Output ordering (normative):** Results are emitted in the following deterministic order:
1. `query_strand`: Forward results before Reverse results.
2. Within each strand: ascending `chrom_id`.
3. Within same chrom_id: ascending `p_min`.
4. Within same p_min: ascending `p_max`.
5. Within same p_max: ascending `read_id`.
6. Within same read_id: ascending SA row order as final tie-break.

**Complexity:** $O(|P| + \text{occ} \cdot s)$ per orientation; $O(|P| + \text{occ} \cdot s)$ total.

---

#### RegionalCount(P, chrom, [a, b])

Counts occurrences of $P$ whose mapped genomic interval $[p_{\min}, p_{\max}]$ intersects $[a, b]$ on chromosome `chrom`, summed over all orientations in $Q(P)$.

**Match condition:**

$$[p_{\min}, p_{\max}] \cap [a, b] \neq \emptyset$$

i.e., $p_{\min} \leq b$ and $p_{\max} \geq a$.

**Complexity — two-tier guarantee (normative):**

- **BASE tier** (required of every compliant implementation): $O(|P| + \text{occ} \cdot s)$ per orientation worst case. Achieved by plain FM backward search followed by per-occurrence locate and window/interval filtering.
- **ENHANCED tier** (holds only when `enable_sarange = true`): $O(|P| + \text{occ}_r \cdot s + |W_r| \cdot \log(|S|/s))$ per orientation, where $\text{occ}_r$ is the number of occurrences whose mapped interval actually intersects $[a, b]$ and $|W_r|$ is the number of windows overlapping $[a, b]$. Achieved by SARange-based range counting over SA samples (§4.3).

Paper 2 (Theory) claims of tight RegionalCount complexity refer to the ENHANCED tier only and MUST state this explicitly. Paper 1 (Systems) may operate in either tier; empirical benchmarks must report which tier is active.

---

#### RegionalExists(P, T_threshold, chrom, [a, b])

Returns `true` if $\text{RegionalCount}(P, \text{chrom}, [a,b]) \geq T_{\text{threshold}}$; `false` otherwise.  
Terminates early when the $T_{\text{threshold}}$-th qualifying hit is found.

**Complexity:** $O(|P| + T_{\text{threshold}} \cdot s)$ best case (early exit); $O(|P| + \text{occ} \cdot s)$ worst case in BASE tier; with SARange, the worst case tightens to $O(|P| + \text{occ}_r \cdot s + |W_r| \log(|S|/s))$. Per orientation; total multiplied by $|Q(P)| \leq 2$.

---

### 4.2 Execution Model

**Output ordering modes (normative):**

Locate results have a formally specified deterministic ordering (per the `Locate(P)` definition above). To avoid forcing all implementations to buffer and fully sort large result sets in memory, two modes are supported and both are compliant:

- **Streaming mode (default).** Results are emitted in SA-row order as they are computed, without sorting. This mode minimises latency to first result and is recommended for large result sets. Order is deterministic (SA rows are enumerated in a fixed order) but does not match the formal ordering.
- **Sorted mode** (opt-in via `--sort-output` or equivalent API flag). Results are fully buffered and sorted per the formal ordering before emission. Required for reproducibility tests and byte-identity verification.

Both modes return the **same multiset of results**. Implementations MUST support sorted mode for testing; they MAY default to streaming mode for production queries. The mode actually used MUST be reported in query output metadata so downstream consumers can interpret the order.

```
// Step 1: Backward search — all queries, per orientation Q ∈ Q(P)
for each Q ∈ Q(P):
    [lo_Q, hi_Q) ← backward_search(Q)

// Step 2a: GlobalCount
GlobalCount: return Σ_{Q ∈ Q(P)} (hi_Q − lo_Q)

// Step 2b: GlobalExists
GlobalExists: return any (hi_Q − lo_Q > 0) over Q ∈ Q(P)

// Step 2c: Locate — for each orientation, for each SA row
for each Q ∈ Q(P) with query_strand label:
    for each row r in [lo_Q, hi_Q):
        if r == sentinel_row: skip
        pos ← locate(r)
        // Note: S[pos] cannot be '#' because P ∈ Σ* (enforced by INVALID_PATTERN at query
        // entry) and the FM-index returns only rows whose suffix begins with Q ∈ Σ*.
        // The earlier v3.0 defensive check `if S[pos] == '#'` has been REMOVED: it
        // required access to raw S at query time, contradicting §2.5 / §8.2 of the
        // Architecture ("raw S never resident at query time"). The check was provably
        // unreachable for any valid query and is elided.
        map ← M_{|P|}(pos) with query_strand
        if streaming mode: emit result (chrom_name, p_min, p_max, read_id, query_strand)
        else: buffer result
if sorted mode:
    sort buffered results per Locate output ordering (§4.1)
    emit sorted results

// Step 2d: RegionalCount — window-filtered, per orientation
chrom_id ← chrom_to_id(chrom)
count ← 0
for each Q ∈ Q(P):
    [lo_Q, hi_Q) ← backward_search(Q)
    for each row r in [lo_Q, hi_Q):
        if r == sentinel_row: skip
        pos ← locate(r)
        window_id ← rank1(B_window, pos) − 1
        W ← windows[window_id]
        if W.chrom_id ≠ chrom_id: skip              // wrong chromosome
        if W.genomic_end < a or W.genomic_start > b: skip   // window doesn't overlap [a,b]
        map ← M_{|P|}(pos) with query_strand
        if map.chrom_id == chrom_id
          and map.p_min ≤ b and map.p_max ≥ a:
            count ← count + 1
return count

// Step 2e: RegionalExists — early-exit version of RegionalCount
count ← 0
for each Q ∈ Q(P):
    for each row r in [lo_Q, hi_Q):
        [same window filter and mapping as RegionalCount]
        if qualifying:
            count ← count + 1
            if count ≥ T_threshold: return true
return false
```

### 4.3 Optional SARange Structure — Tight RegionalCount Bound (Formal Specification)

**Purpose.** Achieve the ENHANCED tier RegionalCount bound $O(|P| + \text{occ}_r \cdot s + |W_r| \cdot \log(|S|/s))$ per orientation, without locating all $\text{occ}$ backward-search occurrences.

**Data structure.** A **wavelet tree over the sampled SA array** $\text{SA\_samples}[0..\lfloor|S|/s\rfloor]$. Equivalently, any 2D range-searchable structure over the pairs $(i, \text{SA\_samples}[i])$ that supports the operation below in $O(\log(|S|/s))$ time is acceptable and remains compliant with the ENHANCED tier. The wavelet tree is specified as the reference structure; alternatives MUST be documented in the index header (`sarange_variant` field, reserved in §7.1 of the Architecture).

**Supported operation.**

```
range_count(lo, hi, l_j, r_j):
    returns: |{ r ∈ [lo, hi)  :  SA_extrapolated[r] ∈ [l_j, r_j] }|
```

where $\text{SA\_extrapolated}[r]$ is the effective SA value. Because the wavelet tree is built over `SA_samples` at step $s$, exact per-row SA values at non-sampled rows are obtained by a single LF-walk from the row down to the nearest sample (cost ≤ $s$); in practice, the structure is built over a bucketed representation that subsumes this walk into the $O(\log(|S|/s))$ query cost.

**Complexity of `range_count`.** $O(\log(|S|/s))$ per call (wavelet tree range-query standard result — Grossi, Gupta, Vitter 2003).

**Integration with the RegionalCount query loop.** For each window $j$ whose `(chrom_id, genomic_start, genomic_end)` overlaps the query's `(chrom, [a, b])`:

```
for each Q ∈ Q(P):
    [lo_Q, hi_Q) ← backward_search(Q)                       // O(|P|)
    for each relevant window j (|W_r| of them):
        c_j ← range_count(lo_Q, hi_Q, windows[j].l, windows[j].r)   // O(log(|S|/s))
        if c_j > 0:
            locate exactly the c_j occurrences that land in window j  // O(c_j · s)
            filter by exact [p_min, p_max] ∩ [a, b] intersection
```

Total per orientation:

$$\underbrace{O(|P|)}_{\text{backward search}} + \underbrace{O(|W_r| \cdot \log(|S|/s))}_{\text{range count over windows}} + \underbrace{O(\text{occ}_r \cdot s)}_{\text{locate only qualifying rows}}$$

which gives the ENHANCED tier bound.

**Space.** $O((|S|/s) \cdot \log(|S|/s))$ bits for the wavelet tree over `SA_samples`.

**Normative status and tier compliance.**

- **BASE tier (enable_sarange=false, default).** Implementations operating in BASE tier MUST NOT claim the ENHANCED bound. Window filtering (via `B_window.rank1` per occurrence) is still used and still improves average-case performance, but the worst case remains $O(|P| + \text{occ} \cdot s)$ per orientation.
- **ENHANCED tier (enable_sarange=true).** Implementations operating in ENHANCED tier MUST store SARange in the `.bsi` and MUST use the integration above. They THEN guarantee the ENHANCED bound.
- **Paper 2 (Theory)** claims of the tight bound require the ENHANCED tier and MUST state this explicitly.
- **Paper 1 (Systems)** operates in either tier; the chosen tier MUST be declared in benchmark headers.

### 4.4 Window Complexity (Unit-Consistent Form)

Under average coverage depth $d$ and window size $T$ in S-character-space, the number of windows overlapping a genomic region of length $L$ reference bases is:

$$|W_r| = O\!\left(\frac{L \cdot d}{T}\right) = O\!\left(\frac{L}{T_{\text{genomic}}}\right), \qquad T_{\text{genomic}} = \frac{T}{d}$$

(See §3.5.1 for the derivation.) The key locality property is that **$|W_r|$ depends on region size $L$ and coverage $d$, not on the total dataset size $|S|$ or number of reads $N$.**

**Key result.** With the SARange structure (§4.3), RegionalCount time depends on region size and coverage, not on dataset size.

**Without SARange (BASE tier):** Window filtering prunes candidates in practice; worst case remains $O(|P| + \text{occ} \cdot s)$ per orientation.

### 4.5 Approximate Pattern Matching — Formal Extension Point

This section addresses Objective 4 ("robust enough for exact and approximate sequence queries"). The v1.0 query model is **exact-only**, but the index structures specified in §3 are **forward-compatible** with a v2.0 approximate-matching extension that:

- Requires **no re-indexing** (the v1.0 `.bsi` remains valid input for v2.0 search)
- Uses the **same FM-index, B_read, B_window, WindowTable, and mapping layer** as exact search
- Adds an algorithmic overlay (k-mismatch via seed-and-extend; k-edit via bidirectional-FM with bounded backtracking)

**Index forward-compatibility requirements (v1.0 normative):** The v1.0 build MUST produce an index satisfying the following extension hooks, even though v1.0 search does not exercise them:

| Hook | Requirement | Used by v2.0 |
|---|---|---|
| **H1 — Bidirectional FM compatibility** | When `BuildConfig.enable_bidirectional = true` (default false in v1.0), an additional FM-index over the reverse of S (`S^R`) is built and stored alongside the forward FM-index. Both share the same SA-IS run via the lemma `SA(S^R)[i] = |S| − 1 − SA(S)[|S|−1−i]`. | Bidirectional search for k-edit matching (Lam et al. 2009; Li 2013) |
| **H2 — Seed length parameter** | The header stores `recommended_seed_length` (uint8, default `⌊|P|/(k+1)⌋` for the smallest expected pattern; v1.0 stores 16 as a default for short-read use). | Seed-and-extend k-mismatch (pigeonhole principle) |
| **H3 — SA range API stability** | Backward-search returning `[lo, hi)` is a stable public API; v2.0 calls it on each seed independently. | k-mismatch via seed enumeration |
| **H4 — Locate stability** | `locate(row)` is a stable public API; v2.0 uses it to position seed hits before extension. | All approximate algorithms |

**v2.0 query semantics (designated, NOT implemented in v1.0):**

- `ApproxLocate(P, k_mismatch)` — returns matches with at most `k_mismatch` Hamming-distance mismatches.
- `ApproxLocate(P, k_edit, BIDIR)` — returns matches with at most `k_edit` Levenshtein-distance edits, requires `enable_bidirectional = true` at build time.
- `ApproxRegionalCount(P, k, chrom, [a, b])` — region-restricted variant.

**v1.0 normative obligation:** Implementations MUST refuse approximate queries with a structured `NOT_IMPLEMENTED_V1` error (rather than silently returning exact-only results) so that callers cannot misinterpret v1.0 results as approximate. The error message MUST identify the v1.0 limitation and point at the v2.0 future work plan.

**Why this matters.** Listing approximate matching as a flat Non-Goal (as v3.0 did) makes Objective 4 unfulfillable. Promoting it to a formal extension point with v1.0 forward-compatibility hooks **ensures** the v1.0 index is the foundation for a v2.0 approximate-matching release without breaking v1.0 commitments — which is exactly what Objective 4 asks for.

---

## 5. Correctness Guarantees

### Theorem 1 — Soundness

**Statement:** Every result reported by any query corresponds to a true occurrence of $Q$ (for some $Q \in Q(P)$) in some $r_i \in \mathcal{R}$.

**Proof dependencies:**

- FM-index correctness: backward search returns only true Q-occurrences (I3)
- Separator correctness: Q ∈ Σ* cannot match across `#` (I4)
- Stream independence: $M_\ell$ fetches CIGAR from S_meta and coordinates from S_map independently (I10)
- CIGAR mapping definition: total, deterministic, BAM-consistent (I5)
- Mapping correctness: $M_\ell$ returns the correct genomic interval for each occurrence
- Sentinel handling: sentinel row is never reported (I11)
- Rank/select correctness: bitvectors yield correct read_id (I8)

### Theorem 2 — Completeness

**Statement:** Every occurrence of $Q$ (for any $Q \in Q(P)$) in any $r_i \in \mathcal{R}$ is discoverable by the query engine.

**Proof dependencies:**

- FM-index completeness: backward search returns all Q-occurrence SA rows (I3)
- Window partition coverage: all positions in S are covered by exactly one window (I7)
- CIGAR mapping definition: total, no occurrence discarded (I6)
- Mapping coverage: $M_\ell$ is defined for all non-separator positions (I6)
- Sentinel handling: sentinel row is skipped correctly (I11)
- Window integrity: each read belongs to exactly one window (I12)

---

## 6. Space Complexity

$$\text{Total space} = |S| \cdot H_k(S) + O(|S|) \text{ bits}$$

where the $O(|S|)$ term collects the FM-index BWT/Occ structure, the bitvectors $B_{\text{read}}$ and $B_{\text{window}}$, the auxiliary metadata streams $S_{\text{qual}}$, $S_{\text{meta}}$, $S_{\text{map}}$ (each $O(N) = O(|S|)$), and the per-read / block-level directories. The compressed sequence stream $S_{\text{seq}}$ achieves $|S| H_k(S) + o(|S|)$ bits individually — this dominates the storage-efficient component — but the **full sealed index** (including everything required for querying) adds a linear-in-$|S|$ overhead, which is standard for compressed full-text indexes (Ferragina & Manzini 2005; Grossi & Vitter 2005).

**Why not $o(|S|)$?** The FM-index BWT/Occ structure alone costs $|S| \cdot H_0(\text{BWT}) + o(|S|)$ bits. For genomic DNA, $H_0(\text{BWT}) \approx \log_2 5 \approx 2.32$ bits/symbol; $H_k(S)$ for $k \in [4, 8]$ is typically in $[1.7, 2.0]$ bits/symbol. The FM-index term cannot be absorbed into the $o(|S|)$ of $S_{\text{seq}}$ — they are two additive linear terms over the same string length. Stating the total as $|S| \cdot H_k(S) + o(|S|)$ (as v3.0 did) is a **mathematical error** that has been corrected in this version.

**Single-BWT optimisation (optional).** An implementation MAY share a single stored BWT between $S_{\text{seq}}$ and the FM-index by retaining the BWT used for FM-index operations and encoding $S_{\text{seq}}$ on top of that same BWT with MTF+RLE+arithmetic metadata stored separately. When this optimisation is active, the $|S| \cdot H_0(\text{BWT})$ term in the table below does NOT appear as a second independent cost; the BWT is counted once. Implementations MUST declare this choice in the `.bsi` header (`shared_bwt` flag reserved in §7.1 of the Architecture). The space bound remains $|S| H_k(S) + O(|S|)$; with the optimisation, the constant in $O(|S|)$ is smaller.

**Detailed component breakdown:**

| Component | Space | Notes |
|---|---|---|
| $S_{\text{seq}}$ (BWT+MTF+RLE+arithmetic) | $|S| H_k(S) + o(|S|)$ bits | dominant storage-efficient term |
| BWT in FM-index (wavelet tree or RRR bitvectors over BWT) | $|S| H_0(\text{BWT}) + o(|S|)$ bits | linear-in-$|S|$ overhead; elided when `shared_bwt=true` |
| C array | $\sigma \cdot 8$ bytes ($\sigma = 6$; negligible) | |
| Occ structure | subsumed in BWT component | |
| SA samples | $(|S|/s) \cdot 8$ bytes | $O(|S|/s)$ |
| ISA samples (optional) | $(|S|/s') \cdot 8$ bytes | $O(|S|/s')$ |
| $B_{\text{read}}$, $B_{\text{window}}$ | $2|S|$ bits + $o(|S|)$ overhead each | linear-in-$|S|$ |
| WindowTable | $O(|S|/T)$ entries × fixed bytes per entry | sublinear in $|S|$ |
| $S_{\text{qual}}$ | $O(N \cdot \text{avg}|Q_i|)$ compressed | $= O(|S|)$ |
| $S_{\text{meta}}$ | $O(N \cdot \text{avg CIGAR+FLAG+tags size})$ compressed | $= O(N) \subseteq O(|S|)$ |
| $S_{\text{map}}$ | $O(N \cdot (\text{sizeof chrom\_id} + \text{sizeof POS}))$ compressed | $= O(N)$ |
| `dir_meta`, `dir_map` (mandatory per-read) | $2N \cdot 12$ bytes | hot-path query directories |
| `dir_seq`, `dir_qual` (permitted block-level) | $2(N/B_{\text{dir}}) \cdot 12$ bytes if block-level, else $2N \cdot 12$ bytes | reconstruction directories |
| SARange (optional) | $O((|S|/s) \log(|S|/s))$ bits | ENHANCED tier only |

**Summary.** Total space = $|S| \cdot H_k(S) + O(|S|)$ bits. The $|S| \cdot H_k(S)$ term is the dominant storage-efficient component (the $S_{\text{seq}}$ stream); the $O(|S|)$ term is the combined linear-in-$|S|$ overhead of the FM-index, bitvectors, auxiliary streams, and directories. This is the standard form for a compressed full-text index carrying additional auxiliary data.

### 6.1 Global Invariants

| ID | Invariant |
|---|---|
| I1 | S encodes all reads with `#` at correct boundaries and `readStarts` correctly records start positions |
| I2 | $B_{\text{read}}$ has a 1-bit at exactly readStarts[i] for all $i \in \{0..N-1\}$, and 0 elsewhere |
| I3 | FM backward search returns all and only SA rows whose suffixes begin with the query pattern |
| I4 | No pattern $P \in \Sigma^*$ matches across `#`; separators cannot be part of any query match (enforced by `P ∈ Σ*` validation + FM-index correctness; no runtime check on S required) |
| I5 | $M_\ell(pos)$ produces intervals matching BAM CIGAR semantics for all non-separator positions |
| I6 | $M_\ell(pos)$ is defined (returns a valid result) for all non-separator positions; no match is discarded |
| I7 | Windows form a disjoint partition of $[0, |S|-1]$; every position belongs to exactly one window |
| I8 | rank/select identities hold for both $B_{\text{read}}$ and $B_{\text{window}}$ (closed-interval convention) |
| I9 | Same BAM inputs + same BuildConfig → bit-identical `.bsi` file (permits parallel SA-IS variant that is bit-identical to reference single-threaded SA-IS, §0.7) |
| I10 | Each stream (S_seq, S_qual, S_meta, S_map) decodes from its own bytes and embedded header only; `dir_meta` and `dir_map` are mandatory per-read; `dir_seq`/`dir_qual` may be block-level |
| I11 | FM-index behaves as if `$` exists; `sentinel_row` is tracked and never produces a reported match |
| I12 | Each read belongs to exactly one window; no read spans a window boundary |
| I13 | BWT Lifecycle — SA and BWT are computed exactly once per build (Architecture §1.1 Stage 4) and shared by S_seq encoding and FM-index construction; no stage recomputes them |
| I14 | No Raw-S Query Access — no query operation accesses raw S or its uncompressed contents at query time; separator detection uses B_read bitvector operations only |
| I15 | Rank API Separation — the half-open BWT-rank is used only within the FM recurrence and LF-mapping; the closed-interval bitvector-rank is used elsewhere; the two MUST NOT be implemented as a shared function |

### 6.2 Theorem-to-Invariant Mapping

| Theorem | Depends on Invariants |
|---|---|
| Theorem 1 (Soundness) | I3, I4, I5, I8, I10, I11, I15 |
| Theorem 2 (Completeness) | I3, I6, I7, I11, I12, I15 |
| Space bound (§6) | I10, I13 |
| Determinism (§0.7) | I9, I13 |
| Query-time S-free (§4.2) | I14, I15 |

---

## 7. Non-Goals (v1.0 — Out of Scope, with Extension Plan)

BAMSI v1.0 explicitly does not address the following. Each is a deliberate scope choice with a documented extension path; reviewers and users should understand both the limitation and the upgrade trajectory.

- **Approximate pattern matching (v1.0 not implemented; v2.0 designed).** v1.0 returns exact substring matches only. The v1.0 index carries forward-compatibility hooks (§4.5: H1–H4) so that a v2.0 release adds k-mismatch and k-edit query overlays **without re-indexing**. v1.0 implementations MUST refuse approximate queries with `NOT_IMPLEMENTED_V1` rather than silently returning exact results.
- **Dynamic updates (insertions/deletions of reads post-build).** The index is immutable; updates require a rebuild. A log-structured overlay is designated future work (v3.0).
- **Variant calling or genotyping.** BAMSI provides substring counts and positions, not variant inference. Pairs naturally with downstream callers (GATK, DeepVariant, bcftools).
- **Full BAM replacement.** BAMSI is a queryable compressed index, not a BAM viewer. Users requiring full record iteration use a BAM reader; users requiring genomically-placed reads + sequences use BAMSI's reconstruction (Table A of §2.3).
- **Multi-sample joint indexing (v1.0 single-manifest; v2.0 designed).** Each `.bsi` indexes one fixed (possibly multi-BAM) input manifest. Joint multi-sample indices with per-sample metadata and per-sample query restriction are future work (v2.0). The v1.0 index format reserves header fields for sample tagging (`sample_count`, `sample_table`) so that v2.0 can introduce per-read sample IDs without breaking the v1.0 format.
- **Reconstruction of unmapped, secondary, or supplementary records (v1.0 primary-only; v2.0 designed).** ℛ contains only primary mapped reads (§0.2). Consequences and trajectory:
  - BAMSI is **well-suited** to motif finding, coverage QC, k-mer frequency analysis, and any task where secondary/supplementary reads are conventionally filtered before counting.
  - BAMSI is **NOT** a drop-in replacement for tools analysing structural variants (SVs), chimeric reads, or long-read-specific alignments that use supplementary alignments as a primary signal.
  - Paper 1 (Systems) MUST include an empirical calibration: for 10 representative motifs on a standard human WGS BAM, report `GlobalCount_BAMSI(P)` vs. `samtools view -c` over all alignments (including secondary/supplementary). This quantifies the practical gap.
  - Supplementary/secondary indexing is designated future work (v2.0) — likely as parallel auxiliary streams indexed by `(primary_read_id, supp_index)`.
- **Quality-aware queries (v1.0 two-phase; v2.0 designed).** $S_{\text{qual}}$ is compressed and reconstruction-recoverable, but quality scores are not used in the v1.0 query model:
  - Exact substring matching on base calls has clean formal properties (soundness, completeness) that quality-weighted semantics would complicate.
  - The intended v1.0 interface is **two-phase**: `Locate(P)` returns matching `read_id`s; the caller retrieves $Q_i$ via `dir_qual[read_id]` and applies quality-based post-filtering externally.
  - Future work (v2.0) will extend `RegionalCount` and friends to accept a minimum average quality threshold, composing as a post-filter without breaking exact-match guarantees.
- **Positional/genomic-deduplication query semantics.** Strand-complete counts are read-sequence-level, not genomic-position-deduplicated (§0.6 counting semantics). Users requiring deduplication post-process `Locate` output.

---

## 8. Research and System Guarantees

### Paper 1 — Systems Contribution

**Required benchmark battery** (Objective 5):

- **Tools compared.** BAMSI vs. CRAM (samtools 1.21+ with reference), Genozip (latest stable), `samtools view` for query latency baseline, NGC (Roguski & Deorowicz 2014) for the BAM-compression-only comparison. All tools run in their default lossless modes; lossy modes (e.g., Genozip's `--optimize`) reported separately and clearly labelled.
- **Datasets — public, real-world, mandatory.**
  - **NA12878** Illumina HiSeq X 30× WGS BAM from the 1000 Genomes Project (or HG002 from GIAB) — short-read benchmark.
  - **HG002** PacBio HiFi 30× from GIAB — long-read PacBio benchmark.
  - **HG002** Oxford Nanopore (PromethION) ~30× from GIAB — long-read ONT benchmark.
  - At least one **exome** BAM (Agilent SureSelect or equivalent capture) from a public source (e.g., NA12878 exome).
  - At least one **multi-sample / cohort** dataset (e.g., a 10-sample subset of 1000 Genomes) demonstrating scalability — even though v1.0 indexes per-manifest, the per-sample build and aggregate metrics MUST be reported.
- **Metrics, per (tool, dataset).**
  - Compression ratio: bits per base pair, separately for sequence, quality, and total.
  - Build (compress) wall-clock and peak RSS, single-threaded and 8/16-thread settings.
  - Decompression to BAM wall-clock and peak RSS (for tools that support it).
  - **Query latency (BAMSI only vs. `samtools view` baseline):** GlobalCount, GlobalExists, Locate, RegionalCount, RegionalExists for at least 100 patterns spanning length 8–32, region sizes 1 kb / 1 Mb / chromosome / genome.
  - Memory footprint at query time (resident set excluding kernel cache).
- **Tier and lossy-mode disclosure.** Each result row in the benchmark tables MUST report: BAMSI tier (BASE / ENHANCED), `is_lossless` flag, codec choices (`qual_codec_id`, `meta_codec_id`, `map_codec_id`, `entropy_order_k`), and `BuildConfig.shared_bwt`.
- **Reproducibility.** All benchmarks MUST be reproducible from a public Docker/Singularity container (or conda-pinned environment) with a single `make benchmarks` command; the container image hash MUST appear in the paper.
- **Sub-experiments (additionally required).**
  - Entropy-order sensitivity analysis: `k ∈ {4,5,6,7,8}` × three datasets (§2.4).
  - Motif-scope calibration: `GlobalCount_BAMSI(P)` vs. `samtools view -c` (all alignments) on 10 representative motifs (§7).
  - Soft-clip interval characterisation: width of `[p_min, p_max]` vs. precise sub-read mapping on CIGAR edge cases (§1.4.3).
  - Window-pruning benefit: RegionalCount latency vs. region size at fixed `T` (validates §4.4).
  - Directory granularity trade-off: per-read vs. block-level `dir_seq`/`dir_qual` (validates §2.3 / §2.10).
  - Stream-codec ablation: TYPED_SPLIT vs. ZSTD_FALLBACK on S_meta; DELTA_RANGE vs. RAW on S_map; per-cycle-range vs. ZSTD on S_qual (validates §2.7–§2.9).

**Performance targets (informational, not contractual):** On the NA12878 30× BAM (~150 GB BAM, ~75 GB compressed by Genozip), the v1.0 reference implementation should aim for total `.bsi` size within 1.10× of Genozip's output (acknowledging the FM-index overhead) and median `GlobalCount` latency under 50 ms on commodity hardware. These targets are aspirational; falling short does not violate any contract clause but should be analysed in Paper 1.

### Paper 2 — Theory Contribution

- Formal proofs of Theorems 1 (Soundness) and 2 (Completeness) with invariant citations.
- Entropy bound: $|S_{\text{seq}}| = |S| H_k(S) + o(|S|)$ via the BWT–MTF–RLE–arithmetic pipeline, parametric in $k$.
- Query complexity per orientation:
  - **BASE tier (all implementations):** GlobalCount $O(|P|)$, GlobalExists $O(|P|)$, Locate $O(|P| + \text{occ} \cdot s)$, RegionalCount $O(|P| + \text{occ} \cdot s)$, RegionalExists $O(|P| + T_{\text{threshold}} \cdot s)$ best / $O(|P| + \text{occ} \cdot s)$ worst.
  - **ENHANCED tier (SARange-enabled):** RegionalCount $O(|P| + \text{occ}_r \cdot s + |W_r| \log(|S|/s))$, RegionalExists correspondingly tightened. All Paper 2 "tight" complexity claims require and MUST state the ENHANCED tier.
- **Overall space bound:** $|S| H_k(S) + O(|S|)$ bits total, where the $O(|S|)$ term collects the FM-index and auxiliary streams (§6). This replaces the earlier (v3.0) misleading `$o(|S|)$' summary.
- **Rank-convention lemma:** the half-open FM rank and the closed-interval bitvector rank differ by at most 1 at any boundary; the formula `read_id = rank1(B_read, pos) − 1` is exact under the closed-interval convention (§3.3).

---

## 9. Operational Guarantees (Clinical / Production Workflows)

This section addresses Objective 2 (clinical-discovery workflows) by enumerating the operational properties an implementation MUST surface. These are not new algorithmic guarantees — they are obligations on the **interface and observability** of the system.

### 9.1 Determinism and Re-Verification

Every `.bsi` is **bit-identically reproducible** from the same `(BAM input set, BuildConfig)` (I9). A clinical or audit workflow can therefore:

- Re-derive any `.bsi` from primary inputs and verify byte-equality with the deployed file.
- Detect tampering: any modification of the `.bsi` or its inputs causes re-build to produce a different byte-stream and fail the comparison.
- Cite a specific `(.bsi SHA-256, source_manifest_hash, ordering_hash, bamsi_version)` tuple in clinical reports.

### 9.2 Provenance Surface

The `bamsi info <file>.bsi` command (§10.2) MUST emit, in a human-readable and machine-parseable form:

- `bamsi_version`, `format_version` (.bsi)
- `source_manifest_hash`, `ordering_hash` (32-byte SHA-256 each, hex)
- `is_lossless` (boolean) — and, if false, which streams are lossy and the parameters used
- All codec IDs and parameters: `entropy_order_k`, `qual_codec_id` + `qual_lossy_bins`, `meta_codec_id`, `map_codec_id`, `sample_step_s`, `sample_step_s_prime`, `window_size_T`, `enable_sarange`, `shared_bwt`, `enable_bidirectional`, `strand_mode`, `seq_block_size`, `qual_block_size`
- `chrom_name_table` — all chromosome names with their `chrom_id` (so callers can confirm reference compatibility)
- `N_reads`, `S_length`, `N_windows`
- Reference identity: when `reference_based_encoding = 1`, the SHA-256 of the reference FASTA used at build time

### 9.3 Query-Time Reproducibility

Identical query inputs MUST return identical results. Specifically:

- Two invocations of `Locate(P)` in **sorted mode** on the same `.bsi` produce byte-identical output (modulo timestamps in headers).
- In **streaming mode** (§4.2), the multiset of results is identical and the SA-row enumeration order is deterministic; a `--sort-output` flag converts streaming output to sorted output without re-querying.
- Per-orientation results are determinately labelled `Forward` or `Reverse`.
- All numeric outputs (`p_min`, `p_max`, `read_id`) are integer-valued; no floating-point arithmetic appears in any query path.

### 9.4 Error Behaviour

All failures produce structured errors with codes (Architecture §8.1). No silent fallback. Specifically forbidden behaviours:

- Returning fewer results than the true match set (silent under-counting).
- Returning extra results not in the true match set (silent over-counting).
- Continuing after a checksum or hash mismatch.
- Auto-rebuilding a `.bsi` when the input has changed.

When ℛ excludes a record (e.g., supplementary alignment), an audit log entry MAY be emitted; it MUST NOT alter the in-scope query result.

### 9.5 No Implicit Network Access

Query operations MUST NOT make implicit network calls. Reference-based reconstruction (§2.6) MUST surface the reference dependency before fetching it. This requirement supports clinical workflows where data residency is a regulatory constraint.

---

## 10. Distribution and Software Deliverables (Reference Implementation Contract)

This section addresses Objective 6 (open-source software, CLI, documentation, tutorials). It defines the **minimum interface** that a reference implementation MUST expose; alternative implementations are free to add capabilities but MUST honour these.

### 10.1 Licensing and Implementation Language

**Licence.** The reference implementation MUST be released under an **OSI-approved open-source licence**. The recommended licence is **Apache 2.0** (compatible with downstream commercial and clinical use). The licence text MUST appear in `LICENSE` at the source-tree root. A `NOTICE` file lists third-party dependencies and their licences (htslib, libsais or equivalent SA-IS, sdsl-lite or equivalent succinct library).

**Implementation language (normative).** The reference implementation MUST be written in **one** of the following two languages:

- **C++**, conforming to **ISO C++20 or later**.
- **Rust**, **2021 edition or later**.

**No other language is permitted for the v1.0 reference implementation.** Specifically, Python, Java, Go, JavaScript/TypeScript, C#, Julia, OCaml, and any garbage-collected or interpreted language are **excluded** as the implementation language for the reference build.

**Rationale (binding, not advisory):**

1. **Latency targets.** Contract §9.3 query-time reproducibility and the §8 Paper 1 informational target of < 50 ms median `GlobalCount` on commodity hardware require deterministic memory layout, predictable cache behaviour, and the absence of garbage-collection pauses. C++ and Rust provide these; managed-runtime languages do not.
2. **Clinical workflows.** Contract §9 forbids silent fallback and requires structured error behaviour; the language must give the implementer full control over allocation failure paths and lifetime bounds. Both C++ and Rust deliver this; Rust additionally enforces it at compile time.
3. **Direct interoperability with htslib.** htslib is C; the SA-IS and succinct-data-structure libraries (libsais, sdsl-lite) have C/C++ APIs. C++ binds directly; Rust binds via FFI with negligible overhead. Other languages add at least one runtime layer that complicates the §9.5 no-implicit-network and deterministic-layout obligations.
4. **Memory safety where it matters most.** Rust's ownership model eliminates entire classes of CVEs in a security-critical genomic-data tool. C++20 with modern hygiene practices (RAII, smart pointers, bounds-checked spans, sanitisers under §13.4 CI) is acceptable but requires the discipline that Rust enforces.

**Choice between C++ and Rust.** Each project chooses one; mixing is not part of the v1.0 reference. The C++ track and the Rust track are interchangeable from the spec's perspective: both must produce bit-identical `.bsi` files on the same `BuildConfig` (Determinism, I9), expose the CLI surface of §10.2, and satisfy every other clause in §9 and §10. Cross-track byte-identity (a C++ build and a Rust build of the same input producing the same `.bsi`) is **not required** — that would over-constrain the codec implementations. **Within a track**, byte-identity is required across operating systems and architectures (validated in CI per Architecture §13.4).

**Consumer languages (not the reference).** Users MAY write bindings in any FFI-capable language (Python via `cffi` / `ctypes`, Go via `cgo`, Java via JNI, R via `Rcpp` or `extendr`, Julia via `ccall`, etc.) on top of the C ABI surface defined in §10.3. These consumer-side bindings are **not part of the v1.0 reference** and are not maintained by the reference project; they are user territory and the project's API stability commitment (§10.6) makes such bindings practical without coupling them to the reference's language choice.

### 10.2 Command-Line Interface (Minimum Surface)

The reference CLI is invoked as `bamsi <subcommand> [options]`. The following subcommands are normative:

| Subcommand | Purpose | Mandatory flags | Optional flags |
|---|---|---|---|
| `bamsi build` | Build a `.bsi` from one or more BAM inputs | `--input <BAM>...`, `--output <file>.bsi` | `--window-size`, `--sample-step`, `--isa-step`, `--entropy-k`, `--strand`, `--enable-sarange`, `--enable-bidirectional`, `--shared-bwt`, `--seq-block-size`, `--qual-block-size`, `--allow-parallel-sa`, `--threads`, `--reference <FASTA>`, `--qual-codec`, `--meta-codec`, `--map-codec`, `--qual-lossy-bins` |
| `bamsi count` | Run GlobalCount(P) | `--index <file>.bsi`, `--pattern <P>` | `--strand`, `--json` |
| `bamsi exists` | Run GlobalExists(P) | `--index`, `--pattern` | `--strand`, `--json` |
| `bamsi locate` | Run Locate(P) | `--index`, `--pattern` | `--strand`, `--sort-output`, `--output <file>`, `--json`, `--bed` |
| `bamsi region-count` | Run RegionalCount(P, chrom, [a,b]) | `--index`, `--pattern`, `--region <chr>:<a>-<b>` | `--strand`, `--json` |
| `bamsi region-exists` | Run RegionalExists | `--index`, `--pattern`, `--region`, `--threshold` | `--strand`, `--json` |
| `bamsi reconstruct` | Recover (subset of) ℛ from streams | `--index`, `--output <BAM>` | `--streams seq,qual,meta,map` (subset; default all), `--read-ids <list>` |
| `bamsi info` | Print provenance and codec disclosure | `--index` | `--json` |
| `bamsi verify` | Re-compute hashes and section checksums; verify byte-integrity | `--index` | `--strict` |

`--json` produces machine-parseable output for scripted/clinical pipelines. The non-JSON form is human-readable and stable across patch versions.

### 10.3 Library Deliverable and Stable C ABI

The reference implementation provides:

- **Core library in the chosen track.** For the C++ track: `libbamsi.so` / `libbamsi.a` (or platform-equivalent), header `bamsi.hpp`, namespace `bamsi::`. For the Rust track: a Cargo crate `bamsi` published to `crates.io`, with native Rust APIs in idiomatic style.
- **Stable C ABI.** Regardless of track, the reference exports a **stable C ABI** through `bamsi.h`. C++ implementations expose this via `extern "C"` shim functions; Rust implementations expose it via `#[no_mangle] pub extern "C"` functions and `#[repr(C)]` types. This C ABI is the ground truth for any non-reference language binding.

The C ABI surface mirrors the CLI subcommands of §10.2 as a small set of C-callable functions plus a streaming match-iterator type. Function names are `bamsi_<verb>_<noun>` (e.g., `bamsi_index_open`, `bamsi_query_global_count`, `bamsi_iter_next`). The stable C ABI is documented in `docs/api.md` and versioned per §10.6.

**Consumer-language bindings (not part of the v1.0 reference).** Any FFI-capable language MAY produce bindings against the C ABI. The reference project does not ship or maintain such bindings; users and downstream tooling produce them. The project commits to ABI stability (§10.6) so that downstream bindings remain valid across patch releases.

### 10.4 Documentation Deliverables

The reference release MUST ship:

- **README.md** — quick-start (`bamsi build` → `bamsi count`) with an example dataset (downloadable from the project page).
- **`docs/format.md`** — `.bsi` file format specification (mirrors Architecture §7), suitable for independent re-implementers.
- **`docs/algorithms.md`** — algorithmic overview citing Contract §3, §4 and Architecture §4–§6.
- **`docs/cli.md`** — full CLI reference (mirrors §10.2 above with examples).
- **`docs/api.md`** — Stable C ABI reference (the ground truth for all bindings) plus the chosen-track native-API reference (C++ in `bamsi.hpp` or the Rust crate's rustdoc).
- **`docs/clinical.md`** — operational guidance for clinical workflows: provenance verification, audit-log integration, lossy-mode caveats.
- **At least three end-to-end tutorials** (in `docs/tutorials/`):
  1. **Motif counting on NA12878:** download the BAM, build the `.bsi`, count 10 motifs.
  2. **Region-restricted query:** count occurrences of a probe sequence in a 1 Mb region of `chr17`.
  3. **Two-phase quality filtering:** Locate matches, retrieve qualities, post-filter at Q ≥ 30.

### 10.5 Distribution Channels

- **Source:** GitHub (or equivalent) with tagged releases following SemVer (`v1.0.0`, `v1.0.1`, `v1.1.0`, …). Releases include source tarball, signed (GPG) checksums, and changelog.
- **Container:** Docker Hub / GHCR image `bamsi/bamsi:v1.0.0` and `bamsi/bamsi:latest`. Image is reproducibly built from `Dockerfile` in the source tree; image SHA-256 published in release notes.
- **Conda (binary distribution of the CLI and core library):** Bioconda recipe (`conda install -c bioconda bamsi`). Same binary is installed regardless of whether it was built from the C++ track or the Rust track.
- **Track-native package index (Rust track only):** `crates.io` publication of the `bamsi` crate. Not applicable to the C++ track (C++ has no equivalent first-party registry; the `vcpkg` / `conan` recipes are user-territory).
- **Cloud-genomics:** A reference Snakemake / Nextflow workflow (in `workflows/`) demonstrates BAMSI in a typical pipeline alongside `samtools` and `bcftools`.
- **Not distributed by the reference project:** PyPI, npm, CRAN, Maven Central, or any other consumer-language package registry. These are user-territory; the reference exposes a stable C ABI (§10.3) so that downstream maintainers can publish bindings to those registries without coupling them to the reference project.

### 10.6 Versioning, Compatibility, and Support Policy

- **Format versioning.** `.bsi` `format_version` is bumped on any breaking change. The reference implementation supports the current major version and the immediately preceding major version for read (decode); writes always produce the current version. v1.x readers MUST refuse v2.x files with a clear `VERSION_MISMATCH` error.
- **API versioning.** The C ABI follows SemVer: minor versions are ABI-compatible additions; major versions are breaking changes (which are avoided in v1.x).
- **Security disclosure.** A `SECURITY.md` describes the vulnerability-disclosure process. CVEs are tracked.
- **Issue tracker.** Public on GitHub.

---

*End of BAMSI Contract v3.3*
