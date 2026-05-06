# BAMSI v1.0 File Format Specification

## Status
Draft stub for Stage 1.  
Normative target document for Stage 6 completion.

## Purpose
This document specifies the `.bsi` binary file format produced by the BAMSI reference implementation.  
It is intended to be precise enough for an independent re-implementer to parse, validate, and inspect BAMSI files.

## Scope
This specification covers:
- File-level structure
- Header fields
- Stream sections
- FM-index section
- Bitvector section
- Window section
- Directory section
- Footer
- Checksums, versioning, and validation behavior

This document does not yet provide full byte diagrams for every field.  
Those are added as the implementation stabilizes.

## Normative language
The key words **MUST**, **MUST NOT**, **SHOULD**, **SHOULD NOT**, and **MAY** are used in their normative sense.

## Legend
- **Normative**: required for compliance
- **Informative**: explanation, rationale, or implementation guidance
- **Stub**: section exists now but will be expanded later
- **Track-agnostic**: applies to both C and Rust implementations

## Canonical references
- BAMSI Contract v3.3
- BAMSI Architecture v4.3
- BAMSI Execution Plan v2.0

## Format overview
A `.bsi` file is a single sealed binary file written in little-endian form.  
The high-level layout is:

1. Header
2. Sseq section
3. Squal section
4. Smeta section
5. Smap section
6. FM-index section
7. Bitvector section
8. Window section
9. Directory section
10. Footer

## Endianness
All multi-byte integers are little-endian.

## Header
### Purpose
The header records format identity, build provenance, frozen build parameters, chromosome-name mapping, and section-interpretation metadata.

### Required fields
- magic
- formatversion
- bamsiversion
- hostosid
- cpuarchid
- buildtimestamputc
- islossless
- sourcefilecount
- sourcemanifesthash
- orderinghash
- Slength
- Nreads
- Nwindows
- samplesteps
- hasisasamples
- samplestepsprime
- enablesarange
- sarangevariant
- sharedbwt
- enablebidirectional
- recommendedseedlength
- windowsizeT
- entropyorderk
- qualcodecid
- quallossybins
- metacodecid
- mapcodecid
- strandmode
- sentinelrow
- chromcount
- chromnametable
- seqblocksize
- qualblocksize
- allowparallelsa
- referencebasedencoding
- referencesha256
- flags

### Notes
This section will later include:
- exact field order
- exact integer widths
- byte-level layout table
- version-specific constraints

## Stream sections
### Purpose
The four stream sections store the compressed payloads for:
- Sseq
- Squal
- Smeta
- Smap

### Common section shape
Each stream section contains:
- payloadlength
- codecid
- codecmetadata
- payloadbytes
- sectionchecksum

### Notes
Future expansion will define codec metadata layouts for:
- sequence entropy configuration
- quality codec parameters
- metadata codec parameters
- mapping codec parameters

## FM-index section
### Purpose
The FM-index section stores the compressed full-text index over the concatenated read-sequence string.

### Expected contents
- BWT length
- BWT payload or shared-BWT ownership semantics
- C array
- Occ metadata
- Occ payload
- SA samples
- optional ISA samples
- optional SARange payload
- section checksum

### Notes
This section must document shared-BWT behavior explicitly once serialization is frozen.

## Bitvector section
### Purpose
Stores the serialized succinct bitvectors used for read-boundary and window-boundary resolution.

### Expected contents
- Bread
- Bwindow
- sectionchecksum

## Window section
### Purpose
Stores the ordered window table used for regional filtering.

### Per-window fields
- chromid
- l
- r
- firstreadid
- lastreadid
- genomicstart
- genomicend

## Directory section
### Purpose
Stores four directories in fixed order:
- dirseq
- dirqual
- dirmeta
- dirmap

### Normative constraints
- `dirmeta` MUST be per-read
- `dirmap` MUST be per-read
- `dirseq` MAY be per-read or block-level
- `dirqual` MAY be per-read or block-level

### Notes
This section will later define:
- directory header layout
- per-read entry layout
- block-level entry layout
- checksum coverage
- size-accounting examples

## Footer
### Purpose
Closes the file and provides whole-file integrity coverage.

### Expected contents
- globalchecksum
- footermagic

## Validation behavior
A compliant reader/writer must define behavior for at least:
- corrupt/truncated file detection
- checksum mismatch
- ordering hash mismatch
- source manifest mismatch
- version mismatch
- unsupported codec
- lossy reconstruction guardrails

## Versioning policy
- The format version MUST increase on breaking format changes.
- A v1.x reader MUST reject unsupported future major versions with a clear version-mismatch error.
- Patch releases MUST NOT silently alter the binary layout for the same format version.

## Error mapping
This document will later map format-level failures to structured BAMSI error codes, including:
- CORRUPTBSI
- CHECKSUMMISMATCH
- ORDERINGHASHMISMATCH
- MANIFESTMISMATCH
- VERSIONMISMATCH
- STREAMDECODEERROR
- UNSUPPORTEDCODEC

## Compliance notes
This specification is track-agnostic.  
The C and Rust reference implementations may differ internally, but observable `.bsi` format behavior must conform to the same contract.

## Open stub sections
The following sections still need full byte-accurate expansion:
- Header field table
- Stream codec metadata tables
- FM-index serialization details
- Bitvector serialization details
- Directory entry encoding
- Footer checksum coverage examples
- Worked minimal `.bsi` example
