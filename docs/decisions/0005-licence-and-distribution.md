# ADR 0005 – Licence and Distribution

## Status
Accepted

## Context
Contract recommends an OSI-approved licence (Apache 2.0) and Execution Plan lists open-source distribution channels.[file:2]

## Decision
- Licence: Apache License 2.0.
- Files:
  - `LICENSE` containing the Apache 2.0 text.
  - `NOTICE` listing third-party dependencies and their licences.

- Distribution channels for v1.0.0:
  - GitHub releases with:
    - Static-linked binaries for Linux x86_64.
    - Source tarball.
  - Docker image hosted on Docker Hub or GHCR.

- Stretch goals:
  - Bioconda recipe after v1.0.0.
  - Additional package managers only if time and resources permit.

## Rationale
- Apache 2.0 is permissive and compatible with clinical and commercial use.
- GitHub + Docker are accessible to both research and production users with minimal overhead.

## Consequences
- Stage 1 will add the licence/governance files to the repo.[file:2]
- The v1.0.0 paper will refer to the Apache 2.0 licence and the release/distribution channels.
