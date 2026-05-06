# Contributing to BAMSI

Thanks for your interest in contributing to BAMSI. A repository-level contributing guide helps contributors know how to open issues, structure pull requests, and follow project conventions.[web:603][web:608]

## Scope

BAMSI is a research and engineering project for lossless BAM file compression and exact query support over compressed data. Contributions are welcome in code, tests, documentation, benchmarks, tutorials, and issue reporting.

## Before contributing

Before starting work, please read these files first:

- `PROJECT.md`
- `docs/decisions/`
- `docs/cli.md`
- `docs/api.md`
- `SECURITY.md`

If the intended change affects behavior, file format, architecture, or benchmarks, open an issue or discussion first. Small fixes such as typo corrections, broken links, or obvious documentation errors may go straight to a pull request.[web:602]

## Types of contributions

Contributions are especially welcome in these areas:

- Build system and repository scaffolding
- Tests: unit, integration, and synthetic
- Documentation and tutorials
- CLI polish and error handling
- Determinism and reproducibility
- Benchmark harness and reporting
- Bug reports with clear reproduction steps

## Reporting issues

When opening an issue, include:

- A clear title
- Expected behavior
- Actual behavior
- Steps to reproduce
- OS, compiler, and dependency details
- If relevant, a small input example or command snippet

Structured issue reports make it easier to reproduce and fix problems.[web:604]

## Pull request process

Please follow this workflow:

1. Create a feature branch from `main`.
2. Make one focused change per pull request.
3. Add or update tests for behavior changes.
4. Update documentation if user-visible behavior changes.
5. Reference the related issue or decision record in the PR description.
6. Open the PR early as draft if feedback is needed.[web:602][web:608]

## Commit guidance

Use short, descriptive commit messages. Keep each commit focused and avoid mixing unrelated refactors with functional changes, because small and scoped changes are easier to review and maintain.[web:606]

Recommended commit style examples:

- `build: add benchmark scaffold`
- `docs: add ADR 0004 dataset details`
- `test: add passing integration stub`
- `cli: split version command handler`

## Coding expectations

For C++ contributions:

- Target the repository’s selected C++ standard and toolchain
- Prefer clear, deterministic behavior over cleverness
- Keep functions small and single-purpose
- Avoid introducing hidden global state
- Preserve cross-platform build compatibility
- Keep public interfaces stable unless the change is explicitly planned

Follow the repository’s build, lint, and test rules before requesting review.[web:608]

## Testing expectations

Every code change should include the appropriate level of validation:

- Unit tests for isolated logic
- Integration tests for end-to-end behavior visible at module boundaries
- Synthetic tests for correctness against controlled inputs

If a change fixes a bug, include a regression test when possible. Clear testing expectations are a standard part of effective contribution guidelines.[web:608][web:603]

## Documentation expectations

Update documentation when changing:

- CLI flags or output
- Public APIs
- File format behavior
- Build instructions
- Benchmark assumptions
- Determinism or reproducibility rules

Docs changes do not need to wait for a major feature branch. Small documentation pull requests are welcome.

## Review expectations

Maintainers may request:

- Smaller scope
- Clearer tests
- Better commit hygiene
- More precise documentation
- Alignment with an ADR or execution-plan stage

Please keep discussion technical, direct, and respectful.

## Security issues

Do not report security vulnerabilities in public issues. Follow the process in `SECURITY.md` instead.[web:571][web:574]

## License

By contributing to this repository, you agree that your contributions are provided under the project license terms in `LICENSE` and accompanied notices in `NOTICE`.[web:566][web:573]
