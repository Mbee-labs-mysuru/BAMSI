# Determinism Failure Runbook

## Purpose

This runbook defines the response protocol when cross-platform determinism fails.

A determinism failure means the same input BAM set and the same BuildConfig did not produce a bit-identical `.bsi` output across supported runners.[file:481][file:482]

## Trigger

This runbook applies when:

- `determinism.yml` fails in CI
- A local reproducibility check produces different `.bsi` bytes for the same input and config
- A rebuild on another supported machine produces different output bytes for the same input and config

## Severity

Determinism failure is **P0**.

No pull request with a determinism failure may merge.[file:481]

## Immediate response

1. Stop further merges if the failing change is not isolated yet.
2. Identify the failing commit or pull request.
3. Reproduce the failure with a small synthetic input first.
4. Compare outputs by checksum and binary diff.
5. Record the incident in the relevant issue or PR thread.

## Standard reproduction input

Use a synthetic 10K-read input as the first reproduction target.[file:481]

Record:

- commit hash
- compiler version
- OS and architecture
- dependency versions
- exact command line
- SHA-256 of each output `.bsi`

## Common causes

Investigate these first, in this order, because the Execution Plan identifies them as the most common practical causes of determinism breaks:[file:481]

1. Parallel iterator collection order
2. Hash-map iteration order
3. Floating-point arithmetic in serialized values
4. Timestamps embedded in headers
5. Uninitialized memory
6. Environment-dependent default thread pool size

## Fix rules

Apply one of these fixes:

- Replace non-deterministic iteration with stable sorted iteration
- Pin ordering explicitly before serialization
- Remove or zero non-deterministic timestamp fields in reproducible builds
- Initialize all serialized memory explicitly
- Pin thread count or seed/order where required
- Replace unstable container traversal with deterministic traversal

## If failure is on a pull request

1. The PR cannot merge.[file:481]
2. The author bisects the failing change using the synthetic 10K-read input.[file:481]
3. The fix must be validated by rerunning the determinism check.
4. The PR description must document the root cause and the fix.

## If failure is on `main`

1. Revert immediately.[file:481]
2. Reproduce the failure under controlled conditions before re-applying any fix.[file:481]
3. Open a P0 incident issue with:
   - failing commit
   - affected platforms
   - reproduction commands
   - output checksums
   - suspected root cause
   - corrective action
4. Do not re-merge until determinism is restored.

## Escalation rule

If three consecutive determinism breaks occur in one week, escalate to the tech lead.

This indicates a structural non-determinism source and requires a design-level fix rather than patch-level fixes.[file:481]

## Evidence to capture

For every determinism incident, preserve:

- input dataset identifier
- BuildConfig used
- compiler and toolchain versions
- dependency lock state
- output SHA-256 values
- diff artifact or mismatch report
- root cause summary
- linked fix commit

## Exit criteria

A determinism incident is closed only when:

- the failure is reproduced
- the root cause is identified
- the fix is merged
- CI determinism checks pass again
- the incident record includes the final cause and corrective action
