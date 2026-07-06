# Security Policy

## Supported Versions

| Version | Supported          |
|---------|--------------------|
| 1.0.x   | :white_check_mark: |
| < 1.0   | :x:                |

## Reporting a Vulnerability

If you discover a security vulnerability in BAMSI, please report it
responsibly:

1. **Do NOT open a public GitHub issue** for security vulnerabilities.
2. Email security reports to: **bamsix-security@example.com**
3. Include:
   - Description of the vulnerability
   - Steps to reproduce
   - Potential impact assessment
   - Suggested fix (if any)

## Response Timeline

- **Acknowledgement**: Within 48 hours of report receipt
- **Triage**: Within 5 business days
- **Fix timeline**: Critical vulnerabilities within 14 days;
  moderate within 30 days; low within 90 days

## Scope

Security vulnerabilities include but are not limited to:
- Buffer overflows in the .bsi reader/writer
- Integer overflows in SA/FM-index operations
- Arbitrary code execution via malformed .bsi files
- Information disclosure through timing side-channels
- Denial of service via crafted inputs

## CVE Tracking

CVEs will be tracked and published via GitHub Security Advisories.

## Clinical Workflow Note

BAMSI is designed for use in clinical genomics workflows (Contract §9).
Security vulnerabilities that affect data integrity, provenance
verification, or determinism guarantees are treated as **critical**
regardless of exploitability.
