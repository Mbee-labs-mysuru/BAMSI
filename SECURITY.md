# Security Policy

## Reporting a Vulnerability

Please do **NOT** report security issues via public GitHub issues, discussions, or other public forums.

**Email vulnerabilities to**: `abhishek.dp263@gmail.com`

Include:
- BAMSI version affected
- Environment details (OS, compiler, dependencies)
- Steps to reproduce
- Impact assessment

## Response SLA

- **Critical** (RCE, data loss): 24 hours
- **High** (arbitrary code execution, memory corruption): 72 hours  
- **Medium** (denial of service): 5 business days
- **Low**: Next release cycle

## Disclosure

We follow [Responsible Disclosure](https://en.wikipedia.org/wiki/Responsible_disclosure).

1. Acknowledge report within 24 hours
2. Assess impact and validity within 72 hours
3. Develop patch privately
4. Release fix + CVE (if applicable)
5. Public disclosure 90 days after fix (or sooner if critical)

## Supported Versions

| Version | Supported     | End-of-Life |
|---------|---------------|-------------|
| 1.0.x   | ✅ Yes        | Dec 2027    |
| 0.x.x   | ❌ No         | ---         |
