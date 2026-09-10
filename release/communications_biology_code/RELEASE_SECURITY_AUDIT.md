# Release security audit

**Scope:** non-hidden files under `release/communications_biology_code/`, excluding this audit and `SHA256SUMS.tsv`.
**Status:** PASS
**Release candidate:** `v1.0-submission-rc1` (local label only)

| Check | Status | Notes |
|---|---|---|
| Raw GWAS summary statistics | PASS | No raw GWAS files are included in the release candidate. |
| LD reference files | PASS | No LAVA, UK Biobank or other LD matrices are included. |
| Restricted/participant-level data | PASS | No individual-level data are included. |
| Credentials and secrets | PASS | Credential/token scan completed. |
| Absolute personal paths | PASS | Path scan completed. |
| Private-key filenames | PASS | Filename scan completed. |
| Upload/push status | PASS | No upload, push or archive deposition performed. |

## Findings

No suspicious absolute paths, credential patterns or private-key filenames were found in the scanned release files.
