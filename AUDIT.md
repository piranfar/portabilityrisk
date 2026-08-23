# Publication-safety audit

**Run 2026-08-23, before this repository was created.** Re-runnable; the scripts are
`tools/public_release/audit_public_repo.py` and `scope_audit_public.py` in the working repository.

## Result

| category | verdict | hits |
|---|---|---:|
| public IP address | **clean** | 0 |
| cloud resource identifier | **clean** | 0 |
| instance / volume name | **clean** | 0 |
| absolute machine path | **clean** | 0 |
| credential material | **clean** | 0 |
| email address | **clean** | 0 |
| adjudication artefact, by filename | **clean** | 0 |
| adjudication artefact, by content shape | **clean** | 0 |
| other-paper code file | **clean** | 0 |
| contract-gate leakage | **clean** | 0 |
| other-paper mention | review — 44, all boundary statements | 44 |

**No blocking finding.**

## Infrastructure redaction

An earlier iteration of this repository retained the NMIS compute instance's name and IP address
on the grounds that they appear in a digest-anchored record. **That was the wrong call for a
public tree**, and it has been corrected: the canonical artefacts keep their infrastructure
fields and their digests *in the private repository and release*, and this repository carries
**redacted derivatives** instead.

Four files are derivatives:

| public file | what was replaced |
|---|---|
| `protocols/NMIS_FROZEN_PROTOCOL_V1.json` | instance name and address in the governance block |
| `code/validation/nmis_freeze.py` | the same governance string |
| `validation/NMV3_RESULT_STATUS.md` | instance name |
| `protocols/environments/isescan_v1.lock` | absolute environment prefix |

Each carries an in-band marker and appears in `REDACTION_MANIFEST.json` with **both** its own
SHA-256 and the canonical private SHA-256, so the difference is documented rather than hidden.

**The derivative does not and must not claim the canonical digest.** The frozen NMIS protocol
digest quoted throughout the manuscript — `5438045a3b73d1233…` — refers to the **private
canonical**. The public derivative hashes to something else by design. That is the honest way to
publish a record whose provenance value depends on a field that must not be published, and the
alternative — silently editing the canonical so the hash still "matches" — would be provenance
theatre.

No scientific field, value, threshold, protocol term or digest reference was altered by redaction.

## Sibling-paper mentions — 44, retained

Every one is a boundary statement: *"No PortabilityEvent or PlasmidCall artefact was read"*,
`"plasmidcall_predicted_location": {"status": "RESERVED, EMPTY"}`, and a claim registered for
withdrawal precisely because it depends on another paper's tool. Naming a sibling project in order
to record that it was excluded is exactly what should be public.

Five scripts that carried PlasmidCall **integration design** rather than a boundary statement —
whole sections titled "Relationship to PlasmidCall", a figure of a future integration — were
**removed** from this repository. They generate a superseded internal report, not the manuscript,
and are Paper 3 planning material.

## Rules that were wrong, and were fixed

Recorded because a scan that is quietly wrong is worse than no scan.

- **Credential material** first fired on `url += "?api_key=" + _KEY`, where `_KEY` is
  `os.environ.get("NCBI_API_KEY")`. That is the correct pattern. A rule that cries wolf on correct
  code gets ignored the day it matters; it now requires a literal value.
- **Adjudication-sensitive** first fired 113 times on the NM-V1 scripts *naming the files they
  write*. Class F is about files, not vocabulary — the scripts must be able to say "unblinding
  key". The rule now matches filenames and the column headers only a real key or casebook would
  carry, and finds none.
- **The infrastructure sweep** first flagged 153 "Windows paths" that were the `s:` of `https:`.
- **A glob is not a defect**: `schemas/*.schema.json` is content, and the asterisk rule now says so.

## What the first build got wrong

It contained six `pe_v3_*.py` scripts belonging to Paper 2, copied by a directory glob that was too
broad, and it retained the instance IP. Both are why the audit runs against the **assembled tree**
rather than against the file list that was meant to be assembled.

## Public history

This repository was **deleted and recreated** from a clean sanitised tree. The IP address was
present in both commits of the previous repository, and removing it only from the latest commit
would have left it reachable in history. There is no earlier history to inspect: the first commit
here is the sanitised one.

## Scope

This audit covers disclosure safety: credentials, personal identifiers, infrastructure detail,
blinding integrity and cross-project boundaries. **It does not review the science.** For that see
`validation/`, and in particular the claim status matrix.
