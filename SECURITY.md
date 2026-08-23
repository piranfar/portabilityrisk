# Security policy

## Scope

This repository contains a data contract, JSON Schemas, an ingestion builder,
tests, and small derived tables. It has no network service, no authentication,
and no runtime that accepts untrusted input from a third party.

The realistic risks are therefore not the usual ones. In order of how likely
they are to actually matter here:

1. **Disclosure** — a credential, a private local path, patient information, or
   licensed source material committed by accident.
2. **Supply-chain** — a malicious change to the builder or the tests that makes
   validation pass when it should fail.
3. **Code execution** — the builder reads files from a path you supply and
   parses TSV, JSON and XML with the standard library. If you point it at a
   hostile source repository, you are executing your own decision to do so.

## Reporting a vulnerability

Please report privately first. Use GitHub's **Report a vulnerability** button
under the Security tab of this repository, which opens a private advisory.

Do not open a public issue for anything that involves a leaked credential or
personal data — that turns a private problem into a public one.

Please include what you found, where, and how to reproduce it. If you are
reporting exposed material, **describe it rather than pasting the value**: the
file path, the commit, and the kind of secret is enough to act on.

Expect an acknowledgement within about a week. This is a research project
maintained by one person, not a product with an on-call rotation, so please
calibrate expectations accordingly.

## Reporting exposed data

If you find any of the following in this repository or its history, it is a
defect and we want to know:

- an API key, token, password, or private key;
- a private filesystem path identifying a person or institution;
- patient-identifying information of any kind;
- a publisher PDF, DOCX, supplement, or full-text XML that should not have been
  redistributed;
- sequence data or a model checkpoint.

None of these should be present. The repository is built around referencing
sources by citation, hash, and locator rather than distributing them, and
`.gitignore` blocks all of the above categories by default.

## What we do about it

For a leaked credential: rotate first, then remove. Removal from the working
tree is not sufficient — history rewriting or repository replacement will be
required, and the credential must be treated as compromised regardless.

For licensed material committed in error: remove and, if it was ever public,
treat redistribution as having occurred and notify the rights holder if they
ask.

For a validation defect that lets an inadmissible record through: fix the rule,
add a negative test that reproduces the defect, and re-run the affected build.
Any derived artefact produced under the broken rule is re-issued rather than
patched in place.

## Handling of credentials in this project

The builder needs no credentials. Where upstream tooling required them, they
were read from environment variables and never written to a file, a log, or a
receipt. No `.env` file, key file, or credential store is committed or expected.

If you are running the upstream acquisition tooling yourself, keep keys in your
environment or a secret manager. Do not put them in a config file inside a
checkout, however private the checkout is.

## Data-handling limitations

This project holds bacterial isolate measurements and public genome accessions.
It holds no patient data and none was accessed. If you believe any row could
identify a person — for example through an unusually specific combination of
isolate metadata — report it privately and it will be removed or generalised.

## Not a clinical system

This repository must not be used for clinical decision-making. It has no
regulatory clearance and no clinical validation. A security report about its
suitability for patient care is out of scope; see `MODEL_CARD.md`, which states
that non-claim formally.
