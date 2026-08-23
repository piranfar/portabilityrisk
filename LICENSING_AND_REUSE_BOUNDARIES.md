# Licensing and reuse boundaries — PortabilityRisk (Paper 1)

**Dated 2026-08-23.** Decisions by Vahhab Piranfar, repository owner, recorded verbatim in
`PORTABILITYRISK_OWNER_DECISIONS_LICENCE_FUNDING_AND_DEPOSIT_V1.json`. Nothing here was chosen by
the analysis agent.

Four different things carry four different terms. Conflating them is the most likely way to get
this wrong.

| what | licence | where it lives |
|---|---|---|
| manuscript, supplementary information, figures | **CC-BY-NC-ND-4.0** | preprint server, public repository |
| NM-V5 occurrence-level dataset and its metadata | **CC-BY-NC-ND-4.0** | **Zenodo only** |
| analysis code | **Apache-2.0** (unchanged) | both repositories |
| restricted private artefacts | **no licence — not distributed at all** | private repository only |

---

## 1. Manuscript and preprint — CC BY-NC-ND 4.0

**Permits:** copying, distributing and sharing the manuscript in unmodified form, for
non-commercial purposes, with attribution.

**Restricts:** commercial use, and distribution of adapted or derivative versions. Translations,
remixes, and reformatted or excerpted redistributions are derivatives.

**Does not restrict:** quoting under fair dealing or fair use, reading, citing, or the
peer-review process.

## 2. Dataset — CC BY-NC-ND 4.0

Applies to the occurrence-level dataset deposited in Zenodo and to the dictionaries, schemas,
table catalogue, provenance map and checksums deposited beside it.

**Reproducibility implications of the ND term, stated as consequences, not as objections.**
The owner has made this decision; these are the facts a reader should know, and they belong in the
record rather than in an argument:

- **Verification is unaffected.** Downloading the dataset, re-running the published checks,
  re-deriving the denominators and reporting a disagreement are all reading, not adapting. The
  entire verification path this project was built for remains open.
- **Redistributing a modified copy is not permitted** without separate permission. A reader who
  filters, reshapes, joins or subsets the table for their own analysis may do so privately, but may
  not publish that derived table under this licence.
- **Pooling into a larger public resource requires permission.** Aggregators that redistribute
  harmonised copies of many datasets cannot take this one without asking.
- **Commercial reuse requires permission**, including inside a commercial product or service.

Requests for either are a matter for the owner, not for this document.

**Attribution requirement.** Cite the Zenodo record with its DOI and version, and the manuscript.
Provenance fields must travel with the values wherever they go.

## 3. Code — Apache-2.0

Unchanged. Permissive: use, modification, distribution and commercial use are all permitted,
subject to the licence's attribution and notice terms and its patent grant.

**The code licence does not extend to the data.** Running Apache-2.0 code over CC-BY-NC-ND data
does not place the output under Apache-2.0.

## 4. Third-party material — governed by its original providers

**No licence is applied to any of the following by this project, retrospectively or otherwise.**

- **Source genome data.** NCBI assemblies, BioSample and BioProject records are used under NCBI's
  terms. The dataset deposited here contains **accessions and computed values derived from them**,
  not redistributed sequence.
- **Third-party software and databases.** AMRFinderPlus, MOB-suite, ISEScan, HMMER, BLAST,
  FragGeneScan, CONJScan/MacSyFinder and their reference databases retain their own licences. The
  versions used are recorded in `PORTABILITYRISK_SOFTWARE_AND_DATABASE_VERSIONS_V1.tsv`.

Where a value in the dataset is a direct restatement of a third-party record — an accession, a
molecule designation — the licence on this dataset governs **this compilation**, not the
underlying record.

## 5. Restricted private artefacts — not distributed

Unblinding keys, blinded casebooks, adjudication instruments and completed workbooks, SSH keys and
any credential, and the restricted recovery bundle. These are **not published anywhere**, not to
GitHub and not to Zenodo, and no licence question arises because no distribution occurs.

The **scoring rubric** is published, because an audit whose rules are secret cannot be checked. The
keys that would unblind it are not.

## 6. Public visibility is not a licence

A file being readable in a public repository grants nothing beyond what its licence grants.
Anything in the public repository not explicitly covered by one of the three licences above — for
example a file inherited from a third party, or a value restating a third-party record — remains
under its own terms. Absence of a licence header is not an implied grant.

The reverse also holds: **the dataset licence is not changed by the dataset being downloadable.**

## 7. Applying this in practice

| you want to | manuscript | dataset | code |
|---|---|---|---|
| read, cite, quote | yes | yes | yes |
| re-run the published verification | n/a | **yes** | yes |
| share an unmodified copy, non-commercially, with attribution | yes | yes | yes |
| publish a modified or derived version | **ask** | **ask** | yes |
| use commercially | **ask** | **ask** | yes |
| redistribute inside an aggregated public resource | **ask** | **ask** | yes |

"Ask" means contact the owner. It does not mean prohibited.

## 8. Funding

Reproduce verbatim wherever a funding statement is required:

> This research received no specific grant from any funding agency in the public, commercial or
> not-for-profit sectors. Cloud-computing resources were supported through promotional credits
> provided by Oracle Cloud Infrastructure.

Oracle promotional credits are **not a grant**. They must not be entered in a funder-award field,
and no grant number exists to cite.
