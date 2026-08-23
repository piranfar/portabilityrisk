"""REDACTED PUBLIC DERIVATIVE. Canonical private artefact: audit/ingest/assay_aware_emergence/v2/nm_validation/nmis_freeze.py, SHA-256 f5a65fc6c1731b5d875264a12a6ab4229c81af49e766c6d9294139b21cfc2329. Infrastructure identifiers replaced by [REDACTED:...]; no scientific content altered."""
"""NMIS design freeze -- written before ISEScan is run on any new block.

NMIS is owner-approved STRUCTURAL CORROBORATION and sensitivity analysis. It is not a rescue
analysis, does not reclassify anything, and must not alter class B or any frozen NM-DIST result.
"""
import argparse, csv, datetime, hashlib, json, os, sys

VERSION = "nmis_freeze_v1.0.0"
SEED = 20260822
L = 10000
LAND = [1000, 2000, 5000, 10000]

INPUTS = [
 ("audit/data/derived/pr_context/out/shared_context_blocks.tsv", "the 21,955 frozen blocks"),
 ("audit/data/derived/pr_context/out/arg_neighbourhood_windows.tsv",
  "occurrence to block key, ARG coordinates, topology, wrap and truncation flags"),
 ("docs/nature_microbiology/nmdist_occurrence_block_distances.tsv",
  "NM-DIST occurrence table: species group, BioProject, block weights, censoring"),
 ("docs/nature_microbiology/NMDIST_FROZEN_PROTOCOL_V1.json", "the NM-DIST frozen design"),
 ("docs/nature_microbiology/NMDIST_RESULT_RECEIPT_V1.json", "the NM-DIST result being corroborated"),
]


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", required=True)
    ap.add_argument("--fasta-manifest", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    if os.path.exists(a.out):
        print("REFUSING: %s exists" % a.out); sys.exit(1)

    dig = {}
    for rel, role in INPUTS:
        p = os.path.join(a.repo, rel)
        if not os.path.exists(p):
            print("REFUSING: missing %s" % rel); sys.exit(1)
        dig[rel] = {"sha256": sha(p), "role": role}

    man = list(csv.DictReader(open(a.fasta_manifest, encoding="utf-8"), delimiter="\t"))
    agg = hashlib.sha256("|".join("%s:%s" % (r["filename"], r["sha256"])
                                  for r in man).encode()).hexdigest()

    Q = {
 "protocol": "NMIS_STRUCTURAL_IS_RECONSTRUCTION", "version": "1.0.0", "builder": VERSION,
 "frozen_utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
 "frozen_before_isescan_was_run_on_any_new_block": True,

 "status": "OWNER-APPROVED STRENGTHENING ANALYSIS. Structural corroboration and sensitivity "
           "only. NOT a rescue analysis. No additional biological computation is required to "
           "support the already frozen Paper-1 claims.",

 "primary_question": "Does the A. baumannii short-range spatial signature reported by NM-DIST "
                     "remain when the endpoint is restricted to structurally resolved COMPLETE "
                     "insertion sequences carrying a complete transposase ORF, bilateral terminal "
                     "inverted repeats and valid element boundaries?",

 "MUST_NOT_CHANGE": {
   "primary_occurrences": 74349, "class_assignments": "A-E unchanged", "class_B": 16303,
   "nmdist_homology_marker_results": "unchanged",
   "in_block_out_of_window_censored": 111,
   "occurrence_level_mge_positivity_pct": 46.39,
   "block_level_mge_positivity_pct": 30.14,
   "statement": "NMIS is structural corroboration and sensitivity analysis, not reclassification"},

 "input_scope": {
   "blocks": 21955,
   "census": "COMPLETE - includes the 1,283 blocks previously sampled in NM-V1, so the whole "
             "output is produced under one uniform environment",
   "new_sequence_retrieval": "NOT AUTHORISED",
   "window_limit": "no analysis extends beyond the frozen +/-10 kb occurrence windows",
   "fasta_aggregate_sha256": agg, "fasta_files": len(man),
   "boundary_rule": "if an IS boundary extends outside available sequence: preserve the call as "
                    "partial/boundary-limited; do NOT infer the missing boundary; do NOT retrieve "
                    "flanking sequence; do NOT exclude the block"},

 "tool": {
   "name": "ISEScan", "version": "1.7.3",
   "conda_package": "isescan-1.7.3-h7b50bb2_0.tar.bz2",
   "conda_md5": "22cf959a984f2cd2d96606209f11a058",
   "executable_sha256": "a6601aab75963dbde04d0100925711fb0f4d601fe85fa536558a21d3d5e16967",
   "dependencies_matched_to_NM_V1": {
     "hmmer": "3.3.2-hdbdd923_4 (md5 805f24f8e8109e5d3d16cfad6799af46)",
     "blast": "2.17.0-h66d330f_0 (md5 405ce6d52eba06fcd48197ae1eb8f5a9)",
     "fraggenescan": "1.32-h7b50bb2_1 (md5 af437ea81e5d02a58f4904bdaa622f6e)",
     "biopython": "1.88-py314h5bd0f2a_0 (md5 da61549005699072c0a8ff45b110da8c)",
     "python": "3.14.7"},
   "environment_lock": "env/isescan_v1.lock",
   "models": {"clusters.faa.hmm": {"bytes": 40472487, "sha256_prefix": "f5de1bf92059e39f"},
              "clusters.single.faa": {"bytes": 100094, "sha256_prefix": "da7b2b67f57394cc"}},
   "command_line": "isescan.py --seqfile <block.fna> --output isescan_out --nthread 1",
   "identical_to_NM_V1": True,
   "DEPENDENCY_DIFFERENCE_FOUND_AND_CORRECTED": {
     "what": "the default conda solve returned HMMER 3.4; NM-V1 used HMMER 3.3.2",
     "why_it_matters": "HMMER drives ISEScan transposase-ORF detection, so it is a "
                       "scientifically relevant dependency",
     "action": "the environment was rebuilt pinned to hmmer=3.3.2, reproducing all five "
               "NM-V1 packages exactly; the isescan.py digest is byte-identical to NM-V1's",
     "declared_here": "so the correction cannot be mistaken for a silent substitution"},
   "parsing_rule": "ISEScan .tsv output is parsed BY ITS OWN HEADER, never by a hard-coded "
                   "column list. NM-V1 error register entry 11 records a defect caused by "
                   "assuming 22 columns when ISEScan emits 24."},

 "structural_definitions": {
   "complete_structural_IS": "an ISEScan element with type == 'c' AND a complete transposase ORF "
                             "(orfLen > 0 with both orfBegin and orfEnd resolved) AND bilateral "
                             "terminal inverted repeats (irLen > 0 with both flanks, "
                             "start1/end1 and start2/end2, resolved)",
   "complete_transposase_ORF": "orfBegin and orfEnd both present and orfLen > 0",
   "bilateral_TIR_requirement": "both inverted-repeat flanks resolved: start1, end1, start2 and "
                                "end2 all present, with irLen > 0",
   "partial_element": "any ISEScan element that is not complete_structural_IS: type == 'p', or "
                      "type == 'c' lacking a resolved ORF or lacking bilateral TIRs",
   "boundary_limited_element": "an element whose interval touches the edge of the available "
                               "window sequence (isBegin <= 1 or isEnd >= sequence length). "
                               "Recorded as partial/boundary-limited and never excluded.",
   "no_structurally_resolved_IS": "the block yielded no ISEScan element meeting "
                                  "complete_structural_IS",
   "tool_failure": "ISEScan exited non-zero or produced no parseable output; recorded, never "
                   "silently dropped or substituted"},

 "distance_and_membership": {
   "occurrence_window_membership": "an element counts for an ARG occurrence only if it lies "
                                   "within that occurrence's own +/-10 kb window, exactly as in "
                                   "NM-DIST. Elements elsewhere in a larger shared block do not "
                                   "count.",
   "distance": "0 when the ARG interval and the element interval overlap, otherwise the bp gap "
               "between nearest boundaries; topology-aware on circular replicons",
   "right_censoring": "an occurrence with no qualifying complete structural IS within 10 kb is "
                      "RIGHT-CENSORED at 10,000 bp; no distance beyond 10 kb is estimated",
   "twenty_kb": "NOT EVALUABLE"},

 "primary_unit_and_weighting": {
   "estimand": "block-balanced occurrence distribution, identical to NM-DIST",
   "rule": "each block contributes total weight 1; an occurrence in a block of m occurrences "
           "receives weight 1/m"},

 "uncertainty": {"method": "BioProject cluster bootstrap", "resamples": 2000, "seed": SEED,
                 "intervals": "percentile 95 percent",
                 "note": "blocks nest within BioProjects, verified in NM-DIST"},

 "groups_and_contrasts": {
   "groups": ["A. baumannii", "P. aeruginosa", "Klebsiella group"],
   "definitions": "identical to NM-DIST and NM-V4/NM-V4C",
   "primary_contrasts": {"N1": "A. baumannii versus Klebsiella group",
                          "N2": "A. baumannii versus P. aeruginosa"},
   "primary_landmarks_bp": LAND,
   "multiplicity": "Holm across the two primary contrasts at each landmark"},

 "gates": {
   "SUCCESS": "for BOTH primary contrasts, the difference in F(d) is positive with a "
              "BioProject-bootstrap 95 percent CI excluding zero at BOTH 1 kb and 2 kb after Holm "
              "correction",
   "REVISE": "direction preserved in both contrasts but at least one of the 1 kb or 2 kb "
             "intervals includes zero",
   "FAIL": "the direction reverses at any primary landmark with a CI excluding zero",
   "interpretation_of_a_smaller_subset": "a smaller structurally complete subset is EXPECTED and "
                                         "is not a failure"},

 "interpretation_rules": {
   "prohibited": ["treating HMM-only calls as false positives",
                  "treating partial ISEScan calls as biologically absent",
                  "treating structurally complete elements as proof of current transposition",
                  "any HGT, transfer, mobilization or phenotype claim",
                  "any distance beyond 10 kb"],
   "required_separate_reporting": ["structurally complete", "partial/boundary-limited",
                                    "no structurally resolved IS", "tool failure"]},

 "execution_gates": {
   "verify_fasta_hashes_locally": "done before transfer",
   "two_sided_manifest": True,
   "benchmark_blocks": 100,
   "benchmark_selection": "deterministic: the 100 blocks with the smallest SHA-256 of block_id",
   "proceed_automatically_if": "projected runtime <= 10 hours AND new disk use <= 20 GB",
   "chunking": "deterministic, resumable, per-block .done markers",
   "failures": "every failure preserved; no silent exclusion or substitution"},

 "stop_conditions": [
   "the 21,955 FASTA set does not verify byte-for-byte on both sides",
   "ISEScan settings cannot be reproduced as recorded for NM-V1",
   "projected runtime exceeds 10 hours or projected new disk use exceeds 20 GB",
   "any frozen NM-DIST or class-B quantity would change",
   "a post-outcome change to endpoint, group, threshold or definition becomes necessary"],

 "governance": {"server": "[REDACTED:INSTANCE_AND_ADDRESS]",
                "prohibited_hosts": "the two prohibited hosts must not be probed or accessed",
                "prohibited": ["commit", "push", "merge", "branch deletion",
                               "manuscript rewrite", "private-release rebuild",
                               "PlasmidCall", "Paper 2", "MicContext",
                               "new biological sequence retrieval"]},
 "input_digests": dig,
    }
    json.dump(Q, open(a.out, "w", encoding="utf-8", newline="\n"), indent=2)
    print("%s\n  wrote %s" % (VERSION, a.out))
    print("  FROZEN PROTOCOL SHA-256: %s" % sha(a.out))
    print("  FASTA aggregate: %s (%d files)" % (agg, len(man)))


if __name__ == "__main__":
    main()
