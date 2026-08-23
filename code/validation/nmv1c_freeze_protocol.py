"""NM-V1C corrected confirmatory audit -- protocol freeze, before any case is drawn.

The classifier under validation is the ORIGINAL frozen rule engine, unchanged. No threshold is
altered: revising F3 in light of the failed audit would convert validation into model
development, and the 19/20 partial-only finding is deferred to a future revised classifier.

The single corrected defect is instrument contamination: Method X, the HMM path under test, was
visible to the adjudicator while the ground-truth engine was barred from using it. Method X is
removed from the instrument entirely.
"""
import argparse, collections, csv, datetime, hashlib, json, os, sys

VERSION = "nmv1c_freeze_protocol_v1.0.0"
ENGINE_SHA = "ed5db383bb0afe1a1a8433886d6666fe72c324975de99c6763a37824d51c2bee"
GT_SHA = "1beecaa39048f4df52a3235f2dbc538056af9adc6912ee057cbac8ca55b85897"
SEED = 20260822
ALLOC = {"A": 25, "B": 34, "D": 1, "E": 60}


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True)
    ap.add_argument("--pool", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    D = a.dir
    if sha(os.path.join(D, "NMV1_RULE_ENGINE_FROZEN.json")) != ENGINE_SHA:
        print("REFUSING: engine digest mismatch"); sys.exit(1)
    if sha(os.path.join(D, "NMV1_RULE_BASED_GROUND_TRUTH.tsv")) != GT_SHA:
        print("REFUSING: ground truth digest mismatch"); sys.exit(1)
    if os.path.exists(a.out):
        print("REFUSING: %s exists" % a.out); sys.exit(1)
    P = json.load(open(a.pool, encoding="utf-8"))
    gt = {r["block_id"]: r for r in csv.DictReader(
        open(os.path.join(D, "NMV1_RULE_BASED_GROUND_TRUTH.tsv"), encoding="utf-8"),
        delimiter="\t")}
    avail = P["available"]
    pool = collections.Counter(gt[b]["rule_id"] for b in avail)

    Q = {
     "protocol": "NMV1C_CORRECTED_CONFIRMATORY_AUDIT",
     "version": "1.0.0", "builder": VERSION,
     "frozen_utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
     "frozen_before_any_case_was_drawn": True,
     "frozen_before_any_outcome_was_viewed": True,

     "classifier_under_validation": {
       "engine": "NMV1_RULE_ENGINE_FROZEN.json", "sha256": ENGINE_SHA,
       "unchanged": True,
       "F3_threshold_unchanged": True,
       "why": "revising F3 now, in light of the failed audit, would convert validation into "
              "model development. The 19/20 partial-only finding is recorded as a "
              "developmental finding for a future revised classifier and is not acted on here."},

     "defect_being_corrected": {
       "defect": "Method X, the HMM path under test, was displayed to the adjudicator while "
                 "the ground-truth engine was prohibited from using it",
       "evidence_of_determinism": "of 47 engine-QUIESCENT cases, the 38 with a Method X marker "
                                  "drew 0 QUIESCENT calls and the 9 without one drew 9 of 9 "
                                  "QUIESCENT calls",
       "correction": "Method X is removed from the instrument in every form"},

     "sampling_frame": {
       "frame": "the frozen 1,283-block NM-V1 sample; the frame itself is NOT extended",
       "exclusion": "every block delivered in any previous package - V1, V2, R2 or R3",
       "previously_shown": len(P["excluded"]),
       "unused_population": len(avail),
       "exclusion_components": {
         "mandatory_adjudication_all_packages": 495, "V1_QC": 80, "V2_QC": 100,
         "R2_audit": 120, "R3_audit": 120,
         "note": "components overlap; the union is the figure above"}},

     "STRUCTURAL_LIMITATION_STATED_BEFORE_DRAWING": {
       "finding": "the unused population contains ZERO NON_EVALUABLE blocks and ZERO "
                  "integron-exclusive (rule C) blocks",
       "cause": "all 114 indeterminate blocks and all 8 rule-C blocks were routed to MANDATORY "
                "adjudication by the frozen rules, so every one of them was delivered in the "
                "V1/V2 packages",
       "consequence": "this confirmatory audit can validate the MOBILE versus QUIESCENT "
                      "discrimination ONLY. It cannot validate the NON_EVALUABLE state, and it "
                      "cannot produce an integron-specific estimate.",
       "why_the_frame_was_not_extended": "sampling beyond the frozen 1,283 blocks would require "
                                         "running ISEScan and IntegronFinder on new sequence, "
                                         "extending the frozen frame. That is a larger deviation "
                                         "than the stated scope limitation.",
       "what_this_does_and_does_not_test": {
         "tested": "whether the engine's mobile/quiescent call agrees with an independent "
                   "blinded expert seeing the same evidence - the discrimination class B rests on",
         "not_tested": "the NON_EVALUABLE boundary, which is exactly where the deferred "
                       "partial-element question lives"}},

     "selection": {
       "seed": SEED,
       "algorithm": "candidates sorted by SHA-256 of '<block_id>|<seed>|<stratum>', first n "
                    "taken; strata processed in the fixed order A, B, D, E so selections are "
                    "mutually exclusive",
       "strata_are_frozen_machine_categories": True,
       "strata_never_exposed_to_the_adjudicator": True,
       "no_case_selected_on_apparent_ease_evidence_or_likely_agreement": True,
       "allocation": ALLOC,
       "pool_by_stratum": {k: pool.get(k, 0) for k in ("A", "B", "C", "D", "E",
                                                       "F1", "F2", "F3", "F4")},
       "sample_size": sum(ALLOC.values()),
       "balance": "MOBILE 60 (A 25 + B 34 + D 1) and QUIESCENT 60 (E 60), chosen for equal "
                  "power in the two states that can be sampled"},

     "weighting": {
       "method": "inverse probability of selection, weight = pool_size / allocated, per stratum",
       "retained": True,
       "why": "strata are sampled at different fractions, so unweighted estimates are biased"},

     "primary_metric": {
       "mapping": {"MOBILE": ["A", "B", "C", "D"], "QUIESCENT": ["E"],
                   "NON_EVALUABLE": ["F1", "F2", "F3", "F4"]},
       "metric": "design-weighted agreement between the adjudicator's three-state call and the "
                 "rule-engine state",
       "ci": "stratified bootstrap, 2000 resamples, seed 20260822",
       "unweighted_reported_as_secondary": True},

     "gates": {
       "SUCCESS": "design-weighted agreement >= 0.90 AND bootstrap lower bound >= 0.80",
       "REVISE": "agreement 0.80 to 0.899, OR lower bound 0.65 to 0.799",
       "FAIL": "agreement < 0.80, OR systematic disagreement within a sampled state",
       "thresholds_carried_over_unchanged_from_the_original_frozen_gate": True,
       "state_specific_agreement_reported": True,
       "scope_of_failure": "a failure confined to one state revises that state only"},

     "indeterminate_handling": {
       "adjudicator_may_still_choose_NON_EVALUABLE": True,
       "reason_code_required_for_NON_EVALUABLE": True,
       "scoring": "a NON_EVALUABLE call against a MOBILE or QUIESCENT machine label counts as "
                  "a disagreement in the primary metric, and is additionally reported "
                  "separately by reason code, since no NON_EVALUABLE ground truth exists in "
                  "this sample",
       "never_forced": "the adjudicator is never required to pick MOBILE or QUIESCENT"},

     "instrument": {
       "form": "single-file offline HTML application",
       "shows": ["Method Y structural IS evidence", "Method Z integron evidence", "topology",
                 "boundary and truncation status", "block-relative feature coordinates",
                 "an anonymous target interval", "independent tool completion status"],
       "removes_completely": ["Method X feature tracks", "Method X counts",
                              "HMM-derived coordinates",
                              "transposase or integrase hits originating from the tested HMM path",
                              "any hidden Method X text in HTML, JavaScript, metadata or exports"],
       "withholds": ["species", "gene identity", "BioProject", "accession", "stratum",
                     "machine classification", "expected label", "rule-engine label"],
       "choices": ["MOBILE", "QUIESCENT", "NON_EVALUABLE"],
       "definitions": {
         "MOBILE": "independent structural evidence supports a mobile-element context",
         "QUIESCENT": "independent tools completed normally and no credible mobile-element "
                      "evidence is present",
         "NON_EVALUABLE": "independent evidence is incomplete, ambiguous, boundary-truncated, "
                          "technically failed or insufficient to establish or exclude mobile "
                          "context"},
       "reason": "optional for MOBILE and QUIESCENT; REQUIRED for NON_EVALUABLE using fixed "
                 "reason codes",
       "reason_codes": ["INCOMPLETE_ELEMENT", "AMBIGUOUS_EVIDENCE", "BOUNDARY_TRUNCATED",
                        "TOOL_FAILURE", "INSUFFICIENT_EVIDENCE"],
       "features": ["keyboard shortcuts", "progress indicator", "local autosave",
                    "export to TSV and JSON"]},

     "tool_versions_pinned": {
       "isescan": "1.7.3", "integron_finder": "2.0.6", "hmmer": "3.3.2",
       "prodigal": "2.6.3", "infernal": "1.1.4", "blast": "2.17.0",
       "note": "no tool is rerun; the frozen outputs of the original run are reused"},
     "input_digests": {},

     "prior_audit_disposition": {
       "registered_result": "FAIL, 62/120, three-state agreement 0.5167",
       "status": "retained as an instrument-failure diagnostic; not erased, not passed, not "
                 "rescored, and not the definitive biological validation",
       "no_adjudication_reused": True},

     "still_prohibited_until_the_completed_file_is_returned_and_frozen": [
       "opening the confirmatory key", "sensitivity, specificity, PPV, NPV, balanced accuracy",
       "any NM-V1 gate verdict"],
    }
    for f in ("NMV1_RULE_ENGINE_FROZEN.json", "NMV1_RULE_BASED_GROUND_TRUTH.tsv",
              "nmv1_block_evidence_table.tsv"):
        Q["input_digests"][f] = sha(os.path.join(D, f))
    json.dump(Q, open(a.out, "w", encoding="utf-8", newline="\n"), indent=2)
    print("%s\n  wrote %s" % (VERSION, a.out))
    print("  FROZEN PROTOCOL SHA-256: %s" % sha(a.out))
    print("  unused population %d | allocation %s | n=%d"
          % (len(avail), ALLOC, sum(ALLOC.values())))
    print("  STRUCTURAL LIMITATION: NON_EVALUABLE and integron strata are EMPTY in the")
    print("  unused population; this audit covers MOBILE vs QUIESCENT only.")


if __name__ == "__main__":
    main()
