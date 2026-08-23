"""NM-V1 deterministic structural decision engine -- frozen before application.

Ground truth is derived from REFERENCE STRUCTURAL EVIDENCE ONLY. The original HMM label is not
an input to any rule; it is used solely to define one audit stratum and is never shown to the
adjudicator. Every rule uses evidence already named in the frozen rubric.

Rules are evaluated in the fixed order below and are mutually exclusive: the first matching
rule wins, so exactly one label is assigned to every block.
"""
import argparse, datetime, hashlib, json, os, sys

VERSION = "nmv1_freeze_rule_engine_v1.0.0"
DESIGN_SHA = "c2aea6cb583c24b997ab376861acc600295e94d59bfa0ef2d55cf2bcc424bb20"
RULES_SHA = "e454873a89d1d56b20adb9ae157f20224076966daa0ed34bcbe44318063260f6"


def sha256_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


E = {
 "engine": "NMV1_STRUCTURAL_DECISION_ENGINE",
 "version": "1.0.0", "builder": VERSION,
 "frozen_before_application": True,
 "frozen_before_any_aggregate_performance": True,
 "nmv1_design_sha256": DESIGN_SHA,
 "groundtruth_rules_sha256": RULES_SHA,

 "principle": ("Ground truth is constructed from reference structural evidence only. The "
               "original HMM annotation is NOT an input to any rule. Using the method under "
               "test to build the truth it is judged against would be circular."),

 "hmm_label_usage": {
   "as_rule_input": "PROHIBITED",
   "as_audit_stratum_definition": "PERMITTED - selection only, never label assignment",
   "shown_to_adjudicator": "PROHIBITED until ground truth and the audit are frozen"},

 "derived_evidence_predicates": {
   "IS_complete_n": "count of ISEScan elements with type = c",
   "IS_bilateral_TIR_n": "count of type = c elements having BOTH terminal inverted-repeat "
                         "coordinate pairs (start1/end1 and start2/end2) with irLen > 0",
   "IS_complete_orf_n": "count of type = c elements with orfLen > 0",
   "IS_strong": "IS_bilateral_TIR_n > 0 AND IS_complete_orf_n > 0",
   "IS_partial_only": "isescan_partial > 0 AND IS_complete_n == 0",
   "INT_complete": "if_complete > 0",
   "INT_intI": "an element annotated intI inside an integron of type complete",
   "INT_attC": "an element of type_elt attC inside an integron of type complete",
   "INT_strong": "INT_complete AND INT_intI AND INT_attC",
   "INT_incomplete_only": "(if_calin + if_in0) > 0 AND NOT INT_complete",
   "no_reference_evidence": "isescan_n == 0 AND (if_complete + if_calin + if_in0) == 0",
   "boundary_problem": "truncated == yes OR wrapped_circular == yes",
   "tool_problem": "tool_status != ok"},

 "rule_order": [
  {"id": "F1", "label": "biologically_indeterminate",
   "condition": "tool_problem",
   "reason": "missing or failed tool output; the evidence base is incomplete"},
  {"id": "F2", "label": "biologically_indeterminate",
   "condition": "boundary_problem",
   "reason": "truncated or origin-wrapped block; an element may extend beyond the retrieved "
             "sequence and boundaries cannot be resolved"},
  {"id": "A", "label": "multiple_MGE_evidence_supported",
   "condition": "IS_strong AND INT_strong",
   "reason": "complete IS with resolved bilateral TIRs AND complete integron with integrase "
             "plus attC"},
  {"id": "B", "label": "IS_associated_supported",
   "condition": "IS_strong AND NOT INT_strong",
   "reason": "at least one complete IS, complete transposase ORF, resolved bilateral TIRs, "
             "no truncation or boundary ambiguity"},
  {"id": "C", "label": "integron_associated_supported",
   "condition": "INT_strong AND NOT IS_strong",
   "reason": "complete integron structure with integrase plus attC, no tool failure or "
             "boundary ambiguity"},
  {"id": "D", "label": "chromosomal_mobile_supported",
   "condition": "(IS_complete_n > 0 OR INT_complete) AND NOT (IS_strong OR INT_strong)",
   "reason": "structurally supported mobile-element evidence at complete level that does not "
             "meet the exclusive IS or integron bar, so it cannot be assigned to one class"},
  {"id": "F3", "label": "biologically_indeterminate",
   "condition": "IS_partial_only OR INT_incomplete_only",
   "reason": "partial-only element or incomplete integron; insufficient to establish a class"},
  {"id": "E", "label": "chromosomal_quiescent_supported",
   "condition": "no_reference_evidence AND NOT tool_problem AND NOT boundary_problem",
   "reason": "no structural MGE evidence from either reference path, both tools completed, "
             "sequence evaluable, no boundary warning"},
  {"id": "F4", "label": "biologically_indeterminate",
   "condition": "otherwise",
   "reason": "conflicting evidence the frozen rules cannot resolve"}],

 "properties": {
   "deterministic": True, "mutually_exclusive": True, "total": True,
   "first_match_wins": True,
   "every_block_receives_exactly_one_label": True},

 "outputs_per_block": ["rule_based_label", "rule_id", "evidence fields supporting the rule",
                       "evaluable_status (evaluable | indeterminate)"],

 "expert_audit": {
   "purpose": "a blinded human check that the rule engine agrees with expert biological "
              "judgement; it replaces exhaustive manual review, it does not replace scrutiny",
   "max_cases": 120,
   "strata": [
     {"id": "S_IS", "target": 20, "definition": "rule_id == B (complete-IS calls)"},
     {"id": "S_INT", "target": 20, "definition": "rule_id == C (complete-integron calls)"},
     {"id": "S_MULTI", "target": 20, "definition": "rule_id == A (multiple-evidence calls)"},
     {"id": "S_QUIET", "target": 20, "definition": "rule_id == E (quiescent calls)"},
     {"id": "S_DISC", "target": 20,
      "definition": "HMM class disagrees with the rule-based mobile/quiescent reading AND the "
                    "rule label is NOT indeterminate (structurally resolved). The HMM label is "
                    "used here for SELECTION only and is never shown."},
     {"id": "S_INDET", "target": 20,
      "definition": "rule_id in F1, F2, F3, F4 (indeterminate, boundary or tool failure)"}],
   "shortfall_rule": "if a stratum has fewer than its target, take all of it and reallocate "
                     "the remainder proportionally across the remaining strata, in the fixed "
                     "order above, before any human decision is seen",
   "seed": 20260821,
   "selection_rule": "deterministic: candidates sorted by SHA-256 of "
                     "'<block_id>|<seed>|<stratum_id>', first n taken; strata processed in the "
                     "fixed order so selections are mutually exclusive",
   "tokenised_and_shuffled": True,
   "withheld_from_adjudicator": ["species", "ARG identity", "BioProject", "accession",
                                 "original HMM class", "tool names", "rule-based label",
                                 "rule id", "aggregate results", "sampling stratum",
                                 "audit stratum"],
   "provided_to_adjudicator": ["evidence-rich feature maps identical in format to package V2",
                               "the fixed seven-outcome rubric"]},

 "audit_gate": {
   "primary_metric": "overall agreement between the adjudicator's outcome and the rule-based "
                     "label across the audited cases, with a 95 per cent CI",
   "ci_method": "Wilson",
   "SUCCESS": "agreement >= 0.90 AND lower 95 per cent CI bound >= 0.80",
   "REVISE": "agreement 0.80 to 0.899, OR lower CI bound 0.65 to 0.799",
   "FAIL": "agreement < 0.80, OR evidence of systematic disagreement within a major class",
   "class_specific_agreement_reported": True,
   "scope_of_failure": "a failure confined to one evidence class revises THAT CLASS ONLY and "
                       "does not invalidate unaffected classes",
   "indeterminate_handling": "indeterminate cases are scored for agreement on INDETERMINATE "
                             "STATUS itself and are never forced into a mobile or quiescent "
                             "label",
   "gate_frozen_before_review": True},

 "prohibited_until_audit_is_frozen": [
   "sensitivity, specificity, PPV, NPV, balanced accuracy",
   "NM-V1 gate verdicts",
   "opening either unblinding key"],
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--design", required=True)
    ap.add_argument("--rules", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    if sha256_file(a.design) != DESIGN_SHA or sha256_file(a.rules) != RULES_SHA:
        print("REFUSING: frozen design or rules digest mismatch"); sys.exit(1)
    if os.path.exists(a.out):
        print("REFUSING: %s exists; a frozen engine is never overwritten" % a.out); sys.exit(1)
    E["frozen_utc"] = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    json.dump(E, open(a.out, "w", encoding="utf-8", newline="\n"), indent=2)
    print("%s\n  wrote %s" % (VERSION, a.out))
    print("  FROZEN ENGINE SHA-256: %s" % sha256_file(a.out))
    print("  rules in order: %s" % ", ".join(r["id"] for r in E["rule_order"]))
    print("  audit strata  : %s  (max %d cases)"
          % (", ".join(s["id"] for s in E["expert_audit"]["strata"]),
             E["expert_audit"]["max_cases"]))
    print("  HMM label as rule input: %s" % E["hmm_label_usage"]["as_rule_input"])


if __name__ == "__main__":
    main()
