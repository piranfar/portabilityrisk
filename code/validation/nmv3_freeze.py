"""NM-V3 tool and database robustness -- design freeze, before any robustness result is computed.

Written after the read-only feasibility and input-integrity audit and before any alternative
classification is calculated. The marker-basis inventory (which observed marker combination
underlies each frozen C/D/E assignment) WAS computed before this freeze, because the owner's
instruction ordered that dependency analysis ahead of the freeze. That exposure is declared in
the document itself: no gate threshold here is derived from those counts. Every threshold is
either carried over verbatim from NM-0 or justified on external grounds stated inline.

Nothing is overwritten. The frozen 74,349-occurrence primary dataset and the A-E classification
are read-only inputs to this module.
"""
import argparse, csv, datetime, hashlib, json, os, sys

VERSION = "nmv3_freeze_v1.0.0"

# inputs whose digests are recorded at freeze time
INPUTS = [
    ("audit/data/derived/pr_context/out/frozen_portability_class_definitions.json",
     "frozen class definitions; pins mob_typer 3.1.9 and the database digest"),
    ("audit/data/derived/pr_context/out/FROZEN_PORTABILITY_CONTEXT_PROTOCOL_V1.json",
     "frozen context protocol"),
    ("audit/data/derived/pr_context/out/determinant_portability_classes.tsv",
     "the frozen 74,349-occurrence A-E classification; READ ONLY"),
    ("audit/data/derived/pr_context/out/plasmid_mobility_annotation.tsv",
     "per-replicon MOB-suite markers; the layer under test"),
    ("audit/ingest/assay_aware_emergence/v2/pr_context/pr_context_mobility_classes.py",
     "the classification code whose decision points are perturbed in Arm R"),
    ("docs/nature_microbiology/NM0_VALIDATION_PROTOCOL_V1.json",
     "NM-0 module definition; note status 1.0.0-PROPOSED"),
    ("docs/nature_microbiology/NM0_CLAIM_EVIDENCE_GAP_MATRIX.tsv",
     "C09 and C10 gates, verbatim"),
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
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    if os.path.exists(a.out):
        print("REFUSING: %s exists; a freeze never overwrites one" % a.out); sys.exit(1)

    digests = {}
    for rel, role in INPUTS:
        p = os.path.join(a.repo, rel)
        if not os.path.exists(p):
            print("REFUSING: missing input %s" % rel); sys.exit(1)
        digests[rel] = {"sha256": sha(p), "role": role}

    # the deterministic subset, fixed here and never re-drawn
    mob = os.path.join(a.repo, "audit/data/derived/pr_context/out/plasmid_mobility_annotation.tsv")
    acc = sorted({r["replicon_accession"] for r in csv.DictReader(
        open(mob, encoding="utf-8"), delimiter="\t")})
    ranked = sorted(acc, key=lambda x: hashlib.sha256(x.encode()).hexdigest())
    subset = ranked[:1000]
    subset_hash = hashlib.sha256("|".join(subset).encode()).hexdigest()

    Q = {
 "protocol": "NMV3_TOOL_AND_DATABASE_ROBUSTNESS",
 "version": "1.0.0", "builder": VERSION,
 "frozen_utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
 "frozen_before_any_robustness_result_was_computed": True,

 "DECLARED_PRE_FREEZE_EXPOSURE": {
   "what_was_seen": "the marker-basis inventory: for each frozen C/D/E assignment, which of "
                    "relaxase / mpf / oriT was observed. Counts were reported to the owner "
                    "before this freeze because the owner's instruction ordered that dependency "
                    "analysis (item 3) ahead of the freeze (item 4).",
   "what_was_NOT_computed": "no alternative classification, no transition matrix, no headline "
                            "recomputation, no concordance statistic, no gate outcome",
   "consequence_accepted": "no threshold in this document may be justified by those counts. "
                           "Every threshold below is either carried over verbatim from NM-0 or "
                           "justified on grounds stated inline that do not reference them.",
   "auditable": "the marker-basis inventory is reproducible from the frozen inputs at any time; "
                "a reader can confirm that no gate value coincides with a natural break in it"},

 "STATUS_CORRECTION": {
   "finding": "NM0_VALIDATION_PROTOCOL_V1.json self-declares version 1.0.0-PROPOSED, status "
              "'PROPOSED - OWNER APPROVAL REQUIRED', executable_in_nm0 false.",
   "consequence": "NM-0's NM-V3 panel is a proposal, not an approved contract. This document "
                  "treats NM-0's numeric bounds as binding anyway, because the owner's "
                  "instruction adopted them and because adopting a pre-existing bound is "
                  "strictly safer than inventing one now.",
   "what_is_genuinely_frozen": ["FROZEN_PORTABILITY_CONTEXT_PROTOCOL_V1.json",
                                "frozen_portability_class_definitions.json",
                                "PORTABILITYRISK_MASTER_REPORT_V1.md (immutable)",
                                "PUBLICATION_SCOPE_BOUNDARIES_V1.md (binding)"],
   "gap_recorded": "FROZEN_PORTABILITY_CONTEXT_PROTOCOL_V1.json never names MOB-suite and pins "
                   "no mobility-tool digest. The entire mobility toolchain pin lives in "
                   "frozen_portability_class_definitions.json and in the classification code. "
                   "No stop condition anywhere gates the MOB-suite binary or database digest."},

 "question": "Which C/D/E assignments are invariant to mobility-tool version, mobility-database "
             "version, marker scheme and decision rule, and which are not, so that a future "
             "database update cannot silently invalidate a published figure.",

 "hard_constraints_carried_over": {
   "denominator_74349_may_not_change": True,
   "A_E_classification_is_read_only": True,
   "A_versus_B_is_out_of_scope": "governed by NM-V1, which PASSED on 2026-08-22",
   "location_layer_bound": "0 pct - any change to a documented molecule designation is a STOP "
                           "CONDITION, not a sensitivity",
   "alternative_marker_arm_is_concordance_only":
     "never reported as sensitivity, specificity or accuracy against truth; a homology-only "
     "detection is never called a false positive",
   "no_marker_detected_remains_a_statement_about_the_database": True,
   "no_PlasmidCall_artefact_is_read_or_run": True,
   "no_Paper_2_result_enters_this_module": True,
   "database_change_and_algorithm_change_reported_separately_never_pooled": True},

 "FEASIBILITY_FINDING_THAT_CHANGES_THE_PANEL": {
   "NM0_specified": "re-call with a LATER database and a later tool version",
   "fact_established_2026_08_22": "MOB-suite 3.1.9 (2024-06-04) is the NEWEST release that "
     "exists (GitHub releases API and PyPI, both checked 2026-08-22). The Zenodo concept record "
     "10.5281/zenodo.3785612 holds exactly two database versions: v3.1.8 (record 10304948, "
     "2023-12-08, 450.9 MB) which is the one the frozen run used, and v2.0.0 (record 3785613, "
     "2020-05-04, 449.3 MB).",
   "consequence": "the 'later version' axis is NOT EXECUTABLE. It is replaced by an EARLIER "
     "version axis plus an INDEPENDENT IMPLEMENTATION axis, which is the stronger test.",
   "honest_limitation": "a backward test bounds version sensitivity; it cannot prove stability "
     "against a future release that does not yet exist. C10 must be worded accordingly."},

 "arms": {
   "R_rule_perturbation": {
     "requires_installation": False, "requires_download": False,
     "input": "plasmid_mobility_annotation.tsv, the frozen per-replicon marker calls",
     "what_varies": "only the decision rule mapping observed markers to C/D/E. The marker "
                    "calls themselves are held fixed, so this arm isolates rule dependence "
                    "from detection dependence.",
     "variants_declared_before_running": {
       "R0": "frozen rule. E = relaxase AND mpf; D = (relaxase OR oriT) AND NOT E; "
             "C = neither relaxase nor oriT",
       "R1": "mpf counts as mobilisation-consistent evidence: D = (relaxase OR oriT OR mpf) "
             "AND NOT E. Rationale: a mate-pair-formation system is transfer machinery; "
             "treating it as no evidence at all is one defensible reading, not the only one.",
       "R2": "conjugation requires corroboration: E = relaxase AND mpf AND oriT. Rationale: "
             "the strictest reading of 'combined evidence consistent with conjugation'.",
       "R3": "oriT alone is insufficient for D: D = relaxase AND NOT E; an oriT-only replicon "
             "falls to C. Rationale: oriT is a short sequence feature and the weakest single "
             "line of evidence in the scheme.",
       "R4": "mpf without relaxase is conjugation-consistent: E = mpf AND (relaxase OR oriT). "
             "Rationale: the most permissive defensible reading of class E."},
     "note": "R1-R4 are alternative readings of the SAME frozen wording, chosen to bracket it "
             "from both directions. Two loosen the rule and two tighten it, so the arm cannot "
             "be one-sided by construction."},
   "T_tool_version": {
     "requires_installation": True, "requires_download": "~0 (database reused)",
     "comparison": "mob_typer 3.1.8 + database v3.1.8   VERSUS   mob_typer 3.1.9 + database "
                   "v3.1.8 (the frozen baseline)",
     "isolates": "algorithm behaviour, database held constant",
     "known_confound_declared_in_advance": "the 3.1.9 changelog states it fixed contig-ordering "
       "effects on reproducibility present in 3.1.8. Disagreement attributable to that fix is "
       "expected, is not evidence of instability in the frozen result, and must be reported as "
       "a fix rather than as noise."},
   "D_database_version": {
     "requires_installation": True, "requires_download": "450.9 MB -> ~1.6 GB unpacked",
     "comparison": "mob_typer 3.1.9 + database v2.0.0   VERSUS   mob_typer 3.1.9 + database "
                   "v3.1.8 (the frozen baseline)",
     "isolates": "database content, tool held constant",
     "compatibility_risk_declared_in_advance":
       "a 2020 database may not load in a 3.x tool, and MOB-cluster identifiers changed between "
       "major versions. A 20-replicon smoke test runs FIRST. If it fails, this arm is reported "
       "as NOT EXECUTABLE with the error, and is not silently dropped.",
     "if_not_executable": "C10's database limb is reported as UNTESTABLE WITH AVAILABLE "
                          "ARTEFACTS, never as passed"},
   "I_independent_scheme": {
     "requires_installation": True, "requires_download": "~100 MB",
     "comparison": "an independent relaxase/MPF scheme (CONJscan via MacSyFinder) and an "
                   "independent oriT scheme, against the frozen MOB-suite calls",
     "isolates": "marker-scheme dependence, i.e. whether the biology or the database drives the "
                 "call",
     "reporting": "CONCORDANCE ONLY - a 4x4 table plus Cohen's kappa, exactly as the MGE "
                  "concordance was reported. Never accuracy against truth."}},

 "sampling": {
   "primary_registered_unit": "ARG-bearing plasmid replicon",
   "target_population": "the 6,621 ARG-bearing plasmid replicons in the frozen run",
   "subset_rule": "sort accessions LC_ALL=C, hash each with SHA-256, take the 1,000 with the "
                  "smallest digest (NM-0 rule, applied verbatim)",
   "n_primary": 1000,
   "subset_token_hash": subset_hash,
   "frozen_before_running": True,
   "prespecified_secondary_census": {
     "n": 6621,
     "condition": "run only if the 1,000-replicon primary completes within the reported time "
                  "envelope; declared HERE, before any result, so it can never be a post hoc "
                  "extension chosen because the primary was unfavourable",
     "why_permitted": "all 6,621 FASTA are already local; a census removes sampling error "
                      "entirely and is strictly stronger than a subset"},
   "arm_R_uses_all_6621": "Arm R is arithmetic on already-computed marker calls, so it is run "
                          "as a census over all 6,621 replicons and all 39,209 plasmid-side "
                          "occurrences. No sampling is involved and none is justified."},

 "units_and_denominators_fixed_now": {
   "why_this_section_exists": "NM-0 wrote its 5 pct bound over 'plasmids' while C10's "
                              "denominator field lists BOTH 6,621 plasmids AND 39,209 "
                              "plasmid-borne occurrences. The two are not interchangeable "
                              "because marker-positive plasmids carry more ARGs. Fixing this "
                              "after seeing transitions would be threshold selection.",
   "PRIMARY": "replicons: the proportion of the 6,621 (or of the 1,000 subset, per arm) whose "
              "mobility_category changes. This is the unit NM-0's bound is written over.",
   "CO_PRIMARY_reported_always": "occurrences: the proportion of the 39,209 plasmid-side "
                                 "occurrences whose portability_class changes.",
   "both_reported_in_every_arm": True,
   "the_5_pct_bound_is_applied_to": "the PRIMARY replicon proportion, as a POINT ESTIMATE, "
                                    "because NM-0 wrote a bound and not an interval",
   "intervals_reported_alongside": "replicon-clustered and BioProject-clustered bootstrap, "
                                   "2,000 resamples, seed 20260822, reported for information "
                                   "and never used to move the gate",
   "unresolved_mobility_handling": "mobility_unresolved is a fourth category. Any transition "
     "INTO or OUT OF it counts as a category change. Rationale: the conservative choice; the "
     "alternative would let a tool that simply fails to produce a record score as agreement. "
     "The frozen run has ZERO unresolved records, so this rule can only ever count against a "
     "comparator, never in its favour."},

 "uncertainty": {
   "requirement": "replicon- and BioProject-aware, never occurrence-only",
   "why": "plasmid-side occurrences cluster at 5-6.6 per replicon across 6,621 replicons and "
          "1,346 BioProjects; effective counts by 1/HHI are replicon 3,699 and BioProject 168. "
          "Occurrence-level intervals would be materially overconfident.",
   "method": "cluster bootstrap resampling whole replicons, and separately whole BioProjects, "
             "2,000 resamples, seed 20260822, recomputing every proportion inside each resample",
   "leave_one_out": "leave-one-BioProject-out influence on every headline proportion, reporting "
                    "the maximum relative change, matching the NM-V2 precedent"},

 "GATES": {
   "why_this_section_exists": "NM-V3 is the only NM-0 module with NO gates object. NM-V1 and "
     "NM-V4 each carry success/revise/failure. The triad below is constructed by mapping NM-0's "
     "OWN already-assigned consequences (stop conditions S4 and S7, and the C09/C10 gate text) "
     "onto the triad. No new numeric threshold is invented.",
   "PRIMARY_GATE_C10": {
     "source": "NM0 acceptable_disagreement_bounds and MATRIX C10 field 14, verbatim: '5 pct or "
               "fewer of plasmids change mobility category between database versions, and no "
               "headline direction reverses'",
     "SUCCESS": "PRIMARY replicon category-change proportion <= 5.0 pct in every executable "
                "arm AND no headline direction reverses AND the location layer is unchanged",
     "REVISE": "the change exceeds 5.0 pct in any executable arm, or a headline direction "
               "reverses. Consequence is NOT open: NM-0 S4 already fixed it - 'C10 fails; "
               "classes C/D/E reported as version-dependent and every class figure "
               "version-stamped', the A/B versus C/D/E split remains the stable headline, and "
               "the destination does not change provided the version dependence is disclosed.",
     "FAIL": "ANY change to the documented location layer, or any result that would alter the "
             "74,349 denominator. This is a hard stop under NM-0 S7: the analysis is refused, "
             "not the claim.",
     "no_threshold_invented": "5.0 pct and 0 pct are both carried over verbatim. The triad "
                              "adds no third number."},
   "SECONDARY_GATE_C09": {
     "source": "MATRIX C09 field 14, verbatim: 'the conjugative against marker-negative "
               "difference in three-or-more-drug-class carriage keeps the same direction and a "
               "CI excluding 0 under lineage-clustered resampling and under the alternative "
               "marker rule'",
     "PASS": "direction preserved AND CI excludes 0, under BOTH replicon/BioProject-clustered "
             "resampling AND every Arm R variant",
     "FAIL": "direction reverses or the CI includes 0 -> C09 demoted to Tier 3 descriptive and "
             "removed from the abstract, exactly as MATRIX C09 field 15 prescribes",
     "baseline_being_tested": "conjugative 61.54 pct vs mobilizable 42.94 pct vs "
                              "marker-negative 41.68 pct carrying >=3 drug classes; 1,455 "
                              "high-concern architectures (21.98 pct)"},
   "combination_rule_fixed_now": "C10 and C09 are scored INDEPENDENTLY and either may fail "
     "while the other passes. NM-0 left this unstated; combining them would let a strong result "
     "on one mask a failure on the other.",
   "concentration_of_change": {
     "prespecified_question": "are changes concentrated among low-evidence replicons or do they "
                              "move high-confidence calls",
     "definition_fixed_before_running": "a HIGH-CONFIDENCE call is a replicon with two or more "
       "independent marker types supporting its class (for E: relaxase AND mpf, which is every "
       "E by definition, so E is stratified further by presence of oriT corroboration). A "
       "LOW-EVIDENCE call is a replicon whose class rests on exactly one marker type.",
     "reported_as": "the category-change proportion within each stratum separately, so a low "
                    "overall figure cannot hide movement among high-confidence calls"}},

 "invariance_to_be_DEMONSTRATED_not_assumed": {
   "finding_from_the_constraint_audit": "no A. baumannii claim, species comparison, host/source "
     "comparison or determinant-family enrichment result in the blueprint rests on the C/D/E "
     "split; they rest on the plasmid/chromosome aggregate, which is arithmetically invariant "
     "to any C/D/E reshuffle.",
   "requirement": "this invariance is DEMONSTRATED arithmetically under every Arm R variant and "
     "every executable tool/database arm, not asserted. If any of these quantities moves, that "
     "is itself a finding and means the audit's independence claim was wrong.",
   "quantities_checked": ["occurrence-weighted plasmid share 52.736 pct",
                          "genome-weighted plasmid share 36.586 pct",
                          "block-weighted chromosomal MGE positivity 30.14 pct",
                          "occurrence-weighted chromosomal MGE positivity 46.39 pct",
                          "A. baumannii 14.0 pct plasmid / 80.9 pct chromosomal-MGE",
                          "K. pneumoniae 67.2 pct / 36.7 pct",
                          "the NM-V2 baseline block-weighted ratio 2.4577",
                          "the NM-V4 ab_residual 2.1126"],
   "expected": "all invariant; any movement is a FAIL under the location-layer stop condition "
               "or a defect in this module"},

 "explicitly_out_of_scope": {
   "axes_4_and_5_of_the_NM0_panel": "PARTIAL/INTERNAL_STOP exclusions (-0.32 pp) and plus-scope "
     "inclusion (-6.55 pp) are ALREADY MEASURED and are cited, not re-run. NM-0 lists both as "
     "sensitivity_only.",
   "prespecified_sensitivity_sets_S1_S2_S3": "no code path classifies them; NM-V3 does not "
     "create one. Recorded as a gap, not silently skipped.",
   "AMRFinderPlus_axis": "the determinant-calling axis changes which occurrences are PRIMARY "
     "and therefore the denominator, which NM-0 S7 makes a hard stop. Deferred with its "
     "rationale, not attempted.",
   "A_versus_B": "NM-V1, already passed",
   "PlasmidCall": "prohibited"},

 "outputs": ["NMV3_TRANSITION_MATRICES.tsv - frozen class -> alternative class, per arm",
             "NMV3_HEADLINE_INVARIANCE.tsv", "NMV3_CONCORDANCE.tsv",
             "NMV3_RESULT_RECEIPT.json", "NMV3_RESULT_STATUS.md"],

 "input_digests": digests,
 "deterministic_subset_first_20_by_digest": subset[:20],
    }

    json.dump(Q, open(a.out, "w", encoding="utf-8", newline="\n"), indent=2)
    print("%s\n  wrote %s" % (VERSION, a.out))
    print("  FROZEN DESIGN SHA-256: %s" % sha(a.out))
    print("  subset: 1,000 of %d replicons | subset hash %s" % (len(acc), subset_hash[:32]))


if __name__ == "__main__":
    main()
