"""PR-CONTEXT step 1 -- freeze the scientific questions and every rule, BEFORE any aggregate.

At the moment this runs, the determinant-to-replicon join has NOT been performed. No
plasmid/chromosome proportion, no enrichment statistic and no portability class exists. What
IS known is structural: the row counts by AMRFinderPlus category, which identifier fields
parse, and that the join is possible. Those are properties of the file format, not results,
and they were established by the read-only audit that this protocol cites.

The one number this protocol must get right before anything else is the DENOMINATOR. The
receipt figure of 184,538 is every AMRFinderPlus row -- metal, biocide, heat, virulence and
point mutations included. Treating it as an acquired-ARG count would inflate the denominator
by a factor of about 2.2 and would make every proportion in the manuscript wrong.
"""
import argparse, datetime, hashlib, json, os, sys

VERSION = "pr_context_freeze_protocol_v1.0.0"


def sha256_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


P = {
 "protocol": "FROZEN_PORTABILITY_CONTEXT_PROTOCOL",
 "version": "1.0.0",
 "programme": "PortabilityRisk",
 "independent_of": "PlasmidCall. No PlasmidCall artefact is read, run, validated or modified.",

 # ------------------------------------------------------------------ inputs
 "inputs": {
   "amrfinder_outputs": {
     "path": "audit/data/derived/g_c_processing/amrfinder",
     "n_tables": 7216, "n_rows_total": 184538,
     "aggregate_sha256":
       "b891eb821a08357bef74d003c5d6e691ebd8b9477fcde4d374270b10e271a7df",
     "note": "Frozen gate G-C output. NOT rerun. Read only."},
   "worklist": {"file": "determinant_processing_worklist_g_c_v2.tsv",
                "sha256": "2d5e6289af5c3aeb98680b542c86a660f740f5c5f1de3c0c1fec52534879dfe3"},
   "toolchain_that_produced_them": {
     "name": "AMRFinderPlus", "version": "4.2.7", "database_version": "2026-08-07.1",
     "database_dir_sha256":
       "e42fcd2cdd0f2ce493da790b8ebc68f33068f1efb62afd0530fa643c92567152"}},

 # ------------------------------------------------------------------ denominator
 "denominator_correction": {
   "observed_total_rows": 184538,
   "breakdown": {"AMR/AMR": 85507, "AMR/POINT": 9391, "AMR/POINT_DISRUPT": 2033,
                 "STRESS/METAL": 65355, "STRESS/BIOCIDE": 5544, "STRESS/HEAT": 5297,
                 "STRESS/ACID": 3, "VIRULENCE/VIRULENCE": 11408},
   "statement":
     "184,538 is every AMRFinderPlus row. It is NOT an acquired-ARG count and is never "
     "reported as one. 85,507 rows are Type=AMR and Subtype=AMR."},

 "acquired_arg_definition": {
   "primary_rule":
     "Type == AMR AND Subtype == AMR AND Scope == core AND Class != EFFLUX",
   "primary_n_expected": 74349,
   "why_this_rule":
     "It is the project's already-frozen dev_parity_core_rule, used for PlasmidCall's ARG "
     "domain. Reusing it keeps the two programmes commensurable and means the definition "
     "was not chosen for this analysis. Scope=core selects genes whose primary function is "
     "resistance; EFFLUX is excluded because efflux pumps in the plus set are "
     "predominantly chromosomal core machinery, so including them would load the "
     "chromosomal arm with genes that were never acquired and bias every proportion.",
   "prespecified_sensitivity_sets": {
     "S1_core_plus_nonefflux": "Type=AMR AND Subtype=AMR AND Class != EFFLUX (n 76,383)",
     "S2_all_amr_amr": "Type=AMR AND Subtype=AMR, EFFLUX included (n 85,507)",
     "S3_complete_calls_only":
       "primary set minus Method containing PARTIAL and minus INTERNAL_STOP"},
   "excluded_and_why": {
     "AMR/POINT and AMR/POINT_DISRUPT":
       "point mutations, not acquired genes; and G_C_AMENDMENT_001 records that they were "
       "only searched where an AMRFinderPlus organism flag existed, so their denominator is "
       "2,542 genomes short. Never counted as acquired ARGs and never pooled with them.",
     "STRESS/*": "metal, biocide, heat and acid tolerance are not antimicrobial resistance",
     "VIRULENCE/*": "not resistance"}},

 # ------------------------------------------------------------------ units
 "units": {
   "primary_biological_unit":
     "one acquired ARG OCCURRENCE = one qualifying AMRFinderPlus row, identified by "
     "(assembly_accession, sequence_accession, gene_start, gene_end, strand, gene_symbol)",
   "genomic_clustering_unit": "assembly accession; BioSample asserted 1:1 and verified",
   "replicon_unit": "unique sequence accession WITHIN a versioned assembly",
   "primary_denominator":
     "acquired ARG occurrences with directly evidenced, unambiguous replicon origin",
   "missingness_denominator": "ALL acquired ARG occurrences, resolved or not",
   "primary_analysis_population": "the frozen 7,216 complete genomes"},

 "handling_rules": {
   "repeated_genes_same_replicon":
     "every occurrence is preserved as its own row. A prespecified sensitivity collapses to "
     "distinct (replicon, gene_symbol) to show the effect of copy number.",
   "duplicated_assemblies_or_biosamples":
     "the G-B frame already deduplicated under rules F1-F4 and only dedup_status=unique "
     "entered the worklist. Re-verified here; any duplicate found is reported, not dropped.",
   "ambiguous_sequence_identifiers":
     "a contig id matching more than one replicon record, or none, is labelled "
     "identifier_ambiguous or identifier_unmatched and RETAINED in the missingness "
     "denominator. Never silently dropped.",
   "pseudogenes_and_partial_genes":
     "Method containing PARTIAL, and INTERNAL_STOP, are retained in the primary set and "
     "flagged. Removing them after seeing results would be selection on the outcome. "
     "Sensitivity S3 excludes them.",
   "minimum_metadata":
     "an occurrence needs only assembly accession, sequence accession and coordinates to "
     "enter the primary analysis. Geography, host, isolation source and collection date are "
     "reported as COVERAGE and never imputed; claims that need them are restricted to the "
     "subset that has them."},

 # ------------------------------------------------------------------ questions
 "frozen_questions": {
   "Q1": "Among directly evidenced acquired ARG occurrences, what proportion are on "
         "plasmids versus chromosomes?",
   "Q2": "Which resistance gene families show significant plasmid enrichment?",
   "Q3": "Which determinants occur in BOTH plasmid and chromosome contexts?",
   "Q4": "Among plasmid-associated ARG occurrences, what proportion are on replicons with "
         "genomic evidence of mobilization or conjugation?",
   "Q5": "Among chromosomal ARG occurrences, what proportion have nearby MGE signatures?",
   "Q6": "How do these distributions differ by organism, determinant class and available "
         "epidemiological metadata?",
   "Q7": "Does a context-aware portability classification distinguish patterns invisible to "
         "gene-presence analysis alone?"},

 # ------------------------------------------------------------------ evidence
 "evidence_rules": {
   "direct_closed_replicon_requires_all_of": [
     "the genome is in the frozen complete-genome cohort",
     "the AMRFinderPlus contig id maps EXACTLY to one documented replicon record",
     "the replicon is explicitly designated Chromosome or Plasmid by NCBI",
     "the gene coordinates lie within that exact sequence's length"],
   "categories": ["direct_chromosome", "direct_plasmid", "direct_other",
                  "replicon_unclassified", "identifier_ambiguous", "identifier_unmatched",
                  "coordinates_missing"],
   "prohibited_inferences": [
     "inferring chromosome because no plasmid annotation was found",
     "inferring plasmid from a gene name, an Inc marker, a sequence length or an "
     "AMRFinderPlus label",
     "imputing a location for any unresolved occurrence"],
   "authority": "NCBI assigned_molecule_location_type from the assembly's sequence report. "
                "Nothing is predicted."},

 # ------------------------------------------------------------------ neighbourhoods
 "neighbourhood_rules": {
   "primary_window_bp": 10000,
   "sensitivity_windows_bp": [5000, 20000],
   "circular_wrap":
     "when a replicon is circular the window WRAPS and is recorded as wrapped, not "
     "truncated. Truncation is recorded only for linear replicons or where topology is "
     "unknown.",
   "topology_source": "NCBI nuccore esummary; recorded per replicon, never assumed",
   "overlapping_neighbourhoods":
     "each ARG keeps its own record AND the shared context block is recorded separately, so "
     "co-located ARGs are never counted as independent contexts"},

 # ------------------------------------------------------------------ statistics
 "statistics": {
   "effect_size": "odds ratio, plasmid vs chromosome, per gene family",
   "interval": "95% confidence interval",
   "test": "two-sided Fisher exact for 2x2 family-vs-rest tables",
   "multiple_testing": "Benjamini-Hochberg FDR, q < 0.05, applied across gene families "
                       "tested in the primary analysis",
   "minimum_family_size_for_testing": 20,
   "non_independence": [
     "genome-level aggregation: each genome contributes at most one event per family",
     "cluster bootstrap by BioProject study group, 2000 resamples",
     "cluster bootstrap by genome, 2000 resamples"],
   "non_independence_statement":
     "multiple genes from one genome are NOT independent replicates. Every headline effect "
     "is reported with a cluster-adjusted interval alongside the naive one.",
   "prespecified_strata": ["organism group (genus)", "gene family", "drug class",
                           "collection geography where present",
                           "collection period where present",
                           "clinical vs non-clinical source where present"],
   "sparse_metadata_rule":
     "if a stratum's metadata coverage is below 50% of the primary denominator, coverage is "
     "reported and the claim is restricted to the covered subset, never generalised."},

 "prohibited_claims": [
   "describing portability classes A-E as measured transmission probability",
   "claiming experimentally demonstrated mobility from sequence annotation",
   "calling a plasmid non-mobilizable because a marker database returned no match",
   "calling a plasmid conjugative from a single isolated transposase",
   "reporting 184,538 as an acquired-ARG count",
   "pooling point mutations with acquired genes",
   "overwriting a direct replicon assignment with any model prediction"],

 "future_model_layer": {
   "reserved_field": "plasmidcall_predicted_location",
   "status": "SCHEMA ONLY. Not run, not populated, in this task.",
   "rule": "a future prediction is a SEPARATE evidence layer. It may never overwrite a "
           "direct assignment, and may never form the target label for the model that "
           "produced it."},

 "stop_conditions": [
   "cumulative new download exceeds 100 GB without owner approval",
   "the AMRFinderPlus aggregate digest does not match the frozen value",
   "acquired-ARG count under the primary rule differs from 74,349",
   "main and independent verifier disagree on any principal count or headline result"],
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    if os.path.exists(a.out):
        print("REFUSING: %s exists; a frozen protocol is never overwritten" % a.out)
        sys.exit(1)
    P["created_utc"] = datetime.datetime.now(
        datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    P["created_before_any_join_or_aggregate"] = True
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    with open(a.out, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(P, fh, indent=2)
    print("%s\nwrote %s\nFROZEN PROTOCOL SHA-256: %s" % (VERSION, a.out, sha256_file(a.out)))
    print("primary acquired-ARG denominator declared: %d" % P["acquired_arg_definition"]
          ["primary_n_expected"])


if __name__ == "__main__":
    main()
