"""NM-V4 freeze -- fix the design BEFORE any outcome is computed.

The protocol's freeze order is binding:
  1) fix the inclusion floor
  2) enumerate the confirmation species
  3) freeze and hash the index definition
  4) only then compute the index on the confirmation set

This script performs steps 1-3 and nothing else. It reads only EXPOSURE quantities -- genome
counts and chromosomal ARG occurrence counts -- and deliberately computes neither the plasmid
fraction nor the chromosomal MGE fraction, because those are the two components of the
outcome. A reader can verify from this file alone that the confirmation species were chosen
without reference to the answer.

Discovery species are Acinetobacter baumannii and Klebsiella pneumoniae: the two from which
the discordance principle was articulated. They are excluded from the fit by construction.
"""
import argparse, collections, csv, datetime, hashlib, json, os, sys

VERSION = "nmv4_freeze_v1.0.0"
PROTOCOL_SHA = "c968fb6d16a528a64d064d6f8bbac745390804df4ad0897c09b12de84ca3fbff"
DISCOVERY = ["Acinetobacter baumannii", "Klebsiella pneumoniae"]
FLOOR_GENOMES = 40
FLOOR_CHROM_OCC = 200


def sha256_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", required=True)
    ap.add_argument("--protocol", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    if sha256_file(a.protocol) != PROTOCOL_SHA:
        print("REFUSING: protocol digest mismatch"); sys.exit(1)
    if os.path.exists(a.out):
        print("REFUSING: %s exists; a frozen design is never overwritten" % a.out); sys.exit(1)
    P = json.load(open(a.protocol, encoding="utf-8"))
    if P["status"] != "PROPOSED - OWNER APPROVAL REQUIRED":
        print("REFUSING: unexpected protocol status"); sys.exit(1)
    print("%s | protocol v%s verified %s" % (VERSION, P["version"], PROTOCOL_SHA[:16]))

    O = os.path.join(a.repo, "audit/data/derived/pr_context/out")
    gpath = os.path.join(O, "genome_level_summary.tsv")
    npath = os.path.join(O, "arg_mge_neighbourhood.tsv")

    # ---- EXPOSURE ONLY: genome counts and chromosomal ARG occurrence counts ----
    genomes = collections.Counter()
    bioprojects = collections.defaultdict(set)
    for r in csv.DictReader(open(gpath, encoding="utf-8"), delimiter="\t"):
        genomes[r["organism"]] += 1
        bioprojects[r["organism"]].add(r["bioproject_accession"])
    chrom_occ = collections.Counter()
    for r in csv.DictReader(open(npath, encoding="utf-8"), delimiter="\t"):
        chrom_occ[r["organism_harmonized"]] += 1

    allsp = sorted(set(genomes) | set(chrom_occ))
    eligible = [s for s in allsp
                if genomes.get(s, 0) >= FLOOR_GENOMES
                and chrom_occ.get(s, 0) >= FLOOR_CHROM_OCC]
    confirmation = [s for s in eligible if s not in DISCOVERY]

    print("\n  species observed              : %d" % len(allsp))
    print("  floor: >=%d genomes AND >=%d chromosomal ARG occurrences"
          % (FLOOR_GENOMES, FLOOR_CHROM_OCC))
    print("  species meeting the floor     : %d" % len(eligible))
    print("  discovery species (excluded)  : %s" % ", ".join(DISCOVERY))
    print("  CONFIRMATION SET              : %d species" % len(confirmation))
    print()
    print("  %-34s %8s %8s %8s" % ("species", "genomes", "chromOcc", "bioproj"))
    for s in eligible:
        tag = "  [DISCOVERY]" if s in DISCOVERY else ""
        print("  %-34s %8d %8d %8d%s"
              % (s, genomes.get(s, 0), chrom_occ.get(s, 0), len(bioprojects.get(s, set())), tag))

    D = {
      "design": "NM_V4_FROZEN_DESIGN",
      "version": "1.0.0",
      "builder": VERSION,
      "frozen_utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
      "authorised_by": "owner, D-NM1 approved 2026-08-21",
      "protocol_sha256": PROTOCOL_SHA,
      "frozen_before_any_outcome_was_computed": True,
      "what_this_file_deliberately_does_not_contain":
        "the plasmid-borne fraction and the chromosomal MGE fraction. Those are the two "
        "components of the outcome. Only genome counts and chromosomal ARG occurrence counts "
        "were read to build this design.",

      "inclusion_floor": {"min_genomes": FLOOR_GENOMES,
                          "min_chromosomal_arg_occurrences": FLOOR_CHROM_OCC,
                          "source": "genome_level_summary.tsv and arg_mge_neighbourhood.tsv"},
      "discovery_species": DISCOVERY,
      "confirmation_species": confirmation,
      "n_confirmation_species": len(confirmation),
      "exposure_counts": {s: {"genomes": genomes.get(s, 0),
                              "chromosomal_arg_occurrences": chrom_occ.get(s, 0),
                              "bioprojects": len(bioprojects.get(s, set()))}
                          for s in eligible},

      "index_definition": {
        "per_species_plasmid_fraction_P":
          "occurrence-weighted: n_plasmid_args / n_arg_occurrences, from genome_level_summary.tsv",
        "per_species_chromosomal_MGE_fraction_M":
          "BLOCK-weighted: blocks carrying >=1 MGE feature / all context blocks on that "
          "species chromosomes. Block-weighted, not occurrence-weighted, because "
          "marker-positive blocks carry 2.46 ARGs against 1.23 and occurrence weighting would "
          "manufacture discordance.",
        "discordance_index_D": "D = logit(M) - logit(P)",
        "logit_guard": "proportions are clamped to [0.5/n, 1 - 0.5/n] with n the species "
                       "denominator, a Haldane-style guard applied identically to every "
                       "species and declared here before any value is seen"},

      "tests": {
        "T1_general_principle":
          "leave-one-species-out within the CONFIRMATION set: regress logit(M) on logit(P) "
          "over confirmation species minus one, predict the held-out species, and test whether "
          "residuals are systematically larger than sampling noise. If P alone predicted M, "
          "residuals would be null.",
        "T2_discovery_species_as_holdout":
          "fit logit(M) ~ logit(P) on the CONFIRMATION set only -- species never used to "
          "articulate the principle -- then predict Acinetobacter baumannii from its plasmid "
          "fraction alone and report the residual. This is the non-circular confirmation: the "
          "relationship is learned without the discovery species and then applied to it.",
        "T3_low_plasmid_control":
          "Pseudomonas aeruginosa sits in the confirmation set and has a plasmid fraction close "
          "to Acinetobacter baumannii. If D were a mechanical consequence of a low plasmid "
          "fraction the two would score alike. They must not."},

      "uncertainty": {
        "method": "BioProject-clustered bootstrap within species: resample BioProjects with "
                  "replacement inside each species, recompute P and M from the resampled "
                  "genomes, refit, 2000 resamples",
        "seed": 20260821,
        "why": "genomes within a BioProject are not independent"},

      "prespecified_direction":
        "a species with a low plasmid fraction AND an AbaR-like chromosomal architecture is "
        "predicted to have POSITIVE D. A low plasmid fraction alone is NOT predicted to "
        "produce positive D.",

      "gates": {
        "success": "the T2 residual for Acinetobacter baumannii is positive with a "
                   "BioProject-clustered bootstrap CI excluding 0, AND T3 holds "
                   "(Pseudomonas aeruginosa does not score like Acinetobacter baumannii), "
                   "AND T1 shows residuals systematically non-zero",
        "revise": "the residual is in the predicted direction but its CI includes 0. Report "
                  "the principle as supported-but-underpowered, name the number of species "
                  "that would be needed, and demote it from the abstract to the results",
        "failure": "the residual is null or reverses direction, or Pseudomonas aeruginosa "
                   "scores like Acinetobacter baumannii. Withdraw the general principle, "
                   "demote to Tier 3, retain the Acinetobacter baumannii observation as a "
                   "species-level finding, and change destination"},

      "input_digests": {"genome_level_summary.tsv": sha256_file(gpath),
                        "arg_mge_neighbourhood.tsv": sha256_file(npath),
                        "shared_context_blocks.tsv": sha256_file(
                            os.path.join(O, "shared_context_blocks.tsv")),
                        "mge_feature_inventory.tsv": sha256_file(
                            os.path.join(O, "mge_feature_inventory.tsv"))},
    }
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    json.dump(D, open(a.out, "w", encoding="utf-8", newline="\n"), indent=2)
    print("\n  wrote %s" % a.out)
    print("  FROZEN DESIGN SHA-256: %s" % sha256_file(a.out))


if __name__ == "__main__":
    main()
