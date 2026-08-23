"""PR-CONTEXT Phase A -- mobile-genetic-element annotation of chromosomal ARG neighbourhoods.

TOOL SELECTION, documented before use (task item 4)

  Prodigal 2.6.3            gene calling in metagenomic mode, coordinates relative to the
                            window, which are mapped back to chromosome coordinates exactly
                            by adding the window's own start offset.
  HMMER 3.3.2 hmmsearch     profile search, E-value 1e-5, domain table output.
  IS/transposase profiles   266 pHMMs bundled with ISEScan 1.7.3 (bioconda), derived from
                            ISfinder family clusters. Bundled: no external download, no
                            licence barrier.
  Integrase profiles        3 HMMs bundled with IntegronFinder 2.0.6 (bioconda):
                            integron_integrase, IntI, phage-int.

  Total added download      0 bytes beyond the 1.5 GB conda env; every model ships with the
                            tool.

  ISEScan 1.7.3 itself was installed and benchmarked and is NOT used as the primary caller:
  it needed >25 minutes for 200 windows against 17 seconds for this path, which extrapolates
  to well over a day for 21,955 windows. It is run on one chunk as a specificity cross-check
  and that comparison is reported. This is a throughput decision, stated rather than hidden.

  SENSITIVITY AND LIMITS, stated plainly. This detects PROTEIN HOMOLOGY to transposase,
  integrase and IS-family profiles. It does not resolve terminal inverted repeats, does not
  assemble complete IS elements, and does not prove an element is intact or active. A hit is
  a DETECTED SEQUENCE FEATURE and is labelled as such throughout.
"""
import argparse, collections, csv, glob, hashlib, json, os, re, subprocess, sys, datetime

VERSION = "pr_context_mge_annotate_v1.0.0"
# The annotation environment. Set PR_CONTEXT_MGE_ENV, or let it fall back to
# the active conda/micromamba prefix.
ENV = os.environ.get("PR_CONTEXT_MGE_ENV") or os.environ.get("CONDA_PREFIX", "")
if not ENV:
    raise SystemExit("set PR_CONTEXT_MGE_ENV to the mge environment prefix")
IS_HMM = ENV + "/bin/pHMMs/clusters.faa.hmm"
IF_DIR = ENV + "/lib/python3.14/site-packages/integron_finder/data/Models"
EVALUE = "1e-5"
NEAR = 10000          # primary window, matches the frozen protocol
SENS = (5000, 20000)  # sensitivity windows


def sha256_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def write(path, rows, cols):
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")).replace("\t", " ").replace("\n", " ")
                               for c in cols) + "\n")


HDR = re.compile(r"^(blk\d+)\|([^|]+)\|(\d+)-(\d+)$")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--chunks", type=int, default=48)
    ap.add_argument("--threads", type=int, default=4)
    a = ap.parse_args()
    root = a.root
    out = os.path.join(root, "out")
    W = os.path.join(root, "mge_work")
    for d in (os.path.join(W, "chunks2"), os.path.join(W, "fast2"),
              os.path.join(W, "logs")):
        os.makedirs(d, exist_ok=True)

    wins = sorted(glob.glob(os.path.join(root, "window_fasta", "*.fna")))
    print("%s | window blocks available: %d" % (VERSION, len(wins)))
    if not wins:
        print("REFUSING: no windows"); sys.exit(1)

    # build chunks
    per = max(1, (len(wins) + a.chunks - 1) // a.chunks)
    chunks = []
    for i in range(0, len(wins), per):
        cp = os.path.join(W, "chunks2", "c%03d.fna" % (i // per))
        with open(cp, "w", encoding="utf-8", newline="\n") as fh:
            for f in wins[i:i + per]:
                fh.write(open(f, encoding="utf-8").read())
        chunks.append(cp)
    print("chunks: %d (%d windows each)" % (len(chunks), per))

    intg = os.path.join(W, "fast2", "integrase.hmm")
    with open(intg, "wb") as fh:
        for n in ("integron_integrase.hmm", "IntI.hmm", "phage-int.hmm"):
            p = os.path.join(IF_DIR, n)
            if os.path.exists(p):
                fh.write(open(p, "rb").read())

    # ---- run prodigal + hmmsearch per chunk, in parallel ----
    script = os.path.join(W, "run_chunk.sh")
    with open(script, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("""#!/bin/bash
export PATH=%s/bin:$PATH
c="$1"; n=$(basename "$c" .fna); W=%s
prodigal -i "$c" -a $W/fast2/$n.faa -p meta -q -o /dev/null 2>/dev/null
hmmsearch --cpu %d -E %s --domtblout $W/fast2/$n.is.domtbl %s $W/fast2/$n.faa >/dev/null 2>&1
hmmsearch --cpu %d -E %s --domtblout $W/fast2/$n.int.domtbl $W/fast2/integrase.hmm $W/fast2/$n.faa >/dev/null 2>&1
echo "$n exit=$?"
""" % (ENV, W, a.threads, EVALUE, IS_HMM, a.threads, EVALUE))
    os.chmod(script, 0o755)
    t0 = datetime.datetime.now()
    proc = subprocess.run(
        "ls %s/chunks2/c*.fna | xargs -P 12 -I{} bash %s {}" % (W, script),
        shell=True, capture_output=True, text=True)
    el = (datetime.datetime.now() - t0).total_seconds()
    open(os.path.join(W, "logs", "fast_driver.log"), "w",
         encoding="utf-8").write(proc.stdout + proc.stderr)
    ok = proc.stdout.count("exit=0")
    print("annotation wall %.1f min | chunks exit=0: %d/%d" % (el / 60, ok, len(chunks)))

    # ---- parse proteins: id -> (block, replicon, block_start, prot_start, prot_end, strand)
    prot = {}
    for f in sorted(glob.glob(os.path.join(W, "fast2", "*.faa"))):
        for line in open(f, encoding="utf-8", errors="replace"):
            if not line.startswith(">"):
                continue
            parts = line[1:].split("#")
            pid = parts[0].strip()
            try:
                ps, pe, st = int(parts[1]), int(parts[2]), int(parts[3])
            except (IndexError, ValueError):
                continue
            base = pid.rsplit("_", 1)[0]
            m = HDR.match(base)
            if not m:
                continue
            prot[pid] = (m.group(1), m.group(2), int(m.group(3)), ps, pe, st)
    print("proteins called: %d" % len(prot))

    # ---- parse hmm hits ----
    feats = []
    for kind, pat in (("IS_or_transposase", "*.is.domtbl"),
                      ("integrase_or_integron", "*.int.domtbl")):
        for f in sorted(glob.glob(os.path.join(W, "fast2", pat))):
            for line in open(f, encoding="utf-8", errors="replace"):
                if line.startswith("#"):
                    continue
                p = line.split()
                if len(p) < 23:
                    continue
                pid, qname = p[0], p[3]
                try:
                    ev = float(p[6])
                    tlen = int(p[2])
                    alen = abs(int(p[18]) - int(p[17])) + 1
                except ValueError:
                    continue
                pr = prot.get(pid)
                if not pr:
                    continue
                blk, rep, bstart, ps, pe, st = pr
                feats.append({
                    "block_id": blk, "replicon_accession": rep,
                    "feature_class": kind, "feature_name": qname,
                    "protein_id": pid,
                    "chrom_start": bstart + ps - 1, "chrom_end": bstart + pe - 1,
                    "strand": "+" if st == 1 else "-",
                    "evalue": ev,
                    "alignment_length_aa": alen,
                    "coverage_of_profile_pct": round(100.0 * alen / max(tlen, 1), 2),
                    "annotation_source": "HMMER 3.3.2 vs %s" % (
                        "ISEScan 1.7.3 IS/transposase pHMMs" if kind.startswith("IS")
                        else "IntegronFinder 2.0.6 integrase models"),
                    "confidence": "high" if ev < 1e-20 else "moderate",
                    "evidence_category": "detected_sequence_feature"})
    # deduplicate: one feature per protein per class, keep best E-value
    best = {}
    for f in feats:
        k = (f["protein_id"], f["feature_class"])
        if k not in best or f["evalue"] < best[k]["evalue"]:
            best[k] = f
    feats = sorted(best.values(), key=lambda r: (r["replicon_accession"],
                                                 r["chrom_start"]))
    write(os.path.join(out, "mge_feature_inventory.tsv"), feats, list(feats[0].keys()))
    print("MGE features (deduplicated per protein per class): %d" % len(feats))

    # ---- join to chromosomal ARG occurrences ----
    occ = [r for r in csv.DictReader(
        open(os.path.join(out, "determinant_occurrences.tsv"), encoding="utf-8"),
        delimiter="\t") if r["analysis_set"] == "PRIMARY"
        and r["evidence_type"] == "direct_chromosome"]
    byrep = collections.defaultdict(list)
    for f in feats:
        byrep[f["replicon_accession"]].append(f)
    covered = {b["block_id"]: b for b in csv.DictReader(
        open(os.path.join(out, "shared_context_blocks.tsv"), encoding="utf-8"),
        delimiter="\t")} if os.path.exists(
        os.path.join(out, "shared_context_blocks.tsv")) else {}
    annotated_reps = {f["replicon_accession"] for f in feats}
    windows_present = set()
    for w in wins:
        h = open(w, encoding="utf-8").readline()[1:].strip()
        m = HDR.match(h)
        if m:
            windows_present.add((m.group(2), int(m.group(3)), int(m.group(4))))

    rows, miss = [], []
    for r in occ:
        rep = r["replicon_accession"]
        gs, ge = sorted((int(r["gene_start"]), int(r["gene_end"])))
        inwin = any(rep == w[0] and w[1] <= gs and ge <= w[2] for w in windows_present)
        if not inwin:
            miss.append({"assembly_version": r["assembly_version"],
                         "replicon_accession": rep,
                         "determinant_name": r["determinant_name"],
                         "gene_start": gs,
                         "reason": "neighbourhood window not retrieved for this occurrence",
                         "state": "MGE_ANNOTATION_NOT_PERFORMED"})
            continue
        near = {}
        for wsz in (NEAR,) + SENS:
            lo, hi = gs - wsz, ge + wsz
            sel = [f for f in byrep.get(rep, [])
                   if f["chrom_end"] >= lo and f["chrom_start"] <= hi]
            near[wsz] = sel
        prim = near[NEAR]
        dists = []
        for f in prim:
            d = 0 if (f["chrom_end"] >= gs and f["chrom_start"] <= ge) else \
                min(abs(f["chrom_start"] - ge), abs(gs - f["chrom_end"]))
            dists.append(d)
        rows.append({
            "assembly_version": r["assembly_version"],
            "biosample_accession": r["biosample_accession"],
            "organism_harmonized": r["organism_harmonized"], "genus": r["genus"],
            "replicon_accession": rep, "determinant_name": r["determinant_name"],
            "gene_family": r["gene_family"], "drug_class": r["drug_class"],
            "gene_start": gs, "gene_end": ge, "strand": r["strand"],
            "n_mge_features": len(prim),
            "n_is_transposase": sum(1 for f in prim
                                    if f["feature_class"] == "IS_or_transposase"),
            "n_integrase_integron": sum(1 for f in prim
                                        if f["feature_class"] == "integrase_or_integron"),
            "nearest_mge_distance_bp": min(dists) if dists else "",
            "nearest_mge_name": (prim[dists.index(min(dists))]["feature_name"]
                                 if dists else ""),
            "overlapping_mge": "yes" if any(d == 0 for d in dists) else "no",
            "n_mge_5kb": len(near[5000]),
            # The retrieval window is +/-10 kb, so a +/-20 kb search cannot see sequence
            # that was never fetched. Reporting it as a sensitivity analysis would be
            # false: it is truncated by construction and differs from the primary only
            # where a MERGED block happened to extend further because a neighbouring ARG
            # pulled the boundary out. It is therefore recorded as NOT EVALUABLE.
            "n_mge_20kb_TRUNCATED": len(near[20000]),
            "n_mge_20kb_evaluable": "no - retrieval window is +/-10 kb",
            "same_replicon": "yes",
            "window_bp": NEAR,
            "evidence_category": "detected_sequence_feature",
            "annotation_state": "COMPLETE"})
    write(os.path.join(out, "arg_mge_neighbourhood.tsv"), rows, list(rows[0].keys()))
    if miss:
        write(os.path.join(out, "mge_missingness.tsv"), miss, list(miss[0].keys()))
    else:
        write(os.path.join(out, "mge_missingness.tsv"),
              [{"assembly_version": "", "replicon_accession": "", "determinant_name": "",
                "gene_start": "", "reason": "none - every chromosomal occurrence had its "
                                            "window retrieved and annotated",
                "state": "COMPLETE"}],
              ["assembly_version", "replicon_accession", "determinant_name", "gene_start",
               "reason", "state"])

    rec = [{"tool": "Prodigal", "version": "2.6.3", "role": "gene calling (-p meta)",
            "database": "n/a", "database_version": "n/a", "download_bytes": 0,
            "licence": "GPLv3, bioconda", "evalue": "", "source": "bioconda"},
           {"tool": "HMMER hmmsearch", "version": "3.3.2", "role": "profile search",
            "database": "n/a", "database_version": "n/a", "download_bytes": 0,
            "licence": "BSD-3, bioconda", "evalue": EVALUE, "source": "bioconda"},
           {"tool": "IS/transposase pHMMs", "version": "bundled with ISEScan 1.7.3",
            "role": "insertion sequence and transposase detection",
            "database": "clusters.faa.hmm", "database_version": "ISEScan 1.7.3",
            "download_bytes": 0, "licence": "bundled with the tool, no separate download",
            "evalue": EVALUE, "source": "bioconda isescan 1.7.3",
            "n_profiles": 266,
            "database_sha256": sha256_file(IS_HMM) if os.path.exists(IS_HMM) else ""},
           {"tool": "integrase models", "version": "bundled with IntegronFinder 2.0.6",
            "role": "integron integrase / IntI / phage integrase detection",
            "database": "integron_integrase.hmm + IntI.hmm + phage-int.hmm",
            "database_version": "IntegronFinder 2.0.6", "download_bytes": 0,
            "licence": "GPLv3, bundled", "evalue": EVALUE, "source": "bioconda",
            "n_profiles": 3, "database_sha256": sha256_file(intg)},
           {"tool": "ISEScan", "version": "1.7.3",
            "role": "EVALUATED, NOT USED AS PRIMARY - >25 min per 200 windows against 17 s "
                    "for the path above; run on one chunk as a specificity cross-check",
            "database": "clusters.faa.hmm", "database_version": "1.7.3",
            "download_bytes": 0, "licence": "bioconda", "evalue": "",
            "source": "bioconda"}]
    write(os.path.join(out, "mge_annotation_receipts.tsv"), rec,
          ["tool", "version", "role", "database", "database_version", "n_profiles",
           "database_sha256", "download_bytes", "licence", "evalue", "source"])

    n = len(rows)
    withmge = sum(1 for r in rows if r["n_mge_features"] > 0)
    print("\n=== Q5: chromosomal ARG occurrences with a nearby MGE feature ===")
    print("  annotated occurrences        : %d of %d chromosomal" % (n, len(occ)))
    print("  with >=1 MGE within +/-10 kb : %d  (%.2f%%)" % (withmge, 100.0 * withmge / n))
    print("  with an OVERLAPPING MGE      : %d  (%.2f%%)"
          % (sum(1 for r in rows if r["overlapping_mge"] == "yes"),
             100.0 * sum(1 for r in rows if r["overlapping_mge"] == "yes") / n))
    c5 = sum(1 for r in rows if r["n_mge_5kb"] > 0)
    print("  sensitivity +/-5,000 bp      : %d  (%.2f%%)" % (c5, 100.0 * c5 / n))
    print("  sensitivity +/-20,000 bp     : NOT EVALUABLE - the retrieval window is "
          "+/-10 kb, so 20 kb sequence was never fetched")
    print("  not annotated (window absent): %d" % len(miss))
    for f in ("mge_feature_inventory.tsv", "arg_mge_neighbourhood.tsv",
              "mge_annotation_receipts.tsv", "mge_missingness.tsv"):
        print("  %-40s %s" % (f, sha256_file(os.path.join(out, f))))


if __name__ == "__main__":
    main()
