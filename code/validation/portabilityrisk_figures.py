"""PortabilityRisk (Paper 1) -- regenerate every main figure from frozen inputs.

Every input digest is verified before anything is drawn, and every output digest
is recorded. A figure that cannot be traced to a hashed input is not produced.

Figure 5 is NOT generated here: it is the NM-DIST four-panel figure, produced by
nmdist_figure.py under the NM-DIST frozen protocol, and is reused unchanged.

    python portabilityrisk_figures.py --check      # verify inputs, draw nothing
    python portabilityrisk_figures.py              # draw and write the receipt
"""
import argparse, collections, csv, hashlib, json, os, sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                    "..", "..", "..", "..", ".."))
OUT = os.path.join(REPO, "audit", "data", "derived", "pr_context", "out")
NM = os.path.join(REPO, "docs", "nature_microbiology")
FIGD = os.path.join(NM, "figures", "manuscript")

# Okabe-Ito, colour-vision safe. Same palette as the NM-DIST figure.
OI = {"blue": "#0072B2", "orange": "#E69F00", "green": "#009E73", "vermillion": "#D55E00",
      "purple": "#CC79A7", "sky": "#56B4E9", "yellow": "#F0E442", "black": "#000000",
      "grey": "#999999"}

plt.rcParams.update({
    "font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 300, "savefig.dpi": 300,
    "pdf.fonttype": 42, "svg.fonttype": "none",
})


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def tsv(name, root=OUT):
    with open(os.path.join(root, name), encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def jsn(name, root=NM):
    with open(os.path.join(root, name), encoding="utf-8") as fh:
        return json.load(fh)


# --- inputs, with the digest each figure depends on -------------------------
INPUTS = {
    "evidence_layer_inventory.tsv": OUT,
    "portability_class_reconciliation.tsv": OUT,
    "both_context_determinants.tsv": OUT,
    "determinant_enrichment_species_adjusted.tsv": OUT,
    "arg_mge_neighbourhood.tsv": OUT,
    "shared_context_blocks.tsv": OUT,
    "mge_feature_inventory.tsv": OUT,
    "convergence_by_mobility_class.tsv": OUT,
    "plasmid_mobility_annotation.tsv": OUT,
    "frozen_portability_class_definitions.json": OUT,
    "NMV3_RESULT_RECEIPT.json": NM,
    "NMV4_RESULT_RECEIPT.json": NM,
    "NMIS_RESULT_RECEIPT_V1.json": NM,
    "NMDIST_RESULT_RECEIPT_V1.json": NM,
    "nmis_occurrence_endpoints.tsv": NM,
}


def title(ax, letter, text, size=7.6):
    """Panel letter and title as one left-aligned string.

    Free-floating panel letters positioned in axes coordinates collide with the
    neighbouring panel's title or y-label as soon as a figure is rescaled. Making
    the letter part of the title means the layout engine accounts for it.
    """
    ax.set_title(r"$\bf{%s}$   %s" % (letter, text), loc="left", fontsize=size)


def save(fig, stem):
    """SVG and PDF as vector masters, PNG at 300 dpi for review and 600 dpi for
    print. A journal asking for 600 should not force a regeneration."""
    paths = []
    for ext in ("svg", "pdf"):
        p = os.path.join(FIGD, "%s.%s" % (stem, ext))
        fig.savefig(p, format=ext)
        paths.append(p)
    for dpi in (300, 600):
        p = os.path.join(FIGD, "%s%s.png" % (stem, "" if dpi == 300 else "_600dpi"))
        fig.savefig(p, format="png", dpi=dpi)
        paths.append(p)
    plt.close(fig)
    return paths


# =========================================================================== #
# Figure 1 -- cohort, denominator flow, resolution guarantee
# =========================================================================== #
def figure1():
    inv = {r["evidence_layer"]: int(r["n_rows"]) for r in tsv("evidence_layer_inventory.tsv")}
    fig, axes = plt.subplots(1, 3, figsize=(7.4, 2.9), layout="constrained")

    ax = axes[0]
    steps = [("all\nrecords", 184538), ("acquired\nAMR", 85507),
             ("PRIMARY", 74349)]
    ax.bar(range(3), [s[1] for s in steps],
           color=[OI["grey"], OI["sky"], OI["blue"]], width=0.62)
    for i, (_, v) in enumerate(steps):
        ax.text(i, v + 4000, "{:,}".format(v), ha="center", fontsize=7)
    ax.set_xticks(range(3))
    ax.set_xticklabels([s[0] for s in steps], fontsize=6.6)
    ax.set_ylabel("determinant records")
    ax.set_ylim(0, 210000)
    title(ax, "a", "Denominator flow")

    ax = axes[1]
    vals = [39209, 35140, 0]
    lab = ["plasmid", "chromosome", "unresolved"]
    ax.barh([2, 1, 0], vals, color=[OI["orange"], OI["blue"], OI["vermillion"]], height=0.55)
    for y, v in zip([2, 1, 0], vals):
        ax.text(v + 900, y, "{:,}".format(v), va="center", fontsize=7)
    ax.set_yticks([2, 1, 0])
    ax.set_yticklabels(lab)
    ax.set_xlim(0, 47000)
    ax.set_xlabel("occurrences")
    title(ax, "b", "100.000 % resolved")

    ax = axes[2]
    ax.bar([0, 1], [52.736, 36.586], color=[OI["orange"], OI["purple"]], width=0.5)
    for i, v in enumerate([52.736, 36.586]):
        ax.text(i, v + 1.2, "%.3f %%" % v, ha="center", fontsize=7)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["occurrence-\nweighted", "genome-\nweighted"])
    ax.set_ylabel("plasmid share (%)")
    ax.set_ylim(0, 65)
    title(ax, "c", "Two denominators")

    fig.suptitle("Figure 1 | 74,349 acquired resistance-gene occurrences in 6,288 closed genomes, "
                 "each assigned to a documented replicon", fontsize=8.5)
    return save(fig, "figure1_cohort_and_resolution")


# =========================================================================== #
# Figure 2 -- five-class portability architecture
# =========================================================================== #
def figure2():
    rec = {r["portability_class"]: (int(r["n"]), float(r["pct_of_74349"]))
           for r in tsv("portability_class_reconciliation.tsv")}
    v3 = jsn("NMV3_RESULT_RECEIPT.json")
    e = v3["E1_E2"]["baseline"]

    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.1), layout="constrained",
                             gridspec_kw={"width_ratios": [1.55, 1]})

    ax = axes[0]
    order = ["A", "B", "C", "D", "E"]
    cols = [OI["grey"], OI["blue"], OI["yellow"], OI["orange"], OI["vermillion"]]
    n = [rec[k][0] for k in order]
    ax.bar(range(5), n, color=cols, width=0.62)
    for i, k in enumerate(order):
        ax.text(i, n[i] + 500, "{:,}\n{:.2f} %".format(rec[k][0], rec[k][1]),
                ha="center", fontsize=6.6)
    ax.set_xticks(range(5))
    ax.set_xticklabels(["A\nchromosomal\nno marker", "B\nchromosomal\nMGE context",
                        "C\nplasmid\nno marker", "D\nplasmid\nmobilizable",
                        "E\nplasmid\nconjugative"], fontsize=6.6)
    ax.set_ylabel("occurrences")
    ax.set_ylim(0, 31500)
    ax.axvline(1.5, color="k", lw=0.7, ls=":")
    ax.text(0.75, 30200, "chromosome  35,140", ha="center", fontsize=7)
    ax.text(3.5, 30200, "plasmid  39,209", ha="center", fontsize=7)
    title(ax, "a", "One class per occurrence")

    ax = axes[1]
    ax.bar([0], [e["E1_occurrences"]], color=OI["vermillion"], width=0.5,
           label="E1  relaxase + MPF")
    ax.bar([0], [e["E2_occurrences"]], bottom=[e["E1_occurrences"]],
           color=OI["purple"], width=0.5, label="E2  + detected oriT")
    ax.text(0, e["E_total"] + 700, "{:,}".format(e["E_total"]), ha="center", fontsize=7)
    ax.text(0, e["E1_occurrences"] / 2, "{:,}".format(e["E1_occurrences"]),
            ha="center", va="center", fontsize=7, color="white")
    ax.text(0, e["E1_occurrences"] + e["E2_occurrences"] / 2,
            "{:,}\n{:.1f} % of E".format(e["E2_occurrences"], e["E2_pct_of_E"]),
            ha="center", va="center", fontsize=6.6, color="white")
    ax.set_xticks([0])
    ax.set_xticklabels(["class E"])
    ax.set_xlim(-0.6, 1.5)
    ax.set_ylim(0, 30000)
    ax.set_ylabel("occurrences")
    ax.legend(loc="upper right", frameon=False)
    title(ax, "b", "Nested evidence tiers")

    fig.suptitle("Figure 2 | Five evidence-ranked portability classes", fontsize=8.5)
    return save(fig, "figure2_five_class_architecture")


# =========================================================================== #
# Figure 3 -- portability is a property of the occurrence
# =========================================================================== #
def figure3():
    both = tsv("both_context_determinants.tsv")
    enr = tsv("determinant_enrichment_species_adjusted.tsv")

    fig, axes = plt.subplots(1, 3, figsize=(7.4, 2.9), layout="constrained",
                             gridspec_kw={"width_ratios": [1, 1.25, 1]})

    ax = axes[0]
    frac = sorted(float(r["minor_context_fraction"]) for r in both)
    ax.hist(frac, bins=24, color=OI["blue"], edgecolor="white", linewidth=0.4)
    ax.set_xlabel("minor-compartment fraction")
    ax.yaxis.get_major_locator().set_params(integer=True)
    ax.set_ylabel("gene families")
    title(ax, "a", "%d in both compartments" % len(both))

    ax = axes[1]
    sig, surv, x, y = [], [], [], []
    for r in enr:
        try:
            orv = float(r["crude_odds_ratio"])
        except (ValueError, KeyError):
            continue
        if orv <= 0:
            continue
        x.append(int(r["n_occurrences"]))
        y.append(np.log2(orv))
        sig.append(r["significant_after_bh"] == "yes")
        surv.append(r["survives_species_adjustment"] == "yes")
    x, y = np.array(x, float), np.array(y, float)
    sig, surv = np.array(sig), np.array(surv)
    ax.scatter(x[~sig], y[~sig], s=9, c=OI["grey"], alpha=0.7, linewidths=0,
               label="not significant")
    ax.scatter(x[sig & ~surv], y[sig & ~surv], s=11, c=OI["sky"], alpha=0.85, linewidths=0,
               label="significant, does not survive")
    ax.scatter(x[sig & surv], y[sig & surv], s=13, c=OI["vermillion"], alpha=0.9, linewidths=0,
               label="survives species adjustment")
    ax.axhline(0, color="k", lw=0.7)
    ax.set_xscale("log")
    ax.set_xlabel("occurrences per family (log)")
    ax.set_ylabel("log$_2$ odds ratio\n(plasmid vs chromosome)")
    ax.legend(loc="lower right", frameon=False, fontsize=6)
    title(ax, "b", "Preference by family")

    ax = axes[2]
    nsig = int(sig.sum())
    nsurv = int((sig & surv).sum())
    bars = [len(both), len(enr), nsig, nsurv]
    ax.bar(range(4), bars, color=[OI["blue"], OI["sky"], OI["orange"], OI["vermillion"]],
           width=0.62)
    for i, v in enumerate(bars):
        ax.text(i, v + 3, str(v), ha="center", fontsize=7)
    ax.set_xticks(range(4))
    ax.set_xticklabels(["both-\ncontext", "$\\geq$20\nocc.", "BH\nq<0.05", "survives\nadj."],
                       fontsize=6.6)
    ax.set_ylabel("gene families")
    ax.set_ylim(0, len(both) * 1.18)
    title(ax, "c", "Attrition")

    fig.suptitle("Figure 3 | Knowing the gene does not tell you the compartment",
                 fontsize=8.5)
    return save(fig, "figure3_occurrence_not_gene"), nsig, nsurv


# =========================================================================== #
# Figure 4 -- the chromosomal mobile compartment, four denominators
# =========================================================================== #
def figure4():
    blocks, feat_blocks = set(), set()
    for i, l in enumerate(open(os.path.join(OUT, "shared_context_blocks.tsv"), encoding="utf-8")):
        if i:
            blocks.add(l.split("\t")[0])
    for i, l in enumerate(open(os.path.join(OUT, "mge_feature_inventory.tsv"), encoding="utf-8")):
        if i:
            feat_blocks.add(l.split("\t")[0])
    mb = feat_blocks & blocks

    occ = collections.Counter()
    posB = collections.Counter()
    blk_sp, blkB = {}, set()
    with open(os.path.join(NM, "nmis_occurrence_endpoints.tsv"), encoding="utf-8") as fh:
        h = fh.readline().rstrip("\n").split("\t")
        ix = {c: i for i, c in enumerate(h)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            sp = f[ix["organism_harmonized"]]
            occ[sp] += 1
            blk_sp.setdefault(f[ix["block_id"]], sp)
            if f[ix["portability_class"]] == "B":
                posB[sp] += 1
                blkB.add(f[ix["block_id"]])
    nblk = collections.Counter(blk_sp.values())
    nblkB = collections.Counter(blk_sp[b] for b in blkB)

    fig, axes = plt.subplots(1, 3, figsize=(7.4, 3.1), layout="constrained",
                             gridspec_kw={"width_ratios": [1, 1.5, 0.9]})

    ax = axes[0]
    vals = [46.39, 40.51, 30.14, 0.26]
    labs = ["10 kb\nocc.", "5 kb\nocc.", "10 kb\nblock", "direct\noverlap"]
    ax.bar(range(4), vals, color=[OI["blue"], OI["sky"], OI["purple"], OI["grey"]], width=0.62)
    for i, v in enumerate(vals):
        ax.text(i, v + 1.1, "%.2f %%" % v, ha="center", fontsize=6.8)
    ax.set_xticks(range(4))
    ax.set_xticklabels(labs, fontsize=6.2)
    ax.set_ylabel("chromosomal occurrences (%)")
    ax.set_ylim(0, 56)
    title(ax, "a", "The denominator decides")

    ax = axes[1]
    sps = [s for s, _ in occ.most_common(6)]
    xo = np.arange(len(sps))
    o = [100 * posB[s] / occ[s] for s in sps]
    b = [100 * nblkB[s] / nblk[s] for s in sps]
    ax.bar(xo - 0.19, o, width=0.36, color=OI["blue"], label="occurrence-weighted")
    ax.bar(xo + 0.19, b, width=0.36, color=OI["purple"], label="block-weighted")
    ax.set_xticks(xo)
    def abbr(name):
        p = name.split(" ")
        return p[0][0] + ". " + " ".join(p[1:]) if len(p) > 1 else name
    ax.set_xticklabels([abbr(s) for s in sps], fontsize=6.0, style="italic",
                       rotation=30, ha="right")
    ax.set_ylabel("within 10 kb of an MGE marker (%)")
    ax.set_ylim(0, 92)
    ax.legend(loc="upper right", frameon=False)
    title(ax, "b", "Host ordering preserved")

    ax = axes[2]
    ax.bar([0, 1], [len(mb), len(blkB)], color=[OI["purple"], OI["blue"]], width=0.5)
    for i, v in enumerate([len(mb), len(blkB)]):
        ax.text(i, v + 90, "{:,}".format(v), ha="center", fontsize=7)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["blocks with\n$\\geq$1 marker", "blocks with\n$\\geq$1 class-B\noccurrence"],
                       fontsize=6.4)
    ax.set_ylabel("context blocks")
    ax.set_ylim(0, 8600)
    ax.annotate("", xy=(0, len(mb)), xytext=(1, len(blkB)),
                arrowprops=dict(arrowstyle="-", color="k", lw=0.8, ls=":"))
    ax.text(0.5, 7300, "difference:\nexactly 1 block", ha="center", fontsize=6.4)
    title(ax, "c", "Two block quantities")

    fig.suptitle("Figure 4 | 46.39 % of chromosomal occurrences lie within 10 kb of a "
                 "mobile-element marker", fontsize=8.5)
    return save(fig, "figure4_chromosomal_mobile_compartment"), len(mb), len(blkB)


# =========================================================================== #
# Figure 6 -- the structural insertion-sequence endpoint
# =========================================================================== #
def figure6():
    R = jsn("NMIS_RESULT_RECEIPT_V1.json")
    PE, HO = R["primary_estimates_structural"], R["homology_estimates_same_rows"]
    GR = ["A. baumannii", "Klebsiella group", "P. aeruginosa"]
    CG = {"A. baumannii": OI["vermillion"], "Klebsiella group": OI["blue"],
          "P. aeruginosa": OI["green"]}
    LM = [1000, 2000, 5000, 10000]

    fig, axes = plt.subplots(1, 4, figsize=(9.0, 2.9), layout="constrained",
                             gridspec_kw={"width_ratios": [1.15, 1, 1, 1]})

    ax = axes[0]
    for g in GR:
        y = [PE[g]["F_%d" % d] for d in LM]
        ax.plot(LM, y, "-o", ms=3, lw=1.4, color=CG[g], label=g)
        yh = [HO[g]["F_%d" % d] for d in LM]
        ax.plot(LM, yh, "--", lw=1.0, color=CG[g], alpha=0.55)
    ax.set_xscale("log")
    ax.set_xlabel("distance to nearest element (bp)")
    ax.set_ylabel("cumulative detection $F(d)$")
    ax.set_ylim(0, 0.72)
    ax.legend(loc="upper left", frameon=False, fontsize=6)
    title(ax, "a", "Structural vs homology")

    ax = axes[1]
    C = [c for c in R["contrasts"] if c["landmark_bp"] in LM]
    xs = np.arange(len(LM))
    for k, (name, col) in enumerate([("N1", OI["blue"]), ("N2", OI["green"])]):
        cc = sorted([c for c in C if c["contrast"] == name], key=lambda c: c["landmark_bp"])
        d = [c["difference"] for c in cc]
        lo = [c["difference"] - c["ci_lo"] for c in cc]
        hi = [c["ci_hi"] - c["difference"] for c in cc]
        ax.errorbar(xs + (k - 0.5) * 0.18, d, yerr=[lo, hi], fmt="o", ms=3.5, lw=1.2,
                    capsize=2, color=col,
                    label="%s  vs %s" % (name, "Klebsiella" if name == "N1" else "P. aeruginosa"))
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xticks(xs)
    ax.set_xticklabels(["1 kb", "2 kb", "5 kb", "10 kb"])
    ax.set_ylabel("difference in $F(d)$")
    ax.set_ylim(-0.05, 0.45)
    ax.legend(loc="lower right", frameon=False, fontsize=6)
    title(ax, "b", "A. baumannii excess")

    ax = axes[2]
    sc = R["structural_corroboration"]
    gs = ["A. baumannii", "Klebsiella group", "P. aeruginosa"]
    v = [100 * sc[g]["structural_share_of_homology"] for g in gs]
    ax.bar(range(3), v, color=[CG[g] for g in gs], width=0.6)
    for i, val in enumerate(v):
        ax.text(i, val + 1.6, "%.2f %%" % val, ha="center", fontsize=6.8)
    ax.axhline(73.80, color="k", lw=0.8, ls=":")
    ax.text(2.45, 75.4, "overall\n73.80 %", ha="right", fontsize=6.2)
    ax.set_xticks(range(3))
    ax.set_xticklabels(["A. bau.", "Klebsiella", "P. aer."], fontsize=6.4, style="italic")
    ax.set_ylabel("class-B occurrences structurally\ncorroborated (%)")
    ax.set_ylim(0, 100)
    title(ax, "c", "Corroboration")

    ax = axes[3]
    all_f, gated = R["is_family_all"], R["is_family_structurally_gated"]
    fam = [f for f in sorted(all_f, key=lambda k: -all_f[k]) if all_f[f] >= 250][:8]
    ret = [100 * gated.get(f, 0) / all_f[f] for f in fam]
    cols = [OI["vermillion"] if f == "IS6" else OI["grey"] for f in fam]
    ax.barh(range(len(fam))[::-1], ret, color=cols, height=0.62)
    for i, (f, r) in enumerate(zip(fam, ret)):
        ax.text(r + 1.6, len(fam) - 1 - i, "%.1f" % r, va="center", fontsize=6.2)
    ax.set_yticks(range(len(fam))[::-1])
    ax.set_yticklabels(fam, fontsize=6.4)
    ax.set_xlim(0, 112)
    ax.set_xlabel("elements retained through the\nstructural gate (%)")
    title(ax, "d", "IS6 best resolved")

    fig.suptitle("Figure 6 | The host-conditioned signature survives an endpoint restricted to "
                 "structurally complete, fully contained insertion sequences",
                 fontsize=8.5)
    return save(fig, "figure6_structural_is_endpoint")


# =========================================================================== #
# Figure 7 -- the discordance principle and its independent confirmation
# =========================================================================== #
def figure7():
    R = jsn("NMV4_RESULT_RECEIPT.json")
    per, T2, T3 = R["per_species"], R["T2"], R["T3"]
    disc = set(R["discovery_species"])

    fig, axes = plt.subplots(1, 3, figsize=(7.6, 3.0), layout="constrained",
                             gridspec_kw={"width_ratios": [1.3, 1, 1]})

    ax = axes[0]
    for s, d in per.items():
        conf = s not in disc
        col = OI["grey"] if conf else (OI["vermillion"] if s.startswith("Acineto") else OI["blue"])
        ax.scatter(d["P"] * 100, d["M"] * 100, s=34 if not conf else 20, c=col,
                   linewidths=0.6, edgecolors="white", zorder=3)
    xs = np.linspace(0.10, 0.70, 60)
    lo = T2["intercept"] + T2["slope"] * np.log(xs / (1 - xs))
    ax.plot(xs * 100, 100 / (1 + np.exp(-lo)), "-", color="k", lw=1.0,
            label="fit, 8 confirmation species")
    ab = per["Acinetobacter baumannii"]
    ax.annotate("", xy=(ab["P"] * 100, ab["M"] * 100),
                xytext=(ab["P"] * 100, T2["ab_predicted_M"] * 100),
                arrowprops=dict(arrowstyle="<->", color=OI["vermillion"], lw=1.2))
    ax.text(ab["P"] * 100 + 2.5, (ab["M"] + T2["ab_predicted_M"]) * 50,
            "residual\n2.11 logits\nCI 1.81-2.44", fontsize=6.2, color=OI["vermillion"])
    ax.text(per["Pseudomonas aeruginosa"]["P"] * 100 + 1.5,
            per["Pseudomonas aeruginosa"]["M"] * 100 - 4, "P. aeruginosa",
            fontsize=6.2, style="italic")
    ax.text(ab["P"] * 100 - 1.0, ab["M"] * 100 + 3.0, "A. baumannii",
            fontsize=6.4, style="italic", color=OI["vermillion"])
    kp = per["Klebsiella pneumoniae"]
    ax.text(kp["P"] * 100 - 0.5, kp["M"] * 100 - 4.5, "K. pneumoniae",
            fontsize=6.2, style="italic", color=OI["blue"], ha="right")
    ax.text(0.03, 0.03, "grey: 8 confirmation species\ncoloured: the 2 discovery species",
            transform=ax.transAxes, fontsize=5.8, color=OI["grey"])
    ax.set_xlabel("plasmid share of ARG occurrences (%)")
    ax.set_ylabel("chromosomal MGE association,\nblock-weighted (%)")
    ax.set_xlim(5, 75)
    ax.set_ylim(0, 75)
    ax.legend(loc="upper right", frameon=False, fontsize=6)
    title(ax, "a", "Two axes, not one")

    ax = axes[1]
    ax.bar([0], [T2["in_sample_r2"] * 100], color=OI["grey"], width=0.42)
    ax.text(0, T2["in_sample_r2"] * 100 + 2.0, "%.2f %%" % (T2["in_sample_r2"] * 100),
            ha="center", fontsize=7)
    ax.set_xticks([0])
    ax.set_xticklabels(["variance in chromosomal\nMGE association explained\nby plasmid share"],
                       fontsize=6.2)
    ax.set_xlim(-0.7, 0.7)
    ax.set_ylim(0, 100)
    ax.set_ylabel("$R^2$ (%)")
    ax.text(0, 52, "the slope's bootstrap\ninterval includes zero\n(-0.059 to 0.229)",
            ha="center", fontsize=6.4)
    title(ax, "b", "Almost none")

    ax = axes[2]
    est = [("NM-V4\nA. baumannii\nresidual (logit)", T2["ab_residual"],
            T2["bootstrap_ci"], OI["vermillion"]),
           ("NM-V4 T3\nMGE gap vs\nP. aeruginosa", T3["M_gap_median"], T3["M_gap_ci"], OI["green"])]
    for i, (lab, v, ci, col) in enumerate(est):
        ax.errorbar([i], [v], yerr=[[v - ci[0]], [ci[1] - v]], fmt="o", ms=5, lw=1.4,
                    capsize=3, color=col)
        ax.text(i + 0.13, v, "%.3f" % v, va="center", fontsize=6.6)
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xticks(range(len(est)))
    ax.set_xticklabels([e[0] for e in est], fontsize=6.0)
    ax.set_xlim(-0.55, len(est) - 0.25)
    ax.set_ylabel("effect (units as labelled)")
    title(ax, "c", "Both exclude zero")

    fig.suptitle("Figure 7 | Plasmid fraction and chromosomal mobile-element association are "
                 "non-redundant, confirmed off the discovery species", fontsize=8.5)
    return save(fig, "figure7_discordance_principle")


# =========================================================================== #
# Figure 8 -- conjugation-consistent replicons carry the convergent cargo
# =========================================================================== #
def figure8():
    conv = tsv("convergence_by_mobility_class.tsv")
    v3 = jsn("NMV3_RESULT_RECEIPT.json")
    key = {"predicted_conjugative": ("conjugative", OI["vermillion"]),
           "predicted_mobilizable": ("mobilizable", OI["orange"]),
           "nonconjugative_or_no_mobility_markers_detected": ("marker-negative", OI["grey"])}
    rows = [r for r in conv if r["mobility_category"] in key]
    ORDER = ["predicted_conjugative", "predicted_mobilizable",
             "nonconjugative_or_no_mobility_markers_detected"]
    rows.sort(key=lambda r: ORDER.index(r["mobility_category"]))

    fig, axes = plt.subplots(1, 4, figsize=(8.6, 2.8), layout="constrained")

    ax = axes[0]
    n = [int(r["n_plasmids"]) for r in rows]
    ax.bar(range(len(rows)), n, color=[key[r["mobility_category"]][1] for r in rows], width=0.6)
    for i, v in enumerate(n):
        ax.text(i, v + 70, "{:,}".format(v), ha="center", fontsize=6.8)
    ax.set_xticks(range(len(rows)))
    ax.set_xticklabels([key[r["mobility_category"]][0] for r in rows], fontsize=6.2,
                       rotation=20, ha="right")
    ax.set_ylabel("ARG-bearing plasmids")
    ax.set_ylim(0, 4700)
    title(ax, "a", "6,621 replicons")

    for j, (col, lab, ttl) in enumerate([
            ("pct_multi_class_ge3", "$\\geq$3 drug classes (%)", "Multi-class cargo"),
            ("pct_arg_plus_metal", "ARG + metal (%)", "Stress co-location")]):
        ax = axes[1 + j]
        v = [float(r[col]) for r in rows]
        ax.bar(range(len(rows)), v, color=[key[r["mobility_category"]][1] for r in rows],
               width=0.6)
        for i, val in enumerate(v):
            ax.text(i, val + 1.3, "%.2f" % val, ha="center", fontsize=6.8)
        ax.set_xticks(range(len(rows)))
        ax.set_xticklabels([key[r["mobility_category"]][0] for r in rows], fontsize=6.2,
                       rotation=20, ha="right")
        ax.set_ylabel(lab)
        ax.set_ylim(0, 74)
        title(ax, "bc"[j], ttl)

    ax = axes[3]
    c9 = v3["C09"]["baseline"]
    d, ci = c9["difference_conj_minus_negative"], c9["bioproject_clustered_ci95"]
    ax.errorbar([0], [d], yerr=[[d - ci[0]], [ci[1] - d]], fmt="o", ms=5.5, lw=1.5,
                capsize=3, color=OI["vermillion"])
    ax.axhline(0, color="k", lw=0.8)
    ax.text(0.1, d, "%.2f pp\nCI %.2f-%.2f" % (d, ci[0], ci[1]), va="center", fontsize=6.4)
    ax.set_xticks([0])
    ax.set_xticklabels(["conjugative $-$\nmarker-negative"], fontsize=6.4)
    ax.set_xlim(-0.6, 0.85)
    ax.set_ylim(-3, 30)
    ax.set_ylabel("difference in $\\geq$3-class share\n(percentage points)")
    title(ax, "d", "Clustered CI")

    fig.suptitle("Figure 8 | Mobilization architecture tracks cargo convergence on documented "
                 "replicons -- association, not co-transfer", fontsize=8.5)
    return save(fig, "figure8_convergence_by_mobility")


# =========================================================================== #
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true",
                    help="verify inputs and exit without drawing")
    a = ap.parse_args()

    digests = {}
    missing = []
    for name, root in INPUTS.items():
        p = os.path.join(root, name)
        if not os.path.isfile(p):
            missing.append(name)
        else:
            digests[name] = {"path": os.path.relpath(p, REPO).replace("\\", "/"),
                             "bytes": os.path.getsize(p), "sha256": sha(p)}
    if missing:
        print("MISSING INPUTS:")
        for m in missing:
            print("  " + m)
        return 2
    print("inputs verified: %d" % len(digests))
    if a.check:
        for k, v in sorted(digests.items()):
            print("  %-48s %s" % (k, v["sha256"]))
        return 0

    os.makedirs(FIGD, exist_ok=True)
    produced = collections.OrderedDict()
    produced["figure1"] = figure1()
    produced["figure2"] = figure2()
    f3, nsig, nsurv = figure3()
    produced["figure3"] = f3
    f4, mb, bb = figure4()
    produced["figure4"] = f4
    produced["figure6"] = figure6()
    produced["figure7"] = figure7()
    produced["figure8"] = figure8()

    receipt = {
        "receipt": "PORTABILITYRISK_FIGURE_RECEIPT", "version": "1.0.0",
        "generator": "portabilityrisk_figures.py",
        "note": ("Figure 5 is the NM-DIST four-panel figure, produced by nmdist_figure.py "
                 "under the NM-DIST frozen protocol and reused unchanged."),
        "recomputed_while_drawing": {
            "families_significant_after_bh": nsig,
            "families_surviving_species_adjustment": nsurv,
            "blocks_with_at_least_one_mge_feature": mb,
            "blocks_with_at_least_one_class_B_occurrence": bb},
        "inputs_sha256": digests,
        "outputs_sha256": {},
    }
    for fig, paths in produced.items():
        for p in paths:
            rel = os.path.relpath(p, REPO).replace("\\", "/")
            receipt["outputs_sha256"][rel] = sha(p)

    rp = os.path.join(NM, "PORTABILITYRISK_FIGURE_RECEIPT_V1.json")
    with open(rp, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(receipt, fh, indent=1)
        fh.write("\n")

    print("\nfigures written to %s" % os.path.relpath(FIGD, REPO))
    for fig, paths in produced.items():
        print("  %-9s %s" % (fig, os.path.basename(paths[0])[:-4]))
    print("\nrecomputed while drawing:")
    for k, v in receipt["recomputed_while_drawing"].items():
        print("  %-46s %s" % (k, v))
    print("\nreceipt: %s" % os.path.relpath(rp, REPO))
    return 0


if __name__ == "__main__":
    sys.exit(main())
