"""NM-DIST Phase 4 -- one manuscript figure, four panels."""
import argparse, csv, json, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# Okabe-Ito, colour-blind safe
COL = {"A. baumannii": "#D55E00", "P. aeruginosa": "#0072B2", "Klebsiella group": "#009E73"}
LBL = {"A. baumannii": "$\\it{A.\\ baumannii}$",
       "P. aeruginosa": "$\\it{P.\\ aeruginosa}$",
       "Klebsiella group": "$\\it{Klebsiella}$ group"}
LAND = [1000, 2000, 5000, 10000]
G = ["A. baumannii", "P. aeruginosa", "Klebsiella group"]


def rd(p):
    return list(csv.DictReader(open(p, encoding="utf-8"), delimiter="\t"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    cur = json.load(open(os.path.join(a.dir, "nmdist_curves.json"), encoding="utf-8"))
    con = rd(os.path.join(a.dir, "nmdist_species_contrasts.tsv"))
    sens = rd(os.path.join(a.dir, "nmdist_sensitivity_results.tsv"))

    plt.rcParams.update({"font.size": 8, "axes.linewidth": 0.7, "xtick.major.width": 0.7,
                         "ytick.major.width": 0.7, "font.family": "DejaVu Sans"})
    fig, ax = plt.subplots(2, 2, figsize=(7.2, 5.6))

    # ---- A: cumulative detection curves with BioProject-bootstrap bands ----
    A = ax[0][0]
    for g in G:
        gd = cur["curves"][g]["grid"]; f = cur["curves"][g]["F"]
        bg = cur["bands"][g]["grid"]; lo = cur["bands"][g]["lo"]; hi = cur["bands"][g]["hi"]
        A.fill_between([x / 1000 for x in bg], lo, hi, color=COL[g], alpha=0.18, linewidth=0)
        A.plot([x / 1000 for x in gd], f, color=COL[g], lw=1.6, label=LBL[g])
    A.set_xlim(0, 10); A.set_ylim(0, 0.75)
    A.set_xlabel("distance from ARG to nearest MGE marker (kb)")
    A.set_ylabel("weighted cumulative detection $F(d)$")
    A.set_title("a   Cumulative MGE detection, block-balanced", loc="left", fontsize=8.5,
                fontweight="bold")
    A.legend(frameon=False, fontsize=7, loc="center right", bbox_to_anchor=(1.0, 0.45))
    A.grid(alpha=0.25, lw=0.5)

    # ---- B: contrasts at landmarks ----
    Bx = ax[0][1]
    off = {"P1": -0.16, "P2": 0.16}
    nm = {"P1": "$\\it{A.\\ baumannii}$ − $\\it{Klebsiella}$",
          "P2": "$\\it{A.\\ baumannii}$ − $\\it{P.\\ aeruginosa}$"}
    mk = {"P1": "o", "P2": "s"}
    for tag in ("P1", "P2"):
        xs, ys, lo, hi = [], [], [], []
        for i, t in enumerate(LAND):
            c = [c for c in con if c["contrast"] == tag and c["landmark_bp"] == str(t)][0]
            xs.append(i + off[tag]); ys.append(float(c["difference"]))
            lo.append(float(c["difference"]) - float(c["ci_lo"]))
            hi.append(float(c["ci_hi"]) - float(c["difference"]))
        Bx.errorbar(xs, ys, yerr=[lo, hi], fmt=mk[tag], ms=4, lw=1.2, capsize=2.5,
                    color="#000000" if tag == "P1" else "#666666", label=nm[tag])
    Bx.axhline(0, color="#999", lw=0.8, ls="--")
    Bx.set_xticks(range(4)); Bx.set_xticklabels(["1 kb", "2 kb", "5 kb", "10 kb"])
    Bx.set_ylim(0, 0.62)
    Bx.set_ylabel("difference in $F(d)$")
    Bx.set_title("b   Species contrasts at frozen landmarks", loc="left", fontsize=8.5,
                 fontweight="bold")
    Bx.legend(frameon=False, fontsize=6.5, loc="lower left")
    Bx.grid(alpha=0.25, lw=0.5, axis="y")
    Bx.text(0.98, 0.95, "Holm-corrected, all $p$ = 0.001", transform=Bx.transAxes,
            ha="right", va="top", fontsize=6, color="#444")

    # ---- C: marker-type decomposition ----
    Cx = ax[1][0]
    S = {(r["sensitivity"], r["group"]): r for r in sens}
    w = 0.36
    for i, g in enumerate(G):
        Cx.bar(i - w / 2, float(S[("S6", g)]["F_10000"]), w, color=COL[g], edgecolor="black",
               lw=0.6)
        Cx.bar(i + w / 2, float(S[("SEC_INT", g)]["F_10000"]), w, color=COL[g], edgecolor="black",
               lw=0.6, hatch="////", alpha=0.55)
    Cx.set_xticks(range(3)); Cx.set_xticklabels([LBL[g] for g in G], fontsize=7)
    Cx.set_ylabel("$F$(10 kb)")
    Cx.set_ylim(0, 0.72)
    Cx.set_title("c   Marker-type decomposition", loc="left", fontsize=8.5, fontweight="bold")
    Cx.legend(handles=[Line2D([], [], color="#555", lw=6, label="IS / transposase"),
                       Line2D([], [], color="#555", lw=6, alpha=0.55, label="integrase / integron")],
              frameon=False, fontsize=6.5, loc="upper right")
    Cx.grid(alpha=0.25, lw=0.5, axis="y")

    # ---- D: occurrence-weighted vs block-balanced ----
    Dx = ax[1][1]
    for i, g in enumerate(G):
        Dx.bar(i - w / 2, float(S[("S1", g)]["F_10000"]), w, color=COL[g], edgecolor="black",
               lw=0.6)
        Dx.bar(i + w / 2, float(S[("PRIMARY", g)]["F_10000"]), w, color=COL[g], edgecolor="black",
               lw=0.6, hatch="....", alpha=0.6)
    Dx.set_xticks(range(3)); Dx.set_xticklabels([LBL[g] for g in G], fontsize=7)
    Dx.set_ylabel("$F$(10 kb)")
    Dx.set_ylim(0, 0.9)
    Dx.set_title("d   Weighting sensitivity", loc="left", fontsize=8.5, fontweight="bold")
    Dx.legend(handles=[Line2D([], [], color="#555", lw=6, label="occurrence-weighted (S1)"),
                       Line2D([], [], color="#555", lw=6, alpha=0.6,
                              label="block-balanced (primary)")],
              frameon=False, fontsize=6.5, loc="upper right")
    Dx.grid(alpha=0.25, lw=0.5, axis="y")

    fig.tight_layout(pad=1.1)
    base = os.path.join(a.outdir, "nmdist_figure")
    for ext, kw in (("pdf", {}), ("svg", {}), ("png", {"dpi": 300})):
        fig.savefig(base + "." + ext, **kw)
        print("  wrote %s.%s" % (os.path.basename(base), ext))
    plt.close(fig)


if __name__ == "__main__":
    main()
