"""Redraw Extended Data Fig. 5 deterministically from the receipts.

The submitted version of this figure was drawn by a generative model from a
written specification. Every value in it was checked, but a reviewer cannot see
that check, and a figure whose provenance is "a model drew it" is a liability
under a journal's AI policy. This replaces it with a figure that a script draws
from the receipts, so the same command reproduces the same file.

No value is typed into this file. Each is read from the artefact that owns it,
and the run aborts if any of them is missing.
"""
import os as _os


def _dir(var):
    """Resolve an input directory from the environment.

    The private copy of this script carried an absolute path on the
    author's machine. The public copy does not substitute a plausible
    path, because a path that looks right and is wrong is worse than one
    that is obviously missing. Set the variable, or the run stops here.
    """
    v = _os.environ.get(var)
    if not v:
        raise SystemExit(
            "set %s to the directory it names; see README" % var)
    return v.rstrip("/\\") + "/"


import hashlib
import io
import json
import os
import sys
import textwrap

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

NM = _dir("PORTABILITYRISK_REPO_DIR") + "docs/nature_microbiology/"
OUT = NM + "figures/submission_v4/"
STEM = "Extended_Data_Fig_5_registration_and_verification"

NAVY, GREY, LGREY = "#1F3064", "#5A5A5A", "#8A8A8A"
BLUE, GREEN, PURPLE, RED, AMBER = "#DCE7F5", "#DCEFE4", "#E8E0F0", "#C8322B", "#FBE7C6"
EDGE = {"blue": "#9DB9DC", "green": "#9CCFB4", "purple": "#BFAAD6", "amber": "#E0B36B"}

plt.rcParams.update({"font.size": 7.5, "font.family": "DejaVu Sans",
                     "pdf.fonttype": 42, "svg.fonttype": "none"})


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def jload(n):
    p = NM + n
    if not os.path.isfile(p):
        raise SystemExit("missing receipt: %s" % n)
    return json.load(io.open(p, encoding="utf-8"))


FIG_W_IN = 7.2


def wrap(text, w_axes, fontsize, pad=0.024):
    """Break `text` so it fits inside a card of width `w_axes`.

    Text that runs past the card edge is the failure this figure had on its first
    render, so width is computed rather than eyeballed: at DejaVu Sans a character
    is about 0.55 em wide, and the card's usable width in points is known.
    """
    usable_pt = (w_axes - pad) * FIG_W_IN * 72.0
    n = max(8, int(usable_pt / (fontsize * 0.55)))
    return textwrap.fill(text, width=n)


FIG_H_IN = 5.6


def _lines(txt, w, fs):
    return wrap(txt, w, fs).count("\n") + 1


def card(ax, x, y, w, h, fc, ec, title, detail, sub, bold_sub=False):
    """Draw a card, flowing the three text levels top-down.

    The first version anchored the title to the top and the sub-detail to the
    bottom, so a detail line that wrapped to two lines collided with the
    sub-detail. Here each block is placed below the one above it, using the line
    count that the wrap actually produced.
    """
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.004,rounding_size=0.008",
                                linewidth=0.7, facecolor=fc, edgecolor=ec, zorder=3))
    ts, ds, ss = 7.3, 5.9, 5.6
    while _lines(title, w, ts) > 1 and ts > 5.9:
        ts -= 0.2
    pt_to_ax = 1.0 / (FIG_H_IN * 72.0)
    cur = y + h - 0.016
    for txt, fs, col, wt in ((title, ts, NAVY, "bold"),
                             (detail, ds, GREY, "normal"),
                             (sub, ss, LGREY, "bold" if bold_sub else "normal")):
        body = wrap(txt, w, fs)
        ax.text(x + 0.012, cur, body, fontsize=fs, color=col, weight=wt,
                zorder=4, va="top", linespacing=1.22)
        cur -= (_lines(txt, w, fs) * fs * 1.30 + fs * 0.55) * pt_to_ax


def padlock(ax, cx, cy, s=0.011):
    """A padlock drawn from shapes. The emoji is not in the embedded font."""
    ax.add_patch(plt.Rectangle((cx - s, cy - s * 0.95), 2 * s, s * 1.5,
                               facecolor="white", edgecolor="none", zorder=5))
    ax.add_patch(matplotlib.patches.Arc((cx, cy + s * 0.55), s * 1.25, s * 1.35,
                                        theta1=0, theta2=180, lw=1.5, color="white",
                                        zorder=5))


def arrow(ax, x1, y1, x2, y2, color=LGREY):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>",
                                 mutation_scale=7, linewidth=0.8, color=color, zorder=2))


def stage(ax, n, y, title):
    ax.add_patch(plt.Circle((0.026, y), 0.014, color=NAVY, zorder=4))
    ax.text(0.026, y, str(n), fontsize=6.8, color="white", ha="center", va="center",
            weight="bold", zorder=5)
    ax.text(0.052, y, title, fontsize=8.6, weight="bold", color=NAVY, va="center")
    ax.plot([0.026, 0.985], [y - 0.024, y - 0.024], color="#CCCCCC", lw=0.7, zorder=1)


def main():
    prim = jload("NM_C1_4_NULL_RESULTS_V1.json")
    gate = jload("NM_C1_6_GATE_RECEIPT_V1.json")
    ver = jload("NM_C1_INDEPENDENT_VERIFICATION_V1.json")
    cen = jload("NMBG_C1_5_VERIFICATION_REPORT_V1.json")
    am = jload("NM_BACKGROUND_ENRICHMENT_AMENDMENT_001.json")

    man = NM + "NMBG_CHROMOSOME_ACCESSION_MANIFEST_FROZEN_V1.tsv"
    rows = sum(1 for _ in io.open(man, encoding="utf-8")) - 1
    sc = cen["spot_checks"]
    n_spot = sc["random_set"]["n"] + sc["zero_yield_set"]["n"]
    n_ok = sc["random_set"]["agree"] + sc["zero_yield_set"]["agree"]
    brute = [c for c in ver["checks"] if "brute-force" in c["check"]][0]["detail"]

    V = {"elements": "{:,}".format(prim["elements_total"]),
         "complete": "{:,}".format(prim["elements_complete_structural"]),
         "chrom": "{:,}".format(rows),
         "B": "{:,}".format(prim["n_permutations"]),
         "seed": str(prim["seed"]),
         "spot": "%d of %d" % (n_ok, n_spot),
         "checks": "%d of %d" % (ver["checks_passed"], ver["checks_total"]),
         "brute": brute.split()[0],
         "verdicts": str(len(gate["verdict_is_one_of"])),
         "sens": str(len(am["required_sensitivities"])),
         "gates": str(len(gate["core_gates"])),
         "verdict": gate["verdict"]}

    fig = plt.figure(figsize=(7.2, 5.6))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    W, H = 0.295, 0.130
    xs = [0.026, 0.352, 0.678]

    # ---- stage 1 --------------------------------------------------------
    stage(ax, 1, 0.968, "Design frozen before any outcome")
    y = 0.806
    card(ax, xs[0], y, W, H, BLUE, EDGE["blue"], "Protocol registered and hashed",
         "null model · permutations · seed · weightings", "analysis aborts on a hash mismatch")
    card(ax, xs[1], y, W, H, BLUE, EDGE["blue"], "Thresholds set in advance",
         "%s verdicts named before any result" % V["verdicts"],
         "no threshold may be added or moved later")
    card(ax, xs[2], y, W, H, BLUE, EDGE["blue"], "Sensitivity panel declared",
         "%s analyses fixed in the protocol" % V["sens"],
         "including the IS family named ahead of time")
    for i in (0, 1):
        arrow(ax, xs[i] + W + 0.004, y + H / 2, xs[i + 1] - 0.004, y + H / 2)

    # ---- stage 2 --------------------------------------------------------
    stage(ax, 2, 0.762, "Genome-wide census, verified independently")
    y = 0.600
    card(ax, xs[0], y, W, H, GREEN, EDGE["green"], "ISEScan 1.7.3, pinned environment",
         "HMMER 3.3.2 · BLAST+ 2.17.0 · FragGeneScan 1.32", "one chromosome at a time")
    card(ax, xs[1], y, W, H, PURPLE, EDGE["purple"], "%s structural elements" % V["elements"],
         "%s complete · %s chromosomes" % (V["complete"], V["chrom"]),
         "one status per chromosome, none dropped")
    card(ax, xs[2], y, W, H, GREEN, EDGE["green"], "%s independent re-runs" % n_spot,
         "seeded selection written before the first re-run",
         "%s reproduced the recorded count" % V["spot"], bold_sub=True)
    for i in (0, 1):
        arrow(ax, xs[i] + W + 0.004, y + H / 2, xs[i + 1] - 0.004, y + H / 2)

    # ---- stage 3: the sealed lane --------------------------------------
    stage(ax, 3, 0.556, "Outcome sealed until the cohort was fixed")
    y = 0.386
    w3, h3 = 0.225, 0.124
    x3 = [0.026, 0.268, 0.510]
    card(ax, x3[0], y, w3, h3, BLUE, EDGE["blue"], "Cohort built",
         "admission rules applied to every occurrence", "no outcome column read")
    card(ax, x3[1], y, w3, h3, BLUE, EDGE["blue"], "Cohort receipt on disk",
         "counts · species · families · exclusions", "written before any estimate exists")
    card(ax, x3[2], y, w3, h3, BLUE, EDGE["blue"], "Null model run",
         "relocation within each chromosome only",
         "B = %s · seed %s" % (V["B"], V["seed"]))
    for i in (0, 1):
        arrow(ax, x3[i] + w3 + 0.004, y + h3 / 2, x3[i + 1] - 0.004, y + h3 / 2, "#7FA3D0")
    ax.text(0.026, y + h3 + 0.012, "COHORT LANE", fontsize=6.0, color="#3A6BA5",
            weight="bold")

    ax.add_patch(FancyBboxPatch((0.752, y), 0.233, h3,
                                boxstyle="round,pad=0.004,rounding_size=0.008",
                                linewidth=0, facecolor=RED, zorder=3))
    padlock(ax, 0.868, y + h3 / 2 + 0.020)
    ax.text(0.868, y + h3 / 2 - 0.022, "OUTCOME SEALED", fontsize=8.0, color="white",
            ha="center", va="center", weight="bold", zorder=4)

    ax.plot([0.026, 0.985], [y - 0.030, y - 0.030], color=RED, lw=1.0, ls=(0, (4, 3)),
            zorder=2)
    ax.text(0.505, y - 0.045, "outcome inaccessible while the cohort was being fixed",
            fontsize=5.9, color=GREY, style="italic", ha="center", va="top")

    y2 = 0.176
    ax.text(0.026, y2 + h3 + 0.012, "OUTCOME LANE", fontsize=6.0, color="#1F6B3A",
            weight="bold")
    card(ax, x3[0], y2, w3, h3, GREEN, EDGE["green"], "Observed distances held unread",
         "no species outcome computed", "no threshold adjustable")
    card(ax, 0.510, y2, 0.235, h3, AMBER, EDGE["amber"], "Outcome read",
         "only after the cohort receipt was written", "gates then applied mechanically")
    arrow(ax, 0.868, y - 0.034, 0.700, y2 + h3 + 0.004)

    # ---- stage 4 --------------------------------------------------------
    stage(ax, 4, 0.152, "Gates and independent verification")
    y3 = 0.006
    hh = 0.108
    card(ax, xs[0], y3, W, hh, BLUE, EDGE["blue"], "%s gates applied mechanically" % V["gates"],
         "enrichment · comparators · balancing · clusters",
         "thresholds unchanged from registration")
    card(ax, xs[1], y3, W, hh, GREEN, EDGE["green"], "Separate code path",
         "every mechanic re-implemented brute force",
         "%s checks pass · %s distances" % (V["checks"], V["brute"]), bold_sub=True)
    card(ax, xs[2], y3, W, hh, PURPLE, EDGE["purple"], "Verdict",
         "one of %s options registered in advance" % V["verdicts"],
         "thresholds were never revisited")
    for i in (0, 1):
        arrow(ax, xs[i] + W + 0.004, y3 + hh / 2, xs[i + 1] - 0.004, y3 + hh / 2)

    ax.add_patch(FancyBboxPatch((0.026, -0.052), 0.959, 0.040,
                                boxstyle="round,pad=0.003,rounding_size=0.006",
                                linewidth=0.7, facecolor="#EEF3FB",
                                edgecolor=EDGE["blue"], zorder=3, clip_on=False))
    ax.text(0.505, -0.026, V["verdict"], fontsize=7.6, weight="bold", color=NAVY,
            ha="center", va="center", zorder=4, clip_on=False)
    ax.set_ylim(-0.062, 1.0)

    for ext, dpi in ((".pdf", 300), (".png", 300)):
        fig.savefig(OUT + STEM + ext, dpi=dpi, bbox_inches="tight",
                    facecolor="white")
    fig.savefig(OUT + STEM + "_600dpi.png", dpi=600, bbox_inches="tight",
                facecolor="white")
    plt.close(fig)

    # ---- verify the values survived the render -------------------------
    import pymupdf
    doc = pymupdf.open(OUT + STEM + ".pdf")
    text = " ".join(" ".join(p.get_text().split()) for p in doc)
    MUST = [V["elements"], V["complete"], V["chrom"], V["B"], V["seed"], V["spot"],
            V["checks"], V["verdict"], "ISEScan 1.7.3", "HMMER 3.3.2", "BLAST+ 2.17.0",
            "FragGeneScan 1.32", "OUTCOME SEALED", "COHORT LANE", "OUTCOME LANE"]
    miss = [m for m in MUST if m not in text]
    print("values present in the rendered PDF: %d of %d%s"
          % (len(MUST) - len(miss), len(MUST), "" if not miss else "  MISSING %s" % miss))
    if miss:
        raise SystemExit("a value did not survive the render")
    for ext in (".pdf", ".png", "_600dpi.png"):
        p = OUT + STEM + ext
        print("  %-14s %8.0f KB  %s" % (ext, os.path.getsize(p) / 1024, sha(p)[:20] + "..."))


if __name__ == "__main__":
    main()
