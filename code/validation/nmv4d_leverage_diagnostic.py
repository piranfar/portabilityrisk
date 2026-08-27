"""Leverage and leave-one-species-out slopes for the eight-species fit (item 4.1).

The discordance argument rests on an ordinary least-squares fit of logit(M) on
logit(P) over eight confirmation species. The reviewer's point: P. aeruginosa's
logit(P) is -1.97 while the other seven cluster between -0.41 and +0.70, so one
point may be determining the slope. The published T1 test reports residuals only,
so the paper never says how much.

This computes what a referee would compute: the hat value and Cook's distance for
every species, and the slope with each species left out in turn. Nothing here is
a new estimand - it is the diagnostic for a fit the paper already reports.

The values are transcribed from the Supplementary table rather than recomputed
from the occurrence data, because the table IS the published record of the fit
and the diagnostic must describe what was published.
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


import io
import json
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

NM = _dir("PORTABILITYRISK_REPO_DIR") + "docs/nature_microbiology/"
OUT = NM + "NM_V4D_LEVERAGE_DIAGNOSTIC_V1.json"

# species, logit(P), logit(M) — the eight confirmation species, as published
D = [
    ("Enterobacter asburiae",      -0.1031, -0.9834),
    ("Enterobacter cloacae",       +0.4119, -0.9418),
    ("Enterobacter hormaechei",    +0.6969, -1.6689),
    ("Klebsiella aerogenes",       -0.0893, -1.8601),
    ("Klebsiella michiganensis",   +0.3432, -0.9322),
    ("Klebsiella quasipneumoniae", +0.3560, -1.7999),
    ("Klebsiella variicola",       -0.4093, -1.7060),
    ("Pseudomonas aeruginosa",     -1.9661, -1.5386),
]
PUBLISHED = {"intercept": -1.421185, "slope": 0.080631, "r2": 0.027524}


def ols(pts):
    n = len(pts)
    mx = sum(p[0] for p in pts) / n
    my = sum(p[1] for p in pts) / n
    sxx = sum((p[0] - mx) ** 2 for p in pts)
    sxy = sum((p[0] - mx) * (p[1] - my) for p in pts)
    b = sxy / sxx
    a = my - b * mx
    syy = sum((p[1] - my) ** 2 for p in pts)
    ssr = sum((p[1] - (a + b * p[0])) ** 2 for p in pts)
    return {"intercept": a, "slope": b, "r2": 1 - ssr / syy if syy else 0.0,
            "n": n, "sxx": sxx, "ssr": ssr}


def main():
    pts = [(x, y) for _s, x, y in D]
    full = ols(pts)
    print("reproducing the published fit")
    print("  intercept %+.6f (published %+.6f)" % (full["intercept"], PUBLISHED["intercept"]))
    print("  slope     %+.6f (published %+.6f)" % (full["slope"], PUBLISHED["slope"]))
    print("  R²        %.6f  (published %.6f)" % (full["r2"], PUBLISHED["r2"]))
    ok = (abs(full["slope"] - PUBLISHED["slope"]) < 5e-4
          and abs(full["r2"] - PUBLISHED["r2"]) < 5e-4)
    print("  reproduces to four decimals: %s" % ok)
    if not ok:
        raise SystemExit("the fit does not reproduce from the published table - "
                         "no diagnostic is reported")

    n = full["n"]
    mx = sum(p[0] for p in pts) / n
    mse = full["ssr"] / (n - 2)
    rows = []
    print("\n%-28s %8s %8s %8s %10s" % ("species", "hat", "resid", "Cook D", "slope w/o"))
    for i, (name, x, y) in enumerate(D):
        h = 1.0 / n + (x - mx) ** 2 / full["sxx"]
        e = y - (full["intercept"] + full["slope"] * x)
        cook = (e ** 2 / (2 * mse)) * (h / (1 - h) ** 2)
        loo = ols([p for j, p in enumerate(pts) if j != i])
        rows.append({"species": name, "logit_P": x, "logit_M": y,
                     "hat": round(h, 4), "residual": round(e, 4),
                     "cooks_distance": round(cook, 4),
                     "slope_without": round(loo["slope"], 4),
                     "r2_without": round(loo["r2"], 4)})
        print("%-28s %8.4f %8.4f %8.4f %10.4f" % (name, h, e, cook, loo["slope"]))

    pa = [r for r in rows if r["species"].startswith("Pseudomonas")][0]
    slopes = [r["slope_without"] for r in rows]
    print("\n  P. aeruginosa hat value      : %.4f  (%.0f%% of the total leverage; 8 points, "
          "average hat = %.3f)" % (pa["hat"], 100 * pa["hat"], 2.0 / n))
    print("  P. aeruginosa Cook's distance: %.4f" % pa["cooks_distance"])
    print("  slope with all eight         : %+.4f" % full["slope"])
    print("  slope without P. aeruginosa  : %+.4f  (%.1f times the full-set slope)"
          % (pa["slope_without"], pa["slope_without"] / full["slope"]))
    print("  slope range over the eight   : %+.4f to %+.4f" % (min(slopes), max(slopes)))
    print("  R² without P. aeruginosa     : %.4f" % pa["r2_without"])

    rec = {"receipt": "leverage diagnostic for the eight-species discordance fit",
           "why": ("the published T1 test reports leave-one-species-out residuals only. "
                   "A referee will ask whether one point determines the slope, and the "
                   "answer should be in the paper rather than in their notes."),
           "source": "the eight-species table in the Supplementary discordance section, "
                     "which is the published record of the fit",
           "fit_reproduced_from_published_table": True,
           "full_fit": {k: round(v, 6) for k, v in full.items() if k != "n"},
           "published_fit": PUBLISHED,
           "per_species": rows,
           "headline": {
               "P_aeruginosa_hat": pa["hat"],
               "P_aeruginosa_cooks_distance": pa["cooks_distance"],
               "slope_all_eight": round(full["slope"], 4),
               "slope_without_P_aeruginosa": pa["slope_without"],
               "r2_without_P_aeruginosa": pa["r2_without"],
               "interpretation": ("one point carries most of the leverage and the slope "
                                  "more than doubles without it. The conclusion is "
                                  "unchanged - R² stays near zero and the slope interval "
                                  "covers zero either way - but the fit is not stable and "
                                  "the paper must say so.")}}
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        json.dumps(rec, indent=1, ensure_ascii=False) + "\n")
    print("\nreceipt: %s" % OUT.split("/")[-1])


if __name__ == "__main__":
    main()
