"""Renumber the figures for the V5 Resource, to first-citation order.

The Resource restructure moved sections, so the figures are no longer cited in
numerical order: main [1, 2, 5, 3, 4, 6] and Extended Data [2, 3, 4, 1, 5]. The
numbers follow the text, not the other way round, so the figures are renumbered.

The mapping is DERIVED from where the manuscript first cites each figure, not
chosen. Substitution happens in one pass through a placeholder so that a figure
already renumbered is never renumbered again - mapping 5 to 3 and then 3 to 4
one after the other would silently destroy both.

The V4 figure directory is left untouched. V5 gets its own.
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
import re
import shutil
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

NM = _dir("PORTABILITYRISK_REPO_DIR") + "docs/nature_microbiology/"
MS = NM + "PORTABILITYRISK_MANUSCRIPT_V5_RESOURCE.md"
SI = NM + "PORTABILITYRISK_SUPPLEMENTARY_INFORMATION_V5.md"
SRC = NM + "figures/submission_v4/"
DST = NM + "figures/submission_v5/"
REC = NM + "PORTABILITYRISK_SUBMISSION_FIGURE_RECEIPT_V7.json"

STEMS = {
    ("main", 1): "Figure_1_matched_family_host_contrast",
    ("main", 2): "Figure_2_occurrence_not_gene",
    ("main", 3): "Figure_3_spatial_and_structural",
    ("main", 4): "Figure_4_genome_background_null",
    ("main", 5): "Figure_5_five_class_architecture",
    ("main", 6): "Figure_6_plasmid_fraction_non_redundant",
    ("ed", 1): "Extended_Data_Fig_1_cohort_and_denominators",
    ("ed", 2): "Extended_Data_Fig_2_marker_type_and_IS_family",
    ("ed", 3): "Extended_Data_Fig_3_cargo_convergence",
    ("ed", 4): "Extended_Data_Fig_4_four_denominators",
    ("ed", 5): "Extended_Data_Fig_5_registration_and_verification",
}
EXTS = (".pdf", ".png", "_600dpi.png")
ED_RX = re.compile(r"Extended\s+Data\s+Fig\.\s*(\d)")
MAIN_RX = re.compile(r"(?<!Data )Fig\.\s*(\d)")


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def first_order(text, kind):
    """Figure numbers in order of first citation, legends excluded."""
    body = text[:text.index("## Figure legends")] if "## Figure legends" in text else text
    if kind == "ed":
        seq = [int(m.group(1)) for m in ED_RX.finditer(body)]
    else:
        seq = [int(m.group(1)) for m in MAIN_RX.finditer(ED_RX.sub("#" * 22, body))]
    out = []
    for x in seq:
        if x not in out:
            out.append(x)
    return out


def main():
    t = io.open(MS, encoding="utf-8").read()
    s = io.open(SI, encoding="utf-8").read()

    maps = {}
    for kind in ("main", "ed"):
        order = first_order(t, kind)
        if sorted(order) != list(range(1, len(STEMS_of(kind)) + 1)):
            raise SystemExit("%s: cited %s but %d figures exist"
                             % (kind, order, len(STEMS_of(kind))))
        maps[kind] = {old: new for new, old in enumerate(order, 1)}
        moved = [(o, n) for o, n in sorted(maps[kind].items()) if o != n]
        print("  %-4s citation order %s -> %s"
              % (kind, order, ", ".join("%d->%d" % p for p in moved) or "already ascending"))

    def rewrite(text):
        # one pass, via a placeholder, so a renumbered figure is never re-hit
        text = ED_RX.sub(lambda m: "\x01ED%d\x01" % maps["ed"][int(m.group(1))], text)
        text = MAIN_RX.sub(lambda m: "\x01MN%d\x01" % maps["main"][int(m.group(1))], text)
        text = re.sub(r"\x01ED(\d)\x01", r"Extended Data Fig. \1", text)
        text = re.sub(r"\x01MN(\d)\x01", r"Fig. \1", text)
        return text

    # legends carry their own numbers and must move with them
    t2, s2 = rewrite(t), rewrite(s)

    # reorder the legend blocks so the section reads 1..n
    def reorder_legends(text, pat, n):
        blocks = {}
        for m in re.finditer(pat, text):
            end = text.find("\n\n", m.start())
            end = len(text) if end < 0 else end
            blocks[int(m.group(1))] = (m.start(), end, text[m.start():end])
        if len(blocks) != n:
            return text
        spans = sorted((v[0], v[1]) for v in blocks.values())
        ordered = "\n\n".join(blocks[k][2] for k in sorted(blocks))
        return text[:spans[0][0]] + ordered + text[spans[-1][1]:]

    t2 = reorder_legends(t2, r"\*\*Fig\. (\d) \|", 6)
    t2 = reorder_legends(t2, r"\*\*Extended Data Fig\. (\d) \|", 5)

    for kind in ("main", "ed"):
        o = first_order(t2, kind)
        print("  %-4s after renumbering: %s  %s"
              % (kind, o, "ascending" if o == sorted(o) else "STILL BROKEN"))
        if o != sorted(o):
            raise SystemExit("renumbering did not produce ascending order")
    leg = re.findall(r"\*\*Fig\. (\d) \|", t2)
    edleg = re.findall(r"\*\*Extended Data Fig\. (\d) \|", t2)
    print("  legend order: main %s | extended %s" % (leg, edleg))
    if leg != sorted(leg) or edleg != sorted(edleg):
        raise SystemExit("legend blocks are not in order")

    # ---- the files ------------------------------------------------------
    os.makedirs(DST, exist_ok=True)
    out = {}
    for (kind, old), stem in STEMS.items():
        new = maps[kind][old]
        prefix = "Figure_%d_" % new if kind == "main" else "Extended_Data_Fig_%d_" % new
        tail = stem.split("_", 2 if kind == "main" else 4)[-1]
        for e in EXTS:
            src = SRC + stem + e
            if not os.path.isfile(src):
                continue
            dst = DST + prefix + tail + e
            shutil.copyfile(src, dst)
            out["docs/nature_microbiology/figures/submission_v5/" + prefix + tail + e] = sha(dst)
    print("\n  copied %d files into figures/submission_v5/" % len(out))

    io.open(MS, "w", encoding="utf-8", newline="\n").write(t2)
    io.open(SI, "w", encoding="utf-8", newline="\n").write(s2)

    rec = {"receipt": "submission figures, v7 (V5 Resource numbering)",
           "why": ("the Resource restructure changed the order in which figures are "
                   "first cited, so they are renumbered to follow the text. The "
                   "mapping is derived from the manuscript, not chosen."),
           "map": {k: {str(o): n for o, n in sorted(v.items())} for k, v in maps.items()},
           "v4_figures_untouched": True,
           "extended_data_caveat": ("Nature Microbiology lists Extended Data for "
                                    "Articles and Brief Communications but not for "
                                    "Resources. Whether Extended Data is permitted "
                                    "here must be confirmed with the editorial "
                                    "office; if not, these five become Supplementary "
                                    "Figures."),
           "outputs_sha256": out}
    io.open(REC, "w", encoding="utf-8", newline="\n").write(
        json.dumps(rec, indent=1, ensure_ascii=False) + "\n")
    print("  manuscript : %s" % sha(MS))
    print("  supplement : %s" % sha(SI))
    print("  receipt    : %s" % os.path.basename(REC))


def STEMS_of(kind):
    return [k for k in STEMS if k[0] == kind]


if __name__ == "__main__":
    main()
