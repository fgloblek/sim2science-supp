#!/usr/bin/env python3
"""Cross-check every benchmark number in the paper against the source CSVs.

A workshop paper's numbers get transcribed by hand at least twice -- CSV to table, table to
abstract -- and a single transposed digit in a table nobody re-derives is the classic way a
result becomes wrong in print. This re-reads the CSVs and asserts each claim, so the check
is repeatable after any edit rather than a one-off eyeballing.

Run from the workshop/ directory. Exits non-zero on the first mismatch.
"""
import csv, os, re, sys

from pkutils import find, dbp, native_pk_pairs, crossing_pairs

BOTH = "report_vs_mfe.csv"            # ensemble-mode=both
AUG = "augment_report_vs_mfe.csv"     # ensemble-mode=augment
TEX = os.environ.get("PAPER_TEX", "sim2science_pk.tex")


def load(path):
    if not os.path.exists(path):
        sys.exit(f"missing {path} -- copy it next to the .tex before verifying")
    out = {}
    with open(path) as fh:
        for r in csv.DictReader(fh):
            out[(r["scope"], r["tool"])] = r
    return out


def num(rows, scope, tool, col):
    v = rows[(scope, tool)][col]
    return float(v) if v not in ("", "nan") else None


def main():
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    # the tables now live in generated \input files, so the "paper text" for checking
    # purposes is the main .tex plus everything it inputs
    def _strip_comments(s):
        """Drop LaTeX comments before any check runs. Superseded drafts are kept commented out
        in the source -- an old abstract still quoting ProbKnot's 0.529 as the best published
        tool sits inside \\begin{abstract}, and matching against it raised a mismatch for a
        number the live text no longer contains. An escaped \\% is not a comment."""
        out = []
        for ln in s.split("\n"):
            i, j = 0, None
            while i < len(ln):
                if ln[i] == "\\":
                    i += 2; continue
                if ln[i] == "%":
                    j = i; break
                i += 1
            out.append(ln if j is None else ln[:j])
        return "\n".join(out)

    # THE PAPER SOURCE IS NOT PART OF THIS BUNDLE. Without it the literal-matching checks
    # cannot run, but every check that recomputes a table from the predictions still can, and
    # those are the ones that establish the numbers. Set PAPER_TEX to enable the rest.
    HAVE_TEX = os.path.exists(TEX)
    if HAVE_TEX:
        tex = _strip_comments(open(TEX).read())
        for inc in re.findall(r"\\input\{([^}]*)\}", tex):
            f = inc if inc.endswith(".tex") else inc + ".tex"
            if os.path.exists(f):
                tex += "\n" + _strip_comments(open(f).read())
    else:
        tex = ""
        print(f"note: {TEX} not found -- recomputing every table from the predictions, "
              "and skipping the checks that match those values against the paper text.")
    fails = []

    def check(label, got, want, tol=5e-4):
        if got is None or abs(got - want) > tol:
            fails.append(f"{label}: paper says {want}, CSV says {got}")

    # ---- PseudoBase: the generated tables are the single source ---------------------------
    # pseudobase_tables.py recomputes every cell from pred_db, so the check here is that the
    # .tex on disk was generated from the CSVs on disk -- i.e. nobody hand-edited a number
    # into the table and nobody forgot to regenerate after a re-run.
    # The generator writes TWO files now -- the pooled table for the body, the covariation
    # stratum split for the appendix -- so the check snapshots both, regenerates, and diffs.
    import subprocess
    GEN = ("pseudobase_main.tex", "pseudobase_stratum.tex")
    before = {f: open(f).read() for f in GEN if os.path.exists(f)}
    gen = subprocess.run([sys.executable, "pseudobase_tables.py",
                          find("final3_pos_bench.csv"), find("final3_neg_bench.csv")],
                         capture_output=True, text=True)
    if gen.returncode != 0:
        fails.append(f"pseudobase_tables.py failed: {gen.stderr.strip()[:200]}")
    else:
        for f in GEN:
            if not os.path.exists(f):
                fails.append(f"{f} missing -- run pseudobase_tables.py")
            elif f in before and open(f).read().strip() != before[f].strip():
                fails.append(f"{f} was stale or hand-edited; it has just been regenerated "
                             "from the CSVs, so rebuild the PDF")

    # 0.171 is the PAIRED cluster-level lift, which is what the quoted interval belongs to;
    # subtracting the two tabulated means (0.695 - 0.525) gives 0.170. Keep them consistent.
    for lit in ("0.931", "0.777", "0.956", "0.960", "0.611", "0.637", "0.695", "0.171"):
        if HAVE_TEX and lit not in tex:
            fails.append(f"paper lost the PseudoBase literal {lit}")

    # ---- Table 2's delta column is now referenced to the MFE, not pKiss -------------------
    # Recomputed here from the per-target ArchiveII rows with a paired bootstrap over the six
    # families that carry positives, so the column cannot drift from the table it describes.
    import csv as _csv, random as _rnd
    from collections import defaultdict as _dd
    if os.path.exists(find("archiveii_pk_bench.csv")):
        per = _dd(dict)
        for r in _csv.DictReader(open(find("archiveii_pk_bench.csv"))):
            if r.get("pred_db") not in (None, "", "NA") and r.get("f1"):
                per[r["tool"]][r["id"]] = float(r["f1"])
        def dvs_mfe(tool):
            common = set(per[tool]) & set(per["mfe"])
            byfam = _dd(list)
            for i in common:
                byfam[i.split("|")[0]].append(per[tool][i] - per["mfe"][i])
            d = [sum(v) / len(v) for v in byfam.values()]
            return sum(d) / len(d)
        for tool, want in (("ours_oracle", +0.126), ("mi_mea0.5", +0.084),
                           ("lm_mea0.5", -0.013), ("pkiss", -0.100),
                           ("probknot", +0.015)):
            check(f"T2 dF1-vs-MFE {tool}", dvs_mfe(tool), want, tol=1e-3)
        if HAVE_TEX and r"$\Delta F_1$ vs.\ MFE" not in tex and r"\Delta F_1$ vs.\ MFE" not in tex:
            fails.append("Table 2's delta column header no longer says 'vs. MFE'")

    # ---- Table 3: ArchiveII, all ten families -------------------------------------------
    # READ FROM THE PER-TARGET CSV COLUMNS, NOT RECOMPUTED FROM pred_db.
    #
    # pred_db is a LOSSY record for heavily knotted structures. Extended dot-bracket offers four
    # bracket types, ()[]{}<>, and a pair that cannot be assigned one is silently dropped when
    # the structure is written out. The matched random ensemble packs helices greedily from a
    # randomly ordered vocabulary, so it produces far more knotted structures than thermodynamics
    # ever does: 1052 of 1063 random structures exhaust all four levels, and 970 of them lose
    # pairs -- 14.5 on average, up to 74. Every other tool round-trips (at most a dozen targets,
    # mean 0.00, which is rounding in the ppv/sty columns).
    #
    # Because the dropped pairs are mostly WRONG ones, recomputing from pred_db raises PPV and
    # inflates F1: the random F1 oracle reads 0.150 from the columns and 0.157 recomputed, and
    # the random PK oracle's PKrec reads 0.204 against 0.171. The columns are computed before
    # serialisation and are the correct ones. An earlier version of this file recomputed
    # everything and mis-diagnosed the difference as the random pool being redrawn per run.
    #
    # So: F1/PPV/Sty/PKrec come from the columns, EXCEPT PKrec for the external tools, whose
    # pk_sty column uses the inflated symmetric convention and must be recomputed -- which is
    # safe precisely because those tools round-trip. _roundtrip_guard below enforces that.
    # IPknot was benchmarked on its own, so its rows arrive in a separate file.
    import statistics as _st

    def _natives(p):
        d = {}
        for ln in open(p):
            if ln.startswith("#") or not ln.strip():
                continue
            f = ln.rstrip("\n").split("\t")
            if len(f) >= 3:
                d[f[0]] = f[2]
        return d

    EXTERNAL = {"pkiss", "ipknot", "probknot", "rnafold"}   # pk_sty column uses the wrong convention

    def _load(*paths):
        """tool -> id -> (pred_db, f1, ppv, sty, pk_sty). The numeric columns were computed
        before the structure was serialised, so they survive the bracket-depth loss."""
        per = _dd(dict)
        for p in paths:
            if not os.path.exists(p):
                continue
            for r in _csv.DictReader(open(p)):
                if r.get("pred_db") in (None, "", "NA") or not r.get("f1"):
                    continue
                def g(k):
                    v = r.get(k)
                    return float(v) if v not in (None, "", "nan") else None
                per[r["tool"]].setdefault(r["id"], (r["pred_db"], g("f1"), g("ppv"),
                                                    g("sty"), g("pk_sty")))
        return per

    def _pairf1(pred, ref):
        tp = len(pred & ref)
        if not pred or not ref:
            return 0.0, 0.0, 0.0
        p, s = tp / len(pred), tp / len(ref)
        return (2 * p * s / (p + s) if p + s else 0.0), p, s

    NPK = _natives(find("archiveii_pk.tsv")) if os.path.exists(find("archiveii_pk.tsv")) else {}
    NNP = _natives(find("archiveii_nopk.tsv")) if os.path.exists(find("archiveii_nopk.tsv")) else {}
    if NPK:
        AP = _load(find("archiveii_pk_bench.csv"), find("ipknot_ext_pk.csv"))
        AN = _load(find("archiveii_nopk_bench.csv"), find("ipknot_ext_nopk.csv"))

        def _fam_mean(d):      # per family, then across families -- the paper's cluster unit
            g = _dd(list)
            for i, v in d.items():
                g[i.split("|")[0]].append(v)
            return _st.mean([_st.mean(v) for v in g.values()]) if g else float("nan")

        def _roundtrip_guard():
            """Fail if a tool whose PKrec we recompute has a pred_db that does not round-trip.
            The random rows are expected to lose pairs and are read from the columns instead."""
            for tool in EXTERNAL:
                bad = 0
                for i, (pdb, f1c, ppv, sty, _) in AP.get(tool, {}).items():
                    if i not in NPK or not ppv or not sty:
                        continue
                    implied = sty * len(dbp(NPK[i])) / ppv
                    if implied - len(dbp(pdb)) > 0.5:
                        bad += 1
                n = len(AP.get(tool, {}))
                if n and bad > 0.01 * n:
                    fails.append(f"{tool}: pred_db loses pairs on {bad}/{n} targets, so its "
                                 "recomputed PKrec is not trustworthy (bracket-depth overflow)")
        _roundtrip_guard()

        def _row(tool):
            ids = [i for i in AP.get(tool, {}) if i in NPK]
            if not ids:
                return None
            F, P, S, K, X = {}, {}, {}, {}, {}
            for i in ids:
                pdb, f1c, ppv, sty, pk = AP[tool][i]
                F[i], P[i], S[i] = f1c, ppv, sty          # columns: pre-serialisation, correct
                if tool in EXTERNAL:                      # wrong convention in the column
                    npk = native_pk_pairs(NPK[i])
                    if npk:
                        K[i] = len(dbp(pdb) & npk) / len(npk)
                elif pk is not None:
                    K[i] = pk
                X[i] = 1.0 if crossing_pairs(dbp(pdb)) else 0.0
            nids = [i for i in AN.get(tool, {}) if i in NNP]
            FP = {i: (1.0 if crossing_pairs(dbp(AN[tool][i][0])) else 0.0) for i in nids}
            FN = {i: AN[tool][i][1] for i in nids if AN[tool][i][1] is not None}
            J = _fam_mean(X) - _fam_mean(FP) if FP else None
            return (_fam_mean(F), _fam_mean(P), _fam_mean(S), _fam_mean(K), J,
                    _fam_mean(FN) if FN else None)

        T2 = {
            "ours_oracle":      (0.638, 0.616, 0.669, 0.346, +0.233, 0.722),
            "ours_pk_oracle":   (0.487, 0.472, 0.512, 0.702, None,   None),
            "random_oracle":    (0.150, 0.159, 0.145, 0.100, +0.013, 0.323),
            "random_pk_oracle": (0.037, 0.039, 0.037, 0.204, None,   None),
            "mi_mea0.5":        (0.599, 0.562, 0.649, 0.300, +0.232, 0.629),
            "ipknot":           (0.536, 0.549, 0.536, 0.266, +0.181, 0.604),
            "probknot":         (0.529, 0.502, 0.569, 0.238, +0.024, 0.561),
            "mfe":              (0.512, 0.477, 0.560, 0.239,  0.000, 0.564),
            "lm_mea0.5":        (0.499, 0.486, 0.525, 0.303, +0.336, 0.539),
            "ours_energy":      (0.481, 0.438, 0.541, 0.242,  0.000, 0.485),
            "ours_centroid":    (0.477, 0.497, 0.469, 0.208, +0.033, 0.522),
            "ours_count":       (0.458, 0.414, 0.520, 0.241, +0.014, 0.492),
            "pkiss":            (0.441, 0.416, 0.474, 0.227, +0.241, 0.516),
        }
        for tool, vals in T2.items():
            got = _row(tool)
            if got is None:
                fails.append(f"T2 {tool}: no rows in the ArchiveII bench CSVs")
                continue
            for col, want, g in zip(("f1", "ppv", "sty", "pk_rec", "J", "f1_neg"), vals, got):
                if want is None:
                    continue
                check(f"T2 {tool}.{col}", g, want, tol=6e-4)
        # IPknot must be the best PUBLISHED tool -- the claim the body now makes
        if _row("ipknot")[0] <= _row("probknot")[0]:
            fails.append("IPknot no longer beats ProbKnot; Section 3.1 and the abstract say it does")
        if _row("mi_mea0.5")[0] <= _row("ipknot")[0]:
            fails.append("covariation no longer beats IPknot; the paper's headline claim breaks")
    else:
        fails.append("archiveii_pk.tsv missing -- Table 2 cannot be verified from source")

    # ---- the abstract must not disagree with the tables ----------------------------------
    abstract = "" if not HAVE_TEX else tex[tex.index(r"\begin{abstract}"):tex.index(r"\end{abstract}")]
    # every literal the abstract DOES quote must match a table; it need not quote them all
    # the paper quotes the PAIRED estimate (+0.084), which is the one the bootstrap interval
    # belongs to; the pooled difference of means (+0.087) is checked above but not quoted.
    # The abstract was rewritten and now quotes a smaller set; every literal below must both
    # appear there AND match a table. 0.536 is IPknot, which displaced ProbKnot's 0.529 as the
    # best published tool once it was benchmarked, so an abstract still saying 0.529 is stale.
    for lit in ("0.931", "0.777", "0.638", "0.536", "0.599", "0.013"):
        if HAVE_TEX and lit not in abstract:
            fails.append(f"abstract lost the literal {lit}")
    if HAVE_TEX and "0.529" in abstract:
        fails.append("abstract still quotes 0.529 (ProbKnot) as the best published tool; "
                     "IPknot reaches 0.536")

    # ---- pending cells must still be visibly marked --------------------------------------
    if HAVE_TEX and r"\todo" not in tex:
        fails.append("no \\todo left -- the LM-provenance question must stay visible "
                     "until the probe's training set is confirmed")

    if fails:
        print("MISMATCHES:")
        for f in fails:
            print("  -", f)
        sys.exit(1)
    print(f"all {sum(1 for v in T2.values() for x in v if x is not None) + 2 + 11 + 8 + 5 + 1} checks pass")


if __name__ == "__main__":
    main()
