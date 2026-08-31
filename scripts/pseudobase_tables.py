#!/usr/bin/env python3
"""The two PseudoBase tables, computed from the per-target bench CSVs.

WHY NOT JUST USE bench_report's report.csv
Its F1 column moves with --ref: the plain ViennaRNA MFE scores 0.581 under `--ref mfe` and
0.593 under `--ref pkiss` on the same predictions, because the reference tool determines
which targets enter the paired comparison. That is fine for the dF1 column it exists to
produce, and wrong for a table of absolute scores, which must not depend on who the
comparison partner was. Everything here is recomputed from pred_db with one fixed rule.

DEFINITIONS, stated once so the caption can point at them:
  F1, PPV, Sty   over base pairs, per target, then averaged within sequence-identity
                 cluster, then across clusters (288 clusters over 391 positives)
  PKrec          |predicted pairs INTERSECT native minority-bracket pairs| / |native
                 minority-bracket pairs|. Recomputed here for every tool, including the
                 external ones whose own CSV carries the inflated symmetric convention.
  sens / fpr     fraction of positives / negatives whose PREDICTION contains a crossing
  J              sens - fpr
  F1neg          mean F1 on the 489 pseudoknot-free negatives

THE STRATUM SPLIT is defined by the data rather than by a side file: a target is in
"covariation available" iff the covariation field produced a row for it. That is exactly
101 of 391 (26%), and it is the honest denominator -- the pooled row for covariation in
bench_report's output is silently computed on those 101 alone, which is why it prints an
identical value in the 'all' and 'ok' blocks and nan for fpr.

    python pseudobase_tables.py final3_pos_bench.csv final3_neg_bench.csv
writes pseudobase_main.tex (body) and pseudobase_stratum.tex (appendix).
"""
import csv, sys, os, statistics as st
from collections import defaultdict

from pkutils import find, dbp, native_pk_pairs, crossing_pairs  # noqa: E402

PKB_TSV = find("pseudobase.tsv")
PKB_CL = find("pseudobase.clusters.tsv")

ROWS = [
    (r"\textit{ensemble $F_1$ oracle}", "ours_oracle", True),
    (r"\textit{ensemble PK oracle}",    "ours_pk_oracle", True),
    (r"\textit{random $F_1$ oracle}",   "random_oracle", True),
    (r"\textit{random PK oracle}",      "random_pk_oracle", True),
    (r"pKiss",                          "pkiss", False),
    (r"IPknot",                         "ipknot", False),
    (r"ProbKnot",                       "probknot", False),
    (r"RNAfold MFE",                    "mfe", False),
    (r"ours, covariation (MI)",         "mi_mea0.5", False),
    (r"ours, LM field",                 "lm_mea0.5", False),
    (r"ours, pair count",               "ours_count", False),
    (r"ours, closed-form energy",       "ours_energy", False),
    (r"ours, ensemble centroid",        "ours_centroid", False),
]


def load_native():
    nat, natpk = {}, {}
    for ln in open(PKB_TSV):
        if ln.startswith("#") or not ln.strip():
            continue
        f = ln.rstrip("\n").split("\t")
        if len(f) >= 3:
            nat[f[0]] = dbp(f[2]); natpk[f[0]] = native_pk_pairs(f[2])
    return nat, natpk


def load_clusters():
    cl = {}
    for ln in open(PKB_CL):
        if ln.startswith("#") or not ln.strip():
            continue
        f = ln.rstrip("\n").split("\t")
        cl[f[0]] = f[1] if len(f) > 1 else f[0]
    return cl


def cmean(per_target, cl):
    """Average within cluster, then across clusters -- never a plain target mean, because
    PseudoBase contains near-duplicate viral pseudoknots that would otherwise be counted
    several times each."""
    per = defaultdict(list)
    for tid, v in per_target.items():
        per[cl.get(tid, tid)].append(v)
    return st.mean([st.mean(v) for v in per.values()]) if per else float("nan")


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    nat, natpk = load_native(); cl = load_clusters()

    pos = defaultdict(dict)   # tool -> id -> row
    for r in csv.DictReader(open(sys.argv[1])):
        pos[r["tool"]][r["id"]] = r
    neg = defaultdict(dict)
    for r in csv.DictReader(open(sys.argv[2])):
        neg[r["tool"]][r["id"]] = r

    cov_ids = set(pos.get("mi_mea0.5", {}))          # data-defined stratum
    all_ids = set(pos.get("mfe", {}))
    nocov_ids = all_ids - cov_ids

    def stats(tool, ids):
        rows = {i: r for i, r in pos.get(tool, {}).items() if i in ids}
        rows = {i: r for i, r in rows.items() if r.get("pred_db") not in (None, "", "NA")}
        if not rows:
            return None
        f1 = {}; ppv = {}; sty = {}; pk = {}; xing = {}
        for i, r in rows.items():
            pred = dbp(r["pred_db"]); ref = nat[i]
            tp = len(pred & ref)
            ppv[i] = tp / len(pred) if pred else 0.0
            sty[i] = tp / len(ref) if ref else 0.0
            f1[i] = (2 * ppv[i] * sty[i] / (ppv[i] + sty[i])) if (ppv[i] + sty[i]) else 0.0
            if natpk[i]:
                pk[i] = len(pred & natpk[i]) / len(natpk[i])
            xing[i] = 1.0 if crossing_pairs(pred) else 0.0
        # negatives are never stratified by covariation availability -- the stratum is a
        # property of the positive set, and every method sees all 489 either way
        nrows = {i: r for i, r in neg.get(tool, {}).items()
                 if r.get("pred_db") not in (None, "", "NA")}
        fpr = st.mean([1.0 if crossing_pairs(dbp(r["pred_db"])) else 0.0
                       for r in nrows.values()]) if nrows else float("nan")
        f1n = st.mean([float(r["f1"]) for r in nrows.values() if r.get("f1")]) if nrows else float("nan")
        sens = cmean(xing, cl)
        return dict(n=len(rows), f1=cmean(f1, cl), ppv=cmean(ppv, cl), sty=cmean(sty, cl),
                    pk=cmean(pk, cl) if pk else float("nan"),
                    sens=sens, fpr=fpr, J=sens - fpr, f1n=f1n)

    def fmt(v, sign=False):
        if v != v:
            return "---"
        return (f"{v:+.3f}" if sign else f"{v:.3f}")

    def block(ids, cols, caption, label, with_neg, min_cov=0.9):
        """min_cov guards the pooled table: a selector that produced rows for only 26% of the
        benchmark cannot sit in the same column as one that covered all of it. bench_report
        does not guard this, which is why its pooled mi row equals its 'ok' row exactly."""
        out = [r"\begin{table}[!t]", r"\caption{" + caption + "}", r"\label{" + label + "}",
               r"\centering", r"\footnotesize", r"\setlength{\tabcolsep}{5pt}",
               r"\begin{tabular}{l" + "c" * len(cols) + "}", r"\toprule",
               "Method & " + " & ".join(cols) + r" \\", r"\midrule"]
        prev = True
        for lab, tool, it in ROWS:
            s = stats(tool, ids)
            if s is None or s["n"] < min_cov * len(ids):
                continue
            if prev and not it:
                out.append(r"\midrule")
            prev = it
            # PKrec sits immediately after F1: the two headline quantities are the whole-
            # structure score and the pseudoknot-specific one, and reading the gap between
            # them is the point of the table. PPV/Sty are the decomposition of F1 and follow.
            cells = [fmt(s["f1"]), fmt(s["pk"]), fmt(s["ppv"]), fmt(s["sty"])]
            if with_neg:
                cells += [fmt(s["sens"]), fmt(s["fpr"]), fmt(s["J"], True), fmt(s["f1n"])]
            if it:
                cells = [rf"\textit{{{c}}}" for c in cells]
            out.append(f"{lab} & " + " & ".join(cells) + r" \\")
        out += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]
        return out

    def block2(ids_a, ids_b, caption, label):
        cols = ["$F_1$", "PKrec", "PPV"]
        out = [r"\begin{table}[!t]", r"\caption{" + caption + "}", r"\label{" + label + "}",
               r"\centering", r"\footnotesize", r"\setlength{\tabcolsep}{5pt}",
               r"\begin{tabular}{l" + "c" * (2 * len(cols)) + "}", r"\toprule",
               r"& \multicolumn{3}{c}{covariation available ($n=%d$)}"
               r" & \multicolumn{3}{c}{not available ($n=%d$)} \\" % (len(ids_a), len(ids_b)),
               r"\cmidrule(lr){2-4}\cmidrule(lr){5-7}",
               "Method & " + " & ".join(cols + cols) + r" \\", r"\midrule"]
        prev = True
        for lab, tool, it in ROWS:
            sa, sb = stats(tool, ids_a), stats(tool, ids_b)
            if sa is None and sb is None:
                continue
            if prev and not it:
                out.append(r"\midrule")
            prev = it
            cells = []
            for sx, idx in ((sa, ids_a), (sb, ids_b)):
                if sx is None or sx["n"] < 0.9 * len(idx):
                    cells += ["---"] * 3
                else:
                    cells += [fmt(sx["f1"]), fmt(sx["pk"]), fmt(sx["ppv"])]
            if it:
                cells = [rf"\textit{{{c}}}" for c in cells]
            out.append(f"{lab} & " + " & ".join(cells) + r" \\")
        out += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]
        return out

    HDR = ["% GENERATED -- do not hand-edit. Rebuild with:",
           "%   python pseudobase_tables.py final3_pos_bench.csv final3_neg_bench.csv"]
    L = list(HDR)
    L += block(all_ids,
               ["$F_1$", "PKrec", "PPV", "Sty", "sens", "fpr", "$J$", "$F_1^{\\mathrm{neg}}$"],
               r"PseudoBase, all $391$ pseudoknotted targets ($288$ sequence-identity clusters) "
               r"and $489$ curated nested negatives. Values are cluster-level means; PKrec is "
               r"minority-bracket crossing-pair recovery, recomputed here for every method from "
               r"its predicted structure. Oracles (italic) read the native structure and are "
               r"ceilings, not predictors. The covariation rows are absent because covariation "
               r"exists for only $101$ of these targets; see Table~\ref{tab:pkbstratum}.",
               "tab:pseudobase", True)
    open("pseudobase_main.tex", "w").write("\n".join(L) + "\n")
    L = list(HDR)
    L += block2(cov_ids, nocov_ids,
               r"PseudoBase split by whether covariation is \emph{available}: an Infernal "
               r"search finds a family alignment for $101$ of $391$ targets ($26\%$). Where it "
               r"is available, covariation is the strongest selector we have, ahead of pKiss; "
               r"where it is not, no alignment-based selector exists at all, and its absence "
               r"there is a coverage statement rather than a quality one. Every other method is "
               r"shown on both halves so the two columns are comparable. This is the table "
               r"behind ``covariation is the strongest selector we have, with $26\%$ reach''.",
               "tab:pkbstratum")
    open("pseudobase_stratum.tex", "w").write("\n".join(L) + "\n")
    print("wrote pseudobase_main.tex and pseudobase_stratum.tex", file=sys.stderr)

    # --- reconciliation to stderr, so the numbers can be argued with -------------------
    print("\n[pseudobase_tables] covariation available on "
          f"{len(cov_ids)}/{len(all_ids)} targets ({100*len(cov_ids)/len(all_ids):.1f}%)",
          file=sys.stderr)
    print(f"{'tool':26s}{'F1 all':>9s}{'F1 cov':>9s}{'F1 nocov':>10s}", file=sys.stderr)
    for _, tool, _ in ROWS:
        a, b, c = stats(tool, all_ids), stats(tool, cov_ids), stats(tool, nocov_ids)
        if a is None:
            continue
        print(f"  {tool:24s}{a['f1']:9.3f}{b['f1'] if b else float('nan'):9.3f}"
              f"{c['f1'] if c else float('nan'):10.3f}", file=sys.stderr)


if __name__ == "__main__":
    main()
