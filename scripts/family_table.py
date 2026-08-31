#!/usr/bin/env python3
"""Per-family ArchiveII table: F1 and crossing-pair recovery, with n and median length.

WHY A SEPARATE SCRIPT
bench_report prints an F1-by-family block to the console but writes only the pooled 'all'
scope to report.csv, and it never breaks PKrec out by family. The pooled numbers are also
dominated by two families -- tmRNA (459) and RNaseP (417) are 82% of the positives -- so the
headline row answers a different question from "does this generalise across families".

INPUT is the per-target bench CSV, which carries one row per (target, tool):
    archiveii_cov/pk_bench.csv        id,n,tool,ppv,sty,f1,pk_sty,pred_db
Family is the part of the id before the first '|' (make_archiveii writes "<fam>|<basename>").

Rows with an empty or NA prediction are dropped per (tool, family) rather than scored as
zero, and the count of survivors is printed, because a tool that TIMED OUT on the long
families -- pKiss does -- would otherwise look merely bad rather than absent.

    python family_table.py archiveii_pk_bench.csv archiveii_pk.tsv [extra.csv ...] > family_table.tex

Extra CSVs are read after the first and merged on (id, tool), which is how a tool run on its
own lands here -- IPknot was benchmarked separately, and the bench CSVs are kept read-only, so
it is passed as ipknot_ext_pk.csv rather than appended to them.

PKrec IS RECOMPUTED HERE, NOT READ FROM THE CSV. The pk_sty column is written by two different
producers that do not agree. multi_select.py uses native_pk_pairs -- the MINORITY-BRACKET arm,
which is the paper's definition. bench_pk_tools.py, which produces every external tool's row,
uses split_nc() instead -- the geometric crossing set, containing BOTH arms of every
pseudoknot. A nested folder is then credited for recovering the primary stem, which is exactly
what thermodynamics already finds, so the external tools' column was inflated and, worse,
inflated by a different amount per family. Measured on ArchiveII: ProbKnot's telomerase cell
read 0.421 where the true value is 0.000, and its RNaseP cell 0.379 against 0.095.
bench_report.py recomputes from pred_db and so was never affected, which is why Table 2 is
correct and this table was not. Taking pred_db plus the native TSV as the only inputs removes
the discrepancy at the source.
"""
import csv, os, sys, statistics as st
from collections import defaultdict

from pkutils import dbp, native_pk_pairs  # noqa: E402

# Families are ordered by size; only these six contain any pseudoknotted target. 5S, SRP and
# tRNA are pseudoknot-free throughout and so have no positives to score.
ORDER = ["tmRNA", "RNaseP", "grp1", "16s", "telomerase", "23s"]

# Row label -> tool key in the CSV. Oracles first, then published tools, then ours.
ROWS = [
    (r"\textit{ensemble $F_1$ oracle}", "ours_oracle", True),
    (r"\textit{ensemble PK oracle}",    "ours_pk_oracle", True),
    (r"\textit{random $F_1$ oracle}",   "random_oracle", True),
    (r"ours, covariation (MI)",         "mi_mea0.5", False),
    (r"IPknot",                         "ipknot", False),
    (r"ProbKnot",                       "probknot", False),
    (r"RNAfold MFE",                    "mfe", False),
    (r"ours, LM field (zero-shot)",     "lm_mea0.5", False),
    (r"ours, closed-form energy",       "ours_energy", False),
    (r"pKiss",                          "pkiss", False),
]


def load_natives(path):
    """id -> native dot-bracket, for the strict PKrec recomputation."""
    nat = {}
    for ln in open(path):
        if ln.startswith("#") or not ln.strip():
            continue
        f = ln.rstrip("\n").split("\t")
        if len(f) >= 3:
            nat[f[0]] = f[2]
    return nat


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    nat = load_natives(sys.argv[2])
    natpk = {i: native_pk_pairs(db) for i, db in nat.items()}
    f1 = defaultdict(list)      # (tool, fam) -> [f1]
    pk = defaultdict(list)      # (tool, fam) -> [pk_sty]
    lens = defaultdict(list)    # fam -> [n]   (one entry per target, from any single tool)
    seen = set()
    # COVARIATION COVERAGE. A target has covariation iff the mi field produced a row for it,
    # i.e. Infernal assigned a family with enough seed sequences. On PseudoBase this is only
    # 26% and is the headline limitation; per family here it says whether covariation's
    # showing on that family rests on the whole family or a fraction of it.
    ids_all = defaultdict(set); ids_mi = defaultdict(set)

    seen_rows = set()
    for path in [sys.argv[1]] + sys.argv[3:]:
      with open(path) as fh:
        for r in csv.DictReader(fh):
            fam = r["id"].split("|")[0]
            if not r.get("f1"):
                continue
            tool = r["tool"]
            if (r["id"], tool) in seen_rows:      # first file wins on a duplicate
                continue
            seen_rows.add((r["id"], tool))
            f1[(tool, fam)].append(float(r["f1"]))
            # strict minority-bracket PKrec, from pred_db -- never the CSV's pk_sty
            if r.get("pred_db") not in (None, "", "NA") and natpk.get(r["id"]):
                npk = natpk[r["id"]]
                pk[(tool, fam)].append(len(dbp(r["pred_db"]) & npk) / len(npk))
            ids_all[fam].add(r["id"])
            if tool == "mi_mea0.5":
                ids_mi[fam].add(r["id"])
            if (r["id"], fam) not in seen:
                seen.add((r["id"], fam)); lens[fam].append(int(r["n"]))

    fams = [f for f in ORDER if f in lens] + sorted(set(lens) - set(ORDER))

    out = []
    out.append(r"\begin{table}[!t]")
    out.append(r"\caption{ArchiveII by family, pseudoknotted targets only. Each cell is "
               r"$F_1$ / PKrec (minority-bracket crossing-pair recovery). tmRNA and RNaseP "
               r"are $82\%$ of the positives, so the pooled row of "
               r"Table~\ref{tab:archiveii} is largely those two; this table is what to read "
               r"when asking whether a method generalises across families. A dash marks a "
               r"tool with no prediction on that family. \emph{Covariation avail.} is the "
               r"fraction of the family's targets for which Infernal found a family with at "
               r"least five seed sequences; on PseudoBase the same figure is $26\%$ overall, "
               r"so ArchiveII is the regime where alignments are plentiful and covariation is "
               r"correspondingly favoured.}")
    out.append(r"\label{tab:byfamily}")
    out.append(r"\centering")
    out.append(r"\scriptsize")
    out.append(r"\setlength{\tabcolsep}{4pt}")
    out.append(r"\begin{tabular}{l" + "c" * len(fams) + "}")
    out.append(r"\toprule")
    out.append("Method & " + " & ".join(fams) + r" \\")
    out.append("$n$ & " + " & ".join(str(len(lens[f])) for f in fams) + r" \\")
    out.append("median nt & " + " & ".join(str(int(st.median(lens[f]))) for f in fams) + r" \\")
    out.append("covariation avail. & " + " & ".join(
        (f"{100*len(ids_mi[f])/len(ids_all[f]):.0f}\\%" if ids_all[f] else "---")
        for f in fams) + r" \\")
    out.append(r"\midrule")
    prev_italic = True
    for label, tool, italic in ROWS:
        if prev_italic and not italic:
            out.append(r"\midrule")
        prev_italic = italic
        cells = []
        for f in fams:
            a, b = f1.get((tool, f)), pk.get((tool, f))
            if not a:
                cells.append("---"); continue
            va = sum(a) / len(a)
            vb = (sum(b) / len(b)) if b else None
            cell = f"{va:.3f}" + (f" / {vb:.3f}" if vb is not None else " / ---")
            cells.append(rf"\textit{{{cell}}}" if italic else cell)
        out.append(f"{label} & " + " & ".join(cells) + r" \\")
    out.append(r"\bottomrule")
    out.append(r"\end{tabular}")
    out.append(r"\end{table}")
    print("\n".join(out))

    # coverage report to stderr so it never lands in the .tex
    print("\n[family_table] targets scored per tool:", file=sys.stderr)
    for _, tool, _ in ROWS:
        got = {f: len(f1.get((tool, f), [])) for f in fams}
        tot = sum(got.values())
        flag = "" if tot >= 0.95 * sum(len(lens[f]) for f in fams) else "   <-- INCOMPLETE"
        print(f"   {tool:16s} {tot:5d}  {got}{flag}", file=sys.stderr)


if __name__ == "__main__":
    main()
