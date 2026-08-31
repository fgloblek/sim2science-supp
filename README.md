# Supplementary material

The data, the per-target predictions of every method, and the scripts that regenerate every
table in the paper. Each number can therefore be re-derived from the predictions.

## What is here

    data/     benchmark TSVs -- id, sequence, native dot-bracket
              archiveii_pk.tsv (1063 pseudoknotted) / archiveii_nopk.tsv (2881 free)
              pseudobase.tsv (391) with its sequence-identity clusters
              nested_negatives.tsv (489 curated pseudoknot-free controls)
    results/  one row per (target, tool): id, n, tool, ppv, sty, f1, pk_sty, pred_db
    scripts/  the table generators, plus pkutils.py (dot-bracket helpers)

## Running it

    cd scripts
    python pseudobase_tables.py ../results/final3_pos_bench.csv ../results/final3_neg_bench.csv
    python family_table.py ../results/archiveii_pk_bench.csv ../data/archiveii_pk.tsv \
                           ../results/ipknot_ext_pk.csv > family_table.tex
    python verify_numbers.py

Python 3 and NumPy only. `verify_numbers.py` recomputes every benchmark number in the paper
from `pred_db` and fails loudly on any mismatch. Point `PAPER_TEX` at the paper source to
additionally check that the text still quotes the values the data give; without it, those
literal-matching checks are skipped and the recomputations still run.

## What is NOT here

The generation pipeline. Building the candidate pools requires the closed-form pseudoknot free
energy of Appendix B, which is the subject of a separate manuscript in preparation and is not
included, so the pools cannot be regenerated from this release. Everything downstream of the
pools -- every selector, metric and table -- is fully reproducible from what is here.

## Two conventions worth knowing

**PKrec** is recovery of the native pairs carrying a *minority* bracket, and it is recomputed
here from `pred_db` for every method. Do not read the `pk_sty` column: rows written by the
external-tool harness use the geometric crossing set, which contains both arms of a pseudoknot
and inflates the score, while our own rows use the minority-bracket convention. The two are not
comparable, which is why nothing downstream trusts that column.

**Cluster-level means.** Nothing is a plain target average. PseudoBase averages within
sequence-identity cluster and then across the 288 clusters; ArchiveII uses the curated family as
the cluster. Bootstrap intervals resample that same unit.

## Provenance and licence

The TSVs are derived from PseudoBase (van Batenburg et al., 2001) and ArchiveII (Sloma and
Mathews, 2016); both are public, and we redistribute derived forms for reproducibility only,
under the terms of the original sources. The derived data files are offered under CC BY 4.0 and
the scripts under the MIT licence.
