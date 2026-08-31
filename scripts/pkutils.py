#!/usr/bin/env python3
"""Structure helpers shared by the release scripts.

These three functions are lifted verbatim from the project's covariation_select.py so that this
supplementary bundle stands alone: the full generation pipeline depends on a closed-form energy
module that is not part of this release, and importing it here would fail.

  dbp              dot-bracket string -> set of (i, j) pairs, all bracket types plus letter pairs
  native_pk_pairs  the pairs of a NATIVE structure carrying a minority bracket, i.e. the arm a
                   nested folder cannot represent alongside the primary one. PKrec is recovery
                   of exactly this set, and it is recomputed from pred_db for every method --
                   never read from the pk_sty column, which two different producers write under
                   two different conventions.
  crossing_pairs   the pairs of a PREDICTED structure that cross some other pair in it
"""

_OPEN = {"(": ")", "[": "]", "{": "}", "<": ">"}; _CLOSE = {v: k for k, v in _OPEN.items()}


def dbp(db):
    st = {o: [] for o in _OPEN}; pr = set()
    for i, c in enumerate(db):
        if c in _OPEN: st[c].append(i)
        elif c in _CLOSE and st[_CLOSE[c]]: pr.add((st[_CLOSE[c]].pop(), i))
    return pr
_BRACKETS = [("(", ")"), ("[", "]"), ("{", "}"), ("<", ">"),
             ("A", "a"), ("B", "b"), ("C", "c"), ("D", "d")]


def native_pk_pairs(db):
    """The genuine pseudoknot pairs of a NATIVE dot-bracket: those annotated with a secondary
    bracket type ([, {, <) or with letter pairs (Aa, Bb, ...). '(' is the primary/nested
    annotation. Recovery is then |predicted_pairs & native_pk_pairs| / |native_pk_pairs|,
    which is what 'did the method find the pseudoknot' actually means."""
    st = {o: [] for o in _OPEN}; lets = {}; sec = set()
    for i, c in enumerate(db):
        if c in _OPEN: st[c].append(i)
        elif c in _CLOSE and st[_CLOSE[c]]:
            p = (st[_CLOSE[c]].pop(), i)
            if _CLOSE[c] != "(": sec.add(p)
        elif c.isalpha():
            u = c.upper()
            if c.isupper(): lets.setdefault(u, []).append(i)
            elif lets.get(u): sec.add((lets[u].pop(), i))
    return sec

# ---------------- Infernal calls ----------------


def cross(pr, k): (i, j), (a, b) = pr, k; return i < a < j < b or a < i < b < j


def crossing_pairs(pairset):
    """Subset of pairs that cross at least one other pair. This is the SYMMETRIC rule, identical
    to split_nc() in bench_pk_tools.py, so numbers are comparable to the existing benchmark
    tables -- but note it is inflated ~5.2x on this dataset (1334 vs 255 pairs over the 132
    natives), because a single tertiary crossing pair drags in every nested stem pair it crosses.
    For tRNA it reports ~10 'crossing' pairs where only 1 is a genuine pseudoknot pair.
    Use native_pk_pairs() below for the metric that answers 'did we find the pseudoknot'."""
    P = list(pairset); out = set()
    for a in range(len(P)):
        for b in range(a + 1, len(P)):
            if cross(P[a], P[b]): out.add(P[a]); out.add(P[b])
    return out


# ---------------------------------------------------------------------------------------------
# PATH RESOLUTION. In the working repository the CSVs, TSVs and scripts sit in three different
# directories; in this bundle they are data/, results/ and scripts/. Rather than hard-code
# either layout, look for a file by basename in the places it could reasonably be.
import os as _os

_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
_SEARCH = [_os.path.join(_ROOT, "data"), _os.path.join(_ROOT, "results"), _ROOT, "."]


def find(name):
    """Absolute path to `name`, searched across the bundle's directories. Returns the bare name
    unchanged if nothing matches, so a caller's own error message is the one the user sees."""
    if _os.path.isabs(name) and _os.path.exists(name):
        return name
    for d in _SEARCH:
        p = _os.path.join(d, _os.path.basename(name))
        if _os.path.exists(p):
            return p
    return name
