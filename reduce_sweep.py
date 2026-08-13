#!/usr/bin/env python3
"""Collapse the full sweep to one row per run.

Reads all_results_full.csv (28 columns, one row per timestep) and writes
run_summary.csv (one row per run). Everything the phase decomposition needs:
pre-shift baseline, the trough, the recovery, and the x statistics that test
the band-width and relocation predictions.

Usage:  python3 reduce_sweep.py
"""
import csv, sys
from collections import defaultdict

IN, OUT = "all_results_combined.csv", "run_summary.csv"

COLS = ("slope disp mutsd pert_treat patch_sd ac draw batch pert_value pert_name rep t n "
        "u_mean u_sd u_skew u_kurt mal_mean mal_sd mal_skew mal_kurt "
        "nn_mean nn_sd nn_skew nn_kurt x_mean x_sd x_skew x_kurt").split()
IDX = {c: i for i, c in enumerate(COLS)}
# batch must be part of the key: the two batches reuse replicate numbers,
# so (rep) alone would silently merge pairs of independent runs.
KEY = ["slope", "disp", "mutsd", "pert_treat", "patch_sd", "ac", "draw", "batch", "rep"]


def f(v):
    try:
        return float(v)
    except (ValueError, TypeError):
        return float("nan")


def median(xs):
    xs = sorted(x for x in xs if x == x)
    if not xs:
        return float("nan")
    m = len(xs) // 2
    return xs[m] if len(xs) % 2 else 0.5 * (xs[m - 1] + xs[m])


runs = defaultdict(list)
with open(IN) as fh:
    r = csv.reader(fh)
    next(r)
    for row in r:
        if len(row) != 29:
            continue
        runs[tuple(row[IDX[k]] for k in KEY)].append(row)

print(f"{len(runs):,} runs read", file=sys.stderr)

out = open(OUT, "w", newline="")
w = csv.writer(out)
w.writerow(KEY + [
    "base_n", "base_usd", "base_xsd", "base_xmean",
    "trough_n", "trough_t", "trough_usd", "trough_xsd", "trough_xmean",
    "peak_n", "peak_t", "final_n", "final_t",
    "final_usd", "final_xsd", "final_xmean", "final_mal",
    "n_steps", "survived"])

for k, rows in runs.items():
    rows.sort(key=lambda z: f(z[IDX["t"]]))
    pre = [z for z in rows if f(z[IDX["pert_value"]]) == 0]
    post = [z for z in rows if f(z[IDX["pert_value"]]) > 0]
    if not pre or not post:
        continue
    g = lambda z, c: f(z[IDX[c]])
    base_n = median([g(z, "n") for z in pre])
    tro = min(post, key=lambda z: g(z, "n"))
    after = [z for z in post if g(z, "t") >= g(tro, "t")]
    pk = max(after, key=lambda z: g(z, "n"))
    fin = after[-1]
    w.writerow([*k,
        base_n, median([g(z, "u_sd") for z in pre]),
        median([g(z, "x_sd") for z in pre]), median([g(z, "x_mean") for z in pre]),
        g(tro, "n"), g(tro, "t"), g(tro, "u_sd"), g(tro, "x_sd"), g(tro, "x_mean"),
        g(pk, "n"), g(pk, "t"), g(fin, "n"), g(fin, "t"),
        g(fin, "u_sd"), g(fin, "x_sd"), g(fin, "x_mean"), g(fin, "mal_mean"),
        len(rows), int(g(fin, "n") > 10)])

out.close()
print(f"wrote {OUT}", file=sys.stderr)
