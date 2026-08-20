# New simulation dataset

I've rerun the simulations from scratch. This is a new dataset — not an extension of
the earlier runs, and not directly comparable to them.

## Design

**30 runs per treatment combination**, across 630 combinations — 18,796 runs recorded of 18,900 launched (0.6% of populations went extinct before the perturbation and have no post-shift record).

Replication has two levels, because two different things can vary:

- **Landscape.** A setting like "autocorrelation 4, patch SD 3" doesn't specify one
  surface — it specifies a recipe. Run that recipe with different random seeds and you
  get different terrain with the same statistical character. I generated **3 landscapes**
  per condition. This is the `draw` column.
- **Population dynamics.** On each landscape, **10 independent simulation runs**, each
  with its own seed governing founder placement, births and deaths. The `batch` column
  is constant (1) in this version.

3 landscapes x 10 runs = **30 runs per combination**.

The two levels answer different questions. Three landscapes tell you whether an effect
is a property of the autocorrelation level or of one particular surface — in the earlier
design each treatment used a single landscape file, so those two were indistinguishable.
Eighteen runs per landscape tell you whether an effect is the treatment or demographic
chance.

At patch SD 0 the surface is flat, so all three draws are the same terrain. Those 30 runs
are still independent (different simulation seeds) — there is simply no landscape
variation to sample.

| Factor | Levels |
|---|---|
| Gradient slope | 0.8, 1.0, 1.2 |
| Dispersal SD | 1.5, 3, 6 |
| Mutation SD | 0 (standing variation only), 0.75 |
| Perturbation magnitude | 9, 10, 11, 12, 13 |
| Landscape heterogeneity (patch SD) | 0, 1, 2 |
| Spatial autocorrelation | 0, 2, 4 |

Autocorrelation is undefined when patch SD is 0 (a flat surface has no spatial structure
to correlate), so the landscape factor has **7 conditions**: flat, plus patch SD 1 and 2
each crossed with autocorrelation 0, 2 and 4.

Measured lag-1 spatial correlation of the generated landscapes: **0.00 / 0.74 / 0.996**
for autocorrelation 0 / 2 / 4. Worth quoting these rather than the exponent — they're
verifiable directly from the deposited landscape files.

## Fixed parameters

700 founders, initialised perfectly adapted, on a 100 x 100 landscape. Perturbation at
t = 250; runs continue to t = 900 (650 steps post-shift); state recorded every 10 steps.
Birth rate 0.54 constant. Carrying capacity `K = 10 * exp(-(u - optimum)^2 / 3^2)`.
Death rate = local competition / K. Clonal reproduction, natal dispersal only.

## Differences from the earlier runs

- **Dispersal** now has three levels (1.5, 3, 6) rather than two. Two levels can't
  detect a non-monotonic response, and there is one: persistence peaks at 3.
- **Mutation SD** is 0.75 rather than 1.5. A pilot showed 0.75 captures the full effect
  at half the mutational magnitude.
- **Patch SD** is 0, 1, 2. A separate viability map (supplement) shows uncorrelated
  terrain becomes uninhabitable between SD 2 and 3, so 2 is the highest level at which
  every autocorrelation treatment supports a population.
- **Three independent landscape draws** per condition rather than one, so landscape
  identity is no longer confounded with treatment level.
- **Perturbation at t = 250.** Most conditions equilibrate by t = 100; the slowest
  (patch SD 2, autocorrelation 0) needs until roughly t = 400, so t = 250 is a
  compromise between a settled baseline and compute. The post-shift window is 650 steps.

## Output format

One row per timestep per run, rather than one row per individual. Columns:

```
slope, disp, mutsd, pert_treat, patch_sd, ac, draw, batch,
pert_value, pert_name, rep, t, n,
u_mean, u_sd, u_skew, u_kurt,
mal_mean, mal_sd, mal_skew, mal_kurt,
nn_mean, nn_sd, nn_skew, nn_kurt,
x_mean, x_sd, x_skew, x_kurt
```

Trait, maladaptation, nearest-neighbour distance and x-position are summarised in the
simulator as mean, SD, skewness and excess kurtosis. Maladaptation is
`|slope * x + patch - u|`.

A run is identified by **rep and batch together** — the two batches reuse replicate
numbers 0–9 and 0–7 respectively, with independent seeds.

Skewness is `m3 / m2^1.5` and kurtosis is excess (`m4 / m2^2 - 3`), which differ slightly
from the `e1071` defaults. Don't mix these with values computed in R.

## Analysis

`analyse_sweep.R` does all of this. Base R, no packages, about a minute to run.

**Outcomes.** From each run's abundance trajectory: the pre-perturbation baseline
(median of the pre-shift record), the trough after the shift, and the peak after the
trough. From those:

| | |
|---|---|
| Resistance | trough / baseline — how far it fell |
| Recovery magnitude | (peak − trough) / baseline — how far it came back |
| Recovery relative to baseline | peak / baseline — whether it regained its starting size |
| Recovery time | peak time − trough time |
| Persistence | survived or not |

The four phase outcomes are fitted on surviving runs; persistence uses all runs.
Resistance and recovery correlate at only +0.29, so they are not mechanically linked
and an opposing pattern between them would be a real one rather than a definitional
artefact.

**From simulation output to `run_summary.csv`.** Two steps.

`combine_results.sh` concatenates the 3,780 per-cell files and parses the treatment
levels out of each filename, so `S0.8_D1.5_M0_P9_sd3_ac4_r2_b2.csv` contributes
slope 0.8, dispersal 1.5, mutation 0, severity 9, patch SD 3, autocorrelation 4,
draw 2, batch 2. The result is one row per timestep with the design attached.

`reduce_sweep.py` then collapses each run to a single row. It groups by the eight-column
key (six treatments plus draw, batch and rep), sorts by time, and pulls out the
pre-perturbation baseline as the *median* of all pre-shift timesteps — not a single
timepoint, because abundance fluctuates by a few percent between steps and a
single-timepoint baseline inherits all of that — then the minimum after the shift, the
maximum at or after that minimum, and the final timestep. Trait SD, x SD and x mean are
recorded at each of those four points. Runs with no post-shift record are dropped; that
is the 4.1%.

The phase metrics are computed in `analyse_sweep.R` rather than during the reduction, so
they stay visible in the analysis instead of being baked into the data.

**Models.** Linear models for the phase outcomes, logistic for persistence, all with
the same predictors: slope, dispersal, mutation, severity, patch SD and autocorrelation,
plus slope x severity and autocorrelation x severity.

Every predictor enters as a factor rather than a linear term, so nothing is assumed to
be monotonic — which matters, because dispersal is not: persistence peaks at the middle
level.

Standard errors are clustered by parameter cell, since runs sharing a cell are not
independent. Landscape terms (autocorrelation, patch SD) are fitted on patchy landscapes
only, because autocorrelation is undefined on a flat surface — including the flat runs
would compare three identical treatments. Effects are tested jointly across all levels
of a factor rather than coefficient by coefficient.

**Robustness.** Three checks on the main result: whether it holds in each of the three
landscape draws separately, whether it holds in each of the two batches separately, and
whether it survives dropping the thinnest treatment combinations.

**What is not here.** No time-to-extinction model. A Cox model was fitted and the
proportional-hazards assumption was rejected for four of the five predictors, so the
hazard ratios would be misleading. Persistence is modelled as a binary outcome instead.

## Files

| | |
|---|---|
| `NEW_SIMULATIONS.md` | this note |
| `run_summary_v2.csv` | one row per run: baseline, trough, peak, endpoint. What the analysis uses. |
| `all_results_v2.csv` | every timestep of every run, with treatment columns |
| `analyse_sweep.R` | analysis and figures, base R only, no packages |
| `combine_results.sh` | concatenates the per-cell output and parses treatment levels from the filenames |
| `reduce_sweep.py` | collapses `all_results_combined.csv` to one row per run |
| `models_table.csv` | joint p-values, each parameter against each outcome |
| `effects_table.csv` | means by treatment level for each outcome |
| `fig1` – `fig5` | phase signatures, comparison matrix, dispersal, slope x severity, draw robustness |
| `simulation/` | C++ source, `Makefile`, `run_design.sh`, `make_landscapes.py`, the sweep logs, and `landscapes/` with the 27 landscape files |

The per-cell raw output (3,780 CSVs) isn't here — all of it is in
`all_results_combined.csv`. Ask if you want it.

## One thing to know about the design

The highest heterogeneity level (patch SD 2) with zero autocorrelation equilibrates at
a much lower abundance than the other conditions — around 230 against ~820 on flat
ground — and takes until roughly t = 400 to settle fully, so it is perturbed slightly
above its long-run equilibrium. Every other condition equilibrates well before t = 250.
Pre-perturbation extinction is 0.6% overall.
