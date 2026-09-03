# RUNBOOK — one section per campaign

Every command below is run from the root of this package, with the environment
of `INSTALL.md` active. No runner takes arguments except where stated. Each
writes a fresh `runs/<family>-<UTC>/` directory holding `progress.json`,
`log.txt` and `RESULTS.md`; the directories already in `runs/` are the
author's originals and are never overwritten by a new run.

Wall times are the **measured** `elapsed_s` field of each original run, on the
machine described in `INSTALL.md`. They are records, not estimates.

## The campaign register

31 directories, 27 completed and 4 not. Both sets appear in
full below. These counts appear here and nowhere else in this file.

| campaign directory | runner | status | wall time | started (UTC) | phases |
|---|---|---|---|---|---|
| `runs/cost-equivalence-20260805T143329Z` | `sim/run_cost_equivalence.py` | done | 5m18s | 2026-08-05T14:33:29Z | 6 |
| `runs/exit-stability-20260805T080943Z` | `sim/run_exit_stability.py` | done | 3m24s | 2026-08-05T08:09:43Z | 6 |
| `runs/frozen-exit-20260807T064845Z` | `sim/run_frozen_exit.py` | done | 3m46s | 2026-08-07T06:48:45Z | 6 |
| `runs/frozen-exit-20260807T065358Z` | `sim/run_frozen_exit.py` | done | 3m44s | 2026-08-07T06:53:58Z | 6 |
| `runs/frozen-exit-20260807T071638Z` | `sim/run_frozen_exit.py` | done | 5m14s | 2026-08-07T07:16:38Z | 6 |
| `runs/heston-ads-lattice-20260805T061703Z` | `sim/run_heston_ads_tables.py` | done | 29m37s | 2026-08-05T06:17:03Z | 4 |
| `runs/heston-ads-lattice-20260805T065744Z` | `sim/run_heston_ads_tables.py` | running | 1m03s | 2026-08-05T06:57:44Z | 4 |
| `runs/heston-ads-lattice-20260805T070506Z` | `sim/run_heston_ads_tables.py` | done | 57m12s | 2026-08-05T07:05:06Z | 4 |
| `runs/heston-lattice-20260804T143759Z` | `sim/run_heston_lattice.py` | done | 31s | 2026-08-04T14:37:59Z | 4 |
| `runs/heston-lattice-20260804T143916Z` | `sim/run_heston_lattice.py` | done | 30s | 2026-08-04T14:39:16Z | 4 |
| `runs/heston-lattice-20260804T145350Z` | `sim/run_heston_lattice.py` | failed | 2m41s | 2026-08-04T14:53:50Z | 5 |
| `runs/heston-lattice-20260804T145710Z` | `sim/run_heston_lattice.py` | done | 2m38s | 2026-08-04T14:57:10Z | 5 |
| `runs/heston-lattice-order-20260805T081819Z` | `sim/run_heston_lattice_order.py` | done | 6m21s | 2026-08-05T08:18:19Z | 2 |
| `runs/heston-lattice-order-20260805T082718Z` | `sim/run_heston_lattice_order.py` | done | 5m39s | 2026-08-05T08:27:18Z | 2 |
| `runs/mc-vs-tree-sweep-20260804T132550Z` | `sim/run_mc_vs_tree_sweep.py` | done | 2m18s | 2026-08-04T13:25:50Z | 5 |
| `runs/nm-bound-20260805T114838Z` | `sim/run_nm_bound.py` | done | 11m35s | 2026-08-05T11:48:38Z | 5 |
| `runs/rheston-american-anchor-20260805T055638Z` | `sim/run_rheston_american_anchor.py` | done | 40s | 2026-08-05T05:56:38Z | 2 |
| `runs/rheston-tables-20260805T055048Z` | `sim/run_rheston_tables.py` | done | 5m39s | 2026-08-05T05:50:48Z | 5 |
| `runs/rough-bergomi-ladder-20260805T084706Z` | `sim/run_rough_bergomi_ladder.py` | done | 6m27s | 2026-08-05T08:47:06Z | 4 |
| `runs/rough-bergomi-ladder-20260805T085623Z` | `sim/run_rough_bergomi_ladder.py` | done | 5m52s | 2026-08-05T08:56:23Z | 4 |
| `runs/rough-heston-20260804T141639Z` | `sim/run_rough_heston.py` | failed | 5s | 2026-08-04T14:16:39Z | 5 |
| `runs/rough-heston-20260804T142039Z` | `sim/run_rough_heston.py` | done | 46s | 2026-08-04T14:20:39Z | 5 |
| `runs/route-aprime-20260804T114010Z` | `sim/run_route_aprime.py` | done | 13m04s | 2026-08-04T11:40:10Z | 4 |
| `runs/route-b-20260804T123921Z` | `sim/run_route_b.py` | done | 42s | 2026-08-04T12:39:21Z | 6 |
| `runs/route-b-20260804T124409Z` | `sim/run_route_b.py` | done | 40s | 2026-08-04T12:44:09Z | 7 |
| `runs/route-b-20260804T125245Z` | `sim/run_route_b.py` | done | 45s | 2026-08-04T12:52:45Z | 7 |
| `runs/route-b-lattice-20260805T071120Z` | `sim/run_route_b_lattice.py` | done | 3m43s | 2026-08-05T07:11:20Z | 6 |
| `runs/route-b-scaling-20260804T130612Z` | `sim/run_route_b_scaling.py` | done | 1m32s | 2026-08-04T13:06:12Z | 4 |
| `runs/routeb-compare-20260805T092052Z` | `sim/run_routeb_compare.py` | running | 4m05s | 2026-08-05T09:20:52Z | 4 |
| `runs/routeb-compare-20260805T094411Z` | `sim/run_routeb_compare.py` | done | 19m10s | 2026-08-05T09:44:11Z | 4 |
| `runs/validation-20260804T111030Z` | `sim/run_validation.py` | done | 3s | 2026-08-04T11:10:30Z | 4 |

Sum of the wall times of the completed campaigns: **3.28 h** (11823 s). The two longest are `runs/heston-ads-lattice-20260805T070506Z` at 57m12s and `runs/heston-ads-lattice-20260805T061703Z` at 29m37s.
Re-running everything is an afternoon's work on one machine, not a coffee
break, and the `make quick` path in `README.md` exists so you can check the
environment before committing to it.

### The campaigns that did not complete

Kept because they are part of the record, not because they are usable. Do not
read a table out of these four.

| campaign directory | status | phase reached | what it means |
|---|---|---|---|
| `runs/heston-ads-lattice-20260805T065744Z` | running | `1. European put and call: our lattice vs ADS's tree vs analytic` | interrupted, so `progress.json` still reads `running`; superseded by the next run in the same family |
| `runs/heston-lattice-20260804T145350Z` | failed | `the contrast with the rough regime` | the runner raised and stopped; `log.txt` holds the traceback |
| `runs/rough-heston-20260804T141639Z` | failed | `the reference's own accuracy` | the runner raised and stopped; `log.txt` holds the traceback |
| `runs/routeb-compare-20260805T092052Z` | running | `1. H = 0.05: prices` | interrupted, so `progress.json` still reads `running`; superseded by the next run in the same family |

## The runners

### `sim/run_rheston_american_anchor.py`

```
run_rheston_american_anchor — how good is our American estimator, really?

The American table of the rough Heston document has no reference column, because
for the rough model none exists.  That leaves the reader with no way to judge the
Longstaff--Schwartz numbers.  This run supplies the missing calibration by
running the SAME estimator at H = 0.5, where the model is classical Heston and
```

```bash
python3 sim/run_rheston_american_anchor.py
```

Campaign directories it produced (1):

| directory | status | wall time |
|---|---|---|
| `runs/rheston-american-anchor-20260805T055638Z` | done | 40s |

Files in the most recent of them: `RESULTS.md`, `log.txt`, `progress.json`

### `sim/run_rough_bergomi_ladder.py`

```
run_rough_bergomi_ladder — the paper's lattice IN THE ROUGH REGIME.

The H = 0.5 tables validate the construction against an exact price, but H = 0.5
is not rough, and it is the rough case that the paper is about.  This run prices
the genuinely rough model with the paper's own lattice and puts it next to a
trustworthy Monte-Carlo, across H, with the error and the Monte-Carlo band on the
```

```bash
python3 sim/run_rough_bergomi_ladder.py
```

Campaign directories it produced (2):

| directory | status | wall time |
|---|---|---|
| `runs/rough-bergomi-ladder-20260805T084706Z` | done | 6m27s |
| `runs/rough-bergomi-ladder-20260805T085623Z` | done | 5m52s |

Files in the most recent of them: `RESULTS.md`, `log.txt`, `progress.json`

### `sim/run_heston_lattice_order.py`

```
run_heston_lattice_order — the lattice's convergence order, measured against a
reference that is precise enough to measure it.

Why this is a separate run.  The main lattice table compares against the
analytical Heston column ADS PUBLISH, which is the right reference for comparing
two lattices on equal terms -- but it is printed to four decimals, so an error
```

```bash
python3 sim/run_heston_lattice_order.py
```

Campaign directories it produced (2):

| directory | status | wall time |
|---|---|---|
| `runs/heston-lattice-order-20260805T081819Z` | done | 6m21s |
| `runs/heston-lattice-order-20260805T082718Z` | done | 5m39s |

Files in the most recent of them: `RESULTS.md`, `log.txt`, `progress.json`

### `sim/run_heston_ads_tables.py`

```
run_heston_ads_tables — THE PAPER'S OWN LATTICE on the ADS2014 grid.

The rough Heston document has Fourier and Monte-Carlo columns but no lattice
column, which is the wrong shape for a paper about a lattice.  This run supplies
it, and the place it can be supplied is `H = 0.5`.
```

```bash
python3 sim/run_heston_ads_tables.py
```

Campaign directories it produced (3):

| directory | status | wall time |
|---|---|---|
| `runs/heston-ads-lattice-20260805T061703Z` | done | 29m37s |
| `runs/heston-ads-lattice-20260805T065744Z` | running | 1m03s |
| `runs/heston-ads-lattice-20260805T070506Z` | done | 57m12s |

Files in the most recent of them: `RESULTS.md`, `log.txt`, `progress.json`

### `sim/run_route_b_lattice.py`

```
run_route_b_lattice — items (B2)/(B3): Route B built as a LATTICE, and the
telescope of Section 9.5 measured term by term.

Everything Route B had until now was kernel approximation: ||K - K^m||, the
covariance surface, m*(n), a cost exponent.  No lattice was ever constructed.
Two things follow from building one.
```

```bash
python3 sim/run_route_b_lattice.py
```

Campaign directories it produced (1):

| directory | status | wall time |
|---|---|---|
| `runs/route-b-lattice-20260805T071120Z` | done | 3m43s |

Files in the most recent of them: `RESULTS.md`, `log.txt`, `progress.json`

### `sim/run_route_b_scaling.py`

```
run_route_b_scaling — preliminary 2: does the number of factors m stay bounded
as the time grid is refined?

The first Route B run (runs/route-b-*) answered "how many factors" at n = 64 and
n = 256 and got small numbers.  But the threshold the lift has to beat is the
scheme's OWN discretisation error, and that shrinks as n grows, so m must grow
```

```bash
python3 sim/run_route_b_scaling.py
```

Campaign directories it produced (1):

| directory | status | wall time |
|---|---|---|
| `runs/route-b-scaling-20260804T130612Z` | done | 1m32s |

Files in the most recent of them: `RESULTS.md`, `log.txt`, `progress.json`

### `sim/run_mc_vs_tree_sweep.py`

```
run_mc_vs_tree_sweep — the tree against a trustworthy Monte-Carlo, swept over the
model parameters, so that the error can be read as a function of them.

Model: rough Bergomi (the paper's Section "Rough Bergomi", Example 4.14), which
is what the code implements.  Rough Heston is a different model -- a fractional
CIR variance -- and would be a new implementation, not a parameter change.
```

```bash
python3 sim/run_mc_vs_tree_sweep.py
```

Campaign directories it produced (1):

| directory | status | wall time |
|---|---|---|
| `runs/mc-vs-tree-sweep-20260804T132550Z` | done | 2m18s |

Files in the most recent of them: `RESULTS.md`, `log.txt`, `progress.json`

### `sim/run_cost_equivalence.py`

```
run_cost_equivalence — how long a lattice price takes, and how many Monte-Carlo
paths (and how much time) it takes to reach the SAME accuracy.

Six phases.

  setup       the Monte-Carlo's one-off cost, timed cold (an nfine x nfine
```

```bash
python3 sim/run_cost_equivalence.py
```

Campaign directories it produced (1):

| directory | status | wall time |
|---|---|---|
| `runs/cost-equivalence-20260805T143329Z` | done | 5m18s |

Files in the most recent of them: `RESULTS.md`, `cost_equivalence.json`, `log.txt`, `progress.json`

### `sim/run_routeb_compare.py`

```
run_routeb_compare — Route B vs the one-step lattice, priced, across H.

Part V established the obstruction: the one-step recombining rough lattice has a
driver whose variance is `2H n^{1-2H}` times the truth, so refining the grid
drives the price away from the reference -- unboundedly.  This run implements
Route B (the Markovian lift) as an actual lattice and shows what it changes.
```

```bash
python3 sim/run_routeb_compare.py
```

Campaign directories it produced (2):

| directory | status | wall time |
|---|---|---|
| `runs/routeb-compare-20260805T092052Z` | running | 4m05s |
| `runs/routeb-compare-20260805T094411Z` | done | 19m10s |

Files in the most recent of them: `RESULTS.md`, `log.txt`, `progress.json`

### `sim/run_rheston_tables.py`

```
run_rheston_tables — rough Heston priced on the Akyildirim--Dolinsky--Soner grid.

This run is deliberately KEPT ASIDE from the paper.  The paper's model class needs
an autonomous driver, and rough Heston is outside it (Section 10.8), so nothing
here is a lattice result.  What it is: a reference table set for the rough Heston
model itself, laid out like the numerical section of
```

```bash
python3 sim/run_rheston_tables.py
```

Campaign directories it produced (1):

| directory | status | wall time |
|---|---|---|
| `runs/rheston-tables-20260805T055048Z` | done | 5m39s |

Files in the most recent of them: `RESULTS.md`, `log.txt`, `progress.json`

### `sim/run_exit_stability.py`

```
run_exit_stability — Lemma 9.6 (freezing the clock): the step the paper left as
an argument, now proved, and here measured.

The gap.  To make the coupled chain Markov, the conditional law of
(Delta W_j, Delta lambda_j) must be a function of the LATTICE state.  Freezing
sigma_y and mu_y at the node perturbs the driver's coefficients, and one then has
```

```bash
python3 sim/run_exit_stability.py
```

Campaign directories it produced (1):

| directory | status | wall time |
|---|---|---|
| `runs/exit-stability-20260805T080943Z` | done | 3m24s |

Files in the most recent of them: `RESULTS.md`, `log.txt`, `progress.json`

### `sim/run_heston_lattice.py`

```
run_heston_lattice — the lattice against a reference containing no Monte-Carlo.

This closes, as far as it can be closed, the question the independent reference
was built for: does the LATTICE agree with a non-Monte-Carlo price?

It cannot be asked in the rough regime, and the reason is structural rather than
```

```bash
python3 sim/run_heston_lattice.py
```

Campaign directories it produced (4):

| directory | status | wall time |
|---|---|---|
| `runs/heston-lattice-20260804T143759Z` | done | 31s |
| `runs/heston-lattice-20260804T143916Z` | done | 30s |
| `runs/heston-lattice-20260804T145350Z` | failed | 2m41s |
| `runs/heston-lattice-20260804T145710Z` | done | 2m38s |

Files in the most recent of them: `RESULTS.md`, `log.txt`, `progress.json`

### `sim/run_route_aprime.py`

```
run_route_aprime — Level 4 redone with Route A' (paper §9.2), against the
four-point scheme of §4.

Four phases.

  1  (V8)   eta = 0: Route A' against the exact CRR reference, with the
```

```bash
python3 sim/run_route_aprime.py
```

Campaign directories it produced (1):

| directory | status | wall time |
|---|---|---|
| `runs/route-aprime-20260804T114010Z` | done | 13m04s |

Files in the most recent of them: `RESULTS.md`, `log.txt`, `progress.json`

### `sim/run_rough_heston.py`

```
run_rough_heston — item (N6): the independent reference, and what it costs.

Every price comparison in this project so far has set our lattice against our own
Monte-Carlo, so a shared misreading of the model would cancel and stay invisible.
Rough Heston breaks that: its characteristic function solves a fractional Riccati
equation and the price follows by Fourier inversion, with no Monte-Carlo anywhere
```

```bash
python3 sim/run_rough_heston.py
```

Campaign directories it produced (2):

| directory | status | wall time |
|---|---|---|
| `runs/rough-heston-20260804T141639Z` | failed | 5s |
| `runs/rough-heston-20260804T142039Z` | done | 46s |

Files in the most recent of them: `RESULTS.md`, `log.txt`, `progress.json`

### `sim/run_frozen_exit.py`

```
run_frozen_exit -- driver for the Route F falsifiers FF1-FF4, FF6, FF7.

    python3 sim/run_frozen_exit.py            # full run (~7 min)
    python3 sim/run_frozen_exit.py --quick    # smoke run (~1 min)

Writes runs/frozen-exit-<UTC>/{progress.json, log.txt, RESULTS.md} through
```

```bash
python3 sim/run_frozen_exit.py [--quick]
```

`--quick` is a smoke run, about 1 min against about 7 min for the full run;
the runner's own docstring states both.

Campaign directories it produced (3):

| directory | status | wall time |
|---|---|---|
| `runs/frozen-exit-20260807T064845Z` | done | 3m46s |
| `runs/frozen-exit-20260807T065358Z` | done | 3m44s |
| `runs/frozen-exit-20260807T071638Z` | done | 5m14s |

Files in the most recent of them: `RESULTS.md`, `log.txt`, `progress.json`

### `sim/run_validation.py`

```
run_validation — the validation stack of paper §10.3, on the simplest example
that still exercises the whole pipeline.

EXAMPLE.  Fractional Brownian driver, constant price volatility:

    Y = W  (sigma_y = 1, mu_y = 0),      K(u) = u^h,  h = H - 1/2
```

```bash
python3 sim/run_validation.py
```

Campaign directories it produced (1):

| directory | status | wall time |
|---|---|---|
| `runs/validation-20260804T111030Z` | done | 3s |

Files in the most recent of them: `RESULTS.md`, `log.txt`, `progress.json`

### `sim/run_nm_bound.py`

```
run_nm_bound — the n<->m link, measured: how many lift factors a grid of n steps
needs, across H, across tolerances, and across the three constructions.

The question.  Route B replaces the single divergent counter of the one-step
lattice by m Ornstein--Uhlenbeck factors.  Part VI showed that this converts an
UNBOUNDED error into a bounded one, but left the obvious practical question
```

```bash
python3 sim/run_nm_bound.py
```

Campaign directories it produced (1):

| directory | status | wall time |
|---|---|---|
| `runs/nm-bound-20260805T114838Z` | done | 11m35s |

Files in the most recent of them: `RESULTS.md`, `log.txt`, `nm_law.json`, `progress.json`

### `sim/run_route_b.py`

```
run_route_b — item (N4)/(B1): quantify the Markovian lift BEFORE proving it.

The question that decides whether Route B is usable: how many factors m does the
lift need, and what does that cost?

Seven phases.
```

```bash
python3 sim/run_route_b.py
```

Campaign directories it produced (3):

| directory | status | wall time |
|---|---|---|
| `runs/route-b-20260804T123921Z` | done | 42s |
| `runs/route-b-20260804T124409Z` | done | 40s |
| `runs/route-b-20260804T125245Z` | done | 45s |

Files in the most recent of them: `RESULTS.md`, `log.txt`, `progress.json`

