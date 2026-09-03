# PROVENANCE — how a number is traced to the code that made it, and where that trace stops

This file exists because the honest answer is not "everything is pinned".

## The mechanism

`sim/progress.py` stamps two fields into every campaign's `progress.json`:

- `algo_version` — the semantic version of the simulation code, from
  `sim/_algo_version.py`;
- `code_sha` — a short SHA over all of `sim/*.py` and `sim/cpp/*.cpp`.

A campaign therefore knows which code produced it, and a changed algorithm
changes the fingerprint. Frozen snapshots live in `sim/versions/vX.Y.Z/` and are
never overwritten, so a fingerprint can be resolved to actual source.

## The limitation, enumerated

The stamping mechanism was added on 2026-08-05 at 09:47:23Z. Campaigns started
before that moment predate it and carry no fingerprint.

Of the 31 campaign directories in `runs/`, **5 carry a recorded
fingerprint** and **26 do not**. Both sets are listed below in full. The
counts appear here and nowhere else in this package; no sentence elsewhere
summarises these tables.

### Campaigns that carry a fingerprint (5)

| campaign | `algo_version` | `code_sha` |
|---|---|---|
| `runs/cost-equivalence-20260805T143329Z` | 1.3.2 | `457dcce1119a` |
| `runs/frozen-exit-20260807T064845Z` | 1.3.2 | `2df95f980738` |
| `runs/frozen-exit-20260807T065358Z` | 1.3.2 | `98144881e333` |
| `runs/frozen-exit-20260807T071638Z` | 1.3.2 | `75a8d3e2f30d` |
| `runs/nm-bound-20260805T114838Z` | 1.1.1 | `6b5cf0db492c` |

### Campaigns that carry no fingerprint (26)

Each is attributable to snapshot `sim/versions/v1.0.0` (`code_sha
e5614e8a4638`) **by inference from its start timestamp, not by a recorded
fingerprint.**

| campaign |
|---|
| `runs/exit-stability-20260805T080943Z` |
| `runs/heston-ads-lattice-20260805T061703Z` |
| `runs/heston-ads-lattice-20260805T065744Z` |
| `runs/heston-ads-lattice-20260805T070506Z` |
| `runs/heston-lattice-20260804T143759Z` |
| `runs/heston-lattice-20260804T143916Z` |
| `runs/heston-lattice-20260804T145350Z` |
| `runs/heston-lattice-20260804T145710Z` |
| `runs/heston-lattice-order-20260805T081819Z` |
| `runs/heston-lattice-order-20260805T082718Z` |
| `runs/mc-vs-tree-sweep-20260804T132550Z` |
| `runs/rheston-american-anchor-20260805T055638Z` |
| `runs/rheston-tables-20260805T055048Z` |
| `runs/rough-bergomi-ladder-20260805T084706Z` |
| `runs/rough-bergomi-ladder-20260805T085623Z` |
| `runs/rough-heston-20260804T141639Z` |
| `runs/rough-heston-20260804T142039Z` |
| `runs/route-aprime-20260804T114010Z` |
| `runs/route-b-20260804T123921Z` |
| `runs/route-b-20260804T124409Z` |
| `runs/route-b-20260804T125245Z` |
| `runs/route-b-lattice-20260805T071120Z` |
| `runs/route-b-scaling-20260804T130612Z` |
| `runs/routeb-compare-20260805T092052Z` |
| `runs/routeb-compare-20260805T094411Z` |
| `runs/validation-20260804T111030Z` |

That snapshot was verified byte-identical to the live `sim/` tree at the import
commit for four files, named here rather than described: `run_validation.py`,
`run_route_b.py`, `mc_reference.py`, `route_b_lattice.py`. The other files in the
snapshot were not individually compared.

What this means for a reader: for the campaigns in the second table you can
establish which code *almost certainly* ran, from the timestamps and the
snapshot, but you cannot verify it from the campaign's own record. Treat a
re-run of those as a re-measurement, not a bit-level replication.

## How to trace one number

1. Find the campaign directory in `runs/` whose name matches the family
   (`route-b-*`, `heston-ads-lattice-*`, …).
2. Read its `progress.json` for the configuration and, where present, the
   fingerprint.
3. Read its `RESULTS.md` for the tables as produced.
4. `log.txt` holds the run's own narration, including anything that went wrong.

`RUNBOOK.md` maps campaign families to the runner that produces them.

## What is not here

The manuscript, the proofs, the project's findings register and its review
record are not part of this package. A number's *interpretation* lived in that
register; what you have here is the number and the code.
