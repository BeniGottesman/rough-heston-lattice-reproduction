# Recombining lattices for rough volatility — reproduction package

Code, campaign results and instructions for reproducing the numerical work of
*Recombining lattices for rough Heston*.

Everything here was run by the author on a single machine. Nothing in this
repository is required to read the paper; it exists so that a reader who wants
to re-derive a number can do so without asking anyone.

## What this reproduces, and what it does not

**It reproduces** every numerical campaign behind the paper: the two-dimensional
Heston lattice, the Markovian-lift construction and its lattice, the
rough-Bergomi route, the rough-Heston Fourier and Monte-Carlo pricers, the
n↔m link, and the cost-equivalence measurement. Each campaign's own output —
`progress.json`, `RESULTS.md`, `log.txt` — is included as it was produced, so a
re-run can be compared against the original rather than against a retyped table.

**It does not** contain the manuscript, the proofs, or the project's research
record. Those are not public.

**It does not claim bit-identical re-runs on a different machine.** The C++
kernels are built with `-march=native`, and the reproduction status of each
campaign is stated per campaign in `RUNBOOK.md` — measured, not assumed. Read
`PROVENANCE.md` before treating any number here as pinned to a specific commit:
one important limitation is recorded there.

## Layout

| path | what it is |
|---|---|
| `sim/*.py` | the live simulation code, 33 files |
| `sim/cpp/` | the two C++17 kernels and their own README |
| `sim/versions/v1.0.0 … v1.3.2` | frozen snapshots of the code, never overwritten |
| `verify/` | three symbolic verifications of the paper's algebra (sympy) |
| `runs/` | 31 campaign directories, exactly as produced |
| `INSTALL.md` | the environment, with the versions actually used |
| `RUNBOOK.md` | one section per campaign: the command, its cost, its output |
| `PROVENANCE.md` | how a number is traced back to the code that made it, and the limits of that trace |
| `MANIFEST.md` | every file in this package with its sha256 |

## Thirty seconds, to check the environment works

The symbolic verifications need only Python and sympy, and finish in seconds:

```bash
python3 verify/verify_claims.py
```

On the author's machine this prints `88/88 passed` in **22 s** (Python 3.9.6,
sympy 1.14.0, Apple M-series). It checks the algebraic content of the paper —
moment equations, non-negativity of the nine-point law, the rate exponents —
symbolically, as identities of rational functions. It is not a proof assistant:
the quantifier structure is supplied by hand. See the header of that file.

If that passes, the Python side of the environment is correct. `make verify`
runs all three, each measured on the author's machine on 2026-09-03:

| script | result | wall time |
|---|---|---|
| `verify/verify_claims.py` | `88/88 passed` | 22 s |
| `verify/verify_lotc_moments.py` | `18/18 passed` | 3 s |
| `verify/verify_mc_band.py` | `all checks passed` | 1 s |

`make quick` — the three verifications plus both C++ builds — took **30.4 s**
end to end, run from a freshly assembled copy of this package on 2026-09-03.
All 240 exported Python files byte-compile under Python 3.9.6.

## The languages, and why

Python (about 12 100 lines) does everything. C++17 (about 880 lines, two files)
appears twice, both times for a stated reason:

- `sim/cpp/rheston_mc.cpp` — the rough-Heston Monte-Carlo. The Volterra
  convolution costs `O(N²)` per path: at `N = 800` steps and 200 000 paths that
  is `6.4 × 10¹⁰` multiplications, which NumPy does not do in a reasonable time.
  No dependencies, not even BLAS — the 8×8 least-squares solve of the
  Longstaff–Schwartz step is a hand-written Cholesky, so the binary cannot drift
  with a library version.
- `sim/cpp/heston_lattice.cpp` — a mirror of `sim/heston_lattice.py`, kept as a
  cross-check. The two agree bit-for-bit; the C++ is not an alternative
  implementation but a witness.

## Author

Benjamin Gottesman Berdah.

## Licence

MIT, see `LICENSE`.
