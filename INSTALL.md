# INSTALL — the environment, as measured

Every version below was read off the machine that produced the campaigns in
`runs/`, not copied from a specification.

## The machine the results came from

| | |
|---|---|
| OS | macOS 26.6.2 |
| architecture | `arm64` (Apple M-series) |
| Python | **3.9.6** (the system interpreter) |
| C++ compiler | **Apple clang 21.0.0** (`clang-2100.1.1.101`), target `arm64-apple-darwin25.6.0` |

Nothing here needs macOS. The Python side is portable as written; the C++ side
needs a C++17 compiler and `-pthread`. See the caveat on `-march=native` below.

## Python

Four packages, and only four:

```
numpy==2.0.2
scipy==1.13.1
mpmath==1.3.0
sympy==1.14.0
```

`environment/requirements.lock` pins these exactly; `REQUIREMENTS.txt` is the
looser declaration the project used (`scipy>=1.13`). The lock file is what was
actually installed.

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r environment/requirements.lock
```

Conda users: `environment/environment.yml` declares the same four.

No plotting library is required to run the campaigns. Figures in the paper were
produced separately; the campaigns write text and JSON.

## C++

Two independent binaries, no build system, no dependencies:

```bash
mkdir -p sim/cpp/build
clang++ -O3 -march=native -std=c++17 -pthread -o sim/cpp/build/rheston_mc     sim/cpp/rheston_mc.cpp
clang++ -O3 -march=native -std=c++17 -pthread -o sim/cpp/build/heston_lattice sim/cpp/heston_lattice.cpp
```

Measured on the machine above: **2.0 s** and **1.4 s** respectively, no
warnings, binaries about 56 kB each. `g++` should work identically — it has not
been tested here, and this file will say so until it has been.

### The `-march=native` caveat, stated rather than hidden

`-march=native` lets the compiler use the host's instruction set, including
fused multiply-add. Floating-point results can therefore differ in the last
bits between two machines with different CPUs, even with identical source. This
is the intended trade: the Monte-Carlo exists because it must be fast.

If you want a build that is more likely to match across machines, drop
`-march=native`. Expect it to be slower, and expect the last digits to move
either way — neither build is "the correct one", and the paper's claims are
stated to a precision that tolerates this. Where a claim depends on a specific
digit, `RUNBOOK.md` says so for that campaign.

## Checking the install

```bash
python3 verify/verify_claims.py            # 88/88 passed, ~21 s — Python side
sim/cpp/build/rheston_mc     < /dev/null   # C++ side
sim/cpp/build/heston_lattice < /dev/null
```

Each binary on empty input prints its CSV header row and exits 0 — measured, not
assumed. `rheston_mc` prints:

```
id,eu_put,eu_put_se,eu_call,eu_call_se,am_put,am_put_se,am_call,am_call_se,am_put_insample,am_call_insample,mean_ST,mean_VT,neg_hits,ex_dates,t_eu,t_am
```

and `heston_lattice`:

```
id,walk,n,eu_put,eu_call,am_put,am_call,grid,offsets,driver_states,violations,max_var_err,max_mean_err,max_abs_drift,v_at_upper,feller,secs
```

The C++ kernels read a CSV on stdin and write a CSV on stdout; they are driven
by `sim/run_rheston_tables.py` and `sim/run_rheston_american_anchor.py`. Their
column layout is documented in `sim/cpp/README.md`. Do not call them by hand for
anything meant to end up in a document.
