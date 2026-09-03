# Reproduction package — the two entry points.
#
#   make quick   the symbolic verifications and the C++ build: ~30 s, no campaign
#   make cpp     just the two C++ binaries
#   make verify  just the three symbolic verifications
#
# There is deliberately no `make all`. Re-running every campaign is a decision
# with a cost, and RUNBOOK.md states that cost per campaign so you can choose.

CXX      ?= clang++
CXXFLAGS ?= -O3 -march=native -std=c++17 -pthread
PY       ?= python3

.PHONY: quick cpp verify clean

quick: cpp verify

cpp: sim/cpp/build/rheston_mc sim/cpp/build/heston_lattice

sim/cpp/build/rheston_mc: sim/cpp/rheston_mc.cpp
	@mkdir -p sim/cpp/build
	$(CXX) $(CXXFLAGS) -o $@ $<

sim/cpp/build/heston_lattice: sim/cpp/heston_lattice.cpp
	@mkdir -p sim/cpp/build
	$(CXX) $(CXXFLAGS) -o $@ $<

# Measured on the author's machine: 22 s, 3 s, 1 s.
verify:
	$(PY) verify/verify_claims.py
	$(PY) verify/verify_lotc_moments.py
	$(PY) verify/verify_mc_band.py

clean:
	rm -rf sim/cpp/build
