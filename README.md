# Positional scoring rules and Kemeny rank aggregation: computations

Supplementary material for the paper

> Y. M. Chee, *Characterising Positional Scoring Rules with Bounded Kemeny Approximation*.

Theorem and section numbers below refer to that paper.

## Contents

| File | Purpose |
|---|---|
| `generate_certificates.py` | Solves the linear programs of Theorem 6.3 for Borda's rule with the HiGHS solver in SciPy and writes rational certificates and integer extremal profiles to `certificates/`. |
| `verify.py` | Exact verifier for Theorem 7.1. Uses only the Python standard library. |
| `certificates/borda_n<n>.json` | Certificates and extremal profile for `n` = 2, ..., 7 candidates. |
| `checks/check_three_candidates.py` | Symbolic check of the multipliers in the proof of Theorem 4.2 (needs SymPy). |
| `checks/check_lemma_average.py` | Exact randomised check of Lemma 3.1 in rational arithmetic. |
| `checks/check_gap_lower.py` | Exact check of the construction in the proof of Theorem 3.4. |
| `checks/conjecture.c` | Exhaustive check of Conjecture 5.5 for all rankings of at most 10 candidates. |

The checks of Lemma 3.1 and Theorem 3.4 test the general proofs on examples. The proofs do not depend on them.

## Conventions

A ranking of the candidates `0, ..., n-1` is stored as a *position map*: `ranking[c]` is the position of candidate `c`, and position `0` is first. Borda totals are sums of positions, so a ranking is a Borda ranking when it orders the candidates by nondecreasing Borda total.

Each file `certificates/borda_n<n>.json` contains

* `lambda`: the value of the approximation ratio `R_n` of Borda's rule, as a fraction;
* `witness`: an integer profile, as a list of rankings with their numbers of voters, whose Borda rankings include the identity and whose quotient of the Kemeny score of the identity by the optimum equals `lambda`;
* `certificates`: for each of the `n!` rankings `rho`, multipliers `mu_0, ..., mu_{n-2}` as fractions satisfying the inequality of Lemma 4.1 with `lambda`.

## Reproducing the results

Requirements: Python 3.9 or later. `generate_certificates.py` also needs NumPy and SciPy, and `check_three_candidates.py` needs SymPy. `verify.py` and the other checks use only the standard library.

```
python3 verify.py                     # checks n = 2, ..., 7 (about a minute)
python3 checks/check_three_candidates.py
python3 checks/check_lemma_average.py 20000
python3 checks/check_gap_lower.py
cc -O2 -o conjecture checks/conjecture.c
for n in 3 4 5 6 7 8 9 10; do ./conjecture $n; done
```

To regenerate the certificates (about ten minutes for `n = 7`):

```
pip install numpy scipy
python3 generate_certificates.py 2 3 4 5 6 7
```

The solver is used only to find the certificates and profiles. `verify.py` checks them in exact integer and rational arithmetic, so the results do not depend on floating-point tolerances.

## Expected output

`verify.py` reports `R_n` = 1, 3/2, 5/3, 75/41, 2, 64/31 for `n` = 2, ..., 7, with 5912 certificates and 25,935,016 inequalities in total. `conjecture` reports the largest value `(n-3)/(n-1)` for `4 <= n <= 10`.

## Licence

MIT. See `LICENSE`.
