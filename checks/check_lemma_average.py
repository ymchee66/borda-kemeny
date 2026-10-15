"""Exact randomized check of Lemma 3.1 in rational arithmetic.

Notation as in the paper. For random scoring vectors x_1 > x_2 > ... > x_n with rational gaps, random
distributions of J on {1, ..., n} and every s in {1, ..., n}, it checks
    |xbar^{-1}(E[x_J]) - s| <= A(x) * E|J - s|,
where xbar is the piecewise-linear interpolation of x and A(x) is the gap irregularity.
Uses only the Python standard library.
Usage: python3 check_lemma_average.py [trials] [seed]
"""
import random
import sys
from fractions import Fraction as Fr


def gap_irregularity(x):
    """A(x) = max over i < j of (x_i - x_j) / ((j - i) * min_{i <= r < j} delta_r), with x[1..n]."""
    n = len(x) - 1
    return max(Fr(x[i] - x[j]) / ((j - i) * min(x[r] - x[r + 1] for r in range(i, j)))
               for i in range(1, n + 1) for j in range(i + 1, n + 1))


def xbar_inverse(x, v):
    """The point u in [1, n] with xbar(u) = v."""
    n = len(x) - 1
    for r in range(1, n):
        if x[r + 1] <= v <= x[r]:
            return r + (x[r] - v) / (x[r] - x[r + 1])
    raise ValueError(v)


rng = random.Random(int(sys.argv[2]) if len(sys.argv) > 2 else 1)
trials = int(sys.argv[1]) if len(sys.argv) > 1 else 20000
worst = Fr(0)
for _ in range(trials):
    n = rng.randint(2, 8)
    gaps = [Fr(rng.randint(1, 30), rng.randint(1, 10)) for _ in range(n - 1)]
    x = [None] + [sum(gaps[r - 1:], Fr(0)) for r in range(1, n + 1)]   # x[n] = 0
    w = [rng.randint(0, 5) for _ in range(n)]
    if sum(w) == 0:
        w[rng.randrange(n)] = 1
    p = {j: Fr(w[j - 1], sum(w)) for j in range(1, n + 1)}
    a = xbar_inverse(x, sum(p[j] * x[j] for j in p))
    A = gap_irregularity(x)
    for s in range(1, n + 1):
        lhs = abs(a - s)
        rhs = A * sum(p[j] * abs(j - s) for j in p)
        assert lhs <= rhs, (x, p, s)
        if rhs > 0:
            worst = max(worst, lhs / rhs)
print(f"{trials} random instances checked exactly; largest observed ratio of the two sides: {float(worst):.6f}")
