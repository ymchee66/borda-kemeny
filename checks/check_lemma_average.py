"""Exact randomised check of Lemma 3.1 in rational arithmetic.

For random strictly increasing rational sequences y_0 < ... < y_{n-1}, random distributions of J on
{0,...,n-1} and every s, it checks |y^{-1}(E[y_J]) - s| <= A(y) * E|J - s|, where y is extended by
linear interpolation. Uses only the Python standard library.
Usage: python3 check_lemma_average.py [trials]
"""
import random, sys
from fractions import Fraction as Fr


def A(y):
    n = len(y)
    return max(Fr(y[j] - y[i], 1) / ((j - i) * min(y[r + 1] - y[r] for r in range(i, j)))
               for i in range(n) for j in range(i + 1, n))


def yinv(y, v):
    for r in range(len(y) - 1):
        if y[r] <= v <= y[r + 1]:
            return r + (v - y[r]) / (y[r + 1] - y[r])
    raise ValueError


rng = random.Random(int(sys.argv[2]) if len(sys.argv) > 2 else 1)
trials = int(sys.argv[1]) if len(sys.argv) > 1 else 20000
worst = Fr(0)
for _ in range(trials):
    n = rng.randint(2, 8)
    gaps = [Fr(rng.randint(1, 30), rng.randint(1, 10)) for _ in range(n - 1)]
    y = [Fr(0)]
    for d in gaps:
        y.append(y[-1] + d)
    w = [rng.randint(0, 5) for _ in range(n)]
    if sum(w) == 0:
        w[rng.randrange(n)] = 1
    tot = sum(w)
    p = [Fr(v, tot) for v in w]
    a = yinv(y, sum(pi * yi for pi, yi in zip(p, y)))
    Ay = A(y)
    for s in range(n):
        lhs = abs(a - s)
        rhs = Ay * sum(pi * abs(j - s) for j, pi in enumerate(p))
        assert lhs <= rhs, (y, p, s)
        if rhs > 0:
            worst = max(worst, lhs / rhs)
print(f"{trials} random instances checked exactly; largest observed ratio of the two sides: {float(worst):.6f}")
