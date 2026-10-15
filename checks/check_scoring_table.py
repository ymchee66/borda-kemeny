"""Exact check of Corollary 3.6 (Tables 1 and 2).

Notation as in the paper: x = x_1, ..., x_n with x_1 >= ... >= x_n. For each scoring vector the script
computes in rational arithmetic
  * the gap irregularity A(x) from its definition,
  * the best lower bound 1 + ((x_i - x_j)/delta_r - 1)/(j - i) of Theorem 3.4 over all i <= r < j,
  * the upper bound 1 + 4A(x) of Theorem 3.3,
and compares them with the entries of Table 2. Vectors with a zero gap are reported as infinite.
For Borda, Dowdall and x_r = 2^{-r} it checks the formulas of Table 2 for 3 <= n <= 12.
Uses only the Python standard library.
"""
from fractions import Fraction as Fr

INF = None


def bounds(x):
    """Return (A(x), lower bound, upper bound) for x given as a list x_1, ..., x_n."""
    n = len(x)
    X = [None] + [Fr(v) for v in x]
    d = {r: X[r] - X[r + 1] for r in range(1, n)}
    if min(d.values()) == 0:
        return INF, INF, INF
    A = max((X[i] - X[j]) / ((j - i) * min(d[r] for r in range(i, j)))
            for i in range(1, n + 1) for j in range(i + 1, n + 1))
    low = max(1 + ((X[i] - X[j]) / d[r] - 1) / (j - i)
              for i in range(1, n + 1) for j in range(i + 1, n + 1) for r in range(i, j))
    return A, low, 1 + 4 * A


TABLE1 = {
    "Eurovision Song Contest": [12, 10, 8, 7, 6, 5, 4, 3, 2, 1],
    "Formula One": [25, 18, 15, 12, 10, 8, 6, 4, 2, 1],
    "MotoGP": [25, 20, 16, 13, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1],
    "FIS World Cup": [100, 80, 60, 50, 45, 40, 36, 32, 29, 26, 24, 22, 20, 18, 16] + list(range(15, 0, -1)),
    "MLB Most Valuable Player": [14, 9, 8, 7, 6, 5, 4, 3, 2, 1],
    "NBA Most Valuable Player": [10, 7, 5, 3, 1],
}
TABLE2 = {  # (name, n): (A, lower, upper)
    ("Eurovision Song Contest", 10): ("5/3", "7/3", "23/3"), ("Eurovision Song Contest", 11): ("5/3", "7/3", "23/3"),
    ("Formula One", 10): ("8/3", "32/9", "35/3"), ("Formula One", 11): ("8/3", "32/9", "35/3"),
    ("MotoGP", 15): ("3", "19/5", "13"), ("MotoGP", 16): ("3", "19/5", "13"),
    ("FIS World Cup", 30): ("17/3", "33/5", "71/3"), ("FIS World Cup", 31): ("17/3", "33/5", "71/3"),
    ("MLB Most Valuable Player", 10): ("3", "7/2", "13"), ("MLB Most Valuable Player", 11): ("3", "7/2", "13"),
    ("NBA Most Valuable Player", 5): ("5/4", "15/8", "6"), ("NBA Most Valuable Player", 6): ("2", "14/5", "9"),
}

for name, pts in TABLE1.items():
    k = len(pts)
    for n in (k, k + 1, k + 2, k + 5):
        x = pts + [0] * (n - k)
        A, low, up = bounds(x)
        if n >= k + 2:
            assert A is INF, (name, n)
            print(f"{name:26s} n = {n:2d}: zero gap, R = infinity")
        else:
            exp = tuple(Fr(v) for v in TABLE2[(name, n)])
            assert (A, low, up) == exp, (name, n, A, low, up)
            print(f"{name:26s} n = {n:2d}: A = {A}, {low} <= R <= {up}")

for n in range(3, 13):
    for k in range(1, n):
        assert bounds([1] * k + [0] * (n - k))[0] is INF          # k-approval
    assert bounds([n - r for r in range(1, n + 1)]) == (1, 2 - Fr(1, n - 1), 5)
    assert bounds([Fr(1, r) for r in range(1, n + 1)]) == (n - 1, n - Fr(1, n - 1), 4 * n - 3)
    assert bounds([Fr(1, 2 ** r) for r in range(1, n + 1)]) == (
        Fr(2 ** (n - 1) - 1, n - 1), 1 + Fr(2 ** (n - 1) - 2, n - 1), 1 + Fr(4 * (2 ** (n - 1) - 1), n - 1))
print("Borda, k-approval, Dowdall and 2^{-r} agree with Table 2 for 3 <= n <= 12")
print("all entries of Table 2 confirmed")
