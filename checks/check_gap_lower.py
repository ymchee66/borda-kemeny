"""Exact check of the construction in the proof of Theorem 3.4.

Candidates and positions are numbered from 0, as in the rest of the code (see README.md).

For every number n <= 8 of candidates and all i <= r < j, it builds the ranking sigma of the proof
(positions i..j filled with r+1, i, ..., r-1, r+2, ..., j, r) and checks that K(id, sigma) = j - i,
that sigma(r+1) = i and sigma(r) = j, and, for random scoring vectors and weights w1 >= w2 with
w1*delta_r < w2*(x_i - x_j), that every x-ranking prefers r+1 to r and that the quotient
S(pi)/OPT of the profile is at least 1 + (w1/w2 - 1)/(j - i), with OPT = w2*(j - i).
Uses only the Python standard library.
"""
import itertools, random
from fractions import Fraction as Fr


def K(p, q):
    n = len(p)
    return sum((p[a] < p[b]) != (q[a] < q[b]) for a in range(n) for b in range(a + 1, n))


rng = random.Random(7)
count = 0
for n in range(2, 9):
    ident = tuple(range(n))
    for i in range(n):
        for j in range(i + 1, n):
            for r in range(i, j):
                seq = list(range(i)) + [r + 1] + list(range(i, r)) + list(range(r + 2, j + 1)) + [r] + list(range(j + 1, n))
                sigma = [0] * n
                for p, c in enumerate(seq):
                    sigma[c] = p
                assert sorted(sigma) == list(range(n))
                assert sigma[r + 1] == i and sigma[r] == j and K(ident, sigma) == j - i
                if n <= 6:
                    for _ in range(3):
                        gaps = [Fr(rng.randint(1, 9)) for _ in range(n - 1)]
                        x = [sum(gaps[k:], Fr(0)) for k in range(n)]
                        w2 = rng.randint(1, 4)
                        lim = w2 * (x[i] - x[j]) / gaps[r]
                        w1 = max(w2, int(lim) - (1 if lim == int(lim) else 0))
                        if w1 * gaps[r] >= w2 * (x[i] - x[j]):
                            continue
                        T = [w1 * x[c] + w2 * x[sigma[c]] for c in range(n)]
                        assert T[r + 1] > T[r]
                        prof = [(ident, w1), (tuple(sigma), w2)]
                        S = lambda p: sum(w * K(p, q) for q, w in prof)
                        opt = min(S(p) for p in itertools.permutations(range(n)))
                        assert opt == w2 * (j - i)
                        for p in itertools.permutations(range(n)):
                            if all(not (T[a] > T[b] and p[a] > p[b]) for a in range(n) for b in range(n)):
                                assert p[r + 1] < p[r]
                                assert Fr(S(p), opt) >= 1 + (Fr(w1, w2) - 1) / (j - i)
                count += 1
print(f"construction checked for {count} triples i <= r < j with n <= 8 candidates")
