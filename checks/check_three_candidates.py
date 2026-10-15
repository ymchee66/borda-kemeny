"""Symbolic check of the certificates in the proof of Theorem 4.2 (three candidates).

Notation as in the paper: candidates and positions are 1, 2, 3, and a ranking sigma is written
sigma(1) sigma(2) sigma(3). The scoring vector is gamma+1, 1, 0 with gamma >= 1, so y_1, y_2, y_3 is
0, gamma, gamma+1 and lambda = 1 + gamma/2. For each comparison ranking rho the multipliers mu_1, mu_2
of the paper are substituted, the residual
    Delta = lambda*K(rho,sigma) + mu_1*(y[sigma(1)]-y[sigma(2)]) + mu_2*(y[sigma(2)]-y[sigma(3)]) - K(iota,sigma)
is multiplied by the positive factor of the paper and written as a polynomial in h = gamma - 1.
The script checks that every coefficient is nonnegative and prints the columns of the table in
Appendix A. Requires sympy.
"""
import itertools
import sympy as sp

g, h = sp.symbols("gamma h", positive=True)
y = {1: 0, 2: g, 3: g + 1}
lam = 1 + g / 2
order = ["123", "132", "213", "231", "312", "321"]
mult = {"123": (0, 0), "132": (0, 1), "213": (1 / g, 0),
        "231": (2 / (2 * g + 1), 4 / (2 * g + 1)), "312": (4 / (g + 2), 2 / (g + 2)),
        "321": (3 / (g + 1), 3 / (g + 1))}
factor = {"123": 2, "132": 2, "213": 2 * g, "231": 2 * (2 * g + 1), "312": 2 * (g + 2), "321": 2 * (g + 1)}


def pos(word):
    """Position map of a ranking written sigma(1) sigma(2) sigma(3)."""
    return {c: int(word[c - 1]) for c in (1, 2, 3)}


def K(p, q):
    return sum((p[a] < p[b]) != (q[a] < q[b]) for a, b in itertools.combinations((1, 2, 3), 2))


iota = pos("123")
ok = True
for rho in order:
    m1, m2 = mult[rho]
    R = pos(rho)
    column = []
    for word in order:
        S = pos(word)
        res = lam * K(R, S) + m1 * (y[S[1]] - y[S[2]]) + m2 * (y[S[2]] - y[S[3]]) - K(iota, S)
        poly = sp.expand(sp.cancel(res * factor[rho]).subs(g, h + 1))
        if poly != 0:
            ok &= all(c >= 0 for c in sp.Poly(poly, h).coeffs())
        column.append(f"{word}: {poly}")
    print(f"rho = {rho}, factor {factor[rho]}:  " + ",  ".join(column))
print("all residuals have nonnegative coefficients:", ok)
assert ok
