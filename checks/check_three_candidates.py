"""Symbolic check of the certificates in the proof of Theorem 4.2 (three candidates).

Scoring vector gamma+1, 1, 0 with gamma >= 1, so y = (0, gamma, gamma+1) and lambda = 1 + gamma/2.
Rankings are position maps written sigma(0) sigma(1) sigma(2). For each comparison ranking rho
the multipliers mu_0, mu_1 of the paper are substituted, the residual
    lambda*K(rho,sigma) + mu_0*(y[sigma(0)]-y[sigma(1)]) + mu_1*(y[sigma(1)]-y[sigma(2)]) - K(id,sigma)
is multiplied by the stated positive factor and written as a polynomial in h = gamma - 1.
The script checks that every coefficient is nonnegative and prints the lists displayed in the paper.
Requires sympy.
"""
import itertools
import sympy as sp

g, h = sp.symbols("gamma h", positive=True)
y = [0, g, g + 1]
lam = 1 + g / 2
order = [(0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0)]
mult = {(0, 1, 2): (0, 0), (0, 2, 1): (0, 1), (1, 0, 2): (1 / g, 0),
        (1, 2, 0): (2 / (2 * g + 1), 4 / (2 * g + 1)), (2, 0, 1): (4 / (g + 2), 2 / (g + 2)),
        (2, 1, 0): (3 / (g + 1), 3 / (g + 1))}
factor = {(0, 1, 2): 2, (0, 2, 1): 2, (1, 0, 2): 2 * g, (1, 2, 0): 2 * (2 * g + 1),
          (2, 0, 1): 2 * (g + 2), (2, 1, 0): 2 * (g + 1)}


def K(p, q):
    return sum((p[a] < p[b]) != (q[a] < q[b]) for a, b in itertools.combinations(range(3), 2))


ok = True
for rho in order:
    m0, m1 = mult[rho]
    row = []
    for s in order:
        res = lam * K(rho, s) + m0 * (y[s[0]] - y[s[1]]) + m1 * (y[s[1]] - y[s[2]]) - K((0, 1, 2), s)
        poly = sp.expand(sp.cancel(res * factor[rho]).subs(g, h + 1))
        ok &= all(c >= 0 for c in sp.Poly(poly, h).coeffs()) if poly != 0 else True
        row.append(poly)
    print("".join(map(str, rho)), row)
print("all residuals have nonnegative coefficients:", ok)
assert ok
