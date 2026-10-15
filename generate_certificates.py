"""Generate exact certificates for the approximation ratio of Borda's rule on n candidates.

For every ranking rho of {0,...,n-1} the linear program of Theorem 6.2 is solved numerically
with the HiGHS solver in SciPy. From the numerical solutions we build
  * a rational number lam (the largest optimum),
  * for every ranking rho, rational multipliers mu_0,...,mu_{n-2} >= 0 (mu_1,...,mu_{n-1} in the paper) satisfying the
    certificate inequality of Lemma 4.1 with lam,
  * an integer profile with the identity as a Borda ranking and quotient lam.
The results are written to certificates/borda_n<n>.json. Correctness does not depend on this
program: run verify.py, which uses exact arithmetic only.

Rankings are stored as position maps: rho[c] is the position of candidate c.
Usage: python3 generate_certificates.py 2 3 4 5 6 7
"""
import itertools, json, math, sys
from fractions import Fraction as Fr
import numpy as np
from scipy.optimize import linprog


def setup(n):
    seqs = list(itertools.permutations(range(n)))
    pos = np.argsort(np.array(seqs), axis=1)          # pos[s][c] = position of candidate c
    a, b = np.triu_indices(n, 1)
    L = (pos[:, a] < pos[:, b]).astype(np.int64)
    return pos, L


def kendall_row(L, i):
    return (L[i][None, :] != L).sum(1)


def solve_exact(rows, rhs, m):
    M = [[Fr(x) for x in r] + [Fr(v)] for r, v in zip(rows, rhs)]
    piv, R = [], 0
    for col in range(m):
        p = next((i for i in range(R, len(M)) if M[i][col] != 0), None)
        if p is None:
            continue
        M[R], M[p] = M[p], M[R]
        for i in range(len(M)):
            if i != R and M[i][col] != 0:
                f = M[i][col] / M[R][col]
                M[i] = [x - f * y for x, y in zip(M[i], M[R])]
        piv.append(col)
        R += 1
    if R < m:
        return None
    y = [Fr(0)] * m
    for i, col in enumerate(piv):
        y[col] = M[i][-1] / M[i][col]
    return y


def certificate_ok(lam, mu, d, c, A):
    D = lam.denominator * math.lcm(*[v.denominator for v in mu], 1)
    mi = [int(v * D) for v in mu]
    lhs = int(lam * D) * d.astype(object) + sum(mi[k] * A[k].astype(object) for k in range(len(mu))) - c.astype(object) * D
    return min(v >= 0 for v in mu) and min(lhs) >= 0


def run(n):
    pos, L = setup(n)
    N = len(pos)
    ident = int(np.where((pos == np.arange(n)).all(1))[0][0])
    c = kendall_row(L, ident)
    A = (pos[:, :-1] - pos[:, 1:]).T                    # A[r][s] = s(r) - s(r+1)
    best = (-1.0, None, None)
    for rho in range(N):
        r = linprog(-c, A_ub=A, b_ub=np.zeros(n - 1), A_eq=kendall_row(L, rho)[None, :], b_eq=[1.0],
                    bounds=(0, None), method="highs")
        if r.status == 0 and -r.fun > best[0]:
            best = (-r.fun, rho, r.x)
    val, rho_star, p = best
    # integer witness from the support of an optimal vertex
    supp = np.nonzero(p > 1e-9)[0]
    witness, lam = None, None
    for D in range(1, 200000):
        w = np.round(p[supp] / p[supp].min() * D).astype(np.int64)
        beta = (pos[supp] * w[:, None]).sum(0)
        if not np.all(beta[:-1] <= beta[1:]):
            continue
        S = np.array([(L[s][None, :] != L[supp]).sum(1) @ w for s in range(N)])
        q = Fr(int(S[ident]), int(S.min()))
        if abs(float(q) - val) < 1e-9:
            witness, lam = [(pos[s].tolist(), int(wt)) for s, wt in zip(supp, w)], q
            break
    if witness is None:
        raise RuntimeError("no integer witness found")
    certs = []
    for rho in range(N):
        d = kendall_row(L, rho)
        rhs = (float(lam) * d - c).astype(float)
        r = linprog(np.ones(n - 1), A_ub=-A.T.astype(float), b_ub=rhs, bounds=(0, None), method="highs")
        y = r.x
        mu = None
        if n == 2 or r.status != 0:
            y = np.zeros(n - 1)
        slack = float(lam) * d + A.T @ y - c
        rows, rr = [], []
        for s in np.argsort(slack):
            if slack[s] > 1e-7:
                break
            rows.append([int(A[k, s]) for k in range(n - 1)])
            rr.append(Fr(int(c[s])) - lam * int(d[s]))
        for k in range(n - 1):
            if y[k] < 1e-9:
                e = [0] * (n - 1); e[k] = 1
                rows.append(e); rr.append(Fr(0))
        cand = solve_exact(rows, rr, n - 1) if rows else None
        if cand is not None and certificate_ok(lam, cand, d, c, A):
            mu = cand
        else:
            for den in [1, 2, 3, 4, 6, 12, 60, 840, 27720]:
                t = [Fr(round(v * den), den) for v in y]
                if certificate_ok(lam, t, d, c, A):
                    mu = t
                    break
        if mu is None:
            raise RuntimeError(f"no certificate for rho={pos[rho].tolist()}")
        certs.append({"rho": pos[rho].tolist(), "mu": [str(v) for v in mu]})
    out = {"n": n, "lambda": str(lam),
           "witness": [{"ranking": r, "weight": w} for r, w in witness],
           "certificates": certs}
    with open(f"certificates/borda_n{n}.json", "w") as f:
        json.dump(out, f)
    print(f"n={n}: lambda={lam}, {len(certs)} certificates written")


if __name__ == "__main__":
    for arg in sys.argv[1:] or ["2", "3", "4", "5", "6", "7"]:
        run(int(arg))
