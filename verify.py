"""Exact verifier for the certificates of Theorem 7.1. Uses only the Python standard library.

For each file certificates/borda_n<n>.json it checks that
  1. there is exactly one certificate for each of the n! rankings rho,
  2. all multipliers are nonnegative,
  3. every certificate satisfies, for every ranking sigma,
        lam*K(rho,sigma) + sum_r mu_r*(sigma(r) - sigma(r+1)) >= K(id,sigma),
  4. the identity is a Borda ranking of the witness profile, and
  5. the Kemeny score of the identity divided by the optimum of the witness profile,
     computed over all n! rankings, equals lam.
Rankings are position maps: ranking[c] is the position of candidate c.
Usage: python3 verify.py [n ...]
"""
import itertools, json, sys
from fractions import Fraction as Fr
from math import gcd


def orientation(p):
    """Bit mask with one bit per pair a < b, set when p prefers a to b."""
    n = len(p)
    m, bit = 0, 1
    for a in range(n):
        for b in range(a + 1, n):
            if p[a] < p[b]:
                m |= bit
            bit <<= 1
    return m


def kendall(p, q):
    return bin(orientation(p) ^ orientation(q)).count("1")


def verify(n):
    data = json.load(open(f"certificates/borda_n{n}.json"))
    assert data["n"] == n
    lam = Fr(data["lambda"])
    rankings = [tuple(p) for p in itertools.permutations(range(n))]
    ident = tuple(range(n))
    mask = {s: orientation(s) for s in rankings}
    base = [bin(mask[s] ^ mask[ident]).count("1") for s in rankings]
    diffs = [[s[r] - s[r + 1] for r in range(n - 1)] for s in rankings]
    seen = set()
    ineq = 0
    for cert in data["certificates"]:
        rho = tuple(cert["rho"])
        mu = [Fr(v) for v in cert["mu"]]
        assert rho in mask and rho not in seen, "duplicate or invalid rho"
        seen.add(rho)
        assert len(mu) == n - 1 and all(v >= 0 for v in mu), "negative multiplier"
        # clear denominators: D*lam, D*mu_r and D are integers
        D = lam.denominator
        for v in mu:
            D = D * v.denominator // gcd(D, v.denominator)
        L = int(lam * D)
        M = [int(v * D) for v in mu]
        mr = mask[rho]
        for s, k0, df in zip(rankings, base, diffs):
            lhs = L * bin(mr ^ mask[s]).count("1") + sum(M[r] * df[r] for r in range(n - 1))
            assert lhs >= D * k0, f"certificate fails: rho={rho}, sigma={s}"
            ineq += 1
    assert seen == set(rankings), "missing certificates"
    wit = [(tuple(w["ranking"]), w["weight"]) for w in data["witness"]]
    assert all(wt > 0 for _, wt in wit)
    beta = [sum(wt * s[c] for s, wt in wit) for c in range(n)]
    assert all(beta[c] <= beta[c + 1] for c in range(n - 1)), "identity is not a Borda ranking"
    score = lambda p: sum(wt * kendall(p, s) for s, wt in wit)
    opt = min(score(p) for p in rankings)
    assert opt > 0 and Fr(score(ident), opt) == lam, "witness quotient differs from lambda"
    print(f"n={n}: R_n = {lam} verified ({len(seen)} certificates, {ineq} inequalities, "
          f"witness with {len(wit)} distinct rankings, Kemeny score {score(ident)}, optimum {opt})")
    return len(seen), ineq


if __name__ == "__main__":
    ns = [int(a) for a in sys.argv[1:]] or [2, 3, 4, 5, 6, 7]
    tc = ti = 0
    for n in ns:
        c, i = verify(n)
        tc += c; ti += i
    print(f"total: {tc} certificates, {ti} inequalities")
