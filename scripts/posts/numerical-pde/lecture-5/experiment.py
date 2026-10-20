"""Experiment for the Lecture 5 post: central vs upwind differences, and how to measure an order.

Problem: -eps u'' + u' = 1 on (0, 1) with u(0) = u(1) = 0, eps = 0.02.
Exact solution: u(x) = x - (exp((x - 1)/eps) - exp(-1/eps)) / (1 - exp(-1/eps)).

Part 1: pure diffusion -eps u'' = 1. The centered scheme is exact at the nodes.
Part 2: small grids. The limit eps -> 0 with h = 1/5 (exact fractions), and the
        actual solutions for eps = 0.02 with n = 5, 10, 20, 40, 80.
Part 3: convergence table for central and upwind, h = 2^-2 ... 2^-12,
        in three norms: max, discrete L2 (with h), and plain l2 (without h).

Standard library only (tridiagonal systems by the Thomas algorithm).
Run: python scripts/posts/numerical-pde/lecture-5/experiment.py   (writes experiment.json next to this file)
"""

import json
import math
from fractions import Fraction
from pathlib import Path

EPS = 0.02


def exact(x, eps=EPS):
    return x - (math.exp((x - 1) / eps) - math.exp(-1 / eps)) / (1 - math.exp(-1 / eps))


def thomas(a, b, c, d):
    """Solve a tridiagonal system. a: lower, b: diagonal, c: upper (lists of length m)."""
    m = len(b)
    cp, dp = [0.0] * m, [0.0] * m
    cp[0], dp[0] = c[0] / b[0], d[0] / b[0]
    for i in range(1, m):
        den = b[i] - a[i] * cp[i - 1]
        cp[i] = c[i] / den if i < m - 1 else 0.0
        dp[i] = (d[i] - a[i] * dp[i - 1]) / den
    x = [0.0] * m
    x[-1] = dp[-1]
    for i in range(m - 2, -1, -1):
        x[i] = dp[i] - cp[i] * x[i + 1]
    return x


def rows(n, eps, scheme):
    """Coefficients of row i: lo*U[i-1] + di*U[i] + up*U[i+1] = 1."""
    h = 1.0 / n
    lo, di, up = -eps / h**2, 2 * eps / h**2, -eps / h**2
    if scheme == "central":
        lo, up = lo - 1 / (2 * h), up + 1 / (2 * h)
    elif scheme == "upwind":
        lo, di = lo - 1 / h, di + 1 / h
    return lo, di, up


def solve(n, eps, scheme):
    lo, di, up = rows(n, eps, scheme)
    m = n - 1
    U = thomas([lo] * m, [di] * m, [up] * m, [1.0] * m)
    # residual check, boundary values are 0
    V = [0.0] + U + [0.0]
    res = max(abs(lo * V[i - 1] + di * V[i] + up * V[i + 1] - 1.0) for i in range(1, n))
    return U, res


def errors(n, U, u):
    h = 1.0 / n
    E = [u(i * h) - U[i - 1] for i in range(1, n)]
    s = sum(e * e for e in E)
    return {"max": max(abs(e) for e in E), "l2h": math.sqrt(h * s), "l2": math.sqrt(s)}


out = {"eps": EPS}

# Part 1: pure diffusion, exact solution x(1 - x)/(2 eps)
part1 = []
for k in range(2, 13):
    n = 2**k
    U, res = solve(n, EPS, "diffusion")
    e = errors(n, U, lambda x: x * (1 - x) / (2 * EPS))
    part1.append({"n": n, "max": e["max"], "peak": max(U), "residual": res})
out["part1"] = part1

# Part 2a: eps -> 0 limit with h = 1/5, central rows (U[i+1] - U[i-1])/(2h) = 1, exact fractions.
# Unknowns U1..U4. Odd count of unknowns would make this limit system singular.
h = Fraction(1, 5)
U2 = 2 * h  # row 1: U2 - U0 = 2h
U3 = -2 * h  # row 4: U5 - U3 = 2h
U1 = U3 - 2 * h  # row 2: U3 - U1 = 2h
U4 = U2 + 2 * h  # row 3: U4 - U2 = 2h
limit = [U1, U2, U3, U4]
out["part2_limit"] = [str(v) for v in limit]

# Part 2b: eps = 0.02 on small grids
part2 = []
for n in (5, 10, 20, 25, 40, 50, 80):
    row = {"n": n, "x": [i / n for i in range(1, n)], "exact": [exact(i / n) for i in range(1, n)]}
    for scheme in ("central", "upwind"):
        U, res = solve(n, EPS, scheme)
        row[scheme] = U
        row[scheme + "_min"] = min(U)
        row[scheme + "_residual"] = res
        # number of sign changes of U[i+1] - U[i]: a monotone-then-drop curve has 1
        d = [U[i + 1] - U[i] for i in range(len(U) - 1)]
        row[scheme + "_turns"] = sum(1 for i in range(len(d) - 1) if d[i] * d[i + 1] < 0)
        row[scheme + "_max_err"] = errors(n, U, exact)["max"]
    part2.append(row)
out["part2"] = part2

# Part 3: convergence table
ks = list(range(2, 13))
part3 = {"k": ks}
for scheme in ("central", "upwind"):
    table = []
    for k in ks:
        n = 2**k
        U, res = solve(n, EPS, scheme)
        e = errors(n, U, exact)
        e["n"] = n
        e["residual"] = res
        table.append(e)
    for norm in ("max", "l2h", "l2"):
        for j in range(1, len(table)):
            table[j]["order_" + norm] = math.log(table[j - 1][norm] / table[j][norm]) / math.log(2)
    part3[scheme] = table
out["part3"] = part3

Path(__file__).with_name("experiment.json").write_text(json.dumps(out, indent=1), encoding="utf-8")

if __name__ == "__main__":
    print("Part 1 (pure diffusion): max error, residual")
    for r in part1:
        print(f"  n={r['n']:5d}  {r['max']:.2e}  peak {r['peak']:.6f}  res {r['residual']:.1e}")
    print("Part 2a (eps -> 0, h = 1/5):", [str(v) for v in limit], " u = x at nodes: 1/5 2/5 3/5 4/5")
    print("Part 2b (eps = 0.02):")
    for r in part2:
        print(f"  n={r['n']:3d}  central min {r['central_min']:+.4f} turns {r['central_turns']:2d} err {r['central_max_err']:.3e} | "
              f"upwind min {r['upwind_min']:+.4f} turns {r['upwind_turns']} err {r['upwind_max_err']:.3e} | res {max(r['central_residual'], r['upwind_residual']):.1e}")
    r = part2[0]
    print("  n=5 central:", [f"{v:.4f}" for v in r["central"]], "upwind:", [f"{v:.4f}" for v in r["upwind"]], "exact:", [f"{v:.4f}" for v in r["exact"]])
    print("Part 3:")
    for scheme in ("central", "upwind"):
        print(" ", scheme)
        for r in part3[scheme]:
            o = lambda nm: f"{r['order_' + nm]:5.2f}" if "order_max" in r else "   - "
            print(f"   n={r['n']:5d}  max {r['max']:.3e} ({o('max')})  l2h {r['l2h']:.3e} ({o('l2h')})  l2 {r['l2']:.3e} ({o('l2')})  res {r['residual']:.1e}")
