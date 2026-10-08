"""Experiment for the post "Numerical PDEs, Lecture 2".

Part 1: finds finite difference weights by the method of undetermined
coefficients, with exact fractions, and checks that the weights add up to 0.

Part 2: solves the two-point boundary value problem u'' = f on (0, 1),
u(0) = alpha, u(1) = beta, with the 3-point difference and the Thomas
algorithm, for n = 4, 8, ..., 4096. Two test solutions:
  smooth: u(x) = sin(3x) + x^2 + 1
  rough:  u(x) = x^(5/2)   (u'''' is unbounded near x = 0)
For each n it records the largest global error max|U_j - u(x_j)| and the
largest local truncation error max|tau_j|.
Writes the numbers to bvp.json next to this file.

Python standard library only.

Run: python scripts/posts/numerical-pde/lecture-2/bvp.py
"""

import json
import math
from fractions import Fraction
from pathlib import Path


# --- Part 1: undetermined coefficients ---

def weights(offsets, m):
    """Weights w_j with sum_j w_j u(x + s_j h) = h^m u^(m)(x) + higher-order terms.

    Matches Taylor coefficients of order 0 .. len(offsets)-1 and solves the
    linear system exactly with Gaussian elimination on fractions.
    """
    n = len(offsets)
    rows = []
    for k in range(n):
        row = [Fraction(s) ** k / math.factorial(k) for s in offsets]
        rows.append(row + [Fraction(1 if k == m else 0)])
    for col in range(n):
        pivot = next(r for r in range(col, n) if rows[r][col] != 0)
        rows[col], rows[pivot] = rows[pivot], rows[col]
        for r in range(n):
            if r != col and rows[r][col] != 0:
                factor = rows[r][col] / rows[col][col]
                rows[r] = [a - factor * b for a, b in zip(rows[r], rows[col])]
    return [rows[i][n] / rows[i][i] for i in range(n)]


def first_leftover(offsets, w, m):
    """First Taylor order k > m whose term does not cancel, and its coefficient."""
    k = len(offsets)
    while True:
        c = sum(wj * Fraction(s) ** k for wj, s in zip(w, offsets)) / math.factorial(k)
        if c != 0:
            return k, c
        k += 1


STENCILS = [
    ("one-sided, 3 points", [0, -1, -2], 1),
    ("one-sided, 4 points", [0, -1, -2, -3], 1),
    ("central, 2 points", [-1, 1], 1),
    ("central, 4 points", [-2, -1, 1, 2], 1),
    ("second derivative, 3 points", [-1, 0, 1], 2),
    ("third derivative, 4 points", [-1, 0, 1, 2], 3),
    ("third derivative, central 4 points", [-2, -1, 1, 2], 3),
]

stencils = []
for name, offsets, m in STENCILS:
    w = weights(offsets, m)
    k, c = first_leftover(offsets, w, m)
    stencils.append({
        "name": name,
        "offsets": offsets,
        "derivative": m,
        "weights": [str(x) for x in w],
        "weight_sum": str(sum(w)),
        "order": k - m,  # leftover h^k term, divided by h^m
        "leading": str(c),
    })


# --- Part 2: two-point boundary value problem ---

def thomas(lower, diag, upper, rhs):
    """Solves a tridiagonal system. lower[0] and upper[-1] are not used."""
    n = len(diag)
    c = [0.0] * n
    d = [0.0] * n
    c[0] = upper[0] / diag[0]
    d[0] = rhs[0] / diag[0]
    for i in range(1, n):
        m = diag[i] - lower[i] * c[i - 1]
        c[i] = upper[i] / m if i < n - 1 else 0.0
        d[i] = (rhs[i] - lower[i] * d[i - 1]) / m
    x = [0.0] * n
    x[-1] = d[-1]
    for i in range(n - 2, -1, -1):
        x[i] = d[i] - c[i] * x[i + 1]
    return x


def solve_bvp(f, alpha, beta, n):
    """3-point scheme for u'' = f, u(0) = alpha, u(1) = beta. Returns interior U_1..U_{n-1}."""
    h = 1.0 / n
    xs = [j * h for j in range(1, n)]
    rhs = [f(x) for x in xs]
    rhs[0] -= alpha / h**2   # known boundary values move to the right side
    rhs[-1] -= beta / h**2
    m = n - 1
    lower = [1 / h**2] * m
    diag = [-2 / h**2] * m
    upper = [1 / h**2] * m
    return xs, thomas(lower, diag, upper, rhs)


CASES = {
    "smooth": {
        "u": lambda x: math.sin(3 * x) + x * x + 1,
        "f": lambda x: -9 * math.sin(3 * x) + 2,
        "max_u4": 81.0,  # u'''' = 81 sin(3x), and 3x reaches pi/2 inside [0, 1]
    },
    "rough": {
        "u": lambda x: x ** 2.5,
        "f": lambda x: 3.75 * x ** 0.5,
        "max_u4": None,
    },
}

ks = list(range(2, 13))
results = {}
for name, case in CASES.items():
    u, f = case["u"], case["f"]
    rows = []
    for k in ks:
        n = 2**k
        h = 1.0 / n
        xs, U = solve_bvp(f, u(0.0), u(1.0), n)
        global_err = max(abs(Uj - u(x)) for Uj, x in zip(U, xs))
        taus = [(u(x - h) - 2 * u(x) + u(x + h)) / h**2 - f(x) for x in xs]
        lte = max(abs(t) for t in taus)
        worst = max(range(len(taus)), key=lambda j: abs(taus[j]))
        row = {"n": n, "h": h, "global": global_err, "lte": lte, "lte_at": xs[worst]}
        if case["max_u4"]:
            row["lte_bound"] = h * h / 12 * case["max_u4"]
        rows.append(row)
    for i in range(1, len(rows)):
        rows[i]["global_order"] = math.log2(rows[i - 1]["global"] / rows[i]["global"])
        rows[i]["lte_order"] = math.log2(rows[i - 1]["lte"] / rows[i]["lte"])
    results[name] = rows

out = {"k": ks, "stencils": stencils, "bvp": results}
Path(__file__).with_name("bvp.json").write_text(json.dumps(out, indent=1), encoding="utf-8")

for s in stencils:
    print(f"{s['name']:30s} w={s['weights']} sum={s['weight_sum']} order={s['order']} leading={s['leading']}")
for name, rows in results.items():
    print(name)
    for r in rows:
        print(f"  n={r['n']:5d} global={r['global']:.3e} ({r.get('global_order', 0):.2f})"
              f"  lte={r['lte']:.3e} ({r.get('lte_order', 0):.2f}) at x={r['lte_at']:.4f}"
              + (f"  bound={r['lte_bound']:.3e}" if 'lte_bound' in r else ""))
