"""Experiment for the post "Numerical PDEs, Lecture 3".

All parts use the test solution u(x) = exp(2x) on [0, 1], so
u'' = f = 4 exp(2x), u'(0) = 2, u'(1) = 2 e^2, u(0) = 1, u(1) = e^2.

Part 1 (Dirichlet at both ends): checks the error equation A E = -tau.
The theory says E_j is close to h^2 * g(x_j) with
  g(x) = (1/12) * (-u''(x) + u''(0) + x (u''(1) - u''(0))).

Part 2 (Neumann at x = 0, Dirichlet at x = 1): three ways to discretize
u'(0) = sigma, for n = 4, 8, ..., 4096:
  one-sided 1st order:  (U_1 - U_0)/h = sigma
  ghost node:           (U_1 - U_0)/h = sigma + (h/2) f(x_0)
  one-sided 2nd order:  (-3 U_0 + 4 U_1 - U_2)/(2h) = sigma
The third row has three unknowns. We remove U_2 with the interior equation
at j = 1, which gives (U_1 - U_0)/h = sigma + (h/2) f(x_1). So every system
stays tridiagonal and the Thomas algorithm solves it.

Part 3 (pure Neumann at both ends): every row of the matrix adds up to 0,
so A times the constant vector is 0. A solution exists only if the
trapezoid sum of f equals sigma_1 - sigma_0. We measure the mismatch, fix
it by shifting f by a constant, pin U_0 = u(0), and measure the error.

Writes the numbers to neumann.json next to this file.

Python standard library only.

Run: python scripts/posts/numerical-pde/lecture-3/neumann.py
"""

import json
import math
from pathlib import Path


def u(x):
    return math.exp(2 * x)


def f(x):
    return 4 * math.exp(2 * x)


SIGMA0, SIGMA1 = 2.0, 2 * math.exp(2)
ALPHA, BETA = 1.0, math.exp(2)


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


def order(errors):
    return [math.log2(a / b) for a, b in zip(errors, errors[1:])]


KS = list(range(2, 13))  # n = 2^k

# --- Part 1: Dirichlet at both ends, error equation ---

def g(x):
    return (-f(x) + f(0) + x * (f(1) - f(0))) / 12


part1 = []
for k in KS:
    n = 2**k
    h = 1 / n
    xs = [j * h for j in range(n + 1)]
    m = n - 1
    rhs = [f(xs[j]) for j in range(1, n)]
    rhs[0] -= ALPHA / h**2
    rhs[-1] -= BETA / h**2
    U = thomas([1 / h**2] * m, [-2 / h**2] * m, [1 / h**2] * m, rhs)
    E = [U[j - 1] - u(xs[j]) for j in range(1, n)]
    emax = max(abs(e) for e in E)
    shape_gap = max(abs(E[j - 1] / h**2 - g(xs[j])) for j in range(1, n))
    gmax = max(abs(g(xs[j])) for j in range(1, n))
    part1.append({"n": n, "global": emax, "max_E_over_h2": emax / h**2,
                  "max_g": gmax, "shape_gap": shape_gap})


# --- Part 2: Neumann at x = 0, three boundary rows ---

def solve_neumann(n, method):
    """Unknowns U_0..U_{n-1}; U_n = BETA is moved to the right-hand side.

    The first row is scaled to (U_1 - U_0)/h = s (coefficients -1, 1 times 1/h),
    and the interior rows keep the 1/h^2 scaling.
    """
    h = 1 / n
    xs = [j * h for j in range(n + 1)]
    if method == "one-sided, 1st order":
        s = SIGMA0
    elif method == "ghost node":
        s = SIGMA0 + h / 2 * f(xs[0])
    else:  # one-sided, 2nd order, with U_2 removed using row j = 1
        s = SIGMA0 + h / 2 * f(xs[1])
    lower = [0.0] + [1 / h**2] * (n - 1)
    diag = [-1 / h] + [-2 / h**2] * (n - 1)
    upper = [1 / h] + [1 / h**2] * (n - 1)
    rhs = [s] + [f(xs[j]) for j in range(1, n)]
    rhs[-1] -= BETA / h**2
    U = thomas(lower, diag, upper, rhs) + [BETA]
    return xs, U


METHODS = ["one-sided, 1st order", "ghost node", "one-sided, 2nd order"]
part2 = {name: [] for name in METHODS}
for k in KS:
    n = 2**k
    for name in METHODS:
        xs, U = solve_neumann(n, name)
        E = [U[j] - u(xs[j]) for j in range(n + 1)]
        part2[name].append({
            "n": n,
            "global": max(abs(e) for e in E),
            "at_0": E[0],
            "at_half": E[n // 2],
            "at_0_over_h": E[0] * n,
        })

# Check the one-sided 2nd-order row directly: plug U into (-3U_0 + 4U_1 - U_2)/(2h).
xs, U = solve_neumann(64, "one-sided, 2nd order")
h = 1 / 64
row_check = (-3 * U[0] + 4 * U[1] - U[2]) / (2 * h) - SIGMA0


# --- Part 3: pure Neumann ---

def trapezoid(values, h):
    return h * (values[0] / 2 + sum(values[1:-1]) + values[-1] / 2)


part3 = []
for k in KS:
    n = 2**k
    h = 1 / n
    xs = [j * h for j in range(n + 1)]
    F = [f(x) for x in xs]
    mismatch = trapezoid(F, h) - (SIGMA1 - SIGMA0)

    def march(Fv):
        # Rows 0..n-1 hold exactly; returns U and the residual of the last row.
        U = [u(0.0), u(0.0) + h * (SIGMA0 + h / 2 * Fv[0])]
        for j in range(1, n):
            U.append(2 * U[j] - U[j - 1] + h**2 * Fv[j])
        last = (U[n - 1] - U[n]) / h - (-SIGMA1 + h / 2 * Fv[n])
        return U, last

    _, last_raw = march(F)
    U, last_fixed = march([v - mismatch for v in F])
    err = max(abs(U[j] - u(xs[j])) for j in range(n + 1))
    part3.append({"n": n, "mismatch": mismatch, "last_row_raw": last_raw,
                  "last_row_fixed": last_fixed, "global_fixed": err})

# Row sums of the pure Neumann matrix (scaled rows) for n = 4.
n = 4
h = 1 / n
rows = [[-1 / h, 1 / h, 0, 0, 0]]
for j in range(1, n):
    r = [0.0] * (n + 1)
    r[j - 1], r[j], r[j + 1] = 1 / h**2, -2 / h**2, 1 / h**2
    rows.append(r)
rows.append([0, 0, 0, 1 / h, -1 / h])
row_sums = [sum(r) for r in rows]

out = {
    "k": KS,
    "part1": part1,
    "part1_order": order([r["global"] for r in part1]),
    "part2": part2,
    "part2_order": {name: order([r["global"] for r in part2[name]]) for name in METHODS},
    "row_check_2nd_order": row_check,
    "part3": part3,
    "part3_mismatch_order": order([abs(r["mismatch"]) for r in part3]),
    "part3_order": order([r["global_fixed"] for r in part3]),
    "row_sums_n4": row_sums,
}
Path(__file__).with_name("neumann.json").write_text(json.dumps(out, indent=1), encoding="utf-8")

if __name__ == "__main__":
    for r in part1:
        print("P1", r)
    print("P1 order", [round(o, 3) for o in out["part1_order"]])
    for name in METHODS:
        for r in part2[name]:
            print("P2", name, r)
        print("P2 order", name, [round(o, 3) for o in out["part2_order"][name]])
    print("row check", row_check)
    for r in part3:
        print("P3", r)
    print("P3 mismatch order", [round(o, 3) for o in out["part3_mismatch_order"]])
    print("P3 order", [round(o, 3) for o in out["part3_order"]])
    print("row sums", row_sums)
