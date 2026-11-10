"""Experiment for the Numerical PDEs Lecture 11 post: von Neumann analysis.

Python standard library only. Writes experiment.json next to this file.

Part 1: the grid mode W_j = exp(i j h xi) is an eigenvector of the central difference D0.
        Check the eigenvalue (i/h) sin(xi h) and its error for a few points per wavelength.
Part 2: measure the amplification factor g(xi) of FTCS, Backward Euler and Crank-Nicolson
        for u_t = u_xx by doing ONE step on a periodic grid from a single Fourier mode.
        Compare with the formulas.
Part 3: the same measurement for two schemes for u_t + a u_x = 0 (FTCS and Lax-Friedrichs),
        and a long run from a smooth bump.
Part 4: 2D Crank-Nicolson: condition number of I - (k/2) L_h from the eigenvalue formula,
        with the largest eigenvalue checked by power iteration.

Run: python scripts/posts/numerical-pde/lecture-11/experiment.py
"""

import cmath
import json
import math
from pathlib import Path

OUT = {}

# ---------------- Part 1: eigenvalue of D0 ----------------
p1 = {"check": [], "table": []}
N = 40
h = 1 / N
for m in (1, 3, 7):
    xi = 2 * math.pi * m  # periodic on [0, 1]
    W = [cmath.exp(1j * j * h * xi) for j in range(N)]
    lam = 1j / h * math.sin(xi * h)
    dev = max(abs((W[(j + 1) % N] - W[j - 1]) / (2 * h) - lam * W[j]) for j in range(N)) / abs(lam)
    p1["check"].append({"m": m, "rel_dev": dev})
for M in (64, 32, 16, 8, 6, 4, 3, 2):  # grid points per wavelength
    th = 2 * math.pi / M  # xi h
    ratio = math.sin(th) / th  # discrete eigenvalue / exact eigenvalue i xi
    p1["table"].append({"M": M, "xih": th, "ratio": ratio, "rel_err": 1 - ratio, "taylor": th * th / 6})
OUT["part1"] = p1


# ---------------- small dense complex solver ----------------
def solve(Amat, b):
    n = len(b)
    A = [row[:] + [b[i]] for i, row in enumerate(Amat)]
    for c in range(n):
        piv = max(range(c, n), key=lambda r: abs(A[r][c]))
        A[c], A[piv] = A[piv], A[c]
        for r in range(c + 1, n):
            f = A[r][c] / A[c][c]
            if f != 0:
                for k in range(c, n + 1):
                    A[r][k] -= f * A[c][k]
    x = [0j] * n
    for r in range(n - 1, -1, -1):
        x[r] = (A[r][n] - sum(A[r][k] * x[k] for k in range(r + 1, n))) / A[r][r]
    return x


def periodic_matrix(n, left, mid, right):
    M = [[0j] * n for _ in range(n)]
    for j in range(n):
        M[j][(j - 1) % n] += left
        M[j][j] += mid
        M[j][(j + 1) % n] += right
    return M


def step_heat(U, r, scheme):
    n = len(U)
    lap = [U[j - 1] - 2 * U[j] + U[(j + 1) % n] for j in range(n)]  # h^2 * u_xx
    if scheme == "FTCS":
        return [U[j] + r * lap[j] for j in range(n)]
    if scheme == "BE":
        return solve(periodic_matrix(n, -r, 1 + 2 * r, -r), U)
    if scheme == "CN":
        rhs = [U[j] + r / 2 * lap[j] for j in range(n)]
        return solve(periodic_matrix(n, -r / 2, 1 + r, -r / 2), rhs)


def g_formula(scheme, r, th):
    c = math.cos(th) - 1
    if scheme == "FTCS":
        return 1 + 2 * r * c
    if scheme == "BE":
        return 1 / (1 - 2 * r * c)
    if scheme == "CN":
        return (1 + r * c) / (1 - r * c)


# ---------------- Part 2: measured g for the heat equation ----------------
N2 = 32
h2 = 1 / N2
p2 = {"N": N2, "runs": {}}
cases = [("FTCS", 0.3), ("FTCS", 0.5), ("FTCS", 0.6), ("BE", 6.0), ("CN", 0.6), ("CN", 6.0), ("CN", 60.0)]
for scheme, r in cases:
    rows = []
    worst_spread = 0.0
    worst_formula = 0.0
    for m in range(0, N2 // 2 + 1):  # xi h = 2 pi m / N from 0 to pi
        th = 2 * math.pi * m / N2
        U = [cmath.exp(1j * j * th) for j in range(N2)]
        V = step_heat(U, r, scheme)
        ratios = [V[j] / U[j] for j in range(N2)]
        g = ratios[0]
        worst_spread = max(worst_spread, max(abs(q - g) for q in ratios))
        worst_formula = max(worst_formula, abs(g - g_formula(scheme, r, th)))
        rows.append({"m": m, "xih": th, "g_re": g.real, "g_im": g.imag})
    gmax = max(abs(complex(q["g_re"], q["g_im"])) for q in rows)
    p2["runs"][f"{scheme} r={r:g}"] = {
        "scheme": scheme, "r": r, "rows": rows, "max_abs_g": gmax,
        "ratio_spread": worst_spread, "formula_dev": worst_formula,
    }
OUT["part2"] = p2


# ---------------- Part 3: advection, FTCS and Lax-Friedrichs ----------------
def step_adv(U, nu, scheme):
    n = len(U)
    if scheme == "FTCS":
        return [U[j] - nu / 2 * (U[(j + 1) % n] - U[j - 1]) for j in range(n)]
    if scheme == "LF":
        return [(U[j - 1] + U[(j + 1) % n]) / 2 - nu / 2 * (U[(j + 1) % n] - U[j - 1]) for j in range(n)]


p3 = {"N": N2, "g": {}, "runs": {}}
for scheme in ("FTCS", "LF"):
    for nu in (0.3, 0.9, 1.0, 1.1):
        gmax = 0.0
        for m in range(0, N2 // 2 + 1):
            th = 2 * math.pi * m / N2
            U = [cmath.exp(1j * j * th) for j in range(N2)]
            V = step_adv(U, nu, scheme)
            gmax = max(gmax, abs(V[0] / U[0]))
        p3["g"][f"{scheme} nu={nu:g}"] = gmax

# long run: smooth bump on a periodic grid with 100 points
N3 = 100
x = [j / N3 for j in range(N3)]
bump = [math.exp(-100 * (xj - 0.5) ** 2) for xj in x]
for scheme, nu in (("FTCS", 0.3), ("FTCS", 0.9), ("LF", 0.9), ("LF", 1.1)):
    U = bump[:]
    sizes = [max(abs(u) for u in U)]
    for n in range(300):
        U = step_adv(U, nu, scheme)
        sizes.append(max(abs(u) for u in U))
    p3["runs"][f"{scheme} nu={nu:g}"] = {"sizes": sizes}
OUT["part3"] = p3


# ---------------- Part 4: 2D Crank-Nicolson ----------------
def lap2(u, n, h):
    """5-point Laplacian on the n x n interior grid, zero Dirichlet boundary. u is a flat list."""
    out = [0.0] * (n * n)
    s = 1 / (h * h)
    for i in range(n):
        for j in range(n):
            c = i * n + j
            v = -4 * u[c]
            if i > 0: v += u[c - n]
            if i < n - 1: v += u[c + n]
            if j > 0: v += u[c - 1]
            if j < n - 1: v += u[c + 1]
            out[c] = s * v
    return out


def dot(a, b):
    return sum(p * q for p, q in zip(a, b))


def power_max(apply, n, iters=3000):
    v = [math.sin(0.37 * c + 1.0) + 0.1 for c in range(n)]
    lam = 0.0
    for _ in range(iters):
        w = apply(v)
        lam = math.sqrt(dot(w, w) / dot(v, v))
        s = math.sqrt(dot(w, w))
        v = [wi / s for wi in w]
    return lam


p4 = {"cond": []}
n4 = 31
h4 = 1 / (n4 + 1)
NN = n4 * n4
lam_lap = lambda p: (2 / h4 ** 2) * (1 - math.cos(p * math.pi * h4))  # 1D, positive
lmin1, lmax1 = lam_lap(1), lam_lap(n4)
kappa_lap = (2 * lmax1) / (2 * lmin1)
p4["h"] = h4
lap_max_measured = power_max(lambda u: lap2(u, n4, h4), NN, 1500)
p4["laplacian"] = {"kappa_formula": kappa_lap, "lmax_formula": 2 * lmax1, "lmax_power": lap_max_measured}
for label, k in (("k = 4h", 4 * h4), ("k = h", h4), ("k = h/4", h4 / 4), ("k = h^2", h4 * h4)):
    big = 1 + k / 2 * 2 * lmax1
    small = 1 + k / 2 * 2 * lmin1
    op = lambda u, k=k: [ui - k / 2 * li for ui, li in zip(u, lap2(u, n4, h4))]
    big_measured = power_max(op, NN, 1500)
    p4["cond"].append({"label": label, "k": k, "kappa": big / small, "lmax_formula": big,
                       "lmax_power": big_measured, "k_over_h2": k / h4 ** 2})

OUT["part4"] = p4

Path(__file__).with_name("experiment.json").write_text(json.dumps(OUT, indent=1), encoding="utf-8")

# ---------------- print a summary ----------------
print("Part 1 check:", p1["check"])
for row in p1["table"]:
    print(f"  M={row['M']:3d} xih={row['xih']:.4f} ratio={row['ratio']:.5f} err={row['rel_err']:.3e} taylor={row['taylor']:.3e}")
print("Part 2:")
for key, run in p2["runs"].items():
    g_at = {rw["m"]: rw["g_re"] for rw in run["rows"]}
    print(f"  {key:12s} max|g|={run['max_abs_g']:.4f} spread={run['ratio_spread']:.1e} formula={run['formula_dev']:.1e} "
          f"g(pi/4)={g_at[4]:.4f} g(pi/2)={g_at[8]:.4f} g(pi)={g_at[16]:.4f}")
print("Part 3 max|g|:", {k: round(v, 5) for k, v in p3["g"].items()})
for key, run in p3["runs"].items():
    s = run["sizes"]
    print(f"  {key:10s} n=0 {s[0]:.4f} n=50 {s[50]:.4g} n=100 {s[100]:.4g} n=200 {s[200]:.4g} n=300 {s[300]:.4g} "
          f"rate(last100)={(s[300] / s[200]) ** 0.01:.5f}")
print("Part 4 laplacian:", p4["laplacian"])
for c in p4["cond"]:
    print("  ", c)
