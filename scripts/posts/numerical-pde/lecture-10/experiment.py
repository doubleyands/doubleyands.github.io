"""Experiment for the Numerical PDEs Lecture 10 post: stability of the heat equation.

Python standard library only. Writes experiment.json next to this file.

Problem: u_t = u_xx on [0, 1], u = 0 at both ends, grid x_j = j h, h = 1/(N+1).
A = (1/h^2) tridiag(1, -2, 1), the N x N second-difference matrix.

Part 1: check the eigenpairs lambda_p = (2/h^2)(cos(p pi h) - 1), u^p_j = sin(p pi j h).
Part 2: 2-norm of A and A^{-1} by power iteration (no formula used), condition number.
Part 3: FTCS with r = k/h^2 = 0.45, 0.50, 0.55 from a hat-shaped start.
Part 4: Crank-Nicolson and Backward Euler with large time steps, same hat start.

Run: python scripts/posts/numerical-pde/lecture-10/experiment.py
"""

import json
import math
from pathlib import Path


def apply_A(u, h):
    n = len(u)
    out = []
    for j in range(n):
        left = u[j - 1] if j > 0 else 0.0
        right = u[j + 1] if j < n - 1 else 0.0
        out.append((left - 2 * u[j] + right) / (h * h))
    return out


def solve_tridiag(a, b, c, d):
    """Thomas algorithm: a = sub, b = main, c = super diagonal (constants), d = right side."""
    n = len(d)
    cp, dp = [0.0] * n, [0.0] * n
    cp[0] = c / b
    dp[0] = d[0] / b
    for i in range(1, n):
        m = b - a * cp[i - 1]
        cp[i] = c / m
        dp[i] = (d[i] - a * dp[i - 1]) / m
    x = [0.0] * n
    x[-1] = dp[-1]
    for i in range(n - 2, -1, -1):
        x[i] = dp[i] - cp[i] * x[i + 1]
    return x


def norm2(v):
    return math.sqrt(sum(x * x for x in v))


def lam(p, h):
    return 2 / h**2 * (math.cos(p * math.pi * h) - 1)


out = {}

# ---------- Part 1: eigenpairs ----------
N = 15
h = 1 / (N + 1)
part1 = []
worst = 0.0
for p in range(1, N + 1):
    u = [math.sin(p * math.pi * j * h) for j in range(1, N + 1)]
    Au = apply_A(u, h)
    l = lam(p, h)
    res = norm2([Au[i] - l * u[i] for i in range(N)]) / (abs(l) * norm2(u))
    worst = max(worst, res)
    part1.append({"p": p, "lambda": l, "continuous": -(p * math.pi) ** 2, "rel_residual": res})
out["part1"] = {"N": N, "rows": part1, "worst_rel_residual": worst}
print("Part 1 (N=15, h=1/16)")
for r in part1:
    if r["p"] in (1, 2, 4, 8, 12, 15):
        print(f"  p={r['p']:2d} lambda={r['lambda']:10.2f}  -(p pi)^2={r['continuous']:10.2f}")
print(f"  worst relative residual {worst:.1e}")

# ---------- Part 2: norms by power iteration ----------
part2 = []
print("Part 2")
for N in (7, 15, 31, 63, 127, 255):
    h = 1 / (N + 1)
    # start vector with all modes present
    v = [1.0 + 0.3 * math.cos(1.7 * j) for j in range(1, N + 1)]
    # ||A||_2: power iteration (A symmetric)
    x = v[:]
    big = 0.0
    for _ in range(20000):
        y = apply_A(x, h)
        nb = norm2(y)
        x = [t / nb for t in y]
        if abs(nb - big) < 1e-13 * nb:
            big = nb
            break
        big = nb
    # ||A^{-1}||_2: inverse iteration with the Thomas algorithm
    x = v[:]
    small = 0.0
    for _ in range(2000):
        y = solve_tridiag(1 / h**2, -2 / h**2, 1 / h**2, x)
        ns = norm2(y)
        x = [t / ns for t in y]
        if abs(ns - small) < 1e-14 * ns:
            small = ns
            break
        small = ns
    kappa = big * small
    row = {"N": N, "h": h, "normA": big, "normAinv": small, "kappa": kappa,
           "formula_kappa": abs(lam(N, h)) / abs(lam(1, h)), "approx": 4 / (math.pi**2 * h**2)}
    part2.append(row)
    print(f"  h=1/{N+1:4d} ||A||={big:12.1f} (4/h^2={4/h**2:10.0f}) ||A^-1||={small:.5f} "
          f"kappa={kappa:10.1f} formula={row['formula_kappa']:10.1f} 4/(pi^2h^2)={row['approx']:10.1f}")
out["part2"] = part2

# ---------- Part 3: FTCS ----------
N = 24
h = 1 / (N + 1)
xs = [j * h for j in range(1, N + 1)]
hat = [min(x, 1 - x) * 2 for x in xs]  # peak 1 at x = 0.5
steps = 200
part3 = {"N": N, "steps": steps, "runs": {}}
print("Part 3 (FTCS, h=1/25)")
for r in (0.45, 0.50, 0.55):
    k = r * h * h
    U = hat[:]
    sizes = [max(abs(t) for t in U)]
    for n in range(steps):
        AU = apply_A(U, h)
        U = [U[i] + k * AU[i] for i in range(N)]
        sizes.append(max(abs(t) for t in U))
    worst_factor = max(abs(1 + k * lam(p, h)) for p in range(1, N + 1))
    worst_odd = max(abs(1 + k * lam(p, h)) for p in range(1, N + 1, 2))
    part3["runs"][f"{r:.2f}"] = {"k": k, "sizes": sizes, "max_factor": worst_factor, "max_factor_odd": worst_odd}
    late = (sizes[200] / sizes[150]) ** (1 / 50)
    print(f"  r={r:.2f} k={k:.6f} max|U| n=50 {sizes[50]:.3e} n=100 {sizes[100]:.3e} n=200 {sizes[200]:.3e}"
          f"  growth/step (n 150..200) {late:.4f}  max|1+k lam| all p {worst_factor:.4f} odd p {worst_odd:.4f}")
out["part3"] = part3


# ---------- Part 4: CN and BE with large steps ----------
def exact_hat(x, t, terms=4000):
    # Fourier sine series of the hat 2 min(x, 1-x): b_m = 8 sin(m pi/2)/(m pi)^2
    s = 0.0
    for m in range(1, terms, 2):
        s += 8 * math.sin(m * math.pi / 2) / (m * math.pi) ** 2 * math.exp(-(m * math.pi) ** 2 * t) * math.sin(m * math.pi * x)
    return s


N = 49
h = 1 / (N + 1)
xs = [j * h for j in range(1, N + 1)]
hat = [min(x, 1 - x) * 2 for x in xs]
T = 0.1
ex = [exact_hat(x, T) for x in xs]
part4 = []
print("Part 4 (h=1/50, T=0.1, hat start)")
for nsteps in (1000, 100, 20, 10, 5, 2):
    k = T / nsteps
    r = k / h**2
    res = {"steps": nsteps, "k": k, "r": r}
    for name in ("CN", "BE"):
        U = hat[:]
        for _ in range(nsteps):
            if name == "CN":
                AU = apply_A(U, h)
                rhs = [U[i] + 0.5 * k * AU[i] for i in range(N)]
                U = solve_tridiag(-0.5 * k / h**2, 1 + k / h**2, -0.5 * k / h**2, rhs)
            else:
                U = solve_tridiag(-k / h**2, 1 + 2 * k / h**2, -k / h**2, U)
        res[name] = max(abs(U[i] - ex[i]) for i in range(N))
    zN = k * lam(N, h)
    z1 = k * lam(1, h)
    res["g_CN_top"] = (1 + zN / 2) / (1 - zN / 2)
    res["g_CN_low"] = (1 + z1 / 2) / (1 - z1 / 2)
    res["g_BE_top"] = 1 / (1 - zN)
    part4.append(res)
    print(f"  steps={nsteps:5d} k={k:.4f} r={r:8.1f}  err CN {res['CN']:.2e}  BE {res['BE']:.2e}"
          f"  g_CN(p=N) {res['g_CN_top']:+.4f}  g_BE(p=N) {res['g_BE_top']:.4f}  g_CN(p=1) {res['g_CN_low']:.4f}")
print(f"  exact max at T: {max(ex):.5f}")
out["part4"] = {"N": N, "T": T, "rows": part4, "exact_peak": max(ex)}

Path(__file__).with_name("experiment.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
