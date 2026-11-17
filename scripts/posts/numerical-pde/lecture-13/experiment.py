"""Experiment for the Lecture 13 post: von Neumann analysis for advection, and a first finite volume scheme.

Python standard library only (cmath, math). Periodic grids, a = 1 unless stated.

Part 1: measure the amplification factor g of upwind, Lax-Friedrichs and Lax-Wendroff by one step
        on a single grid mode, and compare with the closed forms and with g = 1 + k*lambda.
        Leapfrog: both roots have |g| = 1 for |nu| < 1, and a double root at nu = 1 makes |U| grow.
Part 2: Lax-Friedrichs is consistent only if its artificial diffusion h^2/(2k) goes to zero.
Part 3: finite volume form. Constant beta: the diffusive flux gives FTCS. Variable beta: the conservative
        (flux) form keeps the total mass to rounding, a non-conservative form does not.

Run: python scripts/posts/numerical-pde/lecture-13/experiment.py   (writes experiment.json next to this file)
"""

import cmath
import json
import math
from pathlib import Path

PI = math.pi
out = {}


# ---------- Part 1: amplification factors ----------

def step_upwind(U, nu):
    m = len(U)
    return [U[j] - nu * (U[j] - U[j - 1]) for j in range(m)]


def step_lf(U, nu):
    m = len(U)
    return [0.5 * (U[j - 1] + U[(j + 1) % m]) - 0.5 * nu * (U[(j + 1) % m] - U[j - 1]) for j in range(m)]


def step_lw(U, nu):
    m = len(U)
    return [U[j] - 0.5 * nu * (U[(j + 1) % m] - U[j - 1]) + 0.5 * nu * nu * (U[j - 1] - 2 * U[j] + U[(j + 1) % m])
            for j in range(m)]


STEPS = {"upwind": step_upwind, "Lax–Friedrichs": step_lf, "Lax–Wendroff": step_lw}


def g_formula(name, nu, th):
    if name == "upwind":
        return (1 - nu) + nu * cmath.exp(-1j * th)
    if name == "Lax–Friedrichs":
        return math.cos(th) - 1j * nu * math.sin(th)
    return 1 - 1j * nu * math.sin(th) + nu * nu * (math.cos(th) - 1)


def g2_closed(name, nu, th):
    if name == "upwind":
        return 1 - 2 * nu * (1 - nu) * (1 - math.cos(th))
    if name == "Lax–Friedrichs":
        return 1 - (1 - nu * nu) * math.sin(th) ** 2
    return 1 - nu * nu * (1 - nu * nu) * (1 - math.cos(th)) ** 2


def c_of(name, nu):
    # c = 2 k eps / h^2 from Lecture 12: LF 1, upwind nu, LW nu^2
    return {"upwind": nu, "Lax–Friedrichs": 1.0, "Lax–Wendroff": nu * nu}[name]


def measure_g(name, nu, M):
    """One step on each mode e^{i theta j}, theta = 2 pi p / M, p = 0..M/2. Returns thetas, g, max spread."""
    thetas, gs, spread = [], [], 0.0
    for p in range(M // 2 + 1):
        th = 2 * PI * p / M
        U0 = [cmath.exp(1j * th * j) for j in range(M)]
        U1 = STEPS[name](U0, nu)
        ratios = [U1[j] / U0[j] for j in range(M)]
        g = ratios[0]
        spread = max(spread, max(abs(r - g) for r in ratios))
        thetas.append(th)
        gs.append(g)
    return thetas, gs, spread


M1 = 64
nu1 = 0.6
part1 = {"M": M1, "nu": nu1, "schemes": {}}
for name in STEPS:
    th, gs, spread = measure_g(name, nu1, M1)
    err_formula = max(abs(g - g_formula(name, nu1, t)) for g, t in zip(gs, th))
    err_closed = max(abs(abs(g) ** 2 - g2_closed(name, nu1, t)) for g, t in zip(gs, th))
    c = c_of(name, nu1)
    err_klam = max(abs(g - (1 + (-1j * nu1 * math.sin(t) - c * (1 - math.cos(t))))) for g, t in zip(gs, th))
    picks = {}
    for label, p in (("pi/4", M1 // 8), ("pi/2", M1 // 4), ("pi", M1 // 2)):
        picks[label] = abs(gs[p])
    part1["schemes"][name] = {
        "theta": th, "abs_g": [abs(g) for g in gs],
        "spread": spread, "err_formula": err_formula, "err_closed": err_closed, "err_1_plus_klambda": err_klam,
        "abs_g_at": picks,
    }
    print(f"P1 {name:15s} nu={nu1}: spread {spread:.1e}, |g-formula| {err_formula:.1e}, "
          f"||g|^2-closed| {err_closed:.1e}, |g-(1+k lam)| {err_klam:.1e}, |g| at pi/4,pi/2,pi = "
          + ", ".join(f"{v:.4f}" for v in picks.values()))

# small-theta damping per step: 1 - |g| at theta = 2 pi / 64 (64 points per wavelength)
damp = {name: 1 - part1["schemes"][name]["abs_g"][1] for name in STEPS}
part1["damping_64ppw"] = damp
print("P1 damping per step at 64 points per wavelength:", {k: f"{v:.2e}" for k, v in damp.items()})

# max |g| over theta for several nu
table = {}
for name in STEPS:
    row = {}
    for nu in (0.5, 1.0, 1.1, -0.3):
        _, gs, _ = measure_g(name, nu, M1)
        row[str(nu)] = max(abs(g) for g in gs)
    table[name] = row
    print(f"P1 max|g| {name:15s}", {k: f"{v:.4f}" for k, v in row.items()})
part1["max_abs_g"] = table

# leapfrog roots and growth
lf_root_dev = 0.0
for p in range(M1 // 2 + 1):
    th = 2 * PI * p / M1
    s = 0.9 * math.sin(th)
    r = cmath.sqrt(1 - s * s)
    for g in (-1j * s + r, -1j * s - r):
        lf_root_dev = max(lf_root_dev, abs(abs(g) - 1))
part1["leapfrog_root_dev_nu09"] = lf_root_dev
# leapfrog |g| (larger root) at nu1 on the same thetas, for the chart
lf_abs = []
for t in part1["schemes"]["upwind"]["theta"]:
    s = nu1 * math.sin(t)
    r = cmath.sqrt(1 - s * s)
    lf_abs.append(max(abs(-1j * s + r), abs(-1j * s - r)))
part1["leapfrog_abs_g"] = lf_abs
print(f"P1 leapfrog nu=0.9: max ||g|-1| over both roots {lf_root_dev:.1e}")


def leapfrog_growth(nu, M=8, n_list=(10, 100, 1000)):
    th = PI / 2  # mode with sin(theta) = 1, M = 8 so p = 2
    U_old = [cmath.exp(1j * th * j) for j in range(M)]
    U = list(U_old)  # crude start U^1 = U^0, so the start is not an exact eigen-solution
    res = {}
    for n in range(1, max(n_list) + 1):
        U_new = [U_old[j] - nu * (U[(j + 1) % M] - U[j - 1]) for j in range(M)]
        U_old, U = U, U_new
        if n + 1 in n_list:
            res[n + 1] = max(abs(u) for u in U)
    return res


lg = {str(nu): leapfrog_growth(nu) for nu in (0.9, 1.0)}
part1["leapfrog_growth"] = lg
print("P1 leapfrog max|U| at n = 10, 100, 1000:", {k: [f"{v:.3f}" for v in d.values()] for k, d in lg.items()})
out["part1"] = part1


# ---------- Part 2: LF consistency depends on the k-h relation ----------

def run_lf(M, k, T):
    h = 1.0 / M
    nu = k / h
    n = round(T / k)
    U = [math.sin(2 * PI * j * h) for j in range(M)]
    for _ in range(n):
        U = step_lf(U, nu)
    t = n * k
    return h, t, U


T2 = 0.25
part2 = {"T": T2, "cases": {}}
cases = {
    "k = h/2": (lambda h: h / 2, (20, 40, 80, 160)),
    "k = 4h^2": (lambda h: 4 * h * h, (20, 40, 80, 160)),
    "k = 4h^3": (lambda h: 4 * h ** 3, (10, 20, 40, 80)),
}
for label, (kfun, Ms) in cases.items():
    rows = []
    for M in Ms:
        h = 1.0 / M
        k = kfun(h)
        n = round(T2 / k)
        k = T2 / n
        h, t, U = run_lf(M, k, T2)
        eps = h * h / (2 * k)
        adv = max(abs(U[j] - math.sin(2 * PI * (j * h - t))) for j in range(M))
        ad = max(abs(U[j] - math.exp(-eps * 4 * PI * PI * t) * math.sin(2 * PI * (j * h - t))) for j in range(M))
        amp = max(abs(u) for u in U)
        rows.append({"M": M, "steps": n, "eps": eps, "err_advection": adv, "err_adv_diff": ad, "max_U": amp})
        print(f"P2 {label:9s} M={M:4d} steps={n:7d} eps={eps:.4g} err vs advection {adv:.3e}, "
              f"vs u_t+u_x=eps u_xx {ad:.3e}, max|U| {amp:.4f}")
    part2["cases"][label] = rows
out["part2"] = part2


# ---------- Part 3: finite volume, conservative vs non-conservative ----------

def beta(x):
    return 1 + 0.5 * math.sin(2 * PI * x)


def dbeta(x):
    return 0.5 * 2 * PI * math.cos(2 * PI * x)


def q0(x):
    return math.exp(-60 * (x - 0.4) ** 2)


# 3a: constant beta, the flux form equals FTCS
M = 50
h = 1.0 / M
b = 0.7
k = 0.4 * h * h / b
Q = [q0((i + 0.5) * h) for i in range(M)]
F = [-b * (Q[i] - Q[i - 1]) / h for i in range(M)]  # F[i] = flux at x_{i-1/2}
fv = [Q[i] - k / h * (F[(i + 1) % M] - F[i]) for i in range(M)]
ftcs = [Q[i] + b * k / (h * h) * (Q[i - 1] - 2 * Q[i] + Q[(i + 1) % M]) for i in range(M)]
d3a = max(abs(x - y) for x, y in zip(fv, ftcs))
print(f"P3a constant beta: |FV flux step - FTCS step| max {d3a:.1e}")
# upwind flux a*Q_{i-1} gives the upwind scheme
a, nu = 1.0, 0.6
Fa = [a * Q[i - 1] for i in range(M)]
fv_up = [Q[i] - nu / a * (Fa[(i + 1) % M] - Fa[i]) for i in range(M)]
d3a_up = max(abs(x - y) for x, y in zip(fv_up, step_upwind(Q, nu)))
print(f"P3a flux a*Q_(i-1): |FV step - upwind step| max {d3a_up:.1e}")


def run_diffusion(M, T, conservative):
    h = 1.0 / M
    k = 0.25 * h * h / 1.5
    n = round(T / k)
    k = T / n
    xc = [(i + 0.5) * h for i in range(M)]
    Q = [q0(x) for x in xc]
    mass0 = h * sum(Q)
    bf = [beta(i * h) for i in range(M)]  # beta at x_{i-1/2} = i h
    bc = [beta(x) for x in xc]
    dbc = [dbeta(x) for x in xc]
    for _ in range(n):
        if conservative:
            F = [-bf[i] * (Q[i] - Q[i - 1]) / h for i in range(M)]
            Q = [Q[i] - k / h * (F[(i + 1) % M] - F[i]) for i in range(M)]
        else:
            Q = [Q[i] + k * (bc[i] * (Q[i - 1] - 2 * Q[i] + Q[(i + 1) % M]) / (h * h)
                             + dbc[i] * (Q[(i + 1) % M] - Q[i - 1]) / (2 * h)) for i in range(M)]
    return n, mass0, h * sum(Q), Q


T3 = 0.02
part3 = {"T": T3, "const_beta_diff": d3a, "upwind_flux_diff": d3a_up, "rows": []}
sols = {}
for M in (20, 40, 80, 160):
    n, m0, mc, Qc = run_diffusion(M, T3, True)
    _, _, mn, Qn = run_diffusion(M, T3, False)
    gap = max(abs(x - y) for x, y in zip(Qc, Qn))
    row = {"M": M, "steps": n, "mass0": m0, "rel_change_conservative": (mc - m0) / m0,
           "rel_change_nonconservative": (mn - m0) / m0, "max_gap": gap}
    part3["rows"].append(row)
    print(f"P3b M={M:4d} steps={n:6d} mass0 {m0:.6f}  rel change: conservative {row['rel_change_conservative']:+.1e}, "
          f"non-conservative {row['rel_change_nonconservative']:+.2e}; max gap between the two {gap:.2e}")
out["part3"] = part3

Path(__file__).with_name("experiment.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
