"""Experiment for the Numerical PDEs Lecture 12 post: advection and artificial diffusion.

Python standard library only (cmath for complex numbers). Writes experiment.json next to this file.
Equation: u_t + a u_x = 0 on [0, 1) with periodic boundaries.

Part 1: the central-difference MOL matrix A is skew-symmetric. Check A u_p = lambda_p u_p
        for the Fourier eigenvectors, and that every lambda_p is purely imaginary.
Part 2: Lax-Friedrichs, written as FTCS for u_t + a u_x = eps u_xx with eps = h^2/(2k),
        gives the same numbers. Then the eigenvalue ellipse: for a fixed Courant number nu,
        count the scaled eigenvalues k*lambda_p outside the Forward Euler disk |1 + z| <= 1
        as the diffusion c = 2 k eps / h^2 changes. Theory: none outside iff nu^2 <= c <= 1.
Part 3: one full trip around the domain with Lax-Friedrichs, upwind, Lax-Wendroff and
        leapfrog, from a smooth bump and a box. Also the order of convergence on a sine wave.
Part 4: Forward Euler (FTCS) with k = h^2, and upwind used the wrong way (a < 0).

Run: python scripts/posts/numerical-pde/lecture-12/experiment.py
"""

import cmath
import json
import math
from pathlib import Path

OUT = {}


def central_matrix_apply(U, a, h):
    n = len(U)
    return [-a / (2 * h) * (U[(j + 1) % n] - U[(j - 1) % n]) for j in range(n)]


# ---------------- Part 1: skew-symmetric matrix, imaginary eigenvalues ----------------
M = 40
h = 1 / M
a = 1.0
# dense A to check A^T = -A
A = [[0.0] * M for _ in range(M)]
for j in range(M):
    A[j][(j + 1) % M] = -a / (2 * h)
    A[j][(j - 1) % M] = a / (2 * h)
skew = max(abs(A[i][j] + A[j][i]) for i in range(M) for j in range(M))
res = 0.0
max_real = 0.0
lams = []
for p in range(M):
    lam_formula = -1j * a / h * math.sin(2 * math.pi * p * h)
    u = [cmath.exp(2j * math.pi * p * j * h) for j in range(M)]
    Au = [sum(A[j][l] * u[l] for l in range(M)) for j in range(M)]
    # Rayleigh quotient gives the eigenvalue without the formula
    lam_rq = sum(Au[j] * u[j].conjugate() for j in range(M)) / M
    res = max(res, max(abs(Au[j] - lam_formula * u[j]) for j in range(M)) / (a / h))
    max_real = max(max_real, abs(lam_rq.real) / (a / h))
    lams.append(lam_rq.imag)
OUT["part1"] = {
    "M": M,
    "skew_error": skew,
    "eig_residual_rel": res,
    "max_real_part_rel": max_real,
    "imag_min_times_h": min(lams) * h,
    "imag_max_times_h": max(lams) * h,
}
# Forward Euler amplification on these eigenvalues
fe = {}
for nu in (0.1, 0.5, 0.9):
    k = nu * h / a
    fe[str(nu)] = max(abs(1 + k * (-1j * a / h * math.sin(2 * math.pi * p * h))) for p in range(M))
OUT["part1"]["fe_max_abs"] = fe

# ---------------- Part 2: LF = FTCS + diffusion; ellipse ----------------
M2 = 40
h2 = 1 / M2
nu = 0.6
k2 = nu * h2 / a
eps = h2 * h2 / (2 * k2)
U = [math.sin(2 * math.pi * j * h2) + 0.3 * math.cos(6 * math.pi * j * h2) ** 3 for j in range(M2)]
lf = [0.5 * (U[(j - 1) % M2] + U[(j + 1) % M2]) - a * k2 / (2 * h2) * (U[(j + 1) % M2] - U[(j - 1) % M2]) for j in range(M2)]
ftcs_diff = [
    U[j]
    - a * k2 / (2 * h2) * (U[(j + 1) % M2] - U[(j - 1) % M2])
    + k2 * eps / h2**2 * (U[(j - 1) % M2] - 2 * U[j] + U[(j + 1) % M2])
    for j in range(M2)
]
OUT["part2"] = {"M": M2, "nu": nu, "lf_vs_ftcs_diffusion": max(abs(x - y) for x, y in zip(lf, ftcs_diff))}
rows = []
for name, c in (("FTCS", 0.0), ("", 0.15), ("", 0.30), ("Lax–Wendroff", nu * nu), ("upwind", nu), ("", 0.8), ("Lax–Friedrichs", 1.0), ("", 1.1), ("", 1.3)):
    outside = 0
    worst = 0.0
    for p in range(M2):
        th = 2 * math.pi * p * h2
        z = complex(-c + c * math.cos(th), -nu * math.sin(th))
        g = abs(1 + z)
        worst = max(worst, g)
        if g > 1 + 1e-12:
            outside += 1
    rows.append({"name": name, "c": c, "outside": outside, "max_abs": worst})
OUT["part2"]["ellipse"] = rows

# ---------------- Part 3: schemes over one full trip ----------------


def step_lf(U, nu):
    n = len(U)
    return [0.5 * (U[(j - 1) % n] + U[(j + 1) % n]) - nu / 2 * (U[(j + 1) % n] - U[(j - 1) % n]) for j in range(n)]


def step_up(U, nu):
    n = len(U)
    return [U[j] - nu * (U[j] - U[(j - 1) % n]) for j in range(n)]


def step_lw(U, nu):
    n = len(U)
    return [
        U[j] - nu / 2 * (U[(j + 1) % n] - U[(j - 1) % n]) + nu * nu / 2 * (U[(j - 1) % n] - 2 * U[j] + U[(j + 1) % n])
        for j in range(n)
    ]


def run(scheme, U0, nu, steps):
    if scheme == "leapfrog":
        prev = U0[:]
        cur = step_lw(U0, nu)  # one Lax-Wendroff step to start (second order)
        for _ in range(steps - 1):
            n = len(cur)
            nxt = [prev[j] - nu * (cur[(j + 1) % n] - cur[(j - 1) % n]) for j in range(n)]
            prev, cur = cur, nxt
        return cur
    f = {"Lax–Friedrichs": step_lf, "upwind": step_up, "Lax–Wendroff": step_lw}[scheme]
    U = U0[:]
    for _ in range(steps):
        U = f(U, nu)
    return U


def bump_box(x):
    bump = math.exp(-300 * (x - 0.25) ** 2)
    box = 1.0 if 0.55 <= x < 0.8 else 0.0
    return bump + box


SCHEMES = ["Lax–Friedrichs", "upwind", "Lax–Wendroff", "leapfrog"]
M3 = 120
h3 = 1 / M3
nu3 = 0.5
steps3 = round(1 / (nu3 * h3))  # one full trip: T = 1
xs = [j * h3 for j in range(M3)]
U0 = [bump_box(x) for x in xs]
p3 = {"M": M3, "nu": nu3, "steps": steps3, "x": xs, "exact": U0, "runs": {}}
bump_idx = [j for j, x in enumerate(xs) if x < 0.5]
box_idx = [j for j, x in enumerate(xs) if x >= 0.5]
for s in SCHEMES:
    U = run(s, U0, nu3, steps3)
    p3["runs"][s] = {
        "U": U,
        "bump_peak": max(U[j] for j in bump_idx),
        "box_max": max(U[j] for j in box_idx),
        "box_min": min(U[j] for j in box_idx),
        "global_min": min(U),
        "err_max": max(abs(u - e) for u, e in zip(U, U0)),
        "mass": sum(U) * h3,
    }
p3["mass_exact"] = sum(U0) * h3

# order of convergence on a smooth wave, one trip, nu = 0.5
conv = {s: [] for s in SCHEMES}
for Mc in (25, 50, 100, 200, 400):
    hc = 1 / Mc
    stepsc = round(1 / (nu3 * hc))
    Uc0 = [math.sin(2 * math.pi * j * hc) for j in range(Mc)]
    for s in SCHEMES:
        U = run(s, Uc0, nu3, stepsc)
        conv[s].append({"M": Mc, "err": max(abs(u - e) for u, e in zip(U, Uc0))})
for s in SCHEMES:
    r = conv[s]
    for i in range(1, len(r)):
        r[i]["order"] = math.log2(r[i - 1]["err"] / r[i]["err"])
p3["convergence"] = conv
OUT["part3"] = p3

# ---------------- Part 4: two ways to fail (or not) ----------------
p4 = {}
# (a) FTCS with k = h^2: stable in the Lax-Richtmyer sense, bound e^{a^2 T/2}
ftcs = []
for Mf in (20, 40, 80):
    hf = 1 / Mf
    kf = hf * hf
    nuf = a * kf / hf
    steps = round(1 / kf)
    U = [math.sin(2 * math.pi * j * hf) for j in range(Mf)]
    for _ in range(steps):
        n = len(U)
        U = [U[j] - nuf / 2 * (U[(j + 1) % n] - U[(j - 1) % n]) for j in range(n)]
    gmax = max(abs(1 - 1j * nuf * math.sin(2 * math.pi * p * hf)) for p in range(Mf))
    ftcs.append({"M": Mf, "steps": steps, "max_U": max(abs(u) for u in U), "gmax_pow": gmax**steps, "bound": math.exp(a * a / 2)})
p4["ftcs_k_h2"] = ftcs
# and FTCS with fixed nu = 0.5 for comparison: one full trip from the bump and box
fixed = []
for Mf in (20, 40, 80):
    hf = 1 / Mf
    steps = round(1 / (0.5 * hf))
    U = [bump_box(j * hf) for j in range(Mf)]
    for _ in range(steps):
        U = [U[j] - 0.25 * (U[(j + 1) % Mf] - U[(j - 1) % Mf]) for j in range(Mf)]
    gmax = max(abs(1 - 0.5j * math.sin(2 * math.pi * p * hf)) for p in range(Mf))
    fixed.append({"M": Mf, "steps": steps, "max_U": max(abs(u) for u in U), "gmax_pow": gmax**steps})
p4["ftcs_nu05"] = fixed
# (b) upwind formula (U_j - U_{j-1}) used with a = -1: nu = -0.5 -> eps = a h / 2 < 0
M_wrong = 40
h_wrong = 1 / M_wrong
nuw = -0.5
Uw = [math.exp(-100 * (j * h_wrong - 0.5) ** 2) for j in range(M_wrong)]
wrong = []
for n in range(1, 81):
    Uw = [Uw[j] - nuw * (Uw[j] - Uw[(j - 1) % M_wrong]) for j in range(M_wrong)]
    if n % 20 == 0:
        wrong.append({"n": n, "max_U": max(abs(u) for u in Uw)})
p4["wrong_upwind"] = wrong
p4["wrong_upwind_gmax"] = max(abs(1 - nuw * (1 - cmath.exp(-2j * math.pi * p * h_wrong))) for p in range(M_wrong))
OUT["part4"] = p4

Path(__file__).with_name("experiment.json").write_text(json.dumps(OUT, indent=1), encoding="utf-8")

# short report
print(json.dumps(OUT["part1"], indent=1))
print(OUT["part2"]["lf_vs_ftcs_diffusion"])
for r in OUT["part2"]["ellipse"]:
    print(r)
for s in SCHEMES:
    d = p3["runs"][s]
    print(s, {k: round(v, 4) for k, v in d.items() if k != "U"})
print("mass exact", p3["mass_exact"])
for s in SCHEMES:
    print(s, [(r["M"], f'{r["err"]:.3e}', round(r.get("order", 0), 3)) for r in conv[s]])
print(json.dumps(p4, indent=1))
