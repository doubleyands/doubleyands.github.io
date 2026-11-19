"""Experiment for the Lecture 14 post: characteristics, the CFL condition, and numerical fluxes.

Python standard library only (math). Periodic grids on [0, 1], a = 1 unless stated.

Part 1: the CFL condition with upwind. A cosine bump moves with Courant number 1, 0.5 and 1.05.
        Also the domain of dependence: change the initial data near the foot of the characteristic
        and see whether the numerical value at (X, T) notices.
Part 2: four numerical fluxes in conservation form (central, Lax-Friedrichs, Richtmyer, upwind).
        Checks that they reduce to FTCS, LF, LW and upwind for f(q) = a q, and runs them.
Part 3: upwind as a convex combination: the maximum principle for 0 <= nu <= 1, and its loss at nu = 1.05.
Part 4 (my own addition): a hyperbolic system (real eigenvalues) and a non-hyperbolic one
        (complex eigenvalues), both with Lax-Friedrichs.

Run: python scripts/posts/numerical-pde/lecture-14/experiment.py   (writes experiment.json next to this file)
"""

import json
import math
from pathlib import Path

PI = math.pi
out = {}


def bump(x, c=0.2, w=0.1):
    """Cosine bump of height 1 and half-width w, centered at c, periodic on [0, 1]."""
    d = (x - c + 0.5) % 1.0 - 0.5
    return 0.5 * (1 + math.cos(PI * d / w)) if abs(d) < w else 0.0


def cells(m):
    h = 1.0 / m
    return h, [(i + 0.5) * h for i in range(m)]


# ---------- numerical fluxes, conservation form ----------
# F[i] is the flux through the left wall x_{i-1/2} of cell i, built from Q[i-1] and Q[i].

def flux_central(Q, f, r):
    return [0.5 * (f(Q[i - 1]) + f(Q[i])) for i in range(len(Q))]


def flux_lf(Q, f, r):
    # r = dt/dx
    return [0.5 * (f(Q[i - 1]) + f(Q[i])) - 0.5 / r * (Q[i] - Q[i - 1]) for i in range(len(Q))]


def flux_richtmyer(Q, f, r):
    F = []
    for i in range(len(Q)):
        half = 0.5 * (Q[i - 1] + Q[i]) - 0.5 * r * (f(Q[i]) - f(Q[i - 1]))
        F.append(f(half))
    return F


def flux_upwind(Q, f, r, a=1.0):
    return [a * Q[i - 1] if a > 0 else a * Q[i] for i in range(len(Q))]


def fv_step(Q, flux, f, r):
    F = flux(Q, f, r)
    m = len(Q)
    return [Q[i] - r * (F[(i + 1) % m] - F[i]) for i in range(m)]


FLUXES = {"central": flux_central, "Lax–Friedrichs": flux_lf, "Richtmyer": flux_richtmyer, "upwind": flux_upwind}


# node forms we already know (Lectures 11-13), for the identity checks
def ftcs(U, nu):
    m = len(U)
    return [U[j] - 0.5 * nu * (U[(j + 1) % m] - U[j - 1]) for j in range(m)]


def lf(U, nu):
    m = len(U)
    return [0.5 * (U[j - 1] + U[(j + 1) % m]) - 0.5 * nu * (U[(j + 1) % m] - U[j - 1]) for j in range(m)]


def lw(U, nu):
    m = len(U)
    return [U[j] - 0.5 * nu * (U[(j + 1) % m] - U[j - 1]) + 0.5 * nu * nu * (U[j - 1] - 2 * U[j] + U[(j + 1) % m])
            for j in range(m)]


def upwind(U, nu):
    return [nu * U[j - 1] + (1 - nu) * U[j] for j in range(len(U))]


def maxdiff(A, B):
    return max(abs(a - b) for a, b in zip(A, B))


# ---------- Part 1: CFL with upwind ----------
p1 = {}
M = 100
h, xs = cells(M)
U0 = [bump(x) for x in xs]
runs = {}
T_list = [0.5, 0.7, 1.0, 2.0, 3.0]
for nu in (1.0, 0.5, 1.05):
    k = nu * h
    U = U0[:]
    t, n = 0.0, 0
    rec = []
    snaps = {}
    for T in T_list:
        nT = round(T / k)
        while n < nT:
            U = upwind(U, nu)
            n += 1
        t = n * k
        exact = [bump(x - t) for x in xs]
        rec.append({"T": T, "steps": n, "t": t, "max": max(U), "min": min(U),
                    "err": maxdiff(U, exact)})
        snaps[str(T)] = U[:]
    runs[str(nu)] = {"k": k, "records": rec, "profile": snaps["0.7"]}
p1["upwind_runs"] = runs
p1["x"] = xs
p1["profile_T"] = 0.7
p1["profile_t"] = {nu: runs[nu]["records"][1]["t"] for nu in runs}  # nu = 1.05 lands at round(T/k)*k
p1["exact_profile"] = [bump(x - 0.7) for x in xs]

# domain of dependence: perturb the initial data near the foot X - T, look at U(X, T)
dod = []
X, T = 0.8, 0.4
foot = X - T
for m in (50, 100, 200, 400):
    hm, xm = cells(m)
    j = min(range(m), key=lambda i: abs(xm[i] - X))
    for nu in (0.8, 1.25):
        k = nu * hm
        n = round(T / k)
        base = [0.0] * m
        pert = [bump(x, c=foot + (xm[j] - X), w=0.025) for x in xm]  # bump at the foot of the grid point x_j
        Ub, Up = base[:], pert[:]
        for _ in range(n):
            Ub = upwind(Ub, nu)
            Up = upwind(Up, nu)
        # leftmost initial point that U_j^n depends on: x_{j-n}
        dod.append({"cells": m, "nu": nu, "steps": n, "change_at_X": Up[j] - Ub[j],
                    "numerical_dod_left": xm[j] - n * hm, "foot": xm[j] - n * k,
                    "exact_change": 1.0})
p1["domain_of_dependence"] = dod
out["part1"] = p1

# ---------- Part 2: numerical fluxes ----------
p2 = {}
a = 1.0
lin = lambda q: a * q
M = 100
h, xs = cells(M)
nu = 0.5
Q0 = [math.sin(2 * PI * x) + 0.3 * math.cos(6 * PI * x) for x in xs]
checks = {}
for name, node in (("central", ftcs), ("Lax–Friedrichs", lf), ("Richtmyer", lw), ("upwind", upwind)):
    checks[name] = maxdiff(fv_step(Q0, FLUXES[name], lin, nu), node(Q0, nu))
p2["identity_check_one_step"] = checks

# run each flux on the bump, one lap (200 steps at nu = 0.5), record max|Q|, error, mass
U0 = [bump(x) for x in xs]
mass0 = h * sum(U0)
laps = {}
for nu in (0.5, 0.1):
    res = {}
    nsteps = round(1.0 / (nu * h))
    for name in FLUXES:
        Q = U0[:]
        rec = []
        for n in range(1, nsteps + 1):
            Q = fv_step(Q, FLUXES[name], lin, nu)
            if n in (nsteps // 4, nsteps // 2, nsteps):
                rec.append({"step": n, "t": n * nu * h, "maxabs": max(abs(q) for q in Q)})
        res[name] = {"records": rec, "err_one_lap": maxdiff(Q, U0), "mass_change": h * sum(Q) - mass0,
                     "min": min(Q), "max": max(Q)}
    laps[str(nu)] = {"steps": nsteps, "schemes": res}
p2["one_lap"] = laps

# central flux: growth factor per step on the sawtooth-ish worst mode is sqrt(1 + nu^2 sin^2)
p2["central_max_g"] = {str(nu): math.sqrt(1 + nu * nu) for nu in (0.5, 0.1)}

# LF artificial diffusion beta = dx^2/(2 dt) for fixed nu, refining
p2["lf_beta"] = [{"cells": m, "beta": (1 / m) ** 2 / (2 * 0.5 / m)} for m in (50, 100, 200, 400)]
out["part2"] = p2

# ---------- Part 3: maximum principle for upwind ----------
p3 = []
M = 100
h, xs = cells(M)
box = [1.0 if 0.2 <= x < 0.45 else 0.0 for x in xs]
for nu in (0.5, 0.9, 1.0, 1.05):
    U = box[:]
    first = upwind(U, nu)
    hi, lo = max(first), min(first)
    U = first
    ever_hi, ever_lo = hi, lo
    snaps = {}
    for n in range(2, 201):
        U = upwind(U, nu)
        ever_hi, ever_lo = max(ever_hi, max(U)), min(ever_lo, min(U))
        if n in (20, 50, 100, 200):
            snaps[str(n)] = {"max": max(U), "min": min(U)}
    p3.append({"nu": nu, "coeffs": [nu, 1 - nu], "after_1": {"max": hi, "min": lo},
               "over_200": {"max": ever_hi, "min": ever_lo}, "snapshots": snaps})
out["part3"] = p3

# ---------- Part 4: hyperbolic vs non-hyperbolic system, Lax-Friedrichs ----------
# u_t + A u_x = 0, u = (u1, u2). LF: U_j^{n+1} = (U_{j-1} + U_{j+1})/2 - (k/2h) A (U_{j+1} - U_{j-1})


def eig2(A):
    (p, q), (s, t) = A
    tr, det = p + t, p * t - q * s
    disc = tr * tr / 4 - det
    if disc >= 0:
        return [tr / 2 + math.sqrt(disc), tr / 2 - math.sqrt(disc)], True
    return [complex(tr / 2, math.sqrt(-disc)), complex(tr / 2, -math.sqrt(-disc))], False


def lf_system(U1, U2, A, r):
    m = len(U1)
    (p, q), (s, t) = A
    N1, N2 = [], []
    for j in range(m):
        jp, jm = (j + 1) % m, j - 1
        d1, d2 = U1[jp] - U1[jm], U2[jp] - U2[jm]
        N1.append(0.5 * (U1[jm] + U1[jp]) - 0.5 * r * (p * d1 + q * d2))
        N2.append(0.5 * (U2[jm] + U2[jp]) - 0.5 * r * (s * d1 + t * d2))
    return N1, N2


systems = {"A": [[0.0, 1.0], [4.0, 0.0]], "B": [[0.0, -1.0], [1.0, 0.0]]}
p4 = {}
T4 = 0.2
for key, A in systems.items():
    ev, real = eig2(A)
    smax = max(abs(e) for e in ev)
    rows = []
    for m in (25, 50, 100, 200, 400):
        h, xs = cells(m)
        r = 0.8 / smax  # dt/dx so that max |speed| * dt/dx = 0.8
        U1 = [bump(x, c=0.5, w=0.1) for x in xs]
        U2 = [0.0] * m
        n = round(T4 / (r * h))
        for _ in range(n):
            U1, U2 = lf_system(U1, U2, A, r)
        t = n * r * h
        row = {"cells": m, "steps": n, "t": t,
               "max": max(max(abs(u) for u in U1), max(abs(u) for u in U2))}
        if real:
            # A: eigenvalues +2 and -2 with eigenvectors (1, 2) and (1, -2), so the bump splits in two
            ex1 = [0.5 * (bump(x - 2 * t, 0.5, 0.1) + bump(x + 2 * t, 0.5, 0.1)) for x in xs]
            ex2 = [bump(x - 2 * t, 0.5, 0.1) - bump(x + 2 * t, 0.5, 0.1) for x in xs]
            row["err"] = max(maxdiff(U1, ex1), maxdiff(U2, ex2))
            row["max_u1"] = max(U1)
        rows.append(row)
    p4[key] = {"A": A, "eigenvalues": [str(e) for e in ev], "real": real, "runs": rows}
out["part4"] = p4

Path(__file__).with_name("experiment.json").write_text(json.dumps(out, indent=1, ensure_ascii=False), encoding="utf-8")


if __name__ == "__main__":
    import pprint
    s = json.loads(json.dumps(out))
    s["part1"].pop("x"); s["part1"].pop("exact_profile")
    for v in s["part1"]["upwind_runs"].values():
        v.pop("profile")
    pprint.pprint(s, width=140)
