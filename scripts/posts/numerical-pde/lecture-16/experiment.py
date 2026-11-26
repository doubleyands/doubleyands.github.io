"""Experiment for the Lecture 16 post: MC and superbee limiters, the flux-limiter form, Harten's theorem.

Model problem q_t + q_x = 0 (speed a = 1) on the periodic domain [0, 1), cell averages Q_i.
One step of the slope method (Lecture 15) is
    Q_i <- Q_i - nu (Q_i - Q_{i-1}) - (nu/2)(1 - nu) dx (s_i - s_{i-1}),   nu = dt / dx.

Part 1: the two ways to write MC agree, and the slope form equals the flux-limiter form
        written with theta and phi(theta).
Part 2: Harten's coefficients C_{i-1} for random data (limiters vs Lax-Wendroff).
Part 3: convergence on a smooth wave, one period.
Part 4: a smooth bump plus a box, one and five periods.

Standard library only. Run: python scripts/posts/numerical-pde/lecture-16/experiment.py
Writes experiment.json next to this file.
"""

import json
import math
import random
from pathlib import Path


def minmod(*args):
    if all(a > 0 for a in args):
        return min(args)
    if all(a < 0 for a in args):
        return max(args)
    return 0.0


def maxmod(a, b):
    if a * b <= 0:
        return 0.0
    return a if abs(a) > abs(b) else b


def slope(back, fwd, kind):
    """Slope in a cell from the two one-sided slopes (left = back, right = fwd)."""
    if kind == "upwind":
        return 0.0
    if kind == "fromm":
        return 0.5 * (back + fwd)
    if kind == "beam-warming":
        return back
    if kind == "lax-wendroff":
        return fwd
    if kind == "minmod":
        return minmod(back, fwd)
    if kind == "mc":
        return minmod(0.5 * (back + fwd), minmod(2 * back, 2 * fwd))
    if kind == "superbee":
        return maxmod(minmod(fwd, 2 * back), minmod(2 * fwd, back))
    raise ValueError(kind)


def slopes(Q, dx, kind):
    m = len(Q)
    return [slope((Q[i] - Q[i - 1]) / dx, (Q[(i + 1) % m] - Q[i]) / dx, kind) for i in range(m)]


def step(Q, nu, dx, kind):
    s = slopes(Q, dx, kind)
    c = 0.5 * nu * (1 - nu) * dx
    return [Q[i] - nu * (Q[i] - Q[i - 1]) - c * (s[i] - s[i - 1]) for i in range(len(Q))]


def phi_of(theta, kind):
    """phi(theta) of the same method: the limited jump divided by the plain jump.
    Computed from the slope rule with back = theta, fwd = 1 (slope rules are homogeneous)."""
    return slope(theta, 1.0, kind)


def step_flux(Q, nu, kind):
    """Flux-limiter form: F_{i-1/2} = a Q_{i-1} + (a/2)(1 - nu) phi(theta_{i-1/2}) dQ_{i-1/2}, a = 1, in units of dx."""
    m = len(Q)
    dQ = [Q[i] - Q[i - 1] for i in range(m)]  # dQ[i] = Q_i - Q_{i-1} = jump at i-1/2
    delta = []
    for i in range(m):
        # jump at the wall i+1/2 is dQ[i+1]; upwind neighbour jump is dQ[i]
        d, dup = dQ[(i + 1) % m], dQ[i]
        if d == 0.0:
            delta.append(slope(dup, 0.0, kind))  # limit of phi(theta) * d as d -> 0
        else:
            delta.append(phi_of(dup / d, kind) * d)
    # delta[i] = limited jump used by cell i (= dx * sigma_i)
    flux = [Q[i] + 0.5 * (1 - nu) * delta[i] for i in range(m)]  # flux at i+1/2 (times dt/dx below)
    return [Q[i] - nu * (flux[i] - flux[i - 1]) for i in range(m)]


def tv(Q):
    return sum(abs(Q[i] - Q[i - 1]) for i in range(len(Q)))


def cell_avg(F, m):
    dx = 1.0 / m
    return [(F((i + 1) * dx) - F(i * dx)) / dx for i in range(m)]


LIMITERS = ["minmod", "mc", "superbee"]
KINDS = ["upwind", "fromm", "beam-warming", "lax-wendroff"] + LIMITERS
rng = random.Random(16)
out = {}

# ---------- Part 1 ----------
worst = 0.0
for _ in range(200000):
    a, b, c = (rng.uniform(-1, 1) for _ in range(3))
    worst = max(worst, abs(minmod(a, minmod(b, c)) - minmod(a, b, c)))
p1 = {"mc_nested_vs_three": worst, "slope_vs_flux": {}}
m = 64
dx = 1.0 / m
for trial in range(3):
    Q = [rng.uniform(-1, 1) for _ in range(m)]
    for nu in (0.3, 0.8):
        for kind in KINDS:
            d = max(abs(x - y) for x, y in zip(step(Q, nu, dx, kind), step_flux(Q, nu, kind)))
            p1["slope_vs_flux"][kind] = max(p1["slope_vs_flux"].get(kind, 0.0), d)
out["part1"] = p1

# ---------- Part 2: Harten coefficients ----------
# For a > 0 the flux-limiter update is Q_i - C_{i-1} (Q_i - Q_{i-1}) with D_i = 0 and
# C_{i-1} = nu [1 + (1/2)(1 - nu)(phi(theta_{i+1/2}) / theta_{i+1/2} - phi(theta_{i-1/2}))].
nu = 0.7
p2 = {"nu": nu, "C_range": {}, "C_outside_01": {}}
thetas = [rng.uniform(-3, 3) for _ in range(200000)]
for kind in ["upwind", "lax-wendroff", "beam-warming", "fromm"] + LIMITERS:
    lo, hi, bad, n = math.inf, -math.inf, 0, 0
    for j in range(0, len(thetas) - 1, 2):
        t1, t2 = thetas[j], thetas[j + 1]  # t1 = theta_{i+1/2}, t2 = theta_{i-1/2}
        if abs(t1) < 1e-3:
            continue
        C = nu * (1 + 0.5 * (1 - nu) * (phi_of(t1, kind) / t1 - phi_of(t2, kind)))
        lo, hi = min(lo, C), max(hi, C)
        n += 1
        bad += (C < 0 or C > 1)
    p2["C_range"][kind] = [lo, hi]
    p2["C_outside_01"][kind] = bad / n
out["part2"] = p2

# ---------- Part 3: smooth wave ----------
p3 = {"grids": [20, 40, 80, 160, 320], "l1": {}, "max": {}, "l1_order": {}, "max_order": {}}
for kind in ["upwind", "fromm"] + LIMITERS:
    p3["l1"][kind], p3["max"][kind] = [], []
    for m in p3["grids"]:
        dx = 1.0 / m
        Q = cell_avg(lambda x: -math.cos(2 * math.pi * x) / (2 * math.pi), m)
        exact = list(Q)
        nsteps = round(1.0 / (0.6 * dx))
        nu = 1.0 / (nsteps * dx)
        for _ in range(nsteps):
            Q = step(Q, nu, dx, kind)
        err = [abs(a - b) for a, b in zip(Q, exact)]
        p3["l1"][kind].append(sum(err) * dx)
        p3["max"][kind].append(max(err))
    g = len(p3["grids"]) - 1
    p3["l1_order"][kind] = [math.log2(p3["l1"][kind][j] / p3["l1"][kind][j + 1]) for j in range(g)]
    p3["max_order"][kind] = [math.log2(p3["max"][kind][j] / p3["max"][kind][j + 1]) for j in range(g)]
out["part3"] = p3


# ---------- Part 4: bump + box ----------
def F_shape(x):
    c, w = 0.25, 0.12
    if x <= c - w:
        b = 0.0
    elif x >= c + w:
        b = w
    else:
        b = 0.5 * (x - (c - w)) + w / (2 * math.pi) * math.sin(math.pi * (x - c) / w)
    return b + min(max(x - 0.55, 0.0), 0.25)


m = 200
dx = 1.0 / m
nsteps = round(1.0 / (0.6 * dx))
nu = 1.0 / (nsteps * dx)
Qinit = cell_avg(F_shape, m)
xs = [(i + 0.5) * dx for i in range(m)]
p4 = {"cells": m, "nu": nu, "steps_per_period": nsteps, "tv0": tv(Qinit), "x": xs, "exact": Qinit, "runs": {}}
for kind in ["upwind", "fromm"] + LIMITERS:
    Q = list(Qinit)
    tv_up = 0
    res = {}
    for period in range(1, 6):
        for _ in range(nsteps):
            Qn = step(Q, nu, dx, kind)
            if tv(Qn) > tv(Q) + 1e-12:
                tv_up += 1
            Q = Qn
        if period in (1, 5):
            bump = [Q[i] for i in range(m) if xs[i] < 0.5]
            top = max(bump)
            box = [i for i in range(m) if xs[i] >= 0.5]
            res[f"p{period}"] = {
                "profile": Q,
                "max": max(Q), "min": min(Q), "bump_peak": top,
                # cells near the bump top within 0.02 of the peak: a flat top shows up here
                "flat_top_cells": sum(1 for v in bump if v > top - 0.02),
                # cells in the box region that are neither near 0 nor near 1: width of the smeared edges
                "edge_cells": sum(1 for i in box if 0.05 < Q[i] < 0.95),
                "tv": tv(Q), "tv_up_steps": tv_up,
                "l1": sum(abs(a - b) for a, b in zip(Q, Qinit)) * dx,
                "l1_bump": sum(abs(Q[i] - Qinit[i]) for i in range(m) if xs[i] < 0.5) * dx,
                "l1_box": sum(abs(Q[i] - Qinit[i]) for i in range(m) if xs[i] >= 0.5) * dx,
                "top_cells_0005": sum(1 for v in bump if v > top - 0.005),
                "mass_change": abs(sum(Q) - sum(Qinit)) * dx,
            }
    p4["runs"][kind] = res
bump0 = [Qinit[i] for i in range(m) if xs[i] < 0.5]
p4["exact_flat_top_cells"] = sum(1 for v in bump0 if v > max(bump0) - 0.02)
p4["exact_top_cells_0005"] = sum(1 for v in bump0 if v > max(bump0) - 0.005)
p4["exact_edge_cells"] = sum(1 for i in range(m) if xs[i] >= 0.5 and 0.05 < Qinit[i] < 0.95)
out["part4"] = p4

Path(__file__).with_name("experiment.json").write_text(json.dumps(out, indent=1), encoding="utf-8")

print("Part 1", p1["mc_nested_vs_three"], {k: f"{v:.1e}" for k, v in p1["slope_vs_flux"].items()})
print("Part 2", {k: [round(v[0], 3), round(v[1], 3)] for k, v in p2["C_range"].items()}, {k: round(v, 3) for k, v in p2["C_outside_01"].items()})
print("Part 3")
for k in p3["l1"]:
    print(f"  {k:9s} L1", [f"{e:.2e}" for e in p3["l1"][k]], [f"{o:.2f}" for o in p3["l1_order"][k]])
    print(f"  {'':9s} max", [f"{e:.2e}" for e in p3["max"][k]], [f"{o:.2f}" for o in p3["max_order"][k]])
print("Part 4 nu", nu, "steps", nsteps, "tv0", p4["tv0"], "exact flat", p4["exact_flat_top_cells"], "edge", p4["exact_edge_cells"])
for k, r in p4["runs"].items():
    for p, v in r.items():
        print(f"  {k:9s} {p}", {a: (round(b, 4) if isinstance(b, float) else b) for a, b in v.items() if a != "profile"})
