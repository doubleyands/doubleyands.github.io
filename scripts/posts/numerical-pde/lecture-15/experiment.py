"""Experiment for the Lecture 15 post: Godunov's method (REA) with piecewise linear slopes.

Model problem q_t + u q_x = 0 with u = 1 on the periodic domain [0, 1), cell averages Q_i.
One step of the slope method is
    Q_i <- Q_i - nu (Q_i - Q_{i-1}) - (nu/2)(1 - nu) dx (s_i - s_{i-1}),   nu = u dt / dx,
with slope s_i from one of: zero (upwind), centered (Fromm), upwind (Beam-Warming),
downwind (Lax-Wendroff), minmod.

Part 1: the formula equals a literal REA step (reconstruct, shift exactly, average), and the
        downwind / upwind slopes give the known Lax-Wendroff and Beam-Warming formulas.
Part 2: one step from a single jump (1 -> 0), nu = 0.25.
Part 3: convergence on a smooth wave, one period.
Part 4: a smooth bump plus a box, one period: max, min, total variation, L1 error.

Standard library only. Run: python scripts/posts/numerical-pde/lecture-15/experiment.py
Writes experiment.json next to this file.
"""

import json
import math
from pathlib import Path


def minmod(a, b):
    if a * b <= 0:
        return 0.0
    return a if abs(a) < abs(b) else b


def slopes(Q, dx, kind):
    m = len(Q)
    s = [0.0] * m
    if kind == "upwind":
        return s
    for i in range(m):
        back = (Q[i] - Q[i - 1]) / dx
        fwd = (Q[(i + 1) % m] - Q[i]) / dx
        if kind == "fromm":
            s[i] = 0.5 * (back + fwd)
        elif kind == "beam-warming":
            s[i] = back
        elif kind == "lax-wendroff":
            s[i] = fwd
        elif kind == "minmod":
            s[i] = minmod(back, fwd)
        else:
            raise ValueError(kind)
    return s


def step(Q, nu, dx, kind):
    s = slopes(Q, dx, kind)
    c = 0.5 * nu * (1 - nu) * dx
    return [Q[i] - nu * (Q[i] - Q[i - 1]) - c * (s[i] - s[i - 1]) for i in range(len(Q))]


def rea_step(Q, nu, dx, kind):
    """Literal REA: line in each cell, shift right by nu*dx, average over each cell."""
    s = slopes(Q, dx, kind)
    m = len(Q)
    shift = nu * dx
    out = []
    for i in range(m):
        # the new cell i covers old points [x_{i-1/2} - shift, x_{i+1/2} - shift] (relative to cell centers)
        # part from cell i-1: local coordinate xi in [dx/2 - shift, dx/2]
        a, b = dx / 2 - shift, dx / 2
        part_left = Q[i - 1] * (b - a) + s[i - 1] * (b * b - a * a) / 2
        # part from cell i: local coordinate in [-dx/2, dx/2 - shift]
        a, b = -dx / 2, dx / 2 - shift
        part_mid = Q[i] * (b - a) + s[i] * (b * b - a * a) / 2
        out.append((part_left + part_mid) / dx)
    return out


def tv(Q):
    return sum(abs(Q[i] - Q[i - 1]) for i in range(len(Q)))


def cell_avg(F, m):
    """Exact cell averages of a function with antiderivative F on m cells."""
    dx = 1.0 / m
    return [(F((i + 1) * dx) - F(i * dx)) / dx for i in range(m)]


KINDS = ["upwind", "fromm", "beam-warming", "lax-wendroff", "minmod"]
out = {}

# ---------- Part 1 ----------
m = 50
dx = 1.0 / m
Q0 = cell_avg(lambda x: -math.cos(2 * math.pi * x) / (2 * math.pi) + 0.4 * math.sin(6 * math.pi * x) / (6 * math.pi)
              + (x - 0.3 if x > 0.3 else 0.0) - (x - 0.55 if x > 0.55 else 0.0), m)
p1 = {"cells": m, "rea_vs_formula": {}}
for nu in (0.3, 0.7):
    for kind in KINDS:
        a = step(Q0, nu, dx, kind)
        b = rea_step(Q0, nu, dx, kind)
        p1["rea_vs_formula"][f"{kind} nu={nu}"] = max(abs(x - y) for x, y in zip(a, b))
nu = 0.7
lw = [Q0[i] - nu / 2 * (Q0[(i + 1) % m] - Q0[i - 1]) + nu * nu / 2 * (Q0[(i + 1) % m] - 2 * Q0[i] + Q0[i - 1]) for i in range(m)]
bw = [Q0[i] - nu / 2 * (3 * Q0[i] - 4 * Q0[i - 1] + Q0[i - 2]) + nu * nu / 2 * (Q0[i] - 2 * Q0[i - 1] + Q0[i - 2]) for i in range(m)]
p1["lw_diff"] = max(abs(x - y) for x, y in zip(step(Q0, nu, dx, "lax-wendroff"), lw))
p1["bw_diff"] = max(abs(x - y) for x, y in zip(step(Q0, nu, dx, "beam-warming"), bw))
p1["mass_change"] = {k: abs(sum(step(Q0, nu, dx, k)) - sum(Q0)) * dx for k in KINDS}
out["part1"] = p1

# ---------- Part 2: one step from a jump ----------
m = 12
dx = 1.0 / m
Qj = [1.0] * 6 + [0.0] * 6   # periodic, so there is also a 0 -> 1 jump at x = 0
nu = 0.25
p2 = {"nu": nu, "cells": {}, "tv_before": tv(Qj)}
# exact: shift by a quarter cell; the cell right of the 1->0 jump gets 0.25
for kind in KINDS:
    Qn = step(Qj, nu, dx, kind)
    p2["cells"][kind] = {"i-1": Qn[4], "i": Qn[5], "i+1": Qn[6], "i+2": Qn[7], "tv": tv(Qn),
                         "max": max(Qn), "min": min(Qn)}
out["part2"] = p2

# ---------- Part 3: convergence on a smooth wave ----------
p3 = {"nu": 0.6, "grids": [20, 40, 80, 160, 320], "l1": {}, "max": {}}
for kind in KINDS:
    p3["l1"][kind] = []
    p3["max"][kind] = []
    for m in p3["grids"]:
        dx = 1.0 / m
        F = lambda x: -math.cos(2 * math.pi * x) / (2 * math.pi)
        Q = cell_avg(F, m)
        exact = list(Q)
        nsteps = round(1.0 / (0.6 * dx))
        nu = 1.0 / (nsteps * dx)  # land exactly on t = 1
        for _ in range(nsteps):
            Q = step(Q, nu, dx, kind)
        err = [abs(a - b) for a, b in zip(Q, exact)]
        p3["l1"][kind].append(sum(err) * dx)
        p3["max"][kind].append(max(err))
    p3["l1_order"] = p3.get("l1_order", {})
    p3["l1_order"][kind] = [math.log2(p3["l1"][kind][j] / p3["l1"][kind][j + 1]) for j in range(len(p3["grids"]) - 1)]
    p3["max_order"] = p3.get("max_order", {})
    p3["max_order"][kind] = [math.log2(p3["max"][kind][j] / p3["max"][kind][j + 1]) for j in range(len(p3["grids"]) - 1)]
out["part3"] = p3


# ---------- Part 4: bump + box, one period ----------
def F_shape(x):
    # antiderivative of: cosine bump (height 1, center 0.25, half-width 0.12) + box 1 on [0.55, 0.8)
    c, w = 0.25, 0.12
    if x <= c - w:
        b = 0.0
    elif x >= c + w:
        b = w
    else:
        b = 0.5 * (x - (c - w)) + w / (2 * math.pi) * math.sin(math.pi * (x - c) / w)
    box = min(max(x - 0.55, 0.0), 0.25)
    return b + box


m = 200
dx = 1.0 / m
nu = 0.6
nsteps = round(1.0 / (nu * dx))
nu = 1.0 / (nsteps * dx)
Qinit = cell_avg(F_shape, m)
p4 = {"cells": m, "nu": nu, "steps": nsteps, "tv0": tv(Qinit), "x": [(i + 0.5) * dx for i in range(m)],
      "exact": Qinit, "runs": {}}
for kind in KINDS:
    Q = list(Qinit)
    tvs = [tv(Q)]
    tv_up = 0  # steps where TV went up by more than rounding
    for _ in range(nsteps):
        Qn = step(Q, nu, dx, kind)
        if tv(Qn) > tv(Q) + 1e-12:
            tv_up += 1
        Q = Qn
        tvs.append(tv(Q))
    bump = [Q[i] for i in range(m) if (i + 0.5) * dx < 0.5]
    boxr = [Q[i] for i in range(m) if (i + 0.5) * dx >= 0.5]
    p4["runs"][kind] = {
        "profile": Q,
        "max": max(Q), "min": min(Q),
        "bump_peak": max(bump), "box_max": max(boxr), "box_min": min(boxr),
        "tv_end": tvs[-1], "tv_max": max(tvs), "tv_up_steps": tv_up,
        "l1": sum(abs(a - b) for a, b in zip(Q, Qinit)) * dx,
        "mass_change": abs(sum(Q) - sum(Qinit)) * dx,
    }
out["part4"] = p4

Path(__file__).with_name("experiment.json").write_text(json.dumps(out, indent=1), encoding="utf-8")

# ---------- print summary ----------
print("Part 1", {k: f"{v:.1e}" for k, v in p1["rea_vs_formula"].items()})
print("  LW diff", p1["lw_diff"], "BW diff", p1["bw_diff"], "mass", {k: f"{v:.1e}" for k, v in p1["mass_change"].items()})
print("Part 2 tv before", p2["tv_before"])
for k, v in p2["cells"].items():
    print(f"  {k:13s}", {a: round(b, 6) for a, b in v.items()})
print("Part 3")
for k in KINDS:
    print(f"  {k:13s} L1", [f"{e:.2e}" for e in p3["l1"][k]], [f"{o:.2f}" for o in p3["l1_order"][k]])
    print(f"  {'':13s} max", [f"{e:.2e}" for e in p3["max"][k]], [f"{o:.2f}" for o in p3["max_order"][k]])
print("Part 4 nu", nu, "steps", nsteps, "tv0", p4["tv0"])
for k in KINDS:
    r = p4["runs"][k]
    print(f"  {k:13s}", {a: (round(b, 4) if isinstance(b, float) else b) for a, b in r.items() if a != "profile"})
