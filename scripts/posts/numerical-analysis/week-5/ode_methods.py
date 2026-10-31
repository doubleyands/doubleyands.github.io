"""Experiments for Week 5: Euler, midpoint and trapezoidal methods for y' = f(x, y).

Python standard library only. Writes ode_methods.json next to this file.

A. Order of convergence on the logistic equation y' = y - y^2, y(0) = 0.1, x in [0, 4].
   Exact solution y = 1 / (1 + 9 e^(-x)). Max error over the grid, and the Euler error bound.
B. One step vs. many steps for the trapezoidal method (local O(h^3), global O(h^2)).
C. Picard iteration for the implicit trapezoidal step: measured contraction vs. h K / 2.
D. An ill-conditioned problem: y' = 50 (y - cos x) - sin x, y(0) = 1, exact y = cos x.

Run: python scripts/posts/numerical-analysis/week-5/ode_methods.py
"""

import json
import math
from pathlib import Path

A, B, Y0 = 0.0, 4.0, 0.1


def f(x, y):
    return y - y * y


def exact(x):
    return 1.0 / (1.0 + 9.0 * math.exp(-x))


def euler(n):
    h = (B - A) / n
    ys = [Y0]
    for k in range(n):
        ys.append(ys[-1] + h * f(A + k * h, ys[-1]))
    return ys


def midpoint(n):
    """Two-step method; y_1 comes from one Euler step."""
    h = (B - A) / n
    ys = [Y0, Y0 + h * f(A, Y0)]
    for k in range(1, n):
        ys.append(ys[k - 1] + 2 * h * f(A + k * h, ys[k]))
    return ys


def picard_step(x, y, h, guess, tol=1e-15, max_iter=200):
    """Solve z = y + h/2 (f(x, y) + f(x + h, z)) by simple (Picard) iteration."""
    known = y + 0.5 * h * f(x, y)
    z = guess
    history = [z]
    for it in range(1, max_iter + 1):
        z_new = known + 0.5 * h * f(x + h, z)
        history.append(z_new)
        if abs(z_new - z) <= tol * max(1.0, abs(z_new)):
            return z_new, it, history
        z = z_new
    return z, max_iter, history


def trapezoidal(n):
    h = (B - A) / n
    ys = [Y0]
    iters = 0
    for k in range(n):
        x = A + k * h
        guess = ys[-1] + h * f(x, ys[-1])  # Euler step as the first guess
        z, it, _ = picard_step(x, ys[-1], h, guess)
        ys.append(z)
        iters += it
    return ys, iters / n


def max_error(ys, n):
    h = (B - A) / n
    return max(abs(exact(A + k * h) - y) for k, y in enumerate(ys))


# ----- A: orders -----
# Euler bound: C h with C = (e^{(b-a)K} - 1) / (2K) * max|y''|, e_0 = 0.
# K = max |df/dy| = max |1 - 2y| over the range of y in [0.1, exact(4)].
K = max(abs(1 - 2 * Y0), abs(1 - 2 * exact(B)))
grid = [A + (B - A) * i / 40000 for i in range(40001)]


def ypp(x):
    y = exact(x)
    return (1 - 2 * y) * (y - y * y)  # y'' = f_y * f


M2 = max(abs(ypp(x)) for x in grid)
C_bound = (math.exp((B - A) * K) - 1) / (2 * K) * M2

part_a = []
for k in range(8):
    n = 20 * 2**k
    h = (B - A) / n
    trap, avg_it = trapezoidal(n)
    part_a.append({
        "n": n,
        "h": h,
        "euler": max_error(euler(n), n),
        "midpoint": max_error(midpoint(n), n),
        "trapezoidal": max_error(trap, n),
        "euler_bound": C_bound * h,
        "picard_iterations": avg_it,
    })
for prev, row in zip(part_a, part_a[1:]):
    for key in ("euler", "midpoint", "trapezoidal"):
        row[key + "_ratio"] = prev[key] / row[key]

# ----- B: one step vs. the whole interval (trapezoidal) -----
part_b = []
for h in (0.4, 0.2, 0.1, 0.05, 0.025):
    z, _, _ = picard_step(A, Y0, h, Y0)
    n = round((B - A) / h)
    trap, _ = trapezoidal(n)
    part_b.append({
        "h": h,
        "one_step": exact(A + h) - z,
        "global": max_error(trap, n),
    })
for prev, row in zip(part_b, part_b[1:]):
    row["one_step_ratio"] = prev["one_step"] / row["one_step"]
    row["global_ratio"] = prev["global"] / row["global"]

# ----- C: Picard contraction for one step -----
part_c = []
for h in (0.4, 0.1):
    x, y = A, Y0
    z, it, hist = picard_step(x, y, h, y)  # start from y_n
    errs = [abs(z - v) for v in hist[:6]]
    ratios = [errs[i + 1] / errs[i] for i in range(len(errs) - 1)]
    lipschitz_here = abs(1 - 2 * z)
    part_c.append({
        "h": h,
        "errors": errs,
        "ratios": ratios,
        "h_fy_over_2": h * lipschitz_here / 2,
        "iterations": it,
    })
# where Picard fails: h K / 2 > 1 on the stiff problem of part D
lam = 50.0


def g(x, y):
    return lam * (y - math.cos(x)) - math.sin(x)


def picard_g(x, y, h, n_iter):
    known = y + 0.5 * h * g(x, y)
    z = y
    out = []
    for _ in range(n_iter):
        z = known + 0.5 * h * g(x + h, z)
        out.append(z)
    return out


diverge = picard_g(0.0, 1.0, 0.1, 6)  # h * 50 / 2 = 2.5 > 1
true_next = None
# exact trapezoidal step is linear here: solve z directly
h = 0.1
x0, y0 = 0.0, 1.0
true_next = (y0 + 0.5 * h * g(x0, y0) + 0.5 * h * (-lam * math.cos(x0 + h) - math.sin(x0 + h))) / (1 - 0.5 * h * lam)
part_c_fail = {
    "h": h,
    "h_K_over_2": h * lam / 2,
    "trapezoidal_root": true_next,
    "iterates": diverge,
}

# ----- D: ill-conditioned problem -----
part_d = {"perturbation": [], "euler": []}
for eps in (1e-16, 1e-10):
    for x in (0.2, 0.5, 1.0):
        part_d["perturbation"].append({"eps": eps, "x": x, "gap": eps * math.exp(lam * x)})
for n in (1000, 10000, 100000):
    h = 1.0 / n
    y = 1.0
    row = {"n": n, "h": h}
    marks = {round(0.2 * n): "x0.2", round(0.5 * n): "x0.5", n: "x1.0"}
    for k in range(n):
        y = y + h * g(k * h, y)
        if k + 1 in marks:
            row[marks[k + 1]] = abs(y - math.cos((k + 1) * h))
    row["x1.0_over_h"] = row["x1.0"] / h
    part_d["euler"].append(row)
# Euler bound for this problem: K = 50, max|y''| = max|cos x| = 1, e_0 = 0
part_d["euler_bound_constant"] = (math.exp(lam) - 1) / (2 * lam)

out = {
    "K": K,
    "max_ypp": M2,
    "euler_bound_constant": C_bound,
    "part_a": part_a,
    "part_b": part_b,
    "part_c": part_c,
    "part_c_fail": part_c_fail,
    "part_d": part_d,
}
Path(__file__).with_name("ode_methods.json").write_text(json.dumps(out, indent=2), encoding="utf-8")

print(f"K = {K:.4f}, max|y''| = {M2:.5f}, Euler bound C = {C_bound:.4f}")
print("A:  n      h        Euler    ratio   midpoint  ratio   trapez   ratio   bound    picard")
for r in part_a:
    print(f"  {r['n']:5d} {r['h']:.5f}  {r['euler']:.3e} {r.get('euler_ratio', 0):5.2f}  {r['midpoint']:.3e} {r.get('midpoint_ratio', 0):5.2f}  "
          f"{r['trapezoidal']:.3e} {r.get('trapezoidal_ratio', 0):5.2f}  {r['euler_bound']:.3e} {r['picard_iterations']:.2f}")
print("B:  h      one step   ratio   global    ratio")
for r in part_b:
    print(f"  {r['h']:.3f}  {r['one_step']:.3e} {r.get('one_step_ratio', 0):5.2f}  {r['global']:.3e} {r.get('global_ratio', 0):5.2f}")
for r in part_c:
    print("C:", r["h"], "errors", [f"{e:.2e}" for e in r["errors"]], "ratios", [f"{q:.4f}" for q in r["ratios"]],
          "h|f_y|/2", f"{r['h_fy_over_2']:.4f}", "iters", r["iterations"])
print("C fail:", part_c_fail["h_K_over_2"], part_c_fail["trapezoidal_root"], [f"{v:.3e}" for v in part_c_fail["iterates"]])
for r in part_d["perturbation"]:
    print(f"D eps={r['eps']:.0e} x={r['x']}: gap {r['gap']:.2e}")
for r in part_d["euler"]:
    print("D Euler", r)
print("D Euler bound constant", f"{part_d['euler_bound_constant']:.3e}")
