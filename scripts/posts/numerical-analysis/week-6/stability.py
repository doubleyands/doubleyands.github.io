"""Experiments for "Numerical Analysis, Week 6" (stability of ODE methods).

Python standard library only (complex numbers from the built-in type).

Part A. Stiff linear problem y' = LAM (y - sin x) + cos x, y(0) = 0, LAM = -1000,
        exact y = sin x on [0, 1]. Max grid error of forward Euler, RK4,
        backward Euler and the trapezoidal method for many step sizes h.
Part B. RK2 and RK4 on y' = y cos x, y(0) = 1 (exact e^{sin x}) on [0, 2]: order check.
Part C. An unstable but consistent (order 3) two-step method,
        y_{n+1} = -4 y_n + 5 y_{n-1} + h (4 f_n + 2 f_{n-1}),
        rho(r) = r^2 + 4r - 5 = (r - 1)(r + 5). Run on y' = -y with exact start values.
Part D. Milne's method on y' = -y: characteristic roots for small negative z,
        boundary locus z(theta), and the error growth over a long interval.
Part E. Amplification factors on y' = lam y, and the real/imaginary stability limits
        of RK2 and RK4 (found by bisection).

Run: python scripts/posts/numerical-analysis/week-6/stability.py
Writes stability.json next to this file.
"""

import cmath
import json
import math
from pathlib import Path

OUT = {}

# ---------------------------------------------------------------- Part A
LAM = -1000.0


def f_stiff(x, y):
    return LAM * (y - math.sin(x)) + math.cos(x)


def g(x):
    """Forcing term: y' = LAM y + g(x)."""
    return -LAM * math.sin(x) + math.cos(x)


def run_stiff(method, n):
    h = 1.0 / n
    y, worst = 0.0, 0.0
    for k in range(n):
        x = k * h
        if method == "forward_euler":
            y = y + h * f_stiff(x, y)
        elif method == "rk4":
            k1 = f_stiff(x, y)
            k2 = f_stiff(x + h / 2, y + h / 2 * k1)
            k3 = f_stiff(x + h / 2, y + h / 2 * k2)
            k4 = f_stiff(x + h, y + h * k3)
            y = y + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        elif method == "backward_euler":  # linear, so solve for y_{n+1} exactly
            y = (y + h * g(x + h)) / (1 - h * LAM)
        elif method == "trapezoidal":
            y = ((1 + h * LAM / 2) * y + h / 2 * (g(x) + g(x + h))) / (1 - h * LAM / 2)
        if not math.isfinite(y):
            return float("inf")
        worst = max(worst, abs(y - math.sin(x + h)))
    return worst


STEPS = [50, 100, 200, 300, 340, 380, 450, 490, 510, 700, 1000, 2000, 4000, 8000]
part_a = []
for n in STEPS:
    row = {"n": n, "h": 1.0 / n, "z": LAM / n}
    for m in ("forward_euler", "rk4", "backward_euler", "trapezoidal"):
        row[m] = run_stiff(m, n)
    part_a.append(row)
OUT["part_a"] = part_a

# ---------------------------------------------------------------- Part B


def f_b(x, y):
    return y * math.cos(x)


def rk(n, stages):
    h = 2.0 / n
    y, worst = 1.0, 0.0
    for k in range(n):
        x = k * h
        if stages == 2:
            y = y + h * f_b(x + h / 2, y + h / 2 * f_b(x, y))
        else:
            k1 = f_b(x, y)
            k2 = f_b(x + h / 2, y + h / 2 * k1)
            k3 = f_b(x + h / 2, y + h / 2 * k2)
            k4 = f_b(x + h, y + h * k3)
            y = y + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        worst = max(worst, abs(y - math.exp(math.sin(x + h))))
    return worst


OUT["part_b"] = [{"n": n, "h": 2.0 / n, "rk2": rk(n, 2), "rk4": rk(n, 4)} for n in (10, 20, 40, 80, 160)]

# ---------------------------------------------------------------- Part C


def unstable(n, xend=1.0):
    h = xend / n
    ys = [1.0, math.exp(-h)]  # exact start values
    for k in range(1, n):
        fn, fm = -ys[k], -ys[k - 1]
        ys.append(-4 * ys[k] + 5 * ys[k - 1] + h * (4 * fn + 2 * fm))
    return abs(ys[n] - math.exp(-xend))


def one_step_error(h):
    """Local error: start from exact y(0), y(h); compare y_2 with y(2h)."""
    y2 = -4 * math.exp(-h) + 5 * 1.0 + h * (4 * -math.exp(-h) + 2 * -1.0)
    return y2 - math.exp(-2 * h)


OUT["part_c"] = {
    "global": [{"n": n, "h": 1.0 / n, "err_at_1": unstable(n)} for n in (5, 10, 20, 40, 80)],
    "local": [{"h": h, "err": one_step_error(h)} for h in (0.1, 0.05, 0.025, 0.0125)],
    "rho_roots": [1, -5],
}

# ---------------------------------------------------------------- Part D


def quad_roots(a, b, c):
    d = cmath.sqrt(b * b - 4 * a * c)
    return (-b + d) / (2 * a), (-b - d) / (2 * a)


def milne_roots(z):
    # (1 - z/3) r^2 - (4z/3) r - (1 + z/3) = 0
    return quad_roots(1 - z / 3, -4 * z / 3, -(1 + z / 3))


def milne_locus(theta):
    r = cmath.exp(1j * theta)
    return (r * r - 1) / ((r * r + 4 * r + 1) / 3)


def milne_run(h, xend):
    n = round(xend / h)
    ys = [1.0, math.exp(-h)]
    for k in range(1, n):
        # implicit, but linear for y' = -y: solve for y_{n+1}
        rhs = ys[k - 1] + h / 3 * (-ys[k - 1] - 4 * ys[k])
        ys.append(rhs / (1 + h / 3))
    return [(k * h, abs(ys[k] - math.exp(-k * h)), ys[k], math.exp(-k * h)) for k in range(0, n + 1)]


roots = []
for z in (-0.5, -0.1, -0.01, -0.001):
    r1, r2 = milne_roots(z)
    roots.append({"z": z, "abs_roots": sorted([abs(r1), abs(r2)], reverse=True), "roots": [[r1.real, r1.imag], [r2.real, r2.imag]]})
locus = [{"theta_deg": d, "z": [milne_locus(math.radians(d)).real, milne_locus(math.radians(d)).imag]} for d in range(0, 360, 30)]
run = milne_run(0.1, 30.0)
OUT["part_d"] = {
    "roots": roots,
    "locus": locus,
    "locus_max_imag": max(abs(milne_locus(2 * math.pi * k / 3600).imag) for k in range(3600)),
    "locus_max_real": max(abs(milne_locus(2 * math.pi * k / 3600).real) for k in range(3600)),
    "run_h0.1": [{"x": x, "err": e, "y": y, "exact": ex} for (x, e, y, ex) in run if abs(x - round(x)) < 1e-9 and round(x) % 5 == 0],
}

# ---------------------------------------------------------------- Part E


def P_rk2(z):
    return 1 + z + z * z / 2


def P_rk4(z):
    return 1 + z + z**2 / 2 + z**3 / 6 + z**4 / 24


def bisect(fun, a, b, it=200):
    fa = fun(a)
    for _ in range(it):
        m = (a + b) / 2
        fm = fun(m)
        if (fm > 0) == (fa > 0):
            a, fa = m, fm
        else:
            b = m
    return (a + b) / 2


# Real interval: |P(-s)| = 1 for s > 0
rk4_real = bisect(lambda s: abs(P_rk4(-s)) - 1, 2.0, 3.0)
rk2_real = bisect(lambda s: abs(P_rk2(-s)) - 1, 1.0, 3.0)
rk4_imag = bisect(lambda t: abs(P_rk4(1j * t)) - 1, 1.0, 3.5)

factors = []
lam = -30.0
for h in (0.01, 0.05, 0.25):
    z = h * lam
    factors.append({
        "h": h, "z": z, "exact": math.exp(z), "forward_euler": 1 + z,
        "backward_euler": 1 / (1 - z), "trapezoidal": (1 + z / 2) / (1 - z / 2),
        "rk2": P_rk2(z), "rk4": P_rk4(z),
    })

# Trapezoidal ratio on the imaginary axis and in the right half plane
trap = lambda z: abs((1 + z / 2) / (1 - z / 2))
OUT["part_e"] = {
    "lam": lam,
    "factors": factors,
    "rk4_real_limit": -rk4_real,
    "rk2_real_limit": -rk2_real,
    "rk4_imag_limit": rk4_imag,
    "two_sqrt2": 2 * math.sqrt(2),
    "trap_on_imag_axis": [trap(1j * t) for t in (0.5, 5.0, 50.0)],
    "trap_right": trap(0.1 + 0j),
    "trap_left_far": trap(-1e6 + 0j),
    "fe_limit_stiff_h": 2 / abs(LAM),
    "rk4_limit_stiff_h": rk4_real / abs(LAM),
}

Path(__file__).with_name("stability.json").write_text(json.dumps(OUT, indent=1), encoding="utf-8")

if __name__ == "__main__":
    print(json.dumps(OUT, indent=1))
