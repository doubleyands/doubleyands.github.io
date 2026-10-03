"""Experiment for the post "Numerical analysis, week 1".

Solves x = cos(x) with four methods and records the error after every step.
Writes the numbers to rootfinding.json next to this file.

Run: python scripts/posts/numerical-analysis/rootfinding.py
"""

import json
import math
from pathlib import Path

ALPHA = 0.7390851332151607  # the root of x = cos(x), to double precision
STEPS = 40
TINY = 1e-15  # below this the error is round-off, so the series stops


def f(x):
    return x - math.cos(x)


def df(x):
    return 1 + math.sin(x)


# --- the four methods (this part is shown in the post) ---

def bisection(a, b, steps):
    xs = []
    for _ in range(steps + 1):
        c = (a + b) / 2
        xs.append(c)
        if f(a) * f(c) <= 0:
            b = c
        else:
            a = c
    return xs


def fixed_point(x, steps):  # g(x) = cos(x)
    xs = [x]
    for _ in range(steps):
        x = math.cos(x)
        xs.append(x)
    return xs


def newton(x, steps):
    xs = [x]
    for _ in range(steps):
        x = x - f(x) / df(x)
        xs.append(x)
    return xs


def secant(x0, x1, steps):
    xs = [x0, x1]
    for _ in range(steps - 1):
        if f(x1) == f(x0):  # already converged; the next step would divide by 0
            break
        x0, x1 = x1, x1 - f(x1) * (x1 - x0) / (f(x1) - f(x0))
        xs.append(x1)
    return xs


# --- measurements ---

def first_below(errs, level):
    return next((n for n, e in enumerate(errs) if e < level), None)


def observed_order(errs):
    """p from three consecutive errors, using the last triple still above round-off."""
    n = max(n for n in range(2, len(errs)) if errs[n] > 1e-14)
    e0, e1, e2 = errs[n - 2], errs[n - 1], errs[n]
    return math.log(e2 / e1) / math.log(e1 / e0)


LONG = 200  # long enough for every method to reach round-off
runs = {
    "Bisection": bisection(0.0, 1.0, LONG),
    "Fixed point": fixed_point(1.0, LONG),
    "Secant": secant(0.0, 1.0, LONG),
    "Newton": newton(1.0, LONG),
}
raw = {name: [abs(x - ALPHA) for x in xs] for name, xs in runs.items()}

series = {}
summary = {}
for name, e in raw.items():
    floor_at = first_below(e, TINY)
    shown = e[: min(STEPS + 1, floor_at if floor_at is not None else STEPS + 1)]
    # floor_at is drawn as a drop to the bottom of the chart: the error is now round-off
    series[name] = {"errors": shown, "round_off_at": floor_at if floor_at is not None and floor_at <= STEPS else None}
    summary[name] = {"steps_to_1e-10": first_below(e, 1e-10)}

# Linear methods: the error shrinks by a fixed factor per step
summary["Bisection"]["average_factor"] = (raw["Bisection"][40] / raw["Bisection"][0]) ** (1 / 40)
summary["Fixed point"]["average_factor"] = raw["Fixed point"][31] / raw["Fixed point"][30]
# Faster methods: order p from the last clean triple
for name in ("Secant", "Newton"):
    summary[name]["observed_order"] = observed_order(raw[name])

# Fixed point: e(n+1)/e(n) should approach |g'(alpha)| = sin(alpha)
fp_ratio = raw["Fixed point"][31] / raw["Fixed point"][30]
# Newton: e(n+1)/e(n)^2 should approach |f''(alpha) / (2 f'(alpha))| = cos(alpha) / (2 (1 + sin(alpha)))
nw_ratio = raw["Newton"][3] / raw["Newton"][2] ** 2

result = {
    "alpha": ALPHA,
    "steps": STEPS,
    "series": series,
    "summary": summary,
    "fixed_point_ratio": {"measured": fp_ratio, "theory": math.sin(ALPHA)},
    "newton_ratio": {"measured": nw_ratio, "theory": math.cos(ALPHA) / (2 * (1 + math.sin(ALPHA)))},
}
Path(__file__).with_name("rootfinding.json").write_text(json.dumps(result, indent=1))

for name, s in summary.items():
    extra = " ".join(f"{k}={v:.4f}" for k, v in s.items() if k != "steps_to_1e-10")
    print(f"{name:12s} steps to 1e-10: {s['steps_to_1e-10']!s:>4}   {extra}   round-off at: {series[name]['round_off_at']}")
print(f"fixed point e(31)/e(30): {fp_ratio:.4f}   theory sin(alpha) = {math.sin(ALPHA):.4f}")
print(f"newton e(3)/e(2)^2:      {nw_ratio:.4f}   theory = {result['newton_ratio']['theory']:.4f}")
