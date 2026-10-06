"""Experiment for the post "Numerical PDEs, Lecture 1".

Measures the error of four finite difference formulas for u(x) = exp(sin x)
at x = 1, with h = 2^-1, 2^-2, ..., 2^-30. Also evaluates the stability
example u'' = u, u(0) = 2, u'(0) = a on growing intervals.
Writes the numbers to differences.json next to this file.

Run: python scripts/posts/numerical-pde/lecture-1/differences.py
"""

import json
import math
from pathlib import Path

X0 = 1.0
K_MAX = 30


def u(x):
    return math.exp(math.sin(x))


def du(x):
    return math.cos(x) * math.exp(math.sin(x))


def d2u(x):
    return (math.cos(x) ** 2 - math.sin(x)) * math.exp(math.sin(x))


# --- the four formulas (this part is shown in the post) ---

def forward(u, x, h):
    return (u(x + h) - u(x)) / h


def central(u, x, h):
    return (u(x + h) - u(x - h)) / (2 * h)


def fourth_order(u, x, h):
    return (u(x - 2*h) - 8*u(x - h) + 8*u(x + h) - u(x + 2*h)) / (12 * h)


def second_central(u, x, h):
    return (u(x - h) - 2*u(x) + u(x + h)) / h**2

# --- end of the shown part ---


FORMULAS = {
    "Forward": (forward, du),
    "Central": (central, du),
    "Fourth order": (fourth_order, du),
    "Second difference": (second_central, d2u),
}


def observed_order(e_coarse, e_fine):
    """Order p from two errors with h halved: e_coarse / e_fine = 2^p."""
    return math.log2(e_coarse / e_fine)


ks = list(range(1, K_MAX + 1))
series = {}
for name, (formula, exact) in FORMULAS.items():
    errors = [abs(formula(u, X0, 2.0**-k) - exact(X0)) for k in ks]
    best = min(range(len(ks)), key=lambda i: errors[i])
    series[name] = {
        "errors": errors,
        "orders": [None] + [observed_order(errors[i - 1], errors[i]) for i in range(1, len(ks))],
        "best_k": ks[best],
        "best_error": errors[best],
    }

# Stability example: u'' = u, u(0) = 2, u'(0) = -2 + eps. Unperturbed solution: 2 exp(-x).
# The change in the solution is eps * sinh(x), so the constant C on [0, L] is sinh(L).
stability = []
eps = 1e-6
for L in (1, 5, 10, 20, 40):
    grid = [L * i / 1000 for i in range(1001)]
    exact = [2 * math.exp(-x) for x in grid]
    perturbed = [eps / 2 * math.exp(x) + (2 - eps / 2) * math.exp(-x) for x in grid]
    gap = max(abs(p - q) for p, q in zip(perturbed, exact))
    stability.append({"L": L, "max_change": gap, "C": gap / eps, "sinh_L": math.sinh(L)})

out = {"x0": X0, "k": ks, "series": series, "stability": {"eps": eps, "rows": stability}}
Path(__file__).with_name("differences.json").write_text(json.dumps(out, indent=1), encoding="utf-8")

for name, s in series.items():
    print(f"{name:18s} best k={s['best_k']:2d} err={s['best_error']:.2e}")
    print("   orders k=2..12:", " ".join(f"{p:.2f}" for p in s["orders"][1:12]))
    print("   errors:", " ".join(f"{e:.1e}" for e in s["errors"]))
for row in stability:
    print(row)
