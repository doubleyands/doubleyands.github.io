"""Experiment for the post "Numerical Analysis, Week 2".

Interpolates two functions on [-1, 1] with polynomials of degree n = 2, 4, ..., 32,
using Newton's divided difference form, and records the largest error on a
fine grid of 4001 points.

  bump(x) = 1 / (1 + 16 x^2)   smooth, but has poles at x = +-i/4
  exp(x)  = e^x                 smooth everywhere

Two kinds of nodes:
  equal:     x_i = -1 + 2i/n
  clustered: x_i = -cos(i pi / n)   (more nodes near the two ends)

For e^x with equal nodes it also computes the two sides of the error formula
  |E(x)| = |W(x)| e^xi / (n+1)!,   e^-1 <= e^xi <= e,
so the measured error must lie between max|W| / (e (n+1)!) and e max|W| / (n+1)!.

It also records, for both kinds of nodes, the amplification factor
  max over x of sum_i |l_i(x)|   (l_i = Lagrange basis functions).
If every data value has an error of size d (for example round-off), the
interpolant can be off by up to this factor times d.

As a check, every interpolant is also evaluated with the Lagrange formula, and
the largest difference between the two forms is recorded.

Python standard library only. Writes interpolation.json next to this file.

Run: python scripts/posts/numerical-analysis/week-2/interpolation.py
"""

import json
import math
from pathlib import Path

DEGREES = list(range(2, 33, 2))
GRID = [-1 + 2 * k / 4000 for k in range(4001)]


def bump(x):
    return 1 / (1 + 16 * x * x)


# --- divided differences and Newton form (this part is shown in the post) ---

def divided_differences(xs, ys):
    """Returns f[x0], f[x0,x1], ..., f[x0,...,xn] (the top edge of the table)."""
    col = list(ys)
    top = [col[0]]
    for k in range(1, len(xs)):
        # f[x_i..x_{i+k}] = (f[x_{i+1}..x_{i+k}] - f[x_i..x_{i+k-1}]) / (x_{i+k} - x_i)
        col = [(col[i + 1] - col[i]) / (xs[i + k] - xs[i]) for i in range(len(col) - 1)]
        top.append(col[0])
    return top


def newton_eval(xs, coef, x):
    """p(x) = c0 + (x-x0)(c1 + (x-x1)(c2 + ...)), evaluated from the inside out."""
    p = coef[-1]
    for k in range(len(coef) - 2, -1, -1):
        p = coef[k] + (x - xs[k]) * p
    return p


# --- the rest of the experiment ---

def lagrange_eval(xs, ys, x):
    total = 0.0
    for i, xi in enumerate(xs):
        li = 1.0
        for j, xj in enumerate(xs):
            if j != i:
                li *= (x - xj) / (xi - xj)
        total += ys[i] * li
    return total


def nodes(kind, n):
    if kind == "equal":
        return [-1 + 2 * i / n for i in range(n + 1)]
    return [-math.cos(i * math.pi / n) for i in range(n + 1)]


def run(f, kind):
    errors, form_gaps, where = [], [], []
    for n in DEGREES:
        xs = nodes(kind, n)
        ys = [f(x) for x in xs]
        coef = divided_differences(xs, ys)
        worst, at, gap = 0.0, 0.0, 0.0
        for x in GRID:
            p = newton_eval(xs, coef, x)
            e = abs(f(x) - p)
            if e > worst:
                worst, at = e, x
            gap = max(gap, abs(p - lagrange_eval(xs, ys, x)))
        errors.append(worst)
        where.append(abs(at))
        form_gaps.append(gap)
    return {"errors": errors, "worst_at_abs_x": where, "newton_vs_lagrange": form_gaps}


def amplification(kind):
    out = []
    for n in DEGREES:
        xs = nodes(kind, n)
        worst = 0.0
        for x in GRID:
            total = 0.0
            for i, xi in enumerate(xs):
                li = 1.0
                for j, xj in enumerate(xs):
                    if j != i:
                        li *= (x - xj) / (xi - xj)
                total += abs(li)
            worst = max(worst, total)
        out.append(worst)
    return out


def exp_bounds():
    low, high = [], []
    for n in DEGREES:
        xs = nodes("equal", n)
        w = max(abs(math.prod(x - xi for xi in xs)) for x in GRID)
        low.append(w / (math.e * math.factorial(n + 1)))
        high.append(w * math.e / math.factorial(n + 1))
    return low, high


def main():
    series = {
        "bump, equal nodes": run(bump, "equal"),
        "bump, clustered nodes": run(bump, "clustered"),
        "exp, equal nodes": run(math.exp, "equal"),
    }
    low, high = exp_bounds()
    amp = {"equal": amplification("equal"), "clustered": amplification("clustered")}
    out = {"degrees": DEGREES, "series": series, "exp_bound_low": low, "exp_bound_high": high,
           "amplification": amp}
    Path(__file__).with_name("interpolation.json").write_text(json.dumps(out, indent=1), encoding="utf-8")

    print("n   bump/equal  |x|worst  bump/clustered  exp/equal   exp low     exp high")
    for i, n in enumerate(DEGREES):
        print(f"{n:<3} {series['bump, equal nodes']['errors'][i]:10.3e}  "
              f"{series['bump, equal nodes']['worst_at_abs_x'][i]:.4f}  "
              f"{series['bump, clustered nodes']['errors'][i]:12.3e}    "
              f"{series['exp, equal nodes']['errors'][i]:10.3e}  {low[i]:10.3e}  {high[i]:10.3e}")
    print("n   amplification equal / clustered")
    for i, n in enumerate(DEGREES):
        print(f"{n:<3} {amp['equal'][i]:10.3e}  {amp['clustered'][i]:8.3f}")
    for name, s in series.items():
        print(name, "max Newton-Lagrange gap:", f"{max(s['newton_vs_lagrange']):.1e}")


if __name__ == "__main__":
    main()
