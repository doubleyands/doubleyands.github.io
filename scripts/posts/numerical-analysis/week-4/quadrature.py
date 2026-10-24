"""Experiment for the post "Numerical Analysis, Week 4".

Part A. Convergence of four rules on a smooth integrand

  I = integral of 1 / (1 + x) over [0, 1] = ln 2,

with N = 2^k + 1 function values, k = 1..10:
  - composite trapezoidal rule, n = N - 1 subintervals,
  - corrected trapezoidal rule (adds -h^2/12 (f'(1) - f'(0)); f' is exact),
  - composite Simpson rule on the same N points,
  - Gauss-Legendre rule with N points (nodes by Newton's method on the
    Legendre recurrence, weights 2 / ((1 - x^2) P_N'(x)^2)).
  It also records the ratio of the asymptotic error estimate
  -h^2/12 (f'(1) - f'(0)) to the true trapezoidal error.

Part B. Equally spaced vs Gauss nodes on 1 / (1 + 16 x^2) over [-1, 1]
  (exact value atan(4) / 2). Closed Newton-Cotes weights w_j = integral of
  the Lagrange basis l_j are computed exactly with fractions. Compared with
  Gauss-Legendre and composite Simpson on the same number of points.

Part C. A singular integrand: integral of e^x / sqrt(x) over [0, 1].
  Gauss-Legendre directly on it, and after x = u^2 (integrand 2 e^(u^2)).
  Exact value: 2 * sum_k 1 / (k! (2k + 1)).

Part D. Small checks: Simpson and 2-point Gauss on x^k, k = 0..5.

Part E. The piece near the singularity, done by hand: the integral of
  e^x log(x) over [0, eps], eps = 0.1, from the Taylor series of e^x and
  the exact integral of x^n log(x). Partial sums with 1, 2, 3 terms are
  compared with the sum of 30 terms.

Python standard library only. Writes quadrature.json next to this file.
Run: python scripts/posts/numerical-analysis/week-4/quadrature.py
"""

import json
import math
from fractions import Fraction
from pathlib import Path


# ---------- rules ----------

def trapezoid(f, a, b, n):
    h = (b - a) / n
    return h * (0.5 * f(a) + sum(f(a + j * h) for j in range(1, n)) + 0.5 * f(b))


def simpson(f, a, b, n):
    """Composite Simpson with n (even) subintervals, n + 1 points."""
    h = (b - a) / n
    s = f(a) + f(b)
    s += 4 * sum(f(a + j * h) for j in range(1, n, 2))
    s += 2 * sum(f(a + j * h) for j in range(2, n, 2))
    return h / 3 * s


def legendre(n, x):
    """P_n(x) and P_n'(x) by the three-term recurrence."""
    p0, p1 = 1.0, x
    if n == 0:
        return 1.0, 0.0
    for k in range(1, n):
        p0, p1 = p1, ((2 * k + 1) * x * p1 - k * p0) / (k + 1)
    dp = n * (x * p1 - p0) / (x * x - 1)
    return p1, dp


_gauss_cache = {}


def gauss_nodes(n):
    """Nodes and weights of the n-point Gauss-Legendre rule on [-1, 1]."""
    if n in _gauss_cache:
        return _gauss_cache[n]
    nodes, weights = [], []
    for i in range(1, n + 1):
        x = math.cos(math.pi * (i - 0.25) / (n + 0.5))  # standard first guess
        for _ in range(100):
            p, dp = legendre(n, x)
            dx = p / dp
            x -= dx
            if abs(dx) < 1e-16:
                break
        p, dp = legendre(n, x)
        nodes.append(x)
        weights.append(2 / ((1 - x * x) * dp * dp))
    _gauss_cache[n] = (nodes, weights)
    return nodes, weights


def gauss(f, a, b, n):
    nodes, weights = gauss_nodes(n)
    mid, half = (a + b) / 2, (b - a) / 2
    return half * sum(w * f(mid + half * x) for x, w in zip(nodes, weights))


def newton_cotes_weights(n):
    """Exact closed Newton-Cotes weights for n + 1 equal points on [-1, 1]."""
    xs = [Fraction(-1) + Fraction(2 * j, n) for j in range(n + 1)]
    weights = []
    for j in range(n + 1):
        # l_j as a coefficient list (lowest power first)
        poly = [Fraction(1)]
        denom = Fraction(1)
        for m in range(n + 1):
            if m == j:
                continue
            poly = [Fraction(0)] + poly  # times x
            for k in range(len(poly) - 1):
                poly[k] -= xs[m] * poly[k + 1]
            denom *= xs[j] - xs[m]
        integral = sum(c * Fraction(1 - (-1) ** (k + 1), k + 1) for k, c in enumerate(poly))
        weights.append(integral / denom)
    return [float(x) for x in xs], weights


# ---------- Part A ----------

f_a = lambda x: 1 / (1 + x)
df_a = lambda x: -1 / (1 + x) ** 2
exact_a = math.log(2)

part_a = []
for k in range(1, 11):
    n = 2 ** k
    h = 1 / n
    t = trapezoid(f_a, 0.0, 1.0, n)
    est = -h * h / 12 * (df_a(1.0) - df_a(0.0))
    ct = t + est
    s = simpson(f_a, 0.0, 1.0, n)
    g = gauss(f_a, 0.0, 1.0, n + 1)
    part_a.append({
        "points": n + 1,
        "trapezoid": exact_a - t,
        "corrected": exact_a - ct,
        "simpson": exact_a - s,
        "gauss": exact_a - g,
        "estimate_over_error": est / (exact_a - t),
    })
for i in range(1, len(part_a)):
    for key in ("trapezoid", "corrected", "simpson"):
        now = part_a[i][key]
        part_a[i][key + "_ratio"] = part_a[i - 1][key] / now if now != 0 else None

# ---------- Part B ----------

f_b = lambda x: 1 / (1 + 16 * x * x)
exact_b = math.atan(4) / 2
part_b = []
for points in (3, 5, 7, 9, 11, 13, 15, 17, 21):
    n = points - 1
    xs, ws = newton_cotes_weights(n)
    nc = sum(float(w) * f_b(x) for x, w in zip(xs, ws))
    row = {
        "points": points,
        "newton_cotes": exact_b - nc,
        "min_weight": float(min(ws)),
        "gauss": exact_b - gauss(f_b, -1.0, 1.0, points),
        "simpson": exact_b - simpson(f_b, -1.0, 1.0, n) if n % 2 == 0 else None,
    }
    part_b.append(row)

# ---------- Part C ----------

exact_c = 2 * sum(1 / (math.factorial(k) * (2 * k + 1)) for k in range(30))
part_c = []
for points in (2, 4, 8, 16, 32, 64):
    direct = gauss(lambda x: math.exp(x) / math.sqrt(x), 0.0, 1.0, points)
    changed = gauss(lambda u: 2 * math.exp(u * u), 0.0, 1.0, points)
    part_c.append({"points": points, "direct": exact_c - direct, "substituted": exact_c - changed})

# ---------- Part D ----------

part_d = []
g2 = 1 / math.sqrt(3)
for k in range(6):
    exact = Fraction(1, k + 1) * (1 - (-1) ** (k + 1))  # integral over [-1, 1]
    simp = Fraction(1, 3) * ((-1) ** k + 4 * (0 ** k) + 1)
    gs = g2 ** k + (-g2) ** k
    part_d.append({"k": k, "exact": str(exact), "simpson": str(simp), "gauss2": gs})

# ---------- Part E ----------

eps = 0.1
terms = [eps ** (n + 1) / (math.factorial(n) * (n + 1)) * (math.log(eps) - 1 / (n + 1)) for n in range(30)]
i1 = sum(terms)
part_e = {"eps": eps, "value": i1, "partial_errors": [sum(terms[:m]) - i1 for m in (1, 2, 3)]}

out = {
    "exact_a": exact_a,
    "part_a": part_a,
    "exact_b": exact_b,
    "part_b": part_b,
    "exact_c": exact_c,
    "part_c": part_c,
    "part_d": part_d,
    "part_e": part_e,
}
Path(__file__).with_name("quadrature.json").write_text(json.dumps(out, indent=1), encoding="utf-8")

print("Part A  (ln 2 =", exact_a, ")")
print(f"{'N':>5} {'trap':>10} {'ratio':>6} {'CT':>10} {'ratio':>6} {'Simpson':>10} {'ratio':>6} {'Gauss':>10} {'est/err':>8}")
for r in part_a:
    print(f"{r['points']:>5} {r['trapezoid']:>10.2e} {r.get('trapezoid_ratio') or 0:>6.2f} {r['corrected']:>10.2e} "
          f"{r.get('corrected_ratio') or 0:>6.2f} {r['simpson']:>10.2e} {r.get('simpson_ratio') or 0:>6.2f} {r['gauss']:>10.2e} {r['estimate_over_error']:>8.5f}")
print("\nPart B  (atan(4)/2 =", exact_b, ")")
for r in part_b:
    s = f"{r['simpson']:.2e}" if r["simpson"] is not None else "-"
    print(f"{r['points']:>3}  NC {r['newton_cotes']:>10.3e}  min w {r['min_weight']:>10.3e}  Gauss {r['gauss']:>10.3e}  Simpson {s}")
print("\nPart C  (exact =", exact_c, ")")
for r in part_c:
    print(f"{r['points']:>3}  direct {r['direct']:>10.3e}  after x=u^2 {r['substituted']:>10.3e}")
print("\nPart D")
for r in part_d:
    print(r)
print("\nPart E", part_e)
