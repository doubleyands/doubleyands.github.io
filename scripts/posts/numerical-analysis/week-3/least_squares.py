"""Experiment for the post "Numerical Analysis, Week 3".

Finds the best L2 polynomial approximation r_n of degree n = 0..30 to

  f(x) = 1 / (2 + x)   on [-1, 1],  weight w = 1,

in two ways that give the same polynomial in exact arithmetic:

  monomial:  solve the normal equations  sum_j a_j (x^i, x^j) = (f, x^i)
             with Gaussian elimination (partial pivoting).
             (x^i, x^j) = 2 / (i + j + 1) if i + j is even, else 0 (exact).
  shifted:   the same normal equations, but in the variable t = (x + 1) / 2
             on [0, 1]. There (t^i, t^j) = 1 / (i + j + 1): the Hilbert matrix.
             The best polynomial is the same; only the basis changes.
  Legendre:  a_j = (f, phi_j) with the orthonormal Legendre polynomials
             phi_j = sqrt((2j+1)/2) P_j, built by the three-term recurrence.

For each n it records
  - the L2 error ||f - r_n||_2 of both, measured with a 200-point Gauss rule,
  - the L2 error predicted from the Legendre coefficients alone:
      ||f - r_n||^2 = ||f||^2 - sum_{j<=n} a_j^2,  ||f||^2 = 2/3 (exact),
  - the sum sum_{j<=n} a_j^2 (Bessel: it never passes 2/3),
  - the condition numbers of the two monomial matrices in the infinity norm,
    ||A|| ||A^-1||, with A^-1 from the same elimination (rough when huge).

It also compares the cubic Taylor polynomial at 0 with the cubic least squares
polynomial (largest error on 4001 points).

Inner products (f, g) are computed with Gauss-Legendre quadrature. f is smooth
on [-1, 1] (its pole is at x = -2), so 200 points give full double precision.

Python standard library only. Writes least_squares.json next to this file.

Run: python scripts/posts/numerical-analysis/week-3/least_squares.py
"""

import json
import math
from pathlib import Path

N_MAX = 30
NORM_F_SQ = 2.0 / 3.0  # integral of (2+x)^-2 over [-1, 1] = 1 - 1/3


def f(x):
    return 1.0 / (2.0 + x)


def legendre_all(n, x):
    """P_0(x), ..., P_n(x) by (k+1) P_{k+1} = (2k+1) x P_k - k P_{k-1}."""
    p = [1.0, x]
    for k in range(1, n):
        p.append(((2 * k + 1) * x * p[k] - k * p[k - 1]) / (k + 1))
    return p[: n + 1]


def gauss_legendre(m):
    """Nodes and weights of the m-point Gauss-Legendre rule on [-1, 1]."""
    nodes, weights = [], []
    for i in range(m):
        x = math.cos(math.pi * (i + 0.75) / (m + 0.5))
        for _ in range(100):
            p0, p1 = 1.0, x
            for k in range(1, m):
                p0, p1 = p1, ((2 * k + 1) * x * p1 - k * p0) / (k + 1)
            dp = m * (x * p1 - p0) / (x * x - 1.0)  # P_m'(x)
            step = p1 / dp
            x -= step
            if abs(step) < 1e-16:
                break
        nodes.append(x)
        weights.append(2.0 / ((1.0 - x * x) * dp * dp))
    return nodes, weights


def phi_all(n, x):
    return [math.sqrt((2 * j + 1) / 2.0) * p for j, p in enumerate(legendre_all(max(n, 1), x)[: n + 1])]


def solve(a, b):
    """Gaussian elimination with partial pivoting. Returns x with a x = b."""
    m = len(b)
    a = [row[:] + [b[i]] for i, row in enumerate(a)]
    for k in range(m):
        piv = max(range(k, m), key=lambda r: abs(a[r][k]))
        a[k], a[piv] = a[piv], a[k]
        for r in range(k + 1, m):
            t = a[r][k] / a[k][k]
            for c in range(k, m + 1):
                a[r][c] -= t * a[k][c]
    x = [0.0] * m
    for k in range(m - 1, -1, -1):
        x[k] = (a[k][m] - sum(a[k][c] * x[c] for c in range(k + 1, m))) / a[k][k]
    return x


def cond_inf(a):
    m = len(a)
    cols = [solve(a, [1.0 if i == j else 0.0 for i in range(m)]) for j in range(m)]
    inv_norm = max(sum(abs(cols[j][i]) for j in range(m)) for i in range(m))
    return max(sum(abs(v) for v in row) for row in a) * inv_norm


def poly_eval(coef, x):
    p = 0.0
    for c in reversed(coef):
        p = p * x + c
    return p


def main():
    gx, gw = gauss_legendre(200)
    fx = [f(x) for x in gx]
    phis = [phi_all(N_MAX, x) for x in gx]
    grid = [-1.0 + 2.0 * i / 4000 for i in range(4001)]

    def l2(values):
        return math.sqrt(sum(w * v * v for w, v in zip(gw, values)))

    # Legendre coefficients a_j = (f, phi_j)
    leg = [sum(w * fv * ph[j] for w, fv, ph in zip(gw, fx, phis)) for j in range(N_MAX + 1)]
    # (f, x^i) for the normal equations
    rhs_all = [sum(w * fv * x**i for w, fv, x in zip(gw, fx, gx)) for i in range(N_MAX + 1)]
    # (f, t^i) on [0, 1] with t = (x + 1) / 2, dt = dx / 2
    tx = [(x + 1.0) / 2.0 for x in gx]
    rhs_t = [sum(w * fv * t**i for w, fv, t in zip(gw, fx, tx)) / 2.0 for i in range(N_MAX + 1)]

    rows = []
    for n in range(N_MAX + 1):
        r_leg = [sum(leg[j] * ph[j] for j in range(n + 1)) for ph in phis]
        err_leg = l2([a - b for a, b in zip(fx, r_leg)])
        bessel = sum(a * a for a in leg[: n + 1])
        predicted = math.sqrt(max(NORM_F_SQ - bessel, 0.0))

        mat = [[(2.0 / (i + j + 1) if (i + j) % 2 == 0 else 0.0) for j in range(n + 1)] for i in range(n + 1)]
        mono = solve(mat, rhs_all[: n + 1])
        err_mono = l2([fv - poly_eval(mono, x) for fv, x in zip(fx, gx)])
        hilbert = [[1.0 / (i + j + 1) for j in range(n + 1)] for i in range(n + 1)]
        shifted = solve(hilbert, rhs_t[: n + 1])
        err_shift = l2([fv - poly_eval(shifted, t) for fv, t in zip(fx, tx)])
        rows.append(
            {
                "n": n,
                "legendre_coef": leg[n],
                "bessel_sum": bessel,
                "l2_legendre": err_leg,
                "l2_predicted": predicted,
                "l2_monomial": err_mono,
                "cond_monomial": cond_inf(mat),
                "l2_shifted": err_shift,
                "cond_hilbert": cond_inf(hilbert),
            }
        )
        print(f"n={n:2d} a_n={leg[n]: .3e} sum={bessel:.15f} L2 leg={err_leg:.3e} pred={predicted:.3e} "
              f"mono={err_mono:.3e} cond={rows[-1]['cond_monomial']:.2e} "
              f"shift={err_shift:.3e} condH={rows[-1]['cond_hilbert']:.2e}")

    # Monomial coefficients of r_n for small n: how they move when n grows
    mono_coefs = {}
    for n in range(4):
        mat = [[(2.0 / (i + j + 1) if (i + j) % 2 == 0 else 0.0) for j in range(n + 1)] for i in range(n + 1)]
        mono_coefs[n] = solve(mat, rhs_all[: n + 1])
        print(f"r_{n} monomial coefficients:", [round(c, 6) for c in mono_coefs[n]])

    # Cubic Taylor at 0 vs cubic least squares: 1/(2+x) = sum (-1)^k x^k / 2^(k+1)
    taylor = [(-1) ** k / 2 ** (k + 1) for k in range(4)]
    taylor_err = [abs(f(x) - poly_eval(taylor, x)) for x in grid]
    ls_err = [abs(f(x) - poly_eval(mono_coefs[3], x)) for x in grid]
    cubic = {
        "taylor_max": max(taylor_err),
        "taylor_argmax": grid[taylor_err.index(max(taylor_err))],
        "taylor_at_0.5": abs(f(0.5) - poly_eval(taylor, 0.5)),
        "ls_max": max(ls_err),
        "ls_argmax": grid[ls_err.index(max(ls_err))],
        "ls_sign_changes": sum(
            1 for a, b in zip(grid[:-1], grid[1:])
            if (f(a) - poly_eval(mono_coefs[3], a)) * (f(b) - poly_eval(mono_coefs[3], b)) < 0
        ),
    }
    print("cubic:", cubic)

    out = {
        "function": "1/(2+x) on [-1, 1], w = 1",
        "norm_f_sq": NORM_F_SQ,
        "rows": rows,
        "monomial_coefficients": {str(k): v for k, v in mono_coefs.items()},
        "cubic": cubic,
    }
    Path(__file__).with_name("least_squares.json").write_text(json.dumps(out, indent=1), encoding="utf-8")


if __name__ == "__main__":
    main()
