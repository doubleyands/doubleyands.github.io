"""Experiment for the post "Numerical PDEs, Lecture 4".

Part 1 (2D Laplace, Dirichlet): 5-point and 9-point stencils on the unit
square. Exact solution u = exp(x) sin(y), which is harmonic (f = 0).
Then a Poisson problem with u = exp(x + y), f = 2 exp(x + y), where the
9-point stencil is used with the plain right-hand side f.

Also applies the 9-point weights to the harmonic polynomials Re (x + iy)^k
and Im (x + iy)^k at the origin with h = 1, in exact integer arithmetic.
A zero result for every k <= 7 means the 9-point error is O(h^6) for a
harmonic u.

Part 2 (orderings): builds the 5-point matrix for a 4 x 4 interior grid in
natural row order and in red-black order, checks that the red-red and
black-black blocks are diagonal, and checks that both orders give the same
grid values.

Part 3 (2D Neumann): u = exp(x + y) again, Neumann conditions on x = 0 and
x = 1 (du/dn with the outward normal), Dirichlet on y = 0 and y = 1. The
boundary row at each Neumann node is either a ghost node row or a one-sided
first-order row.

Part 4 (advection-diffusion): the exact solution of
  -eps u'' + u' = 1 on (0, 1), u(0) = u(1) = 0,
  u(x) = x - (exp((x - 1)/eps) - exp(-1/eps)) / (1 - exp(-1/eps)),
for several eps. Measures u(1/2), the peak, and the width of the layer.

Uses numpy for the 2D linear systems (dense solves, n up to 64).

Writes the numbers to experiment.json next to this file.

Run: python scripts/posts/numerical-pde/lecture-4/experiment.py
"""

import json
import math
from pathlib import Path

import numpy as np

NS = [4, 8, 16, 32, 64]


def solve_dirichlet(n, u_exact, f, stencil):
    """Dirichlet problem on (0,1)^2, interior points in natural row order."""
    h = 1.0 / n
    m = n - 1
    idx = lambda i, j: (i - 1) + m * (j - 1)
    A = np.zeros((m * m, m * m))
    b = np.zeros(m * m)
    if stencil == 5:
        weights = {(0, 0): -4.0, (1, 0): 1.0, (-1, 0): 1.0, (0, 1): 1.0, (0, -1): 1.0}
        scale = 1.0 / h**2
    else:
        weights = {(0, 0): -20.0}
        for d in [(1, 0), (-1, 0), (0, 1), (0, -1)]:
            weights[d] = 4.0
        for d in [(1, 1), (1, -1), (-1, 1), (-1, -1)]:
            weights[d] = 1.0
        scale = 1.0 / (6 * h**2)
    for j in range(1, n):
        for i in range(1, n):
            row = idx(i, j)
            b[row] = f(i * h, j * h)
            for (di, dj), w in weights.items():
                ii, jj = i + di, j + dj
                if 1 <= ii <= m and 1 <= jj <= m:
                    A[row, idx(ii, jj)] += w * scale
                else:
                    b[row] -= w * scale * u_exact(ii * h, jj * h)
    U = np.linalg.solve(A, b)
    err = max(abs(U[idx(i, j)] - u_exact(i * h, j * h)) for j in range(1, n) for i in range(1, n))
    return err


def order(errs):
    return [None] + [math.log2(errs[k - 1] / errs[k]) for k in range(1, len(errs))]


def part1():
    lap_u = lambda x, y: math.exp(x) * math.sin(y)
    lap_f = lambda x, y: 0.0
    poi_u = lambda x, y: math.exp(x + y)
    poi_f = lambda x, y: 2 * math.exp(x + y)
    out = {}
    for name, u, f in [("laplace", lap_u, lap_f), ("poisson", poi_u, poi_f)]:
        for s in (5, 9):
            errs = [solve_dirichlet(n, u, f, s) for n in NS]
            out[f"{name}_{s}"] = {"err": errs, "order": order(errs)}
    out["harmonic_check"] = harmonic_check()
    return out


def harmonic_check():
    weights = {(0, 0): -20}
    for d in [(1, 0), (-1, 0), (0, 1), (0, -1)]:
        weights[d] = 4
    for d in [(1, 1), (1, -1), (-1, 1), (-1, -1)]:
        weights[d] = 1

    def power(z, k):
        r = (1, 0)
        for _ in range(k):
            r = (r[0] * z[0] - r[1] * z[1], r[0] * z[1] + r[1] * z[0])
        return r

    out = []
    for k in range(10):
        re = sum(w * power(d, k)[0] for d, w in weights.items())
        im = sum(w * power(d, k)[1] for d, w in weights.items())
        out.append({"k": k, "re": re, "im": im})
    return out


def part2():
    n, m = 5, 4
    h = 1.0 / n
    pts = [(i, j) for j in range(1, n) for i in range(1, n)]
    natural = {p: k for k, p in enumerate(pts)}
    red = [p for p in pts if (p[0] + p[1]) % 2 == 0]
    black = [p for p in pts if (p[0] + p[1]) % 2 == 1]
    rb = {p: k for k, p in enumerate(red + black)}

    def build(order_map):
        N = len(order_map)
        A = np.zeros((N, N))
        b = np.zeros(N)
        u = lambda x, y: math.exp(x + y)
        for (i, j), row in order_map.items():
            b[row] = 2 * math.exp((i + j) * h)
            A[row, row] = -4 / h**2
            for di, dj in [(1, 0), (-1, 0), (0, 1), (0, -1)]:
                q = (i + di, j + dj)
                if q in order_map:
                    A[row, order_map[q]] = 1 / h**2
                else:
                    b[row] -= u(q[0] * h, q[1] * h) / h**2
        return A, b

    An, bn = build(natural)
    Ar, br = build(rb)
    nr = len(red)
    D_red, D_black = Ar[:nr, :nr], Ar[nr:, nr:]
    off = lambda M: int(np.count_nonzero(M - np.diag(np.diag(M))))
    Un = np.linalg.solve(An, bn)
    Ur = np.linalg.solve(Ar, br)
    diff = max(abs(Un[natural[p]] - Ur[rb[p]]) for p in pts)
    H = Ar[:nr, nr:]
    return {
        "unknowns": len(pts),
        "red": nr,
        "black": len(black),
        "offdiag_nonzeros_red_block": off(D_red),
        "offdiag_nonzeros_black_block": off(D_black),
        "diag_red": sorted(set(np.diag(D_red) * h**2)),
        "H_symmetric_pair": bool(np.array_equal(Ar[nr:, :nr], H.T)),
        "H_nonzeros_per_row": sorted(set(int(c) for c in np.count_nonzero(H, axis=1))),
        "nonzeros_natural": int(np.count_nonzero(An)),
        "nonzeros_redblack": int(np.count_nonzero(Ar)),
        "max_diff_between_orders": float(diff),
        "red_red_neighbors": sum(
            1 for (i, j) in red for di, dj in [(1, 0), (-1, 0), (0, 1), (0, -1)] if (i + di, j + dj) in red
        ),
    }


def solve_neumann(n, kind):
    """u = exp(x+y). Unknowns u_{i,j}, i = 0..n, j = 1..n-1 (natural row order)."""
    h = 1.0 / n
    u = lambda x, y: math.exp(x + y)
    f = lambda x, y: 2 * math.exp(x + y)
    g_left = lambda y: -math.exp(y)          # du/dn at x = 0, n = (-1, 0)
    g_right = lambda y: math.exp(1 + y)      # du/dn at x = 1, n = (1, 0)
    nx, ny = n + 1, n - 1
    idx = lambda i, j: i + nx * (j - 1)
    N = nx * ny
    A = np.zeros((N, N))
    b = np.zeros(N)

    def add(row, i, j, w):
        if 1 <= j <= n - 1:
            A[row, idx(i, j)] += w
        else:
            b[row] -= w * u(i * h, j * h)  # Dirichlet value on y = 0 or y = 1

    for j in range(1, n):
        y = j * h
        for i in range(0, n + 1):
            row = idx(i, j)
            x = i * h
            if 1 <= i <= n - 1:
                b[row] = f(x, y)
                add(row, i, j, -4 / h**2)
                for di, dj in [(1, 0), (-1, 0), (0, 1), (0, -1)]:
                    add(row, i + di, j + dj, 1 / h**2)
            elif kind == "ghost":
                # 5-point row at the boundary node, ghost value removed with
                # the central difference of du/dn
                inner = 1 if i == 0 else n - 1
                g = g_left(y) if i == 0 else g_right(y)
                b[row] = f(x, y) - 2 * g / h
                add(row, i, j, -4 / h**2)
                add(row, inner, j, 2 / h**2)
                add(row, i, j + 1, 1 / h**2)
                add(row, i, j - 1, 1 / h**2)
            else:
                # (u_boundary - u_inner)/h = du/dn
                inner = 1 if i == 0 else n - 1
                g = g_left(y) if i == 0 else g_right(y)
                b[row] = g
                A[row, idx(i, j)] += 1 / h
                A[row, idx(inner, j)] -= 1 / h
    U = np.linalg.solve(A, b)
    err = max(abs(U[idx(i, j)] - u(i * h, j * h)) for j in range(1, n) for i in range(0, n + 1))
    return err, N


def part3():
    out = {}
    for kind in ("ghost", "one-sided"):
        res = [solve_neumann(n, kind) for n in NS]
        errs = [r[0] for r in res]
        out[kind] = {"err": errs, "order": order(errs), "unknowns": [r[1] for r in res]}
    return out


def u_ad(x, eps):
    return x - (math.exp((x - 1) / eps) - math.exp(-1 / eps)) / (1 - math.exp(-1 / eps))


def part4():
    rows = []
    for eps in [10.0, 1.0, 0.1, 0.01, 0.001]:
        # Peak: u'(x) = 0  <=>  exp((x - 1)/eps) = eps (1 - exp(-1/eps))
        x_star = 1 + eps * math.log(eps * (1 - math.exp(-1 / eps)))
        peak = u_ad(x_star, eps)
        row = {
            "eps": eps,
            "u_half": u_ad(0.5, eps),
            "x_star": x_star,
            "peak": peak,
            "diffusion_only_peak": 1 / (8 * eps),
        }
        if eps <= 0.1:
            # Layer width: the point where u - x first reaches -0.01, from the left
            lo, hi = 0.0, 1.0
            for _ in range(200):
                mid = 0.5 * (lo + hi)
                if u_ad(mid, eps) - mid > -0.01:
                    lo = mid
                else:
                    hi = mid
            row["width"] = 1 - lo
            row["width_over_eps"] = (1 - lo) / eps
            row["peak_approx"] = 1 - eps * math.log(1 / eps) - eps
        rows.append(row)
    # u(1/2) as eps goes down: check it never decreases
    sweep = [10 ** (1 - 4 * k / 400) for k in range(401)]
    halves = [u_ad(0.5, e) for e in sweep]
    monotone = all(b >= a for a, b in zip(halves, halves[1:]))
    curves = {}
    xs = [k / 200 for k in range(201)]
    for eps in [1.0, 0.1, 0.03, 0.01]:
        curves[str(eps)] = [u_ad(x, eps) for x in xs]
    return {
        "table": rows,
        "u_half_monotone": monotone,
        "u_half_min": min(halves),
        "u_half_max": max(halves),
        "xs": xs,
        "curves": curves,
    }


if __name__ == "__main__":
    data = {"n": NS, "part1": part1(), "part2": part2(), "part3": part3(), "part4": part4()}
    Path(__file__).with_name("experiment.json").write_text(json.dumps(data, indent=1), encoding="utf-8")
    p = data
    print([(r["k"], r["re"], r["im"]) for r in p["part1"]["harmonic_check"]])
    for k, v in p["part1"].items():
        if k == "harmonic_check":
            continue
        print(k, ["%.3e" % e for e in v["err"]], ["%.2f" % o for o in v["order"][1:]])
    print(json.dumps(p["part2"], indent=1))
    for k, v in p["part3"].items():
        print(k, ["%.3e" % e for e in v["err"]], ["%.3f" % o for o in v["order"][1:]], v["unknowns"])
    for r in p["part4"]["table"]:
        print({k: (round(v, 6) if isinstance(v, float) else v) for k, v in r.items()})
    print(p["part4"]["u_half_monotone"], p["part4"]["u_half_min"], p["part4"]["u_half_max"])
