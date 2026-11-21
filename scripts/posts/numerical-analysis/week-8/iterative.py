"""Week 8 experiments: condition numbers and iterative methods for Ax = b.

Python standard library only (fractions for the exact Hilbert inverse). Matrices are lists of rows.

Part A  Condition numbers: cond_2 of the Hilbert matrix H_n, and a vector x_hat with a tiny
        relative residual but a large relative error.
Part B  The bound 1/cond <= (rel. error)/(rel. residual) <= cond, tested with many random
        perturbations of b for a nearly singular 2x2 matrix.
Part C  Residual correction with C = D^{-1} (= Jacobi): the error ratio, the estimate of c that
        does not need x, and rho(I - CA).
Part D  Jacobi, Gauss-Seidel, SOR (optimal w) and CG on the 1D Poisson matrix tridiag(-1, 2, -1):
        error history for n = 31, iteration counts for several n, and cond growing like n^2.
Part E  A nonsymmetric matrix: CG on A fails, CG on the normal equations A^t A x = A^t b works,
        and cond_2(A^t A) = cond_2(A)^2.

Run: python scripts/posts/numerical-analysis/week-8/iterative.py   (writes iterative.json next to this file)
"""

import json
import math
import random
from fractions import Fraction
from pathlib import Path


# ---------- helpers ----------

def matvec(A, x):
    return [sum(a * b for a, b in zip(row, x)) for row in A]


def transpose(A):
    return [list(col) for col in zip(*A)]


def matmul(A, B):
    Bt = transpose(B)
    return [[sum(a * b for a, b in zip(row, col)) for col in Bt] for row in A]


def dot(x, y):
    return sum(a * b for a, b in zip(x, y))


def norm(x):
    return math.sqrt(dot(x, x))


def sub(x, y):
    return [a - b for a, b in zip(x, y)]


def lam_max_sym(A, iters=3000):
    """Largest eigenvalue of a symmetric positive (semi)definite matrix, by power iteration."""
    n = len(A)
    v = [1.0 + 0.1 * i for i in range(n)]
    lam = 0.0
    for _ in range(iters):
        w = matvec(A, v)
        lam_new = norm(w) / norm(v)
        v = [a / norm(w) for a in w]
        if abs(lam_new - lam) <= 1e-15 * lam_new:
            lam = lam_new
            break
        lam = lam_new
    return lam, v


def inverse_exact(A):
    """Gauss-Jordan in exact rational arithmetic."""
    n = len(A)
    M = [[Fraction(a) for a in row] + [Fraction(int(i == j)) for j in range(n)] for i, row in enumerate(A)]
    for k in range(n):
        p = next(i for i in range(k, n) if M[i][k] != 0)
        M[k], M[p] = M[p], M[k]
        piv = M[k][k]
        M[k] = [a / piv for a in M[k]]
        for i in range(n):
            if i != k and M[i][k] != 0:
                f = M[i][k]
                M[i] = [a - f * b for a, b in zip(M[i], M[k])]
    return [row[n:] for row in M]


def hilbert(n):
    return [[Fraction(1, i + j + 1) for j in range(n)] for i in range(n)]


def cond2_sym_pd(A_exact):
    """cond_2 = lam_max(A) * lam_max(A^{-1}) for SPD A, with A^{-1} computed exactly."""
    A = [[float(a) for a in row] for row in A_exact]
    Ainv = [[float(a) for a in row] for row in inverse_exact(A_exact)]
    lmax, _ = lam_max_sym(A)
    linv, vmin = lam_max_sym(Ainv)
    return lmax, 1.0 / linv, lmax * linv, vmin


# ---------- Part A ----------

def part_a():
    rows = []
    for n in range(2, 13):
        lmax, lmin, c, _ = cond2_sym_pd(hilbert(n))
        rows.append({"n": n, "lam_max": lmax, "lam_min": lmin, "cond2": c})
    # tiny residual, large error: H_8, exact x = (1, ..., 1), x_hat = x + (eigenvector of lam_min)
    n = 8
    H = [[float(a) for a in row] for row in hilbert(n)]
    lmax, lmin, c, v = cond2_sym_pd(hilbert(n))
    x = [1.0] * n
    b = matvec(H, x)
    scale = 0.5 * norm(x) / norm(v)
    x_hat = [a + scale * w for a, w in zip(x, v)]
    r = sub(b, matvec(H, x_hat))
    demo = {
        "n": n,
        "cond2": c,
        "rel_residual": norm(r) / norm(b),
        "rel_error": norm(sub(x, x_hat)) / norm(x),
    }
    demo["ratio"] = demo["rel_error"] / demo["rel_residual"]
    return {"hilbert": rows, "demo": demo}


# ---------- Part B ----------

def solve2(A, b):
    (a, bb), (c, d) = A
    det = a * d - bb * c
    return [(d * b[0] - bb * b[1]) / det, (-c * b[0] + a * b[1]) / det]


def part_b():
    A = [[6.0, 5.0], [7.0, 6.0]]  # det = 1, nearly singular in the sense of cond
    lmax, _ = lam_max_sym(matmul(transpose(A), A))
    Ainv = [[6.0, -5.0], [-7.0, 6.0]]
    linv, _ = lam_max_sym(matmul(transpose(Ainv), Ainv))
    cond = math.sqrt(lmax * linv)
    rng = random.Random(8)
    cases = []
    for x in ([1.0, 1.0], [1.0, -1.0]):
        b = matvec(A, x)
        ratios = []
        for _ in range(20000):
            t = rng.uniform(0, 2 * math.pi)
            db = [1e-3 * math.cos(t), 1e-3 * math.sin(t)]
            x2 = solve2(A, [b[0] + db[0], b[1] + db[1]])
            rel_err = norm(sub(x2, x)) / norm(x)
            rel_res = norm(db) / norm(b)
            ratios.append(rel_err / rel_res)
        cases.append({"x": x, "b": b, "ratio_min": min(ratios), "ratio_max": max(ratios)})
    return {"A": A, "cond2": cond, "inv_cond": 1.0 / cond, "samples": 20000, "cases": cases}


# ---------- Part C ----------

def part_c():
    A = [[5.0, 2.0, 0.0], [1.0, 6.0, 2.0], [0.0, 3.0, 7.0]]
    x = [1.0, -1.0, 2.0]
    b = matvec(A, x)
    D = [A[i][i] for i in range(3)]
    # iteration matrix I - CA with C = D^{-1}; its eigenvalues are computed from powers (rho = lim ||M^k||^{1/k})
    M = [[(1.0 if i == j else 0.0) - A[i][j] / D[i] for j in range(3)] for i in range(3)]
    P = [[float(i == j) for j in range(3)] for i in range(3)]
    for _ in range(400):
        P = matmul(P, M)
    rho_est = (max(abs(a) for row in P for a in row)) ** (1 / 400)
    # exact: M = -[[0, .4, 0], [1/6, 0, 1/3], [0, 3/7, 0]] -> eigenvalues 0, +-sqrt(.4/6 + 1/7)
    rho_exact = math.sqrt(0.4 / 6 + 1.0 / 7.0)
    norm_M = math.sqrt(lam_max_sym(matmul(transpose(M), M))[0])
    xs = [[0.0, 0.0, 0.0]]
    for _ in range(25):
        xm = xs[-1]
        r = sub(b, matvec(A, xm))
        xs.append([a + ri / d for a, ri, d in zip(xm, r, D)])
    rows = []
    for m in range(len(xs) - 2):
        err = norm(sub(x, xs[m]))
        err_next = norm(sub(x, xs[m + 1]))
        c_est = norm(sub(xs[m + 2], xs[m + 1])) / norm(sub(xs[m + 1], xs[m]))
        two = norm(sub(x, xs[m + 2])) / err
        rows.append({"m": m, "error": err, "error_ratio": err_next / err, "two_step_ratio": two, "c_estimate": c_est})
    return {"A": A, "x": x, "norm2_M": norm_M, "rho_exact": rho_exact, "rho_from_powers": rho_est, "rows": rows}


# ---------- Part D ----------

def poisson(n):
    return [[2.0 if i == j else (-1.0 if abs(i - j) == 1 else 0.0) for j in range(n)] for i in range(n)]


def tri_matvec(x):
    n = len(x)
    return [2 * x[i] - (x[i - 1] if i > 0 else 0.0) - (x[i + 1] if i < n - 1 else 0.0) for i in range(n)]


def run_jacobi(b, x_true, tol, max_it):
    n = len(b)
    x = [0.0] * n
    hist = [norm(sub(x, x_true)) / norm(x_true)]
    while hist[-1] > tol and len(hist) <= max_it:
        x = [(b[i] + (x[i - 1] if i > 0 else 0.0) + (x[i + 1] if i < n - 1 else 0.0)) / 2 for i in range(n)]
        hist.append(norm(sub(x, x_true)) / norm(x_true))
    return hist


def run_sor(b, x_true, w, tol, max_it):
    """w = 1 is Gauss-Seidel."""
    n = len(b)
    x = [0.0] * n
    hist = [norm(sub(x, x_true)) / norm(x_true)]
    while hist[-1] > tol and len(hist) <= max_it:
        for i in range(n):
            gs = (b[i] + (x[i - 1] if i > 0 else 0.0) + (x[i + 1] if i < n - 1 else 0.0)) / 2
            x[i] = w * gs + (1 - w) * x[i]
        hist.append(norm(sub(x, x_true)) / norm(x_true))
    return hist


def run_cg(mv, b, x_true, tol, max_it):
    n = len(b)
    x = [0.0] * n
    r = b[:]
    p = r[:]
    rr = dot(r, r)
    hist = [norm(sub(x, x_true)) / norm(x_true)]
    while hist[-1] > tol and len(hist) <= max_it:
        Ap = mv(p)
        alpha = rr / dot(p, Ap)
        x = [a + alpha * q for a, q in zip(x, p)]
        r = sub(b, mv(x))
        rr_new = dot(r, r)
        beta = rr_new / rr
        p = [a + beta * q for a, q in zip(r, p)]
        rr = rr_new
        hist.append(norm(sub(x, x_true)) / norm(x_true))
    return hist


def poisson_problem(n):
    # exact solution chosen first, with random entries, so it has a part along every eigenvector
    rng = random.Random(n)
    x_true = [rng.uniform(-1, 1) for _ in range(n)]
    return x_true, tri_matvec(x_true)


def part_d():
    tol = 1e-8
    counts = []
    for n in (8, 16, 32, 64):
        x_true, b = poisson_problem(n)
        h = 1.0 / (n + 1)
        rho_j = math.cos(math.pi * h)
        w = 2 / (1 + math.sqrt(1 - rho_j**2))
        lmax, lmin = 4 * math.sin(math.pi * n * h / 2) ** 2, 4 * math.sin(math.pi * h / 2) ** 2
        counts.append({
            "n": n,
            "cond2": lmax / lmin,
            "rho_jacobi": rho_j,
            "rho_gs": rho_j**2,
            "w_opt": w,
            "rho_sor": w - 1,
            "jacobi": len(run_jacobi(b, x_true, tol, 200000)) - 1,
            "gs": len(run_sor(b, x_true, 1.0, tol, 200000)) - 1,
            "sor": len(run_sor(b, x_true, w, tol, 200000)) - 1,
            "cg": len(run_cg(tri_matvec, b, x_true, tol, 10 * n)) - 1,
        })
    # check the closed-form eigenvalues against power iteration once
    lmax_check, _ = lam_max_sym(poisson(16))
    # SOR with other w, n = 32
    n = 32
    x_true, b = poisson_problem(n)
    w_scan = []
    for w in (0.5, 1.0, 1.5, 1.7, 1.8, 1.83, 1.9, 1.95, 1.99, 2.05):
        hist = run_sor(b, x_true, w, tol, 20000)
        w_scan.append({"w": w, "iterations": len(hist) - 1 if hist[-1] <= tol else None, "final_error": hist[-1]})
    # error history for the chart, n = 32
    rho_j = math.cos(math.pi / (n + 1))
    w = 2 / (1 + math.sqrt(1 - rho_j**2))
    hist = {
        "jacobi": run_jacobi(b, x_true, 1e-10, 600)[:601],
        "gs": run_sor(b, x_true, 1.0, 1e-10, 600)[:601],
        "sor": run_sor(b, x_true, w, 1e-10, 600)[:601],
        "cg": run_cg(tri_matvec, b, x_true, 1e-10, 600),
    }
    return {"tol": tol, "counts": counts, "lam_max_check_n16": lmax_check, "w_scan_n32": w_scan, "history_n32": hist}


# ---------- Part E ----------

def part_e():
    A = [[4.0, 1.0, 0.0, 0.0], [-3.0, 4.0, 1.0, 0.0], [0.0, -3.0, 4.0, 1.0], [0.0, 0.0, -3.0, 4.0]]
    x_true = [1.0, 2.0, -1.0, 1.0]
    b = matvec(A, x_true)
    At = transpose(A)
    N = matmul(At, A)
    hist_a = run_cg(lambda v: matvec(A, v), b, x_true, 1e-12, 40)
    hist_n = run_cg(lambda v: matvec(N, v), matvec(At, b), x_true, 1e-12, 40)
    Ninv = [[float(a) for a in row] for row in inverse_exact(N)]
    lmax, _ = lam_max_sym(N)
    linv, _ = lam_max_sym(Ninv)
    cond_n = lmax * linv
    return {
        "A": A,
        "cg_on_A_error": hist_a,
        "cg_on_normal_error": hist_n,
        "cond2_A": math.sqrt(cond_n),
        "cond2_AtA": cond_n,
    }


if __name__ == "__main__":
    result = {"part_a": part_a(), "part_b": part_b(), "part_c": part_c(), "part_d": part_d(), "part_e": part_e()}
    Path(__file__).with_name("iterative.json").write_text(json.dumps(result, indent=1), encoding="utf-8")
    a, bb, c, d, e = (result[k] for k in ("part_a", "part_b", "part_c", "part_d", "part_e"))
    print("A", [(r["n"], f'{r["cond2"]:.3e}') for r in a["hilbert"]])
    print("A demo", a["demo"])
    print("B", bb["cond2"], bb["inv_cond"], bb["cases"])
    print("C rho", c["rho_exact"], c["rho_from_powers"], "norm", c["norm2_M"])
    for r in c["rows"][:12]:
        print("  ", r)
    for r in d["counts"]:
        print("D", r)
    print("D lam check", d["lam_max_check_n16"], 4 * math.sin(math.pi * 16 / 17 / 2) ** 2)
    for r in d["w_scan_n32"]:
        print("  w", r)
    print("D hist lens", {k: len(v) for k, v in d["history_n32"].items()}, d["history_n32"]["cg"][-3:])
    print("E cond", e["cond2_A"], e["cond2_AtA"])
    print("E cg A", [f"{v:.2e}" for v in e["cg_on_A_error"][:12]], len(e["cg_on_A_error"]))
    print("E cg N", [f"{v:.2e}" for v in e["cg_on_normal_error"]])
