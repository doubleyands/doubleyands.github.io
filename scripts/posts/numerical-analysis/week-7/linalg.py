"""Week 7 experiments: matrix powers, the Neumann series, Gaussian elimination.

Python standard library only. Small dense matrices are lists of rows.

Part A  ||A^m||_2 for four matrices: rho(A) < 1 means A^m -> 0, but the norm can grow first.
Part B  Neumann series I + C + C^2 + ... -> (I - C)^{-1}, and the bound 1 / (1 - ||C||).
Part C  Operation count of Gaussian elimination (LU) against 2n^3/3, and the cost of one
        extra right-hand side (forward + back substitution) against 2n^2.
Part D  Pivoting: elimination in 4-digit arithmetic, with and without row exchanges,
        and a tiny pivot in double precision.
Part E  Residual vs. error for a nearly singular matrix.

Run: python scripts/posts/numerical-analysis/week-7/linalg.py   (writes linalg.json next to this file)
"""

import json
import math
import random
from pathlib import Path


# ---------- small matrix helpers ----------

def matmul(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]


def matvec(A, x):
    return [sum(a * b for a, b in zip(row, x)) for row in A]


def transpose(A):
    return [list(r) for r in zip(*A)]


def eye(n):
    return [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]


def sub(A, B):
    return [[a - b for a, b in zip(ra, rb)] for ra, rb in zip(A, B)]


def norm2_vec(x):
    return math.sqrt(sum(v * v for v in x))


def norm2(A, iters=500):
    """||A||_2 = sqrt(largest eigenvalue of A^T A), by power iteration on A^T A."""
    M = matmul(transpose(A), A)
    n = len(M)
    v = [1.0 / math.sqrt(n) + 0.01 * i for i in range(n)]
    lam = 0.0
    for _ in range(iters):
        w = matvec(M, v)
        s = norm2_vec(w)
        if s == 0.0:
            return 0.0
        v = [t / s for t in w]
        lam = s
    return math.sqrt(lam)


# ---------- Part A: matrix powers ----------

def part_a():
    s2 = math.sqrt(2)
    mats = {
        # symmetric, eigenvalues 0.5 and 0.5 +- 0.3*sqrt(2): rho = 0.924
        "symmetric": [[0.5, 0.3, 0.0], [0.3, 0.5, 0.3], [0.0, 0.3, 0.5]],
        # 4x4 Jordan block with lambda = 0.8: rho = 0.8
        "jordan": [[0.8, 1, 0, 0], [0, 0.8, 1, 0], [0, 0, 0.8, 1], [0, 0, 0, 0.8]],
        # upper triangular, eigenvalues 0.5 and 0.6, one big entry: rho = 0.6
        "triangular": [[0.5, 4.0], [0.0, 0.6]],
        # symmetric, eigenvalues 1.05 and 0.75: rho = 1.05
        "growing": [[0.9, 0.15], [0.15, 0.9]],
    }
    rho = {"symmetric": 0.5 + 0.3 * s2, "jordan": 0.8, "triangular": 0.6, "growing": 1.05}
    ms = list(range(0, 121))
    out = {"m": ms, "rho": rho, "norm1": {}, "series": {}}
    for name, A in mats.items():
        P = eye(len(A))
        vals = []
        for m in ms:
            if m > 0:
                P = matmul(P, A)
            vals.append(norm2(P) if m > 0 else 1.0)
        out["series"][name] = vals
        out["norm1"][name] = vals[1]
    return out


# ---------- Part B: Neumann series ----------

def inverse(A):
    """Gauss-Jordan with partial pivoting (for checking only)."""
    n = len(A)
    M = [list(map(float, A[i])) + eye(n)[i] for i in range(n)]
    for k in range(n):
        p = max(range(k, n), key=lambda i: abs(M[i][k]))
        M[k], M[p] = M[p], M[k]
        piv = M[k][k]
        M[k] = [v / piv for v in M[k]]
        for i in range(n):
            if i != k:
                f = M[i][k]
                M[i] = [a - f * b for a, b in zip(M[i], M[k])]
    return [row[n:] for row in M]


def part_b():
    cases = {
        "small_norm": [[0.2, 0.1, 0.3], [0.1, 0.4, 0.1], [0.0, 0.2, 0.3]],
        "big_norm": [[0.4, 3.0], [0.0, 0.5]],  # rho = 0.5 but ||.||_2 > 1
    }
    out = {}
    for name, C in cases.items():
        n = len(C)
        exact = inverse(sub(eye(n), C))
        c = norm2(C)
        S = eye(n)
        P = eye(n)
        rows = []
        for m in range(1, 41):
            P = matmul(P, C)
            S = [[a + b for a, b in zip(ra, rb)] for ra, rb in zip(S, P)]
            if m in (1, 2, 5, 10, 20, 40):
                rows.append({"m": m, "error": norm2(sub(exact, S)), "tail_bound": c ** (m + 1) / (1 - c) if c < 1 else None})
        out[name] = {
            "norm2_C": c,
            "norm2_inverse": norm2(exact),
            "bound": 1 / (1 - c) if c < 1 else None,
            "exact_inverse": exact,
            "rows": rows,
        }
    return out


# ---------- Part C: operation count ----------

def lu_count(A):
    """Plain elimination without row exchanges; returns L, U and the number of
    multiplications, divisions, additions and subtractions."""
    n = len(A)
    U = [row[:] for row in A]
    L = eye(n)
    ops = 0
    for k in range(n - 1):
        for i in range(k + 1, n):
            m = U[i][k] / U[k][k]
            ops += 1
            L[i][k] = m
            U[i][k] = 0.0
            for j in range(k + 1, n):
                U[i][j] -= m * U[k][j]
                ops += 2
    return L, U, ops


def solve_count(L, U, b):
    n = len(b)
    ops = 0
    y = [0.0] * n
    for i in range(n):
        s = b[i]
        for j in range(i):
            s -= L[i][j] * y[j]
            ops += 2
        y[i] = s  # L has ones on the diagonal
    x = [0.0] * n
    for i in reversed(range(n)):
        s = y[i]
        for j in range(i + 1, n):
            s -= U[i][j] * x[j]
            ops += 2
        x[i] = s / U[i][i]
        ops += 1
    return x, ops


def part_c():
    rng = random.Random(7)
    rows = []
    for n in (10, 20, 50, 100, 200):
        # diagonally dominant random matrix: no pivoting needed
        A = [[rng.uniform(-1, 1) + (n if i == j else 0) for j in range(n)] for i in range(n)]
        x_true = [rng.uniform(-1, 1) for _ in range(n)]
        b = matvec(A, x_true)
        L, U, ops_lu = lu_count(A)
        x, ops_solve = solve_count(L, U, b)
        err = max(abs(a - c) for a, c in zip(x, x_true))
        rows.append({
            "n": n,
            "lu_ops": ops_lu,
            "lu_ratio": ops_lu / (2 * n ** 3 / 3),
            "solve_ops": ops_solve,
            "solve_ratio": ops_solve / (2 * n ** 2),
            "max_error": err,
            "rhs_1000_reuse": ops_lu + 1000 * ops_solve,
            "rhs_1000_redo": 1000 * (ops_lu + ops_solve),
        })
    # small check: L U = A for a 3x3 example
    A3 = [[4.0, 2.0, 1.0], [2.0, 5.0, 2.0], [1.0, 2.0, 6.0]]
    L3, U3, _ = lu_count(A3)
    det = U3[0][0] * U3[1][1] * U3[2][2]
    return {"rows": rows, "example": {"A": A3, "L": L3, "U": U3, "det": det, "LU_minus_A": max(abs(v) for r in sub(matmul(L3, U3), A3) for v in r)}}


# ---------- Part D: pivoting ----------

def fl(v, digits=4):
    """Round to `digits` significant digits (simulates short arithmetic)."""
    if v == 0.0:
        return 0.0
    e = math.floor(math.log10(abs(v)))
    q = 10.0 ** (digits - 1 - e)
    return round(v * q) / q


def gauss(A, b, pivot, r=lambda v: v):
    n = len(A)
    M = [[r(v) for v in A[i]] + [r(b[i])] for i in range(n)]
    swaps = []
    for k in range(n - 1):
        if pivot:
            p = max(range(k, n), key=lambda i: abs(M[i][k]))
            if p != k:
                M[k], M[p] = M[p], M[k]
                swaps.append((k + 1, p + 1))
        for i in range(k + 1, n):
            m = r(M[i][k] / M[k][k])
            M[i][k] = 0.0
            for j in range(k + 1, n + 1):
                M[i][j] = r(M[i][j] - r(m * M[k][j]))
    x = [0.0] * n
    for i in reversed(range(n)):
        s = M[i][n]
        for j in range(i + 1, n):
            s = r(s - r(M[i][j] * x[j]))
        x[i] = r(s / M[i][i])
    return x, swaps


def part_d():
    # 4-digit arithmetic on a 3x3 system whose first pivot is small
    A = [[0.003, 1.2, 2.5], [1.6, -0.8, 1.1], [2.1, 1.7, -0.9]]
    x_true = [1.0, 2.0, 3.0]
    b = [sum(a * x for a, x in zip(row, x_true)) for row in A]
    no_piv, _ = gauss(A, b, pivot=False, r=fl)
    piv, swaps = gauss(A, b, pivot=True, r=fl)
    exact_dp, _ = gauss(A, b, pivot=True)
    # double precision with a tiny pivot
    eps = 1e-17
    A2 = [[eps, 1.0], [1.0, 2.0]]
    b2 = [1.0, 3.0]
    # exact solution: x1 = 1 / (1 - 2 eps), x2 = (1 - eps) / (1 - 2 eps)  ~ (1, 1)
    tiny_no, _ = gauss(A2, b2, pivot=False)
    tiny_piv, _ = gauss(A2, b2, pivot=True)
    return {
        "A": A, "b": b, "x_true": x_true, "double_precision": exact_dp,
        "four_digit_no_pivot": no_piv, "four_digit_pivot": piv, "swaps": swaps,
        "max_error_no_pivot": max(abs(a - c) for a, c in zip(no_piv, x_true)),
        "max_error_pivot": max(abs(a - c) for a, c in zip(piv, x_true)),
        "first_multipliers_no_pivot": [A[1][0] / A[0][0], A[2][0] / A[0][0]],
        "tiny": {"eps": eps, "no_pivot": tiny_no, "pivot": tiny_piv},
    }


# ---------- Part E: residual vs error ----------

def part_e():
    A = [[2.0, 3.0], [4.0, 6.0001]]
    x = [1.0, 1.0]
    b = matvec(A, x)
    x_a = [2.5, 0.0]
    r = [bi - v for bi, v in zip(b, matvec(A, x_a))]
    e = [xi - ai for xi, ai in zip(x, x_a)]
    Ae = matvec(A, e)
    # eigenvalues of A (2x2)
    tr = A[0][0] + A[1][1]
    det = A[0][0] * A[1][1] - A[0][1] * A[1][0]
    disc = math.sqrt(tr * tr / 4 - det)
    return {
        "A": A, "b": b, "x": x, "x_A": x_a, "residual": r, "error": e, "A_times_error": Ae,
        "rel_residual": norm2_vec(r) / norm2_vec(b), "rel_error": norm2_vec(e) / norm2_vec(x),
        "det": det, "eigenvalues": [tr / 2 + disc, tr / 2 - disc],
    }


if __name__ == "__main__":
    result = {"part_a": part_a(), "part_b": part_b(), "part_c": part_c(), "part_d": part_d(), "part_e": part_e()}
    Path(__file__).with_name("linalg.json").write_text(json.dumps(result, indent=1), encoding="utf-8")

    a = result["part_a"]
    for name, vals in a["series"].items():
        k = max(range(len(vals)), key=lambda i: vals[i])
        print(f"A {name:10s} rho={a['rho'][name]:.4f} ||A||={vals[1]:.4f} peak m={k} {vals[k]:.4g}  "
              + "  ".join(f"m={m}:{vals[m]:.3g}" for m in (5, 10, 20, 40, 80, 120)))
    for name, d in result["part_b"].items():
        print("B", name, f"||C||={d['norm2_C']:.4f} ||inv||={d['norm2_inverse']:.4f} bound={d['bound']}",
              [(r["m"], f"{r['error']:.2e}", None if r["tail_bound"] is None else f"{r['tail_bound']:.2e}") for r in d["rows"]])
    for r in result["part_c"]["rows"]:
        print("C", r)
    print("C example", result["part_c"]["example"])
    d = result["part_d"]
    print("D b", d["b"], "dp", d["double_precision"])
    print("D no pivot", d["four_digit_no_pivot"], d["max_error_no_pivot"], "mult", d["first_multipliers_no_pivot"])
    print("D pivot", d["four_digit_pivot"], d["max_error_pivot"], d["swaps"])
    print("D tiny", d["tiny"])
    print("E", result["part_e"])
