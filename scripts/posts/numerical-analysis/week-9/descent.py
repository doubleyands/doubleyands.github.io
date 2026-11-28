"""Experiments for "Numerical Analysis, Week 9": gradient methods and the construction of CG.

Python standard library only. Writes descent.json next to this file.

A. Gradient method with a fixed step alpha on a 4x4 SPD matrix with eigenvalues 1, 2.5, 5, 8:
   observed rate against gamma = max |1 - alpha*lambda_j|.
B. Iteration counts to reduce the error by eps = 1e-6, for diagonal matrices with cond = 10 ... 1e4:
   fixed alpha = 1/lambda_max (worst start and random start), against cond*log(1/eps);
   steepest descent with line search and CG for comparison. Also the error history for cond = 100.
C. Steepest descent in 2D: the descent interval 0 < alpha < 2 g.g / g.Ag, and the
   new gradient being orthogonal to the last direction (zigzag).
D. CG built from the "minimization" formulas on a 6x6 SPD matrix: the two forms of alpha and beta,
   orthogonality of residuals, A-orthogonality of directions, and the Krylov span lemma.

Run: python scripts/posts/numerical-analysis/week-9/descent.py
"""

import json
import math
import random
from pathlib import Path

# ---------- small linear algebra ----------

def dot(u, v):
    return sum(a * b for a, b in zip(u, v))


def norm(u):
    return math.sqrt(dot(u, u))


def matvec(A, v):
    return [dot(row, v) for row in A]


def add(u, v, s=1.0):
    return [a + s * b for a, b in zip(u, v)]


def sym_eig(A, sweeps=60):
    """Eigenvalues of a symmetric matrix by cyclic Jacobi rotations."""
    n = len(A)
    M = [row[:] for row in A]
    for _ in range(sweeps):
        off = sum(M[i][j] ** 2 for i in range(n) for j in range(n) if i != j)
        if off < 1e-30:
            break
        for p in range(n):
            for q in range(p + 1, n):
                if abs(M[p][q]) < 1e-300:
                    continue
                th = (M[q][q] - M[p][p]) / (2 * M[p][q])
                t = (1 if th >= 0 else -1) / (abs(th) + math.sqrt(th * th + 1))
                c = 1 / math.sqrt(t * t + 1)
                s = t * c
                for k in range(n):
                    mkp, mkq = M[k][p], M[k][q]
                    M[k][p], M[k][q] = c * mkp - s * mkq, s * mkp + c * mkq
                for k in range(n):
                    mpk, mqk = M[p][k], M[q][k]
                    M[p][k], M[q][k] = c * mpk - s * mqk, s * mpk + c * mqk
    return sorted(M[i][i] for i in range(n))


def householder_spd(eigs, v):
    """A = Q diag(eigs) Q with the reflection Q = I - 2 v v^T / v^T v (Q = Q^T = Q^-1)."""
    n = len(eigs)
    vv = dot(v, v)
    Q = [[(1.0 if i == j else 0.0) - 2 * v[i] * v[j] / vv for j in range(n)] for i in range(n)]
    return [[sum(Q[i][k] * eigs[k] * Q[k][j] for k in range(n)) for j in range(n)] for i in range(n)]


def sci(x, d=3):
    return float(f"{x:.{d}g}")


out = {}

# ---------- A. fixed-step gradient method ----------
eigs_a = [1.0, 2.5, 5.0, 8.0]
A4 = householder_spd(eigs_a, [1.0, 2.0, -1.0, 1.0])
x4 = [1.0, -2.0, 0.5, 3.0]
b4 = matvec(A4, x4)
check_eigs = sym_eig(A4)
part_a = {"matrix": [[round(a, 6) for a in row] for row in A4], "eigs_check": [round(e, 12) for e in check_eigs],
          "x": x4, "rows": []}
for name, alpha in [("1/lmax", 1 / 8), ("2/(lmin+lmax)", 2 / 9), ("2/lmax", 2 / 8), ("0.26", 0.26)]:
    gamma = max(abs(1 - alpha * l) for l in eigs_a)
    xi = [0.0] * 4
    errs = [norm(add(x4, xi, -1))]
    for k in range(200):
        g = add(matvec(A4, xi), b4, -1)
        xi = add(xi, g, -alpha)
        errs.append(norm(add(x4, xi, -1)))
    rate = (errs[100] / errs[80]) ** (1 / 20)
    part_a["rows"].append({"alpha": name, "alpha_value": alpha, "gamma": gamma, "e0": errs[0],
                           "e10": errs[10], "e50": errs[50], "e100": errs[100], "observed_rate": rate})
out["part_a"] = part_a

# ---------- B. iteration counts against the condition number ----------
EPS = 1e-6
N_B = 100
rng = random.Random(9)
e_rand = [rng.uniform(-1, 1) for _ in range(N_B)]


def diag_problem(kappa):
    return [1 + (kappa - 1) * i / (N_B - 1) for i in range(N_B)]


def run_fixed(lam, e0, alpha, eps, kmax, hist=None):
    e = e0[:]
    n0 = norm(e)
    for k in range(1, kmax + 1):
        e = [ei * (1 - alpha * li) for ei, li in zip(e, lam)]
        r = norm(e) / n0
        if hist is not None:
            hist.append(r)
        if r <= eps:
            return k
    return None


def run_sd(lam, e0, eps, kmax, hist=None):
    """Steepest descent with exact line search. Error e = x - xi, gradient g = A xi - b = -A e."""
    e = e0[:]
    n0 = norm(e)
    for k in range(1, kmax + 1):
        g = [-li * ei for li, ei in zip(lam, e)]
        Ag = [li * gi for li, gi in zip(lam, g)]
        alpha = dot(g, g) / dot(g, Ag)
        e = [ei + alpha * gi for ei, gi in zip(e, g)]  # xi_new = xi - alpha g  ->  e_new = e + alpha g
        r = norm(e) / n0
        if hist is not None:
            hist.append(r)
        if r <= eps:
            return k
    return None


def run_cg(lam, e0, eps, kmax, hist=None):
    """CG with the week-9 sign: r = A xi - b = -A e, d0 = -r0."""
    e = e0[:]
    n0 = norm(e)
    r = [-li * ei for li, ei in zip(lam, e)]
    d = [-ri for ri in r]
    for k in range(1, kmax + 1):
        Ad = [li * di for li, di in zip(lam, d)]
        alpha = -dot(r, d) / dot(d, Ad)
        e = [ei - alpha * di for ei, di in zip(e, d)]
        r = add(r, Ad, alpha)
        beta = dot(r, Ad) / dot(d, Ad)
        d = add([-ri for ri in r], d, beta)
        rr = norm(e) / n0
        if hist is not None:
            hist.append(rr)
        if rr <= eps:
            return k
    return None


rows_b = []
for kappa in [10, 100, 1000, 10000]:
    lam = diag_problem(kappa)
    worst = [1.0] + [0.0] * (N_B - 1)
    exact_min = math.ceil(math.log(EPS) / math.log(1 - 1 / kappa))
    rows_b.append({
        "kappa": kappa,
        "fixed_worst": run_fixed(lam, worst, 1 / kappa, EPS, 10**6),
        "formula_min": exact_min,
        "kappa_log": kappa * math.log(1 / EPS),
        "fixed_random": run_fixed(lam, e_rand, 1 / kappa, EPS, 10**6),
        "fixed_best_random": run_fixed(lam, e_rand, 2 / (1 + kappa), EPS, 10**6),
        "sd_random": run_sd(lam, e_rand, EPS, 10**6),
        "cg_random": run_cg(lam, e_rand, EPS, 10**4),
        "sqrt_kappa": math.sqrt(kappa),
    })
lam100 = diag_problem(100)
H = 400
hist = {"fixed_1_over_lmax": [1.0], "fixed_best": [1.0], "sd": [1.0], "cg": [1.0]}
run_fixed(lam100, e_rand, 1 / 100, 0.0, H, hist["fixed_1_over_lmax"])
run_fixed(lam100, e_rand, 2 / 101, 0.0, H, hist["fixed_best"])
run_sd(lam100, e_rand, 0.0, H, hist["sd"])
run_cg(lam100, e_rand, 1e-14, H, hist["cg"])
out["part_b"] = {"n": N_B, "eps": EPS, "rows": rows_b, "history_kappa100": hist}

# ---------- C. steepest descent in 2D ----------
th = 0.5
co, si = math.cos(th), math.sin(th)
l1, l2 = 1.0, 10.0
A2 = [[l1 * co * co + l2 * si * si, (l1 - l2) * co * si], [(l1 - l2) * co * si, l1 * si * si + l2 * co * co]]
x2 = [2.0, -1.0]
b2 = matvec(A2, x2)
f = lambda v: 0.5 * dot(v, matvec(A2, v)) - dot(b2, v)
fstar = f(x2)
xi = [-3.0, 2.0]
g = add(matvec(A2, xi), b2, -1)
alpha_star = dot(g, g) / dot(g, matvec(A2, g))
scan = []
for frac in [0.25, 0.5, 1.0, 1.5, 1.9, 2.0, 2.1, 3.0]:
    a = frac * alpha_star
    scan.append({"alpha_over_opt": frac, "f_change": f(add(xi, g, -a)) - f(xi)})
def sd_steps(start, count=12):
    xi = start[:]
    rows, prev_d = [], None
    ea = lambda v: math.sqrt(dot(add(v, x2, -1), matvec(A2, add(v, x2, -1))))
    for k in range(count):
        g = add(matvec(A2, xi), b2, -1)
        d = [-gi for gi in g]
        alpha = -dot(g, d) / dot(d, matvec(A2, d))
        new = add(xi, d, alpha)
        gnew = add(matvec(A2, new), b2, -1)
        rows.append({
            "k": k, "f_minus_fstar": f(xi) - fstar, "energy_err": ea(xi), "ratio": ea(new) / ea(xi),
            "gnew_dot_d_rel": abs(dot(gnew, d)) / (norm(gnew) * norm(d)),
            "angle_to_prev_deg": (math.degrees(math.acos(max(-1, min(1, dot(d, prev_d) / (norm(d) * norm(prev_d))))))
                                  if prev_d else None),
        })
        prev_d = d
        xi = new
    return rows


v1, v2 = [co, si], [-si, co]  # eigenvectors for l1, l2
starts = {
    "easy": [-3.0, 2.0],
    # error e0 = l2 v1 + l1 v2: the worst start for steepest descent
    "worst": add(x2, add([l2 * c for c in v1], [l1 * c for c in v2]), -1),
}
steps = {name: sd_steps(st) for name, st in starts.items()}
out["part_c"] = {"A": A2, "eigs": [l1, l2], "x": x2, "starts": starts, "alpha_opt_first": alpha_star,
                 "alpha_limit_first": 2 * alpha_star, "scan": scan, "steps": steps,
                 "bound": (l2 / l1 - 1) / (l2 / l1 + 1)}

# ---------- D. CG from the minimization formulas ----------
n = 6
A6 = [[1 / (1 + abs(i - j)) for j in range(n)] for i in range(n)]
b6 = [1.0, 0.0, 2.0, -1.0, 0.0, 1.0]
eig6 = sym_eig(A6)
xi = [0.0] * n
r = add(matvec(A6, xi), b6, -1)  # r = A xi - b
d = [-ri for ri in r]
R, D, alphas, alphas2, betas, betas2, fvals = [r], [d], [], [], [], [], []
fA = lambda v: 0.5 * dot(v, matvec(A6, v)) - dot(b6, v)
fvals.append(fA(xi))
for k in range(n):
    Ad = matvec(A6, d)
    alpha = -dot(r, d) / dot(d, Ad)
    alphas.append(alpha)
    alphas2.append(dot(r, r) / dot(d, Ad))
    xi = add(xi, d, alpha)
    rn = add(matvec(A6, xi), b6, -1)
    fvals.append(fA(xi))
    if k == n - 1:
        R.append(rn)
        break
    beta = dot(rn, Ad) / dot(d, Ad)
    betas.append(beta)
    betas2.append(dot(rn, rn) / dot(r, r))
    d = add([-x for x in rn], d, beta)
    r = rn
    R.append(r)
    D.append(d)
max_rr = max(abs(dot(R[i], R[j])) / (norm(R[i]) * norm(R[j])) for i in range(n) for j in range(n) if i != j)
max_dad = max(abs(dot(D[i], matvec(A6, D[j]))) / math.sqrt(dot(D[i], matvec(A6, D[i])) * dot(D[j], matvec(A6, D[j])))
              for i in range(n) for j in range(n) if i != j)
max_rd = max(abs(dot(R[k + 1], D[k])) / (norm(R[k + 1]) * norm(D[k])) for k in range(n - 1))
# Krylov lemma: distance of d^k from span{r0, A r0, ..., A^k r0} (orthonormal basis by Gram-Schmidt, twice)
basis, kry, v = [], [], R[0][:]
for k in range(n):
    w = v[:]
    for _ in range(2):
        for q in basis:
            w = add(w, q, -dot(w, q))
    basis.append([wi / norm(w) for wi in w])
    dk = D[k]
    res = dk[:]
    for q in basis:
        res = add(res, q, -dot(res, q))
    kry.append(norm(res) / norm(dk))
    v = matvec(A6, v)
x6 = xi
final_res = norm(add(matvec(A6, x6), b6, -1)) / norm(b6)
out["part_d"] = {
    "A_rule": "a_ij = 1/(1+|i-j|), n = 6", "b": b6, "eigs": eig6, "cond2": eig6[-1] / eig6[0],
    "alpha_minimization": alphas, "alpha_week8": alphas2, "beta_minimization": betas, "beta_week8": betas2,
    "max_rel_alpha_diff": max(abs(a - c) / abs(a) for a, c in zip(alphas, alphas2)),
    "max_rel_beta_diff": max(abs(a - c) / abs(a) for a, c in zip(betas, betas2)),
    "max_rel_r_dot_r": max_rr, "max_rel_d_A_d": max_dad, "max_rel_rnext_dot_d": max_rd,
    "krylov_distance": kry, "f_values": fvals, "x_after_6": x6, "rel_residual_after_6": final_res,
}

Path(__file__).with_name("descent.json").write_text(json.dumps(out, indent=1), encoding="utf-8")

# ---------- short report ----------
print("A. eigs check", [round(e, 10) for e in check_eigs])
for row in part_a["rows"]:
    print(f"  alpha {row['alpha']:>14} = {row['alpha_value']:.4f}  gamma {row['gamma']:.4f}  "
          f"e10 {row['e10']:.3g}  e50 {row['e50']:.3g}  e100 {row['e100']:.3g}  rate {row['observed_rate']:.4f}")
print("B.")
for row in rows_b:
    print("  ", {k: (round(v, 1) if isinstance(v, float) else v) for k, v in row.items()})
print("  hist k=400:", {k: f"{v[-1]:.3g}" for k, v in hist.items()}, "cg len", len(hist["cg"]) - 1)
pc = out["part_c"]
print("C. alpha* %.5f, limit %.5f, bound %.4f" % (pc["alpha_opt_first"], pc["alpha_limit_first"], pc["bound"]))
for s in scan:
    print("   ", s)
for name, rows in steps.items():
    print("   start", name, [round(c, 4) for c in starts[name]])
    for s in rows[:4] + rows[-1:]:
        print("   ", {k: (f"{v:.4g}" if isinstance(v, float) else v) for k, v in s.items()})
pd = out["part_d"]
print("D. eigs", [round(e, 4) for e in eig6], "cond", round(pd["cond2"], 4))
for key in ["alpha_minimization", "beta_minimization", "krylov_distance", "f_values"]:
    print("  ", key, [f"{v:.4g}" for v in pd[key]])
for key in ["max_rel_alpha_diff", "max_rel_beta_diff", "max_rel_r_dot_r", "max_rel_d_A_d", "max_rel_rnext_dot_d",
            "rel_residual_after_6"]:
    print("  ", key, f"{pd[key]:.3g}")
print("   x", [round(v, 6) for v in x6])
