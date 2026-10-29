"""Experiment for the Lecture 8 post: convergence of Forward Euler, zero-stability, absolute stability.

Part 1: the convergence proof on a nonlinear problem. u' = cos(u), u(0) = 0, T = 2,
        exact solution u = 2 arctan(tanh(t/2)). f is Lipschitz with L = 1.
        Forward Euler error at the grid times against the proven bound
        e^{LT} T max|tau|, with tau computed from the exact solution.
        RK2 on the same problem, and the Lipschitz constant of its Phi
        (largest difference quotient on a grid of u values) against L + dt L^2 / 2.
Part 2: the linear problem u' = lam (u - sin t) + cos t, u(0) = 0, exact u = sin t,
        T = 2, with lam = -5 and lam = -1000. Forward Euler for many dt = 2/N
        around the limit dt = 2/1000. Error = max over the grid times.
        Also the proven bound T e^{|lam| T} max|tau| (as log10, it overflows a float).
Part 3: Forward Euler written as one lower triangular system A U = F for u' = lam u + g.
        A^{-1}[i][j] = dt (1 + dt lam)^{i-j} for i >= j, so the infinity norm is
        dt * sum_{m=0}^{N-1} |1 + dt lam|^m (largest row = last row).
        Compared with the bound T e^{|lam| T}. Reported as log10.

Standard library only.
Run: python scripts/posts/numerical-pde/lecture-8/experiment.py   (writes experiment.json next to this file)
"""

import json
import math
from pathlib import Path

T = 2.0


# ---------- part 1: nonlinear problem ----------
def f1(u):
    return math.cos(u)


def exact1(t):
    return 2 * math.atan(math.tanh(t / 2))


def fe_run(N):
    k = T / N
    u, err, tau = 0.0, 0.0, 0.0
    for n in range(N):
        t = n * k
        tau = max(tau, abs((exact1(t + k) - exact1(t)) / k - f1(exact1(t))))
        u = u + k * f1(u)
        err = max(err, abs(u - exact1(t + k)))
    return err, tau


def rk2_phi(u, k):
    return f1(u + 0.5 * k * f1(u))


def rk2_run(N):
    k = T / N
    u, err = 0.0, 0.0
    for n in range(N):
        u = u + k * rk2_phi(u, k)
        err = max(err, abs(u - exact1((n + 1) * k)))
    return err


def rk2_lipschitz(k):
    # largest |Phi(u) - Phi(v)| / |u - v| over close pairs on [-4, 4]
    best, h = 0.0, 1e-4
    for i in range(80001):
        u = -4 + i * 1e-4
        best = max(best, abs(rk2_phi(u + h, k) - rk2_phi(u, k)) / h)
    return best


part1 = {"N": [], "dt": [], "fe_err": [], "tau": [], "bound": [], "rk2_err": []}
L = 1.0
for N in (10, 20, 40, 80, 160, 320, 640, 1280):
    err, tau = fe_run(N)
    part1["N"].append(N)
    part1["dt"].append(T / N)
    part1["fe_err"].append(err)
    part1["tau"].append(tau)
    part1["bound"].append(math.exp(L * T) * T * tau)
    part1["rk2_err"].append(rk2_run(N))
part1["rk2_lip"] = {str(k): {"measured": rk2_lipschitz(k), "bound": L + 0.5 * k * L * L} for k in (0.5, 0.2, 0.05)}


# ---------- part 2: linear problem ----------
def run_linear(lam, N, final=False):
    k = T / N
    u, err = 0.0, 0.0
    for n in range(N):
        t = n * k
        u = u + k * (lam * (u - math.sin(t)) + math.cos(t))
        err = max(err, abs(u - math.sin(t + k)))
    return u if final else err


def tau_linear(lam, N):
    k = T / N
    worst = 0.0
    for n in range(N):
        t = n * k
        worst = max(worst, abs((math.sin(t + k) - math.sin(t)) / k - math.cos(t)))
    return worst  # the lam term vanishes on the exact solution


LAMS = (-5.0, -1000.0)
NS = (20, 40, 80, 160, 320, 640, 800, 900, 950, 990, 1000, 1010, 1050, 1100, 1280, 2560, 5120, 10240, 20480)
part2 = {"N": list(NS), "dt": [T / N for N in NS]}
for lam in LAMS:
    key = str(int(lam))
    part2[key] = [run_linear(lam, N) for N in NS]
    part2[key + "_final"] = [run_linear(lam, N, final=True) for N in NS]
    part2[key + "_z"] = [lam * T / N for N in NS]
    part2[key + "_log10_bound"] = [
        math.log10(T) + abs(lam) * T / math.log(10) + math.log10(tau_linear(lam, N)) for N in NS
    ]


# ---------- part 3: the triangular system ----------
def log10_inv_norm(lam, N):
    k = T / N
    g = abs(1 + k * lam)
    if g == 1:
        return math.log10(k * N)
    if g < 1:
        return math.log10(k * (1 - g**N) / (1 - g))
    # g > 1: k (g^N - 1) / (g - 1), in logs to avoid overflow
    return math.log10(k) + N * math.log10(g) + math.log10(1 - g ** (-N)) - math.log10(g - 1)


part3 = {"N": list(NS)}
for lam in LAMS:
    key = str(int(lam))
    part3[key] = [log10_inv_norm(lam, N) for N in NS]
    part3[key + "_log10_bound"] = math.log10(T) + abs(lam) * T / math.log(10)

out = {"part1": part1, "part2": part2, "part3": part3}
Path(__file__).with_name("experiment.json").write_text(json.dumps(out, indent=1), encoding="utf-8")

# ---------- print ----------
print("Part 1: u' = cos u, Forward Euler and RK2")
for i, N in enumerate(part1["N"]):
    e, b = part1["fe_err"][i], part1["bound"][i]
    o = "" if i == 0 else f"order {math.log2(part1['fe_err'][i-1] / e):.2f}  rk2 order {math.log2(part1['rk2_err'][i-1] / part1['rk2_err'][i]):.2f}"
    print(f"dt=1/{N // 2:<5} FE err {e:.3e}  tau {part1['tau'][i]:.3e}  bound {b:.3e}  bound/err {b / e:.2f}  err/dt {e / (T / N):.4f}  RK2 {part1['rk2_err'][i]:.3e}  {o}")
for k, v in part1["rk2_lip"].items():
    print(f"RK2 Lipschitz dt={k}: measured {v['measured']:.4f}  bound {v['bound']:.4f}")
print("Part 2/3: linear problem")
for i, N in enumerate(NS):
    row = f"N={N:<6} dt={T / N:.6f}"
    for lam in LAMS:
        key = str(int(lam))
        row += f" | lam={key}: z={part2[key + '_z'][i]:.3f} |1+z|={abs(1 + part2[key + '_z'][i]):.3f} err {part2[key][i]:.3e} bound 1e{part2[key + '_log10_bound'][i]:.1f} ||Ainv|| 1e{part3[key][i]:.2f}"
    print(row)
print(f"lam=-1000, N=800: U^N = {part2['-1000_final'][NS.index(800)]:.4e}, sin 2 = {math.sin(2):.6f}")
for lam in LAMS:
    print(f"lam={lam}: T e^(|lam|T) = 1e{part3[str(int(lam)) + '_log10_bound']:.2f}")
