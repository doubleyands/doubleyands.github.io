"""Experiment for the Lecture 7 post: Runge-Kutta, multistep methods, Gibbs, stability.

Test problem (parts 2-4): the logistic equation u' = u (1 - u), u(0) = 0.2,
exact solution u = 1 / (1 + 4 e^{-t}). We integrate to T = 2.

Part 1: one step on u' = -2u with k = 0.5 (Forward Euler, RK2, RK4, exact);
        RK4 with f = f(t) only is Simpson's rule; the RK2 and RK4 growth factors
        for u' = lambda u against the Taylor polynomials of e^z.
Part 2: global error at T = 2 for Forward Euler, RK2, RK4, leapfrog (RK2 start),
        AB2 (RK2 start), AB4 (RK4 start), with k = 1/10 ... 1/5120 (halving).
Part 3: starting values. Leapfrog started with Forward Euler / RK2 / exact value,
        AB4 started with Forward Euler / RK4.
Part 4: cost in evaluations of f. RK4 uses 4 per step, AB4 uses 1 per step
        (plus 12 for the three RK4 start steps).
Part 5: Gibbs phenomenon. Partial sums of the square-wave Fourier series.
Part 6: stability and conditioning of the 1D Dirichlet matrix A (u'' on (0,1)).
        2-norm of A^{-1} and condition number from the known eigenvalues,
        infinity-norm of A^{-1} by solving A w = 1 (A^{-1} has one sign).

Standard library only.
Run: python scripts/posts/numerical-pde/lecture-7/experiment.py   (writes experiment.json next to this file)
"""

import json
import math
from pathlib import Path

U0, T = 0.2, 2.0


def f(u, t):
    return u * (1 - u)


def exact(t):
    return 1 / (1 + 4 * math.exp(-t))


# ---------- one-step methods ----------
def fe_step(f, u, t, k):
    return u + k * f(u, t)


def rk2_step(f, u, t, k):
    ustar = u + k / 2 * f(u, t)
    return u + k * f(ustar, t + k / 2)


def rk4_step(f, u, t, k):
    f1 = f(u, t)
    f2 = f(u + k / 2 * f1, t + k / 2)
    f3 = f(u + k / 2 * f2, t + k / 2)
    f4 = f(u + k * f3, t + k)
    return u + k / 6 * (f1 + 2 * f2 + 2 * f3 + f4)


def one_step_run(step, N):
    k = T / N
    u = U0
    for n in range(N):
        u = step(f, u, n * k, k)
    return u


# ---------- multistep methods ----------
def leapfrog_run(N, start):
    k = T / N
    prev = U0
    if start == "Forward Euler":
        cur = fe_step(f, U0, 0.0, k)
    elif start == "RK2":
        cur = rk2_step(f, U0, 0.0, k)
    else:
        cur = exact(k)
    for n in range(1, N):
        prev, cur = cur, prev + 2 * k * f(cur, n * k)
    return cur


AB = {
    2: ([-1, 3], 2),
    4: ([-9, 37, -59, 55], 24),
}


def ab_run(r, N, start_step):
    """Adams-Bashforth with r steps. The first r-1 values come from start_step."""
    k = T / N
    beta, den = AB[r]
    us = [U0]
    for n in range(r - 1):
        us.append(start_step(f, us[-1], n * k, k))
    fs = [f(u, n * k) for n, u in enumerate(us)]
    for n in range(r - 1, N):
        new = us[-1] + k / den * sum(b * fv for b, fv in zip(beta, fs[-r:]))
        us.append(new)
        fs.append(f(new, (n + 1) * k))
    assert len(us) == N + 1
    return us[-1]


def orders(errs):
    return [None] + [math.log2(a / b) for a, b in zip(errs, errs[1:])]


out = {}

# ---------- Part 1 ----------
lam, k = -2.0, 0.5
g = lambda u, t: lam * u
fe = fe_step(g, 1.0, 0.0, k)
rk2 = rk2_step(g, 1.0, 0.0, k)
rk4 = rk4_step(g, 1.0, 0.0, k)
ex = math.exp(lam * k)
print("Part 1: u' = -2u, k = 0.5, one step from 1")
print(f"  FE {fe:.6f}  RK2 {rk2:.6f}  RK4 {rk4:.6f}  exact {ex:.6f}")
h = lambda u, t: math.exp(t)
rk4_quad = rk4_step(h, 0.0, 0.0, 1.0)
simpson = (1 / 6) * (math.exp(0) + 4 * math.exp(0.5) + math.exp(1))
print(f"  f = e^t on [0,1]: RK4 {rk4_quad:.15f}  Simpson {simpson:.15f}  exact {math.e - 1:.15f}")
worst = 0.0
for z in [-3.0, -1.0, -0.3, 0.1, 0.7, 1.5]:
    g2 = lambda u, t, z=z: z * u
    r2 = rk2_step(g2, 1.0, 0.0, 1.0) - (1 + z + z * z / 2)
    r4 = rk4_step(g2, 1.0, 0.0, 1.0) - (1 + z + z**2 / 2 + z**3 / 6 + z**4 / 24)
    worst = max(worst, abs(r2), abs(r4))
print(f"  growth factor vs Taylor polynomial, worst difference {worst:.1e}")
out["part1"] = {"FE": fe, "RK2": rk2, "RK4": rk4, "exact": ex,
                "rk4_quad": rk4_quad, "simpson": simpson, "quad_exact": math.e - 1, "poly_diff": worst}

# ---------- Part 2 ----------
Ns = [10 * 2**j for j in range(10)]
uT = exact(T)
methods = {
    "Forward Euler": lambda N: one_step_run(fe_step, N),
    "RK2": lambda N: one_step_run(rk2_step, N),
    "RK4": lambda N: one_step_run(rk4_step, N),
    "Leapfrog": lambda N: leapfrog_run(N, "RK2"),
    "AB2": lambda N: ab_run(2, N, rk2_step),
    "AB4": lambda N: ab_run(4, N, rk4_step),
}
part2 = {"N": Ns, "exact": uT}
print(f"\nPart 2: logistic, u(2) = {uT:.15f}")
for name, run in methods.items():
    errs = [abs(run(N) - uT) for N in Ns]
    part2[name] = errs
    o = orders(errs)
    print(f"  {name:14s} " + " ".join(f"{e:.2e}" for e in errs))
    print(f"  {'order':14s} " + " ".join("   -    " if x is None else f"{x:8.2f}" for x in o))
out["part2"] = part2

# ---------- Part 3 ----------
part3 = {"N": Ns}
print("\nPart 3: starting values")
for label, run in {
    "Leapfrog, FE start": lambda N: leapfrog_run(N, "Forward Euler"),
    "Leapfrog, RK2 start": lambda N: leapfrog_run(N, "RK2"),
    "Leapfrog, exact start": lambda N: leapfrog_run(N, "exact"),
    "AB4, FE start": lambda N: ab_run(4, N, fe_step),
    "AB4, RK4 start": lambda N: ab_run(4, N, rk4_step),
}.items():
    errs = [abs(run(N) - uT) for N in Ns]
    part3[label] = errs
    o = orders(errs)
    print(f"  {label:22s} " + " ".join(f"{e:.2e}" for e in errs))
    print(f"  {'order':22s} " + " ".join("   -    " if x is None else f"{x:8.2f}" for x in o))
out["part3"] = part3

# ---------- Part 4 ----------
print("\nPart 4: evaluations of f, counted by wrapping f")
base_f = f


class Counter:
    def __init__(self):
        self.n = 0

    def __call__(self, u, t):
        self.n += 1
        return base_f(u, t)


def count_evals(kind, N):
    global f
    c = Counter()
    f = c
    try:
        val = one_step_run(rk4_step, N) if kind == "RK4" else ab_run(4, N, rk4_step)
    finally:
        f = base_f
    return c.n, abs(val - uT)


counts = []
for N in Ns:
    n4, e4 = count_evals("RK4", N)
    nab, eab = count_evals("AB4", N)
    counts.append({"N": N, "RK4 evals": n4, "RK4 error": e4, "AB4 evals": nab, "AB4 error": eab})
    print(f"  N={N:5d}  RK4 {n4:6d} evals {e4:.2e}   AB4 {nab:6d} evals {eab:.2e}")
out["part4"] = counts

# ---------- Part 5 ----------
print("\nPart 5: Gibbs phenomenon, square wave = (4/pi) sum_{k odd} sin(kx)/k")


def partial(x, terms):
    return 4 / math.pi * sum(math.sin((2 * j + 1) * x) / (2 * j + 1) for j in range(terms))


part5 = []
for terms in [4, 16, 64, 256]:
    # first maximum after the jump at x = 0 lies near pi / (2 * terms); search finely around it
    guess = math.pi / (2 * terms)
    best_x, best = guess, partial(guess, terms)
    lo, hi = 0.2 * guess, 3 * guess
    for _ in range(4):
        step = (hi - lo) / 400
        xs = [lo + i * step for i in range(401)]
        vals = [partial(x, terms) for x in xs]
        i = max(range(len(vals)), key=vals.__getitem__)
        best_x, best = xs[i], vals[i]
        lo, hi = best_x - step, best_x + step
    # error at a fixed point away from the jump
    away = abs(partial(math.pi / 2, terms) - 1)
    part5.append({"terms": terms, "max": best, "x_max": best_x, "pi_over_2terms": guess, "error_at_pi_2": away})
    print(f"  {terms:4d} terms: max {best:.6f} at x = {best_x:.6f} (pi/(2n) = {guess:.6f}), |S - 1| at pi/2 = {away:.2e}")
si_pi = sum((-1) ** j * math.pi ** (2 * j + 1) / ((2 * j + 1) * math.factorial(2 * j + 1)) for j in range(40))
print(f"  limit (2/pi) Si(pi) = {2 / math.pi * si_pi:.6f}")
out["part5"] = {"rows": part5, "limit": 2 / math.pi * si_pi}

# ---------- Part 6 ----------
print("\nPart 6: 1D Dirichlet matrix A = (1/h^2) tridiag(1, -2, 1)")


def thomas_ones(m, h):
    """Solve A w = 1 for A = (1/h^2) tridiag(1, -2, 1) of size m."""
    a = [1 / h**2] * m
    b = [-2 / h**2] * m
    c = [1 / h**2] * m
    d = [1.0] * m
    for i in range(1, m):
        w = a[i] / b[i - 1]
        b[i] -= w * c[i - 1]
        d[i] -= w * d[i - 1]
    x = [0.0] * m
    x[-1] = d[-1] / b[-1]
    for i in range(m - 2, -1, -1):
        x[i] = (d[i] - c[i] * x[i + 1]) / b[i]
    return x


part6 = []
for n in [8, 64, 512, 4096]:
    h = 1 / n
    m = n - 1
    lam = [2 / h**2 * (math.cos(p * math.pi * h) - 1) for p in range(1, m + 1)]
    small, large = min(abs(l) for l in lam), max(abs(l) for l in lam)
    w = thomas_ones(m, h)
    inf_norm = max(abs(v) for v in w)  # A^{-1} has one sign, so ||A^{-1}||_inf = max |A^{-1} 1|
    row = {"n": n, "inv2": 1 / small, "cond2": large / small, "invinf": inf_norm}
    part6.append(row)
    print(f"  h = 1/{n:4d}: ||A^-1||_2 = {1/small:.6f}  cond_2 = {large/small:.4e}  ||A^-1||_inf = {inf_norm:.6f}")
print(f"  1/pi^2 = {1/math.pi**2:.6f}, 4/(pi^2 h^2) at h=1/4096: {4*4096**2/math.pi**2:.4e}")
out["part6"] = part6

Path(__file__).with_name("experiment.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
