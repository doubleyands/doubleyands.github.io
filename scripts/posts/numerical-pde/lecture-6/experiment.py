"""Experiment for the Lecture 6 post: the first time-stepping schemes for u' = f(u, t).

Test problem: u' = 1 + u^2, u(0) = 0, exact solution u = tan t (blows up at pi/2).
We integrate to T = 1.

Part 1: the existence/uniqueness bound for this problem against the real blow-up time.
Part 2: global error at T = 1 for Forward Euler, Backward Euler, trapezoidal, and
        the 2nd-order Taylor method, with k = 1/10 ... 1/10240 (halving).
Part 3: one-step (local) error of Forward Euler from the exact value u(0.5),
        against the global error. Local ~ k^2, global ~ k.
Part 4: Newton iterations per step for Backward Euler and trapezoidal.
Part 5: a preview of stability. u' = -50 (u - cos t) - sin t, u(0) = 1, exact u = cos t,
        to T = 2. Forward and Backward Euler with k = 0.01, 0.025, 0.05.

Standard library only.
Run: python scripts/posts/numerical-pde/lecture-6/experiment.py   (writes experiment.json next to this file)
"""

import json
import math
from pathlib import Path

TOL = 1e-13


def f(u, t):
    return 1 + u * u


def fu(u, t):
    return 2 * u


def ft(u, t):
    return 0.0


def newton(g, dg, guess):
    """Solve g(U) = 0 by Newton's method. Returns (root, number of iterations)."""
    U = guess
    for it in range(1, 51):
        step = g(U) / dg(U)
        U -= step
        if abs(step) <= TOL * max(1.0, abs(U)):
            return U, it
    raise RuntimeError("Newton did not converge")


def forward_euler(U, t, k, f=f, fu=fu):
    return U + k * f(U, t), 0


def backward_euler(U, t, k, f=f, fu=fu):
    t1 = t + k
    return newton(lambda V: V - k * f(V, t1) - U, lambda V: 1 - k * fu(V, t1), U)


def trapezoidal(U, t, k, f=f, fu=fu):
    t1 = t + k
    known = U + 0.5 * k * f(U, t)
    return newton(lambda V: V - 0.5 * k * f(V, t1) - known, lambda V: 1 - 0.5 * k * fu(V, t1), U)


def taylor2(U, t, k, f=f, fu=fu):
    F = f(U, t)
    return U + k * F + 0.5 * k * k * (fu(U, t) * F + ft(U, t)), 0


SCHEMES = {
    "Forward Euler": forward_euler,
    "Backward Euler": backward_euler,
    "Trapezoidal": trapezoidal,
    "Taylor 2": taylor2,
}


def run(step, k, N, U0=0.0, **kw):
    U, t, iters = U0, 0.0, []
    for n in range(N):
        U, it = step(U, t, k, **kw)
        t = (n + 1) * k
        iters.append(it)
        if not math.isfinite(U) or abs(U) > 1e300:
            return math.inf, iters
    return U, iters


def order(e1, e2):
    return math.log2(e1 / e2)


out = {}

# Part 1: guaranteed time a / S with S = max |f| = 1 + a^2 on |u| <= a, best a.
best = max(((a / (1 + a * a), a) for a in (i / 1000 for i in range(1, 5001))))
out["part1"] = {"best_a": best[1], "guaranteed_T": best[0], "blowup_T": math.pi / 2}
print(f"Part 1: best a = {best[1]:.3f}, guaranteed T = {best[0]:.4f}, blow-up at pi/2 = {math.pi / 2:.4f}")

# Part 2: global error at T = 1
T = 1.0
exact = math.tan(T)
Ns = [10 * 2**j for j in range(11)]
part2 = {"N": Ns, "exact": exact}
print(f"\nPart 2: error at T = 1 (u(1) = tan 1 = {exact:.15f})")
for name, step in SCHEMES.items():
    errs = []
    for N in Ns:
        U, _ = run(step, T / N, N)
        errs.append(abs(U - exact))
    part2[name] = errs
    orders = ["  -  "] + [f"{order(errs[i - 1], errs[i]):5.2f}" for i in range(1, len(errs))]
    print(f"  {name:15s}", " ".join(f"{e:.2e}" for e in errs))
    print(f"  {'order':15s}", " ".join(f"{o:>8s}" for o in orders))
# signed errors for k = 1/160 to show the sign
signed = {}
for name, step in SCHEMES.items():
    U, _ = run(step, T / 160, 160)
    signed[name] = U - exact
part2["signed_N160"] = signed
print("  signed error U - u at N = 160:", {n: f"{v:+.3e}" for n, v in signed.items()})
out["part2"] = part2

# Part 3: Forward Euler, one step from the exact value at t0 = 0.5
t0 = 0.5
u0 = math.tan(t0)
u2 = (1 + u0 * u0) * 2 * u0  # u'' = f_u f
part3 = []
print(f"\nPart 3: Forward Euler, one step from u({t0}) exact; u''({t0}) = {u2:.6f}")
for N in Ns:
    k = T / N
    one = u0 + k * f(u0, t0) - math.tan(t0 + k)
    part3.append({"k": k, "local": abs(one), "ratio": abs(one) / (k * k / 2 * u2), "global": part2["Forward Euler"][Ns.index(N)]})
for i, r in enumerate(part3):
    lo = "" if i == 0 else f" order {order(part3[i - 1]['local'], r['local']):.2f}"
    print(f"  k = 1/{Ns[i]:5d}: local {r['local']:.3e} (/(k^2/2 u'') = {r['ratio']:.4f}){lo}, global {r['global']:.3e}")
out["part3"] = part3

# Part 4: Newton iterations
part4 = {}
print("\nPart 4: Newton iterations per step (min / mean / max)")
for name in ("Backward Euler", "Trapezoidal"):
    rows = []
    for N in (10, 160, 2560):
        _, it = run(SCHEMES[name], T / N, N)
        rows.append({"N": N, "min": min(it), "mean": sum(it) / len(it), "max": max(it)})
        print(f"  {name:15s} k = 1/{N:5d}: {min(it)} / {sum(it) / len(it):.2f} / {max(it)}")
    part4[name] = rows
out["part4"] = part4

# Part 5: stiff linear problem
LAM = -50.0


def g5(u, t):
    return LAM * (u - math.cos(t)) - math.sin(t)


def g5u(u, t):
    return LAM


part5 = []
print(f"\nPart 5: u' = {LAM:g}(u - cos t) - sin t, u(0) = 1, error at T = 2")
for k in (0.01, 0.025, 0.05):
    N = round(2 / k)
    assert abs(N * k - 2) < 1e-12
    row = {"k": k, "lam_k": LAM * k}
    for name in ("Forward Euler", "Backward Euler"):
        U, _ = run(SCHEMES[name], k, N, U0=1.0, f=g5, fu=g5u)
        row[name] = abs(U - math.cos(2.0))
    part5.append(row)
    print(f"  k = {k}: lambda k = {LAM * k:g}, FE error {row['Forward Euler']:.3e}, BE error {row['Backward Euler']:.3e}")
out["part5"] = part5

Path(__file__).with_name("experiment.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
