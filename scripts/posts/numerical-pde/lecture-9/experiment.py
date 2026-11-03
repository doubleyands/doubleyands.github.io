"""Experiment for the Lecture 9 post: stability regions, the root condition, and the heat equation.

Part 1: complex eigenvalues. The oscillator x' = y, y' = -x (eigenvalues +i and -i),
        x(0) = 1, y(0) = 0, so the exact amplitude sqrt(x^2 + y^2) stays 1.
        Forward Euler, Backward Euler, trapezoid and leapfrog with dt = 0.1 up to
        t = 20*pi (10 turns). Amplitude over time, and the theory
        |P(i dt)|^n for the one-step methods.
Part 2: the root condition for leapfrog on u' = i u (z = i dt), dt = 0.5, 1.0, 1.1.
        At dt = 1 the two roots meet at eta = i (a double root on the unit circle),
        so |U^n| grows linearly. Above 1, one root leaves the circle.
Part 3: leapfrog on u' = -u (real negative z). Zero-stable, so for a fixed final time
        it converges with order 2, but the second root eta ~ -(1 + dt) makes
        the error grow like e^t for long times.
Part 4: boundary locus z(theta) = rho(e^{i theta}) / sigma(e^{i theta}) for the four schemes,
        and the "pick a point and test" step with the largest root modulus.
Part 5: Backward Euler on u' = u with dt = 3 (z = 3, in the right half plane).
Part 6: the heat equation u_t = u_xx on [0, 1] with exact solution e^{-pi^2 t} sin(pi x).
        (a) truncation error of FTCS and Crank-Nicolson, by putting the exact solution
            into each scheme; FTCS against the formula (k/2 - h^2/12) u_xxxx, and the
            special ratio k = h^2/6 where the leading term vanishes.
        (b) Crank-Nicolson global error at t = 0.1 with k = h (tridiagonal solve each step).

Standard library only (cmath for complex numbers).
Run: python scripts/posts/numerical-pde/lecture-9/experiment.py   (writes experiment.json next to this file)
"""

import cmath
import json
import math
from pathlib import Path

out = {}


# ---------- part 1: oscillator ----------
def oscillator(method, dt, n_steps, every):
    x, y = 1.0, 0.0
    amps = [1.0]
    if method == "leapfrog":
        # second starting value from the exact solution, so only the scheme is tested
        xp, yp = x, y
        x, y = math.cos(dt), -math.sin(dt)
    for n in range(1, n_steps + 1):
        if method == "FE":
            x, y = x + dt * y, y - dt * x
        elif method == "BE":
            # (I - dt A) U^{n+1} = U^n with A = [[0, 1], [-1, 0]]
            d = 1 + dt * dt
            x, y = (x + dt * y) / d, (y - dt * x) / d
        elif method == "trap":
            a = dt / 2
            rx, ry = x + a * y, y - a * x
            d = 1 + a * a
            x, y = (rx + a * ry) / d, (ry - a * rx) / d
        elif method == "leapfrog":
            if n >= 2:
                xn, yn = xp + 2 * dt * y, yp - 2 * dt * x
                xp, yp = x, y
                x, y = xn, yn
        if n % every == 0:
            amps.append(math.hypot(x, y))
    return amps


dt1 = 0.1
turns = 10
N1 = round(turns * 2 * math.pi / dt1)  # 628 steps, t_end = 62.8
every = 4
p1 = {"dt": dt1, "N": N1, "every": every, "t": [n * dt1 for n in range(0, N1 + 1, every)]}
for m in ("FE", "BE", "trap", "leapfrog"):
    p1[m] = oscillator(m, dt1, N1, every)
p1["theory_end"] = {
    "FE": (1 + dt1**2) ** (N1 / 2),
    "BE": (1 + dt1**2) ** (-N1 / 2),
    "trap": 1.0,
}
p1["end"] = {m: p1[m][-1] for m in ("FE", "BE", "trap", "leapfrog")}
p1["leapfrog_range"] = [min(p1["leapfrog"]), max(p1["leapfrog"])]
p1["t_end"] = N1 * dt1
out["part1"] = p1


# ---------- part 2: leapfrog root condition on u' = i u ----------
def leapfrog_complex(lam, dt, N, U1=None):
    z = lam * dt
    Up, U = 1 + 0j, (cmath.exp(z) if U1 is None else U1)
    mags = [1.0, abs(U)]
    for _ in range(2, N + 1):
        Up, U = U, Up + 2 * z * U
        mags.append(abs(U))
    return mags


def roots_leapfrog(z):
    # eta^2 - 2 z eta - 1 = 0
    s = cmath.sqrt(z * z + 1)
    return z + s, z - s


p2 = {}
for dt in (0.5, 1.0, 1.1):
    z = 1j * dt
    r1, r2 = roots_leapfrog(z)
    rows = {}
    for N in (10, 100, 1000):
        mags = leapfrog_complex(1j, dt, N)
        rows[N] = {"final": mags[-1], "max": max(mags)}
    p2[str(dt)] = {"roots": [[r1.real, r1.imag], [r2.real, r2.imag]],
                   "root_mods": [abs(r1), abs(r2)], "runs": rows}
out["part2"] = p2


# ---------- part 3: leapfrog on u' = -u ----------
def leapfrog_real(dt, N):
    Up, U = 1.0, math.exp(-dt)
    errs = [0.0, 0.0]
    for n in range(2, N + 1):
        Up, U = U, Up - 2 * dt * U
        errs.append(abs(U - math.exp(-n * dt)))
    return errs


p3 = {"fixed_T": [], "long": []}
for N in (20, 40, 80, 160, 320, 640):
    dt = 2.0 / N
    e = leapfrog_real(dt, N)[-1]
    p3["fixed_T"].append({"dt_den": N // 2, "err": e})
for i in range(1, len(p3["fixed_T"])):
    a, b = p3["fixed_T"][i - 1]["err"], p3["fixed_T"][i]["err"]
    p3["fixed_T"][i]["order"] = math.log2(a / b)
for dt in (0.1, 0.01):
    for Tend in (2, 10, 20, 30):
        N = round(Tend / dt)
        e = leapfrog_real(dt, N)
        p3["long"].append({"dt": dt, "T": Tend, "err": e[-1], "exact": math.exp(-Tend),
                           "rel": e[-1] / math.exp(-Tend)})
    r2 = -dt - math.sqrt(dt * dt + 1)
    p3.setdefault("second_root", {})[str(dt)] = r2
out["part3"] = p3


# ---------- part 4: boundary locus ----------
SCHEMES = {
    # rho(eta), sigma(eta) as coefficient lists from eta^0 upward
    "FE": ([-1, 1], [1, 0]),
    "BE": ([-1, 1], [0, 1]),
    "trap": ([-1, 1], [0.5, 0.5]),
    "leapfrog": ([-1, 0, 1], [0, 2, 0]),
}


def poly(c, x):
    return sum(ck * x**k for k, ck in enumerate(c))


def max_root_mod(rho, sig, z):
    c = [r - z * s for r, s in zip(rho, sig)]
    while abs(c[-1]) < 1e-15:
        c.pop()
    if len(c) == 1:
        return math.inf
    if len(c) == 2:
        return abs(-c[0] / c[1])
    a, b, cc = c[2], c[1], c[0]
    s = cmath.sqrt(b * b - 4 * a * cc)
    return max(abs((-b + s) / (2 * a)), abs((-b - s) / (2 * a)))


p4 = {}
for name, (rho, sig) in SCHEMES.items():
    pts = []
    for j in range(12):
        th = 2 * math.pi * j / 12
        e = cmath.exp(1j * th)
        s = poly(sig, e)
        if abs(s) < 1e-12:
            pts.append(None)
            continue
        z = poly(rho, e) / s
        pts.append([round(z.real, 6), round(z.imag, 6)])
    tests = {lab: max_root_mod(rho, sig, zz) for lab, zz in
             (("-1", -1 + 0j), ("+3", 3 + 0j), ("0.5i", 0.5j), ("-0.1", -0.1 + 0j))}
    p4[name] = {"locus12": pts, "tests": tests}
out["part4"] = p4


# ---------- part 5: Backward Euler in the right half plane ----------
dt5 = 3.0
U, steps = 1.0, []
for n in range(4):
    U = U / (1 - dt5)
    steps.append({"n": n + 1, "BE": U, "exact": math.exp(dt5 * (n + 1))})
out["part5"] = {"dt": dt5, "P": 1 / (1 - dt5), "exact_factor": math.exp(dt5), "steps": steps}


# ---------- part 6: heat equation ----------
PI2 = math.pi**2


def uex(x, t):
    return math.exp(-PI2 * t) * math.sin(math.pi * x)


def lte(h_den, k, t0=0.05):
    h = 1.0 / h_den
    ft, cn = 0.0, 0.0
    ft_formula = 0.0
    for i in range(1, h_den):
        x = i * h
        u0 = [uex(x + d * h, t0) for d in (-1, 0, 1)]
        u1 = [uex(x + d * h, t0 + k) for d in (-1, 0, 1)]
        dtq = (u1[1] - u0[1]) / k
        d2_0 = (u0[0] - 2 * u0[1] + u0[2]) / h**2
        d2_1 = (u1[0] - 2 * u1[1] + u1[2]) / h**2
        ft = max(ft, abs(dtq - d2_0))
        cn = max(cn, abs(dtq - 0.5 * (d2_0 + d2_1)))
        uxxxx = PI2 * PI2 * u0[1]
        ft_formula = max(ft_formula, abs((k / 2 - h * h / 12) * uxxxx))
    return ft, cn, ft_formula


p6 = {"lte_k_eq_h": [], "lte_k_h2_6": [], "lte_k_h2_4": [], "cn_global": []}
for hd in (10, 20, 40, 80, 160):
    h = 1.0 / hd
    ft, cn, form = lte(hd, h)
    p6["lte_k_eq_h"].append({"h_den": hd, "FTCS": ft, "CN": cn, "FTCS_formula": form})
    ft, cn, form = lte(hd, h * h / 6)
    p6["lte_k_h2_6"].append({"h_den": hd, "FTCS": ft, "CN": cn, "FTCS_formula": form})
    ft, cn, form = lte(hd, h * h / 4)
    p6["lte_k_h2_4"].append({"h_den": hd, "FTCS": ft, "FTCS_formula": form})


def thomas(a, b, c, d):
    n = len(d)
    cp, dp = [0.0] * n, [0.0] * n
    cp[0], dp[0] = c[0] / b[0], d[0] / b[0]
    for i in range(1, n):
        m = b[i] - a[i] * cp[i - 1]
        cp[i] = c[i] / m if i < n - 1 else 0.0
        dp[i] = (d[i] - a[i] * dp[i - 1]) / m
    x = [0.0] * n
    x[-1] = dp[-1]
    for i in range(n - 2, -1, -1):
        x[i] = dp[i] - cp[i] * x[i + 1]
    return x


def crank_nicolson(hd, Tend=0.1):
    h = 1.0 / hd
    k = h
    N = round(Tend / k)
    assert abs(N * k - Tend) < 1e-12
    r = k / (2 * h * h)
    m = hd - 1
    U = [uex((i + 1) * h, 0.0) for i in range(m)]
    a, b, c = [-r] * m, [1 + 2 * r] * m, [-r] * m
    for _ in range(N):
        # boundary values are 0 here
        d = [r * (U[i - 1] if i > 0 else 0.0) + (1 - 2 * r) * U[i] + r * (U[i + 1] if i < m - 1 else 0.0)
             for i in range(m)]
        U = thomas(a, b, c, d)
    err = max(abs(U[i] - uex((i + 1) * h, Tend)) for i in range(m))
    return err, r, N


for hd in (10, 20, 40, 80, 160, 320):
    err, r, N = crank_nicolson(hd)
    p6["cn_global"].append({"h_den": hd, "err": err, "r": r, "N": N})
for key in ("lte_k_eq_h", "lte_k_h2_6", "lte_k_h2_4", "cn_global"):
    rows = p6[key]
    for i in range(1, len(rows)):
        for col in ("FTCS", "CN", "err"):
            if col in rows[i]:
                rows[i][col + "_order"] = math.log2(rows[i - 1][col] / rows[i][col])
out["part6"] = p6

Path(__file__).with_name("experiment.json").write_text(json.dumps(out, indent=1), encoding="utf-8")


def show(obj):
    print(json.dumps(obj, indent=1)[:6000])


if __name__ == "__main__":
    for key in ("part2", "part3", "part4", "part5", "part6"):
        print("=====", key)
        show(out[key])
    print("===== part1")
    show({k: v for k, v in out["part1"].items() if k in ("theory_end", "end", "leapfrog_range", "t_end", "N")})
