"""Problem 30: Allen-Cahn equation (sharp moving interfaces).

    u_t - eps u_xx + 5 u^3 - 5 u = 0,  x in [-1, 1] periodic,  t in [0, 1],  eps = 1e-4,
    u(x, 0) = x^2 cos(pi x).
    No closed form: the reference is a Fourier spectral solution with ETDRK4 time stepping
    (Kassam and Trefethen, SIAM J. Sci. Comput. 26, 2005), 4096 modes, step 5e-5.
    u0 is periodic but u0' is not (-2 and 2 at x = -1 and 1), so the modes converge only
    algebraically. Against 8192 modes: relative L2 difference 2e-7, max 1.5e-5 (near the kink,
    t = 0.01). Halving the step changes u by 3e-11.
"""
import math
import torch
from pinn import grad, uniform, face, grid

EPS = 1e-4
DOMAIN = [(-1.0, 1.0), (0.0, 1.0)]  # x, t
OUT = 1
TRAIN = dict(adam=40_000)  # the hardest problem gets twice the Adam steps
SHIFT = torch.tensor([2.0, 0.0])


def pde(net, x):
    u = net(x)[:, 0]
    du = grad(u, x)
    u_xx = grad(du[:, 0], x)[:, 0]
    return du[:, 1] - EPS * u_xx + 5 * u ** 3 - 5 * u


def ic(net, x):
    s = x[:, 0]
    return net(x)[:, 0] - s ** 2 * torch.cos(math.pi * s)


def periodic(net, x):  # x on the face x = -1; u and u_x must match x = 1 at the same t
    xr = x + SHIFT.to(x.device)
    u_l, u_r = net(x)[:, 0], net(xr)[:, 0]
    return torch.stack([u_l - u_r, grad(u_l, x)[:, 0] - grad(u_r, xr)[:, 0]], dim=1)


TERMS = {"pde": pde, "ic": ic, "periodic": periodic}


def points(g):
    return {"pde": uniform(g, 16384, DOMAIN), "ic": face(g, 512, DOMAIN, 1, 0), "periodic": face(g, 256, DOMAIN, 0, 0)}


def spectral(n, h, save_every):
    """u on the grid x_j = -1 + 2j/n at t = 0, save_every*h, ... , 1. Returns (n_t, n)."""
    x = -1 + 2 * torch.arange(n) / n
    v = torch.fft.fft(x ** 2 * torch.cos(math.pi * x))
    k = math.pi * torch.fft.fftfreq(n, d=1 / n)  # wavenumbers for period 2
    L = -EPS * k ** 2
    E, E2 = torch.exp(h * L), torch.exp(h * L / 2)
    M = 32  # contour points for the phi functions
    r = torch.exp(1j * math.pi * (torch.arange(1, M + 1) - 0.5) / M)
    LR = h * L[:, None] + r[None, :]
    Q = h * ((torch.exp(LR / 2) - 1) / LR).mean(1).real
    f1 = h * ((-4 - LR + torch.exp(LR) * (4 - 3 * LR + LR ** 2)) / LR ** 3).mean(1).real
    f2 = h * ((2 + LR + torch.exp(LR) * (-2 + LR)) / LR ** 3).mean(1).real
    f3 = h * ((-4 - 3 * LR - LR ** 2 + torch.exp(LR) * (4 - LR)) / LR ** 3).mean(1).real

    def N(v):
        u = torch.fft.ifft(v).real
        return torch.fft.fft(5 * u - 5 * u ** 3)

    out = [torch.fft.ifft(v).real]
    for step in range(1, round(1 / h) + 1):
        Nv = N(v)
        a = E2 * v + Q * Nv
        Na = N(a)
        b = E2 * v + Q * Na
        Nb = N(b)
        c = E2 * a + Q * (2 * Nb - Nv)
        Nc = N(c)
        v = E * v + Nv * f1 + 2 * (Na + Nb) * f2 + Nc * f3
        if step % save_every == 0:
            out.append(torch.fft.ifft(v).real)
    return torch.stack(out)


def reference():
    sol = spectral(4096, 5e-5, 200)[:, ::16]  # 256 points in x, 101 times (every 0.01)
    x, shape = grid(-1 + 2 * torch.arange(256) / 256, torch.linspace(0, 1, 101))
    return x, shape, sol.T.reshape(-1, 1)


def check():
    coarse = spectral(4096, 5e-5, 200)[:, ::16]
    fine = spectral(8192, 5e-5, 200)[:, ::32]
    diff = float((coarse - fine).norm() / fine.norm())
    assert diff < 1e-6, diff
