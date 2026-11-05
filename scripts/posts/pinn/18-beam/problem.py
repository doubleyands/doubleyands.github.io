"""Problem 18: Euler-Bernoulli cantilever beam under a uniform load (a 4th-order ODE).

    u'''' = 1 on (0, 1)  (EI = 1, load q = 1),
    clamped at x = 0: u = u' = 0;  free at x = 1: u'' = u''' = 0 (no moment, no shear).
    Exact: u = x^2 (6 - 4 x + x^2) / 24, tip deflection 1/8.
"""
import torch
from pinn import grad, uniform, grid

DOMAIN = [(0.0, 1.0)]
OUT = 1


def exact(x):
    s = x[:, :1]
    return s ** 2 * (6 - 4 * s + s ** 2) / 24


def derivs(net, x, n):
    out = [net(x)[:, 0]]
    for _ in range(n):
        out.append(grad(out[-1], x)[:, 0])
    return out


def pde(net, x):
    return derivs(net, x, 4)[4] - 1.0


def clamp(net, x):
    u, u_x = derivs(net, x, 1)
    return torch.stack([u, u_x], dim=1)


def free(net, x):
    _, _, u_xx, u_xxx = derivs(net, x, 3)
    return torch.stack([u_xx, u_xxx], dim=1)


TERMS = {"pde": pde, "clamp": clamp, "free": free}


def points(g):
    return {"pde": uniform(g, 1024, DOMAIN), "clamp": torch.zeros(1, 1), "free": torch.ones(1, 1)}


def reference():
    x, shape = grid(torch.linspace(0, 1, 501))
    return x, shape, exact(x)
