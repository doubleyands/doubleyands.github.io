"""Problem 7: 1D wave equation.

    u_tt = c^2 u_xx,  x in [0, 1],  t in [0, 1],  c = 2,  u(0, t) = u(1, t) = 0,
    u(x, 0) = sin(pi x) + 0.5 sin(4 pi x),  u_t(x, 0) = 0.
    Exact: u = sin(pi x) cos(2 pi t) + 0.5 sin(4 pi x) cos(8 pi t).
"""
import math
import torch
from pinn import grad, uniform, face, faces, grid

C = 2.0
P = math.pi
DOMAIN = [(0.0, 1.0), (0.0, 1.0)]  # x, t
OUT = 1


def exact(x):
    s, t = x[:, :1], x[:, 1:2]
    return torch.sin(P * s) * torch.cos(C * P * t) + 0.5 * torch.sin(4 * P * s) * torch.cos(4 * C * P * t)


def pde(net, x):
    du = grad(net(x)[:, 0], x)
    u_xx = grad(du[:, 0], x)[:, 0]
    u_tt = grad(du[:, 1], x)[:, 1]
    return u_tt - C ** 2 * u_xx


def ic(net, x):
    return net(x)[:, 0] - exact(x)[:, 0]


def ic_t(net, x):
    return grad(net(x)[:, 0], x)[:, 1]


def bc(net, x):
    return net(x)[:, 0]


TERMS = {"pde": pde, "ic": ic, "ic_t": ic_t, "bc": bc}


def points(g):
    return {"pde": uniform(g, 8192, DOMAIN), "ic": face(g, 256, DOMAIN, 1, 0),
            "ic_t": face(g, 256, DOMAIN, 1, 0), "bc": faces(g, 256, DOMAIN, [0])}


def reference():
    x, shape = grid(torch.linspace(0, 1, 101), torch.linspace(0, 1, 101))
    return x, shape, exact(x)
