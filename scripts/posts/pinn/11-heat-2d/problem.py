"""Problem 11: 2D heat equation.

    u_t = a (u_xx + u_yy) on (0, 1)^2,  t in [0, 1],  a = 0.1,  u = 0 on the boundary.
    Exact: u = exp(-2 a pi^2 t) sin(pi x) sin(pi y) + 0.5 exp(-8 a pi^2 t) sin(2 pi x) sin(2 pi y).
"""
import math
import torch
from pinn import grad, uniform, face, faces, grid

A = 0.1
P = math.pi
DOMAIN = [(0.0, 1.0), (0.0, 1.0), (0.0, 1.0)]  # x, y, t
OUT = 1


def exact(x):
    s, y, t = x[:, :1], x[:, 1:2], x[:, 2:3]
    return (torch.exp(-2 * A * P ** 2 * t) * torch.sin(P * s) * torch.sin(P * y)
            + 0.5 * torch.exp(-8 * A * P ** 2 * t) * torch.sin(2 * P * s) * torch.sin(2 * P * y))


def pde(net, x):
    du = grad(net(x)[:, 0], x)
    lap = grad(du[:, 0], x)[:, 0] + grad(du[:, 1], x)[:, 1]
    return du[:, 2] - A * lap


def ic(net, x):
    return net(x)[:, 0] - exact(x)[:, 0]


def bc(net, x):
    return net(x)[:, 0]


TERMS = {"pde": pde, "ic": ic, "bc": bc}


def points(g):
    return {"pde": uniform(g, 8192, DOMAIN), "ic": face(g, 1024, DOMAIN, 2, 0),
            "bc": faces(g, 1024, DOMAIN, [0, 1])}


def reference():
    x, shape = grid(torch.linspace(0, 1, 41), torch.linspace(0, 1, 41), torch.linspace(0, 1, 11))
    return x, shape, exact(x)
