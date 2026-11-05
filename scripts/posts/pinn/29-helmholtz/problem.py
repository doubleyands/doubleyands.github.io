"""Problem 29: 2D Helmholtz equation.

    u_xx + u_yy + k^2 u = q on (-1, 1)^2,  u = 0 on the boundary,  k = 1.
    Exact: u = sin(a1 pi x) sin(a2 pi y),  a1 = 1, a2 = 4, so q = (k^2 - (a1^2 + a2^2) pi^2) u.
"""
import math
import torch
from pinn import grad, uniform, faces, grid

A1, A2, KW = 1.0, 4.0, 1.0
P = math.pi
DOMAIN = [(-1.0, 1.0), (-1.0, 1.0)]
OUT = 1


def exact(x):
    return (torch.sin(A1 * P * x[:, 0]) * torch.sin(A2 * P * x[:, 1]))[:, None]


def pde(net, x):
    u = net(x)[:, 0]
    du = grad(u, x)
    lap = grad(du[:, 0], x)[:, 0] + grad(du[:, 1], x)[:, 1]
    q = (KW ** 2 - (A1 ** 2 + A2 ** 2) * P ** 2) * exact(x)[:, 0]
    return lap + KW ** 2 * u - q


def bc(net, x):
    return net(x)[:, 0]


TERMS = {"pde": pde, "bc": bc}


def points(g):
    return {"pde": uniform(g, 8192, DOMAIN), "bc": faces(g, 1024, DOMAIN, [0, 1])}


def reference():
    x, shape = grid(torch.linspace(-1, 1, 101), torch.linspace(-1, 1, 101))
    return x, shape, exact(x)
