"""Problem 9: 2D Poisson equation on the unit square.

    -(u_xx + u_yy) = f on (0, 1)^2,  u = 0 on the boundary.
    Exact: u = sin(pi x) sin(pi y) + 0.1 sin(4 pi x) sin(4 pi y), and f = -Laplacian(u).
"""
import math
import torch
from pinn import grad, uniform, faces, grid

P = math.pi
DOMAIN = [(0.0, 1.0), (0.0, 1.0)]
OUT = 1


def exact(x):
    s, y = x[:, :1], x[:, 1:2]
    return torch.sin(P * s) * torch.sin(P * y) + 0.1 * torch.sin(4 * P * s) * torch.sin(4 * P * y)


def f(x):
    s, y = x[:, 0], x[:, 1]
    return (2 * P ** 2 * torch.sin(P * s) * torch.sin(P * y)
            + 0.1 * 32 * P ** 2 * torch.sin(4 * P * s) * torch.sin(4 * P * y))


def pde(net, x):
    du = grad(net(x)[:, 0], x)
    lap = grad(du[:, 0], x)[:, 0] + grad(du[:, 1], x)[:, 1]
    return -lap - f(x)


def bc(net, x):
    return net(x)[:, 0]


TERMS = {"pde": pde, "bc": bc}


def points(g):
    return {"pde": uniform(g, 4096, DOMAIN), "bc": faces(g, 512, DOMAIN, [0, 1])}


def reference():
    x, shape = grid(torch.linspace(0, 1, 101), torch.linspace(0, 1, 101))
    return x, shape, exact(x)
