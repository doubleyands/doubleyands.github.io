"""Problem 27: Poisson equation on an annulus (a curved domain with a hole).

    -(u_xx + u_yy) = f on 0.4 < r < 1,  u given on both circles.
    Exact: u = cos(pi x) sin(pi y) + (x^2 + y^2) / 4, so f = 2 pi^2 cos(pi x) sin(pi y) - 1.
"""
import math
import torch
from pinn import grad, uniform, grid

R1, R2 = 0.4, 1.0
P = math.pi
BOX = [(-1.0, 1.0), (-1.0, 1.0)]
DOMAIN = BOX
OUT = 1


def inside(x):
    r2 = (x ** 2).sum(dim=1)
    return (r2 >= R1 ** 2) & (r2 <= R2 ** 2)


def exact(x):
    s, y = x[:, :1], x[:, 1:2]
    return torch.cos(P * s) * torch.sin(P * y) + (s ** 2 + y ** 2) / 4


def f(x):
    return 2 * P ** 2 * torch.cos(P * x[:, 0]) * torch.sin(P * x[:, 1]) - 1


def pde(net, x):
    du = grad(net(x)[:, 0], x)
    return -(grad(du[:, 0], x)[:, 0] + grad(du[:, 1], x)[:, 1]) - f(x)


def bc(net, x):
    return net(x)[:, 0] - exact(x)[:, 0]


TERMS = {"pde": pde, "bc": bc}


def circle(g, n, r):
    theta = 2 * math.pi * torch.rand(n, generator=g)
    return torch.stack([r * torch.cos(theta), r * torch.sin(theta)], dim=1)


def points(g):
    x = uniform(g, 8192, BOX)
    return {"pde": x[inside(x)][:4096], "bc": torch.cat([circle(g, 384, R2), circle(g, 160, R1)])}


def reference():
    x, _ = grid(torch.linspace(-1, 1, 201), torch.linspace(-1, 1, 201))
    x = x[inside(x)]  # the points of the 201 x 201 grid that lie in the annulus
    return x, (len(x),), exact(x)
