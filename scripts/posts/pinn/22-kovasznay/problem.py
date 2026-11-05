"""Problem 22: Kovasznay flow, an exact solution of the steady Navier-Stokes equations.

    (u . grad) u + grad p - (1/Re) Laplacian(u) = 0,  div u = 0  on (-0.5, 1) x (-0.5, 1.5),  Re = 20.
    Exact: u = 1 - exp(l x) cos(2 pi y),  v = l / (2 pi) exp(l x) sin(2 pi y),  p = (1 - exp(2 l x)) / 2,
    with l = Re/2 - sqrt(Re^2/4 + 4 pi^2). Velocity on the boundary from the exact solution.
    Outputs: u, v, p. p is defined up to a constant, so its error is measured after removing the mean.
"""
import math
import torch
from pinn import grad, uniform, faces, grid

RE = 20.0
LAM = RE / 2 - math.sqrt(RE ** 2 / 4 + 4 * math.pi ** 2)
P2 = 2 * math.pi
DOMAIN = [(-0.5, 1.0), (-0.5, 1.5)]
OUT = 3
MEAN_FREE = [2]


def exact(x):
    s, y = x[:, :1], x[:, 1:2]
    e = torch.exp(LAM * s)
    return torch.cat([1 - e * torch.cos(P2 * y), LAM / P2 * e * torch.sin(P2 * y), (1 - e ** 2) / 2], dim=1)


# ponytail: mom_x and mom_y each recompute these derivatives (about 2x work); share them if NS time matters
def fields(net, x):
    q = net(x)
    u, v = q[:, 0], q[:, 1]
    du, dv, dp = grad(u, x), grad(v, x), grad(q[:, 2], x)
    lap_u = grad(du[:, 0], x)[:, 0] + grad(du[:, 1], x)[:, 1]
    lap_v = grad(dv[:, 0], x)[:, 0] + grad(dv[:, 1], x)[:, 1]
    return u, v, du, dv, dp, lap_u, lap_v


def mom_x(net, x):
    u, v, du, dv, dp, lap_u, lap_v = fields(net, x)
    return u * du[:, 0] + v * du[:, 1] + dp[:, 0] - lap_u / RE


def mom_y(net, x):
    u, v, du, dv, dp, lap_u, lap_v = fields(net, x)
    return u * dv[:, 0] + v * dv[:, 1] + dp[:, 1] - lap_v / RE


def cont(net, x):
    q = net(x)
    return grad(q[:, 0], x)[:, 0] + grad(q[:, 1], x)[:, 1]


def bc(net, x):
    return net(x)[:, :2] - exact(x)[:, :2]


TERMS = {"mom_x": mom_x, "mom_y": mom_y, "cont": cont, "bc": bc}


def points(g):
    x = uniform(g, 8192, DOMAIN)
    return {"mom_x": x, "mom_y": x.clone(), "cont": x.clone(), "bc": faces(g, 1024, DOMAIN, [0, 1])}


def reference():
    x, shape = grid(torch.linspace(-0.5, 1, 76), torch.linspace(-0.5, 1.5, 101))
    return x, shape, exact(x)
