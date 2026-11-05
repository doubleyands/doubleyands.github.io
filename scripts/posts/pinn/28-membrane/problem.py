"""Problem 28: 2D wave equation, a vibrating square membrane.

    u_tt = u_xx + u_yy on (0, 1)^2,  t in [0, 1],  u = 0 on the boundary,
    u(x, y, 0) = sin(pi x) sin(pi y) + 0.5 sin(2 pi x) sin(pi y),  u_t(x, y, 0) = 0.
    Exact: u = sin(pi x) sin(pi y) cos(sqrt(2) pi t) + 0.5 sin(2 pi x) sin(pi y) cos(sqrt(5) pi t).
"""
import math
import torch
from pinn import grad, uniform, face, faces, grid

P = math.pi
DOMAIN = [(0.0, 1.0), (0.0, 1.0), (0.0, 1.0)]  # x, y, t
OUT = 1


def exact(x):
    s, y, t = x[:, :1], x[:, 1:2], x[:, 2:3]
    return (torch.sin(P * s) * torch.sin(P * y) * torch.cos(math.sqrt(2) * P * t)
            + 0.5 * torch.sin(2 * P * s) * torch.sin(P * y) * torch.cos(math.sqrt(5) * P * t))


def pde(net, x):
    du = grad(net(x)[:, 0], x)
    lap = grad(du[:, 0], x)[:, 0] + grad(du[:, 1], x)[:, 1]
    return grad(du[:, 2], x)[:, 2] - lap


def ic(net, x):
    return net(x)[:, 0] - exact(x)[:, 0]


def ic_t(net, x):
    return grad(net(x)[:, 0], x)[:, 2]


def bc(net, x):
    return net(x)[:, 0]


TERMS = {"pde": pde, "ic": ic, "ic_t": ic_t, "bc": bc}


def points(g):
    return {"pde": uniform(g, 16384, DOMAIN), "ic": face(g, 1024, DOMAIN, 2, 0),
            "ic_t": face(g, 1024, DOMAIN, 2, 0), "bc": faces(g, 1024, DOMAIN, [0, 1])}


def reference():
    x, shape = grid(torch.linspace(0, 1, 41), torch.linspace(0, 1, 41), torch.linspace(0, 1, 11))
    return x, shape, exact(x)
