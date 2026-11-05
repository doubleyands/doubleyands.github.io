"""Problem 12: nonlinear Klein-Gordon equation.

    u_tt - u_xx + u^3 = f,  x in [0, 1],  t in [0, 1].
    Exact: u = x cos(5 pi t) + (x t)^3, and f follows from it. Initial and boundary values
    come from the exact solution.
"""
import math
import torch
from pinn import grad, uniform, face, faces, grid

P = math.pi
DOMAIN = [(0.0, 1.0), (0.0, 1.0)]  # x, t
OUT = 1


def exact(x):
    s, t = x[:, :1], x[:, 1:2]
    return s * torch.cos(5 * P * t) + (s * t) ** 3


def f(x):
    s, t = x[:, 0], x[:, 1]
    u = s * torch.cos(5 * P * t) + (s * t) ** 3
    return -25 * P ** 2 * s * torch.cos(5 * P * t) + 6 * s ** 3 * t - 6 * s * t ** 3 + u ** 3


def pde(net, x):
    u = net(x)[:, 0]
    du = grad(u, x)
    u_xx = grad(du[:, 0], x)[:, 0]
    u_tt = grad(du[:, 1], x)[:, 1]
    return u_tt - u_xx + u ** 3 - f(x)


def data(net, x):
    return net(x)[:, 0] - exact(x)[:, 0]


def ic_t(net, x):  # u_t(x, 0) = 0
    return grad(net(x)[:, 0], x)[:, 1]


TERMS = {"pde": pde, "ic": data, "ic_t": ic_t, "bc": data}


def points(g):
    return {"pde": uniform(g, 8192, DOMAIN), "ic": face(g, 256, DOMAIN, 1, 0),
            "ic_t": face(g, 256, DOMAIN, 1, 0), "bc": faces(g, 256, DOMAIN, [0])}


def reference():
    x, shape = grid(torch.linspace(0, 1, 101), torch.linspace(0, 1, 101))
    return x, shape, exact(x)
