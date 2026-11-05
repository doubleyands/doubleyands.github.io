"""Problem 5: 1D heat equation.

    u_t = a u_xx,  x in [0, 1],  t in [0, 1],  a = 0.1,  u(0, t) = u(1, t) = 0,
    u(x, 0) = sin(pi x) + 0.5 sin(3 pi x).
    Exact: u = exp(-a pi^2 t) sin(pi x) + 0.5 exp(-9 a pi^2 t) sin(3 pi x).
"""
import math
import torch
from pinn import grad, uniform, face, faces, grid

A = 0.1
P = math.pi
DOMAIN = [(0.0, 1.0), (0.0, 1.0)]  # x, t
OUT = 1


def exact(x):
    s, t = x[:, :1], x[:, 1:2]
    return torch.exp(-A * P ** 2 * t) * torch.sin(P * s) + 0.5 * torch.exp(-9 * A * P ** 2 * t) * torch.sin(3 * P * s)


def pde(net, x):
    u = net(x)[:, 0]
    du = grad(u, x)
    u_xx = grad(du[:, 0], x)[:, 0]
    return du[:, 1] - A * u_xx


def ic(net, x):
    return net(x)[:, 0] - exact(x)[:, 0]


def bc(net, x):
    return net(x)[:, 0]


TERMS = {"pde": pde, "ic": ic, "bc": bc}


def points(g):
    return {"pde": uniform(g, 4096, DOMAIN), "ic": face(g, 256, DOMAIN, 1, 0), "bc": faces(g, 256, DOMAIN, [0])}


def reference():
    x, shape = grid(torch.linspace(0, 1, 101), torch.linspace(0, 1, 101))
    return x, shape, exact(x)
