"""Problem 2: 1D Poisson equation with two scales.

    -u'' = f on (0, 1),  u(0) = u(1) = 0.
    Exact: u = sin(2 pi x) + 0.1 sin(10 pi x), and f = -u''.
"""
import math
import torch
from pinn import grad, uniform, grid

DOMAIN = [(0.0, 1.0)]
OUT = 1
P = math.pi


def exact(x):
    x = x[:, :1]
    return torch.sin(2 * P * x) + 0.1 * torch.sin(10 * P * x)


def f(x):
    x = x[:, 0]
    return (2 * P) ** 2 * torch.sin(2 * P * x) + 0.1 * (10 * P) ** 2 * torch.sin(10 * P * x)


def pde(net, x):
    u = net(x)[:, 0]
    u_xx = grad(grad(u, x)[:, 0], x)[:, 0]
    return -u_xx - f(x)


def bc(net, x):
    return net(x)[:, 0]


TERMS = {"pde": pde, "bc": bc}


def points(g):
    return {"pde": uniform(g, 1024, DOMAIN), "bc": torch.tensor([[0.0], [1.0]])}


def reference():
    x, shape = grid(torch.linspace(0, 1, 501))
    return x, shape, exact(x)
