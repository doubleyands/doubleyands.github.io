"""Problem 1: exponential decay.

    u'(t) = -k u,  u(0) = 1,  t in [0, 3],  k = 2.    Exact: u = exp(-k t).
"""
import torch
from pinn import grad, uniform, grid

K = 2.0
DOMAIN = [(0.0, 3.0)]
OUT = 1


def exact(x):
    return torch.exp(-K * x[:, :1])


def pde(net, x):
    u = net(x)[:, 0]
    return grad(u, x)[:, 0] + K * u


def ic(net, x):
    return net(x)[:, 0] - 1.0


TERMS = {"pde": pde, "ic": ic}


def points(g):
    return {"pde": uniform(g, 512, DOMAIN), "ic": torch.zeros(1, 1)}


def reference():
    x, shape = grid(torch.linspace(0, 3, 301))
    return x, shape, exact(x)
