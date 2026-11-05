"""Problem 4: steady advection-diffusion with a boundary layer.

    -eps u'' + u' = 1 on (0, 1),  u(0) = u(1) = 0,  eps = 0.01.
    Exact: u = x - (exp((x - 1) / eps) - exp(-1 / eps)) / (1 - exp(-1 / eps)).
    The layer at x = 1 has width about eps.
"""
import math
import torch
from pinn import grad, uniform, grid

EPS = 0.01
DOMAIN = [(0.0, 1.0)]
OUT = 1


def exact(x):
    x = x[:, :1]
    e = math.exp(-1 / EPS)
    return x - (torch.exp((x - 1) / EPS) - e) / (1 - e)


def pde(net, x):
    u = net(x)[:, 0]
    u_x = grad(u, x)[:, 0]
    u_xx = grad(u_x, x)[:, 0]
    return -EPS * u_xx + u_x - 1.0


def bc(net, x):
    return net(x)[:, 0]


TERMS = {"pde": pde, "bc": bc}


def points(g):
    return {"pde": uniform(g, 1024, DOMAIN), "bc": torch.tensor([[0.0], [1.0]])}


def reference():
    x, shape = grid(torch.linspace(0, 1, 1001))
    return x, shape, exact(x)
