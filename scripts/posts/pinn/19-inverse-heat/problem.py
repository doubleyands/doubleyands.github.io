"""Problem 19: inverse problem, find the diffusivity of the 1D heat equation from data.

    u_t = kappa u_xx,  x in [0, 1],  t in [0, 1],  kappa unknown (true value 0.1).
    Data: u at 512 random points inside the domain, taken from the exact solution of Problem 5
    (no noise). No initial or boundary conditions are given.
    kappa = exp(log_kappa), log_kappa starts at 0 (kappa = 1, ten times too large).
"""
import math
import torch
from pinn import grad, uniform, grid

KAPPA = 0.1
P = math.pi
DOMAIN = [(0.0, 1.0), (0.0, 1.0)]  # x, t
OUT = 1
PARAMS = {"log_kappa": (0.0, math.log(KAPPA))}


def exact(x):
    s, t = x[:, :1], x[:, 1:2]
    return (torch.exp(-KAPPA * P ** 2 * t) * torch.sin(P * s)
            + 0.5 * torch.exp(-9 * KAPPA * P ** 2 * t) * torch.sin(3 * P * s))


def pde(net, x):
    du = grad(net(x)[:, 0], x)
    u_xx = grad(du[:, 0], x)[:, 0]
    return du[:, 1] - torch.exp(net.extra["log_kappa"]) * u_xx


def data(net, x):
    return net(x)[:, 0] - exact(x)[:, 0]


TERMS = {"pde": pde, "data": data}


def points(g):
    return {"pde": uniform(g, 4096, DOMAIN), "data": uniform(g, 512, DOMAIN)}


def reference():
    x, shape = grid(torch.linspace(0, 1, 101), torch.linspace(0, 1, 101))
    return x, shape, exact(x)
