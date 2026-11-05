"""Problem 8: Fisher-KPP reaction-diffusion, a travelling front.

    u_t = u_xx + u (1 - u),  x in [-10, 20],  t in [0, 5].
    Exact (Ablowitz-Zeppetella): u = (1 + exp(x / sqrt(6) - 5 t / 6))^(-2), speed 5 / sqrt(6).
    The initial and boundary values come from the exact solution.
"""
import math
import torch
from pinn import grad, uniform, face, faces, grid

DOMAIN = [(-10.0, 20.0), (0.0, 5.0)]  # x, t
OUT = 1


def exact(x):
    s, t = x[:, :1], x[:, 1:2]
    return (1 + torch.exp(s / math.sqrt(6) - 5 * t / 6)) ** -2


def pde(net, x):
    u = net(x)[:, 0]
    du = grad(u, x)
    u_xx = grad(du[:, 0], x)[:, 0]
    return du[:, 1] - u_xx - u * (1 - u)


def data(net, x):  # initial and boundary values
    return net(x)[:, 0] - exact(x)[:, 0]


TERMS = {"pde": pde, "ic": data, "bc": data}


def points(g):
    return {"pde": uniform(g, 8192, DOMAIN), "ic": face(g, 512, DOMAIN, 1, 0), "bc": faces(g, 256, DOMAIN, [0])}


def reference():
    x, shape = grid(torch.linspace(-10, 20, 301), torch.linspace(0, 5, 101))
    return x, shape, exact(x)
