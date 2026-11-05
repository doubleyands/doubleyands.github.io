"""Problem 17: diffusion with a variable coefficient.

    -div(k grad u) = f on (0, 1)^2,  u = 0 on the boundary,
    k = 1 + 0.9 sin(2 pi x) sin(2 pi y)  (it varies from 0.1 to 1.9).
    Exact: u = sin(pi x) sin(pi y), so f = 2 pi^2 k u - grad k . grad u.
"""
import math
import torch
from pinn import grad, uniform, faces, grid

P = math.pi
DOMAIN = [(0.0, 1.0), (0.0, 1.0)]
OUT = 1


def k(x):
    return 1 + 0.9 * torch.sin(2 * P * x[:, 0]) * torch.sin(2 * P * x[:, 1])


def exact(x):
    return (torch.sin(P * x[:, 0]) * torch.sin(P * x[:, 1]))[:, None]


def f(x):
    s, y = x[:, 0], x[:, 1]
    u = torch.sin(P * s) * torch.sin(P * y)
    k_x = 0.9 * 2 * P * torch.cos(2 * P * s) * torch.sin(2 * P * y)
    k_y = 0.9 * 2 * P * torch.sin(2 * P * s) * torch.cos(2 * P * y)
    u_x = P * torch.cos(P * s) * torch.sin(P * y)
    u_y = P * torch.sin(P * s) * torch.cos(P * y)
    return 2 * P ** 2 * k(x) * u - (k_x * u_x + k_y * u_y)


def pde(net, x):
    du = grad(net(x)[:, 0], x)
    kx = k(x)
    div = grad(kx * du[:, 0], x)[:, 0] + grad(kx * du[:, 1], x)[:, 1]
    return -div - f(x)


def bc(net, x):
    return net(x)[:, 0]


TERMS = {"pde": pde, "bc": bc}


def points(g):
    return {"pde": uniform(g, 4096, DOMAIN), "bc": faces(g, 512, DOMAIN, [0, 1])}


def reference():
    x, shape = grid(torch.linspace(0, 1, 101), torch.linspace(0, 1, 101))
    return x, shape, exact(x)
