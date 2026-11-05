"""Problem 13: Korteweg-de Vries equation, a single soliton.

    u_t + 6 u u_x + u_xxx = 0,  x in [-10, 10],  t in [0, 2].
    Exact: u = c/2 sech^2(sqrt(c)/2 (x - c t - x0)),  c = 4,  x0 = -4 (it moves from x = -4 to 4).
    Initial and boundary values come from the exact solution.
"""
import math
import torch
from pinn import grad, uniform, face, faces, grid

C, X0 = 4.0, -4.0
DOMAIN = [(-10.0, 10.0), (0.0, 2.0)]  # x, t
OUT = 1


def exact(x):
    s, t = x[:, :1], x[:, 1:2]
    return C / 2 / torch.cosh(math.sqrt(C) / 2 * (s - C * t - X0)) ** 2


def pde(net, x):
    u = net(x)[:, 0]
    du = grad(u, x)
    u_xx = grad(du[:, 0], x)[:, 0]
    u_xxx = grad(u_xx, x)[:, 0]
    return du[:, 1] + 6 * u * du[:, 0] + u_xxx


def data(net, x):
    return net(x)[:, 0] - exact(x)[:, 0]


TERMS = {"pde": pde, "ic": data, "bc": data}


def points(g):
    return {"pde": uniform(g, 8192, DOMAIN), "ic": face(g, 512, DOMAIN, 1, 0), "bc": faces(g, 256, DOMAIN, [0])}


def reference():
    x, shape = grid(torch.linspace(-10, 10, 201), torch.linspace(0, 2, 101))
    return x, shape, exact(x)
