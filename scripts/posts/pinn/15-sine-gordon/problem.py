"""Problem 15: sine-Gordon equation, a breather.

    u_tt - u_xx + sin(u) = 0,  x in [-10, 10],  t in [0, 8].
    Exact: u = 4 arctan( sqrt(1 - w^2) / w * sin(w t) / cosh(sqrt(1 - w^2) x) ),  w = 0.8
    (a bound state that oscillates in place with period 2 pi / w = 7.85).
    Initial values (u and u_t) and boundary values come from the exact solution.
"""
import math
import torch
from pinn import grad, uniform, face, faces, grid

W = 0.8
K = math.sqrt(1 - W ** 2)
DOMAIN = [(-10.0, 10.0), (0.0, 8.0)]  # x, t
OUT = 1


def exact(x):
    s, t = x[:, :1], x[:, 1:2]
    return 4 * torch.atan(K / W * torch.sin(W * t) / torch.cosh(K * s))


def u_t_exact(x):  # d/dt of the exact solution at t = 0: 4 K / cosh(K x)
    return 4 * K / torch.cosh(K * x[:, 0])


def pde(net, x):
    u = net(x)[:, 0]
    du = grad(u, x)
    u_xx = grad(du[:, 0], x)[:, 0]
    u_tt = grad(du[:, 1], x)[:, 1]
    return u_tt - u_xx + torch.sin(u)


def data(net, x):
    return net(x)[:, 0] - exact(x)[:, 0]


def ic_t(net, x):
    return grad(net(x)[:, 0], x)[:, 1] - u_t_exact(x)


TERMS = {"pde": pde, "ic": data, "ic_t": ic_t, "bc": data}


def points(g):
    return {"pde": uniform(g, 8192, DOMAIN), "ic": face(g, 512, DOMAIN, 1, 0),
            "ic_t": face(g, 512, DOMAIN, 1, 0), "bc": faces(g, 256, DOMAIN, [0])}


def reference():
    x, shape = grid(torch.linspace(-10, 10, 201), torch.linspace(0, 8, 161))
    return x, shape, exact(x)
