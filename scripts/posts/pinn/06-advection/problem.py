"""Problem 6: 1D advection with a high speed (a known hard case for PINNs).

    u_t + b u_x = 0,  x in [0, 2 pi],  t in [0, 1],  b = 30,  u(x, 0) = sin(x),  periodic in x.
    Exact: u = sin(x - b t).
"""
import math
import torch
from pinn import grad, uniform, face, grid

B = 30.0
L = 2 * math.pi
DOMAIN = [(0.0, L), (0.0, 1.0)]  # x, t
OUT = 1
SHIFT = torch.tensor([L, 0.0])


def exact(x):
    return torch.sin(x[:, :1] - B * x[:, 1:2])


def pde(net, x):
    du = grad(net(x)[:, 0], x)
    return du[:, 1] + B * du[:, 0]


def ic(net, x):
    return net(x)[:, 0] - torch.sin(x[:, 0])


def periodic(net, x):  # x on the face x = 0; compare with the same t on x = 2 pi
    return net(x)[:, 0] - net(x + SHIFT.to(x.device))[:, 0]


TERMS = {"pde": pde, "ic": ic, "periodic": periodic}


def points(g):
    return {"pde": uniform(g, 8192, DOMAIN), "ic": face(g, 256, DOMAIN, 1, 0), "periodic": face(g, 256, DOMAIN, 0, 0)}


def reference():
    x, shape = grid(torch.linspace(0, L, 201), torch.linspace(0, 1, 101))
    return x, shape, exact(x)
