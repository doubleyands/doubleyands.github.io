"""Problem 10: Laplace equation on an L-shaped domain (a corner singularity).

    u_xx + u_yy = 0 on [-1, 1]^2 without the quadrant (0, 1] x [-1, 0).
    Exact: u = r^(2/3) sin(2 theta / 3), theta in [0, 2 pi). The gradient blows up at the
    corner (0, 0). u = 0 on the two edges that meet there; the other edges take the exact values.
"""
import math
import torch
from pinn import grad, uniform, face, grid

BOX = [(-1.0, 1.0), (-1.0, 1.0)]
DOMAIN = BOX
OUT = 1


def inside(x):
    return ~((x[:, 0] > 0) & (x[:, 1] < 0))


def exact(x):
    s, y = x[:, :1], x[:, 1:2]
    r = torch.sqrt(s ** 2 + y ** 2)
    theta = torch.remainder(torch.atan2(y, s), 2 * math.pi)
    return r ** (2 / 3) * torch.sin(2 * theta / 3)


def pde(net, x):
    du = grad(net(x)[:, 0], x)
    return grad(du[:, 0], x)[:, 0] + grad(du[:, 1], x)[:, 1]


def bc(net, x):
    return net(x)[:, 0] - exact(x)[:, 0]


TERMS = {"pde": pde, "bc": bc}


def points(g):
    x = uniform(g, 8192, BOX)
    outer = [face(g, 128, BOX, d, s) for d in (0, 1) for s in (0, 1)]
    inner = [face(g, 64, [(0.0, 1.0), (0.0, 0.0)], 1, 0),    # y = 0, x in [0, 1]
             face(g, 64, [(0.0, 0.0), (-1.0, 0.0)], 0, 0)]   # x = 0, y in [-1, 0]
    b = torch.cat(outer + inner)
    return {"pde": x[inside(x)][:6144], "bc": b[inside(b)]}


def reference():
    x, _ = grid(torch.linspace(-1, 1, 201), torch.linspace(-1, 1, 201))
    x = x[inside(x)]  # the points of the 201 x 201 grid that lie in the domain
    return x, (len(x),), exact(x)
