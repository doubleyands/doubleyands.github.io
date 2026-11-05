"""Problem 21: Stokes flow in a channel (Poiseuille flow).

    -grad p + mu Laplacian(u) = 0,  div u = 0  on (0, 2) x (0, 1),  mu = 0.1.
    No slip on the walls y = 0 and y = 1, a parabolic inlet at x = 0, p = 0 at the outlet x = 2.
    Exact: u = 4 y (1 - y), v = 0, p = 8 mu (2 - x).
    Outputs: u, v, p. The exact v is zero, so its error is reported as an RMS value.
"""
import torch
from pinn import grad, uniform, face, faces, grid

MU = 0.1
DOMAIN = [(0.0, 2.0), (0.0, 1.0)]
OUT = 3


def exact(x):
    s, y = x[:, :1], x[:, 1:2]
    return torch.cat([4 * y * (1 - y), torch.zeros_like(y), 8 * MU * (2 - s)], dim=1)


def fields(net, x):
    q = net(x)
    du, dv, dp = grad(q[:, 0], x), grad(q[:, 1], x), grad(q[:, 2], x)
    lap_u = grad(du[:, 0], x)[:, 0] + grad(du[:, 1], x)[:, 1]
    lap_v = grad(dv[:, 0], x)[:, 0] + grad(dv[:, 1], x)[:, 1]
    return du, dv, dp, lap_u, lap_v


def mom_x(net, x):
    du, dv, dp, lap_u, lap_v = fields(net, x)
    return -dp[:, 0] + MU * lap_u


def mom_y(net, x):
    du, dv, dp, lap_u, lap_v = fields(net, x)
    return -dp[:, 1] + MU * lap_v


def cont(net, x):
    q = net(x)
    return grad(q[:, 0], x)[:, 0] + grad(q[:, 1], x)[:, 1]


def velocity(net, x):  # walls and inlet: u and v from the exact solution
    return net(x)[:, :2] - exact(x)[:, :2]


def outlet(net, x):
    return net(x)[:, 2]


TERMS = {"mom_x": mom_x, "mom_y": mom_y, "cont": cont, "wall": velocity, "inlet": velocity, "outlet": outlet}


def points(g):
    x = uniform(g, 4096, DOMAIN)
    return {"mom_x": x, "mom_y": x.clone(), "cont": x.clone(), "wall": faces(g, 512, DOMAIN, [1]),
            "inlet": face(g, 256, DOMAIN, 0, 0), "outlet": face(g, 256, DOMAIN, 0, 1)}


def reference():
    x, shape = grid(torch.linspace(0, 2, 101), torch.linspace(0, 1, 51))
    return x, shape, exact(x)
