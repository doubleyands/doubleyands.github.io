"""Problem 23: Taylor-Green vortex, an exact solution of the unsteady 2D Navier-Stokes equations.

    u_t + (u . grad) u + grad p - nu Laplacian(u) = 0,  div u = 0  on (0, 2 pi)^2,  t in [0, 2],  nu = 0.1.
    Exact: u = sin x cos y F,  v = -cos x sin y F,  p = (cos 2x + cos 2y) F^2 / 4,  F = exp(-2 nu t).
    Initial and boundary velocity from the exact solution.
    Outputs: u, v, p. p is defined up to a constant, so its error is measured after removing the mean.
"""
import math
import torch
from pinn import grad, uniform, face, faces, grid

NU = 0.1
L = 2 * math.pi
DOMAIN = [(0.0, L), (0.0, L), (0.0, 2.0)]  # x, y, t
OUT = 3
MEAN_FREE = [2]


def exact(x):
    s, y, t = x[:, :1], x[:, 1:2], x[:, 2:3]
    F = torch.exp(-2 * NU * t)
    return torch.cat([torch.sin(s) * torch.cos(y) * F, -torch.cos(s) * torch.sin(y) * F,
                      (torch.cos(2 * s) + torch.cos(2 * y)) * F ** 2 / 4], dim=1)


# ponytail: mom_x and mom_y each recompute these derivatives (about 2x work); share them if NS time matters
def fields(net, x):
    q = net(x)
    u, v = q[:, 0], q[:, 1]
    du, dv, dp = grad(u, x), grad(v, x), grad(q[:, 2], x)
    lap_u = grad(du[:, 0], x)[:, 0] + grad(du[:, 1], x)[:, 1]
    lap_v = grad(dv[:, 0], x)[:, 0] + grad(dv[:, 1], x)[:, 1]
    return u, v, du, dv, dp, lap_u, lap_v


def mom_x(net, x):
    u, v, du, dv, dp, lap_u, lap_v = fields(net, x)
    return du[:, 2] + u * du[:, 0] + v * du[:, 1] + dp[:, 0] - NU * lap_u


def mom_y(net, x):
    u, v, du, dv, dp, lap_u, lap_v = fields(net, x)
    return dv[:, 2] + u * dv[:, 0] + v * dv[:, 1] + dp[:, 1] - NU * lap_v


def cont(net, x):
    q = net(x)
    return grad(q[:, 0], x)[:, 0] + grad(q[:, 1], x)[:, 1]


def velocity(net, x):
    return net(x)[:, :2] - exact(x)[:, :2]


TERMS = {"mom_x": mom_x, "mom_y": mom_y, "cont": cont, "ic": velocity, "bc": velocity}


def points(g):
    x = uniform(g, 16384, DOMAIN)
    return {"mom_x": x, "mom_y": x.clone(), "cont": x.clone(),
            "ic": face(g, 2048, DOMAIN, 2, 0), "bc": faces(g, 2048, DOMAIN, [0, 1])}


def reference():
    x, shape = grid(torch.linspace(0, L, 41), torch.linspace(0, L, 41), torch.linspace(0, 2, 11))
    return x, shape, exact(x)
