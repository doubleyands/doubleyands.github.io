"""Problem 24: 2D linear elasticity (Navier-Cauchy equations, displacement form).

    mu Laplacian(u) + (lam + mu) grad(div u) + f = 0 on (0, 1)^2,  u = 0 on the boundary,
    lam = 1, mu = 0.5.
    Exact: u1 = sin(pi x) sin(pi y),  u2 = 0.5 sin(2 pi x) sin(pi y); the body force f follows from it.
    Outputs: the two displacement components.
"""
import math
import torch
from pinn import grad, uniform, faces, grid

LAM, MU = 1.0, 0.5
P = math.pi
DOMAIN = [(0.0, 1.0), (0.0, 1.0)]
OUT = 2


def exact(x):
    s, y = x[:, :1], x[:, 1:2]
    return torch.cat([torch.sin(P * s) * torch.sin(P * y), 0.5 * torch.sin(2 * P * s) * torch.sin(P * y)], dim=1)


def force(x):
    s, y = x[:, 0], x[:, 1]
    u1 = torch.sin(P * s) * torch.sin(P * y)
    u2 = 0.5 * torch.sin(2 * P * s) * torch.sin(P * y)
    div_x = -P ** 2 * torch.sin(P * s) * torch.sin(P * y) + P ** 2 * torch.cos(2 * P * s) * torch.cos(P * y)
    div_y = P ** 2 * torch.cos(P * s) * torch.cos(P * y) - 0.5 * P ** 2 * torch.sin(2 * P * s) * torch.sin(P * y)
    f1 = -(MU * (-2 * P ** 2) * u1 + (LAM + MU) * div_x)
    f2 = -(MU * (-5 * P ** 2) * u2 + (LAM + MU) * div_y)
    return f1, f2


def navier(net, x):
    q = net(x)
    d1, d2 = grad(q[:, 0], x), grad(q[:, 1], x)
    lap1 = grad(d1[:, 0], x)[:, 0] + grad(d1[:, 1], x)[:, 1]
    lap2 = grad(d2[:, 0], x)[:, 0] + grad(d2[:, 1], x)[:, 1]
    div = grad(d1[:, 0] + d2[:, 1], x)
    f1, f2 = force(x)
    return MU * lap1 + (LAM + MU) * div[:, 0] + f1, MU * lap2 + (LAM + MU) * div[:, 1] + f2


def eq_x(net, x):
    return navier(net, x)[0]


def eq_y(net, x):
    return navier(net, x)[1]


def bc(net, x):
    return net(x)


TERMS = {"eq_x": eq_x, "eq_y": eq_y, "bc": bc}


def points(g):
    x = uniform(g, 4096, DOMAIN)
    return {"eq_x": x, "eq_y": x.clone(), "bc": faces(g, 512, DOMAIN, [0, 1])}


def reference():
    x, shape = grid(torch.linspace(0, 1, 101), torch.linspace(0, 1, 101))
    return x, shape, exact(x)
