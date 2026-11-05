"""Problem 16: eikonal equation, first-arrival travel time from a point source.

    c(x, y)^2 |grad T|^2 = 1 on (0, 1)^2,  T(s) = 0,  c = 1 + y (speed grows linearly with y),
    source s = (0.5, 0.25).
    Exact (linear speed): T = arccosh(1 + |x - s|^2 / (2 c(s) c(x))) for gradient 1 of c.
    T is not smooth at the source (a cone tip), and only the source value is given.
"""
import torch
from pinn import grad, uniform, grid

S = (0.5, 0.25)
DOMAIN = [(0.0, 1.0), (0.0, 1.0)]
OUT = 1


def speed(x):
    return 1 + x[:, 1]


def exact(x):
    d2 = (x[:, 0] - S[0]) ** 2 + (x[:, 1] - S[1]) ** 2
    return torch.acosh(1 + d2 / (2 * (1 + S[1]) * speed(x)))[:, None]


def pde(net, x):
    dT = grad(net(x)[:, 0], x)
    return speed(x) ** 2 * (dT ** 2).sum(dim=1) - 1


def source(net, x):
    return net(x)[:, 0]


TERMS = {"pde": pde, "source": source}


def points(g):
    return {"pde": uniform(g, 4096, DOMAIN), "source": torch.tensor([list(S)])}


def reference():
    x, shape = grid(torch.linspace(0, 1, 101), torch.linspace(0, 1, 101))
    return x, shape, exact(x)
