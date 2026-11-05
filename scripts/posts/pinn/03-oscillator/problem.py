"""Problem 3: damped harmonic oscillator.

    u'' + 2 delta u' + w0^2 u = 0,  u(0) = 1,  u'(0) = 0,  t in [0, 1],  delta = 2,  w0 = 20.
    Exact (under-damped): u = exp(-delta t) (cos(w t) + delta / w sin(w t)),  w = sqrt(w0^2 - delta^2).
"""
import math
import torch
from pinn import grad, uniform, grid

DELTA, W0 = 2.0, 20.0
W = math.sqrt(W0 ** 2 - DELTA ** 2)
DOMAIN = [(0.0, 1.0)]
OUT = 1


def exact(x):
    t = x[:, :1]
    return torch.exp(-DELTA * t) * (torch.cos(W * t) + DELTA / W * torch.sin(W * t))


def pde(net, x):
    u = net(x)[:, 0]
    u_t = grad(u, x)[:, 0]
    u_tt = grad(u_t, x)[:, 0]
    return u_tt + 2 * DELTA * u_t + W0 ** 2 * u


def ic(net, x):
    return net(x)[:, 0] - 1.0


def ic_t(net, x):
    return grad(net(x)[:, 0], x)[:, 0]


TERMS = {"pde": pde, "ic": ic, "ic_t": ic_t}


def points(g):
    return {"pde": uniform(g, 1024, DOMAIN), "ic": torch.zeros(1, 1), "ic_t": torch.zeros(1, 1)}


def reference():
    x, shape = grid(torch.linspace(0, 1, 501))
    return x, shape, exact(x)
