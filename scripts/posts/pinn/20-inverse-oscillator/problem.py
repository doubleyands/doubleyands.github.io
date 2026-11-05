"""Problem 20: inverse problem, find the damping and the stiffness of an oscillator from data.

    u'' + mu u' + k u = 0,  t in [0, 1],  mu and k unknown (true values 4 and 400, Problem 3).
    Data: u at 64 random times from the exact solution (no noise). No initial conditions.
    mu = exp(log_mu), k = exp(log_k); they start at mu = 1 and k = 100.
"""
import math
import torch
from pinn import grad, uniform, grid

MU, K = 4.0, 400.0
DELTA = MU / 2
W = math.sqrt(K - DELTA ** 2)
DOMAIN = [(0.0, 1.0)]
OUT = 1
PARAMS = {"log_mu": (0.0, math.log(MU)), "log_k": (math.log(100.0), math.log(K))}


def exact(x):
    t = x[:, :1]
    return torch.exp(-DELTA * t) * (torch.cos(W * t) + DELTA / W * torch.sin(W * t))


def pde(net, x):
    u = net(x)[:, 0]
    u_t = grad(u, x)[:, 0]
    u_tt = grad(u_t, x)[:, 0]
    return u_tt + torch.exp(net.extra["log_mu"]) * u_t + torch.exp(net.extra["log_k"]) * u


def data(net, x):
    return net(x)[:, 0] - exact(x)[:, 0]


TERMS = {"pde": pde, "data": data}


def points(g):
    return {"pde": uniform(g, 1024, DOMAIN), "data": uniform(g, 64, DOMAIN)}


def reference():
    x, shape = grid(torch.linspace(0, 1, 501))
    return x, shape, exact(x)
