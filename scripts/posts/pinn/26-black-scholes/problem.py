"""Problem 26: Black-Scholes equation for a European call option.

    In time to maturity tau:  V_tau = 0.5 s^2 S^2 V_SS + r S V_S - r V,  S in [0, 3],  tau in [0, 1],
    strike K = 1, volatility s = 0.3, rate r = 0.05.
    V(S, 0) = max(S - K, 0) (the payoff has a kink at S = K).
    Exact: the Black-Scholes formula V = S N(d1) - K exp(-r tau) N(d2). The boundary values at
    S = 0 and S = 3 come from it.
"""
import math
import torch
from pinn import grad, uniform, face, faces, grid

K, SIG, R = 1.0, 0.3, 0.05
DOMAIN = [(0.0, 3.0), (0.0, 1.0)]  # S, tau
OUT = 1


def normal_cdf(z):
    return 0.5 * (1 + torch.erf(z / math.sqrt(2)))


def exact(x):
    S, tau = x[:, 0], x[:, 1]
    ts, Ss = tau.clamp_min(1e-12), S.clamp_min(1e-12)  # keep both branches of the where finite
    sq = SIG * torch.sqrt(ts)
    d1 = (torch.log(Ss / K) + (R + SIG ** 2 / 2) * ts) / sq
    V = S * normal_cdf(d1) - K * torch.exp(-R * ts) * normal_cdf(d1 - sq)
    return torch.where(tau > 0, V, torch.clamp(S - K, min=0))[:, None]


def pde(net, x):
    V = net(x)[:, 0]
    dV = grad(V, x)
    V_SS = grad(dV[:, 0], x)[:, 0]
    S = x[:, 0]
    return dV[:, 1] - (0.5 * SIG ** 2 * S ** 2 * V_SS + R * S * dV[:, 0] - R * V)


def data(net, x):
    return net(x)[:, 0] - exact(x)[:, 0]


TERMS = {"pde": pde, "ic": data, "bc": data}


def points(g):
    return {"pde": uniform(g, 8192, DOMAIN), "ic": face(g, 512, DOMAIN, 1, 0), "bc": faces(g, 256, DOMAIN, [0])}


def reference():
    x, shape = grid(torch.linspace(0, 3, 151), torch.linspace(0, 1, 101))
    return x, shape, exact(x)
