"""Problem 14: nonlinear Schrodinger equation, a moving bright soliton (complex-valued).

    i h_t + 0.5 h_xx + |h|^2 h = 0,  x in [-6, 6],  t in [0, 4].
    Exact: h = a sech(a (x - x0 - v t)) exp(i (v x + (a^2 - v^2) t / 2)),  a = 1.5,  v = 1,  x0 = -2.
    The network has two outputs, the real part u and the imaginary part w of h = u + i w:
        -w_t + 0.5 u_xx + (u^2 + w^2) u = 0,   u_t + 0.5 w_xx + (u^2 + w^2) w = 0.
    Initial and boundary values come from the exact solution.
"""
import torch
from pinn import grad, uniform, face, faces, grid

A, V, X0 = 1.5, 1.0, -2.0
DOMAIN = [(-6.0, 6.0), (0.0, 4.0)]  # x, t
OUT = 2


def exact(x):
    s, t = x[:, :1], x[:, 1:2]
    amp = A / torch.cosh(A * (s - X0 - V * t))
    phase = V * s + (A ** 2 - V ** 2) * t / 2
    return torch.cat([amp * torch.cos(phase), amp * torch.sin(phase)], dim=1)


def parts(net, x):
    h = net(x)
    u, w = h[:, 0], h[:, 1]
    du, dw = grad(u, x), grad(w, x)
    u_xx, w_xx = grad(du[:, 0], x)[:, 0], grad(dw[:, 0], x)[:, 0]
    return u, w, du, dw, u_xx, w_xx


def pde_re(net, x):
    u, w, du, dw, u_xx, w_xx = parts(net, x)
    return -dw[:, 1] + 0.5 * u_xx + (u ** 2 + w ** 2) * u


def pde_im(net, x):
    u, w, du, dw, u_xx, w_xx = parts(net, x)
    return du[:, 1] + 0.5 * w_xx + (u ** 2 + w ** 2) * w


def data(net, x):
    return net(x) - exact(x)


TERMS = {"pde_re": pde_re, "pde_im": pde_im, "ic": data, "bc": data}


def points(g):
    x = uniform(g, 8192, DOMAIN)
    return {"pde_re": x, "pde_im": x.clone(), "ic": face(g, 512, DOMAIN, 1, 0), "bc": faces(g, 256, DOMAIN, [0])}


def reference():
    x, shape = grid(torch.linspace(-6, 6, 241), torch.linspace(0, 4, 101))
    return x, shape, exact(x)
