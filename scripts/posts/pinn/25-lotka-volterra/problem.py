"""Problem 25: Lotka-Volterra predator-prey system (two coupled nonlinear ODEs).

    x' = a x - b x y,  y' = d x y - c y,  x(0) = y(0) = 1,  t in [0, 12],
    a = 2/3, b = 4/3, c = 1, d = 1.
    No closed form: the reference is classical RK4 with step 1e-3. The orbit is closed, and
    V = d x - c ln x + b y - a ln y stays constant along it (used in the check).
    Outputs: x (prey) and y (predator).
"""
import math
import torch
from pinn import grad, uniform, grid

A, B, C, D = 2 / 3, 4 / 3, 1.0, 1.0
T_END = 12.0
DOMAIN = [(0.0, T_END)]
OUT = 2


def pde_x(net, t):
    q = net(t)
    return grad(q[:, 0], t)[:, 0] - (A * q[:, 0] - B * q[:, 0] * q[:, 1])


def pde_y(net, t):
    q = net(t)
    return grad(q[:, 1], t)[:, 0] - (D * q[:, 0] * q[:, 1] - C * q[:, 1])


def ic(net, t):
    return net(t) - 1.0


TERMS = {"pde_x": pde_x, "pde_y": pde_y, "ic": ic}


def points(g):
    t = uniform(g, 2048, DOMAIN)
    return {"pde_x": t, "pde_y": t.clone(), "ic": torch.zeros(1, 1)}


def rk4(h, every):
    """States at t = 0, every*h, 2*every*h, ... up to T_END."""
    def f(x, y):
        return A * x - B * x * y, D * x * y - C * y
    x, y = 1.0, 1.0
    out = [(x, y)]
    for n in range(1, round(T_END / h) + 1):
        k1 = f(x, y)
        k2 = f(x + h / 2 * k1[0], y + h / 2 * k1[1])
        k3 = f(x + h / 2 * k2[0], y + h / 2 * k2[1])
        k4 = f(x + h * k3[0], y + h * k3[1])
        x += h / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
        y += h / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
        if n % every == 0:
            out.append((x, y))
    return torch.tensor(out)


def reference():
    x, shape = grid(torch.linspace(0, T_END, 601))
    return x, shape, rk4(1e-3, 20)


def check():
    u = rk4(1e-3, 20)
    V = D * u[:, 0] - C * torch.log(u[:, 0]) + B * u[:, 1] - A * torch.log(u[:, 1])
    assert float(V.max() - V.min()) < 1e-10, float(V.max() - V.min())
    assert float((u - rk4(5e-4, 40)).abs().max()) < 1e-10  # halving the step changes nothing
