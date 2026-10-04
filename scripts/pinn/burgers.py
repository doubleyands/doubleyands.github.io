"""Trains a physics-informed neural network (PINN) on Burgers' equation and
records how its prediction changes during training, for the About page.

    u_t + u u_x = nu u_xx,   x in [-1, 1], t in [0, 1],   nu = 0.01 / pi
    u(x, 0) = -sin(pi x),    u(-1, t) = u(1, t) = 0

The network only sees the equation, the initial condition and the boundary
conditions. It never sees the solution. The exact solution (Cole-Hopf) is
used only to measure the error.

Step 1 (train): python scripts/pinn/burgers.py train
    Writes scripts/pinn/out/run.pt (parameters at saved steps and the loss).
Step 2 (export): python scripts/pinn/burgers.py export
    Picks the frames the page shows and writes public/pinn/burgers.json.

Needs numpy and torch (CPU is enough).
"""

import json
import math
import os
import sys
import time
from pathlib import Path

# NumPy (MKL) and PyTorch each bring an OpenMP runtime; on Windows the two clash
# and the process aborts. Keeping MKL single-threaded avoids loading the second one.
os.environ.setdefault("MKL_THREADING_LAYER", "SEQUENTIAL")

import numpy as np
import torch
from torch import nn

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).with_name("out")
RUN = OUT / "run.pt"
EXPORT = ROOT / "public" / "pinn" / "burgers.json"

NU = 0.01 / math.pi
T_SLICE = 0.75  # the time shown on the page
LAYERS = [2] + [20] * 8 + [1]
N_F, N_0, N_B = 10_000, 256, 256  # collocation, initial and boundary points
STEPS = 20_000
LR, LR_END = 1e-3, 5e-5  # Adam learning rate, decayed exponentially to LR_END
DENSE = 200  # keep every step up to here (the prediction changes fast early on)
SAVE_EVERY = 10  # then every SAVE_EVERY steps
FRAMES = 160  # frames the page plays, spaced evenly in log(step)
SEED = 0


# --- exact solution (Cole-Hopf transform, trapezoid rule over eta) ---

trapezoid = getattr(np, "trapezoid", None) or np.trapz  # renamed in NumPy 2

def exact(t, x, n=8001):
    """u(t, x) for t > 0, as a 1-D array over x."""
    s = math.sqrt(4 * NU * t)
    eta = np.linspace(-10 * s, 10 * s, n)
    y = x[:, None] - eta[None, :]
    log_w = -np.cos(np.pi * y) / (2 * np.pi * NU) - eta[None, :] ** 2 / (4 * NU * t)
    w = np.exp(log_w - log_w.max(axis=1, keepdims=True))
    return -trapezoid(np.sin(np.pi * y) * w, eta, axis=1) / trapezoid(w, eta, axis=1)


# --- the network ---

class MLP(nn.Module):
    def __init__(self):
        super().__init__()
        self.linears = nn.ModuleList(nn.Linear(a, b) for a, b in zip(LAYERS[:-1], LAYERS[1:]))
        for layer in self.linears:
            nn.init.xavier_normal_(layer.weight)
            nn.init.zeros_(layer.bias)

    def forward(self, x, t):
        h = torch.cat([x, 2 * t - 1], dim=1)  # both inputs in [-1, 1]
        for layer in self.linears[:-1]:
            h = torch.tanh(layer(h))
        return self.linears[-1](h)


def flat(model):
    return torch.cat([p.detach().reshape(-1) for p in model.parameters()]).clone()


def load_flat(model, vector):
    torch.nn.utils.vector_to_parameters(vector, model.parameters())


def train():
    torch.manual_seed(SEED)
    rng = np.random.default_rng(SEED)
    tensor = lambda a: torch.tensor(a, dtype=torch.float32).reshape(-1, 1)

    x_f = tensor(rng.uniform(-1, 1, N_F)).requires_grad_(True)
    t_f = tensor(rng.uniform(0, 1, N_F)).requires_grad_(True)
    x_0 = tensor(rng.uniform(-1, 1, N_0))
    u_0 = -torch.sin(math.pi * x_0)
    t_b = tensor(rng.uniform(0, 1, N_B))
    x_b = tensor(np.where(rng.random(N_B) < 0.5, -1.0, 1.0))

    model = MLP()
    optimizer = torch.optim.Adam(model.parameters(), lr=LR)
    scheduler = torch.optim.lr_scheduler.ExponentialLR(optimizer, (LR_END / LR) ** (1 / STEPS))

    def loss_fn():
        u = model(x_f, t_f)
        u_x, u_t = torch.autograd.grad(u.sum(), (x_f, t_f), create_graph=True)
        u_xx = torch.autograd.grad(u_x.sum(), x_f, create_graph=True)[0]
        residual = u_t + u * u_x - NU * u_xx
        return (
            residual.pow(2).mean()
            + (model(x_0, torch.zeros_like(x_0)) - u_0).pow(2).mean()
            + model(x_b, t_b).pow(2).mean()
        )

    saved, params, losses = [0], [flat(model)], []
    start = time.time()
    for step in range(1, STEPS + 1):
        optimizer.zero_grad()
        loss = loss_fn()
        loss.backward()
        optimizer.step()
        scheduler.step()
        losses.append(loss.item())
        if step <= DENSE or step % SAVE_EVERY == 0:
            saved.append(step)
            params.append(flat(model))
        if step % 1000 == 0:
            print(f"step {step:6d}  loss {loss.item():.3e}  {time.time() - start:.0f} s", flush=True)

    OUT.mkdir(exist_ok=True)
    torch.save({"steps": torch.tensor(saved), "params": torch.stack(params), "losses": torch.tensor(losses)}, RUN)
    print(f"saved {RUN}")


def export():
    # The exact values first: on Windows, NumPy's OpenMP runtime must start before PyTorch's
    x = np.linspace(-1, 1, 201)
    u_exact = exact(T_SLICE, x)

    # Error grid over the whole domain (t > 0; t = 0 is the initial condition)
    gx = np.linspace(-1, 1, 201)
    gt = np.linspace(0.01, 1, 100)
    grid_exact = np.stack([exact(t, gx) for t in gt])

    run = torch.load(RUN, weights_only=True)
    saved, params, losses = run["steps"].numpy(), run["params"], run["losses"].numpy()
    model = MLP()
    GX, GT = np.meshgrid(gx, gt)
    gx_t = torch.tensor(GX.reshape(-1, 1), dtype=torch.float32)
    gt_t = torch.tensor(GT.reshape(-1, 1), dtype=torch.float32)
    x_t = torch.tensor(x.reshape(-1, 1), dtype=torch.float32)
    t_t = torch.full_like(x_t, T_SLICE)

    def evaluate(k):
        load_flat(model, params[k])
        with torch.no_grad():
            u = model(x_t, t_t).numpy().ravel()
            g = model(gx_t, gt_t).numpy().reshape(GX.shape)
        error = np.linalg.norm(g - grid_exact) / np.linalg.norm(grid_exact)
        return u, error

    # Frames at steps spaced evenly in log(step): the prediction changes fast at first
    # and slowly later, and the page shows the real step number for every frame.
    targets = np.concatenate([[0], np.geomspace(1, STEPS, FRAMES - 1)])
    picks = sorted({int(np.abs(saved - target).argmin()) for target in targets})

    frames = []
    for k in picks:
        u, error = evaluate(k)
        step = int(saved[k])
        frames.append({
            "step": step,
            "loss": float(losses[max(step - 1, 0)]),
            "error": float(error),
            "u": [round(float(v), 3) for v in u],
        })

    # Moments the page describes, measured on the t = 0.75 slice. Each is null if it
    # does not happen in this run, and the page then leaves out that sentence.
    smooth = np.abs(x) >= 0.2  # away from the shock at x = 0
    mid = (x[1:] + x[:-1]) / 2
    near = np.abs(mid) < 0.3

    def shock(u):
        # Position and steepness of the steepest descent near x = 0
        slope = np.where(near, np.diff(u) / np.diff(x), np.inf)
        i = int(np.argmin(slope))
        return float(mid[i]), float(-slope[i])

    exact_slope = shock(u_exact)[1]
    rows = []
    for f in frames:
        u = np.array(f["u"])
        position, slope = shock(u)
        smooth_error = np.linalg.norm((u - u_exact)[smooth]) / np.linalg.norm(u_exact[smooth])
        rows.append((f["step"], smooth_error, position, slope))

    first = lambda test: next((step for step, *row in rows if test(*row)), None)
    smooth_learned = first(lambda e, x0, k: e < 0.1)
    # The shock "slips" when a formed shock (at least half as steep as the exact one)
    # sits more than 0.02 away from x = 0
    off = [step for step, e, x0, k in rows if k > exact_slope / 2 and abs(x0) > 0.02]
    shock_off = {"from": off[0], "to": off[-1]} if off else None
    shock_back = next((step for step, *_ in rows if shock_off and step > shock_off["to"]), None)
    under = next((f["step"] for f in frames if f["error"] < 0.01), None)

    result = {
        "nu": NU,
        "t": T_SLICE,
        "steps": STEPS,
        "x": [round(float(v), 3) for v in x],
        "exact": [round(float(v), 4) for v in u_exact],
        "events": {
            "smoothLearned": smooth_learned,
            "shockOff": shock_off,
            "shockBack": shock_back,
            "underOnePercent": under,
        },
        "frames": frames,
    }
    EXPORT.parent.mkdir(parents=True, exist_ok=True)
    EXPORT.write_text(json.dumps(result, separators=(",", ":")) + "\n", encoding="utf-8")
    print(f"{len(frames)} frames, final error {frames[-1]['error']:.3e}, events {result['events']}")
    print(f"wrote {EXPORT} ({EXPORT.stat().st_size / 1024:.0f} KB)")


if __name__ == "__main__":
    {"train": train, "export": export}[sys.argv[1] if len(sys.argv) > 1 else "train"]()
