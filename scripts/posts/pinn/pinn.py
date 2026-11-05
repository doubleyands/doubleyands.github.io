"""Shared code for the PINN problem posts. Each problem lives in NN-name/problem.py.

Every problem uses the same setup, so the posts compare like with like:

  network    random Fourier features (Gaussian, scale sigma) + MLP, tanh, Xavier init
  weighting  ntk  - one weight per loss term, from the trace of the NTK
                    (Wang, Yu, Perdikaris, J. Comput. Phys. 449, 110768, 2022)
             rba  - one weight per point, residual-based attention
                    (Anagnostopoulos et al., Comput. Methods Appl. Mech. Eng. 421, 116805, 2024)
             brdr - one weight per point, balanced residual decay rate
                    (Chen, Howard, Stinis, J. Comput. Phys. 542, 114226, 2025)
  optimizer  Adam with exponential decay, then L-BFGS with the weights frozen
  precision  float64

Usage (from this folder):
  python pinn.py check                           # exact solutions satisfy each problem; weighting math
  python pinn.py run 01-decay                    # 3 methods x 5 seeds, skips finished runs
  python pinn.py run 01-decay --method rba --seed 0
  python pinn.py run all --shard 3/16            # every 16th job from job 3 (parallel workers)
Results go to NN-name/results/<method>-s<seed>.json and NN-name/results/reference.json.

A problem module defines:
  DOMAIN      [(lo, hi), ...], one pair per input, time last
  OUT         number of outputs
  TERMS       {name: f(net, x) -> residual}, one entry per loss term
  points(g)   {name: x}, the fixed training points of each term (torch.Generator g)
  reference() (x, shape, u): test points, grid shape for plots, reference solution
  optional    exact(x), PARAMS {name: (initial, true)} for inverse problems,
              NET {...}, TRAIN {...}, MEAN_FREE [output indices], check()
Needs only torch.
"""

import argparse
import importlib.util
import json
import math
import os
import sys
import time
from pathlib import Path

import torch
from torch import nn

sys.modules.setdefault("pinn", sys.modules[__name__])  # problems import this file as `pinn`
torch.set_default_dtype(torch.float64)

HERE = Path(__file__).resolve().parent
METHODS = ("ntk", "rba", "brdr")
SEEDS = range(5)
NET = dict(width=128, depth=4, features=128, sigma=1.0)
TRAIN = dict(adam=20_000, lbfgs=2_000, lr=1e-3, decay=0.9, decay_every=2_000, log_every=100)


# --- helpers for problem files ---

def grad(y, x):
    """dy/dx for y of shape (N,) where row i depends only on x[i]. Returns (N, d)."""
    return torch.autograd.grad(y.sum(), x, create_graph=True)[0]


def uniform(g, n, box):
    """n points uniform in a box [(lo, hi), ...]."""
    lo = torch.tensor([a for a, _ in box])
    hi = torch.tensor([b for _, b in box])
    return lo + (hi - lo) * torch.rand(n, len(box), generator=g)


def face(g, n, box, dim, side):
    """n points uniform on the face x[dim] = lo (side 0) or hi (side 1) of a box."""
    x = uniform(g, n, box)
    x[:, dim] = box[dim][side]
    return x


def faces(g, n, box, dims):
    """n points on the faces x[d] = lo or hi for d in dims, spread evenly over the faces."""
    parts = [face(g, n // (2 * len(dims)), box, d, s) for d in dims for s in (0, 1)]
    return torch.cat(parts)


def grid(*axes):
    """Tensor-product grid ('ij' order). Returns (N, d) points and the grid shape."""
    mesh = torch.meshgrid(*[torch.as_tensor(a, dtype=torch.float64) for a in axes], indexing="ij")
    return torch.stack([m.reshape(-1) for m in mesh], dim=1), tuple(len(a) for a in axes)


# --- the network ---

class Net(nn.Module):
    """Inputs scaled to [-1, 1], then random Fourier features [cos(zB), sin(zB)], then an MLP."""

    def __init__(self, domain, out, width, depth, features, sigma):
        super().__init__()
        self.register_buffer("lo", torch.tensor([a for a, _ in domain]))
        self.register_buffer("hi", torch.tensor([b for _, b in domain]))
        self.register_buffer("B", sigma * torch.randn(len(domain), features))
        sizes = [2 * features] + [width] * depth + [out]
        self.layers = nn.ModuleList(nn.Linear(a, b) for a, b in zip(sizes[:-1], sizes[1:]))
        for layer in self.layers:
            nn.init.xavier_normal_(layer.weight)
            nn.init.zeros_(layer.bias)
        self.extra = nn.ParameterDict()  # unknown coefficients of an inverse problem

    def forward(self, x):
        z = (2 * (x - self.lo) / (self.hi - self.lo) - 1) @ self.B
        h = torch.cat([torch.cos(z), torch.sin(z)], dim=1)
        for layer in self.layers[:-1]:
            h = torch.tanh(layer(h))
        return self.layers[-1](h)


class ExactNet(nn.Module):
    """Stands in for the network with the exact solution, to check a problem's residuals."""

    def __init__(self, prob):
        super().__init__()
        self.exact = prob.exact
        self.extra = nn.ParameterDict({k: nn.Parameter(torch.tensor(float(true)))
                                       for k, (_, true) in getattr(prob, "PARAMS", {}).items()})

    def forward(self, x):
        return self.exact(x)


# --- loss weighting ---

class NTKWeights:
    """One weight per term: w_k = sum_j m_j / m_k, where m_k is the mean diagonal of the NTK
    block of term k, i.e. the mean of |d r_i / d theta|^2 over its points. Estimated every
    `every` Adam steps on `n` random points per term, then smoothed (moving average)."""

    def __init__(self, every=1000, n=128, momentum=0.9):
        self.every, self.n, self.momentum = every, n, momentum
        self.w = None

    def update(self, step, net, terms, pts, res):
        if step % self.every:
            return
        params = [p for p in net.parameters() if p.requires_grad]
        m = {}
        for k, x in pts.items():
            # ponytail: n points per term estimate the trace; raise n if the weights look noisy
            idx = torch.randperm(len(x), device=x.device)[: self.n]
            r = terms[k](net, x[idx]).reshape(-1)
            total = 0.0
            for i in range(len(r)):
                gs = torch.autograd.grad(r[i], params, retain_graph=True, allow_unused=True)
                total += sum(float((g * g).sum()) for g in gs if g is not None)
            m[k] = total / len(r)
        s = sum(m.values())
        new = {k: s / max(v, 1e-300) for k, v in m.items()}
        a = self.momentum
        self.w = new if self.w is None else {k: a * self.w[k] + (1 - a) * new[k] for k in new}

    def loss(self, res):
        return sum(self.w[k] * r.pow(2).mean() for k, r in res.items())

    def summary(self):
        return {k: round(v, 6) for k, v in self.w.items()}


class RBAWeights:
    """One weight per point: lambda_i <- gamma lambda_i + eta |r_i| / max_j |r_j|, per term,
    starting from 0. Loss = sum over terms of mean((lambda_i r_i)^2). Bounded by eta / (1 - gamma)."""

    def __init__(self, gamma=0.999, eta=0.01):
        self.gamma, self.eta = gamma, eta
        self.w = None

    def update(self, step, net, terms, pts, res):
        if self.w is None:
            self.w = {k: torch.zeros_like(r.detach()) for k, r in res.items()}
        for k, r in res.items():
            a = r.detach().abs()
            self.w[k] = self.gamma * self.w[k] + self.eta * a / a.max().clamp_min(1e-300)

    def loss(self, res):
        return sum((self.w[k] * r).pow(2).mean() for k, r in res.items())

    def summary(self):
        return {k: [round(float(v.mean()), 6), round(float(v.max()), 6)] for k, v in self.w.items()}


class BRDRWeights:
    """One weight per point, from the inverse residual decay rate
    irdr_i = r_i^2 / sqrt(EMA(r_i^4) + eps), w_i <- beta_w w_i + (1 - beta_w) irdr_i / mean(irdr),
    with the mean over the points of all terms together. Starts at 1.
    Loss = sum over terms of mean(w_i r_i^2)."""

    def __init__(self, beta_c=0.999, beta_w=0.999, eps=1e-14):
        self.beta_c, self.beta_w, self.eps = beta_c, beta_w, eps
        self.w = self.r4 = None
        self.n = 0
        # ponytail: the paper's global scale factor s is left out. Adam is invariant to a
        # (slowly changing) loss scale, so s matters only for plain gradient descent.

    def update(self, step, net, terms, pts, res):
        if self.w is None:
            self.w = {k: torch.ones_like(r.detach()) for k, r in res.items()}
            self.r4 = {k: torch.zeros_like(r.detach()) for k, r in res.items()}
        self.n += 1
        irdr = {}
        for k, r in res.items():
            r2 = r.detach().pow(2)
            self.r4[k] = self.beta_c * self.r4[k] + (1 - self.beta_c) * r2 * r2
            r4 = self.r4[k] / (1 - self.beta_c ** self.n)  # bias correction
            irdr[k] = r2 / torch.sqrt(r4 + self.eps)
        mean = torch.cat(list(irdr.values())).mean().clamp_min(1e-300)
        for k in res:
            self.w[k] = self.beta_w * self.w[k] + (1 - self.beta_w) * irdr[k] / mean

    def loss(self, res):
        return sum((self.w[k] * r.pow(2)).mean() for k, r in res.items())

    def summary(self):
        return {k: [round(float(v.mean()), 6), round(float(v.max()), 6)] for k, v in self.w.items()}


WEIGHTS = {"ntk": NTKWeights, "rba": RBAWeights, "brdr": BRDRWeights}


# --- training ---

def load(name):
    path = HERE / name / "problem.py"
    spec = importlib.util.spec_from_file_location(f"problem_{name.replace('-', '_')}", path)
    prob = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(prob)
    return prob


def problems():
    return sorted(p.parent.name for p in HERE.glob("[0-9][0-9]-*/problem.py"))


def short(v, digits=7):
    return float(f"{v:.{digits}g}")


def relative_errors(net, x, u, mean_free):
    with torch.no_grad():
        p = net(x)
    out = []
    for j in range(u.shape[1]):
        a, b = p[:, j], u[:, j]
        if j in mean_free:  # e.g. pressure, defined up to a constant
            a, b = a - a.mean(), b - b.mean()
        scale = b.norm() if b.norm() > 0 else math.sqrt(len(b))  # zero reference: RMS error
        out.append(float((a - b).norm() / scale))
    return out


def train(name, method, seed, device, adam=None, lbfgs=None, out_dir=None):
    prob = load(name)
    cfg = {**TRAIN, **getattr(prob, "TRAIN", {})}
    if adam is not None:
        cfg["adam"] = adam
    if lbfgs is not None:
        cfg["lbfgs"] = lbfgs
    out_dir = Path(out_dir) if out_dir else HERE / name / "results"
    out_dir.mkdir(parents=True, exist_ok=True)

    torch.manual_seed(seed)
    g = torch.Generator().manual_seed(seed)
    pts = {k: v.to(device).requires_grad_(True) for k, v in prob.points(g).items()}
    assert pts.keys() == prob.TERMS.keys(), f"{name}: points() and TERMS have different keys"
    x_ref, shape, u_ref = prob.reference()
    x_ref, u_ref = x_ref.to(device), u_ref.to(device)
    mean_free = getattr(prob, "MEAN_FREE", [])
    write_reference(out_dir, x_ref, shape, u_ref)

    net = Net(prob.DOMAIN, prob.OUT, **{**NET, **getattr(prob, "NET", {})})
    for k, (init, _) in getattr(prob, "PARAMS", {}).items():
        net.extra[k] = nn.Parameter(torch.tensor(float(init)))
    net = net.to(device)
    params = list(net.parameters())
    weights = WEIGHTS[method]()

    def residuals():
        return {k: prob.TERMS[k](net, x).reshape(-1) for k, x in pts.items()}

    history, start, status = [], time.time(), "ok"

    def log(step, phase, loss, res):
        history.append({
            "step": step, "phase": phase, "loss": short(loss),
            "mse": {k: short(float(r.detach().pow(2).mean())) for k, r in res.items()},
            "error": [short(e) for e in relative_errors(net, x_ref, u_ref, mean_free)],
            "weights": weights.summary(),
            "params": {k: short(float(v)) for k, v in net.extra.items()},
            "time": round(time.time() - start, 1),
        })

    # Adam
    opt = torch.optim.Adam(params, lr=cfg["lr"])
    sched = torch.optim.lr_scheduler.ExponentialLR(opt, cfg["decay"] ** (1 / cfg["decay_every"]))
    for step in range(cfg["adam"]):
        res = residuals()
        weights.update(step, net, prob.TERMS, pts, res)
        loss = weights.loss(res)
        if not torch.isfinite(loss):
            status = f"adam diverged at step {step}"
            break
        if step % cfg["log_every"] == 0:
            log(step, "adam", float(loss.detach()), res)
        opt.zero_grad()
        loss.backward()
        opt.step()
        sched.step()

    # L-BFGS on the same points, weights frozen at their last Adam values
    chunk = cfg["log_every"]
    lb = torch.optim.LBFGS(params, lr=1.0, max_iter=chunk, max_eval=50 * chunk, history_size=50,
                           tolerance_grad=1e-14, tolerance_change=0.0, line_search_fn="strong_wolfe")

    def closure():
        lb.zero_grad()
        loss = weights.loss(residuals())
        loss.backward()
        return loss

    done = 0
    while status == "ok" and done < cfg["lbfgs"]:
        backup = [p.detach().clone() for p in params]
        lb.step(closure)
        res = residuals()
        loss = float(weights.loss(res).detach())
        if not math.isfinite(loss):
            with torch.no_grad():
                for p, b in zip(params, backup):
                    p.copy_(b)
            status = f"lbfgs diverged after {done} iterations, kept the last finite state"
            res = residuals()
            loss = float(weights.loss(res).detach())
        n_iter = lb.state[params[0]]["n_iter"]
        log(cfg["adam"] + n_iter, "lbfgs", loss, res)
        if n_iter - done < chunk:  # stopped early: converged or no progress
            break
        done = n_iter

    result = {
        "problem": name, "method": method, "seed": seed, "status": status,
        "net": {**NET, **getattr(prob, "NET", {})}, "train": cfg,
        "points": {k: len(x) for k, x in pts.items()},
        "error": history[-1]["error"] if history else None,
        "params": {k: float(v) for k, v in net.extra.items()},
        "true_params": {k: t for k, (_, t) in getattr(prob, "PARAMS", {}).items()},
        "time": round(time.time() - start, 1), "device": str(device),
        "history": history,
    }
    if seed == 0:  # the prediction and the point weights, for the figures
        with torch.no_grad():
            result["prediction"] = [short(v, 6) for v in net(x_ref).cpu().reshape(-1).tolist()]
        if method != "ntk":
            result["point_weights"] = {
                k: {"x": [[short(c, 5) for c in row] for row in pts[k][:2000].detach().cpu().tolist()],
                    "w": [short(v, 5) for v in weights.w[k][:2000].cpu().tolist()]}
                for k in pts if len(weights.w[k]) == len(pts[k])
            }
    write_json(out_dir / f"{method}-s{seed}.json", result)
    return result


def write_json(path, data):
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(data, separators=(",", ":")), encoding="utf-8")
    os.replace(tmp, path)  # never leave a half-written result behind


def write_reference(out_dir, x, shape, u):
    path = out_dir / "reference.json"
    if not path.exists():
        write_json(path, {"shape": list(shape),
                          "x": [[short(c, 7) for c in row] for row in x.cpu().tolist()],
                          "u": [[short(c, 7) for c in row] for row in u.cpu().tolist()]})


# --- checks ---

def check():
    torch.manual_seed(0)
    # RBA: a point that always has the largest residual tends to the bound eta / (1 - gamma) = 10
    rba = RBAWeights()
    r = {"pde": torch.tensor([1.0, 2.0, 4.0])}
    for _ in range(20_000):
        rba.update(0, None, None, None, r)
    assert abs(float(rba.w["pde"][2]) - 10) < 1e-6 and abs(float(rba.w["pde"][0]) - 2.5) < 1e-6

    # BRDR: weights keep mean 1, and the point whose residual decays slowest gets the most weight
    brdr = BRDRWeights()
    for n in range(5000):
        brdr.update(n, None, None, None, {"a": torch.tensor([0.999 ** n, 0.99 ** n]),
                                          "b": torch.tensor([0.995 ** n])})
    w = torch.cat([brdr.w["a"], brdr.w["b"]])
    assert abs(float(w.mean()) - 1) < 1e-9 and w[0] > w[2] > w[1], w

    # NTK: a term that is 10 x another has a 100 x larger NTK trace, so 1/100 of its weight
    net = Net([(0.0, 1.0)], 1, 16, 2, 8, 1.0)
    x = torch.rand(64, 1)
    terms = {"a": lambda n, x: n(x)[:, 0], "b": lambda n, x: 10 * n(x)[:, 0]}
    ntk = NTKWeights(n=64)
    ntk.update(0, net, terms, {"a": x, "b": x}, None)
    assert abs(ntk.w["a"] / ntk.w["b"] - 100) < 1e-6, ntk.w

    # Every problem: the exact solution makes every residual (numerically) zero
    for name in problems():
        prob = load(name)
        pts = {k: v.requires_grad_(True) for k, v in prob.points(torch.Generator().manual_seed(0)).items()}
        assert pts.keys() == prob.TERMS.keys(), name
        assert all(v.ndim == 2 and v.shape[1] == len(prob.DOMAIN) for v in pts.values()), name
        x_ref, shape, u_ref = prob.reference()
        assert u_ref.shape == (math.prod(shape), prob.OUT) and torch.isfinite(u_ref).all(), name
        worst = {}
        if hasattr(prob, "exact"):
            net = ExactNet(prob)
            for k, x in pts.items():
                worst[k] = float(prob.TERMS[k](net, x).detach().abs().max())
            assert max(worst.values()) < 1e-6, (name, worst)
        if hasattr(prob, "check"):
            prob.check()
        sizes = {k: len(v) for k, v in pts.items()}
        print(f"ok  {name:28s} points {sizes}  max |residual of exact| "
              + (f"{max(worst.values()):.1e}" if worst else "(numerical reference)"))
    print("all checks passed")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("command", choices=["check", "run"])
    ap.add_argument("problem", nargs="?", default="all")
    ap.add_argument("--method", choices=METHODS)
    ap.add_argument("--seed", type=int)
    ap.add_argument("--shard", default="0/1", help="i/n: run every n-th job starting at job i")
    ap.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    ap.add_argument("--adam", type=int, help="override the number of Adam steps (smoke tests)")
    ap.add_argument("--lbfgs", type=int, help="override the number of L-BFGS iterations")
    ap.add_argument("--out", help="write results here instead of NN-name/results (smoke tests)")
    a = ap.parse_args()
    if a.command == "check":
        return check()

    names = problems() if a.problem == "all" else [a.problem]
    methods = [a.method] if a.method else METHODS
    seeds = [a.seed] if a.seed is not None else SEEDS
    jobs = [(n, m, s) for n in names for m in methods for s in seeds]
    i, n = map(int, a.shard.split("/"))
    for name, method, seed in jobs[i::n]:
        path = (Path(a.out) if a.out else HERE / name / "results") / f"{method}-s{seed}.json"
        if path.exists():
            continue
        r = train(name, method, seed, a.device, a.adam, a.lbfgs, a.out)
        err = ", ".join(f"{e:.2e}" for e in r["error"]) if r["error"] else "-"
        print(f"{name} {method} s{seed}: error {err}  {r['time']:.0f} s  {r['status']}", flush=True)


if __name__ == "__main__":
    main()
