"""Numbers and charts for the Problem 5 post, from results/*.json (standard library only).

Prints the summary table and writes four HTML blocks (legend, SVG, hover script, data
table, no blank lines) that sit directly inside the Markdown post:
  fig-solution.html  the exact solution u(x, t) at four times
  fig-training.html  relative L2 error during training (median of 5 seeds) per weighting
  fig-error.html     RMS over x of the error u - exact at each time t (seed 0)
  fig-weights.html   RBA and BRDR weights of the equation points (seed 0), divided by their
                     mean, averaged in bins of t
Copied from 04-boundary-layer/figure.py; the parts marked PROBLEM are changed.

Run: python scripts/posts/pinn/05-heat/figure.py
"""

import json
import math
import statistics as st
from html import escape
from pathlib import Path

HERE = Path(__file__).parent
RES = HERE / "results"
METHODS = ["ntk", "rba", "brdr"]
LABEL = {"ntk": "NTK", "rba": "RBA", "brdr": "BRDR"}
COLOR = {"ntk": "#2a78d6", "rba": "#eb6834", "brdr": "#1baf7a"}  # fixed for the whole series
# PROBLEM: four gray series for the solution chart, one per time (darker = earlier)
TIMES = [0.0, 0.1, 0.3, 1.0]
for k, (t, c) in enumerate(zip(TIMES, ["#1d1d1f", "#6e6e73", "#8e8e93", "#aeaeb2"])):
    LABEL[f"t{k}"] = f"t = {t:g}"
    COLOR[f"t{k}"] = c
DASH = {}
SEEDS = range(5)

runs = {m: [json.loads((RES / f"{m}-s{s}.json").read_text(encoding="utf-8")) for s in SEEDS] for m in METHODS}
ref = json.loads((RES / "reference.json").read_text(encoding="utf-8"))
X = ref["x"]  # (x, t) pairs, 'ij' order: x outer, t inner
u = [row[0] for row in ref["u"]]
nx, nt = ref["shape"]
xs_grid = [X[i * nt][0] for i in range(nx)]
ts_grid = [X[j][1] for j in range(nt)]
u_rms = math.sqrt(st.mean(v * v for v in u))


def at(field, i, j):  # value at x index i, t index j
    return field[i * nt + j]


SUP = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")


def decade(k):
    return {0: "1", 1: "10"}.get(k, "10" + str(k).translate(SUP))


def sci(v):
    if v == 0:
        return "0"
    e = math.floor(math.log10(abs(v)))
    return f"{v / 10 ** e:.2f}×10{str(e).translate(SUP)}" if not -2 <= e <= 2 else f"{v:.3g}"


def rms(v):
    v = list(v)
    return math.sqrt(sum(a * a for a in v) / len(v))


# --- summary table ---
def adam_end(r):
    return [h for h in r["history"] if h["phase"] == "adam"][-1]


def first_below(r, level):
    return next((h["step"] for h in r["history"] if h["step"] > 0 and h["error"][0] < level), None)


def weight_ratio(r, term):
    """Mean weight of a term / mean weight of the equation term, at the end (NTK: one weight per term)."""
    w = r["history"][-1]["weights"]
    mean = (lambda v: v[0]) if isinstance(w["pde"], list) else (lambda v: v)
    return mean(w[term]) / mean(w["pde"])


summary = {}
for m in METHODS:
    fin = [r["error"][0] for r in runs[m]]
    ad = [adam_end(r)["error"][0] for r in runs[m]]
    summary[m] = {"mean": st.mean(fin), "std": st.stdev(fin), "min": min(fin), "max": max(fin),
                  "adam_median": st.median(ad), "adam_min": min(ad), "adam_max": max(ad),
                  "status": sorted({r["status"] for r in runs[m]}),
                  # PROBLEM: how well the initial and boundary conditions are fitted, and their weights
                  "final": fin, "adam": ad,
                  "first_below_0.01": [first_below(r, 0.01) for r in runs[m]],
                  "first_below_0.001": [first_below(r, 0.001) for r in runs[m]],
                  "rms_residual_final": [{k: math.sqrt(v) for k, v in r["history"][-1]["mse"].items()} for r in runs[m]],
                  "abs_rms_final": [e * u_rms for e in fin],  # relative L2 times the RMS of u
                  "weight_ratio_ic": [weight_ratio(r, "ic") for r in runs[m]],
                  "weight_ratio_bc": [weight_ratio(r, "bc") for r in runs[m]]}
    s = summary[m]
    print(f"{LABEL[m]:5s} final {s['mean']:.2e} ± {s['std']:.1e} (min {s['min']:.2e}, max {s['max']:.2e})"
          f"   after Adam: median {s['adam_median']:.2e} (min {s['adam_min']:.2e}, max {s['adam_max']:.2e})"
          f"  L-BFGS gain (medians) {s['adam_median'] / st.median(fin):.0f}x  {s['status']}")
    for i in SEEDS:
        rr = s["rms_residual_final"][i]
        print(f"      s{i} final {s['final'][i]:.2e}  adam {s['adam'][i]:.2e}  error<0.01 at {s['first_below_0.01'][i]}"
              f"  <0.001 at {s['first_below_0.001'][i]}  residual rms pde {rr['pde']:.1e} ic {rr['ic']:.1e} bc {rr['bc']:.1e}"
              f"  error rms {s['abs_rms_final'][i]:.1e}  weight ic/pde {s['weight_ratio_ic'][i]:.3f} bc/pde {s['weight_ratio_bc'][i]:.3f}")
print(f"RMS of u on the test grid: {u_rms:.4f}")
(HERE / "summary.json").write_text(json.dumps(summary, indent=1), encoding="utf-8")


# --- one line chart, same layout and hover as the earlier posts ---
def line_chart(xs, series, *, logy, y_range, x_ticks, y_ticks, x_title, y_title, aria, caption,
               heads, fmt, vlines=(), table_every=1):
    W, H = 680, 400
    L, R, T, B = 72, 560, 30, 340
    x0, x1 = xs[0], xs[-1]
    y0, y1 = (math.log10(y_range[0]), math.log10(y_range[1])) if logy else y_range

    def px(x):
        return L + (R - L) * (x - x0) / (x1 - x0)

    def py(v):
        v = math.log10(max(v, 1e-300)) if logy else v
        v = min(max(v, y0), y1)
        return T + (B - T) * (y1 - v) / (y1 - y0)

    svg = [f'<svg viewBox="0 0 {W} {H}" role="img" tabindex="0" aria-label="{escape(aria)}" class="viz-plot">']
    for v, text in y_ticks:
        y = py(v)
        svg.append(f'<line x1="{L}" x2="{R}" y1="{y:.1f}" y2="{y:.1f}" class="viz-grid"/>')
        svg.append(f'<text x="{L - 10}" y="{y + 4:.1f}" text-anchor="end" class="viz-tick">{text}</text>')
    for v, text in x_ticks:
        svg.append(f'<text x="{px(v):.1f}" y="{B + 22}" text-anchor="middle" class="viz-tick">{text}</text>')
    for v, text in vlines:
        x = px(v)
        svg.append(f'<line x1="{x:.1f}" x2="{x:.1f}" y1="{T}" y2="{B}" class="viz-leader" stroke-dasharray="4 4"/>')
        svg.append(f'<text x="{x + 6:.1f}" y="{T + 14}" class="viz-label">{escape(text)}</text>')
    svg.append(f'<text x="{(L + R) / 2:.1f}" y="{B + 48}" text-anchor="middle" class="viz-axis-title">{escape(x_title)}</text>')
    svg.append(f'<text x="{L}" y="{T - 16}" class="viz-axis-title">{escape(y_title)}</text>')
    names = list(series)
    for name in names:
        ys = series[name]
        path = " ".join(f"{'M' if i == 0 else 'L'}{px(x):.1f},{py(y):.1f}" for i, (x, y) in enumerate(zip(xs, ys)))
        svg.append(f'<path d="{path}" fill="none" stroke="{COLOR[name]}" stroke-width="2" '
                   f'stroke-linejoin="round" stroke-linecap="round"'
                   + (f' stroke-dasharray="{DASH[name]}"' if name in DASH else "") + "/>")
    # direct labels at the right end, pushed apart so they never overlap
    ys_label = sorted((py(series[n][-1]) + 4, n) for n in names)
    for i in range(1, len(ys_label)):
        if ys_label[i][0] - ys_label[i - 1][0] < 15:
            ys_label[i] = (ys_label[i - 1][0] + 15, ys_label[i][1])
    for y, name in ys_label:
        svg.append(f'<text x="{R + 10}" y="{y:.1f}" class="viz-label">{LABEL[name]}</text>')
    svg.append(f'<line x1="0" x2="0" y1="{T}" y2="{B}" class="viz-crosshair" visibility="hidden"/>')
    svg.append("</svg>")

    rows = [[fmt(series[n][i]) for n in names] for i in range(len(xs))]
    legend = "".join(f'<span class="viz-key"><span class="viz-swatch" style="background:{COLOR[n]}"></span>{LABEL[n]}</span>'
                     for n in names)
    table_head = "<tr><th></th>" + "".join(f"<th>{LABEL[n]}</th>" for n in names) + "</tr>"
    table_body = "".join(f"<tr><td>{escape(heads[i])}</td>" + "".join(f"<td>{v}</td>" for v in rows[i]) + "</tr>"
                         for i in range(0, len(xs), table_every))
    script = (
        "<script>(() => {"
        "const fig = document.currentScript.closest('figure');"
        "const svg = fig.querySelector('svg');"
        "const line = svg.querySelector('.viz-crosshair');"
        "const tip = fig.querySelector('.viz-tooltip');"
        f"const names = {json.dumps([LABEL[n] for n in names])}; const colors = {json.dumps([COLOR[n] for n in names])};"
        f"const rows = {json.dumps(rows)}; const heads = {json.dumps(heads)};"
        f"const xs = {json.dumps([round(px(x), 1) for x in xs])}; const steps = {len(xs) - 1}, width = {W};"
        "let current = -1;"
        "const show = (n) => {"
        " n = Math.max(0, Math.min(steps, n)); current = n;"
        " const x = xs[n];"
        " line.setAttribute('x1', x); line.setAttribute('x2', x); line.setAttribute('visibility', 'visible');"
        " tip.replaceChildren();"
        " const head = document.createElement('div'); head.className = 'viz-tooltip-head'; head.textContent = heads[n]; tip.append(head);"
        " names.forEach((name, i) => {"
        "  const row = document.createElement('div'); row.className = 'viz-tooltip-row';"
        "  const key = document.createElement('span'); key.className = 'viz-swatch'; key.style.background = colors[i];"
        "  const value = document.createElement('strong'); value.textContent = rows[n][i];"
        "  const label = document.createElement('span'); label.textContent = name;"
        "  row.append(key, value, label); tip.append(row);"
        " });"
        " const box = svg.getBoundingClientRect(), frame = fig.getBoundingClientRect();"
        " const xp = box.left - frame.left + x * box.width / width;"
        " tip.hidden = false;"
        " const room = frame.width - tip.offsetWidth - 8;"
        " tip.style.left = Math.max(0, Math.min(room, xp + 12)) + 'px';"
        " tip.style.top = (box.top - frame.top + 8) + 'px';"
        "};"
        "const hide = () => { line.setAttribute('visibility', 'hidden'); tip.hidden = true; current = -1; };"
        "svg.addEventListener('pointermove', (event) => {"
        " const box = svg.getBoundingClientRect();"
        " const x = (event.clientX - box.left) * width / box.width;"
        " let best = 0; xs.forEach((v, i) => { if (Math.abs(v - x) < Math.abs(xs[best] - x)) best = i; }); show(best);"
        "});"
        "svg.addEventListener('pointerleave', hide);"
        "svg.addEventListener('focus', () => show(current < 0 ? 0 : current));"
        "svg.addEventListener('blur', hide);"
        "svg.addEventListener('keydown', (event) => {"
        " if (event.key === 'ArrowRight') { show(current + 1); event.preventDefault(); }"
        " if (event.key === 'ArrowLeft') { show(current - 1); event.preventDefault(); }"
        " if (event.key === 'Escape') hide();"
        "});"
        "})();</script>"
    )
    parts = ['<figure class="viz">', f'<div class="viz-legend">{legend}</div>', '<div class="viz-frame">', "".join(svg),
             '<div class="viz-tooltip" hidden></div>', "</div>", f"<figcaption>{escape(caption)}</figcaption>",
             f'<details class="viz-table"><summary>Show the numbers</summary><table><thead>{table_head}</thead>'
             f"<tbody>{table_body}</tbody></table></details>", script, "</figure>"]
    return "\n".join(parts) + "\n"


def log_ticks(lo, hi):
    return [(10.0 ** k, decade(k)) for k in range(round(math.log10(lo)), round(math.log10(hi)) + 1)]




# --- 1. the exact solution at four times ---  PROBLEM: heat equation
X_TICKS = [(v, f"{v:g}") for v in (0, 0.2, 0.4, 0.6, 0.8, 1)]
t_index = [ts_grid.index(t) for t in TIMES]
sol = {f"t{k}": [at(u, i, j) for i in range(nx)] for k, j in enumerate(t_index)}
html = line_chart(
    xs_grid, sol, logy=False, y_range=(0, 1.25),
    x_ticks=X_TICKS, y_ticks=[(v, f"{v:g}") for v in (0, 0.25, 0.5, 0.75, 1, 1.25)], x_title="x",
    y_title="Exact solution u(x, t) at four times",
    aria="The exact solution of the heat equation at t = 0, 0.1, 0.3 and 1. At t = 0 it is sin(pi x) + "
         "0.5 sin(3 pi x), with two peaks of about 1.08 near x = 0.2 and 0.8 and a dip to 0.5 at x = 0.5. By "
         "t = 0.3 the dip is almost gone, and at t = 1 the solution is a single sine bump of height about 0.37.",
    caption="The exact solution at four times (darker is earlier). The fast part sin(3πx) decays 9 times faster "
            "than the slow part sin(πx), so the solution soon looks like one bump. The predictions of all three "
            "methods (seed 0) are not drawn: on the whole test grid, the largest gap to the exact solution is "
            "below 6×10⁻⁴. Hover or use the arrow keys to read the values.",
    heads=[f"x = {v:.2f}" for v in xs_grid], fmt=lambda v: f"{v:.4f}", table_every=10,
)
(HERE / "fig-solution.html").write_text(html, encoding="utf-8")
gap = {m: max(abs(p - e) for p, e in zip(runs[m][0]["prediction"], u)) for m in METHODS}
where = {m: X[max(range(len(u)), key=lambda i: abs(runs[m][0]["prediction"][i] - u[i]))] for m in METHODS}
print("largest gap, seed 0:", {m: f"{v:.2e} at (x, t) = {where[m]}" for m, v in gap.items()})
assert max(gap.values()) < 6e-4
print("peak of u at t = 0:", f"{max(sol['t0']):.3f}", " at t = 1:", f"{max(sol['t3']):.3f}",
      " u(0.5, 0):", f"{sol['t0'][xs_grid.index(0.5)]:.3f}")

# --- 2. training curves: median over seeds at each logged step ---
steps = [h["step"] for h in runs["ntk"][0]["history"]]
for m in METHODS:
    for r in runs[m]:
        assert [h["step"] for h in r["history"]] == steps, "all runs must log the same steps"
raw = {m: [st.median(r["history"][i]["error"][0] for r in runs[m]) for i in range(len(steps))] for m in METHODS}
adam_steps = runs["ntk"][0]["train"]["adam"]
n_adam = sum(1 for s in steps if s < adam_steps)
# Adam makes the error jump up and down; a rolling median over 11 logs (about 1,000 steps) shows the trend.
# The L-BFGS part is shown as it is.
med = {m: [st.median(v[max(0, i - 5):min(n_adam, i + 6)]) if i < n_adam else v[i] for i in range(len(v))]
       for m, v in raw.items()}
for k in (0, 10, 50, 100, 150, n_adam - 1, n_adam, len(steps) - 1):
    print(f"step {steps[k]}: smoothed median", {m: f"{med[m][k]:.1e}" for m in METHODS})
lo = 10.0 ** math.floor(math.log10(min(min(v) for v in med.values())))
hi = 10.0 ** math.ceil(math.log10(max(max(v) for v in med.values())))
html = line_chart(
    steps, med, logy=True, y_range=(lo, hi),
    x_ticks=[(s, f"{s // 1000}k" if s else "0") for s in range(0, steps[-1] + 1, 4000)],
    y_ticks=log_ticks(lo, hi), x_title="Training step (Adam, then L-BFGS iterations)",
    y_title="Relative L2 error, median of 5 seeds",
    vlines=[(adam_steps, "L-BFGS")],
    aria="Median relative L2 error of the PINN for the 1D heat equation during training, on a log scale, for "
         "NTK, RBA and BRDR weighting. In the Adam phase, NTK falls to about 5e-4 near step 3,400 and then rises to between 1e-3 and 7e-3, RBA falls from "
         "about 3e-2 to 3e-4, and BRDR falls from about 1e-2 to 3e-4. In the L-BFGS phase, NTK drops to about "
         "3e-5 and RBA to about 4e-5, but BRDR stays near 2e-4.",
    caption="Median over 5 seeds of the relative L2 error on the 101 × 101 test grid, logged every 100 steps. In "
            "the Adam part, the line is a rolling median over about 1,000 steps, because the raw error jumps up "
            "and down from one log to the next. The dashed line marks the switch to L-BFGS. Hover over the chart, "
            "or focus it and use the arrow keys, to read the values.",
    heads=[f"step {s:,}" for s in steps], fmt=sci, table_every=10,
)
(HERE / "fig-training.html").write_text(html, encoding="utf-8")

# --- 3. error of seed 0 over time: RMS over x at each t ---  PROBLEM: space-time error as a curve in t
err_t = {m: [rms(at(runs[m][0]["prediction"], i, j) - at(u, i, j) for i in range(nx)) for j in range(nt)]
         for m in METHODS}
for m in METHODS:
    e = [p - q for p, q in zip(runs[m][0]["prediction"], u)]
    sides = [e[i * nt + j] for i in (0, nx - 1) for j in range(nt)]
    inner = [e[i * nt + j] for i in range(1, nx - 1) for j in range(1, nt)]
    print(f"{LABEL[m]} seed 0: RMS error at t = 0 {err_t[m][0]:.1e}, at t = 0.5 {err_t[m][ts_grid.index(0.5)]:.1e}, "
          f"at t = 1 {err_t[m][-1]:.1e}, max over t {max(err_t[m]):.1e} at t = {ts_grid[err_t[m].index(max(err_t[m]))]}, "
          f"min over t {min(err_t[m]):.1e}; on x = 0, 1: {rms(sides):.1e}; inside: {rms(inner):.1e} "
          f"(signed mean {st.mean(inner):.1e})")
lo = 10.0 ** math.floor(math.log10(min(min(v) for v in err_t.values())))
hi = 10.0 ** math.ceil(math.log10(max(max(v) for v in err_t.values())))
T_TICKS = [(v, f"{v:g}") for v in (0, 0.2, 0.4, 0.6, 0.8, 1)]
html = line_chart(
    ts_grid, err_t, logy=True, y_range=(lo, hi), x_ticks=T_TICKS,
    y_ticks=log_ticks(lo, hi), x_title="t", y_title="RMS over x of the error uθ(x, t) − u(x, t), seed 0",
    aria="Error of the trained PINN at each time t from 0 to 1, as the RMS over the 101 x values, on a log scale, "
         "for NTK, RBA and BRDR weighting, seed 0. NTK and RBA stay between about 8e-6 and 5e-5 for all t. BRDR "
         "is about 3e-4 at t = 0, up to 1.6e-4 before t = 0.1, and between about 3e-5 and 9e-5 after t = 0.1.",
    caption="The error of the final network (seed 0) at each of the 101 test times, as the RMS over the 101 test "
            "points in x. The first value (t = 0) is the error in the initial condition. Hover or use the arrow "
            "keys to read the values.",
    heads=[f"t = {v:.2f}" for v in ts_grid], fmt=sci, table_every=10,
)
(HERE / "fig-error.html").write_text(html, encoding="utf-8")

# --- 4. point weights of RBA and BRDR over t, seed 0 ---  PROBLEM: 2D points, so mean in bins of t
W_METHODS = ["rba", "brdr"]
BIN = 0.05
centers = [BIN * (k + 0.5) for k in range(round(1 / BIN))]
wser = {}
for m in W_METHODS:
    p = runs[m][0]["point_weights"]["pde"]
    mean_w = st.mean(p["w"])
    bins = [[] for _ in centers]
    for (_, t), w in zip(p["x"], p["w"]):
        bins[min(int(t / BIN), len(bins) - 1)].append(w / mean_w)
    assert all(bins), "every bin needs points"
    wser[m] = [st.mean(b) for b in bins]
    print(f"{LABEL[m]} seed 0, {len(p['w'])} equation points: weight / mean by t bin",
          [f"{c:.3f}:{v:.2f}" for c, v in zip(centers, wser[m])], "points per bin", [len(b) for b in bins][:4], "...")
    for a, b in [(0, 0.02), (0.02, 0.05), (0.05, 0.1), (0.1, 0.2), (0.2, 0.5), (0.5, 1.01)]:
        ws = [w / mean_w for (_, t), w in zip(p["x"], p["w"]) if a <= t < b]
        print(f"   t in [{a}, {b}): {len(ws)} points, mean {st.mean(ws):.2f}")
    for term in ("ic", "bc"):
        print(f"   {term} weights / mean equation weight: {st.mean(runs[m][0]['point_weights'][term]['w']) / mean_w:.3f}")
hi_w = 3
html = line_chart(
    centers, wser, logy=False, y_range=(0, hi_w), x_ticks=T_TICKS,
    y_ticks=[(v, f"{v:g}") for v in (0, 1, 2, 3)], x_title="t",
    y_title="Mean weight of the equation points ÷ mean weight, seed 0",
    aria="Point weights of the equation residual at the end of training, seed 0, divided by their mean and "
         "averaged in bins of width 0.05 in t. RBA weights are about 2.3 near t = 0 and fall to about 0.75 for "
         "late times. BRDR weights stay between about 0.8 and 1.4 at all times.",
    caption="The weights of the equation points after training (seed 0), divided by the mean weight. The points "
            "are random in (x, t), so each value is the mean over the points in a band of width 0.05 in t (about "
            "100 points per band). NTK has one weight for the whole term, so it is not shown. Hover or use the "
            "arrow keys to read the values.",
    heads=[f"t in [{c - BIN / 2:.2f}, {c + BIN / 2:.2f})" for c in centers], fmt=lambda v: f"{v:.2f}",
)
(HERE / "fig-weights.html").write_text(html, encoding="utf-8")

print("wrote fig-solution.html, fig-training.html, fig-error.html, fig-weights.html")
