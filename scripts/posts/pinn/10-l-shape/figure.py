"""Numbers and charts for the Problem 10 post, from results/*.json (standard library only).

Prints the summary table and writes four HTML blocks (no blank lines) that sit directly
inside the Markdown post:
  fig-solution.html  color maps of the exact solution and the NTK prediction (seed 0)
  fig-training.html  relative L2 error during training (median of 5 seeds) per weighting
  fig-error.html     color maps of the error |u - exact| (seed 0), log scale, one per weighting
  fig-corner.html    RMS error against the distance to the corner (0, 0) (seed 0)
The line chart is copied from 09-poisson-2d/figure.py; the color maps come from ../heatmap.py,
which masks the points outside the L-shape. The parts marked PROBLEM are changed.

Run: python scripts/posts/pinn/10-l-shape/figure.py
"""

import json
import math
import statistics as st
import sys
from html import escape
from pathlib import Path

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.parent))
from heatmap import heatmap  # noqa: E402

RES = HERE / "results"
METHODS = ["ntk", "rba", "brdr"]
LABEL = {"ntk": "NTK", "rba": "RBA", "brdr": "BRDR"}
COLOR = {"ntk": "#2a78d6", "rba": "#eb6834", "brdr": "#1baf7a"}  # fixed for the whole series
DASH = {}
SEEDS = range(5)

runs = {m: [json.loads((RES / f"{m}-s{s}.json").read_text(encoding="utf-8")) for s in SEEDS] for m in METHODS}
ref = json.loads((RES / "reference.json").read_text(encoding="utf-8"))
X = ref["x"]  # PROBLEM: only the points of the 201 x 201 grid that lie in the L-shape
u = [row[0] for row in ref["u"]]
G = [round(-1 + 0.01 * k, 2) for k in range(201)]  # the full grid, 'ij' order: x outer, y inner
index = {(round(x, 2) + 0.0, round(y, 2) + 0.0): n for n, (x, y) in enumerate(X)}
assert len(index) == len(X)
u_rms = math.sqrt(st.mean(v * v for v in u))
pred = {m: [p[0] if isinstance(p, list) else p for p in runs[m][0]["prediction"]] for m in METHODS}
corner = [math.hypot(x, y) for x, y in X]  # distance to the re-entrant corner (0, 0)
OUTLINE = [(-1, -1), (0, -1), (0, 0), (1, 0), (1, 1), (-1, 1)]


def full(field):
    """A field on the L-shape points -> the full grid, None outside."""
    return [None if (i, j) not in index else field[index[(i, j)]] for i in G for j in G]



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


def weight_ratio(r, term, h=-1):
    """Mean weight of a term / mean weight of the equation term (NTK: one weight per term)."""
    w = r["history"][h]["weights"]
    mean = (lambda v: v[0]) if isinstance(w["pde"], list) else (lambda v: v)
    return mean(w[term]) / mean(w["pde"])


summary = {}
CHECK = [0, 1000, 5000]  # PROBLEM: the start, early, and the middle of Adam
TERMS = ("pde", "bc")
for m in METHODS:
    fin = [r["error"][0] for r in runs[m]]
    ad = [adam_end(r)["error"][0] for r in runs[m]]
    summary[m] = {"mean": st.mean(fin), "std": st.stdev(fin), "min": min(fin), "max": max(fin),
                  "adam_median": st.median(ad), "adam_min": min(ad), "adam_max": max(ad),
                  "status": sorted({r["status"] for r in runs[m]}),
                  "final": fin, "adam": ad,
                  "first_below_0.1": [first_below(r, 0.1) for r in runs[m]],
                  "first_below_0.01": [first_below(r, 0.01) for r in runs[m]],
                  "first_below_0.001": [first_below(r, 0.001) for r in runs[m]],
                  "mse_at": {s: [{"error": h["error"][0], **h["mse"]} for r in runs[m]
                                 for h in r["history"] if h["step"] == s] for s in CHECK},
                  "rms_residual_final": [{k: math.sqrt(v) for k, v in r["history"][-1]["mse"].items()} for r in runs[m]],
                  "rms_residual_adam": [{k: math.sqrt(v) for k, v in adam_end(r)["mse"].items()} for r in runs[m]],
                  "weight_ratio": {k: [weight_ratio(r, k) for r in runs[m]] for k in TERMS[1:]},
                  "weight_ratio_1000": {k: [weight_ratio(r, k, 10) for r in runs[m]] for k in TERMS[1:]}}
    s = summary[m]
    print(f"{LABEL[m]:5s} final {s['mean']:.2e} ± {s['std']:.1e} (min {s['min']:.2e}, max {s['max']:.2e})"
          f"   after Adam: median {s['adam_median']:.2e} (min {s['adam_min']:.2e}, max {s['adam_max']:.2e})"
          f"  L-BFGS gain (medians) {s['adam_median'] / st.median(fin):.1f}x  {s['status']}")
    for i in SEEDS:
        rr = s["rms_residual_final"][i]
        print(f"      s{i} final {s['final'][i]:.2e}  adam {s['adam'][i]:.2e}  error<0.1 at {s['first_below_0.1'][i]}"
              f"  <0.01 at {s['first_below_0.01'][i]}  <0.001 at {s['first_below_0.001'][i]}  residual rms "
              + "  ".join(f"{k} {rr[k]:.1e}" for k in TERMS)
              + "  weight/pde " + "  ".join(f"{k} {s['weight_ratio'][k][i]:.3g}" for k in TERMS[1:]))
    for key in ("rms_residual_adam", "rms_residual_final"):
        for k in TERMS:
            v = [rr[k] for rr in s[key]]
            print(f"      {key} {k}: {min(v):.1e}-{max(v):.1e}")
    for key in ("weight_ratio_1000", "weight_ratio"):
        for k in TERMS[1:]:
            v = s[key][k]
            print(f"      {key} {k}/pde: {min(v):.3g}-{max(v):.3g}")
    for step in CHECK:
        q = s["mse_at"][step]
        print(f"      step {step}: error {min(a['error'] for a in q):.3g}-{max(a['error'] for a in q):.3g}  mse " +
              "  ".join(f"{k} {min(a[k] for a in q):.1e}-{max(a[k] for a in q):.1e}" for k in TERMS))
print(f"RMS of u on the test grid: {u_rms:.4f}, min {min(u):.4f}, max {max(u):.4f}")
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


ARIA_TRAINING = (  # PROBLEM: written after reading the printed numbers
    "Median relative L2 error of the PINN for the Laplace equation on an L-shape during training, on a log scale, "
    "for NTK, RBA and BRDR weighting. NTK is the lowest during all of Adam and ends Adam near 2.9e-3. BRDR stays "
    "1.3 to 1.8 times above NTK and ends Adam near 3.4e-3. RBA is the highest and ends Adam near 6.8e-3. In the "
    "L-BFGS part, NTK drops to about 7.3e-4 and BRDR to about 9.1e-4, while RBA ends near 5.5e-3.")
ARIA_ERROR = (
    "Three color maps of the error of seed 0 on the L-shaped domain, one per weighting, on a shared log color "
    "scale from about 2e-5 to 5e-2. The NTK and BRDR maps are light or medium on most of the domain, with one "
    "dark spot at the inner corner (0, 0). The RBA map is dark on almost the whole domain. Thin light lines mark "
    "where the error changes sign.")
ARIA_CORNER = (
    "RMS error of seed 0 against the distance r to the inner corner, on a log scale, for NTK, RBA and BRDR. Near "
    "the corner, all three are about 1e-2. NTK falls to 9.0e-4 near r = 0.1 and 1.1e-4 near r = 1, and BRDR falls "
    "in a similar way. RBA stays high: 4.9e-3 near r = 0.5 and 2.5e-3 near r = 1.")


# --- where the error is, seed 0 ---  PROBLEM: distance to the corner, and the edges
def on_edge(x, y):
    return (abs(abs(x) - 1) < 1e-9 or abs(abs(y) - 1) < 1e-9 or (abs(x) < 1e-9 and y <= 1e-9)
            or (abs(y) < 1e-9 and x >= -1e-9))


edge = [on_edge(x, y) for x, y in X]
print(f"test points {len(X)}, on the edges {sum(edge)}, within 0.1 of the corner {sum(r < 0.1 for r in corner)} "
      f"({sum(r < 0.1 for r in corner) / len(X):.4f}); u max {max(u):.4f} at r = {corner[u.index(max(u))]:.3f}")
for m in METHODS:
    e = [p - q for p, q in zip(pred[m], u)]
    sq = sum(v * v for v in e)
    k = max(range(len(e)), key=lambda i: abs(e[i]))
    print(f"{LABEL[m]} seed 0: largest |error| {abs(e[k]):.2e} at (x, y) = ({X[k][0]:.2f}, {X[k][1]:.2f}), "
          f"signed {e[k]:.2e}; RMS on the edges {rms(v for v, b in zip(e, edge) if b):.1e}; "
          f"RMS where r >= 0.2 {rms(v for v, r in zip(e, corner) if r >= 0.2):.1e}; share of squared error "
          + ", ".join(f"r < {c}: {sum(v * v for v, r in zip(e, corner) if r < c) / sq:.2f}" for c in (0.05, 0.1, 0.2, 0.5))
          + f"; share of points with |error| < 1e-3 {sum(abs(v) < 1e-3 for v in e) / len(e):.2f}")
    near = [(r, v) for r, v in zip(corner, e) if 0 < r < 0.1]
    print(f"      mean signed error for 0 < r < 0.1: {st.mean(v for _, v in near):.2e}, share > 0 "
          f"{sum(v > 0 for _, v in near) / len(near):.2f}; error at the corner {e[index[(0.0, 0.0)]]:.2e}")

# --- 1. exact solution and NTK prediction ---  PROBLEM: u >= 0, so the blue map from 0 to max u
umax = max(u)
html = heatmap(
    [("Exact u", full(u)), ("NTK prediction, seed 0", full(pred["ntk"]))], xs=G, ys=G, cmap="blue",
    vrange=(0, umax), bar_ticks=[(v, f"{v:g}") for v in (0, 0.25, 0.5, 0.75, 1, 1.25)], bar_title="u",
    axis_ticks=(-1, 0, 1), outline=OUTLINE, hover_every=4,
    aria="Two color maps on the L-shaped domain: the exact solution and the NTK prediction for seed 0. Both look "
         "the same. The value is 0 on the two edges that meet at the corner (0, 0), and it grows toward the "
         "top-left corner, where it reaches about 1.26. Near (0, 0) the color changes very fast.",
    caption="The exact solution u = r^(2/3) sin(2θ/3) (left) and the NTK prediction of seed 0 (right), on the "
            "points of the 201 × 201 test grid that lie in the L-shape. Darker blue is a larger value. The missing "
            "square at the bottom right is outside the domain. Hover, or focus the chart and use the arrow keys, "
            "to read the values.",
    fmt=lambda v: f"{v:.4f}",
)
(HERE / "fig-solution.html").write_text(html, encoding="utf-8")

# --- 2. training curves: median over seeds at each logged step ---
steps = [h["step"] for h in runs["ntk"][0]["history"]]
for m in METHODS:
    for r in runs[m]:
        assert [h["step"] for h in r["history"]] == steps, "all runs must log the same steps"
raw = {m: [st.median(r["history"][i]["error"][0] for r in runs[m]) for i in range(len(steps))] for m in METHODS}
adam_steps = runs["ntk"][0]["train"]["adam"]
n_adam = sum(1 for s in steps if s < adam_steps)
# Adam makes the error jump up and down; a rolling median over 11 logs (about 1,000 steps) shows the trend.
med = {m: [st.median(v[max(0, i - 5):min(n_adam, i + 6)]) if i < n_adam else v[i] for i in range(len(v))]
       for m, v in raw.items()}
for k in (0, 2, 5, 10, 20, 30, 40, 50, 60, 80, 100, 120, 150, 190, n_adam - 1, n_adam, n_adam + 5, len(steps) - 1):
    print(f"step {steps[k]}: smoothed median", {m: f"{med[m][k]:.2e}" for m in METHODS})
lo = 10.0 ** math.floor(math.log10(min(min(v) for v in med.values())))
hi = 10.0 ** math.ceil(math.log10(max(max(v) for v in med.values())))
html = line_chart(
    steps, med, logy=True, y_range=(lo, hi),
    x_ticks=[(s, f"{s // 1000}k" if s else "0") for s in range(0, steps[-1] + 1, 4000)],
    y_ticks=log_ticks(lo, hi), x_title="Training step (Adam, then L-BFGS iterations)",
    y_title="Relative L2 error, median of 5 seeds",
    vlines=[(adam_steps, "L-BFGS")],
    aria=ARIA_TRAINING,
    caption="Median over 5 seeds of the relative L2 error on the test points, logged every 100 steps. In the Adam "
            "part, the line is a rolling median over about 1,000 steps, because the raw error jumps up and down "
            "from one log to the next. The dashed line marks the switch to L-BFGS. Hover over the chart, or focus "
            "it and use the arrow keys, to read the values.",
    heads=[f"step {s:,}" for s in steps], fmt=sci, table_every=10,
)
(HERE / "fig-training.html").write_text(html, encoding="utf-8")

# --- 3. error maps, seed 0, orange, log scale shared by the three panels ---
err = {m: [abs(p - q) for p, q in zip(pred[m], u)] for m in METHODS}
pool = sorted(v for m in METHODS for v in err[m])
e_lo = pool[len(pool) // 20]  # 5th percentile
e_hi = pool[-1]
print(f"error color scale: {e_lo:.1e} (5th percentile) to {e_hi:.1e}")
k0, k1 = math.ceil(math.log10(e_lo)), math.floor(math.log10(e_hi))
html = heatmap(
    [(LABEL[m], full(err[m])) for m in METHODS], xs=G, ys=G, cmap="orange", log=True, vrange=(e_lo, e_hi),
    bar_ticks=[(10.0 ** k, decade(k)) for k in range(k0, k1 + 1)], bar_title="|uθ − u|, seed 0 (log scale)",
    axis_ticks=(-1, 0, 1), outline=OUTLINE, hover_every=8,
    aria=ARIA_ERROR,
    caption="The error |uθ − u| of the final network (seed 0) on the test points, for each weighting. All three "
            "maps use the same log color scale; its low end is the 5th percentile of all values, and smaller "
            "errors get the lightest color. Darker is a larger error. Hover, or focus the chart and use the arrow "
            "keys, to read the values.",
    fmt=sci,
)
(HERE / "fig-error.html").write_text(html, encoding="utf-8")

# --- 4. RMS error against the distance to the corner, seed 0 ---  PROBLEM: bins of width 0.02
width = 0.02
nb = math.ceil(max(corner) / width)
bins = [[n for n, r in enumerate(corner) if b * width <= r < (b + 1) * width] for b in range(nb)]
keep = [b for b in range(nb) if bins[b]]
centers = [round((b + 0.5) * width, 3) for b in keep]
by_r = {m: [rms(err[m][n] for n in bins[b]) for b in keep] for m in METHODS}
for m in METHODS:
    print(LABEL[m], "seed 0 RMS error in the bins at r = 0.01, 0.03, 0.05, 0.11, 0.21, 0.51, 1.01, 1.41:",
          [f"{by_r[m][centers.index(c)]:.1e}" for c in (0.01, 0.03, 0.05, 0.11, 0.21, 0.51, 1.01, 1.41) if c in centers])
lo = 10.0 ** math.floor(math.log10(min(min(v) for v in by_r.values())))
hi = 10.0 ** math.ceil(math.log10(max(max(v) for v in by_r.values())))
html = line_chart(
    centers, by_r, logy=True, y_range=(lo, hi), x_ticks=[(v, f"{v:g}") for v in (0.01, 0.4, 0.8, 1.2)],
    y_ticks=log_ticks(lo, hi), x_title="Distance r to the corner (0, 0)", y_title="RMS of |uθ − u|, seed 0",
    aria=ARIA_CORNER,
    caption="The error of the final network (seed 0), as the RMS over the test points in rings of width 0.02 "
            "around the corner (0, 0). The first ring, 0 ≤ r < 0.02, holds the corner and its nearest neighbors. "
            "Hover or use the arrow keys to read the values.",
    heads=[f"r = {c - width / 2:.2f} to {c + width / 2:.2f}" for c in centers], fmt=sci, table_every=5,
)
(HERE / "fig-corner.html").write_text(html, encoding="utf-8")

# --- point weights of RBA and BRDR, seed 0: only printed ---
for m in ["rba", "brdr"]:
    pw = runs[m][0]["point_weights"]
    mean_w = st.mean(pw["pde"]["w"])
    edges = [0, 0.05, 0.1, 0.2, 0.4, 0.7, 1.5]
    rr = [math.hypot(*q) for q in pw["pde"]["x"]]
    means = [st.mean([w / mean_w for r, w in zip(rr, pw["pde"]["w"]) if a <= r < b] or [float("nan")])
             for a, b in zip(edges, edges[1:])]
    counts = [sum(a <= r < b for r in rr) for a, b in zip(edges, edges[1:])]
    br = [math.hypot(*q) for q in pw["bc"]["x"]]
    near = [w for r, w in zip(br, pw["bc"]["w"]) if r < 0.1]
    print(f"{LABEL[m]} seed 0 equation weight / mean, rings {edges}: {[f'{v:.2f}' for v in means]} (points {counts})")
    print(f"{LABEL[m]} seed 0 edge weight / equation mean: all {st.mean(pw['bc']['w']) / mean_w:.3g}, within 0.1 of "
          f"the corner {st.mean(near) / mean_w:.3g} ({len(near)} points); equation weight min-max / mean "
          f"{min(pw['pde']['w']) / mean_w:.2f}-{max(pw['pde']['w']) / mean_w:.2f}")

assert "TBD" not in (ARIA_TRAINING + ARIA_ERROR + ARIA_CORNER), "write the aria texts"
print("wrote fig-solution.html, fig-training.html, fig-error.html, fig-corner.html")
