"""Numbers and charts for the Problem 4 post, from results/*.json (standard library only).

Prints the summary table and writes four HTML blocks (legend, SVG, hover script, data
table, no blank lines) that sit directly inside the Markdown post:
  fig-solution.html  the exact solution and the outer solution u = x
  fig-training.html  relative L2 error during training (median of 5 seeds) per weighting
  fig-error.html     pointwise error |u - exact| over x (seed 0)
  fig-weights.html   RBA and BRDR weights of the equation points (seed 0), divided by their mean
Copied from 03-oscillator/figure.py; the parts marked PROBLEM are changed.

Run: python scripts/posts/pinn/04-boundary-layer/figure.py
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
# PROBLEM: two gray series for the solution chart (the exact solution is gray and dashed in the series)
LABEL.update({"exact": "u", "outer": "u = x"})
COLOR.update({"exact": "#6e6e73", "outer": "#aeaeb2"})
DASH = {"exact": "6 4"}
SEEDS = range(5)

runs = {m: [json.loads((RES / f"{m}-s{s}.json").read_text(encoding="utf-8")) for s in SEEDS] for m in METHODS}
ref = json.loads((RES / "reference.json").read_text(encoding="utf-8"))
x = [row[0] for row in ref["x"]]
u = [row[0] for row in ref["u"]]
u_rms = math.sqrt(st.mean(v * v for v in u))

SUP = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")


def decade(k):
    return {0: "1", 1: "10"}.get(k, "10" + str(k).translate(SUP))


def sci(v):
    if v == 0:
        return "0"
    e = math.floor(math.log10(abs(v)))
    return f"{v / 10 ** e:.2f}×10{str(e).translate(SUP)}" if not -2 <= e <= 2 else f"{v:.3g}"


# --- summary table ---
def adam_end(r):
    return [h for h in r["history"] if h["phase"] == "adam"][-1]


def first_below(r, level):
    return next((h["step"] for h in r["history"] if h["step"] > 0 and h["error"][0] < level), None)


summary = {}
for m in METHODS:
    fin = [r["error"][0] for r in runs[m]]
    ad = [adam_end(r)["error"][0] for r in runs[m]]
    summary[m] = {"mean": st.mean(fin), "std": st.stdev(fin), "min": min(fin), "max": max(fin),
                  "adam_median": st.median(ad), "adam_min": min(ad), "adam_max": max(ad),
                  "status": sorted({r["status"] for r in runs[m]}),
                  # PROBLEM: the end of the plateau, and the boundary error against the error size
                  "final": fin, "adam": ad,
                  "first_below_0.1": [first_below(r, 0.1) for r in runs[m]],
                  "first_below_0.01": [first_below(r, 0.01) for r in runs[m]],
                  "plateau_mse": [r["history"][10]["mse"] for r in runs[m]],  # step 1,000
                  "bc_rms_final": [math.sqrt(r["history"][-1]["mse"]["bc"]) for r in runs[m]],
                  "abs_rms_final": [e * u_rms for e in fin],  # relative L2 times the RMS of u
                  "weights_end": [r["history"][-1]["weights"] for r in runs[m]]}
    s = summary[m]
    print(f"{LABEL[m]:5s} final {s['mean']:.2e} ± {s['std']:.1e} (min {s['min']:.2e}, max {s['max']:.2e})"
          f"   after Adam: median {s['adam_median']:.2e} (min {s['adam_min']:.2e}, max {s['adam_max']:.2e})"
          f"  L-BFGS gain (medians) {s['adam_median'] / st.median(fin):.0f}x  {s['status']}")
    for i in SEEDS:
        print(f"      s{i} final {s['final'][i]:.2e}  adam {s['adam'][i]:.2e}  error<0.1 at {s['first_below_0.1'][i]}"
              f"  <0.01 at {s['first_below_0.01'][i]}  step-1000 mse {s['plateau_mse'][i]}"
              f"  final bc rms {s['bc_rms_final'][i]:.1e} vs error rms {s['abs_rms_final'][i]:.1e}  w {s['weights_end'][i]}")
print(f"RMS of u on the test grid: {u_rms:.4f}")
print(f"relative L2 error of the outer solution u = x: "
      f"{math.sqrt(sum((a - b) ** 2 for a, b in zip(x, u))) / math.sqrt(sum(b * b for b in u)):.3f}")
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


# --- 1. the exact solution and the outer solution ---  PROBLEM: boundary layer
X_TICKS = [(v, f"{v:g}") for v in (0, 0.2, 0.4, 0.6, 0.8, 1)]
html = line_chart(
    x, {"outer": x, "exact": u}, logy=False, y_range=(0, 1),
    x_ticks=X_TICKS, y_ticks=[(v, f"{v:g}") for v in (0, 0.25, 0.5, 0.75, 1)], x_title="x",
    y_title="Exact solution u(x) and the line u = x",
    aria="The exact solution on [0, 1], dashed: it follows the line u = x from 0 almost to 1, and then drops "
         "from about 0.99 to 0 in a thin layer of width about 0.01 next to x = 1.",
    caption="The exact solution (dashed) and the line u = x. Away from x = 1 the two are the same. In the last "
            "few hundredths, the solution drops to 0. The predictions of all three methods (seed 0) lie on the "
            "dashed line (the largest gap is below 4×10⁻⁵), so they are not drawn. Hover or use the arrow keys to "
            "read the values.",
    heads=[f"x = {v:.3f}" for v in x], fmt=lambda v: f"{v:.4f}", table_every=50,
)
(HERE / "fig-solution.html").write_text(html, encoding="utf-8")
gap = {m: max(abs(p - e) for p, e in zip(runs[m][0]["prediction"], u)) for m in METHODS}
print("largest gap, seed 0:", {m: f"{v:.2e}" for m, v in gap.items()})
assert max(gap.values()) < 4e-5

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
lo = 10.0 ** math.floor(math.log10(min(min(v) for v in med.values())))
hi = 10.0 ** math.ceil(math.log10(max(max(v) for v in med.values())))
html = line_chart(
    steps, med, logy=True, y_range=(lo, hi),
    x_ticks=[(s, f"{s // 1000}k" if s else "0") for s in range(0, steps[-1] + 1, 4000)],
    y_ticks=log_ticks(lo, hi), x_title="Training step (Adam, then L-BFGS iterations)",
    y_title="Relative L2 error, median of 5 seeds",
    vlines=[(adam_steps, "L-BFGS")],
    aria="Median relative L2 error of the PINN for the boundary layer problem during training, on a log scale, "
         "for NTK, RBA and BRDR weighting. All three stay near 0.9 for the first few thousand steps and then drop "
         "by about two orders of magnitude, NTK first. BRDR jumps back up to about 0.3 near step 8,500. From step "
         "9,600 to the end of Adam, BRDR is lowest, near 2e-4 to 1e-3, and NTK and RBA are near 3e-3 to 1e-2. In the L-BFGS phase all three drop, to about 1.5e-5 (NTK), 7e-5 (RBA) "
         "and 3e-5 (BRDR).",
    caption="Median over 5 seeds of the relative L2 error on 1,001 test points, logged every 100 steps. In the "
            "Adam part, the line is a rolling median over about 1,000 steps: the raw error jumps up and down by up "
            "to two orders of magnitude from one log to the next. The dashed line marks the switch to L-BFGS. "
            "Hover over the chart, or focus it and use the arrow keys, to read the values.",
    heads=[f"step {s:,}" for s in steps], fmt=sci, table_every=10,
)
(HERE / "fig-training.html").write_text(html, encoding="utf-8")

# --- 3. pointwise error of seed 0 over x ---
err = {m: [abs(p - e) for p, e in zip(runs[m][0]["prediction"], u)] for m in METHODS}
for m in METHODS:  # PROBLEM: where the error is, inside or outside the layer
    e2 = [v * v for v in err[m]]
    i_max = max(range(len(x)), key=lambda i: err[m][i])
    print(f"{LABEL[m]} seed 0: max error {err[m][i_max]:.2e} at x = {x[i_max]:.3f}, share of squared error in "
          f"x >= 0.95: {sum(v for xi, v in zip(x, e2) if xi >= 0.95) / sum(e2):.3f}, "
          f"u(0) = {runs[m][0]['prediction'][0]:.1e}, u(1) = {runs[m][0]['prediction'][-1]:.1e}, "
          f"RMS error x < 0.9: {math.sqrt(st.mean(v for xi, v in zip(x, e2) if xi < 0.9)):.1e}")
all_err = sorted(e for v in err.values() for e in v)
lo = 10.0 ** math.floor(math.log10(max(1e-12, all_err[len(all_err) // 20])))  # 5th percentile; deeper dips are cut
hi = 10.0 ** math.ceil(math.log10(all_err[-1]))
html = line_chart(
    x, err, logy=True, y_range=(lo, hi), x_ticks=X_TICKS,
    y_ticks=log_ticks(lo, hi), x_title="x", y_title="Pointwise error |uθ(x) − u(x)|, seed 0",
    aria="Pointwise error of the trained PINN against the exact solution over x from 0 to 1, on a log scale, "
         "for NTK, RBA and BRDR weighting, seed 0. The errors are between about 1e-6 and 3e-5 across the whole "
         "interval, as a wave that changes sign many times. They are not larger inside the layer near x = 1.",
    caption="Error of the final network (seed 0) at the 1,001 test points. The sharp dips are places where the "
            "error changes sign; the deepest ones are cut at the bottom of the chart. Hover or use the arrow keys "
            "to read the values.",
    heads=[f"x = {v:.3f}" for v in x], fmt=sci, table_every=50,
)
(HERE / "fig-error.html").write_text(html, encoding="utf-8")

# --- 4. point weights of RBA and BRDR, seed 0 ---  PROBLEM: new chart
W_METHODS = ["rba", "brdr"]
pw = {}
for m in W_METHODS:
    p = runs[m][0]["point_weights"]["pde"]
    pairs = sorted(zip((q[0] for q in p["x"]), p["w"]))
    mean_w = st.mean(w for _, w in pairs)
    pw[m] = (pairs, mean_w)
xs_w = [a for a, _ in pw["rba"][0]]
assert xs_w == [a for a, _ in pw["brdr"][0]], "seed 0 uses the same points for every method"
wser = {m: [w / pw[m][1] for _, w in pw[m][0]] for m in W_METHODS}
for m in W_METHODS:
    for a, b in [(0, 0.5), (0.5, 0.9), (0.9, 0.95), (0.95, 0.98), (0.98, 0.99), (0.99, 1.0)]:
        ws = [w for xi, w in zip(xs_w, wser[m]) if a <= xi < b]
        print(f"{LABEL[m]} weight / mean, x in [{a}, {b}): {len(ws)} points, mean {st.mean(ws):.2f}, max {max(ws):.2f}")
    bcw = runs[m][0]["point_weights"]["bc"]["w"]
    print(f"{LABEL[m]} boundary weights / mean equation weight: {[round(w / pw[m][1], 3) for w in bcw]}")
print("equation points with x > 0.97 (seed 0):", sum(1 for v in xs_w if v > 0.97))
hi_w = math.ceil(max(max(v) for v in wser.values()))
html = line_chart(
    xs_w, wser, logy=False, y_range=(0, hi_w), x_ticks=X_TICKS,
    y_ticks=[(v, f"{v:g}") for v in range(0, hi_w + 1, 2)], x_title="x",
    y_title="Weight of each equation point ÷ mean weight, seed 0",
    aria="Point weights of the equation residual at the end of training, seed 0, divided by their mean. RBA "
         "weights are near 0.9 up to x = 0.9 and rise inside the layer to about 4 on average for x above 0.99. "
         "BRDR weights are low, about 0.4, on the left half, and highest, about 6 on average, for x between "
         "0.95 and 0.98.",
    caption="The weight of each of the 1,024 equation points after training (seed 0), divided by the mean weight, "
            "with the points sorted by x. NTK has one weight for the whole term, so it is not shown. Hover or use "
            "the arrow keys to read the values.",
    heads=[f"x = {v:.4f}" for v in xs_w], fmt=lambda v: f"{v:.2f}", table_every=32,
)
(HERE / "fig-weights.html").write_text(html, encoding="utf-8")

print("wrote fig-solution.html, fig-training.html, fig-error.html, fig-weights.html")
