"""Numbers and charts for the Problem 1 post, from results/*.json (standard library only).

Prints the summary table and writes three HTML blocks (legend, SVG, hover script, data
table, no blank lines) that sit directly inside the Markdown post:
  fig-training.html  relative L2 error during training (median of 5 seeds) per weighting
  fig-error.html     pointwise error |u - exact| over t (seed 0)
  fig-weights.html   RBA and BRDR point weights over t (seed 0), each divided by its mean
Later posts copy this file and change the parts marked PROBLEM.

Run: python scripts/posts/pinn/01-decay/figure.py
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
SEEDS = range(5)

runs = {m: [json.loads((RES / f"{m}-s{s}.json").read_text(encoding="utf-8")) for s in SEEDS] for m in METHODS}
ref = json.loads((RES / "reference.json").read_text(encoding="utf-8"))

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
    return [h for h in r["history"] if h["phase"] == "adam"][-1]["error"][0]


summary = {}
for m in METHODS:
    fin = [r["error"][0] for r in runs[m]]
    ad = [adam_end(r) for r in runs[m]]
    summary[m] = {"mean": st.mean(fin), "std": st.stdev(fin), "min": min(fin), "max": max(fin),
                  "adam_median": st.median(ad), "adam_min": min(ad), "adam_max": max(ad),
                  "status": sorted({r["status"] for r in runs[m]})}
    s = summary[m]
    print(f"{LABEL[m]:5s} final {s['mean']:.2e} ± {s['std']:.1e} (min {s['min']:.2e}, max {s['max']:.2e})"
          f"   after Adam: median {s['adam_median']:.2e} (min {s['adam_min']:.2e}, max {s['adam_max']:.2e})  {s['status']}")
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
                   f'stroke-linejoin="round" stroke-linecap="round"/>')
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


# --- 1. training curves: median over seeds at each logged step ---
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
    aria="Median relative L2 error of the PINN for u' = -2u during training, on a log scale, for NTK, RBA and BRDR "
         "weighting. All three fall from about 1 to between 1e-4 and 1e-2 during Adam and wobble, then drop to "
         "about 3e-5 to 5e-5 in the L-BFGS phase.",
    caption="Median over 5 seeds of the relative L2 error on 301 test points, logged every 100 steps. In the Adam "
            "part, the line is a rolling median over about 1,000 steps: the raw error jumps up and down by up to "
            "two orders of magnitude from one log to the next. The dashed line marks the switch to L-BFGS. Hover "
            "over the chart, or focus it and use the arrow keys, to read the values.",
    heads=[f"step {s:,}" for s in steps], fmt=sci, table_every=10,
)
(HERE / "fig-training.html").write_text(html, encoding="utf-8")

# --- 2. pointwise error of seed 0 over t ---  PROBLEM: 1D grid, one output
t = [row[0] for row in ref["x"]]
u = [row[0] for row in ref["u"]]
err = {m: [abs(p - e) for p, e in zip(runs[m][0]["prediction"], u)] for m in METHODS}
all_err = sorted(e for v in err.values() for e in v)
lo = 10.0 ** math.floor(math.log10(max(1e-12, all_err[len(all_err) // 20])))  # 5th percentile; deeper dips are cut
hi = 10.0 ** math.ceil(math.log10(all_err[-1]))
html = line_chart(
    t, err, logy=True, y_range=(lo, hi), x_ticks=[(v, f"{v:g}") for v in (0, 0.5, 1, 1.5, 2, 2.5, 3)],
    y_ticks=log_ticks(lo, hi), x_title="t", y_title="Pointwise error |u(t) − exp(−2t)|, seed 0",
    aria="Pointwise error of the trained PINN against the exact solution exp(-2t) over t from 0 to 3, on a log "
         "scale, for NTK, RBA and BRDR weighting, seed 0. The errors stay between about 1e-7 and 1e-5.",
    caption="Error of the final network (seed 0) at the 301 test points. The sharp dips are places where the "
            "error changes sign; the deepest ones are cut at the bottom of the chart. Hover or use the arrow keys "
            "to read the values.",
    heads=[f"t = {v:.2f}" for v in t], fmt=sci, table_every=25,
)
(HERE / "fig-error.html").write_text(html, encoding="utf-8")

# --- 3. point weights of seed 0 over t, each divided by its mean ---  PROBLEM: the pde term
pw = {m: runs[m][0]["point_weights"]["pde"] for m in ("rba", "brdr")}
order = sorted(range(len(pw["rba"]["x"])), key=lambda i: pw["rba"]["x"][i][0])
assert pw["rba"]["x"] == pw["brdr"]["x"], "same seed, same points"
tw = [pw["rba"]["x"][i][0] for i in order]
wn = {m: [pw[m]["w"][i] / st.mean(pw[m]["w"]) for i in order] for m in ("rba", "brdr")}
top = math.ceil(max(max(v) for v in wn.values()))
html = line_chart(
    tw, wn, logy=False, y_range=(0, top), x_ticks=[(v, f"{v:g}") for v in (0, 0.5, 1, 1.5, 2, 2.5, 3)],
    y_ticks=[(v, f"{v:g}") for v in range(0, top + 1, max(1, top // 5))], x_title="t",
    y_title="Weight of each residual point ÷ mean weight, seed 0",
    aria="Point weights of RBA and BRDR at the end of Adam over t, each divided by its mean, seed 0.",
    caption="Weights of the 512 residual points at the end of Adam (they stay fixed during L-BFGS), sorted by t. "
            "1 means an average point. Hover or use the arrow keys to read the values.",
    heads=[f"t = {v:.3f}" for v in tw], fmt=lambda v: f"{v:.2f}", table_every=32,
)
(HERE / "fig-weights.html").write_text(html, encoding="utf-8")
print("wrote fig-training.html, fig-error.html, fig-weights.html")
