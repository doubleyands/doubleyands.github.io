"""Numbers and charts for the Problem 7 post, from results/*.json (standard library only).

Prints the summary table and writes three HTML blocks (legend, SVG, hover script, data
table, no blank lines) that sit directly inside the Markdown post:
  fig-solution.html  the exact solution u(x, t) at three times
  fig-training.html  relative L2 error during training (median of 5 seeds) per weighting
  fig-error.html     RMS over x of the error u - exact at each time t (seed 0)
Copied from 06-advection/figure.py; the parts marked PROBLEM are changed.

Run: python scripts/posts/pinn/07-wave/figure.py
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
# PROBLEM: three gray series for the solution chart, one per time (darker = earlier)
TIMES = [0.0, 0.25, 0.5]
for k, (t, c) in enumerate(zip(TIMES, ["#1d1d1f", "#6e6e73", "#aeaeb2"])):
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
CHECK = [0, 1000, 5000]  # PROBLEM: the start, early, and the middle of Adam
TERMS = ("pde", "ic", "ic_t", "bc")
for m in METHODS:
    fin = [r["error"][0] for r in runs[m]]
    ad = [adam_end(r)["error"][0] for r in runs[m]]
    summary[m] = {"mean": st.mean(fin), "std": st.stdev(fin), "min": min(fin), "max": max(fin),
                  "adam_median": st.median(ad), "adam_min": min(ad), "adam_max": max(ad),
                  "status": sorted({r["status"] for r in runs[m]}),
                  # PROBLEM: when each run gets below some errors, the residual of each term, the weights
                  "final": fin, "adam": ad,
                  "first_below_0.5": [first_below(r, 0.5) for r in runs[m]],
                  "first_below_0.1": [first_below(r, 0.1) for r in runs[m]],
                  "first_below_0.01": [first_below(r, 0.01) for r in runs[m]],
                  "mse_at": {s: [{"error": h["error"][0], **h["mse"]} for r in runs[m]
                                 for h in r["history"] if h["step"] == s] for s in CHECK},
                  "rms_residual_final": [{k: math.sqrt(v) for k, v in r["history"][-1]["mse"].items()} for r in runs[m]],
                  "weight_ratio": {k: [weight_ratio(r, k) for r in runs[m]] for k in TERMS[1:]}}
    s = summary[m]
    print(f"{LABEL[m]:5s} final {s['mean']:.2e} ± {s['std']:.1e} (min {s['min']:.2e}, max {s['max']:.2e})"
          f"   after Adam: median {s['adam_median']:.2e} (min {s['adam_min']:.2e}, max {s['adam_max']:.2e})"
          f"  L-BFGS gain (medians) {s['adam_median'] / st.median(fin):.1f}x  {s['status']}")
    for i in SEEDS:
        rr = s["rms_residual_final"][i]
        print(f"      s{i} final {s['final'][i]:.2e}  adam {s['adam'][i]:.2e}  error<0.5 at {s['first_below_0.5'][i]}"
              f"  <0.1 at {s['first_below_0.1'][i]}  <0.01 at {s['first_below_0.01'][i]}  residual rms "
              + "  ".join(f"{k} {rr[k]:.1e}" for k in TERMS)
              + "  weight/pde " + "  ".join(f"{k} {s['weight_ratio'][k][i]:.3g}" for k in TERMS[1:]))
    for k in TERMS:
        v = [rr[k] for rr in s["rms_residual_final"]]
        print(f"      final residual rms {k}: {min(v):.1e}-{max(v):.1e}")
    for k in TERMS[1:]:
        v = s["weight_ratio"][k]
        print(f"      final weight {k}/pde: {min(v):.3g}-{max(v):.3g}")
    for step in CHECK:
        q = s["mse_at"][step]
        print(f"      step {step}: error {min(a['error'] for a in q):.2f}-{max(a['error'] for a in q):.2f}  mse " +
              "  ".join(f"{k} {min(a[k] for a in q):.1e}-{max(a[k] for a in q):.1e}" for k in TERMS))
print(f"RMS of u on the test grid: {u_rms:.4f}")
# PROBLEM: the relative L2 error of a network that has only the slow part, or only the fast part
fast = [0.5 * math.sin(4 * math.pi * x) * math.cos(8 * math.pi * t) for x, t in X]
print(f"relative L2 error if the fast part is missing: {rms(fast) / u_rms:.3f}, "
      f"if the slow part is missing: {rms(a - b for a, b in zip(u, fast)) / u_rms:.3f}")
print(f"size of u_tt: fast part 0.5 (8 pi)^2 = {0.5 * (8 * math.pi) ** 2:.0f}, slow part (2 pi)^2 = {(2 * math.pi) ** 2:.1f}")
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




ARIA_TRAINING = (  # PROBLEM: checked against the printed numbers
    "Median relative L2 error of the PINN for the 1D wave equation during training, on a log scale, for NTK, RBA "
    "and BRDR weighting. NTK falls from about 0.45 to about 0.1 at step 5,000 and 1.3e-2 at the end of Adam. RBA "
    "and BRDR fall to about 0.4 by step 4,000, stay near 0.3 to 0.4 until step 10,000, and reach about 0.19 and "
    "0.16 at the end of Adam. After L-BFGS, NTK is about 1.7e-3, RBA about 5.4e-2 and BRDR about 4.3e-2.")
ARIA_ERROR = (
    "Error of the trained PINN at each time t from 0 to 1, as the RMS over the 101 x values, on a log scale, for "
    "NTK, RBA and BRDR weighting, seed 0. NTK grows from about 1.5e-4 at t = 0 to about 1e-3 at t = 0.3, then goes up and down and ends at 1.9e-3 at t = 1. RBA and BRDR "
    "start near 1e-2 at t = 0 and grow to about 4.5e-2 at t = 1, with ups and downs in between.")


# --- 1. the exact solution at three times ---  PROBLEM: wave, x in [0, 1]
X_TICKS = [(v, f"{v:g}") for v in (0, 0.2, 0.4, 0.6, 0.8, 1)]
t_index = [min(range(nt), key=lambda j: abs(ts_grid[j] - t)) for t in TIMES]
sol = {f"t{k}": [at(u, i, j) for i in range(nx)] for k, j in enumerate(t_index)}
gap = {m: max(abs(p - e) for p, e in zip(runs[m][0]["prediction"], u)) for m in METHODS}
where = {m: X[max(range(len(u)), key=lambda i: abs(runs[m][0]["prediction"][i] - u[i]))] for m in METHODS}
print("largest gap, seed 0:", {m: f"{v:.2e} at (x, t) = {where[m]}" for m, v in gap.items()})
assert 9e-3 <= gap["ntk"] < 9.6e-3 and 9e-2 <= min(gap["rba"], gap["brdr"]) and max(gap["rba"], gap["brdr"]) < 1.03e-1
html = line_chart(
    xs_grid, sol, logy=False, y_range=(-1.6, 1.6),
    x_ticks=X_TICKS, y_ticks=[(v, f"{v:g}") for v in (-1.5, -1, -0.5, 0, 0.5, 1, 1.5)], x_title="x",
    y_title="Exact solution u(x, t) at three times",
    aria="The exact solution of the wave equation at t = 0, 0.25 and 0.5. At t = 0 it is one slow arch plus four "
         "small fast bumps of height 0.5, with a highest value near 1.4. At t = 0.25 the slow arch is gone and only "
         "the fast part is left. At t = 0.5 the slow arch is upside down and the fast part is back as at t = 0.",
    caption="The exact solution at three times (darker is earlier). The slow part sin(πx) swings with period 1 in t, "
            "and the fast part 0.5 sin(4πx) swings with period 0.25. At t = 0.25 the slow part is zero and only the "
            "fast part is left. The predictions (seed 0) are not drawn: on the whole test grid, the largest gap to "
            "the exact solution is about 9×10⁻³ for NTK and about 0.09 to 0.10 for RBA and BRDR. Hover or use the "
            "arrow keys to read the values.",
    heads=[f"x = {v:.2f}" for v in xs_grid], fmt=lambda v: f"{v:.4f}", table_every=10,
)
(HERE / "fig-solution.html").write_text(html, encoding="utf-8")
print("max of u at t = 0:", f"{max(sol['t0']):.4f}")
print("times drawn:", [ts_grid[j] for j in t_index])

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
for k in (0, 2, 5, 10, 20, 30, 40, 50, 60, 80, 100, 150, 190, n_adam - 1, n_adam, len(steps) - 1):
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
    caption="Median over 5 seeds of the relative L2 error on the 101 × 101 test grid, logged every 100 steps. In "
            "the Adam part, the line is a rolling median over about 1,000 steps, because the raw error jumps up "
            "and down from one log to the next. The dashed line marks the switch to L-BFGS. Hover over the chart, "
            "or focus it and use the arrow keys, to read the values.",
    heads=[f"step {s:,}" for s in steps], fmt=sci, table_every=10,
)
(HERE / "fig-training.html").write_text(html, encoding="utf-8")

# --- 3. error of seed 0 over time: RMS over x at each t ---
err_t = {m: [rms(at(runs[m][0]["prediction"], i, j) - at(u, i, j) for i in range(nx)) for j in range(nt)]
         for m in METHODS}
for m in METHODS:
    e = err_t[m]
    print(LABEL[m], "seed 0, RMS error every 0.05 in t:", [f"{e[j]:.1e}" for j in range(0, nt, 5)])
    print(f"{LABEL[m]} seed 0: RMS error at t = 0 {e[0]:.1e}, 0.1 {e[10]:.1e}, 0.25 {e[25]:.1e}, 0.5 {e[50]:.1e}, "
          f"0.75 {e[75]:.1e}, 1 {e[-1]:.1e}; max {max(e):.2e} at t = {ts_grid[e.index(max(e))]:.2f}, "
          f"min {min(e):.1e} at t = {ts_grid[e.index(min(e))]:.2f}; mean over t < 0.2 {st.mean(e[:20]):.1e}, "
          f"t >= 0.8 {st.mean(e[80:]):.1e}")
lo = 10.0 ** math.floor(math.log10(min(min(v) for v in err_t.values())))
hi = 10.0 ** math.ceil(math.log10(max(max(v) for v in err_t.values())))
T_TICKS = [(v, f"{v:g}") for v in (0, 0.2, 0.4, 0.6, 0.8, 1)]
html = line_chart(
    ts_grid, err_t, logy=True, y_range=(lo, hi), x_ticks=T_TICKS,
    y_ticks=log_ticks(lo, hi), x_title="t", y_title="RMS over x of the error uθ(x, t) − u(x, t), seed 0",
    aria=ARIA_ERROR,
    caption="The error of the final network (seed 0) at each of the 101 test times, as the RMS over the 101 test "
            "points in x. The first value (t = 0) is the error in the initial condition. Hover or use the arrow "
            "keys to read the values.",
    heads=[f"t = {v:.2f}" for v in ts_grid], fmt=sci, table_every=10,
)
(HERE / "fig-error.html").write_text(html, encoding="utf-8")

# --- point weights of RBA and BRDR, seed 0: only printed (no clear pattern, see the post) ---
for m in ["rba", "brdr"]:
    p = runs[m][0]["point_weights"]["pde"]
    mean_w = st.mean(p["w"])
    for name, k, edges in [("t", 1, [0, 0.1, 0.25, 0.5, 0.75, 1.01]), ("x", 0, [0, 0.2, 0.4, 0.6, 0.8, 1.01])]:
        means = [st.mean(w / mean_w for q, w in zip(p["x"], p["w"]) if a <= q[k] < b) for a, b in zip(edges, edges[1:])]
        print(f"{LABEL[m]} seed 0 equation weight / mean, 5 bins in {name}:", [f"{v:.2f}" for v in means])

print("wrote fig-solution.html, fig-training.html, fig-error.html")
