"""Numbers and charts for the Problem 8 post, from results/*.json (standard library only).

Prints the summary table and writes four HTML blocks (legend, SVG, hover script, data
table, no blank lines) that sit directly inside the Markdown post:
  fig-solution.html  the exact solution u(x, t) at three times
  fig-training.html  relative L2 error during training (median of 5 seeds) per weighting
  fig-error.html     RMS over x of the error u - exact at each time t (seed 0), with a gray e^t line
  fig-front.html     the signed error u - exact at the last time t = 5 over x (seed 0)
Copied from 07-wave/figure.py; the parts marked PROBLEM are changed.

Run: python scripts/posts/pinn/08-fisher-kpp/figure.py
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
# PROBLEM: three gray series for the solution chart, one per time (darker = earlier), and a gray e^t line
TIMES = [0.0, 2.5, 5.0]
for k, (t, c) in enumerate(zip(TIMES, ["#1d1d1f", "#6e6e73", "#aeaeb2"])):
    LABEL[f"t{k}"] = f"t = {t:g}"
    COLOR[f"t{k}"] = c
LABEL["exp"], COLOR["exp"] = "∝ eᵗ", "#8e8e93"
DASH = {"exp": "6 5"}
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
TERMS = ("pde", "ic", "bc")
for m in METHODS:
    fin = [r["error"][0] for r in runs[m]]
    ad = [adam_end(r)["error"][0] for r in runs[m]]
    summary[m] = {"mean": st.mean(fin), "std": st.stdev(fin), "min": min(fin), "max": max(fin),
                  "adam_median": st.median(ad), "adam_min": min(ad), "adam_max": max(ad),
                  "status": sorted({r["status"] for r in runs[m]}),
                  # PROBLEM: when each run gets below some errors, the residual of each term, the weights
                  "final": fin, "adam": ad,
                  "first_below_0.1": [first_below(r, 0.1) for r in runs[m]],
                  "first_below_0.03": [first_below(r, 0.03) for r in runs[m]],
                  "first_below_0.01": [first_below(r, 0.01) for r in runs[m]],
                  "mse_at": {s: [{"error": h["error"][0], **h["mse"]} for r in runs[m]
                                 for h in r["history"] if h["step"] == s] for s in CHECK},
                  "rms_residual_final": [{k: math.sqrt(v) for k, v in r["history"][-1]["mse"].items()} for r in runs[m]],
                  "rms_residual_adam": [{k: math.sqrt(v) for k, v in adam_end(r)["mse"].items()} for r in runs[m]],
                  "weight_ratio": {k: [weight_ratio(r, k) for r in runs[m]] for k in TERMS[1:]}}
    s = summary[m]
    print(f"{LABEL[m]:5s} final {s['mean']:.2e} ± {s['std']:.1e} (min {s['min']:.2e}, max {s['max']:.2e})"
          f"   after Adam: median {s['adam_median']:.2e} (min {s['adam_min']:.2e}, max {s['adam_max']:.2e})"
          f"  L-BFGS gain (medians) {s['adam_median'] / st.median(fin):.1f}x  {s['status']}")
    for i in SEEDS:
        rr = s["rms_residual_final"][i]
        print(f"      s{i} final {s['final'][i]:.2e}  adam {s['adam'][i]:.2e}  error<0.1 at {s['first_below_0.1'][i]}"
              f"  <0.03 at {s['first_below_0.03'][i]}  <0.01 at {s['first_below_0.01'][i]}  residual rms "
              + "  ".join(f"{k} {rr[k]:.1e}" for k in TERMS)
              + "  weight/pde " + "  ".join(f"{k} {s['weight_ratio'][k][i]:.3g}" for k in TERMS[1:]))
    for key in ("rms_residual_adam", "rms_residual_final"):
        for k in TERMS:
            v = [rr[k] for rr in s[key]]
            print(f"      {key} {k}: {min(v):.1e}-{max(v):.1e}")
    for k in TERMS[1:]:
        v = s["weight_ratio"][k]
        print(f"      final weight {k}/pde: {min(v):.3g}-{max(v):.3g}")
    for step in CHECK:
        q = s["mse_at"][step]
        print(f"      step {step}: error {min(a['error'] for a in q):.3f}-{max(a['error'] for a in q):.3f}  mse " +
              "  ".join(f"{k} {min(a[k] for a in q):.1e}-{max(a[k] for a in q):.1e}" for k in TERMS))
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





ARIA_TRAINING = (  # PROBLEM: checked against the printed numbers
    "Median relative L2 error of the PINN for the Fisher-KPP equation during training, on a log scale, for NTK, "
    "RBA and BRDR weighting. RBA falls fastest, to about 1.8e-2 by step 5,000, and ends Adam near 6e-3. BRDR falls "
    "to about 4.5e-2 by step 5,000 and ends Adam near 2e-2. NTK falls more slowly and stays near 3e-2 to 6e-2 in the "
    "second half of Adam. After L-BFGS, RBA is about 3.9e-3, NTK about 7.6e-3 and BRDR about 1.8e-2.")
ARIA_ERROR = (
    "Error of the trained PINN at each time t from 0 to 5, as the RMS over the 301 x values, on a log scale, for "
    "NTK, RBA and BRDR weighting, seed 0, with a gray dashed line that grows like e to the t. All three errors "
    "grow by about 60 to 80 times from t = 0.5 to t = 5, close to the gray line.")
ARIA_FRONT = (
    "Signed error of the trained PINN at the last time t = 5, over x from -10 to 20, seed 0, for NTK, RBA and BRDR. "
    "Behind the front, for x below about 5, all errors are close to zero. Ahead of the front, between x = 8 and "
    "x = 19, all three are positive, with a peak near x = 14: about 0.03 for NTK, 0.01 for RBA and 0.045 for BRDR. "
    "At x = 20 the error is close to zero again.")


# --- 1. the exact solution at three times ---  PROBLEM: Fisher-KPP, x in [-10, 20], t in [0, 5]
X_TICKS = [(v, f"{v:g}") for v in (-10, -5, 0, 5, 10, 15, 20)]
t_index = [min(range(nt), key=lambda j: abs(ts_grid[j] - t)) for t in TIMES]
sol = {f"t{k}": [at(u, i, j) for i in range(nx)] for k, j in enumerate(t_index)}
gap = {m: max(abs(p - e) for p, e in zip(runs[m][0]["prediction"], u)) for m in METHODS}
where = {m: X[max(range(len(u)), key=lambda i: abs(runs[m][0]["prediction"][i] - u[i]))] for m in METHODS}
print("largest gap, seed 0:", {m: f"{v:.2e} at (x, t) = {where[m]}" for m, v in gap.items()})
assert 3.1e-2 < gap["ntk"] < 3.3e-2 and 1.0e-2 < gap["rba"] < 1.05e-2 and 4.5e-2 < gap["brdr"] < 4.6e-2
for k, j in enumerate(t_index):  # where u = 1/2 at each drawn time, against the exact front x = sqrt(6) ln(sqrt 2 - 1) + 5 t / sqrt 6
    half = next(xs_grid[i] for i in range(nx) if sol[f"t{k}"][i] < 0.5)
    print(f"t = {ts_grid[j]:g}: u drops below 1/2 at x = {half:.1f} "
          f"(exact {math.sqrt(6) * math.log(math.sqrt(2) - 1) + 5 * ts_grid[j] / math.sqrt(6):.2f})")
html = line_chart(
    xs_grid, sol, logy=False, y_range=(0, 1.05),
    x_ticks=X_TICKS, y_ticks=[(v, f"{v:g}") for v in (0, 0.25, 0.5, 0.75, 1)], x_title="x",
    y_title="Exact solution u(x, t) at three times",
    aria="The exact solution of the Fisher-KPP equation at t = 0, 2.5 and 5. Each is a smooth step from 1 on the "
         "left to 0 on the right. The step has the same shape at each time and moves to the right by about 5 in x "
         "every 2.5 in t.",
    caption="The exact solution at three times (darker is earlier). It is a front: u = 1 behind it, u = 0 ahead of "
            "it. The front keeps its shape and moves right at speed 5/√6 ≈ 2.04. The predictions (seed 0) are not "
            "drawn: on the whole test grid, the largest gap to the exact solution is about 0.03 for NTK, 0.01 for "
            "RBA and 0.045 for BRDR. Hover or use the arrow keys to read the values.",
    heads=[f"x = {v:.1f}" for v in xs_grid], fmt=lambda v: f"{v:.4f}", table_every=10,
)
(HERE / "fig-solution.html").write_text(html, encoding="utf-8")
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
for k in (0, 2, 5, 10, 20, 30, 40, 50, 60, 80, 100, 120, 150, 190, n_adam - 1, n_adam, len(steps) - 1):
    print(f"step {steps[k]}: smoothed median", {m: f"{med[m][k]:.2e}" for m in METHODS})
for m in METHODS:  # PROBLEM: how much the raw error of one run jumps in the second half of Adam
    for r in runs[m][:2]:
        e = [h["error"][0] for h in r["history"] if 10_000 <= h["step"] < adam_steps]
        print(f"{LABEL[m]} s{r['seed']} raw error, steps 10k-20k: min {min(e):.1e}, max {max(e):.1e}")
lo = 10.0 ** math.floor(math.log10(min(min(v) for v in med.values())))
hi = 10.0 ** math.ceil(math.log10(max(max(v) for v in med.values())))
html = line_chart(
    steps, med, logy=True, y_range=(lo, hi),
    x_ticks=[(s, f"{s // 1000}k" if s else "0") for s in range(0, steps[-1] + 1, 4000)],
    y_ticks=log_ticks(lo, hi), x_title="Training step (Adam, then L-BFGS iterations)",
    y_title="Relative L2 error, median of 5 seeds",
    vlines=[(adam_steps, "L-BFGS")],
    aria=ARIA_TRAINING,
    caption="Median over 5 seeds of the relative L2 error on the 301 × 101 test grid, logged every 100 steps. In "
            "the Adam part, the line is a rolling median over about 1,000 steps, because the raw error jumps up "
            "and down from one log to the next. The dashed line marks the switch to L-BFGS. Hover over the chart, "
            "or focus it and use the arrow keys, to read the values.",
    heads=[f"step {s:,}" for s in steps], fmt=sci, table_every=10,
)
(HERE / "fig-training.html").write_text(html, encoding="utf-8")

# --- 3. error of seed 0 over time: RMS over x at each t, with a gray e^t line ---
err_t = {m: [rms(at(runs[m][0]["prediction"], i, j) - at(u, i, j) for i in range(nx)) for j in range(nt)]
         for m in METHODS}
j1 = ts_grid.index(1.0)
for m in METHODS:
    e = err_t[m]
    # least-squares slope of ln(error) against t for t in [1, 5]
    tt, ll = ts_grid[j1:], [math.log(v) for v in e[j1:]]
    tm, lm = st.mean(tt), st.mean(ll)
    rate = sum((a - tm) * (b - lm) for a, b in zip(tt, ll)) / sum((a - tm) ** 2 for a in tt)
    print(LABEL[m], "seed 0, RMS error every 0.5 in t:", [f"{e[j]:.1e}" for j in range(0, nt, 10)])
    print(f"{LABEL[m]} seed 0: growth rate of ln(error) for t in [1, 5]: {rate:.2f}; "
          f"error(5) / error(0.5) = {e[-1] / e[10]:.0f}, error(5) / error(0) = {e[-1] / e[0]:.0f}")
err_t["exp"] = [err_t["ntk"][j1] * math.exp(t - 1) for t in ts_grid]  # e^t through the NTK value at t = 1
lo = 10.0 ** math.floor(math.log10(min(min(v) for v in err_t.values())))
hi = 10.0 ** math.ceil(math.log10(max(max(v) for v in err_t.values())))
T_TICKS = [(v, f"{v:g}") for v in (0, 1, 2, 3, 4, 5)]
html = line_chart(
    ts_grid, err_t, logy=True, y_range=(lo, hi), x_ticks=T_TICKS,
    y_ticks=log_ticks(lo, hi), x_title="t", y_title="RMS over x of the error uθ(x, t) − u(x, t), seed 0",
    aria=ARIA_ERROR,
    caption="The error of the final network (seed 0) at each of the 101 test times, as the RMS over the 301 test "
            "points in x. The gray dashed line grows like eᵗ and passes through the NTK value at t = 1; it is "
            "there only to compare the slope. Hover or use the arrow keys to read the values.",
    heads=[f"t = {v:.2f}" for v in ts_grid], fmt=sci, table_every=10,
)
(HERE / "fig-error.html").write_text(html, encoding="utf-8")

# --- 4. signed error at t = 5 over x, seed 0 ---
DASH.clear()
front = {m: [at(runs[m][0]["prediction"], i, nt - 1) - at(u, i, nt - 1) for i in range(nx)] for m in METHODS}
for m in METHODS:
    e = front[m]
    k = max(range(nx), key=lambda i: abs(e[i]))
    pos = sum(1 for i in range(nx) if 8 <= xs_grid[i] <= 19 and e[i] > 0)
    tot = sum(1 for i in range(nx) if 8 <= xs_grid[i] <= 19)
    print(f"{LABEL[m]} seed 0 at t = 5: largest |error| {e[k]:.2e} at x = {xs_grid[k]:.1f} (u = {at(u, k, nt - 1):.3f}); "
          f"positive at {pos}/{tot} x in [8, 19]; largest |error| for x < 5: {max(abs(e[i]) for i in range(nx) if xs_grid[i] < 5):.1e}")
for m in METHODS:  # the share of the squared error where u < 0.1 (ahead of the front), whole grid
    p = runs[m][0]["prediction"]
    tot = sum((a - b) ** 2 for a, b in zip(p, u))
    print(f"{LABEL[m]} seed 0: share of squared error where u < 0.1: {sum((a - b) ** 2 for a, b in zip(p, u) if b < 0.1) / tot:.2f}, "
          f"0.1 <= u <= 0.9: {sum((a - b) ** 2 for a, b in zip(p, u) if 0.1 <= b <= 0.9) / tot:.2f}, "
          f"t >= 4: {sum((a - b) ** 2 for (x, t), a, b in zip(X, p, u) if t >= 4) / tot:.2f}")
html = line_chart(
    xs_grid, front, logy=False, y_range=(-0.01, 0.05),
    x_ticks=X_TICKS, y_ticks=[(v, f"{v:g}") for v in (-0.01, 0, 0.01, 0.02, 0.03, 0.04, 0.05)], x_title="x",
    y_title="Error uθ(x, 5) − u(x, 5) at the last time, seed 0",
    aria=ARIA_FRONT,
    caption="The signed error of the final network (seed 0) at t = 5, over the 301 test points in x. At this time "
            "the exact front (u = 1/2) is near x = 8. The error is positive ahead of the front, so the predicted "
            "front is a little too far ahead. Hover or use the arrow keys to read the values.",
    heads=[f"x = {v:.1f}" for v in xs_grid], fmt=lambda v: f"{v:.2e}", table_every=10,
)
(HERE / "fig-front.html").write_text(html, encoding="utf-8")

# --- point weights of RBA and BRDR, seed 0: only printed ---
for m in ["rba", "brdr"]:
    p = runs[m][0]["point_weights"]["pde"]
    mean_w = st.mean(p["w"])
    for name, f, edges in [("t", lambda q: q[1], [0, 1, 2, 3, 4, 5.01]),
                           ("x", lambda q: q[0], [-10, -4, 2, 8, 14, 20.01]),
                           ("x - front", lambda q: q[0] - (-2.159 + 5 * q[1] / math.sqrt(6)), [-40, -6, -2, 0, 2, 6, 40])]:
        means = [st.mean([w / mean_w for q, w in zip(p["x"], p["w"]) if a <= f(q) < b] or [float("nan")])
                 for a, b in zip(edges, edges[1:])]
        print(f"{LABEL[m]} seed 0 equation weight / mean, bins in {name} {edges}:", [f"{v:.2f}" for v in means])

print("wrote fig-solution.html, fig-training.html, fig-error.html, fig-front.html")
