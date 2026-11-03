"""Builds the amplitude chart for the post from experiment.json (part 1).

Writes one HTML block (legend, SVG, hover script, data table) with no blank
lines to figure.html, so it can sit directly inside the Markdown post.
Same layout as scripts/posts/numerical-analysis/figure.py.

Run: python scripts/posts/numerical-pde/lecture-9/figure.py   (writes figure.html next to this file)
"""

import json
import math
from html import escape
from pathlib import Path

data = json.loads(Path(__file__).with_name("experiment.json").read_text(encoding="utf-8"))
p1 = data["part1"]
ts = p1["t"]  # every 4th time step, dt = 0.1
ks = list(range(len(ts)))
T_MAX = 64.0

# Categorical slots 1-4 of the validated reference palette, in fixed order
COLORS = {
    "Forward Euler": "#2a78d6",
    "Backward Euler": "#eb6834",
    "Trapezoid": "#1baf7a",
    "Leapfrog": "#eda100",
}
NAMES = list(COLORS)
KEYS = {"Forward Euler": "FE", "Backward Euler": "BE", "Trapezoid": "trap", "Leapfrog": "leapfrog"}
SERIES = {name: p1[KEYS[name]] for name in NAMES}

W, H = 680, 400
LEFT, RIGHT, TOP, BOTTOM = 64, 520, 30, 340
Y_TOP, Y_BOTTOM = math.log10(30), math.log10(0.03)  # log10 of the amplitude


def px(t):
    return LEFT + (RIGHT - LEFT) * t / T_MAX


def py(amp):
    value = min(max(math.log10(amp), Y_BOTTOM), Y_TOP)
    return TOP + (BOTTOM - TOP) * (Y_TOP - value) / (Y_TOP - Y_BOTTOM)


def fmt(value):
    return f"{value:.4f}" if value < 10 else f"{value:.2f}"


svg = []
svg.append(
    f'<svg viewBox="0 0 {W} {H}" role="img" tabindex="0" '
    f'aria-label="Amplitude of the numerical solution of the oscillator over ten turns, with time step 0.1, on a log scale. Forward Euler grows to about 23, Backward Euler shrinks to about 0.04, and the trapezoid and leapfrog methods stay at 1." '
    f'class="viz-plot">'
)
for value, label in ((10, "10"), (3, "3"), (1, "1"), (0.3, "0.3"), (0.1, "0.1")):
    y = py(value)
    svg.append(f'<line x1="{LEFT}" x2="{RIGHT}" y1="{y:.1f}" y2="{y:.1f}" class="viz-grid"/>')
    svg.append(f'<text x="{LEFT - 10}" y="{y + 4:.1f}" text-anchor="end" class="viz-tick">{label}</text>')
for t in (0, 20, 40, 60):
    svg.append(f'<text x="{px(t):.1f}" y="{BOTTOM + 22}" text-anchor="middle" class="viz-tick">{t}</text>')
svg.append(f'<text x="{(LEFT + RIGHT) / 2:.1f}" y="{BOTTOM + 48}" text-anchor="middle" class="viz-axis-title">Time t (one turn = 2π ≈ 6.28)</text>')
svg.append(f'<text x="{LEFT}" y="{TOP - 16}" class="viz-axis-title">Amplitude √(x² + y²), exact value 1</text>')

for name in NAMES:
    errs = SERIES[name]
    points = [(px(t), py(e)) for t, e in zip(ts, errs)]
    path = " ".join(f"{'M' if i == 0 else 'L'}{x:.1f},{y:.1f}" for i, (x, y) in enumerate(points))
    color = COLORS[name]
    svg.append(f'<path d="{path}" fill="none" stroke="{color}" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/>')

# Trapezoid and leapfrog both end at 1, so move their labels apart
LABEL_DY = {"Forward Euler": 4, "Backward Euler": 4, "Trapezoid": -6, "Leapfrog": 14}
# Direct labels at the right end, where the lines are spread out
for name, dy in LABEL_DY.items():
    errs = SERIES[name]
    svg.append(f'<text x="{RIGHT + 10}" y="{py(errs[-1]) + dy:.1f}" class="viz-label">{escape(name)}</text>')

svg.append(f'<line x1="0" x2="0" y1="{TOP}" y2="{BOTTOM}" class="viz-crosshair" visibility="hidden"/>')
svg.append("</svg>")

rows = [[fmt(SERIES[name][i]) for name in NAMES] for i in range(len(ks))]
heads = [f"t = {t:.1f}" for t in ts]
# the table shows the nearest stored time to each full turn
turn_idx = [min(ks, key=lambda i: abs(ts[i] - 2 * math.pi * j)) for j in range(0, 11, 2)]

legend = "".join(
    f'<span class="viz-key"><span class="viz-swatch" style="background:{COLORS[name]}"></span>{escape(name)}</span>'
    for name in NAMES
)
table_head = "<tr><th>t</th>" + "".join(f"<th>{escape(name)}</th>" for name in NAMES) + "</tr>"
table_body = "".join(
    f"<tr><td>{ts[i]:.1f}</td>" + "".join(f"<td>{v}</td>" for v in rows[i]) + "</tr>" for i in turn_idx
)

script = (
    "<script>(() => {"
    "const fig = document.currentScript.closest('figure');"
    "const svg = fig.querySelector('svg');"
    "const line = svg.querySelector('.viz-crosshair');"
    "const tip = fig.querySelector('.viz-tooltip');"
    f"const names = {json.dumps(NAMES)}; const colors = {json.dumps([COLORS[n] for n in NAMES])};"
    f"const rows = {json.dumps(rows)}; const heads = {json.dumps(heads)};"
    f"const xs = {json.dumps([round(px(t), 1) for t in ts])}; const steps = {len(ks) - 1}, width = {W};"
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

parts = [
    '<figure class="viz">',
    f'<div class="viz-legend">{legend}</div>',
    '<div class="viz-frame">',
    "".join(svg),
    '<div class="viz-tooltip" hidden></div>',
    "</div>",
    "<figcaption>The oscillator x′ = y, y′ = −x with x(0) = 1, y(0) = 0 and Δt = 0.1, for ten turns, on a log scale. "
    "The trapezoid and leapfrog lines lie on top of each other at 1. Leapfrog takes its second value from the exact solution. The table shows every second turn. "
    "Hover over the chart, or focus it and use the arrow keys, to read the values.</figcaption>",
    "<details class=\"viz-table\"><summary>Show the numbers</summary>"
    f"<table><thead>{table_head}</thead><tbody>{table_body}</tbody></table></details>",
    script,
    "</figure>",
]
Path(__file__).with_name("figure.html").write_text("\n".join(parts) + "\n", encoding="utf-8")
