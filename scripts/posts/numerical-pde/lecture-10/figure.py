"""Builds the FTCS growth chart for the post from experiment.json (part 3).

Writes one HTML block (legend, SVG, hover script, data table) with no blank
lines to figure.html, so it can sit directly inside the Markdown post.
Same layout as scripts/posts/numerical-analysis/figure.py.

Run: python scripts/posts/numerical-pde/lecture-10/figure.py   (writes figure.html next to this file)
"""

import json
import math
from html import escape
from pathlib import Path

data = json.loads(Path(__file__).with_name("experiment.json").read_text(encoding="utf-8"))
p3 = data["part3"]
EVERY = 2
ts = list(range(0, p3["steps"] + 1, EVERY))  # time step number n
ks = list(range(len(ts)))
T_MAX = float(p3["steps"])

# Categorical slots 1-3 of the validated reference palette, in fixed order
COLORS = {
    "r = 0.45": "#2a78d6",
    "r = 0.50": "#eb6834",
    "r = 0.55": "#1baf7a",
}
NAMES = list(COLORS)
SERIES = {name: p3["runs"][name[4:]]["sizes"][::EVERY] for name in NAMES}

W, H = 680, 400
LEFT, RIGHT, TOP, BOTTOM = 64, 520, 30, 340
Y_TOP, Y_BOTTOM = math.log10(3e10), math.log10(0.03)  # log10 of max |U|


def px(t):
    return LEFT + (RIGHT - LEFT) * t / T_MAX


def py(amp):
    value = min(max(math.log10(amp), Y_BOTTOM), Y_TOP)
    return TOP + (BOTTOM - TOP) * (Y_TOP - value) / (Y_TOP - Y_BOTTOM)


def fmt(value):
    return f"{value:.4f}" if value < 10 else f"{value:.3g}".replace("e+0", "e").replace("e+", "e")


svg = []
svg.append(
    f'<svg viewBox="0 0 {W} {H}" role="img" tabindex="0" '
    f'aria-label="Largest value of the FTCS solution of the heat equation over 200 time steps, on a log scale, for r = k/h squared equal to 0.45, 0.50 and 0.55. The first two decay slowly to about 0.2. The third grows to about 8 billion." '
    f'class="viz-plot">'
)
for value, label in ((1e10, "10¹⁰"), (1e8, "10⁸"), (1e6, "10⁶"), (1e4, "10⁴"), (1e2, "100"), (1, "1"), (0.1, "0.1")):
    y = py(value)
    svg.append(f'<line x1="{LEFT}" x2="{RIGHT}" y1="{y:.1f}" y2="{y:.1f}" class="viz-grid"/>')
    svg.append(f'<text x="{LEFT - 10}" y="{y + 4:.1f}" text-anchor="end" class="viz-tick">{label}</text>')
for t in (0, 50, 100, 150, 200):
    svg.append(f'<text x="{px(t):.1f}" y="{BOTTOM + 22}" text-anchor="middle" class="viz-tick">{t}</text>')
svg.append(f'<text x="{(LEFT + RIGHT) / 2:.1f}" y="{BOTTOM + 48}" text-anchor="middle" class="viz-axis-title">Time step n</text>')
svg.append(f'<text x="{LEFT}" y="{TOP - 16}" class="viz-axis-title">Largest value max |U| on the grid</text>')

for name in NAMES:
    errs = SERIES[name]
    points = [(px(t), py(e)) for t, e in zip(ts, errs)]
    path = " ".join(f"{'M' if i == 0 else 'L'}{x:.1f},{y:.1f}" for i, (x, y) in enumerate(points))
    color = COLORS[name]
    svg.append(f'<path d="{path}" fill="none" stroke="{color}" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/>')

# r = 0.45 and r = 0.50 end close together, so move their labels apart
LABEL_DY = {"r = 0.45": -6, "r = 0.50": 14, "r = 0.55": 4}
# Direct labels at the right end, where the lines are spread out
for name, dy in LABEL_DY.items():
    errs = SERIES[name]
    svg.append(f'<text x="{RIGHT + 10}" y="{py(errs[-1]) + dy:.1f}" class="viz-label">{escape(name)}</text>')

svg.append(f'<line x1="0" x2="0" y1="{TOP}" y2="{BOTTOM}" class="viz-crosshair" visibility="hidden"/>')
svg.append("</svg>")

rows = [[fmt(SERIES[name][i]) for name in NAMES] for i in range(len(ks))]
heads = [f"n = {t}" for t in ts]
# the table shows every 25th step
turn_idx = [i for i in ks if ts[i] % 25 == 0]

legend = "".join(
    f'<span class="viz-key"><span class="viz-swatch" style="background:{COLORS[name]}"></span>{escape(name)}</span>'
    for name in NAMES
)
table_head = "<tr><th>n</th>" + "".join(f"<th>{escape(name)}</th>" for name in NAMES) + "</tr>"
table_body = "".join(
    f"<tr><td>{ts[i]}</td>" + "".join(f"<td>{v}</td>" for v in rows[i]) + "</tr>" for i in turn_idx
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
    "<figcaption>FTCS for u_t = u_xx with h = 1/25 and a hat-shaped start of height 1, on a log scale. "
    "The lines for r = 0.45 and r = 0.50 almost lie on top of each other. The table shows every 25th step. "
    "Hover over the chart, or focus it and use the arrow keys, to read the values.</figcaption>",
    "<details class=\"viz-table\"><summary>Show the numbers</summary>"
    f"<table><thead>{table_head}</thead><tbody>{table_body}</tbody></table></details>",
    script,
    "</figure>",
]
Path(__file__).with_name("figure.html").write_text("\n".join(parts) + "\n", encoding="utf-8")
