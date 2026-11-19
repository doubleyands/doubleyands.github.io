"""Builds the chart for the Lecture 14 post from experiment.json (part 1).

Upwind on u_t + u_x = 0 from a cosine bump, at Courant numbers 1, 0.5 and 1.05, shown at t = 0.7
(nu = 1.05 lands at t = 0.7035, the nearest step), with the exact solution as a gray dashed line.
Writes one HTML block with no blank lines to figure.html, so it can sit directly inside the Markdown
post. Same layout as scripts/posts/numerical-analysis/figure.py.

Run: python scripts/posts/numerical-pde/lecture-14/figure.py   (writes figure.html next to this file)
"""

import json
from html import escape
from pathlib import Path

data = json.loads(Path(__file__).with_name("experiment.json").read_text(encoding="utf-8"))
p1 = data["part1"]
xs_data = p1["x"]

COLORS = {
    "ν = 1": "#2a78d6",
    "ν = 0.5": "#eb6834",
    "ν = 1.05": "#1baf7a",
}
KEYS = {"ν = 1": "1.0", "ν = 0.5": "0.5", "ν = 1.05": "1.05"}
NAMES = list(COLORS)
SERIES = {name: p1["upwind_runs"][KEYS[name]]["profile"] for name in NAMES}
EXACT = p1["exact_profile"]
ks = list(range(len(xs_data)))

W, H = 680, 400
LEFT, RIGHT, TOP, BOTTOM = 56, 560, 30, 340
Y_TOP, Y_BOTTOM = 1.4, -0.6
X_LO, X_HI = 0.5, 1.0
shown = [i for i in ks if X_LO <= xs_data[i] <= X_HI]


def px(x):
    return LEFT + (RIGHT - LEFT) * (x - X_LO) / (X_HI - X_LO)


def py(u):
    return TOP + (BOTTOM - TOP) * (Y_TOP - u) / (Y_TOP - Y_BOTTOM)


def fmt(value):
    return f"{value:.4f}"


svg = []
svg.append(
    f'<svg viewBox="0 0 {W} {H}" role="img" tabindex="0" '
    f'aria-label="Upwind solution at time 0.7, for x from 0.5 to 1. With Courant number 1 the bump is exact. With 0.5 it is lower and wider, with height 0.57. With 1.05 it has sawtooth wiggles from minus 0.46 up to 1.24, near the front of the bump and behind it." '
    f'class="viz-plot">'
)
for value in (1.25, 1.0, 0.75, 0.5, 0.25, 0.0, -0.25, -0.5):
    y = py(value)
    cls = "viz-leader" if value == 0.0 else "viz-grid"
    svg.append(f'<line x1="{LEFT}" x2="{RIGHT}" y1="{y:.1f}" y2="{y:.1f}" class="{cls}"/>')
    svg.append(f'<text x="{LEFT - 10}" y="{y + 4:.1f}" text-anchor="end" class="viz-tick">{value:g}</text>')
for x in (0.5, 0.6, 0.7, 0.8, 0.9, 1.0):
    svg.append(f'<text x="{px(x):.1f}" y="{BOTTOM + 22}" text-anchor="middle" class="viz-tick">{x:g}</text>')
svg.append(f'<text x="{(LEFT + RIGHT) / 2:.1f}" y="{BOTTOM + 48}" text-anchor="middle" class="viz-axis-title">x</text>')
svg.append(f'<text x="{LEFT}" y="{TOP - 16}" class="viz-axis-title">U at t = 0.7, upwind, 100 cells</text>')
exact_path = " ".join(f"{'M' if n == 0 else 'L'}{px(xs_data[i]):.1f},{py(EXACT[i]):.1f}" for n, i in enumerate(shown))
order = ["ν = 1.05", "ν = 0.5", "ν = 1"]
for name in order:
    path = " ".join(f"{'M' if n == 0 else 'L'}{px(xs_data[i]):.1f},{py(SERIES[name][i]):.1f}" for n, i in enumerate(shown))
    svg.append(f'<path d="{path}" fill="none" stroke="{COLORS[name]}" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/>')
svg.append(f'<path d="{exact_path}" fill="none" stroke="#8e8e93" stroke-width="1.5" stroke-dasharray="4 4"/>')
# direct labels: nu = 1 left of its peak, nu = 1.05 left of its wiggles, nu = 0.5 in the right margin
peaks = {}
for name in NAMES:
    i = max(shown, key=lambda j: SERIES[name][j])
    peaks[name] = (xs_data[i], SERIES[name][i])
label_pos = {"ν = 1": (0.78, 0.97), "ν = 1.05": (0.6, 0.47), "ν = 0.5": (1.01, SERIES["ν = 0.5"][shown[-1]] - 0.02)}
for name in NAMES:
    x, y = label_pos[name]
    svg.append(f'<text x="{px(x):.1f}" y="{py(y):.1f}" class="viz-label">{escape(name)}</text>')
svg.append(f'<line x1="0" x2="0" y1="{TOP}" y2="{BOTTOM}" class="viz-crosshair" visibility="hidden"/>')
svg.append("</svg>")

NAMES_TIP = NAMES + ["exact"]
rows = [[fmt(SERIES[name][i]) for name in NAMES] + [fmt(EXACT[i])] for i in shown]
heads = [f"x = {xs_data[i]:.3f}" for i in shown]
ths = [xs_data[i] for i in shown]
COLORS_TIP = [COLORS[n] for n in NAMES] + ["#8e8e93"]
script = (
    "<script>(() => {"
    "const fig = document.currentScript.closest('figure');"
    "const svg = fig.querySelector('svg');"
    "const line = svg.querySelector('.viz-crosshair');"
    "const tip = fig.querySelector('.viz-tooltip');"
    f"const names = {json.dumps(NAMES_TIP)}; const colors = {json.dumps(COLORS_TIP)};"
    f"const rows = {json.dumps(rows)}; const heads = {json.dumps(heads)};"
    f"const xs = {json.dumps([round(px(t), 1) for t in ths])}; const steps = {len(shown) - 1}, width = {W};"
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

legend = "".join(
    f'<span class="viz-key"><span class="viz-swatch" style="background:{COLORS[name]}"></span>{escape(name)}</span>'
    for name in NAMES
) + '<span class="viz-key"><span class="viz-swatch" style="background:#8e8e93"></span>exact</span>'
table_head = "<tr><th>x</th>" + "".join(f"<th>{escape(name)}</th>" for name in NAMES_TIP) + "</tr>"
table_body = "".join(
    f"<tr><td>{xs_data[i]:.3f}</td>" + "".join(f"<td>{v}</td>" for v in rows[n]) + "</tr>"
    for n, i in enumerate(shown) if n % 2 == 0
)
parts = [
    '<figure class="viz">',
    f'<div class="viz-legend">{legend}</div>',
    '<div class="viz-frame">',
    "".join(svg),
    '<div class="viz-tooltip" hidden></div>',
    "</div>",
    "<figcaption>Upwind for u_t + u_x = 0 on a periodic grid with 100 cells, from a cosine bump centered at x = 0.2. "
    "The curves show the right half of the domain at t = 0.7 (for ν = 1.05 the nearest step is t = 0.7035). "
    "With ν = 1 the scheme is exact and hides the dashed exact line. "
    "Hover over the chart, or focus it and use the arrow keys, to read the values. The table shows every second cell.</figcaption>",
    "<details class=\"viz-table\"><summary>Show the numbers</summary>"
    f"<table><thead>{table_head}</thead><tbody>{table_body}</tbody></table></details>",
    script,
    "</figure>",
]
Path(__file__).with_name("figure.html").write_text("\n".join(parts) + "\n", encoding="utf-8")
