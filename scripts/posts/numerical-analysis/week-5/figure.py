"""Builds the error chart for the post from ode_methods.json (Part A).

Writes one HTML block (legend, SVG, hover script, data table) with no blank
lines to figure.html, so it can sit directly inside the Markdown post.
Same layout as scripts/posts/numerical-analysis/figure.py.

Run: python scripts/posts/numerical-analysis/week-5/figure.py   (writes figure.html next to this file)
"""

import json
import math
from html import escape
from pathlib import Path

data = json.loads(Path(__file__).with_name("ode_methods.json").read_text(encoding="utf-8"))
rows_a = data["part_a"]
degrees = [r["n"] for r in rows_a]  # number of steps n on [0, 4]; h = 4 / n

# Categorical slots 1-4 of the validated reference palette, in fixed order
SERIES = [
    ("Euler", [r["euler"] for r in rows_a], "#2a78d6"),
    ("midpoint", [r["midpoint"] for r in rows_a], "#eb6834"),
    ("trapezoidal", [r["trapezoidal"] for r in rows_a], "#1baf7a"),
    ("Euler error bound", [r["euler_bound"] for r in rows_a], "#eda100"),
]
NAMES = [s[0] for s in SERIES]
COLORS = [s[2] for s in SERIES]
DASHED = {"Euler error bound"}

W, H = 680, 400
LEFT, RIGHT, TOP, BOTTOM = 64, 540, 30, 340
Y_TOP, Y_BOTTOM = 0, -8  # log10 of the error
LAST = len(degrees) - 1


def px(i):
    """Index i (n = 20 * 2^i steps) is equally spaced: a log scale in h."""
    return LEFT + (RIGHT - LEFT) * i / LAST


def py(err):
    value = max(math.log10(err), Y_BOTTOM) if err > 0 else Y_BOTTOM
    return TOP + (BOTTOM - TOP) * (Y_TOP - value) / (Y_TOP - Y_BOTTOM)


def sci(value):
    mantissa, exponent = f"{value:.1e}".split("e")
    return f"{mantissa}e{int(exponent)}"


def power_label(exponent):
    if exponent == 0:
        return "1"
    return f'10<tspan dy="-0.5em" font-size="75%">{"−" if exponent < 0 else ""}{abs(exponent)}</tspan>'


svg = []
svg.append(
    f'<svg viewBox="0 0 {W} {H}" role="img" tabindex="0" '
    f'aria-label="Maximum error of three methods for y prime = y minus y squared on [0, 4] against the step size, both on log scales. The Euler error halves each time h halves. The midpoint and trapezoidal errors drop by 4. The dashed Euler error bound stays above the Euler error by a factor of about 11." '
    f'class="viz-plot">'
)
for exponent in range(Y_TOP, Y_BOTTOM - 1, -2):
    y = py(10.0**exponent)
    svg.append(f'<line x1="{LEFT}" x2="{RIGHT}" y1="{y:.1f}" y2="{y:.1f}" class="viz-grid"/>')
    svg.append(f'<text x="{LEFT - 10}" y="{y + 4:.1f}" text-anchor="end" class="viz-tick">{power_label(exponent)}</text>')
for i, n in enumerate(degrees):
    if i % 2 == 0:
        svg.append(f'<text x="{px(i):.1f}" y="{BOTTOM + 22}" text-anchor="middle" class="viz-tick">{4 / n:g}</text>')
svg.append(f'<text x="{(LEFT + RIGHT) / 2:.1f}" y="{BOTTOM + 48}" text-anchor="middle" class="viz-axis-title">Step size h (halved at each step to the right)</text>')
svg.append(f'<text x="{LEFT}" y="{TOP - 16}" class="viz-axis-title">|error|</text>')

# Lines
ends = {}
for name, values, color in SERIES:
    points = []
    for i, v in enumerate(values):
        points.append((px(i), py(v)))
        if v < 1e-15:  # round-off level: stop the line here
            break
    path = " ".join(f"{'M' if i == 0 else 'L'}{x:.1f},{y:.1f}" for i, (x, y) in enumerate(points))
    dash = ' stroke-dasharray="6 4"' if name in DASHED else ""
    svg.append(f'<path d="{path}" fill="none" stroke="{color}" stroke-width="2"{dash} stroke-linejoin="round" stroke-linecap="round"/>')
    x, y = points[-1]
    svg.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="4" fill="{color}" stroke="#ffffff" stroke-width="2"/>')
    ends[name] = (x, y)

# Direct labels
x, y = ends[NAMES[0]]
svg.append(f'<text x="{x + 10:.1f}" y="{y + 4:.1f}" class="viz-label">Euler</text>')
x, y = ends[NAMES[1]]
svg.append(f'<text x="{x + 10:.1f}" y="{y - 4:.1f}" class="viz-label">midpoint</text>')
x, y = ends[NAMES[2]]
svg.append(f'<text x="{x + 10:.1f}" y="{y + 12:.1f}" class="viz-label">trapezoidal</text>')
x, y = ends[NAMES[3]]
svg.append(f'<text x="{x + 10:.1f}" y="{y + 4:.1f}" class="viz-label">bound</text>')

svg.append(f'<line x1="0" x2="0" y1="{TOP}" y2="{BOTTOM}" class="viz-crosshair" visibility="hidden"/>')
svg.append("</svg>")

rows = [[sci(values[i]) for _, values, _ in SERIES] for i in range(len(degrees))]

legend = "".join(
    f'<span class="viz-key"><span class="viz-swatch" style="background:{color}"></span>{escape(name)}</span>'
    for name, _, color in SERIES
)
table_head = "<tr><th>h</th>" + "".join(f"<th>{escape(name)}</th>" for name in NAMES) + "</tr>"
table_body = "".join(
    "<tr><td>" + f"{4 / n:g}" + "</td>" + "".join(f"<td>{v}</td>" for v in row) + "</tr>" for n, row in zip(degrees, rows)
)

script = (
    "<script>(() => {"
    "const fig = document.currentScript.closest('figure');"
    "const svg = fig.querySelector('svg');"
    "const line = svg.querySelector('.viz-crosshair');"
    "const tip = fig.querySelector('.viz-tooltip');"
    f"const names = {json.dumps(NAMES, ensure_ascii=False)}; const colors = {json.dumps(COLORS)};"
    f"const rows = {json.dumps(rows)}; const degrees = {json.dumps(degrees)};"
    f"const left = {LEFT}, right = {RIGHT}, width = {W}, last = degrees.length - 1;"
    "let current = -1;"
    "const show = (k) => {"
    " k = Math.max(0, Math.min(last, k)); current = k;"
    " const x = left + (right - left) * k / last;"
    " line.setAttribute('x1', x); line.setAttribute('x2', x); line.setAttribute('visibility', 'visible');"
    " tip.replaceChildren();"
    " const head = document.createElement('div'); head.className = 'viz-tooltip-head'; head.textContent = 'h = ' + (4 / degrees[k]) + ' (' + degrees[k] + ' steps)'; tip.append(head);"
    " names.forEach((name, i) => {"
    "  const row = document.createElement('div'); row.className = 'viz-tooltip-row';"
    "  const key = document.createElement('span'); key.className = 'viz-swatch'; key.style.background = colors[i];"
    "  const value = document.createElement('strong'); value.textContent = rows[k][i];"
    "  const label = document.createElement('span'); label.textContent = name;"
    "  row.append(key, value, label); tip.append(row);"
    " });"
    " const box = svg.getBoundingClientRect(), frame = fig.getBoundingClientRect();"
    " const xs = box.left - frame.left + x * box.width / width;"
    " tip.hidden = false;"
    " const room = frame.width - tip.offsetWidth - 8;"
    " tip.style.left = Math.max(0, Math.min(room, xs + 12)) + 'px';"
    " tip.style.top = (box.top - frame.top + 8) + 'px';"
    "};"
    "const hide = () => { line.setAttribute('visibility', 'hidden'); tip.hidden = true; current = -1; };"
    "svg.addEventListener('pointermove', (event) => {"
    " const box = svg.getBoundingClientRect();"
    " const x = (event.clientX - box.left) * width / box.width;"
    " show(Math.round((x - left) / (right - left) * last));"
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
    "<figcaption>Maximum error over the grid for y′ = y − y², y(0) = 0.1 on [0, 4], against the step size h, both on log scales. "
    "Each step to the right halves h. The dashed line is the Euler error bound from the convergence proof. "
    "Hover over the chart, or focus it and use the arrow keys, to read the values.</figcaption>",
    "<details class=\"viz-table\"><summary>Show the numbers</summary>"
    f"<table><thead>{table_head}</thead><tbody>{table_body}</tbody></table></details>",
    script,
    "</figure>",
]
Path(__file__).with_name("figure.html").write_text("\n".join(parts) + "\n", encoding="utf-8")
