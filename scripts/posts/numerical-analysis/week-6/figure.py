"""Builds the error chart for the post from stability.json (Part A).

Writes one HTML block (legend, SVG, hover script, data table) with no blank
lines to figure.html, so it can sit directly inside the Markdown post.
Same layout as scripts/posts/numerical-analysis/figure.py.

Run: python scripts/posts/numerical-analysis/week-6/figure.py   (writes figure.html next to this file)
"""

import json
import math
from html import escape
from pathlib import Path

data = json.loads(Path(__file__).with_name("stability.json").read_text(encoding="utf-8"))
rows_a = data["part_a"]
degrees = [r["n"] for r in rows_a]  # number of steps n on [0, 1]; h = 1 / n

# Categorical slots 1-4 of the validated reference palette, in fixed order
SERIES = [
    ("forward Euler", [r["forward_euler"] for r in rows_a], "#2a78d6"),
    ("RK4", [r["rk4"] for r in rows_a], "#eb6834"),
    ("backward Euler", [r["backward_euler"] for r in rows_a], "#1baf7a"),
    ("trapezoidal", [r["trapezoidal"] for r in rows_a], "#eda100"),
]
NAMES = [s[0] for s in SERIES]
COLORS = [s[2] for s in SERIES]
DASHED = set()

W, H = 680, 400
LEFT, RIGHT, TOP, BOTTOM = 64, 540, 30, 340
Y_TOP, Y_BOTTOM = 2, -12  # log10 of the error; larger errors are drawn at the top edge
LAST = len(degrees) - 1


H_MAX, H_MIN = math.log10(1 / degrees[0]), math.log10(1 / degrees[-1])


def px(i):
    """Position of h = 1 / n on a log scale (largest h on the left)."""
    return LEFT + (RIGHT - LEFT) * (H_MAX - math.log10(1 / degrees[i])) / (H_MAX - H_MIN)


def py(err):
    if not math.isfinite(err):
        return TOP
    value = min(max(math.log10(err), Y_BOTTOM), Y_TOP) if err > 0 else Y_BOTTOM
    return TOP + (BOTTOM - TOP) * (Y_TOP - value) / (Y_TOP - Y_BOTTOM)


def sci(value):
    if not math.isfinite(value):
        return "overflow"
    mantissa, exponent = f"{value:.1e}".split("e")
    return f"{mantissa}e{int(exponent)}"


def power_label(exponent):
    if exponent == 0:
        return "1"
    return f'10<tspan dy="-0.5em" font-size="75%">{"−" if exponent < 0 else ""}{abs(exponent)}</tspan>'


svg = []
svg.append(
    f'<svg viewBox="0 0 {W} {H}" role="img" tabindex="0" '
    f'aria-label="Maximum error of four methods for a stiff problem with lambda = minus 1000 on [0, 1] against the step size, both on log scales. Forward Euler is off the top of the chart until h is below 0.002, and RK4 until h is below about 0.0028. Then both errors drop quickly. Backward Euler and the trapezoidal method give small errors for every h." '
    f'class="viz-plot">'
)
for exponent in range(Y_TOP, Y_BOTTOM - 1, -2):
    y = py(10.0**exponent)
    svg.append(f'<line x1="{LEFT}" x2="{RIGHT}" y1="{y:.1f}" y2="{y:.1f}" class="viz-grid"/>')
    svg.append(f'<text x="{LEFT - 10}" y="{y + 4:.1f}" text-anchor="end" class="viz-tick">{power_label(exponent)}</text>')
for n in (50, 100, 200, 1000, 2000, 8000):
    i = degrees.index(n)
    svg.append(f'<text x="{px(i):.1f}" y="{BOTTOM + 22}" text-anchor="middle" class="viz-tick">{1 / n:g}</text>')
for limit in (0.002, 0.002785):
    x = LEFT + (RIGHT - LEFT) * (H_MAX - math.log10(limit)) / (H_MAX - H_MIN)
    svg.append(f'<line x1="{x:.1f}" x2="{x:.1f}" y1="{TOP}" y2="{BOTTOM}" stroke="#86868b" stroke-width="1" stroke-dasharray="3 4"/>')
svg.append(f'<text x="{(LEFT + RIGHT) / 2:.1f}" y="{BOTTOM + 48}" text-anchor="middle" class="viz-axis-title">Step size h (smaller to the right)</text>')
svg.append(f'<text x="{LEFT}" y="{TOP - 16}" class="viz-axis-title">|error|</text>')

# Lines
ends = {}
for name, values, color in SERIES:
    points = []
    for i, v in enumerate(values):
        points.append((px(i), py(v)))
    path = " ".join(f"{'M' if i == 0 else 'L'}{x:.1f},{y:.1f}" for i, (x, y) in enumerate(points))
    dash = ' stroke-dasharray="6 4"' if name in DASHED else ""
    svg.append(f'<path d="{path}" fill="none" stroke="{color}" stroke-width="2"{dash} stroke-linejoin="round" stroke-linecap="round"/>')
    x, y = points[-1]
    svg.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="4" fill="{color}" stroke="#ffffff" stroke-width="2"/>')
    ends[name] = (x, y)

# Direct labels
x, y = ends[NAMES[0]]
svg.append(f'<text x="{x + 10:.1f}" y="{y - 2:.1f}" class="viz-label">FE</text>')
x, y = ends[NAMES[2]]
svg.append(f'<text x="{x + 10:.1f}" y="{y + 12:.1f}" class="viz-label">BE</text>')
k = degrees.index(2000)  # RK4 and trapezoidal end close together: label RK4 mid-line
x, y = px(k), py(SERIES[1][1][k])
svg.append(f'<text x="{x + 6:.1f}" y="{y - 8:.1f}" class="viz-label">RK4</text>')
x, y = ends[NAMES[3]]
svg.append(f'<text x="{x + 10:.1f}" y="{y + 10:.1f}" class="viz-label">trap.</text>')

svg.append(f'<line x1="0" x2="0" y1="{TOP}" y2="{BOTTOM}" class="viz-crosshair" visibility="hidden"/>')
svg.append("</svg>")

rows = [[sci(values[i]) for _, values, _ in SERIES] for i in range(len(degrees))]

legend = "".join(
    f'<span class="viz-key"><span class="viz-swatch" style="background:{color}"></span>{escape(name)}</span>'
    for name, _, color in SERIES
)
table_head = "<tr><th>h</th>" + "".join(f"<th>{escape(name)}</th>" for name in NAMES) + "</tr>"
table_body = "".join(
    "<tr><td>" + f"1/{n}" + "</td>" + "".join(f"<td>{v}</td>" for v in row) + "</tr>" for n, row in zip(degrees, rows)
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
    f"const xs = {json.dumps([round(px(i), 1) for i in range(len(degrees))])};"
    "let current = -1;"
    "const show = (k) => {"
    " k = Math.max(0, Math.min(last, k)); current = k;"
    " const x = xs[k];"
    " line.setAttribute('x1', x); line.setAttribute('x2', x); line.setAttribute('visibility', 'visible');"
    " tip.replaceChildren();"
    " const head = document.createElement('div'); head.className = 'viz-tooltip-head'; head.textContent = 'h = 1/' + degrees[k] + ', z = hλ = ' + (-1000 / degrees[k]).toFixed(2); tip.append(head);"
    " names.forEach((name, i) => {"
    "  const row = document.createElement('div'); row.className = 'viz-tooltip-row';"
    "  const key = document.createElement('span'); key.className = 'viz-swatch'; key.style.background = colors[i];"
    "  const value = document.createElement('strong'); value.textContent = rows[k][i];"
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
    "<figcaption>Maximum error over the grid for y′ = −1000(y − sin x) + cos x, y(0) = 0 on [0, 1], against the step size h, both on log scales. "
    "Errors above 10² are drawn on the top edge (the table has the real values). The two dotted vertical lines are the stability limits h = 2/1000 for forward Euler and h = 2.785/1000 for RK4. "
    "Hover over the chart, or focus it and use the arrow keys, to read the values.</figcaption>",
    "<details class=\"viz-table\"><summary>Show the numbers</summary>"
    f"<table><thead>{table_head}</thead><tbody>{table_body}</tbody></table></details>",
    script,
    "</figure>",
]
Path(__file__).with_name("figure.html").write_text("\n".join(parts) + "\n", encoding="utf-8")
