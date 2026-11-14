"""Builds the chart for the post from linalg.json (Part A: ||A^m||_2 against m).

Writes one HTML block (legend, SVG, hover script, data table) with no blank
lines to figure.html, so it can sit directly inside the Markdown post.
Same layout as scripts/posts/numerical-analysis/figure.py.

Run: python scripts/posts/numerical-analysis/week-7/figure.py   (writes figure.html next to this file)
"""

import json
import math
from html import escape
from pathlib import Path

data = json.loads(Path(__file__).with_name("linalg.json").read_text(encoding="utf-8"))["part_a"]
ms = data["m"]
series = data["series"]

# Categorical slots 1-4 of the validated reference palette, in fixed order
SERIES = [
    ("symmetric, ρ = 0.924", series["symmetric"], "#2a78d6"),
    ("Jordan block, ρ = 0.8", series["jordan"], "#eb6834"),
    ("triangular, ρ = 0.6", series["triangular"], "#1baf7a"),
    ("symmetric, ρ = 1.05", series["growing"], "#eda100"),
]
SHORT = ["symmetric", "Jordan", "triangular", "ρ = 1.05"]
NAMES = [s[0] for s in SERIES]
COLORS = [s[2] for s in SERIES]

W, H = 680, 400
LEFT, RIGHT, TOP, BOTTOM = 64, 560, 30, 340
Y_TOP, Y_BOTTOM = 3, -8  # log10 of the norm; smaller values are drawn on the bottom edge
LAST = len(ms) - 1


def px(i):
    return LEFT + (RIGHT - LEFT) * ms[i] / ms[-1]


def py(v):
    value = min(max(math.log10(v), Y_BOTTOM), Y_TOP) if v > 0 else Y_BOTTOM
    return TOP + (BOTTOM - TOP) * (Y_TOP - value) / (Y_TOP - Y_BOTTOM)


def sci(value):
    if 0.01 <= value < 1000:
        return f"{value:.3g}"
    mantissa, exponent = f"{value:.2e}".split("e")
    return f"{mantissa}e{int(exponent)}"


def power_label(exponent):
    if exponent == 0:
        return "1"
    return f'10<tspan dy="-0.5em" font-size="75%">{"−" if exponent < 0 else ""}{abs(exponent)}</tspan>'


svg = []
svg.append(
    f'<svg viewBox="0 0 {W} {H}" role="img" tabindex="0" '
    f'aria-label="The 2-norm of A to the power m, for m from 0 to 120, on a log scale, for four matrices. The symmetric matrix with spectral radius 0.924 decays steadily. The Jordan block with spectral radius 0.8 first grows to about 33 near m = 14 and then decays. The triangular matrix with spectral radius 0.6 grows a little and then decays fast. The symmetric matrix with spectral radius 1.05 grows without bound." '
    f'class="viz-plot">'
)
for exponent in range(Y_TOP, Y_BOTTOM - 1, -1):
    y = py(10.0**exponent)
    svg.append(f'<line x1="{LEFT}" x2="{RIGHT}" y1="{y:.1f}" y2="{y:.1f}" class="viz-grid"/>')
    svg.append(f'<text x="{LEFT - 10}" y="{y + 4:.1f}" text-anchor="end" class="viz-tick">{power_label(exponent)}</text>')
for m in range(0, 121, 20):
    svg.append(f'<text x="{px(ms.index(m)):.1f}" y="{BOTTOM + 22}" text-anchor="middle" class="viz-tick">{m}</text>')
svg.append(f'<text x="{(LEFT + RIGHT) / 2:.1f}" y="{BOTTOM + 48}" text-anchor="middle" class="viz-axis-title">Power m</text>')
svg.append(f'<text x="{LEFT}" y="{TOP - 16}" class="viz-axis-title">‖Aᵐ‖₂</text>')

ends = {}
for name, values, color in SERIES:
    points = [(px(i), py(v)) for i, v in enumerate(values)]
    path = " ".join(f"{'M' if i == 0 else 'L'}{x:.1f},{y:.1f}" for i, (x, y) in enumerate(points))
    svg.append(f'<path d="{path}" fill="none" stroke="{color}" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/>')
    x, y = points[-1]
    svg.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="4" fill="{color}" stroke="#ffffff" stroke-width="2"/>')
    ends[name] = (x, y)

# Direct labels
x, y = ends[NAMES[0]]
svg.append(f'<text x="{x + 10:.1f}" y="{y + 4:.1f}" class="viz-label">sym.</text>')
x, y = ends[NAMES[1]]
svg.append(f'<text x="{x + 10:.1f}" y="{y + 4:.1f}" class="viz-label">Jordan</text>')
k = ms.index(20)
svg.append(f'<text x="{px(k) + 12:.1f}" y="{py(SERIES[2][1][k]) + 4:.1f}" class="viz-label">triangular</text>')
x, y = ends[NAMES[3]]
svg.append(f'<text x="{x + 10:.1f}" y="{y + 4:.1f}" class="viz-label">ρ &gt; 1</text>')

svg.append(f'<line x1="0" x2="0" y1="{TOP}" y2="{BOTTOM}" class="viz-crosshair" visibility="hidden"/>')
svg.append("</svg>")

rows = [[sci(values[i]) for _, values, _ in SERIES] for i in range(len(ms))]
TABLE_M = [0, 1, 2, 5, 10, 14, 20, 30, 40, 60, 80, 100, 120]

legend = "".join(
    f'<span class="viz-key"><span class="viz-swatch" style="background:{color}"></span>{escape(name)}</span>'
    for name, _, color in SERIES
)
table_head = "<tr><th>m</th>" + "".join(f"<th>{escape(name)}</th>" for name in NAMES) + "</tr>"
table_body = "".join(
    "<tr><td>" + str(m) + "</td>" + "".join(f"<td>{v}</td>" for v in rows[ms.index(m)]) + "</tr>" for m in TABLE_M
)

script = (
    "<script>(() => {"
    "const fig = document.currentScript.closest('figure');"
    "const svg = fig.querySelector('svg');"
    "const line = svg.querySelector('.viz-crosshair');"
    "const tip = fig.querySelector('.viz-tooltip');"
    f"const names = {json.dumps(SHORT, ensure_ascii=False)}; const colors = {json.dumps(COLORS)};"
    f"const rows = {json.dumps(rows)}; const ms = {json.dumps(ms)};"
    f"const width = {W}, last = ms.length - 1;"
    f"const xs = {json.dumps([round(px(i), 1) for i in range(len(ms))])};"
    "let current = -1;"
    "const show = (k) => {"
    " k = Math.max(0, Math.min(last, k)); current = k;"
    " const x = xs[k];"
    " line.setAttribute('x1', x); line.setAttribute('x2', x); line.setAttribute('visibility', 'visible');"
    " tip.replaceChildren();"
    " const head = document.createElement('div'); head.className = 'viz-tooltip-head'; head.textContent = 'm = ' + ms[k]; tip.append(head);"
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
    "<figcaption>The 2-norm of the power A<sup>m</sup> for four matrices, on a log scale. "
    "Values below 10<sup>−8</sup> are drawn on the bottom edge (the table has the real values). "
    "Hover over the chart, or focus it and use the arrow keys, to read the values.</figcaption>",
    "<details class=\"viz-table\"><summary>Show the numbers</summary>"
    f"<table><thead>{table_head}</thead><tbody>{table_body}</tbody></table></details>",
    script,
    "</figure>",
]
Path(__file__).with_name("figure.html").write_text("\n".join(parts) + "\n", encoding="utf-8")
