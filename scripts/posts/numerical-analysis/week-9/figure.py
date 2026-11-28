"""Builds the chart for the post from descent.json (Part B: error against iteration, cond = 100).

Writes one HTML block (legend, SVG, hover script, data table) with no blank
lines to figure.html, so it can sit directly inside the Markdown post.
Same layout as scripts/posts/numerical-analysis/figure.py.

Run: python scripts/posts/numerical-analysis/week-9/figure.py   (writes figure.html next to this file)
"""

import json
import math
from html import escape
from pathlib import Path

data = json.loads(Path(__file__).with_name("descent.json").read_text(encoding="utf-8"))["part_b"]["history_kappa100"]
LAST_IT = 400
ms = list(range(LAST_IT + 1))

# Categorical slots 1-4 of the validated reference palette, in fixed order
SERIES = [
    ("fixed α = 1/λmax", data["fixed_1_over_lmax"], "#2a78d6"),
    ("fixed α = 2/(λmin + λmax)", data["fixed_best"], "#eb6834"),
    ("steepest descent, line search", data["sd"], "#1baf7a"),
    ("conjugate gradient", data["cg"], "#eda100"),
]
SHORT = ["α = 1/λmax", "best fixed α", "steepest descent", "CG"]
NAMES = [s[0] for s in SERIES]
COLORS = [s[2] for s in SERIES]

W, H = 680, 400
LEFT, RIGHT, TOP, BOTTOM = 64, 560, 30, 340
Y_TOP, Y_BOTTOM = 0, -14


def px(i):
    return LEFT + (RIGHT - LEFT) * i / LAST_IT


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
    f'aria-label="Relative error against iteration number, on a log scale, for a 100 by 100 diagonal matrix with condition number 100. The fixed step 1 over lambda max falls slowly and is at about 2e-4 after 400 iterations. The best fixed step reaches about 3e-5, and steepest descent with line search about 4e-6. Conjugate gradient falls much faster and reaches 1e-14 after 72 iterations." '
    f'class="viz-plot">'
)
for exponent in range(Y_TOP, Y_BOTTOM - 1, -2):
    y = py(10.0**exponent)
    svg.append(f'<line x1="{LEFT}" x2="{RIGHT}" y1="{y:.1f}" y2="{y:.1f}" class="viz-grid"/>')
    svg.append(f'<text x="{LEFT - 10}" y="{y + 4:.1f}" text-anchor="end" class="viz-tick">{power_label(exponent)}</text>')
for m in range(0, LAST_IT + 1, 50):
    svg.append(f'<text x="{px(m):.1f}" y="{BOTTOM + 22}" text-anchor="middle" class="viz-tick">{m}</text>')
svg.append(f'<text x="{(LEFT + RIGHT) / 2:.1f}" y="{BOTTOM + 48}" text-anchor="middle" class="viz-axis-title">Iteration k</text>')
svg.append(f'<text x="{LEFT}" y="{TOP - 16}" class="viz-axis-title">Relative error</text>')

ends = {}
for name, values, color in SERIES:
    points = [(px(i), py(v)) for i, v in enumerate(values)]
    path = " ".join(f"{'M' if i == 0 else 'L'}{x:.1f},{y:.1f}" for i, (x, y) in enumerate(points))
    svg.append(f'<path d="{path}" fill="none" stroke="{color}" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/>')
    x, y = points[-1]
    svg.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="4" fill="{color}" stroke="#ffffff" stroke-width="2"/>')
    ends[name] = (x, y)

# Direct labels (the three gradient lines end close together, so spread the labels a little)
labels = ["1/λmax", "best α", "SD", "CG"]
ys = [ends[name][1] for name in NAMES]
order = sorted(range(3), key=lambda i: ys[i])
placed = []
for i in order:
    y = ys[i] if not placed else max(ys[i], placed[-1] + 16)
    placed.append(y)
    svg.append(f'<text x="{ends[NAMES[i]][0] + 10:.1f}" y="{y + 4:.1f}" class="viz-label">{labels[i]}</text>')
x, y = ends[NAMES[3]]
svg.append(f'<text x="{x + 10:.1f}" y="{y + 4:.1f}" class="viz-label">CG</text>')

svg.append(f'<line x1="0" x2="0" y1="{TOP}" y2="{BOTTOM}" class="viz-crosshair" visibility="hidden"/>')
svg.append("</svg>")

rows = [[sci(values[i]) if i < len(values) else "done" for _, values, _ in SERIES] for i in range(len(ms))]
TABLE_M = [0, 1, 5, 10, 20, 30, 40, 50, 60, 72, 100, 150, 200, 300, 400]
legend = "".join(
    f'<span class="viz-key"><span class="viz-swatch" style="background:{color}"></span>{escape(name)}</span>'
    for name, _, color in SERIES
)
table_head = "<tr><th>k</th>" + "".join(f"<th>{escape(name)}</th>" for name in NAMES) + "</tr>"
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
    " const head = document.createElement('div'); head.className = 'viz-tooltip-head'; head.textContent = 'k = ' + ms[k]; tip.append(head);"
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
    "<figcaption>Relative error of four methods for the same 100 × 100 SPD matrix (eigenvalues spread evenly from 1 to 100) and the same random start, on a log scale. "
    "CG stops when the error is below 10<sup>−14</sup>; \"done\" in the table means it has stopped. "
    "Hover over the chart, or focus it and use the arrow keys, to read the values.</figcaption>",
    "<details class=\"viz-table\"><summary>Show the numbers</summary>"
    f"<table><thead>{table_head}</thead><tbody>{table_body}</tbody></table></details>",
    script,
    "</figure>",
]
Path(__file__).with_name("figure.html").write_text("\n".join(parts) + "\n", encoding="utf-8")
