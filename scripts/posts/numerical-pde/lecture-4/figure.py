"""Builds the chart for the post from experiment.json (part 4).

Writes one HTML block (legend, SVG, hover script, data table) with no blank
lines to figure.html, so it can sit directly inside the Markdown post.
Same layout as scripts/posts/numerical-analysis/figure.py.

Run: python scripts/posts/numerical-pde/lecture-4/figure.py   (writes figure.html next to this file)
"""

import json
from html import escape
from pathlib import Path

data = json.loads(Path(__file__).with_name("experiment.json").read_text(encoding="utf-8"))["part4"]
xs_all = data["xs"]
# Every second sample: x = 0, 0.01, ..., 1
PICK = list(range(0, len(xs_all), 2))
XS = [xs_all[k] for k in PICK]

# Categorical slots 1-4 of the validated reference palette, in fixed order
COLORS = {
    "ε = 1": "#2a78d6",
    "ε = 0.1": "#eb6834",
    "ε = 0.03": "#1baf7a",
    "ε = 0.01": "#eda100",
}
KEYS = {"ε = 1": "1.0", "ε = 0.1": "0.1", "ε = 0.03": "0.03", "ε = 0.01": "0.01"}
NAMES = list(COLORS)
SERIES = {name: [data["curves"][KEYS[name]][k] for k in PICK] for name in NAMES}
REF = "#6b6b6b"

W, H = 680, 420
LEFT, RIGHT, TOP, BOTTOM = 64, 540, 30, 360
Y_TOP, Y_BOTTOM = 1.0, -1.0


def px(x):
    return LEFT + (RIGHT - LEFT) * x


def py(v):
    return TOP + (BOTTOM - TOP) * (Y_TOP - v) / (Y_TOP - Y_BOTTOM)


svg = []
svg.append(
    f'<svg viewBox="0 0 {W} {H}" role="img" tabindex="0" '
    f'aria-label="Exact solutions of minus epsilon u double prime plus u prime equals 1 with u(0) = u(1) = 0, for four values of epsilon. '
    f'For epsilon = 1 the curve is a low bump. As epsilon gets smaller, the curve follows the line u = x and drops to zero in a thin layer near x = 1. '
    f'It never comes close to the other candidate line u = x minus 1." '
    f'class="viz-plot">'
)
for v in (1.0, 0.5, 0.0, -0.5, -1.0):
    y = py(v)
    svg.append(f'<line x1="{LEFT}" x2="{RIGHT}" y1="{y:.1f}" y2="{y:.1f}" class="viz-grid"/>')
    label = f"{v:g}".replace("-", "−")
    svg.append(f'<text x="{LEFT - 10}" y="{y + 4:.1f}" text-anchor="end" class="viz-tick">{label}</text>')
for x in (0, 0.25, 0.5, 0.75, 1):
    svg.append(f'<text x="{px(x):.1f}" y="{BOTTOM + 22}" text-anchor="middle" class="viz-tick">{x:g}</text>')
svg.append(f'<text x="{(LEFT + RIGHT) / 2:.1f}" y="{BOTTOM + 48}" text-anchor="middle" class="viz-axis-title">x</text>')
svg.append(f'<text x="{LEFT}" y="{TOP - 16}" class="viz-axis-title">u(x)</text>')

# The two candidate solutions of the reduced equation u' = 1 (dashed)
for v0, v1 in [(0.0, 1.0), (-1.0, 0.0)]:
    svg.append(
        f'<line x1="{px(0):.1f}" y1="{py(v0):.1f}" x2="{px(1):.1f}" y2="{py(v1):.1f}" '
        f'stroke="{REF}" stroke-width="1.5" stroke-dasharray="5 4"/>'
    )
svg.append(f'<text x="{RIGHT + 10}" y="{py(1.0) - 8:.1f}" class="viz-label">u = x (inflow side)</text>')
# Below the lower line, away from the curves
svg.append(f'<text x="{px(0.45):.1f}" y="{py(-0.55) + 22:.1f}" class="viz-label">{escape("u = x − 1 (outflow side)")}</text>')

for name in NAMES:
    vals = SERIES[name]
    path = " ".join(f"{'M' if i == 0 else 'L'}{px(x):.1f},{py(v):.1f}" for i, (x, v) in enumerate(zip(XS, vals)))
    svg.append(f'<path d="{path}" fill="none" stroke="{COLORS[name]}" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/>')

# Direct labels next to the peak of each curve
LABEL_POS = {"ε = 1": (0.5, -14), "ε = 0.1": (0.62, -16), "ε = 0.03": (0.80, -14), "ε = 0.01": (0.86, -24)}
for name, (x, dy) in LABEL_POS.items():
    k = round(x * 100)
    svg.append(f'<text x="{px(x):.1f}" y="{py(SERIES[name][k]) + dy:.1f}" text-anchor="middle" class="viz-label">{escape(name)}</text>')

svg.append(f'<line x1="0" x2="0" y1="{TOP}" y2="{BOTTOM}" class="viz-crosshair" visibility="hidden"/>')
svg.append("</svg>")

rows = [[f"{SERIES[name][i]:.3f}" for name in NAMES] for i in range(len(XS))]
heads = [f"x = {x:.2f}" for x in XS]

legend = "".join(
    f'<span class="viz-key"><span class="viz-swatch" style="background:{COLORS[name]}"></span>{escape(name)}</span>'
    for name in NAMES
)
legend += f'<span class="viz-key"><span class="viz-swatch" style="background:{REF}"></span>Reduced solutions (dashed)</span>'
table_head = "<tr><th>x</th>" + "".join(f"<th>{escape(name)}</th>" for name in NAMES) + "</tr>"
TABLE_ROWS = range(0, len(XS), 10)
table_body = "".join(
    f"<tr><td>{XS[i]:.1f}</td>" + "".join(f"<td>{v}</td>" for v in rows[i]) + "</tr>" for i in TABLE_ROWS
)
# Extra rows inside the layer, where the curves change fast
table_body += "".join(
    f"<tr><td>{XS[i]:.2f}</td>" + "".join(f"<td>{v}</td>" for v in rows[i]) + "</tr>" for i in (95, 98, 99)
)

script = (
    "<script>(() => {"
    "const fig = document.currentScript.closest('figure');"
    "const svg = fig.querySelector('svg');"
    "const line = svg.querySelector('.viz-crosshair');"
    "const tip = fig.querySelector('.viz-tooltip');"
    f"const names = {json.dumps(NAMES, ensure_ascii=False)}; const colors = {json.dumps([COLORS[n] for n in NAMES])};"
    f"const rows = {json.dumps(rows)}; const heads = {json.dumps(heads)};"
    f"const left = {LEFT}, right = {RIGHT}, steps = {len(XS) - 1}, width = {W};"
    "let current = -1;"
    "const show = (n) => {"
    " n = Math.max(0, Math.min(steps, n)); current = n;"
    " const x = left + (right - left) * n / steps;"
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
    " show(Math.round((x - left) / (right - left) * steps));"
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
    "<figcaption>Exact solutions of −εu″ + u′ = 1 on (0, 1) with u(0) = u(1) = 0. "
    "The dashed lines are the two solutions of the reduced equation u′ = 1 that keep one boundary condition each. "
    "Hover over the chart, or focus it and use the arrow keys, to read the values.</figcaption>",
    "<details class=\"viz-table\"><summary>Show the numbers</summary>"
    f"<table><thead>{table_head}</thead><tbody>{table_body}</tbody></table></details>",
    script,
    "</figure>",
]
Path(__file__).with_name("figure.html").write_text("\n".join(parts) + "\n", encoding="utf-8")
