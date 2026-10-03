"""Builds the error chart for the post from rootfinding.json.

Writes one HTML block (legend, SVG, hover script, data table) with no blank
lines to figure.html, so it can sit directly inside the Markdown post.

Run: python scripts/posts/numerical-analysis/figure.py   (writes figure.html next to this file)
"""

import json
import math
from html import escape
from pathlib import Path

data = json.loads(Path(__file__).with_name("rootfinding.json").read_text())
steps = data["steps"]

# Categorical slots 1-4 of the validated reference palette, in fixed order
COLORS = {
    "Bisection": "#2a78d6",
    "Fixed point": "#eb6834",
    "Secant": "#1baf7a",
    "Newton": "#eda100",
}
NAMES = list(COLORS)

W, H = 680, 400
LEFT, RIGHT, TOP, BOTTOM = 64, 600, 30, 340
Y_TOP, Y_BOTTOM = 0, -16  # log10 of the error


def px(n):
    return LEFT + (RIGHT - LEFT) * n / steps


def py(err):
    value = Y_BOTTOM if err is None else max(math.log10(err), Y_BOTTOM)
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
    f'aria-label="Error per step for four methods on a log scale. Newton and the secant method reach round-off within 7 steps; bisection and fixed-point iteration fall at a steady rate." '
    f'class="viz-plot">'
)
# Grid and y-axis
for exponent in range(Y_TOP, Y_BOTTOM - 1, -4):
    y = py(10.0**exponent)
    svg.append(f'<line x1="{LEFT}" x2="{RIGHT}" y1="{y:.1f}" y2="{y:.1f}" class="viz-grid"/>')
    svg.append(f'<text x="{LEFT - 10}" y="{y + 4:.1f}" text-anchor="end" class="viz-tick">{power_label(exponent)}</text>')
# x-axis
for n in range(0, steps + 1, 10):
    svg.append(f'<text x="{px(n):.1f}" y="{BOTTOM + 22}" text-anchor="middle" class="viz-tick">{n}</text>')
svg.append(f'<text x="{(LEFT + RIGHT) / 2:.1f}" y="{BOTTOM + 48}" text-anchor="middle" class="viz-axis-title">Step n</text>')
svg.append(f'<text x="{LEFT}" y="{TOP - 16}" class="viz-axis-title">Error</text>')

# Lines. A drop to the bottom edge means the error went below round-off level.
for name in NAMES:
    s = data["series"][name]
    points = [(px(n), py(e)) for n, e in enumerate(s["errors"])]
    if s["round_off_at"] is not None:
        points.append((px(s["round_off_at"]), py(None)))
    path = " ".join(f"{'M' if i == 0 else 'L'}{x:.1f},{y:.1f}" for i, (x, y) in enumerate(points))
    color = COLORS[name]
    svg.append(f'<path d="{path}" fill="none" stroke="{color}" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/>')
    x, y = points[-1]
    svg.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="4" fill="{color}" stroke="#ffffff" stroke-width="2"/>')

# Direct labels: at the line end where there is room, with a leader line where the ends crowd
for name in ("Bisection", "Fixed point"):
    errs = data["series"][name]["errors"]
    svg.append(f'<text x="{RIGHT + 10}" y="{py(errs[-1]) + 4:.1f}" class="viz-label">{name}</text>')


def leader(name, at_n, label_x, label_y):
    errs = data["series"][name]["errors"]
    lo = math.floor(at_n)
    t = at_n - lo
    log_err = (1 - t) * math.log10(errs[lo]) + t * math.log10(errs[lo + 1])
    x0, y0 = px(at_n), py(10.0**log_err)
    svg.append(f'<line x1="{x0 + 3:.1f}" y1="{y0:.1f}" x2="{label_x - 4:.1f}" y2="{label_y - 4:.1f}" class="viz-leader"/>')
    svg.append(f'<text x="{label_x:.1f}" y="{label_y:.1f}" class="viz-label">{name}</text>')


# Newton sits left of its own line: every spot to its right is crossed by the secant line
svg.append(f'<text x="{px(2.3) - 8:.1f}" y="{py(1e-6) + 4:.1f}" text-anchor="end" class="viz-label">Newton</text>')
leader("Secant", 5.6, px(12), py(1e-11))

svg.append(f'<line x1="0" x2="0" y1="{TOP}" y2="{BOTTOM}" class="viz-crosshair" visibility="hidden"/>')
svg.append("</svg>")

# Hover readout data, also used for the table
rows = []
for n in range(steps + 1):
    row = []
    for name in NAMES:
        s = data["series"][name]
        if n < len(s["errors"]):
            row.append(sci(s["errors"][n]))
        elif s["round_off_at"] is not None and n >= s["round_off_at"]:
            row.append("< 1e-15")
        else:
            row.append("")
    rows.append(row)

legend = "".join(
    f'<span class="viz-key"><span class="viz-swatch" style="background:{COLORS[name]}"></span>{escape(name)}</span>'
    for name in NAMES
)
table_head = "<tr><th>n</th>" + "".join(f"<th>{escape(name)}</th>" for name in NAMES) + "</tr>"
table_body = "".join("<tr><td>" + str(n) + "</td>" + "".join(f"<td>{v}</td>" for v in row) + "</tr>" for n, row in enumerate(rows))

script = (
    "<script>(() => {"
    "const fig = document.currentScript.closest('figure');"
    "const svg = fig.querySelector('svg');"
    "const line = svg.querySelector('.viz-crosshair');"
    "const tip = fig.querySelector('.viz-tooltip');"
    f"const names = {json.dumps(NAMES)}; const colors = {json.dumps([COLORS[n] for n in NAMES])};"
    f"const rows = {json.dumps(rows)};"
    f"const left = {LEFT}, right = {RIGHT}, steps = {steps}, width = {W};"
    "let current = -1;"
    "const show = (n) => {"
    " n = Math.max(0, Math.min(steps, n)); current = n;"
    " const x = left + (right - left) * n / steps;"
    " line.setAttribute('x1', x); line.setAttribute('x2', x); line.setAttribute('visibility', 'visible');"
    " tip.replaceChildren();"
    " const head = document.createElement('div'); head.className = 'viz-tooltip-head'; head.textContent = 'Step ' + n; tip.append(head);"
    " names.forEach((name, i) => {"
    "  if (!rows[n][i]) return;"
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
    "<figcaption>Error |x<sub>n</sub> − α| after each step, on a log scale. "
    "A line that drops to the bottom edge has gone below 10<sup>−15</sup>, the level of round-off. "
    "Hover over the chart, or focus it and use the arrow keys, to read the values.</figcaption>",
    "<details class=\"viz-table\"><summary>Show the numbers</summary>"
    f"<table><thead>{table_head}</thead><tbody>{table_body}</tbody></table></details>",
    script,
    "</figure>",
]
Path(__file__).with_name("figure.html").write_text("\n".join(parts) + "\n", encoding="utf-8")
