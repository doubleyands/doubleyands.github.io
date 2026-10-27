"""Builds the error chart for the post from experiment.json (part 2).

Writes one HTML block (legend, SVG, hover script, data table) with no blank
lines to figure.html, so it can sit directly inside the Markdown post.
Same layout as scripts/posts/numerical-analysis/figure.py.

Run: python scripts/posts/numerical-pde/lecture-7/figure.py   (writes figure.html next to this file)
"""

import json
import math
from html import escape
from pathlib import Path

data = json.loads(Path(__file__).with_name("experiment.json").read_text(encoding="utf-8"))
Ns = data["part2"]["N"]
ks = list(range(len(Ns)))  # position j means k = 2/N[j] = 1/(5 * 2^j)
K_MIN, K_MAX = ks[0], ks[-1]

# Categorical slots 1-4 of the validated reference palette, in fixed order
COLORS = {
    "Forward Euler": "#2a78d6",
    "RK2": "#eb6834",
    "RK4": "#1baf7a",
    "AB4": "#eda100",
}
NAMES = list(COLORS)
SERIES = {name: data["part2"][name] for name in NAMES}
FLOOR = 1e-16

W, H = 680, 400
LEFT, RIGHT, TOP, BOTTOM = 64, 520, 30, 340
Y_TOP, Y_BOTTOM = -1, -16  # log10 of the error


def px(k):
    return LEFT + (RIGHT - LEFT) * (k - K_MIN) / (K_MAX - K_MIN)


def py(err):
    value = min(max(math.log10(max(err, FLOOR)), Y_BOTTOM), Y_TOP)
    return TOP + (BOTTOM - TOP) * (Y_TOP - value) / (Y_TOP - Y_BOTTOM)


def sci(value):
    if value == 0:
        return "0"
    mantissa, exponent = f"{value:.1e}".split("e")
    return f"{mantissa}e{int(exponent)}"


def power_label(base, exponent):
    if exponent == 0:
        return "1"
    return f'{base}<tspan dy="-0.5em" font-size="75%">{"−" if exponent < 0 else ""}{abs(exponent)}</tspan>'


svg = []
svg.append(
    f'<svg viewBox="0 0 {W} {H}" role="img" tabindex="0" '
    f'aria-label="Error at t = 2 for the logistic equation, as the time step gets smaller, on log scales. Forward Euler falls with slope 1, RK2 with slope 2, RK4 and AB4 with slope 4. RK4 and AB4 flatten out near 1e-15, where rounding takes over." '
    f'class="viz-plot">'
)
for exponent in range(-2, Y_BOTTOM - 1, -2):
    y = py(10.0**exponent)
    svg.append(f'<line x1="{LEFT}" x2="{RIGHT}" y1="{y:.1f}" y2="{y:.1f}" class="viz-grid"/>')
    svg.append(f'<text x="{LEFT - 10}" y="{y + 4:.1f}" text-anchor="end" class="viz-tick">{power_label(10, exponent)}</text>')
for j in (0, 3, 6, 9):
    svg.append(f'<text x="{px(j):.1f}" y="{BOTTOM + 22}" text-anchor="middle" class="viz-tick">1/{Ns[j] // 2}</text>')
svg.append(f'<text x="{(LEFT + RIGHT) / 2:.1f}" y="{BOTTOM + 48}" text-anchor="middle" class="viz-axis-title">Time step Δt (smaller to the right)</text>')
svg.append(f'<text x="{LEFT}" y="{TOP - 16}" class="viz-axis-title">Error</text>')

for name in NAMES:
    errs = SERIES[name]
    points = [(px(k), py(e)) for k, e in zip(ks, errs)]
    path = " ".join(f"{'M' if i == 0 else 'L'}{x:.1f},{y:.1f}" for i, (x, y) in enumerate(points))
    color = COLORS[name]
    svg.append(f'<path d="{path}" fill="none" stroke="{color}" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/>')
    for k, e in zip(ks, errs):
        svg.append(f'<circle cx="{px(k):.1f}" cy="{py(e):.1f}" r="3" fill="{color}" stroke="#ffffff" stroke-width="1.5"/>')

# RK4 and AB4 end near the same error, so move their labels apart
LABEL_DY = {"Forward Euler": 4, "RK2": 4, "RK4": 14, "AB4": -6}
# Direct labels at the right end, where the lines are spread out
for name, dy in LABEL_DY.items():
    errs = SERIES[name]
    svg.append(f'<text x="{RIGHT + 10}" y="{py(errs[-1]) + dy:.1f}" class="viz-label">{escape(name)}</text>')

svg.append(f'<line x1="0" x2="0" y1="{TOP}" y2="{BOTTOM}" class="viz-crosshair" visibility="hidden"/>')
svg.append("</svg>")

rows = [[sci(SERIES[name][i]) for name in NAMES] for i in range(len(ks))]
heads = ["Δt = 1/" + str(N // 2) for N in Ns]

legend = "".join(
    f'<span class="viz-key"><span class="viz-swatch" style="background:{COLORS[name]}"></span>{escape(name)}</span>'
    for name in NAMES
)
table_head = "<tr><th>Δt</th>" + "".join(f"<th>{escape(name)}</th>" for name in NAMES) + "</tr>"
table_body = "".join(
    f"<tr><td>1/{N // 2}</td>" + "".join(f"<td>{v}</td>" for v in row) + "</tr>" for N, row in zip(Ns, rows)
)

script = (
    "<script>(() => {"
    "const fig = document.currentScript.closest('figure');"
    "const svg = fig.querySelector('svg');"
    "const line = svg.querySelector('.viz-crosshair');"
    "const tip = fig.querySelector('.viz-tooltip');"
    f"const names = {json.dumps(NAMES)}; const colors = {json.dumps([COLORS[n] for n in NAMES])};"
    f"const rows = {json.dumps(rows)}; const heads = {json.dumps(heads)};"
    f"const left = {LEFT}, right = {RIGHT}, steps = {len(ks) - 1}, width = {W};"
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
    "<figcaption>Error |U<sup>N</sup> − u(2)| on log scales for the logistic equation u′ = u(1 − u), u(0) = 0.2. "
    "AB4 starts with three RK4 steps. Below about 10<sup>−14</sup> the error is rounding. "
    "Hover over the chart, or focus it and use the arrow keys, to read the values.</figcaption>",
    "<details class=\"viz-table\"><summary>Show the numbers</summary>"
    f"<table><thead>{table_head}</thead><tbody>{table_body}</tbody></table></details>",
    script,
    "</figure>",
]
Path(__file__).with_name("figure.html").write_text("\n".join(parts) + "\n", encoding="utf-8")
