"""Builds the error chart for the post from experiment.json (part 2).

Writes one HTML block (legend, SVG, hover script, data table) with no blank
lines to figure.html, so it can sit directly inside the Markdown post.
Same layout as scripts/posts/numerical-analysis/figure.py. The x axis is log(dt)
with the real spacing, because the time steps are not evenly spaced.
Errors above 10^4 are drawn as hollow points on the top edge.

Run: python scripts/posts/numerical-pde/lecture-8/figure.py   (writes figure.html next to this file)
"""

import json
import math
from html import escape
from pathlib import Path

data = json.loads(Path(__file__).with_name("experiment.json").read_text(encoding="utf-8"))
p2 = data["part2"]
Ns, dts = p2["N"], p2["dt"]

# Categorical slots 1-2 of the validated reference palette, in fixed order
COLORS = {
    "λ = −5": "#2a78d6",
    "λ = −1000": "#eb6834",
}
NAMES = list(COLORS)
SERIES = {"λ = −5": p2["-5"], "λ = −1000": p2["-1000"]}
LIMIT = 2 / 1000

W, H = 680, 400
LEFT, RIGHT, TOP, BOTTOM = 64, 520, 30, 340
Y_TOP, Y_BOTTOM = 4, -8  # log10 of the error
X_LEFT, X_RIGHT = math.log10(0.1), math.log10(T_MIN := 2 / 20480)  # larger dt on the left


def px(dt):
    return LEFT + (RIGHT - LEFT) * (math.log10(dt) - X_LEFT) / (X_RIGHT - X_LEFT)


def py(err):
    value = min(max(math.log10(err), Y_BOTTOM), Y_TOP)
    return TOP + (BOTTOM - TOP) * (Y_TOP - value) / (Y_TOP - Y_BOTTOM)


def sci(value):
    mantissa, exponent = f"{value:.1e}".split("e")
    return f"{mantissa}e{int(exponent)}"


def power_label(base, exponent):
    if exponent == 0:
        return "1"
    return f'{base}<tspan dy="-0.5em" font-size="75%">{"−" if exponent < 0 else ""}{abs(exponent)}</tspan>'


svg = []
svg.append(
    f'<svg viewBox="0 0 {W} {H}" role="img" tabindex="0" '
    f'aria-label="Forward Euler error against the time step on log scales, for lambda = -5 and lambda = -1000. Both fall with slope 1 for small steps. For lambda = -1000 the error is huge for every step larger than 2/1000 and drops suddenly to about 1e-6 at that line." '
    f'class="viz-plot">'
)
for exponent in range(Y_TOP, Y_BOTTOM - 1, -2):
    y = py(10.0**exponent)
    svg.append(f'<line x1="{LEFT}" x2="{RIGHT}" y1="{y:.1f}" y2="{y:.1f}" class="viz-grid"/>')
    svg.append(f'<text x="{LEFT - 10}" y="{y + 4:.1f}" text-anchor="end" class="viz-tick">{power_label(10, exponent)}</text>')
for dt, label in ((0.1, "1/10"), (0.01, "1/100"), (0.001, "1/1000"), (0.0001, "1/10000")):
    if px(dt) <= RIGHT + 1:
        svg.append(f'<text x="{px(dt):.1f}" y="{BOTTOM + 22}" text-anchor="middle" class="viz-tick">{label}</text>')
svg.append(f'<text x="{(LEFT + RIGHT) / 2:.1f}" y="{BOTTOM + 48}" text-anchor="middle" class="viz-axis-title">Time step Δt (smaller to the right)</text>')
svg.append(f'<text x="{LEFT}" y="{TOP - 16}" class="viz-axis-title">Max error on [0, 2]</text>')
# stability limit for lambda = -1000
xl = px(LIMIT)
svg.append(f'<line x1="{xl:.1f}" x2="{xl:.1f}" y1="{TOP}" y2="{BOTTOM}" stroke="#86868b" stroke-width="1.5" stroke-dasharray="5 5"/>')
svg.append(f'<text x="{xl - 6:.1f}" y="{BOTTOM - 10}" text-anchor="end" class="viz-tick">Δt = 2/1000</text>')

for name in NAMES:
    errs = SERIES[name]
    color = COLORS[name]
    points = [(px(dt), py(e)) for dt, e in zip(dts, errs)]
    path = " ".join(f"{'M' if i == 0 else 'L'}{x:.1f},{y:.1f}" for i, (x, y) in enumerate(points))
    svg.append(f'<path d="{path}" fill="none" stroke="{color}" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/>')
    for dt, e in zip(dts, errs):
        if math.log10(e) > Y_TOP:
            svg.append(f'<circle cx="{px(dt):.1f}" cy="{py(e):.1f}" r="3.5" fill="#ffffff" stroke="{color}" stroke-width="2"/>')
        else:
            svg.append(f'<circle cx="{px(dt):.1f}" cy="{py(e):.1f}" r="3" fill="{color}" stroke="#ffffff" stroke-width="1.5"/>')

# Direct labels at the right end
LABEL_DY = {"λ = −5": 4, "λ = −1000": 4}
for name, dy in LABEL_DY.items():
    errs = SERIES[name]
    svg.append(f'<text x="{RIGHT + 10}" y="{py(errs[-1]) + dy:.1f}" class="viz-label">{escape(name)}</text>')

svg.append(f'<line x1="0" x2="0" y1="{TOP}" y2="{BOTTOM}" class="viz-crosshair" visibility="hidden"/>')
svg.append("</svg>")

rows = [[sci(SERIES[name][i]) for name in NAMES] for i in range(len(Ns))]
heads = [f"Δt = 2/{N}" for N in Ns]
xs = [round(px(dt), 1) for dt in dts]

legend = "".join(
    f'<span class="viz-key"><span class="viz-swatch" style="background:{COLORS[name]}"></span>{escape(name)}</span>'
    for name in NAMES
)
table_head = "<tr><th>Δt</th>" + "".join(f"<th>{escape(name)}</th>" for name in NAMES) + "</tr>"
table_body = "".join(
    f"<tr><td>2/{N}</td>" + "".join(f"<td>{v}</td>" for v in row) + "</tr>" for N, row in zip(Ns, rows)
)

script = (
    "<script>(() => {"
    "const fig = document.currentScript.closest('figure');"
    "const svg = fig.querySelector('svg');"
    "const line = svg.querySelector('.viz-crosshair');"
    "const tip = fig.querySelector('.viz-tooltip');"
    f"const names = {json.dumps(NAMES)}; const colors = {json.dumps([COLORS[n] for n in NAMES])};"
    f"const rows = {json.dumps(rows)}; const heads = {json.dumps(heads)};"
    f"const xs = {json.dumps(xs)}; const steps = {len(Ns) - 1}, width = {W};"
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
    " let best = 0; xs.forEach((v, i) => { if (Math.abs(v - x) < Math.abs(xs[best] - x)) best = i; });"
    " show(best);"
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
    "<figcaption>Forward Euler on u′ = λ(u − sin t) + cos t, u(0) = 0, up to t = 2, on log scales. "
    "Hollow points on the top edge are errors larger than 10<sup>4</sup> (up to 10<sup>222</sup>; see the table). "
    "The dashed line is the limit Δt = 2/|λ| for λ = −1000. "
    "Hover over the chart, or focus it and use the arrow keys, to read the values.</figcaption>",
    "<details class=\"viz-table\"><summary>Show the numbers</summary>"
    f"<table><thead>{table_head}</thead><tbody>{table_body}</tbody></table></details>",
    script,
    "</figure>",
]
Path(__file__).with_name("figure.html").write_text("\n".join(parts) + "\n", encoding="utf-8")
