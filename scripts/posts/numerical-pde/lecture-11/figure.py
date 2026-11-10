"""Builds the amplification factor chart for the post from experiment.json (part 2).

Writes one HTML block (legend, SVG, hover script, data table) with no blank
lines to figure.html, so it can sit directly inside the Markdown post.
Same layout as scripts/posts/numerical-analysis/figure.py.
The values are the measured g (one step from a single Fourier mode), not the formulas.

Run: python scripts/posts/numerical-pde/lecture-11/figure.py   (writes figure.html next to this file)
"""

import json
import math
from html import escape
from pathlib import Path

data = json.loads(Path(__file__).with_name("experiment.json").read_text(encoding="utf-8"))
runs = data["part2"]["runs"]
N = data["part2"]["N"]

# Categorical slots 1-4 of the validated reference palette, in fixed order
COLORS = {
    "FTCS, r = 0.3": "#2a78d6",
    "FTCS, r = 0.5": "#eb6834",
    "FTCS, r = 0.6": "#1baf7a",
    "Crank–Nicolson, r = 6": "#eda100",
}
KEYS = {
    "FTCS, r = 0.3": "FTCS r=0.3",
    "FTCS, r = 0.5": "FTCS r=0.5",
    "FTCS, r = 0.6": "FTCS r=0.6",
    "Crank–Nicolson, r = 6": "CN r=6",
}
NAMES = list(COLORS)
ms = [row["m"] for row in runs["FTCS r=0.3"]["rows"]]
ths = [2 * math.pi * m / N for m in ms]
ks = list(range(len(ms)))
SERIES = {name: [row["g_re"] for row in runs[KEYS[name]]["rows"]] for name in NAMES}

W, H = 680, 400
LEFT, RIGHT, TOP, BOTTOM = 64, 500, 30, 340
Y_TOP, Y_BOTTOM = 1.2, -1.6


def px(th):
    return LEFT + (RIGHT - LEFT) * th / math.pi


def py(g):
    return TOP + (BOTTOM - TOP) * (Y_TOP - g) / (Y_TOP - Y_BOTTOM)


def fmt(value):
    value = 0.0 if abs(value) < 5e-13 else value
    return f"{value:+.4f}"


svg = []
svg.append(
    f'<svg viewBox="0 0 {W} {H}" role="img" tabindex="0" '
    f'aria-label="Measured amplification factor g against xi h from 0 to pi, for FTCS with r = 0.3, 0.5 and 0.6 and for Crank-Nicolson with r = 6. All curves start at 1. FTCS with r = 0.6 goes below minus 1 near the sawtooth and ends at minus 1.4. FTCS with r = 0.5 ends exactly at minus 1. Crank-Nicolson stays between minus 1 and 1." '
    f'class="viz-plot">'
)
for value in (1.0, 0.5, 0.0, -0.5, -1.0, -1.5):
    y = py(value)
    cls = "viz-leader" if abs(value) == 1.0 else "viz-grid"
    dash = ' stroke-dasharray="5 4"' if abs(value) == 1.0 else ""
    svg.append(f'<line x1="{LEFT}" x2="{RIGHT}" y1="{y:.1f}" y2="{y:.1f}" class="{cls}"{dash}/>')
    label = f"{value:g}".replace("-", "−")
    svg.append(f'<text x="{LEFT - 10}" y="{y + 4:.1f}" text-anchor="end" class="viz-tick">{label}</text>')
for th, label in ((0, "0"), (math.pi / 4, "π/4"), (math.pi / 2, "π/2"), (3 * math.pi / 4, "3π/4"), (math.pi, "π")):
    svg.append(f'<text x="{px(th):.1f}" y="{BOTTOM + 22}" text-anchor="middle" class="viz-tick">{label}</text>')
svg.append(f'<text x="{(LEFT + RIGHT) / 2:.1f}" y="{BOTTOM + 48}" text-anchor="middle" class="viz-axis-title">Wave number times grid spacing, ξh (π is the sawtooth)</text>')
svg.append(f'<text x="{LEFT}" y="{TOP - 16}" class="viz-axis-title">Amplification factor g per step (dashed lines: ±1)</text>')

for name in NAMES:
    points = [(px(t), py(g)) for t, g in zip(ths, SERIES[name])]
    path = " ".join(f"{'M' if i == 0 else 'L'}{x:.1f},{y:.1f}" for i, (x, y) in enumerate(points))
    color = COLORS[name]
    svg.append(f'<path d="{path}" fill="none" stroke="{color}" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/>')
    for x, y in points:
        svg.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="2.5" fill="{color}"/>')

LABEL_DY = {"FTCS, r = 0.3": 4, "FTCS, r = 0.5": 10, "FTCS, r = 0.6": 4, "Crank–Nicolson, r = 6": -6}
for name, dy in LABEL_DY.items():
    svg.append(f'<text x="{RIGHT + 10}" y="{py(SERIES[name][-1]) + dy:.1f}" class="viz-label">{escape(name)}</text>')

svg.append(f'<line x1="0" x2="0" y1="{TOP}" y2="{BOTTOM}" class="viz-crosshair" visibility="hidden"/>')
svg.append("</svg>")

ts = ms
rows = [[fmt(SERIES[name][i]) for name in NAMES] for i in range(len(ks))]
def frac(m):
    """xi h = 2 pi m / N as a reduced multiple of pi."""
    if m == 0:
        return "0"
    g = math.gcd(m, N // 2)
    a, b = m // g, (N // 2) // g
    top = "π" if a == 1 else f"{a}π"
    return top if b == 1 else f"{top}/{b}"


heads = [f"ξh = {frac(m)}" for m in ms]
# the table shows every second wave number
turn_idx = [i for i in ks if ms[i] % 2 == 0]

legend = "".join(
    f'<span class="viz-key"><span class="viz-swatch" style="background:{COLORS[name]}"></span>{escape(name)}</span>'
    for name in NAMES
)
table_head = "<tr><th>ξh</th>" + "".join(f"<th>{escape(name)}</th>" for name in NAMES) + "</tr>"
table_body = "".join(
    f"<tr><td>{heads[i][5:]}</td>" + "".join(f"<td>{v}</td>" for v in rows[i]) + "</tr>" for i in turn_idx
)

script = (
    "<script>(() => {"
    "const fig = document.currentScript.closest('figure');"
    "const svg = fig.querySelector('svg');"
    "const line = svg.querySelector('.viz-crosshair');"
    "const tip = fig.querySelector('.viz-tooltip');"
    f"const names = {json.dumps(NAMES)}; const colors = {json.dumps([COLORS[n] for n in NAMES])};"
    f"const rows = {json.dumps(rows)}; const heads = {json.dumps(heads)};"
    f"const xs = {json.dumps([round(px(t), 1) for t in ths])}; const steps = {len(ks) - 1}, width = {W};"
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
    "<figcaption>The factor g by which one step multiplies the mode e^(ijhξ), for u_t = u_xx on a periodic grid with 32 points. "
    "Each dot is measured from one step of the scheme. All four curves are real. Only FTCS with r = 0.6 leaves the band between −1 and 1. "
    "Hover over the chart, or focus it and use the arrow keys, to read the values.</figcaption>",
    "<details class=\"viz-table\"><summary>Show the numbers</summary>"
    f"<table><thead>{table_head}</thead><tbody>{table_body}</tbody></table></details>",
    script,
    "</figure>",
]
Path(__file__).with_name("figure.html").write_text("\n".join(parts) + "\n", encoding="utf-8")
