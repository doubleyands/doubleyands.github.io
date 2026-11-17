"""Builds the chart for the Lecture 13 post from experiment.json (part 1).

|g(theta)| measured by one step on each grid mode, for upwind, Lax-Friedrichs and Lax-Wendroff at
Courant number 0.6, with leapfrog (both roots, |g| = 1) for comparison. Writes one HTML block with
no blank lines to figure.html, so it can sit directly inside the Markdown post. Same layout as
scripts/posts/numerical-analysis/figure.py.

Run: python scripts/posts/numerical-pde/lecture-13/figure.py   (writes figure.html next to this file)
"""

import json
import math
from html import escape
from pathlib import Path

data = json.loads(Path(__file__).with_name("experiment.json").read_text(encoding="utf-8"))
p1 = data["part1"]
thetas = p1["schemes"]["upwind"]["theta"]

COLORS = {
    "upwind": "#2a78d6",
    "Lax–Friedrichs": "#eb6834",
    "Lax–Wendroff": "#1baf7a",
    "leapfrog": "#eda100",
}
NAMES = list(COLORS)
SERIES = {name: p1["schemes"][name]["abs_g"] for name in NAMES[:3]}
SERIES["leapfrog"] = p1["leapfrog_abs_g"]
ks = list(range(len(thetas)))

W, H = 680, 400
LEFT, RIGHT, TOP, BOTTOM = 56, 540, 30, 340
Y_TOP, Y_BOTTOM = 1.1, 0.0


def px(t):
    return LEFT + (RIGHT - LEFT) * t / math.pi


def py(u):
    return TOP + (BOTTOM - TOP) * (Y_TOP - u) / (Y_TOP - Y_BOTTOM)


def fmt(value):
    return f"{value:.4f}"


svg = []
svg.append(
    f'<svg viewBox="0 0 {W} {H}" role="img" tabindex="0" '
    f'aria-label="Size of the amplification factor g against the grid frequency theta from 0 to pi, at Courant number 0.6. Leapfrog stays at 1. Lax-Wendroff stays near 1 for small theta and drops to 0.28 at pi. Upwind drops to 0.2 at pi. Lax-Friedrichs drops fastest for small theta, reaches 0.6 at pi/2 and comes back to 1 at pi." '
    f'class="viz-plot">'
)
for value in (1.0, 0.75, 0.5, 0.25, 0.0):
    y = py(value)
    cls = "viz-leader" if value == 0.0 else "viz-grid"
    svg.append(f'<line x1="{LEFT}" x2="{RIGHT}" y1="{y:.1f}" y2="{y:.1f}" class="{cls}"/>')
    svg.append(f'<text x="{LEFT - 10}" y="{y + 4:.1f}" text-anchor="end" class="viz-tick">{value:g}</text>')
for t, lab in ((0, "0"), (math.pi / 4, "π/4"), (math.pi / 2, "π/2"), (3 * math.pi / 4, "3π/4"), (math.pi, "π")):
    svg.append(f'<text x="{px(t):.1f}" y="{BOTTOM + 22}" text-anchor="middle" class="viz-tick">{lab}</text>')
svg.append(f'<text x="{(LEFT + RIGHT) / 2:.1f}" y="{BOTTOM + 48}" text-anchor="middle" class="viz-axis-title">Grid frequency θ = ξh</text>')
svg.append(f'<text x="{LEFT}" y="{TOP - 16}" class="viz-axis-title">|g(θ)|, measured by one step, ν = 0.6</text>')
for name in NAMES:
    points = [(px(t), py(u)) for t, u in zip(thetas, SERIES[name])]
    path = " ".join(f"{'M' if i == 0 else 'L'}{x:.1f},{y:.1f}" for i, (x, y) in enumerate(points))
    svg.append(f'<path d="{path}" fill="none" stroke="{COLORS[name]}" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/>')
# direct labels in the right margin at each curve's value at theta = pi (spread so they do not overlap)
ends = sorted(((SERIES[n][-1], n) for n in NAMES), reverse=True)
ys, last = {}, -1e9
for val, n in ends:
    y = max(py(val) + 4, last + 16)
    ys[n] = y
    last = y
for n in NAMES:
    svg.append(f'<text x="{RIGHT + 10}" y="{ys[n]:.1f}" class="viz-label">{escape(n)}</text>')
svg.append(f'<line x1="0" x2="0" y1="{TOP}" y2="{BOTTOM}" class="viz-crosshair" visibility="hidden"/>')
svg.append("</svg>")

NAMES_TIP = NAMES
rows = [[fmt(SERIES[name][i]) for name in NAMES] for i in ks]
heads = [f"θ = {t / math.pi:.4f}π" for t in thetas]
ths = thetas
script = (
    "<script>(() => {"
    "const fig = document.currentScript.closest('figure');"
    "const svg = fig.querySelector('svg');"
    "const line = svg.querySelector('.viz-crosshair');"
    "const tip = fig.querySelector('.viz-tooltip');"
    f"const names = {json.dumps(NAMES_TIP)}; const colors = {json.dumps([COLORS[n] for n in NAMES])};"
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

legend = "".join(
    f'<span class="viz-key"><span class="viz-swatch" style="background:{COLORS[name]}"></span>{escape(name)}</span>'
    for name in NAMES
)
table_head = "<tr><th>θ / π</th>" + "".join(f"<th>{escape(name)}</th>" for name in NAMES_TIP) + "</tr>"
table_body = "".join(
    f"<tr><td>{thetas[i] / math.pi:.4f}</td>" + "".join(f"<td>{v}</td>" for v in rows[i]) + "</tr>" for i in ks if i % 2 == 0
)
parts = [
    '<figure class="viz">',
    f'<div class="viz-legend">{legend}</div>',
    '<div class="viz-frame">',
    "".join(svg),
    '<div class="viz-tooltip" hidden></div>',
    "</div>",
    "<figcaption>|g(θ)| for u_t + u_x = 0 at Courant number ν = 0.6, measured by one step of each scheme on each grid mode of a periodic grid with 64 points. "
    "Leapfrog is computed from its two roots; both have size 1. Lax–Friedrichs comes back to 1 at θ = π and meets the leapfrog line there. "
    "Hover over the chart, or focus it and use the arrow keys, to read the values. The table shows every second mode.</figcaption>",
    "<details class=\"viz-table\"><summary>Show the numbers</summary>"
    f"<table><thead>{table_head}</thead><tbody>{table_body}</tbody></table></details>",
    script,
    "</figure>",
]
Path(__file__).with_name("figure.html").write_text("\n".join(parts) + "\n", encoding="utf-8")
