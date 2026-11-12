"""Builds the chart for the Lecture 12 post from experiment.json (part 3).

The solution after one full trip around the periodic domain for four schemes, with the exact
solution (the initial data) in gray. Writes one HTML block with no blank lines to figure.html,
so it can sit directly inside the Markdown post. Same layout as
scripts/posts/numerical-analysis/figure.py.

Run: python scripts/posts/numerical-pde/lecture-12/figure.py   (writes figure.html next to this file)
"""

import json
from html import escape
from pathlib import Path

data = json.loads(Path(__file__).with_name("experiment.json").read_text(encoding="utf-8"))
p3 = data["part3"]
xs_data = p3["x"]
exact = p3["exact"]

COLORS = {
    "Lax–Friedrichs": "#2a78d6",
    "upwind": "#eb6834",
    "Lax–Wendroff": "#1baf7a",
    "leapfrog": "#eda100",
}
NAMES = list(COLORS)
SERIES = {name: p3["runs"][name]["U"] for name in NAMES}
ks = list(range(len(xs_data)))

W, H = 680, 400
LEFT, RIGHT, TOP, BOTTOM = 56, 540, 30, 340
Y_TOP, Y_BOTTOM = 1.6, -0.6
GRAY = "#86868b"


def px(x):
    return LEFT + (RIGHT - LEFT) * x


def py(u):
    return TOP + (BOTTOM - TOP) * (Y_TOP - u) / (Y_TOP - Y_BOTTOM)


def fmt(value):
    value = 0.0 if abs(value) < 5e-13 else value
    return f"{value:+.3f}"


svg = []
svg.append(
    f'<svg viewBox="0 0 {W} {H}" role="img" tabindex="0" '
    f'aria-label="Solution of u_t + u_x = 0 after one full trip around a periodic domain, for Lax-Friedrichs, upwind, Lax-Wendroff and leapfrog, with the exact solution in gray. The exact solution is a narrow bump at x = 0.25 and a box from 0.55 to 0.8. Lax-Friedrichs and upwind smear both shapes. Lax-Wendroff and leapfrog keep the bump high but overshoot near the box and leave wiggles." '
    f'class="viz-plot">'
)
for value in (1.5, 1.0, 0.5, 0.0, -0.5):
    y = py(value)
    cls = "viz-leader" if value == 0.0 else "viz-grid"
    svg.append(f'<line x1="{LEFT}" x2="{RIGHT}" y1="{y:.1f}" y2="{y:.1f}" class="{cls}"/>')
    label = f"{value:g}".replace("-", "−")
    svg.append(f'<text x="{LEFT - 10}" y="{y + 4:.1f}" text-anchor="end" class="viz-tick">{label}</text>')
for x in (0, 0.25, 0.5, 0.75, 1.0):
    svg.append(f'<text x="{px(x):.1f}" y="{BOTTOM + 22}" text-anchor="middle" class="viz-tick">{x:g}</text>')
svg.append(f'<text x="{(LEFT + RIGHT) / 2:.1f}" y="{BOTTOM + 48}" text-anchor="middle" class="viz-axis-title">Position x</text>')
svg.append(f'<text x="{LEFT}" y="{TOP - 16}" class="viz-axis-title">Solution u after one trip (t = 1); gray: exact</text>')

# exact solution as a gray step line (the box has sharp edges)
pts = [(px(x), py(u)) for x, u in zip(xs_data, exact)]
path = " ".join(f"{'M' if i == 0 else 'L'}{x:.1f},{y:.1f}" for i, (x, y) in enumerate(pts))
svg.append(f'<path d="{path}" fill="none" stroke="{GRAY}" stroke-width="2" stroke-dasharray="5 4" stroke-linejoin="round"/>')
for name in NAMES:
    points = [(px(x), py(u)) for x, u in zip(xs_data, SERIES[name])]
    path = " ".join(f"{'M' if i == 0 else 'L'}{x:.1f},{y:.1f}" for i, (x, y) in enumerate(points))
    svg.append(f'<path d="{path}" fill="none" stroke="{COLORS[name]}" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/>')

# direct labels in the right margin, at the height of each scheme's largest value near the box
for name in NAMES:
    top = p3["runs"][name]["box_max"]
    svg.append(f'<text x="{RIGHT + 10}" y="{py(top) + 4:.1f}" class="viz-label">{escape(name)}</text>')
svg.append(f'<text x="{px(0.3):.1f}" y="{py(1.0) - 6:.1f}" class="viz-label">exact</text>')

svg.append(f'<line x1="0" x2="0" y1="{TOP}" y2="{BOTTOM}" class="viz-crosshair" visibility="hidden"/>')
svg.append("</svg>")

NAMES_TIP = NAMES + ["exact"]
rows = [[fmt(SERIES[name][i]) for name in NAMES] + [fmt(exact[i])] for i in ks]
heads = [f"x = {x:.3f}" for x in xs_data]
ths = xs_data
script = (
    "<script>(() => {"
    "const fig = document.currentScript.closest('figure');"
    "const svg = fig.querySelector('svg');"
    "const line = svg.querySelector('.viz-crosshair');"
    "const tip = fig.querySelector('.viz-tooltip');"
    f"const names = {json.dumps(NAMES_TIP)}; const colors = {json.dumps([COLORS[n] for n in NAMES] + [GRAY])};"
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
) + f'<span class="viz-key"><span class="viz-swatch" style="background:{GRAY}"></span>exact</span>'
table_head = "<tr><th>x</th>" + "".join(f"<th>{escape(name)}</th>" for name in NAMES_TIP) + "</tr>"
table_body = "".join(
    f"<tr><td>{xs_data[i]:.3f}</td>" + "".join(f"<td>{v}</td>" for v in rows[i]) + "</tr>" for i in ks if i % 4 == 0
)
parts = [
    '<figure class="viz">',
    f'<div class="viz-legend">{legend}</div>',
    '<div class="viz-frame">',
    "".join(svg),
    '<div class="viz-tooltip" hidden></div>',
    "</div>",
    "<figcaption>u_t + u_x = 0 on a periodic grid with 120 points and Courant number ν = 0.5, after 240 steps (one full trip, so the exact solution equals the initial data). "
    "The labels on the right sit at each scheme's highest value near the box. "
    "Hover over the chart, or focus it and use the arrow keys, to read the values. The table shows every fourth grid point.</figcaption>",
    "<details class=\"viz-table\"><summary>Show the numbers</summary>"
    f"<table><thead>{table_head}</thead><tbody>{table_body}</tbody></table></details>",
    script,
    "</figure>",
]
Path(__file__).with_name("figure.html").write_text("\n".join(parts) + "\n", encoding="utf-8")
