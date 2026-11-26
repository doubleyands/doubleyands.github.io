"""Builds the chart for the Lecture 16 post from experiment.json (part 4).

Bump + box test after five periods (200 cells, nu = 0.6): Fromm, minmod, MC and superbee slopes,
with the exact cell averages as a gray dashed line.
Writes one HTML block with no blank lines to figure.html, so it can sit directly inside the Markdown
post. Same layout as scripts/posts/numerical-analysis/figure.py.

Run: python scripts/posts/numerical-pde/lecture-16/figure.py   (writes figure.html next to this file)
"""

import json
from html import escape
from pathlib import Path

data = json.loads(Path(__file__).with_name("experiment.json").read_text(encoding="utf-8"))
p4 = data["part4"]
xs_data = p4["x"]

COLORS = {
    "minmod": "#2a78d6",
    "MC": "#eb6834",
    "superbee": "#1baf7a",
    "Fromm": "#eda100",
}
KEYS = {"minmod": "minmod", "MC": "mc", "superbee": "superbee", "Fromm": "fromm"}
NAMES = list(COLORS)
SERIES = {name: p4["runs"][KEYS[name]]["p5"]["profile"] for name in NAMES}
EXACT = p4["exact"]
ks = list(range(len(xs_data)))
LABELS = {"superbee": (0.195, 1.04, "end"), "MC": (0.195, 0.93, "end"), "minmod": (0.2, 0.8, "end"),
          "Fromm": (0.77, 1.12, "start")}

W, H = 680, 400
LEFT, RIGHT, TOP, BOTTOM = 56, 560, 30, 340
Y_TOP, Y_BOTTOM = 1.2, -0.2
X_LO, X_HI = 0.05, 0.95
shown = [i for i in ks if X_LO <= xs_data[i] <= X_HI]


def px(x):
    return LEFT + (RIGHT - LEFT) * (x - X_LO) / (X_HI - X_LO)


def py(u):
    return TOP + (BOTTOM - TOP) * (Y_TOP - u) / (Y_TOP - Y_BOTTOM)


def fmt(value):
    return f"{value:.4f}"


svg = []
svg.append(
    f'<svg viewBox="0 0 {W} {H}" role="img" tabindex="0" '
    f'aria-label="The bump and the box after five periods. Minmod clips the bump to 0.83 and smears the box edges. MC keeps the bump at 0.96. Superbee keeps the bump at 0.99 but makes its top flatter, and has the sharpest box edges. These three stay between 0 and 1. Fromm keeps the bump but wiggles near the box edges, up to 1.08 and down to minus 0.08." '
    f'class="viz-plot">'
)
for value in (1.0, 0.75, 0.5, 0.25, 0.0):
    y = py(value)
    cls = "viz-leader" if value == 0.0 else "viz-grid"
    svg.append(f'<line x1="{LEFT}" x2="{RIGHT}" y1="{y:.1f}" y2="{y:.1f}" class="{cls}"/>')
    svg.append(f'<text x="{LEFT - 10}" y="{y + 4:.1f}" text-anchor="end" class="viz-tick">{value:g}</text>')
for x in (0.1, 0.3, 0.5, 0.7, 0.9):
    svg.append(f'<text x="{px(x):.1f}" y="{BOTTOM + 22}" text-anchor="middle" class="viz-tick">{x:g}</text>')
svg.append(f'<text x="{(LEFT + RIGHT) / 2:.1f}" y="{BOTTOM + 48}" text-anchor="middle" class="viz-axis-title">x</text>')
svg.append(f'<text x="{LEFT}" y="{TOP - 16}" class="viz-axis-title">Q after five periods, 200 cells, ν = 0.6</text>')
exact_path = " ".join(f"{'M' if n == 0 else 'L'}{px(xs_data[i]):.1f},{py(EXACT[i]):.1f}" for n, i in enumerate(shown))
order = ["Fromm", "minmod", "MC", "superbee"]
for name in order:
    path = " ".join(f"{'M' if n == 0 else 'L'}{px(xs_data[i]):.1f},{py(SERIES[name][i]):.1f}" for n, i in enumerate(shown))
    svg.append(f'<path d="{path}" fill="none" stroke="{COLORS[name]}" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/>')
svg.append(f'<path d="{exact_path}" fill="none" stroke="#8e8e93" stroke-width="1.5" stroke-dasharray="4 4"/>')
# direct labels next to a feature of each curve
label_pos = LABELS
for name in NAMES:
    x, y, anchor = label_pos[name]
    svg.append(f'<text x="{px(x):.1f}" y="{py(y):.1f}" text-anchor="{anchor}" class="viz-label">{escape(name)}</text>')
svg.append(f'<line x1="0" x2="0" y1="{TOP}" y2="{BOTTOM}" class="viz-crosshair" visibility="hidden"/>')
svg.append("</svg>")

NAMES_TIP = NAMES + ["exact"]
rows = [[fmt(SERIES[name][i]) for name in NAMES] + [fmt(EXACT[i])] for i in shown]
heads = [f"x = {xs_data[i]:.3f}" for i in shown]
ths = [xs_data[i] for i in shown]
COLORS_TIP = [COLORS[n] for n in NAMES] + ["#8e8e93"]
script = (
    "<script>(() => {"
    "const fig = document.currentScript.closest('figure');"
    "const svg = fig.querySelector('svg');"
    "const line = svg.querySelector('.viz-crosshair');"
    "const tip = fig.querySelector('.viz-tooltip');"
    f"const names = {json.dumps(NAMES_TIP)}; const colors = {json.dumps(COLORS_TIP)};"
    f"const rows = {json.dumps(rows)}; const heads = {json.dumps(heads)};"
    f"const xs = {json.dumps([round(px(t), 1) for t in ths])}; const steps = {len(shown) - 1}, width = {W};"
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
) + '<span class="viz-key"><span class="viz-swatch" style="background:#8e8e93"></span>exact</span>'
table_head = "<tr><th>x</th>" + "".join(f"<th>{escape(name)}</th>" for name in NAMES_TIP) + "</tr>"
table_body = "".join(
    f"<tr><td>{xs_data[i]:.3f}</td>" + "".join(f"<td>{v}</td>" for v in rows[n]) + "</tr>"
    for n, i in enumerate(shown) if n % 2 == 0
)
parts = [
    '<figure class="viz">',
    f'<div class="viz-legend">{legend}</div>',
    '<div class="viz-frame">',
    "".join(svg),
    '<div class="viz-tooltip" hidden></div>',
    "</div>",
    "<figcaption>The slope method for q_t + q_x = 0 on a periodic grid with 200 cells and ν = 0.6, after five full periods (1665 steps). "
    "The dashed gray line is the exact answer (the initial cell averages). "
    "Hover over the chart, or focus it and use the arrow keys, to read the values. The table shows every second cell.</figcaption>",
    "<details class=\"viz-table\"><summary>Show the numbers</summary>"
    f"<table><thead>{table_head}</thead><tbody>{table_body}</tbody></table></details>",
    script,
    "</figure>",
]
Path(__file__).with_name("figure.html").write_text("\n".join(parts) + "\n", encoding="utf-8")
