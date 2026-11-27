"""Color maps of 2D fields for the PINN posts (standard library only).

heatmap(...) returns one HTML block (no blank lines) with a row of panels, a color bar,
a hover/keyboard readout and a coarse data table. Each panel is a small PNG, made here with zlib.

Color rules for the whole series:
  signed field (solution, prediction)  blue - gray - red, gray at 0, range +-max|reference|
  field >= 0                           blue, light to dark
  error |u - reference|, point weight  orange, light to dark, log scale
No rainbow maps. Colors are mixed in OKLab, so the lightness changes evenly.
A field value of None is a point outside the domain (L-shape, annulus): it is transparent,
it has no readout, and `outline` draws the domain edge instead of the square.
Run this file to print the lightness along each map (the check that it is monotone / symmetric).
"""

import base64
import json
import math
import struct
import zlib
from html import escape


# --- colors ---
def _lin(c):
    return c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4


def _gam(c):
    c = min(max(c, 0.0), 1.0)
    return 12.92 * c if c <= 0.0031308 else 1.055 * c ** (1 / 2.4) - 0.055


def _oklab(hexcolor):
    r, g, b = (_lin(int(hexcolor[i:i + 2], 16) / 255) for i in (1, 3, 5))
    l = (0.4122214708 * r + 0.5363325363 * g + 0.0514459929 * b) ** (1 / 3)
    m = (0.2119034982 * r + 0.6806995451 * g + 0.1073969566 * b) ** (1 / 3)
    s = (0.0883024619 * r + 0.2817188376 * g + 0.6299787005 * b) ** (1 / 3)
    return (0.2104542553 * l + 0.7936177850 * m - 0.0040720468 * s,
            1.9779984951 * l - 2.4285922050 * m + 0.4505937099 * s,
            0.0259040371 * l + 0.7827717662 * m - 0.8086757660 * s)


def _rgb(lab):
    L, a, b = lab
    l = (L + 0.3963377774 * a + 0.2158037573 * b) ** 3
    m = (L - 0.1055613458 * a - 0.0638541728 * b) ** 3
    s = (L - 0.0894841775 * a - 1.2914855480 * b) ** 3
    r = 4.0767416621 * l - 3.3077115913 * m + 0.2309699292 * s
    g = -1.2684380046 * l + 2.6097574011 * m - 0.3413193965 * s
    bb = -0.0041960863 * l - 0.7034186147 * m + 1.7076147010 * s
    return tuple(round(255 * _gam(c)) for c in (r, g, bb))


MAPS = {  # anchors, evenly spaced from 0 to 1
    "diverging": ["#174a8c", "#5d92d6", "#e3e3e8", "#d86b4c", "#86200f"],
    "blue": ["#eef3fb", "#5d92d6", "#174a8c"],
    "orange": ["#fdf0e7", "#eb6834", "#7a2c0c"],
}


def colormap(name, n=256):
    """n RGB tuples along the map."""
    labs = [_oklab(c) for c in MAPS[name]]
    out = []
    for k in range(n):
        t = k / (n - 1) * (len(labs) - 1)
        i = min(int(t), len(labs) - 2)
        f = t - i
        out.append(_rgb(tuple(a + (b - a) * f for a, b in zip(labs[i], labs[i + 1]))))
    return out


def _hex(rgb):
    return "#%02x%02x%02x" % rgb


# --- PNG ---
def png(rows):
    """rows: list of lists of RGB (or RGBA) tuples, top row first. Returns a data URI."""
    h, w = len(rows), len(rows[0])
    kind = 6 if len(rows[0][0]) == 4 else 2  # 6: with alpha
    raw = b"".join(b"\x00" + bytes(c for px in row for c in px) for row in rows)

    def chunk(kind, data):
        return struct.pack(">I", len(data)) + kind + data + struct.pack(">I", zlib.crc32(kind + data) & 0xFFFFFFFF)

    data = (b"\x89PNG\r\n\x1a\n" + chunk(b"IHDR", struct.pack(">IIBBBBB", w, h, 8, kind, 0, 0, 0))
            + chunk(b"IDAT", zlib.compress(raw, 9)) + chunk(b"IEND", b""))
    return "data:image/png;base64," + base64.b64encode(data).decode()


# --- the figure ---
def heatmap(panels, *, xs, ys, cmap, vrange, log=False, bar_ticks, bar_title, aria, caption, fmt,
            x_title="x", y_title="y", hover_every=2, table_every=10, axis_ticks=(0, 0.5, 1), outline=None):
    """panels: [(title, field)], field[i * len(ys) + j] is the value at (xs[i], ys[j]) ('ij' order).
    vrange: (lo, hi) of the color scale; values outside are clipped. log: log10 color scale.
    bar_ticks: [(value, text)] on the color bar.
    None in a field: outside the domain. outline: [(x, y)] corners of the domain edge."""
    nx, ny = len(xs), len(ys)
    colors = colormap(cmap)
    lo, hi = (math.log10(vrange[0]), math.log10(vrange[1])) if log else vrange

    def level(v):
        v = math.log10(max(v, 1e-300)) if log else v
        return min(max((v - lo) / (hi - lo), 0.0), 1.0)

    masked = any(v is None for _, field in panels for v in field)
    n = len(panels)
    W, gap, L, T, S = 680, 18, 52, 34, 0  # panel size S from the width
    S = (W - L - 10 - gap * (n - 1)) / n
    bar_y = T + S + 82  # room for the axis text, which is 20px on a phone
    H = bar_y + 44
    x0s = [L + k * (S + gap) for k in range(n)]
    px = lambda x: (x - xs[0]) / (xs[-1] - xs[0]) * S
    py = lambda y: (ys[-1] - y) / (ys[-1] - ys[0]) * S

    svg = [f'<svg viewBox="0 0 {W} {H:.0f}" role="img" tabindex="0" aria-label="{escape(aria)}" class="viz-plot">']
    for k, (title, field) in enumerate(panels):
        if masked:
            rows = [[(0, 0, 0, 0) if field[i * ny + j] is None else colors[round(255 * level(field[i * ny + j]))] + (255,)
                     for i in range(nx)] for j in reversed(range(ny))]
        else:
            rows = [[colors[round(255 * level(field[i * ny + j]))] for i in range(nx)] for j in reversed(range(ny))]
        x0 = x0s[k]
        svg.append(f'<text x="{x0:.1f}" y="{T - 12}" class="viz-axis-title" style="fill: var(--ink-80); font-weight: 600">{escape(title)}</text>')
        svg.append(f'<image href="{png(rows)}" x="{x0:.1f}" y="{T}" width="{S:.1f}" height="{S:.1f}" preserveAspectRatio="none"/>')
        if outline:
            pts = " ".join(f"{x0 + px(a):.1f},{T + py(b):.1f}" for a, b in outline)
            svg.append(f'<polygon points="{pts}" fill="none" stroke="#c7c7cc" stroke-width="1"/>')
        else:
            svg.append(f'<rect x="{x0:.1f}" y="{T}" width="{S:.1f}" height="{S:.1f}" fill="none" stroke="#c7c7cc" stroke-width="1"/>')
        for v in axis_ticks:  # the end ticks hug the panel, so the ticks of two panels never touch
            anchor = "start" if v == axis_ticks[0] else "end" if v == axis_ticks[-1] else "middle"
            svg.append(f'<text x="{x0 + px(v):.1f}" y="{T + S + 20:.1f}" text-anchor="{anchor}" class="viz-tick">{v:g}</text>')
        if k == 0:
            for v in axis_ticks:
                svg.append(f'<text x="{x0 - 6:.1f}" y="{T + py(v) + 4:.1f}" text-anchor="end" class="viz-tick">{v:g}</text>')
    svg.append(f'<text x="{x0s[0] + S / 2:.1f}" y="{T + S + 44:.1f}" text-anchor="middle" class="viz-axis-title">{escape(x_title)}</text>')
    svg.append(f'<text x="4" y="{T + S / 2 + 4:.1f}" class="viz-axis-title">{escape(y_title)}</text>')
    # color bar
    bx0, bx1 = L, W - 10
    strip = [[colors[round(255 * k / 199)] for k in range(200)]]
    svg.append(f'<image href="{png(strip)}" x="{bx0}" y="{bar_y:.1f}" width="{bx1 - bx0}" height="10" preserveAspectRatio="none"/>')
    for v, text in bar_ticks:
        x = bx0 + (bx1 - bx0) * level(v)
        svg.append(f'<line x1="{x:.1f}" x2="{x:.1f}" y1="{bar_y + 10:.1f}" y2="{bar_y + 14:.1f}" stroke="#86868b"/>')
        svg.append(f'<text x="{x:.1f}" y="{bar_y + 32:.1f}" text-anchor="middle" class="viz-tick">{text}</text>')
    svg.append(f'<text x="{bx0}" y="{bar_y - 8:.1f}" class="viz-axis-title">{escape(bar_title)}</text>')
    svg.append(f'<rect width="{S / (nx - 1) * hover_every:.1f}" height="{S / (ny - 1) * hover_every:.1f}" fill="none" '
               'stroke="#1d1d1f" stroke-width="1.5" class="viz-cursor" visibility="hidden"/>')
    svg.append("</svg>")

    hi_ = list(range(0, nx, hover_every))
    hj = list(range(0, ny, hover_every))
    show_ = lambda v: "" if v is None else fmt(v)
    vals = [[[None if field[i * ny + j] is None else fmt(field[i * ny + j]) for j in hj] for i in hi_]
            for _, field in panels]
    titles = [t for t, _ in panels]
    tj = list(range(0, ny, table_every))[::-1]
    ti = list(range(0, nx, table_every))
    tables = "".join(
        f"<p>{escape(t)}</p><table><thead><tr><th>y \\ x</th>" + "".join(f"<th>{xs[i]:g}</th>" for i in ti)
        + "</tr></thead><tbody>" + "".join(
            f"<tr><td>{ys[j]:g}</td>" + "".join(f"<td>{show_(field[i * ny + j])}</td>" for i in ti) + "</tr>" for j in tj)
        + "</tbody></table>" for t, field in panels)
    script = (
        "<script>(() => {"
        "const fig = document.currentScript.closest('figure');"
        "const svg = fig.querySelector('svg'); const cur = svg.querySelector('.viz-cursor');"
        "const tip = fig.querySelector('.viz-tooltip');"
        f"const vals = {json.dumps(vals, separators=(',', ':'), ensure_ascii=False)}; const titles = {json.dumps(titles, ensure_ascii=False)};"
        f"const gx = {json.dumps([round(xs[i], 4) for i in hi_])}, gy = {json.dumps([round(ys[j], 4) for j in hj])};"
        f"const x0s = {json.dumps([round(v, 1) for v in x0s])}, T = {T}, S = {S:.3f}, width = {W};"
        "const cw = S / (gx.length - 1), ch = S / (gy.length - 1);"
        "let p = -1, a = 0, b = 0;"
        "const show = () => {"
        " a = Math.max(0, Math.min(gx.length - 1, a)); b = Math.max(0, Math.min(gy.length - 1, b));"
        " if (vals[p][a][b] === null) { hide(); return; }"  # outside the domain
        " const cx = x0s[p] + a * cw, cy = T + S - b * ch;"
        " cur.setAttribute('x', cx - cw / 2); cur.setAttribute('y', cy - ch / 2); cur.setAttribute('visibility', 'visible');"
        " tip.replaceChildren();"
        " const head = document.createElement('div'); head.className = 'viz-tooltip-head';"
        f" head.textContent = '{escape(x_title)} = ' + gx[a] + ', {escape(y_title)} = ' + gy[b]; tip.append(head);"
        " titles.forEach((name, k) => {"
        "  const row = document.createElement('div'); row.className = 'viz-tooltip-row';"
        "  const value = document.createElement('strong'); value.textContent = vals[k][a][b];"
        "  const label = document.createElement('span'); label.textContent = name;"
        "  row.append(value, label); tip.append(row);"
        " });"
        " const box = svg.getBoundingClientRect(), frame = fig.getBoundingClientRect(), sc = box.width / width;"
        " tip.hidden = false;"
        " const room = frame.width - tip.offsetWidth - 8;"
        " tip.style.left = Math.max(0, Math.min(room, box.left - frame.left + (cx + cw) * sc + 8)) + 'px';"
        " tip.style.top = (box.top - frame.top + (cy + ch) * sc + 8) + 'px';"
        "};"
        "const hide = () => { cur.setAttribute('visibility', 'hidden'); tip.hidden = true; };"
        "svg.addEventListener('pointermove', (event) => {"
        " const box = svg.getBoundingClientRect(), sc = box.width / width;"
        " const x = (event.clientX - box.left) / sc, y = (event.clientY - box.top) / sc;"
        " const k = x0s.findIndex((v) => x >= v && x <= v + S);"
        " if (k < 0 || y < T || y > T + S) { hide(); return; }"
        " p = k; a = Math.round((x - x0s[k]) / cw); b = Math.round((T + S - y) / ch); show();"
        "});"
        "svg.addEventListener('pointerleave', hide);"
        "svg.addEventListener('focus', () => { if (p < 0) { p = 0; a = Math.floor(gx.length / 2); b = Math.floor(gy.length / 2); } show(); });"
        "svg.addEventListener('blur', hide);"
        "svg.addEventListener('keydown', (event) => {"
        " const move = {ArrowRight: [1, 0], ArrowLeft: [-1, 0], ArrowUp: [0, 1], ArrowDown: [0, -1]}[event.key];"
        " if (move) {"  # past the right or left edge, go on to the next or previous panel
        "  a += move[0]; b += move[1];"
        "  if (a >= gx.length && p < titles.length - 1) { p += 1; a = 0; }"
        "  if (a < 0 && p > 0) { p -= 1; a = gx.length - 1; }"
        "  show(); event.preventDefault(); }"
        " if (event.key === 'Escape') hide();"
        "});"
        "})();</script>"
    )
    parts = ['<figure class="viz">', '<div class="viz-frame">', "".join(svg), '<div class="viz-tooltip" hidden></div>',
             "</div>", f"<figcaption>{escape(caption)}</figcaption>",
             f'<details class="viz-table"><summary>Show the numbers</summary>{tables}</details>', script, "</figure>"]
    return "\n".join(parts) + "\n"


if __name__ == "__main__":  # the check: lightness along each map
    for name in MAPS:
        L = [_oklab(_hex(c))[0] for c in colormap(name, 11)]
        print(f"{name:9s} OKLab L:", " ".join(f"{v:.2f}" for v in L))
