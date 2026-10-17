"""Draw one thumbnail per post: a simple line drawing of the post's main idea.

Writes src/thumbnails/<post id>.svg (1200 x 630). The build turns each SVG into
a PNG for the post cards and the link preview (src/lib/thumbnail.ts).
No text in the drawings: it would be too small on a card.

Run from the repository root: python scripts/thumbnails/make.py
"""

import math
from pathlib import Path

W, H = 1200, 630
BG = "#f5f5f7"
INK = "#1d1d1f"
GRAY = "#aeaeb2"
FAINT = "#dcdce1"
BLUE = "#0066cc"
WHITE = "#ffffff"

OUT = Path("src/thumbnails")


# ---------- helpers ----------

def svg(*parts):
    body = "\n".join(parts)
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">\n'
        f'<rect width="{W}" height="{H}" fill="{BG}"/>\n{body}\n</svg>\n'
    )


def polyline(points, color=INK, width=6, dash=None):
    d = "M" + " L".join(f"{x:.1f},{y:.1f}" for x, y in points)
    extra = f' stroke-dasharray="{dash}"' if dash else ""
    return (
        f'<path d="{d}" fill="none" stroke="{color}" stroke-width="{width}" '
        f'stroke-linecap="round" stroke-linejoin="round"{extra}/>'
    )


def line(x1, y1, x2, y2, color=INK, width=6, dash=None):
    return polyline([(x1, y1), (x2, y2)], color, width, dash)


def dot(x, y, r=10, fill=INK, stroke=BG, width=4):
    return f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r}" fill="{fill}" stroke="{stroke}" stroke-width="{width}"/>'


def arrow(x1, y1, x2, y2, color=BLUE, width=5, head=18):
    a = math.atan2(y2 - y1, x2 - x1)
    left = (x2 - head * math.cos(a - 0.45), y2 - head * math.sin(a - 0.45))
    right = (x2 - head * math.cos(a + 0.45), y2 - head * math.sin(a + 0.45))
    return line(x1, y1, x2, y2, color, width) + polyline([left, (x2, y2), right], color, width)


def mapper(x0, x1, y0, y1, left, right, top, bottom):
    """Data box [x0, x1] x [y0, y1] onto the pixel box (y up)."""
    def f(x, y):
        return (
            left + (x - x0) / (x1 - x0) * (right - left),
            bottom - (y - y0) / (y1 - y0) * (bottom - top),
        )
    return f


def sample(fn, a, b, n=400):
    return [(a + (b - a) * k / n, fn(a + (b - a) * k / n)) for k in range(n + 1)]


# ---------- one drawing per post ----------

def rootfinding():
    """Week 1: fixed-point iteration for x = cos x, drawn as a cobweb."""
    P = mapper(0, 1.25, 0, 1.05, 170, 1030, 70, 560)
    parts = [polyline([P(0, 0), P(1.05, 1.05)], GRAY, 4)]
    parts.append(polyline([P(x, y) for x, y in sample(math.cos, 0, 1.25)], INK, 6))
    x = 0.2
    web = [P(x, x)]
    for _ in range(14):
        y = math.cos(x)
        web += [P(x, y), P(y, y)]
        x = y
    parts.append(polyline(web, BLUE, 4))
    root = 0.7390851332
    parts.append(dot(*P(root, root), r=13))
    return svg(*parts)


def outward_normals():
    """Lecture 1: a region, its boundary, and the outward normal vectors."""
    cx, cy, rx, ry = 600, 315, 310, 165

    def r(t):
        s = 1 + 0.10 * math.cos(3 * t) + 0.07 * math.sin(2 * t + 0.6)
        return cx + rx * s * math.cos(t), cy + ry * s * math.sin(t)

    n = 600
    pts = [r(2 * math.pi * k / n) for k in range(n + 1)]
    d = "M" + " L".join(f"{x:.1f},{y:.1f}" for x, y in pts) + " Z"
    parts = [f'<path d="{d}" fill="{WHITE}" stroke="{INK}" stroke-width="6" stroke-linejoin="round"/>']
    for k in range(12):
        t = 2 * math.pi * (k + 0.3) / 12
        x, y = r(t)
        xa, ya = r(t - 1e-4)
        xb, yb = r(t + 1e-4)
        tx, ty = xb - xa, yb - ya
        norm = math.hypot(tx, ty)
        # The curve runs counterclockwise on screen (y down), so (ty, -tx) points out.
        nx, ny = ty / norm, -tx / norm
        parts.append(arrow(x, y, x + 78 * nx, y + 78 * ny))
        parts.append(dot(x, y, r=7, fill=BLUE, stroke=WHITE, width=3))
    return svg(*parts)


def two_point_bvp():
    """Lecture 2: a 1D grid, the solution at the nodes, and the 3-point stencil."""
    P = mapper(0, 1, 0, 1.15, 170, 1030, 70, 470)
    u = lambda x: 0.25 + 0.35 * x + 0.55 * math.sin(math.pi * x)
    m = 10
    nodes = [k / m for k in range(m + 1)]
    base = 545
    parts = [line(150, base, 1050, base, GRAY, 4)]
    for k, x in enumerate(nodes):
        px, py = P(x, u(x))
        parts.append(line(px, py, px, base, FAINT, 3))
        parts.append(line(px, base - 12, px, base + 12, GRAY, 4))
    parts.append(polyline([P(x, y) for x, y in sample(u, 0, 1)], INK, 6))
    stencil = {5, 6, 7}
    for k, x in enumerate(nodes):
        px, py = P(x, u(x))
        if k in stencil:
            parts.append(dot(px, py, r=14, fill=BLUE))
        elif k in (0, m):
            parts.append(dot(px, py, r=12, fill=INK))
        else:
            parts.append(dot(px, py, r=9, fill=WHITE, stroke=INK, width=4))
    xa, xb = P(nodes[5], 0)[0], P(nodes[7], 0)[0]
    parts.append(polyline([(xa, base), (xb, base)], BLUE, 8))
    for k in stencil:
        parts.append(dot(P(nodes[k], 0)[0], base, r=9, fill=BLUE))
    return svg(*parts)


def runge():
    """Week 2: Runge's function and its degree-10 interpolant on equal nodes."""
    f = lambda x: 1 / (1 + 25 * x * x)
    nodes = [-1 + 2 * k / 10 for k in range(11)]

    def p(x):
        total = 0.0
        for i, xi in enumerate(nodes):
            li = 1.0
            for j, xj in enumerate(nodes):
                if j != i:
                    li *= (x - xj) / (xi - xj)
            total += f(xi) * li
        return total

    P = mapper(-1.05, 1.05, -0.45, 2.05, 150, 1050, 60, 570)
    parts = [line(*P(-1.05, 0), *P(1.05, 0), FAINT, 3)]
    parts.append(polyline([P(x, y) for x, y in sample(p, -1, 1, 800)], BLUE, 5))
    parts.append(polyline([P(x, y) for x, y in sample(f, -1, 1)], INK, 6))
    for x in nodes:
        parts.append(dot(*P(x, f(x)), r=10))
    return svg(*parts)


def five_point_stencil():
    """Lecture 3: a 2D grid with boundary nodes and the 5-point stencil."""
    cols, rows, h = 9, 5, 100
    x0, y0 = 200, 115
    parts = []
    for j in range(rows):
        parts.append(line(x0, y0 + j * h, x0 + (cols - 1) * h, y0 + j * h, FAINT, 3))
    for i in range(cols):
        parts.append(line(x0 + i * h, y0, x0 + i * h, y0 + (rows - 1) * h, FAINT, 3))
    ci, cj = 5, 2
    c = (x0 + ci * h, y0 + cj * h)
    for di, dj in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
        parts.append(line(*c, c[0] + di * h, c[1] + dj * h, BLUE, 7))
    stencil = {(ci, cj), (ci - 1, cj), (ci + 1, cj), (ci, cj - 1), (ci, cj + 1)}
    for j in range(rows):
        for i in range(cols):
            x, y = x0 + i * h, y0 + j * h
            if (i, j) == (ci, cj):
                parts.append(dot(x, y, r=20, fill=BLUE))
            elif (i, j) in stencil:
                parts.append(dot(x, y, r=14, fill=BLUE))
            elif i in (0, cols - 1) or j in (0, rows - 1):
                parts.append(dot(x, y, r=10, fill=INK))
            else:
                parts.append(dot(x, y, r=9, fill=WHITE, stroke=GRAY, width=4))
    return svg(*parts)


def boundary_layer():
    """Lecture 4: -eps u'' + u' = 1 has a thin boundary layer at the outflow end."""
    eps = 0.04
    u = lambda x: x - (math.exp((x - 1) / eps) - math.exp(-1 / eps)) / (1 - math.exp(-1 / eps))
    P = mapper(0, 1, 0, 1.0, 190, 1010, 80, 440)
    left, right = P(0, 0)[0], P(1, 0)[0]
    layer = P(1 - 4 * eps, 0)[0]
    parts = [f'<rect x="{layer:.1f}" y="60" width="{right - layer:.1f}" height="400" fill="{FAINT}" opacity="0.6"/>']
    parts.append(line(left, P(0, 0)[1], right, P(0, 0)[1], GRAY, 4))
    parts.append(polyline([P(x, y) for x, y in sample(u, 0, 1, 1000)], INK, 6))
    parts.append(dot(*P(0, 0), r=12))
    parts.append(dot(*P(1, 0), r=12))
    for k in range(5):
        xs = left + 30 + k * 170
        parts.append(arrow(xs, 535, xs + 110, 535, BLUE, 6, 22))
    return svg(*parts)


def cubic_spline():
    """Week 3: a natural cubic spline, one piece per interval, through its knots."""
    xs = [0, 1, 2, 3, 4, 5, 6]
    ys = [0.45, 0.85, 0.25, 0.55, 0.95, 0.3, 0.6]
    n = len(xs) - 1
    hs = [xs[i + 1] - xs[i] for i in range(n)]
    # Natural spline: solve the tridiagonal system for the second derivatives M.
    a = [0.0] * (n + 1)
    b = [1.0] * (n + 1)
    c = [0.0] * (n + 1)
    d = [0.0] * (n + 1)
    for i in range(1, n):
        a[i], b[i], c[i] = hs[i - 1] / 6, (hs[i - 1] + hs[i]) / 3, hs[i] / 6
        d[i] = (ys[i + 1] - ys[i]) / hs[i] - (ys[i] - ys[i - 1]) / hs[i - 1]
    for i in range(1, n + 1):
        w = a[i] / b[i - 1]
        b[i] -= w * c[i - 1]
        d[i] -= w * d[i - 1]
    M = [0.0] * (n + 1)
    M[n] = d[n] / b[n]
    for i in range(n - 1, -1, -1):
        M[i] = (d[i] - c[i] * M[i + 1]) / b[i]

    def s(i, x):
        h = hs[i]
        t0, t1 = xs[i + 1] - x, x - xs[i]
        return (
            M[i] * t0 ** 3 / (6 * h) + M[i + 1] * t1 ** 3 / (6 * h)
            + (ys[i] / h - M[i] * h / 6) * t0 + (ys[i + 1] / h - M[i + 1] * h / 6) * t1
        )

    P = mapper(-0.2, 6.2, 0, 1.1, 150, 1050, 60, 470)
    base = 545
    parts = [line(150, base, 1050, base, GRAY, 4)]
    for i in range(n):
        piece = sample(lambda x, i=i: s(i, x), xs[i], xs[i + 1], 120)
        parts.append(polyline([P(x, y) for x, y in piece], BLUE if i % 2 else INK, 6))
    for x, y in zip(xs, ys):
        px, py = P(x, y)
        parts.append(line(px, py + 18, px, base, FAINT, 3))
        parts.append(line(px, base - 12, px, base + 12, GRAY, 4))
        parts.append(dot(px, py, r=12))
    return svg(*parts)


DRAWINGS = {
    "2026-10-03-numerical-analysis-week-1": rootfinding,
    "2026-10-06-numerical-pde-lecture-1": outward_normals,
    "2026-10-08-numerical-pde-lecture-2": two_point_bvp,
    "2026-10-10-numerical-analysis-week-2": runge,
    "2026-10-13-numerical-pde-lecture-3": five_point_stencil,
    "2026-10-15-numerical-pde-lecture-4": boundary_layer,
    "2026-10-17-numerical-analysis-week-3": cubic_spline,
}


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    for post_id, draw in DRAWINGS.items():
        (OUT / f"{post_id}.svg").write_text(draw(), encoding="utf-8")
        print(f"wrote {OUT / post_id}.svg")
