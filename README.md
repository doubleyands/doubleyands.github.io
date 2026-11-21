# doubleyands.github.io

My research blog: <https://doubleyands.github.io>

I write about numerical analysis, PDEs, and scientific machine learning.

The site uses [Astro](https://astro.build). GitHub Actions builds it and deploys it to GitHub Pages.

## Run locally

You need Node.js 22.12 or later.

```
npm install
npm run dev
```

## Front page animation

The front page shows flow past a cylinder (DFG 2D-2 benchmark, Re = 100). I computed it with a small lattice Boltzmann code in `scripts/flow/`.

```
node scripts/flow/simulate.mjs --fps=60 --until=6
node scripts/flow/validate.mjs
node scripts/flow/render.mjs --from=0 --to=4.5 --frames=271
```

| | This code | Reference |
|---|---|---|
| Strouhal number | 0.301 | 0.295 – 0.305 |
| Max drag coefficient | 3.36 | 3.22 – 3.24 |
| Max lift coefficient | 1.05 | 0.99 – 1.01 |

The grid is 880 × 164, and the values are from t ≥ 8 s. The Strouhal number is in the reference range. The force coefficients are 4–5% high because the cylinder wall is a staircase on the grid.
