---
title: "PINNs, Problem 4: Advection–Diffusion with a Boundary Layer"
description: "Problem 4 of 30 in my PINN series: -0.01 u'' + u' = 1 with a thin layer at x = 1. All 15 runs reach errors between 1e-5 and 1.3e-4, and NTK is best. The error is not in the layer: it is a small wave over the smooth part."
pubDate: 2026-11-13
project: PINN
tags: [PINN, physics-informed neural networks, loss weighting, NTK, RBA, BRDR, boundary layer, advection-diffusion]
---

This is Problem 4 of my PINN series. The setup is the same as in [Problem 1](/blog/2026-11-06-pinn-01-decay/): Fourier features, an MLP, Adam and then L-BFGS, and three ways to weight the loss (NTK, RBA, BRDR), each with 5 seeds. That post explains all of it. Here I only write what is new.

## 1. The problem

A steady advection–diffusion equation in 1D:

$$
-\varepsilon\,u''(x) + u'(x) = 1, \qquad x \in (0, 1),
$$

$$
u(0) = u(1) = 0, \qquad \varepsilon = 0.01.
$$

The exact solution is

$$
u(x) = x - \frac{e^{(x-1)/\varepsilon} - e^{-1/\varepsilon}}{1 - e^{-1/\varepsilon}}.
$$

{{fig-solution}}

Why is this a test for a PINN?

- **A thin layer.** The diffusion term has a small factor $\varepsilon$. Away from $x = 1$, it hardly matters, and the equation is almost $u' = 1$. So the solution is almost the line $u = x$. But the line ends at $u(1) = 1$, and the boundary condition says $u(1) = 0$. The solution drops from about 1 to 0 in a layer of width about $\varepsilon = 0.01$.
- **Few points in the layer.** The 1,024 residual points are random and uniform. For seed 0, only 47 of them are in $x > 0.97$.
- **Two very different regions.** In the layer, $|u''|$ goes up to $1/\varepsilon^2 = 10{,}000$. Outside the layer, $u''$ is almost 0. The network must be smooth on most of the interval and very steep at one end.
- **Information flows from the left.** Outside the layer, $u' \approx 1$, so $u$ is fixed by its value at $x = 0$. This will matter in Section 4.

If a network ignores the layer and learns $u = x$, its relative L2 error is $0.13$. So any error well below $0.1$ means the network has found the layer.

## 2. What is different from Problem 1

- 1,024 random residual points in $(0, 1)$ (a new set for each seed). The boundary term has two points, $x = 0$ and $x = 1$.
- The test grid has 1,001 evenly spaced points in $[0, 1]$.
- Everything else (network, Adam 20,000 steps, L-BFGS 2,000 iterations, float64, seeds 0 to 4) is the same.

## 3. Results

| weighting | final error, mean ± std | min – max | after Adam only, median (min – max) |
|---|---|---|---|
| NTK | (1.44 ± 0.30)×10⁻⁵ | 9.9×10⁻⁶ – 1.78×10⁻⁵ | 1.5×10⁻³ (1.2×10⁻³ – 6.3×10⁻³) |
| RBA | (7.0 ± 4.1)×10⁻⁵ | 2.16×10⁻⁵ – 1.26×10⁻⁴ | 7.6×10⁻³ (2.7×10⁻³ – 3.2×10⁻²) |
| BRDR | (4.1 ± 1.7)×10⁻⁵ | 2.85×10⁻⁵ – 7.1×10⁻⁵ | 2.1×10⁻⁴ (1.1×10⁻⁴ – 3.9×10⁻³) |

All 15 runs finished without problems. **NTK is the best here.** Its worst seed ($1.78 \times 10^{-5}$) is better than the best RBA seed ($2.16 \times 10^{-5}$) and the best BRDR seed ($2.85 \times 10^{-5}$). So this difference is larger than the seed spread. RBA and BRDR overlap, so I cannot rank them.

{{fig-training}}

The training curves have three parts:

- **A plateau.** All three methods start with an error near 0.9 and stay there for a while. At step 1,000, the boundary loss is still 0.03 to 0.13 in all 15 runs, so the boundary values are off by about 0.2 to 0.35. I did not save the network output at these steps, so I cannot show its shape.
- **A sudden drop.** The plateau ends at different times: steps 2,200 to 3,900 for NTK, 1,500 to 7,300 for RBA, and 2,200 to 7,400 for BRDR (the first log with error below 0.1). Then the error falls fast. In 14 of 15 runs, it goes from above 0.1 to below 0.01 within 100 steps.
- **Adam, then L-BFGS.** The BRDR median jumps back up to about 0.3 near step 8,500 and then comes down again. From step 9,600 to the end of Adam, BRDR has the lowest median error. At the end of Adam, its median is $2.1 \times 10^{-4}$, about 7 times lower than NTK. But L-BFGS cuts the NTK and RBA medians by about 100 times, and the BRDR median by only about 7 times. So NTK ends first.

## 4. The error is not in the layer

The layer is the hard part of this problem. So I expected the final error to sit in the layer. It does not. The next chart shows the error at each test point for seed 0.

{{fig-error}}

For all three methods, less than 1% of the squared error is in $x \ge 0.95$ (0.5% for NTK and RBA, 0.3% for BRDR). The largest errors are at $x = 0.30$ (NTK), $x = 0.68$ (RBA) and $x = 0.16$ (BRDR). The networks did learn the layer. The remaining error is in the smooth part.

In the smooth part, the error is a wave. For seed 0, it changes sign 10 to 12 times in $x < 0.95$ for each method. It is not a constant shift: on $x < 0.9$, the mean of the signed error is small next to its RMS (for BRDR, $-3.1 \times 10^{-6}$ against $1.7 \times 10^{-5}$).

There is one more thing about the boundary. Outside the layer, $u' \approx 1$, so a wrong value at $x = 0$ would move the whole line up or down. To be exact: if $u(0)$ is off by $\delta$, the solution changes by $\delta\,(1 - e^{(x-1)/\varepsilon})/(1 - e^{-1/\varepsilon})$, which is almost $\delta$ everywhere except in the layer. So the left boundary matters here. How well does each method fit it? Here is the RMS of the two boundary residuals at the end, next to the RMS error on the test grid (the relative error times the RMS of $u$, 0.564), over the 5 seeds:

| weighting | boundary residual, RMS | error, RMS |
|---|---|---|
| NTK | 3×10⁻⁸ – 8×10⁻⁷ | 6×10⁻⁶ – 1.0×10⁻⁵ |
| RBA | 3×10⁻⁸ – 1.4×10⁻⁷ | 1.2×10⁻⁵ – 7.1×10⁻⁵ |
| BRDR | 2.8×10⁻⁶ – 4.1×10⁻⁵ | 1.6×10⁻⁵ – 4.0×10⁻⁵ |

For NTK and RBA, the boundary residual is at least 10 times smaller than the error in every seed. For BRDR, in 4 of 5 seeds, the two are about the same size (within a factor of 1.7). For seed 0, $u_\theta(0) = -1.7 \times 10^{-5}$ for BRDR, but at most $3 \times 10^{-7}$ in size for NTK and RBA. So BRDR fits the boundary no better than the rest of the solution. But its error is still mostly a wave, not a shift, so I cannot say how much of the BRDR error comes from the boundary.

The weights explain why the BRDR boundary fit is loose. At the end of training, the mean weight of the BRDR boundary points is 0.015 to 0.17 times the mean weight of the equation points (5 seeds). NTK gives the boundary term 1.2 to 3.2 times the weight of the equation term. RBA gives the boundary points 6.4 to 7.9 times the mean equation weight, and RBA squares its weights in the loss, so the boundary counts even more. L-BFGS keeps the Adam weights fixed, so BRDR's L-BFGS hardly pushes on the boundary. I saw the same small boundary weight for BRDR in [Problem 2](/blog/2026-11-09-pinn-02-poisson-1d/).

## 5. Where RBA and BRDR put their weight

RBA and BRDR give each equation point its own weight. Do they find the layer? The next chart shows the weights at the end of training (seed 0), divided by their mean.

{{fig-weights}}

- **RBA** puts its largest weights in the layer. The mean weight is about 0.9 for $x < 0.9$. It rises to 1.9 for $0.95 \le x < 0.98$, 2.7 for $0.98 \le x < 0.99$, and 4.2 for $x \ge 0.99$.
- **BRDR** puts its largest weights just before the layer. The mean weight is 4.2 for $0.9 \le x < 0.95$ and 5.8 for $0.95 \le x < 0.98$, but only 2.1 for $x \ge 0.99$. On the left half ($x < 0.5$) it is 0.37.

So both methods react to the layer, in different ways. But neither one gave a lower final error than NTK, which has one weight for all equation points.

## 6. Summary

- All three weightings solve the boundary layer problem. The final errors are between $1 \times 10^{-5}$ and $1.3 \times 10^{-4}$.
- NTK has the lowest final error, $(1.44 \pm 0.30) \times 10^{-5}$. Its range does not overlap with RBA or BRDR.
- All runs first sit on a plateau near error 0.9 for 1,500 to 7,400 steps. Then the error drops by about two orders of magnitude within about 100 steps.
- BRDR is best at the end of Adam, but L-BFGS helps it least. BRDR gives the boundary points a small weight, and its boundary residual is as large as its error.
- The final error is not in the layer. It is a small wave over the smooth part. RBA and BRDR both put more weight near the layer.

**Problem 5** (on Monday) is the 1D heat equation.

The code is in [scripts/posts/pinn/04-boundary-layer](https://github.com/doubleyands/doubleyands.github.io/tree/main/scripts/posts/pinn/04-boundary-layer): `problem.py` defines this problem, and `figure.py` makes the numbers and charts from the saved results. The shared code is [`pinn.py`](https://github.com/doubleyands/doubleyands.github.io/tree/main/scripts/posts/pinn).
