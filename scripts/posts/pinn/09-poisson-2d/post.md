---
title: "PINNs, Problem 9: The 2D Poisson Equation"
description: "Problem 9 of 30 in my PINN series, and the first one in 2D: the Poisson equation on the unit square with a slow bump and a fast ripple. NTK reaches 2.2e-5, about 28 times better than RBA and BRDR (6.2e-4 and 6.3e-4). For RBA and BRDR, the error sits on the walls."
pubDate: 2026-11-25
project: PINN
tags: [PINN, physics-informed neural networks, loss weighting, NTK, RBA, BRDR, Poisson equation, boundary conditions]
---

This is Problem 9 of my PINN series. The setup is the same as in [Problem 1](/blog/2026-11-06-pinn-01-decay/): Fourier features, an MLP, Adam and then L-BFGS, and three ways to weight the loss (NTK, RBA, BRDR), each with 5 seeds. That post explains all of it. Here I only write what is new.

This is the first problem with two space variables. So it is also the first post with color maps.

## 1. The problem

The 2D Poisson equation on the unit square, with $u = 0$ on the four walls:

$$
-(u_{xx} + u_{yy}) = f, \qquad (x, y) \in (0, 1)^2.
$$

I choose $f$ so that the exact solution is known:

$$
\begin{aligned}
u(x, y) &= \sin \pi x \,\sin \pi y \\
&\quad + 0.1 \sin 4\pi x \,\sin 4\pi y,
\end{aligned}
$$

and $f = -(u_{xx} + u_{yy})$. It is the 2D version of [Problem 2](/blog/2026-11-09-pinn-02-poisson-1d/): a slow bump plus a fast ripple.

{{fig-solution}}

Why can this be hard for a PINN?

- **The fast part.** In $u$, the ripple is only 10% of the bump. But the Laplacian multiplies each part by its frequency squared, summed over the two directions. So in $f$, the ripple part has size $32\pi^2 \cdot 0.1 \approx 31.6$, and the bump part has size $2\pi^2 \approx 19.7$. The ripple is the larger part of $f$.
- **The two terms have very different sizes.** The equation residual has the size of $f$, which goes up to about 49 on the square. The wall residual has the size of $u$, which is at most 1. At the first step, the mean squared equation residual is 280 to 500, and the wall one is 0.09 to 0.28. Problem 2 had the same kind of gap.

## 2. What is different from Problem 1

- There are two loss terms: the equation at 4,096 random points in the square, and the wall condition at 512 points (128 on each wall). Each seed has its own points.
- The test grid has $101 \times 101$ evenly spaced points in $[0, 1]^2$.
- Everything else (network, Adam 20,000 steps, L-BFGS 2,000 iterations, float64, seeds 0 to 4) is the same. I did not change anything for this problem.

## 3. Results

| weighting | final error, mean ± std | min – max | after Adam only, median (min – max) |
|---|---|---|---|
| NTK | (2.25 ± 0.30)×10⁻⁵ | 1.76×10⁻⁵ – 2.46×10⁻⁵ | 2.8×10⁻³ (1.7×10⁻³ – 9.0×10⁻³) |
| RBA | (6.2 ± 2.0)×10⁻⁴ | 4.35×10⁻⁴ – 9.08×10⁻⁴ | 2.1×10⁻³ (1.1×10⁻³ – 9.1×10⁻³) |
| BRDR | (6.3 ± 2.6)×10⁻⁴ | 4.38×10⁻⁴ – 1.06×10⁻³ | 1.0×10⁻³ (9.2×10⁻⁴ – 1.5×10⁻²) |

All 15 runs finished without problems.

- **NTK is the best by far.** Its worst seed ($2.5 \times 10^{-5}$) is about 18 times better than the best RBA or BRDR seed ($4.4 \times 10^{-4}$). The mean is about 28 times lower.
- **RBA and BRDR are within the seed spread** of each other.
- **After Adam only, the order is different.** BRDR has the lowest median and NTK the highest, and all three ranges overlap. Then L-BFGS cuts the median error by 115 times for NTK, but only by 3.6 times for RBA and 1.8 times for BRDR. This is the same story as in Problem 2.

## 4. Training

{{fig-training}}

- **NTK** is the fastest at first. All seeds are below error 0.01 by step 300 and below 0.001 by step 2,600. But the error does not keep falling. In the second half of Adam, the smoothed median is near 3×10⁻³ to 5×10⁻³. The equation residual is still large at the end of Adam (RMS 0.06 to 0.16). L-BFGS then fixes it: the RMS equation residual drops to about $5 \times 10^{-3}$, and the error to about $2 \times 10^{-5}$.
- **RBA and BRDR** are slower at first, but they keep improving during Adam. BRDR is the lowest of the three from about step 10,000 to the end of Adam (smoothed median). But L-BFGS helps them much less. RBA gets below 0.001 only in the L-BFGS part, and BRDR gets there between steps 12,200 and 18,500.

## 5. Where the error is

Here are the final RMS residuals of each term (range over 5 seeds):

| weighting | equation | walls | final error |
|---|---|---|---|
| NTK | 4.9×10⁻³ – 5.4×10⁻³ | 2.4×10⁻⁵ – 3.7×10⁻⁵ | 1.8×10⁻⁵ – 2.5×10⁻⁵ |
| RBA | 3.7×10⁻³ – 4.8×10⁻³ | 5.2×10⁻⁴ – 1.0×10⁻³ | 4.3×10⁻⁴ – 9.1×10⁻⁴ |
| BRDR | 3.3×10⁻³ – 4.4×10⁻³ | 5.5×10⁻⁴ – 1.2×10⁻³ | 4.4×10⁻⁴ – 1.1×10⁻³ |

The equation residuals of the three methods are about the same. RBA and BRDR are even a little lower than NTK. The big difference is at the walls: the wall residual of RBA and BRDR is about 15 to 50 times larger than that of NTK. The final error follows the wall residual, not the equation residual.

The color maps show the error of seed 0 on the whole square.

{{fig-error}}

The NTK map is light: the NTK error is below $10^{-5}$ at 89% of the test points. The RBA and BRDR maps are dark almost everywhere, and the darkest parts touch the walls. The thin light lines are where the error changes sign. In the NTK map, it is hard to see where the error is largest. So the next chart shows the RMS error against the distance to the nearest wall.

{{fig-wall}}

- On the walls (distance 0), the RMS error of seed 0 is $2.9 \times 10^{-5}$ for NTK, $1.1 \times 10^{-3}$ for RBA and $7.1 \times 10^{-4}$ for BRDR. These are almost the same as the wall residuals of the same runs. This is expected: on the wall the exact $u$ is 0, so the error there is the wall residual.
- Moving away from the walls, the error falls. For NTK it is flat at about $5 \times 10^{-6}$ from distance 0.1 to the center. For RBA and BRDR it falls by about 17 to 18 times, to $6.4 \times 10^{-5}$ and $4.0 \times 10^{-5}$ at the center.
- The points within 0.1 of a wall are 37% of the test grid. They hold 80% (NTK), 88% (RBA) and 92% (BRDR) of the squared error.

So the error comes from the walls. A wrong value on the wall spreads into the inside, and the Laplacian smooths it out. The equation term cannot fix this well: a smooth error with a small Laplacian has a small equation residual.

## 6. Why NTK fixes the walls

The weights explain it, as in Problem 2 and [Problem 7](/blog/2026-11-20-pinn-07-wave/).

- **NTK** gives the wall term 1,700 to 2,400 times the equation weight at the end (and 340 to 630 times at step 1,000). This fills the large gap between the sizes of the two terms. At the end of Adam, the NTK equation residual is still large, but L-BFGS (with these weights fixed) brings it down and keeps the wall residual at about $3 \times 10^{-5}$.
- **RBA** gives the wall points only 1.3 to 2.1 times the mean equation weight (the loss squares the weight, so 1.6 to 4.4 times in the loss). RBA divides each residual by the largest residual of the same term, so it cannot see that one term is much larger than the other.
- **BRDR** gives the wall points 4.6 to 5.5 times the mean equation weight. This is more than in [Problem 5](/blog/2026-11-16-pinn-05-heat/) (where it was below 1), but it is still far from the gap between the two terms. BRDR sets the weights by how fast each residual decays, not by its size.

The point weights inside the square (seed 0) show no strong pattern. For RBA, the equation weights averaged in bands of wall distance stay between 0.96 and 1.07 times their mean. For BRDR they are 1.42 times the mean very close to the walls (distance below 0.02) and 0.88 to 1.06 elsewhere. So I do not show a weight map.

## 7. Summary

- On the 2D Poisson equation, NTK reaches $(2.25 \pm 0.30) \times 10^{-5}$. RBA reaches $(6.2 \pm 2.0) \times 10^{-4}$ and BRDR $(6.3 \pm 2.6) \times 10^{-4}$. NTK is about 28 times better, and the ranges do not overlap.
- After Adam alone, the three are close, and NTK is not the best. L-BFGS cuts the NTK error by 115 times, but the others only by 2 to 4 times.
- For RBA and BRDR, the error is set by the wall residual. It is largest on the walls and falls toward the center. The equation residual is about the same for all three.
- NTK gives the wall term about 2,000 times the equation weight, which fills the gap between the sizes of the two terms. RBA and BRDR give it only 1.3 to 5.5 times.

**Problem 10** (on Friday) is the Laplace equation on an L-shaped domain, with a corner.

The code is in [scripts/posts/pinn/09-poisson-2d](https://github.com/doubleyands/doubleyands.github.io/tree/main/scripts/posts/pinn/09-poisson-2d): `problem.py` defines this problem, and `figure.py` makes the numbers and charts from the saved results. The color maps come from `heatmap.py`, and the shared code is [`pinn.py`](https://github.com/doubleyands/doubleyands.github.io/tree/main/scripts/posts/pinn).
