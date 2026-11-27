---
title: "PINNs, Problem 10: The Laplace Equation on an L-Shaped Domain"
description: "Problem 10 of 30 in my PINN series: the Laplace equation on an L-shape, where the gradient of the solution blows up at one corner. NTK reaches 6.5e-4 and BRDR 1.1e-3, close to each other; RBA is worse at 4.9e-3. For NTK, two thirds of the squared error sits within 0.1 of the corner."
pubDate: 2026-11-27
project: PINN
tags: [PINN, physics-informed neural networks, loss weighting, NTK, RBA, BRDR, Laplace equation, corner singularity]
---

This is Problem 10 of my PINN series. The setup is the same as in [Problem 1](/blog/2026-11-06-pinn-01-decay/): Fourier features, an MLP, Adam and then L-BFGS, and three ways to weight the loss (NTK, RBA, BRDR), each with 5 seeds. That post explains all of it. Here I only write what is new. The color maps work as in [Problem 9](/blog/2026-11-25-pinn-09-poisson-2d/).

Problem 9 had a smooth solution on a square. This time the domain has a corner that points inward, and the solution is not smooth there.

## 1. The problem

The domain is the square $[-1, 1]^2$ without its bottom-right quarter. So it looks like the letter L. The inner corner is at $(0, 0)$. On this domain, I solve the Laplace equation

$$
u_{xx} + u_{yy} = 0.
$$

In polar coordinates $(r, \theta)$ around the corner, with $\theta$ from $0$ to $3\pi/2$, the exact solution is

$$
u = r^{2/3} \sin\!\left(\tfrac{2\theta}{3}\right).
$$

On the two edges that meet at the corner ($\theta = 0$ and $\theta = 3\pi/2$), $u = 0$. On the other four edges, I give the exact values. This is a classic test problem for numerical methods.

{{fig-solution}}

Why is this hard?

- **The gradient blows up at the corner.** The size of $\nabla u$ grows like $r^{-1/3}$, and the second derivatives grow like $r^{-4/3}$. They are infinite at $r = 0$. A neural network is smooth, so it cannot follow this shape exactly. It has to round off the corner.
- **The hard part is very small.** The trouble is in a tiny region. Only 0.76% of the test points are within 0.1 of the corner. So the error near the corner can be large while the rest is good.

The two loss terms have similar sizes here. At the first step, the mean squared equation residual is 2.5 to 7.7, and the edge one is 0.18 to 0.62. This is very different from Problem 9, where the gap was about 1,000 times.

## 2. What is different from Problem 1

- The equation term uses 6,144 random points in the L-shape. The edge term uses about 500 random points on the six edges (504 to 531, depending on the seed). Each seed has its own points.
- The test points are the points of a $201 \times 201$ grid on $[-1, 1]^2$ that lie in the L-shape: 30,301 points. This includes the corner itself.
- Everything else (network, Adam 20,000 steps, L-BFGS 2,000 iterations, float64, seeds 0 to 4) is the same. I did not change anything for this problem.

## 3. Results

| weighting | final error, mean ± std | min – max | after Adam only, median (min – max) |
|---|---|---|---|
| NTK | (6.5 ± 2.1)×10⁻⁴ | 3.35×10⁻⁴ – 8.37×10⁻⁴ | 3.7×10⁻³ (1.7×10⁻³ – 1.0×10⁻²) |
| RBA | (4.9 ± 2.1)×10⁻³ | 1.61×10⁻³ – 6.89×10⁻³ | 6.9×10⁻³ (5.7×10⁻³ – 1.1×10⁻²) |
| BRDR | (1.11 ± 0.45)×10⁻³ | 8.21×10⁻⁴ – 1.88×10⁻³ | 5.3×10⁻³ (2.4×10⁻³ – 9.4×10⁻³) |

All 15 runs finished without problems.

- **NTK and BRDR are close.** NTK has the lower mean, but the best BRDR seed ($8.21 \times 10^{-4}$) is just below the worst NTK seed ($8.37 \times 10^{-4}$). So the difference is within the seed spread.
- **RBA is the worst.** Every RBA seed is worse than every NTK seed. The mean is about 7.5 times higher. Its best seed ($1.61 \times 10^{-3}$) overlaps with the worst BRDR seed ($1.88 \times 10^{-3}$).
- **L-BFGS helps NTK and BRDR, not RBA.** It cuts the median error by 5.1 times for NTK and 5.8 times for BRDR, but only by 1.3 times for RBA.
- All three are much worse than in Problem 9. There, NTK reached $2.3 \times 10^{-5}$. Here it is about 29 times higher. The corner is the reason, as the next sections show.

## 4. Training

{{fig-training}}

- **NTK** is the lowest of the three during all of Adam (smoothed median). All seeds are below error 0.01 by step 3,100. But it falls slowly after that, and ends Adam near $3 \times 10^{-3}$.
- **BRDR** follows NTK from above. Its smoothed median is 1.3 to 1.8 times that of NTK for most of Adam. In the L-BFGS part, both fall to below $10^{-3}$ (median).
- **RBA** is the slowest. Its smoothed median is still about $7 \times 10^{-3}$ at the end of Adam, and no RBA seed gets below $10^{-3}$, even with L-BFGS.

## 5. Where the error is

The color maps show the error of seed 0 on the whole L-shape.

{{fig-error}}

For NTK and BRDR, most of the map is light or medium (errors around $10^{-4}$ to $10^{-3}$), and there is one dark spot at the corner $(0, 0)$. The RBA map is dark on almost the whole domain (errors around $10^{-3}$ to $10^{-2}$). The thin light lines are where the error changes sign.

The next chart makes this clear. It shows the RMS error of seed 0 in thin rings around the corner.

{{fig-corner}}

- **At the corner, all three are about the same.** In the first ring ($r < 0.02$), the RMS error is $1.0 \times 10^{-2}$ for NTK, $1.1 \times 10^{-2}$ for RBA and $1.4 \times 10^{-2}$ for BRDR. The largest error is at the corner point itself, where $u = 0$. There, the networks give $0.025$ to $0.031$ instead of $0$. They round off the sharp tip.
- **Away from the corner, NTK and BRDR fall fast.** For NTK, the RMS error falls from $1.0 \times 10^{-2}$ in the first ring to $9.0 \times 10^{-4}$ near $r = 0.1$ and $1.1 \times 10^{-4}$ near $r = 1$. BRDR falls in a similar way.
- **RBA stays high.** Its RMS error is still $4.9 \times 10^{-3}$ near $r = 0.5$, and $2.5 \times 10^{-3}$ near $r = 1$. That is about 20 times higher than NTK.
- For $r > 1.2$, the error goes up again for all three. These rings only cover a small area near the far corner $(-1, 1)$, where two outer edges meet.

Here are the shares of the squared error of seed 0 that sit within 0.1 of the corner (only 0.76% of the test points):

| weighting | share within r < 0.1 | RMS error where r ≥ 0.2 | RMS error on the edges |
|---|---|---|---|
| NTK | 68% | 1.6×10⁻⁴ | 1.4×10⁻³ |
| RBA | 8% | 3.6×10⁻³ | 5.9×10⁻³ |
| BRDR | 47% | 4.7×10⁻⁴ | 2.4×10⁻³ |

So for NTK, two thirds of the total error comes from a tiny disk around the corner. The corner error is about the same for all three methods. What makes RBA worse is the error on the rest of the domain.

## 6. Why RBA is worse away from the corner

The edge residuals point to the same reason as in Problem 9. Here are the final RMS residuals of each term (range over 5 seeds):

| weighting | equation | edges | final error |
|---|---|---|---|
| NTK | 3.2×10⁻³ – 7.7×10⁻³ | 3.5×10⁻⁴ – 1.6×10⁻³ | 3.4×10⁻⁴ – 8.4×10⁻⁴ |
| RBA | 6.0×10⁻³ – 1.2×10⁻² | 2.1×10⁻³ – 6.6×10⁻³ | 1.6×10⁻³ – 6.9×10⁻³ |
| BRDR | 2.3×10⁻³ – 4.4×10⁻³ | 1.0×10⁻³ – 2.8×10⁻³ | 8.2×10⁻⁴ – 1.9×10⁻³ |

For the same seed, the RBA edge residual is 3.9 to 9.4 times larger than that of NTK. A wrong value on the edge spreads into the inside, so the RBA error is large far from the corner too. RBA also has the largest equation residual.

The weights:

- **NTK** gives the edge term 54 to 100 times the equation weight at the end. This is much less than in Problem 9 (about 2,000 times), because the two terms have similar sizes here. But it is still enough to hold the edges well.
- **RBA** gives the edge points only 1.7 to 2.0 times the mean equation weight. Its point weights also move toward the corner. For seed 0, the equation points within 0.1 of the corner get about 4 times the mean weight, and the points far away (r > 0.7) get 0.79 times. The edge points within 0.1 of the corner get about 7 times. RBA grows a weight where the residual stays large, and near the corner the residual cannot go to zero. I think this pulls the training toward a region it cannot fix, but I did not test this.
- **BRDR** gives the edge points 1.8 to 2.7 times the mean equation weight. For seed 0, its equation points within 0.1 of the corner get about 2.8 times the mean weight. But the edge points within 0.1 of the corner get 32 times the mean equation weight. Its error away from the corner is close to NTK.

## 7. Summary

- On the L-shaped Laplace problem, NTK reaches $(6.5 \pm 2.1) \times 10^{-4}$ and BRDR $(1.11 \pm 0.45) \times 10^{-3}$. Their ranges just touch, so the difference is within the seed spread. RBA reaches $(4.9 \pm 2.1) \times 10^{-3}$, and every RBA seed is worse than every NTK seed.
- The corner limits all three. A smooth network cannot follow $r^{2/3}$, so it rounds off the tip. At the corner point, the error is 0.025 to 0.031 for all three.
- For NTK, 68% of the squared error is within 0.1 of the corner. RBA is worse because of the rest of the domain: its edge residual is 4 to 9 times larger than that of NTK.
- The two loss terms have similar sizes here, so NTK needs a smaller edge weight than in Problem 9 (54 to 100 times instead of about 2,000).

**Problem 11** (on Monday) is the 2D heat equation.

The code is in [scripts/posts/pinn/10-l-shape](https://github.com/doubleyands/doubleyands.github.io/tree/main/scripts/posts/pinn/10-l-shape): `problem.py` defines this problem, and `figure.py` makes the numbers and charts from the saved results. The color maps come from `heatmap.py`, which now also hides the points outside the domain, and the shared code is [`pinn.py`](https://github.com/doubleyands/doubleyands.github.io/tree/main/scripts/posts/pinn).
