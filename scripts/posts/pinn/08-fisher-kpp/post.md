---
title: "PINNs, Problem 8: The Fisher–KPP Equation"
description: "Problem 8 of 30 in my PINN series: a front that moves into an empty region. All three weightings fit the equation to about 1e-4, but the error ahead of the front grows like e^t, so the final relative error is 2.5e-3 to 2.4e-2. RBA and NTK are within the seed spread, and BRDR is worse."
pubDate: 2026-11-23
project: PINN
tags: [PINN, physics-informed neural networks, loss weighting, NTK, RBA, BRDR, Fisher-KPP, reaction-diffusion]
---

This is Problem 8 of my PINN series. The setup is the same as in [Problem 1](/blog/2026-11-06-pinn-01-decay/): Fourier features, an MLP, Adam and then L-BFGS, and three ways to weight the loss (NTK, RBA, BRDR), each with 5 seeds. That post explains all of it. Here I only write what is new.

## 1. The problem

The Fisher–KPP equation adds a growth term to the heat equation:

$$
u_t = u_{xx} + u\,(1 - u), \qquad x \in (-10, 20),\; t \in (0, 5].
$$

You can read $u$ as the density of a population. The term $u_{xx}$ spreads it out, and the term $u(1 - u)$ makes it grow until it reaches $u = 1$. So $u = 1$ is a stable state, and $u = 0$ is an unstable one: a tiny amount of $u$ grows.

This equation has a known exact solution. It is a front that keeps its shape and moves right with speed $5/\sqrt{6} \approx 2.04$:

$$
u(x, t) = \Big(1 + e^{\,x/\sqrt{6} \,-\, 5t/6}\Big)^{-2}.
$$

The start value $u(x, 0)$ and the values at the two ends $x = -10$ and $x = 20$ come from this formula.

{{fig-solution}}

Why can this be hard for a PINN?

- **Long time window.** The front moves about 10 in $x$ during the 5 time units. This is a longer time window than in the earlier problems.
- **An unstable state ahead of the front.** Where $u$ is close to 0, a small error grows. For a small error $e$ near $u = 0$, the equation is about $e_t = e_{xx} + e$, so $e$ can grow like $e^t$. Over 5 time units, that is a factor of up to $e^5 \approx 150$. Behind the front, near $u = 1$, the same term is $-e$, so errors die out there.

## 2. What is different from Problem 1

- There are three loss terms: the equation at 8,192 random points in $(-10, 20) \times (0, 5)$, the start value at 512 points on $t = 0$, and the two ends at 256 points (128 on each end). Each seed has its own points.
- The test grid has $301 \times 101$ evenly spaced points in $[-10, 20] \times [0, 5]$.
- Everything else (network, Adam 20,000 steps, L-BFGS 2,000 iterations, float64, seeds 0 to 4) is the same. I did not change anything for this problem. The inputs are scaled to $[-1, 1]$ as in every problem, so the long domain does not change the network.

## 3. Results

| weighting | final error, mean ± std | min – max | after Adam only, median (min – max) |
|---|---|---|---|
| NTK | (6.7 ± 2.5)×10⁻³ | 3.36×10⁻³ – 9.01×10⁻³ | 2.7×10⁻² (1.5×10⁻² – 3.3×10⁻²) |
| RBA | (4.1 ± 1.2)×10⁻³ | 2.51×10⁻³ – 5.39×10⁻³ | 5.5×10⁻³ (2.6×10⁻³ – 7.7×10⁻³) |
| BRDR | (1.65 ± 0.53)×10⁻² | 1.04×10⁻² – 2.36×10⁻² | 1.9×10⁻² (1.1×10⁻² – 2.5×10⁻²) |

All 15 runs finished without problems.

- **RBA has the lowest mean, but the NTK and RBA ranges overlap**, so the difference between those two is within the seed spread.
- **BRDR is worse than both.** Its best seed ($1.04 \times 10^{-2}$) is worse than the worst NTK seed ($9.0 \times 10^{-3}$) and the worst RBA seed ($5.4 \times 10^{-3}$). The gap to NTK is small, but the ranges do not overlap.

L-BFGS helps less than in most earlier problems. It cuts the median error by 3.6 times for NTK, 1.4 times for RBA and 1.1 times for BRDR.

For seed 0, the largest gap to the exact solution on the whole test grid is about 0.03 for NTK, 0.01 for RBA and 0.045 for BRDR. In all three, this point is at the last time $t = 5$, ahead of the front.

## 4. Training

{{fig-training}}

- **RBA** is the fastest. All seeds are below error 0.1 by step 600 and below 0.03 by step 3,700. During Adam, it is the lowest of the three almost all the time.
- **BRDR** gets below 0.1 between steps 800 and 1,800 and below 0.03 between steps 6,300 and 14,100. Four of the five BRDR seeds never get below 0.01, not even with L-BFGS.
- **NTK** is the slowest during Adam. It stays near 3×10⁻² to 6×10⁻² in the second half of Adam, and its raw error jumps a lot from one log to the next (for seed 0, between 0.028 and 0.12 in steps 10,000 to 20,000). L-BFGS then brings it down to $3 \times 10^{-3}$ to $9 \times 10^{-3}$.

The weights explain why NTK behaves differently here than in [Problem 7](/blog/2026-11-20-pinn-07-wave/). In this problem, the three terms have similar sizes, so NTK has little to fix. At the end, the NTK weights of the start term and the end term are only 0.75 to 0.98 times the equation weight. In Problem 7 they were 400 to about 30,000.

## 5. Small residuals, larger error

Here are the final RMS residuals of each term (range over 5 seeds):

| weighting | equation | start value | two ends | final error |
|---|---|---|---|---|
| NTK | 1.6×10⁻⁴ – 1.8×10⁻⁴ | 5.7×10⁻⁵ – 6.8×10⁻⁵ | 3.9×10⁻⁵ – 4.4×10⁻⁵ | 3.4×10⁻³ – 9.0×10⁻³ |
| RBA | 1.6×10⁻⁴ – 2.5×10⁻⁴ | 2.4×10⁻⁵ – 5.8×10⁻⁵ | 1.8×10⁻⁵ – 4.7×10⁻⁵ | 2.5×10⁻³ – 5.4×10⁻³ |
| BRDR | 1.3×10⁻⁴ – 3.1×10⁻⁴ | 1.0×10⁻⁴ – 2.3×10⁻⁴ | 4.2×10⁻⁵ – 1.2×10⁻⁴ | 1.0×10⁻² – 2.4×10⁻² |

All residuals are small, about $10^{-4}$ or less, and the equation residuals of the three methods are about the same. But the error is much larger. The RMS of $u$ on the test grid is 0.60, so the RMS error is 0.60 times the relative error in the table. This RMS error is 9 to 56 times the RMS equation residual of the same run. For comparison, the same ratio was 0.03 to 0.2 in [Problem 5](/blog/2026-11-16-pinn-05-heat/) and 0.007 to 1.3 in Problem 7.

The next chart shows where this gap comes from.

{{fig-error}}

For seed 0, the error at $t = 0$ is only $2.5 \times 10^{-5}$ (RBA) to $1.3 \times 10^{-4}$ (BRDR). Then it grows steadily. From $t = 0.5$ to $t = 5$ it grows by 58 to 77 times. If I fit a straight line to $\ln(\text{error})$ for $t$ from 1 to 5, the slope is 0.89 to 0.93. The pure $e^t$ growth near $u = 0$ has slope 1. So the error grows at almost the rate that the equation allows ahead of the front.

The last chart shows the error at $t = 5$ over $x$, with its sign.

{{fig-front}}

- **Behind the front** (for $x < 5$), the error is tiny: at most $10^{-4}$ for NTK and RBA, and $7.5 \times 10^{-4}$ for BRDR. Near $u = 1$, errors die out.
- **Ahead of the front** (from $x = 8$ to $x = 19$), the error of all three is positive at every test point. So the predicted $u$ is a little too large there, and the predicted front leads a little. The peak is near $x = 13$ to $15$, where the exact $u$ is only 0.02 to 0.06.
- At the right end $x = 20$, the error is near zero again, because the loss fixes the value there.

Over the whole test grid, 83% to 97% of the squared error (seed 0) is at points where $u < 0.1$, and 83% to 84% is at $t \ge 4$.

So in this problem, the PINN solves the equation well everywhere. The error comes from the equation itself: ahead of the front, any small mistake grows over time. The method with the smallest start error also has the smallest final error. For seed 0, the start error is RBA $2.5 \times 10^{-5}$ < NTK $6.7 \times 10^{-5}$ < BRDR $1.3 \times 10^{-4}$, and the final error is in the same order. BRDR gives the start term and the end term only 0.40 to 0.62 and 0.21 to 0.32 times the mean equation weight, and its start residual ($1.0 \times 10^{-4}$ to $2.3 \times 10^{-4}$) is larger than for the other two ($2.4 \times 10^{-5}$ to $6.8 \times 10^{-5}$). I think this is part of why BRDR is worse. But the saved results do not let me split the final error into a part from the start value and a part from the equation residual, so I cannot say how large this part is.

The point weights of RBA and BRDR (seed 0) show no strong pattern. Averaged in bands in $t$, in $x$, or in the distance to the front, the equation weights stay between 0.56 and 1.64 times their mean. The largest BRDR weights are at early times ($t < 1$, 1.64 times the mean). So I do not show a weight chart.

## 6. Summary

- On the Fisher–KPP front, RBA reaches $(4.1 \pm 1.2) \times 10^{-3}$ and NTK $(6.7 \pm 2.5) \times 10^{-3}$. These overlap. BRDR reaches $(1.65 \pm 0.53) \times 10^{-2}$, which is worse than both, and the ranges do not overlap.
- All three fit the equation to about $10^{-4}$, but the RMS error is 9 to 56 times larger. The error is almost all ahead of the front and at late times.
- Ahead of the front, $u = 0$ is unstable, and a small error grows about like $e^t$. The measured growth rate is 0.89 to 0.93, close to 1.
- Here the terms have similar sizes, so the NTK weights stay near 1, and NTK is not better than the point-wise methods.

**Problem 9** (on Wednesday) is the 2D Poisson equation.

The code is in [scripts/posts/pinn/08-fisher-kpp](https://github.com/doubleyands/doubleyands.github.io/tree/main/scripts/posts/pinn/08-fisher-kpp): `problem.py` defines this problem, and `figure.py` makes the numbers and charts from the saved results. The shared code is [`pinn.py`](https://github.com/doubleyands/doubleyands.github.io/tree/main/scripts/posts/pinn).
