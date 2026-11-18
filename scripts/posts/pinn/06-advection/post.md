---
title: "PINNs, Problem 6: Advection with a High Speed"
description: "Problem 6 of 30 in my PINN series: u_t + 30 u_x = 0 with periodic boundaries, a known hard case for PINNs. With this setup, all 15 runs end with errors between 1e-3 and 1.1e-2. But NTK and BRDR first stay at errors above 0.5 for thousands of steps."
pubDate: 2026-11-18
project: PINN
tags: [PINN, physics-informed neural networks, loss weighting, NTK, RBA, BRDR, advection equation]
---

This is Problem 6 of my PINN series. The setup is the same as in [Problem 1](/blog/2026-11-06-pinn-01-decay/): Fourier features, an MLP, Adam and then L-BFGS, and three ways to weight the loss (NTK, RBA, BRDR), each with 5 seeds. That post explains all of it. Here I only write what is new.

## 1. The problem

The advection equation moves a shape at a fixed speed $\beta$ and does not change it:

$$
u_t + \beta\,u_x = 0, \qquad x \in (0, 2\pi),\; t \in (0, 1], \qquad \beta = 30.
$$

The start is one sine wave, and the domain is periodic in $x$ (what leaves on the right comes back on the left):

$$
u(x, 0) = \sin(x), \qquad u(0, t) = u(2\pi, t).
$$

The exact solution is the start, moved to the right by $\beta t$:

$$
u(x, t) = \sin(x - 30t).
$$

{{fig-solution}}

The speed is high. By $t = 1$, the wave has moved 30, which is about 4.8 times around the domain. So at a fixed $x$, the solution goes up and down about 4.8 times in $t$.

This exact problem is a well-known hard case. Krishnapriyan et al. ([NeurIPS 2021](https://arxiv.org/abs/2109.01050)) showed that a plain PINN fails on it when $\beta$ is large, for example $\beta = 30$. So I expected a failure here.

Why is it hard? Any function of the form $f(x - 30t)$ makes the equation residual zero. The equation alone does not pick the right $f$. Only the initial condition at $t = 0$ does, and that information must travel along the lines $x - 30t = \text{const}$ to all later times. The network has to learn the start and carry it far.

## 2. What is different from Problem 1

- There are three loss terms: the equation at 8,192 random points in $(0, 2\pi) \times (0, 1)$, the initial condition at 256 points on $t = 0$, and the periodic condition at 256 random times $t$, where the term is $u_\theta(0, t) - u_\theta(2\pi, t)$. Each seed has its own points.
- The test grid has $201 \times 101$ evenly spaced points in $[0, 2\pi] \times [0, 1]$.
- Everything else (network, Adam 20,000 steps, L-BFGS 2,000 iterations, float64, seeds 0 to 4) is the same. In particular, I did not change anything for this problem.

## 3. Results

| weighting | final error, mean ± std | min – max | after Adam only, median (min – max) |
|---|---|---|---|
| NTK | (5.8 ± 2.1)×10⁻³ | 2.33×10⁻³ – 7.60×10⁻³ | 2.8×10⁻² (2.4×10⁻² – 3.1×10⁻²) |
| RBA | (1.85 ± 0.62)×10⁻³ | 1.27×10⁻³ – 2.82×10⁻³ | 1.5×10⁻² (6.8×10⁻³ – 2.4×10⁻²) |
| BRDR | (3.9 ± 4.1)×10⁻³ | 1.05×10⁻³ – 1.09×10⁻² | 2.0×10⁻² (1.7×10⁻² – 2.3×10⁻²) |

All 15 runs finished without problems. **It did not fail.** Every run ends with a relative L2 error between $1.05 \times 10^{-3}$ and $1.09 \times 10^{-2}$. For seed 0, the largest gap to the exact solution on the whole test grid is about $5 \times 10^{-3}$ to $6 \times 10^{-3}$, for a wave of height 1.

RBA has the lowest mean, but the three ranges overlap, so the difference between the methods is within the seed spread. BRDR has one bad seed (seed 2, $1.09 \times 10^{-2}$); its other four seeds are between $1.05 \times 10^{-3}$ and $4.6 \times 10^{-3}$.

The errors are much larger than in Problems 1 to 5, where the final errors were around $10^{-5}$ to $10^{-4}$. L-BFGS also helps less here. It cuts the median error by about 4 times for NTK, 9 times for RBA and 12 times for BRDR. In Problem 5, the gain for NTK was about 66 times.

Why does it work here when a plain PINN fails? My setup is different from the plain PINN in many ways: Fourier features, loss weighting, float64, 20,000 Adam steps, and then L-BFGS. I did not test which part matters, so I cannot say.

## 4. The plateau

{{fig-training}}

The training curves show the hard part. For a while, two of the three methods look like the failure in the paper:

- **RBA** leaves the start quickly. Every seed is below error 0.5 by step 1,000 and below 0.1 by step 1,800.
- **BRDR** sits on a plateau. At step 1,000 the error is still 0.86 to 0.89. It gets below 0.5 between steps 3,300 and 3,800.
- **NTK** sits on the plateau longest. At step 1,000 the error is 0.87 to 0.89, and it gets below 0.5 only between steps 4,700 and 6,300.

These ranges do not overlap, so the order RBA, then BRDR, then NTK is larger than the seed spread. On the plateau, the error starts near 0.9 and goes down only slowly. Then it falls fast. (An error of 1 is what a network that gives zero everywhere would get.)

What does the network fit during the plateau? Here are the mean squared residuals at step 1,000 (range over 5 seeds):

| weighting | error | equation | initial condition | periodic condition |
|---|---|---|---|---|
| NTK | 0.87–0.89 | 0.67–0.92 | 1.3×10⁻⁴ – 1.6×10⁻³ | 2.5×10⁻³ – 4.8×10⁻³ |
| RBA | 0.09–0.28 | 2.0×10⁻³ – 7.5×10⁻³ | 4.8×10⁻⁵ – 2.3×10⁻³ | 7.3×10⁻⁴ – 4.7×10⁻³ |
| BRDR | 0.86–0.89 | 0.039–0.059 | 0.013–0.021 | 0.017–0.024 |

- **NTK** fits the start almost exactly, but its equation residual is still large. This fits its weights: the initial-condition term gets about 140 times the weight of the equation term at step 1,000 (seed 0), and 93 to 177 times at the end.
- **BRDR** has a much smaller equation residual than NTK, but it fits the start less well, and its error is just as large.

So in both cases the network is close to the start at $t = 0$, but far from the solution at later times. The saved results do not keep the network during the plateau, so I cannot show its shape. RBA is the only one that gets both the equation and the start small early.

## 5. The error grows in time

{{fig-error}}

For NTK and RBA (seed 0), the error is about $10^{-4}$ at $t = 0$ and grows to about $2 \times 10^{-3}$ at $t = 1$. BRDR is higher at the start ($4.3 \times 10^{-4}$), stays between about $4 \times 10^{-4}$ and $6 \times 10^{-4}$ up to $t = 0.5$, and then also grows to $1.8 \times 10^{-3}$ at $t = 1$.

The equation explains this. Let $e = u_\theta - u$ be the error and $r$ the equation residual of the network. Then $e_t + 30\,e_x = r$, so along each line $x - 30t = \text{const}$,

$$
\begin{aligned}
e(x, t) ={}& e(x - 30t,\, 0) \\
&+ \int_0^t r\big(x - 30(t - s),\, s\big)\,ds.
\end{aligned}
$$

The error at time $t$ is the error in the start plus the residual added up along the way. Nothing removes it. (The periodic condition also adds its error each time the wave goes around.) This is the opposite of the heat equation in [Problem 5](/blog/2026-11-16-pinn-05-heat/), which damps errors over time. So in advection, a small residual everywhere still gives an error that grows with $t$.

The point weights of RBA and BRDR (seed 0) show no clear pattern this time. Averaged in five bands in $t$ or in $x$, the equation weights stay between 0.76 and 1.19 times their mean. So I do not show a weight chart. Unlike in Problem 5, BRDR gives the conditions more weight than the equation points (initial 0.8 to 1.6 times, periodic 1.6 to 2.6 times the mean equation weight, 5 seeds).

## 6. Summary

- Advection with speed 30 is a known failure case for a plain PINN. With this series' setup, it did not fail: all 15 runs end between $1.05 \times 10^{-3}$ and $1.09 \times 10^{-2}$.
- NTK $(5.8 \pm 2.1) \times 10^{-3}$, RBA $(1.85 \pm 0.62) \times 10^{-3}$ and BRDR $(3.9 \pm 4.1) \times 10^{-3}$ overlap, so the final difference is within the seed spread.
- The difference is in the start. RBA leaves the plateau by step 1,000. BRDR needs about 3,500 steps, and NTK about 5,000.
- The error grows with $t$, because advection carries every error forward and never damps it. The final errors are much larger than in the earlier problems.
- I did not test which part of the setup avoids the failure.

**Problem 7** (on Friday) is the 1D wave equation.

The code is in [scripts/posts/pinn/06-advection](https://github.com/doubleyands/doubleyands.github.io/tree/main/scripts/posts/pinn/06-advection): `problem.py` defines this problem, and `figure.py` makes the numbers and charts from the saved results. The shared code is [`pinn.py`](https://github.com/doubleyands/doubleyands.github.io/tree/main/scripts/posts/pinn).
