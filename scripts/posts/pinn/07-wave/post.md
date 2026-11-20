---
title: "PINNs, Problem 7: The 1D Wave Equation"
description: "Problem 7 of 30 in my PINN series: a string with fixed ends and two modes, one slow and one fast. NTK reaches an error of about 1.7e-3. RBA and BRDR stop at about 4e-2 to 7e-2, because they fit the initial and boundary conditions loosely and the wave carries that error forward."
pubDate: 2026-11-20
project: PINN
tags: [PINN, physics-informed neural networks, loss weighting, NTK, RBA, BRDR, wave equation]
---

This is Problem 7 of my PINN series. The setup is the same as in [Problem 1](/blog/2026-11-06-pinn-01-decay/): Fourier features, an MLP, Adam and then L-BFGS, and three ways to weight the loss (NTK, RBA, BRDR), each with 5 seeds. That post explains all of it. Here I only write what is new.

## 1. The problem

The wave equation describes a string that moves up and down. The ends of the string are fixed:

$$
u_{tt} = c^2\,u_{xx}, \qquad x \in (0, 1),\; t \in (0, 1], \qquad c = 2,
$$

$$
u(0, t) = u(1, t) = 0.
$$

At the start, the string has a given shape and does not move:

$$
u(x, 0) = \sin(\pi x) + 0.5\,\sin(4\pi x), \qquad u_t(x, 0) = 0.
$$

Each sine shape swings on its own, so the exact solution is

$$
\begin{aligned}
u(x, t) ={}& \sin(\pi x)\cos(2\pi t) \\
&+ 0.5\,\sin(4\pi x)\cos(8\pi t).
\end{aligned}
$$

{{fig-solution}}

So the solution has a **slow part** (one arch, period 1 in $t$) and a **fast part** (four small bumps, period 0.25 in $t$). In the time window, the fast part swings four times.

Why can this be hard for a PINN?

- **Two conditions at $t = 0$.** The equation has a second derivative in $t$, so it needs both the start shape $u(x, 0)$ and the start speed $u_t(x, 0)$. These are two loss terms, plus the boundary term.
- **Errors do not fade.** A wave keeps its energy. Like advection in [Problem 6](/blog/2026-11-18-pinn-06-advection/), and unlike heat in [Problem 5](/blog/2026-11-16-pinn-05-heat/), an error in the start or at the ends travels on and does not die out.
- **The terms have very different sizes.** In the fast part, $u_{tt}$ is as large as $0.5\,(8\pi)^2 \approx 316$, while $u$ itself is about 1. So the terms in the equation are hundreds of times larger than the values in the conditions. At the first step, the mean squared equation residual is 130 to 740, and the condition ones are 0.06 to 2.5. [Problem 2](/blog/2026-11-09-pinn-02-poisson-1d/) had the same kind of gap.

## 2. What is different from Problem 1

- There are four loss terms: the equation at 8,192 random points in $(0, 1) \times (0, 1)$, the start shape at 256 points on $t = 0$, the start speed at 256 other points on $t = 0$, and the two fixed ends at 256 points (128 on each end). Each seed has its own points.
- The test grid has $101 \times 101$ evenly spaced points in $[0, 1] \times [0, 1]$.
- Everything else (network, Adam 20,000 steps, L-BFGS 2,000 iterations, float64, seeds 0 to 4) is the same. I did not change anything for this problem.

## 3. Results

| weighting | final error, mean ± std | min – max | after Adam only, median (min – max) |
|---|---|---|---|
| NTK | (1.66 ± 0.29)×10⁻³ | 1.21×10⁻³ – 2.00×10⁻³ | 1.9×10⁻² (1.0×10⁻² – 2.8×10⁻²) |
| RBA | (5.8 ± 1.4)×10⁻² | 4.24×10⁻² – 7.21×10⁻² | 0.19 (0.16 – 0.20) |
| BRDR | (4.23 ± 0.38)×10⁻² | 3.62×10⁻² – 4.55×10⁻² | 0.15 (0.12 – 0.16) |

All 15 runs finished without problems. **NTK is clearly the best here.** Its worst seed ($2.0 \times 10^{-3}$) is about 18 times better than the best seed of RBA or BRDR ($3.6 \times 10^{-2}$). The ranges do not overlap, so this is much larger than the seed spread. RBA and BRDR overlap, so the difference between those two is within the seed spread.

For seed 0, the largest gap to the exact solution on the whole test grid is about $9 \times 10^{-3}$ for NTK and about 0.09 to 0.10 for RBA and BRDR. For a string of height about 1.4, the RBA and BRDR answers are visibly wrong in places.

L-BFGS helps all three, but not by as much as in the easy problems. It cuts the median error by about 11 times for NTK, and by about 3.5 times for RBA and BRDR.

## 4. Training

{{fig-training}}

- **NTK** goes down from the start. All seeds are below error 0.5 by step 300 and below 0.1 between steps 4,300 and 5,900.
- **RBA and BRDR** first fall to about 0.4 and then move down slowly. At step 5,000, the error is 0.37 to 0.39 for RBA and 0.40 to 0.43 for BRDR. No run of these two gets below 0.1 during Adam. They get below 0.1 only in the L-BFGS part.

A side note on the number 0.4. If a network had the slow part exactly right and the fast part missing, its relative L2 error would be 0.447. So the slow phase of RBA and BRDR looks like a network that has not yet learned the fast part. But the saved results do not keep the network during training, so I cannot check this.

## 5. Why NTK wins: the conditions

Here are the final RMS residuals of each term (range over 5 seeds):

| weighting | equation | start shape | start speed | fixed ends |
|---|---|---|---|---|
| NTK | 0.090 – 0.12 | 1.1×10⁻⁴ – 2.0×10⁻⁴ | 4.4×10⁻⁴ – 9.0×10⁻⁴ | 2.8×10⁻⁴ – 5.5×10⁻⁴ |
| RBA | 0.030 – 0.036 | 8.1×10⁻³ – 1.3×10⁻² | 2.1×10⁻³ – 4.3×10⁻³ | 1.2×10⁻² – 2.7×10⁻² |
| BRDR | 0.037 – 0.045 | 7.4×10⁻³ – 1.1×10⁻² | 4.9×10⁻³ – 1.1×10⁻² | 9.8×10⁻³ – 1.3×10⁻² |

This table is the main result of the post.

- RBA and BRDR have a **smaller** equation residual than NTK, about 2 to 4 times smaller.
- But their residuals for the start shape and the fixed ends are **20 to 120 times larger**: about $10^{-2}$, against $10^{-4}$ to $5 \times 10^{-4}$ for NTK. The start speed is 2 to 25 times larger.
- And their mean final error is 25 to 35 times larger.

So in this problem, the error follows the conditions, not the equation. This fits the wave equation. If the start or the ends are wrong by about $10^{-2}$, the solution is wrong by about $10^{-2}$ from the beginning, and the wave carries that error to later times. An equation residual of 0.1 sounds large, but it is small next to the size of $u_{tt}$ (up to about 316 in the fast part).

Why do the conditions get so little attention? Look at the weights at the end, compared with the mean weight of the equation points:

- **NTK** gives the start shape 17,000 to 31,000 times the weight of the equation term, the fixed ends 16,000 to 26,000 times, and the start speed 400 to 830 times. NTK measures how strongly each term pulls on the network, so it sees the gap in size between the terms and makes up for it.
- **RBA** gives the condition points 1.6 to 2.6 times the mean equation weight (the loss uses the square of this weight).
- **BRDR** gives the start shape 13 to 16 times, the fixed ends 22 to 30 times, and the start speed 1.0 to 1.4 times.

RBA scales each point by its residual divided by the largest residual **of the same term**, and BRDR looks at how fast each residual **decays**. Neither one looks at the size of a term next to the other terms, so neither one closes a gap of hundreds of times. This is the same reason as for BRDR in Problem 2, but here it hits RBA too, and the effect on the error is much larger.

## 6. Where the error is

{{fig-error}}

The error at $t = 0$ matches the start-shape residual. For seed 0, it is $1.5 \times 10^{-4}$ for NTK (start-shape residual $1.3 \times 10^{-4}$), $9.6 \times 10^{-3}$ for RBA ($8.3 \times 10^{-3}$), and $1.0 \times 10^{-2}$ for BRDR ($9.7 \times 10^{-3}$).

After that, the error grows with $t$ for all three. For NTK it grows to about $10^{-3}$ by $t = 0.3$, then goes up and down, and ends at $1.9 \times 10^{-3}$ at $t = 1$. For RBA and BRDR it goes up and down and reaches about $4.5 \times 10^{-2}$ at $t = 1$. Nothing damps it, as in Problem 6.

The point weights of RBA and BRDR (seed 0) show no clear pattern. Averaged in five bands in $t$ or in $x$, the equation weights stay between 0.70 and 1.16 times their mean. So I do not show a weight chart.

## 7. Summary

- On the 1D wave equation, NTK reaches $(1.66 \pm 0.29) \times 10^{-3}$. RBA $(5.8 \pm 1.4) \times 10^{-2}$ and BRDR $(4.23 \pm 0.38) \times 10^{-2}$ are more than 18 times worse, and the ranges do not overlap.
- RBA and BRDR fit the equation better than NTK, but they fit the start shape and the fixed ends 20 to 120 times worse.
- In a wave, errors in the start and at the ends travel on and do not fade, so the condition residuals set the final error.
- The terms in the equation are hundreds of times larger than the values in the conditions here. NTK makes up for this with weights of about $10^4$. RBA and BRDR do not look at the size of a term, so they do not.

**Problem 8** (on Monday) is the Fisher–KPP equation.

The code is in [scripts/posts/pinn/07-wave](https://github.com/doubleyands/doubleyands.github.io/tree/main/scripts/posts/pinn/07-wave): `problem.py` defines this problem, and `figure.py` makes the numbers and charts from the saved results. The shared code is [`pinn.py`](https://github.com/doubleyands/doubleyands.github.io/tree/main/scripts/posts/pinn).
