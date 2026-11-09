---
title: "PINNs, Problem 2: 1D Poisson with Two Scales"
description: "Problem 2 of 30 in my PINN series: the 1D Poisson equation -u'' = f with a slow wave and a fast ripple. All 15 runs reach errors between 2e-6 and 4e-4. The worst runs all miss the boundary values, because the equation term cannot see a constant shift."
pubDate: 2026-11-09
project: PINN
tags: [PINN, physics-informed neural networks, loss weighting, NTK, RBA, BRDR, Poisson equation, boundary conditions]
---

This is Problem 2 of my PINN series. The setup is the same as in [Problem 1](/blog/2026-11-06-pinn-01-decay/): Fourier features, an MLP, Adam and then L-BFGS, and three ways to weight the loss (NTK, RBA, BRDR), each with 5 seeds. That post explains all of it. Here I only write what is new.

## 1. The problem

$$
-u''(x) = f(x), \qquad x \in (0, 1), \qquad u(0) = u(1) = 0.
$$

I chose the exact solution first and then computed $f$ from it:

$$
u(x) = \sin(2\pi x) + 0.1\,\sin(10\pi x),
$$

$$
f(x) = 4\pi^2 \sin(2\pi x) + 10\pi^2 \sin(10\pi x).
$$

So the solution has **two scales**: one slow wave and five small, fast ripples on top of it.

{{fig-solution}}

Why is this a test for a PINN? Two reasons.

- **The fast part.** In $u$, the ripple is only 10% of the size of the slow wave. But the second derivative multiplies each part by its frequency squared. So in $f$, the ripple is 2.5 times larger than the slow wave. A network that learns only the slow wave has a large residual.
- **The scales of the two terms are very different.** The equation residual $r_{\text{pde}} = -u_\theta'' - f$ has the size of $f$, which goes up to about 140. The boundary residuals $r_{\text{bc}} = u_\theta(0)$ and $u_\theta(1)$ have the size of $u$. At the first step, the mean squared equation residual is in the thousands, and the boundary one is below 1.

## 2. What is different from Problem 1

- 1,024 random residual points in $(0, 1)$ (a new set for each seed), and the 2 boundary points.
- The test grid has 501 evenly spaced points in $[0, 1]$.
- Everything else (network, Adam 20,000 steps, L-BFGS 2,000 iterations, float64, seeds 0 to 4) is the same.

## 3. Results

| weighting | final error, mean ± std | min – max | after Adam only, median (min – max) |
|---|---|---|---|
| NTK | (2.3 ± 2.2)×10⁻⁵ | 2.1×10⁻⁶ – 5.9×10⁻⁵ | 7.0×10⁻³ (9.5×10⁻⁴ – 8.4×10⁻³) |
| RBA | (5.9 ± 3.8)×10⁻⁵ | 6.5×10⁻⁶ – 9.1×10⁻⁵ | 8.6×10⁻⁴ (9.3×10⁻⁵ – 1.9×10⁻³) |
| BRDR | (1.9 ± 1.7)×10⁻⁴ | 2.25×10⁻⁵ – 4.0×10⁻⁴ | 1.2×10⁻³ (3.5×10⁻⁴ – 4.1×10⁻³) |

All 15 runs finished without problems. NTK has the lowest mean and BRDR the highest, about 8 times higher. But the min–max ranges of all three overlap, so I cannot say that one method is better. The difference is within the spread between seeds. Still, the spread of BRDR has a clear cause. I look at it in Section 4.

{{fig-training}}

The training curves look like the ones in Problem 1:

- **At the end of Adam, NTK is the worst.** Its median error is $7 \times 10^{-3}$, and the point-wise methods are near $10^{-3}$.
- **L-BFGS changes the order.** It cuts the NTK error by a factor of about 300 (median), and NTK ends with the lowest mean.

Did the network learn the fast ripple? Yes. For seed 0, I projected the final error onto $\sin(10\pi x)$. This part of the error is about $2 \times 10^{-6}$ or less for all three methods. The ripple itself has size $0.1$. So the Fourier features with $\sigma = 1$ are enough for this frequency.

## 4. Why some runs miss: the boundary

The next chart shows the error of the final network at each test point (seed 0).

{{fig-error}}

The BRDR curve is almost flat at about $3 \times 10^{-4}$, and it does not go to zero at $x = 0$ and $x = 1$. The network is shifted up by a nearly constant amount, and the boundary values are wrong by that amount. The NTK and RBA curves are much lower and drop toward zero at the boundaries.

The reason is simple. **The equation term cannot see a constant shift.** If we add a constant $c$ (or a line $a + bx$) to $u_\theta$, the second derivative does not change, so $r_{\text{pde}}$ does not change. Only the two boundary points can fix this shift. So the final error depends on how much the loss cares about those two points.

The root mean square of the two boundary residuals at the end of training shows this:

| weighting | boundary weight at the end | boundary error (RMS), 5 seeds |
|---|---|---|
| NTK | λ_bc = 1.1–1.4×10⁴, λ_pde ≈ 1.0 | 7×10⁻⁸ – 4×10⁻⁷ |
| RBA | mean 7.2–7.4, equation points mean 1.9–2.7 | 2.7×10⁻⁶ – 7.0×10⁻⁵ |
| BRDR | mean 0.80–1.19, equation points mean 1.00 | 1.2×10⁻⁵ – 3.4×10⁻⁴ |

- **NTK** gives the boundary term a weight about 10,000 times larger than the equation term. The equation residual contains the second derivative, so it reacts much more strongly to the network parameters. Its NTK is much larger, and NTK balances this. The boundary values end up correct to about $10^{-7}$.
- **RBA** gives the two boundary points a weight about 3 to 4 times larger than an average equation point. RBA squares the weight in the loss, so in the loss this is about 10 times more. This helps, but much less than NTK.
- **BRDR** gives the boundary points about the same weight as an average equation point. BRDR looks at how fast a residual **decays**, not at how large it is. So it does not correct the difference in scale between the two terms. At the end, the mean squared equation residual is $10^{-5}$ to $10^{-4}$, and the boundary one is $10^{-10}$ to $10^{-7}$. The boundary term is a tiny part of the loss, and L-BFGS has no reason to push it further down.

In the three worse BRDR seeds (0, 1, 3), the boundary values are off by $1.4$ to $3.4 \times 10^{-4}$, and the final errors are $1.9$ to $4.0 \times 10^{-4}$. In the two good seeds (2, 4), the boundary values are off by only $1.2$ to $2.3 \times 10^{-5}$, and the final errors are about $2 \times 10^{-5}$. So for BRDR, the boundary error explains most of the final error.

## 5. Summary

- All three weightings solve the two-scale Poisson problem. The final errors are between $2 \times 10^{-6}$ and $4 \times 10^{-4}$. The ranges overlap, so the differences are within the seed spread.
- The fast ripple is not the hard part. The error in the $\sin(10\pi x)$ part is about $2 \times 10^{-6}$ or less.
- The hard part is the boundary. The equation term cannot see a constant shift, so only the 2 boundary points fix it. NTK raises the boundary weight by about $10^4$ and gets the boundary right. BRDR keeps it near 1, and its worst seeds are shifted by about $3 \times 10^{-4}$.
- As in Problem 1, NTK is the worst after Adam, and L-BFGS brings it to the lowest final error.

One simple fix would be a network that satisfies $u(0) = u(1) = 0$ by construction, for example $x(1 - x)\,N_\theta(x)$. I do not use it in this series, because all problems keep the same setup.

**Problem 3** (on Wednesday) is a damped harmonic oscillator.

The code is in [scripts/posts/pinn/02-poisson-1d](https://github.com/doubleyands/doubleyands.github.io/tree/main/scripts/posts/pinn/02-poisson-1d): `problem.py` defines this problem, and `figure.py` makes the numbers and charts from the saved results. The shared code is [`pinn.py`](https://github.com/doubleyands/doubleyands.github.io/tree/main/scripts/posts/pinn).
