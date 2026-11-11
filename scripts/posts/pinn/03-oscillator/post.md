---
title: "PINNs, Problem 3: Damped Harmonic Oscillator"
description: "Problem 3 of 30 in my PINN series: a damped oscillator u'' + 4u' + 400u = 0. All 15 runs end with errors between 5e-6 and 1.4e-4. But with RBA weighting the network stays near the zero solution for thousands of Adam steps, because only one point fixes the size of the solution."
pubDate: 2026-11-11
project: PINN
tags: [PINN, physics-informed neural networks, loss weighting, NTK, RBA, BRDR, ODE, oscillator, initial conditions]
---

This is Problem 3 of my PINN series. The setup is the same as in [Problem 1](/blog/2026-11-06-pinn-01-decay/): Fourier features, an MLP, Adam and then L-BFGS, and three ways to weight the loss (NTK, RBA, BRDR), each with 5 seeds. That post explains all of it. Here I only write what is new.

## 1. The problem

A mass on a spring, with friction:

$$
u''(t) + 2\delta\,u'(t) + \omega_0^2\,u(t) = 0, \qquad t \in (0, 1],
$$

$$
u(0) = 1, \qquad u'(0) = 0, \qquad \delta = 2, \qquad \omega_0 = 20.
$$

The friction is small ($\delta < \omega_0$), so the mass swings back and forth, and the swings get smaller. The exact solution is

$$
u(t) = e^{-\delta t}\left(\cos \omega t + \frac{\delta}{\omega}\sin \omega t\right), \qquad \omega = \sqrt{\omega_0^2 - \delta^2} \approx 19.9.
$$

{{fig-solution}}

Why is this a test for a PINN?

- **It oscillates.** The solution makes about three full swings in $[0, 1]$. In [Problem 1](/blog/2026-11-06-pinn-01-decay/) the solution only decayed.
- **The terms have very different sizes.** The term $\omega_0^2 u$ has a factor 400. At the first step, the mean squared equation residual is about 13,000. The two initial residuals are about 1.
- **Zero is almost a solution.** The equation is linear, and its right side is zero. So $c\,u(t)$ solves the equation for any number $c$, and it also has zero slope at $t = 0$. Even $u = 0$ does this. Only the single point $u(0) = 1$ tells the network the right size $c = 1$.

The last point is the key to this post. If the network is $c\,u(t)$, its relative L2 error is exactly $|1 - c|$, and that is also $|u_\theta(0) - 1|$.

## 2. What is different from Problem 1

- 1,024 random residual points in $(0, 1)$ (a new set for each seed). The two initial conditions are two separate loss terms, one point each: $r_{\text{ic}} = u_\theta(0) - 1$ and $r_{\text{ic},t} = u_\theta'(0)$.
- The test grid has 501 evenly spaced points in $[0, 1]$.
- Everything else (network, Adam 20,000 steps, L-BFGS 2,000 iterations, float64, seeds 0 to 4) is the same.

## 3. Results

| weighting | final error, mean ± std | min – max | after Adam only, median (min – max) |
|---|---|---|---|
| NTK | (6.6 ± 3.4)×10⁻⁵ | 3.6×10⁻⁵ – 1.18×10⁻⁴ | 3.7×10⁻³ (2.4×10⁻³ – 7.3×10⁻³) |
| RBA | (8.3 ± 5.7)×10⁻⁵ | 8.3×10⁻⁶ – 1.35×10⁻⁴ | 1.2×10⁻¹ (9.6×10⁻² – 6.5×10⁻¹) |
| BRDR | (7.3 ± 2.3)×10⁻⁶ | 5.1×10⁻⁶ – 1.06×10⁻⁵ | 3.2×10⁻³ (6.9×10⁻⁵ – 6.1×10⁻²) |

All 15 runs finished without problems. This time one result is clear: **BRDR is better than NTK**. The worst BRDR seed ($1.06 \times 10^{-5}$) is better than the best NTK seed ($3.6 \times 10^{-5}$), so the ranges do not overlap. RBA has the widest range. It overlaps both of the other two, so I cannot rank RBA.

{{fig-training}}

The training curves show a big difference during Adam:

- **NTK and BRDR** get below an error of 0.5 within 500 to 1,900 steps (all 10 runs). At the end of Adam their median error is about $3 \times 10^{-3}$.
- **RBA** stays near an error of 1 for a long time. It gets below 0.5 only after 8,400 to 20,100 steps. One seed gets there only in the L-BFGS phase. At the end of Adam, its median error is still $0.12$.
- **L-BFGS** fixes all of this. It cuts the median error by about 70 times (NTK), 1,000 times (RBA) and 500 times (BRDR). So the final RBA errors are in the same range as NTK.

The next chart shows the error at each test point for seed 0. All three stay below $4 \times 10^{-5}$ at every point, and NTK is the highest.

{{fig-error}}

## 4. Why RBA is slow: the size of the solution

An error near 1 for thousands of steps sounds like a failure. But the network is not wrong in shape. It is wrong in **size**. The next chart shows $|u_\theta(0) - 1|$ during training.

{{fig-ic}}

For RBA, this curve and the error curve are almost the same numbers. At the end of Adam, $|u_\theta(0) - 1|$ and the relative L2 error agree to two digits in all 5 RBA seeds (for example $0.103$ and $0.102$ for seed 0). This is exactly what we expect if the network is $c\,u(t)$ with the wrong $c$. In the first few thousand steps, $c$ is close to 0: the network sits near the zero solution. Then $c$ grows slowly toward 1.

So the speed depends on how much the loss cares about the one point $u(0) = 1$, compared with the 1,024 equation points. Here is the weight of the $u(0)$ term divided by the mean weight of the equation points, as it enters the loss (median of 5 seeds). RBA puts the weight inside the square, $(\lambda r)^2$, so I square its ratio.

| weighting | step 100 | step 1,000 | step 5,000 |
|---|---|---|---|
| NTK | 1.6×10⁵ | 1.5×10⁵ | 1.3×10⁵ |
| RBA | 2.5 | 4.8 | 13 |
| BRDR | 32 | 520 | 140 |

- **NTK** sees from the first step that the equation term reacts much more strongly to the network parameters (it has the factor 400 and a second derivative). It gives $u(0)$ a weight about 150,000 times larger. The size is fixed almost at once.
- **BRDR** raises the $u(0)$ weight quickly. The residual at $u(0)$ decays slowly while the equation residual falls fast, and BRDR gives more weight to residuals that decay slowly. By step 1,000 the ratio is about 500.
- **RBA** can only raise the weight slowly. Each weight starts at 0, grows by at most 0.01 per step, and has an upper limit of 10. For a single point, the update always uses the maximum. So $\lambda_{\text{ic}}$ is about 1 at step 100 and 6.3 at step 1,000, and it is close to the limit 10 by step 5,000. The equation weights grow too, so the ratio stays small. With this small push, Adam moves $c$ from 0 to 1 only very slowly.

Why does the network go to zero first? At the start, the equation term is thousands of times larger than the initial terms. The fastest way to reduce the loss is to make the equation residual small. A small, flat network does that. For RBA, the mean squared equation residual falls from about 13,000 to about 100 in the first 100 steps (seed 0), while $u_\theta(0)$ goes near 0.

L-BFGS repairs the size quickly. This makes sense: the size $c$ is one smooth direction in the loss, and L-BFGS uses curvature information. After L-BFGS, the median $|u_\theta(0) - 1|$ for RBA is $7 \times 10^{-5}$.

One small extra: the RBA weights of the equation points are larger near $t = 0$. In seed 0, the mean weight for $t < 0.1$ is 1.8 times the average, and for $t > 0.6$ it is about 0.7 times.

## 5. Summary

- All three weightings solve the damped oscillator. The final errors are between $5 \times 10^{-6}$ and $1.4 \times 10^{-4}$.
- BRDR has the lowest final error, $(7.3 \pm 2.3) \times 10^{-6}$. Its range does not overlap with NTK, so this difference is larger than the seed spread. RBA overlaps both.
- During Adam, RBA starts near the zero solution, and its error stays above 0.5 for 8,400 to 20,100 steps. The equation alone cannot fix the size of the solution, and RBA gives the single point $u(0)$ only a small weight. NTK and BRDR give it a large weight.
- L-BFGS fixes the size in all RBA runs. Without L-BFGS, RBA would be the clear loser here.

**Problem 4** (on Friday) is a 1D advection–diffusion problem with a boundary layer.

The code is in [scripts/posts/pinn/03-oscillator](https://github.com/doubleyands/doubleyands.github.io/tree/main/scripts/posts/pinn/03-oscillator): `problem.py` defines this problem, and `figure.py` makes the numbers and charts from the saved results. The shared code is [`pinn.py`](https://github.com/doubleyands/doubleyands.github.io/tree/main/scripts/posts/pinn).
