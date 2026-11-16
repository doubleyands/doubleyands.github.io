---
title: "PINNs, Problem 5: The 1D Heat Equation"
description: "Problem 5 of 30 in my PINN series: u_t = 0.1 u_xx on [0, 1] × [0, 1], the first problem with time. NTK and RBA reach errors near 3e-5. BRDR is 5 to 6 times worse, because it fits the initial and boundary conditions loosely."
pubDate: 2026-11-16
project: PINN
tags: [PINN, physics-informed neural networks, loss weighting, NTK, RBA, BRDR, heat equation]
---

This is Problem 5 of my PINN series. The setup is the same as in [Problem 1](/blog/2026-11-06-pinn-01-decay/): Fourier features, an MLP, Adam and then L-BFGS, and three ways to weight the loss (NTK, RBA, BRDR), each with 5 seeds. That post explains all of it. Here I only write what is new.

## 1. The problem

The heat equation in one space dimension:

$$
u_t = a\,u_{xx}, \qquad x \in (0, 1),\; t \in (0, 1], \qquad a = 0.1,
$$

with the walls held at zero and a start made of two sine waves:

$$
u(0, t) = u(1, t) = 0,
$$

$$
u(x, 0) = \sin(\pi x) + 0.5\,\sin(3\pi x).
$$

Each sine wave decays on its own, so the exact solution is

$$
\begin{aligned}
u(x, t) ={}& e^{-a\pi^2 t}\,\sin(\pi x) \\
&+ 0.5\,e^{-9 a\pi^2 t}\,\sin(3\pi x).
\end{aligned}
$$

{{fig-solution}}

This is the first problem in the series with two inputs, $x$ and $t$. The network takes the pair $(x, t)$ and gives one number. Time is just a second input. There are three loss terms: the equation inside the domain, the initial condition at $t = 0$, and the boundary condition at $x = 0$ and $x = 1$.

Why should this be easy for a PINN?

- **The solution is smooth.** There are no layers and no shocks. Both parts are plain sine waves times plain exponentials.
- **The heat equation forgets.** Small errors decay over time, and short waves decay fastest. A wave $\sin(k\pi x)$ decays like $e^{-a k^2 \pi^2 t}$. So the fast part $\sin(3\pi x)$ shrinks 9 times faster than the slow part: by $t = 1$, it is below $10^{-4}$.

What could go wrong? The network does not know the solution at $t = 0$ unless the initial-condition term teaches it. The equation alone has many solutions (any start gives one). So the initial and boundary terms must be fitted well. This will matter in Section 4.

## 2. What is different from Problem 1

- 4,096 random residual points in $(0, 1) \times (0, 1)$, 256 points on the line $t = 0$, and 256 points on the two walls (128 on each). Each seed has its own points.
- The test grid has $101 \times 101$ evenly spaced points in $[0, 1] \times [0, 1]$.
- Everything else (network, Adam 20,000 steps, L-BFGS 2,000 iterations, float64, seeds 0 to 4) is the same.

## 3. Results

| weighting | final error, mean ± std | min – max | after Adam only, median (min – max) |
|---|---|---|---|
| NTK | (3.04 ± 0.87)×10⁻⁵ | 1.95×10⁻⁵ – 4.04×10⁻⁵ | 2.0×10⁻³ (5.9×10⁻⁴ – 7.5×10⁻³) |
| RBA | (3.81 ± 0.83)×10⁻⁵ | 2.95×10⁻⁵ – 5.10×10⁻⁵ | 1.4×10⁻⁴ (6.7×10⁻⁵ – 2.1×10⁻²) |
| BRDR | (1.94 ± 0.48)×10⁻⁴ | 1.37×10⁻⁴ – 2.51×10⁻⁴ | 2.9×10⁻⁴ (1.7×10⁻⁴ – 1.2×10⁻³) |

All 15 runs finished without problems. **NTK and RBA are about equal.** Their ranges overlap, so the difference is within the seed spread. **BRDR is worse.** Its best seed ($1.37 \times 10^{-4}$) is worse than the worst NTK seed ($4.04 \times 10^{-5}$) and the worst RBA seed ($5.10 \times 10^{-5}$). Its mean error is 6.4 times the NTK mean and 5.1 times the RBA mean.

{{fig-training}}

This problem is easy at the start. There is no plateau like in [Problem 4](/blog/2026-11-13-pinn-04-boundary-layer/). All 15 runs get below error 0.01 by step 1,300 (NTK by step 200 in every seed) and below 0.001 by step 2,700.

The two phases tell different stories:

- **Adam.** In the first 5,000 steps, RBA has the highest error most of the time. From step 4,900 to the end of Adam, NTK has the highest error. Its median at the end of Adam is $2.0 \times 10^{-3}$. RBA is lowest at the end of Adam ($1.4 \times 10^{-4}$), except seed 1 ($2.1 \times 10^{-2}$). BRDR is close to RBA ($2.9 \times 10^{-4}$).
- **L-BFGS.** L-BFGS cuts the median error of NTK by about 66 times and RBA by about 4 times. For BRDR, it helps very little: about 1.6 times. So BRDR goes from close to the best at the end of Adam to clearly the worst at the end. I saw the same pattern for NTK in Problems 1 and 2: worst after Adam, then fixed by L-BFGS.

## 4. Why BRDR stops early

Where is the BRDR error? The next chart shows the error of seed 0 at each time $t$, as the RMS over $x$.

{{fig-error}}

For NTK and RBA, the error is between $8 \times 10^{-6}$ and $5 \times 10^{-5}$ at all times. For BRDR, the error is largest at $t = 0$: $2.7 \times 10^{-4}$. So BRDR does not fit the initial condition well. After $t = 0.1$, it is smaller, between $3 \times 10^{-5}$ and $9 \times 10^{-5}$. This fits the "forgetting" in Section 1: the heat equation damps an error in the start. But it does not remove all of it, and the walls also matter. On the two walls $x = 0$ and $x = 1$, the RMS error of seed 0 is $2.1 \times 10^{-4}$ for BRDR, and about $5 \times 10^{-5}$ for NTK and RBA.

The training residuals show the same thing for all 5 seeds. Here is the RMS of each residual at the end of training, and the weight of each condition divided by the weight of the equation term (for RBA and BRDR, the mean weights):

| weighting | equation residual | initial residual | boundary residual | initial weight ÷ equation weight | boundary weight ÷ equation weight |
|---|---|---|---|---|---|
| NTK | 3.3–4.5×10⁻⁴ | 1.7–2.9×10⁻⁵ | 2.2–4.6×10⁻⁵ | 6.7–9.1 | 7.3–10.1 |
| RBA | 3.3–4.5×10⁻⁴ | 1.5–3.6×10⁻⁵ | 3.5–6.2×10⁻⁵ | 3.2–3.8 | 1.9–2.7 |
| BRDR | 3.0–5.8×10⁻⁴ | 1.1–4.8×10⁻⁴ | 1.7–2.9×10⁻⁴ | 0.19–0.38 | 0.17–0.26 |

The equation residual is about the same for all three methods. The difference is in the conditions. For NTK and RBA, the initial and boundary residuals are 7 to 23 times smaller than the equation residual. For BRDR, they are only 1.2 to 2.8 times smaller. So the BRDR initial and boundary residuals are several times larger than for NTK and RBA. And the BRDR weights explain it. BRDR gives the initial and boundary points only 0.17 to 0.38 times the mean weight of the equation points. NTK gives them 7 to 10 times more weight, and RBA 2 to 4 times more (RBA squares its weights in the loss, so they count even more). L-BFGS keeps the Adam weights fixed, so BRDR's L-BFGS hardly pushes on the conditions. This is why it improves so little.

I saw a small BRDR boundary weight in [Problem 2](/blog/2026-11-09-pinn-02-poisson-1d/) and Problem 4 too. In Problem 2, I explained it by size: the equation residual was thousands of times larger than the boundary residual at the start. Here that is not the reason. At step 0, the mean squared residual is 0.64 to 1.42 for the equation, 0.15 to 0.87 for the initial condition, and 0.05 to 0.25 for the boundary (5 seeds). They are not far apart. BRDR sets each weight from how fast that residual is still falling, not from how large it is. So a small residual can keep a small weight, and BRDR has no part that makes the conditions count more.

## 5. Where RBA and BRDR put their weight

RBA and BRDR give each equation point its own weight. The points are random in $(x, t)$, so I group them in bands of width 0.05 in $t$ and show the mean weight of each band (seed 0), divided by the mean of all weights.

{{fig-weights}}

- **RBA** puts more weight at early times. The mean weight is 2.4 for $t < 0.02$, 1.7 for $0.05 \le t < 0.1$, and 0.76 for $t \ge 0.5$. At early times, $u$ changes fastest in time, because the fast wave is still there.
- **BRDR** stays between about 0.8 and 1.4 in every band. It shows no clear pattern in $t$.

## 6. Summary

- All three weightings solve the 1D heat equation. The final errors are between $2 \times 10^{-5}$ and $2.5 \times 10^{-4}$.
- NTK $(3.04 \pm 0.87) \times 10^{-5}$ and RBA $(3.81 \pm 0.83) \times 10^{-5}$ are within the seed spread of each other. BRDR $(1.94 \pm 0.48) \times 10^{-4}$ is worse, and its range does not overlap with them.
- The problem is easy at the start: every run is below error 0.001 by step 2,700, with no plateau.
- BRDR fits the initial and boundary conditions loosely. Its error is largest at $t = 0$, and its weights for the conditions are only 0.17 to 0.38 times the equation weight. With these weights fixed, L-BFGS improves BRDR by only about 1.6 times.
- RBA puts more weight on early times, where the fast wave decays.

**Problem 6** (on Wednesday) is 1D advection with a high speed. It is a known hard case for PINNs.

The code is in [scripts/posts/pinn/05-heat](https://github.com/doubleyands/doubleyands.github.io/tree/main/scripts/posts/pinn/05-heat): `problem.py` defines this problem, and `figure.py` makes the numbers and charts from the saved results. The shared code is [`pinn.py`](https://github.com/doubleyands/doubleyands.github.io/tree/main/scripts/posts/pinn).
