---
title: "PINNs, Problem 1: Exponential Decay"
description: "The first of 30 problems that I solve with physics-informed neural networks. All problems use one setup: Fourier features, Adam and then L-BFGS, and three ways to weight the loss (NTK, RBA, BRDR), each with 5 seeds. Problem 1 is the ODE u' = -2u."
pubDate: 2026-11-06
project: PINN
tags: [PINN, physics-informed neural networks, loss weighting, NTK, RBA, BRDR, ODE]
---

This post starts a series. I will solve **30 equations** with physics-informed neural networks (PINNs), one every Monday, Wednesday and Friday. The problems start easy and get harder: ODEs, then 1D and 2D PDEs, then nonlinear waves, inverse problems and the Navier–Stokes equations.

Every problem uses **the same setup**. Only the equation changes. So the posts compare like with like. In each problem I also compare **three ways to weight the loss**, and I train every case with **5 random seeds**. This first post explains the setup in full. The later posts link back here.

## 1. The problem

The first problem is the simplest one I could think of:

$$
u'(t) = -2\,u(t), \qquad u(0) = 1, \qquad t \in [0, 3].
$$

The exact solution is $u(t) = e^{-2t}$. It falls from $1$ to $e^{-6} \approx 0.0025$. A PINN should solve this easily. So this problem mainly checks the whole pipeline, and it gives a baseline for the harder problems.

## 2. What a PINN does

A PINN is a neural network $u_\theta(t)$ with parameters $\theta$. We never show it the solution. Instead we put the network into the equation and into the initial condition, and we get **residuals**:

$$
r_{\text{pde}}(t_i) = u_\theta'(t_i) + 2\,u_\theta(t_i), \qquad r_{\text{ic}} = u_\theta(0) - 1.
$$

Automatic differentiation gives $u_\theta'$ exactly. The points $t_i$ are 512 random points in $[0, 3]$. Training makes the residuals small. If all residuals are zero, the network solves the problem.

The loss adds up the squared residuals of each **term** (here the equation term and the initial-condition term). The question is how much each term, or each point, should count. This is the job of a **weighting** method.

## 3. The common setup

### The network

- The input is first scaled to $[-1, 1]$. Then it goes through **random Fourier features**: $z \mapsto [\cos(zB), \sin(zB)]$. The matrix $B$ has 128 columns, and its entries come from $N(0, \sigma^2)$ with $\sigma = 1$. With these features, the network learns fast changes more easily.
- Then comes an MLP with 4 hidden layers of 128 units and tanh, with Xavier initialization.
- Everything runs in double precision (float64).

### Three ways to weight the loss

**NTK weights** (Wang, Yu and Perdikaris, *J. Comput. Phys.* 449, 2022). One weight per term. The neural tangent kernel (NTK) tells how fast each term can shrink under gradient descent. Let $m_k$ be the mean of $\lVert \partial r_i / \partial \theta \rVert^2$ over the points of term $k$. Then

$$
\lambda_k = \frac{\sum_j m_j}{m_k}, \qquad L = \sum_k \lambda_k \, \frac{1}{N_k} \sum_i r_{k,i}^2 .
$$

If a term shrinks slowly, it gets a larger weight. I estimate $m_k$ on 128 random points of each term every 1,000 steps, and I smooth the weights with a moving average (factor 0.9).

**RBA weights** (residual-based attention; Anagnostopoulos et al., *Comput. Methods Appl. Mech. Eng.* 421, 2024). One weight per point. Every step, each point adds its residual size to a running sum that slowly forgets:

$$
\lambda_i \leftarrow 0.999\,\lambda_i + 0.01\,\frac{|r_i|}{\max_j |r_j|}, \qquad L = \sum_k \frac{1}{N_k} \sum_i (\lambda_i r_i)^2 .
$$

The weights start at 0, the maximum is taken within each term, and no weight goes above $0.01 / (1 - 0.999) = 10$. Points with large residuals get more attention.

**BRDR weights** (balanced residual decay rate; Chen, Howard and Stinis, *J. Comput. Phys.* 542, 2025). One weight per point. Here the weight depends on how fast a residual **decays**. Each point keeps a moving average of $r_i^4$, and

$$
\text{irdr}_i = \frac{r_i^2}{\sqrt{\overline{r_i^4} + 10^{-14}}}, \qquad
w_i \leftarrow 0.999\,w_i + 0.001\,\frac{\text{irdr}_i}{\operatorname{mean}(\text{irdr})},
$$

$$
L = \sum_k \frac{1}{N_k} \sum_i w_i\, r_i^2 .
$$

If a residual falls slowly, its $\text{irdr}$ stays large, so its point gets more weight. The mean runs over the points of all terms together, so the weights have mean 1. The weights start at 1. The paper also scales the whole loss by a global factor. I leave this factor out, because Adam hardly changes when the whole loss is scaled.

### Optimizer

- **Adam** for 20,000 steps. The learning rate starts at $10^{-3}$ and falls by a factor of 0.9 every 2,000 steps.
- Then **L-BFGS** for 2,000 iterations (strong Wolfe line search, history 50). During L-BFGS the weights stay fixed at their last Adam values, because L-BFGS needs a fixed loss.
- The training points stay fixed during training. Each seed draws its own points and its own network.

### How I measure the error

I use the relative $L^2$ error on a test grid. Here the grid has 301 points in $[0, 3]$:

$$
\text{error} = \frac{\lVert u_\theta - u \rVert_2}{\lVert u \rVert_2}.
$$

Every number below comes from 5 seeds (0 to 4).

## 4. Results

| weighting | final error, mean ± std | min – max | after Adam only, median (min – max) |
|---|---|---|---|
| NTK | (4.9 ± 1.3)×10⁻⁵ | 3.2–6.1 ×10⁻⁵ | 7.2×10⁻³ (6.5×10⁻⁵ – 1.2×10⁻²) |
| RBA | (5.3 ± 1.9)×10⁻⁵ | 3.4–7.6 ×10⁻⁵ | 1.1×10⁻⁴ (5.3×10⁻⁵ – 3.2×10⁻²) |
| BRDR | (2.7 ± 1.1)×10⁻⁵ | 1.6–4.1 ×10⁻⁵ | 4.0×10⁻⁵ (1.6×10⁻⁵ – 8.4×10⁻³) |

All 15 runs finished without problems. The three methods end at almost the same level. BRDR has the lowest mean, but its range overlaps the other two ranges. So for this problem, the difference is within the spread between seeds.

The chart shows how the error changes during training.

{{fig-training}}

Two things stand out.

- **Adam alone is not reliable here.** At the end of Adam, the error depends a lot on the seed. It ranges from $1.6 \times 10^{-5}$ to $3.2 \times 10^{-2}$, a factor of about 2,000. The error also goes up and down during Adam.
- **L-BFGS does the last part of the work.** In the L-BFGS phase, all 15 runs fall to between $1.6 \times 10^{-5}$ and $7.6 \times 10^{-5}$. The weighting matters less than this second phase.

The next chart shows the error of the final network at each test point.

{{fig-error}}

The error stays small over the whole interval. It does not grow toward $t = 3$, even though the solution there is 400 times smaller than at $t = 0$.

## 5. Where the weights go

The two point-wise methods both put **more weight near $t = 0$**. Between $t = 0$ and $0.5$, the average RBA weight is about 2.1 times the overall mean, and the average BRDR weight is about 2.4 times. After $t = 2.5$, both are about half of the mean or less. Near $t = 0$ the solution is large and changes fast, so the residual is largest there.

{{fig-weights}}

The weight of the initial condition shows a clear difference between the methods:

- **NTK** gives the initial condition a weight of 5.8 to 7.2, and the equation term a weight of about 1.2 (all 5 seeds).
- **RBA** gives the single initial-condition point the largest possible weight, 10. With one point, $|r| / \max |r|$ is always 1.
- **BRDR** gives it a small weight, 0.07 to 0.17. Its residual falls fast, so BRDR moves weight to slower points.

So the three methods disagree a lot about the initial condition, and they still reach the same final error. For an easy problem like this one, the weighting is not the bottleneck.

## 6. Summary

- I solved $u' = -2u$ with one fixed PINN setup and three weightings, each with 5 seeds.
- All three reach a relative error of about $3$ to $5 \times 10^{-5}$. The differences are within the seed spread.
- Adam alone leaves errors between $10^{-5}$ and $10^{-2}$, depending on the seed. L-BFGS brings all runs to the same level.
- RBA and BRDR put more weight near $t = 0$. For the initial condition, NTK and RBA raise its weight, but BRDR lowers it.

This is only one easy problem, and one fixed value of $\sigma$. The next problems will test the same setup on harder cases. **Problem 2** (on Monday) is a 1D Poisson equation with two scales.

The code is in [scripts/posts/pinn](https://github.com/doubleyands/doubleyands.github.io/tree/main/scripts/posts/pinn): `pinn.py` has the network, the three weightings and the training loop. `01-decay/problem.py` defines this problem, and `01-decay/figure.py` makes the numbers and charts from the saved results.
