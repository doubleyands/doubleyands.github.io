---
title: "Numerical PDEs, Lecture 5: Upwinding, Measuring the Order of Convergence, and Initial Value Problems"
description: "My notes on the fifth lecture: why the centered difference fails when advection dominates, the upwind fix and its price, how to measure the order of a scheme from numbers, singular perturbations, and the start of ODE initial value problems."
pubDate: 2026-10-20
project: numerical-pde
tags: [numerical PDEs, finite differences, advection-diffusion, convergence, ODEs, notes]
---

This is my summary of the fifth lecture in my numerical PDE course. I write the ideas in my own words and check them with a small experiment.

The lecture has three parts. First it finishes the advection-diffusion story from [Lecture 4](/blog/2026-10-15-numerical-pde-lecture-4/): we watch the centered difference fail, and we fix it with **upwinding**. Then it teaches a very practical skill: how to **measure the order of convergence** of a scheme from your own numbers. Finally, it closes the chapter on space and starts a new one, **ODE initial value problems**, the first step toward time-dependent PDEs.

## 1. Watching the centered difference fail

### A boring problem first

Take pure diffusion, $-\varepsilon u'' = 1$ on $(0, 1)$ with $u(0) = u(1) = 0$. The exact solution is

$$
u(x) = \frac{x(1 - x)}{2\varepsilon}.
$$

Use the centered difference from [Lecture 2](/blog/2026-10-08-numerical-pde-lecture-2/), even with only a few grid points. The numerical solution is **exact at the grid points**, for every $h$ and every $\varepsilon$. The reason is the truncation error. In Lecture 2 we found $\tau = \frac{h^2}{12} u'''' + \cdots$, and every later term has a higher derivative too. Here $u$ is a quadratic, so $u'''' = 0$ and $\tau = 0$ exactly. Then $AE = -\tau = 0$ gives $E = 0$ ([Lecture 3](/blog/2026-10-13-numerical-pde-lecture-3/)). Nothing goes wrong, so this problem teaches us nothing. We need a harder one.

### Now add advection

Take the advection-diffusion problem from Lecture 4:

$$
-\varepsilon u'' + u' = 1, \qquad u(0) = u(1) = 0.
$$

To keep second order, we use the centered difference for $u'$ too:

$$
-\varepsilon\,\frac{U_{i-1} - 2U_i + U_{i+1}}{h^2} + \frac{U_{i+1} - U_{i-1}}{2h} = 1.
$$

Here is my own small example. Take $h = 1/5$, so there are four unknowns $U_1, \dots, U_4$, and let $\varepsilon \to 0$. The diffusion term disappears and each row becomes $U_{i+1} - U_{i-1} = 2h = 2/5$, with $U_0 = U_5 = 0$:

- Row 1 sees only the **left** boundary: $U_2 = U_0 + 2/5 = 2/5$. Then row 3 gives $U_4 = U_2 + 2/5 = 4/5$.
- Row 4 sees only the **right** boundary: $U_3 = U_5 - 2/5 = -2/5$. Then row 2 gives $U_1 = U_3 - 2/5 = -4/5$.

So the numerical solution is

$$
(U_1, U_2, U_3, U_4) = \left(-\tfrac45,\ \tfrac25,\ -\tfrac25,\ \tfrac45\right).
$$

The exact solution is close to $u = x$, which is positive everywhere. But $U_1$ and $U_3$ are negative, and the values jump up and down. Look closer. The even nodes give $2/5$ and $4/5$, which is the line $u = x$. The odd nodes give $-4/5 = 1/5 - 1$ and $-2/5 = 3/5 - 1$, which is the line $u = x - 1$. These are exactly the two candidates from Lecture 4. The real solution keeps the **inflow** condition (left), and the ghost solution keeps the **outflow** condition (right).

The centered difference uses the left and the right neighbor with equal weight. So it treats both boundary conditions as equally real. The left value travels to one family of nodes, the right value travels to the other family, and the two families never agree. More grid points do not remove this pattern by themselves. The two families still alternate, and we see an oscillation. This is the **spurious oscillation** from the end of Lecture 4.

### Upwinding

The flow goes from left to right (the term $u'$ has the coefficient $+1$, so the speed is positive). So the real information comes from the left, the **upwind** side. The fix is to use a one-sided difference that looks upwind:

$$
\frac{U_{i+1} - U_{i-1}}{2h} \quad\longrightarrow\quad \frac{U_i - U_{i-1}}{h}.
$$

Now each $U_i$ is decided by its upwind neighbor $U_{i-1}$. In the same limit $\varepsilon \to 0$ with $h = 1/5$, each row becomes $U_i - U_{i-1} = h$, so $U_i = i/5 = x_i$. That is the real solution $u = x$, with no negative values. In this limit the outflow value $U_5 = 0$ is not used at all. With $\varepsilon > 0$ the diffusion term brings it back, and the drop to $0$ happens in the boundary layer.

This idea is one of the most basic ideas in computational fluid dynamics. But it has a **price**: the one-sided difference is only first order. We lose second-order accuracy, and we get the right physical behavior in return.

### Then is the centered difference just wrong?

No. If $h$ is small enough to resolve the boundary layer, the centered difference converges to the exact solution, and it is second order. The experiment below shows this. The problem is the **cost**:

1. In general we do not know in advance how small $h$ must be.
2. In 1D, a few hundred points cost nothing. In 2D and 3D, the number of unknowns grows like $h^{-2}$ or $h^{-3}$, and a very fine grid gets too expensive.
3. Upwinding gives the right qualitative behavior even on a coarse grid. It is not very accurate there, but it is not nonsense.

So the choice is about cost and benefit. We want a good answer for a small price.

## 2. Measuring the order of convergence

Error analysis with Taylor series gives an **upper bound**: $\|E(h)\| \le C h^m$. When $h$ is small, the bound is close to the real error, so

$$
E(h) \approx C h^m.
$$

The question is: what is $m$? Theory gives one answer. Here we find it **from computed numbers**. For this we need a problem with a known exact solution.

### Method A: a convergence table

Halve $h$ again and again, and compute the error each time. Then

$$
\frac{E(h)}{E(h/2)} \approx \frac{C h^m}{C (h/2)^m} = 2^m,
\qquad\text{so}\qquad
m \approx \frac{\ln E(h) - \ln E(h/2)}{\ln 2}.
$$

One value of $m$ needs **two** errors. So the first row of the table has no order. The ratios are never exactly $2^m$, so the measured orders look like $0.97$ or $2.04$, not $1$ or $2$.

### Method B: a log-log plot

Take the logarithm of $E \approx C h^m$:

$$
\ln E(h) = \ln C + m \ln h.
$$

This is a straight line in $\ln h$. Its **slope** is the order $m$, and its value at $\ln h = 0$ (that is, $h = 1$) is $\ln C$. So the plot also gives an estimate of the constant. In practice $h < 1$ and $E < 1$, so both logarithms are negative, and the points sit in the third quadrant of the $(\ln h, \ln E)$ plane.

### Two traps

**The line is not straight at first.** For large $h$ the error does not follow $Ch^m$ yet. The plot can be flat, or the slope can be $0.8$ when we expect $1$. This range of $h$ is called the **pre-asymptotic** range. The true slope appears only after $h$ is small enough. So do not judge the order from the first two or three points.

**There is no slope at all.** There are two common reasons:

1. The code has a bug. This is the most common reason.
2. The scheme is already exact. Then the error is at round-off level for every $h$, and the line is flat. The pure diffusion problem from section 1 is an example.

### Which norm?

The error vector is $E_i = u(x_i) - U_i$ for $i = 1, \dots, n - 1$. Three norms are common:

| Norm | Definition | Comment |
|---|---|---|
| $\ell^1$ | $\sum_i \lvert E_i \rvert$ | plain sum |
| $\ell^2$ | $\left(\sum_i \lvert E_i \rvert^2\right)^{1/2}$ | the standard one, robust |
| $\ell^\infty$ | $\max_i \lvert E_i \rvert$ | the largest pointwise error, the strictest |

For grid functions we usually use the **discrete $L^2$ norm**, with a factor $h$:

$$
\|E\|_{L^2_h} = \Big(h \sum_i \lvert E_i \rvert^2\Big)^{1/2} \approx \Big(\int_0^1 e(x)^2\,dx\Big)^{1/2}.
$$

The factor $h$ makes the sum a Riemann sum: each term is a rectangle with height $\lvert E_i \rvert^2$ and width $h$. So the discrete norm approximates the integral, a kind of average error.

The next point is my own check, not from the lecture. The factor $h$ matters for the **measured order**. Suppose $\lvert E_i \rvert \sim C h^m$ at all $N \sim 1/h$ points. Then

$$
\Big(\sum_i E_i^2\Big)^{1/2} \sim \big(N C^2 h^{2m}\big)^{1/2} = C\,h^{m - 1/2},
\qquad
\Big(h\sum_i E_i^2\Big)^{1/2} \sim C\,h^{m}.
$$

So without the factor $h$, the plot is still a straight line, but its slope is $m - \tfrac12$. A second-order scheme then looks like order $1.5$. In $d$ space dimensions there are about $h^{-d}$ points, and the factor should be $h^d$. The $\ell^\infty$ norm needs no factor.

### Why this matters

- **Code verification.** If the measured slope does not match the theory, something is wrong in the code. This is the most practical tool for testing a numerical scheme.
- **Evidence for a new method.** For a new scheme there may be no theory yet. A convergence plot is numerical evidence that the scheme works, and that it beats other schemes.
- In engineering, codes are also checked on well-known **benchmark problems**, like the lid-driven cavity flow.

## 3. Singular and regular perturbations

There is one simple test: **where is the small parameter $\varepsilon$?**

**Singular perturbation:** $\varepsilon$ multiplies the **highest** derivative, as in

$$
\varepsilon u'' + u' = f.
$$

The highest derivative decides the type of the equation. When $\varepsilon \to 0$, the order drops, and the number of boundary conditions we can use changes. This is where boundary layers come from (Lecture 4).

**Regular perturbation:** $\varepsilon$ multiplies a **lower** derivative, as in

$$
u'' + \varepsilon u' = f.
$$

The highest derivative is untouched. A small $\varepsilon$ gives a small change in the solution, and nothing dramatic happens.

Think of dye in a river again. If the water is almost still and we add a very weak current, the dye cloud just shifts a little. That is a regular perturbation. If the current is strong and diffusion is very weak, a thin layer appears near the outflow end, and the behavior changes completely. That is a singular perturbation. Singular problems are hard for numerical methods because the solution changes very fast in a very narrow region. Flow past an airplane wing (Navier–Stokes with small viscosity) is the famous example.

## 4. Space is done

This lecture also closes the part of the course on space discretization:

- **1D:** the 3-point difference gives a **tridiagonal** matrix.
- **2D:** the 5-point Laplacian gives a **penta-diagonal** matrix, or a **block tridiagonal** matrix with diagonal blocks $T$ and off-diagonal blocks $I$ ([Lecture 3](/blog/2026-10-13-numerical-pde-lecture-3/)).

## 5. ODE initial value problems

Why ODEs now? The goal is time-dependent PDEs: the heat equation, and later hyperbolic equations and conservation laws. Time and space together are too much at once. So we first remove the space part and study the simplest time problem, an ODE. As before, we first understand the equation itself, and then the numerical scheme.

### The standard form

$$
u'(t) = f\big(u(t), t\big) \quad \text{for } t > t_0, \qquad u(t_0) = \eta.
$$

A simple physical example is **Newton's law of cooling**. A cup of coffee sits in a room. The rate of change of the coffee temperature is proportional to the difference between the coffee temperature and the room temperature. Here $u$ is the coffee temperature, and $\eta$ is its temperature at the start.

For special $f$ (linear, separable, exact, Bernoulli, and so on) we can find the solution with calculus. But for a general $f$ we cannot, and this is why we need numerical methods.

If $f$ does not depend on $t$ directly, so $u' = f(u)$, the equation is called **autonomous**.

### Higher order becomes first order

From now on we only study **first-order** equations. This is enough, because a higher-order equation can always be rewritten as a first-order system. Here is my own example:

$$
u''' + u\,u'' = t, \qquad u(0) = \eta_1,\ u'(0) = \eta_2,\ u''(0) = \eta_3.
$$

The equation is third order, so it needs three initial conditions (think of integrating three times, with three constants). Define new unknowns

$$
u_1 = u, \qquad u_2 = u', \qquad u_3 = u''.
$$

Then

$$
\begin{cases}
u_1' = u_2, \\
u_2' = u_3, \\
u_3' = t - u_1 u_3,
\end{cases}
\qquad
\begin{cases}
u_1(0) = \eta_1, \\
u_2(0) = \eta_2, \\
u_3(0) = \eta_3.
\end{cases}
$$

The first two lines just say what the new variables are. The third line is the original equation, where $u'''$ becomes $u_3'$, $u''$ becomes $u_3$, and $u$ becomes $u_1$. A safe way to avoid mistakes: write $u, u', u''$ in a column and match them one by one.

In vector form, with $\underline{u} = (u_1, u_2, u_3)^T$, this is $\underline{u}' = \underline{f}(\underline{u}, t)$. The underline only means "a vector". Later we drop it.

There is a price again: the high derivatives disappear, but we have more unknowns. Still, turning a problem into an equivalent, easier problem is a very useful idea. And a method for a scalar ODE works for a system in the same way. So from now on we focus on the scalar case.

### Linear ODEs and Duhamel's principle

The ODE is **linear** if $f(u, t) = A(t)u + g(t)$, with $A(t) \in \mathbb{R}^{s \times s}$ and $g(t) \in \mathbb{R}^s$. Take a constant matrix $A$.

If $g \equiv 0$ (**homogeneous**), the solution is

$$
u(t) = e^{A(t - t_0)}\,\eta.
$$

If $g \not\equiv 0$, **Duhamel's principle** gives

$$
u(t) = e^{A(t - t_0)}\,\eta + \int_{t_0}^{t} e^{A(t - s)}\, g(s)\, ds.
$$

We can check this quickly.

- **Initial condition.** At $t = t_0$ the integral is over an empty interval, so it is $0$. Also $e^{0} = I$. So $u(t_0) = \eta$.
- **The equation.** Differentiate. The upper limit of the integral gives $e^{0} g(t) = g(t)$ (the fundamental theorem of calculus). Differentiating $e^{A(t - s)}$ inside the integral gives a factor $A$. So

$$
u'(t) = A e^{A(t - t_0)}\eta + g(t) + A\int_{t_0}^{t} e^{A(t - s)} g(s)\,ds = A u(t) + g(t).
$$

The integral has the form of a **convolution** in $t - s$. One way to read it: the source $g(s)$ at each earlier time $s$ is carried forward for a time $t - s$ by $e^{A(t - s)}$, and we add up all these pieces.

## 6. Experiment

The script solves $-\varepsilon u'' + u' = 1$ with $u(0) = u(1) = 0$ and $\varepsilon = 0.02$. The exact solution (Lecture 4) is

$$
u(x) = x - \frac{e^{(x - 1)/\varepsilon} - e^{-1/\varepsilon}}{1 - e^{-1/\varepsilon}}.
$$

It uses the Python standard library only, with the Thomas algorithm for the tridiagonal systems. The residual of every solve is below $2 \times 10^{-10}$.

**Part 1: the limit and a coarse grid.** With exact fractions, the limit $\varepsilon \to 0$ with $h = 1/5$ gives $(-4/5, 2/5, -2/5, 4/5)$, as in section 1. With $\varepsilon = 0.02$ and the same grid:

| $x$ | 0.2 | 0.4 | 0.6 | 0.8 |
|---|---|---|---|---|
| exact | 0.2000 | 0.4000 | 0.6000 | 0.8000 |
| central | −0.0909 | 0.5455 | 0.0909 | 1.2727 |
| upwind | 0.1999 | 0.3993 | 0.5917 | 0.7091 |

The centered solution has a negative value and jumps up and down. The upwind solution has the right shape. With $h = 1/10$, the centered values are still wrong near the layer: $0.7789$ at $x = 0.7$, then $0.6165$ at $x = 0.8$, then $1.3289$ at $x = 0.9$, where the exact values are $0.7000$, $0.8000$, $0.8933$.

**Part 2: convergence.** The chart shows the largest error for $h = 2^{-2}, \dots, 2^{-12}$. It also shows the pure diffusion problem $-\varepsilon u'' = 1$ with the centered scheme.

<figure class="viz">
<div class="viz-legend"><span class="viz-key"><span class="viz-swatch" style="background:#2a78d6"></span>Central</span><span class="viz-key"><span class="viz-swatch" style="background:#eb6834"></span>Upwind</span><span class="viz-key"><span class="viz-swatch" style="background:#1baf7a"></span>Pure diffusion</span></div>
<div class="viz-frame">
<svg viewBox="0 0 680 400" role="img" tabindex="0" aria-label="Largest error as h gets smaller, on log scales. For the advection-diffusion problem with epsilon = 0.02, the central scheme has a bumpy start and then falls with slope 2. The upwind error first grows, then falls with slope 1. For pure diffusion the error stays at round-off level, near 10 to the minus 14, for every h." class="viz-plot"><line x1="64" x2="520" y1="48.2" y2="48.2" class="viz-grid"/><text x="54" y="52.2" text-anchor="end" class="viz-tick">1</text><line x1="64" x2="520" y1="121.2" y2="121.2" class="viz-grid"/><text x="54" y="125.2" text-anchor="end" class="viz-tick">10<tspan dy="-0.5em" font-size="75%">−4</tspan></text><line x1="64" x2="520" y1="194.1" y2="194.1" class="viz-grid"/><text x="54" y="198.1" text-anchor="end" class="viz-tick">10<tspan dy="-0.5em" font-size="75%">−8</tspan></text><line x1="64" x2="520" y1="267.1" y2="267.1" class="viz-grid"/><text x="54" y="271.1" text-anchor="end" class="viz-tick">10<tspan dy="-0.5em" font-size="75%">−12</tspan></text><line x1="64" x2="520" y1="340.0" y2="340.0" class="viz-grid"/><text x="54" y="344.0" text-anchor="end" class="viz-tick">10<tspan dy="-0.5em" font-size="75%">−16</tspan></text><text x="64.0" y="362" text-anchor="middle" class="viz-tick">2<tspan dy="-0.5em" font-size="75%">−2</tspan></text><text x="155.2" y="362" text-anchor="middle" class="viz-tick">2<tspan dy="-0.5em" font-size="75%">−4</tspan></text><text x="246.4" y="362" text-anchor="middle" class="viz-tick">2<tspan dy="-0.5em" font-size="75%">−6</tspan></text><text x="337.6" y="362" text-anchor="middle" class="viz-tick">2<tspan dy="-0.5em" font-size="75%">−8</tspan></text><text x="428.8" y="362" text-anchor="middle" class="viz-tick">2<tspan dy="-0.5em" font-size="75%">−10</tspan></text><text x="520.0" y="362" text-anchor="middle" class="viz-tick">2<tspan dy="-0.5em" font-size="75%">−12</tspan></text><text x="292.0" y="388" text-anchor="middle" class="viz-axis-title">Step size h (smaller to the right)</text><text x="64" y="14" class="viz-axis-title">Error</text><path d="M64.0,45.7 L109.6,53.3 L155.2,58.8 L200.8,67.6 L246.4,79.4 L292.0,90.7 L337.6,101.7 L383.2,112.7 L428.8,123.7 L474.4,134.6 L520.0,145.6" fill="none" stroke="#2a78d6" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><circle cx="64.0" cy="45.7" r="3" fill="#2a78d6" stroke="#ffffff" stroke-width="1.5"/><circle cx="109.6" cy="53.3" r="3" fill="#2a78d6" stroke="#ffffff" stroke-width="1.5"/><circle cx="155.2" cy="58.8" r="3" fill="#2a78d6" stroke="#ffffff" stroke-width="1.5"/><circle cx="200.8" cy="67.6" r="3" fill="#2a78d6" stroke="#ffffff" stroke-width="1.5"/><circle cx="246.4" cy="79.4" r="3" fill="#2a78d6" stroke="#ffffff" stroke-width="1.5"/><circle cx="292.0" cy="90.7" r="3" fill="#2a78d6" stroke="#ffffff" stroke-width="1.5"/><circle cx="337.6" cy="101.7" r="3" fill="#2a78d6" stroke="#ffffff" stroke-width="1.5"/><circle cx="383.2" cy="112.7" r="3" fill="#2a78d6" stroke="#ffffff" stroke-width="1.5"/><circle cx="428.8" cy="123.7" r="3" fill="#2a78d6" stroke="#ffffff" stroke-width="1.5"/><circle cx="474.4" cy="134.6" r="3" fill="#2a78d6" stroke="#ffffff" stroke-width="1.5"/><circle cx="520.0" cy="145.6" r="3" fill="#2a78d6" stroke="#ffffff" stroke-width="1.5"/><path d="M64.0,68.9 L109.6,64.0 L155.2,61.0 L200.8,61.8 L246.4,66.0 L292.0,70.2 L337.6,75.2 L383.2,80.4 L428.8,85.7 L474.4,91.1 L520.0,96.6" fill="none" stroke="#eb6834" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><circle cx="64.0" cy="68.9" r="3" fill="#eb6834" stroke="#ffffff" stroke-width="1.5"/><circle cx="109.6" cy="64.0" r="3" fill="#eb6834" stroke="#ffffff" stroke-width="1.5"/><circle cx="155.2" cy="61.0" r="3" fill="#eb6834" stroke="#ffffff" stroke-width="1.5"/><circle cx="200.8" cy="61.8" r="3" fill="#eb6834" stroke="#ffffff" stroke-width="1.5"/><circle cx="246.4" cy="66.0" r="3" fill="#eb6834" stroke="#ffffff" stroke-width="1.5"/><circle cx="292.0" cy="70.2" r="3" fill="#eb6834" stroke="#ffffff" stroke-width="1.5"/><circle cx="337.6" cy="75.2" r="3" fill="#eb6834" stroke="#ffffff" stroke-width="1.5"/><circle cx="383.2" cy="80.4" r="3" fill="#eb6834" stroke="#ffffff" stroke-width="1.5"/><circle cx="428.8" cy="85.7" r="3" fill="#eb6834" stroke="#ffffff" stroke-width="1.5"/><circle cx="474.4" cy="91.1" r="3" fill="#eb6834" stroke="#ffffff" stroke-width="1.5"/><circle cx="520.0" cy="96.6" r="3" fill="#eb6834" stroke="#ffffff" stroke-width="1.5"/><path d="M64.0,340.0 L109.6,322.7 L155.2,311.7 L200.8,300.3 L246.4,296.9 L292.0,299.0 L337.6,274.4 L383.2,268.4 L428.8,261.7 L474.4,253.0 L520.0,248.2" fill="none" stroke="#1baf7a" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><circle cx="64.0" cy="340.0" r="3" fill="#1baf7a" stroke="#ffffff" stroke-width="1.5"/><circle cx="109.6" cy="322.7" r="3" fill="#1baf7a" stroke="#ffffff" stroke-width="1.5"/><circle cx="155.2" cy="311.7" r="3" fill="#1baf7a" stroke="#ffffff" stroke-width="1.5"/><circle cx="200.8" cy="300.3" r="3" fill="#1baf7a" stroke="#ffffff" stroke-width="1.5"/><circle cx="246.4" cy="296.9" r="3" fill="#1baf7a" stroke="#ffffff" stroke-width="1.5"/><circle cx="292.0" cy="299.0" r="3" fill="#1baf7a" stroke="#ffffff" stroke-width="1.5"/><circle cx="337.6" cy="274.4" r="3" fill="#1baf7a" stroke="#ffffff" stroke-width="1.5"/><circle cx="383.2" cy="268.4" r="3" fill="#1baf7a" stroke="#ffffff" stroke-width="1.5"/><circle cx="428.8" cy="261.7" r="3" fill="#1baf7a" stroke="#ffffff" stroke-width="1.5"/><circle cx="474.4" cy="253.0" r="3" fill="#1baf7a" stroke="#ffffff" stroke-width="1.5"/><circle cx="520.0" cy="248.2" r="3" fill="#1baf7a" stroke="#ffffff" stroke-width="1.5"/><text x="530" y="149.6" class="viz-label">Central</text><text x="530" y="100.6" class="viz-label">Upwind</text><text x="530" y="252.2" class="viz-label">Pure diffusion</text><line x1="0" x2="0" y1="30" y2="340" class="viz-crosshair" visibility="hidden"/></svg>
<div class="viz-tooltip" hidden></div>
</div>
<figcaption>Largest error max |u(x<sub>j</sub>) − U<sub>j</sub>| on log scales. Central and upwind: −εu″ + u′ = 1 with ε = 0.02. Pure diffusion: −εu″ = 1 with the centered scheme. An error of exactly 0 (pure diffusion, h = 1/4) is drawn at 10<sup>−16</sup>. Hover over the chart, or focus it and use the arrow keys, to read the values.</figcaption>
<details class="viz-table"><summary>Show the numbers</summary><table><thead><tr><th>h</th><th>Central</th><th>Upwind</th><th>Pure diffusion</th></tr></thead><tbody><tr><td>1/4</td><td>1.4e0</td><td>7.4e-2</td><td>0</td></tr><tr><td>1/8</td><td>5.2e-1</td><td>1.4e-1</td><td>8.9e-16</td></tr><tr><td>1/16</td><td>2.6e-1</td><td>2.0e-1</td><td>3.6e-15</td></tr><tr><td>1/32</td><td>8.7e-2</td><td>1.8e-1</td><td>1.5e-14</td></tr><tr><td>1/64</td><td>2.0e-2</td><td>1.1e-1</td><td>2.3e-14</td></tr><tr><td>1/128</td><td>4.7e-3</td><td>6.2e-2</td><td>1.8e-14</td></tr><tr><td>1/256</td><td>1.2e-3</td><td>3.3e-2</td><td>3.9e-13</td></tr><tr><td>1/512</td><td>2.9e-4</td><td>1.7e-2</td><td>8.4e-13</td></tr><tr><td>1/1024</td><td>7.3e-5</td><td>8.8e-3</td><td>2.0e-12</td></tr><tr><td>1/2048</td><td>1.8e-5</td><td>4.4e-3</td><td>5.9e-12</td></tr><tr><td>1/4096</td><td>4.6e-6</td><td>2.2e-3</td><td>1.1e-11</td></tr></tbody></table></details>
<script>(() => {const fig = document.currentScript.closest('figure');const svg = fig.querySelector('svg');const line = svg.querySelector('.viz-crosshair');const tip = fig.querySelector('.viz-tooltip');const names = ["Central", "Upwind", "Pure diffusion"]; const colors = ["#2a78d6", "#eb6834", "#1baf7a"];const rows = [["1.4e0", "7.4e-2", "0"], ["5.2e-1", "1.4e-1", "8.9e-16"], ["2.6e-1", "2.0e-1", "3.6e-15"], ["8.7e-2", "1.8e-1", "1.5e-14"], ["2.0e-2", "1.1e-1", "2.3e-14"], ["4.7e-3", "6.2e-2", "1.8e-14"], ["1.2e-3", "3.3e-2", "3.9e-13"], ["2.9e-4", "1.7e-2", "8.4e-13"], ["7.3e-5", "8.8e-3", "2.0e-12"], ["1.8e-5", "4.4e-3", "5.9e-12"], ["4.6e-6", "2.2e-3", "1.1e-11"]]; const heads = ["h = 1/4", "h = 1/8", "h = 1/16", "h = 1/32", "h = 1/64", "h = 1/128", "h = 1/256", "h = 1/512", "h = 1/1024", "h = 1/2048", "h = 1/4096"];const left = 64, right = 520, steps = 10, width = 680;let current = -1;const show = (n) => { n = Math.max(0, Math.min(steps, n)); current = n; const x = left + (right - left) * n / steps; line.setAttribute('x1', x); line.setAttribute('x2', x); line.setAttribute('visibility', 'visible'); tip.replaceChildren(); const head = document.createElement('div'); head.className = 'viz-tooltip-head'; head.textContent = heads[n]; tip.append(head); names.forEach((name, i) => {  const row = document.createElement('div'); row.className = 'viz-tooltip-row';  const key = document.createElement('span'); key.className = 'viz-swatch'; key.style.background = colors[i];  const value = document.createElement('strong'); value.textContent = rows[n][i];  const label = document.createElement('span'); label.textContent = name;  row.append(key, value, label); tip.append(row); }); const box = svg.getBoundingClientRect(), frame = fig.getBoundingClientRect(); const xs = box.left - frame.left + x * box.width / width; tip.hidden = false; const room = frame.width - tip.offsetWidth - 8; tip.style.left = Math.max(0, Math.min(room, xs + 12)) + 'px'; tip.style.top = (box.top - frame.top + 8) + 'px';};const hide = () => { line.setAttribute('visibility', 'hidden'); tip.hidden = true; current = -1; };svg.addEventListener('pointermove', (event) => { const box = svg.getBoundingClientRect(); const x = (event.clientX - box.left) * width / box.width; show(Math.round((x - left) / (right - left) * steps));});svg.addEventListener('pointerleave', hide);svg.addEventListener('focus', () => show(current < 0 ? 0 : current));svg.addEventListener('blur', hide);svg.addEventListener('keydown', (event) => { if (event.key === 'ArrowRight') { show(current + 1); event.preventDefault(); } if (event.key === 'ArrowLeft') { show(current - 1); event.preventDefault(); } if (event.key === 'Escape') hide();});})();</script>
</figure>

| $h$ | central, max | order | upwind, max | order | central, $L^2_h$ | order | central, $\ell^2$ (no $h$) | order |
|---|---|---|---|---|---|---|---|---|
| 1/4 | $1.38$ | — | $7.40 \times 10^{-2}$ | — | $8.42 \times 10^{-1}$ | — | $1.68$ | — |
| 1/8 | $5.25 \times 10^{-1}$ | 1.39 | $1.36 \times 10^{-1}$ | −0.88 | $2.15 \times 10^{-1}$ | 1.97 | $6.09 \times 10^{-1}$ | 1.47 |
| 1/16 | $2.63 \times 10^{-1}$ | 0.99 | $1.98 \times 10^{-1}$ | −0.55 | $6.69 \times 10^{-2}$ | 1.68 | $2.68 \times 10^{-1}$ | 1.18 |
| 1/32 | $8.68 \times 10^{-2}$ | 1.60 | $1.81 \times 10^{-1}$ | 0.14 | $1.62 \times 10^{-2}$ | 2.04 | $9.18 \times 10^{-2}$ | 1.54 |
| 1/64 | $1.96 \times 10^{-2}$ | 2.14 | $1.06 \times 10^{-1}$ | 0.77 | $3.76 \times 10^{-3}$ | 2.11 | $3.01 \times 10^{-2}$ | 1.61 |
| 1/128 | $4.69 \times 10^{-3}$ | 2.07 | $6.21 \times 10^{-2}$ | 0.77 | $9.11 \times 10^{-4}$ | 2.05 | $1.03 \times 10^{-2}$ | 1.55 |
| 1/256 | $1.17 \times 10^{-3}$ | 2.00 | $3.32 \times 10^{-2}$ | 0.90 | $2.26 \times 10^{-4}$ | 2.01 | $3.61 \times 10^{-3}$ | 1.51 |
| 1/512 | $2.93 \times 10^{-4}$ | 2.00 | $1.73 \times 10^{-2}$ | 0.95 | $5.62 \times 10^{-5}$ | 2.00 | $1.27 \times 10^{-3}$ | 1.50 |
| 1/1024 | $7.31 \times 10^{-5}$ | 2.00 | $8.80 \times 10^{-3}$ | 0.97 | $1.41 \times 10^{-5}$ | 2.00 | $4.50 \times 10^{-4}$ | 1.50 |
| 1/2048 | $1.83 \times 10^{-5}$ | 2.00 | $4.45 \times 10^{-3}$ | 0.99 | $3.51 \times 10^{-6}$ | 2.00 | $1.59 \times 10^{-4}$ | 1.50 |
| 1/4096 | $4.57 \times 10^{-6}$ | 2.00 | $2.23 \times 10^{-3}$ | 0.99 | $8.78 \times 10^{-7}$ | 2.00 | $5.62 \times 10^{-5}$ | 1.50 |

What I noticed:

- **Both schemes reach their theoretical order, but only for small $h$.** The centered scheme settles at $2.00$ from $h = 1/256$. The upwind scheme climbs slowly to $0.99$. Before that, both are in the pre-asymptotic range.
- **The upwind error first grows.** From $h = 1/4$ to $1/16$ the upwind error goes up, so the "order" is negative. My reading: on a very coarse grid there is no node inside the layer, so the scheme does not even see the layer. Once some nodes enter the layer, the error first shows up and then starts to shrink. A plot with only these first points would say the scheme diverges.
- **The centered scheme wins once the grid is fine.** At $h = 1/4096$ its error is about $490$ times smaller than the upwind error. On a coarse grid it is the other way around. This is the cost question from section 1.
- **The factor $h$ in the norm changes the slope by $1/2$.** For the centered scheme, the last column settles at $1.50$, while the $L^2_h$ column settles at $2.00$. For the upwind scheme (not in the table) the script gives $0.49$ without the factor and $0.99$ with it. This matches the estimate in section 2.
- **Pure diffusion gives a flat line.** Its error is $0$ at $h = 1/4$ and stays between $10^{-16}$ and $10^{-11}$ for every $h$. It even grows slowly as $h$ gets smaller, because there are more round-off errors in a bigger system. There is no slope to measure, because the scheme is exact at the nodes. This is the second trap from section 2.

The full script and the chart code are in [this site's repository](https://github.com/doubleyands/doubleyands.github.io/tree/main/scripts/posts/numerical-pde/lecture-5).

## Summary

1. When advection dominates, the **centered difference** listens to both boundaries equally. The inflow value and the outflow value travel to alternate nodes, and the solution oscillates.
2. **Upwinding** uses the neighbor on the side the flow comes from. It gives the right physical behavior even on a coarse grid. The price is first-order accuracy.
3. The centered scheme does converge if $h$ resolves the boundary layer, but that can be too expensive, especially in 2D and 3D.
4. To **measure the order**, halve $h$ and use $m \approx \log_2\big(E(h)/E(h/2)\big)$, or read the slope of a log-log plot. Ignore the pre-asymptotic points. A flat line means a bug or an exact scheme.
5. Use the **discrete $L^2$ norm** with the factor $h$ (or $h^d$ in $d$ dimensions). Without it, the measured order drops by $1/2$.
6. A perturbation is **singular** if $\varepsilon$ multiplies the highest derivative, and **regular** otherwise. Singular problems have thin layers.
7. An **ODE initial value problem** is $u' = f(u, t)$, $u(t_0) = \eta$. A higher-order ODE becomes a first-order system. A linear system with constant $A$ is solved by **Duhamel's principle**.
