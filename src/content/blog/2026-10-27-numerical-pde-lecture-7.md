---
title: "Numerical PDEs, Lecture 7: Runge–Kutta Methods, Multistep Methods, and a First Definition of Stability"
description: "My notes on the seventh lecture: the midpoint idea and RK2, Simpson's rule and RK4, why Runge–Kutta methods need no derivatives, linear multistep methods (leapfrog, Adams–Bashforth, Adams–Moulton) and their starting values, the Gibbs phenomenon, and the definition of stability for a linear system."
pubDate: 2026-10-27
project: numerical-pde
tags: [numerical PDEs, ODEs, time stepping, Runge-Kutta, multistep methods, stability, notes]
---

This is my summary of the seventh lecture in my numerical PDE course. I write the ideas in my own words and check them with a small experiment.

[Lecture 6](/blog/2026-10-22-numerical-pde-lecture-6/) ended with one idea: a time-stepping scheme is a rule for the area under $f$ in one step. This lecture goes in three directions from there. We put the point in the **middle** (RK2). We use **three** points (RK4). Or we **reuse old values** instead of making new ones (multistep methods). At the end, the word **stability** gets its first precise definition.

## 1. Recap: a scheme is an area rule

Integrate $u' = f(u, t)$ over one step:

$$
u(t_{n+1}) = u(t_n) + \int_{t_n}^{t_{n+1}} f\big(u(t), t\big)\,dt.
$$

A quadrature rule picks **where** to measure $f$ and **how much weight** to give each value. Lecture 6 gave three choices:

| Points used | Scheme |
|---|---|
| left end | Forward Euler |
| right end | Backward Euler |
| both ends | trapezoidal |

This lecture adds two more: one point in the **middle**, and **three** points (both ends and the middle).

## 2. One point in the middle

The **midpoint rule** measures $f$ only at the center of the step:

$$
U^{n+1} = U^n + \Delta t\, f\big(U^{n+\frac12},\, t_{n+\frac12}\big), \qquad t_{n+\frac12} = \frac{t_n + t_{n+1}}{2}.
$$

It is the simplest Gauss quadrature rule. Why is it good? The left end loses some area and the right end adds too much. The middle point does both: on one half it is too high, on the other half it is too low, and these errors cancel to leading order. This is the same story as the centered difference in [Lecture 2](/blog/2026-10-08-numerical-pde-lecture-2/). There, a balanced choice of two points gave second order. Here, **one** balanced point gives second order.

**But there is a problem.** A finite difference method only computes values at grid times, $U^n$ and $U^{n+1}$. The value $U^{n+1/2}$ does not exist on the grid. The exact solution lives between $t_n$ and $t_{n+1}$, but the moment we discretize, that information is gone. In [Lecture 3](/blog/2026-10-13-numerical-pde-lecture-3/) we met a ghost node outside the domain. This one is a ghost **inside**, between two grid points. The formula above is a good idea, but it is not yet a method, because we cannot compute it.

## 3. RK2: fill the gap with half an Euler step

The fix uses the cheapest tool we have. Take a Forward Euler step of **half** the size to guess the middle value, and then use the slope there for the full step:

$$
U^* = U^n + \frac{\Delta t}{2}\, f(U^n, t_n),
\qquad
U^{n+1} = U^n + \Delta t\, f\big(U^*,\, t_{n+\frac12}\big).
$$

Why $\Delta t/2$? The midpoint is half a step away from $t_n$. We go half a step to look, then come back and take the whole step with what we saw.

Two words help here:

- **One-step** (single-step): the method goes from $t_n$ to $t_{n+1}$ using only $U^n$.
- **Two-stage**: inside that one step, it evaluates $f$ **twice**.

$U^*$ is only a helper. It is a decent guess for $u(t_{n+1/2})$, not a grid value, and we throw it away after the step. Eliminating it gives one line:

$$
U^{n+1} = U^n + \Delta t\, f\Big(U^n + \tfrac{\Delta t}{2} f(U^n, t_n),\; t_{n+\frac12}\Big).
$$

This method is **explicit** and **second order**. I call it RK2 here. A note on names: some books use "midpoint method" for a different scheme, the leapfrog method of section 7. So check which one a book means.

A small example of my own: $u' = -2u$, $u(0) = 1$, one step with $\Delta t = 0.5$.

- Forward Euler: $1 + 0.5 \cdot (-2) = 0$.
- RK2: $U^* = 1 + 0.25 \cdot (-2) = 0.5$, then $U^1 = 1 + 0.5 \cdot (-2 \cdot 0.5) = 0.5$.
- Exact: $e^{-1} \approx 0.368$.

Forward Euler misses by $0.37$, RK2 by $0.13$. The step is large on purpose, so the errors are easy to see.

## 4. Simpson's rule and RK4

**Simpson's rule** uses three points, both ends and the middle. It fits a parabola through them and takes its area. The weights are $1 : 4 : 1$:

$$
\int_{t_n}^{t_{n+1}} f\,dt \approx \frac{\Delta t}{6}\Big(f(t_n) + 4 f(t_{n+\frac12}) + f(t_{n+1})\Big).
$$

The classical **RK4** method (four stages, fourth order) follows this pattern:

$$
\begin{aligned}
Y_1 &= U^n, \\
Y_2 &= U^n + \tfrac12\Delta t\, f(Y_1,\, t_n), \\
Y_3 &= U^n + \tfrac12\Delta t\, f\big(Y_2,\, t_{n+\frac12}\big), \\
Y_4 &= U^n + \Delta t\, f\big(Y_3,\, t_{n+\frac12}\big), \\
U^{n+1} &= U^n + \frac{\Delta t}{6}\Big[f(Y_1, t_n) + 2 f\big(Y_2, t_{n+\frac12}\big) + 2 f\big(Y_3, t_{n+\frac12}\big) + f(Y_4, t_{n+1})\Big].
\end{aligned}
$$

How I read it:

- $Y_2$ and $Y_3$ are **both** guesses at the middle time $t_{n+1/2}$. $Y_2$ uses the slope at $Y_1$, and $Y_3$ uses the slope at $Y_2$.
- $Y_4$ goes a **full** step with the slope at $Y_3$.
- The weights are $1, 2, 2, 1$. The two middle values together have weight $4$, so this is $1 : 4 : 1$, Simpson's rule.
- If $f$ depends on $t$ only (not on $u$), the $Y$'s do not matter, and RK4 **is** Simpson's rule. The experiment checks this.
- It is **explicit**: compute the four lines in order and we are done. No nonlinear solve.

## 5. Why Runge–Kutta? Evaluations instead of derivatives

The Taylor methods of Lecture 6 got stuck. To reach high order they need $u''$, $u'''$, and so on, and each one means differentiating $f$ again with the chain rule. That is a lot of work, and $f$ may not even be smooth enough.

**Runge–Kutta methods take no derivatives at all.** They only **evaluate** $f$, several times. Since $u' = f(u)$, one evaluation of $f$ already gives one derivative of $u$. In RK2 we have $f$ **inside** $f$: $f\big(U^n + \tfrac{\Delta t}{2} f(U^n)\big)$. Nesting two evaluations has the effect of one more derivative.

We can see this exactly on the linear test problem $u' = \lambda u$. Put $z = \lambda \Delta t$. One step gives

$$
\text{RK2:}\quad U^{n+1} = \Big(1 + z + \frac{z^2}{2}\Big) U^n,
\qquad
\text{RK4:}\quad U^{n+1} = \Big(1 + z + \frac{z^2}{2} + \frac{z^3}{6} + \frac{z^4}{24}\Big) U^n.
$$

The exact solution gives $e^{z} U^n$. So RK2 and RK4 match the Taylor series of $e^z$ up to $z^2$ and $z^4$, and we never computed a derivative of $f$.

RK2 and RK4 are two members of a big **family**. There are also implicit Runge–Kutta methods. The whole family can be written in a compact table of coefficients, the **Butcher tableau**. The lecture only gave the main idea, and I stop here too.

## 6. Linear multistep methods

A Runge–Kutta method builds helper values **inside** a step. A multistep method does the opposite. It uses **old grid values** that we already computed.

The general **$r$-step linear multistep method** is

$$
\sum_{j=0}^{r} \alpha_j\, U^{n+j} = \Delta t \sum_{j=0}^{r} \beta_j\, f\big(U^{n+j},\, t_{n+j}\big).
$$

- The left side is a combination of grid values. The newest one is $U^{n+r}$.
- The right side is a weighted sum of $f$ values times $\Delta t$, so it is an **area**.
- $r = 1$ is a one-step method, $r = 2$ uses two old values, and so on.

The test for explicit or implicit is simple. If the newest value $U^{n+r}$ appears inside $f$, we must solve an equation. So the method is **explicit if $\beta_r = 0$**, and implicit otherwise.

### Leapfrog

Use a centered difference in time around $t_n$:

$$
\frac{U^{n+1} - U^{n-1}}{2\Delta t} = f(U^n, t_n)
\qquad\Longrightarrow\qquad
U^{n+1} = U^{n-1} + 2\Delta t\, f(U^n, t_n).
$$

The reference time is $t_n$, and the scheme uses the values before and after it. It is exactly the centered difference, so it is second order. $U^{n+1}$ is not inside $f$, so it is explicit. The new value jumps **over** $U^n$, from $U^{n-1}$ to $U^{n+1}$. That is why it is called **leapfrog**.

A fun special case: take $f(u) = u$ and $\Delta t = 1/2$. Then $U^{n+1} = U^{n-1} + U^n$. This is the **Fibonacci** recursion.

## 7. Multistep methods cannot start by themselves

To get $U^2$ from leapfrog we need $U^0$ **and** $U^1$. The initial condition gives $U^0 = \eta$. Nobody gives us $U^1$.

It is like a second-order ODE, which needs two initial conditions. A one-step method needs only one value, so it is **self-starting**. An $r$-step method needs $r$ values, and we must make all but $U^0$ ourselves, with some other method.

Which method should make the starting values? Explicit and second-order RK2 is a natural choice for leapfrog. A natural worry is that a cheap first-order start, like Forward Euler, spoils the order of the whole method. Here is the rule I found when I checked this (my own reading, confirmed by the experiment below):

- The start is used only **once**. Its error is not added up $N$ times. So what matters is the **one-step** error of the starting method, not its global order.
- One Forward Euler step has error $O(\Delta t^2)$ (Lecture 6, section 11). A method of order $p$ needs starting values with error $O(\Delta t^p)$.
- So for leapfrog ($p = 2$), even a Forward Euler start keeps second order. For a fourth-order method, a Forward Euler start really does destroy the order: it drops to two.

So the worry is real, but it depends on the order of the method. For high-order methods, use a high-order starter, for example RK4.

## 8. The Adams family

**Adams methods** fix the left side: $\alpha_r = 1$, $\alpha_{r-1} = -1$, and all other $\alpha_j = 0$. So

$$
U^{n+r} = U^{n+r-1} + \Delta t \sum_{j=0}^{r} \beta_j\, f^{n+j}, \qquad f^{n+j} = f\big(U^{n+j}, t_{n+j}\big).
$$

In words: **the last value plus an area**. The methods differ only in how they weight the $f$ values. **Adams–Bashforth** (AB) methods are explicit ($\beta_r = 0$): they extrapolate from old $f$ values. **Adams–Moulton** (AM) methods are implicit: they also use the newest $f$.

| $r$ | Adams–Bashforth (explicit) | Adams–Moulton (implicit) |
|---|---|---|
| 1 | $\Delta t\, f^n$ = Forward Euler | $\frac{\Delta t}{2}\big(f^n + f^{n+1}\big)$ = trapezoidal |
| 2 | $\frac{\Delta t}{2}\big(-f^n + 3f^{n+1}\big)$ | $\frac{\Delta t}{12}\big(-f^n + 8f^{n+1} + 5f^{n+2}\big)$ |
| 3 | $\frac{\Delta t}{12}\big(5f^n - 16f^{n+1} + 23f^{n+2}\big)$ | $\frac{\Delta t}{24}\big(f^n - 5f^{n+1} + 19f^{n+2} + 9f^{n+3}\big)$ |
| 4 | $\frac{\Delta t}{24}\big(-9f^n + 37f^{n+1} - 59f^{n+2} + 55f^{n+3}\big)$ | $\frac{\Delta t}{720}\big(-19f^n + 106f^{n+1} - 264f^{n+2} + 646f^{n+3} + 251f^{n+4}\big)$ |

The $r$-step AB method has order $r$. The $r$-step AM method has order $r + 1$. (The trapezoidal method is AM with $r = 1$, and it is second order.)

**A sanity check.** If $f \equiv 1$, the area over one step must be exactly $\Delta t$. So in each row, the coefficients inside the brackets must add up to the denominator. For example, AB4 gives $-9 + 37 - 59 + 55 = 24$, and AM3 gives $1 - 5 + 19 + 9 = 24$. If the sum is anything else, there is a typo. It plays the same role as the "coefficients sum to zero" check for difference formulas in [Lecture 2](/blog/2026-10-08-numerical-pde-lecture-2/). It is a cheap way to catch a typo when you copy one of these tables.

Adams–Moulton methods are implicit, so for nonlinear $f$ every step needs Newton's method, as in Lecture 6.

## 9. One-step or multistep?

| | One-step (Runge–Kutta, ...) | Multistep (leapfrog, Adams, ...) |
|---|---|---|
| Starting | self-starting | needs another method for the first values |
| Changing $\Delta t$ | any time, for example from an error estimate | the formulas assume equal steps |
| A known kink at time $t^*$ | put $t^*$ on the grid, and accuracy is kept | the formula mixes data from both sides of $t^*$ and loses accuracy |
| New $f$ evaluations per step | several (four for RK4) | **one**; the old ones are reused |

The last row is the reason multistep methods exist. RK4 calls $f$ four times per step. AB4 reuses the three $f$ values from earlier steps and computes only one new value. If $f$ is expensive, this matters a lot. The price is in the first three rows: no self-start, no easy step change, and trouble near a non-smooth point.

About the third row: if we know in advance that the solution is not smooth at some time $t^*$, make $t^*$ a grid point. Then each step sees only one smooth side. A one-step method can do this because it can change its step size freely.

## 10. The Gibbs phenomenon

What happens if we approximate a jump with smooth functions and we do **not** know where the jump is? A classic example is the Fourier series of a square wave that is $-1$ on $(-\pi, 0)$ and $+1$ on $(0, \pi)$:

$$
\text{square wave} = \frac{4}{\pi} \sum_{k = 1, 3, 5, \dots} \frac{\sin kx}{k}.
$$

Take more and more terms. Away from the jump the sum converges. But next to the jump it **overshoots**, and the size of the overshoot does not go to zero. The wiggles only move closer to the jump. The peak height tends to about $1.179$, so the overshoot is about $9\%$ of the jump size $2$. The experiment measures this. Going above the true value is called **overshoot**. Going below is **undershoot**.

The lesson from the lecture: **discontinuities are always the hardest part of these problems.** If the solution has one, we should know it and build that knowledge into the scheme. If we ignore it, we meet some version of the Gibbs effect. High-degree polynomial interpolation shows the same problem: linear interpolation is fine near a jump, but high degrees oscillate. (The wiggles of the centered difference in [Lecture 5](/blog/2026-10-20-numerical-pde-lecture-5/) have a different cause, but they look similar: oscillations near a sharp change.)

## 11. A first definition of stability

We saw that some problems force a very small $\Delta t$ and some do not (Lecture 6, part 5 of the experiment). This is the **stability** problem. The lecture introduced the definition with a familiar example: the linear system of the boundary value problem from [Lecture 3](/blog/2026-10-13-numerical-pde-lecture-3/).

Let $U^h$ be the numerical solution and $\hat U^h$ the exact solution at the grid points. Then

$$
A^h U^h = F^h, \qquad A^h \hat U^h = F^h + \tau^h,
$$

where $\tau^h$ is the local truncation error. Subtracting gives an equation for the error $E^h = U^h - \hat U^h$:

$$
A^h E^h = -\tau^h
\qquad\Longrightarrow\qquad
E^h = -(A^h)^{-1} \tau^h
\qquad\Longrightarrow\qquad
\|E^h\| \le \|(A^h)^{-1}\|\, \|\tau^h\|.
$$

Suppose $\|\tau^h\| = O(h^p)$. That is good news, but the other factor can ruin it. If $\|(A^h)^{-1}\|$ is huge, a tiny truncation error can still give a big error. The bound is only a worst case. The error is not always that big. But we cannot hope the worst case never happens. A method should work for every problem in its class, not only on lucky days.

**Definition (stability).** The method is **stable** if there is a constant $C$, independent of $h$, with

$$
\|(A^h)^{-1}\| \le C \qquad \text{for all small enough } h.
$$

Then $\|E^h\| \le C\, \|\tau^h\| = O(h^p)$. The order of the truncation error becomes the order of the global error. In short:

$$
\text{consistency } (\|\tau^h\| \to 0) \;+\; \text{stability} \;\Longrightarrow\; \text{convergence}.
$$

Stability is the main topic for the rest of the course. Next, the same idea will be applied to time-stepping schemes for ODEs.

### Stable is not the same as well-conditioned

The lecture ended with a mention of the **condition number**, $\kappa(A) = \|A\|\,\|A^{-1}\|$. The two ideas are easy to mix up, so here is my own check. The definition of stability uses only $\|A^{-1}\|$. For the 1D Dirichlet matrix $A = \frac{1}{h^2}\,\mathrm{tridiag}(1, -2, 1)$, the eigenvalues are

$$
\lambda_p = \frac{2}{h^2}\big(\cos(p\pi h) - 1\big), \qquad p = 1, \dots, m.
$$

The smallest one in size is close to $-\pi^2$ for every $h$. So $\|A^{-1}\|_2 \approx 1/\pi^2$ stays bounded: the method is **stable**. But the largest one is close to $-4/h^2$, so $\kappa_2(A) \approx 4/(\pi^2 h^2)$ grows without bound. The matrix is **stable but badly conditioned**. Stability controls the discretization error. The condition number controls how much rounding errors can grow when we solve the system. The experiment shows both numbers.

## 12. Experiment

The main test problem is the logistic equation (my own choice):

$$
u' = u(1 - u), \qquad u(0) = 0.2, \qquad u(t) = \frac{1}{1 + 4e^{-t}},
$$

integrated to $T = 2$, where $u(2) \approx 0.648786$. The script uses the Python standard library only. Leapfrog and AB2 get their starting value from RK2, and AB4 gets its three starting values from RK4, unless I say otherwise.

**Part 1: small checks.**

- One step of $u' = -2u$ from $1$ with $\Delta t = 0.5$: Forward Euler $0$, RK2 $0.5$, RK4 $0.375$, exact $0.367879$.
- With $f = e^t$ (no $u$) on $[0, 1]$, one RK4 step gives $1.718861151876593$, and Simpson's rule gives exactly the same digits. The true area is $e - 1 = 1.718281828459045$.
- For $u' = \lambda u$ with several values of $z = \lambda \Delta t$ between $-3$ and $1.5$, the RK2 and RK4 growth factors match the Taylor polynomials of $e^z$ to within $2.2 \times 10^{-16}$, which is rounding.

**Part 2: global error at $T = 2$.** The chart shows four methods for $\Delta t = 1/5, 1/10, \dots, 1/2560$.

<figure class="viz">
<div class="viz-legend"><span class="viz-key"><span class="viz-swatch" style="background:#2a78d6"></span>Forward Euler</span><span class="viz-key"><span class="viz-swatch" style="background:#eb6834"></span>RK2</span><span class="viz-key"><span class="viz-swatch" style="background:#1baf7a"></span>RK4</span><span class="viz-key"><span class="viz-swatch" style="background:#eda100"></span>AB4</span></div>
<div class="viz-frame">
<svg viewBox="0 0 680 400" role="img" tabindex="0" aria-label="Error at t = 2 for the logistic equation, as the time step gets smaller, on log scales. Forward Euler falls with slope 1, RK2 with slope 2, RK4 and AB4 with slope 4. RK4 and AB4 flatten out near 1e-15, where rounding takes over." class="viz-plot"><line x1="64" x2="520" y1="50.7" y2="50.7" class="viz-grid"/><text x="54" y="54.7" text-anchor="end" class="viz-tick">10<tspan dy="-0.5em" font-size="75%">−2</tspan></text><line x1="64" x2="520" y1="92.0" y2="92.0" class="viz-grid"/><text x="54" y="96.0" text-anchor="end" class="viz-tick">10<tspan dy="-0.5em" font-size="75%">−4</tspan></text><line x1="64" x2="520" y1="133.3" y2="133.3" class="viz-grid"/><text x="54" y="137.3" text-anchor="end" class="viz-tick">10<tspan dy="-0.5em" font-size="75%">−6</tspan></text><line x1="64" x2="520" y1="174.7" y2="174.7" class="viz-grid"/><text x="54" y="178.7" text-anchor="end" class="viz-tick">10<tspan dy="-0.5em" font-size="75%">−8</tspan></text><line x1="64" x2="520" y1="216.0" y2="216.0" class="viz-grid"/><text x="54" y="220.0" text-anchor="end" class="viz-tick">10<tspan dy="-0.5em" font-size="75%">−10</tspan></text><line x1="64" x2="520" y1="257.3" y2="257.3" class="viz-grid"/><text x="54" y="261.3" text-anchor="end" class="viz-tick">10<tspan dy="-0.5em" font-size="75%">−12</tspan></text><line x1="64" x2="520" y1="298.7" y2="298.7" class="viz-grid"/><text x="54" y="302.7" text-anchor="end" class="viz-tick">10<tspan dy="-0.5em" font-size="75%">−14</tspan></text><line x1="64" x2="520" y1="340.0" y2="340.0" class="viz-grid"/><text x="54" y="344.0" text-anchor="end" class="viz-tick">10<tspan dy="-0.5em" font-size="75%">−16</tspan></text><text x="64.0" y="362" text-anchor="middle" class="viz-tick">1/5</text><text x="216.0" y="362" text-anchor="middle" class="viz-tick">1/40</text><text x="368.0" y="362" text-anchor="middle" class="viz-tick">1/320</text><text x="520.0" y="362" text-anchor="middle" class="viz-tick">1/2560</text><text x="292.0" y="388" text-anchor="middle" class="viz-axis-title">Time step Δt (smaller to the right)</text><text x="64" y="14" class="viz-axis-title">Error</text><path d="M64.0,51.8 L114.7,58.4 L165.3,64.8 L216.0,71.2 L266.7,77.4 L317.3,83.7 L368.0,89.9 L418.7,96.1 L469.3,102.4 L520.0,108.6" fill="none" stroke="#2a78d6" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><circle cx="64.0" cy="51.8" r="3" fill="#2a78d6" stroke="#ffffff" stroke-width="1.5"/><circle cx="114.7" cy="58.4" r="3" fill="#2a78d6" stroke="#ffffff" stroke-width="1.5"/><circle cx="165.3" cy="64.8" r="3" fill="#2a78d6" stroke="#ffffff" stroke-width="1.5"/><circle cx="216.0" cy="71.2" r="3" fill="#2a78d6" stroke="#ffffff" stroke-width="1.5"/><circle cx="266.7" cy="77.4" r="3" fill="#2a78d6" stroke="#ffffff" stroke-width="1.5"/><circle cx="317.3" cy="83.7" r="3" fill="#2a78d6" stroke="#ffffff" stroke-width="1.5"/><circle cx="368.0" cy="89.9" r="3" fill="#2a78d6" stroke="#ffffff" stroke-width="1.5"/><circle cx="418.7" cy="96.1" r="3" fill="#2a78d6" stroke="#ffffff" stroke-width="1.5"/><circle cx="469.3" cy="102.4" r="3" fill="#2a78d6" stroke="#ffffff" stroke-width="1.5"/><circle cx="520.0" cy="108.6" r="3" fill="#2a78d6" stroke="#ffffff" stroke-width="1.5"/><path d="M64.0,101.4 L114.7,114.7 L165.3,127.5 L216.0,140.1 L266.7,152.6 L317.3,165.1 L368.0,177.5 L418.7,190.0 L469.3,202.4 L520.0,214.9" fill="none" stroke="#eb6834" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><circle cx="64.0" cy="101.4" r="3" fill="#eb6834" stroke="#ffffff" stroke-width="1.5"/><circle cx="114.7" cy="114.7" r="3" fill="#eb6834" stroke="#ffffff" stroke-width="1.5"/><circle cx="165.3" cy="127.5" r="3" fill="#eb6834" stroke="#ffffff" stroke-width="1.5"/><circle cx="216.0" cy="140.1" r="3" fill="#eb6834" stroke="#ffffff" stroke-width="1.5"/><circle cx="266.7" cy="152.6" r="3" fill="#eb6834" stroke="#ffffff" stroke-width="1.5"/><circle cx="317.3" cy="165.1" r="3" fill="#eb6834" stroke="#ffffff" stroke-width="1.5"/><circle cx="368.0" cy="177.5" r="3" fill="#eb6834" stroke="#ffffff" stroke-width="1.5"/><circle cx="418.7" cy="190.0" r="3" fill="#eb6834" stroke="#ffffff" stroke-width="1.5"/><circle cx="469.3" cy="202.4" r="3" fill="#eb6834" stroke="#ffffff" stroke-width="1.5"/><circle cx="520.0" cy="214.9" r="3" fill="#eb6834" stroke="#ffffff" stroke-width="1.5"/><path d="M64.0,135.5 L114.7,160.3 L165.3,185.1 L216.0,209.9 L266.7,234.8 L317.3,259.7 L368.0,284.6 L418.7,312.2 L469.3,329.2 L520.0,316.0" fill="none" stroke="#1baf7a" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><circle cx="64.0" cy="135.5" r="3" fill="#1baf7a" stroke="#ffffff" stroke-width="1.5"/><circle cx="114.7" cy="160.3" r="3" fill="#1baf7a" stroke="#ffffff" stroke-width="1.5"/><circle cx="165.3" cy="185.1" r="3" fill="#1baf7a" stroke="#ffffff" stroke-width="1.5"/><circle cx="216.0" cy="209.9" r="3" fill="#1baf7a" stroke="#ffffff" stroke-width="1.5"/><circle cx="266.7" cy="234.8" r="3" fill="#1baf7a" stroke="#ffffff" stroke-width="1.5"/><circle cx="317.3" cy="259.7" r="3" fill="#1baf7a" stroke="#ffffff" stroke-width="1.5"/><circle cx="368.0" cy="284.6" r="3" fill="#1baf7a" stroke="#ffffff" stroke-width="1.5"/><circle cx="418.7" cy="312.2" r="3" fill="#1baf7a" stroke="#ffffff" stroke-width="1.5"/><circle cx="469.3" cy="329.2" r="3" fill="#1baf7a" stroke="#ffffff" stroke-width="1.5"/><circle cx="520.0" cy="316.0" r="3" fill="#1baf7a" stroke="#ffffff" stroke-width="1.5"/><path d="M64.0,92.3 L114.7,116.6 L165.3,141.5 L216.0,166.5 L266.7,191.4 L317.3,216.4 L368.0,241.3 L418.7,266.2 L469.3,290.9 L520.0,309.5" fill="none" stroke="#eda100" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><circle cx="64.0" cy="92.3" r="3" fill="#eda100" stroke="#ffffff" stroke-width="1.5"/><circle cx="114.7" cy="116.6" r="3" fill="#eda100" stroke="#ffffff" stroke-width="1.5"/><circle cx="165.3" cy="141.5" r="3" fill="#eda100" stroke="#ffffff" stroke-width="1.5"/><circle cx="216.0" cy="166.5" r="3" fill="#eda100" stroke="#ffffff" stroke-width="1.5"/><circle cx="266.7" cy="191.4" r="3" fill="#eda100" stroke="#ffffff" stroke-width="1.5"/><circle cx="317.3" cy="216.4" r="3" fill="#eda100" stroke="#ffffff" stroke-width="1.5"/><circle cx="368.0" cy="241.3" r="3" fill="#eda100" stroke="#ffffff" stroke-width="1.5"/><circle cx="418.7" cy="266.2" r="3" fill="#eda100" stroke="#ffffff" stroke-width="1.5"/><circle cx="469.3" cy="290.9" r="3" fill="#eda100" stroke="#ffffff" stroke-width="1.5"/><circle cx="520.0" cy="309.5" r="3" fill="#eda100" stroke="#ffffff" stroke-width="1.5"/><text x="530" y="112.6" class="viz-label">Forward Euler</text><text x="530" y="218.9" class="viz-label">RK2</text><text x="530" y="330.0" class="viz-label">RK4</text><text x="530" y="303.5" class="viz-label">AB4</text><line x1="0" x2="0" y1="30" y2="340" class="viz-crosshair" visibility="hidden"/></svg>
<div class="viz-tooltip" hidden></div>
</div>
<figcaption>Error |U<sup>N</sup> − u(2)| on log scales for the logistic equation u′ = u(1 − u), u(0) = 0.2. AB4 starts with three RK4 steps. Below about 10<sup>−14</sup> the error is rounding. Hover over the chart, or focus it and use the arrow keys, to read the values.</figcaption>
<details class="viz-table"><summary>Show the numbers</summary><table><thead><tr><th>Δt</th><th>Forward Euler</th><th>RK2</th><th>RK4</th><th>AB4</th></tr></thead><tbody><tr><td>1/5</td><td>8.9e-3</td><td>3.5e-5</td><td>7.8e-7</td><td>9.7e-5</td></tr><tr><td>1/10</td><td>4.2e-3</td><td>8.0e-6</td><td>5.0e-8</td><td>6.5e-6</td></tr><tr><td>1/20</td><td>2.1e-3</td><td>1.9e-6</td><td>3.1e-9</td><td>4.0e-7</td></tr><tr><td>1/40</td><td>1.0e-3</td><td>4.7e-7</td><td>2.0e-10</td><td>2.5e-8</td></tr><tr><td>1/80</td><td>5.1e-4</td><td>1.2e-7</td><td>1.2e-11</td><td>1.5e-9</td></tr><tr><td>1/160</td><td>2.5e-4</td><td>2.9e-8</td><td>7.7e-13</td><td>9.6e-11</td></tr><tr><td>1/320</td><td>1.3e-4</td><td>7.3e-9</td><td>4.8e-14</td><td>6.0e-12</td></tr><tr><td>1/640</td><td>6.3e-5</td><td>1.8e-9</td><td>2.2e-15</td><td>3.7e-13</td></tr><tr><td>1/1280</td><td>3.1e-5</td><td>4.5e-10</td><td>3.3e-16</td><td>2.4e-14</td></tr><tr><td>1/2560</td><td>1.6e-5</td><td>1.1e-10</td><td>1.4e-15</td><td>3.0e-15</td></tr></tbody></table></details>
<script>(() => {const fig = document.currentScript.closest('figure');const svg = fig.querySelector('svg');const line = svg.querySelector('.viz-crosshair');const tip = fig.querySelector('.viz-tooltip');const names = ["Forward Euler", "RK2", "RK4", "AB4"]; const colors = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"];const rows = [["8.9e-3", "3.5e-5", "7.8e-7", "9.7e-5"], ["4.2e-3", "8.0e-6", "5.0e-8", "6.5e-6"], ["2.1e-3", "1.9e-6", "3.1e-9", "4.0e-7"], ["1.0e-3", "4.7e-7", "2.0e-10", "2.5e-8"], ["5.1e-4", "1.2e-7", "1.2e-11", "1.5e-9"], ["2.5e-4", "2.9e-8", "7.7e-13", "9.6e-11"], ["1.3e-4", "7.3e-9", "4.8e-14", "6.0e-12"], ["6.3e-5", "1.8e-9", "2.2e-15", "3.7e-13"], ["3.1e-5", "4.5e-10", "3.3e-16", "2.4e-14"], ["1.6e-5", "1.1e-10", "1.4e-15", "3.0e-15"]]; const heads = ["\u0394t = 1/5", "\u0394t = 1/10", "\u0394t = 1/20", "\u0394t = 1/40", "\u0394t = 1/80", "\u0394t = 1/160", "\u0394t = 1/320", "\u0394t = 1/640", "\u0394t = 1/1280", "\u0394t = 1/2560"];const left = 64, right = 520, steps = 9, width = 680;let current = -1;const show = (n) => { n = Math.max(0, Math.min(steps, n)); current = n; const x = left + (right - left) * n / steps; line.setAttribute('x1', x); line.setAttribute('x2', x); line.setAttribute('visibility', 'visible'); tip.replaceChildren(); const head = document.createElement('div'); head.className = 'viz-tooltip-head'; head.textContent = heads[n]; tip.append(head); names.forEach((name, i) => {  const row = document.createElement('div'); row.className = 'viz-tooltip-row';  const key = document.createElement('span'); key.className = 'viz-swatch'; key.style.background = colors[i];  const value = document.createElement('strong'); value.textContent = rows[n][i];  const label = document.createElement('span'); label.textContent = name;  row.append(key, value, label); tip.append(row); }); const box = svg.getBoundingClientRect(), frame = fig.getBoundingClientRect(); const xs = box.left - frame.left + x * box.width / width; tip.hidden = false; const room = frame.width - tip.offsetWidth - 8; tip.style.left = Math.max(0, Math.min(room, xs + 12)) + 'px'; tip.style.top = (box.top - frame.top + 8) + 'px';};const hide = () => { line.setAttribute('visibility', 'hidden'); tip.hidden = true; current = -1; };svg.addEventListener('pointermove', (event) => { const box = svg.getBoundingClientRect(); const x = (event.clientX - box.left) * width / box.width; show(Math.round((x - left) / (right - left) * steps));});svg.addEventListener('pointerleave', hide);svg.addEventListener('focus', () => show(current < 0 ? 0 : current));svg.addEventListener('blur', hide);svg.addEventListener('keydown', (event) => { if (event.key === 'ArrowRight') { show(current + 1); event.preventDefault(); } if (event.key === 'ArrowLeft') { show(current - 1); event.preventDefault(); } if (event.key === 'Escape') hide();});})();</script>
</figure>

| $\Delta t$ | Forward Euler | RK2 | RK4 | Leapfrog | AB2 | AB4 |
|---|---|---|---|---|---|---|
| 1/5 | $8.86 \times 10^{-3}$ | $3.51 \times 10^{-5}$ | $7.84 \times 10^{-7}$ | $1.04 \times 10^{-3}$ | $2.38 \times 10^{-3}$ | $9.69 \times 10^{-5}$ |
| 1/40 | $1.02 \times 10^{-3}$ | $4.71 \times 10^{-7}$ | $1.97 \times 10^{-10}$ | $1.64 \times 10^{-5}$ | $4.06 \times 10^{-5}$ | $2.49 \times 10^{-8}$ |
| 1/320 | $1.26 \times 10^{-4}$ | $7.26 \times 10^{-9}$ | $4.81 \times 10^{-14}$ | $2.57 \times 10^{-7}$ | $6.41 \times 10^{-7}$ | $5.99 \times 10^{-12}$ |
| 1/2560 | $1.57 \times 10^{-5}$ | $1.13 \times 10^{-10}$ | $1.44 \times 10^{-15}$ | $4.01 \times 10^{-9}$ | $1.00 \times 10^{-8}$ | $3.00 \times 10^{-15}$ |
| order (1/20 → 1/320) | 1.01 | 2.01 | 4.00 | 2.00 | 1.99 | 4.01 |

What I noticed:

- **The orders match the theory**: $1$ for Forward Euler, $2$ for RK2, leapfrog and AB2, $4$ for RK4 and AB4.
- **RK4 hits rounding first.** Its error is about $10^{-15}$ from $\Delta t = 1/640$ on. After that it no longer falls, and the measured "order" jumps around. AB4 gets there at $\Delta t = 1/2560$.
- **Same order, different constants.** With the same $\Delta t$, RK4 is about $125$ times more accurate than AB4 (for example $1.97 \times 10^{-10}$ against $2.49 \times 10^{-8}$ at $\Delta t = 1/40$). Among the second-order methods, RK2 is about $35$ times more accurate than leapfrog on this problem. Order says how fast the error falls, not how small it is.

**Part 3: starting values.**

| $\Delta t$ | Leapfrog, FE start | Leapfrog, RK2 start | Leapfrog, exact start | AB4, FE start | AB4, RK4 start |
|---|---|---|---|---|---|
| 1/5 | $3.63 \times 10^{-4}$ | $1.04 \times 10^{-3}$ | $1.06 \times 10^{-3}$ | $6.65 \times 10^{-3}$ | $9.69 \times 10^{-5}$ |
| 1/40 | $5.61 \times 10^{-6}$ | $1.64 \times 10^{-5}$ | $1.64 \times 10^{-5}$ | $1.25 \times 10^{-4}$ | $2.49 \times 10^{-8}$ |
| 1/320 | $8.77 \times 10^{-8}$ | $2.57 \times 10^{-7}$ | $2.57 \times 10^{-7}$ | $2.00 \times 10^{-6}$ | $5.99 \times 10^{-12}$ |
| 1/2560 | $1.37 \times 10^{-9}$ | $4.01 \times 10^{-9}$ | $4.01 \times 10^{-9}$ | $3.13 \times 10^{-8}$ | $3.00 \times 10^{-15}$ |
| order (1/20 → 1/320) | 2.00 | 2.00 | 2.00 | 1.98 | 4.01 |

- **Leapfrog keeps second order with any of the three starts**, even with Forward Euler. This confirms the rule from section 7: the start is used once, so its one-step error $O(\Delta t^2)$ is enough.
- On this problem the Forward Euler start even gives a **smaller** error, by a factor of about $3$. I do not think this is a real advantage. At $\Delta t = 1/160$ the Forward Euler value $U^1$ is too low, while the leapfrog error itself is positive, so the two partly cancel. The exact start and the RK2 start give almost the same error, so the RK2 start error is tiny.
- **AB4 loses two orders with a Forward Euler start.** The order drops from $4.01$ to $1.98$, and at $\Delta t = 1/2560$ the error is $3.13 \times 10^{-8}$ instead of $3.00 \times 10^{-15}$. Here a cheap start really ruins the method.

**Part 4: cost in evaluations of $f$.** I counted every call to $f$. RK4 uses $4N$. AB4 uses $N + 13$: twelve for the three RK4 start steps, and then one new value per step.

| steps $N$ | RK4: evaluations | RK4: error | AB4: evaluations | AB4: error |
|---|---|---|---|---|
| 20 | 80 | $4.98 \times 10^{-8}$ | 33 | $6.45 \times 10^{-6}$ |
| 40 | 160 | $3.14 \times 10^{-9}$ | 53 | $4.02 \times 10^{-7}$ |
| 80 | 320 | $1.97 \times 10^{-10}$ | 93 | $2.49 \times 10^{-8}$ |
| 160 | 640 | $1.23 \times 10^{-11}$ | 173 | $1.54 \times 10^{-9}$ |
| 320 | 1280 | $7.71 \times 10^{-13}$ | 333 | $9.60 \times 10^{-11}$ |
| 640 | 2560 | $4.81 \times 10^{-14}$ | 653 | $5.99 \times 10^{-12}$ |

Compare rows with about the same work. RK4 with $N = 80$ uses $320$ evaluations and reaches $1.97 \times 10^{-10}$. AB4 with $N = 320$ uses $333$ evaluations and reaches $9.60 \times 10^{-11}$, about half the error. The same happens one row later ($1.23 \times 10^{-11}$ against $5.99 \times 10^{-12}$). So with the same budget of $f$ evaluations, AB4 wins here, even though RK4 is more accurate per step. This is the trade-off from section 9. (For this cheap $f$ the difference does not matter. It matters when one evaluation of $f$ is expensive, for example a whole PDE right-hand side.)

**Part 5: the Gibbs phenomenon.** Partial sums of the square-wave series with $n$ terms ($k = 1, 3, \dots, 2n - 1$):

| terms $n$ | peak height | peak position | error at $x = \pi/2$ |
|---|---|---|---|
| 4 | 1.1842 | 0.39270 | $7.84 \times 10^{-2}$ |
| 16 | 1.1793 | 0.09817 | $1.99 \times 10^{-2}$ |
| 64 | 1.1790 | 0.02454 | $4.97 \times 10^{-3}$ |
| 256 | 1.1790 | 0.00614 | $1.24 \times 10^{-3}$ |

The peak sits exactly at $x = \pi/(2n)$, so it moves toward the jump as $n$ grows. Away from the jump, at $x = \pi/2$, the error falls like $1/n$. But the peak height does not fall. It tends to $\frac{2}{\pi}\int_0^\pi \frac{\sin s}{s}\,ds \approx 1.178980$, an overshoot of about $8.9\%$ of the jump.

**Part 6: stable but badly conditioned.** For the 1D Dirichlet matrix from section 11:

| $h$ | $\lVert A^{-1} \rVert_2$ | $\lVert A^{-1} \rVert_\infty$ | $\kappa_2(A)$ |
|---|---|---|---|
| 1/8 | 0.102633 | 0.125000 | $2.53 \times 10^{1}$ |
| 1/64 | 0.101342 | 0.125000 | $1.66 \times 10^{3}$ |
| 1/512 | 0.101322 | 0.125000 | $1.06 \times 10^{5}$ |
| 1/4096 | 0.101321 | 0.125000 | $6.80 \times 10^{6}$ |

$\|A^{-1}\|_2$ settles at $1/\pi^2 \approx 0.101321$, and $\|A^{-1}\|_\infty$ is exactly $1/8$ for every $h$. (All entries of $A^{-1}$ have one sign, so the $\infty$-norm is the largest entry of $A^{-1}\mathbf{1}$. That vector is the discrete solution of $u'' = 1$ with zero boundary values, whose largest size is $1/8$.) Both norms stay bounded, so the method is stable. Meanwhile $\kappa_2(A)$ grows by a factor of $64$ each time $h$ shrinks by $8$, like $h^{-2}$. At $h = 1/4096$ it is $6.80 \times 10^{6}$, which equals $4/(\pi^2 h^2)$ to three digits.

The full script and the chart code are in [this site's repository](https://github.com/doubleyands/doubleyands.github.io/tree/main/scripts/posts/numerical-pde/lecture-7).

## Summary

1. A time-stepping scheme is an area rule. The **midpoint** rule gives second order with one point, but the middle value is not on the grid.
2. **RK2** fills the gap with a half Euler step. It is explicit, one-step, two-stage, and second order.
3. **RK4** follows Simpson's rule ($1 : 4 : 1$). It is explicit, four-stage, and fourth order. If $f$ depends only on $t$, it **is** Simpson's rule.
4. Runge–Kutta methods reach high order by **nesting evaluations** of $f$, not by differentiating $f$.
5. **Linear multistep methods** reuse old values. They are explicit if $\beta_r = 0$. **Leapfrog** is the centered difference in time.
6. Multistep methods are **not self-starting**. A method of order $p$ needs starting values with error $O(\Delta t^p)$: a Forward Euler start is fine for leapfrog but drops AB4 to second order.
7. **Adams–Bashforth** (explicit, order $r$) and **Adams–Moulton** (implicit, order $r + 1$) add an area to the last value. Their coefficients must pass the $f \equiv 1$ check.
8. Multistep methods need only **one new $f$ evaluation** per step, but they cannot start alone, change $\Delta t$ easily, or handle a known kink well.
9. Smooth approximations of a jump **overshoot** by about $9\%$, however many terms we take (**Gibbs**). Build known discontinuities into the scheme.
10. A method is **stable** if $\|(A^h)^{-1}\|$ is bounded independently of $h$. **Consistency + stability ⇒ convergence**. Stability is not the same as a small condition number.
