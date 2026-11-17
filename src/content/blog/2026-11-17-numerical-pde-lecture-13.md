---
title: "Numerical PDEs, Lecture 13: von Neumann Analysis for Advection, and a First Look at Finite Volumes"
description: "My notes on the thirteenth lecture: the amplification factors of upwind, Lax–Friedrichs, Lax–Wendroff and leapfrog, why the amount of artificial diffusion matters, and the first steps of the finite volume method: integral form, cell averages, numerical flux and conservation form."
pubDate: 2026-11-17
project: numerical-pde
tags: [numerical PDEs, advection, von Neumann analysis, finite volume, conservation laws, notes]
---

This is my summary of the thirteenth lecture in my numerical PDE course. I write the ideas in my own words and check them with a small experiment.

The lecture has two halves. The first half finishes the advection schemes. [Lecture 12](/blog/2026-11-12-numerical-pde-lecture-12/) showed that Lax–Friedrichs, Lax–Wendroff and upwind all add **artificial diffusion**, and found their stability with the eigenvalues of the space matrix. Here we get the same answers again with **von Neumann analysis** ([Lecture 11](/blog/2026-11-10-numerical-pde-lecture-11/)). The second half starts a new topic: the **finite volume method**.

Notation as before: $h$ is the grid spacing, $k$ is the time step, $U^n_j$ is the value at $x_j = jh$ and $t_n = nk$. In exponents, $i = \sqrt{-1}$ and $j$ is the grid index. We solve $u_t + a\,u_x = 0$ with $a > 0$ unless I say otherwise.

## 1. The amount of artificial diffusion matters

Recall the three fixes for FTCS:

| scheme | the idea | artificial diffusion $\epsilon$ |
|---|---|---|
| Lax–Friedrichs | replace $U_j$ by the average of its neighbors | $h^2/(2k)$ |
| upwind | use only the side the flow comes from | $ah/2$ |
| Lax–Wendroff | Taylor expansion in time, with $u_{tt} = a^2u_{xx}$ | $a^2k/2$ |

All three solve, in effect, the **modified equation** $u_t + a\,u_x = \epsilon\,u_{xx}$ instead of the true one. With $\nu = ak/h$, the three values are linked: $\epsilon_{\text{LW}} = \nu\,\epsilon_{\text{up}}$ and $\epsilon_{\text{up}} = \nu\,\epsilon_{\text{LF}}$. Lecture 12 showed that any $\epsilon$ between $a^2k/2$ and $h^2/(2k)$ gives a stable scheme.

So the diffusion is not a bug. Hyperbolic schemes almost always contain some of it. The real question is **how much**. Too little, and the scheme blows up. Too much, and the solution gets flat and wrong.

**Consistency needs $\epsilon \to 0$.** The term $\epsilon\,u_{xx}$ is not in the real equation. So when $h, k \to 0$, we need $\epsilon \to 0$ too. Otherwise the scheme converges to the solution of a **different problem**. For Lax–Friedrichs, $\epsilon = h^2/(2k)$, so this depends on how $k$ and $h$ go to zero together:

- $k/h$ fixed: $\epsilon = ah/(2\nu) \to 0$. Good.
- $k = Ch^2$: $\epsilon = 1/(2C)$ stays fixed. The scheme converges, but to the solution of an advection–diffusion equation.
- $k$ much smaller, like $k = Ch^3$: $\epsilon = 1/(2Ch)$ grows without bound.

This is the opposite of FTCS in [Lecture 12](/blog/2026-11-12-numerical-pde-lecture-12/), where a very small $k$ (like $k = h^2$) was the way to make it converge. Part 2 of the experiment runs all three cases.

## 2. von Neumann analysis for advection

The recipe is the same as for the heat equation in Lecture 11. Put $U^n_j = g(\xi)^n e^{i\xi jh}$ into the scheme, cancel $g^n e^{i\xi jh}$, and ask for $|g(\xi)| \le 1$ for all $\xi$. I write $\theta = \xi h$, which runs over $[-\pi, \pi]$.

All results depend on the **Courant number**

$$
\nu = \frac{ak}{h}.
$$

It is speed times time divided by length, so it has no units. It is one of the most important dimensionless numbers in this subject.

### Upwind

$$
U^{n+1}_j = U^n_j - \nu\big(U^n_j - U^n_{j-1}\big) \quad\Longrightarrow\quad g(\theta) = 1 - \nu\big(1 - e^{-i\theta}\big) = (1 - \nu) + \nu\,e^{-i\theta}.
$$

As $\theta$ moves, $g$ runs around a **circle with center $1 - \nu$ and radius $|\nu|$**. This circle stays inside the unit circle exactly when $0 \le \nu \le 1$:

- **Necessary.** At $\theta = \pi$ (the sawtooth mode), $g = 1 - 2\nu$. We need $|1 - 2\nu| \le 1$, so $0 \le \nu \le 1$.
- **Sufficient.** If $0 \le \nu \le 1$, then $|g| \le |1 - \nu| + |\nu\,e^{-i\theta}| = (1 - \nu) + \nu = 1$.

There is also a one-line formula. Multiply $g$ by its conjugate:

$$
|g|^2 = 1 - 2\nu(1 - \nu)(1 - \cos\theta).
$$

Since $1 - \cos\theta \ge 0$, this is at most $1$ for all $\theta$ exactly when $\nu(1 - \nu) \ge 0$.

### Lax–Friedrichs

$$
U^{n+1}_j = \tfrac12\big(U^n_{j-1} + U^n_{j+1}\big) - \tfrac{\nu}{2}\big(U^n_{j+1} - U^n_{j-1}\big)
\quad\Longrightarrow\quad
g(\theta) = \cos\theta - i\nu\sin\theta,
$$

$$
|g|^2 = \cos^2\theta + \nu^2\sin^2\theta = 1 - (1 - \nu^2)\sin^2\theta.
$$

This is at most $1$ for all $\theta$ exactly when $\nu^2 \le 1$.

### Lax–Wendroff

$$
g(\theta) = 1 - i\nu\sin\theta + \nu^2(\cos\theta - 1), \qquad |g|^2 = 1 - \nu^2(1 - \nu^2)(1 - \cos\theta)^2.
$$

To check the second formula, write $w = 1 - \cos\theta$ and $\sin^2\theta = w(2 - w)$, and expand $(1 - \nu^2w)^2 + \nu^2w(2 - w)$. Again the condition is $|\nu| \le 1$.

### Leapfrog

Leapfrog uses three time levels, so the ansatz gives a **quadratic** equation for $g$:

$$
g^2 + 2i\nu\sin\theta\;g - 1 = 0 \quad\Longrightarrow\quad g = -i\nu\sin\theta \pm \sqrt{1 - \nu^2\sin^2\theta}.
$$

- If $|\nu| < 1$, the square root is real. Then $|g|^2 = \nu^2\sin^2\theta + (1 - \nu^2\sin^2\theta) = 1$ for **both** roots. Nothing grows and nothing is damped: the scheme is nondissipative, as in Lecture 12.
- If $|\nu| = 1$, at $\theta = \pi/2$ the square root is zero and the two roots meet at $g = -i$. A double root of size $1$ makes the solution grow like $n$. This is the root condition from [Lecture 9](/blog/2026-11-03-numerical-pde-lecture-9/). So leapfrog needs the strict inequality $|\nu| < 1$.

### Why the answers match the eigenvalue method

Upwind, Lax–Friedrichs and Lax–Wendroff are all Forward Euler applied to $U' = A_\epsilon U$ ([Lecture 12](/blog/2026-11-12-numerical-pde-lecture-12/)). The Fourier mode $e^{i\theta j}$ is an eigenvector of $A_\epsilon$, so one step multiplies it by exactly

$$
g(\theta) = 1 + k\lambda(\theta).
$$

So the curve of $g$ is the eigenvalue ellipse of Lecture 12, shifted right by $1$. The condition $|g| \le 1$ is the same as $|1 + k\lambda| \le 1$, the Forward Euler disk. von Neumann analysis does the eigenvalue computation and the absolute stability check in one step.

### What the curves say

- For **small $\theta$** (smooth modes), the damping per step is about $\tfrac12(1 - \nu^2)\theta^2$ for Lax–Friedrichs, $\nu(1 - \nu)\tfrac{\theta^2}{2}$ for upwind, and only of order $\theta^4$ for Lax–Wendroff. This is the order of their $\epsilon$: Lax–Friedrichs smears the most.
- At **$\theta = \pi$**, Lax–Friedrichs has $|g| = 1$. It never damps the sawtooth mode. The other two do.
- **Leapfrog** has $|g| = 1$ everywhere.

## 3. The end of finite differences

With this, the finite difference part of the course is done: elliptic, parabolic and hyperbolic problems are all covered. (The ideas will come back.) Two topics are left:

- **Finite volume methods.** These get a lot of time. Hyperbolic problems are hard: their solutions can have **shocks** and other **discontinuities**, and the numerics need care.
- **Finite element methods.** These are more mathematical, and the course only gives a short taste.

The main goal is still to understand finite differences well. The rest builds on them.

## 4. Why finite volumes: conservation laws

A finite difference method computes values at **nodes**. If the scheme is stable, we get the node values, and we are done.

Many physical laws are **conservation laws**: mass, energy and momentum are conserved. We would like the numerical scheme to have the same property. A general finite difference scheme does not guarantee this. Finite volume methods are built so that it holds.

A simple picture is a bank account. The balance now equals the balance before, plus the money that came in, minus the money that went out. No money appears or disappears inside the account.

## 5. The integral form

Let $q(x, t)$ be the **density** of some quantity, for example a pollutant in a pipe. The total amount in the piece of pipe between $x_1$ and $x_2$ is

$$
\int_{x_1}^{x_2} q(x, t)\,dx.
$$

(Density is mass per length, so integrating over length gives mass. An integral is a sum.)

The **flux** $F(t)$ at a fixed point is the rate at which the quantity flows past that point. $F > 0$ means flow to the right, and $F < 0$ means flow to the left. The amount in $[x_1, x_2]$ changes only by what flows in at $x_1$ and out at $x_2$:

$$
\frac{d}{dt}\int_{x_1}^{x_2} q(x, t)\,dx = F_1(t) - F_2(t), \qquad F_i(t) = f\big(q(x_i, t)\big).
$$

This is the **integral form** of the conservation law: in minus out. The function $f$ is the **flux function**. For advection with speed $a$ it is $f(q) = aq$: density (amount per meter) times speed (meters per second) gives amount per second.

## 6. The differential form, and its hidden cost

By the fundamental theorem of calculus, $f(q(x_2, t)) - f(q(x_1, t)) = \int_{x_1}^{x_2} \partial_x f(q)\,dx$. The points $x_1, x_2$ are fixed, so we can move $d/dt$ inside the integral:

$$
\int_{x_1}^{x_2}\Big[\,q_t + f(q)_x\,\Big]\,dx = 0 \quad \text{for all } x_1 < x_2.
$$

If the integral of a continuous function is zero over every interval, the function itself is zero. (If it were positive at one point, it would be positive on a small interval around it, and the integral there would be positive.) So

$$
q_t + f(q)_x = 0.
$$

This is the **differential form**. With $f(q) = aq$ it is our advection equation $q_t + a\,q_x = 0$.

The differential form is more familiar and easier to work with. But on the way we **assumed that $q$ is differentiable**. Hyperbolic problems can have discontinuities, and they move along without changing shape (along the characteristics of Lecture 12). At a jump, $q_x$ does not exist. The integral form still makes sense there, because integrals do not need derivatives. Later in the course this leads to **weak solutions**.

## 7. Cells and cell averages

A finite volume method works with pieces of the domain, not points. Split the line into **cells** (or **volumes**)

$$
C_i = \big(x_{i-1/2},\, x_{i+1/2}\big),
$$

with center $x_i$ and length $\Delta x$. The unknown in each cell is the **cell average**:

$$
Q^n_i \approx \frac{1}{\Delta x}\int_{C_i} q(x, t_n)\,dx.
$$

So $Q^n_i$ is (approximately) the amount in the cell divided by its length. A cell is not alone. It talks to its neighbors $C_{i-1}$ and $C_{i+1}$ through the **interfaces** $x_{i-1/2}$ and $x_{i+1/2}$, like doors. The quantity flows in and out only through these doors.

## 8. From an exact update to a scheme: the numerical flux

Write the integral form for one cell $C_i$:

$$
\frac{d}{dt}\int_{C_i} q(x, t)\,dx = f\big(q(x_{i-1/2}, t)\big) - f\big(q(x_{i+1/2}, t)\big).
$$

Integrate from $t_n$ to $t_{n+1}$ and divide by $\Delta x$:

$$
\frac{1}{\Delta x}\int_{C_i} q(x, t_{n+1})\,dx = \frac{1}{\Delta x}\int_{C_i} q(x, t_n)\,dx - \frac{1}{\Delta x}\left[\int_{t_n}^{t_{n+1}} f\big(q(x_{i+1/2}, t)\big)\,dt - \int_{t_n}^{t_{n+1}} f\big(q(x_{i-1/2}, t)\big)\,dt\right].
$$

So far this is **exact**. Nothing is approximated. Now define the average flux through the door $x_{i-1/2}$ during one step:

$$
F^n_{i-1/2} \approx \frac{1}{\Delta t}\int_{t_n}^{t_{n+1}} f\big(q(x_{i-1/2}, t)\big)\,dt.
$$

Then the update reads

$$
\boxed{\;Q^{n+1}_i = Q^n_i - \frac{\Delta t}{\Delta x}\Big(F^n_{i+1/2} - F^n_{i-1/2}\Big)\;}
$$

On the left we have space averages, and on the right time averages of the flux. It looks like a finite difference scheme, forward in time.

**But we cannot compute the right side exactly.** The integrand contains $q$, the unknown solution. So we must **approximate** the flux.

Hyperbolic problems have a **finite propagation speed**: information moves at a finite speed. In one short step, the flux through $x_{i-1/2}$ can only depend on what is near that door. We do not know in advance whether the flow there goes in or out, so the simplest choice is to let it depend on the two cells on both sides:

$$
F^n_{i-1/2} = \mathcal{F}\big(Q^n_{i-1},\, Q^n_i\big).
$$

The function $\mathcal F$ is the **numerical flux**. The scheme becomes

$$
Q^{n+1}_i = Q^n_i - \frac{\Delta t}{\Delta x}\Big[\mathcal{F}\big(Q^n_i, Q^n_{i+1}\big) - \mathcal{F}\big(Q^n_{i-1}, Q^n_i\big)\Big].
$$

This is the general finite volume scheme. To get $Q^{n+1}_i$ it uses $Q^n_{i-1}, Q^n_i, Q^n_{i+1}$: a **3-point stencil**. The rest of the course is mostly about one question: **how to choose $\mathcal F$**. A good choice gives a good scheme, and a bad choice gives a bad one.

For example, with $f(q) = aq$ and $a > 0$, the flow comes from the left, so take $\mathcal F(Q_{i-1}, Q_i) = aQ_{i-1}$. Then the scheme is exactly the upwind scheme of section 2. The obvious "central" choice, the average $\tfrac12\big(f(Q_{i-1}) + f(Q_i)\big)$, gives back FTCS, which we know is unstable. From section 1 we know the fix: add some artificial diffusion. The next lecture builds numerical fluxes in this way.

## 9. Conservation form and the telescoping sum

A scheme written as $Q^{n+1}_i = Q^n_i - \frac{\Delta t}{\Delta x}\big(F^n_{i+1/2} - F^n_{i-1/2}\big)$ is in **conservation form**. Multiply by $\Delta x$ and add the equations for cells $I$ to $J$:

$$
\Delta x\sum_{i=I}^{J} Q^{n+1}_i = \Delta x\sum_{i=I}^{J} Q^n_i - \Delta t\Big(F^n_{J+1/2} - F^n_{I-1/2}\Big).
$$

Each inner flux appears twice: once as "out" for one cell and once as "in" for its neighbor. They cancel, and only the two outer fluxes remain. This is a **telescoping sum**. The result says: total amount at the new time = total amount before + what came in − what went out. It is the discrete version of the conservation law we started with.

On a periodic domain the two outer fluxes are the same flux, so the total amount stays **exactly** constant, up to rounding.

Dividing by $\Delta t$ gives the finite difference look of the same scheme:

$$
\frac{Q^{n+1}_i - Q^n_i}{\Delta t} + \frac{F^n_{i+1/2} - F^n_{i-1/2}}{\Delta x} = 0.
$$

It is like a central difference in space around the cell center $x_i$, but the values that enter are fluxes at the cell walls.

## 10. A 1D example: diffusion

If the flux depends on $q_x$, the conservation law becomes a diffusion equation. With the flux $f = -\beta(x)\,q_x$ (the quantity flows from high to low density), $q_t + f_x = 0$ is

$$
q_t = \big(\beta(x)\,q_x\big)_x.
$$

A natural numerical flux takes the slope at the door from the two neighboring cells:

$$
\mathcal F(Q_{i-1}, Q_i) = -\beta_{i-1/2}\,\frac{Q_i - Q_{i-1}}{\Delta x}, \qquad \beta_{i-1/2} \approx \beta(x_{i-1/2}).
$$

For constant $\beta$, the finite volume scheme becomes

$$
Q^{n+1}_i = Q^n_i + \frac{\beta\,\Delta t}{\Delta x^2}\big(Q^n_{i-1} - 2Q^n_i + Q^n_{i+1}\big),
$$

which is exactly FTCS for the heat equation ([Lecture 10](/blog/2026-11-05-numerical-pde-lecture-10/)), with the same limit $\Delta t = O(\Delta x^2)$. For Crank–Nicolson, use the average of the fluxes at levels $n$ and $n+1$. So the framework covers the heat equation too. But the main interest of finite volumes is hyperbolic problems.

## 11. What makes a numerical flux good

The basic requirement is **convergence**. By the Lax equivalence theorem ([Lecture 10](/blog/2026-11-05-numerical-pde-lecture-10/)), for a linear problem this means two things:

- **consistency**: the numerical flux approximates the true flux well, locally;
- **stability**: small errors made in each step do not grow too fast later.

For hyperbolic problems, choosing the flux is already the hard part, and the solutions can have discontinuities that move. A key stability tool here is the **CFL condition**, which the next lecture treats properly.

## 12. Experiment

The script uses the Python standard library only (`cmath` and `math`). All grids are periodic.

**Part 1: measured amplification factors.** On a grid with $64$ points, I put each mode $e^{i\theta j}$ ($\theta = 2\pi p/64$, $p = 0, \dots, 32$) into one step of each scheme with $\nu = 0.6$, and divided the result by the input. The ratio is the same at every grid point (largest spread $3.9 \times 10^{-14}$), so each mode really is just multiplied by a number $g$. Against the formulas in section 2:

| | upwind | Lax–Friedrichs | Lax–Wendroff |
|---|---|---|---|
| largest $\lvert g - g_{\text{formula}}\rvert$ | $1.4 \times 10^{-14}$ | $1.8 \times 10^{-14}$ | $1.1 \times 10^{-14}$ |
| largest error of the $\lvert g\rvert^2$ formula | $8.6 \times 10^{-15}$ | $7.3 \times 10^{-15}$ | $1.4 \times 10^{-14}$ |
| largest $\lvert g - (1 + k\lambda)\rvert$ | $1.4 \times 10^{-14}$ | $1.8 \times 10^{-14}$ | $1.1 \times 10^{-14}$ |
| $\lvert g\rvert$ at $\theta = \pi/4$ | $0.9270$ | $0.8246$ | $0.9901$ |
| $\lvert g\rvert$ at $\theta = \pi/2$ | $0.7211$ | $0.6000$ | $0.8773$ |
| $\lvert g\rvert$ at $\theta = \pi$ | $0.2000$ | $1.0000$ | $0.2800$ |

Here $1 + k\lambda$ uses the eigenvalues of $A_\epsilon$ from Lecture 12. So $g = 1 + k\lambda$ holds to rounding.

<figure class="viz">
<div class="viz-legend"><span class="viz-key"><span class="viz-swatch" style="background:#2a78d6"></span>upwind</span><span class="viz-key"><span class="viz-swatch" style="background:#eb6834"></span>Lax–Friedrichs</span><span class="viz-key"><span class="viz-swatch" style="background:#1baf7a"></span>Lax–Wendroff</span><span class="viz-key"><span class="viz-swatch" style="background:#eda100"></span>leapfrog</span></div>
<div class="viz-frame">
<svg viewBox="0 0 680 400" role="img" tabindex="0" aria-label="Size of the amplification factor g against the grid frequency theta from 0 to pi, at Courant number 0.6. Leapfrog stays at 1. Lax-Wendroff stays near 1 for small theta and drops to 0.28 at pi. Upwind drops to 0.2 at pi. Lax-Friedrichs drops fastest for small theta, reaches 0.6 at pi/2 and comes back to 1 at pi." class="viz-plot"><line x1="56" x2="540" y1="58.2" y2="58.2" class="viz-grid"/><text x="46" y="62.2" text-anchor="end" class="viz-tick">1</text><line x1="56" x2="540" y1="128.6" y2="128.6" class="viz-grid"/><text x="46" y="132.6" text-anchor="end" class="viz-tick">0.75</text><line x1="56" x2="540" y1="199.1" y2="199.1" class="viz-grid"/><text x="46" y="203.1" text-anchor="end" class="viz-tick">0.5</text><line x1="56" x2="540" y1="269.5" y2="269.5" class="viz-grid"/><text x="46" y="273.5" text-anchor="end" class="viz-tick">0.25</text><line x1="56" x2="540" y1="340.0" y2="340.0" class="viz-leader"/><text x="46" y="344.0" text-anchor="end" class="viz-tick">0</text><text x="56.0" y="362" text-anchor="middle" class="viz-tick">0</text><text x="177.0" y="362" text-anchor="middle" class="viz-tick">π/4</text><text x="298.0" y="362" text-anchor="middle" class="viz-tick">π/2</text><text x="419.0" y="362" text-anchor="middle" class="viz-tick">3π/4</text><text x="540.0" y="362" text-anchor="middle" class="viz-tick">π</text><text x="298.0" y="388" text-anchor="middle" class="viz-axis-title">Grid frequency θ = ξh</text><text x="56" y="14" class="viz-axis-title">|g(θ)|, measured by one step, ν = 0.6</text><path d="M56.0,58.2 L71.1,58.5 L86.2,59.5 L101.4,61.1 L116.5,63.4 L131.6,66.3 L146.8,69.8 L161.9,74.0 L177.0,78.7 L192.1,84.1 L207.2,90.0 L222.4,96.5 L237.5,103.6 L252.6,111.2 L267.8,119.2 L282.9,127.8 L298.0,136.8 L313.1,146.2 L328.2,156.0 L343.4,166.1 L358.5,176.6 L373.6,187.3 L388.7,198.2 L403.9,209.2 L419.0,220.2 L434.1,231.2 L449.3,242.0 L464.4,252.4 L479.5,262.0 L494.6,270.6 L509.7,277.5 L524.9,282.0 L540.0,283.6" fill="none" stroke="#2a78d6" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><path d="M56.0,58.2 L71.1,59.0 L86.2,61.6 L101.4,65.9 L116.5,71.7 L131.6,79.0 L146.8,87.5 L161.9,97.2 L177.0,107.6 L192.1,118.5 L207.2,129.6 L222.4,140.3 L237.5,150.2 L252.6,158.7 L267.8,165.3 L282.9,169.5 L298.0,170.9 L313.1,169.5 L328.2,165.3 L343.4,158.7 L358.5,150.2 L373.6,140.3 L388.7,129.6 L403.9,118.5 L419.0,107.6 L434.1,97.2 L449.3,87.5 L464.4,79.0 L479.5,71.7 L494.6,65.9 L509.7,61.6 L524.9,59.0 L540.0,58.2" fill="none" stroke="#eb6834" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><path d="M56.0,58.2 L71.1,58.2 L86.2,58.2 L101.4,58.2 L116.5,58.4 L131.6,58.6 L146.8,59.1 L161.9,59.9 L177.0,61.0 L192.1,62.6 L207.2,64.7 L222.4,67.4 L237.5,70.8 L252.6,75.0 L267.8,80.1 L282.9,86.0 L298.0,92.8 L313.1,100.5 L328.2,109.2 L343.4,118.7 L358.5,129.2 L373.6,140.5 L388.7,152.5 L403.9,165.2 L419.0,178.5 L434.1,192.0 L449.3,205.7 L464.4,219.1 L479.5,231.9 L494.6,243.3 L509.7,252.7 L524.9,258.9 L540.0,261.1" fill="none" stroke="#1baf7a" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><path d="M56.0,58.2 L71.1,58.2 L86.2,58.2 L101.4,58.2 L116.5,58.2 L131.6,58.2 L146.8,58.2 L161.9,58.2 L177.0,58.2 L192.1,58.2 L207.2,58.2 L222.4,58.2 L237.5,58.2 L252.6,58.2 L267.8,58.2 L282.9,58.2 L298.0,58.2 L313.1,58.2 L328.2,58.2 L343.4,58.2 L358.5,58.2 L373.6,58.2 L388.7,58.2 L403.9,58.2 L419.0,58.2 L434.1,58.2 L449.3,58.2 L464.4,58.2 L479.5,58.2 L494.6,58.2 L509.7,58.2 L524.9,58.2 L540.0,58.2" fill="none" stroke="#eda100" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><text x="550" y="287.6" class="viz-label">upwind</text><text x="550" y="78.2" class="viz-label">Lax–Friedrichs</text><text x="550" y="265.1" class="viz-label">Lax–Wendroff</text><text x="550" y="62.2" class="viz-label">leapfrog</text><line x1="0" x2="0" y1="30" y2="340" class="viz-crosshair" visibility="hidden"/></svg>
<div class="viz-tooltip" hidden></div>
</div>
<figcaption>|g(θ)| for u_t + u_x = 0 at Courant number ν = 0.6, measured by one step of each scheme on each grid mode of a periodic grid with 64 points. Leapfrog is computed from its two roots; both have size 1. Lax–Friedrichs comes back to 1 at θ = π and meets the leapfrog line there. Hover over the chart, or focus it and use the arrow keys, to read the values. The table shows every second mode.</figcaption>
<details class="viz-table"><summary>Show the numbers</summary><table><thead><tr><th>θ / π</th><th>upwind</th><th>Lax–Friedrichs</th><th>Lax–Wendroff</th><th>leapfrog</th></tr></thead><tbody><tr><td>0.0000</td><td>1.0000</td><td>1.0000</td><td>1.0000</td><td>1.0000</td></tr><tr><td>0.0625</td><td>0.9954</td><td>0.9877</td><td>1.0000</td><td>1.0000</td></tr><tr><td>0.1250</td><td>0.9816</td><td>0.9520</td><td>0.9993</td><td>1.0000</td></tr><tr><td>0.1875</td><td>0.9587</td><td>0.8958</td><td>0.9967</td><td>1.0000</td></tr><tr><td>0.2500</td><td>0.9270</td><td>0.8246</td><td>0.9901</td><td>1.0000</td></tr><tr><td>0.3125</td><td>0.8869</td><td>0.7467</td><td>0.9770</td><td>1.0000</td></tr><tr><td>0.3750</td><td>0.8389</td><td>0.6736</td><td>0.9551</td><td>1.0000</td></tr><tr><td>0.4375</td><td>0.7834</td><td>0.6200</td><td>0.9223</td><td>1.0000</td></tr><tr><td>0.5000</td><td>0.7211</td><td>0.6000</td><td>0.8773</td><td>1.0000</td></tr><tr><td>0.5625</td><td>0.6530</td><td>0.6200</td><td>0.8191</td><td>1.0000</td></tr><tr><td>0.6250</td><td>0.5799</td><td>0.6736</td><td>0.7480</td><td>1.0000</td></tr><tr><td>0.6875</td><td>0.5033</td><td>0.7467</td><td>0.6652</td><td>1.0000</td></tr><tr><td>0.7500</td><td>0.4250</td><td>0.8246</td><td>0.5732</td><td>1.0000</td></tr><tr><td>0.8125</td><td>0.3477</td><td>0.8958</td><td>0.4766</td><td>1.0000</td></tr><tr><td>0.8750</td><td>0.2767</td><td>0.9520</td><td>0.3837</td><td>1.0000</td></tr><tr><td>0.9375</td><td>0.2219</td><td>0.9877</td><td>0.3099</td><td>1.0000</td></tr><tr><td>1.0000</td><td>0.2000</td><td>1.0000</td><td>0.2800</td><td>1.0000</td></tr></tbody></table></details>
<script>(() => {const fig = document.currentScript.closest('figure');const svg = fig.querySelector('svg');const line = svg.querySelector('.viz-crosshair');const tip = fig.querySelector('.viz-tooltip');const names = ["upwind", "Lax\u2013Friedrichs", "Lax\u2013Wendroff", "leapfrog"]; const colors = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"];const rows = [["1.0000", "1.0000", "1.0000", "1.0000"], ["0.9988", "0.9969", "1.0000", "1.0000"], ["0.9954", "0.9877", "1.0000", "1.0000"], ["0.9896", "0.9727", "0.9998", "1.0000"], ["0.9816", "0.9520", "0.9993", "1.0000"], ["0.9712", "0.9262", "0.9984", "1.0000"], ["0.9587", "0.8958", "0.9967", "1.0000"], ["0.9440", "0.8616", "0.9940", "1.0000"], ["0.9270", "0.8246", "0.9901", "1.0000"], ["0.9080", "0.7859", "0.9845", "1.0000"], ["0.8869", "0.7467", "0.9770", "1.0000"], ["0.8639", "0.7087", "0.9673", "1.0000"], ["0.8389", "0.6736", "0.9551", "1.0000"], ["0.8120", "0.6434", "0.9402", "1.0000"], ["0.7834", "0.6200", "0.9223", "1.0000"], ["0.7530", "0.6051", "0.9014", "1.0000"], ["0.7211", "0.6000", "0.8773", "1.0000"], ["0.6877", "0.6051", "0.8498", "1.0000"], ["0.6530", "0.6200", "0.8191", "1.0000"], ["0.6170", "0.6434", "0.7851", "1.0000"], ["0.5799", "0.6736", "0.7480", "1.0000"], ["0.5420", "0.7087", "0.7079", "1.0000"], ["0.5033", "0.7467", "0.6652", "1.0000"], ["0.4642", "0.7859", "0.6201", "1.0000"], ["0.4250", "0.8246", "0.5732", "1.0000"], ["0.3859", "0.8616", "0.5251", "1.0000"], ["0.3477", "0.8958", "0.4766", "1.0000"], ["0.3109", "0.9262", "0.4290", "1.0000"], ["0.2767", "0.9520", "0.3837", "1.0000"], ["0.2463", "0.9727", "0.3430", "1.0000"], ["0.2219", "0.9877", "0.3099", "1.0000"], ["0.2057", "0.9969", "0.2878", "1.0000"], ["0.2000", "1.0000", "0.2800", "1.0000"]]; const heads = ["\u03b8 = 0.0000\u03c0", "\u03b8 = 0.0312\u03c0", "\u03b8 = 0.0625\u03c0", "\u03b8 = 0.0938\u03c0", "\u03b8 = 0.1250\u03c0", "\u03b8 = 0.1562\u03c0", "\u03b8 = 0.1875\u03c0", "\u03b8 = 0.2188\u03c0", "\u03b8 = 0.2500\u03c0", "\u03b8 = 0.2812\u03c0", "\u03b8 = 0.3125\u03c0", "\u03b8 = 0.3437\u03c0", "\u03b8 = 0.3750\u03c0", "\u03b8 = 0.4063\u03c0", "\u03b8 = 0.4375\u03c0", "\u03b8 = 0.4687\u03c0", "\u03b8 = 0.5000\u03c0", "\u03b8 = 0.5312\u03c0", "\u03b8 = 0.5625\u03c0", "\u03b8 = 0.5938\u03c0", "\u03b8 = 0.6250\u03c0", "\u03b8 = 0.6562\u03c0", "\u03b8 = 0.6875\u03c0", "\u03b8 = 0.7188\u03c0", "\u03b8 = 0.7500\u03c0", "\u03b8 = 0.7812\u03c0", "\u03b8 = 0.8125\u03c0", "\u03b8 = 0.8438\u03c0", "\u03b8 = 0.8750\u03c0", "\u03b8 = 0.9062\u03c0", "\u03b8 = 0.9375\u03c0", "\u03b8 = 0.9688\u03c0", "\u03b8 = 1.0000\u03c0"];const xs = [56.0, 71.1, 86.2, 101.4, 116.5, 131.6, 146.8, 161.9, 177.0, 192.1, 207.2, 222.4, 237.5, 252.6, 267.8, 282.9, 298.0, 313.1, 328.2, 343.4, 358.5, 373.6, 388.7, 403.9, 419.0, 434.1, 449.3, 464.4, 479.5, 494.6, 509.7, 524.9, 540.0]; const steps = 32, width = 680;let current = -1;const show = (n) => { n = Math.max(0, Math.min(steps, n)); current = n; const x = xs[n]; line.setAttribute('x1', x); line.setAttribute('x2', x); line.setAttribute('visibility', 'visible'); tip.replaceChildren(); const head = document.createElement('div'); head.className = 'viz-tooltip-head'; head.textContent = heads[n]; tip.append(head); names.forEach((name, i) => {  const row = document.createElement('div'); row.className = 'viz-tooltip-row';  const key = document.createElement('span'); key.className = 'viz-swatch'; key.style.background = colors[i];  const value = document.createElement('strong'); value.textContent = rows[n][i];  const label = document.createElement('span'); label.textContent = name;  row.append(key, value, label); tip.append(row); }); const box = svg.getBoundingClientRect(), frame = fig.getBoundingClientRect(); const xp = box.left - frame.left + x * box.width / width; tip.hidden = false; const room = frame.width - tip.offsetWidth - 8; tip.style.left = Math.max(0, Math.min(room, xp + 12)) + 'px'; tip.style.top = (box.top - frame.top + 8) + 'px';};const hide = () => { line.setAttribute('visibility', 'hidden'); tip.hidden = true; current = -1; };svg.addEventListener('pointermove', (event) => { const box = svg.getBoundingClientRect(); const x = (event.clientX - box.left) * width / box.width; let best = 0; xs.forEach((v, i) => { if (Math.abs(v - x) < Math.abs(xs[best] - x)) best = i; }); show(best);});svg.addEventListener('pointerleave', hide);svg.addEventListener('focus', () => show(current < 0 ? 0 : current));svg.addEventListener('blur', hide);svg.addEventListener('keydown', (event) => { if (event.key === 'ArrowRight') { show(current + 1); event.preventDefault(); } if (event.key === 'ArrowLeft') { show(current - 1); event.preventDefault(); } if (event.key === 'Escape') hide();});})();</script>
</figure>

For a smooth mode with $64$ points per wavelength, the loss per step $1 - |g|$ is $3.08 \times 10^{-3}$ for Lax–Friedrichs, $1.16 \times 10^{-3}$ for upwind, and $2.67 \times 10^{-6}$ for Lax–Wendroff. The first two match $\tfrac12(1 - \nu^2)\theta^2$ and $\tfrac12\nu(1 - \nu)\theta^2$.

The largest $|g|$ over all modes, for other Courant numbers:

| $\nu$ | 0.5 | 1.0 | 1.1 | −0.3 |
|---|---|---|---|---|
| upwind | $1.0000$ | $1.0000$ | $1.2000$ | $1.6000$ |
| Lax–Friedrichs | $1.0000$ | $1.0000$ | $1.1000$ | $1.0000$ |
| Lax–Wendroff | $1.0000$ | $1.0000$ | $1.4200$ | $1.0000$ |

- For $\nu > 1$ all three are unstable. For upwind the worst mode is the sawtooth, with $|1 - 2\nu| = 1.2$.
- For $\nu = -0.3$ (flow to the left), Lax–Friedrichs and Lax–Wendroff are fine, because they only need $|\nu| \le 1$. Upwind with the left-side difference is not: $|g(\pi)| = 1 + 2|\nu| = 1.6$. This is the "wrong side" of Lecture 12.

For leapfrog with $\nu = 0.9$, both roots have $|g| = 1$ to within $1.1 \times 10^{-16}$ for every mode. To see the double root, I ran leapfrog on the mode $\theta = \pi/2$ ($8$ grid points), with the simple start $U^1 = U^0$ (so the start is not exactly one of the two root modes). The largest $|U|$:

| step $n$ | 10 | 100 | 1000 |
|---|---|---|---|
| $\nu = 0.9$ | $2.893$ | $2.508$ | $3.162$ |
| $\nu = 1.0$ | $13.45$ | $140.7$ | $1413.5$ |

With $\nu = 0.9$ the solution stays bounded. (It is larger than $1$ because the bad start excites both roots, but it does not grow.) With $\nu = 1$ it grows like $1.41\,n$: linear growth from the double root.

**Part 2: Lax–Friedrichs is consistent only if $\epsilon \to 0$.** I ran Lax–Friedrichs from $\sin(2\pi x)$ with $a = 1$ to $T = 0.25$, for the three choices of section 1. The exact solution is $\sin(2\pi(x - t))$, so the largest $|U|$ should stay $1$. I also compare with the exact solution of the modified equation $u_t + u_x = \epsilon u_{xx}$, which is $e^{-4\pi^2\epsilon t}\sin(2\pi(x - t))$.

| $k = h/2$ | | | | |
|---|---|---|---|---|
| points | 20 | 40 | 80 | 160 |
| $\epsilon$ | $0.05$ | $0.025$ | $0.0125$ | $0.00625$ |
| largest $\lvert U\rvert$ | $0.689$ | $0.831$ | $0.912$ | $0.955$ |
| error vs. advection | $3.11 \times 10^{-1}$ | $1.69 \times 10^{-1}$ | $8.84 \times 10^{-2}$ | $4.52 \times 10^{-2}$ |

| $k = 4h^2$ | | | | |
|---|---|---|---|---|
| points | 20 | 40 | 80 | 160 |
| $\epsilon$ | $0.125$ | $0.125$ | $0.125$ | $0.125$ |
| largest $\lvert U\rvert$ | $0.300$ | $0.293$ | $0.292$ | $0.291$ |
| error vs. advection | $7.00 \times 10^{-1}$ | $7.07 \times 10^{-1}$ | $7.08 \times 10^{-1}$ | $7.09 \times 10^{-1}$ |
| error vs. $u_t + u_x = 0.125\,u_{xx}$ | $1.78 \times 10^{-2}$ | $4.35 \times 10^{-3}$ | $1.08 \times 10^{-3}$ | $2.70 \times 10^{-4}$ |

| $k = 4h^3$ | | | | |
|---|---|---|---|---|
| points | 10 | 20 | 40 | 80 |
| steps | 62 | 500 | 4000 | 32000 |
| $\epsilon$ | $1.24$ | $2.5$ | $5$ | $10$ |
| largest $\lvert U\rvert$ | $1.96 \times 10^{-6}$ | $1.27 \times 10^{-11}$ | $3.02 \times 10^{-22}$ | $1.24 \times 10^{-43}$ |

- With $k/h$ fixed, the error halves with each refinement: first-order convergence to the right answer.
- With $k = 4h^2$, the error against the advection solution does not go down at all. The scheme converges very nicely (the error drops by $4$ each time) to a different function: the solution of $u_t + u_x = 0.125\,u_{xx}$, whose height is $e^{-\pi^2/8} = 0.291$.
- With $k = 4h^3$, the diffusion grows with refinement and wipes out the wave completely. Smaller steps made the result worse.

**Part 3: conservation form.** First, two small checks of section 8 and section 10 on $50$ cells. With constant $\beta$, one finite volume step with the diffusive flux and one FTCS step differ by at most $1.1 \times 10^{-16}$. With the flux $\mathcal F = aQ_{i-1}$, one finite volume step and one upwind step are identical (difference $0$).

Then a case where the form matters: $q_t = (\beta(x)\,q_x)_x$ with $\beta(x) = 1 + 0.5\sin(2\pi x)$, from the bump $e^{-60(x - 0.4)^2}$, to $T = 0.02$ with $\Delta t = \Delta x^2/6$. I compare two schemes:

- **conservative**: the flux form of section 10, with $\beta$ evaluated at the cell walls;
- **non-conservative**: expand $(\beta q_x)_x = \beta\,q_{xx} + \beta'\,q_x$ and use central differences for $q_{xx}$ and $q_x$, with $\beta$ and $\beta'$ at the cell centers.

Both are consistent and second order in space. The relative change of the total amount $\Delta x\sum_i Q_i$ from $t = 0$ to $T$:

| cells | 20 | 40 | 80 | 160 |
|---|---|---|---|---|
| steps | 48 | 192 | 768 | 3072 |
| conservative | $2.4 \times 10^{-16}$ | $0$ | $-1.2 \times 10^{-16}$ | $1.2 \times 10^{-16}$ |
| non-conservative | $-1.17 \times 10^{-3}$ | $-2.91 \times 10^{-4}$ | $-7.27 \times 10^{-5}$ | $-1.82 \times 10^{-5}$ |
| largest gap between the two solutions | $7.59 \times 10^{-4}$ | $1.93 \times 10^{-4}$ | $4.85 \times 10^{-5}$ | $1.21 \times 10^{-5}$ |

The conservative scheme keeps the total to rounding, at every resolution. The non-conservative one loses mass. The loss goes to zero like $\Delta x^2$, so it is a truncation error, not a bug, but it never vanishes. For a smooth problem like this both schemes converge to the same answer (the gap also drops like $\Delta x^2$). Keeping the total amount exactly is the property that finite volume methods are built for, and it will matter more for problems with shocks later in the course.

The full script and the chart code are in [this site's repository](https://github.com/doubleyands/doubleyands.github.io/tree/main/scripts/posts/numerical-pde/lecture-13).

## Summary

1. Lax–Friedrichs, upwind and Lax–Wendroff all add artificial diffusion. The amount matters: it must be in the stable range, and it must go to zero as $h, k \to 0$, or the scheme solves a different problem.
2. von Neumann analysis with the Courant number $\nu = ak/h$ gives: upwind $g = (1 - \nu) + \nu e^{-i\theta}$, stable for $0 \le \nu \le 1$; Lax–Friedrichs $g = \cos\theta - i\nu\sin\theta$ and Lax–Wendroff, both stable for $|\nu| \le 1$; leapfrog $|g| = 1$ for $|\nu| < 1$, with a double root at $|\nu| = 1$.
3. For these schemes $g = 1 + k\lambda$: von Neumann analysis and the eigenvalue method give the same answer.
4. Finite volume methods start from the **integral form** of a conservation law, which needs no derivatives. The unknowns are **cell averages**.
5. The exact cell update needs the time-averaged flux at the cell walls. Replacing it by a **numerical flux** $\mathcal F(Q_{i-1}, Q_i)$ gives a scheme in **conservation form**, and the telescoping sum makes it conserve the total amount exactly.
6. The choice of $\mathcal F$ decides everything. The upwind flux gives the upwind scheme, and the diffusive flux gives FTCS for the heat equation.
