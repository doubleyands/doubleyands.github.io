---
title: "Numerical PDEs, Lecture 14: Hyperbolic Conservation Laws, the CFL Condition, and Numerical Fluxes"
description: "My notes on the fourteenth lecture: conservation laws in several dimensions, what makes a system hyperbolic, characteristics from the total derivative, domains of dependence, the CFL condition, and four numerical fluxes: central, Lax–Friedrichs, Richtmyer and upwind."
pubDate: 2026-11-19
project: numerical-pde
tags: [numerical PDEs, conservation laws, hyperbolic, CFL condition, finite volume, numerical flux, notes]
---

This is my summary of the fourteenth lecture in my numerical PDE course. I write the ideas in my own words and check them with a small experiment.

[Lecture 13](/blog/2026-11-17-numerical-pde-lecture-13/) started the **finite volume method**: cell averages, the numerical flux, and the conservation form. This lecture does two things. First, it looks at conservation laws more closely: where they come from, what "hyperbolic" means, and how information moves along characteristics. This leads to the **CFL condition**, a simple test that every explicit scheme for a hyperbolic problem must pass. Second, it goes back to the finite volume scheme and tries four numerical fluxes.

Notation as before: $\Delta x$ is the cell size, $\Delta t$ is the time step, and $Q^n_i$ is the average over cell $i$ at time $t_n$. For finite differences I also write $h$ and $k$.

## 1. Conservation laws in several dimensions

Take any fixed region $V$ in space, with density $\rho$ and velocity $\mathbf u$. The mass inside is

$$
m_V(t) = \int_V \rho\,dV.
$$

Mass changes only by flow through the boundary $\partial V$. Only the **normal part** of the velocity matters: if the flow runs along the boundary ($\mathbf n\cdot\mathbf u = 0$), nothing goes in or out. With $\mathbf n$ the outward normal,

$$
\frac{d}{dt}\int_V \rho\,dV = -\int_{\partial V} \rho\,\mathbf n\cdot\mathbf u\,dS.
$$

The surface integral measures what goes **out**, so it gets a minus sign. The divergence theorem ([Lecture 1](/blog/2026-10-06-numerical-pde-lecture-1/)) turns it into a volume integral. The result holds for every region $V$, so the integrand must be zero (the same argument as in Lecture 13):

$$
\frac{\partial\rho}{\partial t} + \nabla\cdot(\rho\,\mathbf u) = 0.
$$

This is the general pattern. For a conserved quantity $u$ with flux function $\mathbf f(u)$,

$$
\frac{\partial u}{\partial t} + \nabla\cdot\mathbf f(u) = 0, \qquad \frac{d}{dt}\int_V u\,dV + \int_{\partial V}\mathbf n\cdot\mathbf f(u)\,dS = 0.
$$

If the density of a fluid is constant, mass conservation becomes $\nabla\cdot\mathbf u = 0$, the **incompressibility condition**. Conservation of mass, momentum and energy together give the equations of fluid dynamics.

In 1D the law is $u_t + f(u)_x = 0$. Two flux functions come up again and again:

- $f(u) = cu$ gives the advection equation $u_t + cu_x = 0$, our model problem.
- $f(u) = u^2/2$ gives **Burgers' equation** $u_t + \left(\tfrac12 u^2\right)_x = 0$. It is nonlinear but simple, and a lot is known about it. So a new scheme is usually tested first on a linear problem, then on Burgers.

If the flux contains a derivative, for example $f = -u_x$, the law becomes the diffusion equation. That is a conservation law too, but it is not **hyperbolic**. Hyperbolic conservation laws describe **transport**: something moves from one place to another.

## 2. Systems, and what "hyperbolic" means

With two unknowns $u, v$ we have a system

$$
u_t + f(u, v)_x = 0, \qquad v_t + g(u, v)_x = 0.
$$

Put them in a vector $\mathbf u = (u, v)$ and write $\mathbf u_t + \mathbf f(\mathbf u)_x = \mathbf 0$. It looks just like the scalar law, in the same way that a system of ODEs looks like one ODE. The chain rule gives $\mathbf f(\mathbf u)_x = D\mathbf f(\mathbf u)\,\mathbf u_x$, where

$$
D\mathbf f(\mathbf u) = \begin{pmatrix} \partial f/\partial u & \partial f/\partial v \\ \partial g/\partial u & \partial g/\partial v \end{pmatrix}
$$

is the **Jacobian**.

**Definition.** The system is **hyperbolic** if the Jacobian is diagonalizable with **real eigenvalues**.

Why real? Look at the linear case $\mathbf u_t + A\mathbf u_x = 0$ with $A = R\Lambda R^{-1}$. Set $\mathbf w = R^{-1}\mathbf u$. Then

$$
\mathbf w_t + \Lambda\,\mathbf w_x = \mathbf 0, \qquad\text{that is}\qquad (w_p)_t + \lambda_p\,(w_p)_x = 0 \ \text{ for each } p.
$$

The system splits into **independent advection equations**. Each one is a wave that moves with speed $\lambda_p$. A speed must be a real number.

Two examples of my own:

- $A = \begin{pmatrix} 0 & 1 \\ 4 & 0 \end{pmatrix}$ has eigenvalues $\pm 2$ with eigenvectors $(1, 2)$ and $(1, -2)$. It is hyperbolic: any initial state splits into one wave moving right with speed $2$ and one moving left with speed $2$. (This is the wave equation $u_{tt} = 4u_{xx}$ written as a first-order system.)
- $B = \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}$ has eigenvalues $\pm i$. It is **not** hyperbolic. Eliminating $v$ gives $u_{tt} + u_{xx} = 0$, which is Laplace's equation, now posed as an initial value problem. A mode $e^{i\xi x}$ grows like $e^{|\xi| t}$, so short waves grow without limit. Part 4 of the experiment shows what this does to a scheme.

## 3. Conservation form and primitive form

For smooth solutions the chain rule gives two ways to write the same 1D law:

$$
u_t + f(u)_x = 0 \quad \text{(conservation form)}, \qquad u_t + f'(u)\,u_x = 0 \quad \text{(primitive form)}.
$$

The primitive form looks like the linear advection equation with speed $f'(u)$. So **$f'(u)$ is the wave speed**. For $f = cu$ it is $c$. For Burgers it is $u$ itself: large values move faster.

The conservation form keeps the flux as it is. It is the continuous partner of the discrete conservation form of Lecture 13,

$$
\frac{Q^{n+1}_i - Q^n_i}{\Delta t} + \frac{F^n_{i+1/2} - F^n_{i-1/2}}{\Delta x} = 0.
$$

For smooth solutions the two forms are the same. They stop being the same when the solution has jumps. That comes later in the course.

## 4. Characteristics, derived

In [Lecture 12](/blog/2026-11-12-numerical-pde-lecture-12/) the solution $u(x, t) = \eta(x - ct)$ of $u_t + cu_x = 0$ was given first and then checked. Here is a way to **find** it.

Follow a point that moves with speed $c$:

$$
X(t) = x_0 + ct, \qquad \dot X(t) = c, \qquad X(0) = x_0.
$$

This is why $c$ is called a speed: it is the time derivative of a position. Now take the total derivative of the solution along this path:

$$
\frac{d}{dt}\,u\big(X(t), t\big) = u_t + u_x\,\dot X = u_t + c\,u_x = 0.
$$

The derivative is zero, so $u$ is **constant along the path**. The value at $(x, t)$ is the value at the start of the line through it, $x_0 = x - ct$. That is exactly $u(x, t) = \eta(x - ct)$. Such a path is a **characteristic**. A good picture: you sit on a small boat and drift with the current. The time goes up, but your position moves with the flow, and the value you see does not change.

The same idea works for a nonlinear law. Use the wave speed from the primitive form: if $\dot X = f'\big(u(X(t), t)\big)$, then again $\frac{d}{dt}u(X(t), t) = u_t + f'(u)\,u_x = 0$. For Burgers this means $\dot X = u$. The characteristics are still straight lines, but now their slopes depend on the data. Burgers' equation will come back later in the course.

**Reading the slopes.** In the $x$–$t$ plane, with $t$ going up, a characteristic $x = x_0 + ct$ has slope $1/c$, not $c$. So a **steeper** line means a **slower** wave. A line that leans left means $c < 0$. This is easy to get backward, so it is worth checking each time.

## 5. Domain of dependence and domain of influence

One characteristic can be read in two directions.

- **Domain of influence.** Start from a point $(x_0, 0)$ on the initial line. Its value is carried forward along the characteristic. The points on that line are the domain of influence of $x_0$.
- **Domain of dependence.** Start from a point $(X, T)$. Which initial values decide $u(X, T)$? For advection, only one: the value at the **foot of the characteristic**, $X - cT$. So the domain of dependence is the single point $\{X - cT\}$.

These are two views of the same thing. In 1D they are simple. In 2D and 3D, and for systems with several speeds, they become regions.

## 6. The CFL condition

Take $c > 0$ and the upwind scheme: Forward Euler in time, and the backward difference in space, because the flow comes from the left:

$$
\frac{U^{n+1}_j - U^n_j}{\Delta t} + c\,\frac{U^n_j - U^n_{j-1}}{\Delta x} = 0.
$$

Now follow the stencil backward. $U^2_3$ comes from $U^1_2$ and $U^1_3$. These come from $U^0_1$, $U^0_2$ and $U^0_3$. In general, $U^n_j$ depends on the initial values at $x_{j-n}, \dots, x_j$. These points form a triangle with its top at $(x_j, t_n)$. This triangle is the **numerical domain of dependence**. Its left edge goes back $\Delta x$ in space for each $\Delta t$ in time, so its slope is $\Delta t/\Delta x$.

The true solution at $(x_j, t_n)$ depends on the foot of the characteristic, $x_j - ct_n$. If the characteristic is steeper than the left edge of the triangle, the foot lies inside. If it is flatter, the foot lies outside, and the scheme **never sees** the value that decides the true solution. Refining the grid with a fixed ratio does not help: the triangle and the characteristic keep the same shapes.

> **The CFL condition.** The numerical domain of dependence must contain the physical domain of dependence.

For upwind this is

$$
\frac{\Delta t}{\Delta x} \le \frac1c \quad\Longleftrightarrow\quad \nu = \frac{c\,\Delta t}{\Delta x} \le 1 \quad\Longleftrightarrow\quad \Delta t \le \frac{\Delta x}{c}.
$$

The number $\nu$ is the **Courant number** (or CFL number) from Lecture 13. The condition is named after Courant, Friedrichs and Lewy. They wrote about it in 1928, before there were computers.

Once you see it, the idea is almost obvious: a scheme must use the **right information**. It also costs almost nothing to check. You only need the stencil and the characteristic speeds. Compare this with the other tools we have: the eigenvalues of the method-of-lines matrix (Lecture 12) and von Neumann analysis (Lecture 13).

Three remarks:

- **The CFL condition is necessary, not sufficient.** A 3-point stencil like FTCS contains the foot whenever $|\nu| \le 1$, but FTCS is unstable for every fixed $\nu$ (Lecture 12). Section 7 meets this again.
- **The time step is the price.** With $\Delta t = \nu\,\Delta x/c$, a larger allowed $\nu$ means fewer steps to reach a given time. A small $\nu$ is safe but expensive. Finding schemes that allow a large CFL number, especially for nonlinear problems, is still a research topic.
- **Ordinary stiff ODEs are a different story.** For $u' = \lambda u$ with $\lambda = -2100$, Forward Euler needs $\Delta t \le 2/2100$ ([Lecture 8](/blog/2026-10-29-numerical-pde-lecture-8/)). This is also a limit on the step size, but it comes from **absolute stability** (the point $k\lambda$ must lie in the stability region), not from a domain of dependence.

## 7. Back to finite volumes: four numerical fluxes

Recall the finite volume scheme from Lecture 13:

$$
Q^{n+1}_i = Q^n_i - \frac{\Delta t}{\Delta x}\Big(F^n_{i+1/2} - F^n_{i-1/2}\Big), \qquad F^n_{i-1/2} = \mathcal F\big(Q^n_{i-1}, Q^n_i\big).
$$

The wall $x_{i-1/2}$ is shared by cells $C_{i-1}$ and $C_i$. Everything depends on the choice of the numerical flux $\mathcal F$.

### Central flux: unstable even under CFL

Without extra information, the natural choice is the average:

$$
\mathcal F(Q_{i-1}, Q_i) = \tfrac12\big[f(Q_{i-1}) + f(Q_i)\big] \quad\Longrightarrow\quad Q^{n+1}_i = Q^n_i - \frac{\Delta t}{2\Delta x}\big[f(Q^n_{i+1}) - f(Q^n_{i-1})\big].
$$

For $f = cq$ this is exactly **FTCS**. Its method-of-lines eigenvalues are purely imaginary, so Forward Euler can never put them inside its stability disk (Lecture 12). The average treats both sides "fairly", but in a hyperbolic problem information comes from one side. This is the standard example of a scheme that passes the CFL test and is still unstable.

### Lax–Friedrichs flux: add numerical diffusion

$$
\mathcal F(Q_{i-1}, Q_i) = \tfrac12\big[f(Q_{i-1}) + f(Q_i)\big] - \frac{\Delta x}{2\Delta t}\big(Q_i - Q_{i-1}\big).
$$

Put this into the scheme. The extra term gives $+\tfrac12\big(Q^n_{i+1} - 2Q^n_i + Q^n_{i-1}\big)$, and the update becomes

$$
Q^{n+1}_i = \tfrac12\big(Q^n_{i-1} + Q^n_{i+1}\big) - \frac{\Delta t}{2\Delta x}\big[f(Q^n_{i+1}) - f(Q^n_{i-1})\big],
$$

the Lax–Friedrichs scheme of Lecture 12. The extra term has the same shape as the diffusive flux $-\beta q_x$ of Lecture 13. So in effect the scheme solves

$$
q_t + f(q)_x = \beta\,q_{xx}, \qquad \beta = \frac{\Delta x^2}{2\Delta t}.
$$

This is the artificial diffusion $\epsilon = h^2/(2k)$ of Lecture 12, now in finite volume language. Numerical diffusion is usually unwelcome. Here it is what makes the scheme stable. But Lax–Friedrichs adds **too much**: a sharp profile gets smeared, and the scheme is only first order.

Is it still consistent? Write $\beta = \frac{\Delta x}{2}\cdot\frac{\Delta x}{\Delta t}$. If $\Delta t/\Delta x$ is fixed, one factor $\Delta x$ is left, so $\beta \to 0$ as the grid is refined. Consistency comes first: without it, the scheme solves a different problem. Consistency plus stability gives convergence, by the Lax equivalence theorem ([Lecture 10](/blog/2026-11-05-numerical-pde-lecture-10/)).

### Richtmyer: a two-step flux

First predict the value at the wall half a step ahead, then evaluate the flux there:

$$
Q^{n+1/2}_{i-1/2} = \tfrac12\big(Q^n_{i-1} + Q^n_i\big) - \frac{\Delta t}{2\Delta x}\big[f(Q^n_i) - f(Q^n_{i-1})\big], \qquad F^n_{i-1/2} = f\big(Q^{n+1/2}_{i-1/2}\big).
$$

The predicted value uses both neighbor cells. For a linear flux $f = aq$ with $\nu = a\Delta t/\Delta x$, put it in:

$$
Q^{n+1}_i = Q^n_i - \frac{\nu}{2}\big(Q^n_{i+1} - Q^n_{i-1}\big) + \frac{\nu^2}{2}\big(Q^n_{i+1} - 2Q^n_i + Q^n_{i-1}\big).
$$

This is **Lax–Wendroff** (Lecture 12). So the Richtmyer method is a way to extend Lax–Wendroff to nonlinear fluxes. It is second order, at the cost of one extra stage.

### Upwind flux: a convex combination

In a hyperbolic problem, information travels along characteristics. A one-sided flux can be the better choice. For $f = aq$:

$$
(a > 0)\qquad F^n_{i-1/2} = a\,Q^n_{i-1} \quad\Longrightarrow\quad Q^{n+1}_i = \nu\,Q^n_{i-1} + (1 - \nu)\,Q^n_i.
$$

The new value is a **weighted average** of two old values. If $0 \le \nu \le 1$, both weights are between $0$ and $1$ and they add up to $1$. That is a **convex combination**, and it gives stability in one line. In fact it gives more: $Q^{n+1}_i$ always lies **between** $Q^n_{i-1}$ and $Q^n_i$. So the maximum can never grow and the minimum can never drop. (This is a **maximum principle**. It is a stronger statement than von Neumann stability.)

If $\nu > 1$, the weight $1 - \nu$ is negative. For example, $\nu = 1.05$ gives the weights $1.05$ and $-0.05$. A negative weight pushes values below the old ones and creates wiggles, and repeated steps make them grow. Part 3 of the experiment shows this.

For $a < 0$ the flow comes from the right. Then $F^n_{i-1/2} = a\,Q^n_i$, the update is $Q^{n+1}_i = Q^n_i - \nu\big(Q^n_{i+1} - Q^n_i\big)$, and the condition is $-1 \le \nu \le 0$.

The next lecture builds a scheme from the upwind idea for general problems.

## 8. Experiment

The script uses the Python standard library only (`math`). All grids are periodic on $[0, 1]$, with speed $1$ unless I say otherwise. The initial data is a cosine bump of height $1$ and half-width $0.1$, centered at $x = 0.2$.

**Part 1: the CFL condition with upwind.** I ran upwind on $100$ cells with three Courant numbers. The largest $|U|$ (the exact value stays at $0.9938$, the largest cell value of the bump):

| time $t$ | 0.5 | 1 | 2 | 3 |
|---|---|---|---|---|
| $\nu = 1$ | $0.9938$ | $0.9938$ | $0.9938$ | $0.9938$ |
| $\nu = 0.5$ | $0.6373$ | $0.4992$ | $0.3742$ | $0.3119$ |
| $\nu = 1.05$ | $1.065$ | $7.24$ | $3.41 \times 10^{4}$ | $2.57 \times 10^{8}$ |

(For $\nu = 1.05$ the step does not divide the times exactly, so these are the nearest steps: $t = 0.504$, $0.9975$, $1.995$, $3.003$.)

- With $\nu = 1$, upwind becomes $U^{n+1}_j = U^n_{j-1}$: it moves the data exactly one cell per step. The error stays below $5.4 \times 10^{-15}$, which is rounding.
- With $\nu = 0.5$, the scheme is stable but the bump gets lower and wider. The numerical diffusion of upwind, $\tfrac{c\Delta x}{2}(1 - \nu)$, is larger for smaller $\nu$. So a small Courant number costs more steps **and** smears more.
- With $\nu = 1.05$, the solution is still almost fine at $t = 0.5$. Then sawtooth wiggles appear and grow by a factor of several thousand per unit of time.

<figure class="viz">
<div class="viz-legend"><span class="viz-key"><span class="viz-swatch" style="background:#2a78d6"></span>ν = 1</span><span class="viz-key"><span class="viz-swatch" style="background:#eb6834"></span>ν = 0.5</span><span class="viz-key"><span class="viz-swatch" style="background:#1baf7a"></span>ν = 1.05</span><span class="viz-key"><span class="viz-swatch" style="background:#8e8e93"></span>exact</span></div>
<div class="viz-frame">
<svg viewBox="0 0 680 400" role="img" tabindex="0" aria-label="Upwind solution at time 0.7, for x from 0.5 to 1. With Courant number 1 the bump is exact. With 0.5 it is lower and wider, with height 0.57. With 1.05 it has sawtooth wiggles from minus 0.46 up to 1.24, near the front of the bump and behind it." class="viz-plot"><line x1="56" x2="560" y1="53.2" y2="53.2" class="viz-grid"/><text x="46" y="57.2" text-anchor="end" class="viz-tick">1.25</text><line x1="56" x2="560" y1="92.0" y2="92.0" class="viz-grid"/><text x="46" y="96.0" text-anchor="end" class="viz-tick">1</text><line x1="56" x2="560" y1="130.8" y2="130.8" class="viz-grid"/><text x="46" y="134.8" text-anchor="end" class="viz-tick">0.75</text><line x1="56" x2="560" y1="169.5" y2="169.5" class="viz-grid"/><text x="46" y="173.5" text-anchor="end" class="viz-tick">0.5</text><line x1="56" x2="560" y1="208.2" y2="208.2" class="viz-grid"/><text x="46" y="212.2" text-anchor="end" class="viz-tick">0.25</text><line x1="56" x2="560" y1="247.0" y2="247.0" class="viz-leader"/><text x="46" y="251.0" text-anchor="end" class="viz-tick">0</text><line x1="56" x2="560" y1="285.8" y2="285.8" class="viz-grid"/><text x="46" y="289.8" text-anchor="end" class="viz-tick">-0.25</text><line x1="56" x2="560" y1="324.5" y2="324.5" class="viz-grid"/><text x="46" y="328.5" text-anchor="end" class="viz-tick">-0.5</text><text x="56.0" y="362" text-anchor="middle" class="viz-tick">0.5</text><text x="156.8" y="362" text-anchor="middle" class="viz-tick">0.6</text><text x="257.6" y="362" text-anchor="middle" class="viz-tick">0.7</text><text x="358.4" y="362" text-anchor="middle" class="viz-tick">0.8</text><text x="459.2" y="362" text-anchor="middle" class="viz-tick">0.9</text><text x="560.0" y="362" text-anchor="middle" class="viz-tick">1</text><text x="308.0" y="388" text-anchor="middle" class="viz-axis-title">x</text><text x="56" y="14" class="viz-axis-title">U at t = 0.7, upwind, 100 cells</text><path d="M61.0,247.0 L71.1,247.0 L81.2,247.0 L91.3,247.0 L101.4,247.0 L111.4,247.0 L121.5,247.0 L131.6,247.0 L141.7,247.0 L151.8,247.0 L161.8,247.0 L171.9,247.0 L182.0,247.0 L192.1,247.0 L202.2,247.0 L212.2,247.0 L222.3,247.0 L232.4,247.1 L242.5,246.6 L252.6,248.7 L262.6,241.2 L272.7,263.2 L282.8,210.3 L292.9,312.6 L303.0,161.0 L313.0,317.7 L323.1,228.8 L333.2,228.4 L343.3,247.4 L353.4,258.7 L363.4,261.3 L373.5,254.9 L383.6,240.2 L393.7,218.5 L403.8,192.0 L413.8,163.3 L423.9,135.3 L434.0,110.5 L444.1,92.0 L454.2,78.5 L464.2,83.5 L474.3,67.9 L484.4,135.6 L494.5,54.9 L504.6,233.0 L514.6,105.0 L524.7,221.9 L534.8,247.0 L544.9,247.0 L555.0,247.0" fill="none" stroke="#1baf7a" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><path d="M61.0,247.0 L71.1,247.0 L81.2,247.0 L91.3,247.0 L101.4,247.0 L111.4,247.0 L121.5,247.0 L131.6,247.0 L141.7,247.0 L151.8,247.0 L161.8,247.0 L171.9,247.0 L182.0,247.0 L192.1,247.0 L202.2,246.9 L212.2,246.9 L222.3,246.8 L232.4,246.6 L242.5,246.3 L252.6,245.9 L262.6,245.4 L272.7,244.5 L282.8,243.4 L292.9,241.8 L303.0,239.7 L313.0,236.9 L323.1,233.5 L333.2,229.2 L343.3,224.2 L353.4,218.3 L363.4,211.7 L373.5,204.5 L383.6,197.0 L393.7,189.3 L403.8,181.8 L413.8,174.8 L423.9,168.8 L434.0,163.9 L444.1,160.4 L454.2,158.7 L464.2,158.7 L474.3,160.4 L484.4,163.9 L494.5,168.8 L504.6,174.8 L514.6,181.8 L524.7,189.3 L534.8,197.0 L544.9,204.5 L555.0,211.7" fill="none" stroke="#eb6834" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><path d="M61.0,247.0 L71.1,247.0 L81.2,247.0 L91.3,247.0 L101.4,247.0 L111.4,247.0 L121.5,247.0 L131.6,247.0 L141.7,247.0 L151.8,247.0 L161.8,247.0 L171.9,247.0 L182.0,247.0 L192.1,247.0 L202.2,247.0 L212.2,247.0 L222.3,247.0 L232.4,247.0 L242.5,247.0 L252.6,247.0 L262.6,247.0 L272.7,247.0 L282.8,247.0 L292.9,247.0 L303.0,247.0 L313.0,247.0 L323.1,247.0 L333.2,247.0 L343.3,247.0 L353.4,247.0 L363.4,246.0 L373.5,238.6 L383.6,224.3 L393.7,204.7 L403.8,181.6 L413.8,157.4 L423.9,134.3 L434.0,114.7 L444.1,100.4 L454.2,93.0 L464.2,93.0 L474.3,100.4 L484.4,114.7 L494.5,134.3 L504.6,157.4 L514.6,181.6 L524.7,204.7 L534.8,224.3 L544.9,238.6 L555.0,246.0" fill="none" stroke="#2a78d6" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><path d="M61.0,247.0 L71.1,247.0 L81.2,247.0 L91.3,247.0 L101.4,247.0 L111.4,247.0 L121.5,247.0 L131.6,247.0 L141.7,247.0 L151.8,247.0 L161.8,247.0 L171.9,247.0 L182.0,247.0 L192.1,247.0 L202.2,247.0 L212.2,247.0 L222.3,247.0 L232.4,247.0 L242.5,247.0 L252.6,247.0 L262.6,247.0 L272.7,247.0 L282.8,247.0 L292.9,247.0 L303.0,247.0 L313.0,247.0 L323.1,247.0 L333.2,247.0 L343.3,247.0 L353.4,247.0 L363.4,246.0 L373.5,238.6 L383.6,224.3 L393.7,204.7 L403.8,181.6 L413.8,157.4 L423.9,134.3 L434.0,114.7 L444.1,100.4 L454.2,93.0 L464.2,93.0 L474.3,100.4 L484.4,114.7 L494.5,134.3 L504.6,157.4 L514.6,181.6 L524.7,204.7 L534.8,224.3 L544.9,238.6 L555.0,246.0" fill="none" stroke="#8e8e93" stroke-width="1.5" stroke-dasharray="4 4"/><text x="338.2" y="96.6" class="viz-label">ν = 1</text><text x="570.1" y="214.8" class="viz-label">ν = 0.5</text><text x="156.8" y="174.1" class="viz-label">ν = 1.05</text><line x1="0" x2="0" y1="30" y2="340" class="viz-crosshair" visibility="hidden"/></svg>
<div class="viz-tooltip" hidden></div>
</div>
<figcaption>Upwind for u_t + u_x = 0 on a periodic grid with 100 cells, from a cosine bump centered at x = 0.2. The curves show the right half of the domain at t = 0.7 (for ν = 1.05 the nearest step is t = 0.7035). With ν = 1 the scheme is exact and hides the dashed exact line. Hover over the chart, or focus it and use the arrow keys, to read the values. The table shows every second cell.</figcaption>
<details class="viz-table"><summary>Show the numbers</summary><table><thead><tr><th>x</th><th>ν = 1</th><th>ν = 0.5</th><th>ν = 1.05</th><th>exact</th></tr></thead><tbody><tr><td>0.505</td><td>0.0000</td><td>0.0000</td><td>-0.0000</td><td>0.0000</td></tr><tr><td>0.525</td><td>0.0000</td><td>0.0000</td><td>-0.0000</td><td>0.0000</td></tr><tr><td>0.545</td><td>0.0000</td><td>0.0000</td><td>-0.0000</td><td>0.0000</td></tr><tr><td>0.565</td><td>0.0000</td><td>0.0000</td><td>-0.0000</td><td>0.0000</td></tr><tr><td>0.585</td><td>0.0000</td><td>0.0000</td><td>-0.0000</td><td>0.0000</td></tr><tr><td>0.605</td><td>0.0000</td><td>0.0000</td><td>-0.0000</td><td>0.0000</td></tr><tr><td>0.625</td><td>0.0000</td><td>0.0002</td><td>-0.0000</td><td>0.0000</td></tr><tr><td>0.645</td><td>0.0000</td><td>0.0005</td><td>0.0000</td><td>0.0000</td></tr><tr><td>0.665</td><td>0.0000</td><td>0.0016</td><td>0.0001</td><td>0.0000</td></tr><tr><td>0.685</td><td>0.0000</td><td>0.0043</td><td>0.0029</td><td>0.0000</td></tr><tr><td>0.705</td><td>0.0000</td><td>0.0105</td><td>0.0375</td><td>0.0000</td></tr><tr><td>0.725</td><td>0.0000</td><td>0.0234</td><td>0.2369</td><td>0.0000</td></tr><tr><td>0.745</td><td>0.0000</td><td>0.0474</td><td>0.5549</td><td>0.0000</td></tr><tr><td>0.765</td><td>0.0000</td><td>0.0874</td><td>0.1174</td><td>0.0000</td></tr><tr><td>0.785</td><td>0.0000</td><td>0.1474</td><td>-0.0026</td><td>0.0000</td></tr><tr><td>0.805</td><td>0.0062</td><td>0.2276</td><td>-0.0923</td><td>0.0062</td></tr><tr><td>0.825</td><td>0.1464</td><td>0.3227</td><td>0.0441</td><td>0.1464</td></tr><tr><td>0.845</td><td>0.4218</td><td>0.4206</td><td>0.3547</td><td>0.4218</td></tr><tr><td>0.865</td><td>0.7270</td><td>0.5048</td><td>0.7207</td><td>0.7270</td></tr><tr><td>0.885</td><td>0.9455</td><td>0.5585</td><td>0.9997</td><td>0.9455</td></tr><tr><td>0.905</td><td>0.9938</td><td>0.5698</td><td>1.0548</td><td>0.9938</td></tr><tr><td>0.925</td><td>0.8536</td><td>0.5363</td><td>0.7190</td><td>0.8536</td></tr><tr><td>0.945</td><td>0.5782</td><td>0.4655</td><td>0.0904</td><td>0.5782</td></tr><tr><td>0.965</td><td>0.2730</td><td>0.3722</td><td>0.1618</td><td>0.2730</td></tr><tr><td>0.985</td><td>0.0545</td><td>0.2739</td><td>0.0000</td><td>0.0545</td></tr></tbody></table></details>
<script>(() => {const fig = document.currentScript.closest('figure');const svg = fig.querySelector('svg');const line = svg.querySelector('.viz-crosshair');const tip = fig.querySelector('.viz-tooltip');const names = ["\u03bd = 1", "\u03bd = 0.5", "\u03bd = 1.05", "exact"]; const colors = ["#2a78d6", "#eb6834", "#1baf7a", "#8e8e93"];const rows = [["0.0000", "0.0000", "-0.0000", "0.0000"], ["0.0000", "0.0000", "0.0000", "0.0000"], ["0.0000", "0.0000", "-0.0000", "0.0000"], ["0.0000", "0.0000", "0.0000", "0.0000"], ["0.0000", "0.0000", "-0.0000", "0.0000"], ["0.0000", "0.0000", "0.0000", "0.0000"], ["0.0000", "0.0000", "-0.0000", "0.0000"], ["0.0000", "0.0000", "0.0000", "0.0000"], ["0.0000", "0.0000", "-0.0000", "0.0000"], ["0.0000", "0.0000", "0.0000", "0.0000"], ["0.0000", "0.0000", "-0.0000", "0.0000"], ["0.0000", "0.0001", "0.0000", "0.0000"], ["0.0000", "0.0002", "-0.0000", "0.0000"], ["0.0000", "0.0003", "-0.0000", "0.0000"], ["0.0000", "0.0005", "0.0000", "0.0000"], ["0.0000", "0.0009", "-0.0000", "0.0000"], ["0.0000", "0.0016", "0.0001", "0.0000"], ["0.0000", "0.0026", "-0.0006", "0.0000"], ["0.0000", "0.0043", "0.0029", "0.0000"], ["0.0000", "0.0068", "-0.0113", "0.0000"], ["0.0000", "0.0105", "0.0375", "0.0000"], ["0.0000", "0.0159", "-0.1043", "0.0000"], ["0.0000", "0.0234", "0.2369", "0.0000"], ["0.0000", "0.0337", "-0.4233", "0.0000"], ["0.0000", "0.0474", "0.5549", "0.0000"], ["0.0000", "0.0651", "-0.4559", "0.0000"], ["0.0000", "0.0874", "0.1174", "0.0000"], ["0.0000", "0.1148", "0.1197", "0.0000"], ["0.0000", "0.1474", "-0.0026", "0.0000"], ["0.0000", "0.1852", "-0.0756", "0.0000"], ["0.0062", "0.2276", "-0.0923", "0.0062"], ["0.0545", "0.2739", "-0.0511", "0.0545"], ["0.1464", "0.3227", "0.0441", "0.1464"], ["0.2730", "0.3722", "0.1840", "0.2730"], ["0.4218", "0.4206", "0.3547", "0.4218"], ["0.5782", "0.4655", "0.5397", "0.5782"], ["0.7270", "0.5048", "0.7207", "0.7270"], ["0.8536", "0.5363", "0.8809", "0.8536"], ["0.9455", "0.5585", "0.9997", "0.9455"], ["0.9938", "0.5698", "1.0869", "0.9938"], ["0.9938", "0.5698", "1.0548", "0.9938"], ["0.9455", "0.5585", "1.1553", "0.9455"], ["0.8536", "0.5363", "0.7190", "0.8536"], ["0.7270", "0.5048", "1.2393", "0.7270"], ["0.5782", "0.4655", "0.0904", "0.5782"], ["0.4218", "0.4206", "0.9162", "0.4218"], ["0.2730", "0.3722", "0.1618", "0.2730"], ["0.1464", "0.3227", "0.0000", "0.1464"], ["0.0545", "0.2739", "0.0000", "0.0545"], ["0.0062", "0.2276", "0.0000", "0.0062"]]; const heads = ["x = 0.505", "x = 0.515", "x = 0.525", "x = 0.535", "x = 0.545", "x = 0.555", "x = 0.565", "x = 0.575", "x = 0.585", "x = 0.595", "x = 0.605", "x = 0.615", "x = 0.625", "x = 0.635", "x = 0.645", "x = 0.655", "x = 0.665", "x = 0.675", "x = 0.685", "x = 0.695", "x = 0.705", "x = 0.715", "x = 0.725", "x = 0.735", "x = 0.745", "x = 0.755", "x = 0.765", "x = 0.775", "x = 0.785", "x = 0.795", "x = 0.805", "x = 0.815", "x = 0.825", "x = 0.835", "x = 0.845", "x = 0.855", "x = 0.865", "x = 0.875", "x = 0.885", "x = 0.895", "x = 0.905", "x = 0.915", "x = 0.925", "x = 0.935", "x = 0.945", "x = 0.955", "x = 0.965", "x = 0.975", "x = 0.985", "x = 0.995"];const xs = [61.0, 71.1, 81.2, 91.3, 101.4, 111.4, 121.5, 131.6, 141.7, 151.8, 161.8, 171.9, 182.0, 192.1, 202.2, 212.2, 222.3, 232.4, 242.5, 252.6, 262.6, 272.7, 282.8, 292.9, 303.0, 313.0, 323.1, 333.2, 343.3, 353.4, 363.4, 373.5, 383.6, 393.7, 403.8, 413.8, 423.9, 434.0, 444.1, 454.2, 464.2, 474.3, 484.4, 494.5, 504.6, 514.6, 524.7, 534.8, 544.9, 555.0]; const steps = 49, width = 680;let current = -1;const show = (n) => { n = Math.max(0, Math.min(steps, n)); current = n; const x = xs[n]; line.setAttribute('x1', x); line.setAttribute('x2', x); line.setAttribute('visibility', 'visible'); tip.replaceChildren(); const head = document.createElement('div'); head.className = 'viz-tooltip-head'; head.textContent = heads[n]; tip.append(head); names.forEach((name, i) => {  const row = document.createElement('div'); row.className = 'viz-tooltip-row';  const key = document.createElement('span'); key.className = 'viz-swatch'; key.style.background = colors[i];  const value = document.createElement('strong'); value.textContent = rows[n][i];  const label = document.createElement('span'); label.textContent = name;  row.append(key, value, label); tip.append(row); }); const box = svg.getBoundingClientRect(), frame = fig.getBoundingClientRect(); const xp = box.left - frame.left + x * box.width / width; tip.hidden = false; const room = frame.width - tip.offsetWidth - 8; tip.style.left = Math.max(0, Math.min(room, xp + 12)) + 'px'; tip.style.top = (box.top - frame.top + 8) + 'px';};const hide = () => { line.setAttribute('visibility', 'hidden'); tip.hidden = true; current = -1; };svg.addEventListener('pointermove', (event) => { const box = svg.getBoundingClientRect(); const x = (event.clientX - box.left) * width / box.width; let best = 0; xs.forEach((v, i) => { if (Math.abs(v - x) < Math.abs(xs[best] - x)) best = i; }); show(best);});svg.addEventListener('pointerleave', hide);svg.addEventListener('focus', () => show(current < 0 ? 0 : current));svg.addEventListener('blur', hide);svg.addEventListener('keydown', (event) => { if (event.key === 'ArrowRight') { show(current + 1); event.preventDefault(); } if (event.key === 'ArrowLeft') { show(current - 1); event.preventDefault(); } if (event.key === 'Escape') hide();});})();</script>
</figure>

**Part 1b: the domain of dependence.** Is it really about information? I took $X = 0.8$ and $T = 0.4$, so the foot of the characteristic is at $0.4$. I ran upwind twice: once from zero initial data, and once from a narrow bump (half-width $0.025$, height $1$) placed exactly at the foot. The scheme is linear, so the difference of the two results at $X$ is the scheme's response to the bump. The exact response is $1$.

| cells | 50 | 100 | 200 | 400 |
|---|---|---|---|---|
| $\nu = 0.8$: left edge of the triangle | $0.290$ | $0.295$ | $0.303$ | $0.299$ |
| $\nu = 0.8$: change of $U$ at $X$ | $0.229$ | $0.333$ | $0.452$ | $0.590$ |
| $\nu = 1.25$: left edge of the triangle | $0.470$ | $0.475$ | $0.483$ | $0.479$ |
| $\nu = 1.25$: change of $U$ at $X$ | $0$ | $0$ | $0$ | $0$ |

(The foot moves a little with the grid, because $X$ is the nearest cell center: $0.390$, $0.395$, $0.4025$, $0.39875$.)

- With $\nu = 0.8$ the foot is inside the triangle. The scheme sees the bump. Its response is smeared, but it moves toward the right value $1$ as the grid is refined.
- With $\nu = 1.25$ the triangle starts at about $0.48$, to the right of the foot. The response is **exactly zero**, on every grid. No refinement can fix this. The scheme computes $U(X, T)$ without the one value it really needs.

**Part 2: four numerical fluxes.** First, one step of each flux in conservation form with $f(q) = q$ and $\nu = 0.5$, compared with the finite difference schemes we already know (data $\sin 2\pi x + 0.3\cos 6\pi x$, $100$ cells). Largest difference: central vs. FTCS $2.2 \times 10^{-16}$, Lax–Friedrichs flux vs. Lax–Friedrichs $2.2 \times 10^{-16}$, Richtmyer vs. Lax–Wendroff $2.2 \times 10^{-16}$, upwind flux vs. upwind $0$. So the algebra of section 7 holds.

Then I moved the bump once around the domain with each flux. The largest $|Q|$ (exact: $0.9938$ at all times):

| $\nu = 0.5$ ($200$ steps) | $t = 0.25$ | $t = 0.5$ | $t = 1$ | error at $t = 1$ |
|---|---|---|---|---|
| central | $1.868$ | $141.6$ | $8.30 \times 10^{6}$ | $8.30 \times 10^{6}$ |
| Lax–Friedrichs | $0.5568$ | $0.4238$ | $0.3120$ | $0.683$ |
| Richtmyer | $0.9947$ | $0.9917$ | $0.9621$ | $0.188$ |
| upwind | $0.7656$ | $0.6373$ | $0.4992$ | $0.495$ |

| $\nu = 0.1$ ($1000$ steps) | $t = 0.25$ | $t = 0.5$ | $t = 1$ | error at $t = 1$ |
|---|---|---|---|---|
| central | $1.080$ | $1.178$ | $1.282$ | $0.702$ |
| Lax–Friedrichs | $0.2468$ | $0.1769$ | $0.1276$ | $0.866$ |
| Richtmyer | $1.005$ | $0.9983$ | $0.9554$ | $0.265$ |
| upwind | $0.6605$ | $0.5213$ | $0.3924$ | $0.603$ |

- Both runs satisfy the CFL condition, but the **central flux grows** in both. A smaller step only slows the growth: the worst mode grows by $\sqrt{1 + \nu^2}$ per step, which is $1.118$ for $\nu = 0.5$ and $1.005$ for $\nu = 0.1$. There are more steps, so it still grows.
- **Lax–Friedrichs** is stable but smears the most, and much more at $\nu = 0.1$. Its diffusion $\beta = \Delta x^2/(2\Delta t)$ grows when $\Delta t$ shrinks.
- **Richtmyer** (= Lax–Wendroff here) keeps the height best. It dips below zero (smallest value $-0.118$ at $\nu = 0.5$): it is not a convex combination.
- **Upwind** stays between $0$ and $1$ but smears.
- All four keep the total amount $\Delta x\sum_i Q_i$ to rounding (change at most $5.6 \times 10^{-17}$; for the central flux $1.8 \times 10^{-11}$ while the values are about $10^{7}$). Conservation form gives conservation, but not stability.

**Part 3: the maximum principle for upwind.** I started from a box ($1$ on $[0.2, 0.45)$, $0$ elsewhere) and ran $200$ steps.

| $\nu$ | weights | after 1 step: max / min | over 200 steps: max / min |
|---|---|---|---|
| $0.5$ | $0.5,\ 0.5$ | $1$ / $0$ | $1$ / $0$ |
| $0.9$ | $0.9,\ 0.1$ | $1$ / $0$ | $1$ / $0$ |
| $1.0$ | $1,\ 0$ | $1$ / $0$ | $1$ / $0$ |
| $1.05$ | $1.05,\ -0.05$ | $1.05$ / $-0.05$ | $1.31 \times 10^{7}$ / $-1.31 \times 10^{7}$ |

For $0 \le \nu \le 1$ no value ever goes above $1$ or below $0$. For $\nu = 1.05$ the very first step already breaks the bounds by exactly $0.05$ at the two jumps, as the weights predict. Then the largest value grows: $2.65$ after $20$ steps, $16.8$ after $50$, $1337$ after $100$.

**Part 4 (my own addition): hyperbolic or not.** I ran Lax–Friedrichs on $\mathbf u_t + A\mathbf u_x = 0$ for the two matrices of section 2, from $u = $ bump (centered at $0.5$, half-width $0.1$) and $v = 0$, to $t = 0.2$. The time step was set so that (largest $|\lambda|$) $\times\,\Delta t/\Delta x = 0.8$.

| cells | 25 | 50 | 100 | 200 | 400 |
|---|---|---|---|---|---|
| $A$ (eigenvalues $\pm 2$): error | $0.513$ | $0.432$ | $0.297$ | $0.179$ | $0.0996$ |
| $B$ (eigenvalues $\pm i$): largest $\lvert U\rvert$ | $1.42$ | $4.86$ | $24.6$ | $800$ | $2.51 \times 10^{7}$ |

(At $25$ and $50$ cells the last step lands at $t = 0.192$.)

- For $A$, the exact solution is two half-bumps moving apart with speeds $\pm 2$, from the eigenvectors $(1, 2)$ and $(1, -2)$. Lax–Friedrichs converges to it. The error ratio per refinement rises toward $2$ ($1.19$, $1.45$, $1.66$, $1.80$), as expected for a first-order scheme.
- For $B$, the result gets **worse** with refinement, and fast. A finer grid can represent shorter waves, and for this problem short waves grow fastest ($e^{|\xi| t}$). This is not a fault of the scheme: the problem itself is not well posed as an initial value problem. Real eigenvalues are what make "waves with speeds" possible.

The full script and the chart code are in [this site's repository](https://github.com/doubleyands/doubleyands.github.io/tree/main/scripts/posts/numerical-pde/lecture-14).

## Summary

1. A conservation law says: the change inside a region equals the flux through its boundary. In differential form, $u_t + \nabla\cdot\mathbf f(u) = 0$.
2. A system $\mathbf u_t + \mathbf f(\mathbf u)_x = 0$ is **hyperbolic** if the Jacobian $D\mathbf f$ is diagonalizable with real eigenvalues. Then it splits into advection equations, and the eigenvalues are wave speeds.
3. In primitive form $u_t + f'(u)\,u_x = 0$, the wave speed is $f'(u)$. Along a characteristic $\dot X = f'(u)$, the total derivative of $u$ is zero, so $u$ is constant there.
4. The **CFL condition**: the numerical domain of dependence must contain the physical one. For upwind this is $0 \le \nu = c\Delta t/\Delta x \le 1$. It is necessary, not sufficient.
5. Four numerical fluxes: **central** (FTCS, unstable), **Lax–Friedrichs** (central plus numerical diffusion $\beta = \Delta x^2/2\Delta t$), **Richtmyer** (two steps; Lax–Wendroff for linear flux), and **upwind** (a convex combination when $0 \le \nu \le 1$, so it obeys a maximum principle).
