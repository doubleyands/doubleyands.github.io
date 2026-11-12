---
title: "Numerical PDEs, Lecture 12: Characteristics, Imaginary Eigenvalues, and Artificial Diffusion"
description: "My notes on the twelfth lecture: characteristics and inflow boundaries for the advection equation, why the central difference gives purely imaginary eigenvalues, why Forward Euler fails and leapfrog is only marginally stable, and how Lax–Friedrichs, Lax–Wendroff and upwind all add artificial diffusion."
pubDate: 2026-11-12
project: numerical-pde
tags: [numerical PDEs, advection, hyperbolic, stability, artificial diffusion, notes]
---

This is my summary of the twelfth lecture in my numerical PDE course. I write the ideas in my own words and check them with a small experiment.

[Lecture 11](/blog/2026-11-10-numerical-pde-lecture-11/) met the advection equation for the first time. It stated two facts without proof: FTCS is unstable for it, and Lax–Friedrichs is stable when $|ak/h| \le 1$. This lecture explains both facts with the eigenvalues of the space matrix, the method of [Lecture 10](/blog/2026-11-05-numerical-pde-lecture-10/). The main message is one idea: **artificial diffusion**. Lax–Friedrichs, Lax–Wendroff and upwind all work because they quietly add a small $u_{xx}$ term.

Notation as before: $h$ is the grid spacing, $k$ is the time step, $U^n_j$ is the value at $x_j = jh$ and $t_n = nk$, and $\nu = ak/h$ is the Courant number. In exponents, $i = \sqrt{-1}$ and $j$ is the grid index.

## 1. Characteristics and the foot

The problem is

$$
u_t + a\,u_x = 0, \qquad u(x, 0) = \eta(x),
$$

with the solution $u(x, t) = \eta(x - at)$. The speed is $a$. If $a > 0$ the shape moves right, and if $a < 0$ it moves left.

The solution depends on $x$ and $t$ only through the combination $x - at$. So it is constant on each line

$$
x - at = x_0 \quad(\text{a constant}).
$$

These lines are the **characteristics**. In the $x$–$t$ plane (with $t$ up), they have slope $1/a$, not $a$, because time is on the vertical axis. All of them are parallel. They differ only in where they start.

To find $u$ at a point $(X, T)$, follow the characteristic back down to $t = 0$. It hits the axis at $x_0 = X - aT$, the **foot of the characteristic**, and the answer is $u(X, T) = \eta(x_0)$.

Compare this with the heat equation. Diffusion spreads a shape out and makes it flatter. Advection does not spread anything. It carries the shape along without any change, like an object that floats on a moving stream.

## 2. Boundary conditions: only at the inflow

In practice we solve on $0 \le x \le 1$. The equation has only one space derivative, so it takes only **one** boundary condition. Two conditions would be too many, and in general no solution would satisfy both.

The condition goes where the flow comes **in**:

| | inflow boundary (give a value) | outflow boundary (give nothing) |
|---|---|---|
| $a > 0$ | $x = 0$: $u(0, t) = g_0(t)$ | $x = 1$ |
| $a < 0$ | $x = 1$ | $x = 0$ |

The value at the outflow boundary is already decided by what flows there from inside. In the characteristic picture: for $a > 0$, a characteristic that reaches the left edge before $t = 0$ takes its value from $g_0$, not from $\eta$.

## 3. Periodic boundaries and a circulant matrix

For analysis, the nicest choice is a **periodic** boundary, $u(0, t) = u(1, t)$. With $h = 1/(m+1)$, the unknowns are $U_1, \dots, U_{m+1}$, and the index wraps around: $U_0 = U_{m+1}$ and $U_{m+2} = U_1$.

Discretize only in space with the central difference (the method of lines):

$$
U_j'(t) = -\frac{a}{2h}\big(U_{j+1}(t) - U_{j-1}(t)\big).
$$

This is $U'(t) = AU(t)$ with

$$
A = -\frac{a}{2h}\begin{bmatrix}0&1&&&-1\\-1&0&1&&\\&\ddots&\ddots&\ddots&\\&&-1&0&1\\1&&&-1&0\end{bmatrix}.
$$

The entries in the corners come from the wrap-around. Every row is the row above, shifted one place to the right. A matrix like this is called **circulant**.

## 4. Skew-symmetric, so the eigenvalues are purely imaginary

Why look at this ODE system at all? Because of **stability**. If we know the eigenvalues $\lambda_p$ of $A$, we split the system into scalar problems $u' = \lambda_p u$ and check if every $z = k\lambda_p$ lies in the stability region of the time stepper. This is what we did for the heat equation in Lecture 10.

Look at $A$: the entry above the diagonal is the negative of the entry below it. So $A^{\mathsf T} = -A$. The matrix is **skew-symmetric**. Two facts from linear algebra:

- A **symmetric** matrix has only **real** eigenvalues. This was the heat equation: all $\lambda_p < 0$.
- A **skew-symmetric** matrix has only **purely imaginary** eigenvalues. This is advection.

The eigenvectors are the grid Fourier modes $u^p_j = e^{2\pi i p j h}$. Putting one in gives the same computation as for $D_0$ in Lecture 11 (with $\xi = 2\pi p$):

$$
\lambda_p = -\frac{ia}{h}\sin(2\pi p h), \qquad p = 1, \dots, m+1.
$$

All of them lie on the imaginary axis, between $-ia/h$ and $ia/h$.

## 5. Forward Euler cannot work

Forward Euler in time plus the central difference in space is FTCS. Its stability region is the disk $|1 + z| \le 1$: center $-1$, radius $1$ ([Lecture 9](/blog/2026-11-03-numerical-pde-lecture-9/)). This disk touches the imaginary axis at **one point only**, $z = 0$. Every other point $z = k\lambda_p$ on the imaginary axis has $|1 + z| = \sqrt{1 + |z|^2} > 1$.

So making $k/h$ smaller does not help. The points move closer to $0$, but they never get inside. For any fixed $k/h$, FTCS is unstable. This is different from the stiff example with $\lambda = -2100$ in [Lecture 8](/blog/2026-10-29-numerical-pde-lecture-8/). That eigenvalue was real, so a small enough $k$ put $k\lambda$ inside the disk. Here that is impossible.

There is one way out, but it is expensive. If $k$ goes to zero much faster than $h$, for example $k = h^2$, then

$$
|1 + k\lambda_p|^2 \le 1 + \Big(\frac{ak}{h}\Big)^2 = 1 + a^2k,
$$

so each step grows by at most a factor $1 + \alpha k$ with $\alpha = a^2/2$. That is Lax–Richtmyer stability from Lecture 10, and the method converges, with $\|U^n\| \le e^{a^2T/2}\|U^0\|$. But $k = h^2$ needs far too many steps. Part 4 of the experiment shows both cases.

## 6. Leapfrog: stable, but just barely

The leapfrog scheme uses a centered difference in time too:

$$
\frac{U^{n+1} - U^{n-1}}{2k} = AU^n \quad\Longrightarrow\quad U^{n+1}_j = U^{n-1}_j - \frac{ak}{h}\big(U^n_{j+1} - U^n_{j-1}\big).
$$

It is explicit, uses three time levels, and is second order in both $k$ and $h$.

Its stability region is the open segment of the imaginary axis between $-i$ and $i$ ([Lecture 9](/blog/2026-11-03-numerical-pde-lecture-9/)). Our $k\lambda_p$ lie exactly on the imaginary axis, with $|k\lambda_p| \le |\nu|$. So leapfrog is stable if

$$
\left|\frac{ak}{h}\right| < 1.
$$

The inequality must be strict. At $z = \pm i$ the two roots of the characteristic polynomial come together as a double root of size $1$, and then the error grows linearly in $n$.

Each stability region is different. Forward Euler has a disk, Backward Euler has the outside of a disk, the trapezoid rule has the left half-plane, and leapfrog has only a segment. It is important not to mix them up.

**Marginally stable.** The leapfrog region has no inside: it is a segment. So every $k\lambda_p$ sits on the edge of the region. Modes neither grow nor decay. For example, at $z = 0.6i$ the characteristic polynomial $\zeta^2 - 2z\zeta - 1$ has roots $\zeta = 0.6i \pm 0.8$, and both have $|\zeta| = 1$ exactly.

- **Good:** the scheme is **nondissipative**. The true advection equation does not lose energy either: the shape only moves.
- **Bad (1):** no mode is damped, but the modes do not all move at the right speed. Near a sharp jump, the fast modes fall behind and leave oscillations.
- **Bad (2):** this kind of stability is very narrow. With variable coefficients or a nonlinear equation, a small change can push it into instability.

## 7. Lax–Friedrichs is FTCS plus artificial diffusion

Lax–Friedrichs replaces $U^n_j$ in FTCS by the average of its neighbors. Why does that save the scheme? One identity explains it:

$$
\frac12\big(U^n_{j-1} + U^n_{j+1}\big) = U^n_j + \frac12\big(U^n_{j-1} - 2U^n_j + U^n_{j+1}\big).
$$

(The two $U^n_j$ on the right cancel.) Put this into Lax–Friedrichs and divide by $k$:

$$
\frac{U^{n+1}_j - U^n_j}{k} + a\,\frac{U^n_{j+1} - U^n_{j-1}}{2h} = \frac{h^2}{2k}\cdot\frac{U^n_{j-1} - 2U^n_j + U^n_{j+1}}{h^2}.
$$

The left side is exactly FTCS. The right side is the 3-point difference for $u_{xx}$. So Lax–Friedrichs is FTCS applied to

$$
u_t + a\,u_x = \epsilon\,u_{xx}, \qquad \epsilon = \frac{h^2}{2k}.
$$

The number $\epsilon$ is not part of the physics. It comes only from the grid. It is called **artificial diffusion** (or **artificial viscosity**). With $k/h$ fixed, $\epsilon = ah/(2\nu) \to 0$ as $h \to 0$, so the scheme is still consistent with the original equation.

This term is what makes the scheme work: it smooths the solution and damps the modes that FTCS lets grow. But if $\epsilon$ is too large, it flattens the solution too much, and the result is no longer physical. We saw the same idea for a steady problem in [Lecture 4](/blog/2026-10-15-numerical-pde-lecture-4/): upwinding was the central difference plus a diffusion term of size $h/2$.

**The matrix view.** Lax–Friedrichs is Forward Euler applied to $U' = A_\epsilon U$, where we discretize the new equation in space in the same way:

$$
A_\epsilon = \underbrace{-\frac{a}{2h}\begin{bmatrix}0&1&&-1\\-1&0&\ddots&\\&\ddots&\ddots&1\\1&&-1&0\end{bmatrix}}_{\text{skew-symmetric}} + \underbrace{\frac{\epsilon}{h^2}\begin{bmatrix}-2&1&&1\\1&-2&\ddots&\\&\ddots&\ddots&1\\1&&1&-2\end{bmatrix}}_{\text{symmetric}}.
$$

## 8. The eigenvalues move left, onto an ellipse

Both matrices in $A_\epsilon$ are circulant, so they have the **same** eigenvectors $e^{2\pi i p j h}$. In general we cannot add the eigenvalues of two matrices, but with shared eigenvectors we can:

$$
\lambda_p = -\frac{ia}{h}\sin(2\pi p h) - \frac{2\epsilon}{h^2}\big(1 - \cos(2\pi p h)\big).
$$

The second term is never positive. So the diffusion moves every eigenvalue to the **left** of the imaginary axis, toward the inside of the Forward Euler disk.

**Why an ellipse.** Multiply by $k$ and write $\theta = 2\pi p h$, $c = 2k\epsilon/h^2$:

$$
k\lambda_p = \underbrace{-c + c\cos\theta}_{x} + i\,\underbrace{(-\nu\sin\theta)}_{y}.
$$

Then $\cos\theta = (x + c)/c$ and $\sin\theta = -y/\nu$, and $\cos^2\theta + \sin^2\theta = 1$ gives

$$
\Big(\frac{x + c}{c}\Big)^2 + \Big(\frac{y}{\nu}\Big)^2 = 1.
$$

This is an ellipse with center $-c$, half-width $c$ along the real axis, and half-height $|\nu|$ along the imaginary axis.

For **Lax–Friedrichs**, $\epsilon = h^2/(2k)$ gives $c = 1$. The ellipse has center $-1$ and half-width $1$, the same as the Forward Euler disk in the real direction. If the half-height $|\nu|$ is at most $1$, the ellipse is inside the disk. This is the CFL condition $|ak/h| \le 1$ from Lecture 11, now with a proof.

**Not too little, not too much.** For which $c$ is the ellipse inside the disk? Write $w = 1 - \cos\theta \in [0, 2]$, so $\sin^2\theta = w(2 - w)$. Then

$$
|1 + k\lambda_p|^2 = (1 - cw)^2 + \nu^2 w(2 - w) = 1 + w\big[\,2(\nu^2 - c) + (c^2 - \nu^2)\,w\,\big].
$$

We need the bracket to be $\le 0$ for all $w$ in $[0, 2]$. It is linear in $w$, so it is enough to check the two ends. Near $w = 0$ we need $\nu^2 \le c$. At $w = 2$ we need $2(c^2 - c) \le 0$, so $c \le 1$. Together:

$$
\nu^2 \le c \le 1 \quad\Longleftrightarrow\quad \frac{a^2k}{2} \le \epsilon \le \frac{h^2}{2k}.
$$

- If $c < \nu^2$ (too little diffusion), the ellipse leaks out of the disk above and below the origin.
- If $c > 1$ (too much diffusion), its left end goes past $-2$ and leaves the disk on the left.
- Lax–Friedrichs is the top end, $c = 1$. The bottom end, $c = \nu^2$, is the smallest diffusion that still gives stability. It is the next scheme.

Part 2 of the experiment counts the eigenvalues outside the disk for several values of $c$.

## 9. Lax–Wendroff: a Taylor method in time

There are several ways to get second order in time, but each has a cost. Leapfrog needs three levels and is nondissipative. Crank–Nicolson is implicit, and for hyperbolic problems that cost is not necessary. A two-stage Runge–Kutta method needs extra storage and care near the boundary. The **Taylor method** gives a two-level explicit scheme.

The idea is to use the equation again, as in the Taylor methods for ODEs ([Lecture 6](/blog/2026-10-22-numerical-pde-lecture-6/)). From $u_t = -a\,u_x$ we get $u_{tt} = -a\,u_{xt} = a^2u_{xx}$. So

$$
u(x, t + k) = u - ka\,u_x + \tfrac12 k^2 a^2\,u_{xx} + O(k^3).
$$

All time derivatives are now space derivatives. Replace them by the standard central differences:

$$
\boxed{\;U^{n+1}_j = U^n_j - \frac{ak}{2h}\big(U^n_{j+1} - U^n_{j-1}\big) + \frac{a^2k^2}{2h^2}\big(U^n_{j-1} - 2U^n_j + U^n_{j+1}\big)\;}
$$

This is the **Lax–Wendroff** scheme.

**The matrix version is wider.** We can also write $U'' = AU' = A^2U$ and take $U^{n+1} = U^n + kAU^n + \frac12k^2A^2U^n$. But applying $A$ twice doubles the stencil. The last term becomes $\frac{a^2k^2}{8h^2}\big(U^n_{j-2} - 2U^n_j + U^n_{j+2}\big)$, a second difference with spacing $2h$. Lax–Wendroff uses the usual 3-point second difference instead. A 5-point stencil would need extra numerical boundary values near the ends.

**Artificial diffusion again.** The last term has the form of $A_\epsilon$ with $\epsilon = a^2k/2$. Then $c = 2k\epsilon/h^2 = \nu^2$, exactly the bottom end of the range in section 8. So

$$
k\lambda_p = -i\nu\sin\theta + \nu^2(\cos\theta - 1),
$$

an ellipse with center $-\nu^2$ and half-axes $\nu^2$ and $|\nu|$. It is stable when $|ak/h| \le 1$, the same condition as Lax–Friedrichs.

## 10. Upwind is artificial diffusion too

The upwind idea from [Lecture 4](/blog/2026-10-15-numerical-pde-lecture-4/): take the information from the side the flow comes from. For $a > 0$, use $u_x \approx (U_j - U_{j-1})/h$:

$$
U^{n+1}_j = U^n_j - \frac{ak}{h}\big(U^n_j - U^n_{j-1}\big).
$$

Add and subtract $U^n_{j+1}$ in the right way, and it becomes

$$
U^{n+1}_j = U^n_j - \frac{ak}{2h}\big(U^n_{j+1} - U^n_{j-1}\big) + \frac{ak}{2h}\big(U^n_{j+1} - 2U^n_j + U^n_{j-1}\big).
$$

Again: central difference plus diffusion, with $\epsilon = ah/2$. Then $c = \nu$, and the condition $\nu^2 \le \nu \le 1$ gives $0 \le ak/h \le 1$.

**The wrong side.** If $a < 0$ and we still use $(U_j - U_{j-1})/h$, then $\epsilon = ah/2 < 0$. That is **negative diffusion**, like running the heat equation backward in time. The eigenvalues move to the right half-plane, and the scheme is unstable for every $k$. For $a < 0$ we must use the other one-sided difference, $(U_{j+1} - U_j)/h$, and then the condition is $-1 \le ak/h \le 0$.

| scheme | artificial diffusion $\epsilon$ | $\epsilon$ in units of $ah/2$ | order | stability |
|---|---|---|---|---|
| FTCS | $0$ | $0$ | 1 in time, 2 in space | never |
| Lax–Friedrichs | $h^2/(2k)$ | $1/\nu$ | 1 | $\lvert\nu\rvert \le 1$ |
| upwind ($a > 0$) | $ah/2$ | $1$ | 1 | $0 \le \nu \le 1$ |
| Lax–Wendroff | $a^2k/2$ | $\nu$ | 2 | $\lvert\nu\rvert \le 1$ |

For $0 < \nu < 1$ the order is clear: Lax–Friedrichs smears the most, then upwind, then Lax–Wendroff. Part 3 of the experiment shows what this means for a real solution.

## 11. Experiment

The script uses the Python standard library only (`cmath` for complex numbers). Everything is on a periodic grid with $a = 1$.

**Part 1: imaginary eigenvalues.** For $40$ grid points I built $A$ as a dense matrix. $A^{\mathsf T} + A$ is exactly zero. For each Fourier mode, $Au^p - \lambda_p u^p$ is at most $2.8 \times 10^{-14}$ times $a/h$. I also computed each eigenvalue from the Rayleigh quotient $\bar u^{\mathsf T}Au/(\bar u^{\mathsf T}u)$, without the formula. The largest real part is $9.4 \times 10^{-18}$ times $a/h$, so the eigenvalues are purely imaginary up to rounding. Their imaginary parts fill $[-a/h, a/h]$.

For Forward Euler, the largest growth factor $\max_p|1 + k\lambda_p|$ is above $1$ for every Courant number:

| $\nu$ | 0.1 | 0.5 | 0.9 |
|---|---|---|---|
| $\max_p\lvert 1 + k\lambda_p\rvert$ | $1.0050$ | $1.1180$ | $1.3454$ |

These are $\sqrt{1 + \nu^2}$, from the mode with $\sin\theta = \pm 1$.

**Part 2: Lax–Friedrichs is FTCS plus diffusion, and the ellipse.** With $40$ points and $\nu = 0.6$, one step of Lax–Friedrichs and one step of FTCS for $u_t + u_x = \epsilon u_{xx}$ with $\epsilon = h^2/(2k)$ differ by at most $2.2 \times 10^{-16}$. They are the same scheme.

Then, for the same $\nu = 0.6$, I counted how many of the $40$ scaled eigenvalues $k\lambda_p$ are outside the Forward Euler disk for different diffusion levels $c$. Theory says none are outside exactly when $\nu^2 = 0.36 \le c \le 1$.

| $c$ | 0 (FTCS) | 0.15 | 0.30 | 0.36 (Lax–Wendroff) | 0.6 (upwind) | 0.8 | 1.0 (Lax–Friedrichs) | 1.1 | 1.3 |
|---|---|---|---|---|---|---|---|---|---|
| outside the disk | 38 | 22 | 12 | 0 | 0 | 0 | 0 | 9 | 15 |
| $\max_p\lvert 1 + k\lambda_p\rvert$ | $1.166$ | $1.063$ | $1.007$ | $1.000$ | $1.000$ | $1.000$ | $1.000$ | $1.200$ | $1.600$ |

- With too little diffusion, the modes near the middle of the ellipse (where $\sin\theta$ is large) are outside. With more diffusion, fewer are outside, but even $c = 0.30$ still leaves $12$.
- From $c = 0.36$ to $c = 1$ all $40$ are inside, as the inequality $\nu^2 \le c \le 1$ says.
- With too much diffusion, the sawtooth mode goes out on the left: at $c = 1.3$ it has $|1 + k\lambda| = |1 - 2c| = 1.6$.
- FTCS has $38$ outside, not $40$: the two modes with $\sin\theta = 0$ ($p = 0$ and the sawtooth) sit at $z = 0$.

**Part 3: one full trip.** The grid has $120$ points, $\nu = 0.5$, and I take $240$ steps, so the shape goes once around the domain and the exact solution is again the initial data. The data is a narrow smooth bump $e^{-300(x - 0.25)^2}$ plus a box of height $1$ on $[0.55, 0.8)$. Leapfrog needs two starting levels, so its first step is one Lax–Wendroff step.

<figure class="viz">
<div class="viz-legend"><span class="viz-key"><span class="viz-swatch" style="background:#2a78d6"></span>Lax–Friedrichs</span><span class="viz-key"><span class="viz-swatch" style="background:#eb6834"></span>upwind</span><span class="viz-key"><span class="viz-swatch" style="background:#1baf7a"></span>Lax–Wendroff</span><span class="viz-key"><span class="viz-swatch" style="background:#eda100"></span>leapfrog</span><span class="viz-key"><span class="viz-swatch" style="background:#86868b"></span>exact</span></div>
<div class="viz-frame">
<svg viewBox="0 0 680 400" role="img" tabindex="0" aria-label="Solution of u_t + u_x = 0 after one full trip around a periodic domain, for Lax-Friedrichs, upwind, Lax-Wendroff and leapfrog, with the exact solution in gray. The exact solution is a narrow bump at x = 0.25 and a box from 0.55 to 0.8. Lax-Friedrichs and upwind smear both shapes. Lax-Wendroff and leapfrog keep the bump high but overshoot near the box and leave wiggles." class="viz-plot"><line x1="56" x2="540" y1="44.1" y2="44.1" class="viz-grid"/><text x="46" y="48.1" text-anchor="end" class="viz-tick">1.5</text><line x1="56" x2="540" y1="114.5" y2="114.5" class="viz-grid"/><text x="46" y="118.5" text-anchor="end" class="viz-tick">1</text><line x1="56" x2="540" y1="185.0" y2="185.0" class="viz-grid"/><text x="46" y="189.0" text-anchor="end" class="viz-tick">0.5</text><line x1="56" x2="540" y1="255.5" y2="255.5" class="viz-leader"/><text x="46" y="259.5" text-anchor="end" class="viz-tick">0</text><line x1="56" x2="540" y1="325.9" y2="325.9" class="viz-grid"/><text x="46" y="329.9" text-anchor="end" class="viz-tick">−0.5</text><text x="56.0" y="362" text-anchor="middle" class="viz-tick">0</text><text x="177.0" y="362" text-anchor="middle" class="viz-tick">0.25</text><text x="298.0" y="362" text-anchor="middle" class="viz-tick">0.5</text><text x="419.0" y="362" text-anchor="middle" class="viz-tick">0.75</text><text x="540.0" y="362" text-anchor="middle" class="viz-tick">1</text><text x="298.0" y="388" text-anchor="middle" class="viz-axis-title">Position x</text><text x="56" y="14" class="viz-axis-title">Solution u after one trip (t = 1); gray: exact</text><path d="M56.0,255.5 L60.0,255.5 L64.1,255.5 L68.1,255.5 L72.1,255.5 L76.2,255.5 L80.2,255.5 L84.2,255.5 L88.3,255.4 L92.3,255.4 L96.3,255.4 L100.4,255.4 L104.4,255.3 L108.4,255.1 L112.5,254.8 L116.5,254.2 L120.5,253.1 L124.6,251.3 L128.6,248.4 L132.6,244.1 L136.7,237.9 L140.7,229.4 L144.7,218.3 L148.8,204.7 L152.8,188.9 L156.8,171.8 L160.9,154.5 L164.9,138.6 L168.9,125.8 L173.0,117.5 L177.0,114.5 L181.0,117.5 L185.1,125.8 L189.1,138.6 L193.1,154.5 L197.2,171.8 L201.2,188.9 L205.2,204.7 L209.3,218.3 L213.3,229.4 L217.3,237.9 L221.4,244.1 L225.4,248.4 L229.4,251.3 L233.5,253.1 L237.5,254.2 L241.5,254.8 L245.6,255.1 L249.6,255.3 L253.6,255.4 L257.7,255.4 L261.7,255.4 L265.7,255.4 L269.8,255.5 L273.8,255.5 L277.8,255.5 L281.9,255.5 L285.9,255.5 L289.9,255.5 L294.0,255.5 L298.0,255.5 L302.0,255.5 L306.1,255.5 L310.1,255.5 L314.1,255.5 L318.2,255.5 L322.2,114.5 L326.2,114.5 L330.3,114.5 L334.3,114.5 L338.3,114.5 L342.4,114.5 L346.4,114.5 L350.4,114.5 L354.5,114.5 L358.5,114.5 L362.5,114.5 L366.6,114.5 L370.6,114.5 L374.6,114.5 L378.7,114.5 L382.7,114.5 L386.7,114.5 L390.8,114.5 L394.8,114.5 L398.8,114.5 L402.9,114.5 L406.9,114.5 L410.9,114.5 L415.0,114.5 L419.0,114.5 L423.0,114.5 L427.1,114.5 L431.1,114.5 L435.1,114.5 L439.2,114.5 L443.2,255.5 L447.2,255.5 L451.3,255.5 L455.3,255.5 L459.3,255.5 L463.4,255.5 L467.4,255.5 L471.4,255.5 L475.5,255.5 L479.5,255.5 L483.5,255.5 L487.6,255.5 L491.6,255.5 L495.6,255.5 L499.7,255.5 L503.7,255.5 L507.7,255.5 L511.8,255.5 L515.8,255.5 L519.8,255.5 L523.9,255.5 L527.9,255.5 L531.9,255.5 L536.0,255.5" fill="none" stroke="#86868b" stroke-width="2" stroke-dasharray="5 4" stroke-linejoin="round"/><path d="M56.0,245.9 L60.0,245.1 L64.1,245.5 L68.1,244.5 L72.1,244.3 L76.2,243.1 L80.2,242.5 L84.2,241.1 L88.3,240.0 L92.3,238.4 L96.3,237.1 L100.4,235.3 L104.4,233.7 L108.4,231.8 L112.5,229.9 L116.5,228.0 L120.5,226.0 L124.6,224.0 L128.6,222.0 L132.6,220.1 L136.7,218.2 L140.7,216.4 L144.7,214.6 L148.8,213.1 L152.8,211.6 L156.8,210.3 L160.9,209.1 L164.9,208.2 L168.9,207.3 L173.0,206.9 L177.0,206.4 L181.0,206.4 L185.1,206.3 L189.1,206.8 L193.1,207.0 L197.2,207.9 L201.2,208.3 L205.2,209.6 L209.3,210.1 L213.3,211.8 L217.3,212.3 L221.4,214.2 L225.4,214.4 L229.4,216.5 L233.5,216.4 L237.5,218.5 L241.5,217.8 L245.6,219.9 L249.6,218.6 L253.6,220.6 L257.7,218.5 L261.7,220.3 L265.7,217.3 L269.8,218.9 L273.8,215.0 L277.8,216.3 L281.9,211.5 L285.9,212.7 L289.9,207.0 L294.0,208.0 L298.0,201.6 L302.0,202.4 L306.1,195.5 L310.1,196.1 L314.1,188.9 L318.2,189.3 L322.2,182.1 L326.2,182.4 L330.3,175.3 L334.3,175.6 L338.3,169.0 L342.4,169.1 L346.4,163.3 L350.4,163.4 L354.5,158.5 L358.5,158.5 L362.5,154.8 L366.6,154.9 L370.6,152.5 L374.6,152.5 L378.7,151.6 L382.7,151.6 L386.7,152.1 L390.8,152.2 L394.8,154.2 L398.8,154.2 L402.9,157.6 L406.9,157.6 L410.9,162.2 L415.0,162.2 L419.0,168.0 L423.0,168.0 L427.1,174.7 L431.1,174.7 L435.1,182.0 L439.2,182.0 L443.2,189.7 L447.2,189.7 L451.3,197.5 L455.3,197.5 L459.3,205.3 L463.4,205.2 L467.4,212.7 L471.4,212.6 L475.5,219.6 L479.5,219.5 L483.5,225.8 L487.6,225.7 L491.6,231.3 L495.6,231.2 L499.7,235.9 L503.7,235.7 L507.7,239.7 L511.8,239.4 L515.8,242.5 L519.8,242.1 L523.9,244.5 L527.9,244.0 L531.9,245.6 L536.0,245.0" fill="none" stroke="#2a78d6" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><path d="M56.0,255.0 L60.0,254.9 L64.1,254.7 L68.1,254.5 L72.1,254.1 L76.2,253.6 L80.2,253.0 L84.2,252.2 L88.3,251.2 L92.3,250.0 L96.3,248.5 L100.4,246.7 L104.4,244.5 L108.4,242.0 L112.5,239.0 L116.5,235.7 L120.5,232.0 L124.6,227.9 L128.6,223.5 L132.6,218.8 L136.7,213.9 L140.7,208.9 L144.7,204.0 L148.8,199.2 L152.8,194.7 L156.8,190.6 L160.9,187.0 L164.9,184.1 L168.9,181.9 L173.0,180.6 L177.0,180.2 L181.0,180.6 L185.1,181.9 L189.1,184.1 L193.1,187.0 L197.2,190.6 L201.2,194.7 L205.2,199.2 L209.3,204.0 L213.3,208.9 L217.3,213.8 L221.4,218.7 L225.4,223.3 L229.4,227.6 L233.5,231.6 L237.5,235.1 L241.5,238.2 L245.6,240.8 L249.6,242.8 L253.6,244.3 L257.7,245.3 L261.7,245.7 L265.7,245.5 L269.8,244.7 L273.8,243.3 L277.8,241.3 L281.9,238.6 L285.9,235.3 L289.9,231.3 L294.0,226.7 L298.0,221.4 L302.0,215.7 L306.1,209.4 L310.1,202.7 L314.1,195.7 L318.2,188.6 L322.2,181.4 L326.2,174.2 L330.3,167.2 L334.3,160.5 L338.3,154.2 L342.4,148.3 L346.4,143.0 L350.4,138.3 L354.5,134.1 L358.5,130.6 L362.5,127.7 L366.6,125.4 L370.6,123.7 L374.6,122.6 L378.7,122.0 L382.7,122.0 L386.7,122.6 L390.8,123.7 L394.8,125.4 L398.8,127.7 L402.9,130.6 L406.9,134.1 L410.9,138.3 L415.0,143.0 L419.0,148.3 L423.0,154.2 L427.1,160.5 L431.1,167.2 L435.1,174.2 L439.2,181.4 L443.2,188.6 L447.2,195.8 L451.3,202.8 L455.3,209.6 L459.3,215.9 L463.4,221.8 L467.4,227.2 L471.4,232.0 L475.5,236.3 L479.5,240.0 L483.5,243.1 L487.6,245.8 L491.6,248.0 L495.6,249.7 L499.7,251.2 L503.7,252.3 L507.7,253.1 L511.8,253.8 L515.8,254.2 L519.8,254.6 L523.9,254.8 L527.9,255.0 L531.9,255.0 L536.0,255.0" fill="none" stroke="#eb6834" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><path d="M56.0,255.4 L60.0,255.5 L64.1,255.6 L68.1,255.6 L72.1,255.6 L76.2,255.5 L80.2,255.3 L84.2,254.9 L88.3,254.6 L92.3,254.4 L96.3,254.6 L100.4,255.4 L104.4,256.8 L108.4,258.7 L112.5,260.7 L116.5,262.3 L120.5,262.6 L124.6,260.7 L128.6,255.9 L132.6,247.7 L136.7,235.8 L140.7,220.8 L144.7,203.4 L148.8,184.8 L152.8,166.4 L156.8,149.6 L160.9,136.0 L164.9,126.6 L168.9,122.0 L173.0,122.0 L177.0,126.2 L181.0,133.7 L185.1,143.8 L189.1,156.0 L193.1,169.6 L197.2,183.9 L201.2,197.8 L205.2,210.1 L209.3,219.9 L213.3,226.9 L217.3,231.8 L221.4,235.9 L225.4,240.7 L229.4,246.5 L233.5,252.8 L237.5,258.1 L241.5,260.7 L245.6,259.6 L249.6,255.0 L253.6,248.6 L257.7,242.7 L261.7,240.0 L265.7,242.0 L269.8,249.1 L273.8,259.8 L277.8,271.7 L281.9,281.6 L285.9,286.9 L289.9,286.0 L294.0,278.5 L298.0,265.4 L302.0,248.4 L306.1,229.0 L310.1,208.9 L314.1,189.1 L318.2,170.8 L322.2,154.8 L326.2,142.2 L330.3,133.5 L334.3,128.4 L338.3,125.8 L342.4,124.0 L346.4,121.5 L350.4,117.7 L354.5,112.9 L358.5,108.9 L362.5,107.2 L366.6,108.9 L370.6,114.0 L374.6,120.9 L378.7,126.9 L382.7,129.8 L386.7,127.8 L390.8,120.8 L394.8,110.1 L398.8,98.5 L402.9,88.6 L406.9,83.2 L410.9,83.9 L415.0,91.1 L419.0,104.1 L423.0,121.6 L427.1,141.5 L431.1,162.1 L435.1,181.6 L439.2,199.1 L443.2,213.9 L447.2,225.8 L451.3,234.9 L455.3,241.7 L459.3,246.5 L463.4,249.8 L467.4,252.0 L471.4,253.4 L475.5,254.2 L479.5,254.8 L483.5,255.1 L487.6,255.3 L491.6,255.3 L495.6,255.4 L499.7,255.4 L503.7,255.4 L507.7,255.5 L511.8,255.5 L515.8,255.5 L519.8,255.5 L523.9,255.4 L527.9,255.4 L531.9,255.4 L536.0,255.4" fill="none" stroke="#1baf7a" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><path d="M56.0,246.4 L60.0,248.0 L64.1,264.4 L68.1,246.3 L72.1,258.9 L76.2,270.2 L80.2,248.0 L84.2,252.3 L88.3,249.9 L92.3,245.3 L96.3,270.9 L100.4,269.0 L104.4,244.5 L108.4,240.2 L112.5,251.7 L116.5,272.9 L120.5,281.0 L124.6,268.4 L128.6,249.6 L132.6,230.1 L136.7,233.5 L140.7,242.1 L144.7,212.1 L148.8,180.8 L152.8,168.4 L156.8,139.1 L160.9,120.9 L164.9,130.1 L168.9,120.8 L173.0,114.9 L177.0,142.6 L181.0,149.2 L185.1,132.4 L189.1,141.5 L193.1,157.2 L197.2,175.3 L201.2,215.3 L205.2,237.9 L209.3,233.0 L213.3,223.2 L217.3,208.2 L221.4,215.4 L225.4,241.5 L229.4,251.3 L233.5,259.0 L237.5,270.5 L241.5,266.2 L245.6,253.7 L249.6,261.2 L253.6,265.1 L257.7,231.5 L261.7,216.7 L265.7,232.0 L269.8,238.9 L273.8,253.6 L277.8,284.1 L281.9,313.0 L285.9,298.8 L289.9,275.0 L294.0,274.5 L298.0,257.2 L302.0,243.8 L306.1,229.8 L310.1,228.9 L314.1,190.4 L318.2,174.3 L322.2,133.2 L326.2,133.6 L330.3,117.4 L334.3,140.4 L338.3,126.9 L342.4,140.8 L346.4,131.9 L350.4,132.2 L354.5,99.0 L358.5,84.8 L362.5,107.2 L366.6,101.9 L370.6,112.0 L374.6,113.1 L378.7,149.2 L382.7,154.5 L386.7,120.6 L390.8,121.6 L394.8,110.0 L398.8,104.6 L402.9,75.8 L406.9,61.6 L410.9,90.2 L415.0,90.5 L419.0,104.1 L423.0,117.9 L427.1,143.7 L431.1,166.7 L435.1,190.7 L439.2,191.0 L443.2,226.9 L447.2,210.8 L451.3,243.7 L455.3,242.8 L459.3,252.0 L463.4,238.5 L467.4,251.2 L471.4,250.7 L475.5,269.7 L479.5,255.5 L483.5,232.1 L487.6,262.6 L491.6,260.4 L495.6,267.3 L499.7,249.4 L503.7,234.5 L507.7,264.6 L511.8,268.2 L515.8,258.1 L519.8,244.8 L523.9,249.3 L527.9,259.0 L531.9,259.1 L536.0,265.6" fill="none" stroke="#eda100" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><text x="550" y="155.6" class="viz-label">Lax–Friedrichs</text><text x="550" y="126.0" class="viz-label">upwind</text><text x="550" y="87.2" class="viz-label">Lax–Wendroff</text><text x="550" y="65.6" class="viz-label">leapfrog</text><text x="201.2" y="108.5" class="viz-label">exact</text><line x1="0" x2="0" y1="30" y2="340" class="viz-crosshair" visibility="hidden"/></svg>
<div class="viz-tooltip" hidden></div>
</div>
<figcaption>u_t + u_x = 0 on a periodic grid with 120 points and Courant number ν = 0.5, after 240 steps (one full trip, so the exact solution equals the initial data). The labels on the right sit at each scheme's highest value near the box. Hover over the chart, or focus it and use the arrow keys, to read the values. The table shows every fourth grid point.</figcaption>
<details class="viz-table"><summary>Show the numbers</summary><table><thead><tr><th>x</th><th>Lax–Friedrichs</th><th>upwind</th><th>Lax–Wendroff</th><th>leapfrog</th><th>exact</th></tr></thead><tbody><tr><td>0.000</td><td>+0.068</td><td>+0.003</td><td>+0.000</td><td>+0.064</td><td>+0.000</td></tr><tr><td>0.033</td><td>+0.079</td><td>+0.010</td><td>-0.001</td><td>-0.025</td><td>+0.000</td></tr><tr><td>0.067</td><td>+0.110</td><td>+0.030</td><td>+0.006</td><td>+0.039</td><td>+0.000</td></tr><tr><td>0.100</td><td>+0.155</td><td>+0.078</td><td>-0.009</td><td>+0.077</td><td>+0.001</td></tr><tr><td>0.133</td><td>+0.209</td><td>+0.167</td><td>-0.051</td><td>-0.182</td><td>+0.017</td></tr><tr><td>0.167</td><td>+0.264</td><td>+0.295</td><td>+0.139</td><td>+0.155</td><td>+0.125</td></tr><tr><td>0.200</td><td>+0.312</td><td>+0.431</td><td>+0.632</td><td>+0.618</td><td>+0.472</td></tr><tr><td>0.233</td><td>+0.341</td><td>+0.522</td><td>+0.947</td><td>+0.955</td><td>+0.920</td></tr><tr><td>0.267</td><td>+0.349</td><td>+0.522</td><td>+0.792</td><td>+0.873</td><td>+0.920</td></tr><tr><td>0.300</td><td>+0.335</td><td>+0.431</td><td>+0.409</td><td>+0.285</td><td>+0.472</td></tr><tr><td>0.333</td><td>+0.306</td><td>+0.295</td><td>+0.168</td><td>+0.335</td><td>+0.125</td></tr><tr><td>0.367</td><td>+0.277</td><td>+0.169</td><td>+0.019</td><td>-0.025</td><td>+0.017</td></tr><tr><td>0.400</td><td>+0.262</td><td>+0.090</td><td>+0.003</td><td>-0.040</td><td>+0.001</td></tr><tr><td>0.433</td><td>+0.271</td><td>+0.071</td><td>+0.095</td><td>+0.166</td><td>+0.000</td></tr><tr><td>0.467</td><td>+0.312</td><td>+0.120</td><td>-0.185</td><td>-0.408</td><td>+0.000</td></tr><tr><td>0.500</td><td>+0.382</td><td>+0.241</td><td>-0.071</td><td>-0.013</td><td>+0.000</td></tr><tr><td>0.533</td><td>+0.472</td><td>+0.424</td><td>+0.471</td><td>+0.462</td><td>+0.000</td></tr><tr><td>0.567</td><td>+0.569</td><td>+0.626</td><td>+0.865</td><td>+0.980</td><td>+1.000</td></tr><tr><td>0.600</td><td>+0.654</td><td>+0.798</td><td>+0.950</td><td>+0.877</td><td>+1.000</td></tr><tr><td>0.633</td><td>+0.714</td><td>+0.907</td><td>+1.052</td><td>+1.052</td><td>+1.000</td></tr><tr><td>0.667</td><td>+0.737</td><td>+0.947</td><td>+0.912</td><td>+0.754</td><td>+1.000</td></tr><tr><td>0.700</td><td>+0.719</td><td>+0.923</td><td>+1.031</td><td>+1.033</td><td>+1.000</td></tr><tr><td>0.733</td><td>+0.661</td><td>+0.832</td><td>+1.218</td><td>+1.173</td><td>+1.000</td></tr><tr><td>0.767</td><td>+0.573</td><td>+0.674</td><td>+0.809</td><td>+0.793</td><td>+1.000</td></tr><tr><td>0.800</td><td>+0.467</td><td>+0.474</td><td>+0.295</td><td>+0.203</td><td>+0.000</td></tr><tr><td>0.833</td><td>+0.356</td><td>+0.281</td><td>+0.064</td><td>+0.025</td><td>+0.000</td></tr><tr><td>0.867</td><td>+0.255</td><td>+0.136</td><td>+0.009</td><td>-0.101</td><td>+0.000</td></tr><tr><td>0.900</td><td>+0.171</td><td>+0.053</td><td>+0.001</td><td>-0.035</td><td>+0.000</td></tr><tr><td>0.933</td><td>+0.112</td><td>+0.017</td><td>+0.000</td><td>-0.065</td><td>+0.000</td></tr><tr><td>0.967</td><td>+0.078</td><td>+0.005</td><td>+0.000</td><td>+0.044</td><td>+0.000</td></tr></tbody></table></details>
<script>(() => {const fig = document.currentScript.closest('figure');const svg = fig.querySelector('svg');const line = svg.querySelector('.viz-crosshair');const tip = fig.querySelector('.viz-tooltip');const names = ["Lax\u2013Friedrichs", "upwind", "Lax\u2013Wendroff", "leapfrog", "exact"]; const colors = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#86868b"];const rows = [["+0.068", "+0.003", "+0.000", "+0.064", "+0.000"], ["+0.073", "+0.004", "-0.000", "+0.053", "+0.000"], ["+0.071", "+0.005", "-0.001", "-0.063", "+0.000"], ["+0.078", "+0.007", "-0.001", "+0.065", "+0.000"], ["+0.079", "+0.010", "-0.001", "-0.025", "+0.000"], ["+0.088", "+0.013", "-0.000", "-0.104", "+0.000"], ["+0.092", "+0.017", "+0.001", "+0.053", "+0.000"], ["+0.102", "+0.023", "+0.004", "+0.022", "+0.000"], ["+0.110", "+0.030", "+0.006", "+0.039", "+0.000"], ["+0.121", "+0.039", "+0.007", "+0.072", "+0.000"], ["+0.131", "+0.049", "+0.006", "-0.110", "+0.000"], ["+0.143", "+0.062", "+0.000", "-0.096", "+0.001"], ["+0.155", "+0.078", "-0.009", "+0.077", "+0.001"], ["+0.168", "+0.096", "-0.023", "+0.109", "+0.002"], ["+0.181", "+0.117", "-0.037", "+0.026", "+0.005"], ["+0.195", "+0.140", "-0.049", "-0.124", "+0.009"], ["+0.209", "+0.167", "-0.051", "-0.182", "+0.017"], ["+0.223", "+0.196", "-0.038", "-0.092", "+0.030"], ["+0.237", "+0.227", "-0.003", "+0.042", "+0.050"], ["+0.251", "+0.260", "+0.055", "+0.180", "+0.080"], ["+0.264", "+0.295", "+0.139", "+0.155", "+0.125"], ["+0.277", "+0.330", "+0.246", "+0.095", "+0.185"], ["+0.290", "+0.365", "+0.369", "+0.307", "+0.264"], ["+0.301", "+0.399", "+0.501", "+0.530", "+0.360"], ["+0.312", "+0.431", "+0.632", "+0.618", "+0.472"], ["+0.321", "+0.461", "+0.751", "+0.825", "+0.594"], ["+0.329", "+0.486", "+0.848", "+0.955", "+0.717"], ["+0.335", "+0.506", "+0.915", "+0.890", "+0.829"], ["+0.341", "+0.522", "+0.947", "+0.955", "+0.920"], ["+0.345", "+0.531", "+0.947", "+0.997", "+0.979"], ["+0.348", "+0.534", "+0.917", "+0.801", "+1.000"], ["+0.348", "+0.531", "+0.864", "+0.754", "+0.979"], ["+0.349", "+0.522", "+0.792", "+0.873", "+0.920"], ["+0.345", "+0.506", "+0.706", "+0.809", "+0.829"], ["+0.344", "+0.486", "+0.610", "+0.698", "+0.717"], ["+0.338", "+0.461", "+0.508", "+0.569", "+0.594"], ["+0.335", "+0.431", "+0.409", "+0.285", "+0.472"], ["+0.325", "+0.399", "+0.322", "+0.125", "+0.360"], ["+0.322", "+0.365", "+0.252", "+0.159", "+0.264"], ["+0.310", "+0.330", "+0.203", "+0.229", "+0.185"], ["+0.306", "+0.295", "+0.168", "+0.335", "+0.125"], ["+0.293", "+0.261", "+0.138", "+0.284", "+0.080"], ["+0.291", "+0.228", "+0.105", "+0.099", "+0.050"], ["+0.277", "+0.197", "+0.064", "+0.030", "+0.030"], ["+0.277", "+0.169", "+0.019", "-0.025", "+0.017"], ["+0.262", "+0.144", "-0.019", "-0.107", "+0.009"], ["+0.267", "+0.122", "-0.037", "-0.076", "+0.005"], ["+0.252", "+0.104", "-0.030", "+0.012", "+0.002"], ["+0.262", "+0.090", "+0.003", "-0.040", "+0.001"], ["+0.248", "+0.079", "+0.049", "-0.069", "+0.001"], ["+0.263", "+0.072", "+0.090", "+0.170", "+0.000"], ["+0.250", "+0.069", "+0.110", "+0.275", "+0.000"], ["+0.271", "+0.071", "+0.095", "+0.166", "+0.000"], ["+0.260", "+0.076", "+0.045", "+0.117", "+0.000"], ["+0.287", "+0.086", "-0.031", "+0.013", "+0.000"], ["+0.278", "+0.100", "-0.115", "-0.203", "+0.000"], ["+0.312", "+0.120", "-0.185", "-0.408", "+0.000"], ["+0.304", "+0.143", "-0.223", "-0.308", "+0.000"], ["+0.344", "+0.171", "-0.216", "-0.139", "+0.000"], ["+0.337", "+0.204", "-0.163", "-0.135", "+0.000"], ["+0.382", "+0.241", "-0.071", "-0.013", "+0.000"], ["+0.377", "+0.282", "+0.050", "+0.083", "+0.000"], ["+0.425", "+0.327", "+0.188", "+0.182", "+0.000"], ["+0.421", "+0.374", "+0.330", "+0.189", "+0.000"], ["+0.472", "+0.424", "+0.471", "+0.462", "+0.000"], ["+0.469", "+0.475", "+0.601", "+0.576", "+0.000"], ["+0.521", "+0.526", "+0.714", "+0.867", "+1.000"], ["+0.519", "+0.577", "+0.803", "+0.865", "+1.000"], ["+0.569", "+0.626", "+0.865", "+0.980", "+1.000"], ["+0.567", "+0.674", "+0.901", "+0.817", "+1.000"], ["+0.614", "+0.719", "+0.920", "+0.912", "+1.000"], ["+0.613", "+0.760", "+0.933", "+0.814", "+1.000"], ["+0.654", "+0.798", "+0.950", "+0.877", "+1.000"], ["+0.654", "+0.832", "+0.978", "+0.874", "+1.000"], ["+0.688", "+0.861", "+1.011", "+1.110", "+1.000"], ["+0.688", "+0.886", "+1.040", "+1.211", "+1.000"], ["+0.714", "+0.907", "+1.052", "+1.052", "+1.000"], ["+0.714", "+0.923", "+1.040", "+1.090", "+1.000"], ["+0.731", "+0.935", "+1.004", "+1.018", "+1.000"], ["+0.730", "+0.943", "+0.955", "+1.010", "+1.000"], ["+0.737", "+0.947", "+0.912", "+0.754", "+1.000"], ["+0.737", "+0.947", "+0.892", "+0.716", "+1.000"], ["+0.733", "+0.943", "+0.906", "+0.957", "+1.000"], ["+0.733", "+0.935", "+0.956", "+0.950", "+1.000"], ["+0.719", "+0.923", "+1.031", "+1.033", "+1.000"], ["+0.719", "+0.907", "+1.114", "+1.070", "+1.000"], ["+0.695", "+0.886", "+1.184", "+1.275", "+1.000"], ["+0.695", "+0.861", "+1.222", "+1.376", "+1.000"], ["+0.661", "+0.832", "+1.218", "+1.173", "+1.000"], ["+0.661", "+0.798", "+1.167", "+1.171", "+1.000"], ["+0.621", "+0.760", "+1.074", "+1.074", "+1.000"], ["+0.621", "+0.719", "+0.950", "+0.976", "+1.000"], ["+0.573", "+0.674", "+0.809", "+0.793", "+1.000"], ["+0.573", "+0.626", "+0.663", "+0.630", "+1.000"], ["+0.521", "+0.577", "+0.524", "+0.460", "+1.000"], ["+0.522", "+0.526", "+0.400", "+0.457", "+1.000"], ["+0.467", "+0.474", "+0.295", "+0.203", "+0.000"], ["+0.467", "+0.423", "+0.211", "+0.317", "+0.000"], ["+0.411", "+0.373", "+0.146", "+0.084", "+0.000"], ["+0.411", "+0.326", "+0.098", "+0.090", "+0.000"], ["+0.356", "+0.281", "+0.064", "+0.025", "+0.000"], ["+0.356", "+0.239", "+0.040", "+0.120", "+0.000"], ["+0.304", "+0.201", "+0.025", "+0.030", "+0.000"], ["+0.304", "+0.166", "+0.015", "+0.034", "+0.000"], ["+0.255", "+0.136", "+0.009", "-0.101", "+0.000"], ["+0.255", "+0.110", "+0.005", "-0.000", "+0.000"], ["+0.210", "+0.088", "+0.003", "+0.166", "+0.000"], ["+0.211", "+0.069", "+0.001", "-0.051", "+0.000"], ["+0.171", "+0.053", "+0.001", "-0.035", "+0.000"], ["+0.172", "+0.041", "+0.000", "-0.084", "+0.000"], ["+0.138", "+0.031", "+0.000", "+0.043", "+0.000"], ["+0.140", "+0.023", "+0.000", "+0.149", "+0.000"], ["+0.112", "+0.017", "+0.000", "-0.065", "+0.000"], ["+0.114", "+0.012", "-0.000", "-0.090", "+0.000"], ["+0.092", "+0.009", "-0.000", "-0.019", "+0.000"], ["+0.095", "+0.006", "-0.000", "+0.075", "+0.000"], ["+0.078", "+0.005", "+0.000", "+0.044", "+0.000"], ["+0.081", "+0.004", "+0.000", "-0.025", "+0.000"], ["+0.070", "+0.003", "+0.000", "-0.026", "+0.000"], ["+0.075", "+0.003", "+0.000", "-0.072", "+0.000"]]; const heads = ["x = 0.000", "x = 0.008", "x = 0.017", "x = 0.025", "x = 0.033", "x = 0.042", "x = 0.050", "x = 0.058", "x = 0.067", "x = 0.075", "x = 0.083", "x = 0.092", "x = 0.100", "x = 0.108", "x = 0.117", "x = 0.125", "x = 0.133", "x = 0.142", "x = 0.150", "x = 0.158", "x = 0.167", "x = 0.175", "x = 0.183", "x = 0.192", "x = 0.200", "x = 0.208", "x = 0.217", "x = 0.225", "x = 0.233", "x = 0.242", "x = 0.250", "x = 0.258", "x = 0.267", "x = 0.275", "x = 0.283", "x = 0.292", "x = 0.300", "x = 0.308", "x = 0.317", "x = 0.325", "x = 0.333", "x = 0.342", "x = 0.350", "x = 0.358", "x = 0.367", "x = 0.375", "x = 0.383", "x = 0.392", "x = 0.400", "x = 0.408", "x = 0.417", "x = 0.425", "x = 0.433", "x = 0.442", "x = 0.450", "x = 0.458", "x = 0.467", "x = 0.475", "x = 0.483", "x = 0.492", "x = 0.500", "x = 0.508", "x = 0.517", "x = 0.525", "x = 0.533", "x = 0.542", "x = 0.550", "x = 0.558", "x = 0.567", "x = 0.575", "x = 0.583", "x = 0.592", "x = 0.600", "x = 0.608", "x = 0.617", "x = 0.625", "x = 0.633", "x = 0.642", "x = 0.650", "x = 0.658", "x = 0.667", "x = 0.675", "x = 0.683", "x = 0.692", "x = 0.700", "x = 0.708", "x = 0.717", "x = 0.725", "x = 0.733", "x = 0.742", "x = 0.750", "x = 0.758", "x = 0.767", "x = 0.775", "x = 0.783", "x = 0.792", "x = 0.800", "x = 0.808", "x = 0.817", "x = 0.825", "x = 0.833", "x = 0.842", "x = 0.850", "x = 0.858", "x = 0.867", "x = 0.875", "x = 0.883", "x = 0.892", "x = 0.900", "x = 0.908", "x = 0.917", "x = 0.925", "x = 0.933", "x = 0.942", "x = 0.950", "x = 0.958", "x = 0.967", "x = 0.975", "x = 0.983", "x = 0.992"];const xs = [56.0, 60.0, 64.1, 68.1, 72.1, 76.2, 80.2, 84.2, 88.3, 92.3, 96.3, 100.4, 104.4, 108.4, 112.5, 116.5, 120.5, 124.6, 128.6, 132.6, 136.7, 140.7, 144.7, 148.8, 152.8, 156.8, 160.9, 164.9, 168.9, 173.0, 177.0, 181.0, 185.1, 189.1, 193.1, 197.2, 201.2, 205.2, 209.3, 213.3, 217.3, 221.4, 225.4, 229.4, 233.5, 237.5, 241.5, 245.6, 249.6, 253.6, 257.7, 261.7, 265.7, 269.8, 273.8, 277.8, 281.9, 285.9, 289.9, 294.0, 298.0, 302.0, 306.1, 310.1, 314.1, 318.2, 322.2, 326.2, 330.3, 334.3, 338.3, 342.4, 346.4, 350.4, 354.5, 358.5, 362.5, 366.6, 370.6, 374.6, 378.7, 382.7, 386.7, 390.8, 394.8, 398.8, 402.9, 406.9, 410.9, 415.0, 419.0, 423.0, 427.1, 431.1, 435.1, 439.2, 443.2, 447.2, 451.3, 455.3, 459.3, 463.4, 467.4, 471.4, 475.5, 479.5, 483.5, 487.6, 491.6, 495.6, 499.7, 503.7, 507.7, 511.8, 515.8, 519.8, 523.9, 527.9, 531.9, 536.0]; const steps = 119, width = 680;let current = -1;const show = (n) => { n = Math.max(0, Math.min(steps, n)); current = n; const x = xs[n]; line.setAttribute('x1', x); line.setAttribute('x2', x); line.setAttribute('visibility', 'visible'); tip.replaceChildren(); const head = document.createElement('div'); head.className = 'viz-tooltip-head'; head.textContent = heads[n]; tip.append(head); names.forEach((name, i) => {  const row = document.createElement('div'); row.className = 'viz-tooltip-row';  const key = document.createElement('span'); key.className = 'viz-swatch'; key.style.background = colors[i];  const value = document.createElement('strong'); value.textContent = rows[n][i];  const label = document.createElement('span'); label.textContent = name;  row.append(key, value, label); tip.append(row); }); const box = svg.getBoundingClientRect(), frame = fig.getBoundingClientRect(); const xp = box.left - frame.left + x * box.width / width; tip.hidden = false; const room = frame.width - tip.offsetWidth - 8; tip.style.left = Math.max(0, Math.min(room, xp + 12)) + 'px'; tip.style.top = (box.top - frame.top + 8) + 'px';};const hide = () => { line.setAttribute('visibility', 'hidden'); tip.hidden = true; current = -1; };svg.addEventListener('pointermove', (event) => { const box = svg.getBoundingClientRect(); const x = (event.clientX - box.left) * width / box.width; let best = 0; xs.forEach((v, i) => { if (Math.abs(v - x) < Math.abs(xs[best] - x)) best = i; }); show(best);});svg.addEventListener('pointerleave', hide);svg.addEventListener('focus', () => show(current < 0 ? 0 : current));svg.addEventListener('blur', hide);svg.addEventListener('keydown', (event) => { if (event.key === 'ArrowRight') { show(current + 1); event.preventDefault(); } if (event.key === 'ArrowLeft') { show(current - 1); event.preventDefault(); } if (event.key === 'Escape') hide();});})();</script>
</figure>

| scheme | bump peak (exact 1) | largest value near the box | smallest value | largest error |
|---|---|---|---|---|
| Lax–Friedrichs | $0.349$ | $0.737$ | $0.068$ | $0.652$ |
| upwind | $0.534$ | $0.947$ | $0.003$ | $0.475$ |
| Lax–Wendroff | $0.947$ | $1.222$ | $-0.223$ | $0.601$ |
| leapfrog | $0.997$ | $1.376$ | $-0.408$ | $0.576$ |

- The first-order schemes smear, in the order of their $\epsilon$. Lax–Friedrichs keeps only $35\%$ of the bump, upwind $53\%$. They stay between $0$ and $1$: no new wiggles.
- The second-order schemes keep the bump almost at full height ($0.947$ and $0.997$). But near the jumps of the box they overshoot (up to $1.22$ and $1.38$) and leave oscillations behind, with negative values the exact solution never has. Leapfrog, with no damping at all, has the largest wiggles.
- So no scheme here wins on both counts. Low order smears, high order oscillates near jumps.
- All four schemes keep the total $h\sum_j U_j$ equal to the exact $0.3523$ to all printed digits. On a periodic grid, each of these updates only moves mass around.

To check the order, I ran one full trip from $\sin(2\pi x)$ with $\nu = 0.5$ and measured the largest error:

| points | 25 | 50 | 100 | 200 | 400 |
|---|---|---|---|---|---|
| Lax–Friedrichs | $6.97 \times 10^{-1}$ | $4.47 \times 10^{-1}$ | $2.56 \times 10^{-1}$ | $1.38 \times 10^{-1}$ | $7.14 \times 10^{-2}$ |
| observed order | | $0.64$ | $0.80$ | $0.90$ | $0.95$ |
| upwind | $3.26 \times 10^{-1}$ | $1.79 \times 10^{-1}$ | $9.40 \times 10^{-2}$ | $4.82 \times 10^{-2}$ | $2.44 \times 10^{-2}$ |
| observed order | | $0.87$ | $0.93$ | $0.97$ | $0.98$ |
| Lax–Wendroff | $4.93 \times 10^{-2}$ | $1.24 \times 10^{-2}$ | $3.10 \times 10^{-3}$ | $7.75 \times 10^{-4}$ | $1.94 \times 10^{-4}$ |
| observed order | | $2.00$ | $2.00$ | $2.00$ | $2.00$ |
| leapfrog | $4.98 \times 10^{-2}$ | $1.24 \times 10^{-2}$ | $3.10 \times 10^{-3}$ | $7.75 \times 10^{-4}$ | $1.94 \times 10^{-4}$ |
| observed order | | $2.00$ | $2.00$ | $2.00$ | $2.00$ |

Lax–Wendroff and leapfrog are second order from the start. Lax–Friedrichs and upwind approach order $1$. (The observed order is computed as in [Lecture 5](/blog/2026-10-20-numerical-pde-lecture-5/).)

Lax–Friedrichs is about three times worse than upwind ($0.0714/0.0244 = 2.9$), although its $\epsilon$ is only twice as large at $\nu = 0.5$. My own check of this: for a smooth mode, one step damps by about $1 - \frac12(1 - \nu^2)\theta^2$ for Lax–Friedrichs and $1 - \frac12\nu(1 - \nu)\theta^2$ for upwind. The ratio of the damping is $(1 + \nu)/\nu = 3$. The reason is that the Forward Euler step itself acts like a negative diffusion of size $a^2k/2$, so the net smoothing is $\epsilon - a^2k/2$. This also explains the lower end $c = \nu^2$ in section 8: there the net diffusion is zero at leading order.

**Part 4: two cases from sections 5 and 10.** First, FTCS from the bump-and-box data for one full trip, with a fixed $\nu = 0.5$:

| points | 20 | 40 | 80 |
|---|---|---|---|
| steps | 40 | 80 | 160 |
| largest $\lvert U\rvert$ at $t = 1$ | $21.4$ | $1.09 \times 10^{3}$ | $4.49 \times 10^{6}$ |
| $\big(\max_p\lvert 1 + k\lambda_p\rvert\big)^{\text{steps}}$ | $86.7$ | $7.52 \times 10^{3}$ | $5.66 \times 10^{7}$ |

Refining the grid makes it worse: more steps, each with growth up to $1.118$. Now FTCS with $k = h^2$ from $\sin(2\pi x)$ to $t = 1$:

| points | 20 | 40 | 80 |
|---|---|---|---|
| steps | 400 | 1600 | 6400 |
| largest $\lvert U\rvert$ at $t = 1$ | $1.043$ | $1.012$ | $1.003$ |
| $\big(\max_p\lvert 1 + k\lambda_p\rvert\big)^{\text{steps}}$ | $1.6477$ | $1.6485$ | $1.6487$ |

The worst growth stays below the bound $e^{a^2T/2} = e^{0.5} = 1.6487$ and gets close to it. The solution itself tends to the exact maximum $1$. So FTCS does converge this way, but $80$ points need $6400$ steps.

Second, the upwind formula $(U_j - U_{j-1})/h$ used with $a = -1$ ($\nu = -0.5$), from a smooth bump on $40$ points:

| step $n$ | 20 | 40 | 60 | 80 |
|---|---|---|---|---|
| largest $\lvert U\rvert$ | $5.65$ | $6.25 \times 10^{3}$ | $1.09 \times 10^{8}$ | $7.56 \times 10^{12}$ |

The sawtooth mode has growth factor $1 + 2|\nu| = 2$ per step. Negative diffusion is a disaster.

The full script and the chart code are in [this site's repository](https://github.com/doubleyands/doubleyands.github.io/tree/main/scripts/posts/numerical-pde/lecture-12).

## Summary

1. For $u_t + au_x = 0$ the solution is constant on the **characteristics** $x - at = x_0$. To find $u(X, T)$, go back to the foot $X - aT$. Give a boundary value only at the **inflow** boundary.
2. With periodic boundaries, the central-difference matrix is circulant and **skew-symmetric**, so its eigenvalues $-\frac{ia}{h}\sin(2\pi ph)$ are purely imaginary.
3. The Forward Euler disk touches the imaginary axis only at $0$, so FTCS is unstable for every fixed $k/h$. Leapfrog is stable for $|ak/h| < 1$, but only marginally: nothing is damped.
4. **Lax–Friedrichs = FTCS + artificial diffusion** with $\epsilon = h^2/(2k)$. The diffusion moves the eigenvalues left, onto an ellipse. The ellipse is inside the disk exactly when $\nu^2 \le c \le 1$, or $a^2k/2 \le \epsilon \le h^2/(2k)$.
5. **Lax–Wendroff** comes from a Taylor expansion with $u_{tt} = a^2u_{xx}$. It is second order and has the smallest stable diffusion, $\epsilon = a^2k/2$. **Upwind** has $\epsilon = ah/2$. Used on the wrong side, it gives negative diffusion and blows up.
6. In practice, the first-order schemes smear sharp features, and the second-order schemes keep them sharp but oscillate near jumps.
