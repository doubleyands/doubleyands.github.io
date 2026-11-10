---
title: "Numerical PDEs, Lecture 11: von Neumann Analysis, the Amplification Factor, and a First Look at Advection"
description: "My notes on the eleventh lecture: Fourier modes as eigenvectors of difference operators, the amplification factor g(ξ), von Neumann stability of FTCS and Crank–Nicolson for the heat equation, the condition number of 2D Crank–Nicolson, and the advection equation with the CFL condition."
pubDate: 2026-11-10
project: numerical-pde
tags: [numerical PDEs, heat equation, advection, stability, Fourier analysis, notes]
---

This is my summary of the eleventh lecture in my numerical PDE course. I write the ideas in my own words and check them with a small experiment.

[Lecture 10](/blog/2026-11-05-numerical-pde-lecture-10/) found the stability limits of the heat equation schemes with the **eigenvalues** of the space matrix. That worked because the eigenvalues were known exactly. For most real problems they are not. At the end of Lecture 10 we started a second tool, **von Neumann analysis**, which uses Fourier modes instead of eigenvectors. This lecture finishes it. It gives the same answers as before (FTCS needs $k/h^2 \le 1/2$, Crank–Nicolson is always stable) without any eigenvalue computation. Then it looks at 2D Crank–Nicolson, and it opens the next topic: the **advection equation** and the **CFL condition**.

Notation as before: $h$ is the grid spacing, $k$ is the time step, $x_j = jh$, and $U^n_j$ is the value at $x_j$ and time $t_n = nk$. In exponents, $i = \sqrt{-1}$ and $j$ is the grid index.

## 1. Reading a scheme from its stencil

A **stencil** is a picture of which grid values one step of a scheme uses.

- **FTCS**: three values at time $n$ ($j-1$, $j$, $j+1$) give one value at time $n+1$. Each new value comes out directly. So the scheme is explicit.
- **Crank–Nicolson**: it averages times $n$ and $n+1$, so three values at time $n+1$ are tied together too. The new values depend on each other, and we must solve a linear system. So the scheme is implicit.

The stencil tells you at once whether a scheme is explicit or implicit, and how far information moves in one step.

Both schemes come from the **method of lines**: discretize only in space and get the ODE system $U'(t) = AU(t)$, then choose a time stepper. We used this in [Lecture 10](/blog/2026-11-05-numerical-pde-lecture-10/). The name comes from the picture: at each grid point $x_j$ there is one line in the time direction, and $U_j(t)$ moves along it.

## 2. Fourier modes are eigenfunctions of the derivative

The key fact behind everything in this lecture is short:

$$
\frac{\partial}{\partial x}\,e^{i\xi x} = i\xi\,e^{i\xi x}.
$$

This has the form $Au = \lambda u$. So $e^{i\xi x}$ is an **eigenfunction** of $\partial_x$, with eigenvalue $i\xi$. The same is true for every linear differential operator with **constant coefficients**: $\partial_x^2 e^{i\xi x} = -\xi^2 e^{i\xi x}$, and so on. With Euler's formula $e^{i\xi x} = \cos\xi x + i\sin\xi x$, sine and cosine series are special cases.

A function with $\int |v|^2\,dx < \infty$ can be written as a mix of these modes with the Fourier transform, and **Parseval's relation** says the transform keeps the 2-norm: $\|\hat v\|_2 = \|v\|_2$. To expand in these modes on the whole domain, we need either no boundary at all or periodic boundaries. These are the assumptions from the end of Lecture 10.

For the heat equation on $[0, 1]$ with zero ends, the exact solution is a sine series in which each coefficient decays on its own. This is beautiful, but it only works in very special cases. The Fourier idea becomes really useful when we apply it to the **scheme**, not to the PDE.

## 3. The grid version: $e^{ijh\xi}$ and the central difference

On the grid, the mode becomes the **grid function** $W_j = e^{ijh\xi}$. A finite difference scheme knows nothing between the grid points. It only relates grid values to each other. So let us apply the central difference $D_0V_j = (V_{j+1} - V_{j-1})/(2h)$ to $W$. Take out the common factor $e^{ijh\xi}$ first:

$$
D_0W_j = \frac{1}{2h}\big(e^{ih\xi} - e^{-ih\xi}\big)e^{ijh\xi} = \frac{i}{h}\sin(h\xi)\,W_j.
$$

In the bracket the cosines cancel (they are even), and only the sine is left. So $W$ is an eigenvector of $D_0$, and the eigenvalue is $\frac{i}{h}\sin(h\xi)$. A general grid function $V$ does not come back as a multiple of itself, but this special one does. In fact $e^{ijh\xi}$ is an eigenvector of **every** difference operator that looks the same at every grid point. The second difference $(V_{j-1} - 2V_j + V_{j+1})/h^2$ is one example, and we need it in section 5.

**Compare with the exact eigenvalue.** Expand the sine:

$$
\frac{i}{h}\sin(h\xi) = i\xi - \frac{i}{6}h^2\xi^3 + O(h^4\xi^5).
$$

It agrees with $i\xi$ up to an $O(h^2)$ term. This is the second-order accuracy of the central difference, seen through the eigenvalue.

**Resolution.** The error term has $\xi^3$ in it. So "second order" alone is not enough. A wave with a large wave number (many wiggles) needs a small $h$. The relative error of the eigenvalue is $1 - \sin(\xi h)/(\xi h) \approx (\xi h)^2/6$. It depends only on $\xi h$, or on the number of grid points per wavelength, $M = 2\pi/(\xi h)$. Part 1 of the experiment gives the numbers.

There is also a hard limit. The grid can only tell apart wave numbers with $|\xi h| \le \pi$, which means more than two grid points per wavelength. A faster wave looks exactly like a slower one at the grid points. For example, on the grid $h = 1/10$, the function $\sin(24\pi x)$ has the same grid values as $\sin(4\pi x)$, because $12j/10$ and $2j/10$ differ by the whole number $j$. This is called **aliasing**. The fastest wave the grid can show is the sawtooth, $\xi h = \pi$.

**Grid norm.** For grid functions we use $\|V\|_2 = \big(h\sum_j |V_j|^2\big)^{1/2}$. The factor $h$ makes the sum imitate an integral. With it, the discrete Fourier transform from [Lecture 10](/blog/2026-11-05-numerical-pde-lecture-10/) satisfies Parseval's relation $\|\hat V\|_2 = \|V\|_2$.

## 4. von Neumann analysis

**The goal.** In Lecture 10, stability meant $\|B^n\| \le M$ for $nk \le T$. A simple way to get this is to show, for every $U^n$,

$$
\|U^{n+1}\|_2 \le (1 + \alpha k)\,\|U^n\|_2
$$

with a constant $\alpha$. Then $\|U^n\|_2 \le (1 + \alpha k)^n\|U^0\|_2 \le e^{\alpha T}\|U^0\|_2$, the same trick as in [Lecture 8](/blog/2026-10-29-numerical-pde-lecture-8/).

Showing this directly is hard: there are many grid values, and each new value mixes its neighbors. Working with the matrix $B$ itself is possible, but $B$ gets bigger every time we refine the grid. So we use Parseval's relation and work on the **Fourier side**, where it is enough to show $\|\hat U^{n+1}\|_2 \le (1 + \alpha k)\|\hat U^n\|_2$.

**On the Fourier side, wave numbers do not mix.** Because $e^{ijh\xi}$ is an eigenvector of the scheme, a one-step scheme acts on each wave number separately:

$$
\hat U^{n+1}(\xi) = g(\xi)\,\hat U^n(\xi).
$$

The number $g(\xi)$ is the **amplification factor** at wave number $\xi$: how much one step multiplies that mode. Then

$$
|g(\xi)| \le 1 + \alpha k \ \text{ for all } \xi \ (\alpha \text{ not depending on } \xi)
\quad\Longrightarrow\quad
\|\hat U^{n+1}\|_2 \le (1 + \alpha k)\|\hat U^n\|_2
\quad\Longrightarrow\quad \text{stable}.
$$

The last arrow is Parseval. In practice we usually check the simpler $|g(\xi)| \le 1$.

**The recipe.**

1. Put $U^n_j = e^{ijh\xi}$ into the scheme. (Strictly, it should be $U^n_j = \hat U^n e^{ijh\xi}$, but the factor cancels.)
2. Write the result as $U^{n+1}_j = g(\xi)\,e^{ijh\xi}$. Take out the common factor $e^{ijh\xi}$ and use Euler's formula.
3. Find the condition on $k$ and $h$ that gives $|g(\xi)| \le 1$ for **every** $\xi$.

This small computation decides the step size limit, and it works for any linear scheme with constant coefficients.

### FTCS for the heat equation

The scheme is

$$
U^{n+1}_j = U^n_j + \frac{k}{h^2}\big(U^n_{j-1} - 2U^n_j + U^n_{j+1}\big).
$$

Put in the mode:

$$
g(\xi)\,e^{ijh\xi} = \Big(1 + \frac{k}{h^2}\big(e^{-i\xi h} - 2 + e^{i\xi h}\big)\Big)e^{ijh\xi}
\quad\Longrightarrow\quad
g(\xi) = 1 + \frac{2k}{h^2}\big(\cos(\xi h) - 1\big).
$$

Since $-1 \le \cos(\xi h) \le 1$, we get $1 - 4k/h^2 \le g(\xi) \le 1$. The top is always fine. The bottom must be at least $-1$:

$$
1 - \frac{4k}{h^2} \ge -1 \quad\Longleftrightarrow\quad \boxed{\;\frac{k}{h^2} \le \frac12.\;}
$$

This is the same limit as in Lecture 10, but we did not need any eigenvalues. The worst mode is $\cos(\xi h) = -1$, the sawtooth. With $\cos\theta - 1 = -2\sin^2(\theta/2)$ the same factor is often written $g = 1 - 4\frac{k}{h^2}\sin^2\frac{\xi h}{2}$.

### Crank–Nicolson for the heat equation

Now $U^{n+1}_j = g\,e^{ijh\xi}$ also appears on the right side:

$$
g = 1 + \frac{k}{2h^2}\big(e^{-i\xi h} - 2 + e^{i\xi h}\big)(1 + g)
\quad\Longrightarrow\quad
g(\xi) = \frac{1 + \frac{k}{h^2}\big(\cos(\xi h) - 1\big)}{1 - \frac{k}{h^2}\big(\cos(\xi h) - 1\big)}.
$$

Write $c = \frac{k}{h^2}(1 - \cos\xi h) \ge 0$. Then $g = (1 - c)/(1 + c)$, and $|1 - c| \le 1 + c$ always. So $|g| \le 1$ for **every** $k$ and $h$: Crank–Nicolson is unconditionally stable. For large $k/h^2$, the sawtooth has $g$ close to $-1$, the same effect we saw in Part 4 of Lecture 10.

### Backward Euler (my own addition)

The lecture treated Backward Euler only with eigenvalues. The recipe gives

$$
g(\xi) = \frac{1}{1 + \frac{2k}{h^2}\big(1 - \cos(\xi h)\big)},
$$

which is between $0$ and $1$ for every $k$. So it is also unconditionally stable, and unlike Crank–Nicolson it damps the sawtooth strongly.

**Why the answers match Lecture 10.** Put $\xi h = p\pi h$. Then the FTCS term $\frac{2k}{h^2}(\cos\xi h - 1)$ is exactly $k\lambda_p$, the scaled eigenvalue from Lecture 10. So for this problem, von Neumann analysis does the eigenvalue computation and the stability region check in one step. Its advantage is that it does not need the eigenvectors of a matrix with boundaries.

Part 2 of the experiment measures $g$ by doing one step of each scheme from a single mode.

## 5. Crank–Nicolson in 2D and the condition number

For $u_t = u_{xx} + u_{yy}$, use the 5-point Laplacian $\nabla_h^2$ ([Lecture 3](/blog/2026-10-13-numerical-pde-lecture-3/)) in space and the trapezoid in time:

$$
U^{n+1}_{ij} = U^n_{ij} + \frac{k}{2}\big(\nabla_h^2U^n_{ij} + \nabla_h^2U^{n+1}_{ij}\big)
\quad\Longrightarrow\quad
\Big(I - \frac{k}{2}\nabla_h^2\Big)U^{n+1} = \Big(I + \frac{k}{2}\nabla_h^2\Big)U^n.
$$

It is implicit, so every step solves a linear system. The matrix $A = I - \frac{k}{2}\nabla_h^2$ has the same sparse pattern as $\nabla_h^2$. (This $A$ is a new matrix, not the 1D matrix of Lecture 10.) With zero boundary values on the unit square, its eigenvalues are

$$
\lambda_{p,q} = 1 - \frac{k}{h^2}\Big[\big(\cos(p\pi h) - 1\big) + \big(\cos(q\pi h) - 1\big)\Big], \qquad p, q = 1, \dots, N.
$$

- The largest is $1 + O(k/h^2)$, and the smallest is $1 + O(k)$.
- So $\kappa(A) = O(k/h^2)$. For the Laplacian alone it was $O(1/h^2)$ ([Lecture 10](/blog/2026-11-05-numerical-pde-lecture-10/)). The identity $I$ makes the diagonal larger, and the time step **reduces** the condition number.

This helps an iterative solver in two ways. A smaller condition number means faster convergence. And there is a very good starting guess: $U^{n+1} = U^n + O(k)$, so we can start from $U^n$, or from the extrapolation $2U^n - U^{n-1}$. A direct solver like Gaussian elimination cannot use a starting guess. For time-dependent problems with thousands of steps, this makes iterative solvers attractive.

## 6. Three tools for stability

| tool | when to use it |
|---|---|
| eigenvalues + stability region (Lecture 10) | the eigenvalues are known (uniform grid, constant coefficients, simple boundary) |
| matrix norm $\lVert B^n\rVert \le M$ (Lecture 10) | the link to convergence: the Lax equivalence theorem |
| **von Neumann analysis** (this lecture) | linear, constant coefficients, no boundary or periodic boundary. This is the one we compute with in practice. |

Designing a scheme is easy. Checking its stability is the real work, and von Neumann analysis is the basic tool for it.

## 7. A first look at the advection equation

The heat equation only smooths things out. The **advection equation** moves them:

$$
u_t + au_x = 0, \qquad u(x, 0) = \eta(x) \quad\Longrightarrow\quad u(x, t) = \eta(x - at).
$$

The initial shape moves with speed $a$, to the right if $a > 0$ and to the left if $a < 0$, and it does not change. The check uses only the chain rule: at $t = 0$ we get $\eta(x)$, and $u_x = \eta'(x - at)$, $u_t = -a\,\eta'(x - at) = -a\,u_x$. The math is short. The meaning (the lines along which the solution is constant, called characteristics) is the topic of the next lecture.

**Two schemes.** The first one copies FTCS from the heat equation: forward in time, central in space.

$$
U^{n+1}_j = U^n_j - \frac{ak}{2h}\big(U^n_{j+1} - U^n_{j-1}\big).
$$

It looks reasonable, but it is **unstable** and we cannot use it. The **Lax–Friedrichs** scheme changes only one thing: it replaces $U^n_j$ by the average of its two neighbors.

$$
U^{n+1}_j = \frac{U^n_{j-1} + U^n_{j+1}}{2} - \frac{ak}{2h}\big(U^n_{j+1} - U^n_{j-1}\big).
$$

This small change makes it work. Lax–Friedrichs is stable when

$$
\left|\frac{ak}{h}\right| \le 1,
$$

the **CFL condition** (Courant–Friedrichs–Lewy). The number $\nu = ak/h$ is called the Courant number. The proofs of these two claims come next lecture. Part 3 of the experiment checks them with numbers.

**An important difference from heat.** CFL allows $k = O(h)$, even for an explicit method. For the heat equation, FTCS needed $k = O(h^2)$. The reason is the order of the derivative. Advection has a first derivative, so the difference formulas have $1/h$. Diffusion has a second derivative, so they have $1/h^2$. So advection problems are **not stiff**. This is why hyperbolic problems like advection are usually solved with explicit methods, and parabolic problems like heat with implicit methods.

The step size $k \approx h/a$ also makes sense on its own. Since $u_t = -a\,u_x$, the solution changes in time $a$ times faster than in space. To resolve time as well as space is resolved by $h$, we need $k$ about $a$ times smaller than $h$.

**Boundaries.** If $a > 0$, the solution moves to the right. So on $0 \le x \le 1$ we give a boundary value only at the **inflow** boundary, $u(0, t) = g_0(t)$. The value at $x = 1$ comes from inside.

## 8. Experiment

The script uses the Python standard library only (`cmath` for complex numbers, and a small dense solver for the periodic implicit schemes). In this post, $r = k/h^2$ for the heat equation and $\nu = ak/h$ for advection.

**Part 1: the eigenvalue of $D_0$.** On a periodic grid with $40$ points, I put $W_j = e^{ijh\xi}$ into $D_0$ for three wave numbers and compared $D_0W$ with $\frac{i}{h}\sin(h\xi)W$ at every point. The largest relative difference is $7.1 \times 10^{-15}$, which is rounding error. Then the ratio of the discrete eigenvalue to the exact one, $\sin(\xi h)/(\xi h)$, for $M$ grid points per wavelength:

| $M$ | 64 | 32 | 16 | 8 | 6 | 4 | 3 | 2 |
|---|---|---|---|---|---|---|---|---|
| $\xi h$ | $0.098$ | $0.196$ | $0.393$ | $0.785$ | $1.047$ | $1.571$ | $2.094$ | $\pi$ |
| ratio | $0.9984$ | $0.9936$ | $0.9745$ | $0.9003$ | $0.8270$ | $0.6366$ | $0.4135$ | $0$ |
| relative error | $1.61 \times 10^{-3}$ | $6.41 \times 10^{-3}$ | $2.55\%$ | $9.97\%$ | $17.3\%$ | $36.3\%$ | $58.7\%$ | $100\%$ |
| $(\xi h)^2/6$ | $1.61 \times 10^{-3}$ | $6.43 \times 10^{-3}$ | $2.57\%$ | $10.3\%$ | $18.3\%$ | $41.1\%$ | $73.1\%$ | $164\%$ |

- Halving $\xi h$ cuts the error by about $4$: second order.
- For fine resolution, the Taylor estimate $(\xi h)^2/6$ is almost exact. For coarse resolution it is too pessimistic, but the real error is large anyway.
- At $M = 2$ (the sawtooth), $D_0$ gives exactly zero: the central difference cannot see this wave at all. To get the slope within about $2.5\%$, we need about $16$ points per wavelength.

**Part 2: the amplification factor, measured.** On a periodic grid with $32$ points, I start from one mode $U_j = e^{ijh\xi}$ with $\xi h = 2\pi m/32$, $m = 0, \dots, 16$, and do **one** step of the scheme. Then I divide $U^{n+1}_j$ by $U^n_j$ at every grid point. If the theory is right, the ratio is the same number at every point, and that number is $g(\xi)$.

- Over all schemes and modes, the ratio differs between grid points by at most $1.5 \times 10^{-14}$. So one step really just multiplies the mode by a number.
- The measured $g$ agrees with the formulas of section 4 to within $5.9 \times 10^{-15}$.

<figure class="viz">
<div class="viz-legend"><span class="viz-key"><span class="viz-swatch" style="background:#2a78d6"></span>FTCS, r = 0.3</span><span class="viz-key"><span class="viz-swatch" style="background:#eb6834"></span>FTCS, r = 0.5</span><span class="viz-key"><span class="viz-swatch" style="background:#1baf7a"></span>FTCS, r = 0.6</span><span class="viz-key"><span class="viz-swatch" style="background:#eda100"></span>Crank–Nicolson, r = 6</span></div>
<div class="viz-frame">
<svg viewBox="0 0 680 400" role="img" tabindex="0" aria-label="Measured amplification factor g against xi h from 0 to pi, for FTCS with r = 0.3, 0.5 and 0.6 and for Crank-Nicolson with r = 6. All curves start at 1. FTCS with r = 0.6 goes below minus 1 near the sawtooth and ends at minus 1.4. FTCS with r = 0.5 ends exactly at minus 1. Crank-Nicolson stays between minus 1 and 1." class="viz-plot"><line x1="64" x2="500" y1="52.1" y2="52.1" class="viz-leader" stroke-dasharray="5 4"/><text x="54" y="56.1" text-anchor="end" class="viz-tick">1</text><line x1="64" x2="500" y1="107.5" y2="107.5" class="viz-grid"/><text x="54" y="111.5" text-anchor="end" class="viz-tick">0.5</text><line x1="64" x2="500" y1="162.9" y2="162.9" class="viz-grid"/><text x="54" y="166.9" text-anchor="end" class="viz-tick">0</text><line x1="64" x2="500" y1="218.2" y2="218.2" class="viz-grid"/><text x="54" y="222.2" text-anchor="end" class="viz-tick">−0.5</text><line x1="64" x2="500" y1="273.6" y2="273.6" class="viz-leader" stroke-dasharray="5 4"/><text x="54" y="277.6" text-anchor="end" class="viz-tick">−1</text><line x1="64" x2="500" y1="328.9" y2="328.9" class="viz-grid"/><text x="54" y="332.9" text-anchor="end" class="viz-tick">−1.5</text><text x="64.0" y="362" text-anchor="middle" class="viz-tick">0</text><text x="173.0" y="362" text-anchor="middle" class="viz-tick">π/4</text><text x="282.0" y="362" text-anchor="middle" class="viz-tick">π/2</text><text x="391.0" y="362" text-anchor="middle" class="viz-tick">3π/4</text><text x="500.0" y="362" text-anchor="middle" class="viz-tick">π</text><text x="282.0" y="388" text-anchor="middle" class="viz-axis-title">Wave number times grid spacing, ξh (π is the sawtooth)</text><text x="64" y="14" class="viz-axis-title">Amplification factor g per step (dashed lines: ±1)</text><path d="M64.0,52.1 L91.2,53.4 L118.5,57.2 L145.8,63.3 L173.0,71.6 L200.2,81.7 L227.5,93.2 L254.8,105.6 L282.0,118.6 L309.2,131.5 L336.5,144.0 L363.7,155.5 L391.0,165.5 L418.3,173.8 L445.5,179.9 L472.8,183.7 L500.0,185.0" fill="none" stroke="#2a78d6" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><circle cx="64.0" cy="52.1" r="2.5" fill="#2a78d6"/><circle cx="91.2" cy="53.4" r="2.5" fill="#2a78d6"/><circle cx="118.5" cy="57.2" r="2.5" fill="#2a78d6"/><circle cx="145.8" cy="63.3" r="2.5" fill="#2a78d6"/><circle cx="173.0" cy="71.6" r="2.5" fill="#2a78d6"/><circle cx="200.2" cy="81.7" r="2.5" fill="#2a78d6"/><circle cx="227.5" cy="93.2" r="2.5" fill="#2a78d6"/><circle cx="254.8" cy="105.6" r="2.5" fill="#2a78d6"/><circle cx="282.0" cy="118.6" r="2.5" fill="#2a78d6"/><circle cx="309.2" cy="131.5" r="2.5" fill="#2a78d6"/><circle cx="336.5" cy="144.0" r="2.5" fill="#2a78d6"/><circle cx="363.7" cy="155.5" r="2.5" fill="#2a78d6"/><circle cx="391.0" cy="165.5" r="2.5" fill="#2a78d6"/><circle cx="418.3" cy="173.8" r="2.5" fill="#2a78d6"/><circle cx="445.5" cy="179.9" r="2.5" fill="#2a78d6"/><circle cx="472.8" cy="183.7" r="2.5" fill="#2a78d6"/><circle cx="500.0" cy="185.0" r="2.5" fill="#2a78d6"/><path d="M64.0,52.1 L91.2,54.3 L118.5,60.6 L145.8,70.8 L173.0,84.6 L200.2,101.3 L227.5,120.5 L254.8,141.3 L282.0,162.9 L309.2,184.5 L336.5,205.2 L363.7,224.4 L391.0,241.1 L418.3,254.9 L445.5,265.1 L472.8,271.4 L500.0,273.6" fill="none" stroke="#eb6834" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><circle cx="64.0" cy="52.1" r="2.5" fill="#eb6834"/><circle cx="91.2" cy="54.3" r="2.5" fill="#eb6834"/><circle cx="118.5" cy="60.6" r="2.5" fill="#eb6834"/><circle cx="145.8" cy="70.8" r="2.5" fill="#eb6834"/><circle cx="173.0" cy="84.6" r="2.5" fill="#eb6834"/><circle cx="200.2" cy="101.3" r="2.5" fill="#eb6834"/><circle cx="227.5" cy="120.5" r="2.5" fill="#eb6834"/><circle cx="254.8" cy="141.3" r="2.5" fill="#eb6834"/><circle cx="282.0" cy="162.9" r="2.5" fill="#eb6834"/><circle cx="309.2" cy="184.5" r="2.5" fill="#eb6834"/><circle cx="336.5" cy="205.2" r="2.5" fill="#eb6834"/><circle cx="363.7" cy="224.4" r="2.5" fill="#eb6834"/><circle cx="391.0" cy="241.1" r="2.5" fill="#eb6834"/><circle cx="418.3" cy="254.9" r="2.5" fill="#eb6834"/><circle cx="445.5" cy="265.1" r="2.5" fill="#eb6834"/><circle cx="472.8" cy="271.4" r="2.5" fill="#eb6834"/><circle cx="500.0" cy="273.6" r="2.5" fill="#eb6834"/><path d="M64.0,52.1 L91.2,54.7 L118.5,62.3 L145.8,74.5 L173.0,91.1 L200.2,111.2 L227.5,134.2 L254.8,159.1 L282.0,185.0 L309.2,210.9 L336.5,235.8 L363.7,258.8 L391.0,278.9 L418.3,295.5 L445.5,307.7 L472.8,315.3 L500.0,317.9" fill="none" stroke="#1baf7a" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><circle cx="64.0" cy="52.1" r="2.5" fill="#1baf7a"/><circle cx="91.2" cy="54.7" r="2.5" fill="#1baf7a"/><circle cx="118.5" cy="62.3" r="2.5" fill="#1baf7a"/><circle cx="145.8" cy="74.5" r="2.5" fill="#1baf7a"/><circle cx="173.0" cy="91.1" r="2.5" fill="#1baf7a"/><circle cx="200.2" cy="111.2" r="2.5" fill="#1baf7a"/><circle cx="227.5" cy="134.2" r="2.5" fill="#1baf7a"/><circle cx="254.8" cy="159.1" r="2.5" fill="#1baf7a"/><circle cx="282.0" cy="185.0" r="2.5" fill="#1baf7a"/><circle cx="309.2" cy="210.9" r="2.5" fill="#1baf7a"/><circle cx="336.5" cy="235.8" r="2.5" fill="#1baf7a"/><circle cx="363.7" cy="258.8" r="2.5" fill="#1baf7a"/><circle cx="391.0" cy="278.9" r="2.5" fill="#1baf7a"/><circle cx="418.3" cy="295.5" r="2.5" fill="#1baf7a"/><circle cx="445.5" cy="307.7" r="2.5" fill="#1baf7a"/><circle cx="472.8" cy="315.3" r="2.5" fill="#1baf7a"/><circle cx="500.0" cy="317.9" r="2.5" fill="#1baf7a"/><path d="M64.0,52.1 L91.2,75.0 L118.5,121.6 L145.8,163.5 L173.0,193.3 L200.2,213.2 L227.5,226.5 L254.8,235.6 L282.0,241.9 L309.2,246.5 L336.5,249.8 L363.7,252.1 L391.0,253.9 L418.3,255.1 L445.5,255.9 L472.8,256.4 L500.0,256.5" fill="none" stroke="#eda100" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><circle cx="64.0" cy="52.1" r="2.5" fill="#eda100"/><circle cx="91.2" cy="75.0" r="2.5" fill="#eda100"/><circle cx="118.5" cy="121.6" r="2.5" fill="#eda100"/><circle cx="145.8" cy="163.5" r="2.5" fill="#eda100"/><circle cx="173.0" cy="193.3" r="2.5" fill="#eda100"/><circle cx="200.2" cy="213.2" r="2.5" fill="#eda100"/><circle cx="227.5" cy="226.5" r="2.5" fill="#eda100"/><circle cx="254.8" cy="235.6" r="2.5" fill="#eda100"/><circle cx="282.0" cy="241.9" r="2.5" fill="#eda100"/><circle cx="309.2" cy="246.5" r="2.5" fill="#eda100"/><circle cx="336.5" cy="249.8" r="2.5" fill="#eda100"/><circle cx="363.7" cy="252.1" r="2.5" fill="#eda100"/><circle cx="391.0" cy="253.9" r="2.5" fill="#eda100"/><circle cx="418.3" cy="255.1" r="2.5" fill="#eda100"/><circle cx="445.5" cy="255.9" r="2.5" fill="#eda100"/><circle cx="472.8" cy="256.4" r="2.5" fill="#eda100"/><circle cx="500.0" cy="256.5" r="2.5" fill="#eda100"/><text x="510" y="189.0" class="viz-label">FTCS, r = 0.3</text><text x="510" y="283.6" class="viz-label">FTCS, r = 0.5</text><text x="510" y="321.9" class="viz-label">FTCS, r = 0.6</text><text x="510" y="250.5" class="viz-label">Crank–Nicolson, r = 6</text><line x1="0" x2="0" y1="30" y2="340" class="viz-crosshair" visibility="hidden"/></svg>
<div class="viz-tooltip" hidden></div>
</div>
<figcaption>The factor g by which one step multiplies the mode e^(ijhξ), for u_t = u_xx on a periodic grid with 32 points. Each dot is measured from one step of the scheme. All four curves are real. Only FTCS with r = 0.6 leaves the band between −1 and 1. Hover over the chart, or focus it and use the arrow keys, to read the values.</figcaption>
<details class="viz-table"><summary>Show the numbers</summary><table><thead><tr><th>ξh</th><th>FTCS, r = 0.3</th><th>FTCS, r = 0.5</th><th>FTCS, r = 0.6</th><th>Crank–Nicolson, r = 6</th></tr></thead><tbody><tr><td>0</td><td>+1.0000</td><td>+1.0000</td><td>+1.0000</td><td>+1.0000</td></tr><tr><td>π/8</td><td>+0.9543</td><td>+0.9239</td><td>+0.9087</td><td>+0.3729</td></tr><tr><td>π/4</td><td>+0.8243</td><td>+0.7071</td><td>+0.6485</td><td>-0.2747</td></tr><tr><td>3π/8</td><td>+0.6296</td><td>+0.3827</td><td>+0.2592</td><td>-0.5748</td></tr><tr><td>π/2</td><td>+0.4000</td><td>+0.0000</td><td>-0.2000</td><td>-0.7143</td></tr><tr><td>5π/8</td><td>+0.1704</td><td>-0.3827</td><td>-0.6592</td><td>-0.7849</td></tr><tr><td>3π/4</td><td>-0.0243</td><td>-0.7071</td><td>-1.0485</td><td>-0.8221</td></tr><tr><td>7π/8</td><td>-0.1543</td><td>-0.9239</td><td>-1.3087</td><td>-0.8406</td></tr><tr><td>π</td><td>-0.2000</td><td>-1.0000</td><td>-1.4000</td><td>-0.8462</td></tr></tbody></table></details>
<script>(() => {const fig = document.currentScript.closest('figure');const svg = fig.querySelector('svg');const line = svg.querySelector('.viz-crosshair');const tip = fig.querySelector('.viz-tooltip');const names = ["FTCS, r = 0.3", "FTCS, r = 0.5", "FTCS, r = 0.6", "Crank\u2013Nicolson, r = 6"]; const colors = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"];const rows = [["+1.0000", "+1.0000", "+1.0000", "+1.0000"], ["+0.9885", "+0.9808", "+0.9769", "+0.7933"], ["+0.9543", "+0.9239", "+0.9087", "+0.3729"], ["+0.8989", "+0.8315", "+0.7978", "-0.0056"], ["+0.8243", "+0.7071", "+0.6485", "-0.2747"], ["+0.7333", "+0.5556", "+0.4667", "-0.4545"], ["+0.6296", "+0.3827", "+0.2592", "-0.5748"], ["+0.5171", "+0.1951", "+0.0341", "-0.6569"], ["+0.4000", "+0.0000", "-0.2000", "-0.7143"], ["+0.2829", "-0.1951", "-0.4341", "-0.7552"], ["+0.1704", "-0.3827", "-0.6592", "-0.7849"], ["+0.0667", "-0.5556", "-0.8667", "-0.8065"], ["-0.0243", "-0.7071", "-1.0485", "-0.8221"], ["-0.0989", "-0.8315", "-1.1978", "-0.8332"], ["-0.1543", "-0.9239", "-1.3087", "-0.8406"], ["-0.1885", "-0.9808", "-1.3769", "-0.8448"], ["-0.2000", "-1.0000", "-1.4000", "-0.8462"]]; const heads = ["\u03beh = 0", "\u03beh = \u03c0/16", "\u03beh = \u03c0/8", "\u03beh = 3\u03c0/16", "\u03beh = \u03c0/4", "\u03beh = 5\u03c0/16", "\u03beh = 3\u03c0/8", "\u03beh = 7\u03c0/16", "\u03beh = \u03c0/2", "\u03beh = 9\u03c0/16", "\u03beh = 5\u03c0/8", "\u03beh = 11\u03c0/16", "\u03beh = 3\u03c0/4", "\u03beh = 13\u03c0/16", "\u03beh = 7\u03c0/8", "\u03beh = 15\u03c0/16", "\u03beh = \u03c0"];const xs = [64.0, 91.2, 118.5, 145.8, 173.0, 200.2, 227.5, 254.8, 282.0, 309.2, 336.5, 363.7, 391.0, 418.3, 445.5, 472.8, 500.0]; const steps = 16, width = 680;let current = -1;const show = (n) => { n = Math.max(0, Math.min(steps, n)); current = n; const x = xs[n]; line.setAttribute('x1', x); line.setAttribute('x2', x); line.setAttribute('visibility', 'visible'); tip.replaceChildren(); const head = document.createElement('div'); head.className = 'viz-tooltip-head'; head.textContent = heads[n]; tip.append(head); names.forEach((name, i) => {  const row = document.createElement('div'); row.className = 'viz-tooltip-row';  const key = document.createElement('span'); key.className = 'viz-swatch'; key.style.background = colors[i];  const value = document.createElement('strong'); value.textContent = rows[n][i];  const label = document.createElement('span'); label.textContent = name;  row.append(key, value, label); tip.append(row); }); const box = svg.getBoundingClientRect(), frame = fig.getBoundingClientRect(); const xp = box.left - frame.left + x * box.width / width; tip.hidden = false; const room = frame.width - tip.offsetWidth - 8; tip.style.left = Math.max(0, Math.min(room, xp + 12)) + 'px'; tip.style.top = (box.top - frame.top + 8) + 'px';};const hide = () => { line.setAttribute('visibility', 'hidden'); tip.hidden = true; current = -1; };svg.addEventListener('pointermove', (event) => { const box = svg.getBoundingClientRect(); const x = (event.clientX - box.left) * width / box.width; let best = 0; xs.forEach((v, i) => { if (Math.abs(v - x) < Math.abs(xs[best] - x)) best = i; }); show(best);});svg.addEventListener('pointerleave', hide);svg.addEventListener('focus', () => show(current < 0 ? 0 : current));svg.addEventListener('blur', hide);svg.addEventListener('keydown', (event) => { if (event.key === 'ArrowRight') { show(current + 1); event.preventDefault(); } if (event.key === 'ArrowLeft') { show(current - 1); event.preventDefault(); } if (event.key === 'Escape') hide();});})();</script>
</figure>

| scheme | $r$ | $g(\pi/4)$ | $g(\pi/2)$ | $g(\pi)$ | largest $\lvert g\rvert$ |
|---|---|---|---|---|---|
| FTCS | 0.3 | $0.8243$ | $0.4000$ | $-0.2000$ | $1.0000$ |
| FTCS | 0.5 | $0.7071$ | $0.0000$ | $-1.0000$ | $1.0000$ |
| FTCS | 0.6 | $0.6485$ | $-0.2000$ | $-1.4000$ | $1.4000$ |
| Crank–Nicolson | 0.6 | $0.7011$ | $0.2500$ | $-0.0909$ | $1.0000$ |
| Crank–Nicolson | 6 | $-0.2747$ | $-0.7143$ | $-0.8462$ | $1.0000$ |
| Crank–Nicolson | 60 | $-0.8923$ | $-0.9672$ | $-0.9835$ | $1.0000$ |
| Backward Euler | 6 | $0.2215$ | $0.0769$ | $0.0400$ | $1.0000$ |

- The largest $|g|$ is always $g(0) = 1$, the constant mode, except for FTCS with $r = 0.6$. There the sawtooth has $g = 1 - 4(0.6) = -1.4$, so it grows by $40\%$ per step and flips sign. At the edge $r = 0.5$, the sawtooth has exactly $g = -1$.
- Crank–Nicolson stays in $[-1, 1]$ for $r = 0.6$, $6$ and $60$. But as $r$ grows, the fast modes get $g$ closer and closer to $-1$. Backward Euler with the same $r = 6$ damps the sawtooth to $0.04$.

**Part 3: advection, FTCS and Lax–Friedrichs.** The same one-step measurement for $u_t + au_x = 0$ on the $32$-point periodic grid gives the largest $|g|$ over all modes:

| $\nu$ | 0.3 | 0.9 | 1.0 | 1.1 |
|---|---|---|---|---|
| FTCS | $1.0440$ | $1.3454$ | $1.4142$ | $1.4866$ |
| Lax–Friedrichs | $1.0000$ | $1.0000$ | $1.0000$ | $1.1000$ |

FTCS is above $1$ for every $\nu$, even small ones. Lax–Friedrichs is at most $1$ exactly when $\nu \le 1$. Then I ran $300$ steps from a smooth bump $e^{-100(x - 1/2)^2}$ on a periodic grid with $100$ points and recorded the largest value:

| run | $n = 0$ | $n = 100$ | $n = 200$ | $n = 300$ | growth per step (last 100 steps) |
|---|---|---|---|---|---|
| FTCS, $\nu = 0.3$ | $1.000$ | $1.097$ | $1.228$ | $1.390$ | $1.0012$ |
| FTCS, $\nu = 0.9$ | $1.000$ | $21.1$ | $4.31 \times 10^{13}$ | $2.39 \times 10^{26}$ | $1.3410$ |
| Lax–Friedrichs, $\nu = 0.9$ | $1.000$ | $0.852$ | $0.754$ | $0.684$ | $0.9990$ |
| Lax–Friedrichs, $\nu = 1.1$ | $1.000$ | $1.286$ | $1.931$ | $6.624$ | $1.0124$ |

- FTCS with $\nu = 0.9$ explodes. In the end it grows by $1.341$ per step, close to the largest $|g| = 1.345$ for this $\nu$ (on the finer grid of this run, the worst mode is not hit exactly).
- FTCS with $\nu = 0.3$ grows only slowly here. The bump is smooth, so it has almost nothing in the fast modes where $|g|$ is largest. But the scheme is still unstable, and it would keep growing. A small growth that never stops is still a failure.
- Lax–Friedrichs with $\nu = 0.9$ stays bounded. The peak slowly shrinks: the scheme smears the bump as it moves. (Why it smears is a question for a later lecture.)
- Lax–Friedrichs with $\nu = 1.1$ grows, and the growth speeds up ($1.29$, $1.93$, then $6.62$). The bump has very little content near the worst mode, so it takes time before the factor $1.1$ per step takes over.

**Part 4: the condition number of 2D Crank–Nicolson.** I use $h = 1/32$ ($31 \times 31$ unknowns). The condition numbers come from the eigenvalue formula of section 5. As a check, I computed the largest eigenvalue of each matrix by power iteration (no formula). The relative difference is below $5 \times 10^{-5}$ in every case.

| matrix | $k/h^2$ | $\kappa$ |
|---|---|---|
| Laplacian $-\nabla_h^2$ | | $414.3$ |
| Crank–Nicolson, $k = 4h$ | $128$ | $229.2$ |
| Crank–Nicolson, $k = h$ | $32$ | $98.4$ |
| Crank–Nicolson, $k = h/4$ | $8$ | $30.6$ |
| Crank–Nicolson, $k = h^2$ | $1$ | $4.94$ |

- Every Crank–Nicolson matrix is better conditioned than the Laplacian, and the smaller $k$, the better.
- For $k/h^2 = 1$ the condition number is below $5$. For very large $k$ it approaches the Laplacian's value, because then $I$ hardly matters.

The full script and the chart code are in [this site's repository](https://github.com/doubleyands/doubleyands.github.io/tree/main/scripts/posts/numerical-pde/lecture-11).

## Summary

1. $e^{i\xi x}$ is an eigenfunction of $\partial_x$, and the grid mode $e^{ijh\xi}$ is an eigenvector of every difference operator that looks the same at every point. For the central difference the eigenvalue is $\frac{i}{h}\sin(h\xi) = i\xi + O(h^2\xi^3)$.
2. The error grows with the wave number. The grid can show only $|\xi h| \le \pi$, and faster waves are aliased onto slower ones.
3. **von Neumann analysis**: put $e^{ijh\xi}$ into the scheme, read off the **amplification factor** $g(\xi)$, and require $|g(\xi)| \le 1$ (or $1 + \alpha k$) for every $\xi$. Parseval's relation turns this into stability.
4. For the heat equation, FTCS has $g = 1 + \frac{2k}{h^2}(\cos\xi h - 1)$ and needs $k/h^2 \le 1/2$. Crank–Nicolson has $|g| \le 1$ for every $k$. These match the eigenvalue results of Lecture 10.
5. In 2D, Crank–Nicolson solves a system with $\kappa = O(k/h^2)$, smaller than the Laplacian's $O(1/h^2)$, and $U^n$ is a good starting guess for an iterative solver.
6. The advection equation $u_t + au_x = 0$ moves its initial shape: $u = \eta(x - at)$. FTCS is unstable for it. Lax–Friedrichs is stable under the **CFL condition** $|ak/h| \le 1$, which allows $k = O(h)$ for an explicit method.
