---
title: "Numerical PDEs, Lecture 10: Stability of the Heat Equation, Condition Numbers, and the Lax Equivalence Theorem"
description: "My notes on the tenth lecture: the eigenvalues of the second-difference matrix, the condition number, why FTCS needs k/h² ≤ 1/2 while Crank–Nicolson is unconditionally stable, Lax–Richtmyer stability and the Lax equivalence theorem, and the idea behind von Neumann analysis."
pubDate: 2026-11-05
project: numerical-pde
tags: [numerical PDEs, heat equation, stability, eigenvalues, linear algebra, notes]
---

This is my summary of the tenth lecture in my numerical PDE course. I write the ideas in my own words and check them with a small experiment.

[Lecture 9](/blog/2026-11-03-numerical-pde-lecture-9/) ended with the heat equation $u_t = u_{xx}$ and two schemes: **FTCS** (Forward Euler in time, central difference in space) and **Crank–Nicolson** (the trapezoid in time). We found their truncation errors, $O(k + h^2)$ and $O(k^2 + h^2)$. So both schemes are consistent. The open question was stability. This lecture answers it with the **eigenvalues** of the space matrix. Then it states the theorem that makes stability so important, and it starts a second tool, **von Neumann analysis**.

I use the same notation as before: $h$ is the grid spacing, $k$ is the time step, and the grid points are $x_j = jh$ with $h = 1/(N+1)$.

## 1. A time-dependent problem is many steady problems

For a steady diffusion problem, we solve one linear system and we are done. An implicit scheme for the heat equation solves a linear system **at every time step**. Crank–Nicolson, for example, solves a tridiagonal system each step. With $100$ steps, we solve $100$ systems. So time-dependent problems cost much more, and the size of the time step matters a lot: if we can take steps four times larger, we solve four times fewer systems.

| | explicit (FTCS) | implicit (Backward Euler, Crank–Nicolson) |
|---|---|---|
| one step | no system to solve, just evaluate | one tridiagonal solve |
| time step | must be small, as we will see | no stability limit |

"No stability limit" does not mean "any step is fine". The step must still be small enough for the scheme to **approximate** the equation. We come back to this in section 6.

## 2. Consistency + stability = convergence

Writing down a scheme is easy. The hard question is whether it is **stable**. A theorem of Lax says: if a (linear) scheme is consistent and stable, then it converges. We already have consistency from Lecture 9. So the only thing left to check is stability. In addition, for a stable scheme, the global error has the same order as the truncation error. Section 7 states the theorem properly and proves the easy direction.

## 3. The eigenvalues of the second-difference matrix

Discretize $u_{xx}$ with the central difference ([Lecture 3](/blog/2026-10-13-numerical-pde-lecture-3/)). With $u = 0$ at both ends, we get the $N \times N$ matrix

$$
A = \frac{1}{h^2}\begin{bmatrix}-2&1&&\\1&-2&1&\\&\ddots&\ddots&\ddots\\&&1&-2\end{bmatrix}.
$$

Its eigenvalues decide almost everything in this lecture. We cannot compute them by hand for a large $N$. But we can **guess them from the continuous problem** and then check the guess.

**The continuous problem.** Solve $u'' = \lambda u$ with $u(0) = u(1) = 0$. The solutions are

$$
u = \sin(p\pi x), \qquad \lambda = -(p\pi)^2, \qquad p = 1, 2, 3, \dots
$$

Differentiate the sine twice: it comes back with a minus sign and a factor $(p\pi)^2$.

**The guess.** The eigenvectors of $A$ are these functions sampled at the grid points:

$$
u^p_j = \sin(p\pi jh), \qquad j = 1, \dots, N.
$$

**The check.** Put $u^p$ into row $j$ of $A$ and use $\sin(a - b) + \sin(a + b) = 2\sin a\cos b$:

$$
(Au^p)_j = \frac{1}{h^2}\Big(\sin(p\pi(j-1)h) - 2\sin(p\pi jh) + \sin(p\pi(j+1)h)\Big) = \frac{2}{h^2}\big(\cos(p\pi h) - 1\big)\sin(p\pi jh).
$$

The first and last rows have only one neighbor. But the sine is zero at $j = 0$ and $j = N + 1$, so the same formula still holds there. So

$$
\boxed{\;\lambda_p = \frac{2}{h^2}\big(\cos(p\pi h) - 1\big), \qquad p = 1, \dots, N.\;}
$$

Some facts follow at once.

- All eigenvalues are **negative real numbers**, because $\cos \le 1$.
- The continuous problem has infinitely many eigenfunctions. The grid only sees the first $N$.
- The one closest to zero: $\cos x - 1 \approx -x^2/2$ with $x = \pi h$ gives $\lambda_1 \approx -\pi^2$. It does not depend on $h$.
- The one farthest from zero: $\cos(N\pi h) \approx -1$ gives $\lambda_N \approx -4/h^2$. It grows as $h \to 0$.
- For a fixed $p$, $\lambda_p = -(p\pi)^2 + O(h^2)$. For large $p$, the discrete values stop at about $-4/h^2$.
- The eigenvector for $p = N$ changes sign at almost every grid point. It is a **sawtooth**, the fastest wiggle the grid can show.

Because $|\lambda_1| \approx \pi^2$ stays away from zero, $\|A^{-1}\|_2 = 1/|\lambda_1| \approx 1/\pi^2$ for every $h$. This is exactly the stability of the steady problem in the sense of [Lecture 7](/blog/2026-10-27-numerical-pde-lecture-7/).

## 4. The condition number

A linear system is **well conditioned** if a small change in the data gives a small change in the answer. We measure this with relative changes. The number that controls it is the **condition number**

$$
\kappa(A) = \|A\|\,\|A^{-1}\| = \frac{\sigma_{\max}(A)}{\sigma_{\min}(A)} \quad \text{(in the 2-norm)}.
$$

- The **singular values** $\sigma$ are the square roots of the eigenvalues of $A^{\mathsf T}A$. This matrix is symmetric and its eigenvalues are real and at least $0$, so the square roots make sense.
- If $A$ is symmetric, $\kappa(A) = |\lambda|_{\max} / |\lambda|_{\min}$.
- Always $\kappa \ge 1$: $1 = \|I\| = \|AA^{-1}\| \le \|A\|\,\|A^{-1}\|$. The best value $1$ belongs to the identity matrix.
- A smaller condition number also helps iterative solvers converge faster.

**A picture.** A $2 \times 2$ system is the crossing point of two lines. A computer stores the coefficients with a small rounding error, so each line becomes a thin band. The computed answer can be anywhere in the region where the two bands overlap. If the lines cross at a large angle, this region is small. If the lines are almost parallel, the region is long and thin, and the computed answer can be far from the true one. The computer did nothing wrong. The problem itself is sensitive. This is why numerical linear algebra cares about more than "$A$ is invertible".

**Our matrix.** $A$ is symmetric, so

$$
\kappa(A) = \frac{|\lambda_N|}{|\lambda_1|} \approx \frac{4/h^2}{\pi^2} = \frac{4}{\pi^2 h^2}.
$$

The order is $1/h^2$, and the constant is $4/\pi^2 \approx 0.4$. Note the two different facts: $\|A^{-1}\|$ stays bounded (the problem is **stable**), but $\kappa(A)$ grows without bound (the problem gets **ill conditioned**). Rounding errors get more and more amplified as $h$ shrinks. Part 2 of the experiment measures both.

## 5. FTCS: put the eigenvalues into the stability region

Here is the key idea. Discretize only in space first. Then the heat equation becomes a system of ODEs,

$$
U'(t) = AU(t).
$$

This is called the **method of lines**. $A$ is diagonalizable, so in the eigenvector basis the system splits into $N$ scalar ODEs,

$$
w_p' = \lambda_p w_p, \qquad p = 1, \dots, N.
$$

These are exactly the test equations of [Lecture 8](/blog/2026-10-29-numerical-pde-lecture-8/). FTCS is Forward Euler applied to this system, so we need $z = k\lambda_p$ in the Forward Euler region for **every** $p$. All $z$ are real and lie in $(-4k/h^2,\ 0)$, and the real part of the Forward Euler region is $[-2, 0]$. So only the farthest eigenvalue matters:

$$
\Big|1 - \frac{4k}{h^2}\Big| \le 1 \quad\Longrightarrow\quad \boxed{\;\frac{k}{h^2} \le \frac12.\;}
$$

(More precisely $|1 + k\lambda_N| \le 1$, and $|\lambda_N| < 4/h^2$, so $k/h^2 \le 1/2$ is enough.)

- Each eigenvalue gives its own condition on $k$, and we need all of them, so we take the intersection. The condition from $\lambda_1 \approx -\pi^2$ is $k \le 2/\pi^2$, which does not depend on $h$ and is no real limit. The strict one comes from $\lambda_N$.
- The edge $k/h^2 = 1/2$ is allowed, but in practice it is safer to stay a bit inside.
- This is a **severe** limit. With $h = 1/100$, we need $k \le 5 \times 10^{-5}$. When we halve $h$, we must cut $k$ by **four**. Tiny steps take a very long time to reach a final time.

**Why the blow-up looks like a sawtooth.** If $k/h^2 > 1/2$, the modes with large $p$ leave the region first. Those are the sawtooth-like eigenvectors. Each step multiplies them by $|1 + k\lambda_p| > 1$, so a sawtooth grows out of the solution. Part 3 of the experiment shows this.

## 6. Crank–Nicolson is unconditionally stable

Crank–Nicolson is the trapezoid in time. The trapezoid is stable on the whole left half plane ([Lecture 9](/blog/2026-11-03-numerical-pde-lecture-9/)). And

$$
z = k\lambda_p = \frac{2k}{h^2}\big(\cos(p\pi h) - 1\big) < 0 \qquad \text{for every } p \text{ and every } k.
$$

So Crank–Nicolson is stable for **any** $k$. We call this **unconditionally stable**.

This is a trade-off. Each implicit step costs more, because we solve a system. In return, there is no limit on the step. Forward Euler is cheap per step, but we pay later with very many steps. Neither one is always better. Explicit methods are simple to code, and some fields use them for diffusion problems on purpose.

**Unconditionally stable does not mean any $k$ is fine.** If $k$ is too large, the scheme is stable but not accurate. In practice we take $k = O(h)$, and the unconditional stability is what allows that. Part 4 of the experiment shows one concrete reason (this part is my own check, not from the lecture). For the sawtooth mode, the Crank–Nicolson factor is

$$
\frac{1 + z/2}{1 - z/2} \to -1 \qquad \text{as } z \to -\infty.
$$

So with a large $k$, the fastest modes, which should die almost at once, barely shrink and flip their sign at every step.

## 7. Lax–Richtmyer stability and the Lax equivalence theorem

All our schemes for the heat equation can be written as

$$
U^{n+1} = B(k)\,U^n + b^n(k),
$$

where $b^n$ holds the boundary data and sources.

| scheme | $B(k)$ | compare: $P(z)$ from Lecture 9 |
|---|---|---|
| Forward Euler (FTCS) | $I + kA$ | $1 + z$ |
| Backward Euler | $(I - kA)^{-1}$ | $1/(1 - z)$ |
| Crank–Nicolson | $(I - \tfrac{k}{2}A)^{-1}(I + \tfrac{k}{2}A)$ | $(1 + z/2)/(1 - z/2)$ |

If we think of $A$ as a number $\lambda$, these are exactly the functions $P(z)$ of the last lecture.

**Matrix norm.** $\|A\| = \max_{x \ne 0} \|Ax\| / \|x\|$: the most that $A$ can stretch a vector. If $A$ is **symmetric**, the 2-norm equals the **spectral radius** $\rho(A) = \max_p |\lambda_p|$. So for symmetric matrices, we can compute the norm from the eigenvalues.

**Definition (Lax–Richtmyer stability).** For every final time $T$ there is a constant $M$ such that

$$
\|B(k)^n\| \le M \qquad \text{for all } k > 0 \text{ and } n \text{ with } kn \le T.
$$

Here $B(k)^n$ is a real power: after $n$ steps, $B$ has been applied $n$ times.

**Theorem (Lax equivalence).** A consistent linear method converges **if and only if** it is Lax–Richtmyer stable.

Proving convergence directly is hard. Proving stability is much easier. The theorem lets us focus on stability alone.

**Proof that stability gives convergence.** Put the exact solution into the scheme. It misses by the truncation error:

$$
u^{n+1} = Bu^n + b^n + k\tau^n.
$$

Subtract the scheme. The $b^n$ cancel, and the error $E^n = U^n - u^n$ satisfies $E^{n+1} = BE^n - k\tau^n$. Repeat this $m$ times:

$$
E^m = B^mE^0 - k\sum_{n=1}^{m} B^{m-n}\tau^{n-1}.
$$

Take norms and use $\|B^j\| \le M$ and $km \le T$:

$$
\|E^m\| \le M\|E^0\| + k\sum_{n=1}^{m} M\,\|\tau^{n-1}\| \le M\|E^0\| + TM \max_{n} \|\tau^{n-1}\|.
$$

Consistency gives $\|\tau\| \to 0$, and a correct start gives $\|E^0\| \to 0$. So $\|E^m\| \to 0$. This is the same structure as the convergence proof for Forward Euler in Lecture 8. There, $(1 + kL)^n \le e^{LT}$ did the job. Here, $\|B^n\| \le M$ does it.

## 8. Computing the norms with eigenvalues

If $Au = \lambda u$, then $(A + 2I)u = (\lambda + 2)u$, $A^2u = \lambda^2 u$, and in general

$$
P(A)\,u = P(\lambda)\,u
$$

for a polynomial or a ratio of polynomials $P$. Applying a function to the matrix applies the same function to each eigenvalue, and the eigenvectors stay the same. Our $B(k)$ is a function of the symmetric matrix $A$, so it is symmetric too, and its 2-norm is its largest eigenvalue in absolute value.

| scheme | eigenvalues of $B(k)$ | $\lVert B(k)\rVert_2 \le 1$ when | result |
|---|---|---|---|
| FTCS | $1 + k\lambda_p$ | $k/h^2 \le 1/2$ | conditionally stable |
| Backward Euler | $1/(1 - k\lambda_p)$ | always | unconditionally stable |
| Crank–Nicolson | $\dfrac{1 + k\lambda_p/2}{1 - k\lambda_p/2}$ | always | unconditionally stable |

For Crank–Nicolson, the denominator looks like a minus, but $\lambda_p < 0$, so it is really $1 + k|\lambda_p|/2$. It is always larger than the numerator in size, so the ratio is at most $1$. If $\|B\| \le 1$, then $\|B^n\| \le 1$, and the scheme is Lax–Richtmyer stable with $M = 1$.

In summary: FTCS is $O(k + h^2)$ and needs $k/h^2 \le 1/2$. Backward Euler is $O(k + h^2)$ with no limit. Crank–Nicolson is $O(k^2 + h^2)$ with no limit. If you can live with the FTCS limit, FTCS is fine. If not, choose an implicit scheme.

## 9. When eigenvalues are out of reach: von Neumann analysis

We were lucky. We could find the eigenvalues exactly because the grid is uniform, the coefficient is constant, and the boundary condition is simple. With a non-uniform grid, an advection term, or a diffusion coefficient that changes in space, the matrix loses this structure and we cannot write down its eigenvalues.

**von Neumann analysis** goes around this problem with Fourier analysis. It assumes:

| assumption | why |
|---|---|
| a linear PDE with **constant coefficients** | so that the Fourier modes do not mix |
| **no boundary** ($-\infty < x < \infty$, the Cauchy problem) or **periodic** boundaries | so that we can expand in Fourier modes on the whole domain |

A **grid function** $V_j$ lives only on the grid points $x_j = jh$. We write it as a combination of the modes $e^{ijh\xi} = e^{ix_j\xi}$, where $\xi$ is the **wave number**. (Here $i$ in the exponent is $\sqrt{-1}$, and $j$ is the grid index.)

$$
V_j = \frac{1}{\sqrt{2\pi}}\int_{-\pi/h}^{\pi/h}\hat V(\xi)\,e^{ijh\xi}\,d\xi, \qquad \hat V(\xi) = \frac{h}{\sqrt{2\pi}}\sum_{j=-\infty}^{\infty}V_j\,e^{-ijh\xi}.
$$

- Only wave numbers in $[-\pi/h, \pi/h]$ are different on the grid. A faster wave looks the same as a slower one at the grid points. The fastest wave the grid can show is the sawtooth again.
- The matching grid 2-norm has a factor $h$: $\|V\|_2 = \big(h\sum_j |V_j|^2\big)^{1/2}$. With this $h$, Parseval's relation $\|\hat V\|_2 = \|V\|_2$ holds. It is the same $h$ as in the discrete $L^2$ norm of [Lecture 5](/blog/2026-10-20-numerical-pde-lecture-5/).

The lecture stopped here. The next step, putting one Fourier mode into a scheme and reading off how much it grows per step, is the topic of the next lecture.

## 10. Experiment

The script uses the Python standard library only. The problem is $u_t = u_{xx}$ on $[0, 1]$ with $u = 0$ at both ends. In this post, $r = k/h^2$. (In Lecture 9, I used $r = k/(2h^2)$ for Crank–Nicolson.)

**Part 1: the eigenpairs.** With $N = 15$ ($h = 1/16$), I put each $u^p$ into $A$ and measured $\|Au^p - \lambda_p u^p\| / (|\lambda_p|\,\|u^p\|)$. The largest value over all $15$ modes is $7.5 \times 10^{-15}$, which is rounding error. So the formula is right. Compared with the continuous values:

| $p$ | 1 | 2 | 4 | 8 | 12 | 15 |
|---|---|---|---|---|---|---|
| $\lambda_p$ | $-9.84$ | $-38.97$ | $-149.96$ | $-512.00$ | $-874.04$ | $-1014.16$ |
| $-(p\pi)^2$ | $-9.87$ | $-39.48$ | $-157.91$ | $-631.65$ | $-1421.22$ | $-2220.66$ |

Small $p$ agree well. Large $p$ do not, and the discrete values stop near $4/h^2 = 1024$.

**Part 2: norms and condition numbers.** I computed $\|A\|_2$ by power iteration and $\|A^{-1}\|_2$ by inverse iteration (each step is one tridiagonal solve). These do not use the eigenvalue formula.

| $h$ | $\lVert A\rVert_2$ | $4/h^2$ | $\lVert A^{-1}\rVert_2$ | $\kappa(A)$ | $\lvert\lambda_N\rvert/\lvert\lambda_1\rvert$ | $4/(\pi^2h^2)$ |
|---|---|---|---|---|---|---|
| 1/8 | $246.3$ | $256$ | $0.10263$ | $25.3$ | $25.3$ | $25.9$ |
| 1/16 | $1014.2$ | $1024$ | $0.10165$ | $103.1$ | $103.1$ | $103.8$ |
| 1/32 | $4086.1$ | $4096$ | $0.10140$ | $414.3$ | $414.3$ | $415.0$ |
| 1/64 | $16374.1$ | $16384$ | $0.10134$ | $1659.4$ | $1659.4$ | $1660.0$ |
| 1/128 | $65526.1$ | $65536$ | $0.10133$ | $6639.5$ | $6639.5$ | $6640.2$ |
| 1/256 | $262134.0$ | $262144$ | $0.10132$ | $26560.1$ | $26560.1$ | $26560.7$ |

- $\|A^{-1}\|_2$ goes to $1/\pi^2 = 0.10132$: the steady problem is stable.
- $\kappa(A)$ grows by $4$ every time $h$ is halved, and matches $4/(\pi^2h^2)$ more and more closely. At $h = 1/256$ it is already above $26{,}000$.
- The measured $\kappa$ agrees with the eigenvalue formula to every printed digit.

**Part 3: the FTCS limit.** I use $h = 1/25$ ($N = 24$), a hat-shaped start $u(x, 0) = 2\min(x, 1 - x)$ with height $1$, and $200$ steps of FTCS with $r = 0.45$, $0.50$ and $0.55$.

<figure class="viz">
<div class="viz-legend"><span class="viz-key"><span class="viz-swatch" style="background:#2a78d6"></span>r = 0.45</span><span class="viz-key"><span class="viz-swatch" style="background:#eb6834"></span>r = 0.50</span><span class="viz-key"><span class="viz-swatch" style="background:#1baf7a"></span>r = 0.55</span></div>
<div class="viz-frame">
<svg viewBox="0 0 680 400" role="img" tabindex="0" aria-label="Largest value of the FTCS solution of the heat equation over 200 time steps, on a log scale, for r = k/h squared equal to 0.45, 0.50 and 0.55. The first two decay slowly to about 0.2. The third grows to about 8 billion." class="viz-plot"><line x1="64" x2="520" y1="42.3" y2="42.3" class="viz-grid"/><text x="54" y="46.3" text-anchor="end" class="viz-tick">10¹⁰</text><line x1="64" x2="520" y1="94.0" y2="94.0" class="viz-grid"/><text x="54" y="98.0" text-anchor="end" class="viz-tick">10⁸</text><line x1="64" x2="520" y1="145.7" y2="145.7" class="viz-grid"/><text x="54" y="149.7" text-anchor="end" class="viz-tick">10⁶</text><line x1="64" x2="520" y1="197.3" y2="197.3" class="viz-grid"/><text x="54" y="201.3" text-anchor="end" class="viz-tick">10⁴</text><line x1="64" x2="520" y1="249.0" y2="249.0" class="viz-grid"/><text x="54" y="253.0" text-anchor="end" class="viz-tick">100</text><line x1="64" x2="520" y1="300.7" y2="300.7" class="viz-grid"/><text x="54" y="304.7" text-anchor="end" class="viz-tick">1</text><line x1="64" x2="520" y1="326.5" y2="326.5" class="viz-grid"/><text x="54" y="330.5" text-anchor="end" class="viz-tick">0.1</text><text x="64.0" y="362" text-anchor="middle" class="viz-tick">0</text><text x="178.0" y="362" text-anchor="middle" class="viz-tick">50</text><text x="292.0" y="362" text-anchor="middle" class="viz-tick">100</text><text x="406.0" y="362" text-anchor="middle" class="viz-tick">150</text><text x="520.0" y="362" text-anchor="middle" class="viz-tick">200</text><text x="292.0" y="388" text-anchor="middle" class="viz-axis-title">Time step n</text><text x="64" y="14" class="viz-axis-title">Largest value max |U| on the grid</text><path d="M64.0,301.1 L68.6,301.8 L73.1,302.2 L77.7,302.5 L82.2,302.8 L86.8,303.1 L91.4,303.4 L95.9,303.6 L100.5,303.8 L105.0,304.0 L109.6,304.3 L114.2,304.5 L118.7,304.7 L123.3,304.9 L127.8,305.0 L132.4,305.2 L137.0,305.4 L141.5,305.6 L146.1,305.8 L150.6,305.9 L155.2,306.1 L159.8,306.3 L164.3,306.5 L168.9,306.6 L173.4,306.8 L178.0,307.0 L182.6,307.1 L187.1,307.3 L191.7,307.5 L196.2,307.6 L200.8,307.8 L205.4,308.0 L209.9,308.1 L214.5,308.3 L219.0,308.5 L223.6,308.6 L228.2,308.8 L232.7,308.9 L237.3,309.1 L241.8,309.3 L246.4,309.4 L251.0,309.6 L255.5,309.7 L260.1,309.9 L264.6,310.1 L269.2,310.2 L273.8,310.4 L278.3,310.6 L282.9,310.7 L287.4,310.9 L292.0,311.0 L296.6,311.2 L301.1,311.4 L305.7,311.5 L310.2,311.7 L314.8,311.8 L319.4,312.0 L323.9,312.2 L328.5,312.3 L333.0,312.5 L337.6,312.6 L342.2,312.8 L346.7,313.0 L351.3,313.1 L355.8,313.3 L360.4,313.4 L365.0,313.6 L369.5,313.8 L374.1,313.9 L378.6,314.1 L383.2,314.2 L387.8,314.4 L392.3,314.6 L396.9,314.7 L401.4,314.9 L406.0,315.0 L410.6,315.2 L415.1,315.3 L419.7,315.5 L424.2,315.7 L428.8,315.8 L433.4,316.0 L437.9,316.1 L442.5,316.3 L447.0,316.5 L451.6,316.6 L456.2,316.8 L460.7,316.9 L465.3,317.1 L469.8,317.3 L474.4,317.4 L479.0,317.6 L483.5,317.7 L488.1,317.9 L492.6,318.1 L497.2,318.2 L501.8,318.4 L506.3,318.5 L510.9,318.7 L515.4,318.9 L520.0,319.0" fill="none" stroke="#2a78d6" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><path d="M64.0,301.1 L68.6,301.8 L73.1,302.3 L77.7,302.6 L82.2,303.0 L86.8,303.3 L91.4,303.5 L95.9,303.8 L100.5,304.0 L105.0,304.3 L109.6,304.5 L114.2,304.7 L118.7,304.9 L123.3,305.1 L127.8,305.3 L132.4,305.5 L137.0,305.7 L141.5,305.9 L146.1,306.1 L150.6,306.3 L155.2,306.5 L159.8,306.7 L164.3,306.9 L168.9,307.1 L173.4,307.3 L178.0,307.4 L182.6,307.6 L187.1,307.8 L191.7,308.0 L196.2,308.2 L200.8,308.3 L205.4,308.5 L209.9,308.7 L214.5,308.9 L219.0,309.1 L223.6,309.2 L228.2,309.4 L232.7,309.6 L237.3,309.8 L241.8,310.0 L246.4,310.1 L251.0,310.3 L255.5,310.5 L260.1,310.7 L264.6,310.9 L269.2,311.0 L273.8,311.2 L278.3,311.4 L282.9,311.6 L287.4,311.7 L292.0,311.9 L296.6,312.1 L301.1,312.3 L305.7,312.5 L310.2,312.6 L314.8,312.8 L319.4,313.0 L323.9,313.2 L328.5,313.3 L333.0,313.5 L337.6,313.7 L342.2,313.9 L346.7,314.1 L351.3,314.2 L355.8,314.4 L360.4,314.6 L365.0,314.8 L369.5,314.9 L374.1,315.1 L378.6,315.3 L383.2,315.5 L387.8,315.7 L392.3,315.8 L396.9,316.0 L401.4,316.2 L406.0,316.4 L410.6,316.5 L415.1,316.7 L419.7,316.9 L424.2,317.1 L428.8,317.3 L433.4,317.4 L437.9,317.6 L442.5,317.8 L447.0,318.0 L451.6,318.1 L456.2,318.3 L460.7,318.5 L465.3,318.7 L469.8,318.9 L474.4,319.0 L479.0,319.2 L483.5,319.4 L488.1,319.6 L492.6,319.7 L497.2,319.9 L501.8,320.1 L506.3,320.3 L510.9,320.5 L515.4,320.6 L520.0,320.8" fill="none" stroke="#eb6834" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><path d="M64.0,301.1 L68.6,301.9 L73.1,302.4 L77.7,302.7 L82.2,303.1 L86.8,303.4 L91.4,303.7 L95.9,303.9 L100.5,304.2 L105.0,304.4 L109.6,304.7 L114.2,304.9 L118.7,305.1 L123.3,305.3 L127.8,305.5 L132.4,305.7 L137.0,305.9 L141.5,305.9 L146.1,305.8 L150.6,305.6 L155.2,305.2 L159.8,304.5 L164.3,303.4 L168.9,302.0 L173.4,300.3 L178.0,298.3 L182.6,296.0 L187.1,293.3 L191.7,290.4 L196.2,287.4 L200.8,284.3 L205.4,281.1 L209.9,277.9 L214.5,274.5 L219.0,271.2 L223.6,267.8 L228.2,264.5 L232.7,261.1 L237.3,257.6 L241.8,254.2 L246.4,250.8 L251.0,247.4 L255.5,244.0 L260.1,240.5 L264.6,237.1 L269.2,233.7 L273.8,230.2 L278.3,226.8 L282.9,223.4 L287.4,219.9 L292.0,216.5 L296.6,213.1 L301.1,209.6 L305.7,206.2 L310.2,202.7 L314.8,199.3 L319.4,195.9 L323.9,192.4 L328.5,189.0 L333.0,185.6 L337.6,182.1 L342.2,178.7 L346.7,175.3 L351.3,171.8 L355.8,168.4 L360.4,165.0 L365.0,161.5 L369.5,158.1 L374.1,154.7 L378.6,151.2 L383.2,147.8 L387.8,144.3 L392.3,140.9 L396.9,137.5 L401.4,134.0 L406.0,130.6 L410.6,127.2 L415.1,123.7 L419.7,120.3 L424.2,116.9 L428.8,113.4 L433.4,110.0 L437.9,106.6 L442.5,103.1 L447.0,99.7 L451.6,96.2 L456.2,92.8 L460.7,89.4 L465.3,85.9 L469.8,82.5 L474.4,79.1 L479.0,75.6 L483.5,72.2 L488.1,68.8 L492.6,65.3 L497.2,61.9 L501.8,58.5 L506.3,55.0 L510.9,51.6 L515.4,48.2 L520.0,44.7" fill="none" stroke="#1baf7a" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><text x="530" y="313.0" class="viz-label">r = 0.45</text><text x="530" y="334.8" class="viz-label">r = 0.50</text><text x="530" y="48.7" class="viz-label">r = 0.55</text><line x1="0" x2="0" y1="30" y2="340" class="viz-crosshair" visibility="hidden"/></svg>
<div class="viz-tooltip" hidden></div>
</div>
<figcaption>FTCS for u_t = u_xx with h = 1/25 and a hat-shaped start of height 1, on a log scale. The lines for r = 0.45 and r = 0.50 almost lie on top of each other. The table shows every 25th step. Hover over the chart, or focus it and use the arrow keys, to read the values.</figcaption>
<details class="viz-table"><summary>Show the numbers</summary><table><thead><tr><th>n</th><th>r = 0.45</th><th>r = 0.50</th><th>r = 0.55</th></tr></thead><tbody><tr><td>0</td><td>0.9600</td><td>0.9600</td><td>0.9600</td></tr><tr><td>50</td><td>0.5696</td><td>0.5465</td><td>1.2351</td></tr><tr><td>100</td><td>0.3967</td><td>0.3664</td><td>1.81e3</td></tr><tr><td>150</td><td>0.2778</td><td>0.2466</td><td>3.83e6</td></tr><tr><td>200</td><td>0.1946</td><td>0.1660</td><td>8.08e9</td></tr></tbody></table></details>
<script>(() => {const fig = document.currentScript.closest('figure');const svg = fig.querySelector('svg');const line = svg.querySelector('.viz-crosshair');const tip = fig.querySelector('.viz-tooltip');const names = ["r = 0.45", "r = 0.50", "r = 0.55"]; const colors = ["#2a78d6", "#eb6834", "#1baf7a"];const rows = [["0.9600", "0.9600", "0.9600"], ["0.9042", "0.9000", "0.8962"], ["0.8713", "0.8650", "0.8593"], ["0.8453", "0.8375", "0.8304"], ["0.8232", "0.8141", "0.8059"], ["0.8035", "0.7933", "0.7841"], ["0.7857", "0.7744", "0.7644"], ["0.7692", "0.7570", "0.7463"], ["0.7538", "0.7408", "0.7295"], ["0.7393", "0.7255", "0.7138"], ["0.7256", "0.7110", "0.6989"], ["0.7126", "0.6973", "0.6850"], ["0.7001", "0.6841", "0.6718"], ["0.6881", "0.6714", "0.6593"], ["0.6766", "0.6593", "0.6477"], ["0.6654", "0.6475", "0.6368"], ["0.6546", "0.6361", "0.6269"], ["0.6442", "0.6251", "0.6266"], ["0.6340", "0.6144", "0.6327"], ["0.6241", "0.6040", "0.6448"], ["0.6145", "0.5938", "0.6649"], ["0.6051", "0.5839", "0.7115"], ["0.5959", "0.5742", "0.7865"], ["0.5870", "0.5648", "0.8911"], ["0.5782", "0.5555", "1.0358"], ["0.5696", "0.5465", "1.2351"], ["0.5612", "0.5376", "1.5149"], ["0.5529", "0.5289", "1.9264"], ["0.5448", "0.5204", "2.4875"], ["0.5369", "0.5121", "3.2520"], ["0.5291", "0.5039", "4.2927"], ["0.5214", "0.4958", "5.7083"], ["0.5138", "0.4879", "7.6334"], ["0.5064", "0.4802", "10.3"], ["0.4991", "0.4725", "13.8"], ["0.4920", "0.4650", "18.6"], ["0.4849", "0.4577", "25.2"], ["0.4780", "0.4504", "34.1"], ["0.4711", "0.4433", "46.2"], ["0.4644", "0.4363", "62.7"], ["0.4578", "0.4294", "85.1"], ["0.4512", "0.4226", "115"], ["0.4448", "0.4160", "157"], ["0.4385", "0.4094", "213"], ["0.4323", "0.4030", "289"], ["0.4261", "0.3966", "392"], ["0.4201", "0.3904", "533"], ["0.4141", "0.3842", "724"], ["0.4082", "0.3782", "983"], ["0.4024", "0.3722", "1.33e3"], ["0.3967", "0.3664", "1.81e3"], ["0.3911", "0.3606", "2.46e3"], ["0.3855", "0.3549", "3.34e3"], ["0.3801", "0.3493", "4.54e3"], ["0.3747", "0.3439", "6.17e3"], ["0.3694", "0.3384", "8.38e3"], ["0.3642", "0.3331", "1.14e4"], ["0.3590", "0.3279", "1.55e4"], ["0.3539", "0.3227", "2.1e4"], ["0.3489", "0.3177", "2.85e4"], ["0.3440", "0.3127", "3.87e4"], ["0.3391", "0.3078", "5.26e4"], ["0.3343", "0.3029", "7.15e4"], ["0.3296", "0.2982", "9.71e4"], ["0.3249", "0.2935", "1.32e5"], ["0.3203", "0.2889", "1.79e5"], ["0.3158", "0.2843", "2.43e5"], ["0.3113", "0.2799", "3.3e5"], ["0.3069", "0.2755", "4.49e5"], ["0.3026", "0.2711", "6.09e5"], ["0.2983", "0.2669", "8.28e5"], ["0.2941", "0.2627", "1.12e6"], ["0.2899", "0.2586", "1.53e6"], ["0.2858", "0.2545", "2.07e6"], ["0.2818", "0.2505", "2.82e6"], ["0.2778", "0.2466", "3.83e6"], ["0.2738", "0.2427", "5.2e6"], ["0.2700", "0.2389", "7.06e6"], ["0.2662", "0.2351", "9.59e6"], ["0.2624", "0.2314", "1.3e7"], ["0.2587", "0.2278", "1.77e7"], ["0.2550", "0.2242", "2.4e7"], ["0.2514", "0.2207", "3.26e7"], ["0.2479", "0.2172", "4.43e7"], ["0.2443", "0.2138", "6.02e7"], ["0.2409", "0.2105", "8.18e7"], ["0.2375", "0.2072", "1.11e8"], ["0.2341", "0.2039", "1.51e8"], ["0.2308", "0.2007", "2.05e8"], ["0.2276", "0.1975", "2.78e8"], ["0.2243", "0.1944", "3.78e8"], ["0.2212", "0.1914", "5.13e8"], ["0.2180", "0.1884", "6.97e8"], ["0.2149", "0.1854", "9.47e8"], ["0.2119", "0.1825", "1.29e9"], ["0.2089", "0.1796", "1.75e9"], ["0.2060", "0.1768", "2.37e9"], ["0.2030", "0.1740", "3.22e9"], ["0.2002", "0.1713", "4.38e9"], ["0.1973", "0.1686", "5.95e9"], ["0.1946", "0.1660", "8.08e9"]]; const heads = ["n = 0", "n = 2", "n = 4", "n = 6", "n = 8", "n = 10", "n = 12", "n = 14", "n = 16", "n = 18", "n = 20", "n = 22", "n = 24", "n = 26", "n = 28", "n = 30", "n = 32", "n = 34", "n = 36", "n = 38", "n = 40", "n = 42", "n = 44", "n = 46", "n = 48", "n = 50", "n = 52", "n = 54", "n = 56", "n = 58", "n = 60", "n = 62", "n = 64", "n = 66", "n = 68", "n = 70", "n = 72", "n = 74", "n = 76", "n = 78", "n = 80", "n = 82", "n = 84", "n = 86", "n = 88", "n = 90", "n = 92", "n = 94", "n = 96", "n = 98", "n = 100", "n = 102", "n = 104", "n = 106", "n = 108", "n = 110", "n = 112", "n = 114", "n = 116", "n = 118", "n = 120", "n = 122", "n = 124", "n = 126", "n = 128", "n = 130", "n = 132", "n = 134", "n = 136", "n = 138", "n = 140", "n = 142", "n = 144", "n = 146", "n = 148", "n = 150", "n = 152", "n = 154", "n = 156", "n = 158", "n = 160", "n = 162", "n = 164", "n = 166", "n = 168", "n = 170", "n = 172", "n = 174", "n = 176", "n = 178", "n = 180", "n = 182", "n = 184", "n = 186", "n = 188", "n = 190", "n = 192", "n = 194", "n = 196", "n = 198", "n = 200"];const xs = [64.0, 68.6, 73.1, 77.7, 82.2, 86.8, 91.4, 95.9, 100.5, 105.0, 109.6, 114.2, 118.7, 123.3, 127.8, 132.4, 137.0, 141.5, 146.1, 150.6, 155.2, 159.8, 164.3, 168.9, 173.4, 178.0, 182.6, 187.1, 191.7, 196.2, 200.8, 205.4, 209.9, 214.5, 219.0, 223.6, 228.2, 232.7, 237.3, 241.8, 246.4, 251.0, 255.5, 260.1, 264.6, 269.2, 273.8, 278.3, 282.9, 287.4, 292.0, 296.6, 301.1, 305.7, 310.2, 314.8, 319.4, 323.9, 328.5, 333.0, 337.6, 342.2, 346.7, 351.3, 355.8, 360.4, 365.0, 369.5, 374.1, 378.6, 383.2, 387.8, 392.3, 396.9, 401.4, 406.0, 410.6, 415.1, 419.7, 424.2, 428.8, 433.4, 437.9, 442.5, 447.0, 451.6, 456.2, 460.7, 465.3, 469.8, 474.4, 479.0, 483.5, 488.1, 492.6, 497.2, 501.8, 506.3, 510.9, 515.4, 520.0]; const steps = 100, width = 680;let current = -1;const show = (n) => { n = Math.max(0, Math.min(steps, n)); current = n; const x = xs[n]; line.setAttribute('x1', x); line.setAttribute('x2', x); line.setAttribute('visibility', 'visible'); tip.replaceChildren(); const head = document.createElement('div'); head.className = 'viz-tooltip-head'; head.textContent = heads[n]; tip.append(head); names.forEach((name, i) => {  const row = document.createElement('div'); row.className = 'viz-tooltip-row';  const key = document.createElement('span'); key.className = 'viz-swatch'; key.style.background = colors[i];  const value = document.createElement('strong'); value.textContent = rows[n][i];  const label = document.createElement('span'); label.textContent = name;  row.append(key, value, label); tip.append(row); }); const box = svg.getBoundingClientRect(), frame = fig.getBoundingClientRect(); const xp = box.left - frame.left + x * box.width / width; tip.hidden = false; const room = frame.width - tip.offsetWidth - 8; tip.style.left = Math.max(0, Math.min(room, xp + 12)) + 'px'; tip.style.top = (box.top - frame.top + 8) + 'px';};const hide = () => { line.setAttribute('visibility', 'hidden'); tip.hidden = true; current = -1; };svg.addEventListener('pointermove', (event) => { const box = svg.getBoundingClientRect(); const x = (event.clientX - box.left) * width / box.width; let best = 0; xs.forEach((v, i) => { if (Math.abs(v - x) < Math.abs(xs[best] - x)) best = i; }); show(best);});svg.addEventListener('pointerleave', hide);svg.addEventListener('focus', () => show(current < 0 ? 0 : current));svg.addEventListener('blur', hide);svg.addEventListener('keydown', (event) => { if (event.key === 'ArrowRight') { show(current + 1); event.preventDefault(); } if (event.key === 'ArrowLeft') { show(current - 1); event.preventDefault(); } if (event.key === 'Escape') hide();});})();</script>
</figure>

| $r$ | $k$ | $\max\lvert U\rvert$ at $n = 50$ | $n = 100$ | $n = 200$ | growth per step (last 50 steps) | largest $\lvert 1 + k\lambda_p\rvert$ (odd $p$) |
|---|---|---|---|---|---|---|
| 0.45 | $0.00072$ | $0.570$ | $0.397$ | $0.195$ | $0.9929$ | $0.9929$ |
| 0.50 | $0.00080$ | $0.547$ | $0.366$ | $0.166$ | $0.9921$ | $0.9921$ |
| 0.55 | $0.00088$ | $1.24$ | $1.81 \times 10^{3}$ | $8.08 \times 10^{9}$ | $1.1654$ | $1.1654$ |

- $r = 0.45$ and $r = 0.50$ decay smoothly, as heat should. At the edge $r = 1/2$, nothing goes wrong.
- Just past the edge, at $r = 0.55$, the solution explodes. It grows by $10^{10}$ in $200$ steps, and the shape is a sawtooth.
- The measured growth per step matches the largest factor $|1 + k\lambda_p|$ to four digits. One detail: the hat is symmetric about $x = 1/2$, so it contains only the odd modes $p = 1, 3, \dots, 23$. The largest factor over **all** modes is $1.1913$ (from $p = 24$), but that mode is not in the start, so the right comparison is over the odd modes. (Rounding errors would excite $p = 24$ too, but from about $10^{-16}$ it does not catch up within $200$ steps.)
- For $r \le 1/2$, the slowest decay per step is $|1 + k\lambda_1| \approx 1 - \pi^2 k$: $0.9929$ and $0.9921$. That is the smooth mode $\sin(\pi x)$, which is what remains at the end.

**Part 4: Crank–Nicolson with large steps** (my own check). I use $h = 1/50$, the same hat start, and the final time $T = 0.1$. The exact solution is the Fourier series of the hat. I compare Crank–Nicolson with Backward Euler and take the largest error over the grid points. $g$ is the factor per step for the sawtooth mode ($p = N$).

| steps | $k$ | $r$ | error, CN | error, BE | $g$ for CN | $g$ for BE |
|---|---|---|---|---|---|---|
| 1000 | $0.0001$ | $0.25$ | $1.98 \times 10^{-4}$ | $3.45 \times 10^{-4}$ | $+0.334$ | $0.500$ |
| 100 | $0.001$ | $2.5$ | $1.95 \times 10^{-4}$ | $1.67 \times 10^{-3}$ | $-0.666$ | $0.091$ |
| 20 | $0.005$ | $12.5$ | $1.93 \times 10^{-3}$ | $7.45 \times 10^{-3}$ | $-0.923$ | $0.020$ |
| 10 | $0.01$ | $25$ | $1.19 \times 10^{-2}$ | $1.45 \times 10^{-2}$ | $-0.961$ | $0.010$ |
| 5 | $0.02$ | $50$ | $3.29 \times 10^{-2}$ | $2.79 \times 10^{-2}$ | $-0.980$ | $0.005$ |
| 2 | $0.05$ | $125$ | $8.07 \times 10^{-2}$ | $6.48 \times 10^{-2}$ | $-0.992$ | $0.002$ |

- Both schemes stay bounded for every $k$, as the theory says. The exact peak at $T = 0.1$ is $0.302$.
- With small steps, the error of Crank–Nicolson stops near $2 \times 10^{-4}$. That is the space error from $h = 1/50$, not the time error.
- With large steps, Crank–Nicolson loses its advantage, and from $5$ steps on it is even **worse** than the first-order Backward Euler. The hat has a corner, so it contains many fast modes. Backward Euler kills them ($g \approx 0$). Crank–Nicolson keeps them alive with $g$ close to $-1$, so they flip sign at every step instead of dying.
- So "unconditionally stable" is about boundedness, not accuracy. This is one reason to keep $k = O(h)$ with Crank–Nicolson.

The full script and the chart code are in [this site's repository](https://github.com/doubleyands/doubleyands.github.io/tree/main/scripts/posts/numerical-pde/lecture-10).

## Summary

1. The second-difference matrix has eigenvectors $\sin(p\pi jh)$ and eigenvalues $\lambda_p = \frac{2}{h^2}(\cos(p\pi h) - 1)$, from about $-\pi^2$ to about $-4/h^2$. We guess them from the continuous problem and check them with a trig identity.
2. The **condition number** $\kappa = \|A\|\,\|A^{-1}\|$ measures how much a system amplifies relative errors. For our matrix, $\|A^{-1}\|$ stays bounded, but $\kappa \approx 4/(\pi^2h^2)$ grows.
3. With the method of lines, $U' = AU$ splits into test equations $w_p' = \lambda_p w_p$. Putting every $k\lambda_p$ into the Forward Euler region gives the FTCS limit $k/h^2 \le 1/2$. Past it, a sawtooth grows.
4. Crank–Nicolson (and Backward Euler) are **unconditionally stable**, because every $k\lambda_p$ is negative. Accuracy still asks for $k = O(h)$.
5. **Lax–Richtmyer stability**: $\|B(k)^n\| \le M$ for $kn \le T$. **Lax equivalence theorem**: for a consistent linear method, convergence and stability are the same thing.
6. When the eigenvalues are out of reach, **von Neumann analysis** uses Fourier modes $e^{ijh\xi}$ instead. It needs constant coefficients and no boundary (or periodic boundaries).
