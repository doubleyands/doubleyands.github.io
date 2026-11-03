---
title: "Numerical PDEs, Lecture 9: Stability Regions, the Root Condition, and the Heat Equation"
description: "My notes on the ninth lecture: why the stability parameter z can be complex, the stability regions of Forward Euler, Backward Euler and the trapezoidal method, the stability polynomial and the root condition for multistep methods, the boundary locus method, and a first look at the heat equation with FTCS and Crank–Nicolson."
pubDate: 2026-11-03
project: numerical-pde
tags: [numerical PDEs, ODEs, time stepping, stability, heat equation, notes]
---

This is my summary of the ninth lecture in my numerical PDE course. I write the ideas in my own words and check them with a small experiment.

[Lecture 8](/blog/2026-10-29-numerical-pde-lecture-8/) split stability into two ideas. **Zero-stability** asks: does the method converge as $\Delta t \to 0$? **Absolute stability** asks: is the method useful for **this** $\Delta t$? In practice we want the largest step that is still safe, because a smaller step costs more. This lecture turns that question into a picture: a **stability region** in the complex plane. It does this for one-step methods and for multistep methods. In the last part, we put space and time together for the first time and meet the **heat equation**.

## 1. Recap: the test equation

Absolute stability uses the test equation $u' = \lambda u$. For Forward Euler,

$$
U^{n+1} = (1 + \lambda\Delta t)\,U^n = (1 + z)\,U^n, \qquad z = \lambda\Delta t.
$$

The step $\Delta t$ and the rate $\lambda$ only appear together, in the product $z$. So we first find which values of $z$ are safe, and then we turn that into a condition on $\Delta t$. For real $z$, Forward Euler needs $-2 \le z \le 0$ (Lecture 8).

## 2. Why $z$ can be complex

Real problems are usually systems. Take $\mathbf{u}' = A\mathbf{u}$ with a matrix $A$. If $A$ can be diagonalized, $A = V\Lambda V^{-1}$, then the new variables $\mathbf{w} = V^{-1}\mathbf{u}$ satisfy

$$
w_j' = \lambda_j\, w_j, \qquad j = 1, 2, \dots
$$

The system breaks into separate test equations, one for each eigenvalue $\lambda_j$ of $A$. A linear method does the same split, so it is safe when **every** $z_j = \lambda_j\Delta t$ is safe. Eigenvalues are roots of $\det(A - \lambda I) = 0$, and these roots can be complex. This is why we write $z$, the usual letter for a complex number.

A simple example is the oscillator $x' = y$, $y' = -x$. Its matrix $\begin{bmatrix}0 & 1\\ -1 & 0\end{bmatrix}$ has eigenvalues $\pm i$, so $z = \pm i\,\Delta t$ lies on the imaginary axis. The solution goes around a circle, and its size $\sqrt{x^2 + y^2}$ never changes.

| Problem | What is $\lambda$? |
|---|---|
| linear, $\mathbf{u}' = A\mathbf{u}$ | the eigenvalues of $A$ |
| nonlinear, $\mathbf{u}' = f(\mathbf{u})$ | first linearize (like a tangent line), then take the eigenvalues of the Jacobian matrix |

## 3. Stability regions of three one-step methods

The **region of absolute stability** is the set of complex $z$ where the method does not let errors grow. For a one-step method, the test equation gives $U^{n+1} = P(z)\,U^n$, and the region is $\{z : |P(z)| \le 1\}$. Write $z = x + iy$. Then $|z| = \sqrt{x^2 + y^2}$ is the distance from $z$ to $0$.

| Method | $P(z)$ | Condition | Region |
|---|---|---|---|
| Forward Euler | $1 + z$ | $(1 + x)^2 + y^2 \le 1$ | **inside** the circle with center $-1$ and radius $1$ |
| Backward Euler | $\dfrac{1}{1 - z}$ | $\lvert 1 - z\rvert \ge 1$ | **outside** the circle with center $+1$ and radius $1$ |
| Trapezoid | $\dfrac{1 + z/2}{1 - z/2}$ | $x \le 0$ | the whole **left half plane** |

(In the Forward Euler row, $|1 + z| \le 1$ has a square root on the left. Squaring both sides gives the same condition, so the table shows the squared form.)

The trapezoid row needs a short computation. The condition is $|1 + z/2| \le |1 - z/2|$. Square both sides and subtract:

$$
\Big(1 + \tfrac{x}{2}\Big)^2 + \tfrac{y^2}{4} - \Big(1 - \tfrac{x}{2}\Big)^2 - \tfrac{y^2}{4} = 2x \le 0.
$$

The $y$ terms cancel completely, and only $x \le 0$ is left.

**Which region is better?** A bigger region gives more freedom to choose $\Delta t$. Backward Euler and the trapezoidal method cover the whole left half plane. So for a decaying problem ($\operatorname{Re}\lambda < 0$), **any** step is safe. Forward Euler has only a small disk. But Forward Euler is still a good method when we know the safe step in advance, as in the Lecture 8 example.

**A caution.** "Bigger is always better" is not quite right. Backward Euler is also stable in most of the **right** half plane. There the exact solution **grows**, but the numerical solution shrinks. So the method hides real growth. (This is my own reading. The experiment shows a case.)

## 4. One-step methods: $P(z)$ should be close to $e^z$

Why should $P(z)$ look like an exponential? Start the exact solution of $u' = \lambda u$ at $u(t_n) = U^n$. Then $u(t) = e^{\lambda(t - t_n)}\,U^n$, so after one step

$$
u(t_{n+1}) = e^{\lambda\Delta t}\,U^n = e^{z}\,U^n.
$$

The exact one-step factor is $e^z$. A consistent method must approximate it.

| Method | $P(z)$ | Type |
|---|---|---|
| Forward Euler | $1 + z$ | explicit: polynomial |
| RK2 | $1 + z + \frac{z^2}{2}$ | explicit: polynomial |
| RK4 | $1 + z + \frac{z^2}{2} + \frac{z^3}{6} + \frac{z^4}{24}$ | explicit: polynomial |
| Backward Euler | $\dfrac{1}{1 - z}$ | implicit: rational function |
| Trapezoid | $\dfrac{1 + z/2}{1 - z/2}$ | implicit: rational function |

For an explicit method, $P$ is a **polynomial**, and it matches the Taylor series of $e^z$ up to the order of the method. This is the same fact as in [Lecture 7](/blog/2026-10-27-numerical-pde-lecture-7/): Runge–Kutta methods copy the Taylor series without derivatives. For an implicit method, we move the unknown $U^{n+1}$ to the left side and divide, so $P$ becomes a **ratio** of polynomials. Approximating $e^z$ by a ratio of polynomials is called **Padé approximation**. Backward Euler and the trapezoidal method are two simple cases.

**How to draw the region.** For a one-step method we do not need to solve anything. Compute $|P(z)|$ on a grid of $z$ values and draw the contour line where $|P(z)| = 1$.

## 5. Multistep methods: the stability polynomial

It helps to count **time levels**. A one-step method uses two levels, $n$ and $n + 1$. A two-step method like leapfrog uses three.

Put $f = \lambda u$ into the linear multistep method of Lecture 7:

$$
\sum_{j=0}^{r}\alpha_j U^{n+j} = \Delta t\sum_{j=0}^{r}\beta_j\,\lambda U^{n+j}
\quad\Longrightarrow\quad
\sum_{j=0}^{r}(\alpha_j - z\beta_j)\,U^{n+j} = 0.
$$

Now there is no single factor $P(z)$. Instead we get the **stability polynomial**

$$
\pi(\eta; z) = \sum_{j=0}^{r}(\alpha_j - z\beta_j)\,\eta^j = \rho(\eta) - z\,\sigma(\eta),
$$

with $\rho(\eta) = \sum_j \alpha_j\eta^j$ and $\sigma(\eta) = \sum_j \beta_j\eta^j$. To see what its roots mean, look at a famous recurrence.

## 6. Fibonacci and the root condition

Take the Fibonacci rule

$$
U^{n+1} = U^n + U^{n-1}.
$$

**Way 1: guess the form of the solution.** For an ODE we guess $e^{\lambda t}$. For a difference equation we guess $U^n = \eta^n$:

$$
\eta^{n+1} = \eta^n + \eta^{n-1} \;\Longrightarrow\; \eta^2 = \eta + 1 \;\Longrightarrow\; \eta = \frac{1 \pm \sqrt5}{2}.
$$

So $U^n = C_1\eta_1^n + C_2\eta_2^n$. The two constants come from the two starting values $U^0$ and $U^1$. This is again the fact from Lecture 7: a multistep method cannot start by itself.

**Way 2: make it a one-step system.**

$$
\begin{bmatrix}U^{n+1}\\ U^{n}\end{bmatrix} =
\begin{bmatrix}1 & 1\\ 1 & 0\end{bmatrix}
\begin{bmatrix}U^{n}\\ U^{n-1}\end{bmatrix},
\qquad
\det\begin{bmatrix}1 - \eta & 1\\ 1 & -\eta\end{bmatrix} = \eta^2 - \eta - 1 = 0.
$$

The same equation appears, and this is not luck. We write a higher-order ODE as a first-order system in the same way. Here a multistep rule becomes a one-step rule for a vector, and the eigenvalues of that matrix are the roots of the characteristic equation.

**The root condition.** A multistep method is absolutely stable at $z$ when all roots $\eta_i$ of $\pi(\eta; z)$ satisfy

$$
|\eta_i| \le 1, \qquad \text{and } |\eta_i| < 1 \text{ if } \eta_i \text{ is a repeated root}.
$$

Why the strict condition for repeated roots? For an ODE, a double root of the characteristic equation gives solutions like $t\,e^{\lambda t}$. A difference equation is the same: a double root gives $n\,\eta^n$. If $|\eta| = 1$, this grows **linearly** in $n$. The experiment shows exactly this.

**Fibonacci is unstable leapfrog.** Leapfrog is $U^{n+1} = U^{n-1} + 2\Delta t\, f(U^n)$. With $f = \lambda u$ it becomes $U^{n+1} = U^{n-1} + 2zU^n$, and its stability polynomial is

$$
\pi(\eta; z) = \eta^2 - 2z\eta - 1.
$$

With $z = 1/2$ this is $\eta^2 - \eta - 1$, the Fibonacci equation. One root is $1.618 > 1$, so this is unstable. (Lecture 7 also met Fibonacci as leapfrog with $\Delta t = 1/2$ and $f(u) = u$.)

**Why the leapfrog region is so small.** The product of the two roots is the constant term, $-1$, so $|\eta_1|\,|\eta_2| = 1$. Both sizes are at most $1$ only if **both are exactly $1$**. Then the two roots are $e^{i\theta}$ and $-e^{-i\theta}$, and their sum $2z = 2i\sin\theta$ is purely imaginary. At $z = \pm i$ the two roots meet ($\eta = \pm i$, a double root on the circle), and the root condition fails. So the region of leapfrog is only the open segment

$$
z \in (-i, i)
$$

on the imaginary axis. **No real $z$ except $0$ is stable.** Even a tiny decay rate makes leapfrog grow in the long run. This will matter soon, because the heat equation gives negative real eigenvalues.

## 7. The boundary locus method

For Forward Euler, $\pi = \eta - (1 + z)$, and for Backward Euler, $\pi = (1 - z)\eta - 1$. We can find their roots by hand. With four time levels, $\pi$ is a cubic. With more levels, we need a computer. Drawing a region by solving for roots at every point would be slow, so here is a better idea.

**At the boundary of the region, one root lies on the unit circle.** So $\eta = e^{i\theta}$ for some angle $\theta$. Put it into $\pi = 0$ and solve for $z$:

$$
\rho(e^{i\theta}) - z\,\sigma(e^{i\theta}) = 0
\quad\Longrightarrow\quad
z(\theta) = \frac{\rho(e^{i\theta})}{\sigma(e^{i\theta})}, \qquad 0 \le \theta \le 2\pi.
$$

Move $\theta$ around the circle and plot the points $z(\theta)$. This draws the candidate boundary curve. A modest number of angles already shows the shape, and more angles make it smooth.

**Check.** Forward Euler has $\rho(\eta) = \eta - 1$ and $\sigma(\eta) = 1$. So $z(\theta) = e^{i\theta} - 1$, the circle with center $-1$ and radius $1$, as in section 3.

**Inside or outside?** The curve alone does not tell us. Forward Euler is stable inside its circle, but Backward Euler is stable outside. So pick one test point on each side of the curve and check the root condition there.

**Faster, but more fragile.** A multistep method that is started correctly can give high order cheaply. But for the Adams families, the stability region gets smaller as the order goes up. Then the step must be small enough for all $z$ to fit inside, and we may not know if they do. The implicit Adams methods have bigger regions than the explicit ones with the same number of steps, but they also shrink with the order.

## 8. The heat equation: space and time together

$$
u_t = \kappa\,u_{xx}\ \ (\kappa = 1),\qquad u(x, 0) = \eta(x),\qquad u(0, t) = g_0(t),\ \ u(1, t) = g_1(t).
$$

This is also called the diffusion equation, a **parabolic** equation. We already have the two tools we need:

- **space**: centered differences ([Lecture 3](/blog/2026-10-13-numerical-pde-lecture-3/) and [Lecture 4](/blog/2026-10-15-numerical-pde-lecture-4/));
- **time**: Forward and Backward Euler, the trapezoidal method, RK4, Adams methods ([Lecture 6](/blog/2026-10-22-numerical-pde-lecture-6/) and Lecture 7).

Designing a scheme is easy: combine the two. The hard question is whether it is **stable**.

From here on I use the notation $U_i^n \approx u(x_i, t_n)$, with $i$ for space and $n$ for time, $h = \Delta x$ and $k = \Delta t$.

**FTCS.** The simplest scheme is

$$
\frac{U_i^{n+1} - U_i^n}{k} = \frac{U_{i-1}^n - 2U_i^n + U_{i+1}^n}{h^2}.
$$

**How to read a scheme.** First find the reference point and the reference time. The space indices $i - 1, i, i + 1$ are centered at $x_i$, so space uses a **centered** difference. The right side uses only time level $n$, so time uses **Forward Euler**. That is the name: **F**orward **T**ime, **C**entered **S**pace.

## 9. Crank–Nicolson: the trapezoid in time

For $u_t = F$, three time schemes sit side by side:

| $\dfrac{U_i^{n+1} - U_i^n}{k} =$ | Name |
|---|---|
| $F_i^n$ | Forward Euler |
| $F_i^{n+1}$ | Backward Euler |
| $\tfrac12\big(F_i^n + F_i^{n+1}\big)$ | **Crank–Nicolson** |

Crank–Nicolson is just the trapezoidal method. The only difference is that there is now a space variable. For the heat equation, put level $n + 1$ on the left and level $n$ on the right:

$$
-rU_{i-1}^{n+1} + (1 + 2r)\,U_i^{n+1} - rU_{i+1}^{n+1} = rU_{i-1}^n + (1 - 2r)\,U_i^n + rU_{i+1}^n,
\qquad r = \frac{k}{2h^2}.
$$

Two quick checks: multiply out $\tfrac{k}{2h^2}$ times the six terms and you get this line. Also, the coefficients on each side add up to $1$, so a constant solution stays constant (away from the boundaries).

**As a matrix.** The left side is a **tridiagonal** matrix with $1 + 2r$ on the diagonal and $-r$ next to it. The boundary values move to the right side. So each time step solves a system with the same structure as the steady problem of Lecture 2. This is the price of an implicit method. In 2D the matrix has five diagonals, and in 3D it gets complicated.

**Accuracy.** Time uses the trapezoidal method (second order) and space uses centered differences (second order), so the truncation error is $O(k^2 + h^2)$.

## 10. The truncation error of FTCS

Put the exact solution into FTCS and expand with Taylor series:

$$
\tau = \Big(u_t + \tfrac12 k\,u_{tt} + \cdots\Big) - \Big(u_{xx} + \tfrac{1}{12}h^2 u_{xxxx} + \cdots\Big).
$$

The terms $u_t$ and $u_{xx}$ cancel because $u_t = u_{xx}$. Use the equation once more: $u_{tt} = (u_t)_t = (u_{xx})_t = (u_t)_{xx} = u_{xxxx}$. So

$$
\tau = \Big(\tfrac12 k - \tfrac{1}{12}h^2\Big)\,u_{xxxx} + O(k^2 + h^4) = O(k + h^2).
$$

**First order in time, second order in space.** First order in time is natural, because time uses Forward Euler. The formula also says something less obvious: if $k = h^2/6$, the leading term vanishes. The experiment checks this.

Is FTCS stable, and for which $k$? That is the topic of the next lecture, and the stability regions from this lecture are the tool.

## 11. Experiment

The script uses the Python standard library only.

**Part 1: complex eigenvalues.** I solve the oscillator $x' = y$, $y' = -x$ with $x(0) = 1$, $y(0) = 0$ and $\Delta t = 0.1$, for ten turns (628 steps, up to $t = 62.8$). The exact size $\sqrt{x^2 + y^2}$ is always $1$. Here $z = \pm 0.1i$, right on the imaginary axis. Leapfrog takes its second value from the exact solution, so only the method itself is tested.

<figure class="viz">
<div class="viz-legend"><span class="viz-key"><span class="viz-swatch" style="background:#2a78d6"></span>Forward Euler</span><span class="viz-key"><span class="viz-swatch" style="background:#eb6834"></span>Backward Euler</span><span class="viz-key"><span class="viz-swatch" style="background:#1baf7a"></span>Trapezoid</span><span class="viz-key"><span class="viz-swatch" style="background:#eda100"></span>Leapfrog</span></div>
<div class="viz-frame">
<svg viewBox="0 0 680 400" role="img" tabindex="0" aria-label="Amplitude of the numerical solution of the oscillator over ten turns, with time step 0.1, on a log scale. Forward Euler grows to about 23, Backward Euler shrinks to about 0.04, and the trapezoid and leapfrog methods stay at 1." class="viz-plot"><line x1="64" x2="520" y1="79.3" y2="79.3" class="viz-grid"/><text x="54" y="83.3" text-anchor="end" class="viz-tick">10</text><line x1="64" x2="520" y1="133.3" y2="133.3" class="viz-grid"/><text x="54" y="137.3" text-anchor="end" class="viz-tick">3</text><line x1="64" x2="520" y1="182.6" y2="182.6" class="viz-grid"/><text x="54" y="186.6" text-anchor="end" class="viz-tick">1</text><line x1="64" x2="520" y1="236.7" y2="236.7" class="viz-grid"/><text x="54" y="240.7" text-anchor="end" class="viz-tick">0.3</text><line x1="64" x2="520" y1="286.0" y2="286.0" class="viz-grid"/><text x="54" y="290.0" text-anchor="end" class="viz-tick">0.1</text><text x="64.0" y="362" text-anchor="middle" class="viz-tick">0</text><text x="206.5" y="362" text-anchor="middle" class="viz-tick">20</text><text x="349.0" y="362" text-anchor="middle" class="viz-tick">40</text><text x="491.5" y="362" text-anchor="middle" class="viz-tick">60</text><text x="292.0" y="388" text-anchor="middle" class="viz-axis-title">Time t (one turn = 2π ≈ 6.28)</text><text x="64" y="14" class="viz-axis-title">Amplitude √(x² + y²), exact value 1</text><path d="M64.0,182.6 L66.8,181.7 L69.7,180.8 L72.5,180.0 L75.4,179.1 L78.2,178.2 L81.1,177.3 L84.0,176.4 L86.8,175.5 L89.7,174.6 L92.5,173.7 L95.3,172.8 L98.2,171.9 L101.1,171.0 L103.9,170.1 L106.8,169.2 L109.6,168.3 L112.5,167.5 L115.3,166.6 L118.2,165.7 L121.0,164.8 L123.8,163.9 L126.7,163.0 L129.6,162.1 L132.4,161.2 L135.2,160.3 L138.1,159.4 L140.9,158.5 L143.8,157.6 L146.7,156.7 L149.5,155.8 L152.4,155.0 L155.2,154.1 L158.1,153.2 L160.9,152.3 L163.8,151.4 L166.6,150.5 L169.4,149.6 L172.3,148.7 L175.2,147.8 L178.0,146.9 L180.9,146.0 L183.7,145.1 L186.6,144.2 L189.4,143.3 L192.2,142.4 L195.1,141.6 L198.0,140.7 L200.8,139.8 L203.7,138.9 L206.5,138.0 L209.4,137.1 L212.2,136.2 L215.1,135.3 L217.9,134.4 L220.8,133.5 L223.6,132.6 L226.5,131.7 L229.3,130.8 L232.2,129.9 L235.0,129.1 L237.9,128.2 L240.7,127.3 L243.6,126.4 L246.4,125.5 L249.2,124.6 L252.1,123.7 L255.0,122.8 L257.8,121.9 L260.6,121.0 L263.5,120.1 L266.4,119.2 L269.2,118.3 L272.1,117.4 L274.9,116.5 L277.8,115.7 L280.6,114.8 L283.5,113.9 L286.3,113.0 L289.1,112.1 L292.0,111.2 L294.9,110.3 L297.7,109.4 L300.6,108.5 L303.4,107.6 L306.2,106.7 L309.1,105.8 L312.0,104.9 L314.8,104.0 L317.6,103.2 L320.5,102.3 L323.3,101.4 L326.2,100.5 L329.1,99.6 L331.9,98.7 L334.8,97.8 L337.6,96.9 L340.5,96.0 L343.3,95.1 L346.2,94.2 L349.0,93.3 L351.9,92.4 L354.7,91.5 L357.6,90.6 L360.4,89.8 L363.2,88.9 L366.1,88.0 L369.0,87.1 L371.8,86.2 L374.7,85.3 L377.5,84.4 L380.4,83.5 L383.2,82.6 L386.1,81.7 L388.9,80.8 L391.8,79.9 L394.6,79.0 L397.5,78.1 L400.3,77.3 L403.2,76.4 L406.0,75.5 L408.9,74.6 L411.7,73.7 L414.6,72.8 L417.4,71.9 L420.2,71.0 L423.1,70.1 L426.0,69.2 L428.8,68.3 L431.7,67.4 L434.5,66.5 L437.4,65.6 L440.2,64.7 L443.1,63.9 L445.9,63.0 L448.8,62.1 L451.6,61.2 L454.5,60.3 L457.3,59.4 L460.2,58.5 L463.0,57.6 L465.9,56.7 L468.7,55.8 L471.6,54.9 L474.4,54.0 L477.2,53.1 L480.1,52.2 L483.0,51.4 L485.8,50.5 L488.7,49.6 L491.5,48.7 L494.4,47.8 L497.2,46.9 L500.1,46.0 L502.9,45.1 L505.8,44.2 L508.6,43.3 L511.5,42.4" fill="none" stroke="#2a78d6" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><path d="M64.0,182.6 L66.8,183.5 L69.7,184.4 L72.5,185.3 L75.4,186.2 L78.2,187.1 L81.1,188.0 L84.0,188.9 L86.8,189.8 L89.7,190.7 L92.5,191.6 L95.3,192.5 L98.2,193.4 L101.1,194.2 L103.9,195.1 L106.8,196.0 L109.6,196.9 L112.5,197.8 L115.3,198.7 L118.2,199.6 L121.0,200.5 L123.8,201.4 L126.7,202.3 L129.6,203.2 L132.4,204.1 L135.2,205.0 L138.1,205.9 L140.9,206.7 L143.8,207.6 L146.7,208.5 L149.5,209.4 L152.4,210.3 L155.2,211.2 L158.1,212.1 L160.9,213.0 L163.8,213.9 L166.6,214.8 L169.4,215.7 L172.3,216.6 L175.2,217.5 L178.0,218.4 L180.9,219.3 L183.7,220.1 L186.6,221.0 L189.4,221.9 L192.2,222.8 L195.1,223.7 L198.0,224.6 L200.8,225.5 L203.7,226.4 L206.5,227.3 L209.4,228.2 L212.2,229.1 L215.1,230.0 L217.9,230.9 L220.8,231.8 L223.6,232.6 L226.5,233.5 L229.3,234.4 L232.2,235.3 L235.0,236.2 L237.9,237.1 L240.7,238.0 L243.6,238.9 L246.4,239.8 L249.2,240.7 L252.1,241.6 L255.0,242.5 L257.8,243.4 L260.6,244.3 L263.5,245.2 L266.4,246.0 L269.2,246.9 L272.1,247.8 L274.9,248.7 L277.8,249.6 L280.6,250.5 L283.5,251.4 L286.3,252.3 L289.1,253.2 L292.0,254.1 L294.9,255.0 L297.7,255.9 L300.6,256.8 L303.4,257.7 L306.2,258.5 L309.1,259.4 L312.0,260.3 L314.8,261.2 L317.6,262.1 L320.5,263.0 L323.3,263.9 L326.2,264.8 L329.1,265.7 L331.9,266.6 L334.8,267.5 L337.6,268.4 L340.5,269.3 L343.3,270.2 L346.2,271.1 L349.0,271.9 L351.9,272.8 L354.7,273.7 L357.6,274.6 L360.4,275.5 L363.2,276.4 L366.1,277.3 L369.0,278.2 L371.8,279.1 L374.7,280.0 L377.5,280.9 L380.4,281.8 L383.2,282.7 L386.1,283.6 L388.9,284.4 L391.8,285.3 L394.6,286.2 L397.5,287.1 L400.3,288.0 L403.2,288.9 L406.0,289.8 L408.9,290.7 L411.7,291.6 L414.6,292.5 L417.4,293.4 L420.2,294.3 L423.1,295.2 L426.0,296.1 L428.8,297.0 L431.7,297.8 L434.5,298.7 L437.4,299.6 L440.2,300.5 L443.1,301.4 L445.9,302.3 L448.8,303.2 L451.6,304.1 L454.5,305.0 L457.3,305.9 L460.2,306.8 L463.0,307.7 L465.9,308.6 L468.7,309.5 L471.6,310.3 L474.4,311.2 L477.2,312.1 L480.1,313.0 L483.0,313.9 L485.8,314.8 L488.7,315.7 L491.5,316.6 L494.4,317.5 L497.2,318.4 L500.1,319.3 L502.9,320.2 L505.8,321.1 L508.6,322.0 L511.5,322.9" fill="none" stroke="#eb6834" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><path d="M64.0,182.6 L66.8,182.6 L69.7,182.6 L72.5,182.6 L75.4,182.6 L78.2,182.6 L81.1,182.6 L84.0,182.6 L86.8,182.6 L89.7,182.6 L92.5,182.6 L95.3,182.6 L98.2,182.6 L101.1,182.6 L103.9,182.6 L106.8,182.6 L109.6,182.6 L112.5,182.6 L115.3,182.6 L118.2,182.6 L121.0,182.6 L123.8,182.6 L126.7,182.6 L129.6,182.6 L132.4,182.6 L135.2,182.6 L138.1,182.6 L140.9,182.6 L143.8,182.6 L146.7,182.6 L149.5,182.6 L152.4,182.6 L155.2,182.6 L158.1,182.6 L160.9,182.6 L163.8,182.6 L166.6,182.6 L169.4,182.6 L172.3,182.6 L175.2,182.6 L178.0,182.6 L180.9,182.6 L183.7,182.6 L186.6,182.6 L189.4,182.6 L192.2,182.6 L195.1,182.6 L198.0,182.6 L200.8,182.6 L203.7,182.6 L206.5,182.6 L209.4,182.6 L212.2,182.6 L215.1,182.6 L217.9,182.6 L220.8,182.6 L223.6,182.6 L226.5,182.6 L229.3,182.6 L232.2,182.6 L235.0,182.6 L237.9,182.6 L240.7,182.6 L243.6,182.6 L246.4,182.6 L249.2,182.6 L252.1,182.6 L255.0,182.6 L257.8,182.6 L260.6,182.6 L263.5,182.6 L266.4,182.6 L269.2,182.6 L272.1,182.6 L274.9,182.6 L277.8,182.6 L280.6,182.6 L283.5,182.6 L286.3,182.6 L289.1,182.6 L292.0,182.6 L294.9,182.6 L297.7,182.6 L300.6,182.6 L303.4,182.6 L306.2,182.6 L309.1,182.6 L312.0,182.6 L314.8,182.6 L317.6,182.6 L320.5,182.6 L323.3,182.6 L326.2,182.6 L329.1,182.6 L331.9,182.6 L334.8,182.6 L337.6,182.6 L340.5,182.6 L343.3,182.6 L346.2,182.6 L349.0,182.6 L351.9,182.6 L354.7,182.6 L357.6,182.6 L360.4,182.6 L363.2,182.6 L366.1,182.6 L369.0,182.6 L371.8,182.6 L374.7,182.6 L377.5,182.6 L380.4,182.6 L383.2,182.6 L386.1,182.6 L388.9,182.6 L391.8,182.6 L394.6,182.6 L397.5,182.6 L400.3,182.6 L403.2,182.6 L406.0,182.6 L408.9,182.6 L411.7,182.6 L414.6,182.6 L417.4,182.6 L420.2,182.6 L423.1,182.6 L426.0,182.6 L428.8,182.6 L431.7,182.6 L434.5,182.6 L437.4,182.6 L440.2,182.6 L443.1,182.6 L445.9,182.6 L448.8,182.6 L451.6,182.6 L454.5,182.6 L457.3,182.6 L460.2,182.6 L463.0,182.6 L465.9,182.6 L468.7,182.6 L471.6,182.6 L474.4,182.6 L477.2,182.6 L480.1,182.6 L483.0,182.6 L485.8,182.6 L488.7,182.6 L491.5,182.6 L494.4,182.6 L497.2,182.6 L500.1,182.6 L502.9,182.6 L505.8,182.6 L508.6,182.6 L511.5,182.6" fill="none" stroke="#1baf7a" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><path d="M64.0,182.6 L66.8,182.6 L69.7,182.6 L72.5,182.6 L75.4,182.6 L78.2,182.6 L81.1,182.6 L84.0,182.6 L86.8,182.6 L89.7,182.6 L92.5,182.6 L95.3,182.6 L98.2,182.6 L101.1,182.6 L103.9,182.6 L106.8,182.6 L109.6,182.6 L112.5,182.6 L115.3,182.6 L118.2,182.6 L121.0,182.6 L123.8,182.6 L126.7,182.6 L129.6,182.6 L132.4,182.6 L135.2,182.6 L138.1,182.6 L140.9,182.6 L143.8,182.6 L146.7,182.6 L149.5,182.6 L152.4,182.6 L155.2,182.6 L158.1,182.6 L160.9,182.6 L163.8,182.6 L166.6,182.6 L169.4,182.6 L172.3,182.6 L175.2,182.6 L178.0,182.6 L180.9,182.6 L183.7,182.6 L186.6,182.6 L189.4,182.6 L192.2,182.6 L195.1,182.6 L198.0,182.6 L200.8,182.6 L203.7,182.6 L206.5,182.6 L209.4,182.6 L212.2,182.6 L215.1,182.6 L217.9,182.6 L220.8,182.6 L223.6,182.6 L226.5,182.6 L229.3,182.6 L232.2,182.6 L235.0,182.6 L237.9,182.6 L240.7,182.6 L243.6,182.6 L246.4,182.6 L249.2,182.6 L252.1,182.6 L255.0,182.6 L257.8,182.6 L260.6,182.6 L263.5,182.6 L266.4,182.6 L269.2,182.6 L272.1,182.6 L274.9,182.6 L277.8,182.6 L280.6,182.6 L283.5,182.6 L286.3,182.6 L289.1,182.6 L292.0,182.6 L294.9,182.6 L297.7,182.6 L300.6,182.6 L303.4,182.6 L306.2,182.6 L309.1,182.6 L312.0,182.6 L314.8,182.6 L317.6,182.6 L320.5,182.6 L323.3,182.6 L326.2,182.6 L329.1,182.6 L331.9,182.6 L334.8,182.6 L337.6,182.6 L340.5,182.6 L343.3,182.6 L346.2,182.6 L349.0,182.6 L351.9,182.6 L354.7,182.6 L357.6,182.6 L360.4,182.6 L363.2,182.6 L366.1,182.6 L369.0,182.6 L371.8,182.6 L374.7,182.6 L377.5,182.6 L380.4,182.6 L383.2,182.6 L386.1,182.6 L388.9,182.6 L391.8,182.6 L394.6,182.6 L397.5,182.6 L400.3,182.6 L403.2,182.6 L406.0,182.6 L408.9,182.6 L411.7,182.6 L414.6,182.6 L417.4,182.6 L420.2,182.6 L423.1,182.6 L426.0,182.6 L428.8,182.6 L431.7,182.6 L434.5,182.6 L437.4,182.6 L440.2,182.6 L443.1,182.6 L445.9,182.6 L448.8,182.6 L451.6,182.6 L454.5,182.6 L457.3,182.6 L460.2,182.6 L463.0,182.6 L465.9,182.6 L468.7,182.6 L471.6,182.6 L474.4,182.6 L477.2,182.6 L480.1,182.6 L483.0,182.6 L485.8,182.6 L488.7,182.6 L491.5,182.6 L494.4,182.6 L497.2,182.6 L500.1,182.6 L502.9,182.6 L505.8,182.6 L508.6,182.6 L511.5,182.6" fill="none" stroke="#eda100" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><text x="530" y="46.4" class="viz-label">Forward Euler</text><text x="530" y="326.9" class="viz-label">Backward Euler</text><text x="530" y="176.6" class="viz-label">Trapezoid</text><text x="530" y="196.6" class="viz-label">Leapfrog</text><line x1="0" x2="0" y1="30" y2="340" class="viz-crosshair" visibility="hidden"/></svg>
<div class="viz-tooltip" hidden></div>
</div>
<figcaption>The oscillator x′ = y, y′ = −x with x(0) = 1, y(0) = 0 and Δt = 0.1, for ten turns, on a log scale. The trapezoid and leapfrog lines lie on top of each other at 1. Leapfrog takes its second value from the exact solution. The table shows every second turn. Hover over the chart, or focus it and use the arrow keys, to read the values.</figcaption>
<details class="viz-table"><summary>Show the numbers</summary><table><thead><tr><th>t</th><th>Forward Euler</th><th>Backward Euler</th><th>Trapezoid</th><th>Leapfrog</th></tr></thead><tbody><tr><td>0.0</td><td>1.0000</td><td>1.0000</td><td>1.0000</td><td>1.0000</td></tr><tr><td>12.4</td><td>1.8532</td><td>0.5396</td><td>1.0000</td><td>1.0000</td></tr><tr><td>25.2</td><td>3.5034</td><td>0.2854</td><td>1.0000</td><td>1.0000</td></tr><tr><td>37.6</td><td>6.4926</td><td>0.1540</td><td>1.0000</td><td>1.0000</td></tr><tr><td>50.4</td><td>12.27</td><td>0.0815</td><td>1.0000</td><td>1.0000</td></tr><tr><td>62.8</td><td>22.75</td><td>0.0440</td><td>1.0000</td><td>1.0000</td></tr></tbody></table></details>
<script>(() => {const fig = document.currentScript.closest('figure');const svg = fig.querySelector('svg');const line = svg.querySelector('.viz-crosshair');const tip = fig.querySelector('.viz-tooltip');const names = ["Forward Euler", "Backward Euler", "Trapezoid", "Leapfrog"]; const colors = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"];const rows = [["1.0000", "1.0000", "1.0000", "1.0000"], ["1.0201", "0.9803", "1.0000", "1.0001"], ["1.0406", "0.9610", "1.0000", "1.0001"], ["1.0615", "0.9420", "1.0000", "1.0001"], ["1.0829", "0.9235", "1.0000", "1.0000"], ["1.1046", "0.9053", "1.0000", "1.0000"], ["1.1268", "0.8874", "1.0000", "0.9999"], ["1.1495", "0.8700", "1.0000", "0.9999"], ["1.1726", "0.8528", "1.0000", "1.0000"], ["1.1961", "0.8360", "1.0000", "1.0001"], ["1.2202", "0.8195", "1.0000", "1.0001"], ["1.2447", "0.8034", "1.0000", "1.0001"], ["1.2697", "0.7876", "1.0000", "1.0000"], ["1.2953", "0.7720", "1.0000", "0.9999"], ["1.3213", "0.7568", "1.0000", "0.9999"], ["1.3478", "0.7419", "1.0000", "1.0000"], ["1.3749", "0.7273", "1.0000", "1.0000"], ["1.4026", "0.7130", "1.0000", "1.0001"], ["1.4308", "0.6989", "1.0000", "1.0001"], ["1.4595", "0.6852", "1.0000", "1.0001"], ["1.4889", "0.6717", "1.0000", "1.0000"], ["1.5188", "0.6584", "1.0000", "0.9999"], ["1.5493", "0.6454", "1.0000", "0.9999"], ["1.5805", "0.6327", "1.0000", "1.0000"], ["1.6122", "0.6203", "1.0000", "1.0000"], ["1.6446", "0.6080", "1.0000", "1.0001"], ["1.6777", "0.5961", "1.0000", "1.0001"], ["1.7114", "0.5843", "1.0000", "1.0000"], ["1.7458", "0.5728", "1.0000", "1.0000"], ["1.7809", "0.5615", "1.0000", "0.9999"], ["1.8167", "0.5504", "1.0000", "0.9999"], ["1.8532", "0.5396", "1.0000", "1.0000"], ["1.8905", "0.5290", "1.0000", "1.0000"], ["1.9285", "0.5185", "1.0000", "1.0001"], ["1.9672", "0.5083", "1.0000", "1.0001"], ["2.0068", "0.4983", "1.0000", "1.0000"], ["2.0471", "0.4885", "1.0000", "1.0000"], ["2.0882", "0.4789", "1.0000", "0.9999"], ["2.1302", "0.4694", "1.0000", "0.9999"], ["2.1730", "0.4602", "1.0000", "1.0000"], ["2.2167", "0.4511", "1.0000", "1.0001"], ["2.2613", "0.4422", "1.0000", "1.0001"], ["2.3067", "0.4335", "1.0000", "1.0001"], ["2.3531", "0.4250", "1.0000", "1.0000"], ["2.4004", "0.4166", "1.0000", "1.0000"], ["2.4486", "0.4084", "1.0000", "0.9999"], ["2.4979", "0.4003", "1.0000", "0.9999"], ["2.5481", "0.3925", "1.0000", "1.0000"], ["2.5993", "0.3847", "1.0000", "1.0001"], ["2.6515", "0.3771", "1.0000", "1.0001"], ["2.7048", "0.3697", "1.0000", "1.0001"], ["2.7592", "0.3624", "1.0000", "1.0000"], ["2.8146", "0.3553", "1.0000", "1.0000"], ["2.8712", "0.3483", "1.0000", "0.9999"], ["2.9289", "0.3414", "1.0000", "0.9999"], ["2.9878", "0.3347", "1.0000", "1.0000"], ["3.0479", "0.3281", "1.0000", "1.0001"], ["3.1091", "0.3216", "1.0000", "1.0001"], ["3.1716", "0.3153", "1.0000", "1.0001"], ["3.2354", "0.3091", "1.0000", "1.0000"], ["3.3004", "0.3030", "1.0000", "0.9999"], ["3.3667", "0.2970", "1.0000", "0.9999"], ["3.4344", "0.2912", "1.0000", "1.0000"], ["3.5034", "0.2854", "1.0000", "1.0000"], ["3.5738", "0.2798", "1.0000", "1.0001"], ["3.6457", "0.2743", "1.0000", "1.0001"], ["3.7190", "0.2689", "1.0000", "1.0001"], ["3.7937", "0.2636", "1.0000", "1.0000"], ["3.8700", "0.2584", "1.0000", "0.9999"], ["3.9477", "0.2533", "1.0000", "0.9999"], ["4.0271", "0.2483", "1.0000", "1.0000"], ["4.1080", "0.2434", "1.0000", "1.0000"], ["4.1906", "0.2386", "1.0000", "1.0001"], ["4.2748", "0.2339", "1.0000", "1.0001"], ["4.3608", "0.2293", "1.0000", "1.0000"], ["4.4484", "0.2248", "1.0000", "1.0000"], ["4.5378", "0.2204", "1.0000", "0.9999"], ["4.6290", "0.2160", "1.0000", "0.9999"], ["4.7221", "0.2118", "1.0000", "1.0000"], ["4.8170", "0.2076", "1.0000", "1.0000"], ["4.9138", "0.2035", "1.0000", "1.0001"], ["5.0126", "0.1995", "1.0000", "1.0001"], ["5.1133", "0.1956", "1.0000", "1.0000"], ["5.2161", "0.1917", "1.0000", "1.0000"], ["5.3210", "0.1879", "1.0000", "0.9999"], ["5.4279", "0.1842", "1.0000", "0.9999"], ["5.5370", "0.1806", "1.0000", "1.0000"], ["5.6483", "0.1770", "1.0000", "1.0000"], ["5.7618", "0.1736", "1.0000", "1.0001"], ["5.8777", "0.1701", "1.0000", "1.0001"], ["5.9958", "0.1668", "1.0000", "1.0000"], ["6.1163", "0.1635", "1.0000", "1.0000"], ["6.2393", "0.1603", "1.0000", "0.9999"], ["6.3647", "0.1571", "1.0000", "0.9999"], ["6.4926", "0.1540", "1.0000", "1.0000"], ["6.6231", "0.1510", "1.0000", "1.0001"], ["6.7562", "0.1480", "1.0000", "1.0001"], ["6.8920", "0.1451", "1.0000", "1.0001"], ["7.0305", "0.1422", "1.0000", "1.0000"], ["7.1719", "0.1394", "1.0000", "1.0000"], ["7.3160", "0.1367", "1.0000", "0.9999"], ["7.4631", "0.1340", "1.0000", "0.9999"], ["7.6131", "0.1314", "1.0000", "1.0000"], ["7.7661", "0.1288", "1.0000", "1.0001"], ["7.9222", "0.1262", "1.0000", "1.0001"], ["8.0814", "0.1237", "1.0000", "1.0001"], ["8.2439", "0.1213", "1.0000", "1.0000"], ["8.4096", "0.1189", "1.0000", "0.9999"], ["8.5786", "0.1166", "1.0000", "0.9999"], ["8.7510", "0.1143", "1.0000", "1.0000"], ["8.9269", "0.1120", "1.0000", "1.0000"], ["9.1064", "0.1098", "1.0000", "1.0001"], ["9.2894", "0.1076", "1.0000", "1.0001"], ["9.4761", "0.1055", "1.0000", "1.0001"], ["9.6666", "0.1034", "1.0000", "1.0000"], ["9.8609", "0.1014", "1.0000", "0.9999"], ["10.06", "0.0994", "1.0000", "0.9999"], ["10.26", "0.0975", "1.0000", "1.0000"], ["10.47", "0.0955", "1.0000", "1.0000"], ["10.68", "0.0937", "1.0000", "1.0001"], ["10.89", "0.0918", "1.0000", "1.0001"], ["11.11", "0.0900", "1.0000", "1.0001"], ["11.33", "0.0882", "1.0000", "1.0000"], ["11.56", "0.0865", "1.0000", "0.9999"], ["11.80", "0.0848", "1.0000", "0.9999"], ["12.03", "0.0831", "1.0000", "1.0000"], ["12.27", "0.0815", "1.0000", "1.0000"], ["12.52", "0.0799", "1.0000", "1.0001"], ["12.77", "0.0783", "1.0000", "1.0001"], ["13.03", "0.0768", "1.0000", "1.0000"], ["13.29", "0.0752", "1.0000", "1.0000"], ["13.56", "0.0738", "1.0000", "0.9999"], ["13.83", "0.0723", "1.0000", "0.9999"], ["14.11", "0.0709", "1.0000", "1.0000"], ["14.39", "0.0695", "1.0000", "1.0000"], ["14.68", "0.0681", "1.0000", "1.0001"], ["14.98", "0.0668", "1.0000", "1.0001"], ["15.28", "0.0655", "1.0000", "1.0000"], ["15.58", "0.0642", "1.0000", "1.0000"], ["15.90", "0.0629", "1.0000", "0.9999"], ["16.22", "0.0617", "1.0000", "0.9999"], ["16.54", "0.0604", "1.0000", "1.0000"], ["16.88", "0.0593", "1.0000", "1.0001"], ["17.22", "0.0581", "1.0000", "1.0001"], ["17.56", "0.0569", "1.0000", "1.0001"], ["17.91", "0.0558", "1.0000", "1.0000"], ["18.27", "0.0547", "1.0000", "1.0000"], ["18.64", "0.0536", "1.0000", "0.9999"], ["19.02", "0.0526", "1.0000", "0.9999"], ["19.40", "0.0516", "1.0000", "1.0000"], ["19.79", "0.0505", "1.0000", "1.0001"], ["20.19", "0.0495", "1.0000", "1.0001"], ["20.59", "0.0486", "1.0000", "1.0001"], ["21.01", "0.0476", "1.0000", "1.0000"], ["21.43", "0.0467", "1.0000", "0.9999"], ["21.86", "0.0457", "1.0000", "0.9999"], ["22.30", "0.0448", "1.0000", "1.0000"], ["22.75", "0.0440", "1.0000", "1.0000"]]; const heads = ["t = 0.0", "t = 0.4", "t = 0.8", "t = 1.2", "t = 1.6", "t = 2.0", "t = 2.4", "t = 2.8", "t = 3.2", "t = 3.6", "t = 4.0", "t = 4.4", "t = 4.8", "t = 5.2", "t = 5.6", "t = 6.0", "t = 6.4", "t = 6.8", "t = 7.2", "t = 7.6", "t = 8.0", "t = 8.4", "t = 8.8", "t = 9.2", "t = 9.6", "t = 10.0", "t = 10.4", "t = 10.8", "t = 11.2", "t = 11.6", "t = 12.0", "t = 12.4", "t = 12.8", "t = 13.2", "t = 13.6", "t = 14.0", "t = 14.4", "t = 14.8", "t = 15.2", "t = 15.6", "t = 16.0", "t = 16.4", "t = 16.8", "t = 17.2", "t = 17.6", "t = 18.0", "t = 18.4", "t = 18.8", "t = 19.2", "t = 19.6", "t = 20.0", "t = 20.4", "t = 20.8", "t = 21.2", "t = 21.6", "t = 22.0", "t = 22.4", "t = 22.8", "t = 23.2", "t = 23.6", "t = 24.0", "t = 24.4", "t = 24.8", "t = 25.2", "t = 25.6", "t = 26.0", "t = 26.4", "t = 26.8", "t = 27.2", "t = 27.6", "t = 28.0", "t = 28.4", "t = 28.8", "t = 29.2", "t = 29.6", "t = 30.0", "t = 30.4", "t = 30.8", "t = 31.2", "t = 31.6", "t = 32.0", "t = 32.4", "t = 32.8", "t = 33.2", "t = 33.6", "t = 34.0", "t = 34.4", "t = 34.8", "t = 35.2", "t = 35.6", "t = 36.0", "t = 36.4", "t = 36.8", "t = 37.2", "t = 37.6", "t = 38.0", "t = 38.4", "t = 38.8", "t = 39.2", "t = 39.6", "t = 40.0", "t = 40.4", "t = 40.8", "t = 41.2", "t = 41.6", "t = 42.0", "t = 42.4", "t = 42.8", "t = 43.2", "t = 43.6", "t = 44.0", "t = 44.4", "t = 44.8", "t = 45.2", "t = 45.6", "t = 46.0", "t = 46.4", "t = 46.8", "t = 47.2", "t = 47.6", "t = 48.0", "t = 48.4", "t = 48.8", "t = 49.2", "t = 49.6", "t = 50.0", "t = 50.4", "t = 50.8", "t = 51.2", "t = 51.6", "t = 52.0", "t = 52.4", "t = 52.8", "t = 53.2", "t = 53.6", "t = 54.0", "t = 54.4", "t = 54.8", "t = 55.2", "t = 55.6", "t = 56.0", "t = 56.4", "t = 56.8", "t = 57.2", "t = 57.6", "t = 58.0", "t = 58.4", "t = 58.8", "t = 59.2", "t = 59.6", "t = 60.0", "t = 60.4", "t = 60.8", "t = 61.2", "t = 61.6", "t = 62.0", "t = 62.4", "t = 62.8"];const xs = [64.0, 66.8, 69.7, 72.5, 75.4, 78.2, 81.1, 84.0, 86.8, 89.7, 92.5, 95.3, 98.2, 101.1, 103.9, 106.8, 109.6, 112.5, 115.3, 118.2, 121.0, 123.8, 126.7, 129.6, 132.4, 135.2, 138.1, 140.9, 143.8, 146.7, 149.5, 152.4, 155.2, 158.1, 160.9, 163.8, 166.6, 169.4, 172.3, 175.2, 178.0, 180.9, 183.7, 186.6, 189.4, 192.2, 195.1, 198.0, 200.8, 203.7, 206.5, 209.4, 212.2, 215.1, 217.9, 220.8, 223.6, 226.5, 229.3, 232.2, 235.0, 237.9, 240.7, 243.6, 246.4, 249.2, 252.1, 255.0, 257.8, 260.6, 263.5, 266.4, 269.2, 272.1, 274.9, 277.8, 280.6, 283.5, 286.3, 289.1, 292.0, 294.9, 297.7, 300.6, 303.4, 306.2, 309.1, 312.0, 314.8, 317.6, 320.5, 323.3, 326.2, 329.1, 331.9, 334.8, 337.6, 340.5, 343.3, 346.2, 349.0, 351.9, 354.7, 357.6, 360.4, 363.2, 366.1, 369.0, 371.8, 374.7, 377.5, 380.4, 383.2, 386.1, 388.9, 391.8, 394.6, 397.5, 400.3, 403.2, 406.0, 408.9, 411.7, 414.6, 417.4, 420.2, 423.1, 426.0, 428.8, 431.7, 434.5, 437.4, 440.2, 443.1, 445.9, 448.8, 451.6, 454.5, 457.3, 460.2, 463.0, 465.9, 468.7, 471.6, 474.4, 477.2, 480.1, 483.0, 485.8, 488.7, 491.5, 494.4, 497.2, 500.1, 502.9, 505.8, 508.6, 511.5]; const steps = 157, width = 680;let current = -1;const show = (n) => { n = Math.max(0, Math.min(steps, n)); current = n; const x = xs[n]; line.setAttribute('x1', x); line.setAttribute('x2', x); line.setAttribute('visibility', 'visible'); tip.replaceChildren(); const head = document.createElement('div'); head.className = 'viz-tooltip-head'; head.textContent = heads[n]; tip.append(head); names.forEach((name, i) => {  const row = document.createElement('div'); row.className = 'viz-tooltip-row';  const key = document.createElement('span'); key.className = 'viz-swatch'; key.style.background = colors[i];  const value = document.createElement('strong'); value.textContent = rows[n][i];  const label = document.createElement('span'); label.textContent = name;  row.append(key, value, label); tip.append(row); }); const box = svg.getBoundingClientRect(), frame = fig.getBoundingClientRect(); const xp = box.left - frame.left + x * box.width / width; tip.hidden = false; const room = frame.width - tip.offsetWidth - 8; tip.style.left = Math.max(0, Math.min(room, xp + 12)) + 'px'; tip.style.top = (box.top - frame.top + 8) + 'px';};const hide = () => { line.setAttribute('visibility', 'hidden'); tip.hidden = true; current = -1; };svg.addEventListener('pointermove', (event) => { const box = svg.getBoundingClientRect(); const x = (event.clientX - box.left) * width / box.width; let best = 0; xs.forEach((v, i) => { if (Math.abs(v - x) < Math.abs(xs[best] - x)) best = i; }); show(best);});svg.addEventListener('pointerleave', hide);svg.addEventListener('focus', () => show(current < 0 ? 0 : current));svg.addEventListener('blur', hide);svg.addEventListener('keydown', (event) => { if (event.key === 'ArrowRight') { show(current + 1); event.preventDefault(); } if (event.key === 'ArrowLeft') { show(current - 1); event.preventDefault(); } if (event.key === 'Escape') hide();});})();</script>
</figure>

| Method | size at $t = 62.8$ | theory $\lvert P(0.1i)\rvert^{628}$ |
|---|---|---|
| Forward Euler | $22.746$ | $(1.01)^{314} = 22.746$ |
| Backward Euler | $0.0440$ | $(1.01)^{-314} = 0.0440$ |
| Trapezoid | $1.000000000000032$ | $1$ |
| Leapfrog | between $0.99992$ and $1.00009$ the whole time | roots on the circle |

- The imaginary axis is **outside** the Forward Euler disk. Each step multiplies the size by $|1 + 0.1i| = \sqrt{1.01} \approx 1.005$. That looks harmless, but after 628 steps the circle has become $23$ times bigger.
- The imaginary axis is **inside** the Backward Euler region, but not on its edge, so the method damps: the circle shrinks to $4\%$ of its size.
- For the trapezoid, $|P(iy)| = 1$ exactly, because the imaginary axis is the boundary of its region. The size stays $1$ up to rounding.
- Leapfrog also keeps the size, because $0.1i$ is inside its segment $(-i, i)$.
- For the three one-step methods, the measured size agrees with $|P(z)|^n$ to all printed digits.

**Part 2: the root condition for leapfrog.** Now take the complex test equation $u' = iu$ with $U^0 = 1$ and $U^1 = e^{i\Delta t}$, so $z = i\Delta t$.

| $\Delta t$ | roots of $\eta^2 - 2z\eta - 1$ | $\lvert\eta\rvert$ | $\lvert U^{10}\rvert$ | $\lvert U^{100}\rvert$ | $\lvert U^{1000}\rvert$ |
|---|---|---|---|---|---|
| 0.5 | $\pm 0.866 + 0.5i$ | $1,\ 1$ | $1.00$ | $1.00$ | $1.00$ |
| 1.0 | $i$ (double) | $1,\ 1$ | $5.43$ | $56.0$ | $563$ |
| 1.1 | $1.558i,\ 0.642i$ | $1.558,\ 0.642$ | $47.7$ | $1.04 \times 10^{19}$ | $2.46 \times 10^{192}$ |

- Two simple roots on the circle: the solution stays bounded (the largest size seen is $1.02$).
- A double root on the circle: the size grows **linearly**, about $0.56\,n$. This is the $n\,\eta^n$ solution, and it is why the root condition asks for $|\eta| < 1$ at a repeated root.
- One root outside the circle: exponential growth, about $1.558^n$.

**Part 3: leapfrog on a decaying problem.** Take $u' = -u$, $u(0) = 1$, so $z = -\Delta t$ is real and negative. Start with the exact $U^1$.

First, a fixed final time $T = 2$. Leapfrog is zero-stable (at $z = 0$ its roots are $\pm 1$, both simple), so it must converge:

| $\Delta t$ | 1/10 | 1/20 | 1/40 | 1/80 | 1/160 | 1/320 |
|---|---|---|---|---|---|---|
| error at $t = 2$ | $9.90 \times 10^{-4}$ | $1.84 \times 10^{-4}$ | $3.74 \times 10^{-5}$ | $8.21 \times 10^{-6}$ | $1.91 \times 10^{-6}$ | $4.59 \times 10^{-7}$ |
| order | | 2.42 | 2.30 | 2.19 | 2.11 | 2.06 |

It converges with order $2$ (the order approaches $2$ from above). Now run longer:

| $\Delta t$ | error at $t = 2$ | $t = 10$ | $t = 20$ | $t = 30$ |
|---|---|---|---|---|
| 0.1 | $9.9 \times 10^{-4}$ | $1.62$ | $3.5 \times 10^{4}$ | $7.6 \times 10^{8}$ |
| 0.01 | $5.1 \times 10^{-6}$ | $1.8 \times 10^{-3}$ | $40$ | $8.8 \times 10^{5}$ |
| exact $u(t)$ | $0.135$ | $4.5 \times 10^{-5}$ | $2.1 \times 10^{-9}$ | $9.4 \times 10^{-14}$ |

The exact solution decays, but the error grows by about $e^{10} \approx 2.2 \times 10^4$ every $10$ time units. The second root is $\eta_2 = -\Delta t - \sqrt{1 + \Delta t^2}$ ($-1.105$ for $\Delta t = 0.1$, $-1.010$ for $\Delta t = 0.01$). Over a time $t$ it grows like $|\eta_2|^{t/\Delta t} \approx e^{t}$, **whatever** $\Delta t$ is. A smaller step only makes the starting seed smaller, so the blow-up comes later. It never goes away. This is "no real $z$ except $0$ is stable" in numbers.

**Part 4: the boundary locus.** I computed $z(\theta) = \rho(e^{i\theta})/\sigma(e^{i\theta})$ at $12$ angles for the four methods. (For Backward Euler, $\rho(\eta) = \eta - 1$ and $\sigma(\eta) = \eta$. For the trapezoid, $\sigma(\eta) = (\eta + 1)/2$. For leapfrog, $\rho(\eta) = \eta^2 - 1$ and $\sigma(\eta) = 2\eta$.)

- Forward Euler: $z(\theta) = e^{i\theta} - 1$, points such as $-1 + i$ and $-2$. The circle around $-1$.
- Backward Euler: $z(\theta) = 1 - e^{-i\theta}$, points such as $1 + i$ and $2$. The circle around $+1$.
- Trapezoid: $z(\theta) = 2i\tan(\theta/2)$, all on the imaginary axis ($0.536i$, $1.155i$, $2i$, ...). At $\theta = \pi$, $\sigma = 0$ and the point goes to infinity. The boundary is the whole imaginary axis.
- Leapfrog: $z(\theta) = i\sin\theta$, which goes back and forth on the segment from $-i$ to $i$. The "curve" has no inside at all, which matches section 6.

Then the test step: the largest root size at a few points.

| test point $z$ | Forward Euler | Backward Euler | Trapezoid | Leapfrog |
|---|---|---|---|---|
| $-1$ | $0$ | $0.5$ | $0.333$ | $2.414$ |
| $-0.1$ | $0.9$ | $0.909$ | $0.905$ | $1.105$ |
| $0.5i$ | $1.118$ | $0.894$ | $1.000$ | $1.000$ |
| $+3$ | $4$ | $0.5$ | $5$ | $6.162$ |

$z = -1$ is inside the Forward Euler circle and $z = +3$ is outside the Backward Euler circle, so each region is on the side where the test point gives a size of at most $1$. The point $+3$ also shows the caution of section 3: Backward Euler calls it stable.

**Part 5: Backward Euler hides growth.** For $u' = u$ with $\Delta t = 3$, $z = 3$ and $P(3) = 1/(1 - 3) = -0.5$. The exact solution grows by $e^3 \approx 20.1$ per step. After four steps the exact value is $1.63 \times 10^5$, and Backward Euler gives $0.0625$, with a sign change at every step. Stable, but completely wrong. A very large step on a growing problem is a bad idea, whatever the region says.

**Part 6: the heat equation.** I use $u_t = u_{xx}$ on $[0, 1]$ with $u = 0$ at both ends and the exact solution $u = e^{-\pi^2 t}\sin(\pi x)$. First, the truncation error: I put the exact solution into each scheme at $t = 0.05$ and take the largest value over the grid points.

| $h$ | FTCS, $k = h$ | ratio to formula | CN, $k = h$ | FTCS, $k = h^2/4$ | FTCS, $k = h^2/6$ |
|---|---|---|---|---|---|
| 1/10 | $2.15$ | 0.73 | $0.272$ | $2.43 \times 10^{-2}$ | $1.08 \times 10^{-4}$ |
| 1/20 | $1.26$ | 0.85 | $8.62 \times 10^{-2}$ | $6.17 \times 10^{-3}$ | $6.78 \times 10^{-6}$ |
| 1/40 | $0.683$ | 0.92 | $2.43 \times 10^{-2}$ | $1.55 \times 10^{-3}$ | $4.24 \times 10^{-7}$ |
| 1/80 | $0.356$ | 0.96 | $6.46 \times 10^{-3}$ | $3.87 \times 10^{-4}$ | $2.65 \times 10^{-8}$ |
| 1/160 | $0.182$ | 0.98 | $1.66 \times 10^{-3}$ | $9.68 \times 10^{-5}$ | $1.67 \times 10^{-9}$ |
| last order | 0.97 | | 1.96 | 2.00 | 3.99 |

- With $k = h$, FTCS is first order: the $\tfrac12 k$ term wins. The ratio to the formula $(\tfrac12 k - \tfrac1{12}h^2)\,u_{xxxx}$ goes to $1$ as $k$ gets small. Crank–Nicolson is second order with the same $k$.
- With $k = h^2/4$, both terms are $O(h^2)$. The measured value agrees with the formula to four digits at $h = 1/160$ ($9.678 \times 10^{-5}$ against $9.679 \times 10^{-5}$).
- With $k = h^2/6$, the leading term cancels and the error drops by $16$ at each halving: **fourth order**, as the formula predicts.

Second, the real error of Crank–Nicolson at $t = 0.1$ with $k = h$. Each step solves the tridiagonal system with the Thomas algorithm.

| $h = k$ | 1/10 | 1/20 | 1/40 | 1/80 | 1/160 | 1/320 |
|---|---|---|---|---|---|---|
| $r = k/(2h^2)$ | 5 | 10 | 20 | 40 | 80 | 160 |
| max error | $2.99 \times 10^{-2}$ | $6.88 \times 10^{-3}$ | $1.69 \times 10^{-3}$ | $4.20 \times 10^{-4}$ | $1.05 \times 10^{-4}$ | $2.62 \times 10^{-5}$ |
| order | | 2.12 | 2.03 | 2.01 | 2.00 | 2.00 |

The global error is second order, like the truncation error. Note that $r$ is as large as $160$ and nothing goes wrong. Why Crank–Nicolson allows such large steps, and what happens to FTCS, is the topic of the next lecture.

The full script and the chart code are in [this site's repository](https://github.com/doubleyands/doubleyands.github.io/tree/main/scripts/posts/numerical-pde/lecture-9).

## Summary

1. Absolute stability is about $z = \lambda\Delta t$. For a system, $\lambda$ runs over the **eigenvalues** (of the Jacobian, for nonlinear problems), so $z$ can be **complex**.
2. A one-step method gives $U^{n+1} = P(z)U^n$, and its **stability region** is $|P(z)| \le 1$: inside a circle for Forward Euler, outside a circle for Backward Euler, the left half plane for the trapezoid.
3. $P(z)$ approximates $e^z$: a Taylor polynomial for explicit methods, a ratio of polynomials (Padé) for implicit ones.
4. A multistep method has a **stability polynomial** $\pi(\eta; z) = \rho(\eta) - z\sigma(\eta)$. The **root condition**: all roots have $|\eta| \le 1$, and repeated roots have $|\eta| < 1$.
5. Leapfrog is stable only on the segment $(-i, i)$. On a decaying problem it converges for fixed $T$, but grows like $e^t$ in the long run.
6. The **boundary locus** $z(\theta) = \rho(e^{i\theta})/\sigma(e^{i\theta})$ draws the boundary. A test point tells which side is stable.
7. The heat equation joins space and time. **FTCS** has truncation error $O(k + h^2)$, and **Crank–Nicolson** (the trapezoid in time) has $O(k^2 + h^2)$ and needs a tridiagonal solve each step.
