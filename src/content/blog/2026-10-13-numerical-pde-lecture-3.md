---
title: "Numerical PDEs, Lecture 3: Error Analysis, Neumann Conditions, and the 5-Point Stencil"
description: "My notes on the third lecture: how the local truncation error controls the global error, three ways to handle a Neumann condition, why pure Neumann problems are singular, and how the 1D scheme becomes a 2D one."
pubDate: 2026-10-13
project: numerical-pde
tags: [numerical PDEs, finite differences, boundary value problems, notes]
---

This is my summary of the third lecture in my numerical PDE course. I write the ideas in my own words and check them with a small experiment.

[Lecture 2](/blog/2026-10-08-numerical-pde-lecture-2/) built the three-point scheme for $u'' = f$ with Dirichlet conditions and found its local truncation error. This lecture finishes the one-dimensional story. First, we connect the local truncation error to the real error. Then we change the boundary condition to a Neumann condition, and see that one careless boundary row can ruin the whole solution. Last, we move to two dimensions.

## 1. Solving the system, and round-off

The scheme gives a linear system $AU = F$ with

$$
A = \frac{1}{h^2} \operatorname{tridiag}(1, -2, 1).
$$

There are many ways to solve it: the Thomas algorithm, general Gaussian elimination, or iterative methods like conjugate gradients. Any good solver is fine.

But a computer works with finite precision. In double precision the machine epsilon is $2^{-52} \approx 2.2 \times 10^{-16}$, and every operation can add an error of about this relative size. One error is tiny. But a large PDE system needs a huge number of operations, and the errors can pile up. They can be controlled, but only if we know how.

**Scaling** helps. If we multiply one row of a linear system by a nonzero number, the solution does not change. But the conditioning of the matrix can change. A good habit is to scale each row so that its coefficients are of size $O(1)$, like the $1, -2, 1$ of the interior rows. We will use this in section 3.

## 2. From local error to global error

We computed $U$. How close is it to the exact solution? We need three objects, and it helps to keep their names apart:

- $u$: the exact solution, a function;
- $U = (U_1, \ldots, U_{n-1})^{\mathsf{T}}$: the computed solution, a vector;
- $\hat{U} = \big(u(x_1), \ldots, u(x_{n-1})\big)^{\mathsf{T}}$: the exact solution sampled at the grid points, also a vector.

**Step 1. The global error** is $E = U - \hat{U}$, so $E_j = U_j - u(x_j)$. The sign does not matter, but we must use one sign consistently.

**Step 2. The local truncation error.** Why do we plug $u(x_j)$ into the scheme, and not $U_j$? The $U_j$ are just numbers at grid points. They have no derivatives. The function $u$ is smooth, so we can Taylor expand it. After we plug it in, the scheme is not satisfied exactly. The leftover is exactly the part of the Taylor series that we threw away, and that is why it is called the *truncation* error. In vector form,

$$
\tau = A\hat{U} - F, \qquad \tau_j = \frac{h^2}{12} u''''(x_j) + O(h^4).
$$

**Step 3. The link.** Apply $A$ to $E$:

$$
AE = AU - A\hat{U} = F - A\hat{U} = -\tau.
$$

So $AE = -\tau$. Row by row, this reads

$$
\frac{E_{j-1} - 2E_j + E_{j+1}}{h^2} = -\tau_j, \qquad E_0 = E_n = 0,
$$

because the boundary values are exact. Look at the left side: it is the same three-point difference, now applied to $E$. So the error solves a discrete version of its own boundary value problem:

$$
e''(x) = -\tau(x), \qquad e(0) = e(1) = 0.
$$

This is the key idea: **the error satisfies the same kind of problem as the solution**, with $-\tau$ as the source. With $\tau(x) \approx \frac{h^2}{12} u''''(x)$ we can integrate twice and use the two boundary conditions:

$$
e(x) \approx \frac{h^2}{12} \Big( -u''(x) + u''(0) + x \big(u''(1) - u''(0)\big) \Big).
$$

Everything in the bracket is bounded if $u$ is smooth, so the global error is $O(h^2)$. The experiment checks this formula, not only its order.

**A note on Lecture 2.** There, a rough solution $u = x^{5/2}$ had $\max|\tau|$ of order only $1/2$, but the global error was still second order. Now we can see why. From $AE = -\tau$ we get $E = -A^{-1}\tau$. The matrix $A$ is sparse, but $A^{-1}$ is dense, so each $\tau_k$ affects every $E_j$. For this matrix, the entries are known exactly:

$$
(A^{-1})_{jk} = h\, G(x_j, x_k), \qquad G(x, s) = \begin{cases} x(s - 1), & x \le s, \\ s(x - 1), & x > s. \end{cases}
$$

(I checked this formula numerically.) The weight $G(x, s)$ goes to zero when $s$ is near the boundary. The large $\tau_1$ sits at $s = h$, so its weight is at most $h \cdot h$, and its effect is only of size $h^2 \cdot h^{1/2}$. A large local error next to a Dirichlet boundary is harmless.

## 3. A Neumann condition at one end

Now we change the left end:

$$
u''(x) = f(x), \quad 0 < x < 1, \qquad u'(0) = \sigma, \quad u(1) = \beta.
$$

At the left we now know the slope, not the value. For heat, this fixes the heat flux through the left end. One end is Dirichlet and the other is Neumann, so strictly this is a *mixed* problem. But the new part is the Neumann end.

**What is unknown?** We know $u'(0)$, but not $u(0)$. So $U_0$ is now an unknown too. We have $n + 1$ unknowns $U_0, \ldots, U_n$. The last one is $U_n = \beta$, which can move to the right-hand side as before. The only real question is how to discretize the first equation, $u'(0) = \sigma$.

### First approach: a one-sided first-order difference

There is no grid point to the left of $x_0$, so we use the point on the right:

$$
\frac{U_1 - U_0}{h} = \sigma.
$$

If we keep the common factor $1/h^2$ in front of the matrix, this row has coefficients $-h, h$, while the other rows have $1, -2, 1$. That is a bad mix of sizes. We scale the first row by dividing it by $h$, as in section 1. Then its coefficients are $-1, 1$ (times $1/h$) and the right side is $\sigma$. The matrix is still tridiagonal, but it is no longer symmetric. That is fine: it is still invertible, and the Thomas algorithm works for nonsymmetric tridiagonal matrices too.

**The problem is accuracy.** Put the exact solution into the first row:

$$
\tau_0 = \frac{u(x_1) - u(x_0)}{h} - \sigma = u'(0) + \frac{h}{2} u''(0) + O(h^2) - \sigma = \frac{h}{2} u''(0) + O(h^2).
$$

This is only first order. All interior rows are second order, but one boundary row is first order. How much does one row matter? A lot. Repeat the argument of section 2: the error now satisfies $e'' \approx 0$ inside, $e(1) = 0$, and the boundary row gives $e'(0) \approx -\frac{h}{2} u''(0)$. So

$$
e(x) \approx \frac{h}{2} u''(0)\, (1 - x).
$$

The error is $O(h)$ **everywhere**, not only near $x = 0$. The low accuracy of one row spreads into the whole domain. This is called the **pollution effect**. Even where the exact solution is very smooth and $\tau$ is tiny, the error is first order. The scheme is *suboptimal*: it gives less than its interior rows could give.

This is not a contradiction with section 2. There, the large $\tau$ was next to a Dirichlet boundary, where the solution is pinned. Here, the error is in the *slope* at the boundary, and nothing pins the solution at $x = 0$. A wrong slope tilts the whole solution.

### Second approach: a ghost node

We want a second-order boundary row. The central difference is second order, but it needs a point to the left of $x_0$. So we invent one: a **ghost node** $x_{-1} = -h$, outside the domain, with an unknown value $U_{-1}$. Then we write two equations at $x_0$:

$$
\frac{U_{-1} - 2U_0 + U_1}{h^2} = f(x_0), \qquad \frac{U_1 - U_{-1}}{2h} = \sigma.
$$

The first one says that the differential equation also holds at the boundary point. The second one is the boundary condition with a central difference, which is second order.

$U_{-1}$ is not a real value, so we remove it. The first equation gives $U_{-1} = 2U_0 - U_1 + h^2 f(x_0)$. Put this into the second one:

$$
\frac{U_1 - U_0}{h} = \sigma + \frac{h}{2} f(x_0).
$$

The left side is the same as in the first approach. The matrix does not change at all. Only the first entry of the right-hand side changes, by $\frac{h}{2} f(x_0)$. This is exactly the term that cancels $\tau_0$, because $f(x_0) = u''(0)$. For almost no extra work, the boundary row becomes second order.

### Third approach: a one-sided second-order difference

We can also stay inside the domain and use one more point. In Lecture 2 we built a three-point one-sided formula for the *right* end. At the left end we mirror it. Mirroring $h \to -h$ also flips the sign:

$$
\frac{-3U_0 + 4U_1 - U_2}{2h} = \sigma.
$$

The weights still add up to zero: $-3 + 4 - 1 = 0$. This row is second order, but it touches three unknowns, so the first row sticks out of the tridiagonal band. And if we want even higher order, a one-sided formula must reach further and further inside.

A small observation of my own: the interior equation at $j = 1$ gives $U_2 = 2U_1 - U_0 + h^2 f(x_1)$. Using it to remove $U_2$ turns the third approach into

$$
\frac{U_1 - U_0}{h} = \sigma + \frac{h}{2} f(x_1).
$$

So it is the ghost node row with $f(x_1)$ in place of $f(x_0)$. The experiment uses this form, so all three systems are tridiagonal.

| approach | order of the row | unknowns in the first row | matrix |
|---|---|---|---|
| one-sided, 1st order | 1 (pollutes the whole solution) | $U_0, U_1$ | tridiagonal |
| ghost node | 2 | $U_0, U_1$ | tridiagonal, only $F$ changes |
| one-sided, 2nd order | 2 | $U_0, U_1, U_2$ | one entry outside the band |

The lesson: if we allow ourselves to step outside the domain for a moment, we get a very clean result.

## 4. Pure Neumann: existence and uniqueness

Now put a Neumann condition at both ends:

$$
u''(x) = f(x), \quad 0 < x < 1, \qquad u'(0) = \sigma_0, \quad u'(1) = \sigma_1.
$$

"Pure" means that the value of $u$ is fixed nowhere. Only slopes are given. This looks harmless, but two things go wrong.

**Uniqueness fails.** If $u$ is a solution, then so is $u + c$ for any constant $c$. Adding a constant does not change $u''$, $u'(0)$, or $u'(1)$. So there are infinitely many solutions.

**Existence can fail.** Integrate the equation over $(0, 1)$:

$$
\int_0^1 f(x)\, dx = \int_0^1 u''(x)\, dx = u'(1) - u'(0) = \sigma_1 - \sigma_0.
$$

If the data do not satisfy this **compatibility condition**, there is no solution at all. For example, with $f(x) = 4e^{2x}$ the left side is $2e^2 - 2$, so we need $\sigma_1 - \sigma_0 = 2e^2 - 2$. With any other $\sigma_1$, nothing works.

Physically, the left side is the total heat made inside the rod per unit time, and the right side is the net heat that leaves through the two ends. In a steady state these must balance. Picture a rod with both ends perfectly insulated ($\sigma_0 = \sigma_1 = 0$) and a heater inside ($f \ne 0$). The heat has nowhere to go, so the temperature rises forever. There is no steady state. On the other hand, with $\sigma_0 = \sigma_1 = 0$ and $f = 0$, every constant is a solution: existence holds, but uniqueness does not.

**The discrete problem behaves the same way.** Use ghost nodes at both ends and scale the boundary rows. Then

$$
\text{row } 0: \frac{U_1 - U_0}{h} = \sigma_0 + \frac{h}{2} f(x_0), \qquad
\text{row } n: \frac{U_{n-1} - U_n}{h} = -\sigma_1 + \frac{h}{2} f(x_n),
$$

with the usual $1, -2, 1$ rows in between. Now use the zero-sum check from Lecture 2: every row of the matrix adds up to zero. That means

$$
A \begin{bmatrix} 1 \\ 1 \\ \vdots \\ 1 \end{bmatrix} = \mathbf{0}.
$$

The constant vector is in the null space, so **$A$ is singular**. With a Dirichlet condition, the known boundary value moved to the right-hand side and broke this zero sum. A pure Neumann problem has nothing to break it.

A singular system has either no solution or infinitely many. Which one depends on $F$. Here is a short calculation of my own for the discrete compatibility condition. Add row $0$, row $n$, and $h$ times each interior row. The left sides telescope to zero. So a solution can exist only if the right sides also add up to zero:

$$
h \Big( \tfrac{1}{2} f(x_0) + f(x_1) + \cdots + f(x_{n-1}) + \tfrac{1}{2} f(x_n) \Big) = \sigma_1 - \sigma_0.
$$

The left side is the **trapezoid rule** for $\int_0^1 f\, dx$. So the discrete condition is the continuous one with the integral replaced by the trapezoid rule. In general the trapezoid rule is not exact, so even if the continuous problem is compatible, the discrete one may not be. The experiment shows this.

So with pure Neumann conditions, we must always check two things: is the problem compatible, and how do we pick one solution out of the infinitely many?

## 5. The steady heat equation in two dimensions

In 2D the steady heat equation is

$$
\nabla \cdot (\kappa \nabla u) = (\kappa u_x)_x + (\kappa u_y)_y = f.
$$

If $\kappa$ is a constant, we divide by it and get the **Poisson equation** $\nabla^2 u = f$. If also $f = 0$ (no source), we get the **Laplace equation** $\nabla^2 u = 0$.

The boundary conditions look the same as before:

$$
u = g \ \text{ on } \partial\Omega \quad \text{(Dirichlet)}, \qquad \frac{\partial u}{\partial \mathbf{n}} = \mathbf{n} \cdot \nabla u = g_N \ \text{ on } \partial\Omega \quad \text{(Neumann)},
$$

where $\mathbf{n}$ is the unit outward normal, so $\mathbf{n} \cdot \nabla u$ is the flux out through the boundary. If $g = 0$ or $g_N = 0$, the condition is *homogeneous*.

We take the unit square $\Omega = (0, 1) \times (0, 1)$. A domain should be an open set, so we write open intervals. To keep things simple: Dirichlet conditions only and $\kappa = 1$. A Neumann condition can be handled with ghost nodes, just as in 1D.

### The 5-point stencil

Use grid points $x_i = i\,\Delta x$, $y_j = j\,\Delta y$ and write $u_{i,j}$ for the approximation of $u(x_i, y_j)$.

The idea is that **2D is just 1D twice**. Fix the row $y = y_j$. Along this row, $u_{xx}$ is a 1D second derivative, and we already know the three-point formula for it. Do the same for $u_{yy}$ along the column $x = x_i$. Then add:

$$
\frac{u_{i-1,j} - 2u_{i,j} + u_{i+1,j}}{(\Delta x)^2} + \frac{u_{i,j-1} - 2u_{i,j} + u_{i,j+1}}{(\Delta y)^2} = f_{i,j}.
$$

With equal spacing $\Delta x = \Delta y = h$:

$$
\frac{u_{i-1,j} + u_{i+1,j} + u_{i,j-1} + u_{i,j+1} - 4u_{i,j}}{h^2} = f_{i,j}.
$$

It uses the center point and its four neighbors (east, west, north, south), so it is called the **5-point stencil**. The center weight is $-4$ because the center appears in both the $x$ part ($-2$) and the $y$ part ($-2$). The weights still add up to zero: $1 + 1 + 1 + 1 - 4 = 0$.

**Accuracy.** Each 1D part has the same Taylor error as before, so

$$
\tau_{i,j} = \frac{h^2}{12} \big( u_{xxxx} + u_{yyyy} \big) + O(h^4).
$$

The error equation $AE = -\tau$ holds again, and the same argument as in 1D gives a global error of $O(h^2)$.

### Ordering the unknowns

In finite differences, the unknowns are always the values at the grid points. In 1D there was one natural order: left to right. In 2D there is no natural order, so we must choose one before we can write the matrix. A random order would be a nightmare.

The standard choice is **lexicographic order** (like words in a dictionary): go along the bottom row from left to right, then move up one row, and repeat. Say we split each side into $5$ intervals. Then the interior points are $u_{i,j}$ with $i, j = 1, \ldots, 4$, so there are $16$ unknowns, and the global number of $u_{i,j}$ is

$$
k = i + 4(j - 1).
$$

So $u_{1,1}, u_{2,1}, u_{3,1}, u_{4,1}$ are numbers $1$ to $4$, then $u_{1,2}$ is number $5$, and so on. Be careful: $i$ (the $x$ index) runs fastest. The point after $u_{1,1}$ is $u_{2,1}$, not $u_{1,2}$. Other orders, such as column by column, are also fine, as long as we know which one we use. The order decides the matrix.

### Filling the matrix

Take one point, find its row, and see which columns its neighbors fall in. Always think of all five terms, and **drop a neighbor that lies on the boundary**: its value is known, so it moves to the right-hand side (or is zero, if the condition is homogeneous).

- A **corner** point like $u_{1,1}$ keeps three terms: itself ($-4$), east, and north. West and south are boundary values.
- An **edge** point like $u_{2,1}$ keeps four: itself, west, east, and north.
- An **inside** point like $u_{2,2}$ keeps all five.

The east and west neighbors are next to each other in the ordering. The north and south neighbors are one whole row of the grid away, which is $4$ positions in our example. The result is a **block tridiagonal** matrix:

$$
A = \frac{1}{h^2}
\begin{bmatrix}
T & I & & \\
I & T & I & \\
& I & T & I \\
& & I & T
\end{bmatrix},
\qquad
T =
\begin{bmatrix}
-4 & 1 & & \\
1 & -4 & 1 & \\
& 1 & -4 & 1 \\
& & 1 & -4
\end{bmatrix},
$$

where $I$ is the $4 \times 4$ identity. Each diagonal block $T$ handles one grid row: the east-west links and the center. It is the 1D matrix with $-4$ instead of $-2$. Each off-diagonal block $I$ handles the north-south links between neighboring grid rows. Only points with the same $i$ are linked, so the block is the identity.

Note the zeros inside the first sub-diagonal: for example, entry $(4, 5)$ is $0$. Point $4$ is $u_{4,1}$ at the right end of the bottom row, and point $5$ is $u_{1,2}$ at the left end of the next row. They are next in the ordering, but they are not neighbors on the grid. The matrix is still sparse: each row has at most $5$ nonzero entries, out of $16$ here, and out of $m^2$ for an $m \times m$ grid of unknowns.

## 6. Experiment

The script uses the test solution $u(x) = e^{2x}$ on $(0, 1)$, so $f = u'' = 4e^{2x}$, $u'(0) = 2$, and $u'(1) = 2e^2$. It uses $n = 4, 8, \ldots, 4096$ intervals. Only the Python standard library is used.

**Part 1: the error formula.** With Dirichlet conditions at both ends, section 2 predicts $E_j \approx h^2 g(x_j)$ with $g(x) = \frac{1}{12}\big(-u''(x) + u''(0) + x(u''(1) - u''(0))\big)$. The measured $\max_j |E_j| / h^2$ goes to $0.50524$, and $\max_j |g(x_j)|$ is also $0.50524$. The largest gap $\max_j |E_j/h^2 - g(x_j)|$ falls from $6.1 \times 10^{-3}$ at $n = 4$ to $3.7 \times 10^{-7}$ at $n = 512$. So the formula predicts not only the size of the error, but its whole shape. For $n \ge 2048$ the gap grows again ($2.4 \times 10^{-4}$ at $n = 4096$), because dividing a tiny error by $h^2$ also magnifies the round-off.

**Part 2: three Neumann rows.** Here $u'(0) = 2$ and $u(1) = e^2$. The chart shows the largest error for each boundary row, with the all-Dirichlet problem as a reference.

<figure class="viz">
<div class="viz-legend"><span class="viz-key"><span class="viz-swatch" style="background:#2a78d6"></span>One-sided, 1st order</span><span class="viz-key"><span class="viz-swatch" style="background:#eb6834"></span>Ghost node</span><span class="viz-key"><span class="viz-swatch" style="background:#1baf7a"></span>One-sided, 2nd order</span><span class="viz-key"><span class="viz-swatch" style="background:#eda100"></span>Dirichlet at both ends</span></div>
<div class="viz-frame">
<svg viewBox="0 0 680 400" role="img" tabindex="0" aria-label="Largest error of the solution as h gets smaller, on log scales, for three ways to handle a Neumann condition at x = 0 and for the all-Dirichlet problem. The first-order boundary row gives a line with slope 1. The other three lines fall with slope 2." class="viz-plot"><line x1="64" x2="520" y1="30.0" y2="30.0" class="viz-grid"/><text x="54" y="34.0" text-anchor="end" class="viz-tick">1</text><line x1="64" x2="520" y1="107.5" y2="107.5" class="viz-grid"/><text x="54" y="111.5" text-anchor="end" class="viz-tick">10<tspan dy="-0.5em" font-size="75%">−2</tspan></text><line x1="64" x2="520" y1="185.0" y2="185.0" class="viz-grid"/><text x="54" y="189.0" text-anchor="end" class="viz-tick">10<tspan dy="-0.5em" font-size="75%">−4</tspan></text><line x1="64" x2="520" y1="262.5" y2="262.5" class="viz-grid"/><text x="54" y="266.5" text-anchor="end" class="viz-tick">10<tspan dy="-0.5em" font-size="75%">−6</tspan></text><line x1="64" x2="520" y1="340.0" y2="340.0" class="viz-grid"/><text x="54" y="344.0" text-anchor="end" class="viz-tick">10<tspan dy="-0.5em" font-size="75%">−8</tspan></text><text x="64.0" y="362" text-anchor="middle" class="viz-tick">2<tspan dy="-0.5em" font-size="75%">−2</tspan></text><text x="155.2" y="362" text-anchor="middle" class="viz-tick">2<tspan dy="-0.5em" font-size="75%">−4</tspan></text><text x="246.4" y="362" text-anchor="middle" class="viz-tick">2<tspan dy="-0.5em" font-size="75%">−6</tspan></text><text x="337.6" y="362" text-anchor="middle" class="viz-tick">2<tspan dy="-0.5em" font-size="75%">−8</tspan></text><text x="428.8" y="362" text-anchor="middle" class="viz-tick">2<tspan dy="-0.5em" font-size="75%">−10</tspan></text><text x="520.0" y="362" text-anchor="middle" class="viz-tick">2<tspan dy="-0.5em" font-size="75%">−12</tspan></text><text x="292.0" y="388" text-anchor="middle" class="viz-axis-title">Step size h (smaller to the right)</text><text x="64" y="14" class="viz-axis-title">Error</text><path d="M64.0,36.7 L109.6,50.6 L155.2,63.6 L200.8,75.9 L246.4,88.0 L292.0,99.8 L337.6,111.6 L383.2,123.3 L428.8,135.0 L474.4,146.6 L520.0,158.3" fill="none" stroke="#2a78d6" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><circle cx="64.0" cy="36.7" r="3" fill="#2a78d6" stroke="#ffffff" stroke-width="1.5"/><circle cx="109.6" cy="50.6" r="3" fill="#2a78d6" stroke="#ffffff" stroke-width="1.5"/><circle cx="155.2" cy="63.6" r="3" fill="#2a78d6" stroke="#ffffff" stroke-width="1.5"/><circle cx="200.8" cy="75.9" r="3" fill="#2a78d6" stroke="#ffffff" stroke-width="1.5"/><circle cx="246.4" cy="88.0" r="3" fill="#2a78d6" stroke="#ffffff" stroke-width="1.5"/><circle cx="292.0" cy="99.8" r="3" fill="#2a78d6" stroke="#ffffff" stroke-width="1.5"/><circle cx="337.6" cy="111.6" r="3" fill="#2a78d6" stroke="#ffffff" stroke-width="1.5"/><circle cx="383.2" cy="123.3" r="3" fill="#2a78d6" stroke="#ffffff" stroke-width="1.5"/><circle cx="428.8" cy="135.0" r="3" fill="#2a78d6" stroke="#ffffff" stroke-width="1.5"/><circle cx="474.4" cy="146.6" r="3" fill="#2a78d6" stroke="#ffffff" stroke-width="1.5"/><circle cx="520.0" cy="158.3" r="3" fill="#2a78d6" stroke="#ffffff" stroke-width="1.5"/><path d="M64.0,59.5 L109.6,82.7 L155.2,106.0 L200.8,129.3 L246.4,152.7 L292.0,176.0 L337.6,199.3 L383.2,222.7 L428.8,246.0 L474.4,269.3 L520.0,292.7" fill="none" stroke="#eb6834" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><circle cx="64.0" cy="59.5" r="3" fill="#eb6834" stroke="#ffffff" stroke-width="1.5"/><circle cx="109.6" cy="82.7" r="3" fill="#eb6834" stroke="#ffffff" stroke-width="1.5"/><circle cx="155.2" cy="106.0" r="3" fill="#eb6834" stroke="#ffffff" stroke-width="1.5"/><circle cx="200.8" cy="129.3" r="3" fill="#eb6834" stroke="#ffffff" stroke-width="1.5"/><circle cx="246.4" cy="152.7" r="3" fill="#eb6834" stroke="#ffffff" stroke-width="1.5"/><circle cx="292.0" cy="176.0" r="3" fill="#eb6834" stroke="#ffffff" stroke-width="1.5"/><circle cx="337.6" cy="199.3" r="3" fill="#eb6834" stroke="#ffffff" stroke-width="1.5"/><circle cx="383.2" cy="222.7" r="3" fill="#eb6834" stroke="#ffffff" stroke-width="1.5"/><circle cx="428.8" cy="246.0" r="3" fill="#eb6834" stroke="#ffffff" stroke-width="1.5"/><circle cx="474.4" cy="269.3" r="3" fill="#eb6834" stroke="#ffffff" stroke-width="1.5"/><circle cx="520.0" cy="292.7" r="3" fill="#eb6834" stroke="#ffffff" stroke-width="1.5"/><path d="M64.0,61.8 L109.6,90.5 L155.2,116.9 L200.8,141.8 L246.4,166.0 L292.0,189.8 L337.6,213.3 L383.2,236.7 L428.8,260.1 L474.4,283.5 L520.0,306.8" fill="none" stroke="#1baf7a" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><circle cx="64.0" cy="61.8" r="3" fill="#1baf7a" stroke="#ffffff" stroke-width="1.5"/><circle cx="109.6" cy="90.5" r="3" fill="#1baf7a" stroke="#ffffff" stroke-width="1.5"/><circle cx="155.2" cy="116.9" r="3" fill="#1baf7a" stroke="#ffffff" stroke-width="1.5"/><circle cx="200.8" cy="141.8" r="3" fill="#1baf7a" stroke="#ffffff" stroke-width="1.5"/><circle cx="246.4" cy="166.0" r="3" fill="#1baf7a" stroke="#ffffff" stroke-width="1.5"/><circle cx="292.0" cy="189.8" r="3" fill="#1baf7a" stroke="#ffffff" stroke-width="1.5"/><circle cx="337.6" cy="213.3" r="3" fill="#1baf7a" stroke="#ffffff" stroke-width="1.5"/><circle cx="383.2" cy="236.7" r="3" fill="#1baf7a" stroke="#ffffff" stroke-width="1.5"/><circle cx="428.8" cy="260.1" r="3" fill="#1baf7a" stroke="#ffffff" stroke-width="1.5"/><circle cx="474.4" cy="283.5" r="3" fill="#1baf7a" stroke="#ffffff" stroke-width="1.5"/><circle cx="520.0" cy="306.8" r="3" fill="#1baf7a" stroke="#ffffff" stroke-width="1.5"/><path d="M64.0,88.8 L109.6,111.7 L155.2,134.8 L200.8,158.2 L246.4,181.5 L292.0,204.8 L337.6,228.1 L383.2,251.5 L428.8,274.8 L474.4,298.1 L520.0,321.4" fill="none" stroke="#eda100" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><circle cx="64.0" cy="88.8" r="3" fill="#eda100" stroke="#ffffff" stroke-width="1.5"/><circle cx="109.6" cy="111.7" r="3" fill="#eda100" stroke="#ffffff" stroke-width="1.5"/><circle cx="155.2" cy="134.8" r="3" fill="#eda100" stroke="#ffffff" stroke-width="1.5"/><circle cx="200.8" cy="158.2" r="3" fill="#eda100" stroke="#ffffff" stroke-width="1.5"/><circle cx="246.4" cy="181.5" r="3" fill="#eda100" stroke="#ffffff" stroke-width="1.5"/><circle cx="292.0" cy="204.8" r="3" fill="#eda100" stroke="#ffffff" stroke-width="1.5"/><circle cx="337.6" cy="228.1" r="3" fill="#eda100" stroke="#ffffff" stroke-width="1.5"/><circle cx="383.2" cy="251.5" r="3" fill="#eda100" stroke="#ffffff" stroke-width="1.5"/><circle cx="428.8" cy="274.8" r="3" fill="#eda100" stroke="#ffffff" stroke-width="1.5"/><circle cx="474.4" cy="298.1" r="3" fill="#eda100" stroke="#ffffff" stroke-width="1.5"/><circle cx="520.0" cy="321.4" r="3" fill="#eda100" stroke="#ffffff" stroke-width="1.5"/><text x="530" y="162.3" class="viz-label">One-sided, 1st order</text><text x="530" y="290.7" class="viz-label">Ghost node</text><text x="530" y="312.8" class="viz-label">One-sided, 2nd order</text><text x="530" y="334.4" class="viz-label">Dirichlet at both ends</text><line x1="0" x2="0" y1="30" y2="340" class="viz-crosshair" visibility="hidden"/></svg>
<div class="viz-tooltip" hidden></div>
</div>
<figcaption>Largest error max |U<sub>j</sub> − u(x<sub>j</sub>)| for u″ = f on (0, 1) with u = e<sup>2x</sup>, on log scales. The first three lines use a Neumann condition u′(0) = 2 with three different boundary rows. The fourth uses Dirichlet conditions at both ends. Hover over the chart, or focus it and use the arrow keys, to read the values.</figcaption>
<details class="viz-table"><summary>Show the numbers</summary><table><thead><tr><th>h</th><th>One-sided, 1st order</th><th>Ghost node</th><th>One-sided, 2nd order</th><th>Dirichlet at both ends</th></tr></thead><tbody><tr><td>1/4</td><td>6.7e-1</td><td>1.7e-1</td><td>1.5e-1</td><td>3.0e-2</td></tr><tr><td>1/8</td><td>2.9e-1</td><td>4.4e-2</td><td>2.7e-2</td><td>7.8e-3</td></tr><tr><td>1/16</td><td>1.4e-1</td><td>1.1e-2</td><td>5.7e-3</td><td>2.0e-3</td></tr><tr><td>1/32</td><td>6.5e-2</td><td>2.7e-3</td><td>1.3e-3</td><td>4.9e-4</td></tr><tr><td>1/64</td><td>3.2e-2</td><td>6.8e-4</td><td>3.1e-4</td><td>1.2e-4</td></tr><tr><td>1/128</td><td>1.6e-2</td><td>1.7e-4</td><td>7.5e-5</td><td>3.1e-5</td></tr><tr><td>1/256</td><td>7.9e-3</td><td>4.3e-5</td><td>1.9e-5</td><td>7.7e-6</td></tr><tr><td>1/512</td><td>3.9e-3</td><td>1.1e-5</td><td>4.6e-6</td><td>1.9e-6</td></tr><tr><td>1/1024</td><td>2.0e-3</td><td>2.7e-6</td><td>1.2e-6</td><td>4.8e-7</td></tr><tr><td>1/2048</td><td>9.8e-4</td><td>6.7e-7</td><td>2.9e-7</td><td>1.2e-7</td></tr><tr><td>1/4096</td><td>4.9e-4</td><td>1.7e-7</td><td>7.2e-8</td><td>3.0e-8</td></tr></tbody></table></details>
<script>(() => {const fig = document.currentScript.closest('figure');const svg = fig.querySelector('svg');const line = svg.querySelector('.viz-crosshair');const tip = fig.querySelector('.viz-tooltip');const names = ["One-sided, 1st order", "Ghost node", "One-sided, 2nd order", "Dirichlet at both ends"]; const colors = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"];const rows = [["6.7e-1", "1.7e-1", "1.5e-1", "3.0e-2"], ["2.9e-1", "4.4e-2", "2.7e-2", "7.8e-3"], ["1.4e-1", "1.1e-2", "5.7e-3", "2.0e-3"], ["6.5e-2", "2.7e-3", "1.3e-3", "4.9e-4"], ["3.2e-2", "6.8e-4", "3.1e-4", "1.2e-4"], ["1.6e-2", "1.7e-4", "7.5e-5", "3.1e-5"], ["7.9e-3", "4.3e-5", "1.9e-5", "7.7e-6"], ["3.9e-3", "1.1e-5", "4.6e-6", "1.9e-6"], ["2.0e-3", "2.7e-6", "1.2e-6", "4.8e-7"], ["9.8e-4", "6.7e-7", "2.9e-7", "1.2e-7"], ["4.9e-4", "1.7e-7", "7.2e-8", "3.0e-8"]]; const heads = ["h = 1/4", "h = 1/8", "h = 1/16", "h = 1/32", "h = 1/64", "h = 1/128", "h = 1/256", "h = 1/512", "h = 1/1024", "h = 1/2048", "h = 1/4096"];const left = 64, right = 520, steps = 10, width = 680;let current = -1;const show = (n) => { n = Math.max(0, Math.min(steps, n)); current = n; const x = left + (right - left) * n / steps; line.setAttribute('x1', x); line.setAttribute('x2', x); line.setAttribute('visibility', 'visible'); tip.replaceChildren(); const head = document.createElement('div'); head.className = 'viz-tooltip-head'; head.textContent = heads[n]; tip.append(head); names.forEach((name, i) => {  const row = document.createElement('div'); row.className = 'viz-tooltip-row';  const key = document.createElement('span'); key.className = 'viz-swatch'; key.style.background = colors[i];  const value = document.createElement('strong'); value.textContent = rows[n][i];  const label = document.createElement('span'); label.textContent = name;  row.append(key, value, label); tip.append(row); }); const box = svg.getBoundingClientRect(), frame = fig.getBoundingClientRect(); const xs = box.left - frame.left + x * box.width / width; tip.hidden = false; const room = frame.width - tip.offsetWidth - 8; tip.style.left = Math.max(0, Math.min(room, xs + 12)) + 'px'; tip.style.top = (box.top - frame.top + 8) + 'px';};const hide = () => { line.setAttribute('visibility', 'hidden'); tip.hidden = true; current = -1; };svg.addEventListener('pointermove', (event) => { const box = svg.getBoundingClientRect(); const x = (event.clientX - box.left) * width / box.width; show(Math.round((x - left) / (right - left) * steps));});svg.addEventListener('pointerleave', hide);svg.addEventListener('focus', () => show(current < 0 ? 0 : current));svg.addEventListener('blur', hide);svg.addEventListener('keydown', (event) => { if (event.key === 'ArrowRight') { show(current + 1); event.preventDefault(); } if (event.key === 'ArrowLeft') { show(current - 1); event.preventDefault(); } if (event.key === 'Escape') hide();});})();</script>
</figure>

| $n$ | one-sided, 1st order | $E_0 / h$ | ghost node | one-sided, 2nd order | Dirichlet (reference) |
|---|---|---|---|---|---|
| 4 | $6.73 \times 10^{-1}$ | $2.69$ | $1.73 \times 10^{-1}$ | $1.51 \times 10^{-1}$ | $3.04 \times 10^{-2}$ |
| 16 | $1.36 \times 10^{-1}$ | $2.17$ | $1.09 \times 10^{-2}$ | $5.73 \times 10^{-3}$ | $1.97 \times 10^{-3}$ |
| 64 | $3.19 \times 10^{-2}$ | $2.04$ | $6.83 \times 10^{-4}$ | $3.09 \times 10^{-4}$ | $1.23 \times 10^{-4}$ |
| 256 | $7.86 \times 10^{-3}$ | $2.01$ | $4.27 \times 10^{-5}$ | $1.86 \times 10^{-5}$ | $7.71 \times 10^{-6}$ |
| 1024 | $1.96 \times 10^{-3}$ | $2.003$ | $2.67 \times 10^{-6}$ | $1.15 \times 10^{-6}$ | $4.82 \times 10^{-7}$ |
| 4096 | $4.88 \times 10^{-4}$ | $2.0007$ | $1.67 \times 10^{-7}$ | $7.18 \times 10^{-8}$ | $3.01 \times 10^{-8}$ |

The measured order from $n = 2048$ to $4096$ is $1.000$ for the first-order row, and $2.000$, $2.001$, and $1.999$ for the other three.

What I noticed:

- **The first-order row gives first order everywhere.** The error at $x = 0$ divided by $h$ goes to $2$, which is $\frac{1}{2} u''(0) = \frac{1}{2} \cdot 4$, as section 3 predicts. At $x = 1/2$ the error is half of that ($2.44 \times 10^{-4}$ at $n = 4096$, against $4.88 \times 10^{-4}$ at $x = 0$), which matches the factor $(1 - x)$. The middle of the domain is far from the bad row, but it is still polluted.
- **Changing one number on the right-hand side fixes it.** The ghost node row has the same matrix as the first-order row. At $n = 4096$ its error is about $2900$ times smaller.
- The one-sided second-order row is also second order. For this $u$ its error is a bit smaller than the ghost node error, but both are about $2$ to $6$ times larger than the all-Dirichlet error. Knowing only the slope at one end costs some accuracy, but not the order.
- As a check, I put the computed solution of the third approach back into its original three-point row $(-3U_0 + 4U_1 - U_2)/(2h) = \sigma$. It holds exactly ($0.0$ at $n = 64$).

**Part 3: pure Neumann.** The data $\sigma_0 = 2$, $\sigma_1 = 2e^2$, $f = 4e^{2x}$ satisfy the continuous compatibility condition exactly. First, for $n = 4$ the five row sums of the matrix are all $0.0$, so $A$ is singular. Next, the discrete compatibility condition from section 4 fails by a small amount:

| $n$ | trapezoid sum $- (\sigma_1 - \sigma_0)$ | error after the fix |
|---|---|---|
| 4 | $2.65 \times 10^{-1}$ | $3.06 \times 10^{-1}$ |
| 16 | $1.66 \times 10^{-2}$ | $1.92 \times 10^{-2}$ |
| 64 | $1.04 \times 10^{-3}$ | $1.20 \times 10^{-3}$ |
| 256 | $6.50 \times 10^{-5}$ | $7.52 \times 10^{-5}$ |
| 1024 | $4.06 \times 10^{-6}$ | $4.70 \times 10^{-6}$ |
| 4096 | $2.54 \times 10^{-7}$ | $2.94 \times 10^{-7}$ |

The mismatch shrinks with order $2.00$. It is the error of the trapezoid rule, about $\frac{h^2}{12}\big(f'(1) - f'(0)\big)$. So, strictly, the discrete system has **no solution** for any $n$, even though the continuous problem has one. To see this, the script fixes $U_0$, solves rows $0$ to $n - 1$ one after another, and checks the last row. The last row is always off by exactly the mismatch (with the opposite sign).

A simple fix is to shift $f$ by a constant equal to the mismatch. A constant has a trapezoid sum of exactly itself on $[0, 1]$, so the shifted data are discretely compatible. After the shift, the last row holds up to round-off (below $3 \times 10^{-11}$). The solution is still unique only up to a constant. To compare with $u$, the script picks the constant by setting $U_0 = u(0)$. The error is then second order ($2.00$), as shown in the table. The shift is $O(h^2)$, so it does not hurt the order.

The full script and the chart code are in [this site's repository](https://github.com/doubleyands/doubleyands.github.io/tree/main/scripts/posts/numerical-pde/lecture-3).

## Summary

1. Round-off is small per operation but can pile up in large systems. **Scale the rows** so the coefficients have similar sizes.
2. The global error and the local truncation error are linked by $AE = -\tau$. **The error solves the same kind of boundary value problem**, with $-\tau$ as the source. For the three-point scheme, this gives a second-order error with a predictable shape.
3. A Neumann condition makes $U_0$ an unknown. A one-sided first-order boundary row makes the error **first order everywhere**: this is the pollution effect.
4. A **ghost node** gives a second-order row with the same matrix. Only the right-hand side changes, by $\frac{h}{2} f(x_0)$. A one-sided second-order row also works, but it leaves the tridiagonal band.
5. With **pure Neumann** conditions, solutions are unique only up to a constant, and they exist only if $\int_0^1 f\, dx = \sigma_1 - \sigma_0$. The matrix is singular, and the discrete condition uses the trapezoid rule in place of the integral.
6. In 2D, the **5-point stencil** is two 1D stencils added together. With lexicographic ordering it gives a sparse **block tridiagonal** matrix with blocks $T$ and $I$.
