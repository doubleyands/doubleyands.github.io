---
title: "Numerical Analysis, Week 2: Nonlinear Systems and Polynomial Interpolation"
description: "Newton's method for systems, a short look at optimization, and polynomial interpolation: Lagrange, divided differences, Runge's problem, and Hermite. Tested on two functions."
pubDate: 2026-10-10
project: numerical-analysis
tags: [numerical analysis, interpolation, newton's method, notes]
---

This is my summary of the second week of a numerical analysis course. As before, I write the ideas in my own words and check one of them with a small experiment. The course follows K. E. Atkinson, *An Introduction to Numerical Analysis*.

The week had two parts. The first part finishes root finding. In [Week 1](/blog/2026-10-03-numerical-analysis-week-1/) we solved one equation in one unknown. Now we have many equations and many unknowns. The second part starts a new topic. We do not know a function. We only know its values at a few points. How do we build a function from that data?

## 1. Systems of nonlinear equations

### Why we need something new

A linear system is easy to handle. We can write it as $A\mathbf{x} = \mathbf{b}$, and there are many good methods for that. Now look at this system:

$$
x^2 + y^2 = 4, \qquad xy = 1.
$$

We cannot write it as $A\mathbf{x} = \mathbf{b}$, because $x^2$ and $xy$ are not linear. So we need one more step first. We replace the problem by a linear one near our current guess. This step is called **linearization**. It is not exact, so we repeat it until we are close enough.

We write a general system of $m$ equations as $\mathbf{f}(\mathbf{x}) = \mathbf{0}$, with $\mathbf{x} \in \mathbb{R}^m$. The true solution is $\boldsymbol{\alpha}$. I use two unknowns below to keep the formulas short. The ideas work for any $m$.

### Approach 1: fixed-point iteration for vectors

As in one dimension, we rewrite the system as $\mathbf{x} = \mathbf{g}(\mathbf{x})$. One simple way is $g_1(x, y) = x + f_1(x, y)$ and $g_2(x, y) = y + f_2(x, y)$. Then we iterate:

$$
\mathbf{x}_{n+1} = \mathbf{g}(\mathbf{x}_n).
$$

To study the error, we apply the mean value theorem to each component. For $g_1$,

$$
\alpha_1 - x_{n+1} = g_1(\boldsymbol{\alpha}) - g_1(\mathbf{x}_n) = \frac{\partial g_1}{\partial x}(\boldsymbol{\xi}^1_n)\,(\alpha_1 - x_n) + \frac{\partial g_1}{\partial y}(\boldsymbol{\xi}^1_n)\,(\alpha_2 - y_n).
$$

Here $\boldsymbol{\xi}^1_n$ is a point on the segment between $\mathbf{x}_n$ and $\boldsymbol{\alpha}$. The same works for $g_2$, with another point $\boldsymbol{\xi}^2_n$. The two points are usually different. Together the two lines give

$$
\boldsymbol{\alpha} - \mathbf{x}_{n+1} = G_n\,(\boldsymbol{\alpha} - \mathbf{x}_n), \qquad
G(\mathbf{x}) = \begin{bmatrix} \partial g_1/\partial x & \partial g_1/\partial y \\ \partial g_2/\partial x & \partial g_2/\partial y \end{bmatrix}.
$$

$G$ is the Jacobian matrix of $\mathbf{g}$. Row 1 of $G_n$ is taken at $\boldsymbol{\xi}^1_n$, and row 2 at $\boldsymbol{\xi}^2_n$. In one dimension the error was multiplied by $g'(\xi_n)$. Here it is multiplied by a matrix.

**Theorem.** Let $D \subset \mathbb{R}^2$ be closed, bounded, and convex, and let $\mathbf{g}$ be continuously differentiable on $D$. Suppose

1. $\mathbf{g}(D) \subset D$, and
2. $\lambda = \max_{\mathbf{x} \in D} \|G(\mathbf{x})\|_\infty < 1$.

Then $\mathbf{x} = \mathbf{g}(\mathbf{x})$ has exactly one solution $\boldsymbol{\alpha}$ in $D$, and the iteration converges to it from any starting point in $D$. Near the end, each step multiplies the error by at most about $\|G(\boldsymbol{\alpha})\|_\infty$.

Here $\|A\|_\infty$ is the largest row sum of $|a_{ij}|$. These are the same two conditions as in one dimension: the map stays inside $D$, and it is a contraction. The size of the system does not matter.

Why do we need these conditions on $D$?

- **Convex**: the mean value theorem uses points on the segment between $\mathbf{x}_n$ and $\boldsymbol{\alpha}$. That segment must stay inside $D$, because we only know the bound on $G$ inside $D$.
- **Closed and bounded**: a continuous function on a closed and bounded set reaches its maximum. So the maximum in condition 2 really exists, and $\lambda$ is a fixed number below $1$.

In practice, condition 2 is very hard to check for a real system. So this method is not used much as a system solver. Newton's method is more common.

### Approach 2: Newton's method for systems

Expand each $f_i$ around the current guess $\mathbf{x}_0$ with Taylor's theorem in two variables:

$$
0 = f_i(\boldsymbol{\alpha}) = f_i(\mathbf{x}_0) + \frac{\partial f_i}{\partial x}(\mathbf{x}_0)\,(\alpha_1 - x_0) + \frac{\partial f_i}{\partial y}(\mathbf{x}_0)\,(\alpha_2 - y_0) + (\text{second-order terms}).
$$

Drop the second-order terms. In matrix form,

$$
\mathbf{0} \approx \mathbf{f}(\mathbf{x}_0) + F(\mathbf{x}_0)\,(\boldsymbol{\alpha} - \mathbf{x}_0).
$$

Here $F$ is the Jacobian matrix of $\mathbf{f}$. Everything here is known except $\boldsymbol{\alpha}$. So this is a linear system for $\boldsymbol{\alpha}$. Its solution is not the true $\boldsymbol{\alpha}$, because we dropped terms. But we hope it is closer than $\mathbf{x}_0$, and we call it $\mathbf{x}_1$. Repeating gives

$$
\mathbf{x}_{n+1} = \mathbf{x}_n - F(\mathbf{x}_n)^{-1}\,\mathbf{f}(\mathbf{x}_n).
$$

The dropped terms are the linearization error. They are small only when $\boldsymbol{\alpha} - \mathbf{x}_n$ is small. So, just like in one dimension, Newton's method for systems needs a good starting guess.

**Do not compute the inverse.** The formula has $F^{-1}$, but we never form it. Instead we solve one linear system for a correction $\boldsymbol{\delta}$:

$$
F(\mathbf{x}_n)\,\boldsymbol{\delta}_{n+1} = -\mathbf{f}(\mathbf{x}_n), \qquad \mathbf{x}_{n+1} = \mathbf{x}_n + \boldsymbol{\delta}_{n+1}.
$$

"Find $A^{-1}$" and "solve $A\mathbf{x} = \mathbf{b}$" are different jobs, and the second one is much cheaper. This matters when there are thousands of unknowns. How fast we can solve a large linear system depends most on how sparse the matrix is, that is, how many entries are zero. A sparse system with millions of unknowns can be easy. A dense one of the same size can be almost impossible. Jacobian matrices are often dense, so each Newton step can be expensive.

When it works, the convergence is still quadratic:

$$
\|\boldsymbol{\alpha} - \mathbf{x}_{n+1}\|_\infty \le B\,\|\boldsymbol{\alpha} - \mathbf{x}_n\|_\infty^2.
$$

**Good and bad points.** Newton's method is simple and works for many problems. If we cannot compute the derivatives, we can replace them with difference quotients. The weak point is cost: we build and solve with a Jacobian at every step. So there are many variants. Two simple ones:

- **Reuse the Jacobian.** Build $F$ once and use it for the next $k$ steps, for example $k = 10$. This saves work. But convergence is not guaranteed, so we use it only when $F$ changes slowly.
- **Shorten the step (damped Newton).** Keep the Newton direction $\mathbf{d}_n = -F(\mathbf{x}_n)^{-1}\mathbf{f}(\mathbf{x}_n)$, but take $\mathbf{x}_{n+1} = \mathbf{x}_n + s\,\mathbf{d}_n$. We choose $s > 0$ so that $\|\mathbf{f}(\mathbf{x}_n + s\,\mathbf{d}_n)\|_2^2$ is as small as possible. $s = 1$ is the plain Newton step. The direction is good near the root, but far from it a full step can jump too far and land somewhere worse. A smaller step makes steady progress.

## 2. A short note on optimization

To find a minimum of a smooth function $f(\mathbf{x})$, we look for points where the gradient is zero:

$$
\nabla f(\mathbf{x}) = \mathbf{0}.
$$

This is a system of nonlinear equations. So we can use Newton's method from above. The Jacobian of $\nabla f$ is the **Hessian** $H$, with $H_{ij} = \partial^2 f / \partial x_i \partial x_j$:

$$
\mathbf{x}_{n+1} = \mathbf{x}_n - H(\mathbf{x}_n)^{-1}\,\nabla f(\mathbf{x}_n).
$$

A zero of the gradient is usually a *local* minimum (or a saddle, or a maximum). It may not be the global minimum. The main message of this section is short: optimization problems can be solved with Newton's method.

## 3. Polynomial interpolation

### The question changes

Until now we knew $f$. Often we do not. We only have data: at $x_0$ we measured $y_0$, at $x_1$ we measured $y_1$, and so on. We want a function that passes through these points. The first choice is a polynomial.

**Theorem.** Given $n+1$ distinct points $x_0, \dots, x_n$ and values $y_0, \dots, y_n$, there is exactly one polynomial $p$ of degree at most $n$ with $p(x_i) = y_i$ for all $i$.

The count makes sense: a polynomial of degree $n$ has $n+1$ coefficients, and we have $n+1$ conditions.

### The direct way and its problems

Write $p(x) = a_0 + a_1 x + \dots + a_n x^n$ and impose the conditions. This gives a linear system with the **Vandermonde matrix**. Its row $i$ is $(1, x_i, x_i^2, \dots, x_i^n)$. This works, but it has two problems:

- **Cost.** We must solve an $(n+1) \times (n+1)$ system. That is fine for 10 points, but not for a million.
- **Ill-conditioning.** If two nodes are close, their rows are almost the same. The matrix is then close to singular, and small errors in the data cause large errors in the coefficients. Vandermonde matrices are known for this. (The condition number measures it. The course will cover it later.)

So we know the polynomial exists, but finding its coefficients this way is expensive. The next two forms avoid the matrix.

### Lagrange form

$$
p_n(x) = \sum_{i=0}^{n} y_i\,\ell_i(x), \qquad \ell_i(x) = \prod_{j \ne i} \frac{x - x_j}{x_i - x_j}.
$$

Each $\ell_i$ is a polynomial of degree $n$. It equals $1$ at its own node $x_i$ and $0$ at every other node. So at $x = x_k$ only one term survives, and $p_n(x_k) = y_k$. No matrix and no inverse are needed.

### The interpolation error

**Theorem.** Let $I$ be an interval that contains $x_0, \dots, x_n$ and $x$, and let $f \in C^{n+1}(I)$. Let $p_n$ interpolate $f$ at the nodes. Then for some $\xi \in I$,

$$
f(x) - p_n(x) = \frac{(x - x_0)(x - x_1)\cdots(x - x_n)}{(n+1)!}\, f^{(n+1)}(\xi).
$$

It looks like the Taylor remainder, but there is a big difference. In Taylor's theorem everything happens near one point. Here every node appears in the product.

**Proof.** Fix $x$ (not a node). Let $W(t) = (t - x_0)\cdots(t - x_n)$, a polynomial of degree $n+1$. Define

$$
\phi(t) = f(t) - p_n(t) - \frac{f(x) - p_n(x)}{W(x)}\,W(t).
$$

Count the zeros of $\phi$. At each node $t = x_i$, both $f - p_n$ and $W$ are zero, so $\phi(x_i) = 0$. That is $n+1$ zeros. At $t = x$, the last term equals $f(x) - p_n(x)$, so $\phi(x) = 0$ too. In total $\phi$ has $n+2$ zeros.

Between two zeros of a smooth function there is a zero of its derivative (Rolle's theorem, a special case of the mean value theorem). So $\phi'$ has at least $n+1$ zeros, $\phi''$ at least $n$, and $\phi^{(n+1)}$ at least one. Call it $\xi$.

Now differentiate $n+1$ times. $p_n$ has degree $n$, so $p_n^{(n+1)} = 0$. $W$ has leading term $t^{n+1}$, so $W^{(n+1)} = (n+1)!$. At $t = \xi$:

$$
0 = f^{(n+1)}(\xi) - \frac{f(x) - p_n(x)}{W(x)}\,(n+1)!,
$$

and this is the formula. $\blacksquare$

**An example.** Interpolate $f(x) = \sin x$ by a straight line on $[x_0, x_1]$ with $h = x_1 - x_0$. For $x$ between the nodes,

$$
|f(x) - p_1(x)| = \frac{(x - x_0)(x_1 - x)}{2}\,|\sin\xi|.
$$

The product $(x - x_0)(x_1 - x)$ is largest at the midpoint. Its value there is $h^2/4$. So the error is at most $\frac{h^2}{8}\max|\sin\xi|$. On $[1.0, 1.1]$ this gives the upper bound $\frac{0.01}{8}\sin 1.1 = 1.114 \times 10^{-3}$. With $\sin 1.0$ instead, we get a lower bound $1.052 \times 10^{-3}$ for the largest error. I measured the largest error on a fine grid: $1.084 \times 10^{-3}$. It is between the two, as it should be.

**Only inside.** The theorem needs $x$ in an interval that contains the nodes. If $x$ is far outside, the product $W(x)$ grows fast, and the formula gives no useful bound. Using $p_n$ outside the data is *extrapolation*, and it is much less safe than *interpolation*.

## 4. Newton's divided differences

### Adding one more point

Suppose we get one new data point. With the Lagrange form, every $\ell_i$ changes, so we must start again. Can we keep the old polynomial and just add a correction?

Let $p_{n-1}$ interpolate at $x_0, \dots, x_{n-1}$, and look for

$$
p_n(x) = p_{n-1}(x) + C(x).
$$

$C$ must have degree $n$. At the old nodes, $C(x_i) = p_n(x_i) - p_{n-1}(x_i) = f(x_i) - f(x_i) = 0$. So $C$ has the $n$ zeros $x_0, \dots, x_{n-1}$, and it must be

$$
C(x) = a_n (x - x_0)\cdots(x - x_{n-1}).
$$

The last condition $p_n(x_n) = f(x_n)$ fixes $a_n$. We call this number the **divided difference** of order $n$:

$$
a_n = f[x_0, x_1, \dots, x_n] = \frac{f(x_n) - p_{n-1}(x_n)}{(x_n - x_0)\cdots(x_n - x_{n-1})}.
$$

Repeating this from one node upward gives **Newton's form**:

$$
p_n(x) = f(x_0) + (x - x_0) f[x_0, x_1] + (x - x_0)(x - x_1) f[x_0, x_1, x_2] + \dots + (x - x_0)\cdots(x - x_{n-1}) f[x_0, \dots, x_n].
$$

### Two properties

**$f[x_0, \dots, x_n]$ is the leading coefficient.** In $p_n = p_{n-1} + C$, only $C$ has an $x^n$ term, and its coefficient is $a_n$. So the divided difference is the coefficient of $x^n$ in the interpolating polynomial. From this we get:

- **Order does not matter.** The interpolating polynomial depends only on the *set* of nodes (it is unique). So its leading coefficient does not change when we reorder the nodes.
- **A recursion.**

$$
f[x_0, \dots, x_n] = \frac{f[x_1, \dots, x_n] - f[x_0, \dots, x_{n-1}]}{x_n - x_0}.
$$

The recursion is the useful part. Using the definition, we would need to build $p_{n-1}$ and evaluate it. With the recursion, each new entry costs one subtraction and one division. So we compute divided differences in a table. We start with the column of values, and each new column comes from the one before:

| $x_i$ | $f(x_i)$ | 1st order | 2nd order | 3rd order |
|---|---|---|---|---|
| $x_0$ | $f(x_0)$ | $f[x_0, x_1]$ | $f[x_0, x_1, x_2]$ | $f[x_0, \dots, x_3]$ |
| $x_1$ | $f(x_1)$ | $f[x_1, x_2]$ | $f[x_1, x_2, x_3]$ | |
| $x_2$ | $f(x_2)$ | $f[x_2, x_3]$ | | |
| $x_3$ | $f(x_3)$ | | | |

Newton's form uses the top row. When a new point arrives, we add one row at the bottom and compute one new diagonal.

**A small example.** Take $f(x) = x^3$ at $0, 1, 2, 3$:

| $x_i$ | $f(x_i)$ | 1st order | 2nd order | 3rd order |
|---|---|---|---|---|
| 0 | 0 | 1 | 3 | 1 |
| 1 | 1 | 7 | 6 | |
| 2 | 8 | 19 | | |
| 3 | 27 | | | |

So $p_3(x) = 0 + x\cdot 1 + x(x-1)\cdot 3 + x(x-1)(x-2)\cdot 1$. Expanding gives $x^3$, as it must, because $f$ is already a cubic. The last entry is $1$, the leading coefficient of $x^3$.

### Every method gives the same polynomial

Let $r_1$ and $r_2$ both interpolate the same data with degree at most $n$, for example one from Lagrange and one from Newton. Then $R = r_1 - r_2$ has degree at most $n$ and $n+1$ zeros $x_0, \dots, x_n$. A nonzero polynomial of degree at most $n$ has at most $n$ zeros. So $R = 0$.

This means the Vandermonde way, the Lagrange form, and Newton's form all give the same polynomial. So they also share the same error formula. We pick the form that is most convenient.

### Divided differences are derivatives

Add one more point $t$ and interpolate at $x_0, \dots, x_n, t$. By Newton's form,

$$
p_{n+1}(x) = p_n(x) + (x - x_0)\cdots(x - x_n)\, f[x_0, \dots, x_n, t].
$$

Put $x = t$. Since $p_{n+1}(t) = f(t)$, we get $f(t) - p_n(t) = (t - x_0)\cdots(t - x_n)\, f[x_0, \dots, x_n, t]$. Compare this with the error formula. The products cancel, and

$$
f[x_0, \dots, x_n, t] = \frac{f^{(n+1)}(\xi)}{(n+1)!}.
$$

In words: a divided difference of $k+1$ points is a $k$-th derivative divided by $k!$, at some point in between. For $k = 1$ this is just the mean value theorem: $f[x_0, x_1] = \frac{f(x_1) - f(x_0)}{x_1 - x_0} = f'(\xi)$. So when the nodes are close together, Newton's form behaves like a Taylor polynomial.

Check with the $x^3$ table. The third-order entry is $1 = f'''/3! = 6/6$. The second-order entries are $f''(\xi)/2 = 3\xi$. They are $3$ and $6$, so $\xi = 1$ for the nodes $0, 1, 2$ and $\xi = 2$ for the nodes $1, 2, 3$. Both points lie inside their node intervals.

## 5. More points do not always help

Interpolate $f$ on $[a, b]$ at $n+1$ equally spaced points. Does $\max_{a \le x \le b} |f(x) - p_n(x)| \to 0$ as $n \to \infty$?

**Not necessarily.** A smooth function can be a counterexample. Take

$$
f(x) = \frac{1}{1 + 16x^2}, \qquad -1 \le x \le 1.
$$

It has derivatives of every order. But with equally spaced nodes, the error near the two ends grows without bound as $n$ grows. This is called the **Runge phenomenon**. My experiment below shows it.

Why at the ends? Inside the interval, each point has data on both sides, and the data holds the polynomial in place. Near an end, there is data on only one side. Nothing holds the polynomial on the outer side, so it swings.

Three ideas to fix it:

1. **Give data outside the interval.** This would help, but often we cannot measure outside.
2. **Put more nodes near the ends.** This reduces the swings. How much it helps depends on how we place the nodes. In my experiment, one common choice removes the problem for this $f$.
3. **Also give derivative values.** If the polynomial must match the slope at the nodes, it has less freedom to swing. This leads to Hermite interpolation.

## 6. Hermite interpolation

Suppose we know both values and slopes at $n$ nodes:

$$
p(x_i) = y_i, \qquad p'(x_i) = y_i', \qquad i = 1, \dots, n.
$$

That is $2n$ conditions, so we look for a polynomial of degree at most $2n - 1$. Such a polynomial exists, and there is an explicit formula built from the Lagrange basis. I skip the formula here, because its derivation is long and not very instructive. The main point is this: if we know derivatives too, there is an interpolation method that uses them.

**Uniqueness.** Suppose $H$ and $G$ both satisfy the $2n$ conditions with degree at most $2n - 1$. Let $R = H - G$. Then $R(x_i) = 0$ and $R'(x_i) = 0$ at every node. So each $x_i$ is a double root, and

$$
R(x) = q(x)\,(x - x_1)^2 \cdots (x - x_n)^2.
$$

If $q \ne 0$, then $R$ has degree at least $2n$. But its degree is at most $2n - 1$. So $q = 0$ and $R = 0$.

This is the same argument as for ordinary interpolation. The only change is that each node now counts twice.

## 7. Experiment: when does a higher degree help?

I interpolated two functions on $[-1, 1]$ with degrees $n = 2, 4, \dots, 32$:

- $f(x) = \dfrac{1}{1 + 16x^2}$, a smooth bump, and
- $f(x) = e^x$.

I used two kinds of nodes:

- **equal**: $x_i = -1 + 2i/n$, and
- **clustered**: $x_i = -\cos(i\pi/n)$. These points are dense near $\pm 1$ and sparse in the middle.

I built every interpolant with the divided difference table and evaluated it in Newton's form. Then I measured the largest error on 4001 equally spaced points in $[-1, 1]$. For $e^x$ I also computed both sides of the error formula. All derivatives of $e^x$ lie between $e^{-1}$ and $e$ on $[-1, 1]$, so the largest error must lie between

$$
\frac{\max|W(x)|}{e\,(n+1)!} \quad\text{and}\quad \frac{e\,\max|W(x)|}{(n+1)!}.
$$

<figure class="viz">
<div class="viz-legend"><span class="viz-key"><span class="viz-swatch" style="background:#2a78d6"></span>1/(1+16x²), equal nodes</span><span class="viz-key"><span class="viz-swatch" style="background:#eb6834"></span>1/(1+16x²), clustered nodes</span><span class="viz-key"><span class="viz-swatch" style="background:#1baf7a"></span>eˣ, equal nodes</span><span class="viz-key"><span class="viz-swatch" style="background:#eda100"></span>eˣ, error formula (upper)</span></div>
<div class="viz-frame">
<svg viewBox="0 0 680 400" role="img" tabindex="0" aria-label="Largest interpolation error against degree, log scale. With equal nodes the error for 1/(1+16x²) grows past 700. With clustered nodes it falls below 10^-3. For e^x the error follows the error formula down to round-off near degree 14, then grows slowly." class="viz-plot"><line x1="64" x2="540" y1="30.0" y2="30.0" class="viz-grid"/><text x="54" y="34.0" text-anchor="end" class="viz-tick">10<tspan dy="-0.5em" font-size="75%">4</tspan></text><line x1="64" x2="540" y1="92.0" y2="92.0" class="viz-grid"/><text x="54" y="96.0" text-anchor="end" class="viz-tick">1</text><line x1="64" x2="540" y1="154.0" y2="154.0" class="viz-grid"/><text x="54" y="158.0" text-anchor="end" class="viz-tick">10<tspan dy="-0.5em" font-size="75%">−4</tspan></text><line x1="64" x2="540" y1="216.0" y2="216.0" class="viz-grid"/><text x="54" y="220.0" text-anchor="end" class="viz-tick">10<tspan dy="-0.5em" font-size="75%">−8</tspan></text><line x1="64" x2="540" y1="278.0" y2="278.0" class="viz-grid"/><text x="54" y="282.0" text-anchor="end" class="viz-tick">10<tspan dy="-0.5em" font-size="75%">−12</tspan></text><line x1="64" x2="540" y1="340.0" y2="340.0" class="viz-grid"/><text x="54" y="344.0" text-anchor="end" class="viz-tick">10<tspan dy="-0.5em" font-size="75%">−16</tspan></text><text x="95.7" y="362" text-anchor="middle" class="viz-tick">4</text><text x="159.2" y="362" text-anchor="middle" class="viz-tick">8</text><text x="222.7" y="362" text-anchor="middle" class="viz-tick">12</text><text x="286.1" y="362" text-anchor="middle" class="viz-tick">16</text><text x="349.6" y="362" text-anchor="middle" class="viz-tick">20</text><text x="413.1" y="362" text-anchor="middle" class="viz-tick">24</text><text x="476.5" y="362" text-anchor="middle" class="viz-tick">28</text><text x="540.0" y="362" text-anchor="middle" class="viz-tick">32</text><text x="302.0" y="388" text-anchor="middle" class="viz-axis-title">Degree n</text><text x="64" y="14" class="viz-axis-title">Largest error on [−1, 1]</text><path d="M64.0,95.7 L95.7,98.4 L127.5,96.8 L159.2,94.1 L190.9,90.9 L222.7,87.4 L254.4,83.8 L286.1,80.0 L317.9,76.2 L349.6,72.3 L381.3,68.3 L413.1,64.3 L444.8,60.2 L476.5,56.1 L508.3,52.0 L540.0,47.8" fill="none" stroke="#2a78d6" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><circle cx="540.0" cy="47.8" r="4" fill="#2a78d6" stroke="#ffffff" stroke-width="2"/><path d="M64.0,95.7 L95.7,98.7 L127.5,102.1 L159.2,105.7 L190.9,109.5 L222.7,113.3 L254.4,116.1 L286.1,119.2 L317.9,122.4 L349.6,125.7 L381.3,129.1 L413.1,132.5 L444.8,135.9 L476.5,139.1 L508.3,142.4 L540.0,145.7" fill="none" stroke="#eb6834" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><circle cx="540.0" cy="145.7" r="4" fill="#eb6834" stroke="#ffffff" stroke-width="2"/><path d="M64.0,109.1 L95.7,137.7 L127.5,169.5 L159.2,204.2 L190.9,241.2 L222.7,280.4 L254.4,314.5 L286.1,306.1 L317.9,300.1 L349.6,303.4 L381.3,278.7 L413.1,267.3 L444.8,265.0 L476.5,261.9 L508.3,245.3 L540.0,238.2" fill="none" stroke="#1baf7a" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><circle cx="540.0" cy="238.2" r="4" fill="#1baf7a" stroke="#ffffff" stroke-width="2"/><path d="M64.0,103.8 L95.7,132.1 L127.5,163.7 L159.2,198.2 L190.9,235.2 L222.7,274.2 L254.4,315.2 L286.1,340.0" fill="none" stroke="#eda100" stroke-width="2" stroke-dasharray="6 4" stroke-linejoin="round" stroke-linecap="round"/><circle cx="286.1" cy="340.0" r="4" fill="#eda100" stroke="#ffffff" stroke-width="2"/><text x="550.0" y="51.8" class="viz-label">equal</text><text x="550.0" y="149.7" class="viz-label">clustered</text><text x="550.0" y="242.2" class="viz-label">eˣ, equal</text><text x="296.1" y="332.0" class="viz-label">formula</text><line x1="0" x2="0" y1="30" y2="340" class="viz-crosshair" visibility="hidden"/></svg>
<div class="viz-tooltip" hidden></div>
</div>
<figcaption>Largest error max |f(x) − p<sub>n</sub>(x)| on [−1, 1], on a log scale. The dashed line is the upper side of the error formula for e<sup>x</sup>; it leaves the chart below 10<sup>−16</sup>. Hover over the chart, or focus it and use the arrow keys, to read the values.</figcaption>
<details class="viz-table"><summary>Show the numbers</summary><table><thead><tr><th>n</th><th>1/(1+16x²), equal nodes</th><th>1/(1+16x²), clustered nodes</th><th>eˣ, equal nodes</th><th>eˣ, error formula (upper)</th></tr></thead><tbody><tr><td>2</td><td>5.7e-1</td><td>5.7e-1</td><td>7.9e-2</td><td>1.7e-1</td></tr><tr><td>4</td><td>3.9e-1</td><td>3.7e-1</td><td>1.1e-3</td><td>2.6e-3</td></tr><tr><td>6</td><td>4.9e-1</td><td>2.2e-1</td><td>1.0e-5</td><td>2.4e-5</td></tr><tr><td>8</td><td>7.3e-1</td><td>1.3e-1</td><td>5.8e-8</td><td>1.4e-7</td></tr><tr><td>10</td><td>1.2e0</td><td>7.5e-2</td><td>2.4e-10</td><td>5.8e-10</td></tr><tr><td>12</td><td>2.0e0</td><td>4.3e-2</td><td>7.0e-13</td><td>1.7e-12</td></tr><tr><td>14</td><td>3.4e0</td><td>2.8e-2</td><td>4.4e-15</td><td>4.0e-15</td></tr><tr><td>16</td><td>5.9e0</td><td>1.8e-2</td><td>1.5e-14</td><td>7.2e-18</td></tr><tr><td>18</td><td>1.0e1</td><td>1.1e-2</td><td>3.7e-14</td><td>1.0e-20</td></tr><tr><td>20</td><td>1.9e1</td><td>6.7e-3</td><td>2.3e-14</td><td>1.2e-23</td></tr><tr><td>22</td><td>3.4e1</td><td>4.0e-3</td><td>9.0e-13</td><td>1.2e-26</td></tr><tr><td>24</td><td>6.2e1</td><td>2.4e-3</td><td>4.9e-12</td><td>1.0e-29</td></tr><tr><td>26</td><td>1.1e2</td><td>1.5e-3</td><td>6.9e-12</td><td>7.6e-33</td></tr><tr><td>28</td><td>2.1e2</td><td>9.1e-4</td><td>1.1e-11</td><td>4.8e-36</td></tr><tr><td>30</td><td>3.8e2</td><td>5.6e-4</td><td>1.3e-10</td><td>2.7e-39</td></tr><tr><td>32</td><td>7.1e2</td><td>3.4e-4</td><td>3.7e-10</td><td>1.3e-42</td></tr></tbody></table></details>
<script>(() => {const fig = document.currentScript.closest('figure');const svg = fig.querySelector('svg');const line = svg.querySelector('.viz-crosshair');const tip = fig.querySelector('.viz-tooltip');const names = ["1/(1+16x²), equal nodes", "1/(1+16x²), clustered nodes", "eˣ, equal nodes", "eˣ, error formula (upper)"]; const colors = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"];const rows = [["5.7e-1", "5.7e-1", "7.9e-2", "1.7e-1"], ["3.9e-1", "3.7e-1", "1.1e-3", "2.6e-3"], ["4.9e-1", "2.2e-1", "1.0e-5", "2.4e-5"], ["7.3e-1", "1.3e-1", "5.8e-8", "1.4e-7"], ["1.2e0", "7.5e-2", "2.4e-10", "5.8e-10"], ["2.0e0", "4.3e-2", "7.0e-13", "1.7e-12"], ["3.4e0", "2.8e-2", "4.4e-15", "4.0e-15"], ["5.9e0", "1.8e-2", "1.5e-14", "7.2e-18"], ["1.0e1", "1.1e-2", "3.7e-14", "1.0e-20"], ["1.9e1", "6.7e-3", "2.3e-14", "1.2e-23"], ["3.4e1", "4.0e-3", "9.0e-13", "1.2e-26"], ["6.2e1", "2.4e-3", "4.9e-12", "1.0e-29"], ["1.1e2", "1.5e-3", "6.9e-12", "7.6e-33"], ["2.1e2", "9.1e-4", "1.1e-11", "4.8e-36"], ["3.8e2", "5.6e-4", "1.3e-10", "2.7e-39"], ["7.1e2", "3.4e-4", "3.7e-10", "1.3e-42"]]; const degrees = [2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32];const left = 64, right = 540, width = 680, last = degrees.length - 1;let current = -1;const show = (k) => { k = Math.max(0, Math.min(last, k)); current = k; const x = left + (right - left) * k / last; line.setAttribute('x1', x); line.setAttribute('x2', x); line.setAttribute('visibility', 'visible'); tip.replaceChildren(); const head = document.createElement('div'); head.className = 'viz-tooltip-head'; head.textContent = 'Degree ' + degrees[k]; tip.append(head); names.forEach((name, i) => {  const row = document.createElement('div'); row.className = 'viz-tooltip-row';  const key = document.createElement('span'); key.className = 'viz-swatch'; key.style.background = colors[i];  const value = document.createElement('strong'); value.textContent = rows[k][i];  const label = document.createElement('span'); label.textContent = name;  row.append(key, value, label); tip.append(row); }); const box = svg.getBoundingClientRect(), frame = fig.getBoundingClientRect(); const xs = box.left - frame.left + x * box.width / width; tip.hidden = false; const room = frame.width - tip.offsetWidth - 8; tip.style.left = Math.max(0, Math.min(room, xs + 12)) + 'px'; tip.style.top = (box.top - frame.top + 8) + 'px';};const hide = () => { line.setAttribute('visibility', 'hidden'); tip.hidden = true; current = -1; };svg.addEventListener('pointermove', (event) => { const box = svg.getBoundingClientRect(); const x = (event.clientX - box.left) * width / box.width; show(Math.round((x - left) / (right - left) * last));});svg.addEventListener('pointerleave', hide);svg.addEventListener('focus', () => show(current < 0 ? 0 : current));svg.addEventListener('blur', hide);svg.addEventListener('keydown', (event) => { if (event.key === 'ArrowRight') { show(current + 1); event.preventDefault(); } if (event.key === 'ArrowLeft') { show(current - 1); event.preventDefault(); } if (event.key === 'Escape') hide();});})();</script>
</figure>

| degree $n$ | bump, equal | bump, clustered | $e^x$, equal | error formula for $e^x$ (lower, upper) |
|---|---|---|---|---|
| 4 | $3.85 \times 10^{-1}$ | $3.71 \times 10^{-1}$ | $1.12 \times 10^{-3}$ | $3.48 \times 10^{-4}$, $2.57 \times 10^{-3}$ |
| 8 | $7.32 \times 10^{-1}$ | $1.31 \times 10^{-1}$ | $5.80 \times 10^{-8}$ | $1.91 \times 10^{-8}$, $1.41 \times 10^{-7}$ |
| 12 | $1.97$ | $4.25 \times 10^{-2}$ | $6.99 \times 10^{-13}$ | $2.37 \times 10^{-13}$, $1.75 \times 10^{-12}$ |
| 20 | $1.88 \times 10^{1}$ | $6.67 \times 10^{-3}$ | $2.29 \times 10^{-14}$ | $1.68 \times 10^{-24}$, $1.24 \times 10^{-23}$ |
| 32 | $7.07 \times 10^{2}$ | $3.42 \times 10^{-4}$ | $3.68 \times 10^{-10}$ | $1.76 \times 10^{-43}$, $1.30 \times 10^{-42}$ |

What I see:

- **The bump with equal nodes gets worse.** The error is smallest at $n = 4$ ($0.385$). After that it grows. At $n = 32$ it is about $707$, but $f$ itself is never larger than $1$. Near the end, each step of $+2$ in degree multiplies the error by about $1.85$. The worst point also moves toward the ends: $|x| = 0.80$ at $n = 4$ and $|x| = 0.986$ at $n = 32$.
- **The same bump with clustered nodes gets better.** The error falls at every step, to $3.4 \times 10^{-4}$ at $n = 32$. Near the end, each step of $+2$ multiplies it by about $0.61$. So the function is not the problem. The equal spacing is.
- **$e^x$ follows the error formula.** Up to $n = 12$, the measured error is always between the two bounds. It is about $0.40$ to $0.45$ times the upper bound.
- **After $n = 14$, round-off takes over.** The formula says the error for $e^x$ should keep falling, to about $10^{-42}$ at $n = 32$. The measured error stops near $4 \times 10^{-15}$ at $n = 14$, and then it *grows* to $3.7 \times 10^{-10}$. The reason is round-off. Each data value $e^{x_i}$ is stored with a tiny relative error, about $10^{-16}$. The interpolant can magnify such data errors by up to $\max_x \sum_i |\ell_i(x)|$. I measured this factor too. For equal nodes it grows fast: $30$ at $n = 10$ and $2.4 \times 10^{7}$ at $n = 32$. For clustered nodes it stays small: $2.4$ at $n = 10$ and $3.2$ at $n = 32$. With equal nodes at $n = 32$, the data errors can grow to about $2.4 \times 10^7 \times e \times 1.1 \times 10^{-16} \approx 7 \times 10^{-9}$. The measured $3.7 \times 10^{-10}$ is below this limit. This is the ill-conditioning from Section 3, seen in numbers.

As a check, I also evaluated every interpolant with the Lagrange form. For the bump with equal nodes, the two forms differ by at most $3.4 \times 10^{-7}$, while the error there is about $700$. For $e^x$ at $n = 32$ they differ by $6.1 \times 10^{-9}$. That is larger than the Newton-form error, but it is still within the round-off limit above. At that degree both forms only show round-off, and their round-off is different.

Here is the core of the code, in Python:

```python
def divided_differences(xs, ys):
    """Returns f[x0], f[x0,x1], ..., f[x0,...,xn] (the top edge of the table)."""
    col = list(ys)
    top = [col[0]]
    for k in range(1, len(xs)):
        # f[x_i..x_{i+k}] = (f[x_{i+1}..x_{i+k}] - f[x_i..x_{i+k-1}]) / (x_{i+k} - x_i)
        col = [(col[i + 1] - col[i]) / (xs[i + k] - xs[i]) for i in range(len(col) - 1)]
        top.append(col[0])
    return top

def newton_eval(xs, coef, x):
    """p(x) = c0 + (x-x0)(c1 + (x-x1)(c2 + ...)), evaluated from the inside out."""
    p = coef[-1]
    for k in range(len(coef) - 2, -1, -1):
        p = coef[k] + (x - xs[k]) * p
    return p
```

Each pass of the loop in `divided_differences` computes one column of the table. `newton_eval` nests the products, so it needs only $n$ multiplications.

The full script, with all the numbers and the chart, is in [this site's repository](https://github.com/doubleyands/doubleyands.github.io/tree/main/scripts/posts/numerical-analysis/week-2).

## Summary

| topic | main idea |
|---|---|
| Fixed-point iteration for systems | Same two conditions as in 1D: stay in $D$, and $\max \lVert G \rVert_\infty < 1$. Hard to check in practice. |
| Newton's method for systems | Linearize and solve $F\boldsymbol{\delta} = -\mathbf{f}$. Quadratic, needs a good start. Never form $F^{-1}$. |
| Optimization | $\nabla f = \mathbf{0}$ is a nonlinear system. Newton's method with the Hessian. |
| Interpolating polynomial | Exactly one of degree $\le n$ through $n+1$ points. Every form gives the same one. |
| Error formula | $\dfrac{W(x)}{(n+1)!} f^{(n+1)}(\xi)$. Valid inside the nodes' interval. |
| Divided differences | Leading coefficient. A recursion gives a cheap table. With $k+1$ points it equals $f^{(k)}(\xi)/k!$. |
| Runge phenomenon | Equal spacing can make the error grow with $n$, even for a smooth $f$. |
| Hermite | Values and slopes: $2n$ conditions, degree $\le 2n-1$, unique. |

The main thing I take from this week: a high-degree polynomial is not automatically a better fit. The error formula has two parts, $f^{(n+1)}$ and $W(x)$. Both can grow with $n$. Where the nodes go controls $W(x)$, and it also controls how much round-off gets magnified. Clustered nodes did much better than equal nodes in both ways.
