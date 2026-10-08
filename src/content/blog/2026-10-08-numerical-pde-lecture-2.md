---
title: "Numerical PDEs, Lecture 2: Building Difference Formulas and the Two-Point Boundary Value Problem"
description: "My notes on the second lecture: what O(h^2) really promises, how to build any difference formula, the heat equation and its boundary conditions, and a first finite difference solver, tested on a smooth and a rough solution."
pubDate: 2026-10-08
project: numerical-pde
tags: [numerical PDEs, finite differences, boundary value problems, notes]
---

This is my summary of the second lecture in my numerical PDE course. I write the ideas in my own words and check them with a small experiment.

[Lecture 1](/blog/2026-10-06-numerical-pde-lecture-1/) ended with the forward, backward, and central differences and their Taylor error. This lecture goes further in two ways. First, we learn to *build* a difference formula for any set of points. Second, we use finite differences to solve a real problem for the first time: the steady heat equation in one dimension.

## 1. What $O(h^2)$ really promises

Lecture 1 showed that the central difference has error $O(h^2)$:

$$
\frac{u(x + h) - u(x - h)}{2h} = u'(x) + \frac{h^2}{6} u'''(\xi)
$$

for some point $\xi$ between $x - h$ and $x + h$. The statement "$O(h^2)$" means: there is a constant $C$, **independent of $h$**, with error $\le C h^2$. Here $C$ is $\frac{1}{6} \max |u'''|$ near $x$.

This hides two assumptions.

1. **$u$ must be differentiable.** At a kink, like $|x|$ at $0$, the formula has nothing to approximate. Away from the kink it is fine.
2. **$u'''$ must be bounded.** If it is not, there is no finite $C$, and we cannot say "$O(h^2)$". If only $u''$ is bounded, the honest claim may be first order.

So the order of a formula is a promise that holds only when the function is smooth enough. The experiment in section 9 tests what happens when this assumption fails.

Also note the **reference point**. A formula approximates the derivative at one chosen point $x$. Whether we use points on the right, on the left, or on both sides is a second choice.

## 2. Building a formula: undetermined coefficients

So far we checked formulas that someone gave us. Now we turn it around: pick the points first, then find the weights.

Say we want $u'(x)$ and we can only use values on the left: $u(x)$, $u(x - h)$, $u(x - 2h)$. This happens at the right end of a domain, where there is no point to the right. Write

$$
Du(x) = a\, u(x) + b\, u(x - h) + c\, u(x - 2h),
$$

where $D$ means "a difference", not a derivative. Expand each value around $x$:

$$
u(x - h) = u - h u' + \frac{h^2}{2} u'' - \frac{h^3}{6} u''' + \cdots, \qquad
u(x - 2h) = u - 2h u' + 2h^2 u'' - \frac{4h^3}{3} u''' + \cdots
$$

(all derivatives at $x$). Collect the terms:

$$
Du(x) = (a + b + c)\, u - (b + 2c)\, h\, u' + \frac{b + 4c}{2}\, h^2 u'' - \frac{b + 8c}{6}\, h^3 u''' + \cdots
$$

We want only $u'$ to survive. So the $u$ term must vanish, the $u'$ term must have coefficient $1$, and the $u''$ term must vanish:

$$
a + b + c = 0, \qquad -(b + 2c)\, h = 1, \qquad b + 4c = 0.
$$

The solution is $a = \frac{3}{2h}$, $b = -\frac{2}{h}$, $c = \frac{1}{2h}$:

$$
Du(x) = \frac{3u(x) - 4u(x - h) + u(x - 2h)}{2h}.
$$

**Why exactly three conditions?** We have three unknowns. With two conditions there are too many solutions. With four there may be none. We kill the lowest terms first, in order. Picking a different set of three conditions would give a worse formula.

**What is the order?** The first term left over is

$$
-\frac{b + 8c}{6}\, h^3 u''' = -\frac{1}{6} \cdot \frac{2}{h} \cdot h^3 u''' = -\frac{h^2}{3}\, u'''.
$$

It is easy to see $h^3$ and say "third order". That is wrong. The weights $b$ and $c$ are not constants: they contain $1/h$. So one power of $h$ is lost, and the formula is **second order**. This is the same thing that happened in Lecture 1, when we divided by $h$ to get the forward difference.

The method works for any set of points and any derivative. Add $u(x - 3h)$ and you get a third-order formula. Mix points on both sides if you like. The algebra gets long by hand, so it is fine to let a computer solve the small linear system. The experiment does exactly that.

## 3. A quick check: the weights add up to zero

Here is a simple test for any difference formula for a derivative. Put in a constant function $u \equiv K$. Its derivative is $0$, so the formula must give $0$. But every value is $K$, so the formula gives $K$ times the sum of the weights. So **the weights must add up to zero.**

For the formula above: $3 - 4 + 1 = 0$. For the central difference: $1 - 1 = 0$. For $u''$: $1 - 2 + 1 = 0$.

A stronger test is to put in $u(x) = x$. A formula for $u'$ must then give exactly $1$. For the formula above, at $x = 0$: $\frac{1}{2h}\big(0 - 4(-h) + (-2h)\big) = 1$. Good.

These tests cannot prove a formula is right. But they catch a wrong sign or a wrong factor in a few seconds.

## 4. More points, higher order

| points used | usual order |
|---|---|
| 2 | 1 |
| 3 | 2 |
| 4 | 3 |
| 5 | 4 |

Each extra point buys one more order, if the point carries useful information. Symmetric formulas often get one order for free. The central difference uses two points but is second order, because the odd terms cancel. In the same way, two points on each side give fourth order, not third (we saw this formula in Lecture 1).

## 5. Higher derivatives by combining differences

A **difference operator** takes a function and returns a function. For example, $D_+$ maps $u$ to the function $x \mapsto (u(x + h) - u(x))/h$. So we can apply one operator after another.

The standard second difference is a difference of differences:

$$
D^2 u(x) = D_+ D_- u(x) = D_- D_+ u(x) = \frac{u(x - h) - 2u(x) + u(x + h)}{h^2}.
$$

The order of $D_+$ and $D_-$ does not matter. Note the $h^2$: each difference divides by $h$ once. This formula uses three points, one on each side, and it is second order. It is the most used formula in this course.

For a third derivative, apply one more difference:

$$
D_+ D^2 u(x) = \frac{u(x + 2h) - 3u(x + h) + 3u(x) - u(x - h)}{h^3} \qquad \text{(first order)},
$$

$$
D_0 D_+ D_- u(x) = \frac{u(x + 2h) - 2u(x + h) + 2u(x - h) - u(x - 2h)}{2h^3} \qquad \text{(second order)},
$$

where $D_0$ is the central difference. Both use four points, but the second one is symmetric around $x$, so it gains an order. Both pass the zero-sum test: $1 - 3 + 3 - 1 = 0$ and $1 - 2 + 2 - 1 = 0$.

## 6. The heat equation and its conditions

Now a real PDE. The heat (or diffusion) equation is

$$
u_t = \nabla \cdot \big(\kappa(\mathbf{x}) \nabla u\big) + f.
$$

Here $u$ is the temperature. The coefficient $\kappa$ describes how easily heat moves through the material. The term $f$ is a heat source (if $f < 0$, a sink). If the material is the same everywhere, $\kappa$ is a constant and comes out of the derivative:

$$
u_t = \kappa \Delta u + f.
$$

**Initial condition.** There is a time derivative, so integrating in time leaves one free function. We fix it with the starting temperature: $u(0, \mathbf{x}) = u^0(\mathbf{x})$.

**Boundary conditions** on $\partial\Omega$:

- **Dirichlet**: the value is given, $u = g_D$.
- **Neumann**: the normal derivative is given, $\partial u / \partial n = \nabla u \cdot \mathbf{n} = g_N$. This is the heat flux through the wall. If $g_N = 0$, no heat passes: the wall is **insulated**.

Think of a metal rod. Hold one end in ice water: that end has a fixed temperature, a Dirichlet condition. Wrap the other end in thick insulation: no heat flows out, a Neumann condition with $g_N = 0$.

On one piece of the boundary we give *one* condition, not both. Giving both the value and the flux on the same piece is too much information, and in general no solution fits both. But different pieces may have different types. This is called a **mixed** boundary condition.

**Steady state.** If nothing depends on time, then $u_t = 0$ and

$$
-\kappa \Delta u = f.
$$

Time is gone, and so is the initial condition.

## 7. The two-point boundary value problem

Three dimensions are hard to analyze and expensive to compute, so we learn the method in 1D first. With $\kappa = 1$:

$$
u''(x) = f(x), \quad 0 < x < 1, \qquad u(0) = \alpha, \quad u(1) = \beta.
$$

This is now an ODE. It is second order, so it needs two conditions, one at each end. It is called a **two-point boundary value problem**.

(A note on signs. The physical form is $-u'' = f$. Many books drop the minus sign because it looks cleaner. That is fine if we are consistent, but it changes the sign of the matrix below.)

**The grid.** Take $n$ intervals, $h = 1/n$, and points $x_j = jh$. The values at $x_0$ and $x_n$ are known from the boundary conditions, so there are $n - 1$ unknowns $U_1, \ldots, U_{n-1}$. For $n = 10$ that is $11$ points and $9$ unknowns. A finite difference method only cares about values at grid points. It says nothing about what happens between them.

**The scheme.** The exact solution satisfies

$$
\underbrace{u''(x_j)}_{= f(x_j)} = \frac{u(x_{j-1}) - 2u(x_j) + u(x_{j+1})}{h^2} + O(h^2).
$$

When we drop the $O(h^2)$ term, the equation is no longer true for $u$. So we give the unknowns a new name, $U_j$, and *define* them by

$$
\frac{U_{j-1} - 2U_j + U_{j+1}}{h^2} = f(x_j), \qquad j = 1, \ldots, n - 1, \qquad U_0 = \alpha,\ U_n = \beta.
$$

Writing $u$ and $U$ with different letters is not pedantry. $u(x_j)$ is the exact value, and $U_j$ is what the computer gives us. Their difference is the error we want to measure.

**The matrix.** In the first equation ($j = 1$), $U_0 = \alpha$ is known. We keep the unknowns on the left and move the known values to the right. The same happens with $U_n = \beta$ in the last equation. The result is $AU = F$ with

$$
A = \frac{1}{h^2}
\begin{bmatrix}
-2 & 1 & & & \\
1 & -2 & 1 & & \\
& \ddots & \ddots & \ddots & \\
& & 1 & -2 & 1 \\
& & & 1 & -2
\end{bmatrix},
\qquad
F =
\begin{bmatrix}
f(x_1) - \alpha / h^2 \\
f(x_2) \\
\vdots \\
f(x_{n-2}) \\
f(x_{n-1}) - \beta / h^2
\end{bmatrix}.
$$

Each equation uses three neighboring values, so this is a **three-point** scheme.

## 8. Properties of the matrix

$A$ is

- **symmetric**: $A^{\mathsf{T}} = A$;
- **tridiagonal**: only the main diagonal and the two next to it are nonzero;
- **sparse**: almost all entries are zero (the opposite is *dense*);
- **negative definite**, with the sign convention $u'' = f$.

The last point needs care. A matrix is positive definite if $\mathbf{v}^{\mathsf{T}} A \mathbf{v} > 0$ for every $\mathbf{v} \ne \mathbf{0}$. Try $\mathbf{v} = (1, 0, \ldots, 0)$: we get $A_{11} = -2/h^2 < 0$. So $A$ is not positive definite. In fact, a short calculation gives

$$
\mathbf{v}^{\mathsf{T}} A \mathbf{v} = -\frac{1}{h^2} \Big( v_1^2 + \sum_{j=1}^{n-2} (v_{j+1} - v_j)^2 + v_{n-1}^2 \Big),
$$

which is negative for every $\mathbf{v} \ne \mathbf{0}$. So $A$ is **negative definite**, and in particular it is invertible: the discrete problem always has exactly one solution. If we write the problem as $-u'' = f$, the matrix becomes $\frac{1}{h^2}\operatorname{tridiag}(-1, 2, -1)$, which is positive definite. The sign convention decides the name.

**Why sparsity matters.** The zeros are known, so we do not store them. We only store three diagonals: about $3n$ numbers instead of $n^2$. For PDEs this is essential, and most modern solvers are built around sparse matrices.

**Local truncation error.** Put the exact solution into the scheme and see how far it is from satisfied:

$$
\tau_j = \frac{u(x_{j-1}) - 2u(x_j) + u(x_{j+1})}{h^2} - f(x_j).
$$

Adding the Taylor expansions of $u(x_j \pm h)$ cancels the odd terms, so the next surviving term has the *fourth* derivative:

$$
\tau_j = \frac{h^2}{12}\, u''''(x_j) + O(h^4).
$$

**Solving the system: the Thomas algorithm.** This is Gaussian elimination for a tridiagonal matrix. One sweep down removes the lower diagonal, and one sweep up finds the unknowns. It takes about $8n$ operations, instead of about $\frac{2}{3}n^3$ for a dense matrix. Here is my version:

```python
def thomas(lower, diag, upper, rhs):
    """Solves a tridiagonal system. lower[0] and upper[-1] are not used."""
    n = len(diag)
    c = [0.0] * n
    d = [0.0] * n
    c[0] = upper[0] / diag[0]
    d[0] = rhs[0] / diag[0]
    for i in range(1, n):
        m = diag[i] - lower[i] * c[i - 1]
        c[i] = upper[i] / m if i < n - 1 else 0.0
        d[i] = (rhs[i] - lower[i] * d[i - 1]) / m
    x = [0.0] * n
    x[-1] = d[-1]
    for i in range(n - 2, -1, -1):
        x[i] = d[i] - c[i] * x[i + 1]
    return x
```

There is no pivoting. For this matrix that is safe, because the elimination never divides by zero. In general, before trusting a solver, we should know whether the matrix is invertible. Otherwise, when the output looks strange, we cannot tell a bug in the code from a problem with no solution.

## 9. Experiment

The script has two parts. Both use only the Python standard library.

**Part 1: undetermined coefficients by computer.** For a set of offsets $s_j$ and a derivative order $m$, the script solves the small linear system from section 2 with exact fractions. Then it finds the first Taylor term that does not cancel, which gives the order.

| formula for | points ($x + s_j h$) | weights (times $1/h^m$) | sum | order | leading error |
|---|---|---|---|---|---|
| $u'$ | $0, -1, -2$ | $\frac{3}{2}, -2, \frac{1}{2}$ | $0$ | 2 | $-\frac{1}{3} h^2 u'''$ |
| $u'$ | $0, -1, -2, -3$ | $\frac{11}{6}, -3, \frac{3}{2}, -\frac{1}{3}$ | $0$ | 3 | $-\frac{1}{4} h^3 u''''$ |
| $u'$ | $-1, 1$ | $-\frac{1}{2}, \frac{1}{2}$ | $0$ | 2 | $\frac{1}{6} h^2 u'''$ |
| $u'$ | $-2, -1, 1, 2$ | $\frac{1}{12}, -\frac{2}{3}, \frac{2}{3}, -\frac{1}{12}$ | $0$ | 4 | $-\frac{1}{30} h^4 u^{(5)}$ |
| $u''$ | $-1, 0, 1$ | $1, -2, 1$ | $0$ | 2 | $\frac{1}{12} h^2 u''''$ |
| $u'''$ | $-1, 0, 1, 2$ | $-1, 3, -3, 1$ | $0$ | 1 | $\frac{1}{2} h\, u''''$ |
| $u'''$ | $-2, -1, 1, 2$ | $-\frac{1}{2}, 1, -1, \frac{1}{2}$ | $0$ | 2 | $\frac{1}{4} h^2 u^{(5)}$ |

The first row is the formula from section 2, with the same error term. The last two rows are $D_+ D^2$ and $D_0 D_+ D_-$ from section 5. Every row passes the zero-sum test, and the symmetric rows are one order higher than the point count suggests.

**Part 2: solving $u'' = f$.** I solved the boundary value problem with the three-point scheme and the Thomas algorithm, for $n = 4, 8, \ldots, 4096$. I chose two exact solutions and computed $f = u''$, $\alpha = u(0)$, $\beta = u(1)$ from them.

- **Smooth:** $u(x) = \sin 3x + x^2 + 1$. Here $u'''' = 81 \sin 3x$, so $\max |u''''| = 81$.
- **Rough:** $u(x) = x^{5/2}$. Here $u'''' = \frac{15}{16} x^{-3/2}$, which is unbounded near $x = 0$. This breaks the assumption from section 1.

For each $n$ I measured the **global error** $\max_j |U_j - u(x_j)|$ and the largest **local truncation error** $\max_j |\tau_j|$.

<figure class="viz">
<div class="viz-legend"><span class="viz-key"><span class="viz-swatch" style="background:#2a78d6"></span>Smooth, global error</span><span class="viz-key"><span class="viz-swatch" style="background:#eb6834"></span>Smooth, max |τ|</span><span class="viz-key"><span class="viz-swatch" style="background:#1baf7a"></span>Rough, global error</span><span class="viz-key"><span class="viz-swatch" style="background:#eda100"></span>Rough, max |τ|</span></div>
<div class="viz-frame">
<svg viewBox="0 0 680 400" role="img" tabindex="0" aria-label="Global error and largest local truncation error of the 3-point scheme as h gets smaller, on log scales. Three lines fall with slope 2. The local truncation error of the rough solution falls with slope one half." class="viz-plot"><line x1="64" x2="520" y1="30.0" y2="30.0" class="viz-grid"/><text x="54" y="34.0" text-anchor="end" class="viz-tick">1</text><line x1="64" x2="520" y1="133.3" y2="133.3" class="viz-grid"/><text x="54" y="137.3" text-anchor="end" class="viz-tick">10<tspan dy="-0.5em" font-size="75%">−3</tspan></text><line x1="64" x2="520" y1="236.7" y2="236.7" class="viz-grid"/><text x="54" y="240.7" text-anchor="end" class="viz-tick">10<tspan dy="-0.5em" font-size="75%">−6</tspan></text><line x1="64" x2="520" y1="340.0" y2="340.0" class="viz-grid"/><text x="54" y="344.0" text-anchor="end" class="viz-tick">10<tspan dy="-0.5em" font-size="75%">−9</tspan></text><text x="64.0" y="362" text-anchor="middle" class="viz-tick">2<tspan dy="-0.5em" font-size="75%">−2</tspan></text><text x="155.2" y="362" text-anchor="middle" class="viz-tick">2<tspan dy="-0.5em" font-size="75%">−4</tspan></text><text x="246.4" y="362" text-anchor="middle" class="viz-tick">2<tspan dy="-0.5em" font-size="75%">−6</tspan></text><text x="337.6" y="362" text-anchor="middle" class="viz-tick">2<tspan dy="-0.5em" font-size="75%">−8</tspan></text><text x="428.8" y="362" text-anchor="middle" class="viz-tick">2<tspan dy="-0.5em" font-size="75%">−10</tspan></text><text x="520.0" y="362" text-anchor="middle" class="viz-tick">2<tspan dy="-0.5em" font-size="75%">−12</tspan></text><text x="292.0" y="388" text-anchor="middle" class="viz-axis-title">Step size h (smaller to the right)</text><text x="64" y="14" class="viz-axis-title">Error</text><path d="M64.0,76.5 L109.6,97.5 L155.2,118.4 L200.8,139.1 L246.4,159.9 L292.0,180.6 L337.6,201.3 L383.2,222.1 L428.8,242.8 L474.4,263.5 L520.0,284.3" fill="none" stroke="#2a78d6" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><circle cx="64.0" cy="76.5" r="3" fill="#2a78d6" stroke="#ffffff" stroke-width="1.5"/><circle cx="109.6" cy="97.5" r="3" fill="#2a78d6" stroke="#ffffff" stroke-width="1.5"/><circle cx="155.2" cy="118.4" r="3" fill="#2a78d6" stroke="#ffffff" stroke-width="1.5"/><circle cx="200.8" cy="139.1" r="3" fill="#2a78d6" stroke="#ffffff" stroke-width="1.5"/><circle cx="246.4" cy="159.9" r="3" fill="#2a78d6" stroke="#ffffff" stroke-width="1.5"/><circle cx="292.0" cy="180.6" r="3" fill="#2a78d6" stroke="#ffffff" stroke-width="1.5"/><circle cx="337.6" cy="201.3" r="3" fill="#2a78d6" stroke="#ffffff" stroke-width="1.5"/><circle cx="383.2" cy="222.1" r="3" fill="#2a78d6" stroke="#ffffff" stroke-width="1.5"/><circle cx="428.8" cy="242.8" r="3" fill="#2a78d6" stroke="#ffffff" stroke-width="1.5"/><circle cx="474.4" cy="263.5" r="3" fill="#2a78d6" stroke="#ffffff" stroke-width="1.5"/><circle cx="520.0" cy="284.3" r="3" fill="#2a78d6" stroke="#ffffff" stroke-width="1.5"/><path d="M64.0,43.2 L109.6,63.8 L155.2,84.4 L200.8,105.1 L246.4,125.9 L292.0,146.6 L337.6,167.3 L383.2,188.1 L428.8,208.8 L474.4,229.5 L520.0,249.7" fill="none" stroke="#eb6834" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><circle cx="64.0" cy="43.2" r="3" fill="#eb6834" stroke="#ffffff" stroke-width="1.5"/><circle cx="109.6" cy="63.8" r="3" fill="#eb6834" stroke="#ffffff" stroke-width="1.5"/><circle cx="155.2" cy="84.4" r="3" fill="#eb6834" stroke="#ffffff" stroke-width="1.5"/><circle cx="200.8" cy="105.1" r="3" fill="#eb6834" stroke="#ffffff" stroke-width="1.5"/><circle cx="246.4" cy="125.9" r="3" fill="#eb6834" stroke="#ffffff" stroke-width="1.5"/><circle cx="292.0" cy="146.6" r="3" fill="#eb6834" stroke="#ffffff" stroke-width="1.5"/><circle cx="337.6" cy="167.3" r="3" fill="#eb6834" stroke="#ffffff" stroke-width="1.5"/><circle cx="383.2" cy="188.1" r="3" fill="#eb6834" stroke="#ffffff" stroke-width="1.5"/><circle cx="428.8" cy="208.8" r="3" fill="#eb6834" stroke="#ffffff" stroke-width="1.5"/><circle cx="474.4" cy="229.5" r="3" fill="#eb6834" stroke="#ffffff" stroke-width="1.5"/><circle cx="520.0" cy="249.7" r="3" fill="#eb6834" stroke="#ffffff" stroke-width="1.5"/><path d="M64.0,118.2 L109.6,136.1 L155.2,154.9 L200.8,174.4 L246.4,194.4 L292.0,214.6 L337.6,234.9 L383.2,255.4 L428.8,275.9 L474.4,296.6 L520.0,317.2" fill="none" stroke="#1baf7a" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><circle cx="64.0" cy="118.2" r="3" fill="#1baf7a" stroke="#ffffff" stroke-width="1.5"/><circle cx="109.6" cy="136.1" r="3" fill="#1baf7a" stroke="#ffffff" stroke-width="1.5"/><circle cx="155.2" cy="154.9" r="3" fill="#1baf7a" stroke="#ffffff" stroke-width="1.5"/><circle cx="200.8" cy="174.4" r="3" fill="#1baf7a" stroke="#ffffff" stroke-width="1.5"/><circle cx="246.4" cy="194.4" r="3" fill="#1baf7a" stroke="#ffffff" stroke-width="1.5"/><circle cx="292.0" cy="214.6" r="3" fill="#1baf7a" stroke="#ffffff" stroke-width="1.5"/><circle cx="337.6" cy="234.9" r="3" fill="#1baf7a" stroke="#ffffff" stroke-width="1.5"/><circle cx="383.2" cy="255.4" r="3" fill="#1baf7a" stroke="#ffffff" stroke-width="1.5"/><circle cx="428.8" cy="275.9" r="3" fill="#1baf7a" stroke="#ffffff" stroke-width="1.5"/><circle cx="474.4" cy="296.6" r="3" fill="#1baf7a" stroke="#ffffff" stroke-width="1.5"/><circle cx="520.0" cy="317.2" r="3" fill="#1baf7a" stroke="#ffffff" stroke-width="1.5"/><path d="M64.0,75.9 L109.6,81.1 L155.2,86.2 L200.8,91.4 L246.4,96.6 L292.0,101.8 L337.6,107.0 L383.2,112.2 L428.8,117.4 L474.4,122.5 L520.0,127.7" fill="none" stroke="#eda100" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><circle cx="64.0" cy="75.9" r="3" fill="#eda100" stroke="#ffffff" stroke-width="1.5"/><circle cx="109.6" cy="81.1" r="3" fill="#eda100" stroke="#ffffff" stroke-width="1.5"/><circle cx="155.2" cy="86.2" r="3" fill="#eda100" stroke="#ffffff" stroke-width="1.5"/><circle cx="200.8" cy="91.4" r="3" fill="#eda100" stroke="#ffffff" stroke-width="1.5"/><circle cx="246.4" cy="96.6" r="3" fill="#eda100" stroke="#ffffff" stroke-width="1.5"/><circle cx="292.0" cy="101.8" r="3" fill="#eda100" stroke="#ffffff" stroke-width="1.5"/><circle cx="337.6" cy="107.0" r="3" fill="#eda100" stroke="#ffffff" stroke-width="1.5"/><circle cx="383.2" cy="112.2" r="3" fill="#eda100" stroke="#ffffff" stroke-width="1.5"/><circle cx="428.8" cy="117.4" r="3" fill="#eda100" stroke="#ffffff" stroke-width="1.5"/><circle cx="474.4" cy="122.5" r="3" fill="#eda100" stroke="#ffffff" stroke-width="1.5"/><circle cx="520.0" cy="127.7" r="3" fill="#eda100" stroke="#ffffff" stroke-width="1.5"/><text x="530" y="288.3" class="viz-label">Smooth, global error</text><text x="530" y="253.7" class="viz-label">Smooth, max |τ|</text><text x="530" y="321.2" class="viz-label">Rough, global error</text><text x="530" y="131.7" class="viz-label">Rough, max |τ|</text><line x1="0" x2="0" y1="30" y2="340" class="viz-crosshair" visibility="hidden"/></svg>
<div class="viz-tooltip" hidden></div>
</div>
<figcaption>The 3-point scheme for u″ = f on (0, 1), on log scales. Smooth: u = sin 3x + x<sup>2</sup> + 1. Rough: u = x<sup>5/2</sup>. τ is the local truncation error. Hover over the chart, or focus it and use the arrow keys, to read the values.</figcaption>
<details class="viz-table"><summary>Show the numbers</summary><table><thead><tr><th>h</th><th>Smooth, global error</th><th>Smooth, max |τ|</th><th>Rough, global error</th><th>Rough, max |τ|</th></tr></thead><tbody><tr><td>1/4</td><td>4.5e-2</td><td>4.1e-1</td><td>2.7e-3</td><td>4.7e-2</td></tr><tr><td>1/8</td><td>1.1e-2</td><td>1.0e-1</td><td>8.3e-4</td><td>3.3e-2</td></tr><tr><td>1/16</td><td>2.7e-3</td><td>2.6e-2</td><td>2.4e-4</td><td>2.3e-2</td></tr><tr><td>1/32</td><td>6.8e-4</td><td>6.6e-3</td><td>6.4e-5</td><td>1.6e-2</td></tr><tr><td>1/64</td><td>1.7e-4</td><td>1.6e-3</td><td>1.7e-5</td><td>1.2e-2</td></tr><tr><td>1/128</td><td>4.2e-5</td><td>4.1e-4</td><td>4.4e-6</td><td>8.2e-3</td></tr><tr><td>1/256</td><td>1.1e-5</td><td>1.0e-4</td><td>1.1e-6</td><td>5.8e-3</td></tr><tr><td>1/512</td><td>2.7e-6</td><td>2.6e-5</td><td>2.9e-7</td><td>4.1e-3</td></tr><tr><td>1/1024</td><td>6.6e-7</td><td>6.4e-6</td><td>7.2e-8</td><td>2.9e-3</td></tr><tr><td>1/2048</td><td>1.7e-7</td><td>1.6e-6</td><td>1.8e-8</td><td>2.1e-3</td></tr><tr><td>1/4096</td><td>4.1e-8</td><td>4.2e-7</td><td>4.6e-9</td><td>1.5e-3</td></tr></tbody></table></details>
<script>(() => {const fig = document.currentScript.closest('figure');const svg = fig.querySelector('svg');const line = svg.querySelector('.viz-crosshair');const tip = fig.querySelector('.viz-tooltip');const names = ["Smooth, global error", "Smooth, max |\u03c4|", "Rough, global error", "Rough, max |\u03c4|"]; const colors = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"];const rows = [["4.5e-2", "4.1e-1", "2.7e-3", "4.7e-2"], ["1.1e-2", "1.0e-1", "8.3e-4", "3.3e-2"], ["2.7e-3", "2.6e-2", "2.4e-4", "2.3e-2"], ["6.8e-4", "6.6e-3", "6.4e-5", "1.6e-2"], ["1.7e-4", "1.6e-3", "1.7e-5", "1.2e-2"], ["4.2e-5", "4.1e-4", "4.4e-6", "8.2e-3"], ["1.1e-5", "1.0e-4", "1.1e-6", "5.8e-3"], ["2.7e-6", "2.6e-5", "2.9e-7", "4.1e-3"], ["6.6e-7", "6.4e-6", "7.2e-8", "2.9e-3"], ["1.7e-7", "1.6e-6", "1.8e-8", "2.1e-3"], ["4.1e-8", "4.2e-7", "4.6e-9", "1.5e-3"]]; const heads = ["h = 1/4", "h = 1/8", "h = 1/16", "h = 1/32", "h = 1/64", "h = 1/128", "h = 1/256", "h = 1/512", "h = 1/1024", "h = 1/2048", "h = 1/4096"];const left = 64, right = 520, steps = 10, width = 680;let current = -1;const show = (n) => { n = Math.max(0, Math.min(steps, n)); current = n; const x = left + (right - left) * n / steps; line.setAttribute('x1', x); line.setAttribute('x2', x); line.setAttribute('visibility', 'visible'); tip.replaceChildren(); const head = document.createElement('div'); head.className = 'viz-tooltip-head'; head.textContent = heads[n]; tip.append(head); names.forEach((name, i) => {  const row = document.createElement('div'); row.className = 'viz-tooltip-row';  const key = document.createElement('span'); key.className = 'viz-swatch'; key.style.background = colors[i];  const value = document.createElement('strong'); value.textContent = rows[n][i];  const label = document.createElement('span'); label.textContent = name;  row.append(key, value, label); tip.append(row); }); const box = svg.getBoundingClientRect(), frame = fig.getBoundingClientRect(); const xs = box.left - frame.left + x * box.width / width; tip.hidden = false; const room = frame.width - tip.offsetWidth - 8; tip.style.left = Math.max(0, Math.min(room, xs + 12)) + 'px'; tip.style.top = (box.top - frame.top + 8) + 'px';};const hide = () => { line.setAttribute('visibility', 'hidden'); tip.hidden = true; current = -1; };svg.addEventListener('pointermove', (event) => { const box = svg.getBoundingClientRect(); const x = (event.clientX - box.left) * width / box.width; show(Math.round((x - left) / (right - left) * steps));});svg.addEventListener('pointerleave', hide);svg.addEventListener('focus', () => show(current < 0 ? 0 : current));svg.addEventListener('blur', hide);svg.addEventListener('keydown', (event) => { if (event.key === 'ArrowRight') { show(current + 1); event.preventDefault(); } if (event.key === 'ArrowLeft') { show(current - 1); event.preventDefault(); } if (event.key === 'Escape') hide();});})();</script>
</figure>

| $n$ | smooth: global error | smooth: $\max\lvert\tau\rvert$ | $\frac{h^2}{12} \cdot 81$ | rough: global error | rough: $\max\lvert\tau\rvert$ |
|---|---|---|---|---|---|
| 4 | $4.47 \times 10^{-2}$ | $4.13 \times 10^{-1}$ | $4.22 \times 10^{-1}$ | $2.75 \times 10^{-3}$ | $4.66 \times 10^{-2}$ |
| 16 | $2.72 \times 10^{-3}$ | $2.63 \times 10^{-2}$ | $2.64 \times 10^{-2}$ | $2.37 \times 10^{-4}$ | $2.33 \times 10^{-2}$ |
| 64 | $1.70 \times 10^{-4}$ | $1.65 \times 10^{-3}$ | $1.65 \times 10^{-3}$ | $1.69 \times 10^{-5}$ | $1.16 \times 10^{-2}$ |
| 256 | $1.06 \times 10^{-5}$ | $1.03 \times 10^{-4}$ | $1.03 \times 10^{-4}$ | $1.12 \times 10^{-6}$ | $5.82 \times 10^{-3}$ |
| 1024 | $6.63 \times 10^{-7}$ | $6.44 \times 10^{-6}$ | $6.44 \times 10^{-6}$ | $7.24 \times 10^{-8}$ | $2.91 \times 10^{-3}$ |
| 4096 | $4.15 \times 10^{-8}$ | $4.19 \times 10^{-7}$ | $4.02 \times 10^{-7}$ | $4.59 \times 10^{-9}$ | $1.46 \times 10^{-3}$ |

The measured order (from $n$ to $2n$):

- smooth, global error: $2.00$ from $n = 32$ on;
- smooth, $\max|\tau|$: $2.00$ from $n = 32$ to $2048$;
- rough, global error: $1.72, 1.81, 1.89, \ldots, 1.99$, slowly rising to $2$;
- rough, $\max|\tau|$: $0.50$ at every step.

What I noticed:

- **The smooth case matches the theory.** The largest $\tau$ agrees with $\frac{h^2}{12} \max|u''''|$ to within $1\%$ for $n \ge 8$. The largest $\tau$ sits near $x = 0.52$, where $3x = \pi/2$ and $|u''''|$ is largest. At $n = 4096$ the measured $\tau$ is $4\%$ too large. That is round-off: computing $\tau$ divides a difference of nearly equal numbers by $h^2$, the same effect we saw in Lecture 1.
- **In the rough case, $\tau$ is only order $1/2$.** The worst point is always $x_1 = h$, next to the boundary, where $u''''$ blows up. A short calculation gives $\tau_1 = (2^{5/2} - 2 - \frac{15}{4})\, h^{1/2} \approx -0.093\, h^{1/2}$. The measured values follow this exactly, for example $0.0466$ at $h = 1/4$. So when the smoothness assumption fails, the local promise of $O(h^2)$ really fails.
- **But the global error is still close to second order.** This surprised me. The large $\tau$ appears only at one or two points next to the boundary, where the solution is pinned by the boundary value. An error made there does not spread much into the rest of the solution. Explaining this properly needs the link between $\tau$ and the global error, which the next part of the course covers. So I only note it here: a large local error does not always mean a large global error.
- **The global error is about $10$ times smaller than $\max|\tau|$** in the smooth case. Again, the link between the two is the next topic.

The full script, with the exact-fraction solver and the chart, is in [this site's repository](https://github.com/doubleyands/doubleyands.github.io/tree/main/scripts/posts/numerical-pde/lecture-2).

## Summary

1. **"$O(h^p)$" is a promise with conditions.** It needs enough smooth derivatives, and a constant $C$ that does not depend on $h$.
2. **Undetermined coefficients** builds a difference formula for any set of points: expand in Taylor series and match terms. Watch out: the weights contain $1/h$, so the order is one less than the leftover power of $h$ suggests.
3. **The weights of a derivative formula add up to zero.** It is a two-second check.
4. More points give higher order, and symmetric formulas get one order for free. $D^2 = D_+ D_-$ is the workhorse.
5. The **heat equation** needs an initial condition and one boundary condition on each piece of the boundary: Dirichlet (value), Neumann (flux), or mixed.
6. The **three-point scheme** for $u'' = f$ gives a symmetric, tridiagonal, sparse, negative definite matrix. Known boundary values move to the right-hand side. The **Thomas algorithm** solves it in $O(n)$ operations.
7. The local truncation error is $\frac{h^2}{12} u''''$. With a smooth solution, both $\tau$ and the global error are second order. With a rough solution, $\tau$ can be much worse than the global error.
