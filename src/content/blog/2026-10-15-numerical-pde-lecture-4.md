---
title: "Numerical PDEs, Lecture 4: More on 2D Matrices, and Advection-Diffusion"
description: "My notes on the fourth lecture: sign conventions and the right-hand side of the 2D system, Neumann conditions in 2D, the 9-point stencil, red-black ordering, and why the direction of flow decides which boundary condition really matters."
pubDate: 2026-10-15
project: numerical-pde
tags: [numerical PDEs, finite differences, advection-diffusion, notes]
---

This is my summary of the fourth lecture in my numerical PDE course. I write the ideas in my own words and check them with a small experiment.

The lecture has two halves. The first half finishes the 2D Poisson problem from [Lecture 3](/blog/2026-10-13-numerical-pde-lecture-3/): signs, the right-hand side, Neumann conditions, a bigger stencil, and a second way to order the unknowns. The second half starts a new topic, the **advection-diffusion equation**. In that half the lecture only covered the physics, not the numerical scheme. So this post is about the physics too, and the scheme comes in the next post.

## 1. The 2D system, once more

### Two sign conventions

The 5-point stencil from Lecture 3 can be written in two ways. For the Laplace equation $\nabla^2 u = 0$:

$$
\frac{u_{i-1,j} + u_{i+1,j} + u_{i,j-1} + u_{i,j+1} - 4u_{i,j}}{h^2} = 0
\qquad\text{or}\qquad
\frac{4u_{i,j} - u_{i-1,j} - u_{i+1,j} - u_{i,j-1} - u_{i,j+1}}{h^2} = 0.
$$

The first form approximates $\nabla^2 u$. Its matrix has $-4$ on the diagonal and $+1$ for the neighbors, and it is **negative definite**. The second form approximates $-\nabla^2 u$. Its matrix has $+4$ and $-1$, and it is **positive definite**. Both are correct. The solution is the same, and only the sign of the matrix (and of the right-hand side) flips. When you mix two references, check which convention each one uses. This is the same issue as in [Lecture 2](/blog/2026-10-08-numerical-pde-lecture-2/), where the 1D matrix for $u'' = f$ was negative definite.

### The right-hand side

Use the setting of Lecture 3: the unit square split into $5$ intervals each way, $16$ interior unknowns in lexicographic order, and Dirichlet data $u = u_B$ on the boundary. In Lecture 3 we dropped every neighbor that lies on the boundary from the matrix. That value is known, so it moves to the right-hand side. With the first sign convention, the row of $u_{i,j}$ is

$$
\frac{(\text{interior neighbors}) - 4u_{i,j}}{h^2} = f_{i,j} - \frac{(\text{boundary neighbors})}{h^2}.
$$

So the right-hand side tells the same story as the matrix:

- A **corner** point like $u_{1,1}$ touches the boundary twice. Its right-hand side gets two values: $u_{0,1} + u_{1,0}$.
- An **edge** point like $u_{2,1}$ touches it once. It gets one value: $u_{2,0}$.
- An **inside** point like $u_{2,2}$ does not touch it. It gets nothing extra.

In Lecture 3 the corner row kept $3$ terms, the edge row $4$, and the inside row $5$. The terms that disappear from the left side are exactly the ones that appear on the right side.

### Penta-diagonal

We called this matrix **block tridiagonal**, because of its blocks $T$ and $I$. If we look at single entries instead of blocks, each row has at most five nonzero entries, on five diagonals: the main diagonal, the two next to it (east and west), and the two that are one grid row away (north and south). So the same matrix is also called **penta-diagonal**. Some of these diagonals have zeros in them. Rows next to the boundary have only three or four nonzero entries, as we saw above.

## 2. Neumann conditions in 2D

Now take the Poisson problem on the unit square with Dirichlet conditions on the bottom and top ($y = 0$ and $y = 1$) and Neumann conditions on the left and right ($x = 0$ and $x = 1$):

$$
\nabla^2 u = f, \qquad u = u_B \text{ on } y = 0, 1, \qquad \frac{\partial u}{\partial \mathbf{n}} = g \text{ on } x = 0, 1.
$$

### Computing $\partial u / \partial \mathbf{n}$

The normal derivative is $\partial u / \partial \mathbf{n} = \mathbf{n} \cdot \nabla u$, where $\mathbf{n}$ is the **unit outward normal** ([Lecture 1](/blog/2026-10-06-numerical-pde-lecture-1/)). On a square the normal is easy, but its sign depends on the side:

| Side | $\mathbf{n}$ | $\partial u / \partial \mathbf{n}$ |
|---|---|---|
| right, $x = 1$ | $(1, 0)$ | $+u_x$ |
| left, $x = 0$ | $(-1, 0)$ | $-u_x$ |
| top, $y = 1$ | $(0, 1)$ | $+u_y$ |
| bottom, $y = 0$ | $(0, -1)$ | $-u_y$ |

Why use one symbol for four different formulas? It is good bookkeeping. One symbol covers every side, and the outward direction makes its meaning clear: it measures how much flows **out** of the domain through the boundary (the flux).

### Discretizing: 1D at every boundary node

On the left side, $\partial u / \partial \mathbf{n} = -u_x$, and $u_x$ along one grid row is just a 1D derivative. So we repeat the 1D **ghost node** trick from [Lecture 3](/blog/2026-10-13-numerical-pde-lecture-3/) at every boundary node, one row at a time. Nothing new is needed.

**What is unknown?** On the Neumann sides we know the slope, but not the value. So the boundary values $u_{0,j}$ and $u_{n,j}$ are unknowns too. On the Dirichlet sides the values $u_{i,0}$ and $u_{i,n}$ are given. With $h = 1/5$, the unknowns are $u_{i,j}$ with $i = 0, \ldots, 5$ and $j = 1, \ldots, 4$: that is $6 \times 4 = 24$ unknowns, against $16$ in the all-Dirichlet case.

**The ghost row.** At the node $(0, j)$, write the 5-point stencil. It needs the value $u_{-1,j}$ at a ghost node outside the square. A central difference of the boundary condition gives

$$
-\frac{u_{1,j} - u_{-1,j}}{2h} = g(0, y_j) \quad\Longrightarrow\quad u_{-1,j} = u_{1,j} + 2h\, g(0, y_j).
$$

Put this into the stencil, and the ghost value is gone:

$$
\frac{2u_{1,j} + u_{0,j-1} + u_{0,j+1} - 4u_{0,j}}{h^2} = f(0, y_j) - \frac{2}{h}\, g(0, y_j).
$$

The right side $x = 1$ gives the same formula with $u_{n-1,j}$ in place of $u_{1,j}$. This works because $g$ is defined with the outward normal.

We could also use a one-sided difference for $u_x$. But that row looks different from all the other rows and spoils the clean pattern of the matrix, and a first-order row would also cost us the second order (the pollution effect of Lecture 3). The central difference with a ghost node keeps the pattern and keeps second order. The experiment checks this in 2D.

## 3. The 9-point stencil

The 5-point stencil uses only east, west, north, and south. We can also use the four diagonal neighbors. The standard **9-point stencil** is

$$
\nabla_9^2 u_{i,j} = \frac{1}{6h^2} \Big( 4\big(u_{i-1,j} + u_{i+1,j} + u_{i,j-1} + u_{i,j+1}\big) + \big(u_{i-1,j-1} + u_{i+1,j-1} + u_{i-1,j+1} + u_{i+1,j+1}\big) - 20\, u_{i,j} \Big).
$$

The weights add up to zero: $4 \cdot 4 + 4 \cdot 1 - 20 = 0$. This is the same quick check as in Lecture 2. A constant has zero derivatives, so any formula for a derivative must give zero for it.

Is it better than the 5-point stencil? In general, no. For the Poisson equation with $f \neq 0$, the 9-point stencil is still **second order**. It only gains for the Laplace equation, $f = 0$. Then the leading error terms vanish, and the lecture said it becomes fourth order. In my experiment it is even **sixth order** for the Laplace equation (section 7).

A wider stencil costs something. It has more nonzero entries per row, and it is harder to handle near the boundary. Higher order needs wider and more complex stencils. So the 5-point stencil is the most popular one: it is very simple, and very strong.

## 4. Ordering the unknowns: natural and red-black

In Lecture 3 we used the **natural row order** (lexicographic): one grid row after another. All the points of one row form one block of the matrix. With $m \times m$ interior points, the matrix is

$$
A = \frac{1}{h^2} \begin{bmatrix} T & I & & \\ I & T & \ddots & \\ & \ddots & \ddots & I \\ & & I & T \end{bmatrix}, \qquad T = \operatorname{tridiag}(1, -4, 1),
$$

with $m$ blocks of size $m \times m$. The $1$s in the blocks $I$ sit $m$ places away from the diagonal, because the north and south neighbors are one whole grid row away in the ordering. In 1D the matrix had all its entries within one place of the diagonal. In 2D this is not possible with any order: every point has two neighbors in the other direction, and in any order some of them must be far away.

**Red-black ordering** is a different idea. Color the grid like a checkerboard: point $(i, j)$ is red if $i + j$ is even and black if $i + j$ is odd. Number all red points first, then all black points. The 5-point stencil links a point only to its east, west, north, and south neighbors, and **all four of them have the other color**. So red points are never linked to red points, and black points are never linked to black points. The matrix becomes

$$
\begin{bmatrix} D_r & H \\ H^{\mathsf{T}} & D_b \end{bmatrix}
\begin{bmatrix} u_{\text{red}} \\ u_{\text{black}} \end{bmatrix}
=
\begin{bmatrix} f_{\text{red}} \\ f_{\text{black}} \end{bmatrix},
\qquad D_r = -\frac{4}{h^2} I, \quad D_b = -\frac{4}{h^2} I,
$$

where $H$ holds the red-black links. The diagonal blocks are **pure diagonal matrices**. (If $m^2$ is odd, the two colors have different numbers of points, so $D_r$ and $D_b$ have different sizes.)

Why is this useful? If we know all the black values, then each red value comes from its own row with one division, and all red values can be updated at the same time. Then the same holds for black. This helps some iterative methods, which we have not covered yet. It is the 2D version of an odd-even ordering in 1D. The 9-point stencil breaks this property, because diagonal neighbors have the same color.

The order changes the matrix, but not the solution. The experiment checks both facts.

### Accuracy in 2D

As in Lecture 3, the local truncation error of the 5-point stencil is

$$
\tau_{i,j} = \frac{h^2}{12} \big( u_{xxxx} + u_{yyyy} \big) + O(h^4),
$$

so the global error is $O(h^2)$. Note the fourth derivatives. The 1D analysis of Lecture 2 also had $u''''$ here. Compared with 1D, there is only one new term, $u_{yyyy}$.

## 5. The advection-diffusion equation

Now the new topic. In one space dimension, the **advection-diffusion equation** is

$$
u_t + \underbrace{a\, u_x}_{\text{advection}} - \underbrace{\kappa\, u_{xx}}_{\text{diffusion}} = \psi, \qquad \kappa > 0.
$$

The lecture spent a lot of time on the role of each term, and this is the main lesson: **before we design a scheme, we should understand what each term does.** A numerical scheme imitates a real process. If it ignores what the process does, it can give numbers with no meaning.

Think of a long, thin river, and let $u$ be the amount of dye in the water.

- **Advection** ($a\,u_x$) carries the dye with the water. If $a > 0$, the water flows in the positive direction, to the right. If $a < 0$, it flows to the left. (Some fields use "advection" for horizontal and "convection" for vertical motion. Here both just mean transport from one place to another.)
- **Diffusion** ($\kappa\, u_{xx}$) spreads the dye even in still water. A drop of dye becomes a wider and lower cloud. A large $\kappa$ means fast spreading and a flat profile. A small $\kappa$ means the dye stays in a tall, narrow peak. This is the flattening effect of the heat equation from Lecture 1.

**Boundary conditions.** Say $a > 0$ on $0 < x < 1$. The water enters at $x = 0$, so it is natural to give the value there:

$$
u(0, t) = \alpha(t) \qquad \text{(the inflow boundary).}
$$

With advection alone, this would be enough: everything inside comes from the left. But diffusion also spreads dye against the flow, so we also need a condition at the outflow boundary:

$$
u(1, t) = \beta(t) \qquad \text{(the outflow boundary).}
$$

**Steady state and one parameter.** Drop the time derivative. With $a > 0$, divide by $a$:

$$
-\varepsilon\, u''(x) + u'(x) = f(x), \qquad \varepsilon = \frac{\kappa}{a}, \quad f = \frac{\psi}{a}.
$$

This trick is called **parameter reduction**. Two parameters $a$ and $\kappa$ become one. And it makes sense: "advection is strong" or "diffusion is strong" only has meaning relative to the other term. Now only the size of $\varepsilon$ matters. A small $\varepsilon$ means advection dominates.

Be careful with signs here. Some references multiply by $-1$ and write $\varepsilon u'' - u' = -f$. This is the same equation. Any form is fine, as long as we use one form consistently.

## 6. The role of $\varepsilon$

### Diffusion only

Start without advection:

$$
-\varepsilon u'' = 1 \text{ on } (0, 1), \qquad u(0) = u(1) = 0.
$$

Integrate twice and use the boundary conditions:

$$
u(x) = \frac{1}{2\varepsilon}\, x(1 - x), \qquad \max u = u\big(\tfrac12\big) = \frac{1}{8\varepsilon}.
$$

This is a symmetric bump. A large $\varepsilon$ gives a low, flat bump. A small $\varepsilon$ gives a tall one. The formula $1/(8\varepsilon)$ is the flattening effect in one line.

### Adding advection

Now add the advection term:

$$
-\varepsilon u'' + u' = 1 \text{ on } (0, 1), \qquad u(0) = u(1) = 0.
$$

The general solution is $u = x + A + B e^{x/\varepsilon}$. The two boundary conditions fix $A$ and $B$, and we get

$$
u(x) = x - \frac{e^{x/\varepsilon} - 1}{e^{1/\varepsilon} - 1}.
$$

(To compute it for small $\varepsilon$, divide the top and bottom by $e^{1/\varepsilon}$ first, or $e^{1/\varepsilon}$ overflows.)

- For a **large** $\varepsilon$, diffusion wins. The solution is almost the symmetric bump, with its top moved a little to the right. For example, with $\varepsilon = 10$ the peak is $0.012498$, against $1/(8\varepsilon) = 0.0125$ for diffusion alone.
- For a **small** $\varepsilon$, advection wins. The solution is close to the line $u = x$ on most of the interval, and then drops sharply to $0$ near $x = 1$.

This is also a first look at **continuous dependence on the data**. If we change $\varepsilon$ a little, the solution changes only a little. For $\varepsilon = 10^7$, the term $u'$ is tiny next to the diffusion term, and the solution is almost exactly the pure diffusion bump.

## 7. Which boundary condition is the real one?

Now let $\varepsilon \to 0$. The equation becomes $u' = 1$, so $u = x + C$. This is a first-order equation, and it can only satisfy **one** boundary condition. But we have two. The problem is over-specified, and we must choose:

$$
u(0) = 0 \;\Rightarrow\; u = x, \qquad\text{or}\qquad u(1) = 0 \;\Rightarrow\; u = x - 1.
$$

Which one is right? **The one on the side the flow comes from.** Here $a > 0$, so the water flows from left to right. The dye in the interval came in at $x = 0$. The value at $x = 1$ is where the water leaves, so it cannot tell us what is coming. The real limit is $u = x$. The other line, $u = x - 1$, solves the equation, but it does not describe the physics.

This is the idea of **upwind** (or upstream). The words come from rivers and wind: to know what is coming to you, look in the direction the flow comes from. If the flow changes direction, the side you must look at changes too.

Two arguments support this answer.

1. **Continuity.** For a large $\varepsilon$, the solution is a positive bump. When we make $\varepsilon$ smaller step by step, the solution changes continuously. It cannot suddenly jump down to the negative line $u = x - 1$. The experiment checks this.
2. **The maximum principle.** This is a strong tool for this kind of equation. The lecture only named it, because it is outside the course.

### But $\varepsilon$ is never zero

However small $\varepsilon$ is, it is not zero. So the equation is still second order, it needs both boundary conditions, and a solution that meets both exists. In this sense, diffusion is always good news: it makes the solution smooth and the problem well-posed. But the two conditions do not play the same role. The inflow condition decides the solution on almost the whole interval. The outflow condition is only needed by the mathematics, and the solution meets it at the last moment.

## 8. Boundary layers

Look again at the solution for small $\varepsilon$, term by term.

- **Inside the interval**, $u \approx x$. The slope is almost constant, so $u''$ is small, and $\varepsilon u''$ can be ignored. The equation is almost $u' = 1$.
- **Near $x = 1$**, the solution must fall from about $1$ to $0$ in a very short distance. There $u''$ is huge, of size $1/\varepsilon^2$, so $\varepsilon u''$ is no longer small. This term is what lets the solution reach the outflow condition.

A **boundary layer** is a small region where the derivatives change a lot. Here it sits at the outflow boundary. It is also called a convective layer, because advection creates it. Its width is $O(\varepsilon)$: near $x = 1$ the correction term behaves like $e^{(x - 1)/\varepsilon}$, which decays over a distance of a few $\varepsilon$.

An equation like $-\varepsilon u'' + u' = f$ with $0 < \varepsilon \ll 1$ is called **singularly perturbed**. A small change (from $\varepsilon = 0$ to a small $\varepsilon > 0$) changes the type of the equation: from first order to second order, with a different number of boundary conditions. This makes it hard for numerical methods. The solution changes very fast in a very thin region, and a grid that does not resolve that region cannot capture it.

### The direction of the flow

The same holds for $-\varepsilon u'' + \beta u' = 1$ with a general $\beta$:

| $\beta$ | shape of the solution |
|---|---|
| $\beta = 0$ | pure diffusion: a symmetric bump |
| $\beta > 0$ | flow to the right: the bump leans right, layer at $x = 1$ |
| $\beta < 0$ | flow to the left: the bump leans left, layer at $x = 0$ |
| $\varepsilon \ll \lvert\beta\rvert$ | advection dominates: a sharp layer |

### The real world: flow past a wing

The steady incompressible Navier-Stokes equations,

$$
-\frac{1}{\mathrm{Re}} \nabla^2 \mathbf{u} + (\mathbf{u} \cdot \nabla)\mathbf{u} + \nabla p = \mathbf{f}, \qquad \nabla \cdot \mathbf{u} = 0,
$$

have the same two parts. The first term is diffusion, and the Reynolds number $\mathrm{Re}$ plays the role of $1/\varepsilon$. The second term is advection. But now the velocity that carries the flow is the solution $\mathbf{u}$ itself, so the equation is nonlinear.

Think of air flowing past an airplane wing. A little away from the wing, the air moves fast and smoothly. On the surface of the wing, the air does not move at all relative to the wing (the no-slip condition, $\mathbf{u} = 0$). In between, the speed changes from zero to very fast in a very thin region. That is a boundary layer, the same thing as in our 1D toy problem. In 2D, and with the nonlinearity, it becomes a very hard problem, both in mathematics and in physics.

### A warning about centered differences

So far the centered difference has always been our friend: it is second order and has no bias to either side. Here that becomes a problem. It treats the information from the left and from the right as equally important. But when advection dominates, only one side carries the real information. A scheme that listens to both sides equally can produce a solution that jumps up and down from point to point. These are called **spurious oscillations**.

The fix and its analysis come in the next lecture. For now, one rule of thumb from this lecture: a first-order equation can have spurious solutions, like $u = x - 1$ above, and adding a second-order (diffusion) term with the right sign picks the physical one.

## 9. Experiment

The script has four parts. Parts 1 to 3 solve 2D problems on the unit square with $n = 4, 8, 16, 32, 64$ intervals each way. They use numpy for the linear systems. Part 4 uses the exact formula from section 6.

**Part 1: 5-point and 9-point stencils.** First, the Laplace equation with Dirichlet data from $u = e^x \sin y$, which is harmonic. Then a Poisson problem with $u = e^{x + y}$ and $f = 2e^{x + y}$. For the 9-point stencil, the right-hand side is just $f$ at the grid point. The table shows the largest error at the grid points.

| $n$ | Laplace, 5-point | Laplace, 9-point | Poisson, 5-point | Poisson, 9-point |
|---|---|---|---|---|
| 4 | $5.79 \times 10^{-4}$ | $4.71 \times 10^{-9}$ | $2.05 \times 10^{-3}$ | $4.32 \times 10^{-3}$ |
| 8 | $1.60 \times 10^{-4}$ | $7.86 \times 10^{-11}$ | $5.44 \times 10^{-4}$ | $1.10 \times 10^{-3}$ |
| 16 | $4.07 \times 10^{-5}$ | $1.24 \times 10^{-12}$ | $1.40 \times 10^{-4}$ | $2.81 \times 10^{-4}$ |
| 32 | $1.02 \times 10^{-5}$ | $2.60 \times 10^{-14}$ | $3.51 \times 10^{-5}$ | $7.02 \times 10^{-5}$ |
| 64 | $2.55 \times 10^{-6}$ | $3.09 \times 10^{-14}$ | $8.79 \times 10^{-6}$ | $1.76 \times 10^{-5}$ |

- The 5-point stencil is second order in both problems (measured order $2.00$ from $n = 32$ to $64$).
- For the Poisson problem, the 9-point stencil is also only second order ($2.00$), and its error is about twice as large as the 5-point error. More points did not help.
- For the Laplace equation, the 9-point stencil is far better. The measured orders are $5.90$, $5.99$, and $5.57$, until the error hits round-off near $10^{-14}$ at $n = 32$. That is sixth order, not fourth.

Why sixth order? This explanation is my own. The script applies the 9-point weights, in exact integer arithmetic, to every harmonic polynomial up to degree $9$ (the real and imaginary parts of $(x + iy)^k$). The result is exactly $0$ for every degree up to $7$. The first nonzero result is at degree $8$. So when we Taylor expand a harmonic $u$, all the terms up to degree $7$ cancel, and the truncation error starts with $h^8 / h^2 = h^6$. For a general $f$, the term $\frac{h^2}{12} \nabla^2 f$ appears in the error, and we are back to second order.

**Part 2: two orderings.** With $h = 1/5$ there are $16$ interior points, $8$ red and $8$ black. In red-black order, the red-red and black-black blocks have $0$ nonzero entries off the diagonal, and every diagonal entry is $-4/h^2$. No red point has a red neighbor. The lower-left block is exactly the transpose of $H$, and the rows of $H$ have $2$, $3$, or $4$ nonzero entries (fewer at the boundary). Both orders give a matrix with the same $64$ nonzero entries, just in different places. Solving both systems gives the same grid values, up to $1.8 \times 10^{-15}$.

**Part 3: Neumann in 2D.** The solution is $u = e^{x + y}$ with $\partial u / \partial \mathbf{n} = g$ on the left and right sides, and Dirichlet data on the bottom and top. At every Neumann node, the boundary row is either the ghost node row from section 2, or the one-sided first-order row $(u_{0,j} - u_{1,j})/h = g(0, y_j)$ (and the mirror image on the right).

| $n$ | unknowns | ghost node | one-sided, 1st order |
|---|---|---|---|
| 4 | 15 | $1.13 \times 10^{-2}$ | $3.06 \times 10^{-1}$ |
| 8 | 63 | $3.20 \times 10^{-3}$ | $1.32 \times 10^{-1}$ |
| 16 | 255 | $8.18 \times 10^{-4}$ | $6.10 \times 10^{-2}$ |
| 32 | 1023 | $2.06 \times 10^{-4}$ | $2.94 \times 10^{-2}$ |
| 64 | 4095 | $5.15 \times 10^{-5}$ | $1.44 \times 10^{-2}$ |

The measured order from $n = 32$ to $64$ is $1.998$ for the ghost node and $1.028$ for the one-sided row. So the 1D result of Lecture 3 carries over to 2D: one first-order boundary row per grid row makes the whole solution first order. At $n = 64$ the ghost node error is about $280$ times smaller. The number of unknowns is $(n + 1)(n - 1)$, because the boundary values on the Neumann sides are unknown.

**Part 4: advection-diffusion.** The exact solution of $-\varepsilon u'' + u' = 1$, $u(0) = u(1) = 0$:

<figure class="viz">
<div class="viz-legend"><span class="viz-key"><span class="viz-swatch" style="background:#2a78d6"></span>ε = 1</span><span class="viz-key"><span class="viz-swatch" style="background:#eb6834"></span>ε = 0.1</span><span class="viz-key"><span class="viz-swatch" style="background:#1baf7a"></span>ε = 0.03</span><span class="viz-key"><span class="viz-swatch" style="background:#eda100"></span>ε = 0.01</span><span class="viz-key"><span class="viz-swatch" style="background:#6b6b6b"></span>Reduced solutions (dashed)</span></div>
<div class="viz-frame">
<svg viewBox="0 0 680 420" role="img" tabindex="0" aria-label="Exact solutions of minus epsilon u double prime plus u prime equals 1 with u(0) = u(1) = 0, for four values of epsilon. For epsilon = 1 the curve is a low bump. As epsilon gets smaller, the curve follows the line u = x and drops to zero in a thin layer near x = 1. It never comes close to the other candidate line u = x minus 1." class="viz-plot"><line x1="64" x2="540" y1="30.0" y2="30.0" class="viz-grid"/><text x="54" y="34.0" text-anchor="end" class="viz-tick">1</text><line x1="64" x2="540" y1="112.5" y2="112.5" class="viz-grid"/><text x="54" y="116.5" text-anchor="end" class="viz-tick">0.5</text><line x1="64" x2="540" y1="195.0" y2="195.0" class="viz-grid"/><text x="54" y="199.0" text-anchor="end" class="viz-tick">0</text><line x1="64" x2="540" y1="277.5" y2="277.5" class="viz-grid"/><text x="54" y="281.5" text-anchor="end" class="viz-tick">−0.5</text><line x1="64" x2="540" y1="360.0" y2="360.0" class="viz-grid"/><text x="54" y="364.0" text-anchor="end" class="viz-tick">−1</text><text x="64.0" y="382" text-anchor="middle" class="viz-tick">0</text><text x="183.0" y="382" text-anchor="middle" class="viz-tick">0.25</text><text x="302.0" y="382" text-anchor="middle" class="viz-tick">0.5</text><text x="421.0" y="382" text-anchor="middle" class="viz-tick">0.75</text><text x="540.0" y="382" text-anchor="middle" class="viz-tick">1</text><text x="302.0" y="408" text-anchor="middle" class="viz-axis-title">x</text><text x="64" y="14" class="viz-axis-title">u(x)</text><line x1="64.0" y1="195.0" x2="540.0" y2="30.0" stroke="#6b6b6b" stroke-width="1.5" stroke-dasharray="5 4"/><line x1="64.0" y1="360.0" x2="540.0" y2="195.0" stroke="#6b6b6b" stroke-width="1.5" stroke-dasharray="5 4"/><text x="550" y="22.0" class="viz-label">u = x (inflow side)</text><text x="278.2" y="307.8" class="viz-label">u = x − 1 (outflow side)</text><path d="M64.0,195.0 L68.8,194.3 L73.5,193.6 L78.3,193.0 L83.0,192.3 L87.8,191.7 L92.6,191.0 L97.3,190.4 L102.1,189.8 L106.8,189.2 L111.6,188.6 L116.4,188.0 L121.1,187.4 L125.9,186.9 L130.6,186.3 L135.4,185.8 L140.2,185.3 L144.9,184.7 L149.7,184.2 L154.4,183.7 L159.2,183.3 L164.0,182.8 L168.7,182.3 L173.5,181.9 L178.2,181.4 L183.0,181.0 L187.8,180.6 L192.5,180.2 L197.3,179.8 L202.0,179.5 L206.8,179.1 L211.6,178.7 L216.3,178.4 L221.1,178.1 L225.8,177.8 L230.6,177.5 L235.4,177.2 L240.1,176.9 L244.9,176.7 L249.6,176.5 L254.4,176.2 L259.2,176.0 L263.9,175.8 L268.7,175.6 L273.4,175.5 L278.2,175.3 L283.0,175.2 L287.7,175.1 L292.5,175.0 L297.2,174.9 L302.0,174.8 L306.8,174.7 L311.5,174.7 L316.3,174.7 L321.0,174.7 L325.8,174.7 L330.6,174.7 L335.3,174.7 L340.1,174.8 L344.8,174.9 L349.6,174.9 L354.4,175.1 L359.1,175.2 L363.9,175.3 L368.6,175.5 L373.4,175.7 L378.2,175.9 L382.9,176.1 L387.7,176.3 L392.4,176.6 L397.2,176.8 L402.0,177.1 L406.7,177.5 L411.5,177.8 L416.2,178.1 L421.0,178.5 L425.8,178.9 L430.5,179.3 L435.3,179.8 L440.0,180.2 L444.8,180.7 L449.6,181.2 L454.3,181.7 L459.1,182.2 L463.8,182.8 L468.6,183.4 L473.4,184.0 L478.1,184.6 L482.9,185.3 L487.6,186.0 L492.4,186.7 L497.2,187.4 L501.9,188.1 L506.7,188.9 L511.4,189.7 L516.2,190.5 L521.0,191.4 L525.7,192.2 L530.5,193.1 L535.2,194.1 L540.0,195.0" fill="none" stroke="#2a78d6" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><path d="M64.0,195.0 L68.8,193.4 L73.5,191.7 L78.3,190.1 L83.0,188.4 L87.8,186.8 L92.6,185.1 L97.3,183.5 L102.1,181.8 L106.8,180.2 L111.6,178.5 L116.4,176.9 L121.1,175.2 L125.9,173.6 L130.6,171.9 L135.4,170.3 L140.2,168.6 L144.9,167.0 L149.7,165.3 L154.4,163.7 L159.2,162.0 L164.0,160.4 L168.7,158.8 L173.5,157.1 L178.2,155.5 L183.0,153.8 L187.8,152.2 L192.5,150.6 L197.3,148.9 L202.0,147.3 L206.8,145.6 L211.6,144.0 L216.3,142.4 L221.1,140.7 L225.8,139.1 L230.6,137.5 L235.4,135.9 L240.1,134.2 L244.9,132.6 L249.6,131.0 L254.4,129.4 L259.2,127.8 L263.9,126.2 L268.7,124.6 L273.4,123.0 L278.2,121.4 L283.0,119.8 L287.7,118.3 L292.5,116.7 L297.2,115.1 L302.0,113.6 L306.8,112.1 L311.5,110.6 L316.3,109.0 L321.0,107.6 L325.8,106.1 L330.6,104.6 L335.3,103.2 L340.1,101.8 L344.8,100.4 L349.6,99.0 L354.4,97.7 L359.1,96.4 L363.9,95.1 L368.6,93.9 L373.4,92.7 L378.2,91.6 L382.9,90.5 L387.7,89.5 L392.4,88.6 L397.2,87.7 L402.0,86.9 L406.7,86.2 L411.5,85.6 L416.2,85.1 L421.0,84.8 L425.8,84.6 L430.5,84.5 L435.3,84.6 L440.0,84.8 L444.8,85.3 L449.6,86.0 L454.3,87.0 L459.1,88.2 L463.8,89.7 L468.6,91.6 L473.4,93.8 L478.1,96.4 L482.9,99.5 L487.6,103.1 L492.4,107.2 L497.2,111.9 L501.9,117.3 L506.7,123.5 L511.4,130.5 L516.2,138.3 L521.0,147.2 L525.7,157.2 L530.5,168.4 L535.2,180.9 L540.0,195.0" fill="none" stroke="#eb6834" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><path d="M64.0,195.0 L68.8,193.4 L73.5,191.7 L78.3,190.1 L83.0,188.4 L87.8,186.8 L92.6,185.1 L97.3,183.5 L102.1,181.8 L106.8,180.2 L111.6,178.5 L116.4,176.9 L121.1,175.2 L125.9,173.6 L130.6,171.9 L135.4,170.3 L140.2,168.6 L144.9,167.0 L149.7,165.3 L154.4,163.7 L159.2,162.0 L164.0,160.4 L168.7,158.7 L173.5,157.1 L178.2,155.4 L183.0,153.8 L187.8,152.1 L192.5,150.5 L197.3,148.8 L202.0,147.2 L206.8,145.5 L211.6,143.9 L216.3,142.2 L221.1,140.6 L225.8,138.9 L230.6,137.3 L235.4,135.6 L240.1,134.0 L244.9,132.3 L249.6,130.7 L254.4,129.0 L259.2,127.4 L263.9,125.7 L268.7,124.1 L273.4,122.4 L278.2,120.8 L283.0,119.1 L287.7,117.5 L292.5,115.8 L297.2,114.2 L302.0,112.5 L306.8,110.9 L311.5,109.2 L316.3,107.6 L321.0,105.9 L325.8,104.3 L330.6,102.6 L335.3,101.0 L340.1,99.3 L344.8,97.7 L349.6,96.0 L354.4,94.4 L359.1,92.7 L363.9,91.1 L368.6,89.4 L373.4,87.8 L378.2,86.1 L382.9,84.5 L387.7,82.8 L392.4,81.2 L397.2,79.5 L402.0,77.9 L406.7,76.2 L411.5,74.6 L416.2,72.9 L421.0,71.3 L425.8,69.7 L430.5,68.0 L435.3,66.4 L440.0,64.8 L444.8,63.2 L449.6,61.6 L454.3,60.1 L459.1,58.6 L463.8,57.2 L468.6,55.9 L473.4,54.7 L478.1,53.6 L482.9,52.8 L487.6,52.4 L492.4,52.4 L497.2,53.1 L501.9,54.7 L506.7,57.6 L511.4,62.2 L516.2,69.4 L521.0,80.1 L525.7,95.7 L530.5,118.0 L535.2,149.9 L540.0,195.0" fill="none" stroke="#1baf7a" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><path d="M64.0,195.0 L68.8,193.3 L73.5,191.7 L78.3,190.0 L83.0,188.4 L87.8,186.8 L92.6,185.1 L97.3,183.4 L102.1,181.8 L106.8,180.2 L111.6,178.5 L116.4,176.8 L121.1,175.2 L125.9,173.6 L130.6,171.9 L135.4,170.2 L140.2,168.6 L144.9,166.9 L149.7,165.3 L154.4,163.7 L159.2,162.0 L164.0,160.3 L168.7,158.7 L173.5,157.1 L178.2,155.4 L183.0,153.8 L187.8,152.1 L192.5,150.4 L197.3,148.8 L202.0,147.1 L206.8,145.5 L211.6,143.8 L216.3,142.2 L221.1,140.5 L225.8,138.9 L230.6,137.2 L235.4,135.6 L240.1,133.9 L244.9,132.3 L249.6,130.6 L254.4,129.0 L259.2,127.4 L263.9,125.7 L268.7,124.1 L273.4,122.4 L278.2,120.8 L283.0,119.1 L287.7,117.5 L292.5,115.8 L297.2,114.2 L302.0,112.5 L306.8,110.8 L311.5,109.2 L316.3,107.5 L321.0,105.9 L325.8,104.2 L330.6,102.6 L335.3,101.0 L340.1,99.3 L344.8,97.7 L349.6,96.0 L354.4,94.4 L359.1,92.7 L363.9,91.1 L368.6,89.4 L373.4,87.8 L378.2,86.1 L382.9,84.5 L387.7,82.8 L392.4,81.2 L397.2,79.5 L402.0,77.9 L406.7,76.2 L411.5,74.6 L416.2,72.9 L421.0,71.3 L425.8,69.6 L430.5,68.0 L435.3,66.3 L440.0,64.7 L444.8,63.0 L449.6,61.4 L454.3,59.7 L459.1,58.1 L463.8,56.4 L468.6,54.8 L473.4,53.1 L478.1,51.5 L482.9,49.8 L487.6,48.2 L492.4,46.5 L497.2,44.9 L501.9,43.3 L506.7,41.7 L511.4,40.3 L516.2,39.4 L521.0,39.6 L525.7,43.2 L530.5,55.6 L535.2,92.4 L540.0,195.0" fill="none" stroke="#eda100" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><text x="302.0" y="160.8" text-anchor="middle" class="viz-label">ε = 1</text><text x="359.1" y="80.4" text-anchor="middle" class="viz-label">ε = 0.1</text><text x="444.8" y="49.2" text-anchor="middle" class="viz-label">ε = 0.03</text><text x="473.4" y="29.1" text-anchor="middle" class="viz-label">ε = 0.01</text><line x1="0" x2="0" y1="30" y2="360" class="viz-crosshair" visibility="hidden"/></svg>
<div class="viz-tooltip" hidden></div>
</div>
<figcaption>Exact solutions of −εu″ + u′ = 1 on (0, 1) with u(0) = u(1) = 0. The dashed lines are the two solutions of the reduced equation u′ = 1 that keep one boundary condition each. Hover over the chart, or focus it and use the arrow keys, to read the values.</figcaption>
<details class="viz-table"><summary>Show the numbers</summary><table><thead><tr><th>x</th><th>ε = 1</th><th>ε = 0.1</th><th>ε = 0.03</th><th>ε = 0.01</th></tr></thead><tbody><tr><td>0.0</td><td>0.000</td><td>0.000</td><td>0.000</td><td>0.000</td></tr><tr><td>0.1</td><td>0.039</td><td>0.100</td><td>0.100</td><td>0.100</td></tr><tr><td>0.2</td><td>0.071</td><td>0.200</td><td>0.200</td><td>0.200</td></tr><tr><td>0.3</td><td>0.096</td><td>0.299</td><td>0.300</td><td>0.300</td></tr><tr><td>0.4</td><td>0.114</td><td>0.398</td><td>0.400</td><td>0.400</td></tr><tr><td>0.5</td><td>0.122</td><td>0.493</td><td>0.500</td><td>0.500</td></tr><tr><td>0.6</td><td>0.122</td><td>0.582</td><td>0.600</td><td>0.600</td></tr><tr><td>0.7</td><td>0.110</td><td>0.650</td><td>0.700</td><td>0.700</td></tr><tr><td>0.8</td><td>0.087</td><td>0.665</td><td>0.799</td><td>0.800</td></tr><tr><td>0.9</td><td>0.051</td><td>0.532</td><td>0.864</td><td>0.900</td></tr><tr><td>1.0</td><td>0.000</td><td>0.000</td><td>0.000</td><td>0.000</td></tr><tr><td>0.95</td><td>0.027</td><td>0.343</td><td>0.761</td><td>0.943</td></tr><tr><td>0.98</td><td>0.011</td><td>0.161</td><td>0.467</td><td>0.845</td></tr><tr><td>0.99</td><td>0.006</td><td>0.085</td><td>0.273</td><td>0.622</td></tr></tbody></table></details>
<script>(() => {const fig = document.currentScript.closest('figure');const svg = fig.querySelector('svg');const line = svg.querySelector('.viz-crosshair');const tip = fig.querySelector('.viz-tooltip');const names = ["ε = 1", "ε = 0.1", "ε = 0.03", "ε = 0.01"]; const colors = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"];const rows = [["0.000", "0.000", "0.000", "0.000"], ["0.004", "0.010", "0.010", "0.010"], ["0.008", "0.020", "0.020", "0.020"], ["0.012", "0.030", "0.030", "0.030"], ["0.016", "0.040", "0.040", "0.040"], ["0.020", "0.050", "0.050", "0.050"], ["0.024", "0.060", "0.060", "0.060"], ["0.028", "0.070", "0.070", "0.070"], ["0.032", "0.080", "0.080", "0.080"], ["0.035", "0.090", "0.090", "0.090"], ["0.039", "0.100", "0.100", "0.100"], ["0.042", "0.110", "0.110", "0.110"], ["0.046", "0.120", "0.120", "0.120"], ["0.049", "0.130", "0.130", "0.130"], ["0.053", "0.140", "0.140", "0.140"], ["0.056", "0.150", "0.150", "0.150"], ["0.059", "0.160", "0.160", "0.160"], ["0.062", "0.170", "0.170", "0.170"], ["0.065", "0.180", "0.180", "0.180"], ["0.068", "0.190", "0.190", "0.190"], ["0.071", "0.200", "0.200", "0.200"], ["0.074", "0.210", "0.210", "0.210"], ["0.077", "0.220", "0.220", "0.220"], ["0.080", "0.230", "0.230", "0.230"], ["0.082", "0.240", "0.240", "0.240"], ["0.085", "0.249", "0.250", "0.250"], ["0.087", "0.259", "0.260", "0.260"], ["0.090", "0.269", "0.270", "0.270"], ["0.092", "0.279", "0.280", "0.280"], ["0.094", "0.289", "0.290", "0.290"], ["0.096", "0.299", "0.300", "0.300"], ["0.098", "0.309", "0.310", "0.310"], ["0.101", "0.319", "0.320", "0.320"], ["0.102", "0.329", "0.330", "0.330"], ["0.104", "0.339", "0.340", "0.340"], ["0.106", "0.349", "0.350", "0.350"], ["0.108", "0.358", "0.360", "0.360"], ["0.109", "0.368", "0.370", "0.370"], ["0.111", "0.378", "0.380", "0.380"], ["0.112", "0.388", "0.390", "0.390"], ["0.114", "0.398", "0.400", "0.400"], ["0.115", "0.407", "0.410", "0.410"], ["0.116", "0.417", "0.420", "0.420"], ["0.117", "0.427", "0.430", "0.430"], ["0.118", "0.436", "0.440", "0.440"], ["0.119", "0.446", "0.450", "0.450"], ["0.120", "0.456", "0.460", "0.460"], ["0.121", "0.465", "0.470", "0.470"], ["0.121", "0.475", "0.480", "0.480"], ["0.122", "0.484", "0.490", "0.490"], ["0.122", "0.493", "0.500", "0.500"], ["0.123", "0.503", "0.510", "0.510"], ["0.123", "0.512", "0.520", "0.520"], ["0.123", "0.521", "0.530", "0.530"], ["0.123", "0.530", "0.540", "0.540"], ["0.123", "0.539", "0.550", "0.550"], ["0.123", "0.548", "0.560", "0.560"], ["0.123", "0.556", "0.570", "0.570"], ["0.123", "0.565", "0.580", "0.580"], ["0.122", "0.573", "0.590", "0.590"], ["0.122", "0.582", "0.600", "0.600"], ["0.121", "0.590", "0.610", "0.610"], ["0.120", "0.598", "0.620", "0.620"], ["0.119", "0.605", "0.630", "0.630"], ["0.118", "0.613", "0.640", "0.640"], ["0.117", "0.620", "0.650", "0.650"], ["0.116", "0.627", "0.660", "0.660"], ["0.115", "0.633", "0.670", "0.670"], ["0.113", "0.639", "0.680", "0.680"], ["0.112", "0.645", "0.690", "0.690"], ["0.110", "0.650", "0.700", "0.700"], ["0.108", "0.655", "0.710", "0.710"], ["0.106", "0.659", "0.720", "0.720"], ["0.104", "0.663", "0.730", "0.730"], ["0.102", "0.666", "0.740", "0.740"], ["0.100", "0.668", "0.750", "0.750"], ["0.098", "0.669", "0.760", "0.760"], ["0.095", "0.670", "0.770", "0.770"], ["0.092", "0.669", "0.779", "0.780"], ["0.090", "0.668", "0.789", "0.790"], ["0.087", "0.665", "0.799", "0.800"], ["0.084", "0.660", "0.808", "0.810"], ["0.081", "0.655", "0.818", "0.820"], ["0.077", "0.647", "0.827", "0.830"], ["0.074", "0.638", "0.835", "0.840"], ["0.070", "0.627", "0.843", "0.850"], ["0.067", "0.613", "0.851", "0.860"], ["0.063", "0.598", "0.857", "0.870"], ["0.059", "0.579", "0.862", "0.880"], ["0.055", "0.557", "0.864", "0.890"], ["0.051", "0.532", "0.864", "0.900"], ["0.046", "0.503", "0.860", "0.910"], ["0.042", "0.471", "0.851", "0.920"], ["0.037", "0.433", "0.833", "0.929"], ["0.032", "0.391", "0.805", "0.938"], ["0.027", "0.343", "0.761", "0.943"], ["0.022", "0.290", "0.696", "0.942"], ["0.017", "0.229", "0.602", "0.920"], ["0.011", "0.161", "0.467", "0.845"], ["0.006", "0.085", "0.273", "0.622"], ["0.000", "0.000", "0.000", "0.000"]]; const heads = ["x = 0.00", "x = 0.01", "x = 0.02", "x = 0.03", "x = 0.04", "x = 0.05", "x = 0.06", "x = 0.07", "x = 0.08", "x = 0.09", "x = 0.10", "x = 0.11", "x = 0.12", "x = 0.13", "x = 0.14", "x = 0.15", "x = 0.16", "x = 0.17", "x = 0.18", "x = 0.19", "x = 0.20", "x = 0.21", "x = 0.22", "x = 0.23", "x = 0.24", "x = 0.25", "x = 0.26", "x = 0.27", "x = 0.28", "x = 0.29", "x = 0.30", "x = 0.31", "x = 0.32", "x = 0.33", "x = 0.34", "x = 0.35", "x = 0.36", "x = 0.37", "x = 0.38", "x = 0.39", "x = 0.40", "x = 0.41", "x = 0.42", "x = 0.43", "x = 0.44", "x = 0.45", "x = 0.46", "x = 0.47", "x = 0.48", "x = 0.49", "x = 0.50", "x = 0.51", "x = 0.52", "x = 0.53", "x = 0.54", "x = 0.55", "x = 0.56", "x = 0.57", "x = 0.58", "x = 0.59", "x = 0.60", "x = 0.61", "x = 0.62", "x = 0.63", "x = 0.64", "x = 0.65", "x = 0.66", "x = 0.67", "x = 0.68", "x = 0.69", "x = 0.70", "x = 0.71", "x = 0.72", "x = 0.73", "x = 0.74", "x = 0.75", "x = 0.76", "x = 0.77", "x = 0.78", "x = 0.79", "x = 0.80", "x = 0.81", "x = 0.82", "x = 0.83", "x = 0.84", "x = 0.85", "x = 0.86", "x = 0.87", "x = 0.88", "x = 0.89", "x = 0.90", "x = 0.91", "x = 0.92", "x = 0.93", "x = 0.94", "x = 0.95", "x = 0.96", "x = 0.97", "x = 0.98", "x = 0.99", "x = 1.00"];const left = 64, right = 540, steps = 100, width = 680;let current = -1;const show = (n) => { n = Math.max(0, Math.min(steps, n)); current = n; const x = left + (right - left) * n / steps; line.setAttribute('x1', x); line.setAttribute('x2', x); line.setAttribute('visibility', 'visible'); tip.replaceChildren(); const head = document.createElement('div'); head.className = 'viz-tooltip-head'; head.textContent = heads[n]; tip.append(head); names.forEach((name, i) => {  const row = document.createElement('div'); row.className = 'viz-tooltip-row';  const key = document.createElement('span'); key.className = 'viz-swatch'; key.style.background = colors[i];  const value = document.createElement('strong'); value.textContent = rows[n][i];  const label = document.createElement('span'); label.textContent = name;  row.append(key, value, label); tip.append(row); }); const box = svg.getBoundingClientRect(), frame = fig.getBoundingClientRect(); const xs = box.left - frame.left + x * box.width / width; tip.hidden = false; const room = frame.width - tip.offsetWidth - 8; tip.style.left = Math.max(0, Math.min(room, xs + 12)) + 'px'; tip.style.top = (box.top - frame.top + 8) + 'px';};const hide = () => { line.setAttribute('visibility', 'hidden'); tip.hidden = true; current = -1; };svg.addEventListener('pointermove', (event) => { const box = svg.getBoundingClientRect(); const x = (event.clientX - box.left) * width / box.width; show(Math.round((x - left) / (right - left) * steps));});svg.addEventListener('pointerleave', hide);svg.addEventListener('focus', () => show(current < 0 ? 0 : current));svg.addEventListener('blur', hide);svg.addEventListener('keydown', (event) => { if (event.key === 'ArrowRight') { show(current + 1); event.preventDefault(); } if (event.key === 'ArrowLeft') { show(current - 1); event.preventDefault(); } if (event.key === 'Escape') hide();});})();</script>
</figure>

| $\varepsilon$ | $u(1/2)$ | peak at $x$ | peak value | $1/(8\varepsilon)$ | layer width / $\varepsilon$ |
|---|---|---|---|---|---|
| 10 | 0.012497 | 0.504 | 0.012498 | 0.0125 | — |
| 1 | 0.122459 | 0.541 | 0.123302 | 0.125 | — |
| 0.1 | 0.493307 | 0.770 | 0.669782 | 1.25 | 4.6007 |
| 0.01 | 0.500000 | 0.954 | 0.943948 | 12.5 | 4.6052 |
| 0.001 | 0.500000 | 0.993 | 0.992092 | 125 | 4.6052 |

What I noticed:

- **For a large $\varepsilon$, advection hardly matters.** With $\varepsilon = 10$ the peak is $0.012498$, and diffusion alone gives $0.0125$. The peak is at $x = 0.504$, just a little right of the middle.
- **For a small $\varepsilon$, the solution follows the inflow line $u = x$.** The value $u(1/2)$ goes to $0.5$, which is the value of $u = x$, not of $u = x - 1$. In a sweep of $401$ values of $\varepsilon$ from $10$ down to $10^{-3}$, $u(1/2)$ never decreased. It rose from $0.0125$ to $0.5$, and it never came near $-0.5$. This is the continuity argument in numbers.
- **The layer width is proportional to $\varepsilon$.** I defined the layer as the region where $u$ differs from $x$ by more than $0.01$. Then $u - x \approx -e^{(x - 1)/\varepsilon}$ predicts a width of $\varepsilon \ln 100 = 4.6052\,\varepsilon$ (my own estimate). The measured width divided by $\varepsilon$ is $4.6052$ for $\varepsilon = 0.01$ and $0.001$. At $\varepsilon = 0.1$ it is $4.6007$, a little off, because $e^{-1/\varepsilon}$ is not yet negligible.
- **The peak moves toward the outflow boundary.** Setting $u' = 0$ gives the peak at $x^* = 1 - \varepsilon \ln(1/\varepsilon)$ and the peak value $1 - \varepsilon \ln(1/\varepsilon) - \varepsilon$, up to terms of size $e^{-1/\varepsilon}$ (also my own estimate). For $\varepsilon = 0.01$ this gives $0.943948$, the same as the exact value to six digits. So the solution goes up to almost $1$ and falls to $0$ within a distance of a few $\varepsilon$.

The full script and the chart code are in [this site's repository](https://github.com/doubleyands/doubleyands.github.io/tree/main/scripts/posts/numerical-pde/lecture-4).

## Summary

1. The 5-point matrix can be negative definite ($\nabla^2$) or positive definite ($-\nabla^2$). The solution is the same. Check the **sign convention** of each reference.
2. Known boundary values move to the **right-hand side**: two values for a corner point, one for an edge point, none inside. The matrix is block tridiagonal, or **penta-diagonal** if we count entries.
3. A **Neumann condition in 2D** is the 1D ghost node at every boundary node. Use the outward normal, and keep the boundary values as unknowns.
4. The **9-point stencil** is second order for Poisson in general. For the Laplace equation it is much more accurate (sixth order in my test). The 5-point stencil is still the most popular.
5. **Red-black ordering** makes the diagonal blocks pure diagonal matrices, because the 5-point stencil only links points of different colors. The order changes the matrix, not the solution.
6. In the **advection-diffusion** equation, advection carries and diffusion spreads. Dividing by $a$ leaves one parameter $\varepsilon = \kappa / a$.
7. When advection dominates, the real information comes from the **upwind** side. The solution follows the inflow condition and meets the outflow condition in a thin **boundary layer** of width $O(\varepsilon)$.
8. A scheme must respect the direction of the flow. The next lecture shows what happens when it does not.
