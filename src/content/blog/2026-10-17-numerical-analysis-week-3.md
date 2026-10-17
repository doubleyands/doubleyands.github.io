---
title: "Numerical Analysis, Week 3: Cubic Splines and Least Squares Approximation"
description: "Cubic splines and their tridiagonal system, then best approximation in the L2 norm: normal equations, orthogonal polynomials, and Bessel's inequality. Tested on 1/(2+x) with three different bases."
pubDate: 2026-10-17
project: numerical-analysis
tags: [numerical analysis, splines, least squares, orthogonal polynomials, notes]
---

This is my summary of the third week of a numerical analysis course. As before, I write the ideas in my own words and check one of them with a small experiment. The course follows K. E. Atkinson, *An Introduction to Numerical Analysis*.

The week had two parts. The first part finishes interpolation. In [Week 2](/blog/2026-10-10-numerical-analysis-week-2/) we saw that one polynomial of high degree can behave badly. Here we cut the interval into pieces and use a low degree on each piece. The second part starts a new question. We stop asking the polynomial to pass through the data. Instead we ask it to be as close as possible to $f$ over the whole interval.

## 1. Piecewise polynomials

### Why pieces

Suppose we have 1,000 data points. One interpolating polynomial through all of them has degree 999. Week 2 showed that this is a bad idea: the error can grow near the ends (Runge's phenomenon), and the linear system is badly conditioned.

So we cut $[a, b]$ into small intervals

$$
a = x_0 < x_1 < \dots < x_n = b,
$$

and we use one polynomial of low degree on each interval $[x_{i-1}, x_i]$. This is a **piecewise polynomial**. When we add more points, we get more pieces, but the degree on each piece stays the same.

The definition does not ask the pieces to join. But a broken curve is not useful, so in practice we ask for continuity. The simplest example joins neighboring data points with straight lines. That is linear interpolation on each piece. It is continuous, but it has a corner at every node. We can also use Hermite interpolation (Week 2) on each piece, if we know the slopes.

### Splines

A **spline of order $m$** on the grid above is a function $s$ such that

- **(P1)** on each $[x_{i-1}, x_i]$, $s$ is a polynomial of degree less than $m$, and
- **(P2)** $s, s', \dots, s^{(m-2)}$ are continuous on all of $[a, b]$.

The most common case is $m = 4$, the **cubic spline**. Each piece has degree at most 3, and $s$, $s'$, and $s''$ are all continuous at every node. This is a strong condition. The joints cannot be seen: there is no jump, no corner, and no sudden change in curvature. And because every piece is only cubic, we do not get the wild swings of a high-degree polynomial.

### Counting conditions

On piece $i$ we write $s(x) = a_i + b_i x + c_i x^2 + d_i x^3$. Let us count.

| | count | where it comes from |
|---|---|---|
| unknowns | $4n$ | four coefficients on each of $n$ pieces |
| interpolation | $n+1$ | $s(x_i) = y_i$ for $i = 0, \dots, n$ |
| smoothness | $3(n-1)$ | at each interior node $x_1, \dots, x_{n-1}$: $s$, $s'$, $s''$ match from both sides |
| total conditions | $4n - 2$ | two short |

Smoothness gives conditions only at the interior nodes. At $x_0$ and $x_n$ there is a piece on one side only, so there is nothing to match. So two conditions are missing. We must add two **end conditions**, and different choices give different kinds of splines.

## 2. Building a cubic spline

Solving for $4n$ coefficients directly is possible, but there is a smarter choice of unknowns. Let

$$
M_i = s''(x_i), \qquad i = 0, 1, \dots, n,
$$

and $h_i = x_{i+1} - x_i$. Below, "the piece on $[x_i, x_{i+1}]$" is the piece to the right of $x_i$.

**Step 1: $s''$ is a straight line on each piece.** $s$ is cubic, so $s''$ is linear. On $[x_i, x_{i+1}]$ it goes through $(x_i, M_i)$ and $(x_{i+1}, M_{i+1})$. The Lagrange form for two points gives it at once:

$$
s''(x) = \frac{(x_{i+1} - x)\,M_i + (x - x_i)\,M_{i+1}}{h_i}.
$$

Both neighboring pieces use the same $M_i$ at $x_i$. So $s''$ is continuous for free.

**Step 2: integrate twice.**

$$
s(x) = \frac{(x_{i+1} - x)^3 M_i + (x - x_i)^3 M_{i+1}}{6 h_i} + C\,(x_{i+1} - x) + D\,(x - x_i).
$$

Two integrations add a linear function. Any linear function can be written as $C(x_{i+1} - x) + D(x - x_i)$. I like this form because the $D$ term is zero at $x_i$ and the $C$ term is zero at $x_{i+1}$. It is the same trick as in the Lagrange form.

**Step 3: use the data.** Put $x = x_i$, then $x = x_{i+1}$:

$$
y_i = \frac{h_i^2}{6} M_i + C h_i, \qquad y_{i+1} = \frac{h_i^2}{6} M_{i+1} + D h_i.
$$

Each equation has one unknown, so $C$ and $D$ follow directly. Now every piece passes through its two data points, so $s$ is continuous too.

**Step 4: ask for a continuous slope.** Only $s'$ is left. At an interior node $x_i$, the slope of the left piece must equal the slope of the right piece. Differentiate the formula from Step 2 on both pieces and put $x = x_i$:

$$
s'(x_i^-) = \frac{h_{i-1}}{6} M_{i-1} + \frac{h_{i-1}}{3} M_i + \frac{y_i - y_{i-1}}{h_{i-1}}, \qquad
s'(x_i^+) = -\frac{h_i}{3} M_i - \frac{h_i}{6} M_{i+1} + \frac{y_{i+1} - y_i}{h_i}.
$$

Set them equal and collect the $M$'s on the left:

$$
\frac{h_{i-1}}{6} M_{i-1} + \frac{h_{i-1} + h_i}{3} M_i + \frac{h_i}{6} M_{i+1} = \frac{y_{i+1} - y_i}{h_i} - \frac{y_i - y_{i-1}}{h_{i-1}}, \qquad i = 1, \dots, n-1.
$$

This is $n - 1$ equations for $n + 1$ unknowns $M_0, \dots, M_n$. Again two are missing, as we counted.

### End conditions

**Clamped spline.** If we know the slopes at the ends, we ask $s'(x_0) = y_0'$ and $s'(x_n) = y_n'$. The slope formulas from Step 4 give two more equations:

$$
\frac{h_0}{3} M_0 + \frac{h_0}{6} M_1 = \frac{y_1 - y_0}{h_0} - y_0', \qquad
\frac{h_{n-1}}{6} M_{n-1} + \frac{h_{n-1}}{3} M_n = y_n' - \frac{y_n - y_{n-1}}{h_{n-1}}.
$$

Now we have a square system $A\mathbf{M} = \mathbf{d}$:

$$
A = \begin{bmatrix}
\frac{h_0}{3} & \frac{h_0}{6} & & & \\
\frac{h_0}{6} & \frac{h_0 + h_1}{3} & \frac{h_1}{6} & & \\
 & \ddots & \ddots & \ddots & \\
 & & \frac{h_{n-2}}{6} & \frac{h_{n-2} + h_{n-1}}{3} & \frac{h_{n-1}}{6} \\
 & & & \frac{h_{n-1}}{6} & \frac{h_{n-1}}{3}
\end{bmatrix}.
$$

The known end slopes $y_0'$ and $y_n'$ go into the first and last entries of $\mathbf{d}$. The result is called the **complete** (or clamped) cubic spline.

**Other end conditions.** If we do not know the end slopes, there are other choices:

- **Natural spline**: $M_0 = M_n = 0$. The curve is straight at both ends.
- **Not-a-knot**: ask $s'''$ to be continuous at $x_1$ and $x_{n-1}$ too. $s'''$ is constant on each piece, so the first two pieces become one cubic, and the last two pieces too. This is the default in MATLAB's `spline` and SciPy's `CubicSpline`.

### Why this system is easy

Look at the matrix $A$. Three things stand out.

- **Tridiagonal.** Each equation only has $M_{i-1}$, $M_i$, $M_{i+1}$. So $A$ has nonzero entries only on the diagonal and next to it. Such a system can be solved in $O(n)$ operations (the Thomas algorithm). Compare this with the dense Vandermonde system from Week 2.
- **Diagonally dominant.** In an interior row, the diagonal is $\frac{h_{i-1} + h_i}{3}$ and the off-diagonal entries add up to $\frac{h_{i-1} + h_i}{6}$. The diagonal is exactly twice as big. In the first and last rows, $\frac{h}{3} > \frac{h}{6}$. This holds for any spacing, equal or not.
- **Symmetric and positive definite.** So the system has exactly one solution, and elimination is stable.

With equal spacing $h_i = h$, divide an interior equation by $h$:

$$
\frac{M_{i-1} + 4 M_i + M_{i+1}}{6} = \frac{y_{i+1} - 2 y_i + y_{i-1}}{h^2}.
$$

The right side is the usual central difference for a second derivative. The left side is a $1{:}4{:}1$ average of the spline's second derivatives. So the system says: "make the curvature of the spline match the second differences of the data."

## 3. From interpolation to approximation

### A different goal

Interpolation asks $p(x_i) = y_i$ exactly. That is not always what we want. Think of data that is almost flat with a little noise. A straight line through the middle is a good description. An interpolating polynomial must hit every point, so it wiggles between them.

So we change the goal. We give up passing through the data. We ask instead that $p$ be **close to $f$ over the whole interval**, where "close" is measured by a norm. This is where the norms from [Week 1](/blog/2026-10-03-numerical-analysis-week-1/) become the main tool.

| | interpolation | approximation |
|---|---|---|
| requirement | $p(x_i) = y_i$ at the nodes | $\lVert f - p \rVert$ as small as possible |
| depends on | the choice of nodes | the choice of norm |
| weak point | can swing at high degree | may not hit the data exactly |

### Weierstrass theorem

**Theorem (Weierstrass).** Let $f$ be continuous on $[a, b]$, and let $\epsilon > 0$. Then there is a polynomial $p$ with

$$
|f(x) - p(x)| \le \epsilon \quad \text{for all } x \in [a, b].
$$

In Week 1 notation, $\|f - p\|_\infty \le \epsilon$. The important word is *any* $\epsilon$. However small we want the error, some polynomial gets there. So approximating continuous functions by polynomials always makes sense.

The theorem needs $f$ to be continuous. A uniform limit of polynomials is always continuous, so a function with a jump cannot be approximated in the $\infty$-norm. In the $L^2$ norm it is different. Near a jump the error can stay large on a tiny interval, but that interval adds almost nothing to an integral. So in $L^2$ even a function with a jump can be approximated as well as we like.

### Taylor polynomials are not the answer

A Taylor polynomial is a natural first try, but it is built around **one point**. Take $f(x) = 1/(2+x)$ on $[-1, 1]$. Its cubic Taylor polynomial at $0$ is

$$
T_3(x) = \frac12 - \frac{x}{4} + \frac{x^2}{8} - \frac{x^3}{16}.
$$

The error is tiny near $0$ (only $0.0016$ at $x = 0.5$), but it reaches $0.0625$ at $x = -1$. The Taylor polynomial spends all its accuracy near the center and leaves the ends poor. We want a polynomial that spreads the error over the whole interval.

## 4. Least squares approximation

### The problem

We use the $L^2$ norm, possibly with a **weight function** $w(x) \ge 0$:

$$
\|g\|_2 = \left( \int_a^b w(x)\,g(x)^2\,dx \right)^{1/2}.
$$

**Problem.** Among all polynomials $r$ of degree at most $n$, find the one that makes $\|f - r\|_2$ smallest.

Two small remarks.

- Minimizing $\|f - r\|_2$ and minimizing $\|f - r\|_2^2$ give the same $r$. The square root is increasing, so it does not change which value is smallest. The square is easier to work with, so we use it.
- The weight squares nothing. It is already non-negative. Only $f - r$ is squared.

We ask two things of the weight:

1. $\int_a^b |x|^k w(x)\,dx$ is finite for every $k \ge 0$. (This allows infinite intervals.)
2. If $g \ge 0$ is continuous and $\int_a^b w\,g\,dx = 0$, then $g \equiv 0$. In words: $w$ does not switch off any part of the interval. If $w$ were $0$ on some subinterval, a bump $g$ inside it would break this rule.

The weight tells us where we care most:

| weight $w(x)$ | interval | focus | orthogonal family |
|---|---|---|---|
| $1$ | $[a, b]$ | the whole interval, evenly | Legendre |
| $1/\sqrt{1 - x^2}$ | $[-1, 1]$ | the two ends | Chebyshev |
| $e^{-x}$ | $[0, \infty)$ | near $0$ | Laguerre |
| $e^{-x^2}$ | $(-\infty, \infty)$ | near $0$, even more | Hermite |

The Chebyshev weight goes to infinity at $x = \pm 1$, but its integral is finite ($\pi$). So it passes condition 1. The last column will make sense in Section 5. (The Hermite family here is not the same thing as Hermite interpolation. They only share a name.)

### A first example

Take $w = 1$, $f(x) = 1/(2+x)$ on $[-1, 1]$, and a line $r(x) = a + bx$. We minimize

$$
F(a, b) = \int_{-1}^{1} \left( \frac{1}{2+x} - a - bx \right)^2 dx.
$$

After integrating, $x$ is gone. $F$ is a function of $a$ and $b$ only. Expand it:

$$
F(a, b) = \int_{-1}^1 f^2\,dx - 2a \int_{-1}^1 f\,dx - 2b \int_{-1}^1 x f\,dx + 2a^2 + \tfrac23 b^2 + 2ab \int_{-1}^1 x\,dx.
$$

The last integral is $0$, so there is no $ab$ term. The coefficients of $a^2$ and $b^2$ are positive. So $F$ is an upward bowl in $a$ and in $b$, and the point where both partial derivatives vanish is the minimum, not a maximum:

$$
\frac{\partial F}{\partial a} = 0 \;\Rightarrow\; a = \frac12 \int_{-1}^1 \frac{dx}{2+x} = \frac{\ln 3}{2} \approx 0.5493,
$$

$$
\frac{\partial F}{\partial b} = 0 \;\Rightarrow\; b = \frac32 \int_{-1}^1 \frac{x\,dx}{2+x} = 3 - 3\ln 3 \approx -0.2958.
$$

The $ab$ term vanished because $\int_{-1}^1 1 \cdot x\,dx = 0$. In the language of Section 5, $1$ and $x$ are **orthogonal** on $[-1, 1]$. Keep this in mind.

### The normal equations

Now the general case with $r(x) = \sum_{j=0}^n a_j x^j$. Let $g = f - r$, so $F = \int w\,g^2\,dx$. Only one term of $g$ contains $a_i$, so $\partial g / \partial a_i = -x^i$. Then

$$
0 = \frac{\partial F}{\partial a_i} = -2 \int_a^b w \left( f - \sum_{j=0}^n a_j x^j \right) x^i\,dx.
$$

Divide by $-2$, take the constants $a_j$ out of the integrals, and move the known part to the right:

$$
\sum_{j=0}^n a_j \int_a^b w\,x^{i+j}\,dx = \int_a^b w\,f\,x^i\,dx, \qquad i = 0, 1, \dots, n.
$$

These are the **normal equations**: $n + 1$ equations for $n + 1$ unknowns. The matrix entry in row $i$, column $j$ is $\int w\,x^{i+j}\,dx$. With $n = 1$, $w = 1$, and $[-1, 1]$, we get back the example above.

The system has exactly one solution. But it has two problems.

- **It is badly conditioned.** On $[0, 1]$ the functions $x^8$, $x^9$, $x^{10}$ look almost the same. So the rows of the matrix are almost the same. With $w = 1$ on $[0, 1]$ the entries are $1/(i + j + 1)$: this is the famous **Hilbert matrix**. Its condition number grows very fast. In my experiment below it is about $2.8 \times 10^4$ for $n = 3$ and $1.2 \times 10^{15}$ for $n = 10$.
- **It is dense.** Every entry is nonzero.

Both problems come from one choice: we wrote $r$ in the basis $1, x, x^2, \dots$. The cure is to choose a better basis.

## 5. Orthogonal polynomials

### A weighted inner product

Define

$$
(f, g) = \int_a^b w(x)\,f(x)\,g(x)\,dx, \qquad \|f\|_2 = \sqrt{(f, f)}.
$$

Is this really an inner product? It is symmetric, and it is linear in each slot because integrals are linear. It is non-negative, because $w \ge 0$ and $u^2 \ge 0$. The only real question is: does $(u, u) = 0$ force $u = 0$? Yes, and this is exactly condition 2 on the weight, used with $g = u^2$. The Cauchy–Schwarz inequality $|(f, g)| \le \|f\|_2 \|g\|_2$ also holds.

We call $f$ and $g$ **orthogonal** if $(f, g) = 0$. If also $\|f\| = \|g\| = 1$, they are **orthonormal**.

Orthogonality always belongs to an inner product. With $w = 1$, the functions $1$ and $x$ are orthogonal on $[-1, 1]$, because $\int_{-1}^1 x\,dx = 0$. But on $[0, 1]$ they are not, because $\int_0^1 x\,dx = \tfrac12$. Change the interval or the weight, and the answer changes.

### Gram–Schmidt

**Theorem.** For any weight as above, there is a sequence of polynomials $\phi_0, \phi_1, \phi_2, \dots$ where $\phi_n$ has degree exactly $n$ and

$$
(\phi_n, \phi_m) = 0 \quad (n \ne m), \qquad \|\phi_n\| = 1.
$$

The proof builds them one at a time.

- Start with a constant $\phi_0 = c$, and choose $c$ so that $\|\phi_0\| = 1$.
- Suppose $\phi_0, \dots, \phi_{n-1}$ are done. Take the new ingredient $x^n$ and remove its parts along the old ones:

$$
\psi_n = x^n - \sum_{j=0}^{n-1} (x^n, \phi_j)\,\phi_j.
$$

  Then $(\psi_n, \phi_k) = (x^n, \phi_k) - (x^n, \phi_k) = 0$ for every $k < n$. Divide by the norm: $\phi_n = \psi_n / \|\psi_n\|$.

This works because $1, x, x^2, \dots$ are linearly independent. No mix of $1$ and $x$ can make $x^2$, so $\psi_n$ is never zero.

**Example.** Take $w = 1$ on $[0, 1]$.

- $\phi_0 = 1$, since $\int_0^1 1\,dx = 1$.
- $\psi_1 = x - (x, 1) = x - \tfrac12$, and $\|\psi_1\|^2 = \int_0^1 (x - \tfrac12)^2 dx = \tfrac1{12}$. So $\phi_1 = \sqrt{3}\,(2x - 1)$.
- $(x^2, \phi_0) = \tfrac13$ and $(x^2, \phi_1) = \sqrt3 \int_0^1 (2x^3 - x^2)\,dx = \tfrac{\sqrt3}{6}$. So $\psi_2 = x^2 - \tfrac13 - \tfrac{\sqrt3}{6} \cdot \sqrt3\,(2x - 1) = x^2 - x + \tfrac16$, and $\|\psi_2\|^2 = \tfrac1{180}$. So $\phi_2 = \sqrt5\,(6x^2 - 6x + 1)$.

If we do the same on $[-1, 1]$, we get the Legendre polynomials (up to scaling).

Gram–Schmidt is not only a theory tool. It is a part of many numerical methods. For example, GMRES, a common iterative solver for large nonsymmetric systems $A\mathbf{x} = \mathbf{b}$, uses a Gram–Schmidt process to build its search directions.

### Three families

Each weight in the table of Section 4 has its own orthogonal family. We do not need closed formulas. Each family has a **three-term recurrence**: start from degree $0$ and $1$, and every next polynomial comes from the two before it.

- **Legendre** ($w = 1$ on $[-1, 1]$): $P_0 = 1$, $P_1 = x$,
  $$P_{n+1}(x) = \frac{2n+1}{n+1}\,x P_n(x) - \frac{n}{n+1}\,P_{n-1}(x).$$
  Then $(P_n, P_n) = \frac{2}{2n+1}$. So they are orthogonal but not orthonormal. $\sqrt{(2n+1)/2}\,P_n$ is orthonormal.
- **Chebyshev** ($w = 1/\sqrt{1 - x^2}$ on $[-1, 1]$): $T_0 = 1$, $T_1 = x$,
  $$T_{n+1}(x) = 2x\,T_n(x) - T_{n-1}(x).$$
  So $T_2 = 2x^2 - 1$ and $T_3 = 4x^3 - 3x$. Here $(T_0, T_0) = \pi$ and $(T_n, T_n) = \pi/2$ for $n \ge 1$. Chebyshev polynomials are very important in approximation theory.
- **Laguerre** ($w = e^{-x}$ on $[0, \infty)$): $L_0 = 1$, $L_1 = 1 - x$,
  $$L_{n+1}(x) = \frac{2n + 1 - x}{n+1}\,L_n(x) - \frac{n}{n+1}\,L_{n-1}(x).$$
  These are already orthonormal.

A quick check with $n = 1$: the Legendre recurrence gives $P_2 = \tfrac32 x^2 - \tfrac12 = \tfrac12(3x^2 - 1)$, and the Laguerre one gives $L_2 = 1 - 2x + \tfrac12 x^2$. Both families can be moved to another interval $[a, b]$ with the linear change $x = \frac{(b - a)t + (a + b)}{2}$.

A side note that connects to Week 2. In the interpolation error formula, the nodes only enter through $W(x) = (x - x_0)\cdots(x - x_n)$. Among all such products, the zeros of $T_{n+1}$ make $\max_{[-1,1]} |W|$ as small as possible. These zeros, $x_k = \cos\frac{(2k+1)\pi}{2(n+1)}$, crowd toward the ends. So they are the best version of "put more nodes near the ends", the cure we used in [Week 2](/blog/2026-10-10-numerical-analysis-week-2/).

### Two properties

**Expansion.** Let $\{\phi_n\}$ be orthogonal, and let $f$ be a polynomial of degree $m$. Then

$$
f = \sum_{k=0}^m \frac{(f, \phi_k)}{(\phi_k, \phi_k)}\,\phi_k.
$$

Proof: $\phi_0, \dots, \phi_m$ have degrees $0$ to $m$, so they can build any polynomial of degree $\le m$. Write $f = \sum \alpha_j \phi_j$ and take the inner product with $\phi_k$. By orthogonality, only $\alpha_k (\phi_k, \phi_k)$ is left. So $\alpha_k = (f, \phi_k) / (\phi_k, \phi_k)$.

A consequence: if $f$ has degree at most $m - 1$, then $(f, \phi_m) = 0$. In words, $\phi_m$ is orthogonal to every polynomial of lower degree.

**Roots.** $\phi_n$ has exactly $n$ distinct real roots, and all of them are inside the open interval $(a, b)$. For Chebyshev: $T_1 = x$ has the root $0$, and $T_2 = 2x^2 - 1$ has the roots $\pm 1/\sqrt2$.

## 6. Least squares with an orthonormal basis

Now let $\{\phi_j\}$ be **orthonormal**, and write $r = \sum_{j=0}^n a_j \phi_j$. We minimize $G(a_0, \dots, a_n) = \|f - r\|_2^2$. There are two ways.

**Way 1: partial derivatives.** As before,

$$
0 = \frac{\partial G}{\partial a_k} = -2(f, \phi_k) + 2 \sum_j a_j\,(\phi_j, \phi_k) = -2(f, \phi_k) + 2 a_k,
$$

because $(\phi_j, \phi_k)$ is $1$ when $j = k$ and $0$ otherwise. So

$$
a_k = (f, \phi_k).
$$

The "matrix" of the normal equations is now the identity. Each coefficient is one inner product. No system to solve.

**Way 2: complete the square.** Expand the norm:

$$
G = (f, f) - 2\sum_j a_j (f, \phi_j) + \sum_j a_j^2
= \|f\|_2^2 - \sum_j (f, \phi_j)^2 + \sum_j \bigl( a_j - (f, \phi_j) \bigr)^2.
$$

The last step adds and subtracts $(f, \phi_j)^2$ in each term. Only the last sum depends on the $a_j$, and it is a sum of squares. So $G$ is smallest exactly when every square is zero: $a_j = (f, \phi_j)$. This gives the same answer, and it also shows that the best approximation **exists and is unique**:

$$
r_n = \sum_{j=0}^n (f, \phi_j)\,\phi_j.
$$

If the $\phi_j$ are only orthogonal, use $a_j = (f, \phi_j)/(\phi_j, \phi_j)$ instead.

Two big gains over the monomial basis:

- **No linear system.** Each coefficient comes from one inner product.
- **Old coefficients stay.** To go from degree $n$ to $n + 1$, add one term: $r_{n+1} = r_n + (f, \phi_{n+1})\,\phi_{n+1}$. This is the same benefit as Newton's divided differences in Week 2. With monomials, every coefficient changes when the degree goes up.

### A picture: projection

Put $a_j = (f, \phi_j)$ into Way 2. The last sum disappears:

$$
\|f - r_n\|_2^2 = \|f\|_2^2 - \sum_{j=0}^n (f, \phi_j)^2 = \|f\|_2^2 - \|r_n\|_2^2.
$$

This is the Pythagorean theorem: $\|f\|^2 = \|r_n\|^2 + \|f - r_n\|^2$. Also, the error is orthogonal to every $\phi_k$ with $k \le n$, because $(f - r_n, \phi_k) = (f, \phi_k) - a_k = 0$. So $r_n$ is the **orthogonal projection** of $f$ onto the polynomials of degree $\le n$. In 3D, the closest point of a plane to a given point is the foot of the perpendicular. It is the same idea here.

### Bessel and Parseval

The left side above is never negative. So we get **Bessel's inequality**:

$$
\|r_n\|_2^2 = \sum_{j=0}^n (f, \phi_j)^2 \le \|f\|_2^2.
$$

The approximation is never "bigger" than $f$. If the interval $[a, b]$ is finite, more is true: $\|f - r_n\|_2 \to 0$ as $n \to \infty$. Then the inequality becomes **Parseval's equality**:

$$
\|f\|_2^2 = \sum_{j=0}^\infty (f, \phi_j)^2.
$$

## 7. Experiment: the same polynomial in three bases

The theory says that the best approximation $r_n$ does not depend on the basis. It is one polynomial. But the *computer* does not see the polynomial. It sees coefficients, and it must compute them. I wanted to see how much the basis matters.

**Setup.** $f(x) = 1/(2+x)$ on $[-1, 1]$, $w = 1$, degree $n = 0, 1, \dots, 30$. I computed $r_n$ in three ways:

1. **Legendre basis**: $a_j = (f, \phi_j)$ with orthonormal Legendre polynomials (built with the recurrence).
2. **Monomials on $[-1, 1]$**: solve the normal equations for $1, x, \dots, x^n$ with Gaussian elimination. The matrix entries $\int_{-1}^1 x^{i+j}\,dx$ are exact.
3. **Monomials on $[0, 1]$**: the same, but in the variable $t = (x + 1)/2 \in [0, 1]$. Now the matrix is the Hilbert matrix. The polynomial is the same; only the basis changes.

All integrals with $f$ use a 200-point Gauss–Legendre rule. $f$ is smooth on $[-1, 1]$ (its pole is at $x = -2$), so this is accurate to full double precision. I measured the $L^2$ error of each result with the same rule. The script uses only the Python standard library and runs in under a second.

<figure class="viz">
<div class="viz-legend"><span class="viz-key"><span class="viz-swatch" style="background:#2a78d6"></span>Legendre basis</span><span class="viz-key"><span class="viz-swatch" style="background:#eb6834"></span>monomials on [−1, 1]</span><span class="viz-key"><span class="viz-swatch" style="background:#1baf7a"></span>monomials on [0, 1] (Hilbert matrix)</span></div>
<div class="viz-frame">
<svg viewBox="0 0 680 400" role="img" tabindex="0" aria-label="L2 error of the best polynomial approximation to 1/(2+x) against degree, log scale. With the Legendre basis the error falls steadily to about 10^-15. With monomials on [-1, 1] it stops near 10^-10 after degree 19. With monomials on [0, 1] it stops near 10^-7 after degree 10." class="viz-plot"><line x1="64" x2="540" y1="30.0" y2="30.0" class="viz-grid"/><text x="54" y="34.0" text-anchor="end" class="viz-tick">1</text><line x1="64" x2="540" y1="107.5" y2="107.5" class="viz-grid"/><text x="54" y="111.5" text-anchor="end" class="viz-tick">10<tspan dy="-0.5em" font-size="75%">−4</tspan></text><line x1="64" x2="540" y1="185.0" y2="185.0" class="viz-grid"/><text x="54" y="189.0" text-anchor="end" class="viz-tick">10<tspan dy="-0.5em" font-size="75%">−8</tspan></text><line x1="64" x2="540" y1="262.5" y2="262.5" class="viz-grid"/><text x="54" y="266.5" text-anchor="end" class="viz-tick">10<tspan dy="-0.5em" font-size="75%">−12</tspan></text><line x1="64" x2="540" y1="340.0" y2="340.0" class="viz-grid"/><text x="54" y="344.0" text-anchor="end" class="viz-tick">10<tspan dy="-0.5em" font-size="75%">−16</tspan></text><text x="64.0" y="362" text-anchor="middle" class="viz-tick">0</text><text x="143.3" y="362" text-anchor="middle" class="viz-tick">5</text><text x="222.7" y="362" text-anchor="middle" class="viz-tick">10</text><text x="302.0" y="362" text-anchor="middle" class="viz-tick">15</text><text x="381.3" y="362" text-anchor="middle" class="viz-tick">20</text><text x="460.7" y="362" text-anchor="middle" class="viz-tick">25</text><text x="540.0" y="362" text-anchor="middle" class="viz-tick">30</text><text x="302.0" y="388" text-anchor="middle" class="viz-axis-title">Degree n</text><text x="64" y="14" class="viz-axis-title">L2 error on [−1, 1]</text><path d="M64.0,41.6 L79.9,52.4 L95.7,63.4 L111.6,74.4 L127.5,85.4 L143.3,96.5 L159.2,107.5 L175.1,118.6 L190.9,129.7 L206.8,140.7 L222.7,151.8 L238.5,162.9 L254.4,174.0 L270.3,185.0 L286.1,196.1 L302.0,207.2 L317.9,218.3 L333.7,229.3 L349.6,240.4 L365.5,251.5 L381.3,262.6 L397.2,273.7 L413.1,284.7 L428.9,295.8 L444.8,306.7 L460.7,316.1 L476.5,319.7 L492.4,320.2 L508.3,320.2 L524.1,320.2 L540.0,320.2" fill="none" stroke="#2a78d6" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><circle cx="540.0" cy="320.2" r="4" fill="#2a78d6" stroke="#ffffff" stroke-width="2"/><path d="M64.0,41.6 L79.9,52.4 L95.7,63.4 L111.6,74.4 L127.5,85.4 L143.3,96.5 L159.2,107.5 L175.1,118.6 L190.9,129.7 L206.8,140.7 L222.7,151.8 L238.5,162.9 L254.4,174.0 L270.3,185.0 L286.1,196.1 L302.0,207.2 L317.9,218.3 L333.7,229.3 L349.6,240.0 L365.5,242.4 L381.3,233.4 L397.2,225.5 L413.1,217.2 L428.9,209.3 L444.8,205.8 L460.7,207.3 L476.5,212.6 L492.4,208.7 L508.3,201.2 L524.1,201.4 L540.0,203.2" fill="none" stroke="#eb6834" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><circle cx="540.0" cy="203.2" r="4" fill="#eb6834" stroke="#ffffff" stroke-width="2"/><path d="M64.0,41.6 L79.9,52.4 L95.7,63.4 L111.6,74.4 L127.5,85.4 L143.3,96.5 L159.2,107.5 L175.1,118.6 L190.9,129.7 L206.8,140.7 L222.7,151.8 L238.5,161.0 L254.4,156.5 L270.3,134.9 L286.1,160.1 L302.0,153.1 L317.9,155.9 L333.7,127.0 L349.6,139.1 L365.5,152.9 L381.3,147.4 L397.2,151.5 L413.1,158.9 L428.9,161.4 L444.8,145.8 L460.7,135.3 L476.5,155.7 L492.4,140.8 L508.3,154.8 L524.1,154.9 L540.0,133.9" fill="none" stroke="#1baf7a" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><circle cx="540.0" cy="133.9" r="4" fill="#1baf7a" stroke="#ffffff" stroke-width="2"/><text x="550.0" y="324.2" class="viz-label">Legendre</text><text x="550.0" y="207.2" class="viz-label">x on [−1, 1]</text><text x="550.0" y="137.9" class="viz-label">t on [0, 1]</text><line x1="0" x2="0" y1="30" y2="340" class="viz-crosshair" visibility="hidden"/></svg>
<div class="viz-tooltip" hidden></div>
</div>
<figcaption>L2 error of the best approximation r<sub>n</sub> to f(x) = 1/(2 + x) on [−1, 1], on a log scale. All three lines compute the same polynomial; only the basis differs. Up to n = 10 the three lines lie on top of each other. Hover over the chart, or focus it and use the arrow keys, to read the values.</figcaption>
<details class="viz-table"><summary>Show the numbers</summary><table><thead><tr><th>n</th><th>Legendre basis</th><th>monomials on [−1, 1]</th><th>monomials on [0, 1] (Hilbert matrix)</th></tr></thead><tbody><tr><td>0</td><td>2.5e-1</td><td>2.5e-1</td><td>2.5e-1</td></tr><tr><td>1</td><td>7.0e-2</td><td>7.0e-2</td><td>7.0e-2</td></tr><tr><td>2</td><td>1.9e-2</td><td>1.9e-2</td><td>1.9e-2</td></tr><tr><td>3</td><td>5.1e-3</td><td>5.1e-3</td><td>5.1e-3</td></tr><tr><td>4</td><td>1.4e-3</td><td>1.4e-3</td><td>1.4e-3</td></tr><tr><td>5</td><td>3.7e-4</td><td>3.7e-4</td><td>3.7e-4</td></tr><tr><td>6</td><td>1.0e-4</td><td>1.0e-4</td><td>1.0e-4</td></tr><tr><td>7</td><td>2.7e-5</td><td>2.7e-5</td><td>2.7e-5</td></tr><tr><td>8</td><td>7.2e-6</td><td>7.2e-6</td><td>7.2e-6</td></tr><tr><td>9</td><td>1.9e-6</td><td>1.9e-6</td><td>1.9e-6</td></tr><tr><td>10</td><td>5.2e-7</td><td>5.2e-7</td><td>5.2e-7</td></tr><tr><td>11</td><td>1.4e-7</td><td>1.4e-7</td><td>1.7e-7</td></tr><tr><td>12</td><td>3.7e-8</td><td>3.7e-8</td><td>3.0e-7</td></tr><tr><td>13</td><td>1.0e-8</td><td>1.0e-8</td><td>3.9e-6</td></tr><tr><td>14</td><td>2.7e-9</td><td>2.7e-9</td><td>1.9e-7</td></tr><tr><td>15</td><td>7.2e-10</td><td>7.2e-10</td><td>4.4e-7</td></tr><tr><td>16</td><td>1.9e-10</td><td>1.9e-10</td><td>3.2e-7</td></tr><tr><td>17</td><td>5.1e-11</td><td>5.1e-11</td><td>9.8e-6</td></tr><tr><td>18</td><td>1.4e-11</td><td>1.4e-11</td><td>2.3e-6</td></tr><tr><td>19</td><td>3.7e-12</td><td>1.1e-11</td><td>4.6e-7</td></tr><tr><td>20</td><td>9.9e-13</td><td>3.2e-11</td><td>8.7e-7</td></tr><tr><td>21</td><td>2.7e-13</td><td>8.1e-11</td><td>5.3e-7</td></tr><tr><td>22</td><td>7.1e-14</td><td>2.2e-10</td><td>2.2e-7</td></tr><tr><td>23</td><td>1.9e-14</td><td>5.6e-10</td><td>1.6e-7</td></tr><tr><td>24</td><td>5.2e-15</td><td>8.4e-10</td><td>1.0e-6</td></tr><tr><td>25</td><td>1.7e-15</td><td>7.1e-10</td><td>3.7e-6</td></tr><tr><td>26</td><td>1.1e-15</td><td>3.8e-10</td><td>3.2e-7</td></tr><tr><td>27</td><td>1.1e-15</td><td>6.0e-10</td><td>1.9e-6</td></tr><tr><td>28</td><td>1.1e-15</td><td>1.5e-9</td><td>3.6e-7</td></tr><tr><td>29</td><td>1.1e-15</td><td>1.4e-9</td><td>3.6e-7</td></tr><tr><td>30</td><td>1.1e-15</td><td>1.2e-9</td><td>4.3e-6</td></tr></tbody></table></details>
<script>(() => {const fig = document.currentScript.closest('figure');const svg = fig.querySelector('svg');const line = svg.querySelector('.viz-crosshair');const tip = fig.querySelector('.viz-tooltip');const names = ["Legendre basis", "monomials on [−1, 1]", "monomials on [0, 1] (Hilbert matrix)"]; const colors = ["#2a78d6", "#eb6834", "#1baf7a"];const rows = [["2.5e-1", "2.5e-1", "2.5e-1"], ["7.0e-2", "7.0e-2", "7.0e-2"], ["1.9e-2", "1.9e-2", "1.9e-2"], ["5.1e-3", "5.1e-3", "5.1e-3"], ["1.4e-3", "1.4e-3", "1.4e-3"], ["3.7e-4", "3.7e-4", "3.7e-4"], ["1.0e-4", "1.0e-4", "1.0e-4"], ["2.7e-5", "2.7e-5", "2.7e-5"], ["7.2e-6", "7.2e-6", "7.2e-6"], ["1.9e-6", "1.9e-6", "1.9e-6"], ["5.2e-7", "5.2e-7", "5.2e-7"], ["1.4e-7", "1.4e-7", "1.7e-7"], ["3.7e-8", "3.7e-8", "3.0e-7"], ["1.0e-8", "1.0e-8", "3.9e-6"], ["2.7e-9", "2.7e-9", "1.9e-7"], ["7.2e-10", "7.2e-10", "4.4e-7"], ["1.9e-10", "1.9e-10", "3.2e-7"], ["5.1e-11", "5.1e-11", "9.8e-6"], ["1.4e-11", "1.4e-11", "2.3e-6"], ["3.7e-12", "1.1e-11", "4.6e-7"], ["9.9e-13", "3.2e-11", "8.7e-7"], ["2.7e-13", "8.1e-11", "5.3e-7"], ["7.1e-14", "2.2e-10", "2.2e-7"], ["1.9e-14", "5.6e-10", "1.6e-7"], ["5.2e-15", "8.4e-10", "1.0e-6"], ["1.7e-15", "7.1e-10", "3.7e-6"], ["1.1e-15", "3.8e-10", "3.2e-7"], ["1.1e-15", "6.0e-10", "1.9e-6"], ["1.1e-15", "1.5e-9", "3.6e-7"], ["1.1e-15", "1.4e-9", "3.6e-7"], ["1.1e-15", "1.2e-9", "4.3e-6"]]; const degrees = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30];const left = 64, right = 540, width = 680, last = degrees.length - 1;let current = -1;const show = (k) => { k = Math.max(0, Math.min(last, k)); current = k; const x = left + (right - left) * k / last; line.setAttribute('x1', x); line.setAttribute('x2', x); line.setAttribute('visibility', 'visible'); tip.replaceChildren(); const head = document.createElement('div'); head.className = 'viz-tooltip-head'; head.textContent = 'Degree ' + degrees[k]; tip.append(head); names.forEach((name, i) => {  const row = document.createElement('div'); row.className = 'viz-tooltip-row';  const key = document.createElement('span'); key.className = 'viz-swatch'; key.style.background = colors[i];  const value = document.createElement('strong'); value.textContent = rows[k][i];  const label = document.createElement('span'); label.textContent = name;  row.append(key, value, label); tip.append(row); }); const box = svg.getBoundingClientRect(), frame = fig.getBoundingClientRect(); const xs = box.left - frame.left + x * box.width / width; tip.hidden = false; const room = frame.width - tip.offsetWidth - 8; tip.style.left = Math.max(0, Math.min(room, xs + 12)) + 'px'; tip.style.top = (box.top - frame.top + 8) + 'px';};const hide = () => { line.setAttribute('visibility', 'hidden'); tip.hidden = true; current = -1; };svg.addEventListener('pointermove', (event) => { const box = svg.getBoundingClientRect(); const x = (event.clientX - box.left) * width / box.width; show(Math.round((x - left) / (right - left) * last));});svg.addEventListener('pointerleave', hide);svg.addEventListener('focus', () => show(current < 0 ? 0 : current));svg.addEventListener('blur', hide);svg.addEventListener('keydown', (event) => { if (event.key === 'ArrowRight') { show(current + 1); event.preventDefault(); } if (event.key === 'ArrowLeft') { show(current - 1); event.preventDefault(); } if (event.key === 'Escape') hide();});})();</script>
</figure>

**What I see.**

- Up to $n = 10$, the three lines lie on top of each other. All three methods give the same polynomial, as the theory says. For example, $n = 1$ gives $0.5493 - 0.2958x$ in every basis, the same line as in Section 4.
- **Legendre** keeps going down by a factor of about $3.7$ per degree. It reaches $1.1 \times 10^{-15}$ at $n = 26$ and stays there. That is round-off level.
- **Monomials on $[0, 1]$** (Hilbert) leave the others at $n = 11$. After that the error does not improve. It jumps around between $1.6 \times 10^{-7}$ and $9.8 \times 10^{-6}$. The condition number of the matrix is $1.2 \times 10^{15}$ at $n = 10$ and $3.8 \times 10^{16}$ at $n = 11$. Beyond that, double precision cannot even compute it reliably.
- **Monomials on $[-1, 1]$** last longer. Their best error is $1.1 \times 10^{-11}$ at $n = 19$. Then the error grows to about $10^{-9}$. Here the condition number is $1.2 \times 10^{14}$ at $n = 19$. On $[-1, 1]$, half of the matrix entries are zero (odd powers), and $x^8$ and $x^{10}$ are a bit more different than on $[0, 1]$. So the trouble comes later, but it still comes.

The condition number is only a worst-case bound. The real loss was smaller than "condition number × $10^{-16}$" in both cases. But it told me correctly which basis would fail first.

**Bessel's inequality in numbers.** $\|f\|_2^2 = \int_{-1}^1 (2+x)^{-2}\,dx = \tfrac23$ exactly. The sum $\sum_{j \le n} (f, \phi_j)^2$ grows toward it from below: $0.6035$ at $n = 0$, $0.66182$ at $n = 1$, $0.666641$ at $n = 3$, and $0.666666666666399$ at $n = 10$. It never passed $\tfrac23$.

I also tried to get the error from Pythagoras, $\|f - r_n\|_2 = \sqrt{\tfrac23 - \sum_j a_j^2}$. This matches the measured error up to about $n = 10$. Then it gets stuck near $1.8 \times 10^{-8}$. The reason is cancellation. When the error is $\varepsilon$, we subtract two numbers near $\tfrac23$ that differ by only $\varepsilon^2$. In double precision we lose everything below about $10^{-16}$ in $\varepsilon^2$, so we cannot see $\varepsilon$ below about $10^{-8}$. The formula is exact, but this way of using it is not.

**Coefficients.** The Legendre coefficients shrink by a fixed ratio, about $0.268$ in size per step, and their sign alternates (for example $|a_{10}/a_9| = 1.856 \times 10^{-6} / 6.915 \times 10^{-6} = 0.268$). This matches $2 - \sqrt3 \approx 0.268$. For a function with a pole at $x = -2$, the theory of analytic functions predicts this ratio. In the monomial basis, the coefficients move when the degree grows: the constant term is $0.5493$ for $n = 0, 1$ but $0.4963$ for $n = 2, 3$, and the $x$ term changes from $-0.2958$ to $-0.2447$ at $n = 3$. In the Legendre basis, a coefficient never changes once it is computed.

**Back to Taylor.** The cubic least squares polynomial has a largest error of $0.0148$, at $x = -1$. The cubic Taylor polynomial had $0.0625$. That is $4.2$ times larger. The least squares error also changes sign 4 times on $[-1, 1]$, so it is spread over the whole interval. The Taylor error is near zero in the middle and large at the ends.

Here is the core of the code, in Python:

```python
def legendre_all(n, x):
    """P_0(x), ..., P_n(x) by (k+1) P_{k+1} = (2k+1) x P_k - k P_{k-1}."""
    p = [1.0, x]
    for k in range(1, n):
        p.append(((2 * k + 1) * x * p[k] - k * p[k - 1]) / (k + 1))
    return p[: n + 1]

# orthonormal basis phi_j = sqrt((2j+1)/2) P_j; a_j = (f, phi_j) by Gauss quadrature
leg = [sum(w * fv * ph[j] for w, fv, ph in zip(gw, fx, phis)) for j in range(N_MAX + 1)]

# monomial basis: normal equations, (x^i, x^j) = 2/(i+j+1) if i+j is even
mat = [[(2.0 / (i + j + 1) if (i + j) % 2 == 0 else 0.0) for j in range(n + 1)] for i in range(n + 1)]
mono = solve(mat, rhs_all[: n + 1])
```

The Legendre part is one loop. The monomial part needs a full linear solve for every $n$.

The full script, with all the numbers and the chart, is in [this site's repository](https://github.com/doubleyands/doubleyands.github.io/tree/main/scripts/posts/numerical-analysis/week-3).

## Summary

| topic | main idea |
|---|---|
| Piecewise polynomials | More points means more pieces, not a higher degree. |
| Cubic spline | Cubic pieces with $s, s', s''$ continuous. $4n$ unknowns, $4n - 2$ conditions, plus 2 end conditions. |
| Spline system | Unknowns $M_i = s''(x_i)$. Tridiagonal, symmetric, diagonally dominant. $O(n)$ to solve. |
| Approximation | Do not pass through the data. Make $\lVert f - p \rVert$ small. Weierstrass: always possible for continuous $f$. |
| Normal equations | $\sum_j a_j (x^i, x^j) = (f, x^i)$. Dense and badly conditioned with monomials. |
| Orthogonal polynomials | Gram–Schmidt on $1, x, x^2, \dots$. Legendre, Chebyshev, Laguerre. Three-term recurrences. $\phi_n$ has $n$ roots inside $(a, b)$. |
| Orthonormal least squares | $a_j = (f, \phi_j)$. No system to solve, and old coefficients stay. |
| Bessel and Parseval | $\sum (f, \phi_j)^2 \le \lVert f \rVert^2$, with equality in the limit on a finite interval. |

The main thing I take from this week: the choice of basis does not change the answer in exact arithmetic, but it decides whether we can compute the answer. In my experiment the same best polynomial was correct to $10^{-15}$ in the Legendre basis and only to about $10^{-7}$ in the monomial basis on $[0, 1]$. Splines teach a similar lesson in another way: choose the unknowns well ($M_i$ instead of $4n$ raw coefficients), and a hard problem becomes a tridiagonal system.
