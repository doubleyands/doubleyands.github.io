---
title: "Numerical Analysis, Week 4: Numerical Integration and a First Look at ODEs"
description: "Trapezoidal and Simpson rules with their error formulas, composite rules, Gauss quadrature and its link to orthogonal polynomials, singular integrals, and the start of initial value problems. Tested on ln 2 and on Runge's function."
pubDate: 2026-10-24
project: numerical-analysis
tags: [numerical analysis, numerical integration, quadrature, gauss quadrature, ODE, notes]
---

This is my summary of the fourth week of a numerical analysis course. As before, I write the ideas in my own words and check some of them with a small experiment. The course follows K. E. Atkinson, *An Introduction to Numerical Analysis*.

The main topic is numerical integration: how to compute

$$
I(f) = \int_a^b f(x)\,dx
$$

when we cannot find an antiderivative. It comes right after interpolation and approximation for a good reason. A computer cannot integrate a general $f$. But it can integrate a polynomial exactly, because $\int_a^b x^k\,dx = \frac{b^{k+1} - a^{k+1}}{k+1}$. So the plan is: replace $f$ by a polynomial, and integrate the polynomial. At the end of the week we also start a new chapter, ordinary differential equations.

## 1. The general form

Almost every integration rule has the same shape:

$$
I_n(f) = \sum_{j=1}^{n} w_j\,f(x_j).
$$

The $x_j$ are the **nodes** and the $w_j$ are the **weights**. A rule is just a choice of nodes and weights. Different choices have different names. A good choice gives a small error with few function values. A bad choice can lose a lot of accuracy.

There are two main ways to make the choice:

1. **Newton–Cotes rules.** Fix equally spaced nodes. Interpolate $f$ at the nodes and integrate the interpolating polynomial exactly. The weights follow.
2. **Gauss rules.** Choose the nodes *and* the weights together, so that the rule is exact for as many polynomials as possible. The nodes do not have to be equally spaced.

In practice, the second way is the more common one in software.

## 2. The trapezoidal rule

### From a line

Take the line $p_1$ through $(a, f(a))$ and $(b, f(b))$:

$$
p_1(x) = \frac{(b - x) f(a) + (x - a) f(b)}{b - a}.
$$

Integrate it exactly. The area under a line is a trapezoid, so

$$
\int_a^b p_1(x)\,dx = \frac{b - a}{2}\bigl(f(a) + f(b)\bigr).
$$

This is the **trapezoidal rule**. In the general form, $x_1 = a$, $x_2 = b$, and $w_1 = w_2 = \frac{b-a}{2}$. Most people guess this formula by intuition. But it is not only a guess: it is the exact integral of the linear interpolant.

### A tool: the mean value theorem for integrals

We need one tool from [Week 1](/blog/2026-10-03-numerical-analysis-week-1/).

**Theorem.** Let $w \ge 0$ be integrable on $[a, b]$, and let $f$ be continuous. Then for some $\xi \in [a, b]$,

$$
\int_a^b w(x)\,f(x)\,dx = f(\xi) \int_a^b w(x)\,dx.
$$

The condition $w \ge 0$ matters. If $w$ changes sign, positive and negative parts can cancel. For example, on $[-1, 1]$ with $w(x) = f(x) = x$, the left side is $\frac23$. But $\int_{-1}^1 x\,dx = 0$, so the right side is $0$ for every $\xi$. This condition will cause a problem for Simpson's rule below.

### The error

From [Week 2](/blog/2026-10-10-numerical-analysis-week-2/), the interpolation error in divided difference form is

$$
f(x) - p_1(x) = (x - a)(x - b)\,f[a, b, x].
$$

Integrate both sides. Here $a$ and $b$ are fixed, and $x$ moves. On $(a, b)$ the factor $(x - a)(x - b)$ is negative everywhere. So $-(x-a)(x-b) \ge 0$ is a valid weight, and we can use the theorem:

$$
E(f) = \int_a^b f\,dx - \frac{b-a}{2}\bigl(f(a) + f(b)\bigr)
= -\int_a^b \bigl[-(x-a)(x-b)\bigr]\,f[a, b, x]\,dx
= -f[a, b, \xi] \int_a^b -(x-a)(x-b)\,dx.
$$

Two facts finish it. First, $\int_a^b -(x-a)(x-b)\,dx = \frac{(b-a)^3}{6}$. Second, from Week 2, $f[a, b, \xi] = \frac{f''(\eta)}{2}$ for some $\eta$. So

$$
E(f) = -\frac{(b - a)^3}{12}\,f''(\eta), \qquad \eta \in [a, b].
$$

The Lagrange form of the error, $\frac{(x-a)(x-b)}{2} f''(\xi_x)$, is the same error. But its $\xi_x$ depends on $x$ in an unknown way, so it is hard to integrate. The divided difference form is easier to analyze.

### The composite rule

The factor $(b-a)^3$ is bad when the interval is long. So we cut $[a, b]$ into $n$ pieces of length $h = (b - a)/n$, with $x_j = a + jh$, and use the rule on each piece. Each inner node belongs to two pieces and gets half a weight from each. So

$$
I_n(f) = h\left(\tfrac12 f_0 + f_1 + f_2 + \dots + f_{n-1} + \tfrac12 f_n\right), \qquad f_j = f(x_j).
$$

The pieces do not have to be equal, but equal pieces give this simple formula.

**Error.** Add up the errors of the pieces:

$$
E_n(f) = -\frac{h^3}{12} \sum_{j=1}^{n} f''(\eta_j) = -\frac{h^3 n}{12}\left[\frac{1}{n}\sum_{j=1}^{n} f''(\eta_j)\right].
$$

The bracket is an average of values of $f''$. It lies between $\min f''$ and $\max f''$. If $f''$ is continuous, the intermediate value theorem gives an $\eta$ with $f''(\eta)$ equal to this average. With $nh = b - a$,

$$
E_n(f) = -\frac{(b - a)\,h^2}{12}\,f''(\eta).
$$

Each piece has an error of size $h^3$, but there are $n = (b-a)/h$ pieces. The sum loses one power of $h$. So the composite rule is $O(h^2)$, not $O(h^3)$. When we halve $h$, the error drops by a factor of $4$.

### Asymptotic error estimate and the corrected rule

We can say more. Divide the error by $h^2$:

$$
\frac{E_n(f)}{h^2} = -\frac{1}{12} \sum_{j=1}^{n} f''(\eta_j)\,h.
$$

The sum on the right is a Riemann sum. As $n \to \infty$ it goes to $\int_a^b f''(x)\,dx = f'(b) - f'(a)$. So for large $n$,

$$
E_n(f) \approx -\frac{h^2}{12}\bigl(f'(b) - f'(a)\bigr) \equiv \tilde E_n(f).
$$

**Definition.** $\tilde E_n(f)$ is an **asymptotic error estimate** of $E_n(f)$ if

$$
\lim_{n\to\infty} \frac{\tilde E_n(f)}{E_n(f)} = 1.
$$

The estimate is not equal to the error for a fixed $n$. But the two behave the same as $n$ grows.

If we know the error so well, we can add it to the rule. This gives the **corrected trapezoidal rule**:

$$
CT_n(f) = I_n(f) + \tilde E_n(f) = h\left(\tfrac12 f_0 + f_1 + \dots + f_{n-1} + \tfrac12 f_n\right) - \frac{h^2}{12}\bigl(f'(b) - f'(a)\bigr).
$$

For smooth $f$, this is $O(h^4)$. The weak point: we must know $f'(a)$ and $f'(b)$. Usually we do not.

## 3. Simpson's rule

### From a parabola

Now use three points: $a$, the midpoint $c = \frac{a+b}{2}$, and $b$. Interpolate with a parabola $p_2$ and integrate it exactly. With $h = \frac{b-a}{2}$, the result is

$$
I_2(f) = \frac{h}{3}\Bigl(f(a) + 4 f(c) + f(b)\Bigr).
$$

This is **Simpson's rule**. The middle point gets weight $4$ and the ends get $1$. The weights add up to $2h = b - a$, which is right, because a rule must integrate $f \equiv 1$ exactly.

### The error needs a trick

As before,

$$
E_2(f) = \int_a^b (x - a)(x - c)(x - b)\,f[a, b, c, x]\,dx.
$$

But now we cannot use the mean value theorem directly. The cubic $(x-a)(x-c)(x-b)$ changes sign at $x = c$. It is positive on $(a, c)$ and negative on $(c, b)$, and because $c$ is the midpoint, the two areas are equal.

The trick is to integrate once. Define

$$
w(x) = \int_a^x (t - a)(t - c)(t - b)\,dt, \qquad \text{so } w'(x) = (x-a)(x-c)(x-b).
$$

Then $w(a) = 0$ (an empty integral), and $w(b) = 0$ (the two equal areas cancel). In between, $w > 0$. To see this, shift to $a = -h$, $c = 0$, $b = h$. Then $(t + h)\,t\,(t - h) = t^3 - h^2 t$, and

$$
w(x) = \frac{(x^2 - h^2)^2}{4}.
$$

It is a square, so it is never negative, and it is zero only at $x = \pm h$.

Now integrate by parts. The boundary term is zero because $w(a) = w(b) = 0$:

$$
E_2(f) = \Bigl[w(x)\,f[a,b,c,x]\Bigr]_a^b - \int_a^b w(x)\,\frac{d}{dx} f[a, b, c, x]\,dx = -\int_a^b w(x)\,f[a, b, c, x, x]\,dx.
$$

Why is the derivative $f[a, b, c, x, x]$? By definition,

$$
\frac{d}{dx} f[a,b,c,x] = \lim_{\delta\to 0} \frac{f[a,b,c,x+\delta] - f[a,b,c,x]}{\delta}.
$$

Divided differences do not depend on the order of their arguments (Week 2). So the fraction is exactly the next divided difference $f[a,b,c,x,x+\delta]$. Let $\delta \to 0$.

Now $w \ge 0$, so the mean value theorem works:

$$
E_2(f) = -f[a, b, c, \xi, \xi] \int_a^b w(x)\,dx = -\frac{f^{(4)}(\eta)}{24} \cdot \frac{4}{15} h^5 = -\frac{h^5}{90}\,f^{(4)}(\eta).
$$

Here $\int_{-h}^{h} \frac{(x^2 - h^2)^2}{4}\,dx = \frac{4}{15}h^5$, and $f[a,b,c,\xi,\xi] = \frac{f^{(4)}(\eta)}{4!}$.

### A gift: exact for cubics

The error has $f^{(4)}$. So if $f$ is a polynomial of degree $3$ or less, the error is **zero**. This is surprising. A parabola cannot match a cubic, so the pointwise error $f - p_2$ is not zero. But its integral is.

Here is a simple reason. Write the cubic around the midpoint. The part of degree $\le 2$ is matched exactly by $p_2$. What is left is a multiple of $(x - c)^3$. This is an odd function around $c$, so its integral over $[a, b]$ is $0$. Simpson's rule also gives $0$ for it: $(a - c)^3 + 4 \cdot 0 + (b - c)^3 = 0$. Symmetry gives one degree for free.

My script checks this on $[-1, 1]$. Simpson's rule gives the exact values $2, 0, \frac23, 0$ for $x^0, \dots, x^3$. For $x^4$ it gives $\frac23$, but the exact value is $\frac25$.

### Composite Simpson

Each Simpson piece needs a midpoint, so we use pairs of subintervals. With $n$ even, $h = (b-a)/n$, and $n/2$ pieces $[x_{2j-2}, x_{2j}]$:

$$
I_n(f) = \frac{h}{3}\Bigl(f_0 + 4f_1 + 2f_2 + 4f_3 + 2f_4 + \dots + 2f_{n-2} + 4f_{n-1} + f_n\Bigr).
$$

Odd nodes are midpoints (weight $4$). Even inner nodes are shared by two pieces (weight $1 + 1 = 2$). The error is

$$
E_n(f) = -\frac{h^5}{90} \sum_{j=1}^{n/2} f^{(4)}(\eta_j) = -\frac{(b - a)\,h^4}{180}\,f^{(4)}(\eta) \approx -\frac{h^4}{180}\bigl(f'''(b) - f'''(a)\bigr).
$$

Be careful here. The formula has $h^5$ inside the sum, but the order is $O(h^4)$, not $O(h^5)$. There are $n/2$ terms, and $n$ grows like $1/h$. The sum eats one power of $h$, just as for the trapezoidal rule.

Compared with the trapezoidal rule, Simpson is clearly better. Compared with the corrected trapezoidal rule, it has the same order, $O(h^4)$, but it needs no derivatives. (The last formula above is an asymptotic error estimate for Simpson. Adding it gives a corrected Simpson rule, but again it needs derivatives at the ends.)

## 4. Newton–Cotes rules

The trapezoidal rule ($p_1$) and Simpson's rule ($p_2$) are the first two members of the **Newton–Cotes** family. Take $n + 1$ equally spaced nodes $x_j = a + jh$, interpolate with $p_n$, and integrate $p_n$ exactly. It does not matter whether we write $p_n$ in Lagrange form or Newton form. The same data give the same polynomial, so the same rule. Higher members exist for every $n$.

They have two weak points:

- **The nodes must be equally spaced.** We saw in Week 2 that a high-degree interpolant on equal nodes can swing wildly near the ends (Runge's phenomenon). The integral of that interpolant inherits the problem. For large $n$, some weights even become negative. Then the rule can fail to converge, even for a nice function. My experiment shows this below.
- **The nodes must include the end points $a$ and $b$.** This is a problem if $f$ is not defined at an end point.

The usual fix for the first problem is a composite rule: keep the degree low and use many pieces. The fix for both is the Gauss rule.

## 5. Gauss quadrature

### The aim: zero error for polynomials

Start with $\int_{-1}^1 f(x)\,dx \approx \sum_{j=1}^n w_j f(x_j)$. Now both the nodes $x_j$ and the weights $w_j$ are unknowns. (Do not mix up the weights $w_j$ with a weight *function* $w(x)$. Even when $w(x) = 1$, the $w_j$ are not $1$ in general.)

The error $E_n(f) = \int_{-1}^1 f\,dx - \sum_j w_j f(x_j)$ is **linear** in $f$. So for a polynomial,

$$
E_n(a_0 + a_1 x + \dots + a_m x^m) = a_0 E_n(1) + a_1 E_n(x) + \dots + a_m E_n(x^m).
$$

The coefficients $a_i$ depend on $f$. But the functions $1, x, \dots, x^m$ are fixed. If

$$
E_n(x^i) = 0, \qquad i = 0, 1, \dots, m,
$$

then the error is zero for **every** polynomial of degree $\le m$. These equations decide the nodes and weights.

### One point

Two unknowns, $w_1$ and $x_1$, so two equations:

$$
E_1(1) = 2 - w_1 = 0, \qquad E_1(x) = 0 - w_1 x_1 = 0.
$$

So $w_1 = 2$ and $x_1 = 0$: the rule is $2 f(0)$, the midpoint rule. This is the natural guess: the length of the interval times the value in the middle.

### Two points

Four unknowns, so four equations, $i = 0, 1, 2, 3$:

$$
w_1 + w_2 = 2, \qquad w_1 x_1 + w_2 x_2 = 0, \qquad w_1 x_1^2 + w_2 x_2^2 = \tfrac23, \qquad w_1 x_1^3 + w_2 x_2^3 = 0.
$$

Solving this by elimination is long. A shorter way: the second equation gives $w_1 x_1 = -w_2 x_2$. Put this into the fourth:

$$
w_1 x_1 \cdot x_1^2 + w_2 x_2 \cdot x_2^2 = w_1 x_1 (x_1^2 - x_2^2) = 0.
$$

If $w_1 x_1 = 0$, then $w_2 x_2 = 0$ too, and the third equation would give $0 = \frac23$. So $x_1^2 = x_2^2$, and since the nodes are different, $x_2 = -x_1$. The second equation becomes $(w_1 - w_2)x_1 = 0$, so $w_1 = w_2 = 1$. The third gives $2x_1^2 = \frac23$. So

$$
\int_{-1}^1 f(x)\,dx \approx f\!\left(-\tfrac{1}{\sqrt3}\right) + f\!\left(\tfrac{1}{\sqrt3}\right).
$$

This explains why the nodes come out symmetric.

### Degree of precision

With $n$ points we have $2n$ unknowns, so we can ask $E_n(x^i) = 0$ for $i = 0, \dots, 2n - 1$. The rule is then exact for all polynomials of degree $\le 2n - 1$. This number is the **degree of precision**. The two-point rule is exact up to degree $3$, the same as Simpson with three points. A three-point Gauss rule is exact up to degree $5$.

My script checks the two-point rule on $x^0, \dots, x^5$. It is exact up to $x^3$ (and for $x^5$, which is odd). For $x^4$ it gives $\frac29$, but the exact value is $\frac25$.

### The nodes are roots of orthogonal polynomials

Solving larger systems by hand quickly becomes painful. Here [Week 3](/blog/2026-10-17-numerical-analysis-week-3/) helps. In the linearity argument above, we used the basis $1, x, x^2, \dots$. We could use an orthogonal basis $\phi_0, \phi_1, \dots$ instead. This leads to a clean answer:

- For a weight function $w(x) \ge 0$ on $(a, b)$, the $n$-point Gauss rule for $\int_a^b w(x) f(x)\,dx$ is unique, and its degree of precision is $2n - 1$.
- Its nodes are the $n$ roots of $\phi_n$, the orthogonal polynomial of degree $n$ for that weight.

In Week 3 we saw that $\phi_n$ has exactly $n$ distinct real roots, all inside $(a, b)$. Now we know what they are good for. There are formulas for the weights and for the error, but in practice we do not compute them by hand. Tables and libraries already have them for each family.

Two examples:

- **Gauss–Legendre** ($w = 1$ on $[-1, 1]$). The Legendre polynomial $P_2 = \frac12(3x^2 - 1)$ has the roots $\pm\frac{1}{\sqrt3}$. These are exactly the two nodes we found above. $P_3 = \frac12(5x^3 - 3x)$ has the roots $0$ and $\pm\sqrt{3/5}$, and the weights are $\frac89$ and $\frac59$.
- **Chebyshev–Gauss** ($w = 1/\sqrt{1 - x^2}$ on $[-1, 1]$). The nodes are the Chebyshev roots $x_j = \cos\frac{(2j - 1)\pi}{2n}$, and all weights are equal: $w_j = \frac{\pi}{n}$. For $n = 1$ the weight is $\pi$, not $2$, because the weight function is different.

### Nodes are inside the interval

A big difference from Newton–Cotes: **Gauss–Legendre nodes are never at the end points.** They are always strictly inside the interval. This matters in practice. If $f$ blows up at an end point, a Gauss rule never evaluates it there. (There are variants that put nodes at the ends on purpose, called Gauss–Lobatto and Gauss–Radau. They are exceptions.)

### Triangles and the affine map

In 2D and 3D, Gauss rules are given on a **reference element**, for example the triangle $\hat T$ with corners $(0,0)$, $(1,0)$, $(0,1)$. A real triangle $T$ has other corners $(x_1, y_1)$, $(x_2, y_2)$, $(x_3, y_3)$. We connect them with an **affine map**

$$
\psi_1(\xi, \eta) = a_1 \xi + b_1 \eta + c_1, \qquad \psi_2(\xi, \eta) = a_2 \xi + b_2 \eta + c_2.
$$

Asking the corners to match gives six equations for six unknowns:

$$
c_1 = x_1,\; a_1 = x_2 - x_1,\; b_1 = x_3 - x_1, \qquad c_2 = y_1,\; a_2 = y_2 - y_1,\; b_2 = y_3 - y_1.
$$

We move the reference nodes with $\psi$. We do not need a new table of weights. But we must multiply by the area ratio, as in any change of variables:

$$
\int_T f\,dx\,dy = |\det J| \int_{\hat T} f(\psi(\xi, \eta))\,d\xi\,d\eta \approx |\det J| \sum_q \hat w_q\,f(\psi(\xi_q, \eta_q)),
$$

where $\det J = a_1 b_2 - b_1 a_2$ is twice the area of $T$.

A small example. Let $T$ have corners $(0, 0)$, $(3, 1)$, $(1, 2)$, so $\det J = 3 \cdot 2 - 1 \cdot 1 = 5$. The simplest rule on $\hat T$ uses the centroid $(\frac13, \frac13)$ with weight $\frac12$ (the area of $\hat T$). The map sends the centroid to $(\frac43, 1)$, the centroid of $T$. For $f = x$, the rule gives $5 \cdot \frac12 \cdot \frac43 = \frac{10}{3}$. This is exact, because the rule is exact for linear functions. Without the factor $5$ we would get $\frac23$.

There is no single "best" set of nodes on a triangle. Many rules with the same number of points exist, with different nodes and weights.

## 6. Singular integrals

The rules above assume a smooth $f$. Two cases need extra care:

1. $f$ has a singularity in the interval, for example $\int_0^1 \log x\,dx$.
2. The interval is infinite.

Adaptive methods can help, but it is usually faster to study the singularity first and treat it directly. There is no single method for every singularity. Here are three tools.

**Tool 1: change of variables.** For smooth $g$,

$$
\int_0^b \frac{g(x)}{\sqrt{x}}\,dx \quad\xrightarrow{\;x = u^2,\; dx = 2u\,du\;}\quad \int_0^{\sqrt b} 2\,g(u^2)\,du.
$$

The $\sqrt{x} = u$ cancels the $u$ in $dx$. The new integrand is smooth, so any standard rule works well. For an infinite interval, $\int_1^\infty \frac{g(x)}{x^p}\,dx$ with $p > 1$, the change $x = u^{-\alpha}$ ($\alpha > 0$) maps $[1, \infty)$ to $(0, 1]$:

$$
\int_1^\infty \frac{g(x)}{x^p}\,dx = \alpha \int_0^1 u^{(p-1)\alpha - 1}\,g(u^{-\alpha})\,du.
$$

Whether the new integrand is nice depends on $g$ and on the choice of $\alpha$. It is an option, not a guarantee.

**Tool 2: Gauss rules with a weight function.** Put the difficult part into $w(x)$ and use the matching Gauss rule:

$$
\int_0^\infty x^\alpha e^{-x} g(x)\,dx, \qquad \int_{-\infty}^\infty e^{-x^2} g(x)\,dx, \qquad \int_0^1 \log\!\left(\tfrac1x\right) g(x)\,dx.
$$

These match the weight functions of the Laguerre (generalized) and Hermite families from Week 3, and a rule with a log weight. The rule never evaluates $g$ at infinity or at the singular point. Even a plain Gauss–Legendre rule never evaluates $f$ at an end point, so it gives a finite answer for $\int_0^1 \log x\,dx$. But it converges slowly, as the experiment shows.

**Tool 3: split and do the bad part by hand.** If $f$ is singular at $0$, write

$$
\int_0^b f\,dx = \int_0^\epsilon f\,dx + \int_\epsilon^b f\,dx.
$$

Use a standard rule on the second part. For the first part, use a series. For example, take $f(x) = e^x \log x$. Expand $e^x$ in its Taylor series and use

$$
\int_0^\epsilon x^n \log x\,dx = \frac{\epsilon^{n+1}}{n+1}\left(\log\epsilon - \frac{1}{n+1}\right),
$$

which follows from integration by parts. (At $x = 0$ the boundary term is $0 \cdot \log 0$, and its limit is $0$.) So

$$
\int_0^\epsilon e^x \log x\,dx = \sum_{n=0}^\infty \frac{\epsilon^{n+1}}{n!\,(n+1)}\left(\log\epsilon - \frac{1}{n+1}\right).
$$

With $\epsilon = 0.1$, the terms shrink fast. In my script, one term has an error of $1.4 \times 10^{-2}$, two terms $4.5 \times 10^{-4}$, and three terms $1.1 \times 10^{-5}$.

## 7. A first look at ODEs

The last part of the week starts a new chapter: numerical methods for ordinary differential equations. The goal is the **initial value problem**

$$
y' = f(x, y), \qquad y(x_0) = y_0,
$$

where $f$ is continuous on some region $D$. A function $u$ is a solution on $[a, b]$ if, for all $x$ in $[a, b]$:

1. $(x, u(x)) \in D$,
2. $u(x_0) = y_0$, and
3. $u'(x)$ exists and $u'(x) = f(x, u(x))$.

### A linear example

Take $y' = \lambda y + g(x)$ with $y(0) = y_0$, where $g$ is continuous. Multiply both sides by $e^{-\lambda x}$. Is this allowed? Multiplying by something can change the problem. For example, if we multiply $h(x) = 0$ by $(x - 2)$, we get a new solution $x = 2$. But $e^{-\lambda x}$ is never zero, so the new equation has exactly the same solutions.

By the product rule, the left side becomes a single derivative:

$$
e^{-\lambda x} y' - \lambda e^{-\lambda x} y = \bigl(e^{-\lambda x} y\bigr)' = e^{-\lambda x} g(x).
$$

Rename the variable to $t$ and integrate from $0$ to $x$. Then $x$ stays as the free variable:

$$
e^{-\lambda x} y(x) - y(0) = \int_0^x e^{-\lambda t} g(t)\,dt
\quad\Longrightarrow\quad
y(x) = y_0\,e^{\lambda x} + \int_0^x e^{\lambda (x - t)}\,g(t)\,dt.
$$

This formula is called **Duhamel's principle**. The meaning is more important than the steps. The solution is decided by two things: the initial value and the forcing. The first term is $y_0$ growing or decaying for time $x$ with no forcing. The second term adds up the forcing: the push $g(t)$ at time $t$ grows or decays by $e^{\lambda(x - t)}$ for the remaining time $x - t$. Change the initial value or the forcing, and the solution changes.

For example, $y' = -y + 1$ with $y(0) = 0$ gives $y = \int_0^x e^{-(x-t)}\,dt = 1 - e^{-x}$.

### A nonlinear example

Take $y' = -y^2$. One solution is $y \equiv 0$. Others are $y = \frac{1}{x + c}$ for any constant $c$. Without an initial condition there are infinitely many solutions. (This alone is not special to nonlinear equations: $y' = -y$ also has the solutions $Ce^{-x}$.) An initial condition picks one: $y(0) = y_0$ gives

$$
y(x) = \frac{y_0}{1 + y_0 x}.
$$

Here is a real nonlinear effect. With $y_0 = -2$, the solution is $y = \frac{-2}{1 - 2x}$, and it blows up at $x = \frac12$. How far the solution exists depends on the initial value. For nonlinear equations we often do not know in advance whether a solution exists, or for how long.

### Direction fields

Even without a formula, we can see how solutions behave. If $y$ is a solution, then $y'$ is its slope. So at each point $(x, y)$ we draw a short line with slope $f(x, y)$. Solutions follow these lines.

- $y' = -y^2$: the slope depends only on $y$. At $y = \pm 1$ the slope is $-1$, at $y = 0$ it is $0$, and at $y = 2$ it is $-4$, much steeper.
- $y' = -y$: the slope again depends only on $y$, so all marks at the same height are parallel. The solutions are $Ce^{-x}$.
- $y' = x^2 + y^2$: on the circle $x^2 + y^2 = r^2$, every slope is $r^2$. On the circle of radius $\frac12$ the slope is $\frac14$ (flat). On radius $1$ it is $1$. On radius $2$ it is $4$ (steep). The curves where the slope is constant are called **isoclines**.

This is a qualitative way to study an equation. It does not give numbers, but it shows the shape of the solutions.

## 8. Experiment: how fast do the rules converge?

The theory gives orders: $O(h^2)$ for the trapezoidal rule, $O(h^4)$ for the corrected trapezoidal and Simpson rules, and degree of precision $2n - 1$ for Gauss. I wanted to see these numbers, and to see Newton–Cotes fail.

**Setup.** The script uses only the Python standard library and runs in under a second.

- **Part A.** $\int_0^1 \frac{dx}{1 + x} = \ln 2$. Four rules, each with $N = 2^k + 1$ function values, $k = 1, \dots, 10$: composite trapezoidal, corrected trapezoidal (with the exact $f'$), composite Simpson, and an $N$-point Gauss–Legendre rule. The Gauss nodes come from Newton's method ([Week 1](/blog/2026-10-03-numerical-analysis-week-1/)) on the Legendre recurrence.
- **Part B.** $\int_{-1}^1 \frac{dx}{1 + 16x^2} = \frac{\arctan 4}{2}$, the Runge-type function from Week 2. Newton–Cotes with one polynomial through all $N$ equal points (weights computed exactly with fractions), against Gauss–Legendre and composite Simpson with the same $N$.
- **Part C.** $\int_0^1 \frac{e^x}{\sqrt{x}}\,dx$. Gauss–Legendre directly, and after the change $x = u^2$.

<figure class="viz">
<div class="viz-legend"><span class="viz-key"><span class="viz-swatch" style="background:#2a78d6"></span>composite trapezoidal</span><span class="viz-key"><span class="viz-swatch" style="background:#eb6834"></span>corrected trapezoidal</span><span class="viz-key"><span class="viz-swatch" style="background:#1baf7a"></span>composite Simpson</span><span class="viz-key"><span class="viz-swatch" style="background:#eda100"></span>Gauss–Legendre</span></div>
<div class="viz-frame">
<svg viewBox="0 0 680 400" role="img" tabindex="0" aria-label="Error of four quadrature rules for the integral of 1/(1+x) on [0, 1] against the number of points, both on log scales. The trapezoidal error falls by 4 each time the points double. The corrected trapezoidal and Simpson errors fall by 16. Gauss-Legendre reaches round-off level with 17 points." class="viz-plot"><line x1="64" x2="540" y1="30.0" y2="30.0" class="viz-grid"/><text x="54" y="34.0" text-anchor="end" class="viz-tick">1</text><line x1="64" x2="540" y1="107.5" y2="107.5" class="viz-grid"/><text x="54" y="111.5" text-anchor="end" class="viz-tick">10<tspan dy="-0.5em" font-size="75%">−4</tspan></text><line x1="64" x2="540" y1="185.0" y2="185.0" class="viz-grid"/><text x="54" y="189.0" text-anchor="end" class="viz-tick">10<tspan dy="-0.5em" font-size="75%">−8</tspan></text><line x1="64" x2="540" y1="262.5" y2="262.5" class="viz-grid"/><text x="54" y="266.5" text-anchor="end" class="viz-tick">10<tspan dy="-0.5em" font-size="75%">−12</tspan></text><line x1="64" x2="540" y1="340.0" y2="340.0" class="viz-grid"/><text x="54" y="344.0" text-anchor="end" class="viz-tick">10<tspan dy="-0.5em" font-size="75%">−16</tspan></text><text x="64.0" y="362" text-anchor="middle" class="viz-tick">3</text><text x="169.8" y="362" text-anchor="middle" class="viz-tick">9</text><text x="275.6" y="362" text-anchor="middle" class="viz-tick">33</text><text x="381.3" y="362" text-anchor="middle" class="viz-tick">129</text><text x="487.1" y="362" text-anchor="middle" class="viz-tick">513</text><text x="540.0" y="362" text-anchor="middle" class="viz-tick">1025</text><text x="302.0" y="388" text-anchor="middle" class="viz-axis-title">Number of function values N</text><text x="64" y="14" class="viz-axis-title">|error|</text><path d="M64.0,65.2 L116.9,76.7 L169.8,88.3 L222.7,100.0 L275.6,111.7 L328.4,123.3 L381.3,135.0 L434.2,146.6 L487.1,158.3 L540.0,170.0" fill="none" stroke="#2a78d6" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><circle cx="540.0" cy="170.0" r="4" fill="#2a78d6" stroke="#ffffff" stroke-width="2"/><path d="M64.0,95.1 L116.9,117.7 L169.8,140.9 L222.7,164.2 L275.6,187.5 L328.4,210.8 L381.3,234.1 L434.2,257.5 L487.1,280.8 L540.0,304.1" fill="none" stroke="#eb6834" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><circle cx="540.0" cy="304.1" r="4" fill="#eb6834" stroke="#ffffff" stroke-width="2"/><path d="M64.0,85.9 L116.9,106.9 L169.8,129.5 L222.7,152.6 L275.6,175.8 L328.4,199.1 L381.3,222.5 L434.2,245.8 L487.1,269.1 L540.0,292.5" fill="none" stroke="#1baf7a" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><circle cx="540.0" cy="292.5" r="4" fill="#1baf7a" stroke="#ffffff" stroke-width="2"/><path d="M64.0,119.0 L116.9,178.1 L169.8,296.6 L222.7,339.1" fill="none" stroke="#eda100" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><circle cx="222.7" cy="339.1" r="4" fill="#eda100" stroke="#ffffff" stroke-width="2"/><text x="550.0" y="174.0" class="viz-label">trapezoidal</text><text x="550.0" y="316.1" class="viz-label">corrected</text><text x="550.0" y="288.5" class="viz-label">Simpson</text><text x="232.7" y="343.1" class="viz-label">Gauss</text><line x1="0" x2="0" y1="30" y2="340" class="viz-crosshair" visibility="hidden"/></svg>
<div class="viz-tooltip" hidden></div>
</div>
<figcaption>Error of four rules for the integral of 1/(1 + x) over [0, 1] (exact value ln 2), against the number of function values N, both on log scales. Each step to the right doubles the number of subintervals. A line stops when it reaches round-off level (below 10<sup>−15</sup>). Hover over the chart, or focus it and use the arrow keys, to read the values.</figcaption>
<details class="viz-table"><summary>Show the numbers</summary><table><thead><tr><th>N</th><th>composite trapezoidal</th><th>corrected trapezoidal</th><th>composite Simpson</th><th>Gauss–Legendre</th></tr></thead><tbody><tr><td>3</td><td>1.5e-2</td><td>4.4e-4</td><td>1.3e-3</td><td>2.5e-5</td></tr><tr><td>5</td><td>3.9e-3</td><td>3.0e-5</td><td>1.1e-4</td><td>2.3e-8</td></tr><tr><td>9</td><td>9.7e-4</td><td>1.9e-6</td><td>7.4e-6</td><td>1.7e-14</td></tr><tr><td>17</td><td>2.4e-4</td><td>1.2e-7</td><td>4.7e-7</td><td>1.1e-16</td></tr><tr><td>33</td><td>6.1e-5</td><td>7.4e-9</td><td>3.0e-8</td><td>1.1e-16</td></tr><tr><td>65</td><td>1.5e-5</td><td>4.7e-10</td><td>1.9e-9</td><td>1.1e-16</td></tr><tr><td>129</td><td>3.8e-6</td><td>2.9e-11</td><td>1.2e-10</td><td>1.1e-16</td></tr><tr><td>257</td><td>9.5e-7</td><td>1.8e-12</td><td>7.3e-12</td><td>2.2e-16</td></tr><tr><td>513</td><td>2.4e-7</td><td>1.1e-13</td><td>4.5e-13</td><td>1.1e-16</td></tr><tr><td>1025</td><td>6.0e-8</td><td>7.1e-15</td><td>2.8e-14</td><td>1.1e-16</td></tr></tbody></table></details>
<script>(() => {const fig = document.currentScript.closest('figure');const svg = fig.querySelector('svg');const line = svg.querySelector('.viz-crosshair');const tip = fig.querySelector('.viz-tooltip');const names = ["composite trapezoidal", "corrected trapezoidal", "composite Simpson", "Gauss–Legendre"]; const colors = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"];const rows = [["1.5e-2", "4.4e-4", "1.3e-3", "2.5e-5"], ["3.9e-3", "3.0e-5", "1.1e-4", "2.3e-8"], ["9.7e-4", "1.9e-6", "7.4e-6", "1.7e-14"], ["2.4e-4", "1.2e-7", "4.7e-7", "1.1e-16"], ["6.1e-5", "7.4e-9", "3.0e-8", "1.1e-16"], ["1.5e-5", "4.7e-10", "1.9e-9", "1.1e-16"], ["3.8e-6", "2.9e-11", "1.2e-10", "1.1e-16"], ["9.5e-7", "1.8e-12", "7.3e-12", "2.2e-16"], ["2.4e-7", "1.1e-13", "4.5e-13", "1.1e-16"], ["6.0e-8", "7.1e-15", "2.8e-14", "1.1e-16"]]; const degrees = [3, 5, 9, 17, 33, 65, 129, 257, 513, 1025];const left = 64, right = 540, width = 680, last = degrees.length - 1;let current = -1;const show = (k) => { k = Math.max(0, Math.min(last, k)); current = k; const x = left + (right - left) * k / last; line.setAttribute('x1', x); line.setAttribute('x2', x); line.setAttribute('visibility', 'visible'); tip.replaceChildren(); const head = document.createElement('div'); head.className = 'viz-tooltip-head'; head.textContent = 'N = ' + degrees[k] + ' points'; tip.append(head); names.forEach((name, i) => {  const row = document.createElement('div'); row.className = 'viz-tooltip-row';  const key = document.createElement('span'); key.className = 'viz-swatch'; key.style.background = colors[i];  const value = document.createElement('strong'); value.textContent = rows[k][i];  const label = document.createElement('span'); label.textContent = name;  row.append(key, value, label); tip.append(row); }); const box = svg.getBoundingClientRect(), frame = fig.getBoundingClientRect(); const xs = box.left - frame.left + x * box.width / width; tip.hidden = false; const room = frame.width - tip.offsetWidth - 8; tip.style.left = Math.max(0, Math.min(room, xs + 12)) + 'px'; tip.style.top = (box.top - frame.top + 8) + 'px';};const hide = () => { line.setAttribute('visibility', 'hidden'); tip.hidden = true; current = -1; };svg.addEventListener('pointermove', (event) => { const box = svg.getBoundingClientRect(); const x = (event.clientX - box.left) * width / box.width; show(Math.round((x - left) / (right - left) * last));});svg.addEventListener('pointerleave', hide);svg.addEventListener('focus', () => show(current < 0 ? 0 : current));svg.addEventListener('blur', hide);svg.addEventListener('keydown', (event) => { if (event.key === 'ArrowRight') { show(current + 1); event.preventDefault(); } if (event.key === 'ArrowLeft') { show(current - 1); event.preventDefault(); } if (event.key === 'Escape') hide();});})();</script>
</figure>

**What I see in Part A.**

| $N$ | trapezoidal | ratio | corrected | ratio | Simpson | ratio | Gauss |
|---|---|---|---|---|---|---|---|
| 3 | $-1.52 \times 10^{-2}$ | | $4.39 \times 10^{-4}$ | | $-1.30 \times 10^{-3}$ | | $2.55 \times 10^{-5}$ |
| 5 | $-3.88 \times 10^{-3}$ | 3.92 | $2.96 \times 10^{-5}$ | 14.8 | $-1.07 \times 10^{-4}$ | 12.2 | $2.27 \times 10^{-8}$ |
| 9 | $-9.75 \times 10^{-4}$ | 3.98 | $1.89 \times 10^{-6}$ | 15.7 | $-7.35 \times 10^{-6}$ | 14.5 | $1.73 \times 10^{-14}$ |
| 33 | $-6.10 \times 10^{-5}$ | 4.00 | $7.45 \times 10^{-9}$ | 16.0 | $-2.97 \times 10^{-8}$ | 15.9 | round-off |
| 1025 | $-5.96 \times 10^{-8}$ | 4.00 | $7.11 \times 10^{-15}$ | 16.0 | $-2.83 \times 10^{-14}$ | 16.1 | round-off |

The ratio is "error in the row before / error in this row" (the rows in between are in the chart).

- **Trapezoidal**: the ratio goes to $4.00$. Doubling $N$ halves $h$, and $2^2 = 4$, so this is $O(h^2)$. The error is negative here because $f'' > 0$, as the formula $-\frac{(b-a)h^2}{12} f''(\eta)$ says.
- **Corrected trapezoidal** and **Simpson**: both ratios go to $16 = 2^4$. So both are $O(h^4)$, not $O(h^5)$. With 1025 points both are near $10^{-14}$, while the trapezoidal rule is still at $6 \times 10^{-8}$.
- **The asymptotic estimate works.** The ratio $\tilde E_n / E_n$ for the trapezoidal rule is $1.029$ at $N = 3$, $1.0019$ at $N = 9$, and $1.00001$ at $N = 129$. It goes to $1$, as the definition asks.
- **Gauss** is in another league. With $9$ points the error is $1.7 \times 10^{-14}$. With $17$ points it is at round-off level ($10^{-16}$). The other rules never reach this, even with $1025$ points.

**Part B: equal nodes break down.**

| $N$ | Newton–Cotes error | smallest NC weight | composite Simpson error | Gauss error |
|---|---|---|---|---|
| 5 | $9.35 \times 10^{-2}$ | $0.156$ | $4.33 \times 10^{-2}$ | $-1.09 \times 10^{-1}$ |
| 9 | $1.78 \times 10^{-1}$ | $-0.320$ | $1.98 \times 10^{-2}$ | $-1.43 \times 10^{-2}$ |
| 13 | $3.29 \times 10^{-1}$ | $-2.78$ | $4.59 \times 10^{-3}$ | $-1.96 \times 10^{-3}$ |
| 17 | $7.38 \times 10^{-1}$ | $-20.9$ | $9.75 \times 10^{-4}$ | $-2.71 \times 10^{-4}$ |
| 21 | $1.86$ | $-180$ | $2.04 \times 10^{-4}$ | $-3.75 \times 10^{-5}$ |

The exact value is $0.6629$. With $21$ equal points, one Newton–Cotes rule has an error of $1.86$: larger than the integral itself. From $9$ points on, some weights are negative, and they get huge. The rule adds and subtracts large numbers, and the swings of the interpolant (Week 2) go straight into the integral. The same $21$ points used as a composite Simpson rule give $2 \times 10^{-4}$. Gauss with $21$ points gives $4 \times 10^{-5}$. So the problem is not the number of points. The problem is one high-degree polynomial on equal nodes.

**Part C: a singular integrand.** The exact value is $2\sum_{k\ge0} \frac{1}{k!\,(2k+1)} = 2.92530\ldots$

| $N$ | Gauss on $e^x/\sqrt{x}$ | Gauss after $x = u^2$ |
|---|---|---|
| 2 | $3.4 \times 10^{-1}$ | $1.7 \times 10^{-2}$ |
| 4 | $1.9 \times 10^{-1}$ | $9.9 \times 10^{-6}$ |
| 8 | $1.0 \times 10^{-1}$ | $2.0 \times 10^{-13}$ |
| 16 | $5.3 \times 10^{-2}$ | round-off |
| 64 | $1.4 \times 10^{-2}$ | round-off |

Gauss–Legendre never touches $x = 0$, so it gives a finite answer. But with the singularity left in, the error only halves when $N$ doubles. That is first order, very slow. After the change of variables the integrand $2e^{u^2}$ is smooth, and $8$ points give $2 \times 10^{-13}$. Studying the singularity first pays off.

Here is the core of the code, in Python:

```python
def legendre(n, x):
    """P_n(x) and P_n'(x) by the three-term recurrence."""
    p0, p1 = 1.0, x
    for k in range(1, n):
        p0, p1 = p1, ((2 * k + 1) * x * p1 - k * p0) / (k + 1)
    return p1, n * (x * p1 - p0) / (x * x - 1)

def gauss_nodes(n):
    nodes, weights = [], []
    for i in range(1, n + 1):
        x = math.cos(math.pi * (i - 0.25) / (n + 0.5))  # first guess
        for _ in range(100):                             # Newton's method
            p, dp = legendre(n, x)
            x -= p / dp
        p, dp = legendre(n, x)
        nodes.append(x)
        weights.append(2 / ((1 - x * x) * dp * dp))
    return nodes, weights
```

The full script, with all the numbers and the chart, is in [this site's repository](https://github.com/doubleyands/doubleyands.github.io/tree/main/scripts/posts/numerical-analysis/week-4).

## Summary

| topic | main idea |
|---|---|
| General form | $I_n(f) = \sum_j w_j f(x_j)$. A rule is a choice of nodes and weights. |
| Trapezoidal rule | Integrate the line through the end points. Error $-\frac{(b-a)^3}{12} f''(\eta)$. Composite: $O(h^2)$. |
| Corrected trapezoidal | Add the asymptotic error $-\frac{h^2}{12}(f'(b) - f'(a))$. $O(h^4)$, but needs $f'$ at the ends. |
| Simpson's rule | Integrate the parabola. A trick with $w(x) \ge 0$ gives the error $-\frac{h^5}{90} f^{(4)}(\eta)$. Exact for cubics. Composite: $O(h^4)$, not $O(h^5)$. |
| Newton–Cotes | Equal nodes, end points included. High degree gives negative weights and no convergence. |
| Gauss quadrature | Choose nodes and weights so $E_n(x^i) = 0$. Degree of precision $2n - 1$. Nodes are roots of orthogonal polynomials, always inside the interval. |
| Singular integrals | Analyze first. Change variables, use a weighted Gauss rule, or split and use a series. |
| ODEs | Initial value problems. Duhamel: initial value plus forcing. Nonlinear solutions can blow up. Direction fields show the shape. |

The main thing I take from this week: integration rules are built from interpolation, but the best rule does not use the obvious nodes. Equal nodes are easy, but a high-degree rule on them can be worse than no rule at all. Letting the nodes move, as Gauss does, gives much more accuracy for the same number of function values. It also links this week to the orthogonal polynomials of Week 3.
