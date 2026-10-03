---
title: "Numerical Analysis, Week 1: Basic Tools and Root Finding"
description: "My notes on the first week: four theorems, norms, error, and four ways to solve f(x) = 0, tested on x = cos x."
pubDate: 2026-10-03
project: numerical-analysis
tags: [numerical analysis, root finding, notes]
---

This semester I am taking a course in numerical analysis. This post is my summary of the first week. I write the ideas in my own words and check them with a small experiment. The course follows K. E. Atkinson, *An Introduction to Numerical Analysis*.

The week had two parts. The first part collects tools we will use again and again to study error. The second part asks one question: how do we solve $f(x) = 0$? Many problems can be written in this form, so root finding is a natural place to start.

## 1. Four theorems we keep using

**Intermediate value theorem.** Let $f$ be continuous on $[a, b]$, and let $m = \inf f$ and $M = \sup f$ on $[a, b]$. Then for every $c \in [m, M]$ there is at least one $x_0 \in [a, b]$ with $f(x_0) = c$.

The theorem says such a point exists. It does not say how many there are, or where. The bisection method below is built on it.

Why inf and sup instead of min and max? A minimum must belong to the set, but an infimum does not have to. On $(0, 1)$, $f(x) = x$ has no minimum, yet its infimum is $0$. Numerical analysis works with limits all the time, so inf and sup are the safer words.

**Mean value theorem.** If $f$ is continuous on $[a, b]$ and differentiable on $(a, b)$, then for some $\xi \in (a, b)$,

$$
\frac{f(b) - f(a)}{b - a} = f'(\xi).
$$

In words: a difference quotient equals a derivative at some point in between. This lets us trade derivatives for differences, and back. It shows up in the secant method, in error estimates for Newton's method, and in the fixed-point theorem below.

**Integral mean value theorem.** If $w(x) \ge 0$ is integrable and $f$ is continuous on $[a, b]$, then for some $\xi \in [a, b]$,

$$
\int_a^b w(x) f(x)\,dx = f(\xi) \int_a^b w(x)\,dx.
$$

It has the same shape as the mean value theorem, with a weight inside the integral. The weight must not change sign.

**Taylor's theorem.** If $f \in C^{n+1}[a, b]$ and $x, x_0 \in [a, b]$, then

$$
f(x) = \sum_{k=0}^{n} \frac{(x - x_0)^k}{k!} f^{(k)}(x_0) + \frac{(x - x_0)^{n+1}}{(n+1)!} f^{(n+1)}(\xi)
$$

for some $\xi$ between $x_0$ and $x$. Most error analysis in this course seems to follow two steps: expand with Taylor's theorem, then ask how large the part we threw away is.

## 2. Norms: a ruler for "how close"

To approximate is to be close to the true answer. To say "close", we need a way to measure distance. A norm does this.

A norm $\|\cdot\|$ maps a vector or a function to a non-negative number, and it follows three rules:

1. $\|u\| \ge 0$, and $\|u\| = 0$ only when $u = 0$.
2. $\|\alpha u\| = |\alpha|\,\|u\|$ for every scalar $\alpha$.
3. $\|u + w\| \le \|u\| + \|w\|$ (the triangle inequality).

Common examples, for vectors and for functions on $[a, b]$:

$$
\|x\|_\infty = \max_i |x_i|, \qquad \|x\|_2 = \Big(\sum_i x_i^2\Big)^{1/2},
$$

$$
\|f\|_\infty = \max_{a \le x \le b} |f(x)|, \qquad \|f\|_{L^2} = \Big(\int_a^b |f(x)|^2\,dx\Big)^{1/2}.
$$

**The square root is not decoration.** Without it, $\int |f|^2\,dx$ breaks rule 2: multiplying $f$ by $\alpha$ multiplies it by $\alpha^2$, not $|\alpha|$. It also breaks rule 3. Take $u = w = f$: the left side is $\int (2f)^2\,dx = 4\int f^2\,dx$, but the right side is only $2\int f^2\,dx$.

**The ruler decides which approximation is better.** Here is an example. On $[0, 1]$, let $g_A$ differ from $f$ by $10^{-3}$ everywhere. Let $g_B$ equal $f$ everywhere except on a tiny interval of length $10^{-8}$, where it is off by $1$.

| | error of $g_A$ | error of $g_B$ | better |
|---|---|---|---|
| $L^\infty$ | $10^{-3}$ | $1$ | $g_A$ |
| $L^2$ | $10^{-3}$ | $\sqrt{10^{-8}} = 10^{-4}$ | $g_B$ |

The two norms give opposite answers. If the bad interval shrinks to a single point, the $L^2$ error of $g_B$ becomes exactly $0$: a single point has measure zero, so it adds nothing to an integral. So choosing the norm is part of the problem. $L^\infty$ looks at the worst point. $L^2$ looks at the average and forgives small bad spots.

## 3. Where error comes from

Let $x$ be the true value and $x_A$ an approximation. Then

$$
\text{error} = x - x_A, \qquad \text{relative error} = \frac{x - x_A}{x}.
$$

A small error is not always a good approximation. If $x = 0.003$ and $x_A = 0.004$, the error is only $0.001$, but the relative error is about $33\%$. Always look at the error together with the scale.

Error comes in through several doors:

- **Modeling.** Simplifications we make when we turn a physical problem into equations.
- **Data.** Measured values are never exact.
- **Blunders.** Mistakes in arithmetic or in code.
- **Machine arithmetic.** A computer keeps a finite number of digits, so it must round or chop. $1/3$ becomes $0.333\ldots3$. Overflow and underflow belong here too.
- **Truncation.** An infinite process, like a series, is cut after finitely many terms.

## 4. How fast does an iteration converge?

Most root-finding methods are iterative. We start from a guess $x_0$ and produce $x_1, x_2, \ldots$ that should approach a root $\alpha$.

A sequence converges to $\alpha$ with **order** $p \ge 1$ if, for some constant $c > 0$,

$$
|\alpha - x_{n+1}| \le c\,|\alpha - x_n|^p, \qquad n \ge 0.
$$

When $p = 1$ (**linear convergence**) we also need $c < 1$. Without it, the inequality says nothing. An error that stays at $3$ forever satisfies $3 \le 1 \cdot 3$, and so does a sequence that jumps between $1$ and $-1$ around a root at $0$. With $c < 1$, every step really shrinks the error. A map that shrinks distances like this is called a contraction, and the idea comes back in the fixed-point theorem.

The gap between $p = 1$ and $p = 2$ is large. Linear convergence multiplies the error by the same factor each step. Quadratic convergence squares it, so the number of correct digits roughly doubles each step. The experiment in section 7 shows this.

## 5. Four methods

### Bisection

Suppose $f$ is continuous on $[a, b]$ and $f(a)\,f(b) < 0$. The signs at the two ends differ, so by the intermediate value theorem there is a root in between.

1. Let $c = (a + b)/2$.
2. If $b - c \le \varepsilon$, accept $c$ as the root and stop.
3. If $f(b)\,f(c) < 0$, set $a = c$. Otherwise set $b = c$. Go back to step 1.

Each step keeps the half that must contain a root. After $n$ halvings the interval has length $(b - a)/2^n$, so its midpoint $c_n$ satisfies

$$
|\alpha - c_n| \le \frac{b - a}{2^{n+1}}.
$$

The bound shrinks by $1/2$ each step: linear convergence with rate $1/2$.

- **Good:** it always converges, and at every step we know an interval that contains a root.
- **Bad:** it is slow. It also works only for a single equation, because there is no "sign change" for a vector.

### Newton's method

Replace $f$ by its tangent line at $x_n$, and take the root of the tangent as the next guess:

$$
x_{n+1} = x_n - \frac{f(x_n)}{f'(x_n)}.
$$

Taylor's theorem gives the same formula and also its error. Expand $f$ around $x_n$ and put in $x = \alpha$, where $f(\alpha) = 0$:

$$
0 = f(x_n) + (\alpha - x_n) f'(x_n) + \frac{(\alpha - x_n)^2}{2} f''(\xi_n).
$$

Divide by $f'(x_n)$ and rearrange:

$$
\alpha - x_{n+1} = -(\alpha - x_n)^2 \, \frac{f''(\xi_n)}{2 f'(x_n)}.
$$

Suppose $f'(\alpha) \ne 0$, and on an interval around $\alpha$ we have $|f''| \le M$ and $|f'| \ge m > 0$. Then

$$
|\alpha - x_{n+1}| \le \frac{M}{2m}\,|\alpha - x_n|^2,
$$

which is quadratic convergence. The constant $M/(2m)$ does not depend on $n$, which is what the definition of order needs.

There is a catch. We dropped the second-order term, and that is only fair when $x_n$ is already close to $\alpha$. From a bad start, the tangent can throw the next guess far away. Newton's method also needs $f'$. For a system of equations it needs the Jacobian and a linear solve at every step, which can be expensive.

### Secant method

If $f'$ is not available, replace it by the slope of the line through the last two points:

$$
x_{n+1} = x_n - f(x_n)\,\frac{x_n - x_{n-1}}{f(x_n) - f(x_{n-1})}.
$$

No derivative is needed. Each step costs only one new function value, because $f(x_{n-1})$ is already known from the step before. Near a simple root its order is $(1 + \sqrt{5})/2 \approx 1.62$: faster than linear, slower than Newton.

### Fixed-point iteration

A **fixed point** of $g$ is a solution of $x = g(x)$. The iteration is simply

$$
x_{n+1} = g(x_n).
$$

Root finding fits this form, since $f(x) = 0$ is the same as $x = x - f(x)$. But one equation can be written as $x = g(x)$ in many ways, and different choices of $g$ can behave very differently. So the real question is when the iteration converges.

**Theorem.** Let $g$ be continuous on $[a, b]$, and suppose

1. $g([a, b]) \subset [a, b]$, and
2. $|g(x) - g(y)| \le \lambda\,|x - y|$ for all $x, y \in [a, b]$, with $0 < \lambda < 1$.

Then $x = g(x)$ has exactly one solution $\alpha$ in $[a, b]$. The iteration converges to $\alpha$ from **any** $x_0 \in [a, b]$, and

$$
|\alpha - x_n| \le \frac{\lambda^n}{1 - \lambda}\,|x_1 - x_0|.
$$

Condition 1 keeps every iterate inside the interval where $g$ is defined. For example, $g(x) = x^2$ on $[1, 2]$ sends $2$ to $4$, and then the next step has nothing to work with. Condition 2 says $g$ is Lipschitz with a constant below $1$, that is, a contraction. Lipschitz is a stronger condition than continuity: $1/x$ is continuous on $(0, 1)$, but $|1/x - 1/y| = |x - y|/(xy)$ and $1/(xy)$ has no upper bound there, so no Lipschitz constant works.

The error bound takes three lines. Because $\alpha = g(\alpha)$ and $x_1 = g(x_0)$,

$$
|\alpha - x_0| \le |\alpha - x_1| + |x_1 - x_0| \le \lambda|\alpha - x_0| + |x_1 - x_0|,
$$

so $|\alpha - x_0| \le |x_1 - x_0|/(1 - \lambda)$. In the same way,

$$
|\alpha - x_n| = |g(\alpha) - g(x_{n-1})| \le \lambda\,|\alpha - x_{n-1}| \le \cdots \le \lambda^n\,|\alpha - x_0|.
$$

The final bound uses only $x_0$ and $x_1$, which we know. So before we even start, we can tell how many steps are enough.

If $g$ is differentiable, the mean value theorem lets us take $\lambda = \max_{[a, b]} |g'(x)|$. It also gives

$$
\lim_{n \to \infty} \frac{\alpha - x_{n+1}}{\alpha - x_n} = g'(\alpha),
$$

so the convergence is linear with rate $|g'(\alpha)|$. Linear convergence can be sped up with extrapolation, which uses the last few iterates to jump ahead. Aitken's method is one such technique.

## 6. Newton's method is a fixed-point iteration

Newton's update has the form $x_{n+1} = g_N(x_n)$ with

$$
g_N(x) = x - \frac{f(x)}{f'(x)}.
$$

So everything about fixed points applies to it. Differentiate:

$$
g_N'(x) = 1 - \frac{f'(x)^2 - f(x) f''(x)}{f'(x)^2} = \frac{f(x)\,f''(x)}{f'(x)^2}.
$$

At a simple root, $f(\alpha) = 0$ and $f'(\alpha) \ne 0$, so

$$
g_N'(\alpha) = 0.
$$

This one fact explains both the strength and the weakness of Newton's method.

- **Why it is fast.** A fixed-point iteration shrinks the error by about $|g'(\alpha)|$ per step. For Newton this factor is $0$, so the linear term disappears. Expanding $g_N$ around $\alpha$ leaves only the second-order term: $\alpha - x_{n+1} = -\tfrac{1}{2} g_N''(\xi_n)(\alpha - x_n)^2$. That is quadratic convergence again, seen from a different side.
- **Why it needs a good start.** $g_N'$ is continuous and equals $0$ at $\alpha$, so there is a small interval $[\alpha - \delta, \alpha + \delta]$ where $|g_N'| \le \lambda < 1$. On that interval, $|g_N(x) - \alpha| \le \lambda\,|x - \alpha| \le \delta$, so both conditions of the theorem hold, and Newton converges from any start inside it. Far from $\alpha$, $|g_N'|$ can be larger than $1$, and the theorem promises nothing.

Compare this with the simplest choice for $x = \cos x$, which is $g(x) = \cos x$. There $|g'(\alpha)| = \sin\alpha \approx 0.674$, so each step removes only about a third of the error.

This also explains a common strategy. Use a slow but safe method, like bisection, to get close to the root. Then switch to Newton's method to finish fast.

## 7. Experiment: solving x = cos x four ways

I solved $x = \cos x$, whose root is $\alpha \approx 0.7390851332$, with all four methods. For bisection, secant, and Newton I used $f(x) = x - \cos x$ and $f'(x) = 1 + \sin x$. For fixed-point iteration I used $g(x) = \cos x$.

- Bisection on $[0, 1]$.
- Fixed-point iteration from $x_0 = 1$.
- Secant method from $x_0 = 0,\ x_1 = 1$.
- Newton's method from $x_0 = 1$.

<figure class="viz">
<div class="viz-legend"><span class="viz-key"><span class="viz-swatch" style="background:#2a78d6"></span>Bisection</span><span class="viz-key"><span class="viz-swatch" style="background:#eb6834"></span>Fixed point</span><span class="viz-key"><span class="viz-swatch" style="background:#1baf7a"></span>Secant</span><span class="viz-key"><span class="viz-swatch" style="background:#eda100"></span>Newton</span></div>
<div class="viz-frame">
<svg viewBox="0 0 680 400" role="img" tabindex="0" aria-label="Error per step for four methods on a log scale. Newton and the secant method reach round-off within 7 steps; bisection and fixed-point iteration fall at a steady rate." class="viz-plot"><line x1="64" x2="600" y1="30.0" y2="30.0" class="viz-grid"/><text x="54" y="34.0" text-anchor="end" class="viz-tick">1</text><line x1="64" x2="600" y1="107.5" y2="107.5" class="viz-grid"/><text x="54" y="111.5" text-anchor="end" class="viz-tick">10<tspan dy="-0.5em" font-size="75%">−4</tspan></text><line x1="64" x2="600" y1="185.0" y2="185.0" class="viz-grid"/><text x="54" y="189.0" text-anchor="end" class="viz-tick">10<tspan dy="-0.5em" font-size="75%">−8</tspan></text><line x1="64" x2="600" y1="262.5" y2="262.5" class="viz-grid"/><text x="54" y="266.5" text-anchor="end" class="viz-tick">10<tspan dy="-0.5em" font-size="75%">−12</tspan></text><line x1="64" x2="600" y1="340.0" y2="340.0" class="viz-grid"/><text x="54" y="344.0" text-anchor="end" class="viz-tick">10<tspan dy="-0.5em" font-size="75%">−16</tspan></text><text x="64.0" y="362" text-anchor="middle" class="viz-tick">0</text><text x="198.0" y="362" text-anchor="middle" class="viz-tick">10</text><text x="332.0" y="362" text-anchor="middle" class="viz-tick">20</text><text x="466.0" y="362" text-anchor="middle" class="viz-tick">30</text><text x="600.0" y="362" text-anchor="middle" class="viz-tick">40</text><text x="332.0" y="388" text-anchor="middle" class="viz-axis-title">Step n</text><text x="64" y="14" class="viz-axis-title">Error</text><path d="M64.0,42.0 L77.4,68.0 L90.8,48.3 L104.2,54.9 L117.6,62.8 L131.0,75.1 L144.4,78.6 L157.8,90.0 L171.2,87.0 L184.6,102.9 L198.0,97.8 L211.4,110.3 L224.8,113.2 L238.2,126.5 L251.6,121.0 L265.0,133.0 L278.4,137.6 L291.8,146.1 L305.2,147.2 L318.6,169.3 L332.0,153.7 L345.4,161.0 L358.8,170.7 L372.2,191.5 L385.6,177.2 L399.0,184.8 L412.4,195.6 L425.8,205.2 L439.2,204.7 L452.6,234.6 L466.0,211.3 L479.4,217.7 L492.8,224.8 L506.2,233.7 L519.6,259.6 L533.0,240.8 L546.4,247.6 L559.8,255.8 L573.2,270.1 L586.6,268.2 L600.0,287.4" fill="none" stroke="#2a78d6" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><circle cx="600.0" cy="287.4" r="4" fill="#2a78d6" stroke="#ffffff" stroke-width="2"/><path d="M64.0,41.3 L77.4,43.6 L90.8,47.9 L104.2,50.8 L117.6,54.5 L131.0,57.6 L144.4,61.1 L157.8,64.3 L171.2,67.7 L184.6,71.0 L198.0,74.3 L211.4,77.6 L224.8,81.0 L238.2,84.3 L251.6,87.6 L265.0,90.9 L278.4,94.3 L291.8,97.6 L305.2,100.9 L318.6,104.2 L332.0,107.6 L345.4,110.9 L358.8,114.2 L372.2,117.5 L385.6,120.9 L399.0,124.2 L412.4,127.5 L425.8,130.8 L439.2,134.2 L452.6,137.5 L466.0,140.8 L479.4,144.1 L492.8,147.5 L506.2,150.8 L519.6,154.1 L533.0,157.4 L546.4,160.8 L559.8,164.1 L573.2,167.4 L586.6,170.7 L600.0,174.1" fill="none" stroke="#eb6834" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><circle cx="600.0" cy="174.1" r="4" fill="#eb6834" stroke="#ffffff" stroke-width="2"/><path d="M64.0,32.5 L77.4,41.3 L90.8,54.6 L104.2,79.5 L117.6,116.5 L131.0,178.7 L144.4,278.0 L157.8,340.0" fill="none" stroke="#1baf7a" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><circle cx="157.8" cy="340.0" r="4" fill="#1baf7a" stroke="#ffffff" stroke-width="2"/><path d="M64.0,41.3 L77.4,67.7 L90.8,118.3 L104.2,219.3 L117.6,340.0" fill="none" stroke="#eda100" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><circle cx="117.6" cy="340.0" r="4" fill="#eda100" stroke="#ffffff" stroke-width="2"/><text x="610" y="291.4" class="viz-label">Bisection</text><text x="610" y="178.1" class="viz-label">Fixed point</text><text x="86.8" y="150.2" text-anchor="end" class="viz-label">Newton</text><line x1="142.0" y1="238.3" x2="220.8" y2="239.1" class="viz-leader"/><text x="224.8" y="243.1" class="viz-label">Secant</text><line x1="0" x2="0" y1="30" y2="340" class="viz-crosshair" visibility="hidden"/></svg>
<div class="viz-tooltip" hidden></div>
</div>
<figcaption>Error |x<sub>n</sub> − α| after each step, on a log scale. A line that drops to the bottom edge has gone below 10<sup>−15</sup>, the level of round-off. Hover over the chart, or focus it and use the arrow keys, to read the values.</figcaption>
<details class="viz-table"><summary>Show the numbers</summary><table><thead><tr><th>n</th><th>Bisection</th><th>Fixed point</th><th>Secant</th><th>Newton</th></tr></thead><tbody><tr><td>0</td><td>2.4e-1</td><td>2.6e-1</td><td>7.4e-1</td><td>2.6e-1</td></tr><tr><td>1</td><td>1.1e-2</td><td>2.0e-1</td><td>2.6e-1</td><td>1.1e-2</td></tr><tr><td>2</td><td>1.1e-1</td><td>1.2e-1</td><td>5.4e-2</td><td>2.8e-5</td></tr><tr><td>3</td><td>5.2e-2</td><td>8.5e-2</td><td>2.8e-3</td><td>1.7e-10</td></tr><tr><td>4</td><td>2.0e-2</td><td>5.4e-2</td><td>3.4e-5</td><td>< 1e-15</td></tr><tr><td>5</td><td>4.7e-3</td><td>3.8e-2</td><td>2.1e-8</td><td>< 1e-15</td></tr><tr><td>6</td><td>3.1e-3</td><td>2.5e-2</td><td>1.6e-13</td><td>< 1e-15</td></tr><tr><td>7</td><td>8.0e-4</td><td>1.7e-2</td><td>< 1e-15</td><td>< 1e-15</td></tr><tr><td>8</td><td>1.1e-3</td><td>1.1e-2</td><td>< 1e-15</td><td>< 1e-15</td></tr><tr><td>9</td><td>1.7e-4</td><td>7.7e-3</td><td>< 1e-15</td><td>< 1e-15</td></tr><tr><td>10</td><td>3.2e-4</td><td>5.2e-3</td><td>< 1e-15</td><td>< 1e-15</td></tr><tr><td>11</td><td>7.1e-5</td><td>3.5e-3</td><td>< 1e-15</td><td>< 1e-15</td></tr><tr><td>12</td><td>5.1e-5</td><td>2.3e-3</td><td>< 1e-15</td><td>< 1e-15</td></tr><tr><td>13</td><td>1.0e-5</td><td>1.6e-3</td><td>< 1e-15</td><td>< 1e-15</td></tr><tr><td>14</td><td>2.0e-5</td><td>1.1e-3</td><td>< 1e-15</td><td>< 1e-15</td></tr><tr><td>15</td><td>4.8e-6</td><td>7.2e-4</td><td>< 1e-15</td><td>< 1e-15</td></tr><tr><td>16</td><td>2.8e-6</td><td>4.8e-4</td><td>< 1e-15</td><td>< 1e-15</td></tr><tr><td>17</td><td>1.0e-6</td><td>3.2e-4</td><td>< 1e-15</td><td>< 1e-15</td></tr><tr><td>18</td><td>8.9e-7</td><td>2.2e-4</td><td>< 1e-15</td><td>< 1e-15</td></tr><tr><td>19</td><td>6.4e-8</td><td>1.5e-4</td><td>< 1e-15</td><td>< 1e-15</td></tr><tr><td>20</td><td>4.1e-7</td><td>9.9e-5</td><td>< 1e-15</td><td>< 1e-15</td></tr><tr><td>21</td><td>1.7e-7</td><td>6.7e-5</td><td>< 1e-15</td><td>< 1e-15</td></tr><tr><td>22</td><td>5.5e-8</td><td>4.5e-5</td><td>< 1e-15</td><td>< 1e-15</td></tr><tr><td>23</td><td>4.6e-9</td><td>3.0e-5</td><td>< 1e-15</td><td>< 1e-15</td></tr><tr><td>24</td><td>2.5e-8</td><td>2.0e-5</td><td>< 1e-15</td><td>< 1e-15</td></tr><tr><td>25</td><td>1.0e-8</td><td>1.4e-5</td><td>< 1e-15</td><td>< 1e-15</td></tr><tr><td>26</td><td>2.8e-9</td><td>9.3e-6</td><td>< 1e-15</td><td>< 1e-15</td></tr><tr><td>27</td><td>9.0e-10</td><td>6.2e-6</td><td>< 1e-15</td><td>< 1e-15</td></tr><tr><td>28</td><td>9.6e-10</td><td>4.2e-6</td><td>< 1e-15</td><td>< 1e-15</td></tr><tr><td>29</td><td>2.8e-11</td><td>2.8e-6</td><td>< 1e-15</td><td>< 1e-15</td></tr><tr><td>30</td><td>4.4e-10</td><td>1.9e-6</td><td>< 1e-15</td><td>< 1e-15</td></tr><tr><td>31</td><td>2.1e-10</td><td>1.3e-6</td><td>< 1e-15</td><td>< 1e-15</td></tr><tr><td>32</td><td>8.9e-11</td><td>8.7e-7</td><td>< 1e-15</td><td>< 1e-15</td></tr><tr><td>33</td><td>3.1e-11</td><td>5.8e-7</td><td>< 1e-15</td><td>< 1e-15</td></tr><tr><td>34</td><td>1.4e-12</td><td>3.9e-7</td><td>< 1e-15</td><td>< 1e-15</td></tr><tr><td>35</td><td>1.3e-11</td><td>2.6e-7</td><td>< 1e-15</td><td>< 1e-15</td></tr><tr><td>36</td><td>5.9e-12</td><td>1.8e-7</td><td>< 1e-15</td><td>< 1e-15</td></tr><tr><td>37</td><td>2.2e-12</td><td>1.2e-7</td><td>< 1e-15</td><td>< 1e-15</td></tr><tr><td>38</td><td>4.0e-13</td><td>8.1e-8</td><td>< 1e-15</td><td>< 1e-15</td></tr><tr><td>39</td><td>5.1e-13</td><td>5.5e-8</td><td>< 1e-15</td><td>< 1e-15</td></tr><tr><td>40</td><td>5.2e-14</td><td>3.7e-8</td><td>< 1e-15</td><td>< 1e-15</td></tr></tbody></table></details>
<script>(() => {const fig = document.currentScript.closest('figure');const svg = fig.querySelector('svg');const line = svg.querySelector('.viz-crosshair');const tip = fig.querySelector('.viz-tooltip');const names = ["Bisection", "Fixed point", "Secant", "Newton"]; const colors = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"];const rows = [["2.4e-1", "2.6e-1", "7.4e-1", "2.6e-1"], ["1.1e-2", "2.0e-1", "2.6e-1", "1.1e-2"], ["1.1e-1", "1.2e-1", "5.4e-2", "2.8e-5"], ["5.2e-2", "8.5e-2", "2.8e-3", "1.7e-10"], ["2.0e-2", "5.4e-2", "3.4e-5", "< 1e-15"], ["4.7e-3", "3.8e-2", "2.1e-8", "< 1e-15"], ["3.1e-3", "2.5e-2", "1.6e-13", "< 1e-15"], ["8.0e-4", "1.7e-2", "< 1e-15", "< 1e-15"], ["1.1e-3", "1.1e-2", "< 1e-15", "< 1e-15"], ["1.7e-4", "7.7e-3", "< 1e-15", "< 1e-15"], ["3.2e-4", "5.2e-3", "< 1e-15", "< 1e-15"], ["7.1e-5", "3.5e-3", "< 1e-15", "< 1e-15"], ["5.1e-5", "2.3e-3", "< 1e-15", "< 1e-15"], ["1.0e-5", "1.6e-3", "< 1e-15", "< 1e-15"], ["2.0e-5", "1.1e-3", "< 1e-15", "< 1e-15"], ["4.8e-6", "7.2e-4", "< 1e-15", "< 1e-15"], ["2.8e-6", "4.8e-4", "< 1e-15", "< 1e-15"], ["1.0e-6", "3.2e-4", "< 1e-15", "< 1e-15"], ["8.9e-7", "2.2e-4", "< 1e-15", "< 1e-15"], ["6.4e-8", "1.5e-4", "< 1e-15", "< 1e-15"], ["4.1e-7", "9.9e-5", "< 1e-15", "< 1e-15"], ["1.7e-7", "6.7e-5", "< 1e-15", "< 1e-15"], ["5.5e-8", "4.5e-5", "< 1e-15", "< 1e-15"], ["4.6e-9", "3.0e-5", "< 1e-15", "< 1e-15"], ["2.5e-8", "2.0e-5", "< 1e-15", "< 1e-15"], ["1.0e-8", "1.4e-5", "< 1e-15", "< 1e-15"], ["2.8e-9", "9.3e-6", "< 1e-15", "< 1e-15"], ["9.0e-10", "6.2e-6", "< 1e-15", "< 1e-15"], ["9.6e-10", "4.2e-6", "< 1e-15", "< 1e-15"], ["2.8e-11", "2.8e-6", "< 1e-15", "< 1e-15"], ["4.4e-10", "1.9e-6", "< 1e-15", "< 1e-15"], ["2.1e-10", "1.3e-6", "< 1e-15", "< 1e-15"], ["8.9e-11", "8.7e-7", "< 1e-15", "< 1e-15"], ["3.1e-11", "5.8e-7", "< 1e-15", "< 1e-15"], ["1.4e-12", "3.9e-7", "< 1e-15", "< 1e-15"], ["1.3e-11", "2.6e-7", "< 1e-15", "< 1e-15"], ["5.9e-12", "1.8e-7", "< 1e-15", "< 1e-15"], ["2.2e-12", "1.2e-7", "< 1e-15", "< 1e-15"], ["4.0e-13", "8.1e-8", "< 1e-15", "< 1e-15"], ["5.1e-13", "5.5e-8", "< 1e-15", "< 1e-15"], ["5.2e-14", "3.7e-8", "< 1e-15", "< 1e-15"]];const left = 64, right = 600, steps = 40, width = 680;let current = -1;const show = (n) => { n = Math.max(0, Math.min(steps, n)); current = n; const x = left + (right - left) * n / steps; line.setAttribute('x1', x); line.setAttribute('x2', x); line.setAttribute('visibility', 'visible'); tip.replaceChildren(); const head = document.createElement('div'); head.className = 'viz-tooltip-head'; head.textContent = 'Step ' + n; tip.append(head); names.forEach((name, i) => {  if (!rows[n][i]) return;  const row = document.createElement('div'); row.className = 'viz-tooltip-row';  const key = document.createElement('span'); key.className = 'viz-swatch'; key.style.background = colors[i];  const value = document.createElement('strong'); value.textContent = rows[n][i];  const label = document.createElement('span'); label.textContent = name;  row.append(key, value, label); tip.append(row); }); const box = svg.getBoundingClientRect(), frame = fig.getBoundingClientRect(); const xs = box.left - frame.left + x * box.width / width; tip.hidden = false; const room = frame.width - tip.offsetWidth - 8; tip.style.left = Math.max(0, Math.min(room, xs + 12)) + 'px'; tip.style.top = (box.top - frame.top + 8) + 'px';};const hide = () => { line.setAttribute('visibility', 'hidden'); tip.hidden = true; current = -1; };svg.addEventListener('pointermove', (event) => { const box = svg.getBoundingClientRect(); const x = (event.clientX - box.left) * width / box.width; show(Math.round((x - left) / (right - left) * steps));});svg.addEventListener('pointerleave', hide);svg.addEventListener('focus', () => show(current < 0 ? 0 : current));svg.addEventListener('blur', hide);svg.addEventListener('keydown', (event) => { if (event.key === 'ArrowRight') { show(current + 1); event.preventDefault(); } if (event.key === 'ArrowLeft') { show(current - 1); event.preventDefault(); } if (event.key === 'Escape') hide();});})();</script>
</figure>

| method | steps until error $< 10^{-10}$ | measured | predicted |
|---|---|---|---|
| Bisection | 29 | error falls by $\approx 0.48$ per step, on average | the bound halves each step |
| Fixed point | 55 | error falls by $0.6736$ per step | $\sin\alpha = 0.6736$ |
| Secant | 6 | order $\approx 1.60$ | order $1.618$ |
| Newton | 4 | order $\approx 2.00$, constant $0.2208$ | order $2$, constant $0.2208$ |

Here $e_n = |x_n - \alpha|$. What I noticed:

- **Bisection** does not improve every single step. Its error goes up and down (for example $0.24$, then $0.011$, then $0.11$). Only the bound halves every time.
- **Fixed-point iteration** shrinks the error by exactly the predicted factor $\sin\alpha$, to four digits.
- **Secant and Newton** reach round-off level in 7 and 4 steps. For Newton, the constant is the limit of $e_{n+1}/e_n^2$. The error formula predicts $|f''(\alpha)|/(2|f'(\alpha)|) = \cos\alpha / (2(1 + \sin\alpha)) \approx 0.2208$, and the measured value matches to four digits.

Here are the four methods as I ran them, in Python:

```python
import math

def f(x):  return x - math.cos(x)
def df(x): return 1 + math.sin(x)

def bisection(a, b, steps):
    xs = []
    for _ in range(steps + 1):
        c = (a + b) / 2
        xs.append(c)
        if f(a) * f(c) <= 0:
            b = c
        else:
            a = c
    return xs

def fixed_point(x, steps):  # g(x) = cos(x)
    xs = [x]
    for _ in range(steps):
        x = math.cos(x)
        xs.append(x)
    return xs

def newton(x, steps):
    xs = [x]
    for _ in range(steps):
        x = x - f(x) / df(x)
        xs.append(x)
    return xs

def secant(x0, x1, steps):
    xs = [x0, x1]
    for _ in range(steps - 1):
        if f(x1) == f(x0):  # already converged
            break
        x0, x1 = x1, x1 - f(x1) * (x1 - x0) / (f(x1) - f(x0))
        xs.append(x1)
    return xs
```

The full script, with the measurements and the chart, is in [this site's repository](https://github.com/doubleyands/doubleyands.github.io/tree/main/scripts/posts/numerical-analysis).

## Summary

| method | order | always converges? | needs $f'$? | cost per step |
|---|---|---|---|---|
| Bisection | 1 (rate $1/2$) | yes, given a sign change | no | one value of $f$ |
| Newton | 2 | no, needs a good start | yes | $f$ and $f'$ |
| Secant | $\approx 1.62$ | no | no | one value of $f$ |
| Fixed point | 1 (rate $\lvert g'(\alpha)\rvert$) | yes, under the two conditions | no | one value of $g$ |

The main thing I take from this week: safety and speed pull in opposite directions. Bisection and contractions are safe but linear. Newton is quadratic but only near the root. And the fixed-point view, $g_N'(\alpha) = 0$, explains both sides of Newton at once.
