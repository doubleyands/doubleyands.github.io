---
title: "Numerical PDEs, Lecture 8: Forward Euler Converges, Zero-Stability, and Absolute Stability"
description: "My notes on the eighth lecture: a full proof that Forward Euler converges, the same proof for any one-step method, what zero-stability does and does not promise, an example that converges but blows up, and the interval of absolute stability."
pubDate: 2026-10-29
project: numerical-pde
tags: [numerical PDEs, ODEs, time stepping, stability, convergence, notes]
---

This is my summary of the eighth lecture in my numerical PDE course. I write the ideas in my own words and check them with a small experiment.

[Lecture 7](/blog/2026-10-27-numerical-pde-lecture-7/) ended with a definition: a method is stable if the inverse of its matrix stays bounded as the grid gets finer. This lecture applies that idea to time stepping. First we **prove** that Forward Euler converges. The same proof works for every one-step method. Then we see that this proof gives a weak promise, called **zero-stability**. An example shows a method that converges and still blows up. This leads to a second, more practical idea: **absolute stability**.

## 1. Recap: why the inverse matters

For a linear finite difference method, the error $E$ and the local truncation error $\tau$ satisfy

$$
A E = -\tau \qquad\Longrightarrow\qquad \|E\| \le \|A^{-1}\|\,\|\tau\|.
$$

The two factors behave differently:

- $\|\tau\|$ is the friendly one. It gets smaller when the grid gets finer. That is what "truncation error" means.
- $\|A^{-1}\|$ can be the dangerous one. A finer grid means a bigger matrix. If $\|A^{-1}\|$ stays bounded independently of the grid, the method is **stable**.

Here $A^{-1}$ is only notation. Nobody computes the inverse in practice. We solve the system.

**Why care so much?** A computer cannot store most numbers exactly. Even $0.1$ has no finite binary expansion ($0.000110011\ldots_2$), so it is cut off and carries a tiny error. A computer only adds, subtracts, multiplies and divides, but it does this billions of times. For a big matrix, these tiny errors pile up. For some matrices the growth is mild. For others it is fast. A good method must work every time, not only sometimes.

## 2. Consistency + stability ⇒ convergence

Three words to keep apart:

- **Consistency**: the truncation error goes to zero as the grid is refined. This is the easy part. A Taylor expansion usually shows it.
- **Convergence**: the numerical solution goes to the exact solution of the **differential** equation. This is what we really want, because we want to solve the ODE or PDE, not the difference equation.
- **Stability**: the bound on $\|A^{-1}\|$ from above.

Convergence is hard to prove directly. The main message of this part of the course is

$$
\text{consistency} + \text{stability} \;\Longrightarrow\; \text{convergence}.
$$

For linear problems even more is true: if the method is consistent, stability and convergence are **equivalent**. This is the **Lax equivalence theorem**. Its precise form comes later in the course.

## 3. A time-stepping scheme is also a linear system

At first, time stepping looks different from the boundary value problems of [Lecture 3](/blog/2026-10-13-numerical-pde-lecture-3/): we march forward and never "solve a system". But we can write it as one. Take $u' = \lambda u + g(t)$ and Forward Euler with step $\Delta t$:

$$
\frac{U^{n+1} - U^n}{\Delta t} - \lambda U^n = g(t_n), \qquad n = 0, 1, \dots, N - 1.
$$

Collect all unknowns $U^1, \dots, U^N$ in one vector. Then all $N$ equations together form $AU = F$, where $A$ is **lower triangular**: $1/\Delta t$ on the diagonal and $-(1 + \Delta t\lambda)/\Delta t$ just below it. (The known $U^0$ moves to the right side.) Solving a lower triangular system from the top down is exactly time stepping, one row at a time.

So the definition of stability from Lecture 7 applies directly. The inverse is easy to write down:

$$
(A^{-1})_{ij} = \Delta t\,(1 + \Delta t\lambda)^{\,i - j} \quad (i \ge j), \qquad
\|A^{-1}\|_\infty = \Delta t \sum_{m=0}^{N-1} |1 + \Delta t\lambda|^m.
$$

Since $|1 + \Delta t\lambda| \le 1 + \Delta t|\lambda| \le e^{\Delta t|\lambda|}$ and $N\Delta t = T$, this gives

$$
\|A^{-1}\|_\infty \le T e^{|\lambda| T},
$$

a bound that does not depend on $\Delta t$. So Forward Euler is stable in the sense of Lecture 7. The next section proves the same thing without matrices, and for nonlinear $f$ too. (Keep the factor $e^{|\lambda| T}$ in mind. It will come back.)

## 4. Proof: Forward Euler converges

**Setting.** Solve $u'(t) = f(u(t))$ with $u(0) = \eta$ up to a fixed final time $T$. Assume $f$ is **Lipschitz** with constant $L$. By [Lecture 6](/blog/2026-10-22-numerical-pde-lecture-6/), the solution then exists and is unique. The scheme is $U^{n+1} = U^n + \Delta t\, f(U^n)$, with $\Delta t = T/N$ and $t_n = n\Delta t$.

**Truncation error.** Put the exact solution into the difference quotient and use $u' = f(u)$:

$$
\tau^n = \frac{u(t_{n+1}) - u(t_n)}{\Delta t} - f(u(t_n)) = \frac{\Delta t}{2}\,u''(t_n) + O(\Delta t^2).
$$

So the method is consistent and first order, if $u''$ is bounded. But this says nothing yet about the real error $E^n = U^n - u(t_n)$.

**Step 1: the error equation.** The exact solution misses the scheme by $\Delta t\,\tau^n$:

$$
u(t_{n+1}) = u(t_n) + \Delta t\, f(u(t_n)) + \Delta t\,\tau^n.
$$

Subtract this from the scheme:

$$
E^{n+1} = E^n + \Delta t\,\big(f(U^n) - f(u(t_n))\big) - \Delta t\,\tau^n.
$$

The next error is the present error plus something. The middle term is nonlinear, so we cannot handle it directly.

**Step 2: Lipschitz removes $f$.** By the Lipschitz condition, $|f(U^n) - f(u(t_n))| \le L|E^n|$. With the triangle inequality,

$$
|E^{n+1}| \le (1 + \Delta t L)\,|E^n| + \Delta t\,|\tau^n|.
$$

Now $f$ is gone and only errors are left. This is the typical role of the Lipschitz condition in proofs: we compare two solutions, their difference appears naturally, and the Lipschitz bound turns a difference of $f$ values into a difference of $u$ values. If $f$ is differentiable, $L$ is a bound on $|f'|$.

**Step 3: unroll it.** If the inequality were an equality, it would be a geometric sequence with ratio $1 + \Delta t L$. Insert the same inequality for $|E^n|$, then for $|E^{n-1}|$, and so on. For $n = 2$:

$$
|E^2| \le (1 + \Delta t L)^2 |E^0| + \Delta t\,(1 + \Delta t L)\,|\tau^0| + \Delta t\,|\tau^1|.
$$

In general,

$$
|E^n| \le (1 + \Delta t L)^n |E^0| + \Delta t \sum_{m=1}^{n} (1 + \Delta t L)^{n-m}\, |\tau^{m-1}|.
$$

An error made early gets multiplied by the ratio more times than an error made late. The factor $(1 + \Delta t L)^{n-m}$ carries an error from time $t_m$ to time $t_n$. This formula is a discrete version of Duhamel's principle from Lecture 6.

**Step 4: fix the final time.** We look at $0 \le n \le N$, so $n\Delta t \le T$. This matters. If $T$ were not fixed, we would be asking about long-time behavior, which is a different question.

**Step 5: the key inequality $1 + x \le e^x$.** As $\Delta t \to 0$, the ratio $1 + \Delta t L$ goes to $1$ but the power $n$ goes to infinity. This looks like the indeterminate form $1^\infty$. One simple inequality settles it:

$$
(1 + \Delta t L)^{n-m} \le e^{(n-m)\Delta t L} \le e^{n \Delta t L} \le e^{LT}.
$$

The factor $e^{LT}$ may look big, but it does **not** depend on $\Delta t$. However small the step, the constant stays the same.

**Step 6: the sum.** There are $n$ terms, each at most $\|\tau\|_\infty = \max_m |\tau^m|$, so

$$
\Delta t \sum_{m=1}^{n} |\tau^{m-1}| \le n\Delta t\,\|\tau\|_\infty \le T\,\|\tau\|_\infty.
$$

**Result.** With an exact initial value, $E^0 = 0$, and

$$
|E^n| \le e^{LT}\, T\, \|\tau\|_\infty = O(\Delta t) \qquad (\Delta t \to 0).
$$

**Forward Euler converges, and it is first order.** If $E^0 \ne 0$, the first term adds $e^{LT}|E^0|$, which is still bounded.

The proof also shows **where stability enters**. Each one-step error $\Delta t\,\tau^{m-1}$ grows by at most the factor $e^{LT}$ on its way to the end. This factor stays bounded as $\Delta t \to 0$. So the global error has the same order as the truncation error.

## 5. The bound is very pessimistic

The bound is true, but it can be far from the real error. The experiment below shows this. On a nonlinear problem with $L = 1$ the bound is about $20$ times too large. That is harmless.

Now take a linear problem $u' = \lambda(u - \sin t) + \cos t$ with $\lambda = -1000$ and $T = 2$. Here $L = |\lambda| = 1000$, so the factor is $e^{LT} = e^{2000} \approx 10^{868}$. The bound promises **nothing useful** until $\Delta t$ is around $10^{-868}$. Yet, as we will see, the real error is $10^{-6}$ or less for every step at or below $\Delta t = 0.002$. A convergence proof can be far from what we need in practice.

## 6. Any one-step method: the same proof

Every one-step method, with any number of stages (RK2, RK4 from Lecture 7), has the form

$$
U^{n+1} = U^n + \Delta t\, \Phi(U^n, t_n, \Delta t).
$$

In the area view of Lecture 6, $\Delta t\,\Phi$ is the method's guess for the area under $f$ in one step. Assume $\Phi$ is continuous in $t$ and $\Delta t$, and Lipschitz in $u$ with some constant $L'$. Then the proof above goes through **word for word**, with $L'$ in place of $L$.

Two conditions to check:

- **Consistency** means $\Phi(u, t, 0) = f(u, t)$: with a zero step, the method uses the true slope.
- **Lipschitz** must hold for $\Phi$, not only for $f$.

**Example: RK2.** For $u' = f(u)$, RK2 has $\Phi(u, \Delta t) = f\big(u + \tfrac12\Delta t\, f(u)\big)$. Consistency is clear: put $\Delta t = 0$. For Lipschitz, use the Lipschitz condition of $f$ twice, once for the outer $f$ and once for the inner $f$:

$$
\begin{aligned}
|\Phi(u) - \Phi(v)| &\le L\,\Big|(u - v) + \tfrac12\Delta t\big(f(u) - f(v)\big)\Big| \\
&\le L\Big(|u - v| + \tfrac12\Delta t\,L\,|u - v|\Big) = \Big(L + \tfrac12\Delta t\,L^2\Big)\,|u - v|.
\end{aligned}
$$

So $L' = L + \tfrac12\Delta t L^2$. Nesting evaluations keeps the Lipschitz property, and the convergence proof applies. The experiment shows that this constant is sharp on my test problem.

## 7. Zero-stability: what the proof really promises

The proof bounds the global error by a sum of one-step errors, so the global error has the same order as the truncation error as $\Delta t \to 0$. This kind of stability is called **zero-stability**. The "zero" refers to the limit $\Delta t \to 0$.

- **What it gives:** if we keep making $\Delta t$ smaller, the method **eventually** converges.
- **What it does not give:** any promise for a **specific** $\Delta t$. Is $\Delta t = 1/1000$ small enough? Zero-stability does not say. For some problems it is not.

By the proof above, every one-step method with a Lipschitz $\Phi$ is zero-stable automatically. Multistep methods need a separate check, which comes in a later lecture.

## 8. An example that converges but blows up

My own test problem:

$$
u'(t) = \lambda\,\big(u - \sin t\big) + \cos t, \qquad u(0) = 0.
$$

It is built so that $u = \sin t$ is a solution: the bracket vanishes and $u' = \cos t$. But is it **the** solution, or just **a** solution? Look at the right side $f(u, t) = \lambda(u - \sin t) + \cos t$. Then

$$
|f(u, t) - f(v, t)| = |\lambda|\,|u - v|,
$$

so $f$ is Lipschitz in $u$ with $L = |\lambda|$, and the solution is unique (Lecture 6). The solution we guessed is the only one. Compare $u' = 3u^{2/3}$ in Lecture 6, which had two solutions because it is not Lipschitz at $0$.

Now take $\lambda = -1000$, $T = 2$, and Forward Euler with $\Delta t = 1/400$. We expect $U^{800} \approx \sin 2 \approx 0.909$. A sine stays between $-1$ and $1$. The experiment gives $U^{800} \approx 1.56 \times 10^{131}$. Note that $1/400$ is not a big step. And the method **does** converge: for small enough $\Delta t$ the error goes to zero, as the theorem says. So "make the step smaller and hope" is not a plan. We want to **know in advance** which steps are safe.

## 9. Absolute stability

To find the safe steps, simplify the problem. The forcing terms do not depend on $u$, so drop them. Keep only the linear part $\lambda u$:

$$
u' = \lambda u, \qquad u(t) = e^{\lambda t} u(0).
$$

This is the **test equation**. We know everything about it. Forward Euler gives

$$
U^{n+1} = U^n + \Delta t\,\lambda U^n = (1 + \lambda\Delta t)\,U^n.
$$

Again a geometric sequence. The ratio $1 + \lambda\Delta t$ is the **amplification factor**: each step multiplies every error by it. If its size is above $1$, say $2$, then after $n$ steps the error has grown by $2^n$. So we need its size to be at most $1$.

Put $z = \lambda\Delta t$. For real $\lambda$:

$$
|1 + z| \le 1 \quad\Longleftrightarrow\quad -2 \le z \le 0.
$$

This interval $[-2, 0]$ is the **interval of absolute stability** of Forward Euler. Note that $\Delta t$ and $\lambda$ never appear alone, only through the product $z$.

For the example, $\lambda = -1000$, so we need $-1000\,\Delta t \ge -2$, that is,

$$
\Delta t \le \frac{2}{1000} = 0.002.
$$

The step $1/400 = 0.0025$ gives $z = -2.5$ and $|1 + z| = 1.5$. Over $800$ steps, errors can grow by $1.5^{800} \approx 10^{141}$. The errors that grow start from the small truncation errors, not from size $1$, so the final value is a bit smaller ($10^{131}$). But this growth is the blow-up we saw. The limit is not luck: the analysis tells us in advance.

For a nonlinear problem we linearize: $\lambda$ is replaced by $\partial f/\partial u$ (or, for systems, the eigenvalues of the Jacobian matrix). In the example, the linear part is exactly the Lipschitz constant: $|\lambda| = L$. When $\lambda$ is complex, the interval becomes a region in the complex plane. That is the topic of the next lecture.

**Two kinds of stability, side by side.**

| | Zero-stability | Absolute stability |
|---|---|---|
| Question | Does the error go to zero as $\Delta t \to 0$? | Does the error stay bounded for **this** $\Delta t$? |
| Depends on | the method only | the method **and** $z = \lambda\Delta t$ |
| Forward Euler | always (one-step method) | only if $-2 \le \lambda\Delta t \le 0$ (real $\lambda$) |

## 10. Experiment

The script uses the Python standard library only.

**Part 1: the proof on a nonlinear problem.** I take $u' = \cos u$, $u(0) = 0$, up to $T = 2$. The exact solution is $u(t) = 2\arctan\big(\tanh(t/2)\big)$, and $f = \cos$ is Lipschitz with $L = 1$. I compute the max error over all grid times, the max truncation error $\|\tau\|_\infty$ from the exact solution, and the bound $e^{LT}\,T\,\|\tau\|_\infty$.

| $\Delta t$ | FE error | $\lVert\tau\rVert_\infty$ | bound | bound / error | RK2 error |
|---|---|---|---|---|---|
| 1/5 | $3.83 \times 10^{-2}$ | $4.99 \times 10^{-2}$ | $0.737$ | 19.2 | $1.59 \times 10^{-3}$ |
| 1/20 | $9.29 \times 10^{-3}$ | $1.25 \times 10^{-2}$ | $0.185$ | 19.9 | $9.07 \times 10^{-5}$ |
| 1/80 | $2.31 \times 10^{-3}$ | $3.13 \times 10^{-3}$ | $4.62 \times 10^{-2}$ | 20.0 | $5.54 \times 10^{-6}$ |
| 1/320 | $5.75 \times 10^{-4}$ | $7.81 \times 10^{-4}$ | $1.15 \times 10^{-2}$ | 20.1 | $3.45 \times 10^{-7}$ |
| 1/640 | $2.88 \times 10^{-4}$ | $3.91 \times 10^{-4}$ | $5.77 \times 10^{-3}$ | 20.1 | $8.61 \times 10^{-8}$ |

- The error is always below the bound, as the proof says. The bound is about $20$ times too large, the same factor for every $\Delta t$. So the bound has the right **order**, just a large constant.
- The measured order of Forward Euler is $1.00$, and error$/\Delta t$ settles at $0.184$. The truncation error is about $\Delta t/4$, because $u'' = -\tfrac12 \sin 2u$ has size at most $\tfrac12$.
- RK2 has order $2.00$ (for example $2.01$ from $1/80$ to $1/160$).
- The Lipschitz constant of the RK2 function $\Phi$, measured as the largest difference quotient on $u \in [-4, 4]$, is $1.2500$, $1.1000$ and $1.0250$ for $\Delta t = 0.5$, $0.2$ and $0.05$. The formula $L + \tfrac12\Delta t L^2$ gives exactly the same numbers. So the bound is sharp here. (It is reached near $u = -\pi/2$, where $f' = 1$.)

**Part 2: converges, but blows up.** The linear problem from section 8 with $\lambda = -5$ and $\lambda = -1000$, $T = 2$, and $\Delta t = 2/N$. The chart shows the max error over the grid times.

<figure class="viz">
<div class="viz-legend"><span class="viz-key"><span class="viz-swatch" style="background:#2a78d6"></span>λ = −5</span><span class="viz-key"><span class="viz-swatch" style="background:#eb6834"></span>λ = −1000</span></div>
<div class="viz-frame">
<svg viewBox="0 0 680 400" role="img" tabindex="0" aria-label="Forward Euler error against the time step on log scales, for lambda = -5 and lambda = -1000. Both fall with slope 1 for small steps. For lambda = -1000 the error is huge for every step larger than 2/1000 and drops suddenly to about 1e-6 at that line." class="viz-plot"><line x1="64" x2="520" y1="30.0" y2="30.0" class="viz-grid"/><text x="54" y="34.0" text-anchor="end" class="viz-tick">10<tspan dy="-0.5em" font-size="75%">4</tspan></text><line x1="64" x2="520" y1="81.7" y2="81.7" class="viz-grid"/><text x="54" y="85.7" text-anchor="end" class="viz-tick">10<tspan dy="-0.5em" font-size="75%">2</tspan></text><line x1="64" x2="520" y1="133.3" y2="133.3" class="viz-grid"/><text x="54" y="137.3" text-anchor="end" class="viz-tick">1</text><line x1="64" x2="520" y1="185.0" y2="185.0" class="viz-grid"/><text x="54" y="189.0" text-anchor="end" class="viz-tick">10<tspan dy="-0.5em" font-size="75%">−2</tspan></text><line x1="64" x2="520" y1="236.7" y2="236.7" class="viz-grid"/><text x="54" y="240.7" text-anchor="end" class="viz-tick">10<tspan dy="-0.5em" font-size="75%">−4</tspan></text><line x1="64" x2="520" y1="288.3" y2="288.3" class="viz-grid"/><text x="54" y="292.3" text-anchor="end" class="viz-tick">10<tspan dy="-0.5em" font-size="75%">−6</tspan></text><line x1="64" x2="520" y1="340.0" y2="340.0" class="viz-grid"/><text x="54" y="344.0" text-anchor="end" class="viz-tick">10<tspan dy="-0.5em" font-size="75%">−8</tspan></text><text x="64.0" y="362" text-anchor="middle" class="viz-tick">1/10</text><text x="215.5" y="362" text-anchor="middle" class="viz-tick">1/100</text><text x="367.0" y="362" text-anchor="middle" class="viz-tick">1/1000</text><text x="518.4" y="362" text-anchor="middle" class="viz-tick">1/10000</text><text x="292.0" y="388" text-anchor="middle" class="viz-axis-title">Time step Δt (smaller to the right)</text><text x="64" y="14" class="viz-axis-title">Max error on [0, 2]</text><line x1="321.4" x2="321.4" y1="30" y2="340" stroke="#86868b" stroke-width="1.5" stroke-dasharray="5 5"/><text x="315.4" y="330" text-anchor="end" class="viz-tick">Δt = 2/1000</text><path d="M64.0,185.1 L109.6,192.9 L155.2,200.7 L200.8,208.5 L246.4,216.3 L292.0,224.1 L306.7,226.6 L314.4,227.9 L318.0,228.5 L320.7,229.0 L321.4,229.1 L322.0,229.2 L324.6,229.7 L327.6,230.2 L337.6,231.9 L383.2,239.7 L428.8,247.4 L474.4,255.2 L520.0,263.0" fill="none" stroke="#2a78d6" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><circle cx="64.0" cy="185.1" r="3" fill="#2a78d6" stroke="#ffffff" stroke-width="1.5"/><circle cx="109.6" cy="192.9" r="3" fill="#2a78d6" stroke="#ffffff" stroke-width="1.5"/><circle cx="155.2" cy="200.7" r="3" fill="#2a78d6" stroke="#ffffff" stroke-width="1.5"/><circle cx="200.8" cy="208.5" r="3" fill="#2a78d6" stroke="#ffffff" stroke-width="1.5"/><circle cx="246.4" cy="216.3" r="3" fill="#2a78d6" stroke="#ffffff" stroke-width="1.5"/><circle cx="292.0" cy="224.1" r="3" fill="#2a78d6" stroke="#ffffff" stroke-width="1.5"/><circle cx="306.7" cy="226.6" r="3" fill="#2a78d6" stroke="#ffffff" stroke-width="1.5"/><circle cx="314.4" cy="227.9" r="3" fill="#2a78d6" stroke="#ffffff" stroke-width="1.5"/><circle cx="318.0" cy="228.5" r="3" fill="#2a78d6" stroke="#ffffff" stroke-width="1.5"/><circle cx="320.7" cy="229.0" r="3" fill="#2a78d6" stroke="#ffffff" stroke-width="1.5"/><circle cx="321.4" cy="229.1" r="3" fill="#2a78d6" stroke="#ffffff" stroke-width="1.5"/><circle cx="322.0" cy="229.2" r="3" fill="#2a78d6" stroke="#ffffff" stroke-width="1.5"/><circle cx="324.6" cy="229.7" r="3" fill="#2a78d6" stroke="#ffffff" stroke-width="1.5"/><circle cx="327.6" cy="230.2" r="3" fill="#2a78d6" stroke="#ffffff" stroke-width="1.5"/><circle cx="337.6" cy="231.9" r="3" fill="#2a78d6" stroke="#ffffff" stroke-width="1.5"/><circle cx="383.2" cy="239.7" r="3" fill="#2a78d6" stroke="#ffffff" stroke-width="1.5"/><circle cx="428.8" cy="247.4" r="3" fill="#2a78d6" stroke="#ffffff" stroke-width="1.5"/><circle cx="474.4" cy="255.2" r="3" fill="#2a78d6" stroke="#ffffff" stroke-width="1.5"/><circle cx="520.0" cy="263.0" r="3" fill="#2a78d6" stroke="#ffffff" stroke-width="1.5"/><path d="M64.0,30.0 L109.6,30.0 L155.2,30.0 L200.8,30.0 L246.4,30.0 L292.0,30.0 L306.7,30.0 L314.4,30.0 L318.0,30.0 L320.7,156.1 L321.4,288.3 L322.0,288.4 L324.6,288.9 L327.6,289.4 L337.6,291.1 L383.2,298.9 L428.8,306.7 L474.4,314.4 L520.0,322.2" fill="none" stroke="#eb6834" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><circle cx="64.0" cy="30.0" r="3.5" fill="#ffffff" stroke="#eb6834" stroke-width="2"/><circle cx="109.6" cy="30.0" r="3.5" fill="#ffffff" stroke="#eb6834" stroke-width="2"/><circle cx="155.2" cy="30.0" r="3.5" fill="#ffffff" stroke="#eb6834" stroke-width="2"/><circle cx="200.8" cy="30.0" r="3.5" fill="#ffffff" stroke="#eb6834" stroke-width="2"/><circle cx="246.4" cy="30.0" r="3.5" fill="#ffffff" stroke="#eb6834" stroke-width="2"/><circle cx="292.0" cy="30.0" r="3.5" fill="#ffffff" stroke="#eb6834" stroke-width="2"/><circle cx="306.7" cy="30.0" r="3.5" fill="#ffffff" stroke="#eb6834" stroke-width="2"/><circle cx="314.4" cy="30.0" r="3.5" fill="#ffffff" stroke="#eb6834" stroke-width="2"/><circle cx="318.0" cy="30.0" r="3.5" fill="#ffffff" stroke="#eb6834" stroke-width="2"/><circle cx="320.7" cy="156.1" r="3" fill="#eb6834" stroke="#ffffff" stroke-width="1.5"/><circle cx="321.4" cy="288.3" r="3" fill="#eb6834" stroke="#ffffff" stroke-width="1.5"/><circle cx="322.0" cy="288.4" r="3" fill="#eb6834" stroke="#ffffff" stroke-width="1.5"/><circle cx="324.6" cy="288.9" r="3" fill="#eb6834" stroke="#ffffff" stroke-width="1.5"/><circle cx="327.6" cy="289.4" r="3" fill="#eb6834" stroke="#ffffff" stroke-width="1.5"/><circle cx="337.6" cy="291.1" r="3" fill="#eb6834" stroke="#ffffff" stroke-width="1.5"/><circle cx="383.2" cy="298.9" r="3" fill="#eb6834" stroke="#ffffff" stroke-width="1.5"/><circle cx="428.8" cy="306.7" r="3" fill="#eb6834" stroke="#ffffff" stroke-width="1.5"/><circle cx="474.4" cy="314.4" r="3" fill="#eb6834" stroke="#ffffff" stroke-width="1.5"/><circle cx="520.0" cy="322.2" r="3" fill="#eb6834" stroke="#ffffff" stroke-width="1.5"/><text x="530" y="267.0" class="viz-label">λ = −5</text><text x="530" y="326.2" class="viz-label">λ = −1000</text><line x1="0" x2="0" y1="30" y2="340" class="viz-crosshair" visibility="hidden"/></svg>
<div class="viz-tooltip" hidden></div>
</div>
<figcaption>Forward Euler on u′ = λ(u − sin t) + cos t, u(0) = 0, up to t = 2, on log scales. Hollow points on the top edge are errors larger than 10<sup>4</sup> (up to 10<sup>222</sup>; see the table). The dashed line is the limit Δt = 2/|λ| for λ = −1000. Hover over the chart, or focus it and use the arrow keys, to read the values.</figcaption>
<details class="viz-table"><summary>Show the numbers</summary><table><thead><tr><th>Δt</th><th>λ = −5</th><th>λ = −1000</th></tr></thead><tbody><tr><td>2/20</td><td>9.9e-3</td><td>1.3e34</td></tr><tr><td>2/40</td><td>4.9e-3</td><td>1.6e61</td></tr><tr><td>2/80</td><td>2.5e-3</td><td>2.4e103</td></tr><tr><td>2/160</td><td>1.2e-3</td><td>1.0e162</td></tr><tr><td>2/320</td><td>6.1e-4</td><td>9.6e221</td></tr><tr><td>2/640</td><td>3.1e-4</td><td>2.1e199</td></tr><tr><td>2/800</td><td>2.5e-4</td><td>1.6e131</td></tr><tr><td>2/900</td><td>2.2e-4</td><td>7.8e68</td></tr><tr><td>2/950</td><td>2.1e-4</td><td>6.2e31</td></tr><tr><td>2/990</td><td>2.0e-4</td><td>1.3e-1</td></tr><tr><td>2/1000</td><td>2.0e-4</td><td>1.0e-6</td></tr><tr><td>2/1010</td><td>1.9e-4</td><td>9.9e-7</td></tr><tr><td>2/1050</td><td>1.9e-4</td><td>9.5e-7</td></tr><tr><td>2/1100</td><td>1.8e-4</td><td>9.1e-7</td></tr><tr><td>2/1280</td><td>1.5e-4</td><td>7.8e-7</td></tr><tr><td>2/2560</td><td>7.7e-5</td><td>3.9e-7</td></tr><tr><td>2/5120</td><td>3.8e-5</td><td>2.0e-7</td></tr><tr><td>2/10240</td><td>1.9e-5</td><td>9.8e-8</td></tr><tr><td>2/20480</td><td>9.6e-6</td><td>4.9e-8</td></tr></tbody></table></details>
<script>(() => {const fig = document.currentScript.closest('figure');const svg = fig.querySelector('svg');const line = svg.querySelector('.viz-crosshair');const tip = fig.querySelector('.viz-tooltip');const names = ["\u03bb = \u22125", "\u03bb = \u22121000"]; const colors = ["#2a78d6", "#eb6834"];const rows = [["9.9e-3", "1.3e34"], ["4.9e-3", "1.6e61"], ["2.5e-3", "2.4e103"], ["1.2e-3", "1.0e162"], ["6.1e-4", "9.6e221"], ["3.1e-4", "2.1e199"], ["2.5e-4", "1.6e131"], ["2.2e-4", "7.8e68"], ["2.1e-4", "6.2e31"], ["2.0e-4", "1.3e-1"], ["2.0e-4", "1.0e-6"], ["1.9e-4", "9.9e-7"], ["1.9e-4", "9.5e-7"], ["1.8e-4", "9.1e-7"], ["1.5e-4", "7.8e-7"], ["7.7e-5", "3.9e-7"], ["3.8e-5", "2.0e-7"], ["1.9e-5", "9.8e-8"], ["9.6e-6", "4.9e-8"]]; const heads = ["\u0394t = 2/20", "\u0394t = 2/40", "\u0394t = 2/80", "\u0394t = 2/160", "\u0394t = 2/320", "\u0394t = 2/640", "\u0394t = 2/800", "\u0394t = 2/900", "\u0394t = 2/950", "\u0394t = 2/990", "\u0394t = 2/1000", "\u0394t = 2/1010", "\u0394t = 2/1050", "\u0394t = 2/1100", "\u0394t = 2/1280", "\u0394t = 2/2560", "\u0394t = 2/5120", "\u0394t = 2/10240", "\u0394t = 2/20480"];const xs = [64.0, 109.6, 155.2, 200.8, 246.4, 292.0, 306.7, 314.4, 318.0, 320.7, 321.4, 322.0, 324.6, 327.6, 337.6, 383.2, 428.8, 474.4, 520.0]; const steps = 18, width = 680;let current = -1;const show = (n) => { n = Math.max(0, Math.min(steps, n)); current = n; const x = xs[n]; line.setAttribute('x1', x); line.setAttribute('x2', x); line.setAttribute('visibility', 'visible'); tip.replaceChildren(); const head = document.createElement('div'); head.className = 'viz-tooltip-head'; head.textContent = heads[n]; tip.append(head); names.forEach((name, i) => {  const row = document.createElement('div'); row.className = 'viz-tooltip-row';  const key = document.createElement('span'); key.className = 'viz-swatch'; key.style.background = colors[i];  const value = document.createElement('strong'); value.textContent = rows[n][i];  const label = document.createElement('span'); label.textContent = name;  row.append(key, value, label); tip.append(row); }); const box = svg.getBoundingClientRect(), frame = fig.getBoundingClientRect(); const xp = box.left - frame.left + x * box.width / width; tip.hidden = false; const room = frame.width - tip.offsetWidth - 8; tip.style.left = Math.max(0, Math.min(room, xp + 12)) + 'px'; tip.style.top = (box.top - frame.top + 8) + 'px';};const hide = () => { line.setAttribute('visibility', 'hidden'); tip.hidden = true; current = -1; };svg.addEventListener('pointermove', (event) => { const box = svg.getBoundingClientRect(); const x = (event.clientX - box.left) * width / box.width; let best = 0; xs.forEach((v, i) => { if (Math.abs(v - x) < Math.abs(xs[best] - x)) best = i; }); show(best);});svg.addEventListener('pointerleave', hide);svg.addEventListener('focus', () => show(current < 0 ? 0 : current));svg.addEventListener('blur', hide);svg.addEventListener('keydown', (event) => { if (event.key === 'ArrowRight') { show(current + 1); event.preventDefault(); } if (event.key === 'ArrowLeft') { show(current - 1); event.preventDefault(); } if (event.key === 'Escape') hide();});})();</script>
</figure>

| $\Delta t$ | $z$ for $\lambda = -1000$ | error, $\lambda = -5$ | error, $\lambda = -1000$ |
|---|---|---|---|
| 2/20 = 0.1 | $-100$ | $9.89 \times 10^{-3}$ | $1.32 \times 10^{34}$ |
| 2/320 | $-6.25$ | $6.13 \times 10^{-4}$ | $9.56 \times 10^{221}$ |
| 2/800 | $-2.5$ | $2.45 \times 10^{-4}$ | $1.56 \times 10^{131}$ |
| 2/950 | $-2.105$ | $2.07 \times 10^{-4}$ | $6.16 \times 10^{31}$ |
| 2/990 | $-2.020$ | $1.98 \times 10^{-4}$ | $1.31 \times 10^{-1}$ |
| 2/1000 = 0.002 | $-2$ | $1.96 \times 10^{-4}$ | $1.00 \times 10^{-6}$ |
| 2/1010 | $-1.980$ | $1.94 \times 10^{-4}$ | $9.90 \times 10^{-7}$ |
| 2/2560 | $-0.781$ | $7.66 \times 10^{-5}$ | $3.91 \times 10^{-7}$ |
| 2/20480 | $-0.098$ | $9.58 \times 10^{-6}$ | $4.88 \times 10^{-8}$ |

What I noticed:

- **Zero-stability** is visible in both lines: for small $\Delta t$, both fall with slope $1$. Even $\lambda = -1000$ converges, as the theorem promises.
- **Absolute stability** is visible in the jump. For $\lambda = -1000$, every step larger than $2/1000$ gives a huge error. A step only $1\%$ larger than the limit ($z = -2.02$) already gives an error of $0.13$. At $z = -2.105$ (about $5\%$ larger), the error is $6 \times 10^{31}$. The sharp edge is exactly where $|1 + z|$ crosses $1$.
- **Inside the limit, $\lambda = -1000$ is more accurate**, by a factor of about $200$. My reading: when $|1 + z| < 1$, old errors are damped quickly, so the error stays near a balance where each step's new error is cancelled by damping. That balance is about $\|\tau\|/|\lambda|$. Here $\|\tau\| \approx \Delta t/2$, so the error is about $\Delta t/2000$. The data agree: $\Delta t = 2/2560$ gives $3.906 \times 10^{-7}$, and $\Delta t/2000 = 3.906 \times 10^{-7}$.
- **The proven bound is useless here.** For $\lambda = -5$, the bound $T e^{|\lambda|T}\|\tau\|_\infty$ at $\Delta t = 2/20480$ is $2.2$, while the real error is $9.6 \times 10^{-6}$. For $\lambda = -1000$ the bound is between $10^{864}$ and $10^{868}$ for every step in the table.

**Part 3: the inverse of the triangular matrix.** From section 3, $\|A^{-1}\|_\infty = \Delta t \sum_{m=0}^{N-1} |1 + \Delta t\lambda|^m$. I compute it in logarithms because it overflows a float.

| $\Delta t$ | $\lambda = -5$ | $\lambda = -1000$ |
|---|---|---|
| 2/20 | $0.200$ | $8.3 \times 10^{36}$ |
| 2/320 | $0.200$ | $4.2 \times 10^{227}$ |
| 2/990 | $0.200$ | $4.0 \times 10^{7}$ |
| 2/1000 | $0.200$ | $2.00$ |
| 2/1010 | $0.200$ | $0.100$ |
| 2/2560 | $0.200$ | $1.00 \times 10^{-3}$ |
| bound $T e^{\lvert\lambda\rvert T}$ | $4.4 \times 10^{4}$ | $\approx 10^{869}$ |

The numbers tell the same story in matrix language. Every value is below the bound $T e^{|\lambda| T}$, so the method is stable in the sense of Lecture 7 for both $\lambda$. But for $\lambda = -1000$ that "bounded" norm can be $10^{227}$. Well inside the limit, the norm is about $1/|\lambda|$ ($0.2$ and $0.001$), which explains the small errors in Part 2. At $z = -2$ exactly, $|1 + z| = 1$, so the sum is just $N\Delta t = T = 2$.

The full script and the chart code are in [this site's repository](https://github.com/doubleyands/doubleyands.github.io/tree/main/scripts/posts/numerical-pde/lecture-8).

## Summary

1. For a linear method, $\|E\| \le \|A^{-1}\|\,\|\tau\|$. **Stability** means $\|A^{-1}\|$ is bounded independently of the grid. **Consistency + stability ⇒ convergence**, and for linear problems the **Lax equivalence theorem** says stability and convergence are equivalent for a consistent method.
2. A time-stepping scheme is a **lower triangular** system, and marching in time is forward substitution.
3. **Forward Euler converges with order 1** if $f$ is Lipschitz. The proof: error equation, Lipschitz removes $f$, unroll the geometric sequence, $1 + x \le e^x$, and a fixed final time $T$. The result is $|E^n| \le e^{LT}\,T\,\|\tau\|_\infty$.
4. The same proof works for **any one-step method** $U^{n+1} = U^n + \Delta t\,\Phi$ with a Lipschitz $\Phi$. For RK2, $L' = L + \tfrac12\Delta t L^2$.
5. **Zero-stability** promises convergence as $\Delta t \to 0$, but nothing for a specific $\Delta t$. The constant $e^{LT}$ can be astronomically large.
6. A problem with $\lambda = -1000$ **converges but blows up** for $\Delta t$ above $2/1000$.
7. **Absolute stability** uses the test equation $u' = \lambda u$. Forward Euler multiplies by $1 + z$ with $z = \lambda\Delta t$, and it is stable for real $\lambda$ when $-2 \le z \le 0$.
