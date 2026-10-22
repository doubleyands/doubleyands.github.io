---
title: "Numerical PDEs, Lecture 6: Well-Posed ODEs, Lipschitz Continuity, and the First Time-Stepping Schemes"
description: "My notes on the sixth lecture: when an ODE has exactly one solution, the Lipschitz condition, blow-up and non-uniqueness, and the first time-stepping schemes (Forward Euler, Backward Euler with Newton's method, the trapezoidal method, Taylor methods), all seen as ways to approximate one integral."
pubDate: 2026-10-22
project: numerical-pde
tags: [numerical PDEs, ODEs, time stepping, Euler method, Newton's method, notes]
---

This is my summary of the sixth lecture in my numerical PDE course. I write the ideas in my own words and check them with a small experiment.

[Lecture 5](/blog/2026-10-20-numerical-pde-lecture-5/) started ODE initial value problems. This lecture has two halves. The first half asks a question we should ask before any computation: **does the problem have a solution, and only one?** The second half builds the first **time-stepping schemes**. At the end, all of these schemes turn out to be one idea: different ways to approximate one integral.

## 1. Duhamel's principle from an integrating factor

In Lecture 5 we checked Duhamel's formula. Here is where it comes from. Take

$$
u' = Au + g(t), \qquad u(t_0) = \eta,
$$

with a constant matrix $A$. Move $Au$ to the left and multiply by $e^{-A(t - t_0)}$. By the product rule, the left side becomes one derivative:

$$
\frac{d}{dt}\Big(e^{-A(t - t_0)}\,u(t)\Big) = e^{-A(t - t_0)}\,g(t).
$$

Integrate from $t_0$ to $t$ and multiply by $e^{A(t - t_0)}$. This gives the formula directly:

$$
u(t) = e^{A(t - t_0)}\,\eta + \int_{t_0}^{t} e^{A(t - s)}\,g(s)\,ds.
$$

Now suppose $g$ depends on $u$ too. The same steps give

$$
u(t) = e^{A(t - t_0)}\,\eta + \int_{t_0}^{t} e^{A(t - s)}\,g\big(u(s)\big)\,ds.
$$

This is no longer a formula, because the unknown $u$ is inside the integral. It is an **integral equation**. So a differential equation can also be read as an integral equation. We come back to this idea in section 12.

## 2. Well-posedness

Before we attack a problem with a computer, we should know that the problem makes sense. Hadamard called a problem **well-posed** if three things hold:

1. **Existence**: a solution exists.
2. **Uniqueness**: there is only one solution.
3. **Continuous dependence**: a small change in the data gives only a small change in the solution.

We already met failures. The pure Neumann problem in [Lecture 3](/blog/2026-10-13-numerical-pde-lecture-3/) breaks uniqueness, because we can add any constant to a solution. Continuous dependence came up in [Lecture 4](/blog/2026-10-15-numerical-pde-lecture-4/). Here we ask when an ODE has properties 1 and 2.

## 3. The Lipschitz condition

Take the standard form

$$
u'(t) = f\big(u(t), t\big), \qquad u(t_0) = \eta.
$$

Look at a region near the initial value:

$$
\mathcal{D} = \{(u, t) : |u - \eta| \le a,\ t_0 \le t \le t_1\}.
$$

We say $f$ is **Lipschitz continuous in $u$** on $\mathcal{D}$ if there is a constant $L \ge 0$ with

$$
|f(u, t) - f(u^*, t)| \le L\,|u - u^*| \qquad \text{for all } (u, t), (u^*, t) \in \mathcal{D}.
$$

Two remarks help me read this.

- The time $t$ is the **same** on both sides. Only the first argument changes. Here $u$ and $u^*$ are just numbers, not solutions.
- Divide by $|u - u^*|$. The left side becomes a difference quotient, a slope. So the condition says: **the slope of $f$ in $u$ is bounded by $L$**. The function can change, but it cannot suddenly shoot up.

**A bounded derivative gives $L$ at once.** If $f_u = \partial f / \partial u$ is bounded on $\mathcal{D}$, take

$$
L = \max_{(u, t) \in \mathcal{D}} |f_u(u, t)|.
$$

The proof is one line. By the mean value theorem (Taylor's theorem stopped at the first term),

$$
f(u, t) - f(u^*, t) = f_u(v, t)\,(u - u^*) \qquad \text{for some } v \text{ between } u \text{ and } u^*,
$$

and $|f_u(v, t)| \le L$. So the three conditions are ordered:

$$
\text{bounded derivative} \;\Longrightarrow\; \text{Lipschitz} \;\Longrightarrow\; \text{continuous}.
$$

Lipschitz continuity sits between continuity and differentiability. The function $|u|$ is an example: it is Lipschitz with $L = 1$, but it has no derivative at $0$.

## 4. Existence and uniqueness theorems

The lecture stated two theorems without proof.

| | Assumption | Conclusion |
|---|---|---|
| **Existence** (Peano) | $f$ is **continuous** on $\lvert u - \eta \rvert \le a$, $\lvert t - t_0 \rvert \le b$, and $\lvert f \rvert \le M$ there | **at least one** solution on $\lvert t - t_0 \rvert \le \alpha$, with $\alpha = \min\{b,\ a/M\}$ |
| **Uniqueness** (Picard–Lindelöf) | $f$ is **Lipschitz** in $u$ on $\mathcal{D}$ | **exactly one** solution, at least up to $T = \min\{t_1,\ t_0 + a/S\}$, with $S = \max_{\mathcal{D}} \lvert f \rvert$ |

A stronger assumption gives a stronger result: continuity gives existence, and Lipschitz continuity adds uniqueness.

Both results are **local**. They only promise a solution "at least up to some time". The time $a/M$ has a simple meaning. The solution leaves the region when it moves a distance $a$ away from $\eta$. Its speed $|u'| = |f|$ is at most $M$. So it needs at least the time $a/M$ to leave.

## 5. Two warnings

### Unique, but it blows up

Take $u' = u^2$ with $u(0) = \eta > 0$. Before solving, read the equation: $u' > 0$, so $u$ grows, and the growth rate is $u$ itself squared. Even $u' = u$ gives exponential growth, and this is faster. Separation of variables gives

$$
u(t) = \frac{1}{\eta^{-1} - t}.
$$

As $t \to 1/\eta$, the solution goes to infinity, and there is no solution after that. The function $f(u) = u^2$ has a bounded derivative on every bounded region, so it is Lipschitz there, and the solution is unique. But **unique does not mean it lives forever**. This is called **finite-time blow-up**.

Here is my own example, which I also use in the experiment: $u' = 1 + u^2$, $u(0) = 0$. The solution is $u = \tan t$, which blows up at $t = \pi/2 \approx 1.571$. What does the theorem promise? On $|u| \le a$ we have $S = 1 + a^2$, so it promises a solution up to $a/(1 + a^2)$. We are free to pick $a$. The best choice is $a = 1$, which gives $1/2$. The real solution lives about three times longer. The theorem says "at least", and it is conservative.

### A solution, but not only one

Take another example of my own: $u' = 3u^{2/3}$ with $u(0) = 0$. Near $u = 0$ the Lipschitz condition fails:

$$
\frac{|f(u) - f(0)|}{|u - 0|} = \frac{3u^{2/3}}{u} = 3u^{-1/3} \;\to\; \infty \quad (u \to 0^+).
$$

The slope is unbounded, so no $L$ works, and the uniqueness theorem does not apply. In fact there are two solutions:

$$
u(t) = 0 \qquad \text{and} \qquad u(t) = t^3.
$$

Separation of variables gives the second one: $u^{-2/3}\,du = 3\,dt$, so $3u^{1/3} = 3t$ and $u = t^3$. Check: $u' = 3t^2$ and $3u^{2/3} = 3t^2$. But the first step divides by $u^{2/3}$, and this loses the solution $u = 0$.

My own note: there are really infinitely many solutions. For any $c \ge 0$, the function that is $0$ up to $t = c$ and $(t - c)^3$ after that is also a solution. It waits as long as it likes, and then it starts to grow. So the Lipschitz condition is really needed for uniqueness.

## 6. Time marching

An initial value problem is different from the boundary value problems of Lectures 2 to 5. A boundary value problem has data at **both ends**, so we solved for all unknowns at once with one matrix. An initial value problem has all its data at the **start**. So we **march forward** one step at a time: from $U^0$ we get $U^1$, from $U^1$ we get $U^2$, and so on.

The time grid is

$$
t_0 < t_1 < \cdots < t_N = t_F, \qquad \Delta t = \frac{t_F - t_0}{N}, \qquad t_n = t_0 + n\,\Delta t,
$$

and $U^n \approx u(t_n)$ with $U^0 = \eta$. As in Lecture 3, capital letters are numerical values and small letters are the exact solution. The index is now a superscript, because later the subscript will be the space index.

## 7. Two Euler methods

### Forward Euler: explicit

Write the ODE at the time $t_n$, $u'(t_n) = f(u(t_n), t_n)$, and replace $u'(t_n)$ by a forward difference:

$$
\frac{U^{n+1} - U^n}{\Delta t} = f(U^n, t_n)
\qquad\Longrightarrow\qquad
U^{n+1} = U^n + \Delta t\, f(U^n, t_n).
$$

If we know $U^n$, we just evaluate the right side. There is no equation to solve. A method like this is called **explicit**. It is a simple recursion, like a sequence defined by a rule. Forward Euler is **first-order** accurate.

### Backward Euler: implicit

Now use a backward difference. Written at $t_{n+1}$ and shifted by one index, this is

$$
\frac{U^{n+1} - U^n}{\Delta t} = f(U^{n+1}, t_{n+1})
\qquad\Longrightarrow\qquad
U^{n+1} = U^n + \Delta t\, f(U^{n+1}, t_{n+1}).
$$

The written form with $(U^n - U^{n-1})/\Delta t = f(U^n, t_n)$ is the same scheme. Only the counting starts at a different index.

A good question for any scheme is: **where is $f$ evaluated?** This is the reference time of the scheme. Forward Euler stands at $t_n$ and looks ahead. Backward Euler stands at $t_{n+1}$ and looks back.

The unknown $U^{n+1}$ now appears on **both sides**. If $f$ is nonlinear, each step needs the solution of a nonlinear equation. A method like this is called **implicit**. We write each step as a root-finding problem:

$$
g(U) = U - \Delta t\, f(U, t_{n+1}) - U^n = 0,
$$

where $U^n$ is already known. Backward Euler is also first-order accurate.

## 8. Newton's method inside each step

To find a root of $g(x) = 0$, Newton's method draws the tangent line at the current guess $x^k$ and takes the point where it meets the axis:

$$
x^{k+1} = x^k - \frac{g(x^k)}{g'(x^k)}, \qquad k = 0, 1, 2, \dots
$$

For Backward Euler, the unknown is $U^{n+1}$, and $g'(U) = 1 - \Delta t\, f_u(U, t_{n+1})$:

$$
U^{n+1, k+1} = U^{n+1, k} - \frac{g(U^{n+1, k})}{g'(U^{n+1, k})}.
$$

There are two superscripts now, and it is easy to mix them up:

- $n + 1$ is the **time level**, the time $t_{n+1}$ we want.
- $k$ is the **Newton iteration count** at that one time level.

We repeat until the result converges. A computer cannot iterate forever, so in practice we stop when the change is below a **tolerance** that we choose. The last iterate is $U^{n+1}$. A natural first guess, and the one I use in the experiment, is the previous value $U^{n+1, 0} = U^n$.

Note that we run Newton's method again at **every** time step: once to get $U^1$, again to get $U^2$, and so on. That is expensive.

## 9. Cheap or expensive? Stability decides

| | Forward Euler | Backward Euler |
|---|---|---|
| cost of one step | **cheap**: only an evaluation | **expensive**: Newton's method |
| allowed $\Delta t$ | must be **small enough** | can be **large** |

So Forward Euler looks like the winner, but it is not that simple. If Forward Euler needs many more steps, the cheap steps add up, and Backward Euler can be faster overall. It depends on the problem. The reason why Forward Euler needs a small $\Delta t$ is **stability**. This is the central topic for the rest of the course. Part 5 of the experiment below gives a first look.

## 10. The trapezoidal method

Put Forward and Backward Euler side by side. The left sides are the same. Only the right sides differ: $f(U^n, t_n)$ or $f(U^{n+1}, t_{n+1})$. So why not take the average?

$$
\frac{U^{n+1} - U^n}{\Delta t} = \frac12\Big(f(U^n, t_n) + f(U^{n+1}, t_{n+1})\Big).
$$

This is the same idea as in space, where the balanced centered difference was second order and the one-sided difference was first order ([Lecture 2](/blog/2026-10-08-numerical-pde-lecture-2/)). The trapezoidal method is **second-order** accurate. It is still implicit, because $U^{n+1}$ is inside $f$, so it still needs Newton's method. In return it has good stability properties.

## 11. Taylor series methods

Expand the exact solution at $t_n$:

$$
u(t_{n+1}) = u(t_n) + \Delta t\, u'(t_n) + \frac{\Delta t^2}{2}\, u''(t_n) + \cdots
$$

The idea is to **use the ODE again and again**. We replace $u'$ by $f$. For $u''$ we differentiate $f(u(t), t)$ in time once more, with the chain rule:

$$
u'' = \frac{d}{dt} f\big(u(t), t\big) = f_u\,\frac{du}{dt} + f_t\,\frac{dt}{dt} = f_u\, f + f_t.
$$

Two points to watch: $d/dt$ is the **total** derivative, because the first argument $u$ also depends on $t$. And $dt/dt = 1$ is usually not written, but it is there. The **second-order Taylor method** is

$$
U^{n+1} = U^n + \Delta t\, f + \frac{\Delta t^2}{2}\big(f_u\, f + f_t\big), \qquad \text{all evaluated at } (U^n, t_n).
$$

It is explicit and second order. But the third-order method needs the time derivative of $f_u f + f_t$, and the product rule and the chain rule pile up. For a complicated $f$ this quickly becomes impractical. The idea is still good: use the differential equation itself to get higher derivatives.

### Is Forward Euler first or second order?

This confused me (and other students too). Forward Euler is the Taylor series cut after the $\Delta t$ term. The first dropped term is $\frac{\Delta t^2}{2} u''$. So is the method first order or second order? The answer: these are two different errors. Here is my own reading.

- **One-step (local) error.** Start one step from the exact value. The result misses the exact value by about $\frac{\Delta t^2}{2} u''$. This is $O(\Delta t^2)$.
- **Global error.** To reach a fixed time $T$ we take $N = T/\Delta t$ steps. Roughly, $N$ errors of size $O(\Delta t^2)$ add up to $O(\Delta t)$. So Forward Euler is **first order**.

The truncation error from Lecture 3 gives the same answer more cleanly. Put the exact solution into the difference form of the scheme:

$$
\tau^n = \frac{u(t_{n+1}) - u(t_n)}{\Delta t} - f\big(u(t_n), t_n\big) = \frac{\Delta t}{2}\, u''(t_n) + O(\Delta t^2).
$$

Dividing by $\Delta t$ lowers the order by one, and the order of $\tau$ matches the global order, as it did for boundary value problems. In the same way, the Taylor method that keeps the terms up to $\Delta t^p$ has a one-step error $O(\Delta t^{p+1})$ and is a method of order $p$. Part 3 of the experiment checks this.

## 12. One integral, three schemes

Integrate the ODE from $t_n$ to $t_{n+1}$. By the fundamental theorem of calculus,

$$
u(t_{n+1}) - u(t_n) = \int_{t_n}^{t_{n+1}} f\big(u(t), t\big)\, dt.
$$

Now every way to approximate this integral gives a scheme. These are the area rules from first-year calculus:

- **Left endpoint**: the integral $\approx \Delta t\, f(U^n, t_n)$. This is **Forward Euler**.
- **Right endpoint**: the integral $\approx \Delta t\, f(U^{n+1}, t_{n+1})$. This is **Backward Euler**.
- **Trapezoid rule**: the integral $\approx \frac{\Delta t}{2}\big(f(U^n, t_n) + f(U^{n+1}, t_{n+1})\big)$. This is the **trapezoidal method**, and its name comes from here.

This view explains a lot at once:

- The trapezoidal method is second order because the trapezoid rule is one order more accurate than a rectangle rule.
- Backward Euler is implicit because the height at the right end, $f(U^{n+1}, t_{n+1})$, is not known yet.
- The "differential equation = integral equation" idea from section 1 comes back as a numerical scheme: the choice of a quadrature rule is the choice of a scheme.

| Scheme | Update | Reference time | Explicit or implicit | Order | Integral rule |
|---|---|---|---|---|---|
| Forward Euler | $U^{n+1} = U^n + \Delta t\, f(U^n)$ | $t_n$ | explicit | 1 | left endpoint |
| Backward Euler | $U^{n+1} = U^n + \Delta t\, f(U^{n+1})$ | $t_{n+1}$ | implicit | 1 | right endpoint |
| Trapezoidal | $U^{n+1} = U^n + \frac{\Delta t}{2}\big(f(U^n) + f(U^{n+1})\big)$ | midway | implicit | 2 | trapezoid rule |
| Taylor (order $p$) | terms up to $\Delta t^p$ | $t_n$ | explicit | $p$ | (not practical) |

## 13. Experiment

The script uses my example from section 5:

$$
u' = 1 + u^2, \qquad u(0) = 0, \qquad u(t) = \tan t,
$$

and integrates it to $T = 1$, where $u(1) = \tan 1 \approx 1.5574$. Here $f_u = 2u$ and $f_t = 0$. It uses the Python standard library only. The implicit schemes use Newton's method with the first guess $U^n$, and stop when the Newton step is below $10^{-13}$ (relative).

**Part 1: the theorem's promise.** A search over $a$ confirms the best choice $a = 1$, which promises a solution up to $t = 0.5$. The real blow-up time is $\pi/2 \approx 1.5708$.

**Part 2: global error.** The chart shows $|U^N - u(1)|$ for $\Delta t = 1/10, 1/20, \dots, 1/10240$.

<figure class="viz">
<div class="viz-legend"><span class="viz-key"><span class="viz-swatch" style="background:#2a78d6"></span>Forward Euler</span><span class="viz-key"><span class="viz-swatch" style="background:#eb6834"></span>Backward Euler</span><span class="viz-key"><span class="viz-swatch" style="background:#1baf7a"></span>Trapezoidal</span><span class="viz-key"><span class="viz-swatch" style="background:#eda100"></span>Taylor 2</span></div>
<div class="viz-frame">
<svg viewBox="0 0 680 400" role="img" tabindex="0" aria-label="Error at t = 1 for u prime = 1 + u squared, as the time step gets smaller, on log scales. Forward and Backward Euler fall with slope 1 and almost lie on top of each other. The trapezoidal method and the second-order Taylor method fall with slope 2, and the trapezoidal error is about half the Taylor error." class="viz-plot"><line x1="64" x2="520" y1="30.0" y2="30.0" class="viz-grid"/><text x="54" y="34.0" text-anchor="end" class="viz-tick">1</text><line x1="64" x2="520" y1="107.5" y2="107.5" class="viz-grid"/><text x="54" y="111.5" text-anchor="end" class="viz-tick">10<tspan dy="-0.5em" font-size="75%">−2</tspan></text><line x1="64" x2="520" y1="185.0" y2="185.0" class="viz-grid"/><text x="54" y="189.0" text-anchor="end" class="viz-tick">10<tspan dy="-0.5em" font-size="75%">−4</tspan></text><line x1="64" x2="520" y1="262.5" y2="262.5" class="viz-grid"/><text x="54" y="266.5" text-anchor="end" class="viz-tick">10<tspan dy="-0.5em" font-size="75%">−6</tspan></text><line x1="64" x2="520" y1="340.0" y2="340.0" class="viz-grid"/><text x="54" y="344.0" text-anchor="end" class="viz-tick">10<tspan dy="-0.5em" font-size="75%">−8</tspan></text><text x="64.0" y="362" text-anchor="middle" class="viz-tick">1/10</text><text x="155.2" y="362" text-anchor="middle" class="viz-tick">1/40</text><text x="246.4" y="362" text-anchor="middle" class="viz-tick">1/160</text><text x="337.6" y="362" text-anchor="middle" class="viz-tick">1/640</text><text x="428.8" y="362" text-anchor="middle" class="viz-tick">1/2560</text><text x="520.0" y="362" text-anchor="middle" class="viz-tick">1/10240</text><text x="292.0" y="388" text-anchor="middle" class="viz-axis-title">Time step Δt (smaller to the right)</text><text x="64" y="14" class="viz-axis-title">Error</text><path d="M64.0,60.7 L109.6,70.3 L155.2,80.8 L200.8,91.9 L246.4,103.2 L292.0,114.7 L337.6,126.3 L383.2,137.9 L428.8,149.5 L474.4,161.2 L520.0,172.8" fill="none" stroke="#2a78d6" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><circle cx="64.0" cy="60.7" r="3" fill="#2a78d6" stroke="#ffffff" stroke-width="1.5"/><circle cx="109.6" cy="70.3" r="3" fill="#2a78d6" stroke="#ffffff" stroke-width="1.5"/><circle cx="155.2" cy="80.8" r="3" fill="#2a78d6" stroke="#ffffff" stroke-width="1.5"/><circle cx="200.8" cy="91.9" r="3" fill="#2a78d6" stroke="#ffffff" stroke-width="1.5"/><circle cx="246.4" cy="103.2" r="3" fill="#2a78d6" stroke="#ffffff" stroke-width="1.5"/><circle cx="292.0" cy="114.7" r="3" fill="#2a78d6" stroke="#ffffff" stroke-width="1.5"/><circle cx="337.6" cy="126.3" r="3" fill="#2a78d6" stroke="#ffffff" stroke-width="1.5"/><circle cx="383.2" cy="137.9" r="3" fill="#2a78d6" stroke="#ffffff" stroke-width="1.5"/><circle cx="428.8" cy="149.5" r="3" fill="#2a78d6" stroke="#ffffff" stroke-width="1.5"/><circle cx="474.4" cy="161.2" r="3" fill="#2a78d6" stroke="#ffffff" stroke-width="1.5"/><circle cx="520.0" cy="172.8" r="3" fill="#2a78d6" stroke="#ffffff" stroke-width="1.5"/><path d="M64.0,48.8 L109.6,64.8 L155.2,78.1 L200.8,90.5 L246.4,102.5 L292.0,114.3 L337.6,126.1 L383.2,137.8 L428.8,149.5 L474.4,161.2 L520.0,172.8" fill="none" stroke="#eb6834" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><circle cx="64.0" cy="48.8" r="3" fill="#eb6834" stroke="#ffffff" stroke-width="1.5"/><circle cx="109.6" cy="64.8" r="3" fill="#eb6834" stroke="#ffffff" stroke-width="1.5"/><circle cx="155.2" cy="78.1" r="3" fill="#eb6834" stroke="#ffffff" stroke-width="1.5"/><circle cx="200.8" cy="90.5" r="3" fill="#eb6834" stroke="#ffffff" stroke-width="1.5"/><circle cx="246.4" cy="102.5" r="3" fill="#eb6834" stroke="#ffffff" stroke-width="1.5"/><circle cx="292.0" cy="114.3" r="3" fill="#eb6834" stroke="#ffffff" stroke-width="1.5"/><circle cx="337.6" cy="126.1" r="3" fill="#eb6834" stroke="#ffffff" stroke-width="1.5"/><circle cx="383.2" cy="137.8" r="3" fill="#eb6834" stroke="#ffffff" stroke-width="1.5"/><circle cx="428.8" cy="149.5" r="3" fill="#eb6834" stroke="#ffffff" stroke-width="1.5"/><circle cx="474.4" cy="161.2" r="3" fill="#eb6834" stroke="#ffffff" stroke-width="1.5"/><circle cx="520.0" cy="172.8" r="3" fill="#eb6834" stroke="#ffffff" stroke-width="1.5"/><path d="M64.0,99.9 L109.6,123.6 L155.2,147.0 L200.8,170.4 L246.4,193.7 L292.0,217.0 L337.6,240.4 L383.2,263.7 L428.8,287.0 L474.4,310.4 L520.0,333.7" fill="none" stroke="#1baf7a" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><circle cx="64.0" cy="99.9" r="3" fill="#1baf7a" stroke="#ffffff" stroke-width="1.5"/><circle cx="109.6" cy="123.6" r="3" fill="#1baf7a" stroke="#ffffff" stroke-width="1.5"/><circle cx="155.2" cy="147.0" r="3" fill="#1baf7a" stroke="#ffffff" stroke-width="1.5"/><circle cx="200.8" cy="170.4" r="3" fill="#1baf7a" stroke="#ffffff" stroke-width="1.5"/><circle cx="246.4" cy="193.7" r="3" fill="#1baf7a" stroke="#ffffff" stroke-width="1.5"/><circle cx="292.0" cy="217.0" r="3" fill="#1baf7a" stroke="#ffffff" stroke-width="1.5"/><circle cx="337.6" cy="240.4" r="3" fill="#1baf7a" stroke="#ffffff" stroke-width="1.5"/><circle cx="383.2" cy="263.7" r="3" fill="#1baf7a" stroke="#ffffff" stroke-width="1.5"/><circle cx="428.8" cy="287.0" r="3" fill="#1baf7a" stroke="#ffffff" stroke-width="1.5"/><circle cx="474.4" cy="310.4" r="3" fill="#1baf7a" stroke="#ffffff" stroke-width="1.5"/><circle cx="520.0" cy="333.7" r="3" fill="#1baf7a" stroke="#ffffff" stroke-width="1.5"/><path d="M64.0,92.4 L109.6,113.9 L155.2,136.3 L200.8,159.1 L246.4,182.3 L292.0,205.5 L337.6,228.8 L383.2,252.1 L428.8,275.4 L474.4,298.7 L520.0,322.0" fill="none" stroke="#eda100" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><circle cx="64.0" cy="92.4" r="3" fill="#eda100" stroke="#ffffff" stroke-width="1.5"/><circle cx="109.6" cy="113.9" r="3" fill="#eda100" stroke="#ffffff" stroke-width="1.5"/><circle cx="155.2" cy="136.3" r="3" fill="#eda100" stroke="#ffffff" stroke-width="1.5"/><circle cx="200.8" cy="159.1" r="3" fill="#eda100" stroke="#ffffff" stroke-width="1.5"/><circle cx="246.4" cy="182.3" r="3" fill="#eda100" stroke="#ffffff" stroke-width="1.5"/><circle cx="292.0" cy="205.5" r="3" fill="#eda100" stroke="#ffffff" stroke-width="1.5"/><circle cx="337.6" cy="228.8" r="3" fill="#eda100" stroke="#ffffff" stroke-width="1.5"/><circle cx="383.2" cy="252.1" r="3" fill="#eda100" stroke="#ffffff" stroke-width="1.5"/><circle cx="428.8" cy="275.4" r="3" fill="#eda100" stroke="#ffffff" stroke-width="1.5"/><circle cx="474.4" cy="298.7" r="3" fill="#eda100" stroke="#ffffff" stroke-width="1.5"/><circle cx="520.0" cy="322.0" r="3" fill="#eda100" stroke="#ffffff" stroke-width="1.5"/><text x="530" y="166.8" class="viz-label">Forward Euler</text><text x="530" y="186.8" class="viz-label">Backward Euler</text><text x="530" y="347.7" class="viz-label">Trapezoidal</text><text x="530" y="316.0" class="viz-label">Taylor 2</text><line x1="0" x2="0" y1="30" y2="340" class="viz-crosshair" visibility="hidden"/></svg>
<div class="viz-tooltip" hidden></div>
</div>
<figcaption>Error |U<sup>N</sup> − u(1)| on log scales for u′ = 1 + u<sup>2</sup>, u(0) = 0, exact solution u = tan t. Forward and Backward Euler almost overlap. Hover over the chart, or focus it and use the arrow keys, to read the values.</figcaption>
<details class="viz-table"><summary>Show the numbers</summary><table><thead><tr><th>Δt</th><th>Forward Euler</th><th>Backward Euler</th><th>Trapezoidal</th><th>Taylor 2</th></tr></thead><tbody><tr><td>1/10</td><td>1.6e-1</td><td>3.3e-1</td><td>1.6e-2</td><td>2.5e-2</td></tr><tr><td>1/20</td><td>9.1e-2</td><td>1.3e-1</td><td>3.8e-3</td><td>6.9e-3</td></tr><tr><td>1/40</td><td>4.9e-2</td><td>5.7e-2</td><td>9.6e-4</td><td>1.8e-3</td></tr><tr><td>1/80</td><td>2.5e-2</td><td>2.7e-2</td><td>2.4e-4</td><td>4.6e-4</td></tr><tr><td>1/160</td><td>1.3e-2</td><td>1.3e-2</td><td>6.0e-5</td><td>1.2e-4</td></tr><tr><td>1/320</td><td>6.5e-3</td><td>6.7e-3</td><td>1.5e-5</td><td>3.0e-5</td></tr><tr><td>1/640</td><td>3.3e-3</td><td>3.3e-3</td><td>3.7e-6</td><td>7.4e-6</td></tr><tr><td>1/1280</td><td>1.6e-3</td><td>1.7e-3</td><td>9.3e-7</td><td>1.9e-6</td></tr><tr><td>1/2560</td><td>8.2e-4</td><td>8.2e-4</td><td>2.3e-7</td><td>4.7e-7</td></tr><tr><td>1/5120</td><td>4.1e-4</td><td>4.1e-4</td><td>5.8e-8</td><td>1.2e-7</td></tr><tr><td>1/10240</td><td>2.1e-4</td><td>2.1e-4</td><td>1.5e-8</td><td>2.9e-8</td></tr></tbody></table></details>
<script>(() => {const fig = document.currentScript.closest('figure');const svg = fig.querySelector('svg');const line = svg.querySelector('.viz-crosshair');const tip = fig.querySelector('.viz-tooltip');const names = ["Forward Euler", "Backward Euler", "Trapezoidal", "Taylor 2"]; const colors = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"];const rows = [["1.6e-1", "3.3e-1", "1.6e-2", "2.5e-2"], ["9.1e-2", "1.3e-1", "3.8e-3", "6.9e-3"], ["4.9e-2", "5.7e-2", "9.6e-4", "1.8e-3"], ["2.5e-2", "2.7e-2", "2.4e-4", "4.6e-4"], ["1.3e-2", "1.3e-2", "6.0e-5", "1.2e-4"], ["6.5e-3", "6.7e-3", "1.5e-5", "3.0e-5"], ["3.3e-3", "3.3e-3", "3.7e-6", "7.4e-6"], ["1.6e-3", "1.7e-3", "9.3e-7", "1.9e-6"], ["8.2e-4", "8.2e-4", "2.3e-7", "4.7e-7"], ["4.1e-4", "4.1e-4", "5.8e-8", "1.2e-7"], ["2.1e-4", "2.1e-4", "1.5e-8", "2.9e-8"]]; const heads = ["\u0394t = 1/10", "\u0394t = 1/20", "\u0394t = 1/40", "\u0394t = 1/80", "\u0394t = 1/160", "\u0394t = 1/320", "\u0394t = 1/640", "\u0394t = 1/1280", "\u0394t = 1/2560", "\u0394t = 1/5120", "\u0394t = 1/10240"];const left = 64, right = 520, steps = 10, width = 680;let current = -1;const show = (n) => { n = Math.max(0, Math.min(steps, n)); current = n; const x = left + (right - left) * n / steps; line.setAttribute('x1', x); line.setAttribute('x2', x); line.setAttribute('visibility', 'visible'); tip.replaceChildren(); const head = document.createElement('div'); head.className = 'viz-tooltip-head'; head.textContent = heads[n]; tip.append(head); names.forEach((name, i) => {  const row = document.createElement('div'); row.className = 'viz-tooltip-row';  const key = document.createElement('span'); key.className = 'viz-swatch'; key.style.background = colors[i];  const value = document.createElement('strong'); value.textContent = rows[n][i];  const label = document.createElement('span'); label.textContent = name;  row.append(key, value, label); tip.append(row); }); const box = svg.getBoundingClientRect(), frame = fig.getBoundingClientRect(); const xs = box.left - frame.left + x * box.width / width; tip.hidden = false; const room = frame.width - tip.offsetWidth - 8; tip.style.left = Math.max(0, Math.min(room, xs + 12)) + 'px'; tip.style.top = (box.top - frame.top + 8) + 'px';};const hide = () => { line.setAttribute('visibility', 'hidden'); tip.hidden = true; current = -1; };svg.addEventListener('pointermove', (event) => { const box = svg.getBoundingClientRect(); const x = (event.clientX - box.left) * width / box.width; show(Math.round((x - left) / (right - left) * steps));});svg.addEventListener('pointerleave', hide);svg.addEventListener('focus', () => show(current < 0 ? 0 : current));svg.addEventListener('blur', hide);svg.addEventListener('keydown', (event) => { if (event.key === 'ArrowRight') { show(current + 1); event.preventDefault(); } if (event.key === 'ArrowLeft') { show(current - 1); event.preventDefault(); } if (event.key === 'Escape') hide();});})();</script>
</figure>

| $\Delta t$ | Forward Euler | order | Backward Euler | order | Trapezoidal | order | Taylor 2 | order |
|---|---|---|---|---|---|---|---|---|
| 1/10 | $1.61 \times 10^{-1}$ | — | $3.26 \times 10^{-1}$ | — | $1.57 \times 10^{-2}$ | — | $2.45 \times 10^{-2}$ | — |
| 1/20 | $9.11 \times 10^{-2}$ | 0.82 | $1.27 \times 10^{-1}$ | 1.36 | $3.84 \times 10^{-3}$ | 2.03 | $6.86 \times 10^{-3}$ | 1.84 |
| 1/40 | $4.88 \times 10^{-2}$ | 0.90 | $5.74 \times 10^{-2}$ | 1.14 | $9.55 \times 10^{-4}$ | 2.01 | $1.81 \times 10^{-3}$ | 1.92 |
| 1/80 | $2.53 \times 10^{-2}$ | 0.95 | $2.75 \times 10^{-2}$ | 1.06 | $2.38 \times 10^{-4}$ | 2.00 | $4.65 \times 10^{-4}$ | 1.96 |
| 1/160 | $1.29 \times 10^{-2}$ | 0.97 | $1.35 \times 10^{-2}$ | 1.03 | $5.96 \times 10^{-5}$ | 2.00 | $1.18 \times 10^{-4}$ | 1.98 |
| 1/320 | $6.52 \times 10^{-3}$ | 0.99 | $6.66 \times 10^{-3}$ | 1.01 | $1.49 \times 10^{-5}$ | 2.00 | $2.96 \times 10^{-5}$ | 1.99 |
| 1/640 | $3.28 \times 10^{-3}$ | 0.99 | $3.31 \times 10^{-3}$ | 1.01 | $3.72 \times 10^{-6}$ | 2.00 | $7.43 \times 10^{-6}$ | 2.00 |
| 1/1280 | $1.64 \times 10^{-3}$ | 1.00 | $1.65 \times 10^{-3}$ | 1.00 | $9.31 \times 10^{-7}$ | 2.00 | $1.86 \times 10^{-6}$ | 2.00 |
| 1/2560 | $8.23 \times 10^{-4}$ | 1.00 | $8.25 \times 10^{-4}$ | 1.00 | $2.33 \times 10^{-7}$ | 2.00 | $4.65 \times 10^{-7}$ | 2.00 |
| 1/5120 | $4.12 \times 10^{-4}$ | 1.00 | $4.12 \times 10^{-4}$ | 1.00 | $5.82 \times 10^{-8}$ | 2.00 | $1.16 \times 10^{-7}$ | 2.00 |
| 1/10240 | $2.06 \times 10^{-4}$ | 1.00 | $2.06 \times 10^{-4}$ | 1.00 | $1.45 \times 10^{-8}$ | 2.00 | $2.91 \times 10^{-8}$ | 2.00 |

What I noticed:

- **The orders match the theory.** Forward and Backward Euler settle at $1.00$, the trapezoidal and Taylor methods at $2.00$. Backward Euler reaches its order from above ($1.36$ first), and Forward Euler from below ($0.82$ first). Both start in the pre-asymptotic range from [Lecture 5](/blog/2026-10-20-numerical-pde-lecture-5/).
- **The two Euler errors have opposite signs.** At $\Delta t = 1/160$, Forward Euler gives $U - u = -1.29 \times 10^{-2}$ (too low) and Backward Euler gives $+1.35 \times 10^{-2}$ (too high). This fits the integral view. The function $f(u(t)) = 1 + \tan^2 t$ grows, so the left endpoint rule underestimates the area and the right endpoint rule overestimates it. The trapezoid sits between them.
- **Two second-order methods, two constants.** At the smallest step the trapezoidal error is $1.45 \times 10^{-8}$ and the Taylor error is $2.91 \times 10^{-8}$. The slope is the same, but the constant $C$ in $E \approx C\,\Delta t^2$ differs by about a factor of two.

**Part 3: one-step error against global error (Forward Euler).** Start one step from the exact value $u(0.5) = \tan 0.5$ and compare with $\tan(0.5 + \Delta t)$:

| $\Delta t$ | one-step error | order | ratio to $\frac{\Delta t^2}{2}u''$ | global error at $T = 1$ | order |
|---|---|---|---|---|---|
| 1/10 | $7.99 \times 10^{-3}$ | — | 1.1263 | $1.61 \times 10^{-1}$ | — |
| 1/20 | $1.88 \times 10^{-3}$ | 2.09 | 1.0604 | $9.11 \times 10^{-2}$ | 0.82 |
| 1/40 | $4.56 \times 10^{-4}$ | 2.04 | 1.0295 | $4.88 \times 10^{-2}$ | 0.90 |
| 1/80 | $1.12 \times 10^{-4}$ | 2.02 | 1.0146 | $2.53 \times 10^{-2}$ | 0.95 |
| 1/160 | $2.79 \times 10^{-5}$ | 2.01 | 1.0073 | $1.29 \times 10^{-2}$ | 0.97 |
| 1/320 | $6.95 \times 10^{-6}$ | 2.01 | 1.0036 | $6.52 \times 10^{-3}$ | 0.99 |
| 1/640 | $1.73 \times 10^{-6}$ | 2.00 | 1.0018 | $3.28 \times 10^{-3}$ | 0.99 |
| 1/1280 | $4.33 \times 10^{-7}$ | 2.00 | 1.0009 | $1.64 \times 10^{-3}$ | 1.00 |
| 1/2560 | $1.08 \times 10^{-7}$ | 2.00 | 1.0005 | $8.23 \times 10^{-4}$ | 1.00 |
| 1/5120 | $2.71 \times 10^{-8}$ | 2.00 | 1.0002 | $4.12 \times 10^{-4}$ | 1.00 |
| 1/10240 | $6.77 \times 10^{-9}$ | 2.00 | 1.0001 | $2.06 \times 10^{-4}$ | 1.00 |

The one-step error has order $2.00$, and its ratio to $\frac{\Delta t^2}{2} u''(0.5)$ goes to $1$ (here $u'' = f_u f = 2u(1 + u^2) \approx 1.4187$). The global error has order $1.00$. So both statements in section 11 are true: the one-step error is second order, and the method is first order.

**Part 4: the price of an implicit step.** With the first guess $U^n$, Newton's method needs 4 to 5 iterations per step at $\Delta t = 1/10$, and exactly 3 per step at $\Delta t = 1/160$ and $1/2560$, for both Backward Euler and the trapezoidal method. A smaller step makes the first guess better, so Newton needs fewer iterations. But each implicit step still costs a few evaluations of $f$ and $f_u$, while an explicit step costs one.

**Part 5: a first look at stability.** This part uses a different, linear problem (my own):

$$
u' = -50\,\big(u - \cos t\big) - \sin t, \qquad u(0) = 1, \qquad u(t) = \cos t.
$$

The solution is smooth and slow. But any deviation from $\cos t$ decays very fast, like $e^{-50 t}$. Here is the error at $T = 2$:

| $\Delta t$ | Forward Euler | Backward Euler |
|---|---|---|
| 0.01 | $4.01 \times 10^{-5}$ | $3.95 \times 10^{-5}$ |
| 0.025 | $1.01 \times 10^{-4}$ | $9.75 \times 10^{-5}$ |
| 0.05 | $5.53 \times 10^{3}$ | $1.91 \times 10^{-4}$ |

For the two smaller steps, both methods are fine and have about the same error. With $\Delta t = 0.05$, Backward Euler is still accurate, but Forward Euler explodes to an error of about $5.5 \times 10^{3}$, even though the exact solution stays between $-1$ and $1$. The step is not too large for **accuracy**: the solution $\cos t$ barely changes over $0.05$. It is too large for **stability**. Why this happens, and how small $\Delta t$ must be, is the subject of the next lectures.

The full script and the chart code are in [this site's repository](https://github.com/doubleyands/doubleyands.github.io/tree/main/scripts/posts/numerical-pde/lecture-6).

## Summary

1. An integrating factor gives **Duhamel's formula**. If the source depends on $u$, the formula becomes an **integral equation**.
2. A problem is **well-posed** if a solution exists, is unique, and depends continuously on the data.
3. $f$ is **Lipschitz** in $u$ if its slope in $u$ is bounded by some $L$. A bounded derivative gives $L = \max |f_u|$.
4. Continuity of $f$ gives **existence**, and the Lipschitz condition adds **uniqueness**. Both results are local.
5. Uniqueness does not prevent **blow-up** ($u' = 1 + u^2$ gives $\tan t$). Without the Lipschitz condition, uniqueness can fail ($u' = 3u^{2/3}$ has $0$ and $t^3$).
6. **Forward Euler** is explicit and first order. **Backward Euler** is implicit and first order, and needs Newton's method at every step. The **trapezoidal method** is implicit and second order. **Taylor methods** reach any order but become impractical.
7. A one-step error of $O(\Delta t^{p+1})$ gives a method of order $p$.
8. All these schemes come from approximating $\int_{t_n}^{t_{n+1}} f\,dt$ with an area rule.
9. Cheap explicit steps can lose to expensive implicit steps when **stability** forces a small $\Delta t$.
