---
title: "Numerical Analysis, Week 5: Euler's Method, Multistep Methods, and the Trapezoidal Method"
description: "Existence, uniqueness and stability of initial value problems, a stable but ill-conditioned example, Euler's method with a full convergence proof, multistep methods, the midpoint and trapezoidal methods, and Picard iteration. Tested on the logistic equation."
pubDate: 2026-10-31
project: numerical-analysis
tags: [numerical analysis, ODE, euler method, multistep methods, trapezoidal method, notes]
---

This is my summary of the fifth week of a numerical analysis course. As before, I write the ideas in my own words and check some of them with a small experiment. The course follows K. E. Atkinson, *An Introduction to Numerical Analysis*.

Last week ended with the start of a new chapter: initial value problems (see [Week 4](/blog/2026-10-24-numerical-analysis-week-4/)). The problem is

$$
y'(x) = f(x, y(x)), \qquad y(x_0) = y_0. \tag{1}
$$

This week has two parts. First, the theory: when does (1) have a solution, is it the only one, and does it change a lot when the data changes a little? Second, the numerical methods: Euler's method, its convergence proof, and two methods of second order.

A word on names. The equation is *ordinary* because there is only one independent variable, $x$. With two or more independent variables and their partial derivatives, it would be a partial differential equation. The condition $y(x_0) = y_0$ is called an *initial* condition even when $x$ is not time.

## 1. Existence and uniqueness

**Theorem.** Let $D = \{(x, y) : |x - x_0| < \alpha,\ |y - y_0| < \beta\}$. Assume that $f(x, y)$ is continuous in $x$ and $y$ on $D$, and that for some $K \ge 0$

$$
|f(x, y_1) - f(x, y_2)| \le K\,|y_1 - y_2| \qquad \text{for all } (x, y_1), (x, y_2) \in D.
$$

Then (1) has a unique solution $y(x)$ on $[x_0 - \gamma,\ x_0 + \gamma]$, where

$$
\gamma = \min\left\{\alpha,\ \frac{\beta}{M}\right\}, \qquad M = \max_{(x, y)\in D} |f(x, y)|.
$$

The inequality is the **Lipschitz condition**, and $K$ is the **Lipschitz constant**. Two points are easy to miss:

- There are two conditions. $f$ must be continuous in both $x$ and $y$. *And* $f$ must be Lipschitz in $y$. The Lipschitz condition is about $y$ only, not $x$.
- The result is **local**. It gives a solution near $x_0$, on an interval of half-width $\gamma$. It does not promise a solution for all $x$.

The theorem is often called the Picard–Lindelöf theorem. Many ODE books start with it.

**Why $\gamma = \min\{\alpha, \beta/M\}$?** The slope of the solution is $|y'| = |f| \le M$. So while $x$ moves by $\beta/M$, $y$ moves by at most $\beta$. This means the solution cannot leave the box $D$ through the top or the bottom before $x$ moves by $\beta/M$. The box also has half-width $\alpha$ in $x$. The smaller of the two is $\gamma$.

**A sufficient condition, not a necessary one.** The theorem says: if the conditions hold, then a unique solution exists. It does not say the opposite. If a condition fails, a solution may still exist. The theorem just cannot promise it. In logic: $A \Rightarrow B$ does not give $\text{not } A \Rightarrow \text{not } B$.

### Checking the Lipschitz condition

If $\partial f / \partial y$ exists, the mean value theorem gives an easy check. For $y_1 \ne y_2$ there is a $\xi$ between them with

$$
f(x, y_1) - f(x, y_2) = \frac{\partial f}{\partial y}(x, \xi)\,(y_1 - y_2).
$$

So we can take $K = \max_D |\partial f / \partial y|$.

**Example.** Take $f(x, y) = x \arctan y$ on the strip $0 \le x \le 2$, $-\infty < y < \infty$. It is continuous. And

$$
\frac{\partial f}{\partial y} = \frac{x}{1 + y^2}, \qquad \left|\frac{\partial f}{\partial y}\right| \le 2,
$$

so $K = 2$. Also $|f| \le 2 \cdot \frac{\pi}{2} = \pi$, so $M = \pi$ and $\gamma = \min\{\alpha, \beta/\pi\}$. The values of $\alpha$ and $\beta$ come from the initial point. Once $(x_0, y_0)$ is given, the box and $\gamma$ follow.

Here the strip has no top or bottom, so $\beta$ can be as large as we like. Then $\beta/\pi$ is large too, and only $\alpha$ limits $\gamma$. Compare this with $y' = -y^2$ from Week 4. There, $\partial f/\partial y = -2y$ grows with $|y|$, so no single $K$ works on a whole strip. And indeed that solution blows up in finite time. "Only local" is a real limit, not just a careful word.

## 2. Stability of the problem

Why do we care about existence and uniqueness? If we know a solution exists in some region, we know where to look for an approximation. If uniqueness fails, a numerical method can give strange, oscillating answers.

The third property is **stability**. A problem is stable if a small change in the data gives a small change in the solution. Consider the perturbed problem

$$
y' = f(x, y) + \delta(x), \qquad y(x_0) = y_0 + \epsilon, \tag{2}
$$

where $\delta$ is continuous.

**Theorem.** Under the same assumptions as the existence theorem, take $\epsilon_0$ small enough. Then for all perturbations with $|\epsilon| \le \epsilon_0$ and $\|\delta\|_\infty \le \epsilon_0$, (2) has a unique solution $y(x; \delta, \epsilon)$ on $[x_0 - \gamma, x_0 + \gamma]$, and

$$
\max_{|x - x_0| \le \gamma} |y(x) - y(x; \delta, \epsilon)| \le k\bigl(|\epsilon| + \gamma\,\|\delta\|_\infty\bigr), \qquad k = \frac{1}{1 - \gamma K}.
$$

In words: the solution depends continuously on the data. If we change the right-hand side or the initial value by a tiny amount, the solution moves by a tiny amount times a constant.

### Three different things

Three ideas sound alike, and it is important to keep them apart:

1. **Stability of the problem.** This is a property of the ODE itself. It is what the theorem above is about.
2. **Stability of a numerical method.** This is a property of a method applied to the problem. We will see one for Euler's method below.
3. **Ill-conditioning.** When we discretize a differential equation, we often get an algebraic system $A\mathbf{x} = \mathbf{b}$. If the condition number of $A$ is large, small errors in $\mathbf{b}$ give large errors in $\mathbf{x}$. A problem can be stable and still lead to an ill-conditioned matrix. Condition numbers come in a later chapter.

### A stable problem that is ill-conditioned

Here is my own example:

$$
y' = 50\,\bigl(y - \cos x\bigr) - \sin x, \qquad y(0) = 1.
$$

Check: $y = \cos x$ gives $y' = -\sin x$ on the left, and $50 \cdot 0 - \sin x$ on the right. So $y(x) = \cos x$ is the solution. Now change the initial value to $1 + \epsilon$. The difference $w = y - \cos x$ solves $w' = 50w$, $w(0) = \epsilon$. So

$$
y(x; \epsilon) = \cos x + \epsilon\,e^{50x}.
$$

On $[0, 1]$ the gap is at most $\epsilon\,e^{50}$. So a constant exists, and in the sense of the theorem the problem is stable. But the constant is $e^{50} \approx 5 \times 10^{21}$. With $\epsilon = 10^{-16}$, about the size of one rounding error in double precision, the gap at $x = 1$ is already $5 \times 10^{5}$. The true solution $\cos x$ is less than $1$ in size. Any method will meet rounding errors, and each one acts like an $\epsilon$. So no method can follow $\cos x$ to $x = 1$ in double precision. A problem like this is called **ill-conditioned**. If we let $x$ go to infinity, no fixed constant works at all.

## 3. Higher order equations become first order systems

The chapter only solves first order equations. This is enough. A higher order equation can always be rewritten as a first order system.

A system of $m$ first order equations is

$$
y_i' = f_i(x, y_1, \dots, y_m), \quad y_i(x_0) = y_{i,0} \quad (i = 1, \dots, m),
$$

or in vector form, $\mathbf{y}' = \mathbf{f}(x, \mathbf{y})$, $\mathbf{y}(x_0) = \mathbf{y}_0$. Now take an equation of order $m$:

$$
y^{(m)} = f\bigl(x, y, y', \dots, y^{(m-1)}\bigr), \qquad y(x_0) = y_0,\ y'(x_0) = y_{0,1},\ \dots,\ y^{(m-1)}(x_0) = y_{0,m-1}.
$$

Set $y_1 = y$, $y_2 = y'$, ..., $y_m = y^{(m-1)}$. Then

$$
y_1' = y_2, \quad y_2' = y_3, \quad \dots, \quad y_{m-1}' = y_m, \quad y_m' = f(x, y_1, \dots, y_m),
$$

with $y_1(x_0) = y_0$, $y_2(x_0) = y_{0,1}$, and so on. An equation of order $m$ gives $m$ unknowns and needs $m$ initial conditions.

**Example.** $y'' + 2y' + 5y = \cos x$, with $y(0) = 1$, $y'(0) = 0$. Write it as $y'' = -5y - 2y' + \cos x$. With $y_1 = y$ and $y_2 = y'$:

$$
\begin{bmatrix} y_1 \\ y_2 \end{bmatrix}' =
\begin{bmatrix} 0 & 1 \\ -5 & -2 \end{bmatrix}
\begin{bmatrix} y_1 \\ y_2 \end{bmatrix} +
\begin{bmatrix} 0 \\ \cos x \end{bmatrix},
\qquad
\begin{bmatrix} y_1(0) \\ y_2(0) \end{bmatrix} = \begin{bmatrix} 1 \\ 0 \end{bmatrix}.
$$

So a good method for first order equations handles every higher order equation too. All methods below work for systems without change: just read $y$ and $f$ as vectors.

There are many ways to discretize a differential equation. For an ODE in one variable, the **finite difference method** is the usual first choice. For a time-dependent PDE, we may choose finite volumes or finite elements in space, but the time direction is still usually done with finite differences. Time is continuous, and it only marches forward. So we build the solution step by step.

## 4. Euler's method

Use a uniform grid $x_j = x_0 + jh$. (A uniform grid is not required, but it is easier to analyze.) Write $y_n$ for the approximation of $y(x_n)$.

**Euler's method:**

$$
y_{n+1} = y_n + h\,f(x_n, y_n), \qquad n = 0, 1, 2, \dots
$$

It is simple. Start from $y_0$, compute $y_1$, then $y_2$, and so on. There are three ways to derive it.

**1. Taylor expansion.** Expand $y(x_{n+1})$ around $x_n$:

$$
y(x_{n+1}) = y(x_n) + h\,y'(x_n) + \frac{h^2}{2}\,y''(\xi_n), \qquad x_n \le \xi_n \le x_{n+1}.
$$

Put $y'(x_n) = f(x_n, y(x_n))$. Then

$$
y(x_{n+1}) = y(x_n) + h\,f(x_n, y(x_n)) + \frac{h^2}{2}\,y''(\xi_n).
$$

We do not know the exact $y$, so we replace $y(x_n)$ by $y_n$ and drop the last term. This is Euler's method. To be exact, it is **forward Euler**: the new value comes from information at the old point.

**2. Numerical differentiation.** Solve the Taylor expansion for the derivative:

$$
y'(x_n) = \frac{y(x_{n+1}) - y(x_n)}{h} - \frac{h}{2}\,y''(\xi_n).
$$

The fraction replaces the derivative. The part we cut off, $-\frac{h}{2} y''(\xi_n)$, is the **truncation error**. The name comes from this: we truncate the formula for the derivative.

**3. Numerical integration.** Integrate $y' = f$ over one step:

$$
y(x_{n+1}) - y(x_n) = \int_{x_n}^{x_{n+1}} f(x, y(x))\,dx.
$$

Approximate the integral by one rectangle. If its height is the value at the left end, $f(x_n, y(x_n))$, we get forward Euler. If its height is the value at the right end, we get

$$
y_{n+1} = y_n + h\,f(x_{n+1}, y_{n+1}).
$$

This is **backward Euler**. The same formula comes from expanding $y(x_n)$ around $x_{n+1}$. Its truncation error is also of size $h^2$ per step. If we use the average of both ends, we get the trapezoidal rule of Week 4. That gives the trapezoidal method in Section 7.

Backward Euler has a catch. The unknown $y_{n+1}$ is on both sides. For $y' = -y^3$, one step is $y_1 = y_0 - h\,y_1^3$, a cubic equation for $y_1$. We must solve an equation at every step.

## 5. Euler's method converges

The big question: does $y_n$ really approach $y(x_n)$?

**Theorem.** Let $f$ be Lipschitz in $y$ with constant $K$, and let the solution $y(x)$ of (1) have a bounded second derivative on $[x_0, b]$. Let $\{y_n\}$ be the Euler solution, and assume $|y(x_0) - y_0| \le c_1 h$ for some $c_1 \ge 0$. Then there is a constant $C \ge 0$, **independent of $h$**, with

$$
\max_{x_0 \le x_n \le b} |y(x_n) - y_n| \le C\,h.
$$

Usually the initial value is exact, so $c_1 = 0$. The theorem says: if we make $h$ smaller, the answer gets better, at a rate proportional to $h$.

**Proof.** Define $\tau_n = \frac{h}{2} y''(\xi_n)$. Then the exact solution and the method satisfy

$$
y(x_{n+1}) = y(x_n) + h\,f(x_n, y(x_n)) + h\tau_n, \qquad y_{n+1} = y_n + h\,f(x_n, y_n).
$$

Subtract, and let $e_n = y(x_n) - y_n$:

$$
e_{n+1} = e_n + h\bigl(f(x_n, y(x_n)) - f(x_n, y_n)\bigr) + h\tau_n. \tag{3}
$$

*Step 1: Lipschitz.* $|f(x_n, y(x_n)) - f(x_n, y_n)| \le K|e_n|$. So

$$
|e_{n+1}| \le (1 + hK)\,|e_n| + h\,|\tau_n|.
$$

*Step 2: repeat.* Use the same bound for $|e_n|$, then for $|e_{n-1}|$, and so on:

$$
|e_{n+1}| \le (1 + hK)^{n+1}|e_0| + h\Bigl(|\tau_n| + (1 + hK)|\tau_{n-1}| + \dots + (1 + hK)^n|\tau_0|\Bigr).
$$

*Step 3: one bound for all $\tau$.* All of $\tau_0, \dots, \tau_n$ appear, not only $\tau_n$. Let $\tilde\tau = \max_{0 \le m \le n} |\tau_m|$. Since $y''$ is bounded,

$$
\tilde\tau \le \frac{h}{2}\,\|y''\|_\infty.
$$

*Step 4: geometric sum.* With $r = hK$,

$$
1 + (1 + r) + \dots + (1 + r)^n = \frac{(1 + r)^{n+1} - 1}{r}.
$$

The factor $h$ in front cancels the $h$ in $r = hK$:

$$
|e_{n+1}| \le (1 + hK)^{n+1}|e_0| + \frac{(1 + hK)^{n+1} - 1}{K}\,\tilde\tau.
$$

*Step 5: remove $h$ from $(1 + hK)^{n+1}$.* This is the key step. The problem is that $n$ also depends on $h$: on a fixed interval, $n = (x_n - x_0)/h$. If $h$ shrinks, $n$ grows. So we need a bound that does not depend on $h$. Use $e^x = 1 + x + \frac{x^2}{2} + \dots \ge 1 + x$ for $x \ge 0$. Then

$$
(1 + hK)^{n+1} \le e^{(n+1)hK} = e^{(x_{n+1} - x_0)K} \le e^{(b - x_0)K}.
$$

The product $nh$ is a length, so $h$ is gone.

*Step 6: put it together.* With $|e_0| \le c_1 h$ and $\tilde\tau \le \frac{h}{2}\|y''\|_\infty$,

$$
|e_{n+1}| \le \left(c_1 e^{(b - x_0)K} + \frac{e^{(b - x_0)K} - 1}{2K}\,\|y''\|_\infty\right) h = C\,h.
$$

The constant $C$ depends on the length of the interval, on $K$, and on $y''$. It does not depend on $h$. $\blacksquare$

**Why must $C$ be independent of $h$?** Suppose $C$ were allowed to be $\frac{1}{h}$ times something. Then $Ch$ is just a constant, and making $h$ small does not help at all. So "$O(h)$" only means something when $C$ has no $h$ inside.

## 6. Euler's method is stable

A numerical method also needs to be stable. Its stability is a separate question from the stability of the problem, and we must check it method by method.

Perturb Euler's method in the same way as (2):

$$
z_{n+1} = z_n + h\bigl(f(x_n, z_n) + \delta(x_n)\bigr), \qquad z_0 = y_0 + \epsilon.
$$

Let $w_n = z_n - y_n$. Then $w_0 = \epsilon$ and

$$
w_{n+1} = w_n + h\bigl(f(x_n, z_n) - f(x_n, y_n)\bigr) + h\,\delta(x_n).
$$

This has exactly the same form as (3). The same steps give

$$
|z_n - y_n| \le e^{(x_n - x_0)K}\,|\epsilon| + \frac{e^{(x_n - x_0)K} - 1}{K}\,\|\delta\|_\infty.
$$

Unlike in the convergence proof, there is no factor $h$ in front. There, $|e_0|$ and $\tilde\tau$ were both of size $h$. Here, $\epsilon$ and $\delta$ do not depend on $h$. (A quick check: with $f \equiv 0$ and $\delta \equiv 0$, we get $z_n - y_n = \epsilon$ for every $n$.)

The constant depends only on the length of the interval and on $K$. It does not blow up as $h \to 0$. So Euler's method is stable.

**Zero-stability and absolute stability.** Both the convergence proof and this stability proof look at the limit $h \to 0$. Stability in this sense is called **zero-stability**. In practice we must choose one fixed $h$, for example $h = 10^{-3}$. Whether the method behaves well for that $h$ is a different question. It is called **absolute stability**, and it comes later in the course.

## 7. Multistep methods

Euler's method uses only $y_n$ to make $y_{n+1}$. But by then we also know $y_{n-1}, y_{n-2}, \dots$ More information might give a better approximation. This is the idea of **multistep methods**:

$$
y_{n+1} = \sum_{j=0}^{p} a_j\,y_{n-j} + h\sum_{j=-1}^{p} b_j\,f(x_{n-j}, y_{n-j}), \qquad n \ge p. \tag{4}
$$

The coefficients $a_j$ and $b_j$ are constants, and $p \ge 0$.

- If $a_p \ne 0$ or $b_p \ne 0$, it is a **$(p+1)$-step method**. It uses $p + 1$ old values.
- Only $y_0$ is given. The starting values $y_1, \dots, y_p$ must come from another method, usually a simple one-step method.
- If $b_{-1} = 0$, the unknown $y_{n+1}$ appears only on the left. The method is **explicit**.
- If $b_{-1} \ne 0$, $y_{n+1}$ is on both sides. The method is **implicit**, and we must solve an equation at each step.

Forward Euler is the case $p = 0$, $a_0 = 1$, $b_{-1} = 0$, $b_0 = 1$. Backward Euler is $p = 0$, $a_0 = 1$, $b_{-1} = 1$, $b_0 = 0$.

Explicit methods look much easier. So why use implicit methods at all? This is a fair question. The answer comes with absolute stability, later in the course.

### Truncation error and consistency

Rearrange (4) and divide by $h$. The left side becomes $\frac{1}{h}\bigl(y_{n+1} - \sum_j a_j y_{n-j}\bigr)$, which plays the role of $y'$. The right side is a linear combination of values $f(x_{n-j}, y_{n-j})$, which are values of $y'$. So we measure how well an exact solution fits the formula. For a differentiable function $y$, define

$$
T_n(y) = y(x_{n+1}) - \left(\sum_{j=0}^{p} a_j\,y(x_{n-j}) + h\sum_{j=-1}^{p} b_j\,y'(x_{n-j})\right), \qquad \tau_n(y) = \frac{T_n(y)}{h}.
$$

For Euler, $T_n = \frac{h^2}{2} y''(\xi_n)$ and $\tau_n = \frac{h}{2}y''(\xi_n)$. This is the same $\tau_n$ as in the proof above.

**Definition.** The method is **consistent** if

$$
\tau(h) = \max_{x_p \le x_n \le b} |\tau_n(y)| \to 0 \qquad \text{as } h \to 0.
$$

Consistency means the part we cut off really vanishes as the grid gets finer.

There are theorems that turn consistency and the order into conditions on the coefficients $a_j$, $b_j$. There is also a theorem that bounds the error of a multistep method by the error of the starting values plus the truncation error. I do not go into their formulas here. The main message is a relation between three words:

- **Consistency:** the truncation error goes to zero.
- **Stability:** small perturbations give small changes.
- **Convergence:** the approximation goes to the exact solution.

**Stable and consistent together give convergent.** If convergence is hard to prove directly, we can prove stability and consistency instead. The Euler proof above is a small example of this: the recursion (3) is the stability part, and $\tau_n \to 0$ is the consistency part.

## 8. The midpoint method

Replace $y'(x_n)$ by a centered difference over two steps. Expand $g$ around $a$:

$$
g(a \pm h) = g(a) \pm h\,g'(a) + \frac{h^2}{2}g''(a) \pm \frac{h^3}{6}g'''(\xi_\pm).
$$

Subtract the two lines. The $g(a)$ and $g''(a)$ terms cancel:

$$
g'(a) = \frac{g(a + h) - g(a - h)}{2h} - \frac{h^2}{12}\bigl(g'''(\xi_+) + g'''(\xi_-)\bigr) = \frac{g(a + h) - g(a - h)}{2h} - \frac{h^2}{6}\,g'''(\xi).
$$

In the last step, the average of two values of a continuous $g'''$ is equal to $g'''(\xi)$ for some $\xi \in [a - h, a + h]$. This is the intermediate value theorem, as in the composite trapezoidal error of Week 4.

Apply this to $y'(x_n) = f(x_n, y(x_n))$:

$$
y(x_{n+1}) = y(x_{n-1}) + 2h\,f(x_n, y(x_n)) + \frac{h^3}{3}\,y'''(\xi_n).
$$

Drop the last term. The result is the **midpoint method**:

$$
y_{n+1} = y_{n-1} + 2h\,f(x_n, y_n), \qquad n \ge 1.
$$

It is explicit, and it is a two-step method. We get $y_1$ from one Euler step and then continue.

**Order.** The term we drop is of size $h^3$. But the order is $O(h^2)$, not $O(h^3)$. The reason: what replaces $y'$ is $\frac{y_{n+1} - y_{n-1}}{2h}$, and its error is $\frac{h^2}{6}y'''$. Dividing by $h$ costs one power. In the language above, $T_n = O(h^3)$ and $\tau_n = T_n/h = O(h^2)$. The order of a method is the order of $\tau$.

Does the Euler start spoil the order? No. One Euler step from the exact $y_0$ has an error of size $h^2$. So the starting error is $O(h^2)$, the same as the truncation error. The experiment below confirms this.

## 9. The trapezoidal method

Integrate $y' = f$ over one step and use the trapezoidal rule with its error from Week 4. The integrand is $f(x, y(x)) = y'(x)$, so its second derivative is $y'''$:

$$
y(x_{n+1}) = y(x_n) + \frac{h}{2}\Bigl(f(x_n, y(x_n)) + f(x_{n+1}, y(x_{n+1}))\Bigr) - \frac{h^3}{12}\,y'''(\xi_n).
$$

Drop the last term. The result is the **trapezoidal method**:

$$
y_{n+1} = y_n + \frac{h}{2}\Bigl(f(x_n, y_n) + f(x_{n+1}, y_{n+1})\Bigr), \qquad n \ge 0. \tag{5}
$$

It is implicit, and it is a one-step method. Again $T_n = O(h^3)$ and the order is $O(h^2)$.

For a PDE $u_t = F(u, x, t)$, the same idea is very famous. It is called the **Crank–Nicolson method**:

$$
\frac{u^{n+1} - u^{n}}{\Delta t} = \frac{1}{2}\Bigl(F(u^n, x, t_n) + F(u^{n+1}, x, t_{n+1})\Bigr).
$$

### Solving the implicit step: Picard iteration

If $f$ is nonlinear in $y$, then (5) is a nonlinear equation for $y_{n+1}$. Any method from Chapter 2 can solve it. The simplest one is a fixed-point iteration, as in [Week 1](/blog/2026-10-03-numerical-analysis-week-1/). Choose a good first guess $y_{n+1}^{(0)}$ (for example $y_n$, or one Euler step), and repeat

$$
y_{n+1}^{(j+1)} = y_n + \frac{h}{2}\Bigl(f(x_n, y_n) + f\bigl(x_{n+1}, y_{n+1}^{(j)}\bigr)\Bigr), \qquad j = 0, 1, \dots
$$

Subtract this from (5). Only the last term differs, and the Lipschitz condition gives

$$
\bigl|y_{n+1} - y_{n+1}^{(j+1)}\bigr| = \frac{h}{2}\Bigl|f(x_{n+1}, y_{n+1}) - f\bigl(x_{n+1}, y_{n+1}^{(j)}\bigr)\Bigr| \le \frac{hK}{2}\,\bigl|y_{n+1} - y_{n+1}^{(j)}\bigr|.
$$

So if $\frac{hK}{2} < 1$, the iterates converge to $y_{n+1}$. For small $h$ this always holds. This is called **Picard iteration**. It is easy to use but slow: the error shrinks by a fixed factor each time (linear convergence). Newton's method would be faster, but it needs $\partial f/\partial y$.

### The error of the trapezoidal method

This proof has the same structure as the Euler proof, with one new term. Subtract (5) from the exact relation. With $e_n = y(x_n) - y_n$ and the Lipschitz condition on both $f$ terms:

$$
|e_{n+1}| \le |e_n| + \frac{h}{2}\bigl(K|e_n| + K|e_{n+1}|\bigr) + \frac{h^3}{12}\|y'''\|_\infty.
$$

The new term is $\frac{hK}{2}|e_{n+1}|$. It comes from the implicit part. Move it to the left:

$$
\left(1 - \frac{hK}{2}\right)|e_{n+1}| \le \left(1 + \frac{hK}{2}\right)|e_n| + \frac{h^3}{12}\|y'''\|_\infty.
$$

For small $h$ we can divide by $1 - \frac{hK}{2}$:

$$
|e_{n+1}| \le \beta\,|e_n| + \tilde\beta\,h^3\,\|y'''\|_\infty, \qquad \beta = \frac{1 + hK/2}{1 - hK/2}, \quad \tilde\beta = \frac{1}{12\,(1 - hK/2)}.
$$

Repeat as before:

$$
|e_{n+1}| \le \beta^{n+1}|e_0| + \bigl(1 + \beta + \dots + \beta^n\bigr)\,\tilde\beta\,h^3\,\|y'''\|_\infty.
$$

Now $\beta - 1 = \frac{hK}{1 - hK/2}$, so the geometric sum is $\frac{(\beta^{n+1} - 1)(1 - hK/2)}{hK}$. Multiply by $\tilde\beta$. The factor $1 - hK/2$ cancels:

$$
|e_{n+1}| \le \beta^{n+1}|e_0| + \frac{\beta^{n+1} - 1}{12K}\,\|y'''\|_\infty\,h^2.
$$

Finally, $\beta = 1 + \frac{hK}{1 - hK/2} \le e^{hK/(1 - hK/2)}$. If $hK \le 1$, this gives $\beta^{n+1} \le e^{2(x_{n+1} - x_0)K}$, a constant without $h$. With $e_0 = 0$, the total error is $O(h^2)$.

So there are two orders here, and it is easy to mix them up. The **local** error, made in one step from the exact value, is $O(h^3)$. The **global** error at a fixed $x$ is $O(h^2)$. About $1/h$ local errors add up, and the sum loses one power of $h$. This is the same story as the composite rules of Week 4.

## 10. The test problem $y' = \lambda y$

At the end of the week we met the standard test problem for every new method:

$$
y' = \lambda y, \qquad y(0) = 1, \qquad y(x) = e^{\lambda x}.
$$

The solution is unique: the only functions whose derivative is a constant multiple of themselves are exponentials, and the initial condition fixes the constant. When someone proposes a new method, they apply it to this problem and report its stability and its order. This makes it easy to compare methods.

The letter $\lambda$ is a hint: it stands for an **eigenvalue**. Recall that for a square $n \times n$ matrix $A$, a number $\lambda$ is an eigenvalue if $A\mathbf{x} = \lambda\mathbf{x}$ for some vector $\mathbf{x} \ne \mathbf{0}$. The eigenvalue can be zero (then $A$ is singular), but the eigenvector cannot be zero. In the direction of $\mathbf{x}$, the matrix acts like multiplication by one number, however large $A$ is. To find eigenvalues we solve $\det(A - \lambda I) = 0$, a polynomial equation of degree $n$. For $n \ge 5$ there is no general formula for the roots, so in practice we compute eigenvalues numerically. The same idea works for operators: $\mathcal{L}u = \lambda u$, for example with $\mathcal{L} = -\nabla^2$.

Why is one scalar equation enough to test a method? For a linear system $\mathbf{y}' = A\mathbf{y}$, we can change to the eigenvectors of $A$ (if $A$ has a full set of them). Then the system splits into one equation $y' = \lambda y$ for each eigenvalue. For a nonlinear equation, a small perturbation $w$ near a solution follows $w' \approx \frac{\partial f}{\partial y}\,w$, which is the same form locally. My ill-conditioned example in Section 2 is the case $\lambda = 50$.

## 11. Experiment: three methods on the logistic equation

I tested the three methods on the logistic equation

$$
y' = y - y^2, \qquad y(0) = 0.1, \qquad x \in [0, 4].
$$

The exact solution is $y(x) = \frac{1}{1 + 9e^{-x}}$. The equation is nonlinear, so the trapezoidal method needs Picard iteration at every step. I started each step from one Euler step and stopped when two iterates agreed to about $10^{-15}$. For each method I measured the maximum error over the grid, with $h = 0.2, 0.1, \dots, 0.2/2^7$.

I also computed the Euler error bound from Section 5. On this solution $y$ stays in $[0.1, 0.86]$, so $K = \max|1 - 2y| = 0.8$. Also $\|y''\|_\infty = 0.0962$, and $c_1 = 0$. This gives $C = \frac{e^{3.2} - 1}{1.6} \cdot 0.0962 = 1.415$.

<figure class="viz">
<div class="viz-legend"><span class="viz-key"><span class="viz-swatch" style="background:#2a78d6"></span>Euler</span><span class="viz-key"><span class="viz-swatch" style="background:#eb6834"></span>midpoint</span><span class="viz-key"><span class="viz-swatch" style="background:#1baf7a"></span>trapezoidal</span><span class="viz-key"><span class="viz-swatch" style="background:#eda100"></span>Euler error bound</span></div>
<div class="viz-frame">
<svg viewBox="0 0 680 400" role="img" tabindex="0" aria-label="Maximum error of three methods for y prime = y minus y squared on [0, 4] against the step size, both on log scales. The Euler error halves each time h halves. The midpoint and trapezoidal errors drop by 4. The dashed Euler error bound stays above the Euler error by a factor of about 11." class="viz-plot"><line x1="64" x2="540" y1="30.0" y2="30.0" class="viz-grid"/><text x="54" y="34.0" text-anchor="end" class="viz-tick">1</text><line x1="64" x2="540" y1="107.5" y2="107.5" class="viz-grid"/><text x="54" y="111.5" text-anchor="end" class="viz-tick">10<tspan dy="-0.5em" font-size="75%">−2</tspan></text><line x1="64" x2="540" y1="185.0" y2="185.0" class="viz-grid"/><text x="54" y="189.0" text-anchor="end" class="viz-tick">10<tspan dy="-0.5em" font-size="75%">−4</tspan></text><line x1="64" x2="540" y1="262.5" y2="262.5" class="viz-grid"/><text x="54" y="266.5" text-anchor="end" class="viz-tick">10<tspan dy="-0.5em" font-size="75%">−6</tspan></text><line x1="64" x2="540" y1="340.0" y2="340.0" class="viz-grid"/><text x="54" y="344.0" text-anchor="end" class="viz-tick">10<tspan dy="-0.5em" font-size="75%">−8</tspan></text><text x="64.0" y="362" text-anchor="middle" class="viz-tick">0.2</text><text x="200.0" y="362" text-anchor="middle" class="viz-tick">0.05</text><text x="336.0" y="362" text-anchor="middle" class="viz-tick">0.0125</text><text x="472.0" y="362" text-anchor="middle" class="viz-tick">0.003125</text><text x="302.0" y="388" text-anchor="middle" class="viz-axis-title">Step size h (halved at each step to the right)</text><text x="64" y="14" class="viz-axis-title">|error|</text><path d="M64.0,91.7 L132.0,103.4 L200.0,115.0 L268.0,126.7 L336.0,138.4 L404.0,150.0 L472.0,161.7 L540.0,173.4" fill="none" stroke="#2a78d6" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><circle cx="540.0" cy="173.4" r="4" fill="#2a78d6" stroke="#ffffff" stroke-width="2"/><path d="M64.0,133.0 L132.0,156.1 L200.0,179.5 L268.0,202.8 L336.0,226.1 L404.0,249.4 L472.0,272.8 L540.0,296.1" fill="none" stroke="#eb6834" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><circle cx="540.0" cy="296.1" r="4" fill="#eb6834" stroke="#ffffff" stroke-width="2"/><path d="M64.0,162.2 L132.0,185.4 L200.0,208.8 L268.0,232.1 L336.0,255.4 L404.0,278.8 L472.0,302.1 L540.0,325.4" fill="none" stroke="#1baf7a" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><circle cx="540.0" cy="325.4" r="4" fill="#1baf7a" stroke="#ffffff" stroke-width="2"/><path d="M64.0,51.2 L132.0,62.9 L200.0,74.6 L268.0,86.2 L336.0,97.9 L404.0,109.6 L472.0,121.2 L540.0,132.9" fill="none" stroke="#eda100" stroke-width="2" stroke-dasharray="6 4" stroke-linejoin="round" stroke-linecap="round"/><circle cx="540.0" cy="132.9" r="4" fill="#eda100" stroke="#ffffff" stroke-width="2"/><text x="550.0" y="177.4" class="viz-label">Euler</text><text x="550.0" y="292.1" class="viz-label">midpoint</text><text x="550.0" y="337.4" class="viz-label">trapezoidal</text><text x="550.0" y="136.9" class="viz-label">bound</text><line x1="0" x2="0" y1="30" y2="340" class="viz-crosshair" visibility="hidden"/></svg>
<div class="viz-tooltip" hidden></div>
</div>
<figcaption>Maximum error over the grid for y′ = y − y², y(0) = 0.1 on [0, 4], against the step size h, both on log scales. Each step to the right halves h. The dashed line is the Euler error bound from the convergence proof. Hover over the chart, or focus it and use the arrow keys, to read the values.</figcaption>
<details class="viz-table"><summary>Show the numbers</summary><table><thead><tr><th>h</th><th>Euler</th><th>midpoint</th><th>trapezoidal</th><th>Euler error bound</th></tr></thead><tbody><tr><td>0.2</td><td>2.6e-2</td><td>2.2e-3</td><td>3.9e-4</td><td>2.8e-1</td></tr><tr><td>0.1</td><td>1.3e-2</td><td>5.6e-4</td><td>9.7e-5</td><td>1.4e-1</td></tr><tr><td>0.05</td><td>6.4e-3</td><td>1.4e-4</td><td>2.4e-5</td><td>7.1e-2</td></tr><tr><td>0.025</td><td>3.2e-3</td><td>3.5e-5</td><td>6.1e-6</td><td>3.5e-2</td></tr><tr><td>0.0125</td><td>1.6e-3</td><td>8.7e-6</td><td>1.5e-6</td><td>1.8e-2</td></tr><tr><td>0.00625</td><td>8.0e-4</td><td>2.2e-6</td><td>3.8e-7</td><td>8.8e-3</td></tr><tr><td>0.003125</td><td>4.0e-4</td><td>5.4e-7</td><td>9.5e-8</td><td>4.4e-3</td></tr><tr><td>0.0015625</td><td>2.0e-4</td><td>1.4e-7</td><td>2.4e-8</td><td>2.2e-3</td></tr></tbody></table></details>
<script>(() => {const fig = document.currentScript.closest('figure');const svg = fig.querySelector('svg');const line = svg.querySelector('.viz-crosshair');const tip = fig.querySelector('.viz-tooltip');const names = ["Euler", "midpoint", "trapezoidal", "Euler error bound"]; const colors = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"];const rows = [["2.6e-2", "2.2e-3", "3.9e-4", "2.8e-1"], ["1.3e-2", "5.6e-4", "9.7e-5", "1.4e-1"], ["6.4e-3", "1.4e-4", "2.4e-5", "7.1e-2"], ["3.2e-3", "3.5e-5", "6.1e-6", "3.5e-2"], ["1.6e-3", "8.7e-6", "1.5e-6", "1.8e-2"], ["8.0e-4", "2.2e-6", "3.8e-7", "8.8e-3"], ["4.0e-4", "5.4e-7", "9.5e-8", "4.4e-3"], ["2.0e-4", "1.4e-7", "2.4e-8", "2.2e-3"]]; const degrees = [20, 40, 80, 160, 320, 640, 1280, 2560];const left = 64, right = 540, width = 680, last = degrees.length - 1;let current = -1;const show = (k) => { k = Math.max(0, Math.min(last, k)); current = k; const x = left + (right - left) * k / last; line.setAttribute('x1', x); line.setAttribute('x2', x); line.setAttribute('visibility', 'visible'); tip.replaceChildren(); const head = document.createElement('div'); head.className = 'viz-tooltip-head'; head.textContent = 'h = ' + (4 / degrees[k]) + ' (' + degrees[k] + ' steps)'; tip.append(head); names.forEach((name, i) => {  const row = document.createElement('div'); row.className = 'viz-tooltip-row';  const key = document.createElement('span'); key.className = 'viz-swatch'; key.style.background = colors[i];  const value = document.createElement('strong'); value.textContent = rows[k][i];  const label = document.createElement('span'); label.textContent = name;  row.append(key, value, label); tip.append(row); }); const box = svg.getBoundingClientRect(), frame = fig.getBoundingClientRect(); const xs = box.left - frame.left + x * box.width / width; tip.hidden = false; const room = frame.width - tip.offsetWidth - 8; tip.style.left = Math.max(0, Math.min(room, xs + 12)) + 'px'; tip.style.top = (box.top - frame.top + 8) + 'px';};const hide = () => { line.setAttribute('visibility', 'hidden'); tip.hidden = true; current = -1; };svg.addEventListener('pointermove', (event) => { const box = svg.getBoundingClientRect(); const x = (event.clientX - box.left) * width / box.width; show(Math.round((x - left) / (right - left) * last));});svg.addEventListener('pointerleave', hide);svg.addEventListener('focus', () => show(current < 0 ? 0 : current));svg.addEventListener('blur', hide);svg.addEventListener('keydown', (event) => { if (event.key === 'ArrowRight') { show(current + 1); event.preventDefault(); } if (event.key === 'ArrowLeft') { show(current - 1); event.preventDefault(); } if (event.key === 'Escape') hide();});})();</script>
</figure>

| $h$ | Euler | ratio | midpoint | ratio | trapezoidal | ratio | Euler bound $Ch$ |
|---|---|---|---|---|---|---|---|
| 0.2 | $2.56 \times 10^{-2}$ | | $2.20 \times 10^{-3}$ | | $3.87 \times 10^{-4}$ | | $2.83 \times 10^{-1}$ |
| 0.1 | $1.28 \times 10^{-2}$ | 2.00 | $5.55 \times 10^{-4}$ | 3.96 | $9.74 \times 10^{-5}$ | 3.98 | $1.42 \times 10^{-1}$ |
| 0.05 | $6.39 \times 10^{-3}$ | 2.00 | $1.39 \times 10^{-4}$ | 4.00 | $2.43 \times 10^{-5}$ | 4.00 | $7.08 \times 10^{-2}$ |
| 0.025 | $3.20 \times 10^{-3}$ | 2.00 | $3.48 \times 10^{-5}$ | 4.00 | $6.08 \times 10^{-6}$ | 4.00 | $3.54 \times 10^{-2}$ |
| 0.0125 | $1.60 \times 10^{-3}$ | 2.00 | $8.69 \times 10^{-6}$ | 4.00 | $1.52 \times 10^{-6}$ | 4.00 | $1.77 \times 10^{-2}$ |
| 0.00625 | $7.98 \times 10^{-4}$ | 2.00 | $2.17 \times 10^{-6}$ | 4.00 | $3.80 \times 10^{-7}$ | 4.00 | $8.85 \times 10^{-3}$ |
| 0.003125 | $3.99 \times 10^{-4}$ | 2.00 | $5.43 \times 10^{-7}$ | 4.00 | $9.51 \times 10^{-8}$ | 4.00 | $4.42 \times 10^{-3}$ |
| 0.0015625 | $2.00 \times 10^{-4}$ | 2.00 | $1.36 \times 10^{-7}$ | 4.00 | $2.38 \times 10^{-8}$ | 4.00 | $2.21 \times 10^{-3}$ |

The ratio is "error of the row above / error of this row". When $h$ halves, the Euler error halves (ratio $2$, first order). The midpoint and trapezoidal errors drop by $4$ (second order). The midpoint method starts with one Euler step, and this does not spoil its order, as Section 8 predicted. At the same $h$ the trapezoidal error is about $6$ times smaller than the midpoint error, but each trapezoidal step needs Picard iterations: on average $10$ per step at $h = 0.2$, and $3.8$ at the smallest $h$, because a smaller $h$ gives a better first guess and a smaller contraction factor.

The Euler bound is always above the true error, as it must be. But it is about $11$ times too large. The proof uses the worst case at every step: it assumes the error grows by $1 + hK$ each time. In this problem $\partial f/\partial y = 1 - 2y$ becomes negative once $y > 0.5$, so the true error is damped in the second half of the interval. The bound cannot see this.

**Part B: local and global error of the trapezoidal method.** The one-step error starts from the exact $y(0)$ and takes one step. The global error is the maximum over $[0, 4]$.

| $h$ | one-step error | ratio | global error | ratio |
|---|---|---|---|---|
| 0.4 | $-2.39 \times 10^{-4}$ | | $1.55 \times 10^{-3}$ | |
| 0.2 | $-2.92 \times 10^{-5}$ | 8.21 | $3.87 \times 10^{-4}$ | 4.01 |
| 0.1 | $-3.56 \times 10^{-6}$ | 8.20 | $9.74 \times 10^{-5}$ | 3.98 |
| 0.05 | $-4.38 \times 10^{-7}$ | 8.12 | $2.43 \times 10^{-5}$ | 4.00 |
| 0.025 | $-5.43 \times 10^{-8}$ | 8.06 | $6.08 \times 10^{-6}$ | 4.00 |

The one-step error drops by $8$ ($O(h^3)$) and the global error by $4$ ($O(h^2)$). This is the "two orders" point of Section 9 in numbers.

**Part C: Picard iteration.** For the first step from $y_0 = 0.1$, starting from the guess $y_1^{(0)} = y_0$, the error shrinks by a constant factor each time:

| $h$ | errors of iterates $0, 1, \dots, 5$ | factor per iteration | $\frac{h}{2}\lvert\partial f/\partial y\rvert$ at the root |
|---|---|---|---|
| 0.4 | $4.2 \times 10^{-2}$, $6.4 \times 10^{-3}$, $9.3 \times 10^{-4}$, $1.3 \times 10^{-4}$, $1.9 \times 10^{-5}$, $2.7 \times 10^{-6}$ | 0.143 | 0.143 |
| 0.1 | $9.4 \times 10^{-3}$, $3.7 \times 10^{-4}$, $1.5 \times 10^{-5}$, $5.7 \times 10^{-7}$, $2.2 \times 10^{-8}$, $8.6 \times 10^{-10}$ | 0.0391 | 0.0391 |

The measured factor matches $\frac{h}{2}|\partial f/\partial y|$ at the root. It is below the worst case $\frac{hK}{2}$ ($0.16$ and $0.04$). This is linear convergence, as predicted.

When $\frac{hK}{2} > 1$, Picard iteration fails. I tried it on the ill-conditioned example of Section 2, where $K = 50$, with $h = 0.1$, so $\frac{hK}{2} = 2.5$. That equation is linear, so I could also solve (5) directly: $y_1 = 0.99500$. The Picard iterates were $1.0075$, $1.0262$, $1.0731$, $1.1903$, $1.4831$, $2.2154$. The distance to $y_1$ grows by $2.5$ every time. The trapezoidal step itself is fine. Only the simple way of solving it breaks.

**Part D: Euler on the ill-conditioned problem.** Back to $y' = 50(y - \cos x) - \sin x$, $y(0) = 1$, exact $y = \cos x$. Here $K = 50$ and $\|y''\|_\infty = 1$, so the Euler bound gives $C = \frac{e^{50} - 1}{100} = 5.2 \times 10^{19}$.

| $h$ | error at $x = 0.2$ | error at $x = 0.5$ | error at $x = 1$ | error at $x = 1$, divided by $h$ |
|---|---|---|---|---|
| $10^{-3}$ | $0.17$ | $3.9 \times 10^{5}$ | $1.5 \times 10^{16}$ | $1.5 \times 10^{19}$ |
| $10^{-4}$ | $0.021$ | $6.8 \times 10^{4}$ | $4.6 \times 10^{15}$ | $4.6 \times 10^{19}$ |
| $10^{-5}$ | $0.0022$ | $7.2 \times 10^{3}$ | $5.1 \times 10^{14}$ | $5.1 \times 10^{19}$ |

Euler still converges: each column shrinks by about $10$ when $h$ shrinks by $10$. The theorem is right. But the constant is about $5 \times 10^{19}$, and here the bound is almost sharp, because the error really does grow like $e^{50x}$. Even $100{,}000$ steps give an error of $5 \times 10^{14}$ at $x = 1$. This is what "stable but ill-conditioned" means in practice.

Here is the core of the code, in Python:

```python
def euler(n):
    h = (B - A) / n
    ys = [Y0]
    for k in range(n):
        ys.append(ys[-1] + h * f(A + k * h, ys[-1]))
    return ys

def midpoint(n):
    """Two-step method; y_1 comes from one Euler step."""
    h = (B - A) / n
    ys = [Y0, Y0 + h * f(A, Y0)]
    for k in range(1, n):
        ys.append(ys[k - 1] + 2 * h * f(A + k * h, ys[k]))
    return ys

def picard_step(x, y, h, guess, tol=1e-15, max_iter=200):
    """Solve z = y + h/2 (f(x, y) + f(x + h, z)) by Picard iteration."""
    known = y + 0.5 * h * f(x, y)
    z = guess
    for it in range(1, max_iter + 1):
        z_new = known + 0.5 * h * f(x + h, z)
        if abs(z_new - z) <= tol * max(1.0, abs(z_new)):
            return z_new, it
        z = z_new
    return z, max_iter
```

The full script, with all the numbers and the chart, is in [this site's repository](https://github.com/doubleyands/doubleyands.github.io/tree/main/scripts/posts/numerical-analysis/week-5).

## Summary

| topic | main idea |
|---|---|
| Existence and uniqueness | $f$ continuous and Lipschitz in $y$ gives a unique solution on $[x_0 - \gamma, x_0 + \gamma]$, $\gamma = \min\{\alpha, \beta/M\}$. Local only. Sufficient, not necessary. |
| Stability of the problem | Small changes in $f$ or $y_0$ give small changes in $y$. Different from the stability of a method, and different from ill-conditioning. |
| Ill-conditioned example | $y' = 50(y - \cos x) - \sin x$: a change $\epsilon$ in $y(0)$ grows like $\epsilon e^{50x}$. |
| Higher order equations | Rewrite as a first order system. Order $m$ gives $m$ unknowns and $m$ initial conditions. |
| Euler's method | $y_{n+1} = y_n + h f(x_n, y_n)$. From Taylor, a difference quotient, or a rectangle. Converges with $O(h)$; the constant has no $h$. Zero-stable. |
| Multistep methods | Use several old values. Explicit if $b_{-1} = 0$, implicit otherwise. Consistent if $\tau(h) \to 0$. Stable and consistent give convergent. |
| Midpoint method | $y_{n+1} = y_{n-1} + 2h f(x_n, y_n)$. Explicit, two-step, $O(h^2)$. |
| Trapezoidal method | $y_{n+1} = y_n + \frac{h}{2}(f_n + f_{n+1})$. Implicit, one-step, $O(h^2)$. Crank–Nicolson for PDEs. Picard iteration converges if $hK/2 < 1$. |
| Local vs. global | Local error $O(h^3)$, global error $O(h^2)$. |
| Test problem | $y' = \lambda y$. Every new method is tested on it. |

The main thing I take from this week: convergence is not one idea but two ideas put together. Consistency says that each step makes only a small error. Stability says that these errors do not grow out of control. The Euler proof shows both in a few lines, and the trapezoidal proof has the same shape. The experiment also shows the limits of a proof. A bound with $e^{(b - x_0)K}$ can be loose by a factor of $11$, or almost sharp, as in the ill-conditioned problem. In both cases it is true, but it does not tell us how big the error really is.
