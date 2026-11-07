---
title: "Numerical Analysis, Week 6: Stability Regions, Runge–Kutta, and the Root Condition"
description: "Stability of ODE methods on the test problem y' = λy, A-stability of the trapezoidal method, Taylor and Runge–Kutta methods, Milne's method, one-step vs. multistep methods, the root condition, stability regions and the boundary locus method. Tested on a stiff problem with λ = −1000."
pubDate: 2026-11-07
project: numerical-analysis
tags: [numerical analysis, ODE, stability, runge-kutta, multistep methods, stiff problems, notes]
---

This is my summary of the sixth week of a numerical analysis course. As before, I write the ideas in my own words and check some of them with a small experiment. The course follows K. E. Atkinson, *An Introduction to Numerical Analysis*.

Last week we built Euler's method, the midpoint method and the trapezoidal method, and we met the test problem $y' = \lambda y$ (see [Week 5](/blog/2026-10-31-numerical-analysis-week-5/)). Last week also left two questions open. Why is Euler's method often slow in practice? And why would anyone use an implicit method, when every step needs an equation solve? This week answers both. The answer is **stability for a fixed step size $h$**.

The plan:

1. Put a method on $y' = \lambda y$ and watch how a small perturbation grows.
2. More one-step methods: Taylor series methods and Runge–Kutta methods.
3. More multistep methods, and a comparison of one-step and multistep methods.
4. The root condition, which decides if a multistep method converges.
5. Stability regions in the complex plane, and how to draw them.

## 1. Stability on the test problem

To study the stability of a method, we always apply it to the test problem

$$
y' = \lambda y, \qquad y(0) = 1, \qquad y(x) = e^{\lambda x}.
$$

Accuracy is easier: we can study it with Taylor expansions or with integration rules. Stability needs this test problem.

The solution is unique. Take two solutions $y_1$ and $y_2$ and let $d = y_1 - y_2$. Then $d' = \lambda d$ and $d(0) = 0$. A function whose derivative is a constant times itself must be $Ce^{\lambda x}$, and $d(0) = 0$ forces $C = 0$. So $y_1 = y_2$. This works because the problem is linear. (It also follows from last week's theorem, since $\lambda y$ is Lipschitz in $y$ with $K = |\lambda|$.)

Be careful with the word "eigenvalue" here. For a matrix, one eigenvalue can have many independent eigenvectors. Here the solution is unique only because it must be an exponential.

### The trapezoidal method with a perturbed start

Take a slightly more general linear problem, $y' = \lambda y + g(x)$ with $y(0) = y_0$. The trapezoidal method gives

$$
\Bigl(1 - \frac{h\lambda}{2}\Bigr) y_{n+1} = \Bigl(1 + \frac{h\lambda}{2}\Bigr) y_n + \frac{h}{2}\bigl(g(x_n) + g(x_{n+1})\bigr).
$$

Now start from $z_0 = y_0 + \epsilon$ instead and call the result $z_n$. Let $w_n = z_n - y_n$. The $g$ terms are the same in both, so they cancel:

$$
w_{n+1} = r\, w_n, \qquad r = \frac{1 + h\lambda/2}{1 - h\lambda/2}, \qquad w_n = r^n \epsilon.
$$

We need $h\lambda \ne 2$, or the denominator is zero. The number $r$ is the **amplification factor**: each step multiplies the perturbation by $r$. A method is stable on this problem if a small change at the start stays small. So we want $|r| \le 1$, and $|r| < 1$ is even better, because then $w_n \to 0$.

### When is $|r| < 1$?

The solution of $y' = \lambda y$ may be complex. For example, $y' = iy$ has the solution $e^{ix}$. So $\lambda$ lives in the complex plane. Write $w = h\lambda/2$. For any complex $w$,

$$
|1 + w|^2 - |1 - w|^2 = 4\,\mathrm{Re}(w).
$$

So $|r| = \frac{|1 + w|}{|1 - w|} < 1$ exactly when $\mathrm{Re}(w) < 0$, that is, when $\mathrm{Re}(\lambda) < 0$. The size of $h$ does not matter. The imaginary part of $\lambda$ does not matter either. If $\lambda$ is purely imaginary, $|r| = 1$: the perturbation keeps its size. It does not grow, but it does not die out.

**A-stability.** A method is called **A-stable** if, for every $h > 0$ and every $\lambda$ with $\mathrm{Re}(\lambda) < 0$, the numerical solution of $y' = \lambda y$ goes to $0$. The trapezoidal method is A-stable.

**Stable is not the same as accurate.** In theory, any $h$ is fine. In practice, $h$ must still be small. Here are the one-step factors for $\lambda = -30$. The exact factor is $e^{h\lambda}$.

| $h$ | $z = h\lambda$ | exact $e^{z}$ | forward Euler $1 + z$ | backward Euler $\frac{1}{1 - z}$ | trapezoidal $r$ |
|---|---|---|---|---|---|
| 0.01 | $-0.3$ | $0.7408$ | $0.7$ | $0.7692$ | $0.7391$ |
| 0.05 | $-1.5$ | $0.2231$ | $-0.5$ | $0.4$ | $0.1429$ |
| 0.25 | $-7.5$ | $0.00055$ | $-6.5$ | $0.1176$ | $-0.5789$ |

At $h = 0.25$ the true solution almost vanishes in one step. The trapezoidal factor is $-0.58$: the numerical solution does not blow up, but it changes sign every step and shrinks slowly. Forward Euler, with factor $-6.5$, explodes.

## 2. One-step methods

The three methods we know are:

| method | formula | type | order |
|---|---|---|---|
| forward Euler | $y_{n+1} = y_n + h f(x_n, y_n)$ | explicit | 1 |
| backward Euler | $y_{n+1} = y_n + h f(x_{n+1}, y_{n+1})$ | implicit | 1 |
| trapezoidal | $y_{n+1} = y_n + \frac{h}{2}\bigl(f(x_n, y_n) + f(x_{n+1}, y_{n+1})\bigr)$ | implicit | 2 |

Each one uses only $y_n$ to build $y_{n+1}$, so each is a **one-step method**.

**Forward and backward Euler from Taylor.** Expand around $x_n$:

$$
y(x_{n+1}) = y(x_n) + h y'(x_n) + \frac{h^2}{2} y''(\xi_n)
\quad\Longrightarrow\quad
\frac{y(x_{n+1}) - y(x_n)}{h} = y'(x_n) + \frac{h}{2} y''(\xi_n).
$$

The difference quotient equals $y'(x_n)$ up to a term $\frac{h}{2}y''$. That term is a constant times $h$, and the constant does not depend on $h$. So forward Euler is first order. For backward Euler, expand around $x_{n+1}$ instead:

$$
y(x_n) = y(x_{n+1}) - h y'(x_{n+1}) + \frac{h^2}{2} y''(\eta_n).
$$

The sign of the $h$ term is negative, so we divide by $-h$. The result is again a first order formula, now for $y'(x_{n+1})$. The name says where we stand: forward Euler steps forward from the known point $x_n$. Backward Euler stands at the unknown point $x_{n+1}$ and looks back.

**Why explicit is easy.** Take $f(x, y) = e^{-xy}$. Forward Euler gives $y_1 = y_0 + h e^{-x_0 y_0}$. Every value on the right is known, and $y_2$ follows from $y_1$ in the same way. Backward Euler gives $y_1 = y_0 + h e^{-x_1 y_1}$. Now $y_1$ sits inside the exponential, so we must solve a nonlinear equation at every step. The trapezoidal method has the same problem.

**The trapezoidal method is the average of the two Euler methods.** Its right-hand side is the mean of the two Euler right-hand sides. It is second order. We can see this from the trapezoidal rule error (Week 5), or from the two Taylor expansions above with one more term each:

$$
\frac{y(x_{n+1}) - y(x_n)}{h} = y'(x_n) + \frac{h}{2} y''(x_n) + \frac{h^2}{6} y'''(x_n) + \cdots
$$

$$
\frac{y(x_{n+1}) - y(x_n)}{h} = y'(x_{n+1}) - \frac{h}{2} y''(x_{n+1}) + \frac{h^2}{6} y'''(x_{n+1}) - \cdots
$$

Add them and divide by $2$. The $y''$ terms appear as $\frac{h}{4}\bigl(y''(x_n) - y''(x_{n+1})\bigr)$. The signs are opposite, and this matters. Since $y''(x_{n+1}) = y''(x_n) + h y'''(\cdot)$, the difference is only $O(h)$, so the whole term is $O(h^2)$. Collecting the $h^2$ terms gives $\bigl(\frac{1}{6} - \frac{1}{4}\bigr)h^2 y''' = -\frac{h^2}{12} y'''$. With a plus sign instead, a term $\frac{h}{2} y''$ would remain and the method would be only first order.

### Taylor series methods

Can a one-step method have a higher order? The most direct idea is to keep more Taylor terms:

$$
y(x_{n+1}) = y(x_n) + h y'(x_n) + \frac{h^2}{2} y''(x_n) + \cdots + \frac{h^p}{p!} y^{(p)}(x_n).
$$

We know $y' = f$. We can get $y''$ by differentiating the ODE. For an autonomous problem $y' = f(y)$, the chain rule gives $y'' = f'(y)\,y' = f f'$. Then the error term starts at $h^2$, and so on for higher orders.

This is a **Taylor series method** of order $p$. It is explicit and it can have any order. But it is **not practical**. When $f$ depends on $x$ too, the derivatives grow fast:

$$
y'' = f_x + f_y f, \qquad y''' = f_{xx} + 2 f_{xy} f + f_{yy} f^2 + f_y (f_x + f_y f).
$$

For a system, these become matrices and tensors. And if $f$ is only given as a computer program, we may not be able to differentiate it at all. Almost nobody uses this method.

## 3. Runge–Kutta methods

The popular high order one-step methods are the **Runge–Kutta methods**. They do not differentiate $f$. Instead they evaluate $f$ at several points inside one step. One step has several **stages**, so these are *multistage* one-step methods.

### RK2

The idea: take half a step first, and build a value at the midpoint $x_{n+1/2} = x_n + h/2$. The midpoint is closer to $x_n$ than $x_{n+1}$ is, so this value is more accurate. Then use the slope at the midpoint for the full step.

- **Stage 1.** One Euler step of size $h/2$: $\;y^* = y_n + \frac{h}{2} f(x_n, y_n)$.
- **Stage 2.** The full step with the midpoint slope: $\;y_{n+1} = y_n + h f\bigl(x_{n+1/2},\, y^*\bigr)$.

Why the midpoint slope? Integrate $y' = f$ over one step: $y(x_{n+1}) - y(x_n) = \int_{x_n}^{x_{n+1}} f\,dx$. If we must approximate this integral with one value of $f$, we can use the left end, the right end, or the midpoint. The midpoint rule is the most accurate of the three: its error is $O(h^3)$ on one step, and the two end-point rules have $O(h^2)$. We do not know $y$ at the midpoint, so we use $y^*$ in its place. This is the **second order Runge–Kutta method**:

$$
y_{n+1} = y_n + h f\Bigl(x_{n+1/2},\ y_n + \frac{h}{2} f(x_n, y_n)\Bigr).
$$

### RK4

The most widely used member is the **fourth order, four-stage Runge–Kutta method**:

$$
\begin{aligned}
Y_1 &= y_n, \\
Y_2 &= y_n + \tfrac{h}{2} f(x_n, Y_1), \\
Y_3 &= y_n + \tfrac{h}{2} f(x_{n+1/2}, Y_2), \\
Y_4 &= y_n + h f(x_{n+1/2}, Y_3), \\
y_{n+1} &= y_n + \tfrac{h}{6}\Bigl[f(x_n, Y_1) + 2 f(x_{n+1/2}, Y_2) + 2 f(x_{n+1/2}, Y_3) + f(x_{n+1}, Y_4)\Bigr].
\end{aligned}
$$

There is a simple way to read it:

- **The last line is Simpson's rule.** Simpson's rule on one step is $\frac{h}{6}\bigl(f_n + 4 f_{n+1/2} + f_{n+1}\bigr)$. RK4 splits the middle weight $4$ into $2 + 2$ and uses two different guesses for the midpoint value.
- **$Y_2$ is a forward Euler guess for the midpoint.** It takes the slope at $x_n$: $y_{n+1/2} \approx y_n + \frac{h}{2} f(x_n, y_n)$.
- **$Y_3$ is like a backward Euler guess for the midpoint.** Backward Euler would be $y_{n+1/2} = y_n + \frac{h}{2} f(x_{n+1/2}, y_{n+1/2})$. We already have a midpoint value $Y_2$, so we put it on the right-hand side. No equation solve is needed.
- **$Y_4$ is a full step** from $y_n$, with the slope from the midpoint guess $Y_3$.

Every stage starts from $y_n$. Only the slope changes.

Is this really fourth order? The check is not hard, but it is long. Assume for simplicity that $f$ depends only on $y$. The midpoint value is $f\bigl(y_n + \frac{h}{2} f(y_n)\bigr)$, an $f$ inside an $f$. Expand it as $f(y_n + \delta) = f + \delta f' + \frac{\delta^2}{2} f'' + \cdots$ with $\delta = \frac{h}{2} f$. Do the same for the other stages, put everything together, and compare with the Taylor series of the exact solution. All terms agree up to $h^4$. The experiment in Section 10 checks the order with numbers.

## 4. More multistep methods

Recall the general multistep method from Week 5:

$$
y_{n+1} = \sum_{j=0}^{p} a_j\, y_{n-j} + h \sum_{j=-1}^{p} b_j\, f(x_{n-j}, y_{n-j}), \qquad n \ge p. \tag{1}
$$

It is explicit if $b_{-1} = 0$, and implicit otherwise. There are two main ways to build methods of high order:

1. **Undetermined coefficients.** Choose the $a_j$ and $b_j$ so that the Taylor expansion matches to the order we want.
2. **Numerical integration.** Integrate the ODE over some interval $[x_{n-r}, x_{n+1}]$:
   $$
   y(x_{n+1}) = y(x_{n-r}) + \int_{x_{n-r}}^{x_{n+1}} f(t, y(t))\,dt,
   $$
   and replace the integrand by a polynomial that interpolates $f$ at the grid points. Each integration rule gives a method.

**Milne's method.** Take $r = 1$ and Simpson's rule on $[x_{n-1}, x_{n+1}]$. The interval has length $2h$, so the weight $\frac{2h}{6}$ becomes $\frac{h}{3}$:

$$
y_{n+1} = y_{n-1} + \frac{h}{3}\Bigl[f(x_{n-1}, y_{n-1}) + 4 f(x_n, y_n) + f(x_{n+1}, y_{n+1})\Bigr].
$$

The Simpson error is $-\frac{h^5}{90} y^{(5)}$, so after dividing by $h$ the method is fourth order. It is implicit.

**Explicit and implicit families.** If we integrate only over the last step, $[x_n, x_{n+1}]$, and interpolate at older points, we get explicit methods. If we also use the new point $x_{n+1}$, we get implicit methods. There are standard families of both kinds for every order, and you can look up their coefficients when you need them. Two members are old friends: the implicit family starts with backward Euler (order 1) and the trapezoidal method (order 2).

One pattern is worth knowing. At the same order, the truncation error constant of the implicit member is usually several times smaller than the constant of the explicit one. This is one reason to use implicit methods. It is not the main reason. Milne's method has an even smaller error constant than the implicit method of the same order, but people still prefer the implicit family. The reason is stability, and Section 8 shows what goes wrong with Milne.

### One-step vs. multistep

| | one-step methods | multistep methods |
|---|---|---|
| starting | Self-starting: $y_0$ is enough. | Need $y_0, \dots, y_p$. The first values must come from another method, usually a one-step method. |
| changing $h$ | Easy at any step, because only $y_n$ is used. Good for adaptive step sizes. | Hard. The old values sit on the old grid, so the step size must stay the same for a while. |
| a point where $y$ is not smooth | Put the point on the grid, and we usually keep the full accuracy. | The formula uses values on both sides of the point, so accuracy can drop. |
| cost per step | Several evaluations of $f$ (RK4 needs $4$), or derivatives of $f$ (Taylor). | Old values of $f$ are stored and reused. Usually only one new $f$ per step. |

The last row matters when $f$ is expensive. A fourth order explicit multistep method needs one new $f$ per step, and RK4 needs four. We choose by the needs of the problem.

## 5. Convergence and the root condition

For a multistep method, the old values are tangled together, so the direct stability analysis of Week 5 becomes hard. Fortunately there is a simple test. First, the three ideas, in words:

- **Stable**: if the starting values $y_0, \dots, y_p$ change by at most $\epsilon$, then all later values change by at most $c\,\epsilon$, with $c$ independent of $h$ (for $h$ small enough).
- **Convergent**: if the starting values tend to the exact values as $h \to 0$, then the maximum error over the whole interval tends to $0$.
- **Consistent**: the truncation error goes to $0$ as $h \to 0$ (Week 5).

### The root condition

Look at the $a_j$ in (1) only, and define

$$
\rho(r) = r^{p+1} - \sum_{j=0}^{p} a_j\, r^{p-j}.
$$

Where does it come from? Let $h \to 0$ in (1). The $f$ terms vanish and we are left with $y_{n+1} = \sum a_j y_{n-j}$. Try $y_n = r^n$, divide by $r^{n-p}$, and we get $\rho(r) = 0$. Consistency includes $\sum a_j = 1$, so $\rho(1) = 0$: one root is always $r_0 = 1$.

$\rho$ has degree $p + 1$, so it has $p + 1$ complex roots $r_0 = 1, r_1, \dots, r_p$, counted with multiplicity. The method satisfies the **root condition** if

1. $|r_j| \le 1$ for all $j$: every root lies in the closed unit disk, and
2. $|r_j| = 1 \Rightarrow \rho'(r_j) \ne 0$: every root on the unit circle is a simple root.

Why these two rules? A small error in the starting values spreads as a combination of the solutions $r_j^n$. To cover the interval we need $N = (b - x_0)/h$ steps.

- If $|r_j| > 1$, then $r_j^n$ grows exponentially. A smaller $h$ means more steps, so the error gets *worse* as $h \to 0$.
- If a root on the circle is a double root, there is also a solution $n\,r_j^n$. It grows like $n \sim 1/h$, which also blows up as $h \to 0$.

**Main theorem.** Assume the method is consistent. Then

$$
\text{convergent} \iff \text{root condition} \iff \text{stable}.
$$

The assumption matters. The root condition only looks at the $a_j$, so it cannot guarantee consistency by itself. The usual short form, "stable + consistent $\iff$ convergent", says the same thing.

### An unstable method

Here is my own example. The explicit two-step method

$$
y_{n+1} = -4 y_n + 5 y_{n-1} + h\bigl(4 f(x_n, y_n) + 2 f(x_{n-1}, y_{n-1})\bigr)
$$

is consistent. In fact it is third order: its one-step error drops by about $16$ when $h$ halves (see the table below). But

$$
\rho(r) = r^2 + 4r - 5 = (r - 1)(r + 5),
$$

and the root $-5$ is far outside the unit circle. So the root condition fails, and the method cannot converge.

Let us see it directly. Solve $y' = 0$, $y(0) = 0$, so the exact solution is $0$. Start with $z_0 = 0$ and $z_1 = \epsilon$. Since $f = 0$, we get $z_{n+1} = -4 z_n + 5 z_{n-1}$. The general solution is $A \cdot 1^n + B(-5)^n$, and the starting values give $B = -\epsilon/6$. So $|z_n|$ grows like $\frac{\epsilon}{6}\,5^n$. With $\epsilon = 10^{-16}$, a single rounding error, $|z_{40}|$ is already about $1.5 \times 10^{11}$.

I also ran it on $y' = -y$, $y(0) = 1$, on $[0, 1]$, with the exact values for $y_0$ and $y_1$:

| $h$ | one-step error | ratio | error at $x = 1$ |
|---|---|---|---|
| 0.2 | | | $3.1 \times 10^{-2}$ |
| 0.1 | $-1.54 \times 10^{-5}$ | | $7.0$ |
| 0.05 | $-1.00 \times 10^{-6}$ | 15.4 | $4.7 \times 10^{6}$ |
| 0.025 | $-6.38 \times 10^{-8}$ | 15.7 | $2.9 \times 10^{19}$ |
| 0.0125 | $-4.03 \times 10^{-9}$ | 15.8 | $1.7 \times 10^{46}$ |

The one-step error is excellent and shrinks like $h^4$ (third order). But the global error explodes, and it explodes faster when $h$ is smaller. Here even the tiny truncation errors feed the root $-5$. A good order is useless without stability.

## 6. Fixed $h$: the characteristic polynomial

The main theorem is about $h \to 0$. In practice we pick one $h$. How small must it be? For this we go back to the test problem. Apply (1) to $y' = \lambda y$, so $f(x_{n-j}, y_{n-j}) = \lambda y_{n-j}$:

$$
(1 - h\lambda\, b_{-1})\, y_{n+1} - \sum_{j=0}^{p} (a_j + h\lambda\, b_j)\, y_{n-j} = 0.
$$

This is a linear difference equation of order $p + 1$. Try $y_n = r^n$ and divide by $r^{n-p}$. With $z = h\lambda$ we get the **characteristic polynomial**

$$
\pi(r) = (1 - z\,b_{-1})\, r^{p+1} - \sum_{j=0}^{p} (a_j + z\, b_j)\, r^{p-j} = \rho(r) - z\,\sigma(r),
$$

where $\sigma(r) = b_{-1} r^{p+1} + \sum_{j=0}^{p} b_j r^{p-j}$. Its roots are the **characteristic roots**. Only the product $z = h\lambda$ matters, not $h$ and $\lambda$ separately.

Compare the two polynomials:

- $\rho(r)$ is $\pi(r)$ at $z = 0$. Its root condition tells us about the limit $h \to 0$. This is **zero-stability**. It is a statement about theory: it says nothing about $\lambda$.
- $\pi(r)$ at a fixed $z$ tells us whether the method behaves well for this $h$ and this $\lambda$. This is **absolute stability**. For a real problem, this is the polynomial to watch.

### An example: Euler on a stiff problem

Here is my own example:

$$
y' = -1000\,\bigl(y - \sin x\bigr) + \cos x, \qquad y(0) = 0.
$$

The exact solution is $y = \sin x$: then the bracket is zero and $y' = \cos x$. The solution is smooth and slow. But near it, a perturbation $w$ follows $w' = -1000\,w$. So $\lambda = -1000$, and every error is multiplied by about $1 + h\lambda = 1 - 1000h$ per forward Euler step.

If $h > 0.002$, then $|1 - 1000h| > 1$ and every small error grows. At $h = 1/450$, forward Euler ends with an error of $5 \times 10^{29}$. At $h = 1/510$ the error is $8 \times 10^{-7}$. The method does converge as $h \to 0$, and the theory is right. But until $h$ is below $2/1000$, the answer is garbage. Section 10 has the full numbers.

This is the answer to last week's question. The exact solution $\sin x$ is gentle, and accuracy alone would allow a large step. Yet stability forces $h < 2/|\lambda|$. A problem like this, where a strongly damped component forces a tiny $h$ on an explicit method, is called **stiff**.

## 7. Regions of absolute stability

For the test problem, perturb the start by $\epsilon$ as in Section 1.

**Forward Euler.** $w_{n+1} = (1 + h\lambda)\,w_n$, so $w_n = (1 + z)^n \epsilon$. The method is absolutely stable if $|1 + z| \le 1$. For real $\lambda < 0$ this means $-2 \le z \le 0$: the **interval of absolute stability** is $[-2, 0]$. In the complex plane, $|1 + z| = 1$ is the circle with center $-1$ and radius $1$. To see which side is stable, test one point: $z = -1$ gives $0 \le 1$. So the **region of absolute stability** is the inside of the disk.

**Backward Euler.** $(1 - h\lambda)\,w_{n+1} = w_n$, so $w_n = (1 - z)^{-n}\epsilon$. We need $|1 - z| \ge 1$. The boundary is the circle with center $1$ and radius $1$. The test point $z = 1$ fails, so the stable region is everything *outside* this disk. It contains the whole left half plane.

**Trapezoidal.** From Section 1, the region is $\mathrm{Re}(z) \le 0$, the closed left half plane. In a way it sits between the two Euler regions, just as the method is their average.

The difference between explicit and implicit is huge. With forward Euler, $h$ must be small enough to put $h\lambda$ in a small disk. With backward Euler or the trapezoidal method, any $h$ is stable when $\mathrm{Re}(\lambda) < 0$. (In practice, $h$ must still be small for accuracy, as the table in Section 1 shows.)

**Definition.** The region of absolute stability of a multistep method is the set of $z \in \mathbb{C}$ for which $\pi(r)$ satisfies the root condition: $|r_j| \le 1$ for all roots, and $|r_j| < 1$ for repeated roots.

For one-step methods we can use the direct perturbation argument, as above. For multistep methods that argument is very hard, so we use $\pi(r)$. For example, forward Euler has $\pi(r) = r - (1 + z)$, with the single root $r = 1 + z$. Backward Euler has $\pi(r) = (1 - z)\,r - 1$, with the root $r = (1 - z)^{-1}$. Both give the regions above.

## 8. Drawing the region: the boundary locus method

For a method with many steps, $\pi(r)$ has high degree. There is no formula for the roots of a general polynomial of degree five or more, so we cannot write the region down in closed form. We need a numerical way to draw it.

**Idea.** If $z$ is on the boundary of the region, then $\pi(r)$ has a root with $|r| = 1$. Such a root can be written as $r = e^{i\theta}$. Then

$$
\pi(e^{i\theta}) = 0 \;\Longrightarrow\; \rho(e^{i\theta}) - z\,\sigma(e^{i\theta}) = 0 \;\Longrightarrow\; z(\theta) = \frac{\rho(e^{i\theta})}{\sigma(e^{i\theta})}.
$$

So every boundary point lies on the curve $z(\theta)$, $0 \le \theta \le 2\pi$. The **boundary locus method** is:

1. Move $\theta$ from $0$ to $2\pi$ and plot the points $z(\theta)$.
2. Connect the points. The region's boundary is part of this curve.
3. Pick one test point on each side. Compute the roots of $\pi(r)$ there and check the root condition. This tells us which side is inside.

For forward Euler, $\rho(r) = r - 1$ and $\sigma(r) = 1$, so $z(\theta) = e^{i\theta} - 1$: the circle around $-1$ again.

### Why nobody uses Milne's method

For Milne's method, $\rho(r) = r^2 - 1$ and $\sigma(r) = \frac{1}{3}(r^2 + 4r + 1)$. The roots of $\rho$ are $\pm 1$, both simple and on the unit circle. So the method passes the root condition: it is zero-stable and it converges.

Now the boundary locus. Divide the top and the bottom by $e^{i\theta}$:

$$
z(\theta) = \frac{3\,(e^{i\theta} - e^{-i\theta})}{e^{i\theta} + 4 + e^{-i\theta}} = \frac{3i\,\sin\theta}{2 + \cos\theta}.
$$

This is purely imaginary for every $\theta$. I computed it at $\theta = 0°, 30°, \dots, 330°$ and found points like $0.523i$, $1.039i$, $1.5i$, $1.732i$. All have real part zero, and the largest is $\sqrt{3}\,i$. So the "boundary" is the segment from $-\sqrt{3}\,i$ to $\sqrt{3}\,i$ on the imaginary axis, and the region has no inside at all.

Check with test points on the negative real axis:

| $z$ | $\lvert r \rvert$ of the two roots of $\pi$ |
|---|---|
| $-0.5$ | $1.178$, $0.606$ |
| $-0.1$ | $1.034$, $0.905$ |
| $-0.01$ | $1.0033$, $0.9900$ |
| $-0.001$ | $1.00033$, $0.99900$ |

One root is always slightly outside the circle. It is the root near $-1$, which comes from the second root of $\rho$. For $y' = \lambda y$ with $\lambda < 0$, no step size is stable. For small $z$ the bad root is close to $-(1 - z/3)$, so its size is about $1 + |z|/3$. A tiny error grows like $e^{-\lambda x/3}$ and changes sign every step.

I ran Milne's method on $y' = -y$ with $h = 0.1$ and exact starting values. The exact solution decays like $e^{-x}$, but the error grows:

| $x$ | 5 | 10 | 15 | 20 | 25 | 30 |
|---|---|---|---|---|---|---|
| error | $1.6 \times 10^{-7}$ | $7.2 \times 10^{-7}$ | $3.8 \times 10^{-6}$ | $2.0 \times 10^{-5}$ | $1.1 \times 10^{-4}$ | $5.7 \times 10^{-4}$ |

The error grows by about $5.3$ every $5$ units of $x$, and $e^{5/3} = 5.29$. At $x = 30$ the exact value is $9 \times 10^{-14}$, but Milne gives $-5.7 \times 10^{-4}$, with the wrong sign. This is why a method with a slightly larger error constant but better stability wins.

## 9. Stability regions of one-step methods

For a one-step method on $y' = \lambda y$, each step is a multiplication:

$$
y_{n+1} = P(z)\, y_n, \qquad z = h\lambda.
$$

The exact solution gives $y(x_{n+1}) = e^{z}\, y(x_n)$, so a consistent method has $P(z) \approx e^z$. The region of absolute stability is simply

$$
\{ z \in \mathbb{C} : |P(z)| \le 1 \},
$$

because then $|y_n| \le |P(z)|^n |y_0|$ stays bounded for all $n$.

**RK4.** Put $f(x, y) = \lambda y$ into the four stages. Each stage becomes a multiple of $y_n$:

$$
Y_2 = \Bigl(1 + \frac{z}{2}\Bigr) y_n, \quad
Y_3 = \Bigl(1 + \frac{z}{2} + \frac{z^2}{4}\Bigr) y_n, \quad
Y_4 = \Bigl(1 + z + \frac{z^2}{2} + \frac{z^3}{4}\Bigr) y_n.
$$

Then

$$
y_{n+1} = y_n + \frac{z}{6}\bigl[\,Y_1 + 2Y_2 + 2Y_3 + Y_4\,\bigr]
= \Bigl(1 + z + \frac{z^2}{2} + \frac{z^3}{6} + \frac{z^4}{24}\Bigr) y_n.
$$

So $P(z)$ is exactly the Taylor polynomial of $e^z$ of degree $4$. This is one more sign of how accurate RK4 is. In the same way, RK2 gives $P(z) = 1 + z + \frac{z^2}{2}$.

To draw $|P(z)| \le 1$, evaluate $|P|$ on a fine grid of points in the complex plane and draw the level curve $|P| = 1$ (a contour plot). Then test one point to find the inside. By bisection on the real axis, I found these limits:

| method | real interval | on the imaginary axis |
|---|---|---|
| forward Euler | $[-2, 0]$ | only $z = 0$ |
| RK2 | $[-2, 0]$ | only $z = 0$ |
| RK4 | $[-2.785, 0]$ | from $-2.828i$ to $2.828i$ ($= \pm 2\sqrt{2}\,i$) |

The RK4 region is larger than the forward Euler disk, but it is still bounded. Every explicit method here has a bounded region. So when $|\lambda|$ is large, $h$ is bound to about $2.8/|\lambda|$ or less, no matter how smooth the solution is.

## 10. Experiment: four methods on a stiff problem

I solved the stiff problem from Section 6,

$$
y' = -1000\,(y - \sin x) + \cos x, \qquad y(0) = 0, \qquad x \in [0, 1],
$$

with forward Euler, RK4, backward Euler and the trapezoidal method, for $h = 1/50$ down to $h = 1/8000$. The problem is linear, so I solved the implicit steps exactly, with no iteration. I measured the maximum error over the grid.

The theory predicts the stability limits: forward Euler needs $z = -1000h \ge -2$, so $h \le 0.002$ ($n \ge 500$ steps). RK4 needs $z \ge -2.785$, so $h \le 0.002785$ ($n \ge 360$ steps). The two implicit methods should be fine for every $h$.

<figure class="viz">
<div class="viz-legend"><span class="viz-key"><span class="viz-swatch" style="background:#2a78d6"></span>forward Euler</span><span class="viz-key"><span class="viz-swatch" style="background:#eb6834"></span>RK4</span><span class="viz-key"><span class="viz-swatch" style="background:#1baf7a"></span>backward Euler</span><span class="viz-key"><span class="viz-swatch" style="background:#eda100"></span>trapezoidal</span></div>
<div class="viz-frame">
<svg viewBox="0 0 680 400" role="img" tabindex="0" aria-label="Maximum error of four methods for a stiff problem with lambda = minus 1000 on [0, 1] against the step size, both on log scales. Forward Euler is off the top of the chart until h is below 0.002, and RK4 until h is below about 0.0028. Then both errors drop quickly. Backward Euler and the trapezoidal method give small errors for every h." class="viz-plot"><line x1="64" x2="540" y1="30.0" y2="30.0" class="viz-grid"/><text x="54" y="34.0" text-anchor="end" class="viz-tick">10<tspan dy="-0.5em" font-size="75%">2</tspan></text><line x1="64" x2="540" y1="74.3" y2="74.3" class="viz-grid"/><text x="54" y="78.3" text-anchor="end" class="viz-tick">1</text><line x1="64" x2="540" y1="118.6" y2="118.6" class="viz-grid"/><text x="54" y="122.6" text-anchor="end" class="viz-tick">10<tspan dy="-0.5em" font-size="75%">−2</tspan></text><line x1="64" x2="540" y1="162.9" y2="162.9" class="viz-grid"/><text x="54" y="166.9" text-anchor="end" class="viz-tick">10<tspan dy="-0.5em" font-size="75%">−4</tspan></text><line x1="64" x2="540" y1="207.1" y2="207.1" class="viz-grid"/><text x="54" y="211.1" text-anchor="end" class="viz-tick">10<tspan dy="-0.5em" font-size="75%">−6</tspan></text><line x1="64" x2="540" y1="251.4" y2="251.4" class="viz-grid"/><text x="54" y="255.4" text-anchor="end" class="viz-tick">10<tspan dy="-0.5em" font-size="75%">−8</tspan></text><line x1="64" x2="540" y1="295.7" y2="295.7" class="viz-grid"/><text x="54" y="299.7" text-anchor="end" class="viz-tick">10<tspan dy="-0.5em" font-size="75%">−10</tspan></text><line x1="64" x2="540" y1="340.0" y2="340.0" class="viz-grid"/><text x="54" y="344.0" text-anchor="end" class="viz-tick">10<tspan dy="-0.5em" font-size="75%">−12</tspan></text><text x="64.0" y="362" text-anchor="middle" class="viz-tick">0.02</text><text x="129.0" y="362" text-anchor="middle" class="viz-tick">0.01</text><text x="194.0" y="362" text-anchor="middle" class="viz-tick">0.005</text><text x="345.0" y="362" text-anchor="middle" class="viz-tick">0.001</text><text x="410.0" y="362" text-anchor="middle" class="viz-tick">0.0005</text><text x="540.0" y="362" text-anchor="middle" class="viz-tick">0.000125</text><line x1="280.0" x2="280.0" y1="30" y2="340" stroke="#86868b" stroke-width="1" stroke-dasharray="3 4"/><line x1="248.9" x2="248.9" y1="30" y2="340" stroke="#86868b" stroke-width="1" stroke-dasharray="3 4"/><text x="302.0" y="388" text-anchor="middle" class="viz-axis-title">Step size h (smaller to the right)</text><text x="64" y="14" class="viz-axis-title">|error|</text><path d="M64.0,30.0 L129.0,30.0 L194.0,30.0 L232.0,30.0 L243.8,30.0 L254.2,30.0 L270.1,30.0 L278.1,95.8 L281.8,209.0 L311.5,212.0 L345.0,215.5 L410.0,222.1 L475.0,228.8 L540.0,235.5" fill="none" stroke="#2a78d6" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><circle cx="540.0" cy="235.5" r="4" fill="#2a78d6" stroke="#ffffff" stroke-width="2"/><path d="M64.0,30.0 L129.0,30.0 L194.0,30.0 L232.0,30.0 L243.8,30.0 L254.2,191.2 L270.1,208.8 L278.1,214.3 L281.8,216.6 L311.5,232.4 L345.0,248.2 L410.0,277.1 L475.0,304.8 L540.0,332.1" fill="none" stroke="#eb6834" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><circle cx="540.0" cy="332.1" r="4" fill="#eb6834" stroke="#ffffff" stroke-width="2"/><path d="M64.0,186.7 L129.0,193.4 L194.0,200.0 L232.0,203.9 L243.8,205.1 L254.2,206.2 L270.1,207.8 L278.1,208.6 L281.8,209.0 L311.5,212.0 L345.0,215.5 L410.0,222.1 L475.0,228.8 L540.0,235.5" fill="none" stroke="#1baf7a" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><circle cx="540.0" cy="235.5" r="4" fill="#1baf7a" stroke="#ffffff" stroke-width="2"/><path d="M64.0,234.1 L129.0,248.3 L194.0,263.1 L232.0,272.2 L243.8,275.0 L254.2,277.6 L270.1,281.6 L278.1,283.7 L281.8,284.5 L311.5,290.6 L345.0,297.5 L410.0,310.8 L475.0,324.1 L540.0,337.5" fill="none" stroke="#eda100" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><circle cx="540.0" cy="337.5" r="4" fill="#eda100" stroke="#ffffff" stroke-width="2"/><text x="550.0" y="233.5" class="viz-label">FE</text><text x="550.0" y="247.5" class="viz-label">BE</text><text x="416.0" y="269.1" class="viz-label">RK4</text><text x="550.0" y="347.5" class="viz-label">trap.</text><line x1="0" x2="0" y1="30" y2="340" class="viz-crosshair" visibility="hidden"/></svg>
<div class="viz-tooltip" hidden></div>
</div>
<figcaption>Maximum error over the grid for y′ = −1000(y − sin x) + cos x, y(0) = 0 on [0, 1], against the step size h, both on log scales. Errors above 10² are drawn on the top edge (the table has the real values). The two dotted vertical lines are the stability limits h = 2/1000 for forward Euler and h = 2.785/1000 for RK4. Hover over the chart, or focus it and use the arrow keys, to read the values.</figcaption>
<details class="viz-table"><summary>Show the numbers</summary><table><thead><tr><th>h</th><th>forward Euler</th><th>RK4</th><th>backward Euler</th><th>trapezoidal</th></tr></thead><tbody><tr><td>1/50</td><td>4.9e56</td><td>2.6e179</td><td>8.4e-6</td><td>6.1e-8</td></tr><tr><td>1/100</td><td>3.1e87</td><td>1.8e238</td><td>4.2e-6</td><td>1.4e-8</td></tr><tr><td>1/200</td><td>4.3e111</td><td>1.0e219</td><td>2.1e-6</td><td>3.0e-9</td></tr><tr><td>1/300</td><td>4.6e100</td><td>2.9e94</td><td>1.4e-6</td><td>1.2e-9</td></tr><tr><td>1/340</td><td>2.5e87</td><td>2.2e27</td><td>1.2e-6</td><td>8.6e-10</td></tr><tr><td>1/380</td><td>1.0e71</td><td>5.3e-6</td><td>1.1e-6</td><td>6.6e-10</td></tr><tr><td>1/450</td><td>4.8e29</td><td>8.4e-7</td><td>9.3e-7</td><td>4.3e-10</td></tr><tr><td>1/490</td><td>1.1e-1</td><td>4.8e-7</td><td>8.6e-7</td><td>3.5e-10</td></tr><tr><td>1/510</td><td>8.2e-7</td><td>3.7e-7</td><td>8.2e-7</td><td>3.2e-10</td></tr><tr><td>1/700</td><td>6.0e-7</td><td>7.2e-8</td><td>6.0e-7</td><td>1.7e-10</td></tr><tr><td>1/1000</td><td>4.2e-7</td><td>1.4e-8</td><td>4.2e-7</td><td>8.3e-11</td></tr><tr><td>1/2000</td><td>2.1e-7</td><td>7.0e-10</td><td>2.1e-7</td><td>2.1e-11</td></tr><tr><td>1/4000</td><td>1.1e-7</td><td>3.9e-11</td><td>1.1e-7</td><td>5.2e-12</td></tr><tr><td>1/8000</td><td>5.3e-8</td><td>2.3e-12</td><td>5.3e-8</td><td>1.3e-12</td></tr></tbody></table></details>
<script>(() => {const fig = document.currentScript.closest('figure');const svg = fig.querySelector('svg');const line = svg.querySelector('.viz-crosshair');const tip = fig.querySelector('.viz-tooltip');const names = ["forward Euler", "RK4", "backward Euler", "trapezoidal"]; const colors = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"];const rows = [["4.9e56", "2.6e179", "8.4e-6", "6.1e-8"], ["3.1e87", "1.8e238", "4.2e-6", "1.4e-8"], ["4.3e111", "1.0e219", "2.1e-6", "3.0e-9"], ["4.6e100", "2.9e94", "1.4e-6", "1.2e-9"], ["2.5e87", "2.2e27", "1.2e-6", "8.6e-10"], ["1.0e71", "5.3e-6", "1.1e-6", "6.6e-10"], ["4.8e29", "8.4e-7", "9.3e-7", "4.3e-10"], ["1.1e-1", "4.8e-7", "8.6e-7", "3.5e-10"], ["8.2e-7", "3.7e-7", "8.2e-7", "3.2e-10"], ["6.0e-7", "7.2e-8", "6.0e-7", "1.7e-10"], ["4.2e-7", "1.4e-8", "4.2e-7", "8.3e-11"], ["2.1e-7", "7.0e-10", "2.1e-7", "2.1e-11"], ["1.1e-7", "3.9e-11", "1.1e-7", "5.2e-12"], ["5.3e-8", "2.3e-12", "5.3e-8", "1.3e-12"]]; const degrees = [50, 100, 200, 300, 340, 380, 450, 490, 510, 700, 1000, 2000, 4000, 8000];const left = 64, right = 540, width = 680, last = degrees.length - 1;const xs = [64.0, 129.0, 194.0, 232.0, 243.8, 254.2, 270.1, 278.1, 281.8, 311.5, 345.0, 410.0, 475.0, 540.0];let current = -1;const show = (k) => { k = Math.max(0, Math.min(last, k)); current = k; const x = xs[k]; line.setAttribute('x1', x); line.setAttribute('x2', x); line.setAttribute('visibility', 'visible'); tip.replaceChildren(); const head = document.createElement('div'); head.className = 'viz-tooltip-head'; head.textContent = 'h = 1/' + degrees[k] + ', z = hλ = ' + (-1000 / degrees[k]).toFixed(2); tip.append(head); names.forEach((name, i) => {  const row = document.createElement('div'); row.className = 'viz-tooltip-row';  const key = document.createElement('span'); key.className = 'viz-swatch'; key.style.background = colors[i];  const value = document.createElement('strong'); value.textContent = rows[k][i];  const label = document.createElement('span'); label.textContent = name;  row.append(key, value, label); tip.append(row); }); const box = svg.getBoundingClientRect(), frame = fig.getBoundingClientRect(); const xp = box.left - frame.left + x * box.width / width; tip.hidden = false; const room = frame.width - tip.offsetWidth - 8; tip.style.left = Math.max(0, Math.min(room, xp + 12)) + 'px'; tip.style.top = (box.top - frame.top + 8) + 'px';};const hide = () => { line.setAttribute('visibility', 'hidden'); tip.hidden = true; current = -1; };svg.addEventListener('pointermove', (event) => { const box = svg.getBoundingClientRect(); const x = (event.clientX - box.left) * width / box.width; let best = 0; xs.forEach((v, i) => { if (Math.abs(v - x) < Math.abs(xs[best] - x)) best = i; }); show(best);});svg.addEventListener('pointerleave', hide);svg.addEventListener('focus', () => show(current < 0 ? 0 : current));svg.addEventListener('blur', hide);svg.addEventListener('keydown', (event) => { if (event.key === 'ArrowRight') { show(current + 1); event.preventDefault(); } if (event.key === 'ArrowLeft') { show(current - 1); event.preventDefault(); } if (event.key === 'Escape') hide();});})();</script>
</figure>

| $n$ | $h$ | $z = h\lambda$ | forward Euler | RK4 | backward Euler | trapezoidal |
|---|---|---|---|---|---|---|
| 50 | 0.02 | $-20$ | $4.9 \times 10^{56}$ | $2.6 \times 10^{179}$ | $8.4 \times 10^{-6}$ | $6.1 \times 10^{-8}$ |
| 100 | 0.01 | $-10$ | $3.1 \times 10^{87}$ | $1.8 \times 10^{238}$ | $4.2 \times 10^{-6}$ | $1.4 \times 10^{-8}$ |
| 200 | 0.005 | $-5$ | $4.3 \times 10^{111}$ | $1.0 \times 10^{219}$ | $2.1 \times 10^{-6}$ | $3.0 \times 10^{-9}$ |
| 300 | 0.00333 | $-3.33$ | $4.6 \times 10^{100}$ | $2.9 \times 10^{94}$ | $1.4 \times 10^{-6}$ | $1.2 \times 10^{-9}$ |
| 340 | 0.00294 | $-2.94$ | $2.5 \times 10^{87}$ | $2.2 \times 10^{27}$ | $1.2 \times 10^{-6}$ | $8.6 \times 10^{-10}$ |
| 380 | 0.00263 | $-2.63$ | $1.0 \times 10^{71}$ | $5.3 \times 10^{-6}$ | $1.1 \times 10^{-6}$ | $6.6 \times 10^{-10}$ |
| 450 | 0.00222 | $-2.22$ | $4.8 \times 10^{29}$ | $8.4 \times 10^{-7}$ | $9.3 \times 10^{-7}$ | $4.3 \times 10^{-10}$ |
| 490 | 0.00204 | $-2.04$ | $1.1 \times 10^{-1}$ | $4.8 \times 10^{-7}$ | $8.6 \times 10^{-7}$ | $3.5 \times 10^{-10}$ |
| 510 | 0.00196 | $-1.96$ | $8.2 \times 10^{-7}$ | $3.7 \times 10^{-7}$ | $8.2 \times 10^{-7}$ | $3.2 \times 10^{-10}$ |
| 700 | 0.00143 | $-1.43$ | $6.0 \times 10^{-7}$ | $7.2 \times 10^{-8}$ | $6.0 \times 10^{-7}$ | $1.7 \times 10^{-10}$ |
| 1000 | 0.001 | $-1$ | $4.2 \times 10^{-7}$ | $1.4 \times 10^{-8}$ | $4.2 \times 10^{-7}$ | $8.3 \times 10^{-11}$ |
| 2000 | 0.0005 | $-0.5$ | $2.1 \times 10^{-7}$ | $7.0 \times 10^{-10}$ | $2.1 \times 10^{-7}$ | $2.1 \times 10^{-11}$ |
| 4000 | 0.00025 | $-0.25$ | $1.1 \times 10^{-7}$ | $3.9 \times 10^{-11}$ | $1.1 \times 10^{-7}$ | $5.2 \times 10^{-12}$ |
| 8000 | 0.000125 | $-0.125$ | $5.3 \times 10^{-8}$ | $2.3 \times 10^{-12}$ | $5.3 \times 10^{-8}$ | $1.3 \times 10^{-12}$ |

What the numbers say:

- **The limits are sharp.** RK4 explodes at $z = -2.94$ and works at $z = -2.63$, on the two sides of $-2.785$. Forward Euler explodes at $z = -2.22$, is still bad at $z = -2.04$ ($|1 + z| = 1.04$, and $1.04^{490}$ is about $10^8$), and works at $z = -1.96$. The jump from $10^{29}$ to $10^{-6}$ happens over a change in $h$ of less than $15\%$.
- **The implicit methods never fail.** Backward Euler and the trapezoidal method give small errors even at $h = 0.02$, ten times beyond the forward Euler limit.
- **Once stable, the errors are small and simple.** For forward and backward Euler the error is close to $\frac{h}{2}\,\frac{\max|y''|}{|\lambda|}$. For example, at $h = 0.001$: $\frac{0.001}{2} \cdot \frac{\sin 1}{1000} = 4.2 \times 10^{-7}$, which is the value in the table. The strong damping kills old errors, so they do not pile up. The trapezoidal error matches $\frac{h^2}{12}\,\frac{\max|y'''|}{|\lambda|}$ in the same way: $\frac{10^{-6}}{12 \cdot 1000} = 8.3 \times 10^{-11}$ at $h = 0.001$.
- **The trapezoidal method beats RK4** at every $h$ in the table, even though RK4 is fourth order and the trapezoidal method is second order. At $z = -1$ or $-0.5$, $z$ is not small, and RK4 is not yet in its $h^4$ regime. From $h = 0.001$ on, the RK4 errors drop by $16$ to $20$ per halving of $h$, so RK4 wins only at even smaller $h$.

**Part B: the order of RK2 and RK4.** On a non-stiff problem, $y' = y\cos x$, $y(0) = 1$ (exact $e^{\sin x}$), on $[0, 2]$:

| $h$ | RK2 error | ratio | RK4 error | ratio |
|---|---|---|---|---|
| 0.2 | $5.00 \times 10^{-3}$ | | $1.85 \times 10^{-5}$ | |
| 0.1 | $1.36 \times 10^{-3}$ | 3.67 | $1.14 \times 10^{-6}$ | 16.2 |
| 0.05 | $3.58 \times 10^{-4}$ | 3.81 | $7.06 \times 10^{-8}$ | 16.2 |
| 0.025 | $9.20 \times 10^{-5}$ | 3.89 | $4.38 \times 10^{-9}$ | 16.1 |
| 0.0125 | $2.33 \times 10^{-5}$ | 3.95 | $2.73 \times 10^{-10}$ | 16.1 |

Halving $h$ divides the RK2 error by about $4$ and the RK4 error by about $16$, as the theory says.

Here is the core of the stiff test, in Python:

```python
LAM = -1000.0

def f(x, y):
    return LAM * (y - math.sin(x)) + math.cos(x)

def g(x):  # y' = LAM * y + g(x)
    return -LAM * math.sin(x) + math.cos(x)

def step(method, x, y, h):
    if method == "forward_euler":
        return y + h * f(x, y)
    if method == "rk4":
        k1 = f(x, y)
        k2 = f(x + h / 2, y + h / 2 * k1)
        k3 = f(x + h / 2, y + h / 2 * k2)
        k4 = f(x + h, y + h * k3)
        return y + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    if method == "backward_euler":  # linear: solve for y_{n+1} exactly
        return (y + h * g(x + h)) / (1 - h * LAM)
    if method == "trapezoidal":
        return ((1 + h * LAM / 2) * y + h / 2 * (g(x) + g(x + h))) / (1 - h * LAM / 2)
```

## 11. So why use an implicit method?

Three answers came up this week:

1. **Smaller error constants.** At the same order, implicit multistep methods usually have a much smaller truncation error constant than explicit ones.
2. **Stability.** This is the main reason. Forward Euler and RK4 have bounded stability regions, so on a stiff problem $h$ must be below about $2/|\lambda|$. Backward Euler and the trapezoidal method are stable in the whole left half plane, so $h$ can be chosen for accuracy alone.
3. **The price.** Every implicit step needs the solution of an equation for $y_{n+1}$, often nonlinear (Picard iteration from Week 5, or Newton's method from Week 1).

On a stiff problem, one expensive implicit step can replace thousands of cheap explicit steps.

## Summary

| topic | main idea |
|---|---|
| Test problem | Apply the method to $y' = \lambda y$ and perturb the start by $\epsilon$. The perturbation is multiplied by a factor every step. |
| Trapezoidal stability | $r = \frac{1 + z/2}{1 - z/2}$, $\lvert r \rvert < 1 \iff \mathrm{Re}(z) < 0$ for every $h$: A-stable. Stable is not the same as accurate. |
| One-step methods | Forward Euler (explicit, order 1), backward Euler (implicit, order 1), trapezoidal (implicit, order 2). Explicit is easy, implicit needs an equation solve. |
| Taylor series methods | Any order, but need derivatives of $f$. Not practical. |
| Runge–Kutta | Several stages per step, no derivatives. RK2: Euler half step + midpoint rule. RK4: Simpson's rule with two midpoint guesses, order 4. |
| Multistep methods | Built by undetermined coefficients or by integration rules. Milne = Simpson's rule on two steps. |
| One-step vs. multistep | Self-starting and easy to change $h$ vs. reuse of old $f$ values. |
| Root condition | Roots of $\rho(r)$ in the closed unit disk, simple on the circle. With consistency: convergent $\iff$ root condition $\iff$ stable. |
| Fixed $h$ | Use $\pi(r) = \rho(r) - z\sigma(r)$, $z = h\lambda$. Zero-stability ($h \to 0$) vs. absolute stability (fixed $h$). |
| Stability regions | Forward Euler: disk $\lvert 1 + z \rvert \le 1$. Backward Euler: outside $\lvert 1 - z \rvert < 1$. Trapezoidal: left half plane. RK4: $\lvert P(z) \rvert \le 1$, real interval $[-2.785, 0]$. |
| Boundary locus | Plot $z(\theta) = \rho(e^{i\theta})/\sigma(e^{i\theta})$, then use a test point. Milne's region has no inside. |
| Stiff problems | A strongly damped component forces explicit methods to take $h \lesssim 2/\lvert\lambda\rvert$. |

The main thing I take from this week: convergence as $h \to 0$ is not enough to make a method useful. The root condition and the main theorem tell us that a method works in the limit. But we always compute with one fixed $h$, and then what matters is where $h\lambda$ falls in the complex plane. In the experiment, a change in $h$ of less than $15\%$ turned an error of $10^{29}$ into an error of $10^{-6}$. Accuracy tells us how good a method is when it works. Stability tells us whether it works at all.

The full script, with all the numbers and the chart, is in [this site's repository](https://github.com/doubleyands/doubleyands.github.io/tree/main/scripts/posts/numerical-analysis/week-6).
