---
title: "Numerical PDEs, Lecture 15: Godunov's Method, Slopes, and Total Variation"
description: "My notes on the fifteenth lecture: Godunov's method as reconstruct–evolve–average, piecewise linear reconstruction with a slope in each cell, the update formula that follows, four classic slopes, the overshoot next to a jump, total variation, TVD schemes, and the minmod slope."
pubDate: 2026-11-24
project: numerical-pde
tags: [numerical PDEs, finite volume, Godunov's method, slope limiter, TVD, minmod, notes]
---

This is my summary of the fifteenth lecture in my numerical PDE course. I write the ideas in my own words and check them with a small experiment.

In [Lecture 13](/blog/2026-11-17-numerical-pde-lecture-13/) and [Lecture 14](/blog/2026-11-19-numerical-pde-lecture-14/) the key word was the **numerical flux**. This lecture changes the key word to the **slope**. The idea is called **Godunov's method**. We do not guess a flux directly. We first rebuild a function inside each cell, then move it **exactly** with the PDE, then take cell averages again. With a constant function in each cell this gives upwind. With a straight line in each cell it gives a whole family of second-order schemes, and the slope of the line decides which one. The price is new wiggles next to jumps, and the rest of the lecture is about how to measure and stop them.

Notation as before: $\Delta x$ is the cell size, $\Delta t$ is the time step, and $Q^n_i$ is the average over cell $C_i = [x_{i-1/2}, x_{i+1/2}]$ at time $t_n$. The cell center is $x_i$. The model problem is advection $q_t + a q_x = 0$ with a constant speed $a > 0$, and $\nu = a\Delta t/\Delta x$ is the Courant number.

## 1. The starting point is exact

Integrate the conservation law $q_t + f(q)_x = 0$ over the cell and over one time step ([Lecture 13, section 8](/blog/2026-11-17-numerical-pde-lecture-13/)):

$$
Q^{n+1}_i = Q^n_i - \frac{1}{\Delta x}\left[\int_{t_n}^{t_{n+1}} f\big(q(x_{i+1/2}, t)\big)\,dt - \int_{t_n}^{t_{n+1}} f\big(q(x_{i-1/2}, t)\big)\,dt\right].
$$

There is no approximation here yet. The only problem is that the two time integrals need the true solution on the cell walls, and we do not know it. Every finite volume scheme is a way to compute these two integrals. Note also that we only need $q$ **on the walls**: a cell talks to its neighbors only through its interfaces.

## 2. Godunov's method: reconstruct, evolve, average

At the start of each time step we only have one number per cell, the cell average. Godunov's method uses three steps, called **REA**:

1. **Reconstruct.** From the cell averages, build a function $\tilde q^n(x)$: a simple polynomial on each cell.
2. **Evolve.** Use $\tilde q^n$ as initial data and solve the PDE **exactly** from $t_n$ to $t_{n+1}$. This gives the flux integrals in section 1.
3. **Average.** Take cell averages of the result at $t_{n+1}$. Now we are back where we started, one step later.

A useful habit: write down what is given and what is wanted. Both are the same kind of object, cell averages. The only place we lose information is step 1. In step 2 we use the PDE itself, and step 3 is just an integral.

The hard part of step 2 is time. The cell walls are fixed, but $t$ runs over the whole interval $[t_n, t_{n+1}]$, and we only know the state at $t_n$. So we need the PDE to tell us what happens at the walls at every time in between. For advection this is easy, as section 4 shows.

## 3. Reconstruct: a line through each cell average

The simplest choice keeps the cell average as a constant on each cell (piecewise constant). The next choice is a straight line on each cell (piecewise linear):

$$
\tilde q^n(x) = Q^n_i + \sigma^n_i\,(x - x_i) \qquad \text{for } x \in C_i.
$$

The number $\sigma^n_i$ is the **slope** in cell $i$. We have not chosen it yet.

Why write the line as $Q_i + \sigma_i(x - x_i)$ and not as $\alpha x + \beta$? Because $x - x_i$ is an odd function about the cell center, so its integral over the cell is zero:

$$
\frac{1}{\Delta x}\int_{C_i} \tilde q^n(x)\,dx = Q^n_i + \sigma^n_i \cdot \frac{1}{\Delta x}\int_{x_i - \Delta x/2}^{x_i + \Delta x/2} (x - x_i)\,dx = Q^n_i.
$$

So for **any** slope, the line passes through the cell average at the cell center, and the average of the reconstruction is exactly the cell average. The reconstruction does not create or destroy anything. For conservation laws this is the important property: the slope can change the shape inside a cell, but never the amount in it.

## 4. Evolve: follow the characteristics back

For $q_t + aq_x = 0$ with data $q(x, t_n) = \tilde q^n(x)$ at time $t_n$, the solution is the data shifted to the right:

$$
q(x, t) = \tilde q^n\big(x - a(t - t_n)\big), \qquad t_n \le t \le t_{n+1}.
$$

The shift is $t - t_n$, not $t$, because we start at $t_n$. A quick check: at $t = t_n$ we must get the data back, and we do.

Now stand on the right wall $x_{i+1/2}$ at some time $t$ in the step. Follow the characteristic back to time $t_n$ ([Lecture 12, section 1](/blog/2026-11-12-numerical-pde-lecture-12/)). Its foot is at $x_{i+1/2} - a(t - t_n)$, a little to the left of the wall. The value on the wall is the value of the reconstruction at the foot.

- On the **right wall** $x_{i+1/2}$ (the outflow wall of cell $i$), the foot lies in cell $i$, so we read the line of cell $i$.
- On the **left wall** $x_{i-1/2}$ (the inflow wall), the foot lies in cell $i-1$, so we read the line of cell $i-1$. It is the same calculation with the index shifted by one.

For the whole step, the foot moves back by at most $a\Delta t$. It stays inside cell $i$ only if

$$
a\Delta t \le \Delta x \quad\Longleftrightarrow\quad \nu \le 1.
$$

This is the **CFL condition** again ([Lecture 14, section 6](/blog/2026-11-19-numerical-pde-lecture-14/)). If $\nu > 1$, the foot moves past cell $i$ into cell $i-1$, and the formula below is wrong. We assume $0 \le \nu \le 1$ from now on.

## 5. The update formula

On the right wall, $x_{i+1/2} - x_i = \Delta x/2$, so

$$
q(x_{i+1/2}, t) = Q^n_i + \sigma^n_i\left(\frac{\Delta x}{2} - a(t - t_n)\right).
$$

This is a line in $t$, so the flux integral is easy:

$$
\int_{t_n}^{t_{n+1}} a\,q(x_{i+1/2}, t)\,dt = a\Delta t\,Q^n_i + \frac{a\Delta t\,\Delta x}{2}\sigma^n_i - \frac{a^2\Delta t^2}{2}\sigma^n_i = \Delta t\left[aQ^n_i + \frac{a}{2}(1 - \nu)\Delta x\,\sigma^n_i\right].
$$

The bracket is a numerical flux:

$$
F^n_{i+1/2} = aQ^n_i + \frac{a}{2}(1 - \nu)\,\Delta x\,\sigma^n_i.
$$

The left wall gives the same expression with $i - 1$. Put both into section 1:

$$
\boxed{\;Q^{n+1}_i = Q^n_i - \nu\big(Q^n_i - Q^n_{i-1}\big) - \frac{\nu}{2}(1 - \nu)\,\Delta x\,\big(\sigma^n_i - \sigma^n_{i-1}\big).\;}
$$

The first two terms are the upwind scheme from [Lecture 14, section 7](/blog/2026-11-19-numerical-pde-lecture-14/). The last term is the correction from the slopes. The scheme is in conservation form, so the total amount is kept exactly for every choice of slope ([Lecture 13, section 9](/blog/2026-11-17-numerical-pde-lecture-13/)).

The formula is fixed. All the freedom is in $\sigma$, and choosing $\sigma$ well needs both mathematics and physical intuition. A new good slope is a new scheme.

## 6. Four classic slopes

The natural choices use the two neighbors of cell $i$:

| name | slope $\sigma^n_i$ | which information | order for smooth $q$ |
|---|---|---|---|
| upwind | $0$ | none (do nothing) | 1 |
| Fromm | $\dfrac{Q^n_{i+1} - Q^n_{i-1}}{2\Delta x}$ | centered slope | 2 |
| Beam–Warming | $\dfrac{Q^n_i - Q^n_{i-1}}{\Delta x}$ | upwind slope (from the left) | 2 |
| Lax–Wendroff | $\dfrac{Q^n_{i+1} - Q^n_i}{\Delta x}$ | downwind slope (from the right) | 2 |

The names are not random. Put the downwind slope into the update formula. With $\Delta x\,\sigma_i = Q_{i+1} - Q_i$, the slope term becomes $-\tfrac{\nu}{2}(1-\nu)(Q_{i+1} - 2Q_i + Q_{i-1})$, and after collecting terms,

$$
Q^{n+1}_i = Q^n_i - \frac{\nu}{2}\big(Q^n_{i+1} - Q^n_{i-1}\big) + \frac{\nu^2}{2}\big(Q^n_{i+1} - 2Q^n_i + Q^n_{i-1}\big).
$$

This is exactly Lax–Wendroff from [Lecture 12, section 9](/blog/2026-11-12-numerical-pde-lecture-12/). In the same way the upwind slope gives the **Beam–Warming** scheme

$$
Q^{n+1}_i = Q^n_i - \frac{\nu}{2}\big(3Q^n_i - 4Q^n_{i-1} + Q^n_{i-2}\big) + \frac{\nu^2}{2}\big(Q^n_i - 2Q^n_{i-1} + Q^n_{i-2}\big),
$$

which uses only points on the upwind side. So "which slope do we reconstruct with" puts many schemes we know into one frame.

About accuracy: a piecewise constant function approximates a smooth function to first order in $\Delta x$, and a piecewise linear one to second order. A Taylor expansion confirms that the last three slopes give second-order schemes for smooth solutions, one order better than upwind. But for hyperbolic problems smooth solutions are not the main concern. Jumps are. For linear advection a jump in the data simply moves along. For nonlinear laws a jump (a **shock**) can form even from smooth data. So the real question is what the slopes do at a jump.

## 7. One step at a jump: overshoot

Take data that jump from $1$ down to $0$ between cell $i$ and cell $i+1$, and do one REA step with the Lax–Wendroff (downwind) slope.

- **Reconstruct.** Only cell $i$ gets a nonzero slope: $\sigma_i = (0 - 1)/\Delta x < 0$. Every other cell has the same value as its right neighbor, so its line is flat. Everything interesting happens next to the jump.
- **Evolve.** Shift the whole picture right by $a\Delta t = \nu\Delta x$.
- **Average.** After the shift, cell $i$ is covered by a piece of cell $i-1$ (value $1$) and by the left part of its own tilted line. That part of the line lies **above $1$**. So the new average of cell $i$ is larger than $1$.

The new value is bigger than any value in the data. This is an **overshoot**. The true solution never goes above $1$; the extra height is created by the scheme. It is the same kind of effect as the Gibbs phenomenon ([Lecture 7](/blog/2026-10-27-numerical-pde-lecture-7/)): a smooth (here, linear) shape tries to fit a jump and goes past it. If nothing stops it, it can grow step by step.

My own numbers for this picture, with $\nu = 1/4$ (experiment part 2). The exact answer: cell $i$ stays $1$, and cell $i+1$ gets $1/4$, because a quarter of it is now covered by the value $1$.

| slope | cell $i$ | cell $i+1$ | cell $i+2$ | new max / min? |
|---|---|---|---|---|
| exact | $1$ | $0.25$ | $0$ | no |
| upwind | $1$ | $0.25$ | $0$ | no |
| Lax–Wendroff | $\mathbf{1.09375}$ | $0.15625$ | $0$ | max |
| Beam–Warming | $1$ | $0.34375$ | $\mathbf{-0.09375}$ | min |
| Fromm | $\mathbf{1.046875}$ | $0.25$ | $\mathbf{-0.046875}$ | both |
| minmod (section 10) | $1$ | $0.25$ | $0$ | no |

All three second-order slopes leave the range $[0, 1]$ after one step. Lax–Wendroff overshoots behind the jump, and Beam–Warming undershoots in front of it. (Upwind is exact here only because the data are a single clean jump.)

## 8. Total variation

How do we stop this? The slope is the problem, so we limit the slope near jumps. This is called a **slope limiter**. A closely related idea is a **flux limiter**: mix a low-order flux and a high-order flux and use less of the high-order one near jumps. The two ideas look different, but in the end they do the same job. We cannot always forbid all growth, but we can control it.

To control the wiggles, we need to measure them with a number. The **total variation** of the cell averages is

$$
TV(Q) = \sum_i \big|Q_i - Q_{i-1}\big|.
$$

It adds up all the changes between neighbors. A step down is also a change, so we add absolute values. The more the cell averages go up and down, the larger $TV$ is.

A small example of my own. The data $0, 0, 3, 3$ have $TV = 3$. The data $0, 1, 2, 3$ also have $TV = 3$. In general, if the data are **monotone**, then $TV$ is just $|\text{last value} - \text{first value}|$, no matter how the climb is split. If a wiggle appears, for example $0, 1, 3.5, 3$, then the trip down adds to the sum: $TV = 1 + 2.5 + 0.5 = 4$. In the jump step of section 7, each second-order slope raises $TV$ near the jump from $1$ to $1.1875$.

## 9. TVD and monotonicity-preserving schemes

**Definition.** A scheme is **total variation diminishing (TVD)** if, for any data,

$$
TV(Q^{n+1}) \le TV(Q^n).
$$

Wiggles raise the total variation, so a TVD scheme cannot create new wiggles. This is a clever and safe test. (It is not the only possible test.) The exact solution of advection just moves the data, so it keeps $TV$ constant. Strictly speaking the name should be "total variation non-increasing", because "not growing" is all we need.

**Definition.** A scheme is **monotonicity-preserving** if monotone data stay monotone:

$$
Q^n_i \ge Q^n_{i+1} \text{ for all } i \quad\Longrightarrow\quad Q^{n+1}_i \ge Q^{n+1}_{i+1} \text{ for all } i.
$$

TVD looks at the whole sum; this property looks at the order of neighbors.

**TVD implies monotonicity-preserving.** Suppose $Q^n$ decreases from a left value $q_L$ to a right value $q_R$. By section 8, $TV(Q^n) = q_L - q_R$. The scheme uses only a few neighbors, so the far left and far right values do not change in one step. If $Q^{n+1}$ were not monotone, it would go up somewhere on the way down, and then $TV(Q^{n+1}) > q_L - q_R = TV(Q^n)$. That contradicts TVD.

**Upwind is TVD.** For $0 \le \nu \le 1$ upwind is the convex combination $Q^{n+1}_i = \nu Q^n_{i-1} + (1 - \nu)Q^n_i$ ([Lecture 14, section 7](/blog/2026-11-19-numerical-pde-lecture-14/)). Subtract the same formula at $i - 1$:

$$
Q^{n+1}_i - Q^{n+1}_{i-1} = \nu\big(Q^n_{i-1} - Q^n_{i-2}\big) + (1 - \nu)\big(Q^n_i - Q^n_{i-1}\big).
$$

Take absolute values, use the triangle inequality, and sum over $i$ (on a periodic grid both sums on the right are $TV(Q^n)$). This gives $TV(Q^{n+1}) \le \nu\,TV(Q^n) + (1 - \nu)\,TV(Q^n) = TV(Q^n)$.

What does TVD give us? If the data are monotone, a TVD scheme cannot make the solution oscillate. It can still **smear** a jump: upwind is TVD and smears a lot. We want sharp jumps, but sharpness tends to bring overshoot and undershoot. Getting both is the hard design problem, and it is still an active research topic.

## 10. The minmod slope

The **minmod** slope (minimum modulus) picks the smaller of the two one-sided slopes in size:

$$
\sigma^n_i = \operatorname{minmod}\!\left(\frac{Q^n_i - Q^n_{i-1}}{\Delta x},\ \frac{Q^n_{i+1} - Q^n_i}{\Delta x}\right), \qquad
\operatorname{minmod}(\alpha, \beta) = \begin{cases} \alpha & \text{if } |\alpha| < |\beta| \text{ and } \alpha\beta > 0, \\ \beta & \text{if } |\beta| \le |\alpha| \text{ and } \alpha\beta > 0, \\ 0 & \text{if } \alpha\beta \le 0. \end{cases}
$$

The sign rule matters. If the two slopes have different signs, then $Q_i$ is a local maximum or minimum, and the slope is set to $0$. A flat line at a peak cannot push a value above the peak, so no new maximum or minimum appears.

At the jump of section 7, cell $i$ has upwind slope $0$ and downwind slope $-1/\Delta x$. Minmod takes $0$, so cell $i$ stays flat, and the step gives exactly the right averages. Away from jumps and peaks, both one-sided slopes are close to $q_x$, so minmod is a good slope there too. Minmod is TVD, and it is second order where the solution is smooth. Because its slope depends on the data, the scheme is **nonlinear**, even for linear advection.

At a smooth peak minmod also sets the slope to zero, and the top of the peak gets clipped a little. The experiment shows how much this costs. Other limiters, less strict than minmod, are the topic of the next lecture.

## 11. Experiment

I implemented the update formula of section 5 with all five slopes, for $q_t + q_x = 0$ on the periodic interval $[0, 1)$. The script uses only the Python standard library.

**Part 1: the formula is REA.** I also wrote a literal REA step: build the lines, shift them by $\nu\Delta x$, and integrate them exactly over each new cell. On $50$ cells with a mix of smooth and jump data, the formula and the literal REA step agree to $4.4 \times 10^{-16}$ for all five slopes and for $\nu = 0.3$ and $0.7$. The downwind slope matches the Lax–Wendroff formula and the upwind slope matches the Beam–Warming formula, both to $4.4 \times 10^{-16}$. The total amount $\Delta x \sum_i Q_i$ does not change at all in these runs.

**Part 2: one step at a jump.** These are the numbers in the table of section 7 ($12$ cells, $\nu = 1/4$). The periodic data also have a jump up from $0$ to $1$ at $x = 0$, so the whole $TV$ goes from $2$ to $2.375$ for Fromm, Beam–Warming and Lax–Wendroff ($0.1875$ at each jump), and stays $2$ for upwind and minmod.

**Part 3: smooth data.** I moved $\sin 2\pi x$ (as exact cell averages) once around the domain with $\nu \approx 0.6$ and compared with the start. Error ratios are written as orders $\log_2(e_{\text{coarse}}/e_{\text{fine}})$.

| slope | $L^1$ error, $20$ cells | $L^1$ error, $320$ cells | $L^1$ order (last) | max error, $320$ cells | max order (last) |
|---|---|---|---|---|---|
| upwind | $2.06 \times 10^{-1}$ | $1.55 \times 10^{-2}$ | $0.99$ | $2.43 \times 10^{-2}$ | $0.99$ |
| Fromm | $5.23 \times 10^{-3}$ | $1.04 \times 10^{-5}$ | $2.01$ | $1.63 \times 10^{-5}$ | $2.01$ |
| Beam–Warming | $3.56 \times 10^{-2}$ | $1.44 \times 10^{-4}$ | $2.01$ | $2.26 \times 10^{-4}$ | $2.01$ |
| Lax–Wendroff | $4.08 \times 10^{-2}$ | $1.64 \times 10^{-4}$ | $2.00$ | $2.58 \times 10^{-4}$ | $2.00$ |
| minmod | $5.35 \times 10^{-2}$ | $4.13 \times 10^{-4}$ | $1.88$ | $3.93 \times 10^{-3}$ | $1.29$ |

- Upwind is first order and the three linear slopes are second order, as section 6 says. Fromm has the smallest error, more than ten times smaller than the one-sided slopes.
- Minmod is close to second order in the $L^1$ norm, but only about $1.3$ in the max norm. The largest error sits at the top and bottom of the sine wave, where minmod sets the slope to zero. This is the clipping from section 10. It affects only a few cells, so the average error is still small.

**Part 4: a bump and a box.** The data are a cosine bump (height $1$, center $0.25$, half-width $0.12$) and a box ($1$ on $[0.55, 0.8)$), on $200$ cells. One full period takes $333$ steps with $\nu = 0.6006$. The initial $TV$ is $3.997$.

| slope | bump peak | max | min | $TV$ at the end | steps where $TV$ grew | $L^1$ error |
|---|---|---|---|---|---|---|
| upwind | $0.751$ | $0.995$ | $0.000$ | $3.467$ | $0$ | $0.110$ |
| Fromm | $0.998$ | $1.067$ | $-0.067$ | $4.488$ | $160$ | $0.021$ |
| Beam–Warming | $0.999$ | $1.237$ | $-0.236$ | $6.161$ | $177$ | $0.053$ |
| Lax–Wendroff | $0.998$ | $1.219$ | $-0.219$ | $5.702$ | $181$ | $0.052$ |
| minmod | $0.941$ | $1.000$ | $0.000$ | $3.882$ | $0$ | $0.035$ |

<figure class="viz">
<div class="viz-legend"><span class="viz-key"><span class="viz-swatch" style="background:#2a78d6"></span>upwind</span><span class="viz-key"><span class="viz-swatch" style="background:#eb6834"></span>Lax–Wendroff</span><span class="viz-key"><span class="viz-swatch" style="background:#1baf7a"></span>Beam–Warming</span><span class="viz-key"><span class="viz-swatch" style="background:#eda100"></span>minmod</span><span class="viz-key"><span class="viz-swatch" style="background:#8e8e93"></span>exact</span></div>
<div class="viz-frame">
<svg viewBox="0 0 680 400" role="img" tabindex="0" aria-label="The box after one period, for x from 0.45 to 0.9. Upwind stays between 0 and 1 but its edges are smeared. Lax-Wendroff wiggles behind each jump, up to 1.22 and down to minus 0.22. Beam-Warming wiggles ahead of each jump, up to 1.24 and down to minus 0.24. Minmod stays between 0 and 1 with sharper edges than upwind." class="viz-plot"><line x1="56" x2="560" y1="55.8" y2="55.8" class="viz-grid"/><text x="46" y="59.8" text-anchor="end" class="viz-tick">1.25</text><line x1="56" x2="560" y1="98.9" y2="98.9" class="viz-grid"/><text x="46" y="102.9" text-anchor="end" class="viz-tick">1</text><line x1="56" x2="560" y1="141.9" y2="141.9" class="viz-grid"/><text x="46" y="145.9" text-anchor="end" class="viz-tick">0.75</text><line x1="56" x2="560" y1="185.0" y2="185.0" class="viz-grid"/><text x="46" y="189.0" text-anchor="end" class="viz-tick">0.5</text><line x1="56" x2="560" y1="228.1" y2="228.1" class="viz-grid"/><text x="46" y="232.1" text-anchor="end" class="viz-tick">0.25</text><line x1="56" x2="560" y1="271.1" y2="271.1" class="viz-leader"/><text x="46" y="275.1" text-anchor="end" class="viz-tick">0</text><line x1="56" x2="560" y1="314.2" y2="314.2" class="viz-grid"/><text x="46" y="318.2" text-anchor="end" class="viz-tick">-0.25</text><text x="112.0" y="362" text-anchor="middle" class="viz-tick">0.5</text><text x="224.0" y="362" text-anchor="middle" class="viz-tick">0.6</text><text x="336.0" y="362" text-anchor="middle" class="viz-tick">0.7</text><text x="448.0" y="362" text-anchor="middle" class="viz-tick">0.8</text><text x="560.0" y="362" text-anchor="middle" class="viz-tick">0.9</text><text x="308.0" y="388" text-anchor="middle" class="viz-axis-title">x</text><text x="56" y="14" class="viz-axis-title">Q after one period, 200 cells, ν = 0.6</text><path d="M58.8,268.1 L64.4,267.4 L70.0,266.5 L75.6,265.3 L81.2,263.8 L86.8,261.9 L92.4,259.7 L98.0,257.0 L103.6,254.0 L109.2,250.4 L114.8,246.3 L120.4,241.7 L126.0,236.6 L131.6,231.0 L137.2,224.9 L142.8,218.3 L148.4,211.4 L154.0,204.2 L159.6,196.7 L165.2,189.1 L170.8,181.4 L176.4,173.8 L182.0,166.3 L187.6,159.0 L193.2,152.0 L198.8,145.4 L204.4,139.2 L210.0,133.5 L215.6,128.3 L221.2,123.7 L226.8,119.5 L232.4,115.9 L238.0,112.7 L243.6,110.0 L249.2,107.8 L254.8,105.9 L260.4,104.4 L266.0,103.1 L271.6,102.1 L277.2,101.4 L282.8,100.8 L288.4,100.3 L294.0,100.0 L299.6,99.8 L305.2,99.8 L310.8,99.8 L316.4,99.9 L322.0,100.1 L327.6,100.4 L333.2,100.9 L338.8,101.5 L344.4,102.3 L350.0,103.3 L355.6,104.6 L361.2,106.1 L366.8,108.0 L372.4,110.3 L378.0,112.9 L383.6,116.0 L389.2,119.6 L394.8,123.7 L400.4,128.3 L406.0,133.4 L411.6,139.0 L417.2,145.1 L422.8,151.7 L428.4,158.6 L434.0,165.8 L439.6,173.3 L445.2,180.9 L450.8,188.6 L456.4,196.2 L462.0,203.7 L467.6,211.0 L473.2,218.0 L478.8,224.6 L484.4,230.8 L490.0,236.5 L495.6,241.7 L501.2,246.3 L506.8,250.5 L512.4,254.1 L518.0,257.3 L523.6,260.0 L529.2,262.2 L534.8,264.1 L540.4,265.7 L546.0,266.9 L551.6,267.9 L557.2,268.7" fill="none" stroke="#2a78d6" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><path d="M58.8,271.6 L64.4,271.5 L70.0,271.3 L75.6,271.1 L81.2,270.9 L86.8,270.8 L92.4,270.8 L98.0,270.7 L103.6,270.7 L109.2,270.4 L114.8,269.9 L120.4,269.0 L126.0,267.7 L131.6,265.6 L137.2,262.6 L142.8,258.4 L148.4,252.5 L154.0,244.4 L159.6,233.8 L165.2,220.2 L170.8,203.6 L176.4,184.1 L182.0,162.1 L187.6,138.8 L193.2,115.4 L198.8,93.8 L204.4,75.9 L210.0,63.6 L215.6,58.1 L221.2,59.9 L226.8,68.1 L232.4,81.0 L238.0,95.5 L243.6,108.7 L249.2,117.6 L254.8,120.5 L260.4,117.4 L266.0,109.6 L271.6,100.0 L277.2,91.6 L282.8,87.0 L288.4,87.2 L294.0,91.4 L299.6,97.6 L305.2,103.1 L310.8,105.9 L316.4,105.4 L322.0,102.2 L327.6,98.3 L333.2,95.4 L338.8,94.7 L344.4,96.2 L350.0,98.6 L355.6,100.7 L361.2,101.4 L366.8,100.6 L372.4,99.2 L378.0,98.1 L383.6,98.0 L389.2,98.8 L394.8,100.3 L400.4,101.7 L406.0,103.1 L411.6,104.6 L417.2,107.1 L422.8,111.1 L428.4,117.2 L434.0,125.6 L439.6,136.5 L445.2,150.0 L450.8,166.5 L456.4,185.8 L462.0,207.7 L467.6,231.1 L473.2,254.6 L478.8,276.3 L484.4,294.1 L490.0,306.4 L495.6,311.8 L501.2,310.1 L506.8,301.9 L512.4,289.1 L518.0,274.5 L523.6,261.3 L529.2,252.4 L534.8,249.5 L540.4,252.6 L546.0,260.4 L551.6,270.0 L557.2,278.4" fill="none" stroke="#1baf7a" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><path d="M58.8,270.2 L64.4,263.1 L70.0,256.9 L75.6,253.5 L81.2,254.5 L86.8,260.2 L92.4,269.9 L98.0,282.0 L103.6,294.0 L109.2,303.7 L114.8,308.9 L120.4,308.2 L126.0,301.3 L131.6,288.7 L137.2,271.3 L142.8,250.7 L148.4,228.6 L154.0,206.6 L159.6,185.7 L165.2,166.9 L170.8,150.7 L176.4,137.4 L182.0,126.8 L187.6,118.7 L193.2,112.7 L198.8,108.2 L204.4,105.0 L210.0,102.6 L215.6,101.0 L221.2,99.9 L226.8,99.6 L232.4,99.6 L238.0,99.9 L243.6,100.0 L249.2,99.7 L254.8,98.8 L260.4,97.8 L266.0,97.0 L271.6,97.0 L277.2,98.0 L282.8,99.9 L288.4,101.9 L294.0,103.1 L299.6,102.6 L305.2,100.3 L310.8,96.7 L316.4,92.9 L322.0,90.5 L327.6,90.6 L333.2,93.9 L338.8,99.8 L344.4,106.9 L350.0,113.1 L355.6,116.5 L361.2,115.5 L366.8,109.8 L372.4,100.1 L378.0,88.0 L383.6,76.0 L389.2,66.3 L394.8,61.2 L400.4,61.8 L406.0,68.7 L411.6,81.3 L417.2,98.7 L422.8,119.3 L428.4,141.4 L434.0,163.5 L439.6,184.4 L445.2,203.1 L450.8,219.2 L456.4,232.5 L462.0,243.2 L467.6,251.4 L473.2,257.5 L478.8,262.0 L484.4,265.1 L490.0,267.2 L495.6,268.7 L501.2,269.6 L506.8,270.2 L512.4,270.6 L518.0,270.8 L523.6,270.9 L529.2,271.0 L534.8,271.1 L540.4,271.1 L546.0,271.1 L551.6,271.1 L557.2,271.1" fill="none" stroke="#eb6834" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><path d="M58.8,271.1 L64.4,271.1 L70.0,271.1 L75.6,271.1 L81.2,271.0 L86.8,271.0 L92.4,270.9 L98.0,270.7 L103.6,270.4 L109.2,270.0 L114.8,269.2 L120.4,268.0 L126.0,266.0 L131.6,263.1 L137.2,258.7 L142.8,252.4 L148.4,243.5 L154.0,231.5 L159.6,215.6 L165.2,195.3 L170.8,172.7 L176.4,153.7 L182.0,138.5 L187.6,126.8 L193.2,118.1 L198.8,111.8 L204.4,107.4 L210.0,104.3 L215.6,102.3 L221.2,101.0 L226.8,100.1 L232.4,99.6 L238.0,99.3 L243.6,99.1 L249.2,99.0 L254.8,99.0 L260.4,98.9 L266.0,98.9 L271.6,98.9 L277.2,98.9 L282.8,98.9 L288.4,98.9 L294.0,98.9 L299.6,98.9 L305.2,98.9 L310.8,98.9 L316.4,98.9 L322.0,98.9 L327.6,98.9 L333.2,98.9 L338.8,98.9 L344.4,98.9 L350.0,98.9 L355.6,98.9 L361.2,99.0 L366.8,99.0 L372.4,99.1 L378.0,99.3 L383.6,99.6 L389.2,100.0 L394.8,100.8 L400.4,102.0 L406.0,104.0 L411.6,106.9 L417.2,111.3 L422.8,117.6 L428.4,126.5 L434.0,138.5 L439.6,154.4 L445.2,174.7 L450.8,197.3 L456.4,216.3 L462.0,231.5 L467.6,243.2 L473.2,251.9 L478.8,258.2 L484.4,262.6 L490.0,265.7 L495.6,267.7 L501.2,269.0 L506.8,269.9 L512.4,270.4 L518.0,270.7 L523.6,270.9 L529.2,271.0 L534.8,271.0 L540.4,271.1 L546.0,271.1 L551.6,271.1 L557.2,271.1" fill="none" stroke="#eda100" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><path d="M58.8,271.1 L64.4,271.1 L70.0,271.1 L75.6,271.1 L81.2,271.1 L86.8,271.1 L92.4,271.1 L98.0,271.1 L103.6,271.1 L109.2,271.1 L114.8,271.1 L120.4,271.1 L126.0,271.1 L131.6,271.1 L137.2,271.1 L142.8,271.1 L148.4,271.1 L154.0,271.1 L159.6,271.1 L165.2,271.1 L170.8,98.9 L176.4,98.9 L182.0,98.9 L187.6,98.9 L193.2,98.9 L198.8,98.9 L204.4,98.9 L210.0,98.9 L215.6,98.9 L221.2,98.9 L226.8,98.9 L232.4,98.9 L238.0,98.9 L243.6,98.9 L249.2,98.9 L254.8,98.9 L260.4,98.9 L266.0,98.9 L271.6,98.9 L277.2,98.9 L282.8,98.9 L288.4,98.9 L294.0,98.9 L299.6,98.9 L305.2,98.9 L310.8,98.9 L316.4,98.9 L322.0,98.9 L327.6,98.9 L333.2,98.9 L338.8,98.9 L344.4,98.9 L350.0,98.9 L355.6,98.9 L361.2,98.9 L366.8,98.9 L372.4,98.9 L378.0,98.9 L383.6,98.9 L389.2,98.9 L394.8,98.9 L400.4,98.9 L406.0,98.9 L411.6,98.9 L417.2,98.9 L422.8,98.9 L428.4,98.9 L434.0,98.9 L439.6,98.9 L445.2,98.9 L450.8,271.1 L456.4,271.1 L462.0,271.1 L467.6,271.1 L473.2,271.1 L478.8,271.1 L484.4,271.1 L490.0,271.1 L495.6,271.1 L501.2,271.1 L506.8,271.1 L512.4,271.1 L518.0,271.1 L523.6,271.1 L529.2,271.1 L534.8,271.1 L540.4,271.1 L546.0,271.1 L551.6,271.1 L557.2,271.1" fill="none" stroke="#8e8e93" stroke-width="1.5" stroke-dasharray="4 4"/><text x="61.6" y="224.6" text-anchor="start" class="viz-label">upwind</text><text x="420.0" y="61.0" text-anchor="start" class="viz-label">Lax–Wendroff</text><text x="196.0" y="61.0" text-anchor="end" class="viz-label">Beam–Warming</text><text x="468.2" y="150.6" text-anchor="start" class="viz-label">minmod</text><line x1="0" x2="0" y1="30" y2="340" class="viz-crosshair" visibility="hidden"/></svg>
<div class="viz-tooltip" hidden></div>
</div>
<figcaption>The slope method for q_t + q_x = 0 on a periodic grid with 200 cells and ν = 0.6, after one full period (333 steps). Only the box part of the data is shown; the dashed gray line is the exact answer (the initial cell averages). Hover over the chart, or focus it and use the arrow keys, to read the values. The table shows every second cell.</figcaption>
<details class="viz-table"><summary>Show the numbers</summary><table><thead><tr><th>x</th><th>upwind</th><th>Lax–Wendroff</th><th>Beam–Warming</th><th>minmod</th><th>exact</th></tr></thead><tbody><tr><td>0.453</td><td>0.0175</td><td>0.0053</td><td>-0.0030</td><td>0.0000</td><td>0.0000</td></tr><tr><td>0.463</td><td>0.0270</td><td>0.0825</td><td>-0.0011</td><td>0.0001</td><td>0.0000</td></tr><tr><td>0.473</td><td>0.0427</td><td>0.0965</td><td>0.0013</td><td>0.0004</td><td>0.0000</td></tr><tr><td>0.482</td><td>0.0663</td><td>0.0069</td><td>0.0020</td><td>0.0012</td><td>0.0000</td></tr><tr><td>0.492</td><td>0.0996</td><td>-0.1331</td><td>0.0026</td><td>0.0038</td><td>0.0000</td></tr><tr><td>0.503</td><td>0.1441</td><td>-0.2191</td><td>0.0069</td><td>0.0111</td><td>0.0000</td></tr><tr><td>0.512</td><td>0.2004</td><td>-0.1755</td><td>0.0199</td><td>0.0296</td><td>0.0000</td></tr><tr><td>0.522</td><td>0.2684</td><td>-0.0010</td><td>0.0492</td><td>0.0721</td><td>0.0000</td></tr><tr><td>0.532</td><td>0.3465</td><td>0.2466</td><td>0.1083</td><td>0.1602</td><td>0.0000</td></tr><tr><td>0.542</td><td>0.4319</td><td>0.4961</td><td>0.2168</td><td>0.3225</td><td>0.0000</td></tr><tr><td>0.552</td><td>0.5208</td><td>0.6991</td><td>0.3919</td><td>0.5712</td><td>1.0000</td></tr><tr><td>0.562</td><td>0.6088</td><td>0.8381</td><td>0.6327</td><td>0.7701</td><td>1.0000</td></tr><tr><td>0.573</td><td>0.6917</td><td>0.9200</td><td>0.9042</td><td>0.8884</td><td>1.0000</td></tr><tr><td>0.583</td><td>0.7659</td><td>0.9645</td><td>1.1333</td><td>0.9508</td><td>1.0000</td></tr><tr><td>0.593</td><td>0.8291</td><td>0.9880</td><td>1.2366</td><td>0.9802</td><td>1.0000</td></tr><tr><td>0.603</td><td>0.8803</td><td>0.9961</td><td>1.1786</td><td>0.9927</td><td>1.0000</td></tr><tr><td>0.613</td><td>0.9196</td><td>0.9940</td><td>1.0194</td><td>0.9975</td><td>1.0000</td></tr><tr><td>0.623</td><td>0.9483</td><td>0.9954</td><td>0.8913</td><td>0.9992</td><td>1.0000</td></tr><tr><td>0.633</td><td>0.9682</td><td>1.0066</td><td>0.8927</td><td>0.9998</td><td>1.0000</td></tr><tr><td>0.642</td><td>0.9812</td><td>1.0111</td><td>0.9937</td><td>0.9999</td><td>1.0000</td></tr><tr><td>0.652</td><td>0.9891</td><td>0.9939</td><td>1.0691</td><td>1.0000</td><td>1.0000</td></tr><tr><td>0.662</td><td>0.9933</td><td>0.9755</td><td>1.0437</td><td>1.0000</td><td>1.0000</td></tr><tr><td>0.672</td><td>0.9949</td><td>0.9916</td><td>0.9756</td><td>1.0000</td><td>1.0000</td></tr><tr><td>0.682</td><td>0.9941</td><td>1.0348</td><td>0.9623</td><td>1.0000</td><td>1.0000</td></tr><tr><td>0.693</td><td>0.9910</td><td>1.0480</td><td>1.0037</td><td>1.0000</td><td>1.0000</td></tr><tr><td>0.703</td><td>0.9847</td><td>0.9947</td><td>1.0241</td><td>1.0000</td><td>1.0000</td></tr><tr><td>0.713</td><td>0.9742</td><td>0.9175</td><td>1.0014</td><td>0.9999</td><td>1.0000</td></tr><tr><td>0.723</td><td>0.9579</td><td>0.9035</td><td>0.9857</td><td>0.9996</td><td>1.0000</td></tr><tr><td>0.733</td><td>0.9340</td><td>0.9931</td><td>0.9981</td><td>0.9988</td><td>1.0000</td></tr><tr><td>0.743</td><td>0.9005</td><td>1.1331</td><td>1.0053</td><td>0.9962</td><td>1.0000</td></tr><tr><td>0.753</td><td>0.8560</td><td>1.2191</td><td>0.9921</td><td>0.9889</td><td>1.0000</td></tr><tr><td>0.763</td><td>0.7996</td><td>1.1755</td><td>0.9757</td><td>0.9704</td><td>1.0000</td></tr><tr><td>0.772</td><td>0.7316</td><td>1.0010</td><td>0.9523</td><td>0.9279</td><td>1.0000</td></tr><tr><td>0.782</td><td>0.6535</td><td>0.7533</td><td>0.8937</td><td>0.8398</td><td>1.0000</td></tr><tr><td>0.792</td><td>0.5681</td><td>0.5038</td><td>0.7818</td><td>0.6775</td><td>1.0000</td></tr><tr><td>0.802</td><td>0.4792</td><td>0.3013</td><td>0.6076</td><td>0.4288</td><td>0.0000</td></tr><tr><td>0.812</td><td>0.3912</td><td>0.1622</td><td>0.3681</td><td>0.2299</td><td>0.0000</td></tr><tr><td>0.823</td><td>0.3083</td><td>0.0790</td><td>0.0957</td><td>0.1116</td><td>0.0000</td></tr><tr><td>0.833</td><td>0.2341</td><td>0.0350</td><td>-0.1336</td><td>0.0492</td><td>0.0000</td></tr><tr><td>0.843</td><td>0.1709</td><td>0.0141</td><td>-0.2364</td><td>0.0198</td><td>0.0000</td></tr><tr><td>0.853</td><td>0.1197</td><td>0.0052</td><td>-0.1786</td><td>0.0073</td><td>0.0000</td></tr><tr><td>0.863</td><td>0.0804</td><td>0.0018</td><td>-0.0195</td><td>0.0025</td><td>0.0000</td></tr><tr><td>0.873</td><td>0.0516</td><td>0.0006</td><td>0.1087</td><td>0.0008</td><td>0.0000</td></tr><tr><td>0.883</td><td>0.0317</td><td>0.0002</td><td>0.1073</td><td>0.0002</td><td>0.0000</td></tr><tr><td>0.893</td><td>0.0186</td><td>0.0000</td><td>0.0062</td><td>0.0001</td><td>0.0000</td></tr></tbody></table></details>
<script>(() => {const fig = document.currentScript.closest('figure');const svg = fig.querySelector('svg');const line = svg.querySelector('.viz-crosshair');const tip = fig.querySelector('.viz-tooltip');const names = ["upwind", "Lax\u2013Wendroff", "Beam\u2013Warming", "minmod", "exact"]; const colors = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#8e8e93"];const rows = [["0.0175", "0.0053", "-0.0030", "0.0000", "0.0000"], ["0.0216", "0.0463", "-0.0024", "0.0001", "0.0000"], ["0.0270", "0.0825", "-0.0011", "0.0001", "0.0000"], ["0.0339", "0.1020", "0.0002", "0.0002", "0.0000"], ["0.0427", "0.0965", "0.0013", "0.0004", "0.0000"], ["0.0534", "0.0634", "0.0018", "0.0007", "0.0000"], ["0.0663", "0.0069", "0.0020", "0.0012", "0.0000"], ["0.0817", "-0.0631", "0.0021", "0.0022", "0.0000"], ["0.0996", "-0.1331", "0.0026", "0.0038", "0.0000"], ["0.1204", "-0.1891", "0.0040", "0.0066", "0.0000"], ["0.1441", "-0.2191", "0.0069", "0.0111", "0.0000"], ["0.1707", "-0.2154", "0.0120", "0.0184", "0.0000"], ["0.2004", "-0.1755", "0.0199", "0.0296", "0.0000"], ["0.2330", "-0.1019", "0.0319", "0.0467", "0.0000"], ["0.2684", "-0.0010", "0.0492", "0.0721", "0.0000"], ["0.3064", "0.1184", "0.0738", "0.1088", "0.0000"], ["0.3465", "0.2466", "0.1083", "0.1602", "0.0000"], ["0.3885", "0.3749", "0.1551", "0.2302", "0.0000"], ["0.4319", "0.4961", "0.2168", "0.3225", "0.0000"], ["0.4762", "0.6052", "0.2955", "0.4404", "0.0000"], ["0.5208", "0.6991", "0.3919", "0.5712", "1.0000"], ["0.5652", "0.7766", "0.5053", "0.6820", "1.0000"], ["0.6088", "0.8381", "0.6327", "0.7701", "1.0000"], ["0.6511", "0.8851", "0.7684", "0.8379", "1.0000"], ["0.6917", "0.9200", "0.9042", "0.8884", "1.0000"], ["0.7300", "0.9457", "1.0296", "0.9250", "1.0000"], ["0.7659", "0.9645", "1.1333", "0.9508", "1.0000"], ["0.7990", "0.9783", "1.2048", "0.9684", "1.0000"], ["0.8291", "0.9880", "1.2366", "0.9802", "1.0000"], ["0.8562", "0.9939", "1.2264", "0.9878", "1.0000"], ["0.8803", "0.9961", "1.1786", "0.9927", "1.0000"], ["0.9014", "0.9956", "1.1041", "0.9957", "1.0000"], ["0.9196", "0.9940", "1.0194", "0.9975", "1.0000"], ["0.9352", "0.9934", "0.9431", "0.9986", "1.0000"], ["0.9483", "0.9954", "0.8913", "0.9992", "1.0000"], ["0.9593", "1.0003", "0.8743", "0.9996", "1.0000"], ["0.9682", "1.0066", "0.8927", "0.9998", "1.0000"], ["0.9754", "1.0112", "0.9379", "0.9999", "1.0000"], ["0.9812", "1.0111", "0.9937", "0.9999", "1.0000"], ["0.9857", "1.0049", "1.0423", "1.0000", "1.0000"], ["0.9891", "0.9939", "1.0691", "1.0000", "1.0000"], ["0.9916", "0.9823", "1.0682", "1.0000", "1.0000"], ["0.9933", "0.9755", "1.0437", "1.0000", "1.0000"], ["0.9944", "0.9782", "1.0078", "1.0000", "1.0000"], ["0.9949", "0.9916", "0.9756", "1.0000", "1.0000"], ["0.9948", "1.0128", "0.9591", "1.0000", "1.0000"], ["0.9941", "1.0348", "0.9623", "1.0000", "1.0000"], ["0.9929", "1.0489", "0.9807", "1.0000", "1.0000"], ["0.9910", "1.0480", "1.0037", "1.0000", "1.0000"], ["0.9883", "1.0290", "1.0203", "1.0000", "1.0000"], ["0.9847", "0.9947", "1.0241", "1.0000", "1.0000"], ["0.9801", "0.9537", "1.0156", "0.9999", "1.0000"], ["0.9742", "0.9175", "1.0014", "0.9999", "1.0000"], ["0.9669", "0.8980", "0.9896", "0.9998", "1.0000"], ["0.9579", "0.9035", "0.9857", "0.9996", "1.0000"], ["0.9470", "0.9366", "0.9899", "0.9993", "1.0000"], ["0.9340", "0.9931", "0.9981", "0.9988", "1.0000"], ["0.9185", "1.0631", "1.0045", "0.9978", "1.0000"], ["0.9005", "1.1331", "1.0053", "0.9962", "1.0000"], ["0.8797", "1.1891", "1.0003", "0.9934", "1.0000"], ["0.8560", "1.2191", "0.9921", "0.9889", "1.0000"], ["0.8293", "1.2154", "0.9836", "0.9816", "1.0000"], ["0.7996", "1.1755", "0.9757", "0.9704", "1.0000"], ["0.7670", "1.1020", "0.9666", "0.9533", "1.0000"], ["0.7316", "1.0010", "0.9523", "0.9279", "1.0000"], ["0.6936", "0.8816", "0.9289", "0.8912", "1.0000"], ["0.6535", "0.7533", "0.8937", "0.8398", "1.0000"], ["0.6115", "0.6249", "0.8450", "0.7698", "1.0000"], ["0.5681", "0.5038", "0.7818", "0.6775", "1.0000"], ["0.5238", "0.3949", "0.7031", "0.5596", "1.0000"], ["0.4792", "0.3013", "0.6076", "0.4288", "0.0000"], ["0.4348", "0.2239", "0.4952", "0.3180", "0.0000"], ["0.3912", "0.1622", "0.3681", "0.2299", "0.0000"], ["0.3489", "0.1146", "0.2320", "0.1621", "0.0000"], ["0.3083", "0.0790", "0.0957", "0.1116", "0.0000"], ["0.2700", "0.0532", "-0.0300", "0.0750", "0.0000"], ["0.2341", "0.0350", "-0.1336", "0.0492", "0.0000"], ["0.2010", "0.0225", "-0.2048", "0.0316", "0.0000"], ["0.1709", "0.0141", "-0.2364", "0.0198", "0.0000"], ["0.1438", "0.0087", "-0.2263", "0.0122", "0.0000"], ["0.1197", "0.0052", "-0.1786", "0.0073", "0.0000"], ["0.0986", "0.0031", "-0.1042", "0.0043", "0.0000"], ["0.0804", "0.0018", "-0.0195", "0.0025", "0.0000"], ["0.0648", "0.0010", "0.0569", "0.0014", "0.0000"], ["0.0516", "0.0006", "0.1087", "0.0008", "0.0000"], ["0.0407", "0.0003", "0.1258", "0.0004", "0.0000"], ["0.0317", "0.0002", "0.1073", "0.0002", "0.0000"], ["0.0244", "0.0001", "0.0621", "0.0001", "0.0000"], ["0.0186", "0.0000", "0.0062", "0.0001", "0.0000"], ["0.0140", "0.0000", "-0.0423", "0.0000", "0.0000"]]; const heads = ["x = 0.453", "x = 0.458", "x = 0.463", "x = 0.468", "x = 0.473", "x = 0.478", "x = 0.482", "x = 0.487", "x = 0.492", "x = 0.497", "x = 0.503", "x = 0.508", "x = 0.512", "x = 0.517", "x = 0.522", "x = 0.527", "x = 0.532", "x = 0.537", "x = 0.542", "x = 0.547", "x = 0.552", "x = 0.557", "x = 0.562", "x = 0.568", "x = 0.573", "x = 0.578", "x = 0.583", "x = 0.588", "x = 0.593", "x = 0.598", "x = 0.603", "x = 0.608", "x = 0.613", "x = 0.618", "x = 0.623", "x = 0.628", "x = 0.633", "x = 0.638", "x = 0.642", "x = 0.647", "x = 0.652", "x = 0.657", "x = 0.662", "x = 0.667", "x = 0.672", "x = 0.677", "x = 0.682", "x = 0.688", "x = 0.693", "x = 0.698", "x = 0.703", "x = 0.708", "x = 0.713", "x = 0.718", "x = 0.723", "x = 0.728", "x = 0.733", "x = 0.738", "x = 0.743", "x = 0.748", "x = 0.753", "x = 0.758", "x = 0.763", "x = 0.768", "x = 0.772", "x = 0.777", "x = 0.782", "x = 0.787", "x = 0.792", "x = 0.797", "x = 0.802", "x = 0.807", "x = 0.812", "x = 0.818", "x = 0.823", "x = 0.828", "x = 0.833", "x = 0.838", "x = 0.843", "x = 0.848", "x = 0.853", "x = 0.858", "x = 0.863", "x = 0.868", "x = 0.873", "x = 0.878", "x = 0.883", "x = 0.888", "x = 0.893", "x = 0.897"];const xs = [58.8, 64.4, 70.0, 75.6, 81.2, 86.8, 92.4, 98.0, 103.6, 109.2, 114.8, 120.4, 126.0, 131.6, 137.2, 142.8, 148.4, 154.0, 159.6, 165.2, 170.8, 176.4, 182.0, 187.6, 193.2, 198.8, 204.4, 210.0, 215.6, 221.2, 226.8, 232.4, 238.0, 243.6, 249.2, 254.8, 260.4, 266.0, 271.6, 277.2, 282.8, 288.4, 294.0, 299.6, 305.2, 310.8, 316.4, 322.0, 327.6, 333.2, 338.8, 344.4, 350.0, 355.6, 361.2, 366.8, 372.4, 378.0, 383.6, 389.2, 394.8, 400.4, 406.0, 411.6, 417.2, 422.8, 428.4, 434.0, 439.6, 445.2, 450.8, 456.4, 462.0, 467.6, 473.2, 478.8, 484.4, 490.0, 495.6, 501.2, 506.8, 512.4, 518.0, 523.6, 529.2, 534.8, 540.4, 546.0, 551.6, 557.2]; const steps = 89, width = 680;let current = -1;const show = (n) => { n = Math.max(0, Math.min(steps, n)); current = n; const x = xs[n]; line.setAttribute('x1', x); line.setAttribute('x2', x); line.setAttribute('visibility', 'visible'); tip.replaceChildren(); const head = document.createElement('div'); head.className = 'viz-tooltip-head'; head.textContent = heads[n]; tip.append(head); names.forEach((name, i) => {  const row = document.createElement('div'); row.className = 'viz-tooltip-row';  const key = document.createElement('span'); key.className = 'viz-swatch'; key.style.background = colors[i];  const value = document.createElement('strong'); value.textContent = rows[n][i];  const label = document.createElement('span'); label.textContent = name;  row.append(key, value, label); tip.append(row); }); const box = svg.getBoundingClientRect(), frame = fig.getBoundingClientRect(); const xp = box.left - frame.left + x * box.width / width; tip.hidden = false; const room = frame.width - tip.offsetWidth - 8; tip.style.left = Math.max(0, Math.min(room, xp + 12)) + 'px'; tip.style.top = (box.top - frame.top + 8) + 'px';};const hide = () => { line.setAttribute('visibility', 'hidden'); tip.hidden = true; current = -1; };svg.addEventListener('pointermove', (event) => { const box = svg.getBoundingClientRect(); const x = (event.clientX - box.left) * width / box.width; let best = 0; xs.forEach((v, i) => { if (Math.abs(v - x) < Math.abs(xs[best] - x)) best = i; }); show(best);});svg.addEventListener('pointerleave', hide);svg.addEventListener('focus', () => show(current < 0 ? 0 : current));svg.addEventListener('blur', hide);svg.addEventListener('keydown', (event) => { if (event.key === 'ArrowRight') { show(current + 1); event.preventDefault(); } if (event.key === 'ArrowLeft') { show(current - 1); event.preventDefault(); } if (event.key === 'Escape') hide();});})();</script>
</figure>

- **Upwind** never leaves $[0, 1]$ and its $TV$ only goes down, but the bump drops to $0.75$ and the box edges are wide.
- **Lax–Wendroff** and **Beam–Warming** keep the bump almost perfectly, but wiggle at the box by more than $0.2$. Lax–Wendroff wiggles **behind** each jump (its largest value is at $x = 0.75$, just behind the front at $0.8$). Beam–Warming wiggles **in front** of each jump (largest value at $x = 0.59$, just in front of the back edge at $0.55$). This is the same pattern as the one step in section 7. Their $TV$ grows by about $40$–$50\%$.
- **Fromm** has the smallest $L^1$ error, and its wiggles are small ($\pm 0.067$), but it is still not TVD: $TV$ grew in $160$ of the $333$ steps.
- **Minmod** stays exactly in $[0, 1]$, and $TV$ never grew in any step. Its box edges are much sharper than upwind's, and its $L^1$ error is about three times smaller. The cost is the bump peak, $0.94$ instead of $1$, because the peak is a local maximum and gets a zero slope.

All five runs keep the total amount exactly (no change at all in floating point).

The full script and the chart code are in [this site's repository](https://github.com/doubleyands/doubleyands.github.io/tree/main/scripts/posts/numerical-pde/lecture-15).

## Summary

1. **Godunov's method** = reconstruct, evolve, average. Start from cell averages, build a function in each cell, solve the PDE exactly for one step, and average again.
2. A line $Q_i + \sigma_i(x - x_i)$ keeps the cell average for any slope $\sigma_i$, because $x - x_i$ integrates to zero over the cell.
3. For advection, following characteristics back gives the update $Q^{n+1}_i = Q^n_i - \nu(Q^n_i - Q^n_{i-1}) - \tfrac{\nu}{2}(1 - \nu)\Delta x(\sigma^n_i - \sigma^n_{i-1})$, valid for $0 \le \nu \le 1$ (the CFL condition).
4. Slope zero gives upwind (first order). Centered, upwind and downwind slopes give Fromm, Beam–Warming and Lax–Wendroff (second order).
5. Near a jump the second-order slopes create **overshoot** and undershoot. **Total variation** $\sum_i |Q_i - Q_{i-1}|$ measures the wiggles; a **TVD** scheme never increases it, so it is also monotonicity-preserving.
6. The **minmod** slope takes the smaller one-sided slope, and zero at a local extremum. It is TVD, sharper than upwind, and second order where the solution is smooth, but it clips smooth peaks a little.
