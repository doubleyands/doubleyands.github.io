---
title: "Numerical PDEs, Lecture 1: Vector Calculus, Model Equations, and Finite Differences"
description: "My notes on the first lecture: gradient, divergence, curl, the classic PDEs, well-posedness, Green's formula, and a first look at finite differences, tested on four difference formulas."
pubDate: 2026-10-06
project: numerical-pde
tags: [numerical PDEs, vector calculus, finite differences, notes]
---

This semester I am also taking a course in numerical methods for partial differential equations (PDEs). This post is my summary of the first lecture. I write the ideas in my own words and check one of them with a small experiment.

The lecture had two parts. The first part reviews vector calculus: the gradient, the divergence, and the curl. Almost every PDE in the course is written with these three operators. The second part looks at PDEs themselves: the classic equations, what makes a problem "well-posed", Green's formula, and the first idea of a numerical method, the finite difference.

## 1. Directional derivative and gradient

In one dimension there are only two directions, left and right. In two or three dimensions there are infinitely many. So we first pick a direction, a unit vector $\mathbf{d}$, and walk along it. On that path the problem is one-dimensional again, and the slope is

$$
\partial_{\mathbf{d}} f(\mathbf{x}) = \lim_{h \to 0} \frac{f(\mathbf{x} + h\mathbf{d}) - f(\mathbf{x})}{h}.
$$

If it is positive, $f$ grows in that direction. If it is negative, $f$ falls.

The **gradient** is a vector, so it needs a direction and a size. Its direction $\mathbf{d}^*$ is the one with the largest slope, and its size is that slope:

$$
\nabla f(\mathbf{x}) = \big(\partial_{\mathbf{d}^*} f(\mathbf{x})\big)\, \mathbf{d}^*, \qquad \partial_{\mathbf{d}^*} f(\mathbf{x}) = \max_{|\mathbf{d}| = 1} \partial_{\mathbf{d}} f(\mathbf{x}).
$$

So the gradient points in the direction of steepest increase. The usual formula, $\nabla f = (\partial_1 f, \partial_2 f, \partial_3 f)$, is easy to compute with. But the definition above tells us what the vector means.

**The gradient is normal to level sets.** A level set is the set of points where $f$ has one fixed value, like a contour line on a map:

$$
\mathcal{L}_\lambda = \{\mathbf{y} : f(\mathbf{y}) = \lambda\}.
$$

If we walk along a contour line, the height does not change, so the slope along it is zero. The gradient points where the slope is largest, so it must be perpendicular to the contour line. Of the two normal directions, it picks the one where $f$ grows.

## 2. Computing the outward normal

To state a boundary condition on a region $\Omega$, we often need the **unit outward normal** $\mathbf{n}$ on the boundary $\partial\Omega$. The gradient gives it directly. Write the region as $\Omega = \{\mathbf{y} : g(\mathbf{y}) < c\}$, so the boundary is the level set $g = c$. Outside, $g$ is larger, so the gradient points outward, and

$$
\mathbf{n} = \frac{\nabla g}{|\nabla g|}.
$$

**Example.** Take the half-plane $3x + 4y < 12$. Here $g(x, y) = 3x + 4y$ and $\nabla g = (3, 4)$. This vector has the right direction, but its length is $5$, not $1$. So

$$
\mathbf{n} = \tfrac{1}{5}(3, 4) = (0.6,\ 0.8).
$$

Forgetting to divide by the length is an easy mistake. "A normal vector" can have any length. "The unit normal" must have length $1$.

For a curved boundary the steps are the same, but $\mathbf{n}$ changes from point to point. A nice fact: the divergence of the unit normal field is the curvature of the level curve,

$$
\kappa = \nabla \cdot \left( \frac{\nabla \phi}{|\nabla \phi|} \right).
$$

For a circle $\phi = x^2 + y^2 - R^2$, the unit normal is $(x, y)/r$ with $r = \sqrt{x^2 + y^2}$, and its divergence in 2D is $1/r$. On the circle this is $1/R$, which is the curvature of a circle of radius $R$.

## 3. Divergence: two definitions

A good habit with any operator: first check what goes in and what comes out.

| operator | notation | input → output | meaning |
|---|---|---|---|
| gradient | $\nabla f$ | scalar → vector | direction of steepest increase |
| divergence | $\nabla \cdot \mathbf{F}$ | vector → scalar | net outflow at a point |
| curl | $\nabla \times \mathbf{F}$ | vector → vector | how much the field rotates |

The divergence has two common definitions.

**For computing.** Differentiate each component in its own variable and add:

$$
\nabla \cdot \mathbf{F} = \partial_1 F_1 + \partial_2 F_2 + \partial_3 F_3.
$$

**For meaning.** Put a small ball $B_r$ of radius $r$ around the point. Measure the net flux of $\mathbf{F}$ out through its surface, divide by the volume of the ball, and shrink the ball:

$$
\operatorname{div} \mathbf{F}(\mathbf{x}) = \lim_{r \to 0} \frac{1}{|B_r|} \int_{\partial B_r} \mathbf{F} \cdot d\mathbf{S}, \qquad |B_r| = \tfrac{4}{3}\pi r^3.
$$

So the divergence is the outflow per unit volume. Think of a flow meter: we draw a closed surface, count what goes in and what comes out, and learn how much is produced inside. A ball is convenient because it is symmetric. A small cube also works well, because its face normals are just the coordinate directions.

Both definitions give the same thing. You can start from either one and prove the other.

## 4. The divergence theorem

If $\mathbf{F}$ is a $C^1$ vector field ($C^1$ means the function and its first derivatives are continuous), then

$$
\int_\Omega \nabla \cdot \mathbf{F} \, d\mathbf{y} = \int_{\partial\Omega} \mathbf{F} \cdot \mathbf{n} \, dS.
$$

The left side adds up the divergence over the whole region. The right side measures the total flow through the boundary. The integral drops by one dimension: a volume integral becomes a surface integral. To know how many people are inside a room, we do not need to look inside. We only need to count who goes in and out at the door.

The normal points outward because we agree that "out" is positive. This is a sign convention, and it comes from physics.

## 5. Curl, circulation, and Stokes' theorem

The curl takes a vector field and returns a vector field:

$$
\nabla \times \mathbf{F} = \begin{vmatrix} \mathbf{e}_1 & \mathbf{e}_2 & \mathbf{e}_3 \\ \partial_1 & \partial_2 & \partial_3 \\ F_1 & F_2 & F_3 \end{vmatrix}.
$$

Its meaning comes from **circulation**, the line integral of $\mathbf{F}$ around a closed curve $C$:

$$
\oint_C \mathbf{F} \cdot d\mathbf{l}.
$$

Take a tiny flat loop, divide its circulation by its area, and shrink it. Turn the loop until this value is largest. That largest value is the size of the curl, and the normal of the loop (by the right-hand rule) is its direction. This is the same idea as the divergence: there we divided an outflow by a volume, and here we divide a circulation by an area.

**Stokes' theorem** has the same shape as the divergence theorem. For a surface $S$ with boundary curve $C$,

$$
\int_S (\nabla \times \mathbf{F}) \cdot d\mathbf{S} = \oint_C \mathbf{F} \cdot d\mathbf{l}.
$$

The curl is central in electromagnetics. In this course the divergence comes up much more often.

## 6. Useful identities

Combining the three operators gives some identities we will use again:

$$
\Delta u = \nabla \cdot (\nabla u),
$$

$$
\nabla \cdot (\varphi \mathbf{V}) = \nabla \varphi \cdot \mathbf{V} + \varphi \, \nabla \cdot \mathbf{V},
$$

$$
\nabla \times (\nabla \times \mathbf{V}) = \nabla (\nabla \cdot \mathbf{V}) - \Delta \mathbf{V},
$$

$$
\nabla \cdot (\nabla \times \mathbf{V}) = 0, \qquad \nabla \times (\nabla \varphi) = \mathbf{0}.
$$

The first one defines the Laplacian: scalar to vector, then vector back to scalar. The second one is the product rule $(fg)' = f'g + fg'$ in several variables. The last line says a pure rotation has no outflow, and a gradient field has no rotation.

Checking input and output types is a quick test for these formulas. In the third identity both sides are vectors. If a formula puts a vector where a scalar should be, it is wrong.

## 7. The Laplacian and the mean value property

In three variables,

$$
\Delta u = \frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2} + \frac{\partial^2 u}{\partial z^2}.
$$

The unknown $u$ means different things in different fields: a potential, a displacement, a groundwater head, a voltage. But the mathematics is the same. If you can solve one of these equations well, you can use it in many fields.

A good way to read the Laplacian is "the difference from the average of the neighbors". In one dimension, Taylor's theorem gives

$$
u''(x) = \frac{u(x + h) + u(x - h) - 2u(x)}{h^2} + O(h^2).
$$

So for small $h$:

- if $u''(x) = 0$, then $u(x) \approx \frac{1}{2}\big(u(x+h) + u(x-h)\big)$: the value equals the average of its neighbors;
- if $u''(x) > 0$, the value is below the average;
- if $u''(x) < 0$, the value is above the average.

In higher dimensions there is an exact version: if $\Delta u = 0$, then $u(\mathbf{x})$ equals the average of $u$ over any sphere centered at $\mathbf{x}$ (inside the domain). This is the **mean value property**. It is a theorem, not just a picture. And the formula above is already the finite difference we will meet in section 13.

## 8. The model equations

Four basic linear PDEs:

| equation | formula | $u$ means |
|---|---|---|
| Laplace | $\Delta u = 0$ | a potential at $(x, y, z)$ |
| Helmholtz | $\Delta u + k^2 u = 0$ | a pressure at $(x, y, z)$ |
| heat | $\partial_t u - \Delta u = 0$ | temperature at time $t$ |
| wave | $\partial_t^2 u - \Delta u = 0$ | displacement at time $t$ |

The Laplace equation is the most basic one. If you keep simplifying a hard problem, you often end up here. A nonzero right-hand side is called a **source** or a **load**, for example the weight of a structure.

**The heat equation smooths.** In 2D, convolving data with a Gaussian kernel,

$$
w(\mathbf{x}, t) = \int \frac{1}{4\pi t}\, e^{-|\mathbf{x} - \mathbf{y}|^2 / (4t)}\, u(\mathbf{y})\, d\mathbf{y},
$$

solves the heat equation with $w(\cdot, 0) = u$. As $t$ grows, the kernel gets wider and the data gets blurrier. This matches the mean value picture: at every moment, each value is pulled toward the average of its neighbors.

**Vector equations.** Some important equations have a vector unknown.

- *Maxwell's equations* for the electric and magnetic fields. In time-harmonic form, every $\partial_t$ becomes $-i\omega$, for example $\nabla \times \mathbf{E} = -\partial_t \mathbf{B}$ becomes $\nabla \times \mathbf{E} = i\omega \mathbf{B}$.
- *Linear elasticity*, for the displacement $\mathbf{u}$ of a solid, with Lamé constants $\lambda, \mu$ (density set to $1$):

$$
\partial_t^2 \mathbf{u} - \mu \Delta \mathbf{u} - (\lambda + \mu) \nabla (\nabla \cdot \mathbf{u}) = \mathbf{0}.
$$

- *Incompressible Navier–Stokes*, for a velocity $\mathbf{u}$ and pressure $p$ (viscosity set to $1$):

$$
\partial_t \mathbf{u} + (\mathbf{u} \cdot \nabla) \mathbf{u} - \Delta \mathbf{u} + \nabla p = \mathbf{0}, \qquad \nabla \cdot \mathbf{u} = 0.
$$

The condition $\nabla \cdot \mathbf{u} = 0$ means the fluid is incompressible: no net outflow at any point.

Roughly, difficulty grows like this: Laplace (no time, scalar), then heat and wave (time, scalar), then elasticity (vector unknown with coupled components), then Maxwell (several coupled vector fields). Helmholtz is a scalar equation, but it becomes hard when $k$ is large. The next section explains why.

## 9. From the wave equation to Helmholtz

Start from the wave equation with speed $c$:

$$
\frac{1}{c^2} \frac{\partial^2 u}{\partial t^2} - \Delta u = 0.
$$

Look for a solution that oscillates in time with one angular frequency $\omega$:

$$
u(\mathbf{x}, t) = v(\mathbf{x})\, e^{-i\omega t}.
$$

Two time derivatives give $(-i\omega)^2 = -\omega^2$. Substitute, and divide by $e^{-i\omega t}$, which is never zero:

$$
\Delta v + k^2 v = 0, \qquad k = \frac{\omega}{c}.
$$

This is the **Helmholtz equation**, and $k$ is the **wave number**. Time has disappeared. A problem in four variables becomes a problem in three, and all the time information is packed into one number $k$.

**Why large $k$ is hard.** The wavelength is $2\pi / k$. A large $k$ means the solution oscillates many times in a small region. To follow those oscillations, a numerical method needs many grid points per wavelength, so the cost grows very fast with $k$. High-frequency problems are still an open challenge.

## 10. Linear and nonlinear

Consider

$$
-\Delta u = f.
$$

If $f = 0$ it is the Laplace equation. If $f = f(\mathbf{x})$ it is the Poisson equation. If $f = f(\mathbf{x}, u)$ depends on the unknown, it is a nonlinear Poisson equation.

A PDE is **linear** if it is a linear combination of $u$ and its derivatives. The coefficients may depend on $x$ and $t$, but not on $u$. Some examples:

- $u_t + b(x, t)\, u_x = 0$ is linear;
- $u_t + (u_x)^2 = 0$ is nonlinear (a derivative is squared);
- $u_t = k(u)\, u_{xx}$ is nonlinear (the coefficient depends on $u$);
- $\nabla \cdot (k(u) \nabla u) = 0$ is nonlinear.

Why does this matter? For a linear equation, the sum of two solutions is again a solution (superposition). So we can solve simple pieces and add them up. For a nonlinear equation this fails. There may be many solutions, no solution, or a solution that blows up in finite time. Navier–Stokes is nonlinear because of the term $(\mathbf{u} \cdot \nabla)\mathbf{u}$.

## 11. Forward and inverse problems

- A **forward problem**: the equation, the coefficients, and the boundary or initial conditions are known. Find $u$.
- An **inverse problem**: measurements of $u$ are known (often only on the boundary). Find the coefficients or material properties inside.

Medical imaging is a good example. We apply currents through electrodes on the skin, measure voltages outside the body, and try to find the tissue properties inside. Finding the properties of rock underground from surface data is another one. In cooking terms: a forward problem predicts the taste from the recipe, and an inverse problem guesses the recipe from the taste.

Inverse problems are hard because the answer may not be unique, or because a tiny change in the data can change the answer a lot. That second issue is the topic of the next section.

## 12. Well-posedness

A problem is **well-posed** if three things hold:

1. **Existence.** There is at least one solution.
2. **Uniqueness.** There is only one solution.
3. **Stability.** The solution depends continuously on the data.

**Existence and uniqueness need the right conditions.** A second-order ODE is integrated twice, so it has two free constants and needs two conditions. These can be *initial conditions* (the value and the slope at one point) or *boundary conditions* (values at two points). But two conditions are not always enough. Here is my own example with $u'' + \pi^2 u = 0$, whose general solution is $u = A \cos \pi x + B \sin \pi x$:

- With $u(0) = 0,\ u(1) = 0$: we get $A = 0$, and $B$ can be anything. There are infinitely many solutions $u = B \sin \pi x$. Uniqueness fails.
- With $u(0) = 0,\ u(1) = 1$: we get $A = 0$, but then $u(1) = B \sin \pi = 0 \ne 1$. There is no solution. Existence fails.

Same equation, very similar data, and two different failures.

**Stability.** Take $u'' = u$ with $u(0) = 2,\ u'(0) = a$. The solution is

$$
u(x) = \frac{2 + a}{2} e^{x} + \frac{2 - a}{2} e^{-x}.
$$

For $a = -2$ the growing mode vanishes and $u = 2e^{-x}$. Now change the data a little, $a = -2 + \varepsilon$. The growing mode comes back:

$$
u_\varepsilon(x) - u(x) = \frac{\varepsilon}{2}\left(e^{x} - e^{-x}\right) = \varepsilon \sinh x.
$$

Stability means a bound of the form

$$
\| u_\varepsilon - u \| \le C\, \|\text{change in the data}\|
$$

with a finite constant $C$. On $[0, L]$ with the max norm, $C = \sinh L$. This is finite for every finite $L$, so the problem is stable there. On $[0, \infty)$ no finite $C$ works, so it is unstable.

Two points are easy to miss. First, "stable" does not mean "$C$ is small". Any finite $C$ is enough. Second, the same equation can be stable on one domain and unstable on another. The table below shows how fast $C$ grows. I perturbed $a$ by $\varepsilon = 10^{-6}$ and measured the largest change in $u$ on $[0, L]$ (1001 grid points):

| $L$ | largest change in $u$ | $C$ = change / $\varepsilon$ | $\sinh L$ |
|---|---|---|---|
| 1 | $1.18 \times 10^{-6}$ | $1.18$ | $1.18$ |
| 5 | $7.42 \times 10^{-5}$ | $74.2$ | $74.2$ |
| 10 | $1.10 \times 10^{-2}$ | $1.10 \times 10^{4}$ | $1.10 \times 10^{4}$ |
| 20 | $2.43 \times 10^{2}$ | $2.43 \times 10^{8}$ | $2.43 \times 10^{8}$ |
| 40 | $1.18 \times 10^{11}$ | $1.18 \times 10^{17}$ | $1.18 \times 10^{17}$ |

On $[0, 20]$, an error in the sixth digit of $a$ already changes $u$ by more than $200$. Real data always has measurement error, so this matters in practice.

## 13. Green's formula

One formula is enough to remember:

$$
\int_\Omega \mathbf{F} \cdot \nabla u \, d\mathbf{x} + \int_\Omega (\nabla \cdot \mathbf{F})\, u \, d\mathbf{x} = \int_{\partial\Omega} (\mathbf{F} \cdot \mathbf{n})\, u \, dS.
$$

It is integration by parts in several variables. The proof is the product rule from section 6 plus the divergence theorem. The idea: we can move a derivative from $u$ to $\mathbf{F}$, and the price is a boundary term.

Everything else follows from it.

- **$u = 1$.** Then $\nabla u = 0$, and we get the divergence theorem back.
- **$\mathbf{F} = \nabla v$.** Then $\nabla \cdot \mathbf{F} = \Delta v$, and we get **Green's first identity**:

$$
\int_\Omega \nabla v \cdot \nabla u \, d\mathbf{x} + \int_\Omega (\Delta v)\, u \, d\mathbf{x} = \int_{\partial\Omega} \frac{\partial v}{\partial n}\, u \, dS.
$$

- **Swap $u$ and $v$ and subtract.** The term $\nabla v \cdot \nabla u$ is symmetric, so it cancels, and we get **Green's second identity**:

$$
\int_\Omega (u\, \Delta v - v\, \Delta u) \, d\mathbf{x} = \int_{\partial\Omega} \left( u \frac{\partial v}{\partial n} - v \frac{\partial u}{\partial n} \right) dS.
$$

Here $\partial v / \partial n = \nabla v \cdot \mathbf{n}$ is the directional derivative in the normal direction. Some books write $\nu$ instead of $n$.

Green's first identity is the starting point of the finite element method. It turns the second derivative $\Delta u$ into a product of two first derivatives, $\nabla u \cdot \nabla v$. That is the core of the weak form.

## 14. Finite differences

A derivative is a limit of difference quotients. A finite difference simply drops the limit and keeps a small $h$:

$$
u'(x) \approx \frac{u(x + h) - u(x)}{h}.
$$

On a grid $x_j$ with spacing $h$ and values $u_j \approx u(x_j)$, the three basic formulas are

$$
\delta^+ u_j = \frac{u_{j+1} - u_j}{h}, \qquad \delta^- u_j = \frac{u_j - u_{j-1}}{h}, \qquad \delta u_j = \frac{u_{j+1} - u_{j-1}}{2h}.
$$

These are the forward, backward, and central differences. The central one divides by $2h$ because the two points are $2h$ apart.

**Taylor's theorem gives the error.**

$$
u(x \pm h) = u(x) \pm h\, u'(x) + \frac{h^2}{2} u''(x) \pm \frac{h^3}{6} u'''(x) + \cdots
$$

Changing $h$ to $-h$ flips the sign of the odd terms. For the forward difference, the first term we drop is $\frac{h}{2} u''(x)$, so the error is $O(h)$. For the central difference, *subtract* the two expansions. The value $u(x)$ and all even terms cancel:

$$
\frac{u(x + h) - u(x - h)}{2h} = u'(x) + \frac{h^2}{6} u'''(x) + \cdots,
$$

so the error is $O(h^2)$. Here $g(h) = O(h^p)$ means $|g(h)| \le C h^p$ for some constant $C$.

To get the second derivative, *add* the two expansions instead. Now the odd terms cancel:

$$
\frac{u(x - h) - 2u(x) + u(x + h)}{h^2} = u''(x) + \frac{h^2}{12} u''''(x) + \cdots.
$$

This is the formula from section 7. Setting it to zero says "the value is the average of the two neighbors", so the discrete Laplace equation inherits the mean value property.

More points give higher order. For example,

$$
u'(x) = \frac{u(x - 2h) - 8u(x - h) + 8u(x + h) - u(x + 2h)}{12h} + O(h^4).
$$

If we halve $h$, an $O(h)$ error halves, an $O(h^2)$ error drops to $1/4$, and an $O(h^4)$ error drops to $1/16$. The central difference uses information from both sides, and it buys one extra order for the same two function values.

But there is no free lunch. Every approximation has a cost, and we choose $h$ to balance accuracy against work. The experiment below shows one more cost: on a computer, $h$ cannot be made as small as we like.

## 15. Experiment: four difference formulas

I tested four formulas on $u(x) = e^{\sin x}$ at $x = 1$, where the exact derivatives are known:

$$
u'(x) = \cos x\; e^{\sin x}, \qquad u''(x) = (\cos^2 x - \sin x)\, e^{\sin x}.
$$

The step sizes were $h = 2^{-1}, 2^{-2}, \ldots, 2^{-30}$. "Second difference" approximates $u''$. The other three approximate $u'$.

```python
def forward(u, x, h):
    return (u(x + h) - u(x)) / h

def central(u, x, h):
    return (u(x + h) - u(x - h)) / (2 * h)

def fourth_order(u, x, h):
    return (u(x - 2*h) - 8*u(x - h) + 8*u(x + h) - u(x + 2*h)) / (12 * h)

def second_central(u, x, h):
    return (u(x - h) - 2*u(x) + u(x + h)) / h**2
```

<figure class="viz">
<div class="viz-legend"><span class="viz-key"><span class="viz-swatch" style="background:#2a78d6"></span>Forward</span><span class="viz-key"><span class="viz-swatch" style="background:#eb6834"></span>Central</span><span class="viz-key"><span class="viz-swatch" style="background:#1baf7a"></span>Fourth order</span><span class="viz-key"><span class="viz-swatch" style="background:#eda100"></span>Second difference</span></div>
<div class="viz-frame">
<svg viewBox="0 0 680 400" role="img" tabindex="0" aria-label="Error of four finite difference formulas as h gets smaller, on log scales. Each error first falls with slope 1, 2, 4, or 2, then rises again because of round-off." class="viz-plot"><line x1="64" x2="540" y1="30.0" y2="30.0" class="viz-grid"/><text x="54" y="34.0" text-anchor="end" class="viz-tick">10<tspan dy="-0.5em" font-size="75%">2</tspan></text><line x1="64" x2="540" y1="107.5" y2="107.5" class="viz-grid"/><text x="54" y="111.5" text-anchor="end" class="viz-tick">10<tspan dy="-0.5em" font-size="75%">−2</tspan></text><line x1="64" x2="540" y1="185.0" y2="185.0" class="viz-grid"/><text x="54" y="189.0" text-anchor="end" class="viz-tick">10<tspan dy="-0.5em" font-size="75%">−6</tspan></text><line x1="64" x2="540" y1="262.5" y2="262.5" class="viz-grid"/><text x="54" y="266.5" text-anchor="end" class="viz-tick">10<tspan dy="-0.5em" font-size="75%">−10</tspan></text><line x1="64" x2="540" y1="340.0" y2="340.0" class="viz-grid"/><text x="54" y="344.0" text-anchor="end" class="viz-tick">10<tspan dy="-0.5em" font-size="75%">−14</tspan></text><text x="64.0" y="362" text-anchor="middle" class="viz-tick">2<tspan dy="-0.5em" font-size="75%">−1</tspan></text><text x="129.7" y="362" text-anchor="middle" class="viz-tick">2<tspan dy="-0.5em" font-size="75%">−5</tspan></text><text x="211.7" y="362" text-anchor="middle" class="viz-tick">2<tspan dy="-0.5em" font-size="75%">−10</tspan></text><text x="293.8" y="362" text-anchor="middle" class="viz-tick">2<tspan dy="-0.5em" font-size="75%">−15</tspan></text><text x="375.9" y="362" text-anchor="middle" class="viz-tick">2<tspan dy="-0.5em" font-size="75%">−20</tspan></text><text x="457.9" y="362" text-anchor="middle" class="viz-tick">2<tspan dy="-0.5em" font-size="75%">−25</tspan></text><text x="540.0" y="362" text-anchor="middle" class="viz-tick">2<tspan dy="-0.5em" font-size="75%">−30</tspan></text><text x="302.0" y="388" text-anchor="middle" class="viz-axis-title">Step size h (smaller to the right)</text><text x="64" y="14" class="viz-axis-title">Error</text><path d="M64.0,75.1 L80.4,82.3 L96.8,89.0 L113.2,95.3 L129.7,101.4 L146.1,107.4 L162.5,113.3 L178.9,119.2 L195.3,125.0 L211.7,130.9 L228.1,136.7 L244.6,142.5 L261.0,148.4 L277.4,154.2 L293.8,160.0 L310.2,165.9 L326.6,171.7 L343.0,177.5 L359.4,183.4 L375.9,189.2 L392.3,195.0 L408.7,200.8 L425.1,206.9 L441.5,212.8 L457.9,226.8 L474.3,225.7 L490.8,225.7 L507.2,225.7 L523.6,225.7 L540.0,225.7" fill="none" stroke="#2a78d6" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><circle cx="457.9" cy="226.8" r="4" fill="#2a78d6" stroke="#ffffff" stroke-width="2"/><path d="M64.0,84.3 L80.4,95.5 L96.8,107.1 L113.2,118.7 L129.7,130.4 L146.1,142.0 L162.5,153.7 L178.9,165.4 L195.3,177.0 L211.7,188.7 L228.1,200.4 L244.6,212.0 L261.0,223.7 L277.4,235.4 L293.8,247.0 L310.2,258.5 L326.6,277.8 L343.0,269.8 L359.4,269.8 L375.9,269.8 L392.3,250.4 L408.7,250.4 L425.1,248.8 L441.5,248.8 L457.9,248.8 L474.3,225.7 L490.8,217.2 L507.2,225.7 L523.6,203.5 L540.0,197.3" fill="none" stroke="#eb6834" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><circle cx="326.6" cy="277.8" r="4" fill="#eb6834" stroke="#ffffff" stroke-width="2"/><path d="M64.0,96.1 L80.4,117.9 L96.8,140.8 L113.2,164.1 L129.7,187.4 L146.1,210.7 L162.5,234.0 L178.9,257.4 L195.3,280.7 L211.7,305.8 L228.1,314.4 L244.6,307.6 L261.0,318.1 L277.4,296.9 L293.8,307.6 L310.2,296.9 L326.6,270.8 L343.0,268.1 L359.4,269.8 L375.9,264.3 L392.3,245.7 L408.7,269.8 L425.1,238.4 L441.5,248.8 L457.9,233.9 L474.3,228.8 L490.8,211.7 L507.2,222.3 L523.6,198.9 L540.0,197.3" fill="none" stroke="#1baf7a" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><circle cx="261.0" cy="318.1" r="4" fill="#1baf7a" stroke="#ffffff" stroke-width="2"/><path d="M64.0,100.4 L80.4,113.0 L96.8,125.0 L113.2,136.7 L129.7,148.4 L146.1,160.1 L162.5,171.7 L178.9,183.4 L195.3,195.1 L211.7,206.7 L228.1,217.0 L244.6,221.6 L261.0,221.6 L277.4,200.3 L293.8,211.5 L310.2,179.4 L326.6,166.0 L343.0,153.9 L359.4,142.2 L375.9,134.5 L392.3,124.1 L408.7,124.1 L425.1,99.9 L441.5,99.9 L457.9,70.9 L474.3,66.7 L490.8,50.0 L507.2,66.7 L523.6,30.0 L540.0,30.0" fill="none" stroke="#eda100" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><circle cx="244.6" cy="221.6" r="4" fill="#eda100" stroke="#ffffff" stroke-width="2"/><text x="550" y="229.7" class="viz-label">Forward</text><text x="550" y="193.3" class="viz-label">Central</text><text x="550" y="211.3" class="viz-label">Fourth order</text><text x="550" y="34.0" class="viz-label">Second difference</text><line x1="0" x2="0" y1="30" y2="340" class="viz-crosshair" visibility="hidden"/></svg>
<div class="viz-tooltip" hidden></div>
</div>
<figcaption>Error of each formula for u(x) = e<sup>sin x</sup> at x = 1, on log scales. The dot marks the smallest error. Values above 10<sup>2</sup> are drawn at the top edge. Hover over the chart, or focus it and use the arrow keys, to read the values.</figcaption>
<details class="viz-table"><summary>Show the numbers</summary><table><thead><tr><th>h</th><th>Forward</th><th>Central</th><th>Fourth order</th><th>Second difference</th></tr></thead><tbody><tr><td>2<sup>−1</sup></td><td>4.7e-1</td><td>1.6e-1</td><td>3.9e-2</td><td>2.3e-2</td></tr><tr><td>2<sup>−2</sup></td><td>2.0e-1</td><td>4.1e-2</td><td>2.9e-3</td><td>5.2e-3</td></tr><tr><td>2<sup>−3</sup></td><td>9.0e-2</td><td>1.1e-2</td><td>1.9e-4</td><td>1.3e-3</td></tr><tr><td>2<sup>−4</sup></td><td>4.2e-2</td><td>2.6e-3</td><td>1.2e-5</td><td>3.1e-4</td></tr><tr><td>2<sup>−5</sup></td><td>2.1e-2</td><td>6.6e-4</td><td>7.5e-7</td><td>7.7e-5</td></tr><tr><td>2<sup>−6</sup></td><td>1.0e-2</td><td>1.6e-4</td><td>4.7e-8</td><td>1.9e-5</td></tr><tr><td>2<sup>−7</sup></td><td>5.0e-3</td><td>4.1e-5</td><td>2.9e-9</td><td>4.8e-6</td></tr><tr><td>2<sup>−8</sup></td><td>2.5e-3</td><td>1.0e-5</td><td>1.8e-10</td><td>1.2e-6</td></tr><tr><td>2<sup>−9</sup></td><td>1.2e-3</td><td>2.6e-6</td><td>1.1e-11</td><td>3.0e-7</td></tr><tr><td>2<sup>−10</sup></td><td>6.2e-4</td><td>6.4e-7</td><td>5.9e-13</td><td>7.6e-8</td></tr><tr><td>2<sup>−11</sup></td><td>3.1e-4</td><td>1.6e-7</td><td>2.1e-13</td><td>2.2e-8</td></tr><tr><td>2<sup>−12</sup></td><td>1.6e-4</td><td>4.0e-8</td><td>4.7e-13</td><td>1.3e-8</td></tr><tr><td>2<sup>−13</sup></td><td>7.8e-5</td><td>1.0e-8</td><td>1.3e-13</td><td>1.3e-8</td></tr><tr><td>2<sup>−14</sup></td><td>3.9e-5</td><td>2.5e-9</td><td>1.7e-12</td><td>1.6e-7</td></tr><tr><td>2<sup>−15</sup></td><td>1.9e-5</td><td>6.3e-10</td><td>4.7e-13</td><td>4.3e-8</td></tr><tr><td>2<sup>−16</sup></td><td>9.7e-6</td><td>1.6e-10</td><td>1.7e-12</td><td>2.0e-6</td></tr><tr><td>2<sup>−17</sup></td><td>4.9e-6</td><td>1.6e-11</td><td>3.7e-11</td><td>9.6e-6</td></tr><tr><td>2<sup>−18</sup></td><td>2.4e-6</td><td>4.2e-11</td><td>5.2e-11</td><td>4.0e-5</td></tr><tr><td>2<sup>−19</sup></td><td>1.2e-6</td><td>4.2e-11</td><td>4.2e-11</td><td>1.6e-4</td></tr><tr><td>2<sup>−20</sup></td><td>6.1e-7</td><td>4.2e-11</td><td>8.1e-11</td><td>4.1e-4</td></tr><tr><td>2<sup>−21</sup></td><td>3.0e-7</td><td>4.2e-10</td><td>7.3e-10</td><td>1.4e-3</td></tr><tr><td>2<sup>−22</sup></td><td>1.5e-7</td><td>4.2e-10</td><td>4.2e-11</td><td>1.4e-3</td></tr><tr><td>2<sup>−23</sup></td><td>7.4e-8</td><td>5.1e-10</td><td>1.7e-9</td><td>2.5e-2</td></tr><tr><td>2<sup>−24</sup></td><td>3.7e-8</td><td>5.1e-10</td><td>5.1e-10</td><td>2.5e-2</td></tr><tr><td>2<sup>−25</sup></td><td>6.9e-9</td><td>5.1e-10</td><td>3.0e-9</td><td>7.7e-1</td></tr><tr><td>2<sup>−26</sup></td><td>8.0e-9</td><td>8.0e-9</td><td>5.5e-9</td><td>1.3e0</td></tr><tr><td>2<sup>−27</sup></td><td>8.0e-9</td><td>2.2e-8</td><td>4.2e-8</td><td>9.3e0</td></tr><tr><td>2<sup>−28</sup></td><td>8.0e-9</td><td>8.0e-9</td><td>1.2e-8</td><td>1.3e0</td></tr><tr><td>2<sup>−29</sup></td><td>8.0e-9</td><td>1.1e-7</td><td>1.9e-7</td><td>1.3e2</td></tr><tr><td>2<sup>−30</sup></td><td>8.0e-9</td><td>2.3e-7</td><td>2.3e-7</td><td>5.1e2</td></tr></tbody></table></details>
<script>(() => {const fig = document.currentScript.closest('figure');const svg = fig.querySelector('svg');const line = svg.querySelector('.viz-crosshair');const tip = fig.querySelector('.viz-tooltip');const names = ["Forward", "Central", "Fourth order", "Second difference"]; const colors = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"];const rows = [["4.7e-1", "1.6e-1", "3.9e-2", "2.3e-2"], ["2.0e-1", "4.1e-2", "2.9e-3", "5.2e-3"], ["9.0e-2", "1.1e-2", "1.9e-4", "1.3e-3"], ["4.2e-2", "2.6e-3", "1.2e-5", "3.1e-4"], ["2.1e-2", "6.6e-4", "7.5e-7", "7.7e-5"], ["1.0e-2", "1.6e-4", "4.7e-8", "1.9e-5"], ["5.0e-3", "4.1e-5", "2.9e-9", "4.8e-6"], ["2.5e-3", "1.0e-5", "1.8e-10", "1.2e-6"], ["1.2e-3", "2.6e-6", "1.1e-11", "3.0e-7"], ["6.2e-4", "6.4e-7", "5.9e-13", "7.6e-8"], ["3.1e-4", "1.6e-7", "2.1e-13", "2.2e-8"], ["1.6e-4", "4.0e-8", "4.7e-13", "1.3e-8"], ["7.8e-5", "1.0e-8", "1.3e-13", "1.3e-8"], ["3.9e-5", "2.5e-9", "1.7e-12", "1.6e-7"], ["1.9e-5", "6.3e-10", "4.7e-13", "4.3e-8"], ["9.7e-6", "1.6e-10", "1.7e-12", "2.0e-6"], ["4.9e-6", "1.6e-11", "3.7e-11", "9.6e-6"], ["2.4e-6", "4.2e-11", "5.2e-11", "4.0e-5"], ["1.2e-6", "4.2e-11", "4.2e-11", "1.6e-4"], ["6.1e-7", "4.2e-11", "8.1e-11", "4.1e-4"], ["3.0e-7", "4.2e-10", "7.3e-10", "1.4e-3"], ["1.5e-7", "4.2e-10", "4.2e-11", "1.4e-3"], ["7.4e-8", "5.1e-10", "1.7e-9", "2.5e-2"], ["3.7e-8", "5.1e-10", "5.1e-10", "2.5e-2"], ["6.9e-9", "5.1e-10", "3.0e-9", "7.7e-1"], ["8.0e-9", "8.0e-9", "5.5e-9", "1.3e0"], ["8.0e-9", "2.2e-8", "4.2e-8", "9.3e0"], ["8.0e-9", "8.0e-9", "1.2e-8", "1.3e0"], ["8.0e-9", "1.1e-7", "1.9e-7", "1.3e2"], ["8.0e-9", "2.3e-7", "2.3e-7", "5.1e2"]]; const heads = ["h = 2^-1", "h = 2^-2", "h = 2^-3", "h = 2^-4", "h = 2^-5", "h = 2^-6", "h = 2^-7", "h = 2^-8", "h = 2^-9", "h = 2^-10", "h = 2^-11", "h = 2^-12", "h = 2^-13", "h = 2^-14", "h = 2^-15", "h = 2^-16", "h = 2^-17", "h = 2^-18", "h = 2^-19", "h = 2^-20", "h = 2^-21", "h = 2^-22", "h = 2^-23", "h = 2^-24", "h = 2^-25", "h = 2^-26", "h = 2^-27", "h = 2^-28", "h = 2^-29", "h = 2^-30"];const left = 64, right = 540, steps = 29, width = 680;let current = -1;const show = (n) => { n = Math.max(0, Math.min(steps, n)); current = n; const x = left + (right - left) * n / steps; line.setAttribute('x1', x); line.setAttribute('x2', x); line.setAttribute('visibility', 'visible'); tip.replaceChildren(); const head = document.createElement('div'); head.className = 'viz-tooltip-head'; head.textContent = heads[n]; tip.append(head); names.forEach((name, i) => {  const row = document.createElement('div'); row.className = 'viz-tooltip-row';  const key = document.createElement('span'); key.className = 'viz-swatch'; key.style.background = colors[i];  const value = document.createElement('strong'); value.textContent = rows[n][i];  const label = document.createElement('span'); label.textContent = name;  row.append(key, value, label); tip.append(row); }); const box = svg.getBoundingClientRect(), frame = fig.getBoundingClientRect(); const xs = box.left - frame.left + x * box.width / width; tip.hidden = false; const room = frame.width - tip.offsetWidth - 8; tip.style.left = Math.max(0, Math.min(room, xs + 12)) + 'px'; tip.style.top = (box.top - frame.top + 8) + 'px';};const hide = () => { line.setAttribute('visibility', 'hidden'); tip.hidden = true; current = -1; };svg.addEventListener('pointermove', (event) => { const box = svg.getBoundingClientRect(); const x = (event.clientX - box.left) * width / box.width; show(Math.round((x - left) / (right - left) * steps));});svg.addEventListener('pointerleave', hide);svg.addEventListener('focus', () => show(current < 0 ? 0 : current));svg.addEventListener('blur', hide);svg.addEventListener('keydown', (event) => { if (event.key === 'ArrowRight') { show(current + 1); event.preventDefault(); } if (event.key === 'ArrowLeft') { show(current - 1); event.preventDefault(); } if (event.key === 'Escape') hide();});})();</script>
</figure>

Each line has two parts.

**Left part: truncation error.** As $h$ shrinks, the error falls along a straight line on the log–log chart. Its slope is the order. I measured the order as $\log_2(e_{h}/e_{h/2})$ for $h = 2^{-4}$ to $2^{-9}$:

| formula | approximates | predicted order | measured order |
|---|---|---|---|
| Forward | $u'$ | 1 | 1.00 to 1.05 |
| Central | $u'$ | 2 | 2.00 |
| Fourth order | $u'$ | 4 | 4.00 |
| Second difference | $u''$ | 2 | 2.00 |

**Right part: round-off error.** A computer stores about 16 significant digits, so each value of $u$ has a tiny error of about $\epsilon \approx 2.2 \times 10^{-16}$. The formulas subtract values that are almost equal and then divide by $h$ (or $h^2$). So the round-off part of the error grows like $\epsilon / h$ (or $\epsilon / h^2$) as $h$ shrinks. When it becomes larger than the truncation error, the total error starts to rise.

The best $h$ is where the two parts are about equal. For an order-$p$ formula with round-off $\epsilon / h^q$, setting $h^p \approx \epsilon / h^q$ gives $h \approx \epsilon^{1/(p+q)}$ and a smallest error of about $\epsilon^{p/(p+q)}$. These are rough estimates, because they ignore the constants in front:

| formula | best $h$, measured | best $h$, estimate | smallest error, measured | smallest error, estimate |
|---|---|---|---|---|
| Forward | $2^{-25} \approx 3.0 \times 10^{-8}$ | $1.5 \times 10^{-8}$ | $6.9 \times 10^{-9}$ | $1 \times 10^{-8}$ |
| Central | $2^{-17} \approx 7.6 \times 10^{-6}$ | $6.1 \times 10^{-6}$ | $1.6 \times 10^{-11}$ | $4 \times 10^{-11}$ |
| Fourth order | $2^{-13} \approx 1.2 \times 10^{-4}$ | $7.4 \times 10^{-4}$ | $1.3 \times 10^{-13}$ | $3 \times 10^{-13}$ |
| Second difference | $2^{-12} \approx 2.4 \times 10^{-4}$ | $1.2 \times 10^{-4}$ | $1.3 \times 10^{-8}$ | $1 \times 10^{-8}$ |

The estimates are close, within a factor of about $6$ for $h$ and $3$ for the error. For the fourth-order formula the error stays near $10^{-13}$ from $h = 2^{-11}$ to $2^{-13}$, so the exact position of the minimum is decided by round-off noise.

What I noticed:

- **Higher order wins, but only to a point.** The fourth-order formula reaches an error of $10^{-13}$ with a fairly large $h$. The forward difference never gets below about $10^{-8}$, however small we make $h$.
- **The second derivative is the most fragile.** Its round-off grows like $\epsilon / h^2$. At $h = 2^{-30}$ its error is about $500$, which is useless. This matters for PDEs, because the discrete Laplacian has exactly this $1/h^2$ form.
- **Smaller $h$ is not always better.** In practice we rarely go near these limits, because PDE grids are much coarser. But it is good to know where the wall is.

The full script, with the stability table from section 12 and the chart, is in [this site's repository](https://github.com/doubleyands/doubleyands.github.io/tree/main/scripts/posts/numerical-pde/lecture-1).

## Summary

1. The **gradient** points in the direction of steepest increase and is normal to level sets. The unit outward normal is $\nabla g / |\nabla g|$. Do not forget to divide by the length.
2. The **divergence** has a computing definition (sum of partial derivatives) and a physical one (outflow per unit volume). The **divergence theorem** turns a volume integral into a boundary integral.
3. $\Delta u = \nabla \cdot \nabla u$. **$\Delta u = 0$ means each value equals the average of its neighbors.** The discrete Laplacian keeps this property.
4. Putting $u = v(\mathbf{x}) e^{-i\omega t}$ into the wave equation gives the **Helmholtz equation** $\Delta v + k^2 v = 0$. Large $k$ is expensive.
5. **Well-posed = existence + uniqueness + stability.** Stability asks whether $C$ is finite, not whether it is small, and the answer can depend on the domain.
6. **One integration-by-parts formula** gives the divergence theorem and both Green's identities. Green's first identity is the root of the finite element method.
7. **Central differences are second order**, because the even terms cancel. Round-off sets a floor, and the second difference hits it first.
