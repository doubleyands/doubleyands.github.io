---
title: "Numerical Analysis, Week 7: Linear Algebra Review and Gaussian Elimination"
description: "Vector spaces, symmetric matrices, similar matrices, the Schur form, the principal axes theorem, the SVD and the Jordan form; vector and matrix norms, the spectral radius, powers of a matrix and the Neumann series; Gaussian elimination as LU, its O(n³) cost, pivoting, Cholesky, and error vs. residual. Checked with five small experiments."
pubDate: 2026-11-14
project: numerical-analysis
tags: [numerical analysis, linear algebra, matrix norms, gaussian elimination, LU decomposition, notes]
---

This is my summary of the seventh week of a numerical analysis course. As before, I write the ideas in my own words and check some of them with a small experiment. The course follows K. E. Atkinson, *An Introduction to Numerical Analysis*.

Until now the course was about one equation at a time: roots, interpolation, integrals and ODEs (see [Week 6](/blog/2026-11-07-numerical-analysis-week-6/)). This week starts a new part. The next chapters are all about one problem:

$$
A\mathbf{x} = \mathbf{b}.
$$

Why is this problem so important? Three reasons:

- When we discretize a PDE, we almost always get a linear system.
- Nonlinear problems lead to linear systems too. Newton's method for a system ([Week 2](/blog/2026-10-10-numerical-analysis-week-2/)) solves one linear system in every step.
- Real problems have millions or billions of unknowns. Then an exact solve is too expensive, and we need good approximations. To build them, we must know the matrix well.

So the week has two parts. First a review of linear algebra, with the facts we will use again and again. Then the first method for $A\mathbf{x} = \mathbf{b}$: Gaussian elimination.

## 1. Vector spaces and bases

A **vector space** $V$ is a set with two operations, vector addition and multiplication by a scalar. The set must be **closed** under both operations:

$$
\mathbf{u}, \mathbf{v} \in V \ \Rightarrow\ \mathbf{u} + \mathbf{v} \in V, \qquad c\,\mathbf{u} \in V.
$$

If a result could leave $V$, we could not continue to compute inside $V$. The scalars $c$ are real or complex numbers. Examples: $\mathbb{R}^n$, $\mathbb{C}^n$, the polynomials of degree at most $n$, and $C[a,b]$, the continuous functions on $[a,b]$.

Vectors $\mathbf{v}_1, \dots, \mathbf{v}_m$ are **linearly independent** if

$$
\alpha_1 \mathbf{v}_1 + \cdots + \alpha_m \mathbf{v}_m = \mathbf{0} \ \Rightarrow\ \alpha_1 = \cdots = \alpha_m = 0.
$$

A **basis** is a set of vectors with two properties:

1. It **spans** $V$: every $\mathbf{u} \in V$ is a linear combination of the basis vectors.
2. It is linearly independent. Because of this, the combination in 1 is unique.

A basis is not unique. For example, $\{(1,0), (0,1)\}$ and $\{(1,2), (2,-1)\}$ are both bases of $\mathbb{R}^2$. For polynomials we can use $\{1, x, x^2, \dots\}$ or the orthogonal polynomials from [Week 3](/blog/2026-10-17-numerical-analysis-week-3/).

The **dimension** is the number of vectors in a basis. The polynomials of degree at most $n$ have dimension $n + 1$. $\mathbb{R}^n$ and $\mathbb{C}^n$ use the same basis $\{\mathbf{e}_1, \dots, \mathbf{e}_n\}$; only the coefficients differ (real or complex). $C[a,b]$ has infinite dimension.

## 2. Matrices: symmetric and skew-symmetric

Matrix multiplication is not commutative: in general $AB \ne BA$. We write $A^t$ for the transpose and $A^* = \overline{A}^{\,t}$ for the conjugate transpose.

Two special kinds of matrix are very important:

- **Symmetric**, $A^t = A$. A real symmetric matrix has only **real** eigenvalues.
- **Skew-symmetric**, $A^t = -A$. A real skew-symmetric matrix has only **purely imaginary** eigenvalues (or zero).

The proof is short. Let $A\mathbf{x} = \lambda \mathbf{x}$ with $\mathbf{x} \ne \mathbf{0}$, and use the complex inner product from Section 4. If $A$ is real symmetric, then $A^* = A$, so

$$
\lambda \|\mathbf{x}\|^2 = (A\mathbf{x}, \mathbf{x}) = (\mathbf{x}, A^*\mathbf{x}) = (\mathbf{x}, A\mathbf{x}) = \bar\lambda \|\mathbf{x}\|^2 .
$$

So $\lambda = \bar\lambda$, and $\lambda$ is real. If $A$ is real skew-symmetric, then $A^* = -A$, and the same lines give $\lambda = -\bar\lambda$, so $\mathrm{Re}\,\lambda = 0$. For example, $\begin{bmatrix} 0 & 3 \\ -3 & 0 \end{bmatrix}$ has eigenvalues $\pm 3i$.

Why do we care? Real eigenvalues can be put in order, like $0.5 < 2 < 7$. Later we will see that the speed of many solvers for $A\mathbf{x} = \mathbf{b}$ depends on where the eigenvalues are, especially the smallest one.

The **rank** of a matrix is the largest number of linearly independent rows. It is also the largest number of linearly independent columns.

## 3. Five ways to say "invertible"

Let $A$ be $n \times n$. These five statements are equivalent:

1. For every $\mathbf{b}$, $A\mathbf{x} = \mathbf{b}$ has exactly one solution.
2. $A\mathbf{x} = \mathbf{0}$ only for $\mathbf{x} = \mathbf{0}$.
3. $A^{-1}$ exists.
4. $\det A \ne 0$.
5. $\mathrm{rank}(A) = n$.

Some of the links are easy to see.

**1 ⇒ 3.** Solve $A\mathbf{x}_k = \mathbf{e}_k$ for each unit vector $\mathbf{e}_k$. Put the solutions side by side as columns. Then $A[\mathbf{x}_1 \cdots \mathbf{x}_n] = [\mathbf{e}_1 \cdots \mathbf{e}_n] = I$, so this matrix is $A^{-1}$.

**4.** The determinant is the product of the eigenvalues, $\det A = \lambda_1 \lambda_2 \cdots \lambda_n$ (this follows from the Schur form in Section 6). So $\det A \ne 0$ means no eigenvalue is zero. If some eigenvalue were zero, there would be an $\mathbf{x} \ne \mathbf{0}$ with $A\mathbf{x} = \mathbf{0}$, and then every multiple $c\,\mathbf{x}$ also solves $A\mathbf{x} = \mathbf{0}$. For example, $\begin{bmatrix} 1 & 1 \\ 1 & 1 \end{bmatrix}$ sends $(1,-1)$, $(50,-50)$ and $(-10^4, 10^4)$ all to zero.

**5 ⇒ 2.** The product $A\mathbf{x}$ is a combination of the columns of $A$:

$$
A\mathbf{x} = x_1 \mathbf{a}_1 + x_2 \mathbf{a}_2 + \cdots + x_n \mathbf{a}_n ,
$$

where $\mathbf{a}_j$ is column $j$. If the $n$ columns are independent, this is zero only when all $x_j = 0$.

## 4. Inner products, norms and orthogonality

In $\mathbb{C}^n$ the inner product has a complex conjugate:

$$
(\mathbf{x}, \mathbf{y}) = x_1 \overline{y_1} + \cdots + x_n \overline{y_n}.
$$

We need the conjugate so that $(\mathbf{x}, \mathbf{x}) \ge 0$. For one entry $x = 1 + 2i$: with the conjugate, $(1 + 2i)(1 - 2i) = 5 \ge 0$. Without it, $(1 + 2i)^2 = -3 + 4i$, which is not even real.

From the inner product we get:

- the norm $\|\mathbf{x}\|_2 = \sqrt{(\mathbf{x}, \mathbf{x})} = \sqrt{|x_1|^2 + \cdots + |x_n|^2}$,
- the Cauchy–Schwarz inequality $|(\mathbf{x}, \mathbf{y})| \le \|\mathbf{x}\|\,\|\mathbf{y}\|$,
- the triangle inequality $\|\mathbf{x} + \mathbf{y}\| \le \|\mathbf{x}\| + \|\mathbf{y}\|$,
- the rule $(A\mathbf{x}, \mathbf{y}) = (\mathbf{x}, A^*\mathbf{y})$.

Two vectors are **orthogonal** if $(\mathbf{x}, \mathbf{y}) = 0$. If they also have norm 1, they are **orthonormal**. A matrix $U$ is **unitary** if $U^*U = UU^* = I$, so $U^{-1} = U^*$. A real unitary matrix is called an **orthogonal matrix**. The columns of a unitary matrix form an orthonormal basis.

Orthogonal bases make coefficients easy. If $\mathbf{u}_1, \dots, \mathbf{u}_n$ is an orthogonal basis and $\mathbf{u} = \alpha_1 \mathbf{u}_1 + \cdots + \alpha_n \mathbf{u}_n$, take the inner product with $\mathbf{u}_j$. All the other terms vanish, and

$$
\alpha_j = \frac{(\mathbf{u}, \mathbf{u}_j)}{\|\mathbf{u}_j\|^2}.
$$

For an orthonormal basis, simply $\alpha_j = (\mathbf{u}, \mathbf{u}_j)$. No linear system is needed.

## 5. Eigenvalues and similar matrices

If $A\mathbf{x} = \lambda \mathbf{x}$ with $\mathbf{x} \ne \mathbf{0}$, then $\lambda$ is an **eigenvalue** and $\mathbf{x}$ is an **eigenvector**. An eigenvalue can be zero, but an eigenvector cannot. Eigenvalues are defined only for square matrices. (For rectangular matrices we use singular values, below.)

The eigenvalues are the roots of the **characteristic polynomial** $f_A(\lambda) = \det(A - \lambda I)$. An $n \times n$ matrix has between $1$ and $n$ distinct eigenvalues (if we allow complex ones), so we must always think about **multiplicity**. For example, $\begin{bmatrix} 3 & 1 \\ 1 & 3 \end{bmatrix}$ has $f_A(\lambda) = (3 - \lambda)^2 - 1 = (\lambda - 2)(\lambda - 4)$, eigenvalues $2$ and $4$, and eigenvectors $(1, -1)$ and $(1, 1)$. The matrix is symmetric, so the eigenvalues are real and the eigenvectors are orthogonal.

### Similar matrices

Two matrices $A$ and $B$ are **similar** if $B = P^{-1}AP$ for some nonsingular $P$.

This idea is the starting point of **preconditioning**. If we multiply both sides of $A\mathbf{x} = \mathbf{b}$ by a nonsingular matrix, the solution does not change. In one dimension: $2x = 6$ and $x = 3$ have the same solution. So, when $A\mathbf{x} = \mathbf{b}$ is hard to solve, we can look for a related system that is easier. For example, solve $(P^{-1}AP)\,\tilde{\mathbf{x}} = P^{-1}\mathbf{b}$. Then $A(P\tilde{\mathbf{x}}) = \mathbf{b}$, so $\mathbf{x} = P\tilde{\mathbf{x}}$.

**Similar matrices have the same characteristic polynomial, so the same eigenvalues.** Here are two proofs.

Through the eigenvalue equation: if $B\mathbf{x} = \lambda\mathbf{x}$, then $P^{-1}AP\mathbf{x} = \lambda\mathbf{x}$, so $A(P\mathbf{x}) = \lambda (P\mathbf{x})$. So $\lambda$ is an eigenvalue of $A$ with eigenvector $P\mathbf{x} \ne \mathbf{0}$. The eigenvalues are the same, but the eigenvectors differ by $P$.

Through the determinant:

$$
\det(B - \lambda I) = \det\bigl(P^{-1}(A - \lambda I)P\bigr) = \det(P^{-1})\,\det(A - \lambda I)\,\det(P) = \det(A - \lambda I),
$$

because $\det(XY) = \det X \det Y$ and $\det(P^{-1})\det(P) = \det(I) = 1$.

## 6. Four canonical forms

These four results tell us what a matrix "really looks like" after a good change of basis.

### The Schur normal form

**Theorem.** For every $A \in \mathbb{C}^{n \times n}$ there is a unitary $U$ such that $T = U^*AU$ is upper triangular.

Since $U^* = U^{-1}$, $A$ and $T$ are similar. So **every square matrix is similar to an upper triangular matrix**. They have the same eigenvalues, and the eigenvalues of a triangular matrix are its diagonal entries:

$$
f_A(\lambda) = f_T(\lambda) = (t_{11} - \lambda) \cdots (t_{nn} - \lambda).
$$

The determinant of a triangular matrix is the product of its diagonal. So $\det A = \det T = \lambda_1 \cdots \lambda_n$. This is the fact we used in Section 3.

### The principal axes theorem

**Theorem.** If $A$ is Hermitian ($A^* = A$), then there is a unitary $U$ with

$$
U^*AU = D = \mathrm{diag}(\lambda_1, \dots, \lambda_n).
$$

The columns of $U$ are eigenvectors of $A$, and they form an orthonormal basis. In words: **a symmetric (Hermitian) matrix always has a full set of orthogonal eigenvectors.** This is the content to remember, more than the name.

The proof for a real symmetric $A$ uses the Schur form. All eigenvalues are real, so we can take a real orthogonal $U$ with $T = U^tAU$ upper triangular. Now transpose:

$$
T^t = (U^tAU)^t = U^tA^tU = U^tAU = T.
$$

$T^t$ is lower triangular and $T$ is upper triangular, and they are equal. Only a diagonal matrix can be both. So $U^tAU = D$, and multiplying by $U$ on the left gives

$$
AU = UD.
$$

Now read this column by column. Multiplying by a diagonal matrix on the **right** scales the **columns**: column $k$ of $UD$ is $d_{kk}\mathbf{u}_k$. (On the left, it would scale the rows.) Column $k$ of $AU$ is $A\mathbf{u}_k$. So

$$
A\mathbf{u}_k = d_{kk}\,\mathbf{u}_k, \qquad k = 1, \dots, n.
$$

The diagonal entries of $D$ are the eigenvalues, and the orthonormal columns of $U$ are the eigenvectors. For the $2 \times 2$ example above, $U = \frac{1}{\sqrt 2}\begin{bmatrix} 1 & 1 \\ -1 & 1 \end{bmatrix}$ and $U^tAU = \mathrm{diag}(2, 4)$.

### The singular value decomposition (SVD)

**Theorem.** If $A$ is $n \times m$, there are unitary matrices $V$ ($n \times n$) and $U$ ($m \times m$) such that

$$
V^*AU = F,
$$

where $F$ is an $n \times m$ "diagonal" matrix: its diagonal holds $\mu_1 \ge \mu_2 \ge \cdots \ge \mu_r > 0$, and all other entries are zero.

The numbers $\mu_1, \dots, \mu_r$ are the **singular values** of $A$. They are real and positive, and $r$ is the rank of $A$. A rectangular matrix has no eigenvalues, so singular values take their place. From $A = VFU^*$ we get $A^*A = U(F^*F)U^*$, so the nonzero eigenvalues of $A^*A$ are $\mu_k^2$. The singular values are the square roots of the eigenvalues of $A^*A$.

### The Jordan canonical form

A **Jordan block** $J_k(\lambda)$ is a $k \times k$ matrix with $\lambda$ on the diagonal, $1$ just above the diagonal, and $0$ everywhere else.

**Theorem.** For every $n \times n$ matrix $A$ there is a nonsingular $P$ (not always unitary) with

$$
P^{-1}AP = \mathrm{diag}\bigl(J_{k_1}(\lambda_1), \dots, J_{k_s}(\lambda_s)\bigr).
$$

The $\lambda_i$ can repeat. Suppose one eigenvalue has a block of size 3 and a block of size 2. Then it is a root of multiplicity 5 of the characteristic polynomial, but it has only 2 independent eigenvectors, one per block. A matrix is diagonalizable only when all its blocks have size 1. A Hermitian matrix is like that, by the principal axes theorem.

## 7. Vector norms and matrix norms

A **vector norm** must satisfy three conditions. They are the same as for function norms:

1. $\|\mathbf{u}\| \ge 0$, and $\|\mathbf{u}\| = 0$ only for $\mathbf{u} = \mathbf{0}$,
2. $\|\alpha\mathbf{u}\| = |\alpha|\,\|\mathbf{u}\|$,
3. $\|\mathbf{u} + \mathbf{v}\| \le \|\mathbf{u}\| + \|\mathbf{v}\|$.

The $p$-norms are the discrete versions of the $L^p$ norms of functions:

$$
\|\mathbf{x}\|_p = \bigl(|x_1|^p + \cdots + |x_n|^p\bigr)^{1/p}, \qquad \|\mathbf{x}\|_\infty = \max_j |x_j| .
$$

For $\mathbf{x} = (3, -4, 0, 1)$: $\|\mathbf{x}\|_1 = 8$, $\|\mathbf{x}\|_2 = \sqrt{26} \approx 5.10$, $\|\mathbf{x}\|_\infty = 4$.

**In finite dimensions, all norms are equivalent:** for any two norms there are constants with

$$
C_1 \|\mathbf{x}\|_p \le \|\mathbf{x}\|_q \le C_2 \|\mathbf{x}\|_p .
$$

So if we prove that something goes to zero in one easy norm, it goes to zero in every norm. For example, $\|\mathbf{x}\|_\infty \le \|\mathbf{x}\|_2 \le \sqrt{n}\,\|\mathbf{x}\|_\infty$. This is **false in infinite dimensions**. On $C[0,1]$, the functions $f_k(x) = x^k$ all have maximum $1$, but $\int_0^1 |f_k|\,dx = \frac{1}{k+1} \to 0$. No constant can bound the max norm by the integral norm.

### Operator norms

A **matrix norm** satisfies the three conditions above and also

$$
\|AB\| \le \|A\|\,\|B\|, \qquad \|A\mathbf{x}\| \le \|A\|\,\|\mathbf{x}\| .
$$

We think of the matrix as an operator that acts on vectors: $A$ maps $\mathbf{x}$ to $A\mathbf{x}$. The **operator norm** measures the largest stretch:

$$
\|A\| = \sup_{\mathbf{x} \ne \mathbf{0}} \frac{\|A\mathbf{x}\|}{\|\mathbf{x}\|} .
$$

Note: $\|A\|_2$ is **not** $\sqrt{\sum_{i,j} |a_{ij}|^2}$. That is the Frobenius norm, which is a matrix norm but not an operator norm. (For the identity it gives $\sqrt{n}$, but every operator norm gives $\|I\| = 1$.)

The definition is not practical: to find the sup, we would have to try every $\mathbf{x}$. One test vector $\mathbf{x}'$ only gives a **lower bound**, $\|A\| \ge \|A\mathbf{x}'\|/\|\mathbf{x}'\|$. Luckily, the common operator norms have closed formulas. For $p = 1$ and $p = \infty$ they are simple sums of absolute values of the entries, and they are easy to look up. For $p = 2$:

$$
\|A\|_2 = \sqrt{\rho(A^*A)} ,
$$

the largest singular value of $A$.

### The spectral radius

The **spectrum** $\sigma(A)$ is the set of all eigenvalues of $A$. The **spectral radius** is

$$
\rho(A) = \max_{\lambda \in \sigma(A)} |\lambda| .
$$

Eigenvalues can be complex, so $\rho(A)$ is the radius of the smallest disk around $0$ that holds all of them.

**Theorem.** For every operator norm, $\rho(A) \le \|A\|$.

Proof: take $A\mathbf{x} = \lambda\mathbf{x}$. Then $|\lambda|\,\|\mathbf{x}\| = \|A\mathbf{x}\| \le \|A\|\,\|\mathbf{x}\|$, so $|\lambda| \le \|A\|$.

**Claim.** If $A$ is real symmetric, then $\|A\|_2 = \rho(A)$.

"$\ge$" is the theorem above. For "$\le$", take any $\mathbf{x}$. Because $A$ is symmetric, it has orthonormal eigenvectors $\mathbf{x}_1, \dots, \mathbf{x}_n$ (principal axes theorem). Order them so that $|\lambda_1| \ge |\lambda_2| \ge \cdots \ge |\lambda_n|$, and write

$$
\mathbf{x} = \alpha_1\mathbf{x}_1 + \cdots + \alpha_n\mathbf{x}_n, \qquad A\mathbf{x} = \alpha_1\lambda_1\mathbf{x}_1 + \cdots + \alpha_n\lambda_n\mathbf{x}_n .
$$

The eigenvectors are orthonormal, so all cross terms vanish:

$$
\frac{\|A\mathbf{x}\|_2^2}{\|\mathbf{x}\|_2^2}
= \frac{\alpha_1^2\lambda_1^2 + \cdots + \alpha_n^2\lambda_n^2}{\alpha_1^2 + \cdots + \alpha_n^2}
= \lambda_1^2 \cdot \frac{\alpha_1^2 + \alpha_2^2 (\lambda_2/\lambda_1)^2 + \cdots + \alpha_n^2 (\lambda_n/\lambda_1)^2}{\alpha_1^2 + \cdots + \alpha_n^2}
\le \lambda_1^2 .
$$

Each ratio $(\lambda_j/\lambda_1)^2$ is at most 1, so the fraction is at most 1. This bound does not depend on $\mathbf{x}$, so the sup is at most $\lambda_1^2$ too.

Where did we use symmetry? Only in the orthogonal eigenvectors, which killed the cross terms. Without symmetry, only $\|A\|_2 \ge \rho(A)$ holds, and the gap can be huge. For $\begin{bmatrix} 0.5 & 4 \\ 0 & 0.6 \end{bmatrix}$, $\rho = 0.6$ but $\|A\|_2 \approx 4.07$. (The same holds for Hermitian matrices: $A^*A = A^2$, so $\|A\|_2 = \sqrt{\rho(A^2)} = \rho(A)$.)

## 8. Powers of a matrix and the Neumann series

**Theorem.** For a square matrix $A$: $A^m \to 0$ as $m \to \infty$ $\iff$ $\rho(A) < 1$.

Why does this matter? In an iterative method for $A\mathbf{x} = \mathbf{b}$, the error is multiplied by the same matrix in every step. The error goes away only if the powers of that matrix go to zero.

The proof uses the Jordan form. If $A = PJP^{-1}$, then $A^m = PJ^mP^{-1}$, so we only need $J^m \to 0$. Look at one Jordan block and split it as $J = D + N$, where $D = \lambda I$ is the diagonal part and $N$ has the ones above the diagonal. Inside one block $D$ and $N$ commute, so the binomial formula works:

$$
J^m = (D + N)^m = \sum_{j=0}^{m} \binom{m}{j} D^{m-j} N^j .
$$

Two things happen:

- Each power of $N$ pushes the band of ones one step further up. For a $k \times k$ block, $N^k = 0$ ($N$ is **nilpotent**). So only the terms $j = 0, \dots, k-1$ survive.
- In each surviving term, $D^{m-j}$ has entries $\lambda^{m-j}$, and $\bigl|\binom{m}{j}\lambda^{m-j}\bigr| \le m^j |\lambda|^{m-j}$. If $|\lambda| < 1$, the geometric decay beats the polynomial growth $m^j$, and the term goes to 0.

Conversely, if $|\lambda| \ge 1$ for some eigenvalue, then $A^m\mathbf{x} = \lambda^m\mathbf{x}$ does not shrink.

Note the factor $m^j$. The powers do go to zero, but they can **grow first**. Experiment A shows this.

### The Neumann series

**Theorem.** If $\rho(A) < 1$, then $I - A$ is invertible and

$$
(I - A)^{-1} = I + A + A^2 + A^3 + \cdots .
$$

This is the matrix version of the geometric series $1 + r + r^2 + \cdots = \frac{1}{1 - r}$ for $|r| < 1$.

It is tempting to "prove" it like this: multiply $(I - A)(I + A + A^2 + \cdots)$ and cancel $A$ with $-A$, $A^2$ with $-A^2$, and so on. But be careful with infinite sums. With numbers and $a = 2$: $(1 + 2 + 4 + \cdots) - (2 + 4 + 8 + \cdots)$ is $\infty - \infty$, and cancelling term by term to get $1$ is wrong. We cannot split $\sum(a_k - b_k)$ into $\sum a_k - \sum b_k$ when the sums diverge.

The correct proof works with a finite sum. The cancellation there is exact:

$$
(I - A)(I + A + \cdots + A^m) = I - A^{m+1}.
$$

Since $\rho(A) < 1$, $A$ has no eigenvalue equal to 1, so $I - A$ is invertible. Then $I + A + \cdots + A^m = (I - A)^{-1}(I - A^{m+1})$, and $A^{m+1} \to 0$ gives the result.

**Theorem.** If $\|A\| < 1$ in some operator norm, then $I - A$ is invertible and

$$
\|(I - A)^{-1}\| \le \frac{1}{1 - \|A\|} .
$$

Proof: $\rho(A) \le \|A\| < 1$, so the Neumann series converges. Then use the triangle inequality, $\|I\| = 1$ and $\|A^k\| \le \|A\|^k$:

$$
\|(I - A)^{-1}\| \le 1 + \|A\| + \|A\|^2 + \cdots = \frac{1}{1 - \|A\|} .
$$

Here we need an operator norm, because the proof uses $\|I\| = 1$.

## 9. Gaussian elimination is LU decomposition

There are two families of methods for $A\mathbf{x} = \mathbf{b}$:

- **Direct methods.** They give the exact solution after a finite number of steps (up to rounding). All of them are based on Gaussian elimination.
- **Iterative methods.** They build a sequence of approximations. These come later.

Gaussian elimination changes $A\mathbf{x} = \mathbf{b}$ into an equivalent system $U\mathbf{x} = \mathbf{g}$, where $U$ is upper triangular, by row operations. Then it solves $U\mathbf{x} = \mathbf{g}$ from the bottom up (**back substitution**). The same row operations applied to $[A \mid I]$ until the left side becomes $I$ give $[I \mid A^{-1}]$ (the Gauss–Jordan method).

The key fact is that elimination secretly computes a factorization.

**Theorem (LU decomposition).** If elimination runs without row exchanges, then $A = LU$. Here $U$ is the upper triangular matrix from elimination, and $L$ is lower triangular with ones on the diagonal. Below the diagonal, $L$ holds the multipliers $m_{ik}$ that elimination used to clear entry $(i, k)$.

A small example:

$$
A = \begin{bmatrix} 4 & 2 & 1 \\ 2 & 5 & 2 \\ 1 & 2 & 6 \end{bmatrix}
= \underbrace{\begin{bmatrix} 1 & 0 & 0 \\ 0.5 & 1 & 0 \\ 0.25 & 0.375 & 1 \end{bmatrix}}_{L}
\underbrace{\begin{bmatrix} 4 & 2 & 1 \\ 0 & 4 & 1.5 \\ 0 & 0 & 5.1875 \end{bmatrix}}_{U} .
$$

**Corollary.** $\det L = 1$, so $\det A = \det U = u_{11} u_{22} \cdots u_{nn}$. Here $\det A = 4 \cdot 4 \cdot 5.1875 = 83$.

Variants of the same process give the QR and Cholesky factorizations. (Cholesky is in Section 11.)

### The cost: $O(n^3)$

Count the work. To clear the first column, each of the $n - 1$ rows below the first needs one multiplier, and then about $2n$ operations (a multiply and a subtract per entry). That is about $2n^2$ for the first column. The next step works on an $(n-1) \times (n-1)$ block, and so on. In total,

$$
\sum_{k=1}^{n} 2k^2 \approx 2 \cdot \frac{n^3}{3} = \frac{2}{3}n^3 .
$$

This is **extremely expensive** for large $n$. A dense system with a million unknowns, which is not big for a PDE, needs about $\frac{2}{3} \times 10^{18}$ operations. That is the scale of the fastest supercomputers today (about $10^{18}$ operations per second). Storing the matrix would also need $10^{12}$ numbers, or 8 TB. So for big problems we need methods that use the structure of the matrix and cost much less, like $O(n \log n)$.

### When LU pays off: one $A$, many $\mathbf{b}$

If $A = LU$, then $A\mathbf{x} = \mathbf{b}$ splits into two triangular systems:

$$
L\mathbf{y} = \mathbf{b}, \qquad U\mathbf{x} = \mathbf{y}.
$$

- $L\mathbf{y} = \mathbf{b}$: the first row gives $y_1$, the second row then gives $y_2$, and so on (**forward substitution**).
- $U\mathbf{x} = \mathbf{y}$: start from the last row (**back substitution**).

Each triangular solve costs about $n^2$ operations, so about $2n^2$ in total. If $A$ stays the same and we have a thousand right-hand sides, we factor once ($\frac{2}{3}n^3$) and then only substitute.

## 10. Pivoting: what it is for

Elimination divides by the pivots $a_{kk}^{(k)}$. If a pivot is zero, elimination stops. If a pivot is very small, the multipliers become huge, and rounding errors grow. **Pivoting** fixes this by exchanging rows (or columns) so that a large entry becomes the pivot.

- **Row pivoting.** At step $k$, pick the row with the largest $|a_{ik}^{(k)}|$ in column $k$ and swap it with row $k$, in both $A$ and $\mathbf{b}$. This only changes the order of the equations, so $\mathbf{x}$ stays the same. Now every multiplier has $|m_{ik}| \le 1$, the entries do not blow up, and less precision is lost. In matrix form the swaps are a permutation matrix $P$, and $PA = LU$.
- **Column pivoting.** Swap columns instead. This changes the order of the unknowns, so we must also swap the entries of $\mathbf{x}$; $\mathbf{b}$ stays the same.

Pivoting is related to preconditioning: we change the system into an equivalent one that is nicer to compute with. Even a simple scaling by the diagonal, $D^{-1}A\mathbf{x} = D^{-1}\mathbf{b}$, can make a big difference. In practice nobody pivots by hand: linear algebra software does it automatically. What matters is to understand its role. Experiment D shows the effect.

## 11. Cholesky for symmetric matrices

If $A$ is symmetric, the factorization can be symmetric too. Suppose $A = LU$ with ones on the diagonal of $L$, and let $D$ be the diagonal of $U$. Then $U = DL^t$ (because $A = A^t$ and the factorization is unique), so

$$
A = LDL^t .
$$

If all pivots are positive, we can take square roots and write $A = (L\sqrt{D})(L\sqrt{D})^t = CC^t$. This is the **Cholesky factorization**. The pivots are all positive exactly when $A$ is **symmetric positive definite**. For the example in Section 9 (pivots $4, 4, 5.1875$):

$$
C = L\sqrt{D} \approx \begin{bmatrix} 2 & 0 & 0 \\ 1 & 2 & 0 \\ 0.5 & 0.75 & 2.278 \end{bmatrix}, \qquad CC^t = A .
$$

Cholesky needs about half the work of LU, because it uses the symmetry.

## 12. Error and residual

Let $\mathbf{x}_A$ be an approximate solution of $A\mathbf{x} = \mathbf{b}$. Define

$$
\mathbf{e} = \mathbf{x} - \mathbf{x}_A \quad (\text{error}), \qquad \mathbf{r} = \mathbf{b} - A\mathbf{x}_A \quad (\text{residual}).
$$

Subtracting $A\mathbf{x}_A$ from $A\mathbf{x} = \mathbf{b}$ gives

$$
A\mathbf{e} = \mathbf{r} .
$$

(Take care with the sign: if we define the residual the other way, as $A\mathbf{x}_A - \mathbf{b}$, we get $A\mathbf{e} = -\mathbf{r}$.)

- The **error** is what we want, but we cannot compute it, because we do not know $\mathbf{x}$.
- The **residual** can always be computed. So we use it to judge an approximation.

But is a small residual a proof of a small error? Only if the matrix is "nice". If $A = I$, the residual equals the error. If $A$ is singular, like $\begin{bmatrix} 1 & 1 \\ 1 & 1 \end{bmatrix}$ in Section 3, then $\mathbf{x}_A = \mathbf{x} + (10^4, -10^4)$ has residual zero and error $10^4$. Experiment E shows the same thing for a matrix that is not singular, only close to singular. How much a matrix can amplify the residual is measured by the **condition number**, the topic of next week.

## 13. Experiments

I wrote one Python script with five small experiments, using only the standard library.

### A. Powers of a matrix: $A^m \to 0$, but not always at once

I computed $\|A^m\|_2$ for $m = 0, \dots, 120$ for four matrices:

- a symmetric $3 \times 3$ matrix with $\rho = 0.924$,
- a $4 \times 4$ Jordan block with $\lambda = 0.8$,
- the triangular matrix $\begin{bmatrix} 0.5 & 4 \\ 0 & 0.6 \end{bmatrix}$ from Section 7, with $\rho = 0.6$,
- a symmetric $2 \times 2$ matrix with eigenvalues $1.05$ and $0.75$, so $\rho = 1.05$.

<figure class="viz">
<div class="viz-legend"><span class="viz-key"><span class="viz-swatch" style="background:#2a78d6"></span>symmetric, ρ = 0.924</span><span class="viz-key"><span class="viz-swatch" style="background:#eb6834"></span>Jordan block, ρ = 0.8</span><span class="viz-key"><span class="viz-swatch" style="background:#1baf7a"></span>triangular, ρ = 0.6</span><span class="viz-key"><span class="viz-swatch" style="background:#eda100"></span>symmetric, ρ = 1.05</span></div>
<div class="viz-frame">
<svg viewBox="0 0 680 400" role="img" tabindex="0" aria-label="The 2-norm of A to the power m, for m from 0 to 120, on a log scale, for four matrices. The symmetric matrix with spectral radius 0.924 decays steadily. The Jordan block with spectral radius 0.8 first grows to about 33 near m = 14 and then decays. The triangular matrix with spectral radius 0.6 grows a little and then decays fast. The symmetric matrix with spectral radius 1.05 grows without bound." class="viz-plot"><line x1="64" x2="560" y1="30.0" y2="30.0" class="viz-grid"/><text x="54" y="34.0" text-anchor="end" class="viz-tick">10<tspan dy="-0.5em" font-size="75%">3</tspan></text><line x1="64" x2="560" y1="58.2" y2="58.2" class="viz-grid"/><text x="54" y="62.2" text-anchor="end" class="viz-tick">10<tspan dy="-0.5em" font-size="75%">2</tspan></text><line x1="64" x2="560" y1="86.4" y2="86.4" class="viz-grid"/><text x="54" y="90.4" text-anchor="end" class="viz-tick">10<tspan dy="-0.5em" font-size="75%">1</tspan></text><line x1="64" x2="560" y1="114.5" y2="114.5" class="viz-grid"/><text x="54" y="118.5" text-anchor="end" class="viz-tick">1</text><line x1="64" x2="560" y1="142.7" y2="142.7" class="viz-grid"/><text x="54" y="146.7" text-anchor="end" class="viz-tick">10<tspan dy="-0.5em" font-size="75%">−1</tspan></text><line x1="64" x2="560" y1="170.9" y2="170.9" class="viz-grid"/><text x="54" y="174.9" text-anchor="end" class="viz-tick">10<tspan dy="-0.5em" font-size="75%">−2</tspan></text><line x1="64" x2="560" y1="199.1" y2="199.1" class="viz-grid"/><text x="54" y="203.1" text-anchor="end" class="viz-tick">10<tspan dy="-0.5em" font-size="75%">−3</tspan></text><line x1="64" x2="560" y1="227.3" y2="227.3" class="viz-grid"/><text x="54" y="231.3" text-anchor="end" class="viz-tick">10<tspan dy="-0.5em" font-size="75%">−4</tspan></text><line x1="64" x2="560" y1="255.5" y2="255.5" class="viz-grid"/><text x="54" y="259.5" text-anchor="end" class="viz-tick">10<tspan dy="-0.5em" font-size="75%">−5</tspan></text><line x1="64" x2="560" y1="283.6" y2="283.6" class="viz-grid"/><text x="54" y="287.6" text-anchor="end" class="viz-tick">10<tspan dy="-0.5em" font-size="75%">−6</tspan></text><line x1="64" x2="560" y1="311.8" y2="311.8" class="viz-grid"/><text x="54" y="315.8" text-anchor="end" class="viz-tick">10<tspan dy="-0.5em" font-size="75%">−7</tspan></text><line x1="64" x2="560" y1="340.0" y2="340.0" class="viz-grid"/><text x="54" y="344.0" text-anchor="end" class="viz-tick">10<tspan dy="-0.5em" font-size="75%">−8</tspan></text><text x="64.0" y="362" text-anchor="middle" class="viz-tick">0</text><text x="146.7" y="362" text-anchor="middle" class="viz-tick">20</text><text x="229.3" y="362" text-anchor="middle" class="viz-tick">40</text><text x="312.0" y="362" text-anchor="middle" class="viz-tick">60</text><text x="394.7" y="362" text-anchor="middle" class="viz-tick">80</text><text x="477.3" y="362" text-anchor="middle" class="viz-tick">100</text><text x="560.0" y="362" text-anchor="middle" class="viz-tick">120</text><text x="312.0" y="388" text-anchor="middle" class="viz-axis-title">Power m</text><text x="64" y="14" class="viz-axis-title">‖Aᵐ‖₂</text><path d="M64.0,114.5 L68.1,115.5 L72.3,116.5 L76.4,117.4 L80.5,118.4 L84.7,119.4 L88.8,120.3 L92.9,121.3 L97.1,122.3 L101.2,123.2 L105.3,124.2 L109.5,125.1 L113.6,126.1 L117.7,127.1 L121.9,128.0 L126.0,129.0 L130.1,130.0 L134.3,130.9 L138.4,131.9 L142.5,132.9 L146.7,133.8 L150.8,134.8 L154.9,135.8 L159.1,136.7 L163.2,137.7 L167.3,138.6 L171.5,139.6 L175.6,140.6 L179.7,141.5 L183.9,142.5 L188.0,143.5 L192.1,144.4 L196.3,145.4 L200.4,146.4 L204.5,147.3 L208.7,148.3 L212.8,149.2 L216.9,150.2 L221.1,151.2 L225.2,152.1 L229.3,153.1 L233.5,154.1 L237.6,155.0 L241.7,156.0 L245.9,157.0 L250.0,157.9 L254.1,158.9 L258.3,159.9 L262.4,160.8 L266.5,161.8 L270.7,162.7 L274.8,163.7 L278.9,164.7 L283.1,165.6 L287.2,166.6 L291.3,167.6 L295.5,168.5 L299.6,169.5 L303.7,170.5 L307.9,171.4 L312.0,172.4 L316.1,173.3 L320.3,174.3 L324.4,175.3 L328.5,176.2 L332.7,177.2 L336.8,178.2 L340.9,179.1 L345.1,180.1 L349.2,181.1 L353.3,182.0 L357.5,183.0 L361.6,183.9 L365.7,184.9 L369.9,185.9 L374.0,186.8 L378.1,187.8 L382.3,188.8 L386.4,189.7 L390.5,190.7 L394.7,191.7 L398.8,192.6 L402.9,193.6 L407.1,194.6 L411.2,195.5 L415.3,196.5 L419.5,197.4 L423.6,198.4 L427.7,199.4 L431.9,200.3 L436.0,201.3 L440.1,202.3 L444.3,203.2 L448.4,204.2 L452.5,205.2 L456.7,206.1 L460.8,207.1 L464.9,208.0 L469.1,209.0 L473.2,210.0 L477.3,210.9 L481.5,211.9 L485.6,212.9 L489.7,213.8 L493.9,214.8 L498.0,215.8 L502.1,216.7 L506.3,217.7 L510.4,218.6 L514.5,219.6 L518.7,220.6 L522.8,221.5 L526.9,222.5 L531.1,223.5 L535.2,224.4 L539.3,225.4 L543.5,226.4 L547.6,227.3 L551.7,228.3 L555.9,229.3 L560.0,230.2" fill="none" stroke="#2a78d6" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><circle cx="560.0" cy="230.2" r="4" fill="#2a78d6" stroke="#ffffff" stroke-width="2"/><path d="M64.0,114.5 L68.1,108.1 L72.3,101.9 L76.4,96.0 L80.5,90.7 L84.7,86.2 L88.8,82.5 L92.9,79.5 L97.1,77.1 L101.2,75.4 L105.3,74.0 L109.5,73.1 L113.6,72.4 L117.7,72.1 L121.9,71.9 L126.0,72.0 L130.1,72.2 L134.3,72.6 L138.4,73.2 L142.5,73.8 L146.7,74.6 L150.8,75.5 L154.9,76.4 L159.1,77.5 L163.2,78.6 L167.3,79.8 L171.5,81.0 L175.6,82.3 L179.7,83.7 L183.9,85.1 L188.0,86.5 L192.1,88.0 L196.3,89.5 L200.4,91.1 L204.5,92.7 L208.7,94.3 L212.8,96.0 L216.9,97.7 L221.1,99.4 L225.2,101.2 L229.3,103.0 L233.5,104.8 L237.6,106.6 L241.7,108.4 L245.9,110.3 L250.0,112.2 L254.1,114.1 L258.3,116.0 L262.4,118.0 L266.5,119.9 L270.7,121.9 L274.8,123.9 L278.9,125.9 L283.1,127.9 L287.2,130.0 L291.3,132.0 L295.5,134.1 L299.6,136.1 L303.7,138.2 L307.9,140.3 L312.0,142.4 L316.1,144.5 L320.3,146.7 L324.4,148.8 L328.5,150.9 L332.7,153.1 L336.8,155.2 L340.9,157.4 L345.1,159.6 L349.2,161.8 L353.3,164.0 L357.5,166.2 L361.6,168.4 L365.7,170.6 L369.9,172.8 L374.0,175.1 L378.1,177.3 L382.3,179.6 L386.4,181.8 L390.5,184.1 L394.7,186.3 L398.8,188.6 L402.9,190.9 L407.1,193.2 L411.2,195.4 L415.3,197.7 L419.5,200.0 L423.6,202.3 L427.7,204.6 L431.9,206.9 L436.0,209.3 L440.1,211.6 L444.3,213.9 L448.4,216.2 L452.5,218.6 L456.7,220.9 L460.8,223.3 L464.9,225.6 L469.1,228.0 L473.2,230.3 L477.3,232.7 L481.5,235.0 L485.6,237.4 L489.7,239.8 L493.9,242.1 L498.0,244.5 L502.1,246.9 L506.3,249.3 L510.4,251.7 L514.5,254.0 L518.7,256.4 L522.8,258.8 L526.9,261.2 L531.1,263.6 L535.2,266.0 L539.3,268.5 L543.5,270.9 L547.6,273.3 L551.7,275.7 L555.9,278.1 L560.0,280.5" fill="none" stroke="#eb6834" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><circle cx="560.0" cy="280.5" r="4" fill="#eb6834" stroke="#ffffff" stroke-width="2"/><path d="M64.0,114.5 L68.1,97.4 L72.3,96.4 L76.4,98.7 L80.5,102.4 L84.7,106.9 L88.8,111.9 L92.9,117.2 L97.1,122.6 L101.2,128.3 L105.3,134.1 L109.5,139.9 L113.6,145.9 L117.7,151.9 L121.9,157.9 L126.0,164.0 L130.1,170.1 L134.3,176.2 L138.4,182.4 L142.5,188.6 L146.7,194.8 L150.8,201.0 L154.9,207.2 L159.1,213.4 L163.2,219.6 L167.3,225.8 L171.5,232.1 L175.6,238.3 L179.7,244.5 L183.9,250.8 L188.0,257.0 L192.1,263.3 L196.3,269.5 L200.4,275.7 L204.5,282.0 L208.7,288.2 L212.8,294.5 L216.9,300.7 L221.1,307.0 L225.2,313.2 L229.3,319.5 L233.5,325.7 L237.6,332.0 L241.7,338.2 L245.9,340.0 L250.0,340.0 L254.1,340.0 L258.3,340.0 L262.4,340.0 L266.5,340.0 L270.7,340.0 L274.8,340.0 L278.9,340.0 L283.1,340.0 L287.2,340.0 L291.3,340.0 L295.5,340.0 L299.6,340.0 L303.7,340.0 L307.9,340.0 L312.0,340.0 L316.1,340.0 L320.3,340.0 L324.4,340.0 L328.5,340.0 L332.7,340.0 L336.8,340.0 L340.9,340.0 L345.1,340.0 L349.2,340.0 L353.3,340.0 L357.5,340.0 L361.6,340.0 L365.7,340.0 L369.9,340.0 L374.0,340.0 L378.1,340.0 L382.3,340.0 L386.4,340.0 L390.5,340.0 L394.7,340.0 L398.8,340.0 L402.9,340.0 L407.1,340.0 L411.2,340.0 L415.3,340.0 L419.5,340.0 L423.6,340.0 L427.7,340.0 L431.9,340.0 L436.0,340.0 L440.1,340.0 L444.3,340.0 L448.4,340.0 L452.5,340.0 L456.7,340.0 L460.8,340.0 L464.9,340.0 L469.1,340.0 L473.2,340.0 L477.3,340.0 L481.5,340.0 L485.6,340.0 L489.7,340.0 L493.9,340.0 L498.0,340.0 L502.1,340.0 L506.3,340.0 L510.4,340.0 L514.5,340.0 L518.7,340.0 L522.8,340.0 L526.9,340.0 L531.1,340.0 L535.2,340.0 L539.3,340.0 L543.5,340.0 L547.6,340.0 L551.7,340.0 L555.9,340.0 L560.0,340.0" fill="none" stroke="#1baf7a" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><circle cx="560.0" cy="340.0" r="4" fill="#1baf7a" stroke="#ffffff" stroke-width="2"/><path d="M64.0,114.5 L68.1,113.9 L72.3,113.4 L76.4,112.8 L80.5,112.2 L84.7,111.6 L88.8,111.0 L92.9,110.4 L97.1,109.8 L101.2,109.2 L105.3,108.6 L109.5,108.0 L113.6,107.4 L117.7,106.8 L121.9,106.2 L126.0,105.6 L130.1,105.0 L134.3,104.4 L138.4,103.8 L142.5,103.2 L146.7,102.6 L150.8,102.0 L154.9,101.4 L159.1,100.8 L163.2,100.2 L167.3,99.6 L171.5,99.0 L175.6,98.4 L179.7,97.8 L183.9,97.2 L188.0,96.6 L192.1,96.0 L196.3,95.4 L200.4,94.8 L204.5,94.2 L208.7,93.6 L212.8,93.0 L216.9,92.5 L221.1,91.9 L225.2,91.3 L229.3,90.7 L233.5,90.1 L237.6,89.5 L241.7,88.9 L245.9,88.3 L250.0,87.7 L254.1,87.1 L258.3,86.5 L262.4,85.9 L266.5,85.3 L270.7,84.7 L274.8,84.1 L278.9,83.5 L283.1,82.9 L287.2,82.3 L291.3,81.7 L295.5,81.1 L299.6,80.5 L303.7,79.9 L307.9,79.3 L312.0,78.7 L316.1,78.1 L320.3,77.5 L324.4,76.9 L328.5,76.3 L332.7,75.7 L336.8,75.1 L340.9,74.5 L345.1,73.9 L349.2,73.3 L353.3,72.7 L357.5,72.1 L361.6,71.6 L365.7,71.0 L369.9,70.4 L374.0,69.8 L378.1,69.2 L382.3,68.6 L386.4,68.0 L390.5,67.4 L394.7,66.8 L398.8,66.2 L402.9,65.6 L407.1,65.0 L411.2,64.4 L415.3,63.8 L419.5,63.2 L423.6,62.6 L427.7,62.0 L431.9,61.4 L436.0,60.8 L440.1,60.2 L444.3,59.6 L448.4,59.0 L452.5,58.4 L456.7,57.8 L460.8,57.2 L464.9,56.6 L469.1,56.0 L473.2,55.4 L477.3,54.8 L481.5,54.2 L485.6,53.6 L489.7,53.0 L493.9,52.4 L498.0,51.8 L502.1,51.2 L506.3,50.7 L510.4,50.1 L514.5,49.5 L518.7,48.9 L522.8,48.3 L526.9,47.7 L531.1,47.1 L535.2,46.5 L539.3,45.9 L543.5,45.3 L547.6,44.7 L551.7,44.1 L555.9,43.5 L560.0,42.9" fill="none" stroke="#eda100" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/><circle cx="560.0" cy="42.9" r="4" fill="#eda100" stroke="#ffffff" stroke-width="2"/><text x="570.0" y="234.2" class="viz-label">sym.</text><text x="570.0" y="284.5" class="viz-label">Jordan</text><text x="158.7" y="198.8" class="viz-label">triangular</text><text x="570.0" y="46.9" class="viz-label">ρ &gt; 1</text><line x1="0" x2="0" y1="30" y2="340" class="viz-crosshair" visibility="hidden"/></svg>
<div class="viz-tooltip" hidden></div>
</div>
<figcaption>The 2-norm of the power A<sup>m</sup> for four matrices, on a log scale. Values below 10<sup>−8</sup> are drawn on the bottom edge (the table has the real values). Hover over the chart, or focus it and use the arrow keys, to read the values.</figcaption>
<details class="viz-table"><summary>Show the numbers</summary><table><thead><tr><th>m</th><th>symmetric, ρ = 0.924</th><th>Jordan block, ρ = 0.8</th><th>triangular, ρ = 0.6</th><th>symmetric, ρ = 1.05</th></tr></thead><tbody><tr><td>0</td><td>1</td><td>1</td><td>1</td><td>1</td></tr><tr><td>1</td><td>0.924</td><td>1.69</td><td>4.07</td><td>1.05</td></tr><tr><td>2</td><td>0.854</td><td>2.81</td><td>4.42</td><td>1.1</td></tr><tr><td>5</td><td>0.674</td><td>10.1</td><td>1.86</td><td>1.28</td></tr><tr><td>10</td><td>0.455</td><td>27.4</td><td>0.203</td><td>1.63</td></tr><tr><td>14</td><td>0.332</td><td>32.5</td><td>0.0289</td><td>1.98</td></tr><tr><td>20</td><td>0.207</td><td>26.1</td><td>1.42e-3</td><td>2.65</td></tr><tr><td>30</td><td>0.0942</td><td>9.89</td><td>8.81e-6</td><td>4.32</td></tr><tr><td>40</td><td>0.0428</td><td>2.58</td><td>5.35e-8</td><td>7.04</td></tr><tr><td>60</td><td>8.87e-3</td><td>0.103</td><td>1.96e-12</td><td>18.7</td></tr><tr><td>80</td><td>1.84e-3</td><td>2.84e-3</td><td>7.15e-17</td><td>49.6</td></tr><tr><td>100</td><td>3.80e-4</td><td>6.44e-5</td><td>2.61e-21</td><td>132</td></tr><tr><td>120</td><td>7.86e-5</td><td>1.29e-6</td><td>9.56e-26</td><td>349</td></tr></tbody></table></details>
<script>(() => {const fig = document.currentScript.closest('figure');const svg = fig.querySelector('svg');const line = svg.querySelector('.viz-crosshair');const tip = fig.querySelector('.viz-tooltip');const names = ["symmetric", "Jordan", "triangular", "ρ = 1.05"]; const colors = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"];const rows = [["1", "1", "1", "1"], ["0.924", "1.69", "4.07", "1.05"], ["0.854", "2.81", "4.42", "1.1"], ["0.79", "4.54", "3.65", "1.16"], ["0.73", "7", "2.69", "1.22"], ["0.674", "10.1", "1.86", "1.28"], ["0.623", "13.8", "1.24", "1.34"], ["0.576", "17.6", "0.808", "1.41"], ["0.533", "21.2", "0.516", "1.48"], ["0.492", "24.6", "0.325", "1.55"], ["0.455", "27.4", "0.203", "1.63"], ["0.42", "29.6", "0.126", "1.71"], ["0.389", "31.2", "0.0773", "1.8"], ["0.359", "32.2", "0.0474", "1.89"], ["0.332", "32.5", "0.0289", "1.98"], ["0.307", "32.3", "0.0176", "2.08"], ["0.284", "31.7", "0.0107", "2.18"], ["0.262", "30.7", "6.47e-3", "2.29"], ["0.242", "29.4", "3.91e-3", "2.41"], ["0.224", "27.8", "2.36e-3", "2.53"], ["0.207", "26.1", "1.42e-3", "2.65"], ["0.191", "24.3", "8.59e-4", "2.79"], ["0.177", "22.5", "5.17e-4", "2.93"], ["0.163", "20.7", "3.11e-4", "3.07"], ["0.151", "18.9", "1.87e-4", "3.23"], ["0.14", "17.2", "1.13e-4", "3.39"], ["0.129", "15.5", "6.77e-5", "3.56"], ["0.119", "13.9", "4.07e-5", "3.73"], ["0.11", "12.5", "2.44e-5", "3.92"], ["0.102", "11.1", "1.47e-5", "4.12"], ["0.0942", "9.89", "8.81e-6", "4.32"], ["0.087", "8.75", "5.29e-6", "4.54"], ["0.0804", "7.72", "3.18e-6", "4.76"], ["0.0743", "6.79", "1.91e-6", "5"], ["0.0687", "5.96", "1.14e-6", "5.25"], ["0.0635", "5.21", "6.87e-7", "5.52"], ["0.0587", "4.55", "4.12e-7", "5.79"], ["0.0543", "3.96", "2.47e-7", "6.08"], ["0.0501", "3.44", "1.48e-7", "6.39"], ["0.0463", "2.98", "8.91e-8", "6.7"], ["0.0428", "2.58", "5.35e-8", "7.04"], ["0.0396", "2.22", "3.21e-8", "7.39"], ["0.0366", "1.91", "1.92e-8", "7.76"], ["0.0338", "1.65", "1.15e-8", "8.15"], ["0.0313", "1.41", "6.93e-9", "8.56"], ["0.0289", "1.21", "4.16e-9", "8.99"], ["0.0267", "1.04", "2.49e-9", "9.43"], ["0.0247", "0.885", "1.50e-9", "9.91"], ["0.0228", "0.755", "8.98e-10", "10.4"], ["0.0211", "0.644", "5.39e-10", "10.9"], ["0.0195", "0.548", "3.23e-10", "11.5"], ["0.018", "0.466", "1.94e-10", "12"], ["0.0166", "0.395", "1.16e-10", "12.6"], ["0.0154", "0.335", "6.99e-11", "13.3"], ["0.0142", "0.284", "4.19e-11", "13.9"], ["0.0131", "0.24", "2.51e-11", "14.6"], ["0.0122", "0.203", "1.51e-11", "15.4"], ["0.0112", "0.171", "9.05e-12", "16.1"], ["0.0104", "0.145", "5.43e-12", "16.9"], ["9.59e-3", "0.122", "3.26e-12", "17.8"], ["8.87e-3", "0.103", "1.96e-12", "18.7"], ["8.20e-3", "0.0863", "1.17e-12", "19.6"], ["7.57e-3", "0.0726", "7.04e-13", "20.6"], ["7.00e-3", "0.061", "4.22e-13", "21.6"], ["6.47e-3", "0.0512", "2.53e-13", "22.7"], ["5.98e-3", "0.0429", "1.52e-13", "23.8"], ["5.53e-3", "0.036", "9.12e-14", "25"], ["5.11e-3", "0.0301", "5.47e-14", "26.3"], ["4.72e-3", "0.0252", "3.28e-14", "27.6"], ["4.36e-3", "0.0211", "1.97e-14", "29"], ["4.03e-3", "0.0176", "1.18e-14", "30.4"], ["3.73e-3", "0.0147", "7.09e-15", "31.9"], ["3.45e-3", "0.0123", "4.26e-15", "33.5"], ["3.19e-3", "0.0102", "2.55e-15", "35.2"], ["2.94e-3", "8.54e-3", "1.53e-15", "37"], ["2.72e-3", "7.12e-3", "9.19e-16", "38.8"], ["2.51e-3", "5.93e-3", "5.52e-16", "40.8"], ["2.32e-3", "4.94e-3", "3.31e-16", "42.8"], ["2.15e-3", "4.11e-3", "1.99e-16", "45"], ["1.99e-3", "3.41e-3", "1.19e-16", "47.2"], ["1.84e-3", "2.84e-3", "7.15e-17", "49.6"], ["1.70e-3", "2.36e-3", "4.29e-17", "52"], ["1.57e-3", "1.96e-3", "2.57e-17", "54.6"], ["1.45e-3", "1.62e-3", "1.54e-17", "57.4"], ["1.34e-3", "1.35e-3", "9.27e-18", "60.2"], ["1.24e-3", "1.12e-3", "5.56e-18", "63.3"], ["1.14e-3", "9.27e-4", "3.34e-18", "66.4"], ["1.06e-3", "7.68e-4", "2.00e-18", "69.7"], ["9.77e-4", "6.36e-4", "1.20e-18", "73.2"], ["9.03e-4", "5.26e-4", "7.21e-19", "76.9"], ["8.35e-4", "4.36e-4", "4.32e-19", "80.7"], ["7.72e-4", "3.60e-4", "2.59e-19", "84.8"], ["7.13e-4", "2.98e-4", "1.56e-19", "89"], ["6.59e-4", "2.46e-4", "9.34e-20", "93.5"], ["6.09e-4", "2.04e-4", "5.60e-20", "98.1"], ["5.63e-4", "1.68e-4", "3.36e-20", "103"], ["5.21e-4", "1.39e-4", "2.02e-20", "108"], ["4.81e-4", "1.15e-4", "1.21e-20", "114"], ["4.45e-4", "9.46e-5", "7.26e-21", "119"], ["4.11e-4", "7.81e-5", "4.36e-21", "125"], ["3.80e-4", "6.44e-5", "2.61e-21", "132"], ["3.51e-4", "5.31e-5", "1.57e-21", "138"], ["3.24e-4", "4.37e-5", "9.41e-22", "145"], ["3.00e-4", "3.60e-5", "5.65e-22", "152"], ["2.77e-4", "2.97e-5", "3.39e-22", "160"], ["2.56e-4", "2.45e-5", "2.03e-22", "168"], ["2.37e-4", "2.01e-5", "1.22e-22", "176"], ["2.19e-4", "1.66e-5", "7.32e-23", "185"], ["2.02e-4", "1.36e-5", "4.39e-23", "194"], ["1.87e-4", "1.12e-5", "2.63e-23", "204"], ["1.73e-4", "9.22e-6", "1.58e-23", "214"], ["1.60e-4", "7.58e-6", "9.48e-24", "225"], ["1.48e-4", "6.23e-6", "5.69e-24", "236"], ["1.36e-4", "5.12e-6", "3.41e-24", "248"], ["1.26e-4", "4.21e-6", "2.05e-24", "260"], ["1.17e-4", "3.46e-6", "1.23e-24", "273"], ["1.08e-4", "2.84e-6", "7.37e-25", "287"], ["9.96e-5", "2.33e-6", "4.42e-25", "301"], ["9.20e-5", "1.91e-6", "2.65e-25", "316"], ["8.51e-5", "1.57e-6", "1.59e-25", "332"], ["7.86e-5", "1.29e-6", "9.56e-26", "349"]]; const ms = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120];const width = 680, last = ms.length - 1;const xs = [64.0, 68.1, 72.3, 76.4, 80.5, 84.7, 88.8, 92.9, 97.1, 101.2, 105.3, 109.5, 113.6, 117.7, 121.9, 126.0, 130.1, 134.3, 138.4, 142.5, 146.7, 150.8, 154.9, 159.1, 163.2, 167.3, 171.5, 175.6, 179.7, 183.9, 188.0, 192.1, 196.3, 200.4, 204.5, 208.7, 212.8, 216.9, 221.1, 225.2, 229.3, 233.5, 237.6, 241.7, 245.9, 250.0, 254.1, 258.3, 262.4, 266.5, 270.7, 274.8, 278.9, 283.1, 287.2, 291.3, 295.5, 299.6, 303.7, 307.9, 312.0, 316.1, 320.3, 324.4, 328.5, 332.7, 336.8, 340.9, 345.1, 349.2, 353.3, 357.5, 361.6, 365.7, 369.9, 374.0, 378.1, 382.3, 386.4, 390.5, 394.7, 398.8, 402.9, 407.1, 411.2, 415.3, 419.5, 423.6, 427.7, 431.9, 436.0, 440.1, 444.3, 448.4, 452.5, 456.7, 460.8, 464.9, 469.1, 473.2, 477.3, 481.5, 485.6, 489.7, 493.9, 498.0, 502.1, 506.3, 510.4, 514.5, 518.7, 522.8, 526.9, 531.1, 535.2, 539.3, 543.5, 547.6, 551.7, 555.9, 560.0];let current = -1;const show = (k) => { k = Math.max(0, Math.min(last, k)); current = k; const x = xs[k]; line.setAttribute('x1', x); line.setAttribute('x2', x); line.setAttribute('visibility', 'visible'); tip.replaceChildren(); const head = document.createElement('div'); head.className = 'viz-tooltip-head'; head.textContent = 'm = ' + ms[k]; tip.append(head); names.forEach((name, i) => {  const row = document.createElement('div'); row.className = 'viz-tooltip-row';  const key = document.createElement('span'); key.className = 'viz-swatch'; key.style.background = colors[i];  const value = document.createElement('strong'); value.textContent = rows[k][i];  const label = document.createElement('span'); label.textContent = name;  row.append(key, value, label); tip.append(row); }); const box = svg.getBoundingClientRect(), frame = fig.getBoundingClientRect(); const xp = box.left - frame.left + x * box.width / width; tip.hidden = false; const room = frame.width - tip.offsetWidth - 8; tip.style.left = Math.max(0, Math.min(room, xp + 12)) + 'px'; tip.style.top = (box.top - frame.top + 8) + 'px';};const hide = () => { line.setAttribute('visibility', 'hidden'); tip.hidden = true; current = -1; };svg.addEventListener('pointermove', (event) => { const box = svg.getBoundingClientRect(); const x = (event.clientX - box.left) * width / box.width; let best = 0; xs.forEach((v, i) => { if (Math.abs(v - x) < Math.abs(xs[best] - x)) best = i; }); show(best);});svg.addEventListener('pointerleave', hide);svg.addEventListener('focus', () => show(current < 0 ? 0 : current));svg.addEventListener('blur', hide);svg.addEventListener('keydown', (event) => { if (event.key === 'ArrowRight') { show(current + 1); event.preventDefault(); } if (event.key === 'ArrowLeft') { show(current - 1); event.preventDefault(); } if (event.key === 'Escape') hide();});})();</script>
</figure>

| $m$ | symmetric, $\rho = 0.924$ | Jordan, $\rho = 0.8$ | triangular, $\rho = 0.6$ | symmetric, $\rho = 1.05$ |
|---|---|---|---|---|
| 1 | 0.924 | 1.69 | 4.07 | 1.05 |
| 5 | 0.674 | 10.1 | 1.86 | 1.28 |
| 14 | 0.332 | **32.5** (peak) | $2.9 \times 10^{-2}$ | 1.98 |
| 20 | 0.207 | 26.1 | $1.4 \times 10^{-3}$ | 2.65 |
| 40 | $4.3 \times 10^{-2}$ | 2.58 | $5.4 \times 10^{-8}$ | 7.04 |
| 80 | $1.8 \times 10^{-3}$ | $2.8 \times 10^{-3}$ | $7.2 \times 10^{-17}$ | 49.6 |
| 120 | $7.9 \times 10^{-5}$ | $1.3 \times 10^{-6}$ | $9.6 \times 10^{-26}$ | 349 |

What the numbers say:

- **Symmetric, $\rho < 1$:** $\|A\|_2 = \rho(A) = 0.924$ exactly, as the claim in Section 7 says, and $\|A^m\|_2 = 0.924^m$. The decay is smooth from the start.
- **Jordan block:** $\rho = 0.8 < 1$, but the norm grows to $32.5$ at $m = 14$ before it decays. This is the $\binom{m}{j}\lambda^{m-j}$ term from Section 8: for $j = 3$ it peaks near $m = 14$. In the end the geometric factor wins.
- **Triangular:** the big off-diagonal entry gives $\|A\|_2 = 4.07$ and a small bump ($4.42$ at $m = 2$), then fast decay at the rate $0.6^m$.
- **$\rho = 1.05$:** the norm grows like $1.05^m$; $1.05^{120} \approx 349$.

So only $\rho(A)$ decides the limit, but the norm decides what happens on the way. For an iteration, a large hump means many wasted steps (and lost digits) before the error starts to shrink.

### B. The Neumann series

For $C = \begin{bmatrix} 0.2 & 0.1 & 0.3 \\ 0.1 & 0.4 & 0.1 \\ 0 & 0.2 & 0.3 \end{bmatrix}$, $\|C\|_2 = 0.598$. The bound gives $\|(I - C)^{-1}\|_2 \le 1/(1 - 0.598) = 2.49$, and the true value is $2.37$. I compared the partial sums $S_m = I + C + \cdots + C^m$ with $(I - C)^{-1}$. The tail of the series is bounded by $\|C\|^{m+1}/(1 - \|C\|)$:

| $m$ | $\lVert (I - C)^{-1} - S_m \rVert_2$ | bound $\lVert C\rVert^{m+1}/(1 - \lVert C\rVert)$ |
|---|---|---|
| 1 | $7.96 \times 10^{-1}$ | $8.90 \times 10^{-1}$ |
| 2 | $4.51 \times 10^{-1}$ | $5.32 \times 10^{-1}$ |
| 5 | $8.11 \times 10^{-2}$ | $1.14 \times 10^{-1}$ |
| 10 | $4.67 \times 10^{-3}$ | $8.71 \times 10^{-3}$ |
| 20 | $1.55 \times 10^{-5}$ | $5.10 \times 10^{-5}$ |
| 40 | $1.70 \times 10^{-10}$ | $1.75 \times 10^{-9}$ |

The error is always below the bound. Then I took $C = \begin{bmatrix} 0.4 & 3 \\ 0 & 0.5 \end{bmatrix}$. Here $\|C\|_2 = 3.07 > 1$, so the second theorem says nothing. But $\rho(C) = 0.5 < 1$, so the first theorem still applies. Indeed the series converges to $(I - C)^{-1} = \begin{bmatrix} 5/3 & 10 \\ 0 & 2 \end{bmatrix}$: the error is $7.0$ at $m = 1$, $2.7 \times 10^{-2}$ at $m = 10$, and $2.7 \times 10^{-11}$ at $m = 40$.

### C. The cost of elimination

I counted every multiplication, division and subtraction in plain elimination (no pivoting, on random diagonally dominant matrices) and in one forward + back substitution:

| $n$ | LU operations | ÷ $\frac{2}{3}n^3$ | one solve | ÷ $2n^2$ | 1000 right sides: factor once | 1000 right sides: redo elimination |
|---|---|---|---|---|---|---|
| 10 | 615 | 0.923 | 190 | 0.950 | $1.9 \times 10^5$ | $8.1 \times 10^5$ |
| 20 | 5,130 | 0.962 | 780 | 0.975 | $7.9 \times 10^5$ | $5.9 \times 10^6$ |
| 50 | 82,075 | 0.985 | 4,950 | 0.990 | $5.0 \times 10^6$ | $8.7 \times 10^7$ |
| 100 | 661,650 | 0.992 | 19,900 | 0.995 | $2.1 \times 10^7$ | $6.8 \times 10^8$ |
| 200 | 5,313,300 | 0.996 | 79,800 | 0.998 | $8.5 \times 10^7$ | $5.4 \times 10^9$ |

The counts approach $\frac{2}{3}n^3$ and $2n^2$, as predicted. Doubling $n$ from 100 to 200 multiplies the LU work by $8.03$. With a thousand right-hand sides and $n = 200$, factoring once is $63$ times cheaper than eliminating again for every $\mathbf{b}$. The solutions were correct to about $10^{-15}$.

### D. Pivoting in short arithmetic

To make rounding visible, I simulated a computer with **4 significant digits**: every result is rounded to 4 digits. The system has a small first pivot:

$$
\begin{bmatrix} 0.003 & 1.2 & 2.5 \\ 1.6 & -0.8 & 1.1 \\ 2.1 & 1.7 & -0.9 \end{bmatrix}\mathbf{x} = \begin{bmatrix} 9.903 \\ 3.3 \\ 2.8 \end{bmatrix}, \qquad \mathbf{x} = (1, 2, 3).
$$

| | $x_1$ | $x_2$ | $x_3$ | max error |
|---|---|---|---|---|
| exact | 1 | 2 | 3 | |
| no pivoting | 0.6667 | 2.463 | 2.778 | 0.46 |
| row pivoting | 0.9995 | 2.001 | 3.001 | 0.001 |

Without pivoting, the first multipliers are $1.6/0.003 \approx 533$ and $2.1/0.003 = 700$. The second and third rows are swamped by 700 times the first row, and their own information is lost in the rounding. With row pivoting, the code swapped rows 1 and 3 first, all multipliers were at most 1, and the error dropped by a factor of about 460.

The same failure happens in normal double precision when the pivot is tiny enough. For $\begin{bmatrix} 10^{-17} & 1 \\ 1 & 2 \end{bmatrix}\mathbf{x} = \begin{bmatrix} 1 \\ 3 \end{bmatrix}$, the solution is very close to $(1, 1)$. Without pivoting, the multiplier is $10^{17}$, and elimination returns $(0, 1)$: $x_1$ is completely wrong. With a row swap it returns $(1, 1)$.

### E. A small residual does not mean a small error

Take a matrix that is close to singular:

$$
A = \begin{bmatrix} 2 & 3 \\ 4 & 6.0001 \end{bmatrix}, \qquad \mathbf{x} = (1, 1), \qquad \mathbf{b} = (5,\ 10.0001).
$$

Its eigenvalues are $8.00$ and $2.5 \times 10^{-5}$, and $\det A = 2 \times 10^{-4}$. The guess $\mathbf{x}_A = (2.5, 0)$ is far from $\mathbf{x}$, but

$$
\mathbf{r} = \mathbf{b} - A\mathbf{x}_A = (0,\ 10^{-4}), \qquad \mathbf{e} = \mathbf{x} - \mathbf{x}_A = (-1.5,\ 1).
$$

The relative residual is $\|\mathbf{r}\|/\|\mathbf{b}\| = 8.9 \times 10^{-6}$, but the relative error is $\|\mathbf{e}\|/\|\mathbf{x}\| = 1.27$. That is, the error is 127% while the residual looks like about 5 correct digits. I also checked $A\mathbf{e} = \mathbf{r}$ numerically. Next week's condition number will explain the factor of about $10^5$ between the two.

Here is the core of the elimination and the substitution, in Python:

```python
def lu(A):
    n = len(A)
    U = [row[:] for row in A]
    L = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    for k in range(n - 1):
        for i in range(k + 1, n):
            m = U[i][k] / U[k][k]          # multiplier, stored in L
            L[i][k] = m
            U[i][k] = 0.0
            for j in range(k + 1, n):
                U[i][j] -= m * U[k][j]
    return L, U

def solve(L, U, b):
    n = len(b)
    y = [0.0] * n
    for i in range(n):                     # forward: L y = b
        y[i] = b[i] - sum(L[i][j] * y[j] for j in range(i))
    x = [0.0] * n
    for i in reversed(range(n)):           # backward: U x = y
        x[i] = (y[i] - sum(U[i][j] * x[j] for j in range(i + 1, n))) / U[i][i]
    return x
```

## Summary

| topic | main idea |
|---|---|
| Why $A\mathbf{x} = \mathbf{b}$ | Discretized PDEs and Newton steps all lead to linear systems, often huge. |
| Basis | Spans the space and is linearly independent, so coefficients are unique. Not unique itself. |
| Symmetric matrices | Real eigenvalues (skew-symmetric: purely imaginary). |
| Invertible | Unique solution $\iff$ only $\mathbf{x} = \mathbf{0}$ solves $A\mathbf{x} = \mathbf{0}$ $\iff$ $A^{-1}$ exists $\iff$ $\det A \ne 0$ $\iff$ rank $n$. |
| Similar matrices | $B = P^{-1}AP$ has the same eigenvalues. Idea behind preconditioning. |
| Schur | Every matrix is unitarily similar to a triangular one, so $\det A = \lambda_1 \cdots \lambda_n$. |
| Principal axes | A Hermitian matrix has an orthonormal basis of eigenvectors: $U^*AU = D$. |
| SVD | $V^*AU = F$, singular values $\mu_1 \ge \cdots \ge \mu_r > 0$, $r$ = rank. |
| Jordan form | Blocks with $\lambda$ on the diagonal and 1 above it. |
| Norms | All norms are equivalent in finite dimensions only. Operator norm = largest stretch. $\|A\|_2 = \sqrt{\rho(A^*A)}$. |
| Spectral radius | $\rho(A) \le \|A\|$ always; $= \|A\|_2$ for symmetric $A$. |
| Powers | $A^m \to 0 \iff \rho(A) < 1$, but $\|A^m\|$ can grow first. |
| Neumann series | $(I - A)^{-1} = \sum A^k$ if $\rho(A) < 1$. $\|(I - A)^{-1}\| \le \frac{1}{1 - \|A\|}$ if $\|A\| < 1$. |
| Gaussian elimination | $= LU$ decomposition. Costs $\frac{2}{3}n^3$; each extra right side costs $2n^2$. |
| Pivoting | Keeps multipliers $\le 1$ in size, so rounding errors do not blow up. |
| Cholesky | $A = CC^t$ for symmetric positive definite $A$. |
| Error vs. residual | $A\mathbf{e} = \mathbf{r}$. The residual is computable, but it can be tiny while the error is large. |

The main thing I take from this week: the eigenvalues tell us a lot, but not everything. The spectral radius decides whether $A^m$ goes to zero, but the hump in Experiment A comes from a matrix that is far from symmetric. The determinant tells us whether $A$ is invertible, but Experiment E shows that "invertible" is not the same as "safe to trust the residual". Next week the condition number will put a number on this.

The full script, with all the numbers and the chart, is in [this site's repository](https://github.com/doubleyands/doubleyands.github.io/tree/main/scripts/posts/numerical-analysis/week-7).
