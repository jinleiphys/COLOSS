# Lagrange-Legendre Method for Scattering (Method 3)

## Overview

Method 3 solves the scattering problem using Lagrange-Legendre basis functions on a finite interval $[0, R]$. Unlike Lagrange-Laguerre basis (which has built-in exponential decay), Lagrange-Legendre is a pure polynomial basis that can properly represent oscillating scattering waves.

**Reference:** D. Baye, Physics Reports 565 (2015) 1-107

---

## 1. Shifted Legendre Mesh

### 1.1 Mesh Points

The mesh points $x_j$ ($j = 1, \ldots, N$) are zeros of the shifted Legendre polynomial on $(0,1)$:

$$P_N(2x_j - 1) = 0$$

These are obtained by mapping Gauss-Legendre points $t_j \in (-1, 1)$ to $(0, 1)$:

$$x_j = \frac{t_j + 1}{2}$$

Physical coordinates: $r_j = R \cdot x_j$

### 1.2 Gauss Weights

The Gauss-Legendre weights on $(-1, 1)$:

$$w_j^{GL} = \frac{2}{(1 - t_j^2)[P_N'(t_j)]^2}$$

Weights on $(0, 1)$:

$$\lambda_j = \frac{w_j^{GL}}{2}$$

---

## 2. Standard Lagrange Basis Functions

The standard (unregularized) Lagrange interpolation functions:

$$L_j(x) = \prod_{k \neq j} \frac{x - x_k}{x_j - x_k}$$

**Key property:** $L_j(x_i) = \delta_{ij}$

---

## 3. x-Regularized Lagrange Basis

For scattering problems, we use basis functions regularized at the origin ($x = 0$) to handle the $1/r^2$ centrifugal singularity:

$$\hat{f}_j(x) = (-1)^{N-j} \sqrt{\frac{1-x_j}{x_j}} \frac{x \, P_N(2x-1)}{x - x_j}$$

This can be written as:

$$\hat{f}_j(x) = \alpha_j \cdot x \cdot L_j(x)$$

where $\alpha_j$ is a normalization constant.

**Key properties:**
- $\hat{f}_j(0) = 0$ (regularization at origin)
- $\hat{f}_j(x_i) = \alpha_j x_j \delta_{ij}$
- No constraint at $x = 1$ (allows oscillating waves at boundary)

---

## 4. Derivative Matrices

### 4.1 First Derivative Matrix D1 (Standard Basis)

**Off-diagonal ($i \neq j$):**

$$D1_{ij} = \frac{(-1)^{i-j}}{x_i - x_j} \sqrt{\frac{\lambda_j}{\lambda_i}}$$

Or equivalently using barycentric weights $w_i = 1/\prod_{k \neq i}(x_i - x_k)$:

$$D1_{ij} = \frac{w_j}{w_i(x_i - x_j)}$$

**Diagonal:**

$$D1_{ii} = \frac{1 - 2x_i}{2x_i(1-x_i)}$$

Or from row-sum property (derivative of constant = 0):

$$D1_{ii} = -\sum_{k \neq i} D1_{ik}$$

### 4.2 Second Derivative Matrix D2 (Standard Basis)

Computed as matrix product:

$$D2 = D1 \times D1$$

**Analytical off-diagonal ($i \neq j$):**

$$D2_{ij} = \frac{(-1)^{i-j}}{(x_i-x_j)^2} \sqrt{\frac{\lambda_j}{\lambda_i}} \left[\frac{1-2x_i}{x_i(1-x_i)} - \frac{2}{x_i-x_j}\right]$$

**Analytical diagonal:**

$$D2_{ii} = \frac{N(N+1)}{3x_i(1-x_i)} + \frac{(1-2x_i)^2 - 1}{4x_i^2(1-x_i)^2}$$

### 4.3 Effective D2 for x-Regularized Basis

For the wave function expanded in x-regularized basis:

$$\phi(x) = \sum_j c_j \hat{f}_j(x)$$

The second derivative becomes:

$$\phi''(x_i) = \sum_j \tilde{D2}_{ij} \, \phi(x_j)$$

where the **effective second derivative matrix** is:

$$\boxed{\tilde{D2}_{ij} = \frac{2 D1_{ij}}{x_j} + \frac{x_i}{x_j} D2_{ij}}$$

**Derivation:**

From $\hat{f}_j(x) = \alpha_j x L_j(x)$:

$$\hat{f}_j''(x) = \alpha_j [2 L_j'(x) + x L_j''(x)]$$

At mesh points with $c_j = \phi(x_j)/(\alpha_j x_j)$:

$$\phi''(x_i) = \sum_j \frac{\phi(x_j)}{x_j} [2 D1_{ij} + x_i D2_{ij}]$$

### 4.4 Scaling to Physical Coordinates

For physical coordinates $r = R \cdot x$:

$$D1^{(r)} = \frac{1}{R} D1^{(x)}, \qquad D2^{(r)} = \frac{1}{R^2} D2^{(x)}$$

The effective D2 in physical coordinates:

$$\tilde{D2}_{ij}^{(r)} = \frac{2 D1_{ij}^{(r)}}{r_j} + \frac{r_i}{r_j} D2_{ij}^{(r)}$$

---

## 5. Scattering Equation

### 5.1 Physical Decomposition

The scattering wave function is decomposed as:

$$\psi_l = F_l(kr) + \phi_l$$

where:
- $F_l(kr)$ = regular Coulomb function (known, oscillating)
- $\phi_l$ = scattered wave (to be solved)

### 5.2 Equation for Scattered Wave

$$\left[\frac{d^2}{dr^2} - \frac{l(l+1)}{r^2} - U(r) + k^2\right] \phi = U(r) \cdot F_l(kr)$$

where $U(r) = 2\mu V_{\text{short}}(r) / \hbar^2$ is the reduced short-range potential.

### 5.3 Matrix Form

Discretizing on mesh points $r_1, \ldots, r_N$:

$$\sum_j M_{ij} \, \phi(r_j) = b_i$$

**Matrix elements (rows 1 to N-1):**

$$M_{ij} = D2_{ij} + \left[k^2 - \frac{l(l+1)}{r_i^2} - U(r_i)\right] \delta_{ij}$$

**Source vector:**

$$b_i = U(r_i) \cdot F_l(k r_i)$$

---

## 6. Boundary Conditions

### 6.1 Origin ($r = 0$)

Automatically satisfied by x-regularization: $\phi(0) = 0$

### 6.2 Outer Boundary ($r = R$): Outgoing Wave

$$\phi'(R) = \gamma_s \cdot \phi(R)$$

where $\gamma_s$ is the logarithmic derivative of the outgoing Hankel function:

$$\gamma_s = \frac{H_l^{+\prime}(kR)}{H_l^+(kR)}$$

with $H_l^+ = G_l + i F_l$ (outgoing Coulomb-Hankel function).

### 6.3 Implementation

**Row N of matrix:** The boundary condition replaces the differential equation:

$$\sum_j \left[D1_{\text{at }R}(j) - \gamma_s \cdot L_j(R)\right] \phi(r_j) = 0$$

where:
- $D1_{\text{at }R}(j) = L_j'(R)$ = derivative of Lagrange function at $r = R$
- $L_j(R) = \prod_{k \neq j} \frac{R - r_k}{r_j - r_k}$ = Lagrange function at $r = R$

---

## 7. Scattering Amplitude Extraction

After solving the linear system for $\phi(r_j)$:

### 7.1 Wave Function at Boundary

$$\phi(R) = \sum_j L_j(R) \cdot \phi(r_j)$$

### 7.2 Scattering Amplitude

From the asymptotic form $\phi \sim f_l k \cdot H_l^+(kr)$:

$$f_l = \frac{\phi(R)}{kH_l^+(kR)}$$

### 7.3 S-Matrix

$$S_l = 1 + 2i kf_l$$

### 7.4 Reaction Cross Section

$$\sigma_{\text{reac}} = \frac{\pi}{k^2} \frac{2J+1}{2S+1} \left(1 - |S_l|^2\right)$$

---

## 8. Method 3 vs Method 4

| Aspect | Method 3 | Method 4 (Corrected) |
|--------|----------|---------------------|
| D1 | Barycentric formula | Same |
| D2 | $D1 \times D1$ | $\tilde{D2}_{ij} = \frac{2D1_{ij}}{r_j} + \frac{r_i}{r_j}D2_{ij}$ |
| Basis | Standard polynomial | x-regularized |
| Origin BC | Explicit or implicit | Built into basis |

**Note:** The original Method 4 implementation incorrectly used Baye's T-matrix formula which assumes regularization by $x(1-x)$ (with Bloch operator). The correct formula for x-only regularization is given above.

---

## 9. Summary of Key Formulas

### Mesh
$$x_j = \frac{t_j + 1}{2}, \quad r_j = R \cdot x_j$$

### Standard D1
$$D1_{ij} = \frac{w_j}{w_i(x_i - x_j)} \quad (i \neq j), \qquad D1_{ii} = -\sum_{k \neq i} D1_{ik}$$

### Standard D2
$$D2 = D1 \times D1$$

### x-Regularized D2
$$\tilde{D2}_{ij} = \frac{2 D1_{ij}}{r_j} + \frac{r_i}{r_j} D2_{ij}$$

### Outgoing Wave BC
$$\gamma_s = \frac{H_l^{+\prime}(kR)}{H_l^+(kR)}, \qquad H_l^+ = G_l + iF_l$$

### Scattering Amplitude
$$f_l = \frac{\phi(R)}{kH_l^+(kR)}, \qquad S_l = 1 + 2ikf_l$$







有库仑的时候direct matching可能有问题！ 是不是因为有库仑力的时候散射振幅的定义有问题，， $$\psi_l = F_l(kr) + \phi_l$$
