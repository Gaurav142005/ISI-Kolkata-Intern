**Optimal LPQ Estimator Without Moore–Penrose Inverse**

---

### 1. Model Setup

We consider the linear regression model:
$y = X\beta + \varepsilon$
where:

* $y \in \mathbb{R}^n$, $X \in \mathbb{R}^{n \times p}$,
* $\mathbb{E}[\varepsilon] = 0$,
* ${Var}(\varepsilon) = \sigma^2 I_n$.

The ordinary least squares (OLS) estimator is:
$b = (X^\top X)^{-1} X^\top y \quad \text{(assuming full column rank)}$

Let $M = I - X(X^\top X)^{-1} X^\top$ denote the residual-maker matrix.

---

### 2. LPQ Estimator Definition

We define the LPQ estimator for the $i$th coefficient as:
$\tilde{\beta}_i = b_i^\top y + \frac{1}{\sigma} y^\top H_i y - \sigma \operatorname{tr}(H_i)$
where:

* $b_i^\top = e_i^\top (X^\top X)^{-1} X^\top$
* $H_i \in \mathbb{R}^{n \times n}$ is symmetric and satisfies $H_i X = 0 \Rightarrow H_i = M H_i M$

---

### 3. Unbiasedness Proof

We compute the expectation:

$$
\mathbb{E}[\tilde{\beta}_i] = \mathbb{E}[b_i^\top y] + \frac{1}{\sigma} \mathbb{E}[y^\top H_i y] - \sigma \operatorname{tr}(H_i)
$$

Since $\mathbb{E}[y] = X\beta$, and $H_i X = 0$:

* $\mathbb{E}[b_i^\top y] = b_i^\top X\beta = \beta_i$
* $\mathbb{E}[y^\top H_i y] = \mathbb{E}[\varepsilon^\top H_i \varepsilon] = \sigma^2 \operatorname{tr}(H_i)$

Thus:
$\mathbb{E}[\tilde{\beta}_i] = \beta_i + \sigma \operatorname{tr}(H_i) - \sigma \operatorname{tr}(H_i) = \beta_i$
$\Rightarrow \tilde{\beta}_i$ is unbiased.

---

### 4. Variance of the LPQ Estimator

Let:

* $\mu = (m_{31}, \dots, m_{3n})^\top$: skewness vector
* $\gamma = (m_{41} - 3, \dots, m_{4n} - 3)^\top$: excess kurtosis vector

Then:
$\operatorname{Var}(\tilde{\beta}_i) = \sigma^2 \left[ b_i^\top b_i + 2\operatorname{tr}(H_i^2) + 3 \sum_j m_{3j} b_{ij} H_{i,jj} + \sum_j \gamma_j H_{i,jj}^2 \right]$

---

### 5. Optimization: Finding the Optimal $H_i$

Let $h_i = \operatorname{diag}(H_i) \in \mathbb{R}^n$.

Define:

$$
A = 2(M \circ M) + \operatorname{diag}(\gamma), \quad c = 3 (b_i \circ \mu)
$$

Then the MSE minimization becomes:
$\min_{h_i} \quad h_i^\top A h_i + c^\top h_i \quad \Rightarrow \boxed{ h_i = -\frac{1}{2} A^{-1} c }$

Construct the optimal:
$H_i = M \cdot \operatorname{diag}(h_i) \cdot M$

---

### 6. Final Optimal LPQ Estimator

$$
\boxed{ \tilde{\beta}_i = b_i^\top y + \frac{1}{\sigma} y^\top M \operatorname{diag}(h_i) M y - \sigma \cdot \operatorname{tr}(M \operatorname{diag}(h_i) M) }
$$

---

### 7. Comparison with Moore–Penrose Approach

* When $X$ is full column rank:

  * $X^+ = (X^\top X)^{-1} X^\top$, so both approaches give **identical results**.
* When $X$ is rank-deficient:

  * Standard inverse fails
  * Moore–Penrose $X^+$ is well-defined and yields a valid LPQ estimator

$\Rightarrow \text{Moore–Penrose estimator is more general and stable}$
