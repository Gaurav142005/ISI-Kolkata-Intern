# Robust Regression Analysis: L4 and LPQ Estimators vs. LS

This repository contains simulation code and supporting material for evaluating and comparing **higher-order norm-based estimators**—specifically **L4 estimators** and Non linear estimators - **LPQ estimators**—with the traditional **Least Squares (LS)** estimator in linear regression settings under conditions of non normality.

## Purpose

The goal of this project is to:
- Explore the performance of **L4** (fourth-order loss) and **LPQ** (power-q norm) regression estimators.
- Compare these with **LS (L2)** estimators across a variety of error distributions.
- Analyze their robustness and efficiency under both **homoskedastic** and **heteroskedastic** settings.
- Evaluate theoretical and empirical efficiency using metrics like MSE and pseudo-R².
- Provide insights into when L4/LPQ methods should be preferred over LS.

## Key Components

### Estimators
- **LS Estimator:** Minimizes squared residuals.
- **LPQ Estimator:** Minimizes the Lᵩ-norm.
- **L4 Estimator:** Minimizes the fourth power of residuals.

### Simulations
- Monte Carlo simulations across 500 iterations.
- Simulated design matrices and controlled error variances.
- Wide range of **error distributions**:
  - Gaussian, Laplace, Uniform, Truncated Normal, Lognormal
  - Pareto (heavy-tailed), Exponential, Gamma, Beta
  - Mixture of Gaussians (GMM)
  - Cauchy (for robustness testing)
  - U-shaped and bell-shaped distributions from a Pearsonian family

### 📈 Evaluation Metrics
- **Mean Squared Error (MSE)** of estimated coefficients
- **Efron's pseudo-R²** for in-sample fit
- **Execution Time** for each estimator
- **Skewness & Kurtosis** of error distributions
- **L4 Efficiency Condition** from theory
