# Maximum Likelihood Estimation and Bootstrap Inference for the Generalized Exponential Distribution

## Overview

This project implements Maximum Likelihood Estimation (MLE) for the Generalized Exponential (GE) distribution using the Newton-Raphson optimization algorithm in R. The project also constructs 95% confidence intervals for the model parameters using both Parametric and Non-Parametric Bootstrap methods.

The work was completed as part of a Statistical Computing course project under **Prof. Debasis Kundu, IIT Kanpur**.

---

## Distribution

The Generalized Exponential distribution is defined by the probability density function

[
f(x;\alpha,\lambda)=\alpha\lambda e^{-\lambda x}(1-e^{-\lambda x})^{\alpha-1},
\quad x>0,; \alpha>0,; \lambda>0
]

where

* (\alpha) = shape parameter
* (\lambda) = rate parameter

---

## Objectives

* Estimate the parameters (\alpha) and (\lambda) using Maximum Likelihood Estimation.
* Implement the Newton-Raphson algorithm from scratch.
* Study convergence behavior under multiple initialization strategies.
* Construct 95% confidence intervals using:

  * Parametric Bootstrap
  * Non-Parametric Bootstrap
* Compare uncertainty estimates obtained from both methods.

---

## Methodology

### 1. Maximum Likelihood Estimation

The log-likelihood function was derived analytically and optimized numerically.

The implementation includes:

* Analytical score vector
* Analytical Hessian matrix
* Newton-Raphson updates
* Backtracking line search
* Hessian regularization
* Positivity constraints on parameters
* Gradient-ascent fallback for numerical stability

---

### 2. Bootstrap Confidence Intervals

#### Parametric Bootstrap

* Generate synthetic datasets from the fitted GE distribution.
* Re-estimate parameters for each sample.
* Obtain percentile-based confidence intervals.

#### Non-Parametric Bootstrap

* Resample observations with replacement from the original dataset.
* Re-estimate parameters.
* Construct percentile-based confidence intervals.

---

## Dataset

Dataset file:

```text
fortdataset.36
```

Sample size:

```text
n = 25
```

---

## Initial Values Used

The Newton-Raphson algorithm was tested using six different starting points:

| α₀ | λ₀   |
| -- | ---- |
| 1  | 1/x̄ |
| 2  | 1/x̄ |
| 3  | 1/x̄ |
| 2  | 0.5  |
| 3  | 1.0  |
| 4  | 1.5  |

All initializations converged to the same solution.

---

## Results

### Maximum Likelihood Estimates

| Parameter | Estimate |
| --------- | -------- |
| α̂        | 3.8798   |
| λ̂        | 1.6035   |

Log-likelihood:

```text
-24.7141
```

Convergence achieved in:

```text
5–9 iterations
```

---

## Bootstrap Confidence Intervals

### Parametric Bootstrap (B = 500)

| Parameter | Lower 95% CI | Upper 95% CI |
| --------- | ------------ | ------------ |
| α         | 2.1927       | 12.8089      |
| λ         | 1.1346       | 2.5488       |

### Non-Parametric Bootstrap (B = 500)

| Parameter | Lower 95% CI | Upper 95% CI |
| --------- | ------------ | ------------ |
| α         | 2.7206       | 8.2768       |
| λ         | 1.1956       | 2.5496       |

---

## Key Observations

* The Newton-Raphson algorithm converged rapidly from all initialization points.
* Shape parameter estimates exhibit noticeable positive skewness.
* Non-parametric confidence intervals for α are narrower than parametric intervals.
* Confidence intervals for λ are highly consistent across both bootstrap methods.

---

## Technologies Used

* R
* Numerical Optimization
* Newton-Raphson Algorithm
* Bootstrap Resampling
* Statistical Inference

---

## Repository Structure

```text
├── data/
│   └── fortdataset.36
├── src/
│   └── mle_bootstrap.R
├── report/
│   └── Statistical_Report.pdf
└── README.md
```

---

## Author

**Yagnesh Bonnada**

B.Tech Undergraduate
Indian Institute of Technology Kanpur

---
