data <- scan("fortdataset.36")

# Log-likelihood (positive, for maximization)
loglik <- function(theta, x){
  alpha  <- theta[1]; lambda <- theta[2]
  if(alpha <= 0 || lambda <= 0) return(-Inf)
  n <- length(x)
  exp_term  <- exp(-lambda * x)
  one_minus <- pmax(1 - exp_term, 1e-10)
  n*log(alpha) + n*log(lambda) - lambda*sum(x) + (alpha-1)*sum(log(one_minus))
}
NR_GE <- function(x, init = c(2, 0.8), tol = 1e-8, max_iter = 500) {
  
  theta <- init
  n     <- length(x)
  iter  <- 0
  grad  <- rep(NA, 2)
  
  for (i in 1:max_iter) {
    
    alpha  <- theta[1]
    lambda <- theta[2]
    
    # ---------- Precompute ----------
    exp_term  <- exp(-lambda * x)
    one_minus <- pmax(1 - exp_term, 1e-10)
    
    A <- log(one_minus)
    B <- (x * exp_term) / one_minus          # x·e^{-λx}/(1-e^{-λx})
    D <- (x^2 * exp_term) / one_minus^2      # x²·e^{-λx}/(1-e^{-λx})²
    
    # ---------- Gradient ----------
    grad <- c(
      n / alpha + sum(A),
      n / lambda - sum(x) + (alpha - 1) * sum(B)
    )
    
    # ---------- Hessian (corrected H22) ----------
    H11 <- -n / alpha^2
    H12 <- sum(B)
    H22 <- -n / lambda^2 - (alpha - 1) * sum(D)   # <-- fixed sign & formula
    
    H <- matrix(c(H11, H12, H12, H22), 2, 2)
    H <- H - diag(1e-6, 2)   # regularize: subtract (H is negative-definite at max)
    
    # ---------- Newton Step ----------
    step <- tryCatch(solve(H, grad), error = function(e) NULL)
    
    if (is.null(step) || any(is.nan(step)) || any(is.infinite(step))) {
      step <- grad / sqrt(sum(grad^2))   # gradient ascent fallback
    }
    
    # ---------- Backtracking Line Search ----------
    step_size <- 1
    ll_curr   <- loglik(theta, x)
    theta_new <- theta
    
    for (j in 1:50) {
      candidate <- theta - step_size * step   # subtract: Newton on negative H
      
      if (all(candidate > 0)) {
        ll_new <- loglik(candidate, x)
        if (ll_new > ll_curr) {
          theta_new <- candidate
          break
        }
      }
      step_size <- step_size / 2
    }
    
    iter <- iter + 1
    
    if (max(abs(theta_new - theta)) < tol) break
    
    theta <- theta_new
  }
  
  # Recompute gradient at final theta for convergence check
  alpha  <- theta[1]; lambda <- theta[2]
  exp_term  <- exp(-lambda * x)
  one_minus <- pmax(1 - exp_term, 1e-10)
  grad_final <- c(
    n / alpha + sum(log(one_minus)),
    n / lambda - sum(x) + (alpha - 1) * sum((x * exp_term) / one_minus)
  )
  
  return(list(
    alpha     = theta[1],
    lambda    = theta[2],
    loglik    = loglik(theta, x),
    Iter      = iter,
    step_size = step_size,
    converged = (max(abs(grad_final)) < 1e-4)
  ))
}

NR_GE(data)

# -------- Initial value selection --------
xbar <- mean(data)

init_list <- list(
  c(1, 1/xbar),
  c(2, 1/xbar),
  c(3, 1/xbar),
  c(2, 0.5),
  c(3, 1),
  c(4, 1.5)
)

results <- lapply(init_list, function(init){
  res <- NR_GE(data, init=init)
  list(
    alpha = res$alpha,
    lambda = res$lambda,
    loglik = res$loglik,
    ITER = res$Iter,
    CONVERGENCE = (res$converged == 1)
  )
})
results_mat <- do.call(rbind, results)
results_mat
best <- results_mat[which.max(results_mat[,"loglik"]), ]
best


#Parametric Boot Strap
set.seed(123)

alpha_hat  <- 3.88
lambda_hat <- 1.60

# Generator (inverse CDF)
rgen <- function(n, alpha, lambda){
  u <- runif(n)
  -log(1 - u^(1/alpha))/lambda
}

B <- 500
boot_par <- matrix(0, B, 2)

for(i in 1:B){
  x_star <- rgen(length(data), alpha_hat, lambda_hat)
  
  res <- NR_GE(x_star, init=c(alpha_hat, lambda_hat))
  
  boot_par[i,] <- c(res$alpha, res$lambda)
}

# 95% CI
ci_par <- apply(boot_par, 2, quantile, c(0.025, 0.975))
ci_par



#Non paramteric Bootstrap
set.seed(123)

boot_np <- matrix(0, B, 2)

for(i in 1:B){
  x_star <- sample(data, replace=TRUE)
  
  res <- NR_GE(x_star, init=c(alpha_hat, lambda_hat))
  
  boot_np[i,] <- c(res$alpha, res$lambda)
}

ci_np <- apply(boot_np, 2, quantile, c(0.025, 0.975))
ci_np


