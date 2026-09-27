# ============================================================
# FUNCTIONAL TRANSFORMATION REGRESSION MODEL
# PART 1: PACKAGES, SETTINGS, AND BASIC FUNCTIONS
# ============================================================

rm(list = ls())

set.seed(20260925)

# ------------------------------------------------------------
# Packages
# ------------------------------------------------------------

required_packages <- c(
  "ggplot2",
  "moments"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(ggplot2)
library(moments)

# ------------------------------------------------------------
# Simulation settings
# ------------------------------------------------------------

R <- 1000

sample_sizes <- c(
  100,
  200,
  400,
  800
)

alpha <- 0.05

B_perm <- 500

grid_x <- seq(
  0.05,
  0.95,
  by = 0.05
)

mu2_epanechnikov <- 1 / 5

# ------------------------------------------------------------
# Epanechnikov kernel
# ------------------------------------------------------------

epanechnikov_kernel <- function(u) {
  
  result <- 0.75 * (1 - u^2)
  
  result[abs(u) > 1] <- 0
  
  result
}

# ------------------------------------------------------------
# Default bandwidth
# ------------------------------------------------------------

default_bandwidth <- function(x) {
  
  n <- length(x)
  
  s <- sd(x, na.rm = TRUE)
  
  iqr_value <- IQR(
    x,
    na.rm = TRUE
  )
  
  scale_value <- min(
    s,
    iqr_value / 1.34
  )
  
  if (
    !is.finite(scale_value) ||
    scale_value <= 0
  ) {
    scale_value <- 1
  }
  
  bandwidth <- 0.9 *
    scale_value *
    n^(-1 / 5)
  
  bandwidth
}

# ------------------------------------------------------------
# Univariate Nadaraya-Watson estimator
# ------------------------------------------------------------

nw_estimate <- function(
    x,
    y,
    x_eval,
    bandwidth) {
  
  n_eval <- length(x_eval)
  
  result <- numeric(n_eval)
  
  valid_y <- is.finite(y)
  
  if (sum(valid_y) == 0) {
    return(rep(NA_real_, n_eval))
  }
  
  for (k in seq_len(n_eval)) {
    
    u <- (
      x - x_eval[k]
    ) / bandwidth
    
    weights <- epanechnikov_kernel(u)
    
    weights[!valid_y] <- 0
    
    denominator <- sum(
      weights,
      na.rm = TRUE
    )
    
    if (
      !is.finite(denominator) ||
      denominator <= 1e-12
    ) {
      
      result[k] <- mean(
        y[valid_y],
        na.rm = TRUE
      )
      
    } else {
      
      result[k] <-
        sum(
          weights * y,
          na.rm = TRUE
        ) /
        denominator
    }
  }
  
  result
}

# ------------------------------------------------------------
# Multivariate product-kernel estimator
# ------------------------------------------------------------

multivariate_kernel_estimate <- function(
    Z,
    y,
    Z_eval,
    bandwidth) {
  
  Z <- as.matrix(Z)
  
  Z_eval <- as.matrix(Z_eval)
  
  n <- nrow(Z)
  
  m <- nrow(Z_eval)
  
  p <- ncol(Z)
  
  if (length(bandwidth) != p) {
    stop(
      "Length of bandwidth must equal ",
      "the number of columns of Z."
    )
  }
  
  result <- numeric(m)
  
  valid_y <- is.finite(y)
  
  if (sum(valid_y) == 0) {
    return(rep(NA_real_, m))
  }
  
  for (k in seq_len(m)) {
    
    weights <- rep(
      1,
      n
    )
    
    for (j in seq_len(p)) {
      
      u <- (
        Z[, j] - Z_eval[k, j]
      ) / bandwidth[j]
      
      weights <- weights *
        epanechnikov_kernel(u)
    }
    
    weights[!valid_y] <- 0
    
    denominator <- sum(
      weights,
      na.rm = TRUE
    )
    
    if (
      !is.finite(denominator) ||
      denominator <= 1e-12
    ) {
      
      result[k] <- mean(
        y[valid_y],
        na.rm = TRUE
      )
      
    } else {
      
      result[k] <-
        sum(
          weights * y,
          na.rm = TRUE
        ) /
        denominator
    }
  }
  
  result
}

# ------------------------------------------------------------
# PAVA transformation
# ------------------------------------------------------------

pava_transformation <- function(
    Y,
    eta) {
  
  ord <- order(eta)
  
  eta_sorted <- eta[ord]
  
  Y_sorted <- Y[ord]
  
  keep <- is.finite(eta_sorted) &
    is.finite(Y_sorted)
  
  eta_sorted <- eta_sorted[keep]
  
  Y_sorted <- Y_sorted[keep]
  
  if (length(Y_sorted) == 0) {
    return(rep(NA_real_, length(Y)))
  }
  
  fit <- stats::isoreg(
    eta_sorted,
    Y_sorted
  )
  
  g_sorted <- fit$yf
  
  g_original <- rep(
    NA_real_,
    length(Y)
  )
  
  original_indices <- ord[keep]
  
  g_original[original_indices] <-
    g_sorted
  
  finite_values <- is.finite(
    g_original
  )
  
  if (
    any(finite_values) &&
    any(!finite_values)
  ) {
    
    g_original[!finite_values] <-
      mean(
        g_original[finite_values],
        na.rm = TRUE
      )
  }
  
  g_original
}

# ------------------------------------------------------------
# Function-performance measures
# ------------------------------------------------------------

calc_function_metrics <- function(
    estimated,
    truth) {
  
  error <- estimated - truth
  
  c(
    Bias = mean(
      error,
      na.rm = TRUE
    ),
    
    AbsBias = mean(
      abs(error),
      na.rm = TRUE
    ),
    
    RMSE = sqrt(
      mean(
        error^2,
        na.rm = TRUE
      )
    ),
    
    ISE = mean(
      error^2,
      na.rm = TRUE
    )
  )
}

# ------------------------------------------------------------
# Check that Part 1 works
# ------------------------------------------------------------

cat("\nPart 1 completed successfully.\n")

cat(
  "Epanechnikov kernel at 0 = ",
  epanechnikov_kernel(0),
  "\n"
)

cat(
  "Default bandwidth example = ",
  default_bandwidth(
    runif(100)
  ),
  "\n"
)

##################################################################################

# ============================================================
# FUNCTIONAL TRANSFORMATION REGRESSION MODEL
# PART 2: DATA-GENERATING FUNCTIONS
# ============================================================

# ------------------------------------------------------------
# Error generator
# ------------------------------------------------------------

generate_error <- function(
    n,
    distribution = c(
      "normal",
      "t5",
      "gamma"
    )) {
  
  distribution <- match.arg(
    distribution
  )
  
  if (distribution == "normal") {
    
    e <- rnorm(n)
    
  } else if (distribution == "t5") {
    
    e <- rt(
      n,
      df = 5
    )
    
    e <- e / sqrt(
      5 / 3
    )
    
  } else {
    
    e <- rgamma(
      n,
      shape = 4,
      rate = 4
    )
    
    e <- e - mean(e)
    
    e_sd <- sd(e)
    
    if (
      is.finite(e_sd) &&
      e_sd > 0
    ) {
      e <- e / e_sd
    }
  }
  
  e
}

# ------------------------------------------------------------
# True component functions
# ------------------------------------------------------------

f1_true_function <- function(x) {
  
  sin(
    2 * pi * x
  )
}

f2_true_function <- function(x) {
  
  x^2
}

h_true_function <- function(
    u,
    v) {
  
  u +
    v +
    0.5 * u * v
}

# ------------------------------------------------------------
# Baseline data-generating mechanism
# ------------------------------------------------------------

generate_baseline <- function(
    n,
    sigma = 0.5,
    response_transformation = "log",
    error_distribution = "normal") {
  
  # Covariates
  X1 <- runif(
    n,
    min = 0,
    max = 1
  )
  
  X2 <- runif(
    n,
    min = 0,
    max = 1
  )
  
  # True nonparametric functions
  f1 <- f1_true_function(X1)
  
  f2 <- f2_true_function(X2)
  
  # Joint transformation
  eta <- h_true_function(
    f1,
    f2
  )
  
  # Random error
  eps <- sigma *
    generate_error(
      n = n,
      distribution = error_distribution
    )
  
  # Latent response
  gY <- eta + eps
  
  # Observed response
  if (
    response_transformation == "log"
  ) {
    
    Y <- exp(gY)
    
  } else if (
    response_transformation == "identity"
  ) {
    
    Y <- gY
    
  } else {
    
    stop(
      "response_transformation must be ",
      "'log' or 'identity'."
    )
  }
  
  data.frame(
    Y = Y,
    X1 = X1,
    X2 = X2,
    eta = eta,
    eps = eps,
    f1_true = f1,
    f2_true = f2
  )
}

# ------------------------------------------------------------
# Identity-response data-generating mechanism
# ------------------------------------------------------------

generate_identity <- function(
    n,
    sigma = 0.5,
    error_distribution = "normal") {
  
  X1 <- runif(
    n,
    0,
    1
  )
  
  X2 <- runif(
    n,
    0,
    1
  )
  
  f1 <- f1_true_function(X1)
  
  f2 <- f2_true_function(X2)
  
  eta <- h_true_function(
    f1,
    f2
  )
  
  eps <- sigma *
    generate_error(
      n = n,
      distribution = error_distribution
    )
  
  Y <- eta + eps
  
  data.frame(
    Y = Y,
    X1 = X1,
    X2 = X2,
    eta = eta,
    eps = eps,
    f1_true = f1,
    f2_true = f2
  )
}

# ------------------------------------------------------------
# Log-response data-generating mechanism
# ------------------------------------------------------------

generate_log_response <- function(
    n,
    sigma = 0.5,
    error_distribution = "normal") {
  
  X1 <- runif(
    n,
    0,
    1
  )
  
  X2 <- runif(
    n,
    0,
    1
  )
  
  f1 <- f1_true_function(X1)
  
  f2 <- f2_true_function(X2)
  
  eta <- h_true_function(
    f1,
    f2
  )
  
  eps <- sigma *
    generate_error(
      n = n,
      distribution = error_distribution
    )
  
  gY <- eta + eps
  
  Y <- exp(gY)
  
  data.frame(
    Y = Y,
    X1 = X1,
    X2 = X2,
    eta = eta,
    eps = eps,
    f1_true = f1,
    f2_true = f2
  )
}

# ------------------------------------------------------------
# Prediction RMSE
# ------------------------------------------------------------

prediction_rmse <- function(
    y_true,
    y_hat) {
  
  sqrt(
    mean(
      (
        y_true - y_hat
      )^2,
      na.rm = TRUE
    )
  )
}

# ------------------------------------------------------------
# Common evaluation grids
# ------------------------------------------------------------

f1_true_grid <- f1_true_function(
  grid_x
)

f2_true_grid <- f2_true_function(
  grid_x
)

grid_h <- expand.grid(
  u = grid_x,
  v = grid_x
)

h_true_grid <- h_true_function(
  grid_h$u,
  grid_h$v
)

# ------------------------------------------------------------
# Basic checks
# ------------------------------------------------------------

test_data <- generate_baseline(
  n = 100,
  sigma = 0.5
)

cat("\nPart 2 completed successfully.\n")

cat(
  "Number of observations: ",
  nrow(test_data),
  "\n"
)

cat(
  "Number of variables: ",
  ncol(test_data),
  "\n"
)

cat(
  "Variables:\n"
)

print(
  names(test_data)
)

cat(
  "\nFirst five observations:\n"
)

print(
  head(
    test_data,
    5
  )
)

##################################################################################

# ============================================================
# FUNCTIONAL TRANSFORMATION REGRESSION MODEL
# PART 3: MODEL FITTING FUNCTION
# ============================================================

fit_functional_transformation <- function(
    dat,
    max_iter = 30,
    tol = 1e-4,
    bandwidth_multiplier = 1) {
  
  # ----------------------------------------------------------
  # Check input
  # ----------------------------------------------------------
  
  if (!is.data.frame(dat)) {
    stop(
      "'dat' must be a data.frame."
    )
  }
  
  if (!("Y" %in% names(dat))) {
    stop(
      "The data must contain a response variable named 'Y'."
    )
  }
  
  # ----------------------------------------------------------
  # Response
  # ----------------------------------------------------------
  
  Y <- dat$Y
  
  n <- length(Y)
  
  if (n == 0) {
    stop(
      "The response variable Y is empty."
    )
  }
  
  # ----------------------------------------------------------
  # Extract covariates X1, X2, ...
  # ----------------------------------------------------------
  
  x_names <- grep(
    "^X[0-9]+$",
    names(dat),
    value = TRUE
  )
  
  if (length(x_names) == 0) {
    stop(
      "No covariates named X1, X2, ... were found."
    )
  }
  
  X <- as.matrix(
    dat[
      ,
      x_names,
      drop = FALSE
    ]
  )
  
  storage.mode(X) <- "numeric"
  
  p <- ncol(X)
  
  # ----------------------------------------------------------
  # Response transformation
  # ----------------------------------------------------------
  
  if (
    all(
      is.finite(Y) &
      Y > 0
    )
  ) {
    
    eta_response <- log(Y)
    
    response_type <- "log"
    
  } else {
    
    eta_response <- Y
    
    response_type <- "identity"
  }
  
  # ----------------------------------------------------------
  # Initial bandwidths
  # ----------------------------------------------------------
  
  bandwidths <- numeric(p)
  
  for (j in seq_len(p)) {
    
    bandwidths[j] <-
      bandwidth_multiplier *
      default_bandwidth(
        X[, j]
      )
    
    if (
      !is.finite(bandwidths[j]) ||
      bandwidths[j] <= 0
    ) {
      
      bandwidths[j] <- 0.1
    }
  }
  
  # ----------------------------------------------------------
  # Initial estimates of f_j
  # ----------------------------------------------------------
  
  f_hat <- matrix(
    0,
    nrow = n,
    ncol = p
  )
  
  for (j in seq_len(p)) {
    
    f_hat[, j] <-
      nw_estimate(
        x = X[, j],
        y = eta_response,
        x_eval = X[, j],
        bandwidth = bandwidths[j]
      )
    
    f_mean <- mean(
      f_hat[, j],
      na.rm = TRUE
    )
    
    if (is.finite(f_mean)) {
      
      f_hat[, j] <-
        f_hat[, j] -
        f_mean
    }
  }
  
  # ----------------------------------------------------------
  # Initial transformed covariates
  # ----------------------------------------------------------
  
  Z <- f_hat
  
  # ----------------------------------------------------------
  # Initial estimate of h
  # ----------------------------------------------------------
  
  initial_value <- mean(
    eta_response,
    na.rm = TRUE
  )
  
  if (!is.finite(initial_value)) {
    
    initial_value <- 0
  }
  
  h_hat <- rep(
    initial_value,
    n
  )
  
  # ----------------------------------------------------------
  # Iteration controls
  # ----------------------------------------------------------
  
  convergence <- Inf
  
  iteration <- 0
  
  objective_history <- numeric(
    max_iter
  )
  
  # ----------------------------------------------------------
  # Main backfitting iteration
  # ----------------------------------------------------------
  
  while (
    iteration < max_iter &&
    convergence > tol
  ) {
    
    iteration <- iteration + 1
    
    f_old <- f_hat
    
    h_old <- h_hat
    
    # --------------------------------------------------------
    # Bandwidth for multivariate h
    # --------------------------------------------------------
    
    h_bandwidth_value <-
      max(
        mean(bandwidths),
        0.05
      )
    
    h_bandwidth <- rep(
      h_bandwidth_value,
      p
    )
    
    # --------------------------------------------------------
    # Estimate h
    # --------------------------------------------------------
    
    h_hat <- multivariate_kernel_estimate(
      Z = Z,
      y = eta_response,
      Z_eval = Z,
      bandwidth = h_bandwidth
    )
    
    # --------------------------------------------------------
    # Update each f_j
    # --------------------------------------------------------
    
    for (j in seq_len(p)) {
      
      partial_residual <-
        eta_response - h_hat
      
      f_new <- nw_estimate(
        x = X[, j],
        y = partial_residual,
        x_eval = X[, j],
        bandwidth = bandwidths[j]
      )
      
      f_mean <- mean(
        f_new,
        na.rm = TRUE
      )
      
      if (is.finite(f_mean)) {
        
        f_new <- f_new - f_mean
      }
      
      f_hat[, j] <- f_new
    }
    
    # --------------------------------------------------------
    # Update transformed covariates
    # --------------------------------------------------------
    
    Z <- f_hat
    
    # --------------------------------------------------------
    # PAVA transformation
    # --------------------------------------------------------
    
    gY_new <- pava_transformation(
      Y = Y,
      eta = h_hat
    )
    
    g_mean <- mean(
      gY_new,
      na.rm = TRUE
    )
    
    if (is.finite(g_mean)) {
      
      gY_new <- gY_new - g_mean
    }
    
    g_sd <- sd(
      gY_new,
      na.rm = TRUE
    )
    
    if (
      is.finite(g_sd) &&
      g_sd > 0
    ) {
      
      gY_new <- gY_new / g_sd
    }
    
    # --------------------------------------------------------
    # Update response scale
    # --------------------------------------------------------
    
    eta_response <- gY_new
    
    # --------------------------------------------------------
    # Re-estimate h
    # --------------------------------------------------------
    
    h_hat <- multivariate_kernel_estimate(
      Z = Z,
      y = eta_response,
      Z_eval = Z,
      bandwidth = h_bandwidth
    )
    
    # --------------------------------------------------------
    # Convergence criterion
    # --------------------------------------------------------
    
    diff_f <- max(
      abs(
        f_hat - f_old
      ),
      na.rm = TRUE
    )
    
    diff_h <- max(
      abs(
        h_hat - h_old
      ),
      na.rm = TRUE
    )
    
    if (
      !is.finite(diff_f)
    ) {
      diff_f <- Inf
    }
    
    if (
      !is.finite(diff_h)
    ) {
      diff_h <- Inf
    }
    
    convergence <- max(
      diff_f,
      diff_h
    )
    
    # --------------------------------------------------------
    # Objective function
    # --------------------------------------------------------
    
    residual_current <-
      eta_response - h_hat
    
    objective_history[iteration] <-
      mean(
        residual_current^2,
        na.rm = TRUE
      )
  }
  
  # ----------------------------------------------------------
  # Construct estimated f_j functions
  # ----------------------------------------------------------
  
  f_functions <- vector(
    "list",
    p
  )
  
  for (j in seq_len(p)) {
    
    x_train <- X[, j]
    
    y_train <- f_hat[, j]
    
    bw_train <- bandwidths[j]
    
    f_functions[[j]] <- local({
      
      x0 <- x_train
      
      y0 <- y_train
      
      bw0 <- bw_train
      
      function(x_new) {
        
        nw_estimate(
          x = x0,
          y = y0,
          x_eval = x_new,
          bandwidth = bw0
        )
      }
    })
  }
  
  # ----------------------------------------------------------
  # Construct estimated h function
  # ----------------------------------------------------------
  
  h_function <- local({
    
    Z_train <- Z
    
    y_train <- eta_response
    
    bw_h <- rep(
      max(
        mean(bandwidths),
        0.05
      ),
      p
    )
    
    function(
    u,
    v = NULL) {
      
      if (is.null(v)) {
        
        Z_new <- as.matrix(u)
        
      } else {
        
        Z_new <- cbind(
          u,
          v
        )
      }
      
      if (
        ncol(Z_new) !=
        ncol(Z_train)
      ) {
        
        stop(
          "The supplied evaluation points ",
          "do not have the correct dimension."
        )
      }
      
      multivariate_kernel_estimate(
        Z = Z_train,
        y = y_train,
        Z_eval = Z_new,
        bandwidth = bw_h
      )
    }
  })
  
  # ----------------------------------------------------------
  # Construct estimated transformation g
  # ----------------------------------------------------------
  
  g_function <- local({
    
    Y_train <- Y
    
    eta_train <- eta_response
    
    function(y_new) {
      
      ord <- order(
        Y_train
      )
      
      x_values <- Y_train[ord]
      
      y_values <- eta_train[ord]
      
      finite_values <-
        is.finite(x_values) &
        is.finite(y_values)
      
      x_values <- x_values[
        finite_values
      ]
      
      y_values <- y_values[
        finite_values
      ]
      
      if (
        length(x_values) < 2
      ) {
        
        return(
          rep(
            mean(
              y_values,
              na.rm = TRUE
            ),
            length(y_new)
          )
        )
      }
      
      approx(
        x = x_values,
        y = y_values,
        xout = y_new,
        rule = 2
      )$y
    }
  })
  
  # ----------------------------------------------------------
  # Final residuals
  # ----------------------------------------------------------
  
  residuals <- eta_response -
    h_hat
  
  # ----------------------------------------------------------
  # Final output
  # ----------------------------------------------------------
  
  result <- list(
    
    f1 = if (p >= 1) {
      f_functions[[1]]
    } else {
      NULL
    },
    
    f2 = if (p >= 2) {
      f_functions[[2]]
    } else {
      NULL
    },
    
    f3 = if (p >= 3) {
      f_functions[[3]]
    } else {
      NULL
    },
    
    f4 = if (p >= 4) {
      f_functions[[4]]
    } else {
      NULL
    },
    
    h = h_function,
    
    g = g_function,
    
    fitted = h_hat,
    
    residuals = residuals,
    
    transformed_covariates = Z,
    
    bandwidths = bandwidths,
    
    iterations = iteration,
    
    convergence = convergence,
    
    objective_history =
      objective_history[
        seq_len(iteration)
      ],
    
    response_type = response_type
  )
  
  return(result)
}

# ============================================================
# TEST PART 3
# ============================================================

test_fit_data <- generate_baseline(
  n = 100,
  sigma = 0.5
)

test_fit <- fit_functional_transformation(
  dat = test_fit_data,
  max_iter = 10,
  tol = 1e-4
)

cat("\nPart 3 completed successfully.\n")

cat(
  "Response type: ",
  test_fit$response_type,
  "\n"
)

cat(
  "Number of iterations: ",
  test_fit$iterations,
  "\n"
)

cat(
  "Convergence value: ",
  test_fit$convergence,
  "\n"
)

cat(
  "Estimated bandwidths:\n"
)

print(
  test_fit$bandwidths
)

cat(
  "\nAvailable fitted components:\n"
)

print(
  names(test_fit)
)
result <- list(
  
  f1 = if (p >= 1) {
    f_functions[[1]]
  } else {
    NULL
  },
  
  f2 = if (p >= 2) {
    f_functions[[2]]
  } else {
    NULL
  },
  
  f3 = if (p >= 3) {
    f_functions[[3]]
  } else {
    NULL
  },
  
  f4 = if (p >= 4) {
    f_functions[[4]]
  } else {
    NULL
  },
  
  h = h_function,
  
  g = g_function,
  
  fitted = h_hat,
  
  residuals = residuals,
  
  transformed_covariates = Z,
  
  bandwidths = bandwidths,
  
  iterations = iteration,
  
  convergence = convergence,
  
  objective_history =
    objective_history[
      seq_len(iteration)
    ],
  
  response_type = response_type
)

return(result)

}

##################################################################################

# ============================================================
# FUNCTIONAL TRANSFORMATION REGRESSION MODEL
# PART 4: MODEL EVALUATION
# ============================================================

# ------------------------------------------------------------
# Generate one test dataset
# ------------------------------------------------------------

evaluation_data <- generate_baseline(
  n = 200,
  sigma = 0.5,
  response_transformation = "log",
  error_distribution = "normal"
)

# ------------------------------------------------------------
# Fit the functional transformation model
# ------------------------------------------------------------

evaluation_fit <- fit_functional_transformation(
  dat = evaluation_data,
  max_iter = 30,
  tol = 1e-4,
  bandwidth_multiplier = 1
)

# ------------------------------------------------------------
# Estimate f1 on the common grid
# ------------------------------------------------------------

f1_hat_grid <- evaluation_fit$f1(
  grid_x
)

# ------------------------------------------------------------
# Estimate f2 on the common grid
# ------------------------------------------------------------

f2_hat_grid <- evaluation_fit$f2(
  grid_x
)

# ------------------------------------------------------------
# Estimate h on the two-dimensional grid
# ------------------------------------------------------------

h_hat_grid <- evaluation_fit$h(
  grid_h$u,
  grid_h$v
)

# ------------------------------------------------------------
# Calculate performance measures
# ------------------------------------------------------------

f1_metrics <- calc_function_metrics(
  estimated = f1_hat_grid,
  truth = f1_true_grid
)

f2_metrics <- calc_function_metrics(
  estimated = f2_hat_grid,
  truth = f2_true_grid
)

h_metrics <- calc_function_metrics(
  estimated = h_hat_grid,
  truth = h_true_grid
)

# ------------------------------------------------------------
# Display results
# ------------------------------------------------------------

cat("\n============================================\n")
cat("PART 4: MODEL EVALUATION\n")
cat("============================================\n")

cat("\nNumber of observations:\n")
print(
  nrow(evaluation_data)
)

cat("\nNumber of iterations:\n")
print(
  evaluation_fit$iterations
)

cat("\nConvergence:\n")
print(
  evaluation_fit$convergence
)

cat("\nBandwidths:\n")
print(
  evaluation_fit$bandwidths
)

cat("\n\nf1 performance:\n")
print(
  f1_metrics
)

cat("\n\nf2 performance:\n")
print(
  f2_metrics
)

cat("\n\nh performance:\n")
print(
  h_metrics
)

# ------------------------------------------------------------
# Combine metrics into one data frame
# ------------------------------------------------------------

evaluation_summary <- data.frame(
  
  component = c(
    "f1",
    "f2",
    "h"
  ),
  
  Bias = c(
    f1_metrics["Bias"],
    f2_metrics["Bias"],
    h_metrics["Bias"]
  ),
  
  AbsBias = c(
    f1_metrics["AbsBias"],
    f2_metrics["AbsBias"],
    h_metrics["AbsBias"]
  ),
  
  RMSE = c(
    f1_metrics["RMSE"],
    f2_metrics["RMSE"],
    h_metrics["RMSE"]
  ),
  
  ISE = c(
    f1_metrics["ISE"],
    f2_metrics["ISE"],
    h_metrics["ISE"]
  )
)

cat("\n\nCombined evaluation summary:\n")

print(
  evaluation_summary
)

# ------------------------------------------------------------
# Check for invalid values
# ------------------------------------------------------------

cat("\n\nChecking estimated functions...\n")

cat(
  "f1 finite values: ",
  sum(is.finite(f1_hat_grid)),
  " / ",
  length(f1_hat_grid),
  "\n"
)

cat(
  "f2 finite values: ",
  sum(is.finite(f2_hat_grid)),
  " / ",
  length(f2_hat_grid),
  "\n"
)

cat(
  "h finite values: ",
  sum(is.finite(h_hat_grid)),
  " / ",
  length(h_hat_grid),
  "\n"
)

# ------------------------------------------------------------
# Simple plots
# ------------------------------------------------------------

plot_data_f1 <- data.frame(
  x = grid_x,
  Truth = f1_true_grid,
  Estimate = f1_hat_grid
)

plot_data_f2 <- data.frame(
  x = grid_x,
  Truth = f2_true_grid,
  Estimate = f2_hat_grid
)

# ------------------------------------------------------------
# Plot f1
# ------------------------------------------------------------

p_f1 <- ggplot(
  plot_data_f1,
  aes(x = x)
) +
  geom_line(
    aes(
      y = Truth,
      linetype = "True"
    ),
    linewidth = 0.8
  ) +
  geom_line(
    aes(
      y = Estimate,
      linetype = "Estimated"
    ),
    linewidth = 0.8
  ) +
  labs(
    title = "True and Estimated f1",
    x = "x",
    y = "f1(x)",
    linetype = ""
  ) +
  theme_bw()

print(
  p_f1
)

# ------------------------------------------------------------
# Plot f2
# ------------------------------------------------------------

p_f2 <- ggplot(
  plot_data_f2,
  aes(x = x)
) +
  geom_line(
    aes(
      y = Truth,
      linetype = "True"
    ),
    linewidth = 0.8
  ) +
  geom_line(
    aes(
      y = Estimate,
      linetype = "Estimated"
    ),
    linewidth = 0.8
  ) +
  labs(
    title = "True and Estimated f2",
    x = "x",
    y = "f2(x)",
    linetype = ""
  ) +
  theme_bw()

print(
  p_f2
)

cat("\nPart 4 completed successfully.\n")

#################################################################################

# ============================================================
# FUNCTIONAL TRANSFORMATION REGRESSION MODEL
# PART 5: FINITE-SAMPLE MONTE CARLO SIMULATION
# ============================================================

# ------------------------------------------------------------
# Storage object
# ------------------------------------------------------------

simulation_results <- list()

# ------------------------------------------------------------
# Start simulation
# ------------------------------------------------------------

simulation_start_time <- Sys.time()

for (n in sample_sizes) {
  
  cat("\n============================================\n")
  cat("Sample size n =", n, "\n")
  cat("============================================\n")
  
  results_n <- vector(
    "list",
    R
  )
  
  for (r in seq_len(R)) {
    
    # --------------------------------------------------------
    # Generate data
    # --------------------------------------------------------
    
    dat <- generate_baseline(
      n = n,
      sigma = 0.5,
      response_transformation = "log",
      error_distribution = "normal"
    )
    
    # --------------------------------------------------------
    # Fit model
    # --------------------------------------------------------
    
    fit <- tryCatch(
      
      fit_functional_transformation(
        dat = dat,
        max_iter = 30,
        tol = 1e-4,
        bandwidth_multiplier = 1
      ),
      
      error = function(e) {
        
        message(
          "Error in replication ",
          r,
          " for n = ",
          n,
          ": ",
          e$message
        )
        
        NULL
      }
    )
    
    # --------------------------------------------------------
    # Handle failed replication
    # --------------------------------------------------------
    
    if (is.null(fit)) {
      
      results_n[[r]] <- NULL
      
      next
    }
    
    # --------------------------------------------------------
    # Estimate f1
    # --------------------------------------------------------
    
    f1_hat <- tryCatch(
      
      fit$f1(grid_x),
      
      error = function(e) {
        
        rep(
          NA_real_,
          length(grid_x)
        )
      }
    )
    
    # --------------------------------------------------------
    # Estimate f2
    # --------------------------------------------------------
    
    f2_hat <- tryCatch(
      
      fit$f2(grid_x),
      
      error = function(e) {
        
        rep(
          NA_real_,
          length(grid_x)
        )
      }
    )
    
    # --------------------------------------------------------
    # Estimate h
    # --------------------------------------------------------
    
    h_hat <- tryCatch(
      
      fit$h(
        grid_h$u,
        grid_h$v
      ),
      
      error = function(e) {
        
        rep(
          NA_real_,
          length(h_true_grid)
        )
      }
    )
    
    # --------------------------------------------------------
    # Calculate performance measures
    # --------------------------------------------------------
    
    f1_metrics <- calc_function_metrics(
      estimated = f1_hat,
      truth = f1_true_grid
    )
    
    f2_metrics <- calc_function_metrics(
      estimated = f2_hat,
      truth = f2_true_grid
    )
    
    h_metrics <- calc_function_metrics(
      estimated = h_hat,
      truth = h_true_grid
    )
    
    # --------------------------------------------------------
    # Store results
    # --------------------------------------------------------
    
    results_n[[r]] <- list(
      
      f1 = f1_metrics,
      
      f2 = f2_metrics,
      
      h = h_metrics,
      
      iterations = fit$iterations,
      
      convergence = fit$convergence
    )
    
    # --------------------------------------------------------
    # Progress message
    # --------------------------------------------------------
    
    if (
      r %% 50 == 0 ||
      r == 1 ||
      r == R
    ) {
      
      cat(
        "n = ",
        n,
        " | replication = ",
        r,
        "/",
        R,
        "\n",
        sep = ""
      )
    }
  }
  
  # ----------------------------------------------------------
  # Store results for this sample size
  # ----------------------------------------------------------
  
  simulation_results[[as.character(n)]] <- results_n
}

# ------------------------------------------------------------
# Simulation time
# ------------------------------------------------------------

simulation_end_time <- Sys.time()

simulation_time <- simulation_end_time -
  simulation_start_time

cat("\n============================================\n")
cat("FINITE-SAMPLE SIMULATION COMPLETED\n")
cat("============================================\n")

print(
  simulation_time
)

# ------------------------------------------------------------
# Number of successful replications
# ------------------------------------------------------------

successful_replications <- data.frame(
  n = sample_sizes,
  successful = NA_integer_,
  failed = NA_integer_
)

for (i in seq_along(sample_sizes)) {
  
  n <- sample_sizes[i]
  
  current_results <-
    simulation_results[
      [as.character(n)]
    ]
  
  successful <- sum(
    !vapply(
      current_results,
      is.null,
      logical(1)
    )
  )
  
  successful_replications$successful[i] <-
    successful
  
  successful_replications$failed[i] <-
    R - successful
}

cat("\nSuccessful and failed replications:\n")

print(
  successful_replications
)

####################################################################################

# ============================================================
# FUNCTIONAL TRANSFORMATION REGRESSION MODEL
# PART 6: FINITE-SAMPLE SUMMARY
# ============================================================

# ------------------------------------------------------------
# Check that simulation results exist
# ------------------------------------------------------------

if (!exists("simulation_results")) {
  
  stop(
    "simulation_results does not exist. ",
    "Please run Part 5 first."
  )
}

# ------------------------------------------------------------
# Empty data frame for extracted results
# ------------------------------------------------------------

rmse_results <- data.frame(
  n = integer(),
  replication = integer(),
  component = character(),
  Bias = numeric(),
  AbsBias = numeric(),
  RMSE = numeric(),
  ISE = numeric(),
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# Extract results
# ------------------------------------------------------------

for (n in sample_sizes) {
  
  current_results <-
    simulation_results[[as.character(n)]]
  
  if (is.null(current_results)) {
    next
  }
  
  for (r in seq_along(current_results)) {
    
    current_replication <-
      current_results[[r]]
    
    if (is.null(current_replication)) {
      next
    }
    
    component_names <- c(
      "f1",
      "f2",
      "h"
    )
    
    for (
      component_name
      in component_names
    ) {
      
      metric <-
        current_replication[
          [component_name]
        ]
      
      if (is.null(metric)) {
        next
      }
      
      rmse_results <- rbind(
        rmse_results,
        
        data.frame(
          n = n,
          
          replication = r,
          
          component =
            component_name,
          
          Bias =
            as.numeric(
              metric["Bias"]
            ),
          
          AbsBias =
            as.numeric(
              metric["AbsBias"]
            ),
          
          RMSE =
            as.numeric(
              metric["RMSE"]
            ),
          
          ISE =
            as.numeric(
              metric["ISE"]
            ),
          
          stringsAsFactors =
            FALSE
        )
      )
    }
  }
}

# ------------------------------------------------------------
# Correct extraction syntax
# ------------------------------------------------------------

# The following line is the correct R syntax:
#
# metric <- current_replication[[component_name]]
#
# If the block above produces a syntax error because of the
# bracket formatting, use the replacement loop below.
# ------------------------------------------------------------

rmse_results <- data.frame(
  n = integer(),
  replication = integer(),
  component = character(),
  Bias = numeric(),
  AbsBias = numeric(),
  RMSE = numeric(),
  ISE = numeric(),
  stringsAsFactors = FALSE
)

for (n in sample_sizes) {
  
  current_results <-
    simulation_results[[as.character(n)]]
  
  if (is.null(current_results)) {
    next
  }
  
  for (r in seq_along(current_results)) {
    
    current_replication <-
      current_results[[r]]
    
    if (is.null(current_replication)) {
      next
    }
    
    for (
      component_name
      in c("f1", "f2", "h")
    ) {
      
      metric <-
        current_replication[[component_name]]
      
      if (is.null(metric)) {
        next
      }
      
      rmse_results <- rbind(
        rmse_results,
        
        data.frame(
          n = n,
          replication = r,
          component = component_name,
          Bias =
            as.numeric(
              metric["Bias"]
            ),
          AbsBias =
            as.numeric(
              metric["AbsBias"]
            ),
          RMSE =
            as.numeric(
              metric["RMSE"]
            ),
          ISE =
            as.numeric(
              metric["ISE"]
            ),
          stringsAsFactors = FALSE
        )
      )
    }
  }
}

# ------------------------------------------------------------
# Check extracted results
# ------------------------------------------------------------

cat("\n============================================\n")
cat("EXTRACTED SIMULATION RESULTS\n")
cat("============================================\n")

cat(
  "\nNumber of rows: ",
  nrow(rmse_results),
  "\n"
)

cat(
  "Number of columns: ",
  ncol(rmse_results),
  "\n"
)

if (nrow(rmse_results) > 0) {
  
  print(
    head(
      rmse_results,
      10
    )
  )
  
} else {
  
  warning(
    "No simulation results were extracted."
  )
}

# ------------------------------------------------------------
# Finite-sample summary
# ------------------------------------------------------------

if (nrow(rmse_results) > 0) {
  
  finite_sample_summary <-
    aggregate(
      cbind(
        Bias,
        AbsBias,
        RMSE,
        ISE
      ) ~ n + component,
      data = rmse_results,
      FUN = mean,
      na.rm = TRUE
    )
  
} else {
  
  finite_sample_summary <- data.frame()
}

# ------------------------------------------------------------
# Display summary
# ------------------------------------------------------------

cat(
  "\n============================================\n"
)

cat(
  "FINITE-SAMPLE SUMMARY\n"
)

cat(
  "============================================\n"
)

print(
  finite_sample_summary
)

# ------------------------------------------------------------
# Save summary
# ------------------------------------------------------------

write.csv(
  finite_sample_summary,
  file = "finite_sample_summary.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# RMSE plot
# ------------------------------------------------------------

if (nrow(rmse_results) > 0) {
  
  p_rmse <- ggplot(
    rmse_results,
    aes(
      x = n,
      y = RMSE,
      linetype = component
    )
  ) +
    geom_point(
      alpha = 0.35
    ) +
    geom_smooth(
      method = "loess",
      se = FALSE
    ) +
    labs(
      title =
        "Finite-Sample RMSE",
      x =
        "Sample size",
      y =
        "RMSE",
      linetype =
        "Component"
    ) +
    theme_bw()
  
  print(
    p_rmse
  )
  
  ggsave(
    filename =
      "RMSE_vs_sample_size.png",
    plot = p_rmse,
    width = 8,
    height = 5,
    dpi = 300
  )
}

# ------------------------------------------------------------
# ISE plot
# ------------------------------------------------------------

if (nrow(rmse_results) > 0) {
  
  p_ise <- ggplot(
    rmse_results,
    aes(
      x = n,
      y = ISE,
      linetype = component
    )
  ) +
    geom_point(
      alpha = 0.35
    ) +
    geom_smooth(
      method = "loess",
      se = FALSE
    ) +
    labs(
      title =
        "Finite-Sample Integrated Squared Error",
      x =
        "Sample size",
      y =
        "ISE",
      linetype =
        "Component"
    ) +
    theme_bw()
  
  print(
    p_ise
  )
  
  ggsave(
    filename =
      "ISE_vs_sample_size.png",
    plot = p_ise,
    width = 8,
    height = 5,
    dpi = 300
  )
}

cat(
  "\nPart 6 completed successfully.\n"
)

# ============================================================
# FUNCTIONAL TRANSFORMATION REGRESSION MODEL
# PART 7: ASYMPTOTIC NORMALITY ANALYSIS
# ============================================================

# ------------------------------------------------------------
# Bias functions
# ------------------------------------------------------------

bias_f1 <- function(
    x,
    h) {
  
  # f1(x) = sin(2*pi*x)
  #
  # f1''(x) = -4*pi^2*sin(2*pi*x)
  #
  # Leading kernel bias:
  # (mu2 / 2) * h^2 * f1''(x)
  
  (
    mu2_epanechnikov / 2
  ) *
    h^2 *
    (
      -4 *
        pi^2 *
        sin(
          2 * pi * x
        )
    )
}


bias_f2 <- function(
    x,
    h) {
  
  # f2(x) = x^2
  #
  # f2''(x) = 2
  #
  # Leading kernel bias:
  # (mu2 / 2) * h^2 * 2
  
  (
    mu2_epanechnikov / 2
  ) *
    h^2 *
    2
}


# ------------------------------------------------------------
# Evaluation points
# ------------------------------------------------------------

asymptotic_x_values <- c(
  0.25,
  0.50,
  0.75
)


# ------------------------------------------------------------
# Storage object
# ------------------------------------------------------------

asymptotic_results <- data.frame(
  n = integer(),
  x = numeric(),
  component = character(),
  estimate = numeric(),
  truth = numeric(),
  bias_approx = numeric(),
  standardized_error = numeric(),
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# Simulation
# ------------------------------------------------------------

for (n in sample_sizes) {
  
  cat(
    "\n--------------------------------------------\n"
  )
  
  cat(
    "Asymptotic analysis: n = ",
    n,
    "\n",
    sep = ""
  )
  
  cat(
    "--------------------------------------------\n"
  )
  
  
  # ----------------------------------------------------------
  # Approximate bandwidth
  # ----------------------------------------------------------
  
  reference_x <- runif(
    n,
    0,
    1
  )
  
  reference_bandwidth <-
    default_bandwidth(
      reference_x
    )
  
  
  # ----------------------------------------------------------
  # Loop over evaluation points
  # ----------------------------------------------------------
  
  for (x0 in asymptotic_x_values) {
    
    # --------------------------------------------------------
    # True values
    # --------------------------------------------------------
    
    true_f1 <-
      f1_true_function(x0)
    
    true_f2 <-
      f2_true_function(x0)
    
    
    # --------------------------------------------------------
    # Leading bias approximations
    # --------------------------------------------------------
    
    bias1 <-
      bias_f1(
        x = x0,
        h = reference_bandwidth
      )
    
    bias2 <-
      bias_f2(
        x = x0,
        h = reference_bandwidth
      )
    
    
    # --------------------------------------------------------
    # Monte Carlo replications
    # --------------------------------------------------------
    
    for (r in seq_len(R)) {
      
      dat <- generate_baseline(
        n = n,
        sigma = 0.5,
        response_transformation = "log",
        error_distribution = "normal"
      )
      
      
      # ------------------------------------------------------
      # Fit model
      # ------------------------------------------------------
      
      fit <- tryCatch(
        
        fit_functional_transformation(
          dat = dat,
          max_iter = 30,
          tol = 1e-4,
          bandwidth_multiplier = 1
        ),
        
        error = function(e) {
          NULL
        }
      )
      
      
      if (is.null(fit)) {
        next
      }
      
      
      # ------------------------------------------------------
      # Estimate f1 at x0
      # ------------------------------------------------------
      
      estimate_f1 <- tryCatch(
        
        as.numeric(
          fit$f1(x0)
        ),
        
        error = function(e) {
          NA_real_
        }
      )
      
      
      # ------------------------------------------------------
      # Estimate f2 at x0
      # ------------------------------------------------------
      
      estimate_f2 <- tryCatch(
        
        as.numeric(
          fit$f2(x0)
        ),
        
        error = function(e) {
          NA_real_
        }
      )
      
      
      # ------------------------------------------------------
      # Standard error scale
      # ------------------------------------------------------
      
      se_scale <- 1 /
        sqrt(
          n *
            reference_bandwidth
        )
      
      
      # ------------------------------------------------------
      # Standardized errors
      # ------------------------------------------------------
      
      standardized_f1 <-
        (
          estimate_f1 -
            true_f1 -
            bias1
        ) /
        se_scale
      
      
      standardized_f2 <-
        (
          estimate_f2 -
            true_f2 -
            bias2
        ) /
        se_scale
      
      
      # ------------------------------------------------------
      # Store f1
      # ------------------------------------------------------
      
      asymptotic_results <-
        rbind(
          asymptotic_results,
          
          data.frame(
            n = n,
            x = x0,
            component = "f1",
            estimate = estimate_f1,
            truth = true_f1,
            bias_approx = bias1,
            standardized_error =
              standardized_f1,
            stringsAsFactors = FALSE
          )
        )
      
      
      # ------------------------------------------------------
      # Store f2
      # ------------------------------------------------------
      
      asymptotic_results <-
        rbind(
          asymptotic_results,
          
          data.frame(
            n = n,
            x = x0,
            component = "f2",
            estimate = estimate_f2,
            truth = true_f2,
            bias_approx = bias2,
            standardized_error =
              standardized_f2,
            stringsAsFactors = FALSE
          )
        )
      
      
      # ------------------------------------------------------
      # Progress
      # ------------------------------------------------------
      
      if (
        r %% 100 == 0 ||
        r == R
      ) {
        
        cat(
          "x = ",
          x0,
          " | replication = ",
          r,
          "/",
          R,
          "\n",
          sep = ""
        )
      }
    }
  }
}


# ============================================================
# Summary of standardized errors
# ============================================================

asymptotic_normality_summary <-
  aggregate(
    standardized_error ~
      n +
      x +
      component,
    data = asymptotic_results,
    FUN = function(z) {
      
      c(
        Mean = mean(
          z,
          na.rm = TRUE
        ),
        
        SD = sd(
          z,
          na.rm = TRUE
        ),
        
        Skewness = moments::skewness(
          z,
          na.rm = TRUE
        ),
        
        Kurtosis = moments::kurtosis(
          z,
          na.rm = TRUE
        )
      )
    }
  )


# ------------------------------------------------------------
# Convert aggregated matrix column to separate columns
# ------------------------------------------------------------

summary_matrix <-
  asymptotic_normality_summary[
    ,
    "standardized_error",
    drop = FALSE
  ]


x <- asymptotic_normality_summary$standardized_error

x <- x[
  is.finite(x)
]

n_x <- length(x)

Mean <- mean(x)

SD <- sd(x)

Skewness <-
  mean(
    (x - Mean)^3
  ) /
  SD^3

Kurtosis <-
  mean(
    (x - Mean)^4
  ) /
  SD^4


summary_matrix <-
  data.frame(
    Mean = Mean,
    SD = SD,
    Skewness = Skewness,
    Kurtosis = Kurtosis
  )

asymptotic_normality_summary_final <-
  data.frame(
    
    asymptotic_normality_summary[
      ,
      c(
        "n",
        "x",
        "component"
      ),
      drop = FALSE
    ],
    
    Mean = summary_matrix$Mean,
    
    SD = summary_matrix$SD,
    
    Skewness = summary_matrix$Skewness,
    
    Kurtosis = summary_matrix$Kurtosis,
    
    row.names = NULL
  )
# ============================================================
# Q-Q plots
# ============================================================

if (
  nrow(asymptotic_results) > 0
) {
  
  for (component_name in c(
    "f1",
    "f2"
  )) {
    
    for (x0 in asymptotic_x_values) {
      
      plot_subset <-
        asymptotic_results[
          asymptotic_results$component ==
            component_name &
            asymptotic_results$x ==
            x0,
          ,
          drop = FALSE
        ]
      
      if (
        nrow(plot_subset) < 10
      ) {
        next
      }
      
      
      qq_plot <- ggplot(
        plot_subset,
        aes(
          sample = standardized_error
        )
      ) +
        
        stat_qq() +
        
        stat_qq_line() +
        
        labs(
          title = paste(
            "Q-Q Plot:",
            component_name,
            "at x =",
            x0
          ),
          x = "Theoretical Quantiles",
          y = "Standardized Error"
        ) +
        
        theme_bw()
      
      
      print(
        qq_plot
      )
      
      
      file_name <- paste0(
        "QQ_",
        component_name,
        "_x_",
        gsub(
          "\\.",
          "_",
          as.character(x0)
        ),
        ".png"
      )
      
      
      ggsave(
        filename = file_name,
        plot = qq_plot,
        width = 6,
        height = 5,
        dpi = 300
      )
    }
  }
}


cat(
  "\nPart 7 completed successfully.\n"
)

###################################################################################

# ============================================================
# FUNCTIONAL TRANSFORMATION REGRESSION MODEL
# PART 8: RESPONSE TRANSFORMATION EXPERIMENT
# ============================================================

# ------------------------------------------------------------
# Settings
# ------------------------------------------------------------

transformation_sample_size <- 400

transformation_types <- c(
  "identity",
  "log"
)

# ------------------------------------------------------------
# Storage
# ------------------------------------------------------------

response_transformation_results <- data.frame(
  transformation = character(),
  replication = integer(),
  f1_RMSE = numeric(),
  f2_RMSE = numeric(),
  h_RMSE = numeric(),
  residual_RMSE = numeric(),
  iterations = numeric(),
  convergence = numeric(),
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# Start experiment
# ------------------------------------------------------------

transformation_start_time <- Sys.time()

for (
  transformation_name
  in transformation_types
) {
  
  cat(
    "\n============================================\n"
  )
  
  cat(
    "Response transformation: ",
    transformation_name,
    "\n",
    sep = ""
  )
  
  cat(
    "============================================\n"
  )
  
  for (r in seq_len(R)) {
    
    # --------------------------------------------------------
    # Generate data
    # --------------------------------------------------------
    
    if (
      transformation_name == "identity"
    ) {
      
      dat <- generate_identity(
        n = transformation_sample_size,
        sigma = 0.5,
        error_distribution = "normal"
      )
      
    } else {
      
      dat <- generate_log_response(
        n = transformation_sample_size,
        sigma = 0.5,
        error_distribution = "normal"
      )
    }
    
    # --------------------------------------------------------
    # Fit model
    # --------------------------------------------------------
    
    fit <- tryCatch(
      
      fit_functional_transformation(
        dat = dat,
        max_iter = 30,
        tol = 1e-4,
        bandwidth_multiplier = 1
      ),
      
      error = function(e) {
        
        message(
          "Transformation = ",
          transformation_name,
          ", replication = ",
          r,
          ": ",
          e$message
        )
        
        NULL
      }
    )
    
    # --------------------------------------------------------
    # Skip failed replication
    # --------------------------------------------------------
    
    if (is.null(fit)) {
      next
    }
    
    # --------------------------------------------------------
    # Estimate f1
    # --------------------------------------------------------
    
    f1_hat <- tryCatch(
      
      fit$f1(
        grid_x
      ),
      
      error = function(e) {
        
        rep(
          NA_real_,
          length(grid_x)
        )
      }
    )
    
    # --------------------------------------------------------
    # Estimate f2
    # --------------------------------------------------------
    
    f2_hat <- tryCatch(
      
      fit$f2(
        grid_x
      ),
      
      error = function(e) {
        
        rep(
          NA_real_,
          length(grid_x)
        )
      }
    )
    
    # --------------------------------------------------------
    # Estimate h
    # --------------------------------------------------------
    
    h_hat <- tryCatch(
      
      fit$h(
        grid_h$u,
        grid_h$v
      ),
      
      error = function(e) {
        
        rep(
          NA_real_,
          length(h_true_grid)
        )
      }
    )
    
    # --------------------------------------------------------
    # Calculate function RMSE
    # --------------------------------------------------------
    
    f1_RMSE <- sqrt(
      mean(
        (
          f1_hat -
            f1_true_grid
        )^2,
        na.rm = TRUE
      )
    )
    
    f2_RMSE <- sqrt(
      mean(
        (
          f2_hat -
            f2_true_grid
        )^2,
        na.rm = TRUE
      )
    )
    
    h_RMSE <- sqrt(
      mean(
        (
          h_hat -
            h_true_grid
        )^2,
        na.rm = TRUE
      )
    )
    
    # --------------------------------------------------------
    # Residual RMSE
    # --------------------------------------------------------
    
    residual_RMSE <- sqrt(
      mean(
        fit$residuals^2,
        na.rm = TRUE
      )
    )
    
    # --------------------------------------------------------
    # Store results
    # --------------------------------------------------------
    
    response_transformation_results <-
      rbind(
        response_transformation_results,
        
        data.frame(
          transformation =
            transformation_name,
          
          replication = r,
          
          f1_RMSE =
            f1_RMSE,
          
          f2_RMSE =
            f2_RMSE,
          
          h_RMSE =
            h_RMSE,
          
          residual_RMSE =
            residual_RMSE,
          
          iterations =
            fit$iterations,
          
          convergence =
            fit$convergence,
          
          stringsAsFactors = FALSE
        )
      )
    
    # --------------------------------------------------------
    # Progress
    # --------------------------------------------------------
    
    if (
      r %% 100 == 0 ||
      r == R
    ) {
      
      cat(
        "Replication ",
        r,
        "/",
        R,
        "\n",
        sep = ""
      )
    }
  }
}

# ------------------------------------------------------------
# End time
# ------------------------------------------------------------

transformation_end_time <- Sys.time()

cat(
  "\nTransformation experiment time:\n"
)

print(
  transformation_end_time -
    transformation_start_time
)

# ============================================================
# Summary
# ============================================================

response_transformation_summary <-
  aggregate(
    cbind(
      f1_RMSE,
      f2_RMSE,
      h_RMSE,
      residual_RMSE,
      iterations,
      convergence
    ) ~ transformation,
    data =
      response_transformation_results,
    FUN = mean,
    na.rm = TRUE
  )

# ------------------------------------------------------------
# Display summary
# ------------------------------------------------------------

cat(
  "\n============================================\n"
)

cat(
  "RESPONSE TRANSFORMATION SUMMARY\n"
)

cat(
  "============================================\n"
)

print(
  response_transformation_summary
)

# ------------------------------------------------------------
# Save summary
# ------------------------------------------------------------

write.csv(
  response_transformation_summary,
  file =
    "response_transformation_summary.csv",
  row.names = FALSE
)

# ============================================================
# Long-format RMSE data for plotting
# ============================================================

response_rmse_plot_data <- rbind(
  
  data.frame(
    transformation =
      response_transformation_results$transformation,
    replication =
      response_transformation_results$replication,
    component = "f1",
    RMSE =
      response_transformation_results$f1_RMSE
  ),
  
  data.frame(
    transformation =
      response_transformation_results$transformation,
    replication =
      response_transformation_results$replication,
    component = "f2",
    RMSE =
      response_transformation_results$f2_RMSE
  ),
  
  data.frame(
    transformation =
      response_transformation_results$transformation,
    replication =
      response_transformation_results$replication,
    component = "h",
    RMSE =
      response_transformation_results$h_RMSE
  )
)

# ============================================================
# Boxplot
# ============================================================

if (
  nrow(response_rmse_plot_data) > 0
) {
  
  p_response <- ggplot(
    response_rmse_plot_data,
    aes(
      x = transformation,
      y = RMSE,
      linetype = component
    )
  ) +
    
    geom_boxplot() +
    
    labs(
      title =
        "Effect of Response Transformation",
      x =
        "Response transformation",
      y =
        "RMSE",
      linetype =
        "Component"
    ) +
    
    theme_bw()
  
  print(
    p_response
  )
  
  ggsave(
    filename =
      "Response_transformation_RMSE.png",
    plot =
      p_response,
    width = 8,
    height = 5,
    dpi = 300
  )
}

cat(
  "\nPart 8 completed successfully.\n"
)

################################################################################

# ============================================================
# FUNCTIONAL TRANSFORMATION REGRESSION MODEL
# PART 9: ERROR-DISTRIBUTION SENSITIVITY
# ============================================================

# ------------------------------------------------------------
# Settings
# ------------------------------------------------------------

error_distributions <- c(
  "normal",
  "t5",
  "gamma"
)

error_sensitivity_n <- 400

# ------------------------------------------------------------
# Storage
# ------------------------------------------------------------

error_distribution_results <- data.frame(
  distribution = character(),
  replication = integer(),
  f1_RMSE = numeric(),
  f2_RMSE = numeric(),
  h_RMSE = numeric(),
  residual_RMSE = numeric(),
  iterations = numeric(),
  convergence = numeric(),
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# Start experiment
# ------------------------------------------------------------

error_start_time <- Sys.time()

for (
  error_name
  in error_distributions
) {
  
  cat(
    "\n============================================\n"
  )
  
  cat(
    "Error distribution: ",
    error_name,
    "\n",
    sep = ""
  )
  
  cat(
    "============================================\n"
  )
  
  for (r in seq_len(R)) {
    
    # --------------------------------------------------------
    # Generate data
    # --------------------------------------------------------
    
    dat <- generate_baseline(
      n = error_sensitivity_n,
      sigma = 0.5,
      response_transformation = "log",
      error_distribution = error_name
    )
    
    # --------------------------------------------------------
    # Fit model
    # --------------------------------------------------------
    
    fit <- tryCatch(
      
      fit_functional_transformation(
        dat = dat,
        max_iter = 30,
        tol = 1e-4,
        bandwidth_multiplier = 1
      ),
      
      error = function(e) {
        
        message(
          "Distribution = ",
          error_name,
          ", replication = ",
          r,
          ": ",
          e$message
        )
        
        NULL
      }
    )
    
    # --------------------------------------------------------
    # Skip failed replication
    # --------------------------------------------------------
    
    if (is.null(fit)) {
      next
    }
    
    # --------------------------------------------------------
    # Estimate f1
    # --------------------------------------------------------
    
    f1_hat <- tryCatch(
      
      fit$f1(
        grid_x
      ),
      
      error = function(e) {
        
        rep(
          NA_real_,
          length(grid_x)
        )
      }
    )
    
    # --------------------------------------------------------
    # Estimate f2
    # --------------------------------------------------------
    
    f2_hat <- tryCatch(
      
      fit$f2(
        grid_x
      ),
      
      error = function(e) {
        
        rep(
          NA_real_,
          length(grid_x)
        )
      }
    )
    
    # --------------------------------------------------------
    # Estimate h
    # --------------------------------------------------------
    
    h_hat <- tryCatch(
      
      fit$h(
        grid_h$u,
        grid_h$v
      ),
      
      error = function(e) {
        
        rep(
          NA_real_,
          length(h_true_grid)
        )
      }
    )
    
    # --------------------------------------------------------
    # RMSE of f1
    # --------------------------------------------------------
    
    f1_RMSE <- sqrt(
      mean(
        (
          f1_hat -
            f1_true_grid
        )^2,
        na.rm = TRUE
      )
    )
    
    # --------------------------------------------------------
    # RMSE of f2
    # --------------------------------------------------------
    
    f2_RMSE <- sqrt(
      mean(
        (
          f2_hat -
            f2_true_grid
        )^2,
        na.rm = TRUE
      )
    )
    
    # --------------------------------------------------------
    # RMSE of h
    # --------------------------------------------------------
    
    h_RMSE <- sqrt(
      mean(
        (
          h_hat -
            h_true_grid
        )^2,
        na.rm = TRUE
      )
    )
    
    # --------------------------------------------------------
    # Residual RMSE
    # --------------------------------------------------------
    
    residual_RMSE <- sqrt(
      mean(
        fit$residuals^2,
        na.rm = TRUE
      )
    )
    
    # --------------------------------------------------------
    # Store results
    # --------------------------------------------------------
    
    error_distribution_results <-
      rbind(
        error_distribution_results,
        
        data.frame(
          distribution =
            error_name,
          
          replication =
            r,
          
          f1_RMSE =
            f1_RMSE,
          
          f2_RMSE =
            f2_RMSE,
          
          h_RMSE =
            h_RMSE,
          
          residual_RMSE =
            residual_RMSE,
          
          iterations =
            fit$iterations,
          
          convergence =
            fit$convergence,
          
          stringsAsFactors =
            FALSE
        )
      )
    
    # --------------------------------------------------------
    # Progress
    # --------------------------------------------------------
    
    if (
      r %% 100 == 0 ||
      r == R
    ) {
      
      cat(
        "Replication ",
        r,
        "/",
        R,
        "\n",
        sep = ""
      )
    }
  }
}

# ------------------------------------------------------------
# End time
# ------------------------------------------------------------

error_end_time <- Sys.time()

cat(
  "\nError-distribution experiment time:\n"
)

print(
  error_end_time -
    error_start_time
)

# ============================================================
# Summary
# ============================================================

error_distribution_summary <-
  aggregate(
    cbind(
      f1_RMSE,
      f2_RMSE,
      h_RMSE,
      residual_RMSE,
      iterations,
      convergence
    ) ~ distribution,
    data =
      error_distribution_results,
    FUN = mean,
    na.rm = TRUE
  )

# ------------------------------------------------------------
# Display
# ------------------------------------------------------------

cat(
  "\n============================================\n"
)

cat(
  "ERROR-DISTRIBUTION SUMMARY\n"
)

cat(
  "============================================\n"
)

print(
  error_distribution_summary
)

# ------------------------------------------------------------
# Save summary
# ------------------------------------------------------------

write.csv(
  error_distribution_summary,
  file =
    "error_distribution_summary.csv",
  row.names = FALSE
)

# ============================================================
# Long-format RMSE data
# ============================================================

error_rmse_plot_data <- rbind(
  
  data.frame(
    distribution =
      error_distribution_results$distribution,
    
    replication =
      error_distribution_results$replication,
    
    component = "f1",
    
    RMSE =
      error_distribution_results$f1_RMSE
  ),
  
  data.frame(
    distribution =
      error_distribution_results$distribution,
    
    replication =
      error_distribution_results$replication,
    
    component = "f2",
    
    RMSE =
      error_distribution_results$f2_RMSE
  ),
  
  data.frame(
    distribution =
      error_distribution_results$distribution,
    
    replication =
      error_distribution_results$replication,
    
    component = "h",
    
    RMSE =
      error_distribution_results$h_RMSE
  )
)

# ============================================================
# Boxplot
# ============================================================

if (
  nrow(error_rmse_plot_data) > 0
) {
  
  p_error <- ggplot(
    error_rmse_plot_data,
    aes(
      x = distribution,
      y = RMSE,
      linetype = component
    )
  ) +
    
    geom_boxplot() +
    
    labs(
      title =
        "Sensitivity to Error Distribution",
      x =
        "Error distribution",
      y =
        "RMSE",
      linetype =
        "Component"
    ) +
    
    theme_bw()
  
  print(
    p_error
  )
  
  ggsave(
    filename =
      "Error_distribution_RMSE.png",
    plot =
      p_error,
    width = 8,
    height = 5,
    dpi = 300
  )
}

cat(
  "\nPart 9 completed successfully.\n"
)

################################################################################

# ============================================================
# FUNCTIONAL TRANSFORMATION REGRESSION MODEL
# PART 10: DIMENSIONALITY STUDY
# ============================================================

# ------------------------------------------------------------
# Generate covariates
# ------------------------------------------------------------

generate_covariates <- function(
    n,
    p) {
  
  X <- matrix(
    runif(
      n * p,
      min = 0,
      max = 1
    ),
    nrow = n,
    ncol = p
  )
  
  colnames(X) <-
    paste0(
      "X",
      seq_len(p)
    )
  
  X
}

# ------------------------------------------------------------
# Transform covariates
# ------------------------------------------------------------

transform_covariates <- function(
    X) {
  
  X <- as.matrix(X)
  
  n <- nrow(X)
  
  p <- ncol(X)
  
  Z <- matrix(
    NA_real_,
    nrow = n,
    ncol = p
  )
  
  # First component
  Z[, 1] <-
    sin(
      2 * pi * X[, 1]
    )
  
  # Second component
  if (p >= 2) {
    
    Z[, 2] <-
      X[, 2]^2
  }
  
  # Third component
  if (p >= 3) {
    
    Z[, 3] <-
      exp(
        X[, 3] - 1
      ) - 1
  }
  
  # Fourth component
  if (p >= 4) {
    
    Z[, 4] <-
      log(
        1 + X[, 4]
      )
  }
  
  colnames(Z) <-
    paste0(
      "Z",
      seq_len(p)
    )
  
  Z
}

# ------------------------------------------------------------
# Generate multivariate h
# ------------------------------------------------------------

generate_h <- function(
    Z) {
  
  Z <- as.matrix(Z)
  
  p <- ncol(Z)
  
  eta <- rowSums(Z)
  
  if (p >= 2) {
    
    eta <-
      eta +
      0.5 *
      Z[, 1] *
      Z[, 2]
  }
  
  if (p >= 3) {
    
    eta <-
      eta +
      0.25 *
      Z[, 1] *
      Z[, 3]
  }
  
  if (p >= 4) {
    
    eta <-
      eta +
      0.20 *
      Z[, 2] *
      Z[, 4]
  }
  
  eta
}

# ------------------------------------------------------------
# Generate dimensionality-study data
# ------------------------------------------------------------

generate_dimension_data <- function(
    n,
    p,
    sigma = 0.5) {
  
  X <- generate_covariates(
    n = n,
    p = p
  )
  
  Z <- transform_covariates(
    X
  )
  
  eta <- generate_h(
    Z
  )
  
  eps <- sigma *
    rnorm(n)
  
  Y <- exp(
    eta + eps
  )
  
  data <- data.frame(
    Y = Y
  )
  
  for (j in seq_len(p)) {
    
    data[[paste0("X", j)]] <- X[, j]
    
  }
  
  data
}

# ------------------------------------------------------------
# Settings
# ------------------------------------------------------------

dimension_values <- c(
  2,
  3,
  4
)

dimension_sample_size <- 400

# ------------------------------------------------------------
# Storage
# ------------------------------------------------------------

dimension_results <- data.frame(
  p = integer(),
  replication = integer(),
  RMSE = numeric(),
  residual_RMSE = numeric(),
  iterations = numeric(),
  convergence = numeric(),
  computation_time = numeric(),
  stringsAsFactors = FALSE
)

# ============================================================
# Dimensionality simulation
# ============================================================

dimension_start_time <- Sys.time()

for (p in dimension_values) {
  
  cat(
    "\n============================================\n"
  )
  
  cat(
    "Dimension p = ",
    p,
    "\n",
    sep = ""
  )
  
  cat(
    "============================================\n"
  )
  
  for (r in seq_len(R)) {
    
    # --------------------------------------------------------
    # Generate data
    # --------------------------------------------------------
    
    dat <- generate_dimension_data(
      n = dimension_sample_size,
      p = p,
      sigma = 0.5
    )
    
    # --------------------------------------------------------
    # Fit model and measure time
    # --------------------------------------------------------
    
    start_time <- Sys.time()
    
    fit <- tryCatch(
      
      fit_functional_transformation(
        dat = dat,
        max_iter = 30,
        tol = 1e-4,
        bandwidth_multiplier = 1
      ),
      
      error = function(e) {
        
        message(
          "p = ",
          p,
          ", replication = ",
          r,
          ": ",
          e$message
        )
        
        NULL
      }
    )
    
    end_time <- Sys.time()
    
    computation_time <-
      as.numeric(
        end_time -
          start_time,
        units = "secs"
      )
    
    # --------------------------------------------------------
    # Skip failed fit
    # --------------------------------------------------------
    
    if (is.null(fit)) {
      next
    }
    
    # --------------------------------------------------------
    # Residual RMSE
    # --------------------------------------------------------
    
    residual_RMSE <- sqrt(
      mean(
        fit$residuals^2,
        na.rm = TRUE
      )
    )
    
    # --------------------------------------------------------
    # A general fitted-value RMSE
    # --------------------------------------------------------
    
    fitted_RMSE <- sqrt(
      mean(
        (
          fit$fitted -
            mean(fit$fitted, na.rm = TRUE)
        )^2,
        na.rm = TRUE
      )
    )
    
    # --------------------------------------------------------
    # Store
    # --------------------------------------------------------
    
    dimension_results <-
      rbind(
        dimension_results,
        
        data.frame(
          p = p,
          replication = r,
          RMSE = fitted_RMSE,
          residual_RMSE =
            residual_RMSE,
          iterations =
            fit$iterations,
          convergence =
            fit$convergence,
          computation_time =
            computation_time,
          stringsAsFactors = FALSE
        )
      )
    
    # --------------------------------------------------------
    # Progress
    # --------------------------------------------------------
    
    if (
      r %% 100 == 0 ||
      r == R
    ) {
      
      cat(
        "Replication ",
        r,
        "/",
        R,
        "\n",
        sep = ""
      )
    }
  }
}

dimension_end_time <- Sys.time()

cat(
  "\nTotal dimensionality-study time:\n"
)

print(
  dimension_end_time -
    dimension_start_time
)

# ============================================================
# Summary
# ============================================================

dimension_summary <-
  aggregate(
    cbind(
      RMSE,
      residual_RMSE,
      iterations,
      convergence,
      computation_time
    ) ~ p,
    data = dimension_results,
    FUN = mean,
    na.rm = TRUE
  )

# ------------------------------------------------------------
# Display
# ------------------------------------------------------------

cat(
  "\n============================================\n"
)

cat(
  "DIMENSIONALITY SUMMARY\n"
)

cat(
  "============================================\n"
)

print(
  dimension_summary
)

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

write.csv(
  dimension_summary,
  file =
    "dimension_summary.csv",
  row.names = FALSE
)

# ============================================================
# Plot computation time
# ============================================================

if (
  nrow(dimension_results) > 0
) {
  
  p_dimension_time <- ggplot(
    dimension_results,
    aes(
      x = factor(p),
      y = computation_time
    )
  ) +
    
    geom_boxplot() +
    
    labs(
      title =
        "Computation Time versus Dimension",
      x =
        "Number of covariates",
      y =
        "Computation time (seconds)"
    ) +
    
    theme_bw()
  
  print(
    p_dimension_time
  )
  
  ggsave(
    filename =
      "Dimension_computation_time.png",
    plot =
      p_dimension_time,
    width = 8,
    height = 5,
    dpi = 300
  )
}

# ============================================================
# Plot residual RMSE
# ============================================================

if (
  nrow(dimension_results) > 0
) {
  
  p_dimension_rmse <- ggplot(
    dimension_results,
    aes(
      x = factor(p),
      y = residual_RMSE
    )
  ) +
    
    geom_boxplot() +
    
    labs(
      title =
        "Residual RMSE versus Dimension",
      x =
        "Number of covariates",
      y =
        "Residual RMSE"
    ) +
    
    theme_bw()
  
  print(
    p_dimension_rmse
  )
  
  ggsave(
    filename =
      "Dimension_residual_RMSE.png",
    plot =
      p_dimension_rmse,
    width = 8,
    height = 5,
    dpi = 300
  )
}

cat(
  "\nPart 10 completed successfully.\n"
)

################################################################################

# ============================================================
# FUNCTIONAL TRANSFORMATION REGRESSION MODEL
# PART 11: MISSPECIFIED MODEL EXPERIMENT
# ============================================================

# ------------------------------------------------------------
# Generate misspecified data
# ------------------------------------------------------------

generate_misspecified <- function(
    n,
    delta = 0.5,
    sigma = 0.5) {
  
  X1 <- runif(
    n,
    min = 0,
    max = 1
  )
  
  X2 <- runif(
    n,
    min = 0,
    max = 1
  )
  
  # True component functions
  f1 <- sin(
    2 * pi * X1
  )
  
  f2 <- X2^2
  
  # Correct structural component
  h0 <-
    f1 +
    f2 +
    0.5 * f1 * f2
  
  # Omitted component
  omitted <-
    delta * X1 * X2
  
  # Random error
  eps <-
    sigma * rnorm(n)
  
  # Latent response
  gY <-
    h0 +
    omitted +
    eps
  
  # Observed response
  Y <- exp(gY)
  
  data.frame(
    Y = Y,
    X1 = X1,
    X2 = X2
  )
}

# ------------------------------------------------------------
# Misspecification settings
# ------------------------------------------------------------

misspecification_values <- c(
  0,
  0.25,
  0.50,
  0.75
)

misspecification_n <- 400

# ------------------------------------------------------------
# Storage
# ------------------------------------------------------------

misspecification_results <- data.frame(
  delta = numeric(),
  replication = integer(),
  f1_RMSE = numeric(),
  f2_RMSE = numeric(),
  h_RMSE = numeric(),
  residual_RMSE = numeric(),
  iterations = numeric(),
  convergence = numeric(),
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# Simulation
# ------------------------------------------------------------

misspecification_start_time <-
  Sys.time()

for (delta in misspecification_values) {
  
  cat(
    "\n============================================\n"
  )
  
  cat(
    "Misspecification delta = ",
    delta,
    "\n",
    sep = ""
  )
  
  cat(
    "============================================\n"
  )
  
  for (r in seq_len(R)) {
    
    # --------------------------------------------------------
    # Generate data
    # --------------------------------------------------------
    
    dat <- generate_misspecified(
      n = misspecification_n,
      delta = delta,
      sigma = 0.5
    )
    
    # --------------------------------------------------------
    # Fit the assumed model
    # --------------------------------------------------------
    
    fit <- tryCatch(
      
      fit_functional_transformation(
        dat = dat,
        max_iter = 30,
        tol = 1e-4,
        bandwidth_multiplier = 1
      ),
      
      error = function(e) {
        
        message(
          "delta = ",
          delta,
          ", replication = ",
          r,
          ": ",
          e$message
        )
        
        NULL
      }
    )
    
    # --------------------------------------------------------
    # Skip failed fits
    # --------------------------------------------------------
    
    if (is.null(fit)) {
      next
    }
    
    # --------------------------------------------------------
    # Evaluate estimated functions
    # --------------------------------------------------------
    
    f1_hat <-
      fit$f1(grid_x)
    
    f2_hat <-
      fit$f2(grid_x)
    
    h_hat <-
      fit$h(
        grid_h$u,
        grid_h$v
      )
    
    # --------------------------------------------------------
    # Calculate function RMSEs
    # --------------------------------------------------------
    
    f1_RMSE <-
      sqrt(
        mean(
          (
            f1_hat -
              f1_true_grid
          )^2,
          na.rm = TRUE
        )
      )
    
    f2_RMSE <-
      sqrt(
        mean(
          (
            f2_hat -
              f2_true_grid
          )^2,
          na.rm = TRUE
        )
      )
    
    h_RMSE <-
      sqrt(
        mean(
          (
            h_hat -
              h_true_grid
          )^2,
          na.rm = TRUE
        )
      )
    
    # --------------------------------------------------------
    # Residual RMSE
    # --------------------------------------------------------
    
    residual_RMSE <-
      sqrt(
        mean(
          fit$residuals^2,
          na.rm = TRUE
        )
      )
    
    # --------------------------------------------------------
    # Store results
    # --------------------------------------------------------
    
    misspecification_results <-
      rbind(
        misspecification_results,
        
        data.frame(
          delta = delta,
          replication = r,
          f1_RMSE = f1_RMSE,
          f2_RMSE = f2_RMSE,
          h_RMSE = h_RMSE,
          residual_RMSE =
            residual_RMSE,
          iterations =
            fit$iterations,
          convergence =
            fit$convergence,
          stringsAsFactors = FALSE
        )
      )
    
    # --------------------------------------------------------
    # Progress
    # --------------------------------------------------------
    
    if (
      r %% 100 == 0 ||
      r == R
    ) {
      
      cat(
        "Replication ",
        r,
        "/",
        R,
        "\n",
        sep = ""
      )
    }
  }
}

misspecification_end_time <-
  Sys.time()

cat(
  "\nTotal misspecification-study time:\n"
)

print(
  misspecification_end_time -
    misspecification_start_time
)

# ============================================================
# Summary
# ============================================================

misspecification_summary <-
  aggregate(
    cbind(
      f1_RMSE,
      f2_RMSE,
      h_RMSE,
      residual_RMSE,
      iterations,
      convergence
    ) ~ delta,
    data =
      misspecification_results,
    FUN = mean,
    na.rm = TRUE
  )

# ------------------------------------------------------------
# Display summary
# ------------------------------------------------------------

cat(
  "\n============================================\n"
)

cat(
  "MISSPECIFICATION SUMMARY\n"
)

cat(
  "============================================\n"
)

print(
  misspecification_summary
)

# ------------------------------------------------------------
# Save summary
# ------------------------------------------------------------

write.csv(
  misspecification_summary,
  file =
    "misspecification_summary.csv",
  row.names = FALSE
)

# ============================================================
# Convert results to long format for plotting
# ============================================================

misspecification_plot_data <-
  rbind(
    
    data.frame(
      delta =
        misspecification_results$delta,
      replication =
        misspecification_results$replication,
      component = "f1",
      RMSE =
        misspecification_results$f1_RMSE
    ),
    
    data.frame(
      delta =
        misspecification_results$delta,
      replication =
        misspecification_results$replication,
      component = "f2",
      RMSE =
        misspecification_results$f2_RMSE
    ),
    
    data.frame(
      delta =
        misspecification_results$delta,
      replication =
        misspecification_results$replication,
      component = "h",
      RMSE =
        misspecification_results$h_RMSE
    )
  )

# ============================================================
# Plot RMSE under misspecification
# ============================================================

if (
  nrow(misspecification_plot_data) > 0
) {
  
  p_misspecification <-
    ggplot(
      misspecification_plot_data,
      aes(
        x = factor(delta),
        y = RMSE,
        linetype = component
      )
    ) +
    
    geom_boxplot(
      aes(
        group =
          interaction(
            delta,
            component
          )
      )
    ) +
    
    labs(
      title =
        "Effect of Model Misspecification",
      x =
        "Misspecification parameter",
      y =
        "RMSE",
      linetype =
        "Component"
    ) +
    
    theme_bw()
  
  print(
    p_misspecification
  )
  
  ggsave(
    filename =
      "Misspecification_RMSE.png",
    plot =
      p_misspecification,
    width = 8,
    height = 5,
    dpi = 300
  )
}

# ============================================================
# Residual RMSE versus misspecification
# ============================================================

if (
  nrow(misspecification_results) > 0
) {
  
  p_residual_misspecification <-
    ggplot(
      misspecification_results,
      aes(
        x = factor(delta),
        y = residual_RMSE
      )
    ) +
    
    geom_boxplot() +
    
    labs(
      title =
        "Residual RMSE under Model Misspecification",
      x =
        "Misspecification parameter",
      y =
        "Residual RMSE"
    ) +
    
    theme_bw()
  
  print(
    p_residual_misspecification
  )
  
  ggsave(
    filename =
      "Misspecification_Residual_RMSE.png",
    plot =
      p_residual_misspecification,
    width = 8,
    height = 5,
    dpi = 300
  )
}

cat(
  "\nPart 11 completed successfully.\n"
)

################################################################################

# ============================================================
# FUNCTIONAL TRANSFORMATION REGRESSION MODEL
# PART 12: KENDALL PERMUTATION INDEPENDENCE TEST
# ============================================================

# ------------------------------------------------------------
# Kendall permutation test
# ------------------------------------------------------------

kendall_permutation_test <- function(
    x,
    y,
    B = 500) {
  
  # Remove missing and non-finite observations
  keep <- is.finite(x) &
    is.finite(y)
  
  x <- x[keep]
  y <- y[keep]
  
  n <- length(x)
  
  # Check sample size
  if (n < 3) {
    
    return(
      list(
        statistic = NA_real_,
        p.value = NA_real_
      )
    )
  }
  
  # Observed Kendall's tau
  observed <-
    suppressWarnings(
      cor(
        x,
        y,
        method = "kendall"
      )
    )
  
  # Permutation statistics
  permuted <- numeric(B)
  
  for (b in seq_len(B)) {
    
    permuted[b] <-
      suppressWarnings(
        cor(
          sample(x),
          y,
          method = "kendall"
        )
      )
  }
  
  # Two-sided permutation p-value
  p_value <-
    (
      1 +
        sum(
          abs(permuted) >=
            abs(observed),
          na.rm = TRUE
        )
    ) /
    (B + 1)
  
  list(
    statistic = observed,
    p.value = p_value,
    permutations = permuted
  )
}

# ============================================================
# Simulation settings
# ============================================================

kendall_n <- 400

# ------------------------------------------------------------
# Storage
# ------------------------------------------------------------

kendall_results <- data.frame(
  replication = integer(),
  statistic = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# ============================================================
# Kendall test simulation
# ============================================================

kendall_start_time <- Sys.time()

for (r in seq_len(R)) {
  
  # ----------------------------------------------------------
  # Generate data
  # ----------------------------------------------------------
  
  dat <- generate_baseline(
    n = kendall_n,
    sigma = 0.5,
    response_transformation = "log",
    error_distribution = "normal"
  )
  
  # ----------------------------------------------------------
  # Fit functional transformation model
  # ----------------------------------------------------------
  
  fit <- tryCatch(
    
    fit_functional_transformation(
      dat = dat,
      max_iter = 30,
      tol = 1e-4,
      bandwidth_multiplier = 1
    ),
    
    error = function(e) {
      
      message(
        "Replication ",
        r,
        ": ",
        e$message
      )
      
      NULL
    }
  )
  
  if (is.null(fit)) {
    next
  }
  
  # ----------------------------------------------------------
  # Estimate the nonparametric component
  # ----------------------------------------------------------
  
  h_hat <-
    fit$h(
      dat$X1,
      dat$X2
    )
  
  # ----------------------------------------------------------
  # Independence test between h_hat and residuals
  # ----------------------------------------------------------
  
  test_result <-
    kendall_permutation_test(
      x = h_hat,
      y = fit$residuals,
      B = B_perm
    )
  
  # ----------------------------------------------------------
  # Store result
  # ----------------------------------------------------------
  
  kendall_results <-
    rbind(
      kendall_results,
      
      data.frame(
        replication = r,
        statistic =
          test_result$statistic,
        p_value =
          test_result$p.value,
        stringsAsFactors = FALSE
      )
    )
  
  # ----------------------------------------------------------
  # Progress
  # ----------------------------------------------------------
  
  if (
    r %% 100 == 0 ||
    r == R
  ) {
    
    cat(
      "Replication ",
      r,
      "/",
      R,
      "\n",
      sep = ""
    )
  }
}

kendall_end_time <- Sys.time()

cat(
  "\nTotal Kendall permutation-test time:\n"
)

print(
  kendall_end_time -
    kendall_start_time
)

# ============================================================
# Summary
# ============================================================

if (
  nrow(kendall_results) > 0
) {
  
  kendall_summary <- data.frame(
    
    sample_size =
      kendall_n,
    
    replications =
      nrow(kendall_results),
    
    mean_tau =
      mean(
        kendall_results$statistic,
        na.rm = TRUE
      ),
    
    sd_tau =
      sd(
        kendall_results$statistic,
        na.rm = TRUE
      ),
    
    mean_p_value =
      mean(
        kendall_results$p_value,
        na.rm = TRUE
      ),
    
    rejection_rate =
      mean(
        kendall_results$p_value <
          alpha,
        na.rm = TRUE
      )
  )
  
} else {
  
  kendall_summary <- data.frame(
    
    sample_size =
      kendall_n,
    
    replications = 0,
    
    mean_tau = NA_real_,
    
    sd_tau = NA_real_,
    
    mean_p_value = NA_real_,
    
    rejection_rate = NA_real_
  )
}

# ------------------------------------------------------------
# Display
# ------------------------------------------------------------

cat(
  "\n============================================\n"
)

cat(
  "KENDALL PERMUTATION TEST SUMMARY\n"
)

cat(
  "============================================\n"
)

print(
  kendall_summary
)

# ------------------------------------------------------------
# Save summary
# ------------------------------------------------------------

write.csv(
  kendall_summary,
  file =
    "kendall_permutation_summary.csv",
  row.names = FALSE
)

write.csv(
  kendall_results,
  file =
    "kendall_permutation_results.csv",
  row.names = FALSE
)

# ============================================================
# Distribution of Kendall's tau
# ============================================================

if (
  nrow(kendall_results) > 0
) {
  
  p_kendall_tau <-
    ggplot(
      kendall_results,
      aes(
        x = statistic
      )
    ) +
    
    geom_histogram(
      bins = 30
    ) +
    
    labs(
      title =
        "Distribution of Kendall's Tau",
      x =
        "Kendall's tau",
      y =
        "Frequency"
    ) +
    
    theme_bw()
  
  print(
    p_kendall_tau
  )
  
  ggsave(
    filename =
      "Kendall_tau_distribution.png",
    plot =
      p_kendall_tau,
    width = 8,
    height = 5,
    dpi = 300
  )
}

# ============================================================
# Distribution of permutation p-values
# ============================================================

if (
  nrow(kendall_results) > 0
) {
  
  p_kendall_pvalues <-
    ggplot(
      kendall_results,
      aes(
        x = p_value
      )
    ) +
    
    geom_histogram(
      bins = 20
    ) +
    
    geom_vline(
      xintercept = alpha,
      linetype = "dashed"
    ) +
    
    labs(
      title =
        "Distribution of Kendall Permutation p-values",
      x =
        "Permutation p-value",
      y =
        "Frequency"
    ) +
    
    theme_bw()
  
  print(
    p_kendall_pvalues
  )
  
  ggsave(
    filename =
      "Kendall_permutation_pvalues.png",
    plot =
      p_kendall_pvalues,
    width = 8,
    height = 5,
    dpi = 300
  )
}

cat(
  "\nPart 12 completed successfully.\n"
)

################################################################################

# ============================================================
# FUNCTIONAL TRANSFORMATION REGRESSION MODEL
# PART 13: TAU-STAR DIAGNOSTIC / PERMUTATION PROCEDURE
# ============================================================

# ------------------------------------------------------------
# Diagnostic tau-star statistic
# ------------------------------------------------------------

tau_star_diagnostic <- function(
    x,
    y) {
  
  # Remove non-finite observations
  keep <- is.finite(x) &
    is.finite(y)
  
  x <- x[keep]
  y <- y[keep]
  
  n <- length(x)
  
  if (n < 4) {
    
    return(NA_real_)
  }
  
  # Pairwise ordering matrices
  dx <- outer(
    x,
    x,
    FUN = "-"
  )
  
  dy <- outer(
    y,
    y,
    FUN = "-"
  )
  
  # Sign matrices
  sx <- sign(dx)
  sy <- sign(dy)
  
  # Diagnostic concordance-discordance measure
  product_matrix <-
    sx * sy
  
  upper_index <-
    upper.tri(
      product_matrix
    )
  
  values <-
    product_matrix[
      upper_index
    ]
  
  if (length(values) == 0) {
    
    return(NA_real_)
  }
  
  mean(
    values,
    na.rm = TRUE
  )
}

# ------------------------------------------------------------
# Permutation test based on the diagnostic statistic
# ------------------------------------------------------------

tau_star_permutation_test <- function(
    x,
    y,
    B = 500) {
  
  keep <- is.finite(x) &
    is.finite(y)
  
  x <- x[keep]
  y <- y[keep]
  
  if (length(x) < 4) {
    
    return(
      list(
        statistic = NA_real_,
        p.value = NA_real_
      )
    )
  }
  
  observed <-
    tau_star_diagnostic(
      x,
      y
    )
  
  permuted <-
    numeric(B)
  
  for (b in seq_len(B)) {
    
    permuted[b] <-
      tau_star_diagnostic(
        x,
        sample(y)
      )
  }
  
  p_value <-
    (
      1 +
        sum(
          abs(permuted) >=
            abs(observed),
          na.rm = TRUE
        )
    ) /
    (B + 1)
  
  list(
    statistic = observed,
    p.value = p_value,
    permutations = permuted
  )
}

# ============================================================
# Settings
# ============================================================

tau_star_n <- 200

tau_star_deltas <- c(
  0,
  0.25,
  0.50,
  0.75
)

# ------------------------------------------------------------
# Storage
# ------------------------------------------------------------

tau_star_results <- data.frame(
  delta = numeric(),
  replication = integer(),
  statistic = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# ============================================================
# Simulation
# ============================================================

tau_star_start_time <-
  Sys.time()

for (delta in tau_star_deltas) {
  
  cat(
    "\n============================================\n"
  )
  
  cat(
    "Tau-star diagnostic: delta = ",
    delta,
    "\n",
    sep = ""
  )
  
  cat(
    "============================================\n"
  )
  
  for (r in seq_len(R)) {
    
    # --------------------------------------------------------
    # Generate data
    # --------------------------------------------------------
    
    dat <-
      generate_misspecified(
        n = tau_star_n,
        delta = delta,
        sigma = 0.5
      )
    
    # --------------------------------------------------------
    # Fit model
    # --------------------------------------------------------
    
    fit <- tryCatch(
      
      fit_functional_transformation(
        dat = dat,
        max_iter = 30,
        tol = 1e-4,
        bandwidth_multiplier = 1
      ),
      
      error = function(e) {
        
        message(
          "delta = ",
          delta,
          ", replication = ",
          r,
          ": ",
          e$message
        )
        
        NULL
      }
    )
    
    if (is.null(fit)) {
      next
    }
    
    # --------------------------------------------------------
    # Estimated nonparametric component
    # --------------------------------------------------------
    
    h_hat <-
      fit$h(
        dat$X1,
        dat$X2
      )
    
    # --------------------------------------------------------
    # Tau-star diagnostic
    # --------------------------------------------------------
    
    test_result <-
      tau_star_permutation_test(
        x = h_hat,
        y = fit$residuals,
        B = B_perm
      )
    
    # --------------------------------------------------------
    # Store results
    # --------------------------------------------------------
    
    tau_star_results <-
      rbind(
        tau_star_results,
        
        data.frame(
          delta = delta,
          replication = r,
          statistic =
            test_result$statistic,
          p_value =
            test_result$p.value,
          stringsAsFactors = FALSE
        )
      )
    
    # --------------------------------------------------------
    # Progress
    # --------------------------------------------------------
    
    if (
      r %% 50 == 0 ||
      r == R
    ) {
      
      cat(
        "Replication ",
        r,
        "/",
        R,
        "\n",
        sep = ""
      )
    }
  }
}

tau_star_end_time <-
  Sys.time()

cat(
  "\nTotal tau-star diagnostic time:\n"
)

print(
  tau_star_end_time -
    tau_star_start_time
)

# ============================================================
# Summary
# ============================================================

if (
  nrow(tau_star_results) > 0
) {
  
  tau_star_summary <-
    aggregate(
      cbind(
        statistic,
        p_value
      ) ~ delta,
      data =
        tau_star_results,
      FUN = mean,
      na.rm = TRUE
    )
  
  rejection_rates <-
    aggregate(
      (
        p_value < alpha
      ) ~ delta,
      data =
        tau_star_results,
      FUN = mean,
      na.rm = TRUE
    )
  
  names(
    rejection_rates
  )[2] <-
    "rejection_rate"
  
  tau_star_summary <-
    merge(
      tau_star_summary,
      rejection_rates,
      by = "delta"
    )
  
} else {
  
  tau_star_summary <- data.frame(
    delta = tau_star_deltas,
    statistic = NA_real_,
    p_value = NA_real_,
    rejection_rate = NA_real_
  )
}

# ------------------------------------------------------------
# Display
# ------------------------------------------------------------

cat(
  "\n============================================\n"
)

cat(
  "TAU-STAR DIAGNOSTIC SUMMARY\n"
)

cat(
  "============================================\n"
)

print(
  tau_star_summary
)

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

write.csv(
  tau_star_summary,
  file =
    "tau_star_diagnostic_summary.csv",
  row.names = FALSE
)

write.csv(
  tau_star_results,
  file =
    "tau_star_diagnostic_results.csv",
  row.names = FALSE
)

# ============================================================
# Plot statistic against misspecification
# ============================================================

if (
  nrow(tau_star_results) > 0
) {
  
  p_tau_star <-
    ggplot(
      tau_star_results,
      aes(
        x = factor(delta),
        y = statistic
      )
    ) +
    
    geom_boxplot() +
    
    labs(
      title =
        "Tau-Star Diagnostic under Misspecification",
      x =
        "Misspecification parameter",
      y =
        "Diagnostic statistic"
    ) +
    
    theme_bw()
  
  print(
    p_tau_star
  )
  
  ggsave(
    filename =
      "Tau_star_diagnostic.png",
    plot =
      p_tau_star,
    width = 8,
    height = 5,
    dpi = 300
  )
}

# ============================================================
# Plot p-values
# ============================================================

if (
  nrow(tau_star_results) > 0
) {
  
  p_tau_star_pvalues <-
    ggplot(
      tau_star_results,
      aes(
        x = factor(delta),
        y = p_value
      )
    ) +
    
    geom_boxplot() +
    
    geom_hline(
      yintercept = alpha,
      linetype = "dashed"
    ) +
    
    labs(
      title =
        "Tau-Star Diagnostic p-values",
      x =
        "Misspecification parameter",
      y =
        "Permutation p-value"
    ) +
    
    theme_bw()
  
  print(
    p_tau_star_pvalues
  )
  
  ggsave(
    filename =
      "Tau_star_diagnostic_pvalues.png",
    plot =
      p_tau_star_pvalues,
    width = 8,
    height = 5,
    dpi = 300
  )
}

cat(
  "\nPart 13 completed successfully.\n"
)

#################################################################################

# ============================================================
# FUNCTIONAL TRANSFORMATION REGRESSION MODEL
# PART 14: DEPENDENCE-DIAGNOSTIC SIMULATION
# ============================================================

# ------------------------------------------------------------
# Settings
# ------------------------------------------------------------

diagnostic_n <- 200

diagnostic_deltas <- c(
  0,
  0.25,
  0.50,
  0.75
)

# ------------------------------------------------------------
# Storage
# ------------------------------------------------------------

dependence_diagnostic_results <- data.frame(
  delta = numeric(),
  replication = integer(),
  kendall_statistic = numeric(),
  kendall_p_value = numeric(),
  tau_star_statistic = numeric(),
  tau_star_p_value = numeric(),
  stringsAsFactors = FALSE
)

# ============================================================
# Simulation
# ============================================================

diagnostic_start_time <-
  Sys.time()

for (delta in diagnostic_deltas) {
  
  cat(
    "\n============================================\n"
  )
  
  cat(
    "Dependence diagnostic: delta = ",
    delta,
    "\n",
    sep = ""
  )
  
  cat(
    "============================================\n"
  )
  
  for (r in seq_len(R)) {
    
    # --------------------------------------------------------
    # Generate data
    # --------------------------------------------------------
    
    dat <-
      generate_misspecified(
        n = diagnostic_n,
        delta = delta,
        sigma = 0.5
      )
    
    # --------------------------------------------------------
    # Fit model
    # --------------------------------------------------------
    
    fit <- tryCatch(
      
      fit_functional_transformation(
        dat = dat,
        max_iter = 30,
        tol = 1e-4,
        bandwidth_multiplier = 1
      ),
      
      error = function(e) {
        
        message(
          "delta = ",
          delta,
          ", replication = ",
          r,
          ": ",
          e$message
        )
        
        NULL
      }
    )
    
    if (is.null(fit)) {
      next
    }
    
    # --------------------------------------------------------
    # Estimated nonparametric component
    # --------------------------------------------------------
    
    h_hat <-
      fit$h(
        dat$X1,
        dat$X2
      )
    
    # --------------------------------------------------------
    # Kendall permutation test
    # --------------------------------------------------------
    
    kendall_result <-
      kendall_permutation_test(
        x = h_hat,
        y = fit$residuals,
        B = B_perm
      )
    
    # --------------------------------------------------------
    # Tau-star diagnostic
    # --------------------------------------------------------
    
    tau_star_result <-
      tau_star_permutation_test(
        x = h_hat,
        y = fit$residuals,
        B = B_perm
      )
    
    # --------------------------------------------------------
    # Store
    # --------------------------------------------------------
    
    dependence_diagnostic_results <-
      rbind(
        dependence_diagnostic_results,
        
        data.frame(
          delta = delta,
          replication = r,
          
          kendall_statistic =
            kendall_result$statistic,
          
          kendall_p_value =
            kendall_result$p.value,
          
          tau_star_statistic =
            tau_star_result$statistic,
          
          tau_star_p_value =
            tau_star_result$p.value,
          
          stringsAsFactors = FALSE
        )
      )
    
    # --------------------------------------------------------
    # Progress
    # --------------------------------------------------------
    
    if (
      r %% 50 == 0 ||
      r == R
    ) {
      
      cat(
        "Replication ",
        r,
        "/",
        R,
        "\n",
        sep = ""
      )
    }
  }
}

diagnostic_end_time <-
  Sys.time()

cat(
  "\nTotal dependence-diagnostic time:\n"
)

print(
  diagnostic_end_time -
    diagnostic_start_time
)

# ============================================================
# Calculate rejection indicators
# ============================================================

if (
  nrow(
    dependence_diagnostic_results
  ) > 0
) {
  
  dependence_diagnostic_results$kendall_reject <-
    dependence_diagnostic_results$kendall_p_value <
    alpha
  
  dependence_diagnostic_results$tau_star_reject <-
    dependence_diagnostic_results$tau_star_p_value <
    alpha
}

# ============================================================
# Summary
# ============================================================

if (
  nrow(
    dependence_diagnostic_results
  ) > 0
) {
  
  dependence_diagnostic_summary <-
    aggregate(
      cbind(
        kendall_statistic,
        kendall_p_value,
        tau_star_statistic,
        tau_star_p_value,
        kendall_reject,
        tau_star_reject
      ) ~ delta,
      data =
        dependence_diagnostic_results,
      FUN = mean,
      na.rm = TRUE
    )
  
  names(
    dependence_diagnostic_summary
  )[names(
    dependence_diagnostic_summary
  ) == "kendall_reject"] <-
    "kendall_rejection_rate"
  
  names(
    dependence_diagnostic_summary
  )[names(
    dependence_diagnostic_summary
  ) == "tau_star_reject"] <-
    "tau_star_rejection_rate"
  
} else {
  
  dependence_diagnostic_summary <-
    data.frame(
      delta = diagnostic_deltas,
      kendall_statistic = NA_real_,
      kendall_p_value = NA_real_,
      tau_star_statistic = NA_real_,
      tau_star_p_value = NA_real_,
      kendall_rejection_rate = NA_real_,
      tau_star_rejection_rate = NA_real_
    )
}

# ============================================================
# Display summary
# ============================================================

cat(
  "\n============================================\n"
)

cat(
  "DEPENDENCE-DIAGNOSTIC SUMMARY\n"
)

cat(
  "============================================\n"
)

print(
  dependence_diagnostic_summary
)

# ============================================================
# Save results
# ============================================================

write.csv(
  dependence_diagnostic_results,
  file =
    "dependence_diagnostic_results.csv",
  row.names = FALSE
)

write.csv(
  dependence_diagnostic_summary,
  file =
    "dependence_diagnostic_summary.csv",
  row.names = FALSE
)

# ============================================================
# Prepare rejection-rate data
# ============================================================

rejection_rate_data <-
  rbind(
    
    data.frame(
      delta =
        dependence_diagnostic_summary$delta,
      
      test = "Kendall",
      
      rejection_rate =
        dependence_diagnostic_summary$
        kendall_rejection_rate
    ),
    
    data.frame(
      delta =
        dependence_diagnostic_summary$delta,
      
      test = "Tau-star diagnostic",
      
      rejection_rate =
        dependence_diagnostic_summary$
        tau_star_rejection_rate
    )
  )

# ============================================================
# Plot rejection rates
# ============================================================

if (
  nrow(rejection_rate_data) > 0
) {
  
  p_rejection_rate <-
    ggplot(
      rejection_rate_data,
      aes(
        x = factor(delta),
        y = rejection_rate,
        linetype = test,
        group = test
      )
    ) +
    
    geom_point() +
    
    geom_line() +
    
    geom_hline(
      yintercept = alpha,
      linetype = "dashed"
    ) +
    
    labs(
      title =
        "Rejection Rates under Increasing Misspecification",
      x =
        "Misspecification parameter",
      y =
        "Rejection rate",
      linetype =
        "Test"
    ) +
    
    theme_bw()
  
  print(
    p_rejection_rate
  )
  
  ggsave(
    filename =
      "Dependence_diagnostic_rejection_rates.png",
    plot =
      p_rejection_rate,
    width = 8,
    height = 5,
    dpi = 300
  )
}

# ============================================================
# Plot p-values
# ============================================================

if (
  nrow(
    dependence_diagnostic_results
  ) > 0
) {
  
  p_diagnostic_pvalues <-
    ggplot(
      dependence_diagnostic_results
    ) +
    
    geom_boxplot(
      aes(
        x = factor(delta),
        y = kendall_p_value,
        group = factor(delta)
      )
    ) +
    
    geom_hline(
      yintercept = alpha,
      linetype = "dashed"
    ) +
    
    labs(
      title =
        "Kendall Permutation p-values",
      x =
        "Misspecification parameter",
      y =
        "p-value"
    ) +
    
    theme_bw()
  
  print(
    p_diagnostic_pvalues
  )
  
  ggsave(
    filename =
      "Dependence_diagnostic_Kendall_pvalues.png",
    plot =
      p_diagnostic_pvalues,
    width = 8,
    height = 5,
    dpi = 300
  )
}

cat(
  "\nPart 14 completed successfully.\n"
)

#################################################################################

# ============================================================
# FUNCTIONAL TRANSFORMATION REGRESSION MODEL
# PART 15: BANDWIDTH SENSITIVITY ANALYSIS
# ============================================================

# ------------------------------------------------------------
# Bandwidth settings
# ------------------------------------------------------------

bandwidth_multipliers <- c(
  0.75,
  1.00,
  1.25
)

bandwidth_n <- 400

# ------------------------------------------------------------
# Storage
# ------------------------------------------------------------

bandwidth_rmse <- data.frame(
  bandwidth_multiplier = numeric(),
  replication = integer(),
  component = character(),
  RMSE = numeric(),
  stringsAsFactors = FALSE
)

# ============================================================
# Bandwidth sensitivity simulation
# ============================================================

bandwidth_start_time <-
  Sys.time()

for (
  multiplier in bandwidth_multipliers
) {
  
  cat(
    "\n============================================\n"
  )
  
  cat(
    "Bandwidth multiplier = ",
    multiplier,
    "\n",
    sep = ""
  )
  
  cat(
    "============================================\n"
  )
  
  for (r in seq_len(R)) {
    
    # --------------------------------------------------------
    # Generate data
    # --------------------------------------------------------
    
    dat <-
      generate_baseline(
        n = bandwidth_n,
        sigma = 0.5,
        response_transformation = "log",
        error_distribution = "normal"
      )
    
    # --------------------------------------------------------
    # Fit model
    # --------------------------------------------------------
    
    fit <- tryCatch(
      
      fit_functional_transformation(
        dat = dat,
        max_iter = 30,
        tol = 1e-4,
        bandwidth_multiplier =
          multiplier
      ),
      
      error = function(e) {
        
        message(
          "Bandwidth = ",
          multiplier,
          ", replication = ",
          r,
          ": ",
          e$message
        )
        
        NULL
      }
    )
    
    if (is.null(fit)) {
      next
    }
    
    # --------------------------------------------------------
    # Evaluate estimated functions
    # --------------------------------------------------------
    
    f1_hat <-
      fit$f1(
        grid_x
      )
    
    f2_hat <-
      fit$f2(
        grid_x
      )
    
    h_hat <-
      fit$h(
        grid_h$u,
        grid_h$v
      )
    
    # --------------------------------------------------------
    # RMSE for f1
    # --------------------------------------------------------
    
    rmse_f1 <-
      sqrt(
        mean(
          (
            f1_hat -
              f1_true_grid
          )^2,
          na.rm = TRUE
        )
      )
    
    # --------------------------------------------------------
    # RMSE for f2
    # --------------------------------------------------------
    
    rmse_f2 <-
      sqrt(
        mean(
          (
            f2_hat -
              f2_true_grid
          )^2,
          na.rm = TRUE
        )
      )
    
    # --------------------------------------------------------
    # RMSE for h
    # --------------------------------------------------------
    
    rmse_h <-
      sqrt(
        mean(
          (
            h_hat -
              h_true_grid
          )^2,
          na.rm = TRUE
        )
      )
    
    # --------------------------------------------------------
    # Store results
    # --------------------------------------------------------
    
    bandwidth_rmse <-
      rbind(
        
        bandwidth_rmse,
        
        data.frame(
          bandwidth_multiplier =
            multiplier,
          replication = r,
          component = "f1",
          RMSE = rmse_f1,
          stringsAsFactors = FALSE
        ),
        
        data.frame(
          bandwidth_multiplier =
            multiplier,
          replication = r,
          component = "f2",
          RMSE = rmse_f2,
          stringsAsFactors = FALSE
        ),
        
        data.frame(
          bandwidth_multiplier =
            multiplier,
          replication = r,
          component = "h",
          RMSE = rmse_h,
          stringsAsFactors = FALSE
        )
      )
    
    # --------------------------------------------------------
    # Progress
    # --------------------------------------------------------
    
    if (
      r %% 100 == 0 ||
      r == R
    ) {
      
      cat(
        "Replication ",
        r,
        "/",
        R,
        "\n",
        sep = ""
      )
    }
  }
}

bandwidth_end_time <-
  Sys.time()

cat(
  "\nTotal bandwidth-sensitivity time:\n"
)

print(
  bandwidth_end_time -
    bandwidth_start_time
)

# ============================================================
# Summary
# ============================================================

bandwidth_summary <-
  aggregate(
    RMSE ~
      bandwidth_multiplier +
      component,
    data =
      bandwidth_rmse,
    FUN = mean,
    na.rm = TRUE
  )

# ------------------------------------------------------------
# Display
# ------------------------------------------------------------

cat(
  "\n============================================\n"
)

cat(
  "BANDWIDTH SENSITIVITY SUMMARY\n"
)

cat(
  "============================================\n"
)

print(
  bandwidth_summary
)

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

write.csv(
  bandwidth_summary,
  file =
    "bandwidth_summary.csv",
  row.names = FALSE
)

write.csv(
  bandwidth_rmse,
  file =
    "bandwidth_rmse_results.csv",
  row.names = FALSE
)

# ============================================================
# Boxplot of RMSE
# ============================================================

if (
  nrow(bandwidth_rmse) > 0
) {
  
  p_bandwidth <-
    ggplot(
      bandwidth_rmse,
      aes(
        x =
          factor(
            bandwidth_multiplier
          ),
        y = RMSE,
        group =
          interaction(
            bandwidth_multiplier,
            component
          ),
        linetype = component
      )
    ) +
    
    geom_boxplot() +
    
    labs(
      title =
        "Bandwidth Sensitivity",
      x =
        "Bandwidth multiplier",
      y =
        "RMSE",
      linetype =
        "Function"
    ) +
    
    theme_bw()
  
  print(
    p_bandwidth
  )
  
  ggsave(
    filename =
      "Bandwidth_sensitivity.png",
    plot =
      p_bandwidth,
    width = 8,
    height = 5,
    dpi = 300
  )
}

# ============================================================
# Mean RMSE plot
# ============================================================

if (
  nrow(bandwidth_summary) > 0
) {
  
  p_bandwidth_mean <-
    ggplot(
      bandwidth_summary,
      aes(
        x =
          bandwidth_multiplier,
        y = RMSE,
        linetype = component,
        group = component
      )
    ) +
    
    geom_point() +
    
    geom_line() +
    
    labs(
      title =
        "Mean RMSE versus Bandwidth Multiplier",
      x =
        "Bandwidth multiplier",
      y =
        "Mean RMSE",
      linetype =
        "Function"
    ) +
    
    theme_bw()
  
  print(
    p_bandwidth_mean
  )
  
  ggsave(
    filename =
      "Bandwidth_mean_RMSE.png",
    plot =
      p_bandwidth_mean,
    width = 8,
    height = 5,
    dpi = 300
  )
}

cat(
  "\nPart 15 completed successfully.\n"
)

###################################################################################

# ============================================================
# FUNCTIONAL TRANSFORMATION REGRESSION MODEL
# PART 16: FINAL SUMMARIES, DIAGNOSTICS, AND EXPORT
# ============================================================

# ============================================================
# 1. Check available simulation objects
# ============================================================

cat(
  "\n============================================\n"
)

cat(
  "FINAL SIMULATION SUMMARY\n"
)

cat(
  "============================================\n\n"
)

cat(
  "Objects available in the workspace:\n\n"
)

objects_to_check <- c(
  "finite_sample_summary",
  "asymptotic_normality_summary",
  "response_transformation_summary",
  "error_distribution_summary",
  "dimension_summary",
  "misspecification_summary",
  "kendall_summary",
  "tau_star_summary",
  "dependence_diagnostic_summary",
  "bandwidth_summary"
)

for (object_name in objects_to_check) {
  
  object_exists <-
    exists(
      object_name,
      inherits = TRUE
    )
  
  cat(
    object_name,
    ": ",
    object_exists,
    "\n",
    sep = ""
  )
}

# ============================================================
# 2. Create final summary directory
# ============================================================

final_directory <-
  "final_simulation_results"

if (
  !dir.exists(
    final_directory
  )
) {
  
  dir.create(
    final_directory
  )
}

# ============================================================
# 3. Export all available summary objects
# ============================================================

if (
  exists(
    "finite_sample_summary"
  )
) {
  
  write.csv(
    finite_sample_summary,
    file =
      file.path(
        final_directory,
        "finite_sample_summary.csv"
      ),
    row.names = FALSE
  )
}

if (
  exists(
    "asymptotic_normality_summary"
  )
) {
  
  write.csv(
    asymptotic_normality_summary,
    file =
      file.path(
        final_directory,
        "asymptotic_normality_summary.csv"
      ),
    row.names = FALSE
  )
}

if (
  exists(
    "response_transformation_summary"
  )
) {
  
  write.csv(
    response_transformation_summary,
    file =
      file.path(
        final_directory,
        "response_transformation_summary.csv"
      ),
    row.names = FALSE
  )
}

if (
  exists(
    "error_distribution_summary"
  )
) {
  
  write.csv(
    error_distribution_summary,
    file =
      file.path(
        final_directory,
        "error_distribution_summary.csv"
      ),
    row.names = FALSE
  )
}

if (
  exists(
    "dimension_summary"
  )
) {
  
  write.csv(
    dimension_summary,
    file =
      file.path(
        final_directory,
        "dimension_summary.csv"
      ),
    row.names = FALSE
  )
}

if (
  exists(
    "misspecification_summary"
  )
) {
  
  write.csv(
    misspecification_summary,
    file =
      file.path(
        final_directory,
        "misspecification_summary.csv"
      ),
    row.names = FALSE
  )
}

if (
  exists(
    "kendall_summary"
  )
) {
  
  write.csv(
    kendall_summary,
    file =
      file.path(
        final_directory,
        "kendall_summary.csv"
      ),
    row.names = FALSE
  )
}

if (
  exists(
    "tau_star_summary"
  )
) {
  
  write.csv(
    tau_star_summary,
    file =
      file.path(
        final_directory,
        "tau_star_summary.csv"
      ),
    row.names = FALSE
  )
}

if (
  exists(
    "dependence_diagnostic_summary"
  )
) {
  
  write.csv(
    dependence_diagnostic_summary,
    file =
      file.path(
        final_directory,
        "dependence_diagnostic_summary.csv"
      ),
    row.names = FALSE
  )
}

if (
  exists(
    "bandwidth_summary"
  )
) {
  
  write.csv(
    bandwidth_summary,
    file =
      file.path(
        final_directory,
        "bandwidth_summary.csv"
      ),
    row.names = FALSE
  )
}

# ============================================================
# 4. Create compact final summary
# ============================================================

final_summary <- data.frame(
  Analysis = character(),
  Status = character(),
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# Finite-sample study
# ------------------------------------------------------------

if (
  exists(
    "finite_sample_summary"
  )
) {
  
  final_summary <-
    rbind(
      final_summary,
      data.frame(
        Analysis =
          "Finite-sample estimation",
        Status =
          "Completed",
        stringsAsFactors = FALSE
      )
    )
}

# ------------------------------------------------------------
# Asymptotic normality
# ------------------------------------------------------------

if (
  exists(
    "asymptotic_normality_summary"
  )
) {
  
  final_summary <-
    rbind(
      final_summary,
      data.frame(
        Analysis =
          "Asymptotic normality",
        Status =
          "Completed",
        stringsAsFactors = FALSE
      )
    )
}

# ------------------------------------------------------------
# Response transformation
# ------------------------------------------------------------

if (
  exists(
    "response_transformation_summary"
  )
) {
  
  final_summary <-
    rbind(
      final_summary,
      data.frame(
        Analysis =
          "Response transformation",
        Status =
          "Completed",
        stringsAsFactors = FALSE
      )
    )
}

# ------------------------------------------------------------
# Error distribution
# ------------------------------------------------------------

if (
  exists(
    "error_distribution_summary"
  )
) {
  
  final_summary <-
    rbind(
      final_summary,
      data.frame(
        Analysis =
          "Error distribution",
        Status =
          "Completed",
        stringsAsFactors = FALSE
      )
    )
}

# ------------------------------------------------------------
# Dimensionality
# ------------------------------------------------------------

if (
  exists(
    "dimension_summary"
  )
) {
  
  final_summary <-
    rbind(
      final_summary,
      data.frame(
        Analysis =
          "Dimensionality",
        Status =
          "Completed",
        stringsAsFactors = FALSE
      )
    )
}

# ------------------------------------------------------------
# Misspecification
# ------------------------------------------------------------

if (
  exists(
    "misspecification_summary"
  )
) {
  
  final_summary <-
    rbind(
      final_summary,
      data.frame(
        Analysis =
          "Model misspecification",
        Status =
          "Completed",
        stringsAsFactors = FALSE
      )
    )
}

# ------------------------------------------------------------
# Kendall test
# ------------------------------------------------------------

if (
  exists(
    "kendall_summary"
  )
) {
  
  final_summary <-
    rbind(
      final_summary,
      data.frame(
        Analysis =
          "Kendall permutation test",
        Status =
          "Completed",
        stringsAsFactors = FALSE
      )
    )
}

# ------------------------------------------------------------
# Tau-star diagnostic
# ------------------------------------------------------------

if (
  exists(
    "tau_star_summary"
  )
) {
  
  final_summary <-
    rbind(
      final_summary,
      data.frame(
        Analysis =
          "Tau-star diagnostic",
        Status =
          "Completed",
        stringsAsFactors = FALSE
      )
    )
}

# ------------------------------------------------------------
# Dependence diagnostic
# ------------------------------------------------------------

if (
  exists(
    "dependence_diagnostic_summary"
  )
) {
  
  final_summary <-
    rbind(
      final_summary,
      data.frame(
        Analysis =
          "Dependence diagnostics",
        Status =
          "Completed",
        stringsAsFactors = FALSE
      )
    )
}

# ------------------------------------------------------------
# Bandwidth sensitivity
# ------------------------------------------------------------

if (
  exists(
    "bandwidth_summary"
  )
) {
  
  final_summary <-
    rbind(
      final_summary,
      data.frame(
        Analysis =
          "Bandwidth sensitivity",
        Status =
          "Completed",
        stringsAsFactors = FALSE
      )
    )
}

# ============================================================
# 5. Display final summary
# ============================================================

cat(
  "\n============================================\n"
)

cat(
  "AVAILABLE ANALYSES\n"
)

cat(
  "============================================\n\n"
)

print(
  final_summary
)

# ============================================================
# 6. Save final summary
# ============================================================

write.csv(
  final_summary,
  file =
    file.path(
      final_directory,
      "final_simulation_summary.csv"
    ),
  row.names = FALSE
)

# ============================================================
# 7. Save session information
# ============================================================

sink(
  file =
    file.path(
      final_directory,
      "sessionInfo.txt"
    )
)

print(
  sessionInfo()
)

sink()

# ============================================================
# 8. Save workspace objects
# ============================================================

save(
  list = ls(),
  file =
    file.path(
      final_directory,
      "functional_transformation_simulation.RData"
    )
)

# ============================================================
# 9. Final completion message
# ============================================================

cat(
  "\n============================================\n"
)

cat(
  "SIMULATION STUDY COMPLETED\n"
)

cat(
  "============================================\n\n"
)

cat(
  "All available summaries have been exported to:\n"
)

cat(
  final_directory,
  "\n\n",
  sep = ""
)

cat(
  "The principal summary file is:\n"
)

cat(
  file.path(
    final_directory,
    "final_simulation_summary.csv"
  ),
  "\n",
  sep = ""
)

cat(
  "\nPart 16 completed successfully.\n"
)

