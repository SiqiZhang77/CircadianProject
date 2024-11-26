# Function to calculate model values
model_function <- function(t, A, gamma, omega, phi, y_shift) {
  A * exp(-gamma * t) * cos(omega * t + phi) + y_shift
}

# Function to simulate data with noise
simulate_data <- function(t, params, sigma) {
  predicted <- model_function(t, params[1], params[2], params[3], params[4], params[5])
  noise <- rnorm(length(t), mean = 0, sd = sigma)
  simulated_data <- predicted + noise
  return(simulated_data)
}

# Function to find peaks
find_peaks <- function(data, minpeakheight, minpeakdistance) {
  findpeaks(data, minpeakheight = minpeakheight, minpeakdistance = minpeakdistance)
}

# Asymmetrical oscillatory model function
asym_model_function <- function(t, A_max, eta_max, A_min, eta_min, omega, phi) {
  A_max_t <- A_max * exp(-eta_max * t)
  A_min_t <- A_min * exp(-eta_min * t)
  result <- 0.5 * (A_max_t + A_min_t) + 0.5 * (A_max_t - A_min_t) * cos(omega * t + phi)
  return(result)
}

# Simulate data for asymmetrical oscillatory model
simulate_asymm_func <- function(t, params, sigma) {
  predicted <- asym_model_function(t, params[1], params[2], params[3], params[4], params[5], params[6])
  noise <- rnorm(length(t), mean = 0, sd = sigma)
  simulated_data <- predicted + noise
  return(simulated_data)
}

# Bayesian variance estimation
bayesian_variance_estimation <- function(observed_variances, n, alpha_prior = 1, beta_prior = 1) {
  weights <- numeric(length(observed_variances))
  for (i in 1:length(observed_variances)) {
    S2_i <- observed_variances[i]
    alpha_post <- alpha_prior + (n - 1) / 2
    beta_post <- beta_prior + S2_i * (n - 1) / 2
    posterior_mean <- alpha_post / beta_post
    weights[i] <- posterior_mean
  }
  return(weights)
}

# Method of Moments for prior parameters
calculate_moments <- function(observed_variances) {
  n <- length(observed_variances)
  sample_mean <- mean(observed_variances)
  sample_variance <- sum((observed_variances - sample_mean)^2) / (n - 1)
  beta_prior <- sample_mean / sample_variance
  alpha_prior <- sample_mean^2 / sample_variance
  return(list(alpha_prior = alpha_prior, beta_prior = beta_prior))
}

# Calculate variance at each time point
calculate_variance <- function(replicates) {
  n <- nrow(replicates)
  mean_vals <- rowMeans(replicates)
  variance <- rowSums((replicates - mean_vals)^2) / (n - 1)
  return(variance)
}

# Weighted nonlinear least squares estimation
weighted_nls <- function(t, replicates, initial_params) {
  variances <- calculate_variance(replicates)
  fit <- nlsLM(
    y ~ model_function(t, A, gamma, omega, phi, y_shift),
    start = initial_params,
    data = list(t = t, y = rowMeans(replicates)),
    weights = 1 / variances
  )
  return(fit)
}

# Simulate and estimate parameters
simulate_and_estimate <- function(true_params, sigma0 = 10) {
  t <- seq(0, 100, by = 0.1)
  y_sim <- simulate_data(t, true_params, sigma0)
  peaks <- find_peaks(y_sim, minpeakheight = 100, minpeakdistance = 15)
  t_peaks <- t[peaks[, 2]]
  y_peaks <- peaks[, 1]
  log_y_peaks <- log(y_peaks)
  fit_linear <- lm(log_y_peaks ~ t_peaks)
  log_A <- coef(fit_linear)[1]
  init_gamma <- -coef(fit_linear)[2]
  init_A <- exp(log_A)
  Delta <- t_peaks[2] - t_peaks[1]
  init_omega <- 2 * pi / Delta
  fit_nls <- try(
    nls(
      y_sim ~ model_function(t, A, gamma, omega, phi, y_shift),
      start = list(A = init_A, gamma = init_gamma, omega = init_omega, phi = 0.1, y_shift = 0),
      data = data.frame(t = t, y_sim = y_sim)
    ),
    silent = TRUE
  )
  if (class(fit_nls) == "try-error") {
    return(NULL)
  }
  param_estimates <- coef(fit_nls)
  conf_intervals <- confint2(fit_nls, level = 0.99, method = 'asymptotic')
  return(list(estimates = param_estimates, conf_intervals = conf_intervals))
}
