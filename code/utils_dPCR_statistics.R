## Simulation of dPCR measurements ----

# Simulation of measurements based on the droplet count model
sim_conc_pre <- function(lambda, n, c, cv, n_replicates = 1, distribution = "gamma", n_samples = 1) {
  if (cv==0) {
    lam <- rep(lambda, n_samples)
  } else {
    if (distribution == "gamma") {
      lam <- rgamma2(n_samples, lambda, cv = cv)
    } else if (distribution == "lnorm") {
      lam <- rlnorm3(n_samples, lambda, cv = cv)
    } else {
      stop("Distribution not supported.")
    }
  }
  p <- 1 - exp(-lam * c)
  concentrations <- sapply(p, function(p) {
    sim_drop <- rbinom(n_replicates, n, p)
    sim_conc <- - log(1 - (sim_drop / n)) * (1/c)
    return(mean(sim_conc))
  })
  return(concentrations)
}

# Simulation of measurements based on the derived CV and non-detection probability
# and using a continuous likelihood
sim_theoretical_conc_pre <- function(lambda, n, c, cv, n_replicates = 1, distribution = "gamma", likelihood = "gamma", n_samples = 1, use_taylor = T) {
  if (use_taylor) {
    conc_cv <- est_cv_pre_taylor(lambda, n, c, cv, n_replicates)
  } else if (distribution == "gamma") {
    conc_cv <- est_cv_pre_gamma(lambda, n, c, cv, n_replicates)
  } else if (distribution == "lnorm") {
    conc_cv <- est_cv_pre_lnorm(lambda, n, c, cv, n_replicates)
  } else {
    stop("Distribution not supported.")
  }
  
  p_nondetect <- est_nondetection_prob_pre(
    lambda = lambda, cv = cv, n = n, c = c, n_replicates = n_replicates, noise_dist = distribution, method = "lambert"
    )
  
  if (likelihood == "gamma") {
    concentrations <- rgamma2(n_samples, mean = lambda, cv = conc_cv)
  } else if (likelihood %in% c("lnorm", "lognormal", "log-normal")) {
    concentrations <- extraDistr::rtnorm(n_samples, mean = lambda, sd = lambda * conc_cv, a = 0)
  } else {
    stop("Likelihood not supported.")
  }
  
  detection <- extraDistr::rbern(length(concentrations), prob = 1-p_nondetect)
  concentrations <- detection * concentrations
  
  return(concentrations)
}

## Coefficient of variation ----

# Simulation of coefficient of variation via the droplet count model
sim_cv_pre <- function(lambda, n, c, cv, n_replicates = 1, distribution = "gamma", n_samples = 10000) {
  if (cv==0) {
    lam <- rep(lambda, n_samples)
  } else {
    if (distribution == "gamma") {
      lam <- rgamma2(n_samples, lambda, cv = cv)
    } else if (distribution %in% c("lnorm", "lognormal", "log-normal")) {
      lam <- rlnorm3(n_samples, lambda, cv = cv)
    } else {
      stop("Distribution not supported.")
    }
  }
  p <- 1 - exp(-lam * c)
  concentrations <- sapply(p, function(p) {
    sim_drop <- rbinom(n_replicates, n, p)
    sim_conc <- - log(1 - (sim_drop / n)) * (1/c)
    return(mean(sim_conc))
  })
  return(sd(concentrations, na.rm = T) / mean(concentrations, na.rm = T))
}

est_cv <- function(lambda, n, c, n_replicates = 1) {
  return(sqrt(exp(lambda*c)-1) / (sqrt(n) * lambda * c) / sqrt(n_replicates))
}

est_cv_pre_gamma <- function(lambda, n, c, cv, n_replicates = 1) {
  if (cv==0) {
    return(est_cv(lambda, n, c, n_replicates = n_replicates))
  } else {
    if (c>=lambda/(cv^2*lambda^2)) {
      # only defined if c<lambda/(cv^2*lambda^2), otherwise MGF of gamma is undefined
      # in other words: maximum concentration is 1/(c*cv^2)
      var_unexpl = NA
    } else {
      # regular case
      var_unexpl <- 1/(c^2*n*n_replicates) * ((1 - lambda * cv^2 * c)^(-1/(cv^2))-1)
    }
    var_expl <- cv^2 * lambda^2
    return(sqrt(var_unexpl + var_expl) / lambda)
  }
}

est_cv_pre_lnorm <- function(lambda, n, c, cv, n_replicates = 1) {
  if (cv==0) {
    return(est_cv(lambda, n, c, n_replicates = n_replicates))
  } else {
    sigma2 <- log(1 + cv^2)
    mu <- log(lambda) - sigma2/2
    var_unexpl <- 1/(c^2*n*n_replicates) * (expL(c, mu, sigma2) - 1)
    var_expl <- cv^2 * lambda^2
    return(sqrt(var_unexpl + var_expl) / lambda)
  }
}

est_cv_pre_taylor <- function(lambda, n, c, cv, n_replicates = 1) {
  if (cv==0) {
    return(est_cv(lambda, n, c, n_replicates = n_replicates))
  } else {
    var_unexpl <- 1/(c^2*n*n_replicates) * (lambda*c + (lambda^2*c^2 + cv^2*lambda^2*c^2)/2)
    var_expl <- cv^2 * lambda^2
    return(sqrt(var_unexpl + var_expl) / lambda)
  }
}

# same as est_cv_pre_taylor, just simplified (does not require a special case)
est_cv_pre_taylor2 <- function(lambda, n, c, cv, n_replicates = 1) {
  return(sqrt(cv^2 + 1/(n*c*n_replicates*lambda) + (1+cv^2)/(2*n*n_replicates)))
}

## Non-detection probability ----

#' Get probability of non-detection given a number of gene copies per reaction in the PCR
#'
#' @return probability of non-detection (zero measurement)
sim_nondetection_prob <- function(lambda, c, n = 20000, n_replicates = 1) {
  nondetection_prob <- dbinom(0, size = n, prob = 1 - exp(-lambda * c))
  # same as -exp(-lambda * droplet_vol * n)
  return(nondetection_prob^n_replicates)
}

#' Get probability of non-detection given a number of gene copies per reaction in the PCR
#'
#' @param gc_per_rxn overall number of gene copies in the reaction
#'
#' @return probability of non-detection (zero measurement)
sim_nondetection_prob_pre <- function(lambda, cv, c, n = 20000, n_replicates = 1, noise_dist = "gamma") {
  if (cv == 0) {
    return(sim_nondetection_prob(lambda = lambda, c = c, n = n, n_replicates = n_replicates))
  }
  if (noise_dist == "gamma") {
    lambdas <- rgamma(10000, shape = 1/cv^2, rate = 1/(lambda * cv^2))
  } else if (noise_dist == "lognormal") {
    lambdas <- rlnorm3(10000, lambda, cv)
  } else {
    stop(paste("Noise dist", noise_dist, "not supported."))
  }
  nondetection_prob <- sum(sapply(lambdas, function(lambda) {
    mean(rbinom(n_replicates, size = n, prob = 1 - exp(-lambda * c))) == 0
  })) / length(lambdas)
  return(nondetection_prob)
}



est_nondetection_prob <- function(lambda, c, n = 20000, n_replicates = 1) {
  nondetection_prob <- exp(-lambda * c * n)
  return(nondetection_prob^n_replicates)
}

est_nondetection_prob_pre <- function(lambda, cv, c, n = 20000, n_replicates = 1, noise_dist = "gamma", method = "lambert") {
  if (cv == 0) {
    return(est_nondetection_prob(lambda, n = n, c = c, n_replicates = n_replicates))
  }
  
  if (lambda == 0) {
    return(1)
  }
  
  if (noise_dist == "gamma") {
    nondetection_prob <- (1 + n * lambda * c * cv^2 * n_replicates)^(-1/cv^2)
    # can also be written as:
    # --> nondetection_prob = exp(log(1 + n * lambda * c * cv^2 * n_replicates) *(-1/cv^2))
    # --> log(nondetection_prob) = -log(1 + n * lambda * c * cv^2 * n_replicates)/cv^2
  } else if (noise_dist %in% c("lnorm", "lognormal", "log-normal")) {
    sigma2 = log(1 + cv^2)
    mu = log(lambda) - sigma2/2
    if (method == "lambert") {
      nondetection_prob <- expL(theta = -c * n * n_replicates, mu = mu, sigma2 = sigma2)
    } else if (method == "integrate") {
      nondetection_prob <- expL_int(theta = -c * n * n_replicates, mu = mu, sigma2 = sigma2)
    }
  } else {
    stop(paste("Noise dist", noise_dist, "not supported."))
  }
  return(nondetection_prob)
}

## Expected value ----
sim_expectation_pre <- function(lambda, n, c, cv, n_replicates = 1, distribution = "gamma", n_samples = 10000) {
  if (distribution == "gamma") {
    lam <- rgamma2(n_samples, lambda, cv = cv)
  } else if (distribution == "lnorm") {
    lam <- rlnorm3(n_samples, lambda, cv = cv)
  } else {
    stop("Distribution not supported.")
  }
  p <- 1 - exp(-lam * c)
  concentrations <- sapply(p, function(p) {
    sim_drop <- rbinom(n_replicates, n, p)
    sim_conc <- - log(1 - (sim_drop / n)) * (1/c)
    return(mean(sim_conc))
  })
  return(mean(concentrations, na.rm = T))
}

## Standard deviation ----
est_sd <- function(lambda, n, c, n_replicates = 1) {
  return(sqrt(exp(lambda*c)-1) / (sqrt(n) * c) / sqrt(n_replicates))
}