## Simulation of dPCR measurements ----

#' Simulate number of valid partitions across dPCR runs
#'
#' @param n_samples Number of different samples
#' @param n_replicates Number of replicates per sample
#' @param max_partitions Maximum number of partitions
#' @param partition_loss_logit_mean Mean proportion of invalid partitions
#'   (logit-level)
#' @param partition_loss_logit_sd Sd of proportion of invalid partitions
#'   (logit-level)
#' @param partition_loss_max Upper bound for partition loss proportion
#'
#' @return Simulated partition numbers, structured as a list with one element
#'   per sample, where each element is a vector giving the partition numbers
#'   across replicates
sim_partitions <- function(n_samples, n_replicates, max_partitions, partition_loss_logit_mean, partition_loss_logit_sd, partition_loss_max) {
  lapply(1:n_samples, function(x) extraDistr::rtbinom(n_replicates, size = max_partitions, prob = 1-plogis(rnorm(n_replicates, mean = partition_loss_logit_mean, sd = partition_loss_logit_sd)), a = max_partitions * (1-partition_loss_max)))
}

#' Simulate concentration measurements based on the dPCR partition count model
#'
#' @param lambda Expected concentration in the sample
#' @param m_partitions Number of partitions (vector of length n_replicates or scalar)
#' @param c Conversion factor, i.e. volume v * scaling factor s
#' @param n_replicates Number of replicates
#' @param cv Pre-PCR coefficient of variation
#' @param distribution Pre-PCR variation distribution
#' @param n_samples How many measurements of the sample do you want to simulate?
#'
#' @return A vector of simulated measurements
sim_conc_pre <- function(lambda, m_partitions, c, cv, n_replicates = 1, distribution = "gamma", n_samples = 1, get_partitions = FALSE) {
  # Ensure m_partitions is a vector of length n_replicates
  if (is.list(m_partitions)) {
    stopifnot(length(m_partitions[[1]]) == n_replicates)
    stopifnot(length(m_partitions) == n_samples)
  } else {
    if (length(m_partitions) == 1) m_partitions <- rep(m_partitions, n_replicates)
  }
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
  sim_res <- mapply(function(p, m_partitions) {
    sim_drop <- mapply(function(mp) rbinom(1, mp, p), m_partitions)
    sim_positive_partitions <- sum(sim_drop)
    sim_total_partitions <- sum(m_partitions)
    sim_conc <- - log(1 - (sim_positive_partitions / sim_total_partitions)) * (1/c)
    return(list(
      concentration = sim_conc,
      positive_partitions = sim_positive_partitions/n_replicates, # average over replicates
      total_partitions = sim_total_partitions/n_replicates # average over replicates
      ))
  }, p = p, m_partitions = m_partitions, SIMPLIFY = FALSE)
  if (get_partitions) {
    return(data.frame(
      concentration = sapply(sim_res, function(x) x$concentration),
      positive_partitions = sapply(sim_res, function(x) x$positive_partitions),
      total_partitions = sapply(sim_res, function(x) x$total_partitions)
    ))
  } else {
    return(sapply(sim_res, function(x) x$concentration))
  }
}

# Simulation of measurements based on the derived CV and non-detection probability
# and using a continuous likelihood
sim_theoretical_conc_pre <- function(lambda, m_partitions, c, cv, n_replicates = 1, distribution = "gamma", likelihood = "gamma", n_samples = 1, use_taylor = T) {
  if (length(m_partitions) == 1) m_partitions <- rep(m_partitions, n_replicates)
  if (use_taylor) {
    conc_cv <- est_cv_pre_taylor(lambda, m_partitions, c, cv, n_replicates)
  } else if (distribution == "gamma") {
    conc_cv <- est_cv_pre_gamma(lambda, m_partitions, c, cv, n_replicates)
  } else if (distribution == "lnorm") {
    conc_cv <- est_cv_pre_lnorm(lambda, m_partitions, c, cv, n_replicates)
  } else {
    stop("Distribution not supported.")
  }
  
  p_nondetect <- est_nondetection_prob_pre(
    lambda = lambda, cv = cv, m_partitions = m_partitions, c = c, n_replicates = n_replicates, noise_dist = distribution, method = "lambert"
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
sim_cv_pre <- function(lambda, m_partitions, c, cv, n_replicates = 1, distribution = "gamma", n_samples = 10000) {
  if (length(m_partitions) == 1) m_partitions <- rep(m_partitions, n_replicates)
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
    sim_drop <- mapply(function(mp) rbinom(1, mp, p), m_partitions)
    sim_conc <- - log(1 - (sum(sim_drop) / sum(m_partitions))) * (1/c)
    return(sim_conc)
  })
  return(sd(concentrations, na.rm = T) / mean(concentrations, na.rm = T))
}

est_cv <- function(lambda, m_partitions, c, n_replicates = 1) {
  if (length(m_partitions) == 1) m_partitions <- rep(m_partitions, n_replicates)
  return(sqrt(exp(lambda*c)-1) / (sqrt(sum(m_partitions)) * lambda * c))
}

est_cv_pre_gamma <- function(lambda, m_partitions, c, cv, n_replicates = 1) {
  if (length(m_partitions) == 1) m_partitions <- rep(m_partitions, n_replicates)
  if (cv==0) {
    return(est_cv(lambda, m_partitions, c, n_replicates = n_replicates))
  } else {
    if (c>=lambda/(cv^2*lambda^2)) {
      # only defined if c<lambda/(cv^2*lambda^2), otherwise MGF of gamma is undefined
      # in other words: maximum concentration is 1/(c*cv^2)
      var_unexpl = NA
    } else {
      # regular case
      var_unexpl <- 1/(c^2*sum(m_partitions)) * ((1 - lambda * cv^2 * c)^(-1/(cv^2))-1)
    }
    var_expl <- cv^2 * lambda^2
    return(sqrt(var_unexpl + var_expl) / lambda)
  }
}

est_cv_pre_lnorm <- function(lambda, m_partitions, c, cv, n_replicates = 1) {
  if (length(m_partitions) == 1) m_partitions <- rep(m_partitions, n_replicates)
  if (cv==0) {
    return(est_cv(lambda, m_partitions, c, n_replicates = n_replicates))
  } else {
    sigma2 <- log(1 + cv^2)
    mu <- log(lambda) - sigma2/2
    var_unexpl <- 1/(c^2*sum(m_partitions)) * (expL(c, mu, sigma2) - 1)
    var_expl <- cv^2 * lambda^2
    return(sqrt(var_unexpl + var_expl) / lambda)
  }
}

est_cv_pre_taylor <- function(lambda, m_partitions, c, cv, n_replicates = 1) {
  if (length(m_partitions) == 1) m_partitions <- rep(m_partitions, n_replicates)
  if (cv==0) {
    return(est_cv(lambda, m_partitions, c, n_replicates = n_replicates))
  } else {
    var_unexpl <- 1/(c^2*sum(m_partitions)) * (lambda*c + (lambda^2*c^2 + cv^2*lambda^2*c^2)/2)
    var_expl <- cv^2 * lambda^2
    return(sqrt(var_unexpl + var_expl) / lambda)
  }
}

# same as est_cv_pre_taylor, just simplified (does not require a special case)
est_cv_pre_taylor2 <- function(lambda, m_partitions, c, cv, n_replicates = 1) {
  if (length(m_partitions) == 1) m_partitions <- rep(m_partitions, n_replicates)
  return(sqrt(cv^2 + 1/(sum(m_partitions)*c*lambda) + (1+cv^2)/(2*sum(m_partitions))))
}

## Non-detection probability ----

#' Get probability of non-detection given a number of gene copies per reaction in the PCR
#'
#' @return probability of non-detection (zero measurement)
sim_nondetection_prob <- function(lambda, c, m_partitions = 20000, n_replicates = 1) {
  if (length(m_partitions) == 1) m_partitions <- rep(m_partitions, n_replicates)
  nondetection_prob <- prod(sapply(m_partitions, function(mp) dbinom(0, size = mp, prob = 1 - exp(-lambda * c))))
  return(nondetection_prob)
}

#' Get probability of non-detection given a number of gene copies per reaction in the PCR
#'
#' @param gc_per_rxn overall number of gene copies in the reaction
#'
#' @return probability of non-detection (zero measurement)
sim_nondetection_prob_pre <- function(lambda, cv, c, m_partitions = 20000, n_replicates = 1, noise_dist = "gamma") {
  if (length(m_partitions) == 1) m_partitions <- rep(m_partitions, n_replicates)
  if (cv == 0) {
    return(sim_nondetection_prob(lambda = lambda, c = c, m_partitions = m_partitions, n_replicates = n_replicates))
  }
  if (noise_dist == "gamma") {
    lambdas <- rgamma(10000, shape = 1/cv^2, rate = 1/(lambda * cv^2))
  } else if (noise_dist == "lognormal") {
    lambdas <- rlnorm3(10000, lambda, cv)
  } else {
    stop(paste("Noise dist", noise_dist, "not supported."))
  }
  nondetection_prob <- mean(sapply(lambdas, function(lambda) {
    prod(sapply(m_partitions, function(mp) rbinom(1, mp, 1 - exp(-lambda * c)) == 0))
  }))
  return(nondetection_prob)
}

est_nondetection_prob <- function(lambda, c, m_partitions = 20000, n_replicates = 1) {
  if (length(m_partitions) == 1) m_partitions <- rep(m_partitions, n_replicates)
  nondetection_prob <- exp(-lambda * c * sum(m_partitions))
  return(nondetection_prob)
}

est_nondetection_prob_pre <- function(lambda, cv, c, m_partitions = 20000, n_replicates = 1, noise_dist = "gamma", method = "lambert") {
  if (length(m_partitions) == 1) m_partitions <- rep(m_partitions, n_replicates)
  if (cv == 0) {
    return(est_nondetection_prob(lambda, m_partitions = m_partitions, c = c, n_replicates = n_replicates))
  }
  if (lambda == 0) {
    return(1)
  }
  total_partitions <- sum(m_partitions)
  if (noise_dist == "gamma") {
    # can also be written as:
    # --> nondetection_prob = exp(log(1 + m * lambda * c * cv^2 * n_replicates) *(-1/cv^2))
    # --> log(nondetection_prob) = -log(1 + m * lambda * c * cv^2 * n_replicates)/cv^2
    nondetection_prob <- (1 + total_partitions * lambda * c * cv^2)^(-1/cv^2)
  } else if (noise_dist %in% c("lnorm", "lognormal", "log-normal")) {
    sigma2 = log(1 + cv^2)
    mu = log(lambda) - sigma2/2
    if (method == "lambert") {
      nondetection_prob <- expL(theta = -c * total_partitions, mu = mu, sigma2 = sigma2)
    } else if (method == "integrate") {
      nondetection_prob <- expL_int(theta = -c * total_partitions, mu = mu, sigma2 = sigma2)
    }
  } else {
    stop(paste("Noise dist", noise_dist, "not supported."))
  }
  return(nondetection_prob)
}

## Expected value ----
sim_expectation_pre <- function(lambda, m_partitions, c, cv, n_replicates = 1, distribution = "gamma", n_samples = 10000) {
  if (length(m_partitions) == 1) m_partitions <- rep(m_partitions, n_replicates)
  if (distribution == "gamma") {
    lam <- rgamma2(n_samples, lambda, cv = cv)
  } else if (distribution == "lnorm") {
    lam <- rlnorm3(n_samples, lambda, cv = cv)
  } else {
    stop("Distribution not supported.")
  }
  p <- 1 - exp(-lam * c)
  concentrations <- sapply(p, function(p) {
    sim_drop <- mapply(function(mp) rbinom(1, mp, p), m_partitions)
    sim_conc <- - log(1 - (sum(sim_drop) / sum(m_partitions))) * (1/c)
    return(sim_conc)
  })
  return(mean(concentrations, na.rm = T))
}

## Standard deviation ----
est_sd <- function(lambda, m_partitions, c, n_replicates = 1) {
  if (length(m_partitions) == 1) m_partitions <- rep(m_partitions, n_replicates)
  return(sqrt(exp(lambda*c)-1) / (sqrt(sum(m_partitions)) * c))
}

## Back-calculation of partitions ----
backcalc_partitions <- function(conc, c = 0.85*0.001, partitions_mean = 20000, partitions_lower = partitions_mean*0.8, partitions_upper = partitions_mean*1.2, tolerance = 0.01) {
  if (conc == 0) return(list(positive = 0, total = NA))
  pre_c <- (1-exp(-conc*c))
  positive_initial <- round(pre_c*partitions_mean)
  positive_mod <- seq(-positive_initial+1, positive_initial, 1)
  positive_mod <- positive_mod[order(abs(positive_mod))]
  
  positive_all <- positive_initial + positive_mod
  total_all <- positive_all / pre_c
  
  matching <- abs(total_all - round(total_all)) < tolerance
  positive_all <- positive_all[matching]
  total_all <- total_all[matching]
  if (total_all[1] > partitions_upper || total_all[1] < partitions_lower) {
    warning(paste0("Back-calculated total partitions (", round(total_all[1],2) ,
                   ") are out of bounds ", 
                   partitions_lower, " to ", partitions_upper, 
                   " for concentration ", conc, "."))
  }
  return(list(positive = positive_all[1], total = round(total_all[1],2)))
}