# Assumptions ----

get_generation_dist <- function(target, shift_mean = 0, shift_sd = 0) {
  switch(target,
         `SARS-N1` = get_discrete_lognormal(
           unit_mean = 3.0 + shift_mean,
           unit_sd = 1.5 + shift_sd,
           maxX = 8 + shift_mean + 2*shift_sd,
           include_zero = F,
         ),
         `SARS-N2` = get_discrete_lognormal(
           unit_mean = 3.0 + shift_mean,
           unit_sd = 1.5 + shift_sd,
           maxX = 8 + shift_mean + 2*shift_sd,
           include_zero = F,
         ),
         `IAV-M` = get_discrete_gamma_shifted(
           gamma_mean = 2.6 + shift_mean,
           gamma_sd = 1.7 + shift_sd,
           maxX = 12 + shift_mean + 2*shift_sd
         ),
         `IBV-M` = get_discrete_gamma_shifted(
           gamma_mean = 2.6 + shift_mean,
           gamma_sd = 1.7 + shift_sd,
           maxX = 12 + shift_mean + 2*shift_sd
         ),
         `RSV-N` =  get_discrete_gamma_shifted(
           gamma_mean = 7.5 + shift_mean,
           gamma_sd = 2.1 + shift_sd,
           maxX = 14 + shift_mean + 2*shift_sd
         ),
         `RSV-A` =  get_discrete_gamma_shifted(
           gamma_mean = 7.5 + shift_mean,
           gamma_sd = 2.1 + shift_sd,
           maxX = 14 + shift_mean + 2*shift_sd
         ),
         `RSV-B` =  get_discrete_gamma_shifted(
           gamma_mean = 7.5 + shift_mean,
           gamma_sd = 2.1 + shift_sd,
           maxX = 14 + shift_mean + 2*shift_sd
         ),
         stop(paste("PCR target", target, "not found."))
  )
}

# This is just a dummy, since we model shedding load distributions indexed by
# date of infection, not by date of symptom onset, i.e. no incubation dist needed
get_incubation_dist <- function(target, shift_mean = 0, shift_sd = 0) {
  switch(target,
         `SARS-N1` = get_discrete_gamma(
           gamma_mean = 0.2 + shift_mean,
           gamma_sd = 0.01 + shift_sd,
           maxX = 1 + shift_mean + 2*shift_sd
         ),
         `SARS-N2` = get_discrete_gamma(
           gamma_mean = 0.2 + shift_mean,
           gamma_sd = 0.01 + shift_sd,
           maxX = 1 + shift_mean + 2*shift_sd
         ),
         `IAV-M` = get_discrete_gamma(
           gamma_mean = 0.2 + shift_mean,
           gamma_sd = 0.01 + shift_sd,
           maxX = 1 + shift_mean + 2*shift_sd
         ),
         `IBV-M` = get_discrete_gamma(
           gamma_mean = 0.2 + shift_mean,
           gamma_sd = 0.01 + shift_sd,
           maxX = 1 + shift_mean + 2*shift_sd
         ),
         `RSV-N` =  get_discrete_gamma(
           gamma_mean = 0.2 + shift_mean,
           gamma_sd = 0.01 + shift_sd,
           maxX = 1 + shift_mean + 2*shift_sd
         ),
         `RSV-A` =  get_discrete_gamma(
           gamma_mean = 0.2 + shift_mean,
           gamma_sd = 0.01 + shift_sd,
           maxX = 1 + shift_mean + 2*shift_sd
         ),
         `RSV-B` =  get_discrete_gamma(
           gamma_mean = 0.2 + shift_mean,
           gamma_sd = 0.01 + shift_sd,
           maxX = 1 + shift_mean + 2*shift_sd
         ),
         stop(paste("PCR target", target, "not found."))
  )
}

get_shedding_dist <- function(target, shift_mean = 0, shift_cv = 0) {
  switch(target,
         `SARS-N1` = get_discrete_gamma(
           gamma_mean = 12.4 + shift_mean,
           gamma_cv = 0.7 + shift_cv
         ),
         `SARS-N2` = get_discrete_gamma(
           gamma_mean = 12.4 + shift_mean,
           gamma_cv = 0.7 + shift_cv
         ),
         `IAV-M` = get_discrete_gamma(
           gamma_mean = 2.5 + shift_mean,
           gamma_cv = 0.34 + shift_cv
         ),
         `IBV-M` = get_discrete_gamma(
           gamma_mean = 2.5 + shift_mean,
           gamma_cv = 0.34 + shift_cv
         ),
         `RSV-N` =  get_discrete_gamma(
           gamma_mean = 6.7 + shift_mean,
           gamma_cv = 0.35 + shift_cv
         ),
         `RSV-A` =  get_discrete_gamma(
           gamma_mean = 6.7 + shift_mean,
           gamma_cv = 0.35 + shift_cv
         ),
         `RSV-B` =  get_discrete_gamma(
           gamma_mean = 6.7 + shift_mean,
           gamma_cv = 0.35 + shift_cv
         ),
         stop(paste("PCR target", target, "not found."))
  )
}

get_shedding_dist_info_fixed <- function(target, shift_mean = 0, shift_cv = 0) {
  info <- list()
  info$estimate <-  FALSE
  info$shedding_reference <- "symptom_onset"
  
  info$shedding_dist <- get_shedding_dist(target, shift_mean, shift_cv)
  
  return(info)
}

get_shedding_dist_info_estimate <- function(target, shift_mean = 0, shift_cv = 0) {
  info <- list()
  info$estimate <- TRUE
  info$shedding_reference <- "infection"
  
  if (target == "SARS-N1" | target == "SARS-N2") {
    info$shedding_dist_mean_prior_mean = c(9.16, 15.69)
    info$shedding_dist_mean_prior_sd = c(0.83, 0.79)
    info$shedding_dist_cv_prior_mean = c(0.87, 0.53)
    info$shedding_dist_cv_prior_sd = c(0.06, 0.04)
    info$shedding_dist_type <- "gamma"
  } else if (target == "IAV-M" | target == "IBV-M") {
    info$shedding_dist_mean_prior_mean = 2.50
    info$shedding_dist_mean_prior_sd = 0.19
    info$shedding_dist_cv_prior_mean = 0.34
    info$shedding_dist_cv_prior_sd = 0.01
    info$shedding_dist_type <- "gamma"
  } else if (target %in% c("RSV-N","RSV-A","RSV-B")) {
    info$shedding_dist_mean_prior_mean <- 6.76
    info$shedding_dist_mean_prior_sd <- 1.17
    info$shedding_dist_cv_prior_mean <- 0.35
    info$shedding_dist_cv_prior_sd <- 0.05
    info$shedding_dist_type <- "gamma"
  } else {
    stop(paste("PCR target", target, "not found."))
  }
  
  return(info)
}

get_load_per_case <- function(wwtp_select, target_select, multiplier = 1, load_per_case_data = NULL) {
  reporting_prop <- switch(target_select,
                           `SARS-N1` = 1,
                           `SARS-N2` = 1,
                           `IAV-M` = 1,
                           `IBV-M` = 1,
                           `RSV-N` = 1,
                           `RSV-A` = 1,
                           `RSV-B` = 1
                           )
  if (is.null(load_per_case_data)) {
    load_per_case_data <- tibble::tribble(
      ~target, ~wwtp, ~load_per_case,
      # SARS-N1
      "SARS-N1", "ARA Altenrhein", 146181818122,
      "SARS-N1", "ARA Basel/Prorheno",  89159331949,
      "SARS-N1", "ARA Buholz", 197376292399,
      "SARS-N1", "ARA Chur", 121038933187,
      "SARS-N1", "ARA Region Bern", 110881842596,
      "SARS-N1", "ARA Schwyz", 195384124819,
      "SARS-N1", "ARA Sensetal",  61443310601,
      "SARS-N1", "ARA Werdhoelzli", 203251833882,
      "SARS-N1", "ARA Zuchwil", 115071113816,
      "SARS-N1", "CDA Lugano", 148745496205,
      "SARS-N1", "STEP Aire", 156555905066,
      "SARS-N1", "STEP Neuchatel", 132906148872,
      "SARS-N1", "STEP Porrentruy",  73737612282,
      "SARS-N1", "STEP Vidy", 253930670773,
      # SARS-N2 (assumed to be 1.6x higher than SARS-N1 in concentration)
      "SARS-N2", "ARA Altenrhein", 199004805702,
      "SARS-N2", "ARA Basel/Prorheno", 100323753439,
      "SARS-N2", "ARA Buholz", 234225747868,
      "SARS-N2", "ARA Chur", 172533729746,
      "SARS-N2", "ARA Region Bern", 134561007811,
      "SARS-N2", "ARA Schwyz", 244418078803,
      "SARS-N2", "ARA Sensetal",  90213918792,
      "SARS-N2", "ARA Werdhoelzli", 272628558994,
      "SARS-N2", "ARA Zuchwil", 133740662133,
      "SARS-N2", "CDA Lugano", 195706209753,
      "SARS-N2", "STEP Aire", 195958806184,
      "SARS-N2", "STEP Neuchatel", 158595205833,
      "SARS-N2", "STEP Porrentruy",  94288209968,
      "SARS-N2", "STEP Vidy", 294109900150,
      # IAV-M
      "IAV-M", "ARA Altenrhein",   5800262087,
      "IAV-M", "ARA Basel/Prorheno",   4695751232,
      "IAV-M", "ARA Buholz",   8594122855,
      "IAV-M", "ARA Chur",   4884244469,
      "IAV-M", "ARA Region Bern",   3144171320,
      "IAV-M", "ARA Schwyz",   9997787656,
      "IAV-M", "ARA Sensetal",   2210166604,
      "IAV-M", "ARA Werdhoelzli",  10012626492,
      "IAV-M", "ARA Zuchwil",   4504546931,
      "IAV-M", "CDA Lugano",   6203485874,
      "IAV-M", "STEP Aire",   6065668819,
      "IAV-M", "STEP Neuchatel",   6780428075,
      "IAV-M", "STEP Porrentruy",   4943481313,
      "IAV-M", "STEP Vidy",   7950390016,
      # IBV-M
      "IBV-M", "ARA Werdhoelzli",  23923706651,
      "IBV-M", "ARA Zuchwil", 40406633534,
      "IBV-M", "STEP Aire", 17610441005,
      "IBV-M", "CDA Lugano",19043336922,
      "IBV-M", "ARA Chur", 8909904201,
      "IBV-M", "ARA Altenrhein", 13304208755,
      "IBV-M", "ARA Sensetal", 5834678236,
      "IBV-M", "STEP Vidy",  24487863356,
      "IBV-M", "STEP Neuchatel",  18132521448,
      "IBV-M", "ARA Buholz", 17696709959,
      "IBV-M", "ARA Basel/Prorheno", 15732238608,
      "IBV-M", "ARA Region Bern", 16317539852,
      "IBV-M", "ARA Schwyz", 17292993096,
      "IBV-M", "STEP Porrentruy", 9843175122,
      # RSV-N
      "RSV-N", "ARA Altenrhein",  14894827704,
      "RSV-N", "ARA Basel/Prorheno",   8903311398,
      "RSV-N", "ARA Buholz",  23218537464,
      "RSV-N", "ARA Chur",  14537185547,
      "RSV-N", "ARA Region Bern",   9878053196,
      "RSV-N", "ARA Schwyz",  25485002695,
      "RSV-N", "ARA Sensetal",   6995582047,
      "RSV-N", "ARA Werdhoelzli",  24070579221,
      "RSV-N", "ARA Zuchwil",  15396036925,
      "RSV-N", "CDA Lugano",  17701561882,
      "RSV-N", "STEP Aire",  13389833020,
      "RSV-N", "STEP Neuchatel",  16631932208,
      "RSV-N", "STEP Porrentruy",  12200209366,
      "RSV-N", "STEP Vidy",  21078217567,
      # RSV A is assumed to be the same as RSV-N
      "RSV-A", "ARA Altenrhein",  14894827704,
      "RSV-A", "ARA Basel/Prorheno",   8903311398,
      "RSV-A", "ARA Buholz",  23218537464,
      "RSV-A", "ARA Chur",  14537185547,
      "RSV-A", "ARA Region Bern",   9878053196,
      "RSV-A", "ARA Schwyz",  25485002695,
      "RSV-A", "ARA Sensetal",   6995582047,
      "RSV-A", "ARA Werdhoelzli",  24070579221,
      "RSV-A", "ARA Zuchwil",  15396036925,
      "RSV-A", "CDA Lugano",  17701561882,
      "RSV-A", "STEP Aire",  13389833020,
      "RSV-A", "STEP Neuchatel",  16631932208,
      "RSV-A", "STEP Porrentruy",  12200209366,
      "RSV-A", "STEP Vidy",  21078217567,
      # RSV B is assumed to be the same as RSV-N
      "RSV-B", "ARA Altenrhein",  14894827704,
      "RSV-B", "ARA Basel/Prorheno",   8903311398,
      "RSV-B", "ARA Buholz",  23218537464,
      "RSV-B", "ARA Chur",  14537185547,
      "RSV-B", "ARA Region Bern",   9878053196,
      "RSV-B", "ARA Schwyz",  25485002695,
      "RSV-B", "ARA Sensetal",   6995582047,
      "RSV-B", "ARA Werdhoelzli",  24070579221,
      "RSV-B", "ARA Zuchwil",  15396036925,
      "RSV-B", "CDA Lugano",  17701561882,
      "RSV-B", "STEP Aire",  13389833020,
      "RSV-B", "STEP Neuchatel",  16631932208,
      "RSV-B", "STEP Porrentruy",  12200209366,
      "RSV-B", "STEP Vidy",  21078217567
    )
  }
  sel <- load_per_case_data |> dplyr::filter(wwtp == wwtp_select, target == target_select)
  if (nrow(sel)==0) {
    stop(paste("Load per case for", wwtp_select, target_select, "not found."))
  } else {
    return(sel$load_per_case[1] * reporting_prop * multiplier)
  }
}