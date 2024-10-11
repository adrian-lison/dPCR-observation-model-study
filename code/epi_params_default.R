# Assumptions ----

get_generation_dist <- function(target, shift_mean = 0, shift_sd = 0) {
  switch(target,
         `SARS-N1` = get_discrete_lognormal(
           unit_mean = 2.976 + shift_mean,
           unit_sd = 1.479 + shift_sd,
           maxX = 8 + shift_mean + 2*shift_sd,
           include_zero = F,
         ),
         `SARS-N2` = get_discrete_lognormal(
           unit_mean = 2.976 + shift_mean,
           unit_sd = 1.479 + shift_sd,
           maxX = 8 + shift_mean + 2*shift_sd,
           include_zero = F,
         ),
         `IAV-M` = get_discrete_gamma_shifted(
           gamma_mean = 2.6 + shift_mean,
           gamma_sd = 1.5 + shift_sd,
           maxX = 12 + shift_mean + 2*shift_sd
         ),
         `IBV-M` = get_discrete_gamma_shifted(
           gamma_mean = 2.6 + shift_mean,
           gamma_sd = 1.5 + shift_sd,
           maxX = 12 + shift_mean + 2*shift_sd
         ),
         `RSV-N` =  get_discrete_gamma_shifted(
           gamma_mean = 7.5 + shift_mean,
           gamma_sd = 2.1 + shift_sd,
           maxX = 14 + shift_mean + 2*shift_sd
         ),
         stop(paste("PCR target", target, "not found."))
  )
}

get_incubation_dist <- function(target, shift_mean = 0, shift_sd = 0) {
  switch(target,
         `SARS-N1` = get_discrete_gamma(
           gamma_mean = 3.485 + shift_mean,
           gamma_sd = 1.195 + shift_sd,
           maxX = 8 + shift_mean + 2*shift_sd
         ),
         `SARS-N2` = get_discrete_gamma(
           gamma_mean = 3.485 + shift_mean,
           gamma_sd = 1.195 + shift_sd,
           maxX = 8 + shift_mean + 2*shift_sd
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
         stop(paste("PCR target", target, "not found."))
  )
}

get_shedding_dist <- function(target, shift_mean = 0, shift_cv = 0) {
  switch(target,
         `SARS-N1` = get_discrete_gamma(
           gamma_mean = 6.731885 + shift_mean,
           gamma_cv = 6.981995/6.731885 + shift_cv
         ),
         `SARS-N2` = get_discrete_gamma(
           gamma_mean = 6.731885 + shift_mean,
           gamma_cv = 6.981995/6.731885 + shift_cv
         ),
         `IAV-M` = get_discrete_gamma(
           gamma_mean = 2.491217 + shift_mean,
           gamma_cv = 1.004283/2.491217 + shift_cv
         ),
         `IBV-M` = get_discrete_gamma(
           gamma_mean = 2.491217 + shift_mean,
           gamma_cv = 1.004283/2.491217 + shift_cv
         ),
         `RSV-N` =  get_discrete_gamma(
           gamma_mean = 6.7672 + shift_mean,
           gamma_cv = 2.281223/6.7672 + shift_cv
         ),
         stop(paste("PCR target", target, "not found."))
  )
}

get_load_per_case <- function(wwtp, target, multiplier = 1) {
  reporting_prop <- switch(target,
                           `SARS-N1` = 1,
                           `SARS-N2` = 1,
                           `IAV-M` = 1,
                           `IBV-M` = 1,
                           `RSV-N` = 1)
  
  dt <- tibble::tribble(
    ~target_, ~wwtp_, ~load_per_case,
    # SARS-N1
    "SARS-N1", "ARA Werdhoelzli", 1e+11,
    "SARS-N1", "ARA Zuchwil", 1e+11,
    "SARS-N1", "STEP Aire", 1e+11,
    "SARS-N1", "CDA Lugano", 1e+11,
    "SARS-N1", "ARA Chur", 1e+11,
    "SARS-N1", "ARA Altenrhein", 1e+11,
    "SARS-N1", "ARA Sensetal", 1e+11,
    "SARS-N1", "STEP Vidy", 1e+11, # no comparison with catchment-case data yet available
    "SARS-N1", "STEP Neuchatel", 1e+11, # no comparison with catchment-case data yet available
    "SARS-N1", "ARA Buholz", 1e+11, # no comparison with catchment-case data yet available
    "SARS-N1", "ARA Basel/Prorheno", 1e+11, # no comparison with catchment-case data yet available
    "SARS-N1", "ARA Region Bern", 1e+11, # no comparison with catchment-case data yet available
    "SARS-N1", "ARA Schwyz", 1e+11, # no comparison with catchment-case data yet available
    "SARS-N1", "STEP Porrentruy", 1e+11, # no comparison with catchment-case data yet available
    # SARS-N2 (assumed to be 1.6x higher than SARS-N1 in concentration)
    "SARS-N2", "ARA Werdhoelzli", 1e+11*1.6,
    "SARS-N2", "ARA Zuchwil", 1e+11*1.6,
    "SARS-N2", "STEP Aire", 1e+11*1.6,
    "SARS-N2", "CDA Lugano", 1e+11*1.6,
    "SARS-N2", "ARA Chur", 1e+11*1.6,
    "SARS-N2", "ARA Altenrhein", 1e+11*1.6,
    "SARS-N2", "ARA Sensetal", 1e+11*1.6,
    "SARS-N2", "STEP Vidy", 1e+11*1.6, # no comparison with catchment-case data yet available
    "SARS-N2", "STEP Neuchatel", 1e+11*1.6, # no comparison with catchment-case data yet available
    "SARS-N2", "ARA Buholz", 1e+11*1.6, # no comparison with catchment-case data yet available
    "SARS-N2", "ARA Basel/Prorheno", 1e+11*1.6, # no comparison with catchment-case data yet available
    "SARS-N2", "ARA Region Bern", 1e+11*1.6, # no comparison with catchment-case data yet available
    "SARS-N2", "ARA Schwyz", 1e+11*1.6, # no comparison with catchment-case data yet available
    "SARS-N2", "STEP Porrentruy", 1e+11*1.6, # no comparison with catchment-case data yet available
    # IAV-M
    "IAV-M", "ARA Werdhoelzli", 1e+11,
    "IAV-M", "ARA Zuchwil", 1e+11,
    "IAV-M", "STEP Aire", 1e+11,
    "IAV-M", "CDA Lugano", 1e+10,
    "IAV-M", "ARA Chur", 1e+11,
    "IAV-M", "ARA Altenrhein", 1e+11,
    "IAV-M", "ARA Sensetal", 1e+10,
    "IAV-M", "STEP Vidy", 1e+11, # no comparison with catchment-case data yet available
    "IAV-M", "STEP Neuchatel", 1e+11, # no comparison with catchment-case data yet available
    "IAV-M", "ARA Buholz", 1e+11, # no comparison with catchment-case data yet available
    "IAV-M", "ARA Basel/Prorheno", 1e+11, # no comparison with catchment-case data yet available
    "IAV-M", "ARA Region Bern", 1e+11, # no comparison with catchment-case data yet available
    "IAV-M", "ARA Schwyz", 1e+11, # no comparison with catchment-case data yet available
    "IAV-M", "STEP Porrentruy", 1e+11, # no comparison with catchment-case data yet available
    # IBV-M
    "IBV-M", "ARA Werdhoelzli", 1e+11,
    "IBV-M", "ARA Zuchwil", 1e+11,
    "IBV-M", "STEP Aire", 1e+11,
    "IBV-M", "CDA Lugano", 1e+11,
    "IBV-M", "ARA Chur", 1e+11,
    "IBV-M", "ARA Altenrhein", 1e+11,
    "IBV-M", "ARA Sensetal", 1e+11,
    "IBV-M", "STEP Vidy", 1e+11, # no comparison with catchment-case data yet available
    "IBV-M", "STEP Neuchatel", 1e+11, # no comparison with catchment-case data yet available
    "IBV-M", "ARA Buholz", 1e+11, # no comparison with catchment-case data yet available
    "IBV-M", "ARA Basel/Prorheno", 1e+11, # no comparison with catchment-case data yet available
    "IBV-M", "ARA Region Bern", 1e+11, # no comparison with catchment-case data yet available
    "IBV-M", "ARA Schwyz", 1e+11, # no comparison with catchment-case data yet available
    "IBV-M", "STEP Porrentruy", 1e+11, # no comparison with catchment-case data yet available
    # RSV-N: all values are assumed to be comparable to IAV:
    # As RSV is not notifiable, we have no confirmed case data for it
    #   - the number of detected cases of RSV in the 2022/2023 peaked around 
    #     400 cases per week in Switzerland
    #   - the confirmed number of cases for influenza peaked around 4000 in the 
    #     same season (it was mostly an IAV wave)
    #   - the wastewater concentrations had trajectories at comparable levels 
    #     for RSV and IAV-M (both peaking between 100-200 gc/mlWW)
    #   - the RSV cases are likely much more underreported than IAV
    #   - if we assume that underreporting is 10 times stronger than for IAV, 
    #     we would get the same load_per_case
    "RSV-N", "ARA Werdhoelzli", 1e+11,
    "RSV-N", "ARA Zuchwil", 1e+11,
    "RSV-N", "STEP Aire", 1e+11,
    "RSV-N", "CDA Lugano", 1e+10,
    "RSV-N", "ARA Chur", 1e+11,
    "RSV-N", "ARA Altenrhein", 1e+11,
    "RSV-N", "ARA Sensetal", 1e+10,
    "RSV-N", "STEP Vidy", 1e+11,
    "RSV-N", "STEP Neuchatel", 1e+11,
    "RSV-N", "ARA Buholz", 1e+11,
    "RSV-N", "ARA Basel/Prorheno", 1e+11,
    "RSV-N", "ARA Region Bern", 1e+11,
    "RSV-N", "ARA Schwyz", 1e+11,
    "RSV-N", "STEP Porrentruy", 1e+11
  )
  sel <- dt |> dplyr::filter(wwtp_ == wwtp, target_ == target)
  if (nrow(sel)==0) {
    stop(paste("Load per case for", wwtp, target, "not found."))
  } else {
    return(sel$load_per_case[1] * reporting_prop * multiplier)
  }
}