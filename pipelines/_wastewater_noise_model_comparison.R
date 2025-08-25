# user settings
dry_run <- readRDS(file.path(targets::tar_path_store(), "settings/dry_run.rds"))
run_cluster <- readRDS(file.path(targets::tar_path_store(), "settings/run_cluster.rds"))

# Define custom functions and other global objects.
# This is where you write source(\"R/functions.R\")
# if you keep your functions in external scripts.
source("code/local_config.R")
source("code/utils_dists.R")

source("code/pipeline/utils_pipeline.R")
source("code/pipeline/functions_pipeline.R")
source("code/pipeline/components_pipeline_dPCR_paper.R")
source("data/assumptions/epi_params_default.R")

source("code/pipeline/run_local_pipeline.R")

if (run_cluster) {
  source("code/utils_euler_slurm.R")
  source("code/pipeline/run_cluster_pipeline.R")
}

# target-specific options
tar_option_set(
  packages = c(
    "dplyr", "tidyr", "readr", "EpiSewer",
    "data.table", "stringr", "targets", "ssh"),
  workspace_on_error = TRUE, controller = crew_controller_local(workers = 4)
  )

## selections ----
selection_targets <- list(
  tar_target(
    wwtp_select,
    c("ARA Werdhoelzli")
  ),
  tar_target(
    assay_select,
    list(c("RESPV4", "RESPV6"))
  ),
  tar_target(
    target_select,
    c("IAV-M")
  ),
  tar_target(
    date_select,
    list(
      # 2022
      c(from = as.Date("2022-09-15"), to = as.Date("2022-12-29")),
      c(from = as.Date("2022-09-15"), to = as.Date("2023-04-15")),
      # 2023
      c(from = as.Date("2023-10-15"), to = as.Date("2024-02-03")),
      c(from = as.Date("2023-10-15"), to = as.Date("2024-04-15")),
      # 2024
      c(from = as.Date("2024-10-15"), to = as.Date("2025-02-07")),
      c(from = as.Date("2024-10-15"), to = as.Date("2025-04-15"))
      )
  ),
  tar_target(
    load_per_case_data,
    NULL
  )
)

## shedding info
get_shedding_dist_info <- get_shedding_dist_info_fixed

## sensitivity analyses of assumptions ----
sensitivity_targets <- list(
  tar_target(
    sensitivity_generation,
    list(
      list(shift_mean = 0, shift_sd = 0)
    )
  ),
  tar_target(
    sensitivity_incubation,
    list(
      list(shift_mean = 0, shift_sd = 0)
    )
  ),
  tar_target(
    sensitivity_shedding,
    list(
      list(shift_mean = 0, shift_cv = 0)
    )
  ),
  tar_target(
    sensitivity_load_per_case,
    list(
      list(multiplier = 1)
    ),
  )
)

## subsampling ----
subsampling_targets <- list(
  tar_target(
    subsampling,
    list(
      list(
        type = "Every day",
        subtype = "Every day",
        subsampling_f = function(dates) {rep(TRUE,length(dates))}
      )
      )
    )
  )

## modeling modules ----
modeling_targets <- list(
  tar_target(
    module_measurements,
    list(
      normal = model_measurements(
        concentrations = concentrations_observe(
          concentration_col = "gc_per_mlww",
          n_averaged_col = "n_reps",
          distribution = "normal"
        ),
        noise = noise_estimate_constant_var(
          cv_prior_mu = 0, cv_prior_sigma = 1, warn = FALSE
          ),
        LOD = LOD_none()
      ),
      lognormal = model_measurements(
        concentrations = concentrations_observe(
          concentration_col = "gc_per_mlww",
          n_averaged_col = "n_reps",
          distribution = "lognormal"
        ),
        noise = noise_estimate(
          cv_prior_mu = 0, cv_prior_sigma = 1
        ),
        LOD = LOD_none()
      ),
      dPCR = model_measurements(
        concentrations = concentrations_observe(
          concentration_col = "gc_per_mlww",
          n_averaged_col = "n_reps",
          distribution = "gamma"
        ),
        noise = noise_estimate_dPCR(
          cv_prior_mu = 0, cv_prior_sigma = 1,
          max_partitions_prior_lower = 10000,
          max_partitions_prior_upper = 50000,
          partition_loss_mean_prior_lower = 0.01,
          partition_loss_mean_prior_upper = 0.3,
          partition_loss_variation_prior_lower = 0.5,
          partition_loss_variation_prior_upper = 2,
          partition_loss_max = 0.5,
          volume_scaled_prior_mu = 1e-5,
          volume_scaled_prior_sigma = 4e-5,
          prePCR_noise_type = "log-normal"
        ),
        LOD = LOD_estimate_dPCR()
      ),
      dPCR_informed = model_measurements(
        concentrations = concentrations_observe(
          concentration_col = "gc_per_mlww",
          n_averaged_col = "n_reps",
          total_partitions_col = "total_droplets",
          distribution = "gamma"
        ),
        noise = noise_estimate_dPCR(
          cv_prior_mu = 0, cv_prior_sigma = 1,
          total_partitions_observe = TRUE,
          volume_scaled_prior_mu = 1.73e-5,
          volume_scaled_prior_sigma = 0,
          prePCR_noise_type = "log-normal"
        ),
        LOD = LOD_estimate_dPCR()
      ),
      dPCR_binomial = model_measurements(
        concentrations = concentrations_observe_partitions(
          concentration_col = "gc_per_mlww",
          n_averaged_col = "n_reps",
          total_partitions_col = "total_droplets",
          positive_partitions_col = "positive_droplets"
        ),
        noise = noise_estimate_dPCR(
          cv_prior_mu = 0, cv_prior_sigma = 1,
          total_partitions_observe = TRUE,
          volume_scaled_prior_mu = 1.73e-5,
          volume_scaled_prior_sigma = 0,
          prePCR_noise_type = "log-normal"
        ),
        LOD = LOD_estimate_dPCR()
      )
    )
  ),
  tar_target(
    module_sampling,
    list(
      model_sampling(
        sample_effects = sample_effects_none()
      )
    )
  ),
  tar_target(
    module_sewage,
    list(
      model_sewage(
        flows = flows_observe(),
        residence_dist = residence_dist_assume(
          residence_dist = c(1)
          )
      )
    )
  ),
  tar_target(
    module_shedding,
    list(
      model_shedding(
        incubation_dist = incubation_dist_assume(),
        shedding_dist = shedding_dist_assume(),
        load_per_case = load_per_case_assume(),
        load_variation = load_variation_none()
      )
    )
  ),
  tar_target(
    module_infections,
    list(
      rw = model_infections(
        generation_dist = generation_dist_assume(),
        R = R_estimate_rw(),
        seeding = seeding_estimate_growth(),
        infection_noise = infection_noise_estimate()
      )
    )
  ),
  tar_target(
    module_forecast,
    list(
     model_forecast(
       horizon = horizon_none()
     )
    )
  ),
  tar_target(
    data_handling_opts,
    list(list(aggregate_data = TRUE, remove_outliers = TRUE))
  )
)

## sampling options ----
### sampler ----
option_targets <- list(
  tar_target(
    fit_opts,
    list(
      set_fit_opts(
        sampler = sampler_stan_mcmc(
          chains = 4, iter_warmup = 1000, iter_sampling = 1000
          )
      )
    )
  ),
  tar_target(
    results_opts,
    list(
      set_results_opts(
        fitted = TRUE, summary_intervals = c(0.5, 0.8, 0.95)
      )
    )
  )
)

### run local vs euler ----
if (run_cluster) {
  run_targets <- euler_targets
  euler_up <- !all(is.na(pingr::ping("euler.ethz.ch", count = 2))) 
  if (euler_up) { try(sync_from_euler(target_path = "pipelines")) }
} else {
  run_targets <- local_run_targets[c("EpiSewer")]
}

# pipeline ----
c(
  selection_targets,
  data_targets,
  subsampling_targets,
  sensitivity_targets,
  modeling_targets,
  option_targets,
  EpiSewer_job_targets,
  run_targets
)
