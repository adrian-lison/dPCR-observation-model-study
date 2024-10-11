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
source("code/epi_params_default.R")

source("code/pipeline/run_local_pipeline.R")

if (run_cluster) {
  source("code/utils_euler_slurm.R")
  source("code/pipeline/run_cluster_pipeline.R")
}

# target-specific options
tar_option_set(
  packages = c(
    "dplyr", "tidyr", "readr", "EpiSewer",
    "data.table", "stringr", "targets"),
  workspace_on_error = TRUE, controller = crew_controller_local(workers = 4)
  )

## selections ----
selection_targets <- list(
  tar_target(
    wwtp_select,
    c("ARA Werdhoelzli", "STEP Aire")
  ),
  tar_target(
    assay_select,
    c("respv4")
  ),
  tar_target(
    target_select,
    c("IAV-M")
  ),
  tar_target(
    date_select,
    list(
      c(from = as.Date("2022-09-15"), to = as.Date("2022-12-15")),
      c(from = as.Date("2022-09-15"), to = as.Date("2022-12-19")),
      c(from = as.Date("2022-09-15"), to = as.Date("2022-12-22")),
      c(from = as.Date("2022-09-15"), to = as.Date("2022-12-26")),
      c(from = as.Date("2022-09-15"), to = as.Date("2022-12-29")),
      c(from = as.Date("2022-09-15"), to = as.Date("2023-05-01"))
      )
  )
)

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
        subsampling_f = function(date) {TRUE}
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
          replicate_col = "replicate_id",
          distribution = "normal"
        ),
        noise = noise_estimate_constant_var(warn = FALSE),
        LOD = LOD_none()
      ),
      lognormal = model_measurements(
        concentrations = concentrations_observe(
          concentration_col = "gc_per_mlww",
          replicate_col = "replicate_id",
          distribution = "lognormal"
        ),
        noise = noise_estimate(),
        LOD = LOD_none()
      ),
      dPCR = model_measurements(
        concentrations = concentrations_observe(
          concentration_col = "gc_per_mlww",
          replicate_col = "replicate_id",
          distribution = "gamma"
        ),
        noise = noise_estimate_dPCR(
          total_partitions_prior_mu = 20000,
          total_partitions_prior_sigma = 5000,
          partition_variation_prior_mu = 0,
          partition_variation_prior_sigma = 0.05,
          volume_scaled_prior_mu = 1e-5,
          volume_scaled_prior_sigma = 4e-5,
          prePCR_noise_type = "log-normal"
        ),
        LOD = LOD_estimate_dPCR()
      ),
      dPCR_informed = model_measurements(
        concentrations = concentrations_observe(
          concentration_col = "gc_per_mlww",
          replicate_col = "replicate_id",
          total_partitions_col = "total_droplets",
          distribution = "gamma"
        ),
        noise = noise_estimate_dPCR(
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
        residence_dist = residence_dist_assume(residence_dist = c(1))
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
        seeding = seeding_estimate_rw(),
        infection_noise = infection_noise_estimate()
      )
    )
  ),
  tar_target(
    module_forecast,
    list(
     model_forecast(
       horizon = horizon_assume(14)
     )
    )
  ),
  tar_target(
    aggregate_data,
    FALSE
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
