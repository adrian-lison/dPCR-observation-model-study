EpiSewer_job_targets <- list(
  tar_target(
    job_EpiSewer,
    {
      generation_dist <- get_generation_dist(
        target = data_selection_EpiSewer$target_select,
        shift_mean = sensitivity_generation[[1]]$shift_mean,
        shift_sd = sensitivity_generation[[1]]$shift_sd
      )
      
      incubation_dist <- get_incubation_dist(
        target = data_selection_EpiSewer$target_select,
        shift_mean = sensitivity_incubation[[1]]$shift_mean,
        shift_sd = sensitivity_incubation[[1]]$shift_sd
      )
      
      shedding_dist <- get_shedding_dist(
        target = data_selection_EpiSewer$target_select,
        shift_mean = sensitivity_shedding[[1]]$shift_mean,
        shift_cv = sensitivity_shedding[[1]]$shift_cv
      )
      
      load_per_case <- get_load_per_case(
        wwtp = data_selection_EpiSewer$wwtp_select,
        target = data_selection_EpiSewer$target_select,
        multiplier = sensitivity_load_per_case[[1]]$multiplier
      )
      
      selection <- c(
        data_selection_EpiSewer,
        list(
          generation_dist = generation_dist,
          sensitivity_generation = sensitivity_generation,
          incubation_dist = incubation_dist,
          sensitivity_incubation = sensitivity_incubation,
          shedding_dist = shedding_dist,
          sensitivity_shedding = sensitivity_shedding,
          load_per_case = load_per_case,
          sensitivity_load_per_case = sensitivity_load_per_case,
          subsampling = subsampling,
          module_measurements = names(module_measurements[1]),
          module_sampling = names(module_sampling[1]),
          module_sewage = names(module_sewage[1]),
          module_shedding = names(module_shedding[1]),
          module_infections = names(module_infections[1]),
          module_forecast = names(module_forecast[1])
        )
      )
      
      if (aggregate_data[1]) {
        selected_measurements <- data_PCR_agg_select |> 
          dplyr::filter(sapply(date,subsampling[[1]]$subsampling_f), !is_outlier)
      } else {
        selected_measurements <- data_PCR_select |> 
          dplyr::filter(sapply(date,subsampling[[1]]$subsampling_f), !is_outlier)
      }
      
      inputs <- list(
        PCR_select = data_PCR_select,
        PCR_agg_select = data_PCR_agg_select,
        flow_select = data_flow_select,
        selected_measurements = selected_measurements
      )
      
      return(tryCatch(
        {
          ww_res <- EpiSewer::EpiSewer(
            data = EpiSewer::sewer_data(
              measurements = selected_measurements,
              flows = data_flow_select
            ),
            assumptions = EpiSewer::sewer_assumptions(
              generation_dist = generation_dist,
              incubation_dist = incubation_dist,
              shedding_dist = shedding_dist,
              shedding_reference = "symptom_onset",
              load_per_case = load_per_case
            ),
            measurements = module_measurements[[1]],
            sampling = module_sampling[[1]],
            sewage = module_sewage[[1]],
            shedding = module_shedding[[1]],
            infections = module_infections[[1]],
            forecast = module_forecast[[1]],
            fit_opts = fit_opts[[1]],
            results_opts = results_opts[[1]],
            run_fit = FALSE
          )
          ww_res$job$job_dir <- file.path(basename(tar_path_store()), "results")
          ww_res$job$job_target_name <- tar_name()
          ww_res$job$job_name <- paste0(basename(tar_path_store()),"_EpiSewer_", digest::digest(
            EpiSewer::get_checksums(ww_res$job), algo = "md5"
          ))
          ww_res$job$selection <- selection
          ww_res$inputs <- inputs
          return(ww_res)
        },
        error = function(e) {
          return(list(
            error = e,
            job = list(selection = selection),
            inputs = c(inputs, list(
              module_measurements = module_measurements,
              module_sampling = module_sampling,
              module_sewage = module_sewage,
              module_shedding = module_shedding,
              module_infections = module_infections,
              module_forecast = module_forecast,
              fit_opts = fit_opts,
              results_opts = results_opts
            ))
          ))
        }
      ))
    },
    pattern = cross(
      # data
      map(
        data_PCR_select,
        data_PCR_agg_select,
        data_flow_select,
        data_selection_EpiSewer
      ),
      # sensitivity analyses
      sensitivity_generation,
      sensitivity_incubation,
      sensitivity_shedding,
      sensitivity_load_per_case,
      # subsampling
      subsampling,
      # modules
      module_measurements,
      module_sampling,
      module_sewage,
      module_shedding,
      module_infections,
      module_forecast,
      # options
      fit_opts,
      results_opts,
      aggregate_data
    ),
    iteration = "list"
  ),
  tar_target(
    jobs_EpiSewer_valid,
    job_EpiSewer[sapply(job_EpiSewer, function(x) !"error" %in% names(x))],
    deployment = "main"
  ),
  tar_target(
    jobs_EpiSewer_invalid,
    job_EpiSewer[sapply(job_EpiSewer, function(x) "error" %in% names(x))],
    deployment = "main"
  ),
  tar_target(
    job_EpiSewer_resultpath,
    target_job_resultpath(jobs_EpiSewer_valid[[1]]),
    pattern = map(jobs_EpiSewer_valid),
    iteration = "list",
    format = "file",
    deployment = "main"
  ),
  tar_target(
    job_EpiSewer_result,
    readRDS(job_EpiSewer_resultpath),
    pattern = map(job_EpiSewer_resultpath),
    iteration = "list"
  ),
  tar_target(
    job_EpiSewer_UpToDate,
    target_job_UpToDate(jobs_EpiSewer_valid[[1]]$job, job_EpiSewer_result),
    pattern = map(jobs_EpiSewer_valid, job_EpiSewer_result),
    iteration = "vector",
    deployment = "main"
  ),
  tar_target(
    job_EpiSewer_failed,
    job_EpiSewer_UpToDate && class(job_EpiSewer_result$summary) == "try-error",
    pattern = map(job_EpiSewer_UpToDate, job_EpiSewer_result),
    iteration = "vector",
    deployment = "main"
  ),
  tar_target(
    job_EpiSewer_state,
    target_job_state(job_EpiSewer_result),
    pattern = map(job_EpiSewer_result),
    iteration = "vector",
    deployment = "main"
  )
)
