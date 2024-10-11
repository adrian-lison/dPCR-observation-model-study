# Data targets ----
data_targets <- list(
  ## overall selection ----
  tar_target(
    data_selection_EpiSewer_all,
    {
      list(
        wwtp_select = wwtp_select,
        assay_select = assay_select,
        target_select = target_select,
        date_select = date_select
      )
    },
    pattern = cross(wwtp_select, assay_select, date_select, target_select),
    iteration = "list"
  ),
  tar_target(
    data_selection_EpiSewer,
    data_selection_EpiSewer_all[!data_PCR_duplicated],
    iteration = "list"
  ),
  ## dPCR measurements ----
  tar_target(
    file_ww_data,
    here::here("data", "ww_data", "dPCR_time_series_IAV_Zurich_Geneva.csv"),
    format = "file"
  ),
  tar_target(
    ww_data,
    readr::read_csv(
      file_ww_data,
      skip = 1,
      col_names = c(
        "wwtp", "date", "target",
        "assay", "gc_per_lww",
        "replicate_id", "dilution", "total_droplets"
      ),
      col_types = readr::cols(
        date = readr::col_date(format = "%Y-%m-%dT%H:%M:%SZ"),
        gc_per_lww = readr::col_double()
      ),
      show_col_types = F
    )
  ),
  tar_target(
    outlier_df_file,
    file.path(here::here("data", "ww_data", "outliers_manual.csv")),
    format = "file"
  ),
  tar_target(
    outlier_df,
    readr::read_csv(outlier_df_file, show_col_types = FALSE)
  ),
  tar_target(
    data_PCR_select_all,
    {
    measurements <- ww_data |> dplyr::filter(
      wwtp == wwtp_select,
      assay == assay_select,
      target == target_select,
      date >= date_select[[1]][["from"]], date <= date_select[[1]][["to"]]
    ) |> 
      dplyr::mutate(gc_per_mlww = gc_per_lww / 1000) |> 
      dplyr::select(target, wwtp, date, gc_per_mlww, dilution, total_droplets) |> 
      dplyr::group_by(target, wwtp, date) |>
      dplyr::mutate(replicate_id = 1:dplyr::n())
    data.table::setDT(measurements)
    measurements[, is_outlier := FALSE]
    measurements <- mark_outlier_manual(
      measurements,
      outlier_df,
      target_col = "target",
      wwtp_col = "wwtp",
      date_col = "date"
    )
    data.table::setDT(measurements)
    measurements[is.na(dilution) & !is.na(gc_per_mlww), dilution := median(measurements$dilution, na.rm = T)]
    measurements[is.na(total_droplets) & !is.na(gc_per_mlww), total_droplets := median(measurements$total_droplets, na.rm = T)]
    return(measurements)
    },
    pattern = cross(wwtp_select, assay_select, date_select, target_select),
    iteration = "list"
  ),
  tar_target(
    data_PCR_duplicated,
    duplicated(t(sapply(data_PCR_select_all, function(x) x |> slice_max(date))))
  ),
  tar_target(
    data_PCR_select,
    data_PCR_select_all[!data_PCR_duplicated],
    iteration = "list"
  ),
  tar_target(
    data_PCR_agg_select,
    {
      key_cols <- c(
        "target",
        "wwtp",
        "date"
      )
      aggregate_replicates_mean(data_PCR_select, key_cols)
    },
    pattern = map(data_PCR_select),
    iteration = "list"
  ),
  tar_target(
    data_PCR_plot,
    EpiSewer::plot_concentration(
      measurements = data_PCR_agg_select, concentration_col = "gc_per_mlww",
      mark_outliers = TRUE, outlier_col = "is_outlier"
    ),
    pattern = map(data_PCR_agg_select),
    iteration = "list"
  ),
  ## flow ----
  tar_target(
    file_flow_data,
    here::here("data", "ww_data", "flow_time_series_Zurich_Geneva.csv"),
    format = "file"
  ),
  tar_target(
    flow_data,
    readr::read_csv(file_flow_data)
  ),
  tar_target(
    data_flow_select_all,
    flow_data |>
      dplyr::filter(
        wwtp == wwtp_select,
        date >= date_select[[1]][["from"]],
        date <= date_select[[1]][["to"]] + 28 # this allows for up to 4 weeks forecast
      ),
    pattern = cross(wwtp_select, assay_select, date_select, target_select),
    iteration = "list"
  ),
  tar_target(
   data_flow_select,
   data_flow_select_all[!data_PCR_duplicated],
   iteration = "list"
  )
)

## EpiSewer_job_targets ----
source("code/pipeline/EpiSewer_job_targets.R")