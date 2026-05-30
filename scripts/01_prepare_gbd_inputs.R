library(readr)
library(dplyr)

root_dir <- normalizePath(Sys.getenv("AGEING_ROOT", unset = getwd()), mustWork = TRUE)
demo_mode <- identical(Sys.getenv("AGEING_DEMO", unset = "FALSE"), "TRUE")
processed_dir <- if (demo_mode) file.path(root_dir, "data", "demo") else file.path(root_dir, "data", "processed")
dir.create(processed_dir, showWarnings = FALSE, recursive = TRUE)

resolve_input_file <- function(directory, file) {
  csv_path <- file.path(directory, file)
  gz_path <- paste0(csv_path, ".gz")

  if (file.exists(csv_path)) {
    csv_path
  } else if (file.exists(gz_path)) {
    gz_path
  } else {
    csv_path
  }
}

expected_files <- tibble::tribble(
  ~file, ~purpose,
  "incidence_age_pattern.csv", "Incidence rates used to screen age-related causes",
  "prevalence_age_pattern.csv", "Prevalence rates used when incidence was unavailable",
  "daly_age_specific.csv", "Age-specific DALY rates and numbers for selected detailed causes",
  "disease_system_map.csv", "Mapping from detailed causes to disease systems",
  "population_history_2021.csv", "Historical GBD 2021 population by location, sex, year and age",
  "population_projection_who.csv", "Projected population denominators used by the projection analysis"
)

required_columns <- c("measure", "location", "sex", "age", "year", "cause", "metric", "val", "lower", "upper")

check_table <- function(path, required = required_columns) {
  dat <- suppressMessages(read_csv(path, show_col_types = FALSE, n_max = 20))
  missing <- setdiff(required, names(dat))
  tibble::tibble(
    file = basename(path),
    exists = TRUE,
    missing_columns = if (length(missing) == 0) "" else paste(missing, collapse = "; ")
  )
}

available <- expected_files %>%
  mutate(path = vapply(file, resolve_input_file, character(1), directory = processed_dir)) %>%
  rowwise() %>%
  mutate(exists = file.exists(path)) %>%
  ungroup()

checks <- lapply(seq_len(nrow(available)), function(i) {
  if (!available$exists[i]) {
    return(tibble::tibble(file = available$file[i], exists = FALSE, missing_columns = "file not found"))
  }
  if (available$file[i] == "disease_system_map.csv") {
    return(check_table(available$path[i], required = c("cause", "category")))
  }
  if (grepl("^population_", available$file[i])) {
    return(check_table(available$path[i], required = c("location", "sex", "year", "age", "val")))
  }
  check_table(available$path[i])
}) %>%
  bind_rows() %>%
  left_join(expected_files, by = "file")

location_group_path <- file.path(processed_dir, "projection_location_groups.csv")
if (file.exists(location_group_path)) {
  location_group_check <- check_table(
    location_group_path,
    required = c("location", "aggregate_location")
  ) %>%
    mutate(purpose = "Optional mapping from projected country locations to aggregate SDI or WHO-region outputs")

  checks <- bind_rows(checks, location_group_check)
}

write_csv(checks, file.path(processed_dir, "input_file_check.csv"))

if (any(!checks$exists) || any(checks$missing_columns != "")) {
  message("Some processed inputs are missing or incomplete. See data/README.md and data/processed/input_file_check.csv.")
} else {
  message("All processed input files are present with the expected columns.")
}
