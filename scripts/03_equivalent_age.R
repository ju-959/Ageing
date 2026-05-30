library(dplyr)
library(readr)
library(stringr)
library(tidyr)

root_dir <- normalizePath(Sys.getenv("AGEING_ROOT", unset = getwd()), mustWork = TRUE)
demo_mode <- identical(Sys.getenv("AGEING_DEMO", unset = "FALSE"), "TRUE")

input_dir <- if (demo_mode) file.path(root_dir, "outputs", "demo") else file.path(root_dir, "outputs", "analysis")
output_dir <- input_dir
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

overall_input_file <- file.path(input_dir, "total_age_related_daly_by_age.csv")
system_input_file <- file.path(input_dir, "system_age_related_daly_by_age.csv")

overall_output_file <- file.path(output_dir, "equivalent_age_overall_1990_2021.csv")
system_output_file <- file.path(output_dir, "equivalent_age_system_1990_2021.csv")
combined_output_file <- file.path(output_dir, "equivalent_age_1990_2021.csv")

age_groups_for_analysis <- c(
  "25 to 29", "30 to 34", "35 to 39", "40 to 44", "45 to 49",
  "50 to 54", "55 to 59", "60 to 64", "65 to 69", "70 to 74",
  "75 to 79", "80 to 84", "85 to 89", "90 to 94"
)

reference_age_groups <- c("60 to 64", "65 to 69")
benchmark_age <- 65

compute_midpoint_age <- function(age_group) {
  values <- as.numeric(str_extract_all(age_group, "\\d+")[[1]])
  mean(values, na.rm = TRUE)
}

interpolate_equivalent_age <- function(prev_rate, next_rate, prev_age, next_age, benchmark) {
  if (is.na(prev_rate) || is.na(next_rate) || is.na(benchmark)) return(NA_real_)
  if (abs(next_rate - prev_rate) < 1e-9) return(mean(c(prev_age, next_age)))
  prev_age + (benchmark - prev_rate) * (next_age - prev_age) / (next_rate - prev_rate)
}

as_display_value <- function(value, lower_boundary, upper_boundary) {
  dplyr::case_when(
    lower_boundary ~ "<25",
    upper_boundary ~ ">94",
    !is.na(value) ~ sprintf("%.2f", value),
    TRUE ~ "Cannot Determine"
  )
}

as_numeric_value <- function(value, lower_boundary, upper_boundary) {
  dplyr::case_when(
    lower_boundary ~ 25,
    upper_boundary ~ 94,
    !is.na(value) ~ value,
    TRUE ~ NA_real_
  )
}

prepare_rate_data <- function(path, analysis_level) {
  dat <- read_csv(path, show_col_types = FALSE) %>%
    filter(
      measure == "DALYs (Disability-Adjusted Life Years)",
      metric == "Rate",
      age %in% age_groups_for_analysis
    )

  if (!"sex" %in% names(dat)) dat$sex <- "Both"
  if (!"category" %in% names(dat)) dat$category <- "Overall"

  dat %>%
    mutate(
      analysis_level = analysis_level,
      category = if_else(is.na(category), "Overall", as.character(category)),
      age = factor(age, levels = age_groups_for_analysis, ordered = TRUE),
      year = as.integer(year),
      val = as.numeric(val),
      lower = as.numeric(lower),
      upper = as.numeric(upper)
    )
}

get_benchmarks <- function(rate_data) {
  benchmarks <- rate_data %>%
    filter(location == "Global", sex == "Both", year == 2021, age %in% reference_age_groups) %>%
    group_by(analysis_level, category) %>%
    summarise(
      benchmark_rate = mean(val, na.rm = TRUE),
      benchmark_lower = mean(lower, na.rm = TRUE),
      benchmark_upper = mean(upper, na.rm = TRUE),
      .groups = "drop"
    )

  if (nrow(benchmarks) == 0 || any(is.na(benchmarks$benchmark_rate))) {
    stop("The Global Both 2021 age-65 benchmark could not be derived.", call. = FALSE)
  }

  benchmarks
}

calculate_equivalent_age <- function(rate_data) {
  benchmarks <- get_benchmarks(rate_data)
  grouping_columns <- c("analysis_level", "category", "location", "sex", "year")

  analysis_units <- rate_data %>%
    distinct(across(all_of(grouping_columns)))

  interpolatable_data <- rate_data %>%
    left_join(benchmarks, by = c("analysis_level", "category")) %>%
    group_by(across(all_of(grouping_columns))) %>%
    arrange(age, .by_group = TRUE) %>%
    mutate(
      prev_group = as.character(age),
      prev_rate = val,
      next_group = lead(as.character(age)),
      next_rate = lead(val)
    ) %>%
    filter(
      !is.na(next_rate),
      (prev_rate <= benchmark_rate & next_rate >= benchmark_rate) |
        (prev_rate >= benchmark_rate & next_rate <= benchmark_rate)
    ) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(
      prev_midpoint_age = vapply(prev_group, compute_midpoint_age, numeric(1)),
      next_midpoint_age = vapply(next_group, compute_midpoint_age, numeric(1)),
      equivalent_age_raw = mapply(
        interpolate_equivalent_age,
        prev_rate,
        next_rate,
        prev_midpoint_age,
        next_midpoint_age,
        benchmark_rate
      ),
      lower_bound_raw = mapply(
        interpolate_equivalent_age,
        prev_rate,
        next_rate,
        prev_midpoint_age,
        next_midpoint_age,
        ifelse(next_rate > prev_rate, benchmark_lower, benchmark_upper)
      ),
      upper_bound_raw = mapply(
        interpolate_equivalent_age,
        prev_rate,
        next_rate,
        prev_midpoint_age,
        next_midpoint_age,
        ifelse(next_rate > prev_rate, benchmark_upper, benchmark_lower)
      )
    ) %>%
    select(
      all_of(grouping_columns), prev_group, next_group, prev_rate, next_rate,
      equivalent_age_raw, lower_bound_raw, upper_bound_raw
    )

  boundary_data <- rate_data %>%
    filter(age %in% c("25 to 29", "90 to 94")) %>%
    select(all_of(grouping_columns), age, val) %>%
    pivot_wider(
      id_cols = all_of(grouping_columns),
      names_from = age,
      values_from = val,
      values_fn = mean
    ) %>%
    rename(
      rate_at_start = `25 to 29`,
      rate_at_end = `90 to 94`
    )

  analysis_units %>%
    left_join(interpolatable_data, by = grouping_columns) %>%
    left_join(boundary_data, by = grouping_columns) %>%
    left_join(benchmarks, by = c("analysis_level", "category")) %>%
    mutate(
      lower_boundary = !is.na(rate_at_start) & rate_at_start > benchmark_rate,
      upper_boundary = !is.na(rate_at_end) & rate_at_end < benchmark_rate,
      anchor_benchmark = location == "Global" & sex == "Both" & year == 2021,
      equivalent_age = if_else(
        anchor_benchmark,
        sprintf("%.2f", benchmark_age),
        as_display_value(equivalent_age_raw, lower_boundary, upper_boundary)
      ),
      lower_bound = if_else(
        anchor_benchmark,
        sprintf("%.2f", benchmark_age),
        as_display_value(lower_bound_raw, lower_boundary, upper_boundary)
      ),
      upper_bound = if_else(
        anchor_benchmark,
        sprintf("%.2f", benchmark_age),
        as_display_value(upper_bound_raw, lower_boundary, upper_boundary)
      ),
      equivalent_age_numeric = if_else(
        anchor_benchmark,
        as.numeric(benchmark_age),
        as_numeric_value(equivalent_age_raw, lower_boundary, upper_boundary)
      ),
      lower_bound_numeric = if_else(
        anchor_benchmark,
        as.numeric(benchmark_age),
        as_numeric_value(lower_bound_raw, lower_boundary, upper_boundary)
      ),
      upper_bound_numeric = if_else(
        anchor_benchmark,
        as.numeric(benchmark_age),
        as_numeric_value(upper_bound_raw, lower_boundary, upper_boundary)
      )
    ) %>%
    select(
      analysis_level, location, sex, category, year,
      benchmark_rate, benchmark_lower, benchmark_upper,
      prev_group, next_group, prev_rate, next_rate,
      equivalent_age, lower_bound, upper_bound,
      equivalent_age_numeric, lower_bound_numeric, upper_bound_numeric
    ) %>%
    arrange(analysis_level, location, sex, category, year)
}

overall_rate_data <- prepare_rate_data(overall_input_file, analysis_level = "overall")
overall_results <- calculate_equivalent_age(overall_rate_data)
write_csv(overall_results, overall_output_file)

all_results <- list(overall_results)

if (file.exists(system_input_file)) {
  system_rate_data <- prepare_rate_data(system_input_file, analysis_level = "system")
  system_results <- calculate_equivalent_age(system_rate_data)
  write_csv(system_results, system_output_file)
  all_results <- c(all_results, list(system_results))
}

write_csv(bind_rows(all_results), combined_output_file)
