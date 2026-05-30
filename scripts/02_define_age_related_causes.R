library(dplyr)
library(readr)
library(stringr)
library(tidyr)

root_dir <- normalizePath(Sys.getenv("AGEING_ROOT", unset = getwd()), mustWork = TRUE)
demo_mode <- identical(Sys.getenv("AGEING_DEMO", unset = "FALSE"), "TRUE")

input_dir <- if (demo_mode) file.path(root_dir, "data", "demo") else file.path(root_dir, "data", "processed")
output_dir <- if (demo_mode) file.path(root_dir, "outputs", "demo") else file.path(root_dir, "outputs", "analysis")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

incidence_input_file <- file.path(input_dir, "incidence_age_pattern.csv")
prevalence_input_file <- file.path(input_dir, "prevalence_age_pattern.csv")
daly_input_file <- file.path(input_dir, "daly_age_specific.csv")
disease_system_map_file <- file.path(input_dir, "disease_system_map.csv")

age_related_disease_output <- file.path(output_dir, "age_related_disease_list.csv")
total_age_related_output <- file.path(output_dir, "total_age_related_daly_by_age.csv")
system_age_related_output <- file.path(output_dir, "system_age_related_daly_by_age.csv")

n_draws <- if (demo_mode) 25 else 1000
set.seed(123)

age_groups_for_definition <- c(
  "25 to 29", "30 to 34", "35 to 39", "40 to 44", "45 to 49",
  "50 to 54", "55 to 59", "60 to 64", "65 to 69", "70 to 74",
  "75 to 79", "80 to 84"
)

age_groups_for_aggregation <- c(
  "25 to 29", "30 to 34", "35 to 39", "40 to 44", "45 to 49",
  "50 to 54", "55 to 59", "60 to 64", "65 to 69", "70 to 74",
  "75 to 79", "80 to 84", "85 to 89", "90 to 94"
)

age_midpoints <- c(
  "25 to 29" = 27, "30 to 34" = 32, "35 to 39" = 37, "40 to 44" = 42,
  "45 to 49" = 47, "50 to 54" = 52, "55 to 59" = 57, "60 to 64" = 62,
  "65 to 69" = 67, "70 to 74" = 72, "75 to 79" = 77, "80 to 84" = 82,
  "85 to 89" = 87, "90 to 94" = 92, "95 plus" = 97
)

simulate_positive_rate <- function(mean_rate, se, n_draws) {
  if (!is.finite(mean_rate) || !is.finite(se) || mean_rate <= 0 || se <= 0) {
    return(rep(pmax(mean_rate, 0), n_draws))
  }
  meanlog <- log(mean_rate^2 / sqrt(se^2 + mean_rate^2))
  sdlog <- sqrt(log(1 + (se^2 / mean_rate^2)))
  rlnorm(n_draws, meanlog = meanlog, sdlog = sdlog)
}

screen_age_related_diseases <- function(input_data) {
  screened_input <- input_data %>%
    mutate(
      age_mid = unname(age_midpoints[age]),
      se = (upper - lower) / (2 * 1.96)
    ) %>%
    filter(
      sex == "Both",
      location == "Global",
      metric == "Rate",
      age %in% age_groups_for_definition
    ) %>%
    arrange(cause, age_mid)

  if (nrow(screened_input) == 0) {
    return(tibble(
      cause = character(),
      prop_beta1_positive = numeric(),
      prop_beta2_positive = numeric(),
      passes_step1 = logical(),
      passes_step2 = logical(),
      age_related = logical()
    ))
  }

  incomplete_causes <- screened_input %>%
    count(cause, name = "n_age_groups") %>%
    filter(n_age_groups != length(age_groups_for_definition))

  if (nrow(incomplete_causes) > 0) {
    stop(
      "Age-pattern screening requires complete 25-29 through 80-84 age panels for every cause.",
      call. = FALSE
    )
  }

  disease_results <- screened_input %>%
    group_by(cause) %>%
    group_modify(~ {
      draw_summary <- lapply(seq_len(n_draws), function(draw_id) {
        simulated_values <- mapply(
          simulate_positive_rate,
          .x$val,
          .x$se,
          MoreArgs = list(n_draws = 1)
        )
        draw_data <- mutate(.x, simulated_rate = as.numeric(simulated_values))
        beta_linear <- coef(lm(simulated_rate ~ age_mid, data = draw_data))[["age_mid"]]
        beta_quadratic <- NA_real_
        if (is.finite(beta_linear) && beta_linear > 0) {
          beta_quadratic <- coef(lm(simulated_rate ~ age_mid + I(age_mid^2), data = draw_data))[["I(age_mid^2)"]]
        }
        tibble(beta_linear = beta_linear, beta_quadratic = beta_quadratic)
      }) %>%
        bind_rows()

      tibble(
        prop_beta1_positive = mean(draw_summary$beta_linear > 0, na.rm = TRUE),
        prop_beta2_positive = mean(draw_summary$beta_quadratic > 0, na.rm = TRUE),
        passes_step1 = prop_beta1_positive >= 0.95,
        passes_step2 = ifelse(is.na(prop_beta2_positive), FALSE, prop_beta2_positive >= 0.95),
        age_related = passes_step1 & passes_step2
      )
    }) %>%
    ungroup()

  disease_results
}

aggregate_with_uncertainty <- function(input_data, grouping_columns, n_draws = 1000) {
  aggregated_draws <- lapply(seq_len(n_draws), function(draw_id) {
    input_data %>%
      mutate(
        se = (upper - lower) / (2 * 1.96),
        draw_value = pmax(0, rnorm(n(), mean = val, sd = se))
      ) %>%
      group_by(across(all_of(grouping_columns))) %>%
      summarise(total_value = sum(draw_value, na.rm = TRUE), .groups = "drop") %>%
      mutate(draw = draw_id)
  }) %>%
    bind_rows()

  aggregated_draws %>%
    group_by(across(all_of(grouping_columns))) %>%
    summarise(
      val = mean(total_value, na.rm = TRUE),
      lower = quantile(total_value, 0.025, na.rm = TRUE),
      upper = quantile(total_value, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
}

incidence_data <- read_csv(incidence_input_file, show_col_types = FALSE) %>%
  filter(measure == "Incidence", year == 2021)

prevalence_data <- read_csv(prevalence_input_file, show_col_types = FALSE) %>%
  filter(measure == "Prevalence", year == 2021)

# Diseases with available incidence data are screened using incidence. Diseases
# without incidence data are screened using prevalence.
incidence_screen <- screen_age_related_diseases(incidence_data)

missing_incidence_causes <- setdiff(
  unique(prevalence_data$cause),
  unique(incidence_data$cause)
)

prevalence_screen <- prevalence_data %>%
  filter(cause %in% missing_incidence_causes) %>%
  screen_age_related_diseases()

age_related_disease_list <- bind_rows(
  incidence_screen %>%
    filter(age_related) %>%
    transmute(cause, source_measure = "Incidence"),
  prevalence_screen %>%
    filter(age_related) %>%
    transmute(cause, source_measure = "Prevalence")
) %>%
  distinct(cause, .keep_all = TRUE) %>%
  arrange(cause)

if (file.exists(disease_system_map_file)) {
  disease_system_map <- read_csv(disease_system_map_file, show_col_types = FALSE) %>%
    select(cause, category)
  age_related_disease_list <- age_related_disease_list %>%
    left_join(disease_system_map, by = "cause")
} else {
  age_related_disease_list <- age_related_disease_list %>%
    mutate(category = NA_character_)
}

write_csv(age_related_disease_list, age_related_disease_output)

daly_data <- read_csv(daly_input_file, show_col_types = FALSE) %>%
  filter(
    measure == "DALYs (Disability-Adjusted Life Years)",
    metric %in% c("Rate", "Number"),
    age %in% age_groups_for_aggregation
  ) %>%
  semi_join(age_related_disease_list, by = "cause")

total_age_related_daly <- aggregate_with_uncertainty(
  input_data = daly_data,
  grouping_columns = c("measure", "metric", "location", "year", "age", "sex"),
  n_draws = n_draws
) %>%
  arrange(measure, metric, location, year, age, sex)

write_csv(total_age_related_daly, total_age_related_output)

if (any(!is.na(age_related_disease_list$category))) {
  system_age_related_daly <- daly_data %>%
    left_join(
      age_related_disease_list %>% select(cause, category),
      by = "cause"
    ) %>%
    filter(!is.na(category)) %>%
    aggregate_with_uncertainty(
      grouping_columns = c("measure", "metric", "location", "year", "age", "sex", "category"),
      n_draws = n_draws
    ) %>%
    arrange(category, measure, metric, location, year, age, sex)

  write_csv(system_age_related_daly, system_age_related_output)
}
