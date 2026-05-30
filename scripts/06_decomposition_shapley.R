library(data.table)
library(dplyr)
library(stringr)
library(tidyr)

# Shapley decomposition with bootstrap confidence intervals.

root_dir <- normalizePath(Sys.getenv("AGEING_ROOT", unset = getwd()), mustWork = TRUE)
demo_mode <- identical(Sys.getenv("AGEING_DEMO", unset = "FALSE"), "TRUE")

input_dir <- if (demo_mode) file.path(root_dir, "outputs", "demo") else file.path(root_dir, "outputs", "analysis")
output_dir <- if (demo_mode) file.path(root_dir, "outputs", "demo") else file.path(root_dir, "outputs", "analysis", "decomposition")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

age_groups <- c(
  "25 to 29", "30 to 34", "35 to 39", "40 to 44", "45 to 49",
  "50 to 54", "55 to 59", "60 to 64", "65 to 69", "70 to 74",
  "75 to 79", "80 to 84", "85 to 89", "90 to 94"
)
reference_ages <- c("60 to 64", "65 to 69")
year0 <- 1990
year1 <- 2021
n_perm <- if (demo_mode) 50 else 1000
n_boot <- if (demo_mode) 25 else 1000
seed <- 2026

mean_age <- function(age_group) {
  mean(as.numeric(str_extract_all(age_group, "\\d+")[[1]]), na.rm = TRUE)
}

age_mid <- setNames(vapply(age_groups, mean_age, numeric(1)), age_groups)

equiv_age_from_curve <- function(curve, reference_rate) {
  y_start <- curve[["25 to 29"]]
  y_end <- curve[["90 to 94"]]
  if (is.na(y_start) || is.na(y_end)) return(NA_real_)
  if (y_start > reference_rate) return(25)
  if (y_end < reference_rate) return(94)

  for (i in seq_len(length(age_groups) - 1)) {
    a1 <- age_groups[i]
    a2 <- age_groups[i + 1]
    y1 <- curve[[a1]]
    y2 <- curve[[a2]]
    if (any(is.na(c(y1, y2)))) next

    crossed <- (y1 <= reference_rate && y2 >= reference_rate) ||
      (y1 >= reference_rate && y2 <= reference_rate)
    if (crossed) {
      if (y2 == y1) return(mean(c(age_mid[[a1]], age_mid[[a2]])))
      return(age_mid[[a1]] + (reference_rate - y1) * (age_mid[[a2]] - age_mid[[a1]]) / (y2 - y1))
    }
  }
  NA_real_
}

shapley_equiv_age <- function(mat0, mat1, reference_rate, n_perm, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  k <- nrow(mat0)
  contribution <- numeric(k)
  names(contribution) <- rownames(mat0)
  base_sum <- colSums(mat0)

  for (b in seq_len(n_perm)) {
    order_k <- sample.int(k, k, replace = FALSE)
    current_sum <- base_sum
    current_e <- equiv_age_from_curve(setNames(current_sum, age_groups), reference_rate)

    for (idx in order_k) {
      current_sum2 <- current_sum + (mat1[idx, ] - mat0[idx, ])
      new_e <- equiv_age_from_curve(setNames(current_sum2, age_groups), reference_rate)
      if (!is.na(current_e) && !is.na(new_e)) {
        contribution[rownames(mat0)[idx]] <- contribution[rownames(mat0)[idx]] + (new_e - current_e)
      }
      current_sum <- current_sum2
      current_e <- new_e
    }
  }
  contribution / n_perm
}

make_matrix <- function(dt, categories, year_value, value_col = "val") {
  wide <- dt %>%
    filter(year == year_value) %>%
    pivot_wider(id_cols = category, names_from = age, values_from = all_of(value_col))

  wide2 <- data.frame(category = categories) %>%
    left_join(wide, by = "category") %>%
    arrange(category)

  mat <- as.matrix(wide2[, age_groups, drop = FALSE])
  mat[is.na(mat)] <- 0
  rownames(mat) <- wide2$category
  mat
}

run_shapley <- function(input_file, output_file, sex_value = "Both") {
  raw <- fread(input_file)
  dt <- raw %>%
    filter(
      measure == "DALYs (Disability-Adjusted Life Years)",
      metric == "Rate",
      location == "Global",
      sex == sex_value,
      year %in% c(year0, year1),
      age %in% age_groups
    ) %>%
    select(year, age, category, val, lower, upper)

  stopifnot(all(age_groups %in% unique(dt$age)))
  stopifnot(all(c(year0, year1) %in% unique(dt$year)))

  categories <- sort(unique(dt$category))
  mat0 <- make_matrix(dt, categories, year0)
  mat1 <- make_matrix(dt, categories, year1)
  reference_rate <- mean(colSums(mat1)[reference_ages], na.rm = TRUE)

  point <- shapley_equiv_age(mat0, mat1, reference_rate, n_perm = n_perm, seed = seed)

  set.seed(seed)
  boot <- matrix(NA_real_, nrow = n_boot, ncol = length(categories))
  colnames(boot) <- categories

  for (b in seq_len(n_boot)) {
    dt_boot <- dt %>%
      rowwise() %>%
      mutate(
        sd_est = (upper - lower) / (2 * 1.96),
        val_boot = rnorm(1, mean = val, sd = sd_est)
      ) %>%
      ungroup() %>%
      select(year, age, category, val = val_boot)

    mat0_boot <- make_matrix(dt_boot, categories, year0)
    mat1_boot <- make_matrix(dt_boot, categories, year1)
    boot[b, ] <- shapley_equiv_age(mat0_boot, mat1_boot, reference_rate, n_perm = n_perm)
  }

  out <- data.frame(
    category = names(point),
    shapley_contribution = as.numeric(point),
    lower_95ci = apply(boot, 2, quantile, probs = 0.025, na.rm = TRUE),
    upper_95ci = apply(boot, 2, quantile, probs = 0.975, na.rm = TRUE),
    se = apply(boot, 2, sd, na.rm = TRUE)
  ) %>%
    mutate(
      ci_width = upper_95ci - lower_95ci,
      significant = (lower_95ci > 0 & upper_95ci > 0) | (lower_95ci < 0 & upper_95ci < 0)
    ) %>%
    arrange(desc(abs(shapley_contribution)))

  fwrite(out, output_file)
}

if (demo_mode) {
  run_shapley(
    file.path(input_dir, "system_age_related_daly_by_age.csv"),
    file.path(output_dir, "Global_EA_Shapley_demo_with_CI.csv"),
    sex_value = "Both"
  )
} else {
  input_file <- file.path(input_dir, "system_age_related_daly_by_age.csv")
  available_sexes <- fread(input_file, select = "sex") %>%
    distinct(sex) %>%
    pull(sex) %>%
    sort()

  for (sex_name in available_sexes) {
    run_shapley(
      input_file,
      file.path(output_dir, paste0(sex_name, "_Global_EA_Shapley_1990_2021_with_CI.csv")),
      sex_value = sex_name
    )
  }
}
