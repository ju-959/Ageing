library(dplyr)
library(readr)

root_dir <- normalizePath(Sys.getenv("AGEING_ROOT", unset = getwd()), mustWork = TRUE)
demo_mode <- identical(Sys.getenv("AGEING_DEMO", unset = "FALSE"), "TRUE")

output_dir <- if (demo_mode) file.path(root_dir, "outputs", "demo") else file.path(root_dir, "outputs", "analysis")
aapc_input_dir <- file.path(output_dir, "aapc_input")
aapc_output_dir <- file.path(output_dir, "aapc")
aapc_work_dir <- Sys.getenv("AGEING_AAPC_WORKDIR", unset = file.path(tempdir(), "ageing_aapc_joinpoint"))

dir.create(aapc_input_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(aapc_output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(aapc_work_dir, showWarnings = FALSE, recursive = TRUE)

historical_equivalent_file <- file.path(output_dir, "equivalent_age_1990_2021.csv")
projected_equivalent_file <- file.path(output_dir, "projected_equivalent_age_2022_2040.csv")

required_files <- c(historical_equivalent_file, projected_equivalent_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required equivalent-age file(s):\n", paste(missing_files, collapse = "\n"), call. = FALSE)
}

resolve_joinpoint_cmd <- function(root_dir) {
  candidates <- c(
    Sys.getenv("JOINPOINT_CMD", unset = ""),
    file.path(root_dir, "tools", "joinpoint", "jpCommand.exe"),
    "C:/Program Files (x86)/Joinpoint Command/jpCommand.exe",
    "C:/Program Files/Joinpoint Command/jpCommand.exe"
  )
  candidates <- candidates[nzchar(candidates)]
  found <- candidates[file.exists(candidates)]
  if (length(found) > 0) normalizePath(found[1], mustWork = TRUE) else NA_character_
}

write_joinpoint_run_options <- function(path, model = "ln", max_joinpoints = 5, n_cores = 2) {
  lines <- c(
    "[Session Options]",
    paste0("Model=", model),
    "Data shift=0",
    "Minimum joinpoints=0",
    paste0("Maximum joinpoints=", max_joinpoints),
    "Method=grid",
    "Min obs end=2",
    "Min obs between=2",
    "Num obs between=0",
    "Model selection method=permutation test",
    "Permutations significance level=0.05",
    "Num permutations=4499",
    "Early stopping=b-value",
    "Run type=provided",
    "Rates per N=100000",
    "Dependent variable type=age-adjusted rate",
    "Het error=constant variance",
    "CI method=parametric",
    paste0("Num cores=", n_cores),
    "Delay type=delay",
    "Autocorr errors=number",
    "Jump model=false",
    "Comparability ratio=false",
    "Include standard analysis=false",
    "Jump location=9999",
    "Comparability ratio value=0",
    "CR variance=0",
    "Joinpoint alpha level=0.05",
    "APC alpha level=0.05",
    "AAPC alpha level=0.05",
    "Jump CR alpha level=0.05",
    "Random number generator seed=7160",
    "empirical quantile seed=10000",
    "empirical quantile seed type=constant",
    "number of resamples=1000",
    "madwd=false",
    "madwd psi value=0"
  )
  writeLines(lines, path, useBytes = TRUE)
}

write_joinpoint_export_options <- function(path) {
  lines <- c(
    "[Export Options]",
    "Models=best fit",
    "Line delimiter=unix",
    "Missing character=period",
    "Field delimiter=comma",
    "By-var format=quoted labels",
    "Output by-group headers=false",
    "Remove JP flags=true",
    "Remove best fit flags=false",
    "All models in same column=false",
    "Include JP estimates=false",
    "Include apcs in data export=true",
    "X Values Precision=9",
    "Y Values Precision=9",
    "Model Statistics Precision=3",
    "Estimated Joinpoints Precision=3",
    "Regression Coefficients Precision=3",
    "Covariance Matrix Precision=3",
    "Correlation Matrix Precision=3",
    "APC Precision=3",
    "AAPC Precision=3",
    "AAPC Segment Ranges Precision=3",
    "P-Value Precision=3",
    "AAPC Full Range=true",
    "AAPC Last Obs=false",
    "Export Bad Cohorts=true",
    "Export Report=true",
    "Export data=true",
    "Export apc=true",
    "Export aapc=true",
    "Export ftest=true",
    "Export pairwise=true",
    "Export jump_cr=true"
  )
  writeLines(lines, path, useBytes = TRUE)
}

write_joinpoint_session <- function(path, dataset_file) {
  lines <- c(
    "[Datafile options]",
    paste0("Datafile name=", dataset_file),
    "File format=DOS/Windows",
    "Field delimiter=tab",
    "Missing character=period",
    "Fields with delimiter in quotes=false",
    "Variable names include=false",
    "",
    "[Joinpoint Session Parameters]",
    "age-adjusted rate=val",
    "age-adjusted rate location=6",
    "independent variable=year",
    "independent variable location=5",
    "by-var1=measure",
    "by-var1 location=1",
    "by-var2=location",
    "by-var2 location=2",
    "by-var3=sex",
    "by-var3 location=3",
    "by-var4=cause",
    "by-var4 location=4"
  )
  writeLines(lines, path, useBytes = TRUE)
}

write_joinpoint_launcher <- function(path, session_file, export_file, run_file, output_file) {
  lines <- c(
    "[Joinpoint Input Files]",
    paste0("Session File=", session_file),
    paste0("Output File=", output_file),
    paste0("Export Options File=", export_file),
    paste0("Run Options File=", run_file)
  )
  writeLines(lines, path, useBytes = TRUE)
}

clean_label <- function(x) {
  gsub(",", "", as.character(x), fixed = TRUE)
}

format_ci <- function(point, lower, upper, digits = 2, sep = ",") {
  paste0(
    round(as.numeric(point), digits),
    "\n(",
    round(as.numeric(lower), digits),
    sep,
    round(as.numeric(upper), digits),
    ")"
  )
}

prepare_aapc_input <- function(data, cause_column = NULL, default_cause = "Age-related disease") {
  if (is.null(cause_column)) {
    data$cause <- default_cause
  } else {
    data$cause <- data[[cause_column]]
  }

  out <- data %>%
    mutate(
      location = clean_label(location),
      sex = as.character(sex),
      year = as.integer(year),
      cause = clean_label(if_else(cause == "Overall", default_cause, as.character(cause))),
      age = "Age-standardized",
      metric = "Rate",
      measure = "DALYs (Disability-Adjusted Life Years)",
      val = as.numeric(equivalent_age_numeric),
      lower = as.numeric(lower_bound_numeric),
      upper = as.numeric(upper_bound_numeric)
    ) %>%
    filter(!is.na(val), val > 0) %>%
    select(any_of(c("aggregation_level", "aggregate_location")),
           location, sex, year, cause, age, metric, measure, val, lower, upper)

  out %>% arrange(cause, location, sex, year)
}

read_joinpoint_table <- function(path, columns) {
  if (!file.exists(path)) {
    stop("Missing Joinpoint output: ", path, call. = FALSE)
  }

  out <- read.table(
    path,
    row.names = NULL,
    header = TRUE,
    sep = ",",
    quote = "\"",
    fill = TRUE,
    check.names = FALSE
  )

  if (ncol(out) >= length(columns)) {
    names(out)[seq_along(columns)] <- columns
  }
  out
}

postprocess_aapc <- function(aapc, digits = 2, sep = ",") {
  if ("Significant_indicator" %in% names(aapc)) {
    aapc$Significant_indicator[aapc$Significant_indicator == 0] <- "No"
    aapc$Significant_indicator[aapc$Significant_indicator == 1] <- "Yes"
  }

  numeric_columns <- intersect(c("AAPC", "AAPC_LCI", "AAPC_UCI", "P.Value"), names(aapc))
  aapc[numeric_columns] <- lapply(aapc[numeric_columns], function(x) suppressWarnings(as.numeric(x)))

  if (all(c("AAPC", "AAPC_LCI", "AAPC_UCI") %in% names(aapc))) {
    aapc$AAPC_95CI <- format_ci(aapc$AAPC, aapc$AAPC_LCI, aapc$AAPC_UCI, digits = digits, sep = sep)
  }

  if ("AAPC.Index" %in% names(aapc)) {
    aapc <- aapc[aapc$AAPC.Index == "Full Range", , drop = FALSE]
  }

  aapc$metric <- "Rate"
  aapc$age <- "Age-standardized"
  aapc
}

run_joinpoint_aapc <- function(data, analysis_name, joinpoint_cmd, joinpoints = 5, digits = 2, sep = ",") {
  work_dir <- file.path(aapc_work_dir, analysis_name)
  dir.create(work_dir, showWarnings = FALSE, recursive = TRUE)

  joinpoint_input <- data %>%
    select(measure, location, sex, cause, year, val) %>%
    arrange(measure, location, sex, cause, year)

  write.table(
    joinpoint_input,
    file.path(work_dir, "joinpoint.csv"),
    row.names = FALSE,
    col.names = FALSE,
    sep = "\t",
    quote = TRUE,
    na = "."
  )

  detected_cores <- parallel::detectCores(logical = TRUE)
  if (is.na(detected_cores)) detected_cores <- 2L
  n_cores <- max(2L, detected_cores - 2L)

  session_file <- file.path(work_dir, "AAPC.Created.Session.ini")
  export_file <- file.path(work_dir, "AAPC.Export.Options.ini")
  run_file <- file.path(work_dir, "AAPC.Run.Options.ini")
  launcher_file <- file.path(work_dir, "AAPC.JPRun.ini")
  output_file <- file.path(work_dir, "AAPC.jpo")

  write_joinpoint_session(session_file, file.path(work_dir, "joinpoint.csv"))
  write_joinpoint_export_options(export_file)
  write_joinpoint_run_options(run_file, model = "ln", max_joinpoints = joinpoints, n_cores = n_cores)
  write_joinpoint_launcher(launcher_file, session_file, export_file, run_file, output_file)

  old_wd <- getwd()
  setwd(work_dir)
  on.exit(setwd(old_wd), add = TRUE)

  if (.Platform$OS.type == "windows") {
    status <- system2(joinpoint_cmd, shQuote(launcher_file), stdout = TRUE, stderr = TRUE)
  } else {
    wine_cmd <- Sys.getenv("WINE", unset = "")
    if (nzchar(wine_cmd)) {
      status <- system2(wine_cmd, c(shQuote(joinpoint_cmd), shQuote(launcher_file)), stdout = TRUE, stderr = TRUE)
    } else {
      status <- system2(joinpoint_cmd, shQuote(launcher_file), stdout = TRUE, stderr = TRUE)
    }
  }

  exit_status <- attr(status, "status")
  if (!is.null(exit_status) && !identical(exit_status, 0L)) {
    stop("Joinpoint failed for ", analysis_name, ".", call. = FALSE)
  }

  error_file <- file.path(work_dir, "AAPC.ErrorFile.txt")
  if (file.exists(error_file)) {
    stop("Joinpoint error for ", analysis_name, ": ", error_file, call. = FALSE)
  }

  aapc_path <- file.path(work_dir, "AAPC.aapcexport.txt")
  if (!file.exists(aapc_path)) {
    candidates <- list.files(work_dir, pattern = "[.]aapcexport[.]txt$", full.names = TRUE)
    if (length(candidates) > 0) aapc_path <- candidates[1]
  }

  aapc <- read_joinpoint_table(
    aapc_path,
    c(
      "measure", "location", "sex", "cause", "joinpoint", "AAPC.Index",
      "Start.Obs", "End.Obs", "AAPC", "AAPC_LCI", "AAPC_UCI",
      "Significant_indicator", "Test.Statistic", "P.Value"
    )
  )

  postprocess_aapc(aapc, digits = digits, sep = sep)
}

run_split_by_cause <- function(data, analysis_name, joinpoint_cmd) {
  causes <- sort(unique(data$cause))
  pieces <- vector("list", length(causes))
  names(pieces) <- causes

  for (cause_name in causes) {
    message("Running AAPC: ", analysis_name, " | ", cause_name)
    pieces[[cause_name]] <- run_joinpoint_aapc(
      data = data %>% filter(cause == cause_name),
      analysis_name = paste0(analysis_name, "_", gsub("[^A-Za-z0-9]+", "_", cause_name)),
      joinpoint_cmd = joinpoint_cmd
    )
  }

  bind_rows(pieces)
}

historical <- read_csv(historical_equivalent_file, show_col_types = FALSE)
projected <- read_csv(projected_equivalent_file, show_col_types = FALSE)

historical_overall <- historical %>%
  filter(analysis_level == "overall") %>%
  prepare_aapc_input(cause_column = NULL, default_cause = "Age-related disease")

historical_system <- historical %>%
  filter(analysis_level == "system") %>%
  prepare_aapc_input(cause_column = "category")

projected_prepared <- projected %>%
  mutate(
    cause = if_else(cause == "Age-related diseases", "Age-related disease", cause)
  ) %>%
  prepare_aapc_input(cause_column = "cause")

projected_overall <- projected_prepared %>% filter(cause == "Age-related disease")
projected_system <- projected_prepared %>% filter(cause != "Age-related disease")

write_csv(historical_overall, file.path(aapc_input_dir, "AAPC_input_Overall_Historical_1990_2021.csv"))
write_csv(historical_system, file.path(aapc_input_dir, "AAPC_input_DiseaseSystems_Historical_1990_2021.csv"))
write_csv(projected_overall, file.path(aapc_input_dir, "AAPC_input_Overall_Projection_2022_2040.csv"))
write_csv(projected_system, file.path(aapc_input_dir, "AAPC_input_DiseaseSystems_Projection_2022_2040.csv"))

joinpoint_cmd <- resolve_joinpoint_cmd(root_dir)
if (is.na(joinpoint_cmd)) {
  stop("Set JOINPOINT_CMD to the full path of jpCommand.exe.", call. = FALSE)
}

historical_overall_aapc <- run_joinpoint_aapc(
  historical_overall,
  analysis_name = "overall_historical",
  joinpoint_cmd = joinpoint_cmd
)

projected_overall_aapc <- run_joinpoint_aapc(
  projected_overall,
  analysis_name = "overall_projection",
  joinpoint_cmd = joinpoint_cmd
)

historical_system_aapc <- run_split_by_cause(
  historical_system,
  analysis_name = "system_historical",
  joinpoint_cmd = joinpoint_cmd
)

projected_system_aapc <- run_split_by_cause(
  projected_system,
  analysis_name = "system_projection",
  joinpoint_cmd = joinpoint_cmd
)

write_csv(historical_overall_aapc, file.path(aapc_output_dir, "AAPC_Overall_Historical_1990_2021.csv"))
write_csv(projected_overall_aapc, file.path(aapc_output_dir, "AAPC_Overall_Projection_2022_2040.csv"))
write_csv(historical_system_aapc, file.path(aapc_output_dir, "AAPC_DiseaseSystems_Historical_1990_2021.csv"))
write_csv(projected_system_aapc, file.path(aapc_output_dir, "AAPC_DiseaseSystems_Projection_2022_2040.csv"))
