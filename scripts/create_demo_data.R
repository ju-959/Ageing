root_dir <- normalizePath(Sys.getenv("AGEING_ROOT", unset = getwd()), mustWork = TRUE)
demo_dir <- file.path(root_dir, "data", "demo")
dir.create(demo_dir, showWarnings = FALSE, recursive = TRUE)

age_definition <- c(
  "25 to 29", "30 to 34", "35 to 39", "40 to 44", "45 to 49",
  "50 to 54", "55 to 59", "60 to 64", "65 to 69", "70 to 74",
  "75 to 79", "80 to 84"
)
age_projection <- c(age_definition, "85 to 89", "90 to 94")
age_all <- c("0 to 4", "5 to 9", "10 to 14", "15 to 19", "20 to 24", age_projection, "95 plus")
mid <- c(
  "0 to 4" = 2, "5 to 9" = 7, "10 to 14" = 12, "15 to 19" = 17, "20 to 24" = 22,
  "25 to 29" = 27, "30 to 34" = 32, "35 to 39" = 37, "40 to 44" = 42,
  "45 to 49" = 47, "50 to 54" = 52, "55 to 59" = 57, "60 to 64" = 62,
  "65 to 69" = 67, "70 to 74" = 72, "75 to 79" = 77, "80 to 84" = 82,
  "85 to 89" = 87, "90 to 94" = 92, "95 plus" = 97
)

make_rate_rows <- function(cause, measure, ages, fn) {
  rows <- lapply(ages, function(age) {
    val <- fn(mid[[age]])
    data.frame(
      measure = measure, location = "Global", sex = "Both", age = age, year = 2021,
      cause = cause, metric = "Rate", val = val, lower = val * 0.98, upper = val * 1.02
    )
  })
  do.call(rbind, rows)
}

demo_causes <- data.frame(
  cause = c(
    "Atrial fibrillation and flutter",
    "Ischemic heart disease",
    "Hypertensive heart disease",
    "Chronic kidney disease due to diabetes mellitus type 2",
    "Chronic kidney disease due to hypertension",
    "Chronic kidney disease due to other and unspecified causes",
    "Aortic aneurysm",
    "Diabetes mellitus type 2"
  ),
  category = c(
    "Cardiovascular diseases",
    "Cardiovascular diseases",
    "Cardiovascular diseases",
    "Diabetes and kidney diseases",
    "Diabetes and kidney diseases",
    "Diabetes and kidney diseases",
    "Cardiovascular diseases",
    "Diabetes and kidney diseases"
  ),
  screen_source = c(
    "Incidence", "Incidence", "Prevalence",
    "Incidence", "Incidence", "Incidence",
    "Prevalence", "Incidence"
  ),
  age_related = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE),
  daly_offset = c(3900, 5400, 3600, 3100, 2700, 2300, 1200, 1000),
  daly_slope = c(64, 82, 60, 54, 50, 46, 12, 10),
  stringsAsFactors = FALSE
)

screen_curve <- function(is_age_related) {
  if (is_age_related) {
    function(x) 12 + 0.35 * (x - 25) + 0.050 * (x - 25)^2
  } else {
    function(x) 65 - 0.20 * (x - 25)
  }
}

incidence <- do.call(
  rbind,
  lapply(which(demo_causes$screen_source == "Incidence"), function(i) {
    make_rate_rows(
      demo_causes$cause[i],
      "Incidence",
      age_definition,
      screen_curve(demo_causes$age_related[i])
    )
  })
)

prevalence <- do.call(
  rbind,
  lapply(which(demo_causes$screen_source == "Prevalence"), function(i) {
    make_rate_rows(
      demo_causes$cause[i],
      "Prevalence",
      age_definition,
      screen_curve(demo_causes$age_related[i])
    )
  })
)

write.csv(incidence, file.path(demo_dir, "incidence_age_pattern.csv"), row.names = FALSE)
write.csv(prevalence, file.path(demo_dir, "prevalence_age_pattern.csv"), row.names = FALSE)

map <- demo_causes[, c("cause", "category")]
write.csv(map, file.path(demo_dir, "disease_system_map.csv"), row.names = FALSE)

locations <- c("Global", "Japan", "Zimbabwe")
years <- 1990:2021
causes <- demo_causes$cause
projection_locations <- locations

population_history_file <- file.path(demo_dir, "population_history_2021.csv")
population_projection_file <- file.path(demo_dir, "population_projection_who.csv")
full_population_history_file <- file.path(root_dir, "data", "processed", "population_history_2021.csv.gz")
full_population_projection_file <- file.path(root_dir, "data", "processed", "population_projection_who.csv.gz")

if (file.exists(full_population_history_file) && file.exists(full_population_projection_file)) {
  full_population_history <- read.csv(full_population_history_file, stringsAsFactors = FALSE)
  full_population_projection <- read.csv(full_population_projection_file, stringsAsFactors = FALSE)

  population_history_excerpt <- full_population_history[
    full_population_history$location %in% projection_locations &
      full_population_history$sex %in% c("Male", "Female") &
      full_population_history$year %in% 1990:2021 &
      full_population_history$age %in% age_all,
  ]

  population_projection_excerpt <- full_population_projection[
    full_population_projection$location %in% projection_locations &
      full_population_projection$sex %in% c("Male", "Female") &
      full_population_projection$year %in% 2022:2041 &
      full_population_projection$age %in% age_all,
  ]

  write.csv(population_history_excerpt, population_history_file, row.names = FALSE)
  write.csv(population_projection_excerpt, population_projection_file, row.names = FALSE)
}

if (!file.exists(population_history_file) || !file.exists(population_projection_file)) {
  stop("Demo population files are missing from data/demo.", call. = FALSE)
}

population_history <- read.csv(population_history_file, stringsAsFactors = FALSE)
population_projection <- read.csv(population_projection_file, stringsAsFactors = FALSE)

required_population_columns <- c("location", "sex", "year", "age", "val")
if (!all(required_population_columns %in% names(population_history)) ||
    !all(required_population_columns %in% names(population_projection))) {
  stop("Demo population files must contain location, sex, year, age, and val columns.", call. = FALSE)
}

demo_population <- function(location, sex, year, age) {
  if (sex == "Both") {
    male <- demo_population(location, "Male", year, age)
    female <- demo_population(location, "Female", year, age)
    return(male + female)
  }

  row <- population_history[
    population_history$location == location &
      population_history$sex == sex &
      population_history$year == year &
      population_history$age == age,
  ]

  if (nrow(row) != 1) {
    stop(
      "Population lookup failed for ",
      paste(location, sex, year, age, sep = " / "),
      call. = FALSE
    )
  }

  row$val
}

validate_projection_population <- function(population_projection) {
  expected <- expand.grid(
    location = projection_locations,
    sex = c("Male", "Female"),
    year = 2022:2041,
    age = age_all,
    stringsAsFactors = FALSE
  )

  observed <- population_projection[, required_population_columns]
  missing_rows <- dplyr::anti_join(expected, observed, by = required_population_columns[-5])

  if (nrow(missing_rows) > 0) {
    stop("Demo projected population file is incomplete.", call. = FALSE)
  }

  invisible(TRUE)
}

validate_projection_population(population_projection)

demo_location_offset <- function(location) {
  dplyr::case_when(
    location == "Japan" ~ -360,
    location == "Zimbabwe" ~ 420,
    TRUE ~ 0
  )
}

demo_sex_offset <- function(sex) {
  if (sex == "Female") {
    -220
  } else if (sex == "Male") {
    220
  } else {
    0
  }
}

demo_rate <- function(cause, location, sex, year, age) {
  cause_info <- demo_causes[demo_causes$cause == cause, ]
  trend <- if (isTRUE(cause_info$age_related)) -10 * (year - 2021) else -3 * (year - 2021)
  pmax(
    25,
    cause_info$daly_offset +
      cause_info$daly_slope * (mid[[age]] - 65) +
      trend + demo_location_offset(location) + demo_sex_offset(sex)
  )
}

metric_rows <- list()
k <- 1
for (cause in causes) {
  for (location in locations) {
    for (sex in c("Male", "Female", "Both")) {
      for (year in years) {
        for (age in age_projection) {
          value_rate <- demo_rate(cause, location, sex, year, age)
          value_number <- value_rate * demo_population(location, sex, year, age) / 100000
          metric_rows[[k]] <- data.frame(
            measure = "DALYs (Disability-Adjusted Life Years)", location = location, sex = sex,
            age = age, year = year, cause = cause, metric = "Rate",
            val = value_rate, lower = value_rate * 0.95, upper = value_rate * 1.05
          )
          k <- k + 1
          metric_rows[[k]] <- data.frame(
            measure = "DALYs (Disability-Adjusted Life Years)", location = location, sex = sex,
            age = age, year = year, cause = cause, metric = "Number",
            val = value_number, lower = value_number * 0.95, upper = value_number * 1.05
          )
          k <- k + 1
        }
      }
    }
  }
}
daly <- do.call(rbind, metric_rows)
write.csv(daly, file.path(demo_dir, "daly_age_specific.csv"), row.names = FALSE)

write.csv(
  data.frame(location = c("Japan", "Zimbabwe"), aggregate_location = "Demo region"),
  file.path(demo_dir, "projection_location_groups.csv"),
  row.names = FALSE
)

message("Demo input files written to: ", demo_dir)
