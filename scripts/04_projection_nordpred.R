library(dplyr)
library(readr)
library(tidyr)

# Nordpred helper ----------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(reshape2)
})

normalise_age_labels <- function(data) {
  data$age <- sub("95\\+", "95 plus", data$age)
  data$age <- sub("95 to 99", "95 plus", data$age)
  data
}

prepare_standard_population <- function(standard_population, gbd_edition) {
  std <- standard_population
  if (gbd_edition %in% c(2021, 2023)) {
    std$std_population[6] <- sum(std$std_population[1:6])
    std[6, 1] <- "0 to 4"
    std <- std[-c(1:5), ]
  } else if (gbd_edition == 2019) {
    std <- std[-c(1:3), ]
    std$std_population[2] <- std$std_population[1] + std$std_population[2]
    std <- std[-1, ]
    std[1, 1] <- "0 to 4"
  }
  std$percent <- std$std_population / sum(std$std_population)
  std
}

datagg <- function(cases, pyr, nagg = 5) {
  aggregate_columns <- function(x, fun) {
    x <- as.matrix(x)
    groups <- ceiling(seq_len(ncol(x)) / nagg)
    out <- sapply(sort(unique(groups)), function(g) {
      apply(x[, groups == g, drop = FALSE], 1, fun, na.rm = TRUE)
    })
    out <- as.data.frame(out)
    rownames(out) <- rownames(x)
    colnames(out) <- tapply(colnames(x), groups, function(y) paste0(y[1], "-", y[length(y)]))
    out
  }

  list(
    cases = aggregate_columns(cases, sum),
    pyr = aggregate_columns(pyr, sum)
  )
}

asrpy <- function(age_specific_rate, data, population, startp, nagg = 5) {
  age_specific_rate <- as.data.frame(age_specific_rate)
  population <- as.data.frame(population)

  if (ncol(age_specific_rate) == ncol(population)) {
    colnames(age_specific_rate) <- colnames(population)
    return(age_specific_rate)
  }

  period_index <- pmin(ceiling(seq_len(ncol(population)) / nagg), ncol(age_specific_rate))
  out <- age_specific_rate[, period_index, drop = FALSE]
  colnames(out) <- colnames(population)
  rownames(out) <- rownames(age_specific_rate)
  out
}

asry <- function(age_specific_rate, population, standpop) {
  rate_matrix <- as.matrix(age_specific_rate)
  population_matrix <- as.matrix(population)
  standpop <- as.numeric(standpop)
  standpop <- standpop / sum(standpop, na.rm = TRUE)

  cases <- rate_matrix * population_matrix / 100000
  data.frame(
    ASR = colSums(rate_matrix * standpop, na.rm = TRUE),
    crude_rate = colSums(cases, na.rm = TRUE) / colSums(population_matrix, na.rm = TRUE) * 100000,
    case = colSums(cases, na.rm = TRUE),
    row.names = colnames(rate_matrix)
  )
}

nordpred.estimate <- function(cases, pyr, noperiod, startestage, linkfunc = "power5") {
  if (nrow(cases) != nrow(pyr)) {
    stop("\"cases\" and \"pyr\" must have the same age groups.", call. = FALSE)
  }
  if (ncol(cases) > ncol(pyr)) {
    stop("\"pyr\" must include all observed periods and future periods.", call. = FALSE)
  }
  if (ncol(pyr) == ncol(cases)) {
    stop("\"pyr\" must include future periods.", call. = FALSE)
  }
  if ((ncol(cases) - noperiod) < 0 || noperiod < 3) {
    stop("\"noperiod\" must be at least 3 and no larger than the observed periods.", call. = FALSE)
  }

  dnoperiods <- ncol(cases)
  dnoagegr <- nrow(cases)
  ageno <- rep(seq_len(dnoagegr), dnoperiods)
  periodno <- sort(rep(seq_len(dnoperiods), dnoagegr))
  cohort <- max(ageno) - ageno + periodno
  y <- c(as.matrix(pyr[, seq_len(dnoperiods), drop = FALSE]))

  apcdata <- data.frame(
    Cases = c(as.matrix(cases)),
    Age = ageno,
    Cohort = cohort,
    Period = periodno,
    y = y
  )
  apcdata <- apcdata[apcdata$Age >= startestage, ]
  apcdata <- apcdata[apcdata$Period > (dnoperiods - noperiod), ]

  if (linkfunc == "power5") {
    y <- apcdata$y
    power5link <- poisson()
    power5link$link <- "0.2 root link Poisson family"
    power5link$linkfun <- function(mu) (mu / y)^0.2
    power5link$linkinv <- function(eta) pmax(.Machine$double.eps, y * eta^5)
    power5link$mu.eta <- function(eta) pmax(.Machine$double.eps, 5 * y * eta^4)
    res.glm <- suppressWarnings(glm(
      Cases ~ as.factor(Age) + Period + as.factor(Period) + as.factor(Cohort) - 1,
      family = power5link,
      data = apcdata
    ))
  } else if (linkfunc == "poisson") {
    res.glm <- suppressWarnings(glm(
      Cases ~ as.factor(Age) + Period + as.factor(Period) + as.factor(Cohort) +
        offset(log(y)) - 1,
      family = poisson(),
      data = apcdata
    ))
  } else {
    stop("Unknown link function.", call. = FALSE)
  }

  pvalue <- 1 - pchisq(res.glm$deviance, res.glm$df.residual)
  mod1 <- suppressWarnings(glm(
    Cases ~ as.factor(Age) + Period + as.factor(Cohort) + offset(log(y)) - 1,
    family = poisson(),
    data = apcdata
  ))
  mod2 <- suppressWarnings(glm(
    Cases ~ as.factor(Age) + Period + I(Period^2) + as.factor(Cohort) + offset(log(y)) - 1,
    family = poisson(),
    data = apcdata
  ))
  pdiff <- suppressWarnings(anova(mod1, mod2, test = "Chisq")$"P(>|Chi|)"[2])
  if (is.null(pdiff)) pdiff <- suppressWarnings(anova(mod1, mod2, test = "Chisq")$"Pr(>Chi)"[2])

  res <- list(
    glm = res.glm,
    cases = cases,
    pyr = pyr,
    noperiod = noperiod,
    gofpvalue = pvalue,
    startestage = startestage,
    suggestionrecent = isTRUE(pdiff < 0.05),
    pvaluerecent = pdiff,
    linkfunc = linkfunc
  )
  class(res) <- "nordpred.estimate"
  res
}

nordpred.prediction <- function(nordpred.estimate.object, startuseage, recent,
                                cuttrend = c(0, 0.25, 0.5, 0.75, 0.75)) {
  cases <- nordpred.estimate.object$cases
  pyr <- nordpred.estimate.object$pyr
  noperiod <- nordpred.estimate.object$noperiod
  nototper <- ncol(pyr)
  noobsper <- ncol(cases)
  nonewpred <- nototper - noobsper
  dnoagegr <- nrow(cases)

  if (length(cuttrend) < nonewpred) {
    cuttrend <- c(cuttrend, rep(tail(cuttrend, 1), nonewpred - length(cuttrend)))
  } else {
    cuttrend <- cuttrend[seq_len(nonewpred)]
  }

  datatable <- matrix(NA_real_, dnoagegr, nototper)
  datatable[, seq_len(noobsper)] <- as.matrix(cases)
  datatable <- as.data.frame(datatable)
  row.names(datatable) <- rownames(cases)
  names(datatable) <- colnames(pyr)

  for (age in seq_len(startuseage - 1)) {
    obsinc <- cases[age, (noobsper - 1):noobsper] / pyr[age, (noobsper - 1):noobsper]
    obsinc[is.na(obsinc)] <- 0
    datatable[age, (noobsper + 1):nototper] <- mean(as.numeric(obsinc)) * pyr[age, (noobsper + 1):nototper]
  }

  for (age in startuseage:dnoagegr) {
    startestage <- nordpred.estimate.object$startestage
    coefficients <- nordpred.estimate.object$glm$coefficients
    coh <- (dnoagegr - startestage) - (age - startestage) + (noperiod + seq_len(nonewpred))
    noages <- dnoagegr - startestage + 1
    driftmp <- cumsum(1 - cuttrend)
    cohfind <- noages + (noperiod - 1) + 1 + (coh - 1)
    maxcoh <- dnoagegr - startuseage + noperiod
    agepar <- as.numeric(coefficients[age - startestage + 1])
    driftfind <- pmatch("Period", names(coefficients))
    driftpar <- as.numeric(coefficients[driftfind])
    cohpar <- rep(NA_real_, length(coh))

    for (i in seq_along(coh)) {
      if (coh[i] < maxcoh) {
        cohpar[i] <- as.numeric(coefficients[cohfind[i]])
      } else {
        cohpar[i] <- as.numeric(coefficients[length(coefficients) - (startuseage - startestage)])
        if (is.na(cohpar[i])) cohpar[i] <- 0
      }
    }

    if (recent) {
      lpfind <- driftfind + noperiod - 2
      driftrecent <- driftpar - as.numeric(coefficients[lpfind])
    }

    if (nordpred.estimate.object$linkfunc == "power5") {
      if (recent) {
        rate <- (agepar + driftpar * noobsper + driftrecent * driftmp + cohpar)^5
      } else {
        rate <- (agepar + driftpar * (noobsper + driftmp) + cohpar)^5
      }
    } else {
      if (recent) {
        rate <- exp(agepar + driftpar * noobsper + driftrecent * driftmp + cohpar)
      } else {
        rate <- exp(agepar + driftpar * (noobsper + driftmp) + cohpar)
      }
    }

    datatable[age, (noobsper + 1):nototper] <- rate * pyr[age, (noobsper + 1):nototper]
  }

  res <- list(
    predictions = datatable,
    pyr = pyr,
    nopred = nonewpred,
    noperiod = noperiod,
    gofpvalue = nordpred.estimate.object$gofpvalue,
    recent = recent,
    pvaluerecent = nordpred.estimate.object$pvaluerecent,
    cuttrend = cuttrend,
    startuseage = startuseage,
    startestage = startestage,
    glm = nordpred.estimate.object$glm
  )
  class(res) <- "nordpred"
  res
}

nordpred <- function(cases, pyr, startestage, startuseage, noperiods = NULL,
                     recent = NULL, cuttrend = c(0, 0.25, 0.5, 0.75, 0.75),
                     linkfunc = "power5") {
  percases <- ncol(cases)
  if (percases < 3) stop("Too few periods in cases.", call. = FALSE)
  if (is.null(noperiods)) noperiods <- min(percases, 6)
  noperiod <- max(sort(noperiods))

  if (is.null(recent)) {
    recent <- nordpred.estimate(cases, pyr, noperiod, startestage, linkfunc)$suggestionrecent
  }

  est <- nordpred.estimate(cases, pyr, noperiod, startestage, linkfunc)
  nordpred.prediction(est, startuseage = startuseage, recent = recent, cuttrend = cuttrend)
}

nordpred.getpred <- function(nordpred.object, incidence = TRUE) {
  datatable <- nordpred.object$predictions
  pyr <- as.data.frame(nordpred.object$pyr)
  if (incidence) {
    res <- (datatable / pyr) * 100000
    res[is.na(res)] <- 0
    res
  } else {
    datatable
  }
}

gbd_nordpred_prediction <- function(data, measure_name, cause_name,
                                  location_name, rei_name = NULL, By_sex = FALSE,
                                  predyear = 2040, full_age_adjusted = FALSE,
                                  pop_predict = "WHO", gbd_edition = 2021,
                                  population_history, population_projection,
                                  standard_population) {
  if (missing(population_history) || missing(population_projection) || missing(standard_population)) {
    stop("population_history, population_projection, and standard_population must be supplied.", call. = FALSE)
  }
  
  
  
  if ('location_name' %in% names(data) | 'location_id' %in% names(data)){
    stop("Rename *_id and *_name columns before running this script.", call. = FALSE)
  }
  
  
  
  
  if (gbd_edition==2023){
    basicyear = 1994
    endyear =2023
    predyear_selects <- c(2028,2033,2038,2043,2048)
    data <- data %>% filter(year %in%  basicyear:endyear)
    if (predyear>predyear_selects[length(predyear_selects)]){
      stop("predyear exceeds the supported range for GBD 2023.", call. = FALSE)
    }
    predyear1 <- predyear
    if (predyear1 <= predyear_selects[1]){
      predyear <- predyear_selects[1]
    } else if (predyear1 <= predyear_selects[2]){
      predyear <- predyear_selects[2]
    } else if (predyear1 <= predyear_selects[3]){
      predyear <- predyear_selects[3]
    } else if (predyear1 <= predyear_selects[4]){
      predyear <- predyear_selects[4]
    } else if (predyear1 <= predyear_selects[5]){
      predyear <- predyear_selects[5]
    }
    
  } else if (gbd_edition==2021) {
    basicyear = 1992
    endyear =2021
    predyear_selects <- c(2026,2031,2036,2041,2046)
    data <- data %>% filter(year %in%  basicyear:endyear)
    if (predyear>predyear_selects[length(predyear_selects)]){
      stop("predyear exceeds the supported range for GBD 2021.", call. = FALSE)
    }
    predyear1 <- predyear
    if (predyear1 <= predyear_selects[1]){
      predyear <- predyear_selects[1]
    } else if (predyear1 <= predyear_selects[2]){
      predyear <- predyear_selects[2]
    } else if (predyear1 <= predyear_selects[3]){
      predyear <- predyear_selects[3]
    } else if (predyear1 <= predyear_selects[4]){
      predyear <- predyear_selects[4]
    } else if (predyear1 <= predyear_selects[5]){
      predyear <- predyear_selects[5]
    }
    
  } else if (gbd_edition==2019) {
    basicyear = 1990
    endyear = 2019
    predyear_selects <- c(2024,2029,2034,2039,2044)
    if (predyear>predyear_selects[length(predyear_selects)]){
      stop("predyear exceeds the supported range for GBD 2019.", call. = FALSE)
    }
    predyear1 <- predyear
    if (predyear1 <= predyear_selects[1]){
      predyear <- predyear_selects[1]
    } else if (predyear1 <= predyear_selects[2]){
      predyear <- predyear_selects[2]
    } else if (predyear1 <= predyear_selects[3]){
      predyear <- predyear_selects[3]
    } else if (predyear1 <= predyear_selects[4]){
      predyear <- predyear_selects[4]
    } else if (predyear1 <= predyear_selects[5]){
      predyear <- predyear_selects[5]
    }
    
  }
  
  
  
  if (pop_predict=='GBD'){
    GBDpredict_loct <- unique(population_projection$location)
  } else if (pop_predict=='WHO'){
    GBDpredict_loct <- unique(population_projection$location)
  }
  
  data <- data %>% dplyr::filter(location %in% GBDpredict_loct)
  
  location_name <- location_name[location_name %in% unique(data$location)]
  
  if ('location_name' %in% names(data) | 'location_id' %in% names(data)){
    stop("Rename *_id and *_name columns before running this script.", call. = FALSE)
  }
  if (!is.null(rei_name)){
    if (is.logical(rei_name)) {
      stop("rei_name must be a label, not a logical value.", call. = FALSE)
    }
  }
  
  
  if (("rei" %in% names(data))==TRUE & length(rei_name)==0){
    stop("Specify rei_name for inputs with a rei column.", call. = FALSE)
  }
  
  if (!'Number' %in% unique(data$metric)) {
    stop("Input data must include metric == 'Number' for age-specific counts.", call. = FALSE)
  }
  
  if (min(data$val)<0) {
    stop("Input data contain val < 0.", call. = FALSE)
  }
  
  if (nrow(data)==0){
    stop("No usable records remain after filtering to metric == 'Number'.", call. = FALSE)
  }
  
  if (length(unique(data$sex)) !=3){
    stop("Input data must include Male, Female, and Both sex categories.", call. = FALSE)
  }
  
  
  if (!is.null(rei_name)){
    if (length(rei_name %in% unique(data$rei)) != sum(rei_name %in% unique(data$rei))){
      stop("One or more rei_name labels are not present in the input data.", call. = FALSE)
    }
  }
  
  if (length(cause_name %in% unique(data$cause)) != sum(cause_name %in% unique(data$cause))){
    stop("One or more cause_name labels are not present in the input data.", call. = FALSE)
  }
  
  if (length(measure_name %in% unique(data$measure)) != sum(measure_name %in% unique(data$measure))){
    stop("One or more measure_name labels are not present in the input data.", call. = FALSE)
  }
  if (length(location_name %in% unique(data$location)) != sum(location_name %in% unique(data$location))){
    stop("One or more location_name labels are not present in the input data.", call. = FALSE)
  }
  
  
  ages_2 <- c("0 to 4","5 to 9","10 to 14", "15 to 19","20 to 24", "25 to 29",
              "30 to 34", "35 to 39", "40 to 44", "45 to 49", "50 to 54", "55 to 59",
              "60 to 64", "65 to 69", "70 to 74", "75 to 79", "80 to 84",  "85 to 89",
              "90 to 94", "95 plus")
  data$location <- sub("Côte d'Ivoire", replacement = "Coted'Ivoire", data$location)
  data$location <- sub("C.te d'Ivoire", replacement = "Coted'Ivoire", data$location)
  
  location_name <- sub("Côte d'Ivoire", replacement = "Coted'Ivoire", location_name)
  location_name <- sub("C.te d'Ivoire", replacement = "Coted'Ivoire", location_name)
  
  
  
  data <- normalise_age_labels(data)
  
  
  dat <- data
  data <- subset(data, metric=='Number')
  dat2 <- data
  data_NA <- subset(data,val == 0)
  #data <- subset(data,val != 0)
  
  
  
  if (full_age_adjusted==T){
    
  } else {
    data_indicate <- dat2 %>% dplyr::filter(!age %in% c('All ages','Age-standardized'))
    unique_age <- unique(data_indicate$age)
    unique_age <- sub('95 plus',replacement = '95 to 99', unique_age)
    unique_age <- matrix(as.numeric(unlist(strsplit(unique(unique_age),split="to"))),
                         ncol=2,byrow=T)[,1]
    startage <- min(unique_age)
    endage <- max(unique_age)
    
    if (endage==95) {
      endage = '95 plus'
    } else {
      endage <-  paste(endage,endage+4,sep = ' to ')
    }
    startage <-  paste(startage,startage+4,sep = ' to ')
    age_zero <- unique(data_NA$age)
    if (length(age_zero)>0) invisible(age_zero)
    
    
  }
  
  
  
  if (pop_predict=='GBD'){
    GBDpredict <- population_projection
  } else if (pop_predict=='WHO'){
    GBDpredict <- population_projection
  }
  
  location_name2 <- unique(GBDpredict$location)
  location_name <- location_name[which(location_name %in% location_name2)]
  
  
  
  if (gbd_edition==2023){
    population <- population_history %>% dplyr::filter(location %in% location_name) %>%
      dplyr::select(location,sex,year,age,val) %>%
      dplyr::filter(age %in% ages_2) %>% dplyr::filter(year %in% 1994:2023)
    
  } else if (gbd_edition==2021) {
    population <- population_history %>% dplyr::filter(location %in% location_name) %>%
      dplyr::select(location,sex,year,age,val) %>%
      dplyr::filter(age %in% ages_2) %>% dplyr::filter(year %in% 1992:2021)
    
  } else if (gbd_edition==2019) {
    population <- population_history %>% dplyr::filter(location %in% location_name) %>%
      dplyr::filter(age %in% ages_2)
    
  }
  
  GBDpredict <- GBDpredict %>% dplyr::filter(year %in% ((max(population$year)+1):predyear),location %in% location_name,age %in% ages_2)
  
  population <- rbind(population,GBDpredict)
  
  std <- prepare_standard_population(standard_population, gbd_edition)
  
  
  
  a = measure_name[1]
  b = cause_name[1]
  d = location_name[1]
  e = rei_name[1]
  if (is.null(rei_name)==T){
    i=0
    for (a in measure_name){
      for (b in cause_name){
        for (d in location_name){
          if (i==0) {
            result <- .gbd_nordpred_prediction(data,data_NA,NULL,ages_2,dat,population,std,
                                             measure_name=a,cause_name=b,
                                             location_name=d,rei_name,
                                             By_sex,predyear,full_age_adjusted,basicyear,endyear,predyear1)
            i = i + 1
          } else {
            temp <- .gbd_nordpred_prediction(data,data_NA,NULL,ages_2,dat,population,std,
                                           measure_name=a,cause_name=b,
                                           location_name=d,rei_name,
                                           By_sex,predyear,full_age_adjusted,basicyear,endyear,predyear1)
            
            result[["ASR_Number"]] <- rbind(result[["ASR_Number"]],temp[["ASR_Number"]])
            result[["age_specific_rate"]] <- rbind(result[["age_specific_rate"]],temp[["age_specific_rate"]])
            result[["age_specific_projection"]] <- rbind(result[["age_specific_projection"]],temp[["age_specific_projection"]])
            i = i + 1
          }
        }
        
      }
    }
  } else {
    i=0
    for (a in measure_name){
      for (b in cause_name){
        for (d in location_name){
          for (e in rei_name){
            if (i==0) {
              result <- .gbd_nordpred_prediction(data,data_NA,NULL,ages_2,dat,population,std,
                                               measure_name=a,cause_name=b,
                                               location_name=d,rei_name=e,
                                               By_sex,predyear,full_age_adjusted,basicyear,endyear,predyear1)
              i = i + 1
            } else {
              temp <- .gbd_nordpred_prediction(data,data_NA,NULL,ages_2,dat,population,std,
                                             measure_name=a,cause_name=b,
                                             location_name=d,rei_name=e,
                                             By_sex,predyear,full_age_adjusted,basicyear,endyear,predyear1)
              
              
              result[["ASR_Number"]] <- rbind(result[["ASR_Number"]],temp[["ASR_Number"]])
              result[["age_specific_rate"]] <- rbind(result[["age_specific_rate"]],temp[["age_specific_rate"]])
              result[["age_specific_projection"]] <- rbind(result[["age_specific_projection"]],temp[["age_specific_projection"]])
              
              i = i + 1
            }
          }
        }
      }
    }
  }
  return(result)
  
}

.gbd_nordpred_prediction <- function(data,data_NA,ages,ages_2,dat,population,std,measure_name,cause_name,
                                   location_name,rei_name,By_sex,predyear,full_age_adjusted,basicyear,endyear,predyear1) {
  
  
  if (is.null(rei_name) == T) {
    Male_data <- subset(data,age != 'Age-standardized' &
                          age != 'All ages' &
                          sex == 'Male' &
                          metric == 'Number' &
                          measure == measure_name &
                          location  == location_name &
                          cause  == cause_name)
    
    Female_data <- subset(data,age != 'Age-standardized' &
                            age != 'All ages' &
                            sex == 'Female' &
                            metric == 'Number' &
                            measure == measure_name &
                            location  == location_name &
                            cause  == cause_name)
    
    Both_data <- subset(data,age != 'Age-standardized' &
                          age != 'All ages' &
                          sex == 'Both' &
                          metric == 'Number' &
                          measure == measure_name &
                          location  == location_name &
                          cause  == cause_name)
    
    data_NA <- subset(data_NA,age != 'Age-standardized' &
                        age != 'All ages' &
                        metric == 'Number' &
                        measure == measure_name &
                        location  == location_name &
                        cause  == cause_name)
  } else {
    Male_data <- subset(data,age != 'Age-standardized' &
                          age != 'All ages' &
                          sex == 'Male' &
                          metric == 'Number' &
                          measure == measure_name &
                          location  == location_name &
                          cause  == cause_name &
                          rei == rei_name)
    
    Female_data <- subset(data,age != 'Age-standardized' &
                            age != 'All ages' &
                            sex == 'Female' &
                            metric == 'Number' &
                            measure == measure_name &
                            location  == location_name &
                            cause  == cause_name &
                            rei == rei_name)
    
    Both_data <- subset(data,age != 'Age-standardized' &
                          age != 'All ages' &
                          sex == 'Both' &
                          metric == 'Number' &
                          measure == measure_name &
                          location  == location_name &
                          cause  == cause_name &
                          rei == rei_name)
    
    data_NA <- subset(data_NA,age != 'Age-standardized' &
                        age != 'All ages' &
                        metric == 'Number' &
                        measure == measure_name &
                        location  == location_name &
                        cause  == cause_name &
                        rei == rei_name)
  }
  data_all <- rbind(Male_data, Female_data, Both_data)
  
  
  age_0 <- unique(data_NA$age)
  
  
  unique_age <- unique(Both_data$age)
  
  
  unique_age <- sub('95 plus',replacement = '95 to 99', unique_age)
  unique_age2 <- matrix(as.numeric(unlist(strsplit(unique(unique_age),split="to"))),
                        ncol=2,byrow=T)[,2]
  unique_age <- matrix(as.numeric(unlist(strsplit(unique(unique_age),split="to"))),
                       ncol=2,byrow=T)[,1]
  
  diff_age <- unique_age2 - unique_age
  
  if (length(unique(diff_age)) !=1 | max(diff_age) !=4 |min(diff_age) !=4){
    stop("Age groups must be in five-year intervals; ASR and all-age records may also be present.")
    
  }
  
  
  
  if (sum(location_name %in% unique(data$location)) !=length(location_name)){
    stop("One or more location_name labels are not present in the input data.", call. = FALSE)
  }
  
  if (sum(measure_name %in% unique(data$measure)) !=length(measure_name)){
    stop("One or more measure_name labels are not present in the input data.", call. = FALSE)
  }
  
  if (sum(cause_name %in% unique(data$cause)) !=length(cause_name)){
    stop("One or more cause_name labels are not present in the input data.", call. = FALSE)
  }
  
  if (is.null(rei_name) & 'rei' %in% names(data)){
    stop("Specify rei_name for inputs with a rei column.", call. = FALSE)
  }
  
  if (!is.null(rei_name)){
    if (sum(rei_name %in% unique(data$rei)) !=length(rei_name)){
      stop("One or more rei_name labels are not present in the input data.", call. = FALSE)
    }
  }
  
  startage <- min(unique_age)
  endage <- max(unique_age)
  ages <- c()
  for (j in seq(from=startage,to=endage,by=5)) {
    ages <- c(ages,paste(j,j+4,sep=' to '))
  }
  if (endage == 95) {
    ages[length(ages)] <- '95 plus' }

  wstand <- std$percent
  std2 <- subset(std,age %in% ages)
  std2$percent <- std2$std_population/sum(std2$std_population)
  wstand_2 <- std2$percent
  
  population_Male <- subset(population, location == location_name & sex == 'Male')
  population_Female <- subset(population, location == location_name & sex == 'Female')
  
  population_Male_n <- reshape2::dcast(data = population_Male %>% as.data.frame() %>% unique(),
                                       age ~ year,
                                       value.var = c("val")) %>% dplyr::mutate(age=factor(age,levels = ages_2,ordered = T)) %>%
    dplyr::arrange(age) %>% dplyr::select(-age) %>% as.data.frame()
  population_Female_n <- reshape2::dcast(data = population_Female %>% as.data.frame() %>% unique(),
                                         age ~ year,
                                         value.var = c("val")) %>% dplyr::mutate(age=factor(age,levels = ages_2,ordered = T)) %>%
    dplyr::arrange(age) %>% dplyr::select(-age) %>% as.data.frame()
  
  if (ncol(population_Male_n)<=30){
    stop(paste0("Projection population data are unavailable for ", location_name, "."), call. = FALSE)
  }
  
  rownames(population_Male_n) <- ages_2
  rownames(population_Female_n) <- ages_2
  
  population_Male_n <- apply(population_Male_n, c(1,2), as.numeric) %>% as.data.frame()
  population_Female_n <- apply(population_Female_n, c(1,2), as.numeric) %>% as.data.frame()
  
  population_Both_n <- population_Female_n + population_Male_n
  
  
  
  Male_data_n <- reshape2::dcast(data = Male_data %>% as.data.frame() %>% unique(),
                                 age ~ year,
                                 value.var = c("val")) %>% dplyr::arrange(age=factor(age,levels = ages,ordered = T))
  
  Female_data_n <- reshape2::dcast(data = Female_data %>% as.data.frame() %>% unique(),
                                   age ~ year,
                                   value.var = c("val"))  %>% dplyr::arrange(age=factor(age,levels = ages,ordered = T))
  
  Both_data_n <- reshape2::dcast(data = Both_data %>% as.data.frame() %>% unique(),
                                 age ~ year,
                                 value.var = c("val"))  %>% dplyr::arrange(age=factor(age,levels = ages,ordered = T))
  ages_add <- setdiff(ages_2,ages)
  if (length(ages_add) > 0) {
    ages_add_data <- matrix(0,nrow = length(ages_add),ncol=ncol(Male_data_n)) %>% as.data.frame()
    names(ages_add_data) <- c("age",basicyear:((basicyear-1)+length(unique(Male_data$year))))
    ages_add_data[,1] <- ages_add
    Male_data_n <- rbind(ages_add_data,Male_data_n)
    Female_data_n <- rbind(ages_add_data,Female_data_n)
    Both_data_n <- rbind(ages_add_data,Both_data_n)
  }
  
  
  Male_data_n <- Male_data_n %>% dplyr::mutate(age=factor(age,levels = ages_2,ordered = T)) %>% dplyr::arrange(age) %>% dplyr::select(-1) %>% as.data.frame()
  Female_data_n <- Female_data_n %>%  dplyr::mutate(age=factor(age,levels = ages_2,ordered = T)) %>% dplyr::arrange(age) %>% dplyr::select(-1) %>% as.data.frame()
  Both_data_n <- Both_data_n %>%  dplyr::mutate(age=factor(age,levels = ages_2,ordered = T)) %>% dplyr::arrange(age) %>% dplyr::select(-1) %>% as.data.frame()
  
  Female_data_g <- datagg(Female_data_n,population_Female_n,nagg=5)[["cases"]]
  population_Female_g <- datagg(Female_data_n,population_Female_n,nagg=5)[["pyr"]]
  
  Both_data_g <- datagg(Both_data_n,population_Both_n,nagg=5)[["cases"]]
  population_Both_g <- datagg(Both_data_n,population_Both_n,nagg=5)[["pyr"]]
  
  Male_data_g <- datagg(Male_data_n,population_Male_n,nagg=5)[["cases"]]
  population_Male_g <- datagg(Male_data_n,population_Male_n,nagg=5)[["pyr"]]
  
  rownames(Male_data_g) <- ages_2
  rownames(population_Male_g) <- ages_2
  
  rownames(Both_data_g) <- ages_2
  rownames(population_Both_g) <- ages_2
  
  rownames(Female_data_g) <- ages_2
  rownames(population_Female_g) <- ages_2
  
  
  
  
  n <- which(matrix(as.numeric(unlist(strsplit(sub('plus','to 100',ages_2),split=" to "))),ncol=2,byrow=T)[,1] %in% startage)
  
  Male_nordpred <- nordpred(cases=Male_data_g,pyr=population_Male_g,startestage=n,startuseage=n,
                           cuttrend = c(0, .25, .5, .75, .75), linkfunc = "power5",
                           recent = NULL)
  
  if (length(age_0)>0){
    for (i in which(ages_2 %in% age_0)) {
      Male_nordpred[["predictions"]][i,] <- matrix(0,nrow=1,ncol=ncol(Male_nordpred[["predictions"]])) %>% as.data.frame()
    }
  }
  Male_age_specific_rate <- nordpred.getpred(Male_nordpred)
  Male_age_specific_rate_n <- asrpy(Male_age_specific_rate,Male_data_n,population_Male_n,startp=(basicyear+length(unique(Male_data$year))),nagg=5)
  Male_age_count_n <- Male_age_specific_rate_n*population_Male_n/100000
  
  
  Female_nordpred <- nordpred(Female_data_g,population_Female_g,startestage=n,startuseage=n,
                             cuttrend = c(0, .25, .5, .75, .75), linkfunc = "power5",
                             recent = NULL)
  
  
  if (length(age_0)>0){
    for (i in which(ages_2 %in% age_0)) {
      Female_nordpred[["predictions"]][i,] <- matrix(0,nrow=1,ncol=ncol(Female_nordpred[["predictions"]])) %>% as.data.frame()
    }
  }
  
  Female_age_specific_rate <- nordpred.getpred(Female_nordpred)
  Female_age_specific_rate_n <- asrpy(Female_age_specific_rate,Female_data_n,population_Female_n,startp=(basicyear+length(unique(Female_data$year))),nagg=5)
  Female_age_count_n <- Female_age_specific_rate_n*population_Female_n/100000
  
  
  
  if (By_sex==F) {
    Both_nordpred <- nordpred(Both_data_g,population_Both_g,startestage=n,startuseage=n,
                             cuttrend = c(0, .25, .5, .75, .75), linkfunc = "power5",
                             recent = NULL)
    
    if (length(age_0)>0){
      for (i in which(ages_2 %in% age_0)) {
        Both_nordpred[["predictions"]][i,] <- matrix(0,nrow=1,ncol=ncol(Both_nordpred[["predictions"]])) %>% as.data.frame()
      }
    }
    
    Both_age_specific_rate <- nordpred.getpred(Both_nordpred)
    Both_age_specific_rate_n <- asrpy(Both_age_specific_rate,Both_data_n,population_Both_n,startp=(basicyear+length(unique(Both_data$year))),nagg=5)
    Both_age_count_n <- Both_age_specific_rate_n*population_Both_n/100000
  } else {
    Both_age_count_n <- Female_age_count_n + Male_age_count_n
    Both_age_specific_rate_n <- Both_age_count_n/population_Both_n*100000
  }
  
  
  
  
  if (full_age_adjusted==T) {
    Male_asr_sum <- asry(Male_age_specific_rate_n, population_Male_n, standpop=wstand) %>% as.data.frame()
    Female_asr_sum <- asry(Female_age_specific_rate_n, population_Female_n, standpop=wstand) %>% as.data.frame()
    Both_asr_sum <- asry(Both_age_specific_rate_n, population_Both_n, standpop=wstand) %>% as.data.frame()
  } else {
    loc_age <- which((ages_2 %in% ages))
    Male_asr_sum <- asry(Male_age_specific_rate_n[loc_age,], population_Male_n[loc_age,], standpop=wstand_2) %>% as.data.frame()
    Female_asr_sum <- asry(Female_age_specific_rate_n[loc_age,], population_Female_n[loc_age,], standpop=wstand_2) %>% as.data.frame()
    Both_asr_sum <- asry(Both_age_specific_rate_n[loc_age,], population_Both_n[loc_age,], standpop=wstand_2) %>% as.data.frame()
  }
  
  Male_asr_sum$year <- rownames(Male_asr_sum) %>% as.numeric()
  Male_asr_sum$sex <- 'Male'
  Male_asr_sum$location <- location_name
  Male_asr_sum$cause <- cause_name
  Male_asr_sum$measure <- measure_name
  Male_age_specific_rate_n$age <- rownames(Male_age_specific_rate_n)
  Male_age_specific_rate_n$sex <- 'Male'
  Male_age_specific_rate_n$location <- location_name
  Male_age_specific_rate_n$cause <- cause_name
  Male_age_specific_rate_n$measure <- measure_name
  Male_age_count_n$age <- rownames(Male_age_count_n)
  Male_age_count_n$sex <- 'Male'
  Male_age_count_n$location <- location_name
  Male_age_count_n$cause <- cause_name
  Male_age_count_n$measure <- measure_name
  
  Female_asr_sum$year <- rownames(Female_asr_sum) %>% as.numeric()
  Female_asr_sum$sex <- 'Female'
  Female_asr_sum$location <- location_name
  Female_asr_sum$cause <- cause_name
  Female_asr_sum$measure <- measure_name
  Female_age_specific_rate_n$age <- rownames(Female_age_specific_rate_n)
  Female_age_specific_rate_n$sex <- 'Female'
  Female_age_specific_rate_n$location <- location_name
  Female_age_specific_rate_n$cause <- cause_name
  Female_age_specific_rate_n$measure <- measure_name
  Female_age_count_n$age <- rownames(Female_age_count_n)
  Female_age_count_n$sex <- 'Female'
  Female_age_count_n$location <- location_name
  Female_age_count_n$cause <- cause_name
  Female_age_count_n$measure <- measure_name
  
  Both_asr_sum$year <- rownames(Both_asr_sum) %>% as.numeric()
  Both_asr_sum$sex <- 'Both'
  Both_asr_sum$location <- location_name
  Both_asr_sum$cause <- cause_name
  Both_asr_sum$measure <- measure_name
  Both_age_specific_rate_n$age <- rownames(Both_age_specific_rate_n)
  Both_age_specific_rate_n$sex <- 'Both'
  Both_age_specific_rate_n$location <- location_name
  Both_age_specific_rate_n$cause <- cause_name
  Both_age_specific_rate_n$measure <- measure_name
  Both_age_count_n$age <- rownames(Both_age_count_n)
  Both_age_count_n$sex <- 'Both'
  Both_age_count_n$location <- location_name
  Both_age_count_n$cause <- cause_name
  Both_age_count_n$measure <- measure_name
  
  
  val_num <- 5
  if(is.null(rei_name)==F) {
    Male_asr_sum$rei <- rei_name
    Male_age_specific_rate_n$rei <- rei_name
    Male_age_count_n$rei <- rei_name
    
    Both_asr_sum$rei <- rei_name
    Both_age_specific_rate_n$rei <- rei_name
    Both_age_count_n$rei <- rei_name
    
    Female_asr_sum$rei <- rei_name
    Female_age_specific_rate_n$rei <- rei_name
    Female_age_count_n$rei <- rei_name
    val_num <- 6
  }
  
  ## Male
  Male_age_count_s <- Male_age_count_n %>%
    tidyr::pivot_longer(1:(ncol(Male_age_count_n)-val_num),
                        names_to = "year",
                        values_to = 'val') %>%
    dplyr::filter(age %in% ages)
  
  Male_age_specific_rate_s <- Male_age_specific_rate_n %>%
    tidyr::pivot_longer(1:(ncol(Male_age_specific_rate_n)-val_num),
                        names_to = "year",
                        values_to = 'val') %>%
    dplyr::filter(age %in% ages)
  
  ## Female
  Female_age_count_s <- Female_age_count_n %>%
    tidyr::pivot_longer(1:(ncol(Female_age_count_n)-val_num),
                        names_to = "year",
                        values_to = 'val') %>%
    dplyr::filter(age %in% ages)
  
  Female_age_specific_rate_s <- Female_age_specific_rate_n %>%
    tidyr::pivot_longer(1:(ncol(Female_age_specific_rate_n)-val_num),
                        names_to = "year",
                        values_to = 'val') %>%
    dplyr::filter(age %in% ages)
  
  ## Both
  Both_age_count_s <- Both_age_count_n %>%
    tidyr::pivot_longer(1:(ncol(Both_age_count_n)-val_num),
                        names_to = "year",
                        values_to = 'val') %>%
    dplyr::filter(age %in% ages)
  
  Both_age_specific_rate_s <- Both_age_specific_rate_n %>%
    tidyr::pivot_longer(1:(ncol(Both_age_specific_rate_n)-val_num),
                        names_to = "year",
                        values_to = 'val') %>%
    dplyr::filter(age %in% ages)
  
  ######## summary
  age_specific_rate <- rbind(Male_age_specific_rate_s,Female_age_specific_rate_s,Both_age_specific_rate_s)
  age_count <- rbind(Male_age_count_s,Female_age_count_s,Both_age_count_s)
  asr_sum <- rbind(Male_asr_sum,Female_asr_sum,Both_asr_sum)
  
  if (is.null(rei_name)==T){
    age_specific_rate <- age_specific_rate %>% dplyr::select(measure,location,cause,sex,age,year,val)
    age_count <- age_count %>% dplyr::select(measure,location,cause,sex,age,year,val)
    asr_sum <- asr_sum %>% dplyr::select(measure,location,cause,sex,year,ASR,crude_rate,case)
  } else {
    age_specific_rate <- age_specific_rate %>% dplyr::select(measure,location,cause,rei,sex,age,year,val)
    age_count <- age_count %>% dplyr::select(measure,location,cause,rei,sex,age,year,val)
    asr_sum <- asr_sum %>% dplyr::select(measure,location,cause,rei,sex,year,ASR,crude_rate,case)
  }
  #### output
  
  asr_sum <- asr_sum %>% dplyr::filter(year <= predyear1)
  age_specific_rate <- age_specific_rate %>% dplyr::filter(year <= predyear1) %>% as.data.frame()
  age_count <- age_count %>% dplyr::filter(year <= predyear1) %>% as.data.frame()
  asr_sum <- asr_sum %>% as.data.frame()
  
  row.names(asr_sum) <- 1:nrow(asr_sum)
  row.names(age_specific_rate) <- 1:nrow(age_specific_rate)
  row.names(age_count) <- 1:nrow(age_count)
  
  asr_sum$year <- as.numeric(asr_sum$year)
  age_specific_rate$year <- as.numeric(age_specific_rate$year)
  age_count$year <- as.numeric(age_count$year)
  
  
  result <- list(ASR_Number = asr_sum,
                 age_specific_rate = age_specific_rate,
                 age_specific_projection = age_count,
                 label = 'nordpred')
  invisible()
  return(result)
}


# Projection analysis ------------------------------------------------------

# User settings ------------------------------------------------------------

root_dir <- normalizePath(Sys.getenv("AGEING_ROOT", unset = getwd()), mustWork = TRUE)
demo_mode <- identical(Sys.getenv("AGEING_DEMO", unset = "FALSE"), "TRUE")

input_dir <- if (demo_mode) file.path(root_dir, "data", "demo") else file.path(root_dir, "data", "processed")
output_dir <- if (demo_mode) file.path(root_dir, "outputs", "demo") else file.path(root_dir, "outputs", "analysis")
upstream_dir <- output_dir

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

overall_projection_input_file <- file.path(upstream_dir, "total_age_related_daly_by_age.csv")
system_projection_input_file <- file.path(upstream_dir, "system_age_related_daly_by_age.csv")
population_history_file <- resolve_input_file(input_dir, "population_history_2021.csv")
population_projection_file <- resolve_input_file(input_dir, "population_projection_who.csv")
location_group_file <- file.path(input_dir, "projection_location_groups.csv")
standard_population_file <- file.path(root_dir, "data", "metadata", "std_GBD2021.csv")
rate_output_file <- file.path(output_dir, "projected_age_specific_rates.csv")
equivalent_age_output_file <- file.path(output_dir, "projected_equivalent_age_2022_2040.csv")

target_measure <- "DALYs (Disability-Adjusted Life Years)"
target_predyear <- 2040
n_simulations <- if (demo_mode) 10 else 300
random_seed <- 123
gbd_edition <- 2021

analysis_age_groups <- c(
  "25 to 29", "30 to 34", "35 to 39", "40 to 44", "45 to 49",
  "50 to 54", "55 to 59", "60 to 64", "65 to 69", "70 to 74",
  "75 to 79", "80 to 84", "85 to 89", "90 to 94"
)

# Helper functions ---------------------------------------------------------

assert_required_columns <- function(data, required_columns) {
  missing_columns <- setdiff(required_columns, names(data))
  if (length(missing_columns) > 0) {
    stop(
      "Input data are missing required columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
}

load_projection_input <- function(overall_path, system_path) {
  overall_data <- read_csv(overall_path, show_col_types = FALSE) %>%
    mutate(cause = "Age-related diseases")

  assert_required_columns(
    overall_data,
    c("measure", "cause", "location", "sex", "year", "age", "metric", "val", "lower", "upper")
  )

  projection_data <- overall_data

  if (file.exists(system_path)) {
    system_data <- read_csv(system_path, show_col_types = FALSE) %>%
      rename(cause = category)

    assert_required_columns(
      system_data,
      c("measure", "cause", "location", "sex", "year", "age", "metric", "val", "lower", "upper")
    )

    projection_data <- bind_rows(projection_data, system_data)
  }

  projection_data %>%
    filter(
      measure == target_measure,
      age %in% analysis_age_groups
    ) %>%
    mutate(
      year = as.integer(year),
      val = as.numeric(val),
      lower = as.numeric(lower),
      upper = as.numeric(upper)
    ) %>%
    mutate(
      val = if_else(metric == "Number", round(pmax(val, 0)), val),
      lower = if_else(metric == "Number", round(pmax(lower, 0)), lower),
      upper = if_else(metric == "Number", round(pmax(upper, 0)), upper)
    )
}

load_population_file <- function(path) {
  population_data <- read_csv(path, show_col_types = FALSE)
  assert_required_columns(population_data, c("location", "sex", "year", "age", "val"))
  population_data
}

load_standard_population <- function(path) {
  standard_population_data <- read_csv(path, show_col_types = FALSE, name_repair = "minimal")
  assert_required_columns(standard_population_data, c("age", "std_population"))

  standard_population_data %>%
    select(age, std_population) %>%
    as.data.frame()
}

load_location_groups <- function(path) {
  if (!file.exists(path)) {
    return(NULL)
  }

  group_data <- read_csv(path, show_col_types = FALSE)
  assert_required_columns(group_data, c("location", "aggregate_location"))

  group_data %>%
    distinct(location, aggregate_location) %>%
    filter(!is.na(location), !is.na(aggregate_location))
}

validate_projection_input <- function(data) {
  if (nrow(data) == 0) {
    stop("The filtered input dataset is empty.", call. = FALSE)
  }

  if (!all(c("Male", "Female", "Both") %in% unique(data$sex))) {
    stop("The input dataset must contain Male, Female, and Both rows.", call. = FALSE)
  }

  if (!"Number" %in% unique(data$metric)) {
    stop("Projection input must include metric == 'Number' rows.", call. = FALSE)
  }

  if (!"Rate" %in% unique(data$metric)) {
    stop("Projection input must include metric == 'Rate' rows for equivalent-age benchmark mapping.", call. = FALSE)
  }

  unexpected_ages <- setdiff(unique(data$age), analysis_age_groups)
  if (length(unexpected_ages) > 0) {
    stop(
      "Unexpected ages are present after filtering: ",
      paste(unexpected_ages, collapse = ", "),
      call. = FALSE
    )
  }

  historical_counts <- data %>%
    filter(metric == "Number") %>%
    distinct(cause, location, sex, year, age)

  incomplete_panels <- historical_counts %>%
    count(cause, location, sex, year, name = "n_age_groups") %>%
    filter(n_age_groups != length(analysis_age_groups))

  if (nrow(incomplete_panels) > 0) {
    stop(
      "Historical Number data do not contain a complete age panel for every location-sex-year combination.",
      call. = FALSE
    )
  }

  invisible(data)
}

prepare_location_dataset <- function(data, location_value) {
  location_data <- data %>%
    filter(location == location_value) %>%
    arrange(metric, sex, year, factor(age, levels = analysis_age_groups))

  if (nrow(location_data) == 0) {
    stop("No data were found for location: ", location_value, call. = FALSE)
  }

  location_data
}

build_age_period_matrix <- function(count_data) {
  count_data %>%
    select(sex, year, age, val) %>%
    mutate(age = factor(age, levels = analysis_age_groups)) %>%
    arrange(sex, year, age) %>%
    pivot_wider(
      id_cols = c(sex, year),
      names_from = age,
      values_from = val
    ) %>%
    arrange(sex, year)
}

run_nordpred_projection <- function(location_data) {
  projection_fit <- NULL
  suppressMessages(capture.output({
    projection_fit <- gbd_nordpred_prediction(
      data = location_data,
      measure_name = target_measure,
      cause_name = unique(location_data$cause),
      location_name = unique(location_data$location),
      rei_name = NULL,
      By_sex = FALSE,
      predyear = target_predyear,
      full_age_adjusted = FALSE,
      pop_predict = "WHO",
      gbd_edition = gbd_edition,
      population_history = population_history,
      population_projection = population_projection,
      standard_population = standard_population
    )
  }))

  historical_end_year <- max(location_data$year[location_data$metric == "Number"], na.rm = TRUE)

  projection_fit$age_specific_rate %>%
    as_tibble() %>%
    filter(
      age %in% analysis_age_groups,
      year > historical_end_year,
      year <= target_predyear
    ) %>%
    transmute(
      cause = unique(location_data$cause),
      location,
      sex,
      age,
      year = as.integer(year),
      point_estimate = as.numeric(val)
    ) %>%
    arrange(sex, year, factor(age, levels = analysis_age_groups))
}

simulate_historical_counts <- function(location_data) {
  simulated_data <- location_data
  index_to_simulate <- which(simulated_data$metric == "Number")

  simulated_data$val[index_to_simulate] <- rpois(
    n = length(index_to_simulate),
    lambda = pmax(simulated_data$val[index_to_simulate], 0)
  )

  simulated_data$lower[index_to_simulate] <- NA_real_
  simulated_data$upper[index_to_simulate] <- NA_real_

  simulated_data
}

estimate_simulation_intervals <- function(location_data, n_simulations, seed) {
  set.seed(seed)

  simulation_results <- vector("list", n_simulations)

  for (i in seq_len(n_simulations)) {
    simulated_input <- simulate_historical_counts(location_data)
    simulation_results[[i]] <- run_nordpred_projection(simulated_input) %>%
      transmute(cause, location, sex, age, year, simulation_value = point_estimate)
  }

  bind_rows(simulation_results) %>%
    group_by(cause, location, sex, age, year) %>%
    summarise(
      point_estimate = mean(simulation_value, na.rm = TRUE),
      lower = quantile(simulation_value, 0.025, na.rm = TRUE),
      upper = quantile(simulation_value, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
}

project_one_location <- function(location_data, n_simulations, seed) {
  count_matrix <- build_age_period_matrix(location_data %>% filter(metric == "Number"))

  if (nrow(count_matrix) == 0) {
    stop("No historical Number rows were available for projection.", call. = FALSE)
  }

  estimate_simulation_intervals(location_data, n_simulations, seed) %>%
    arrange(cause, location, sex, year, factor(age, levels = analysis_age_groups))
}

run_projection_analysis <- function(data, n_simulations, seed) {
  unit_list <- data %>%
    distinct(cause, location) %>%
    arrange(cause, location)

  location_results <- lapply(seq_len(nrow(unit_list)), function(i) {
    location_data <- data %>%
      filter(cause == unit_list$cause[i]) %>%
      prepare_location_dataset(unit_list$location[i])
    project_one_location(
      location_data = location_data,
      n_simulations = n_simulations,
      seed = seed + i - 1
    )
  })

  bind_rows(location_results) %>%
    arrange(cause, location, sex, year, factor(age, levels = analysis_age_groups))
}

aggregate_projected_rates <- function(projected_rates, population_projection, location_groups) {
  if (is.null(location_groups) || nrow(location_groups) == 0) {
    return(projected_rates)
  }

  population_for_weights <- population_projection %>%
    filter(age %in% analysis_age_groups) %>%
    transmute(
      location,
      sex,
      year = as.integer(year),
      age,
      val = as.numeric(val)
    )

  both_population <- population_for_weights %>%
    filter(sex %in% c("Male", "Female")) %>%
    group_by(location, year, age) %>%
    summarise(val = sum(val, na.rm = TRUE), .groups = "drop") %>%
    mutate(sex = "Both")

  aggregation_weights <- bind_rows(population_for_weights, both_population) %>%
    inner_join(location_groups, by = "location") %>%
    transmute(
      location,
      aggregate_location,
      sex,
      year = as.integer(year),
      age,
      population = as.numeric(val)
    )

  aggregated_rates <- projected_rates %>%
    inner_join(aggregation_weights, by = c("location", "sex", "year", "age")) %>%
    group_by(cause, aggregate_location, sex, age, year) %>%
    summarise(
      point_estimate = weighted.mean(point_estimate, population, na.rm = TRUE),
      lower = weighted.mean(lower, population, na.rm = TRUE),
      upper = weighted.mean(upper, population, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    transmute(
      cause,
      location = aggregate_location,
      sex,
      age,
      year,
      point_estimate,
      lower,
      upper
    )

  if (nrow(aggregated_rates) == 0) {
    return(projected_rates)
  }

  aggregated_rates <- aggregated_rates %>%
    anti_join(
      projected_rates %>% distinct(cause, location, sex, age, year),
      by = c("cause", "location", "sex", "age", "year")
    )

  bind_rows(projected_rates, aggregated_rates) %>%
    arrange(cause, location, sex, year, factor(age, levels = analysis_age_groups))
}

compute_midpoint_age <- function(age_group) {
  values <- as.numeric(stringr::str_extract_all(age_group, "\\d+")[[1]])
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

get_projection_benchmarks <- function(historical_data) {
  benchmarks <- historical_data %>%
    filter(
      metric == "Rate",
      location == "Global",
      sex == "Both",
      year == 2021,
      age %in% c("60 to 64", "65 to 69")
    ) %>%
    group_by(cause) %>%
    summarise(
      benchmark_rate = mean(val, na.rm = TRUE),
      benchmark_lower = mean(lower, na.rm = TRUE),
      benchmark_upper = mean(upper, na.rm = TRUE),
      .groups = "drop"
    )

  if (nrow(benchmarks) == 0 || any(is.na(benchmarks$benchmark_rate))) {
    stop("Projection benchmarks require Global Both 2021 Rate rows for each cause.", call. = FALSE)
  }

  benchmarks
}

calculate_projected_equivalent_age <- function(projected_rates, historical_data) {
  benchmarks <- get_projection_benchmarks(historical_data)
  grouping_columns <- c("cause", "location", "sex", "year")

  projected_prepared <- projected_rates %>%
    filter(age %in% analysis_age_groups) %>%
    mutate(age = factor(age, levels = analysis_age_groups, ordered = TRUE)) %>%
    left_join(benchmarks, by = "cause")

  analysis_units <- projected_prepared %>%
    distinct(across(all_of(grouping_columns)))

  interpolatable_data <- projected_prepared %>%
    group_by(across(all_of(grouping_columns))) %>%
    arrange(age, .by_group = TRUE) %>%
    mutate(
      prev_group = as.character(age),
      prev_rate = point_estimate,
      next_group = lead(as.character(age)),
      next_rate = lead(point_estimate)
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
      )
    ) %>%
    select(all_of(grouping_columns), prev_group, next_group, prev_rate, next_rate, equivalent_age_raw)

  interval_long <- projected_prepared %>%
    select(all_of(grouping_columns), age, lower, upper, benchmark_rate) %>%
    pivot_longer(
      cols = c(lower, upper),
      names_to = "interval_type",
      values_to = "interval_rate"
    ) %>%
    group_by(across(all_of(grouping_columns)), interval_type) %>%
    arrange(age, .by_group = TRUE) %>%
    mutate(
      prev_group = as.character(age),
      prev_rate = interval_rate,
      next_group = lead(as.character(age)),
      next_rate = lead(interval_rate)
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
      interval_age = mapply(
        interpolate_equivalent_age,
        prev_rate,
        next_rate,
        prev_midpoint_age,
        next_midpoint_age,
        benchmark_rate
      )
    ) %>%
    select(all_of(grouping_columns), interval_type, interval_age) %>%
    pivot_wider(names_from = interval_type, values_from = interval_age)

  boundary_data <- projected_prepared %>%
    filter(age %in% c("25 to 29", "90 to 94")) %>%
    select(all_of(grouping_columns), age, point_estimate) %>%
    pivot_wider(
      id_cols = all_of(grouping_columns),
      names_from = age,
      values_from = point_estimate,
      values_fn = mean
    ) %>%
    rename(
      rate_at_start = `25 to 29`,
      rate_at_end = `90 to 94`
    )

  analysis_units %>%
    left_join(interpolatable_data, by = grouping_columns) %>%
    left_join(interval_long, by = grouping_columns) %>%
    left_join(boundary_data, by = grouping_columns) %>%
    left_join(benchmarks, by = "cause") %>%
    mutate(
      lower_boundary = !is.na(rate_at_start) & rate_at_start > benchmark_rate,
      upper_boundary = !is.na(rate_at_end) & rate_at_end < benchmark_rate,
      lower_bound_raw = pmin(lower, upper),
      upper_bound_raw = pmax(lower, upper),
      equivalent_age = as_display_value(equivalent_age_raw, lower_boundary, upper_boundary),
      lower_bound = as_display_value(lower_bound_raw, lower_boundary, upper_boundary),
      upper_bound = as_display_value(upper_bound_raw, lower_boundary, upper_boundary),
      equivalent_age_numeric = as_numeric_value(equivalent_age_raw, lower_boundary, upper_boundary),
      lower_bound_numeric = as_numeric_value(lower_bound_raw, lower_boundary, upper_boundary),
      upper_bound_numeric = as_numeric_value(upper_bound_raw, lower_boundary, upper_boundary)
    ) %>%
    select(
      cause, location, sex, year,
      benchmark_rate, benchmark_lower, benchmark_upper,
      prev_group, next_group, prev_rate, next_rate,
      equivalent_age, lower_bound, upper_bound,
      equivalent_age_numeric, lower_bound_numeric, upper_bound_numeric
    ) %>%
    arrange(cause, location, sex, year)
}

# Final execution block ----------------------------------------------------

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

projection_input <- load_projection_input(
  overall_projection_input_file,
  system_projection_input_file
)
population_history <- load_population_file(population_history_file)
population_projection <- load_population_file(population_projection_file)
location_groups <- load_location_groups(location_group_file)
standard_population <- load_standard_population(standard_population_file)

validate_projection_input(projection_input)

projection_output <- run_projection_analysis(
  data = projection_input,
  n_simulations = n_simulations,
  seed = random_seed
)

projection_output <- aggregate_projected_rates(
  projected_rates = projection_output,
  population_projection = population_projection,
  location_groups = location_groups
)

projected_equivalent_age <- calculate_projected_equivalent_age(projection_output, projection_input)

write_csv(projection_output, rate_output_file)
write_csv(projected_equivalent_age, equivalent_age_output_file)
