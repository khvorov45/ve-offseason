# Script to simulate data for the investigation of VE estimates when
# flu activity is low
#
# Arseniy Khvorov
# Created 2019/10/04
# Last edit 2019/11/12

library(tidyverse)
library(furrr) # Like purrr but with parallel support

plan(multiprocess) # Parallel setup

sim_folder <- "sim"

# Functions ===================================================================

# Categorical variable generation (e.g. infection status)
rcat <- function(n, prob) {
  prob <- c(prob, 1 - sum(prob))
  mat <- rmultinom(n, 1, prob)
  out <- rep(0L, n)
  reduce(1:nrow(mat), ~ if_else(mat[.y, ] == 1, .y, .x), .init = out)
}

# Simulation function
simulate_pop <- function(nsam = 1000,
                         pvac = 0.5,
                         pflu = 0.1,
                         ve = 0.5,
                         pnonflu = 0.1,
                         prisk = 1,
                         sens_vac = 1,
                         spec_vac = 1,
                         sens_flu = 1,
                         spec_flu = 1,
                         simseed = NULL) {
  
  if (!is.null(simseed)) set.seed(simseed)
  
  # Work out the vaccination record
  calc_vac_record <- function(.tbl, .key, sens_vac, spec_vac) {
    if (.key$vaccinated) prob <- sens_vac
    else prob <- 1 - spec_vac
    .tbl %>% mutate(vac_record = rbinom(n(), 1, prob) %>% as.logical())
  }
  
  # Work out infection status
  calc_status <- function(.tbl, .key, pflu, ve, pnonflu) {
    if (!.key$at_risk) prob <- c(0, pnonflu)
    else if (.key$vaccinated) prob <- c(pflu * (1 - ve), pnonflu)
    else prob <- c(pflu, pnonflu)
    .tbl %>% mutate(status = rcat(n(), prob))
  }
  
  # Work out the test result
  calc_test_result <- function(.tbl, .key, sens_flu, spec_flu) {
    if (.key$status == "Flu") prob <- sens_flu
    else prob <- 1 - spec_flu
    .tbl %>% mutate(test_result = rbinom(n(), 1, prob) %>% as.logical())
  }
  
  # Simulate data
  tibble(.rows = nsam) %>% 
    mutate(
      vaccinated = rbinom(n(), 1, pvac) %>% as.logical(),
      at_risk = rbinom(n(), 1, prisk) %>% as.logical()
    ) %>%
    group_by(vaccinated) %>%
    group_modify(~ calc_vac_record(.x, .y, sens_vac, spec_vac)) %>%
    group_by(vaccinated, at_risk) %>%
    group_modify(~ calc_status(.x, .y, pflu, ve, pnonflu)) %>%
    ungroup() %>%
    mutate(
      status = recode(status, "3" = "Healthy", "2" = "Nonflu", "1" = "Flu")
    ) %>%
    group_by(status) %>%
    group_modify(~ calc_test_result(.x, .y, sens_flu, spec_flu)) %>%
    ungroup()
}

# Calculates the number of vaccinated/unvaccinated cases/controls
calc_counts <- function(.tbl, type) {
  type <- match.arg(type, c("tn", "cc"))
  if (type == "tn") {
    .tbl <- .tbl %>% 
      mutate(
        case = case_when(
          status == "Healthy" ~ NA,
          test_result ~ TRUE,
          !test_result ~ FALSE
        )
      )
  } else .tbl <- .tbl %>% mutate(case = test_result)
  .tbl %>%
    filter(!is.na(case)) %>%
    count(vac_record, case, name = "n_group") %>%
    mutate(
      study = if_else(type == "tn", "Test negative", "Case control"), 
    )
}

# Calculates vaccine efficacy given the number of 
# vaccinated/unvaccinated cases/controls
calc_ve <- function(.tbl) {
  n_group <- sum(.tbl$n_group)
  .tbl <- .tbl %>%
    mutate(case = if_else(case, "case", "noncase")) %>%
    pivot_wider(names_from = case, values_from = n_group) %>%
    mutate(odds = case / noncase)
  ve = 1 - .tbl$odds[.tbl$vac_record] / .tbl$odds[!.tbl$vac_record]
  if (is.na(ve)) stop("ve calculation failed", call. = FALSE)
  tibble(ve = ve, n_study = n_group)
}

# Summarises a population
# Calculates VE estimates as if a test-negative and a case-control study was
# done on the population
summarise_pop <- function(pop, studies = c("tn", "cc")) {
  map_dfr(studies, ~ calc_counts(pop, .x)) %>%
    group_by(study) %>%
    group_modify(~ calc_ve(.x))
}

# Wraps the simulation and the population summary functions
sim_sum <- function(..., init_seed = NULL, ind = NULL) {
  simseed <- if (!is.null(init_seed)) init_seed + ind else NULL
  one <- summarise_pop(simulate_pop(..., simseed = simseed))
  if (!is.null(simseed)) one$simseed <- simseed
  if (!is.null(ind)) one$simind <- ind
  one
}

# Many parallel simulations for a parameter set
sim_par_set <- function(..., nsim, init_seed = NULL) {
  args <- list(...)
  future_map_dfr(
    1:nsim, 
    function(ind) do.call(sim_sum, c(args, init_seed = init_seed, ind = ind))
  )
}

# Summarise a set of simulations
sum_par_set <- function(set, ...) {
  set_sum <- set %>%
    group_by(study) %>%
    summarise(ve_mean = mean(ve), ve_sd = sd(ve), n_study_mean = mean(n_study))
  pars <- list(...) %>% as_tibble() %>% slice(rep(1, nrow(set_sum)))
  bind_cols(set_sum, pars)
}

# Wraps the simulation and the summary functions
sim_sum_set <- function(..., nsim, name = NULL, init_seed = NULL) {
  res <- sim_par_set(..., nsim = nsim, init_seed = init_seed) %>%
    sum_par_set(...) %>%
    mutate(nsim = nsim)
  if (!is.null(init_seed)) res$init_seed <- init_seed
  if (!is.null(name)) {
    res <- res %>% 
      mutate(name = name) %>% 
      select(name, everything())
  }
  res
}

# Adds a range for parameter combinations (one parameter at a time variation)
add_range <- function(.tbl, var_name, range, combo_name, ...) {
  bind_rows(.tbl, tibble(!!ensym(var_name) := range, name = combo_name, ...))
}

# Reads previous results
read_previous <- function(path) {
  if (!file.exists(path)) return(NULL)
  suppressMessages(read_csv(path))
}

# Gets rid of combos that already were simulated
filter_combos <- function(combos, previous) {
  if (is.null(previous)) return(list(combos = combos, prev = previous))
  combos <- mutate_if(combos, is.numeric, function(vec) round(vec, 2))
  previous <- mutate_if(previous, is.numeric, function(vec) round(vec, 2))
  pars <- names(combos)
  previous <- mutate(previous, id_prev = row_number())
  combos <- mutate(combos, id_comb = row_number())
  common <- inner_join(previous, combos, by = pars) %>% select(id_comb, id_prev)
  list(
    combos = filter(combos, !id_comb %in% common$id_comb) %>% select(-id_comb),
    prev = filter(previous, id_prev %in% common$id_prev) %>% select(-id_prev)
  )
}

# Script ======================================================================

nsim <- 2 # Numer of simulations
res_name <- paste0("result-", nsim, "sims.csv") # Name of csv to save to

# Parameter combinations
combos <- tibble(
  name = "vary-pflu",
  pflu = c(0.005, seq(0.01, 0.1, by = 0.01)), 
  prisk = 1, 
  pnonflu = 0.1,
  sens_vac = 1,
  spec_vac = 1,
  sens_flu = 1,
  spec_flu = 1) %>%
  add_range(
    prisk, seq(0.1, 0.9, by = 0.1), "vary-prisk", 
    pflu = 0.1, 
    pnonflu = 0.1, 
    sens_vac = 1,
    spec_vac = 1,
    sens_flu = 1,
    spec_flu = 1) %>%
  add_range(
    pnonflu, 0.005, "small-sample",
    prisk = 1, 
    pflu = 0.005,
    sens_vac = 1,
    spec_vac = 1,
    sens_flu = 1,
    spec_flu = 1) %>%
  add_range(
    pflu, c(0.005, seq(0.01, 0.1, by = 0.01)), "vary-pflu-miss",
    prisk = 1, 
    pnonflu = 0.1,
    sens_vac = 0.95,
    spec_vac = 0.8,
    sens_flu = 0.86,
    spec_flu = 0.984) %>%
  add_range(
    pnonflu, 0.005, "small-sample-miss", 
    prisk = 1, 
    pflu = 0.005,
    sens_vac = 0.95,
    spec_vac = 0.8,
    sens_flu = 0.86,
    spec_flu = 0.984) %>%
  mutate(init_seed = seq(20191007, by = nsim, length.out = n()))

# Filter out what's been simulated
previous_results <- read_previous(file.path(sim_folder, res_name))
filtered <- filter_combos(combos, previous_results)

# Simualte the rest
res <- pmap_dfr(
  filtered$combos,
  function(name, pflu, prisk, pnonflu, 
           sens_vac, spec_vac, sens_flu, spec_flu, 
           init_seed) sim_sum_set(
    nsam = 1e4,
    pvac = 0.5, pflu = pflu, ve = 0.5, pnonflu = pnonflu, prisk = prisk,
    sens_vac = sens_vac, spec_vac = spec_vac, sens_flu = sens_flu, 
    spec_flu = spec_flu,
    nsim = nsim, name = name, init_seed = init_seed
  )
)

# Bind to previous and save
res_all <- bind_rows(filtered$prev, res)
write_csv(res_all, file.path(sim_folder, paste0("result-", nsim, "sims.csv")))
