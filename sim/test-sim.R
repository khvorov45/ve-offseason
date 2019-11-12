# Tests of simulation functions. 
#
# Cannot be used by itself, requires functions
# to be in the global environment. RStudio's "Run Tests" will also not work
# because it will not look in the global environment apparently.
#
# Arseniy Khvorov
# Created 2019/10/04
# Last edit 2019/11/12

library(testthat)

test_that("rcat works", {
  prop_true_1 <- 0.2
  prop_true_2 <- 0.3
  prop_true_3 <- 1 - prop_true_1 - prop_true_2
  catvar <- rcat(1e6, c(prop_true_1, prop_true_2))
  tol <- 0.01
  walk2(
    1:3, c(prop_true_1, prop_true_2, prop_true_3),
    ~ expect_equal(sum(catvar == .x) / length(catvar), .y, tol = tol)
  )
})

test_that("simulation function works", {
  nsam <- 1e6
  pvac <- 0.4
  pflu <- 0.2
  ve <- 0.6
  pnonflu <- 0.3
  sens_vac <- 0.8
  spec_vac <- 0.9
  sens_flu <- 0.85
  spec_flu <- 0.95
  pop <- simulate_pop(
    nsam = nsam, pvac = pvac, pflu = pflu, ve = ve, pnonflu = pnonflu,
    sens_vac = sens_vac, spec_vac = spec_vac, sens_flu = sens_flu,
    spec_flu = spec_flu
  )
  
  expect_equal(nrow(pop), nsam)
  
  tol <- 0.01
  
  # Vaccinated proportion
  expect_equal(sum(pop$vaccinated) / nsam, pvac, tol = tol)
  
  # Flu risk (and ve)
  pop_sum <- pop %>%
    mutate(case = status == "Flu") %>%
    count(vaccinated, case, name = "n_cell") %>%
    mutate(vaccinated = if_else(vaccinated, "vac", "unvac")) %>%
    pivot_wider(names_from = vaccinated, values_from = n_cell)
  r_vac <- pop_sum$vac[pop_sum$case] / sum(pop_sum$vac)
  r_unvac <- pop_sum$unvac[pop_sum$case] / sum(pop_sum$unvac)
  expect_equal(r_vac, pflu * (1 - ve), tol = tol)
  expect_equal(r_unvac, pflu, tol = tol)
  
  # Vaccination record
  true_vac <- sum(pop$vac_record)
  exp_rec <- sum(pop$vaccinated) * sens_vac +
    sum(!pop$vaccinated) * (1 - spec_vac)
  expect_equal(exp_rec / nsam, sum(pop$vac_record) / nsam, tol = tol)
  
  # Flu test
  true_flu <- sum(pop$status == "Flu")
  exp_pos <- sum(pop$status == "Flu") * sens_flu + 
    sum(pop$status != "Flu") * (1 - spec_flu)
  expect_equal(exp_pos / nsam, sum(pop$test_result) / nsam, tol = tol)
})

test_that("counting works (no missclassification)", {
  nsam <- 1e6
  pvac <- 0.4
  pflu <- 0.2
  ve <- 0.6
  pnonflu <- 0.3
  pop <- simulate_pop(
    nsam = nsam, pvac = pvac, pflu = pflu, ve = ve, pnonflu = pnonflu
  )
  tn <- calc_counts(pop, "tn")
  tol <- 0.01
  p_tn <- pvac * pflu * (1 - ve) + (1 - pvac) * pflu + pnonflu
  expect_equal(
    tn$n_group[tn$vac_record & tn$case] / sum(tn$n_group),
    pvac * pflu * (1 - ve) / p_tn,
    tol = tol
  )
  expect_equal(
    tn$n_group[tn$vac_record & !tn$case] / sum(tn$n_group), 
    pvac * pnonflu / p_tn,
    tol = tol
  )
  expect_equal(
    tn$n_group[!tn$vac_record & tn$case] / sum(tn$n_group), 
    (1 - pvac) * pflu / p_tn,
    tol = tol
  )
  expect_equal(
    tn$n_group[!tn$vac_record & !tn$case] / sum(tn$n_group), 
    (1 - pvac) * pnonflu / p_tn,
    tol = tol
  )
  cc <- calc_counts(pop, "cc")
  expect_equal(
    cc$n_group[cc$vac_record & cc$case] / nsam,
    pvac * pflu * (1 - ve),
    tol = tol
  )
  expect_equal(
    cc$n_group[cc$vac_record & !cc$case] / nsam, 
    pvac * (1 - pflu * (1 - ve)),
    tol = tol
  )
  expect_equal(
    cc$n_group[!cc$vac_record & cc$case] / nsam, 
    (1 - pvac) * pflu,
    tol = tol
  )
  expect_equal(
    cc$n_group[!cc$vac_record & !cc$case] / nsam, 
    (1 - pvac) * (1 - pflu),
    tol = tol
  )
})

test_that("VE calculations work", {
  ve <- 0.6
  pop <- simulate_pop(nsam = 1e6, ve = ve)
  tn <- calc_counts(pop, "tn")
  ve_tn <- calc_ve(tn)
  tol <- 0.02
  expect_equal(ve_tn$ve, ve, tol = tol)
})

test_that("population summary works", {
  pop <- simulate_pop()
  popsum <- summarise_pop(pop, studies = c("tn", "cc"))
  expect_equal(nrow(popsum), 2)
})

