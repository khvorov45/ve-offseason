# Investigation of the behaviour of the expected value of the 
# reciprocal of a binomial random variable.
#
# Arseniy Khvorov
# Created 2019/10/09
# Last edit 2019/11/12

library(tidyverse)
library(scales)
library(ggdark)

binrec_folder <- "binomial-reciprocal"

# Functions ===================================================================

# Generate binomial samples for one combination
gen_bin <- function(pop_size, exp_vac_case, nsim) {
  x <- rbinom(nsim, size = pop_size, prob = exp_vac_case / pop_size)
  if (any(x == 0)) stop("0 x produced", call. = FALSE)
  tibble(
    nsim = nsim, npop = pop_size, 
    expected = exp_vac_case, mean = mean(x),
    mean_rec = 1 / mean(x), rec_mean = mean(1 / x)
  )
}

# Generate samples for many combinations
gen_combos <- function(pars, nsim) {
  map2_dfr(pars$pop_size, pars$exp_vac_case, gen_bin, nsim = nsim)
}

to_long <- function(data, exp_vac_case) {
  data %>%
    pivot_longer(
      c(mean_rec, rec_mean), names_to = "mean_type", values_to = "value",
    ) %>%
    mutate(bias = value - 1 / expected)
}

plot_lines <- function(data, xvar, xlab, ylim, add_label) {
  xvar <- ensym(xvar)
  xbreaks <- unique(pull(data, !!xvar))
  sum_bias <- data %>%
    group_by(mean_type) %>%
    summarise(mean_bias = mean(bias), min_x = min(!!xvar))
  data %>%
    mutate(
      mean_type = recode(
        mean_type,
        "mean_rec" = "Reciprocal of mean", "rec_mean" = "Mean of reciprocal"
      )
    ) %>%
    ggplot(aes(!!xvar, bias, shape = mean_type, lty = mean_type)) +
    dark_theme_bw(verbose = FALSE) +
    theme(
      legend.position = "bottom",
      legend.margin = margin(0, 0, 0, 0, "lines"),
      legend.title = element_blank(),
      legend.box.spacing = unit(0, "lines"),
      panel.grid.minor = element_blank(),
      plot.margin = unit(c(0.5, 1, 0.5, 0.5), "lines")
    ) +
    scale_x_log10(xlab, breaks = xbreaks, labels = comma_format()) +
    scale_y_continuous("Estimated bias", labels = percent_format()) +
    scale_shape_manual(values = c(15, 18)) +
    coord_cartesian(ylim = ylim) +
    geom_line() +
    geom_point(size = 4)
}

# Script ======================================================================

nsim <- 1e6

pars <- list(
  "npsame" = list(
    pop_size = c(200, 50 * 10^(1:5)),
    exp_vac_case = 25
  ),
  "npdiff" = list(
    pop_size = 500,
    exp_vac_case = 15 * 2^(0:5)
  )
)

res <- map(pars, gen_combos, nsim = nsim)

res_long <- map(res, to_long)

plots <- pmap(
  list(
    res_long, 
    c("npop", "expected"), 
    c("Number of trials", "Expected count"),
    list(c(0, 0.0019), c(0, 0.0051)),
    c(2e-4, 7e-4)
  ), 
  plot_lines
)
plots[[1]]
plots[[2]]

walk2(
  file.path(binrec_folder, paste0(c("npop", "np"), ".pdf")),
  plots,
  ggsave_dark,
  dark = FALSE,
  width = 10, height = 5, units = "cm", device = "pdf"
)
