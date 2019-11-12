# Graphs of simulation results
# Arseniy Khvorov
# Created 2019/10/24
# Last edit 2019/11/12

library(readr)
library(dplyr)
library(rlang)
library(ggplot2)
library(scales)
library(ggdark) # devtools::install_github("khvorov45", "ggdark")
library(shadowtext)

graph_folder <- "graph"
sim_folder <- "sim"
sim_filename <- "result-1e+05sims.csv"

# Functions ===================================================================

expected_vac_cases <- function(nsam, 
                               prisk, pvac, pflu, ve, pnonflu,
                               sens_vac, spec_vac,
                               sens_flu, spec_flu) {
  vac_flu_true <- nsam * prisk * pvac * pflu * (1 - ve)
  vac_nonflu_true <- nsam * pvac * pnonflu
  unvac_flu_true <- nsam * prisk * (1 - pvac) * pflu
  unvac_nonflu_true <- nsam * (1 - pvac) * pnonflu
  vac_flu_true * sens_vac * sens_flu + 
    unvac_flu_true * (1 - spec_vac) * sens_flu + 
    vac_nonflu_true * sens_vac * (1 - spec_flu) +
    unvac_nonflu_true * (1 - spec_vac) * (1 - spec_flu)
}

ve_lines <- function(dat, combo_name, var_name,
                     xlab = var_name, 
                     xbreaks = waiver(), 
                     errorwidth = 0.005,
                     ylim = c(0.26, 0.7),
                     textsize = 3) {
  dat %>%
    filter(name == combo_name) %>%
    ggplot(aes(!!sym(var_name), ve_mean)) + 
    dark_theme_bw(verbose = FALSE) +
    theme(panel.grid.minor = element_blank()) +
    scale_x_continuous(xlab, breaks = xbreaks) +
    scale_y_continuous(
      "Estimated vaccine effectiveness",
      breaks = seq(0, 1, 0.05), labels = percent_format(accuracy = 1),
      expand = expand_scale(0)
    ) +
    coord_cartesian(ylim = ylim) +
    geom_hline(aes(yintercept = unique(ve)), lty = "33") +
    geom_errorbar(
      aes(ymin = ve_mean - ve_sd, ymax = ve_mean + ve_sd),
      width = errorwidth, size = 0.025
    ) +
    geom_line() +
    geom_point(shape = 15, size = 3) +
    geom_shadowtext(
      aes(label = signif(expected_vac_cases, 2), y = ve + 0.02),
      col = "black", bg.color = "white",
      size = textsize
    )
}

save_plot <- function(plot, name, folder) {
  ggsave_dark(
    file.path(folder, paste0(name, ".pdf")), plot, dark = FALSE,
    device = "pdf", width = 14, height = 8, units = "cm"
  )
}

# Script ======================================================================

sim_results <- read_csv(file.path(sim_folder, sim_filename))

res_tn <- sim_results %>%
  filter(study == "Test negative") %>%
  mutate(
    expected_vac_cases = expected_vac_cases(
      nsam, prisk, pvac, pflu, ve, pnonflu, 
      sens_vac, spec_vac, sens_flu, spec_flu
    )
  )

plot_pflu <- ve_lines(
  res_tn, "vary-pflu", "pflu", "Probability of infection by influenza",
  errorwidth = 0.003, xbreaks = seq(0, 0.1, 0.01)
)
plot_pflu
save_plot(plot_pflu, "pflu", graph_folder)

plot_prisk <- ve_lines(
  res_tn, "vary-prisk", "prisk", "Proportion at risk",
  errorwidth = 0.03, xbreaks = seq(0, 1, 0.1)
)
plot_prisk
save_plot(plot_prisk, "prisk", graph_folder)

plot_pflu_miss <- ve_lines(
  res_tn, "vary-pflu-miss", "pflu", "Probability of infection by influenza",
  errorwidth = 0.003, xbreaks = seq(0, 0.1, 0.01), ylim = c(0, 0.55),
  textsize = 3
)
plot_pflu_miss
save_plot(plot_pflu_miss, "pflu-miss", graph_folder)
