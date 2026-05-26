library(climniche)

sim <- simulate_climniche(n = 1000, p = 3, rho = 0.5, shift = 0.5)

fit <- fit_climniche(
  current = sim$current,
  future = sim$future_away,
  occupied = sim$occupied,
  center = sim$center,
  sensitivity = sim$sensitivity,
  scale = FALSE
)

report <- climniche_report(fit, species = "simulated species")
print(report)

if (requireNamespace("ggplot2", quietly = TRUE)) {
  plot_climniche_exposure(fit)
  plot_climniche_class_summary(fit)
  plot_climniche_variable_contribution(fit)
}

# write_climniche_report(report, "climniche_report.md")
