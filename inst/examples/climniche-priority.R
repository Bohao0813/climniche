library(climniche)

# Fit one projected climate comparison.
sim <- simulate_climniche(n = 1200, p = 6, seed = 18)
fit <- fit_climniche(
  current = sim[["current"]],
  future = sim[["future_away"]],
  occupied = sim[["occupied"]],
  sensitivity = sim[["sensitivity"]]
)

# Rank current reference cells by reference weight and positive niche shift.
priority <- climniche_priority(fit)
priority

priority_table <- priority[["table"]]
priority_cells <- priority_table[priority_table[["included"]], ]
head(priority_cells[order(priority_cells[["pareto_rank"]]), ])

if (requireNamespace("ggplot2", quietly = TRUE)) {
  plot_climniche_priority(priority, type = "plane")
}
