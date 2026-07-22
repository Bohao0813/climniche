library(climniche)

# Fit one current-to-future climate comparison.
sim <- simulate_climniche(n = 1200, p = 6, seed = 31)
fit <- fit_climniche(
  current = sim[["current"]],
  future = sim[["future_away"]],
  occupied = sim[["occupied"]],
  sensitivity = sim[["sensitivity"]]
)

# Identify the largest climatic contribution at each reference cell.
contribution <- climniche_dominant_contribution(fit, scope = "current")
contribution
contribution[["summary"]]

cell_table <- contribution[["table"]]
head(cell_table[cell_table[["included"]], ])
