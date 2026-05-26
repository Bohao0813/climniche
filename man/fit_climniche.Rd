\name{fit_climniche}
\alias{fit_climniche}
\alias{fit_climniche_raster}
\alias{fit_climniche_terra}
\alias{climniche_table}
\alias{climniche_summary}
\alias{climniche_report}
\alias{print.climniche_report}
\alias{write_climniche_report}
\alias{plot_climniche_map}
\alias{plot_climniche_classes}
\alias{plot_climniche_exposure}
\alias{plot_climniche_class_summary}
\alias{plot_climniche_distribution}
\alias{plot_climniche_report}
\alias{climniche_diagram_data}
\alias{plot_climniche_diagram}
\alias{climniche_showcase_data}
\alias{plot_climniche_showcase}
\alias{niche_metric}
\alias{niche_potential}
\alias{niche_radius}
\alias{boundary_exceedance}
\alias{niche_percentile}
\alias{variable_contribution}
\alias{simulate_climniche}
\alias{plot_variable_contribution}
\alias{summary.climniche_fit}
\alias{plot.climniche_fit}
\title{Assess niche climate exposure}
\description{
Core functions for assessing niche climate exposure: projected climate change
interpreted relative to the current realised climatic niche of a species. The
framework separates climate change amount, niche distance change, composition
change and niche boundary exceedance.
}
\usage{
fit_climniche(current, future, occupied = NULL, cnfa = NULL,
  center = NULL, sensitivity = NULL, A = NULL,
  metric = c("diag", "factor"), boundary = 0.95, scale = TRUE,
  global_mean = NULL, global_sd = NULL, conflict_ratio = 0.25)

fit_climniche_raster(current, future, occupied = NULL,
  occupied_threshold = 0, domain = NULL, domain_threshold = 0, ...)

fit_climniche_terra(current, future, occupied = NULL,
  occupied_threshold = 0, domain = NULL, domain_threshold = 0, ...)

climniche_table(x, scope = c("current", "all"))

climniche_summary(x, scope = c("current", "all"))

climniche_report(x, species = NULL, scope = c("current", "all"),
  top_variables = 5)

\method{print}{climniche_report}(x, ...)

write_climniche_report(report, file)

plot_climniche_map(x, metric = c("niche_distance_change",
  "outside_niche_exceedance", "climate_change_amount", "composition_change",
  "change_alignment"), occupied = NULL, occupied_only = FALSE,
  title = NULL, midpoint = 0)

plot_climniche_classes(x, occupied = NULL, occupied_only = FALSE,
  title = NULL)

plot_climniche_exposure(x, scope = c("current", "all"),
  max_points = 6000, seed = 1, title = NULL)

plot_climniche_class_summary(x, scope = c("current", "all"),
  title = NULL)

plot_climniche_distribution(x, metric = c("niche_distance_change",
  "climate_change_amount", "outside_niche_exceedance", "composition_change"),
  scope = c("current", "all"), title = NULL)

plot_climniche_report(x, scope = c("current", "all"))

climniche_diagram_data(x, scope = c("current", "all"),
  max_arrows = 350, seed = 1)

plot_climniche_diagram(x, scope = c("current", "all"),
  type = c("summary", "sample"), max_arrows = 350, seed = 1,
  show_reference = FALSE, show_hulls = TRUE,
  boundary_shape = c("hull", "circle", "none"),
  show_boundary_label = TRUE, show_points = NULL,
  show_endpoints = FALSE, show_center = FALSE, show_variables = FALSE,
  variable_labels = NULL, title = NULL)

climniche_showcase_data(x, scope = c("current", "all"),
  max_points = 6000, seed = 1, plane_bins = 45,
  boundary_probs = seq(0.50, 0.99, 0.01), top_variables = 6)

plot_climniche_showcase(x, scope = c("current", "all"),
  max_points = 6000, seed = 1, plane_bins = 45,
  boundary_probs = seq(0.50, 0.99, 0.01), top_variables = 6,
  variable_labels = NULL, title = NULL)

niche_metric(sensitivity = NULL, cnfa = NULL, type = c("diag", "factor"))

niche_potential(x, center, A)

niche_radius(psi)

boundary_exceedance(psi_future, boundary_value,
  scale = c("radial", "potential"))

niche_percentile(psi_current, psi_future, occupied)

variable_contribution(current, future, center, A)

simulate_climniche(n = 2000, p = 2, seed = 1, rho = 0,
  prevalence = 0.30, shift = 0.40)

plot_variable_contribution(x, occupied_only = TRUE)

\method{summary}{climniche_fit}(object, ...)

\method{plot}{climniche_fit}(x, type = c("distance", "boundary", "amount",
  "classification"), ...)
}
\arguments{
\item{current}{Numeric matrix or data frame of current environmental values.}
\item{future}{Numeric matrix or data frame of future environmental values, with
the same dimensions and variables as \code{current}.}
\item{occupied}{NULL, logical vector, row indices, or a one layer raster
identifying current occurrence, range, or thresholded SDM cells. If NULL, all
analysed cells are used as the reference set.}
\item{cnfa}{Optional CENFA model object or compatible list.}
\item{center}{Optional realised niche centre in standardised climate space.}
\item{sensitivity}{Optional environmental sensitivity weights.}
\item{A}{Optional niche metric matrix. It must be finite, symmetric and positive
semidefinite.}
\item{metric}{Metric construction when \code{A} is not supplied, or metric to
plot for map and distribution functions.}
\item{boundary}{Quantile defining the empirical realised niche boundary.}
\item{scale}{Logical standardisation flag for \code{fit_climniche()}, or output
scale for \code{boundary_exceedance()}.}
\item{global_mean}{Optional means used for standardisation.}
\item{global_sd}{Optional standard deviations used for standardisation.}
\item{conflict_ratio}{Minimum minority sign contribution share used to mark
mixed variable response as a diagnostic flag. Set to NULL to disable this flag.}
\item{occupied_threshold}{Threshold used for raster \code{occupied} layers.
Continuous suitability rasters are thresholded; values greater than this
threshold are treated as current occurrence cells.}
\item{domain}{Optional one layer raster limiting the cells where exposure is
analysed, for example a marine region, accessible area, management unit or
suitable area predicted by an SDM.}
\item{domain_threshold}{Threshold used when \code{domain} is supplied.}
\item{...}{Additional arguments passed to lower level methods.}
\item{x}{A fitted climniche object, raster object, or numeric matrix depending
on the function.}
\item{object}{A fitted climniche object.}
\item{report}{An object returned by \code{climniche_report()}.}
\item{file}{Output Markdown file.}
\item{species}{Optional species name used in reports.}
\item{top_variables}{Number of variable contributions to show in reports.}
\item{scope}{\code{"current"} for current occurrence, range or thresholded SDM
cells; \code{"all"} for all evaluated cells.}
\item{title}{Optional plot title. Use \code{FALSE} to suppress titles in map
plots.}
\item{midpoint}{Midpoint for the niche distance change colour scale.}
\item{occupied_only}{Logical; summarize or plot current occurrence cells only.}
\item{max_points}{Maximum number of cells to keep as source data for exposure
plane visualisation.}
\item{max_arrows}{Maximum number of current to future arrows to keep in the
niche climate exposure diagram.}
\item{seed}{Random seed.}
\item{show_reference}{Logical; draw the full analysed environmental domain in
the diagram.}
\item{show_hulls}{Logical; draw current and future niche hulls in the diagram.}
\item{boundary_shape}{Boundary display in the diagram. \code{"hull"} draws the
empirical occupied niche polygon, \code{"circle"} draws the constant niche
distance boundary, and \code{"none"} suppresses the boundary.}
\item{show_boundary_label}{Logical; add reference-area explanations below the
exposure-class legend.}
\item{show_points}{Logical; draw future points in the diagram. If NULL, points
are shown only when \code{type = "sample"}.}
\item{show_endpoints}{Logical; draw class mean future positions in the
diagram.}
\item{show_center}{Logical; mark the realised niche centre in the diagram.}
\item{show_variables}{Logical; draw environmental variable directions in the
diagram.}
\item{variable_labels}{Optional named vector replacing variable labels in the
diagram.}
\item{plane_bins}{Number of fixed bins used to summarize the exposure plane.}
\item{boundary_probs}{Boundary quantiles retained for optional sensitivity
diagnostics.}
\item{type}{Metric type or base plot type.}
\item{psi}{Numeric niche potential values.}
\item{psi_current}{Current niche potential.}
\item{psi_future}{Future niche potential.}
\item{boundary_value}{Empirical boundary value in niche potential units.}
\item{p}{Number of environmental variables in simulations.}
\item{n}{Number of simulated cells.}
\item{rho}{Pairwise correlation among simulated environmental variables.}
\item{prevalence}{Proportion of background cells treated as true current
occurrence locations under the virtual niche.}
\item{shift}{Climate change amount imposed in simulated future scenarios.}
}
\value{
\code{fit_climniche()} returns an object of class \code{climniche_fit}. Raster
workflows add a \code{rasters} element with mapped outputs. Report functions
return tables, summaries, Markdown files or an object of class
\code{climniche_report}. Plot functions return ggplot objects, or a patchwork
object when a multi panel report is requested and \pkg{patchwork} is installed.
}
\details{
Niche climate exposure is an exposure interpretation. It does not estimate
population growth, dispersal, adaptation, biotic interactions or conservation
priority by itself. Continuous SDM or suitability rasters are accepted by
thresholding; they are not used as continuous weights unless the user converts
the workflow before fitting.
}
\examples{
sim <- simulate_climniche(n = 200)
fit <- fit_climniche(sim$current, sim$future_away,
  occupied = sim$occupied, center = sim$center,
  sensitivity = sim$sensitivity)
climniche_summary(fit)

report <- climniche_report(fit, species = "simulated species")
print(report)
}
