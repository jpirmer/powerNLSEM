% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.powerNLSEM.R
\name{plot.powerNLSEM}
\alias{plot.powerNLSEM}
\title{plot powerNLSEM object}
\usage{
\method{plot}{powerNLSEM}(
  x,
  test = NULL,
  plot = "power_model",
  power_modeling_method = NULL,
  se = FALSE,
  power_aim = NULL,
  alpha = NULL,
  alpha_power_modeling = NULL,
  min_num_bins = 10,
  defaultgg = FALSE,
  ...
)
}
\arguments{
\item{x}{object of class powerNLSEM}

\item{test}{Should the parameter be tested with a directed hypothesis (onesided) or with an undirected hypothesis (twosided, also equivalent to Wald-Test for single parameter). Default to \code{NULL}, then the same as in fitted \code{powerNLSEM} object in \code{x} is used.}

\item{plot}{Character indicating what type of plot to create. Default to \code{"power_model"}, referencing to the prediction of significant parameters using the model specified in \code{power_modeling_method}.}

\item{power_modeling_method}{Character indicating the power modeling method used. This is only relevant when \code{plot = "power_model"} is used. Default to \code{NULL}, indicating to use the same power modeling method as was used in the \code{powerNLSEM} function.}

\item{se}{Logical indicating to use confidence intervals based on normal approximation using the standard errors. Default to \code{FALSE}.}

\item{power_aim}{Power level to be included into the plot with respective N. If \code{NULL} the same power level as in the \code{powerNLSEM} function will be used. If set to \code{0} no power level and corresponding N will be plotted. Default to \code{NULL}, indicating to use the same power modeling method as was used in the \code{powerNLSEM} function.}

\item{alpha}{Alpha value used for confidence intervals, when \code{se = TRUE}. Default to \code{NULL}, indicating to use the same alpha as was used in the powerNLSEM function. This does not influence the significance decision, although same alpha is used per default.}

\item{alpha_power_modeling}{Type I-error rate for confidence band around predicted power rate. Used to ensure that the computed \code{N} keeps the desired power value (with the given Type I-error rate \code{alpha_power_modeling} divided by 2). If set to 1, no confidence band is used. Default to \code{.05}.}

\item{min_num_bins}{minimal number of bins used for aggregating results. Default to 10.}

\item{defaultgg}{Logical to return default ggplot object. Default to \code{FALSE}, which returns \code{theme_minimal} and other changes in theme.}

\item{...}{Additional arguments passed on to the plot function.}
}
\value{
Returns \code{ggplot} object of the type specified in plot.
}
\description{
plot powerNLSEM object
}
\examples{
\donttest{
# write model in lavaan syntax
model <- "
# measurement models
          X =~ 1*x1 + 0.8*x2 + 0.7*x3
          Y =~ 1*y1 + 0.85*y2 + 0.78*y3
          Z =~ 1*z1 + 0.9*z2 + 0.6*z3

# structural models
          Y ~ 0.3*X + .2*Z +  .2*X:Z

# residual variances
         Y~~.7975*Y
         X~~1*X
         Z~~1*Z

# covariances
         X~~0.5*Z

# measurement error variances
         x1~~.1*x1
         x2~~.2*x2
         x3~~.3*x3
         z1~~.2*z1
         z2~~.3*z2
         z3~~.4*z3
         y1~~.5*y1
         y2~~.4*y2
         y3~~.3*y3
"
# run model-implied simulation-based power estimation
# for the effects: c("Y~X", "Y~Z", "Y~X:Z")
Result_Power <- powerNLSEM(model = model, POI = c("Y~X", "Y~Z", "Y~X:Z"),
                           method = "UPI", search_method = "adaptive",
                           steps = 10, power_modeling_method = "probit",
                           R = 1000, power_aim = .8, alpha = .05,
                           alpha_power_modeling = .05,
                           CORES = 1, seed = 2024)

Result_Power
plot(Result_Power)
}
}
