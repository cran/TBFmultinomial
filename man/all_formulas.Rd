% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all_formulas.R
\name{all_formulas}
\alias{all_formulas}
\title{Formulas of all the candidate models}
\usage{
all_formulas(fullModel, discreteSurv = TRUE)
}
\arguments{
\item{fullModel}{formula of the model including all potential variables}

\item{discreteSurv}{Boolean variable telling us whether a `simple'
multinomial regression is looked for or if the goal is a discrete
survival-time model for multiple modes of failure is needed.}
}
\value{
character vector with all the formulas; the first one will be the
reference model; the last element will be the full model.
}
\description{
This function retrieves the formulas of all the candidate models if the
reference model is the null / baseline model.
}
\examples{
data("VAP_data")
FULL <- outcome ~ ns(day, df = 4) + male + type + SOFA
models <- TBFmultinomial:::all_formulas(fullModel = FULL,
discreteSurv = TRUE)
# models

}
\author{
Rachel Heyard
}
