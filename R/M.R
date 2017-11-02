M <-
function (a,b) {
  # Normalizing constante of the IG(a, b)
  #
  #
  # This function computes the TBF as well as g
  # @param a hyperparameter of the IG(., .)
  # @param b hyperparameter of the IG(., .)
  # @importFrom stats integrate
  # @author Rachel Heyard

  integrand <- function(x) {(x ^ (a - 1)) * exp(- x)}

  return((b^a) / (integrate(integrand, lower = 0, upper = b)$value))
}
