model_priors <-
function(fullModel, discreteSurv = TRUE, modelPrior = 'flat'){
  #' Prior model probability
  #'
  #'
  #' This function computes the prior model probabilities of the candidate
  #' models
  #' @param fullModel formula of the model including all potential variables
  #' @param discreteSurv Boolean var telling us whether a 'simple' multinomial
  #' regression is looked for or if the goal is a discrete survival-time model
  #' for multiple modes of failure is needed.
  #' @param modelPrior what prior should be used on the model space?
  #' \code{modelPrior} should be included in \{\code{'flat','dependent'}\}
  #' where \code{'flat'} means a uniform pior and \code{'dependent'} sets a
  #' multiplicity-corrected model prior on the model space.
  #' @return a numerical vector with the prior model probabilities
  #' @importFrom utils combn
  #' @export
  #' @examples
  #' # the definition of the full model with three potential predictors:
  #' FULL <- outcome ~ ns(day, df = 4) + gender + type + SOFA
  #' # here we define time as a spline with 3 knots
  #'
  #' priors <- model_priors(fullModel = FULL, discreteSurv = TRUE,
  #'                        modelPrior = 'dependent')
  #' @author Rachel Heyard

  candidateModels <- all_formulas(fullModel = fullModel,
                                  discreteSurv = discreteSurv)
  if (modelPrior == 'flat'){
    prior <- rep(1/(length(candidateModels)), length(candidateModels))
  }
  if (modelPrior == 'dependent'){
    # getting out the variables:
    vars <- str_trim(str_split(as.character(fullModel), '\\+')[[3]][-1])

    p <- length(vars) # Number of variables

    # All possible combinations of the variables:
    comb.list <- lapply(1:p, function(i){combn(vars, i) })

    prior <- lapply(1:(length(comb.list) - 1),
                    function(i){rep(1 / ((length(comb.list) + 1) * ncol(comb.list[[i]])),
                                    ncol(comb.list[[i]])) })
    prior <- c(1/(length(comb.list) + 1), unlist(prior), 1/(length(comb.list) + 1))

  }
  return(prior)
}
