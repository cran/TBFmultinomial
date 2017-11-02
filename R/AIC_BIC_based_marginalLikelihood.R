AIC_BIC_based_marginalLikelihood <-
function(fullModel = NULL, candidateModels = NULL,
         data, discreteSurv = TRUE, AIC = TRUE,
         package = 'nnet', maxit = 150, numberCores = 1){
  #' Marginal likelihoods based on AIC or BIC
  #'
  #'
  #' This function computes the marginal likelihoods based on the AIC or on the
  #' BIC, that will later be used to calculate the TBF.
  #' @param fullModel formula of the model including all potential variables
  #' @param candidateModels Instead of defining the full model we can also
  #' specify the candidate models whose deviance statistic and d.o.f should be
  #' computed
  #' @param data the data
  #' @param discreteSurv Boolean variable telling us whether a `simple'
  #' multinomial regression is looked for or if the goal is a discrete
  #' survival-time model for multiple modes of failure is needed.
  #' @param AIC if \code{TRUE}, AIC will be used, else we use BIC
  #' @param package Which package should be used to fit the models; by default
  #' the \code{nnet} package is used; we could also specify to use the package
  #' 'VGAM'
  #' @param maxit Only needs to be specified with package \code{nnet}: maximal
  #' number of iterations
  #' @param numberCores How many cores should be used in parallel?
  #' @return a vector with the marginal likelihoods of all candidate models
  #' @export
  #' @examples
  #' # data extraction:
  #' data("VAP_data")
  #'
  #' # the definition of the full model with three potential predictors:
  #' FULL <- outcome ~ ns(day, df = 4) + gender + type + SOFA
  #' # here the define time as a spline with 3 knots
  #'
  #' # now we can compute the marginal likelihoods based on the AIC f.ex:
  #' mL_AIC <-
  #' AIC_BIC_based_marginalLikelihood(fullModel = FULL,
  #'                                  data = VAP_data,
  #'                                  discreteSurv = TRUE,
  #'                                  AIC = TRUE)
  #' @author Rachel Heyard

  # Define the candidate models:
  if (is.null(candidateModels) & !is.null(fullModel)){
    candidateModels <- all_formulas(fullModel = fullModel,
                                    discreteSurv = discreteSurv)
  } else {
    if (!is.null(candidateModels) & !is.null(fullModel)){
      if (tail(candidateModels, 1) == fullModel){
        candidateModels <- all_formulas(fullModel = fullModel,
                                        discreteSurv = discreteSurv)
      } else stop('Either define fullModel OR candidatModels')
    }
  }

  # list with results (runs in parallel on numberCores cores)
  mod.IC <- mclapply(1:length(candidateModels), function(counter){
    m <- candidateModels[counter]
    model <- if(package == 'VGAM') {
      vglm(formula = m, data = data, family = multinomial(refLevel = 1))
    } else {
      nnet::multinom(formula = m, data = data, maxit = maxit, trace = FALSE)
    }
    if (AIC){
      IC <- AIC(model)
    } else {
      IC <- BIC(model)
    }
    IC
  }, mc.cores = numberCores)

  # The information criterions:
  IC <- unlist(lapply(mod.IC, function(x) x[1]))

  # increment to make sure that the exponential does not deviate
  incr <- mean(IC/2)
  # The marginal likelihoods:
  marginal.likelihoods <- sapply(IC, function(x) exp(- x / 2  + incr) )
  return(marginal.likelihoods)
}
