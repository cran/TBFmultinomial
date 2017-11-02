PMP <-
function(fullModel= NULL, candidateModels = NULL, data = NULL,
         discreteSurv = TRUE, method = 'LEB', prior = 'flat',
         package = 'nnet', maxit = 150, numberCores = 1) {
  #' Posterior model probability
  #'
  #'
  #' This function computes the posterior probability of all candidate models
  #' @param fullModel formula of the model including all potential variables
  #' @param candidateModels Instead of defining the full model we can also
  #' specify the candidate models whose deviance statistic and d.o.f should be
  #' computed
  #' @param data the data frame with all the information
  #' @param discreteSurv Boolean variable telling us whether a 'simple'
  #' multinomial regression is looked for or if the goal is a discrete
  #' survival-time model for multiple modes of failure is needed.
  #' @param method tells us which method for the definition of g should be
  #' used. Possibilities are: \code{LEB}, \code{GEB}, \code{g=n}, \code{hyperG},
  #' \code{ZS}, \code{ZSadapted} and \code{hyperGN}
  #' @param prior should a dependent or a flat prior be used on the model space?
  #'  Only needed if \code{method = `GEB`}.
  #' @param numberCores How many cores should be used in parallel?
  #' @param package Which package should be used to fit the models; by default
  #' the \code{nnet} package is used; we could also specify to use the package
  #' 'VGAM'
  #' @param maxit Only needs to be specified with package \code{nnet}: maximal
  #' number of iterations
  #' @return an object of class \code{TBF.ingredients}
  #' @export
  #' @examples
  #' # extract the data:
  #' data("VAP_data")
  #'
  #' # the definition of the full model with three potential predictors:
  #' FULL <- outcome ~ ns(day, df = 4) + gender + type + SOFA
  #' # here we define time as a spline with 3 knots
  #'
  #' # computation of the posterior model probabilities:
  #' test <- PMP(fullModel = FULL, data = VAP_data,
  #'             discreteSurv = TRUE, maxit = 150)
  #' class(test)
  #'
  #' @author Rachel Heyard

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

  mod.prior <- model_priors(fullModel = fullModel, discreteSurv = discreteSurv,
                            modelPrior = prior)
  if (method %in% c('LEB', 'hyperG', 'GEB', 'ZSadapted', 'g=n', 'hyperGN',
                    'ZS')){
    ingredients.for.TBF <- TBF_ingredients(fullModel = fullModel, data= data,
                                           discreteSurv = discreteSurv,
                                           numberCores = numberCores,
                                           candidateModels = candidateModels,
                                           package = package, maxit = maxit)
    tbf_results <- TBF(ingredients = ingredients.for.TBF,
                       method = method, prior = prior,
                       discreteSurv = discreteSurv)
    TBF <- tbf_results$TBF
    g <- tbf_results$G
    post.mod.prob <- (TBF * mod.prior) / sum(TBF * mod.prior)
  }
  if (method %in% c('AIC', 'BIC')){
    AIC = (method == 'AIC')
    margLike <- AIC_BIC_based_marginalLikelihood(fullModel = fullModel,
                                                 data = data,
                                                 discreteSurv = discreteSurv,
                                                 AIC = AIC)
    post.mod.prob <- (margLike * mod.prior) / sum(margLike * mod.prior)
    g <- NULL
  }
  res <- list(posterior = post.mod.prob, prior = mod.prior,
              model = candidateModels, method = method,
              G = g)
  class(res) <- c('PMP', 'list')
  return(res)
}
