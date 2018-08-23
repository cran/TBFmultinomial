TBF_ingredients <-
function(fullModel = NULL, data, discreteSurv = FALSE, numberCores = 1,
         candidateModels = NULL, package = 'nnet', maxit = 150){
  #' Ingredients to calculate the TBF
  #'
  #'
  #' This function calculates the ingredients needed to
  #' compute the TBFs: like the deviances with their degrees of freedom of the
  #' relevant candidate models.
  #' @param fullModel formula of the model including all potential variables
  #' @param data the data frame with all the information
  #' @param discreteSurv Boolean variable telling us whether a 'simple'
  #' multinomial regression is looked for or if the goal is a discrete
  #' survival-time model for multiple modes of failure is needed.
  #' @param numberCores How many cores should be used in parallel?
  #' @param candidateModels Instead of defining the full model we can also
  #' specify the candidate models whose deviance statistic and d.o.f should be
  #' computed
  #' @param package Which package should be used to fit the models; by default
  #' the \code{nnet} package is used; we could also specify to use the package
  #' 'VGAM'
  #' @param maxit Only needs to be specified with package \code{nnet}: maximal
  #' number of iterations
  #' @return an object of class \code{TBF.ingredients}
  #' @import VGAM splines nnet parallel
  #' @export
  #'
  #' @author Rachel Heyard


  # Define the candidate models if they have not been defined in candidateModes
  if (is.null(candidateModels)){
    candidateModels <- all_formulas(fullModel = fullModel,
                                    discreteSurv = discreteSurv)
  }
  # the null Model:
  model_0 <- if(package == 'nnet'){
    nnet::multinom(formula = candidateModels[1], data = data,
                   maxit = maxit, trace = FALSE)
  } else {
    if (package == 'VGAM'){
      VGAM::vglm(formula = candidateModels[1], data = data,
                 family = multinomial(refLevel = 1))
    }
  }
  # setTxtProgressBar(pb, 1/(length(candidateModels)))

  # list with results (runs in parallel on numberCores cores)
  mod.ds <- mclapply(2:length(candidateModels),function(counter){
    model <- if (package == 'nnet') {
      nnet::multinom(formula = candidateModels[counter], data = data,
                     maxit = maxit, trace = FALSE)
    } else {
        if (package == 'VGAM'){
          VGAM::vglm(formula = candidateModels[counter], data = data,
                     family = multinomial(refLevel = 1))
        }
      }

      # setTxtProgressBar(pb, counter/(length(candidateModels)))
    if (package == 'nnet'){
      c(abs(deviance(model) - deviance(model_0)),
        abs(model_0$edf - model$edf))
    } else {
      if (package == 'VGAM'){
        c(abs(deviance(model) - deviance(model_0)),
          abs(model_0@df.residual - model@df.residual))
      }
    }
  }, mc.cores = numberCores)
  # close(pb)
  obj <- list(model = candidateModels,
              deviance = c(0, sapply(mod.ds, function(x) x[1])), # the deviance of the nullModel is 0
              degreesOFfreedom = c(0, sapply(mod.ds, function(x) x[2])))
  class(obj) <- append(class(obj), 'TBF.ingredients')
  return(obj)

}
