GEBargMax <-
  function (g, z, d, modelPriors, incr) {
    return(sum((g + 1)^(- d / 2) * exp((g / (g + 1))* z / 2 - incr) * modelPriors))
  }

search_GEB <-
  function(fullModel, data, discreteSurv, modelPrior = 'flat', incr){

    ingred <- TBF_ingredients(fullModel = fullModel, data = data,
                              discreteSurv = discreteSurv)
    deviances <- ingred$deviance
    dofs <- ingred$degreesOFfreedom

    mod.prior <- model_priors(fullModel = fullModel, discreteSurv = discreteSurv,
                              modelPrior = modelPrior)

    G <- optimise(function(g) {GEBargMax(g, z = deviances, d = dofs,
                                         modelPriors = mod.prior, incr)},
                  c(0:1000), maximum = TRUE)$maximum

    return(list(GEB = G, TBFingredients = ingred))
  }

helper_forTBF <- function(candMod, deviances, dofs, method,
                          prior = NULL, n = NULL, discreteSurv = NULL){
  # Initialisation of the TBF vector:
  TBF <- rep(1, length(dofs))
  g <- rep(0, length(dofs))

  # Some methods for g estimation need n:

  # In order to make sure that the exponential of the deviance does not
  # deviate we need to substract from the deviance, f.ex. this increment
  increment <-  max(deviances - dofs)/3
  TBF[1] <- exp(- increment) # the TBF of the null / reference model should be
  # 1, since we substract the increment in the
  # exponential we need to normalise the TBF of the
  # reference model!
  if (method == 'g=n'){
    g[1] <- n
  }
  if (method== 'LEB'){
    g[1] <- 0
  }


  if (method == 'GEB'){
    fullModel <- utils::tail(candMod, 1)
    mod.prior <- model_priors(fullModel = fullModel,
                              discreteSurv = discreteSurv,
                              modelPrior = prior)
    Ggeb <- optimise(function(g) {
      GEBargMax(g, z = deviances, d = dofs,
                modelPriors = mod.prior,
                incr = increment)},
      c(0:1000), maximum = TRUE)$maximum
    g[1] <- Ggeb
  }
  for (i in 2:length(dofs)){
    if (method == 'LEB') {
      TBF[i] <- max( (exp( (deviances[i] - dofs[i]) / 2 - increment )) *
                       (deviances[i] / dofs[i])^(- dofs[i] / 2), exp( - increment) )
      g[i] <- max(deviances[i] / dofs[i] - 1, 0)
    }
    if (method == 'hyperG') {
      TBF[i] <- (1 / M(1 + dofs[i] / 2,
                       deviances[i] / 2)) * exp(deviances[i] / 2 - increment)
      g <- NULL
    }
    if (method == 'ZSadapted') {
      TBF[i] <- (M(0.5, (n + 3) / 2) / M(0.5 + dofs[i] / 2,
                                         ((n + 3) + deviances[i]) / 2)) *
        exp(deviances[i]/2 - increment)
      g <- NULL
    }
    if (method == 'g=n') {
      TBF[i] <- ((n + 1)^(- dofs[i] / 2)) * exp((n * deviances[i]) /
                                                  (2 * (n + 1)) - increment)
      g[i] <- n
    }
    if(method == 'GEB'){
      TBF[i] <- (Ggeb + 1)^(- dofs[i] / 2) *
        exp((Ggeb / (Ggeb + 1)) * deviances[i] / 2 - increment)
      g[i] <- Ggeb
    }
    if (method %in% c('hyperGN', 'ZS')){
      # Under the null model:
      margLik_0 <- exp(- deviances[i] / 2 + increment)
      # the approximate marginal likelihoods
      if (method == 'ZS'){
        integrand <- function(g) {
          (g + 1)^(-dofs[i] / 2) * g^(-3 / 2) *
            exp(- (deviances[i]/(g+1) + n/g) / 2)
        }
        margLik_j <- (((n / 2)^(1 / 2)) / gamma(1 / 2)) *
          (integrate(integrand, lower = 0, upper = Inf)$value)
      }
      if (method == 'hyperGN'){
        integrand <- function(g) {
          (g + 1)^(-dofs[i] / 2) * exp(- deviances[i] / (2 * (g + 1))) * (n/((n + g)^2))
        }
        margLik_j <- (integrate(integrand, lower = 0, upper = Inf)$value)
      }
      TBF[i] <- margLik_j / margLik_0
      g <- NULL
    }
  }
  return(list(TBF = TBF, g = g))
}

TBF <- function(ingredients = NULL, fullModel = NULL,
                method = 'LEB', data = NULL, discreteSurv = TRUE,
                prior = NULL, package = 'nnet',
                maxit = 150) {
  #' Test-based Bayes factor
  #'
  #'
  #' This function computes the TBF as well as g
  #' @param ingredients \code{TBF_ingredients_object} ingredients for the TBF
  #' (and g) calculation.
  #' @param fullModel if \code{ingredients} is \code{NULL}, formula of the model
  #'  including all potential variables
  #' @param method tells us which method for the definition of g should be
  #' used. Possibilities are: \code{LEB}, \code{GEB}, \code{g=n}, \code{hyperG},
  #' \code{ZS}, \code{ZSadapted} and \code{hyperGN}
  #' @param data the data frame with all the information. Only needed if
  #' \code{ingredients} is \code{NULL}
  #' @param discreteSurv Boolean variable telling us whether a 'simple'
  #' multinomial regression is looked for or if the goal is a discrete
  #' survival-time model for multiple modes of failure is needed.
  #' @param prior should a dependent or a flat prior be used on the model space?
  #'  Only needed if \code{method = `GEB`}.
  #' @param package Which package should be used to fit the models; by default
  #' the \code{nnet} package is used; we could also specify to use the package
  #' 'VGAM'
  #' @param maxit Only needs to be specified with package \code{nnet}: maximal
  #' number of iterations
  #' @return A list with the TBF and the g (if it is fixed) for all the
  #' candidate models.
  #' @importFrom stats optimise integrate
  #' @export
  #'
  #' @author Rachel Heyard


  if ('TBF.ingredients' %in% class(ingredients)){
    # Get all the ingredients:
    candidateModels <- ingredients$model
    deviances <- ingredients$deviance
    dofs <- ingredients$degreesOFfreedom

    out <- stringr::str_trim(stringr::str_split(candidateModels[1],
                                                '~')[[1]][1])
    n <- sum(table(data[, out])[-1])
    tbf_results <- helper_forTBF(candMod = candidateModels,
                                 deviances = deviances, dofs = dofs,
                                 method = method, prior = prior, n = n,
                                 discreteSurv = discreteSurv)
    TBF <- tbf_results$TBF
    g <- tbf_results$g
  }
  else {
    # Get all the ingredients:
    ingredients.for.TBF <- TBF_ingredients(fullModel = fullModel, data = data,
                                           discreteSurv = discreteSurv,
                                           package = package, maxit = maxit)

    candidateModels <- ingredients.for.TBF$model
    deviances <- ingredients.for.TBF$deviance
    dofs <- ingredients.for.TBF$degreesOFfreedom

    out <- stringr::str_trim(stringr::str_split(candidateModels[1],
                                                '~')[[1]][1])
    n <- sum(table(data[, out])[-1])

    tbf_results <- helper_forTBF(candMod = candidateModels,
                                 deviances = deviances, n = n, dofs = dofs,
                                 method = method, prior = prior,
                                 discreteSurv = discreteSurv)
    TBF <- tbf_results$TBF
    g <- tbf_results$g

  }
  TBF <- list(TBF = TBF, deviance = deviances, dof = dofs,
              method = method, G = g, candidateModels = candidateModels)
  class(TBF) <- append(class(TBF), 'TBF')
  return(TBF)
}
