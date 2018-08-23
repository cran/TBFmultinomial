sample_multinomial <- function (PMP_object,
                                shrink = TRUE, data, which = "MPM",
                                discreteSurv = TRUE) {
  #' Samples from a PMP object
  #'
  #'
  #' This function samples from a specific model inside a \code{PMP} object.
  #' @param PMP_object formula of the model including all potential variables
  #' @param shrink should the coefficients be shrunken towards their prior mean?
  #' @param data the (training) data frame with all the information
  #' @param which which model should be sampled from? either an integer, 'MPM'
  #' or 'MAP'
  #' @param discreteSurv Boolean variable telling us whether a 'simple'
  #' multinomial regression is looked for or if the goal is a discrete
  #' survival-time model for multiple modes of failure is needed.
  #' @return returns an object with the model coefficients and supplementary
  #' information
  #' @export
  #'
  #' @author Rachel Heyard
  modelFormulas <- as.character(PMP_object$model)
  fullModel <- tail(modelFormulas, 1)
  PMP <- PMP_object$posterior
  g <- PMP_object$G
  if (length(g) == 1) {
    g <- rep(g, length(PMP))
  }
  PIP <- postInclusionProb(PMP_object)
  outcome <- str_trim(str_split(as.character(fullModel), "\\~")[[1]][1])
  t <- g/(g + 1)
  if (is.numeric(which)) {
    Mindex <- which
    Mformula <- modelFormulas[Mindex]
    Mmodel <- vglm(formula = Mformula, data = data, family = multinomial(refLevel = 1))
  }
  else {
    if (which == "MAP") {
      Mindex <- which(PMP == max(PMP))
      Mformula <- modelFormulas[Mindex]
      Mmodel <- vglm(formula = Mformula, data = data, family = multinomial(refLevel = 1))
    }
    else {
      MPMvariable <- names(PIP[which(PIP$inclusionProb >= 0.5)])
      Mformula <- paste(outcome, "~", intercept, "+", paste(MPMvariable,
                                                            collapse = "+"), sep = "")
      Mindex <- which(modelFormulas == Mformula)
      Mmodel <- vglm(formula = Mformula, data = data, family = multinomial(refLevel = 1))
    }
  }
  if (discreteSurv) {
    vars <- str_split(str_split(as.character(Mformula), "\\~")[[1]][-1],
                      "\\+")[[1]][-1]
    intercept <- str_split(str_split(as.character(Mformula),
                                     "\\~")[[1]][-1], "\\+")[[1]][1]
  }
  else {
    vars <- strsplit(strsplit(as.character(Mformula), "\\~")[[1]][-1],
                     "\\+")[[1]]
    intercept <- NA
  }

  outcomeSpecific <-
    lapply(1:(length(levels(data[, outcome])) - 1),
           function(i) {which(grepl(paste(":", i, sep = ""),
                                    names(coef(Mmodel))))
           })
  if (length(vars)>0){
    variableIndex <- lapply(1:length(vars), function(i) {
      which(grepl(vars[i], names(coef(Mmodel))))
    })
    names(variableIndex) <- vars
    coefs <- lapply(outcomeSpecific, function(x) {
      coefs <- coef(Mmodel)[x]
      if (shrink){
        toShrink <- which(x %in% unlist(variableIndex))
        coefs[toShrink] <- t[Mindex] * coefs[toShrink]
      }
      return(coefs)
    })
  } else coefs <- lapply(outcomeSpecific, function(x) coef(Mmodel)[x])

  X <- model.matrix(as.formula(Mformula), data)
  predictions <- lapply(coefs, function(x) {
    X %*% x
  })
  names(predictions) <- names(coefs) <- levels(data[, outcome])[-1]
  return(list(Coefficients = coefs, Predictions = predictions,
              shrinkageCoef = t[Mindex], which = which, inclusionProbs = PIP))
}
