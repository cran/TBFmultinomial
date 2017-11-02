postInclusionProb <- function(object) {
  #' Posterior inclusion probability (PIP)
  #'
  #'
  #' This function computes the PIPs of all potential predictors
  #' @param object An object of class \code{PMP}
  #' @return an named vector with all PIPs
  #' @importFrom utils head
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
  #' #computation of the posterior inclusion probabilities:
  #' postInclusionProb(test)
  #'
  #' @author Rachel Heyard
    if (!('PMP' %in% class(object))){
      stop('n needs to be defined')
    }
    posts <- object$posterior
    f <- object$model
    fullModel <- tail(f, 1)
    nullModel <- head(f, 1)

    discreteSurv <- !(str_trim(str_split(nullModel, '\\~')[[1]][2]) == '1')
    vars <- str_trim(str_split(fullModel, '\\+')[[1]][-1])

    inclusionProb <- rep(0, length(vars))  # initialisation of vector with the
                                           # posterior inclusion probabilities
    names(inclusionProb) <- vars
    inclusionProb <- sapply(vars, function(v){
      valold <- 0
      for (i in 1:(length(f))) {
        if (grepl(v, f[i])) {
          valnew <- valold + posts[i]
          valold <- valnew
        }
      }
      return(valnew)
    })

    return(inclusionProb)
  }
