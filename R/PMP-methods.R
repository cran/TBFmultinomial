
##' Convert a PMP object into a data frame
##'
##' This function takes a PMP object an returns a \code{data.frame} summarising
##' the information.
##' @method as.data.frame PMP
##' @param x valid \code{\link{PMP}} object
##' @param ... arguments to be passed to \code{data.frame}
##' @rdname as.data.frame
##' @import methods
##' @importFrom utils tail
##' @return a \code{data.frame} with the posterior and prior probabilities as
##' well as the definition of the models
##' @export
##' @author Rachel Heyard

as.data.frame.PMP <- function(x, ...) {
  if (!('PMP' %in% class(x))){
    stop('The function needs an object of class PMP')
  }
  ## posterior probabilites:
  ret <- data.frame(posterior= x$posterior, ...)
  ## log prior
  ret$logPrior <- log(x$prior)

  ## variables
  vars <- str_trim(str_split(str_split(tail(x$model, 1), '\\~')[[1]][2],
                             '\\+')[[1]][-1])

  ## TRUE / FALSE lines
  for (v in vars)
  {
    ret[[v]] <- sapply(x$model, function(one) grepl(v, one))
  }

  ret <- ret[order(ret$posterior, decreasing = TRUE),]
  return(ret)
}
