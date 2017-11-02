BF <- function(g, t){
  term1 <- 1/sqrt(1+g)
  term2 <- exp(0.5*g/(g+1)*t^2)
  return(term1*term2)
}

logBF <- function(g, t){
  term1 <- -0.5*log(1+g)
  term2 <- 0.5*g/(g+1)*t^2
  return(term1+term2)
}

CSVS <- function(g, model, discreteSurv = TRUE, package = 'nnet'){
  #' Cause-specific variable selection (CSVS)
  #'
  #'
  #' This function performs CSVS given a model fitted using the \code{multinom()}
  #' function of the \code{nnet} package or the \code{vglm()} function of the
  #' \code{VGAM} package.
  #' @param g the estimated g, must be fixed to one value
  #' @param model the model fitted using either \code{nnet} or \code{VGAM}
  #' @param discreteSurv Boolean variable telling us whether a 'simple'
  #' multinomial regression is looked for or if the goal is a discrete
  #' survival-time model for multiple modes of failure is needed.
  #' @param package Which package has been used to fit the model, \code{nnet}
  #' or \code{VGAM}?
  #' @import nnet VGAM
  #' @importFrom utils getFromNamespace
  #' @examples
  #' # data extraction:
  #' data("VAP_data")
  #'
  #' # the definition of the full model with three potential predictors:
  #' FULL <- outcome ~ ns(day, df = 4) + gender + type + SOFA
  #' # here the define time as a spline with 3 knots
  #'
  #' # we first need to fit the multinomial model:
  #' model_full <- multinom(formula = FULL, data = VAP_data,
  #'                        maxit = 150, trace = FALSE)
  #'
  #' G <- 9 # let's suppose g equals to nine
  #'
  #' # then we proceed to CSVS
  #' CSVS_nnet <- CSVS(g = G, model = model_full,
  #'                   discreteSurv = TRUE, package = 'nnet')
  #' @export
  #' @author Rachel Heyard

  sum <- summary(model)
  nbIntercepts <-
    ifelse((discreteSurv & package == 'VGAM'),
           length(sum@extra$colnames.y)*length(sum@assign[[2]]),
           ifelse((!discreteSurv & package == 'VGAM'),
                  length(sum@extra$colnames.y)-1,
                  ifelse(discreteSurv & package == 'nnet',
                         (as.numeric(str_split(attr(sum$terms, 'dataClasses')[2],
                                               pattern = '\\.')[[1]][2]) + 1),
                         1)))
  all_coefs <- if (package == 'VGAM') sum@coef3[ , 1] else sum$coefficients

  intercepts <- if (package == 'VGAM') all_coefs[1:nbIntercepts] else
    all_coefs[ ,(1:nbIntercepts)] # The intercepts
  betas <- if (package == 'VGAM') all_coefs[-(1:nbIntercepts)] else
    all_coefs[ , -(1:nbIntercepts)]
  # the original betas / actually the log odds who's importance need to be evaluated
  sd_betas <- if (package == 'VGAM') sum@coef3[-(1:nbIntercepts), 2] else
    sum$standard.errors[ , -(1:nbIntercepts)]
  # the standard errors of the betas
  standardizedCoefs <- betas/sd_betas # the standardized coefficients to be evaluated
  if (package == 'nnet') {
    betas <- matrix(betas, nrow = 1, byrow = TRUE)[1,]
    all_coefs <- matrix(all_coefs, nrow = 1, byrow = TRUE)[1,]
    standardizedCoefs <- matrix(standardizedCoefs, nrow = 1, byrow = TRUE)[1,]
    names(standardizedCoefs) <- names(betas) <-
      unlist(lapply(colnames(sum$coefficients)[-(1:nbIntercepts)],
                    function(out) paste(out, 1:(length(sum$lev)-1),
                                        sep =  ':')))
    names(all_coefs) <- unlist(lapply(colnames(sum$coefficients),
                                      function(out) paste(out,
                                                          1:(length(sum$lev)-1),
                                                          sep =  ':')))
  }
  logBFs <- logBF(g, standardizedCoefs) # the log Bayes factor of the coefficients
  # Which coefs will be set to 0, and which ones won't

  betas_not0 <- betas[which(logBFs > 0)]
  # the original coefficients that are not set to 0
  betas_tobe0 <- betas[which(logBFs <= 0)]
  betas_0 <- rep(0, length(which(logBFs <= 0)))

  order <- c(1:(ifelse(package == 'VGAM', nbIntercepts,
                       nbIntercepts * (length(sum$lev)-1))), # first the intercepts (that need to be corrected too)
             which(logBFs > 0) + (ifelse(package == 'VGAM', nbIntercepts,
                                         nbIntercepts * (length(sum$lev)-1))), # than the non zero coefficients (that need to be corrected)
             which(logBFs <= 0) + (ifelse(package == 'VGAM', nbIntercepts,
                                          nbIntercepts * (length(sum$lev)-1)))) # the zero coefficients
  if (package == 'nnet') vcov.multinom <- getFromNamespace('vcov.multinom', 'nnet')
  covVarMat <- if (package == 'VGAM') sum@cov.unscaled else
    vcov.multinom(model)
  # The function was originally written for the VGAM package, so we try to
  # rearrange the covariance matrix of nnet to the same structure as the cov
  # matrix of the vgam package
  if(package == 'nnet'){
    covVarMat_new <- matrix(NA, nrow = length(all_coefs),
                            ncol = length(all_coefs))
    colnames(covVarMat_new) <- rownames(covVarMat_new) <-
      names(all_coefs)

    rownames(covVarMat) <- colnames(covVarMat) <-
      unname(sapply(rownames(covVarMat), function(rn){
        ele <- str_split(rn, ':')[[1]]
        paste(ele[2], which(ele[1] == sum$lev[-1]), sep = ':')
        }))
    newOrder <-
      sapply(colnames(covVarMat), function(cn){
        which(colnames(covVarMat_new) == cn)
      })

    covVarMat_new[newOrder, newOrder] <- covVarMat
    covVarMat <- covVarMat_new
  }

  # Reorder the covariance matrix:
  #     to have first the part of the coeffs that need correction and then the part that are set to 0
  covVarMat_ordered <- matrix(NA, nrow = length(order), ncol = length(order))
  J <- 1
  for (j in order){
    I <- 1
    for (i in order){
      covVarMat_ordered[J, I] <- covVarMat[j, i]
      I <- I + 1
    }
    J <- J + 1
  }
  # Matrix reordered

  # The different relevant parts of the covariance matrix:
  if (package == 'VGAM') numberIntercepts <- nbIntercepts else
    numberIntercepts <- nbIntercepts * (length(sum$lev)-1)
  sigma_betaNot0 <- covVarMat_ordered[1:(length(which(logBFs > 0)) + numberIntercepts),
                                      1:(length(which(logBFs > 0)) + numberIntercepts)]
  sigma_beta0 <- covVarMat_ordered[(length(which(logBFs > 0)) + numberIntercepts + 1):length(order),
                                   (length(which(logBFs > 0)) + numberIntercepts + 1):length(order)]
  sigma_betaNot0_beta0 <- covVarMat_ordered[(length(which(logBFs > 0)) + numberIntercepts + 1):length(order),
                                            1:(length(which(logBFs > 0)) + numberIntercepts)]
  sigma_beta0_betaNot0 <- covVarMat_ordered[1:(length(which(logBFs > 0)) + numberIntercepts),
                                            (length(which(logBFs > 0)) + numberIntercepts + 1):length(order)]

  # corrected betas:
  alpha_betas_not0_corrected <- c(intercepts, betas_not0) +
    sigma_beta0_betaNot0%*%solve(sigma_beta0)%*%(betas_0 - betas_tobe0)
  sigma_corr <- sigma_betaNot0 - sigma_beta0_betaNot0 %*%
    solve(sigma_beta0) %*% sigma_betaNot0_beta0
  sigma_correct <- rbind(cbind(sigma_corr,
                               matrix(0, ncol= ncol(sigma_beta0),
                                      nrow = nrow(sigma_corr))),
                         matrix(0, ncol = ncol(sigma_beta0) + ncol(sigma_corr),
                                nrow = nrow(sigma_beta0)))
  new_coefs <- c(alpha_betas_not0_corrected, betas_0)
  if (package == 'nnet'){
    interceptNames <-
      unlist(lapply(colnames(sum$coefficients)[1:nbIntercepts],
                    function(out) paste(out, 1:(length(sum$lev)-1),
                                        sep =  ':')))
  } else interceptNames <- names(intercepts)
  names(new_coefs) <- colnames(sigma_correct) <-
    rownames(sigma_correct) <- c(interceptNames, names(betas_not0),
                                 names(betas_tobe0))

  rownames(sigma_corr) <- colnames(sigma_corr) <-
    c(interceptNames, names(betas_not0))

  # reorder:
  new_coefficients <- rep(NA, length(new_coefs))
  sigma_correction <- matrix(NA, nrow = nrow(sigma_correct),
                             ncol = ncol(sigma_correct))

  for (i in 1:length(new_coefs)){
    who <- which(names(new_coefs) == names(all_coefs)[i])
    new_coefficients[i] <- new_coefs[who]
  }
  names(new_coefficients) <- names(all_coefs)
  for (i in 1:nrow(sigma_correct)){
    who <- which(rownames(sigma_correct) == names(all_coefs)[i])
    for (j in 1:ncol(sigma_correct)){
      whoWho <- which(colnames(sigma_correct) == names(all_coefs)[j])
      sigma_correction[i, j] <- sigma_correct[who, whoWho]
    }
  }
  rownames(sigma_correction) <- colnames(sigma_correction) <- names(all_coefs)


  shrunken_standardized_correctedCoefs <-
    c(new_coefficients[1:numberIntercepts] /
        sqrt(diag(sigma_correction[1:numberIntercepts, 1:numberIntercepts])),
      sqrt((g/(g+1))) * new_coefficients[-(1:numberIntercepts)] /
        sqrt(diag(sigma_correction[(-(1:numberIntercepts)),(-(1:numberIntercepts))])))
  whichNA <- which(is.na(shrunken_standardized_correctedCoefs))
  shrunken_standardized_correctedCoefs[whichNA] <- 0
  ret <- list(correctedCoefs = new_coefficients, correctedSD = sigma_correction,
              oldCoefs = all_coefs, old_standardized_coefs = all_coefs / sqrt(diag(covVarMat)),
              old_shrunken_standardized_coefs = sqrt((g/(g+1))) * all_coefs / sqrt(diag(covVarMat)),
              shrunken_correctedCoefs = c(new_coefficients[1:numberIntercepts],
                                          (g/(g+1)) * new_coefficients[-(1:numberIntercepts)]),
              shrunken_oldCoefs = c(all_coefs[1:numberIntercepts],
                                    (g/(g+1)) * all_coefs[-(1:numberIntercepts)]),
              shrunken_correctedSD = (g/(g+1)) * sigma_correction,
              shrunken_standardized_correctedCoefs = shrunken_standardized_correctedCoefs)
  class(ret) <- c('CSVS', class(ret))
  return(ret)

}

