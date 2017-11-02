plot_CSVS <- function(CSVSobject, namesVar = NULL, shrunken = FALSE,
                      standardized = FALSE, numberIntercepts, ...){
  ##' Plot a CSVS object
  ##'
  ##' @param CSVSobject valid \code{\link{CSVS}} object
  ##' @param namesVar names of the variables
  ##' @param shrunken should the coefficients be shrunken?
  ##' @param standardized should the coefficients be standardized?
  ##' @param numberIntercepts how many cause-specific intercepts are in the
  ##' model for each outcome
  ##' @param ... parameters for plot
  ##' @importFrom plotrix color2D.matplot
  ##' @importFrom graphics par axis
  ##' @export
  ##' @author Rachel Heyard

  corrected <- if (shrunken & standardized){
    CSVSobject$shrunken_standardized_correctedCoefs
  } else{
      if (shrunken & !standardized) {
        CSVSobject$shrunken_correctedCoefs
      } else{
        if (!shrunken & !standardized){
          CSVSobject$correctedCoefs
        } else stop('standardized and non shrunken coefs not possible')
      }
  }

  uncorrected <- if (shrunken & standardized) {
    CSVSobject$old_shrunken_standardized_coefs
  } else {
    if (shrunken & !standardized){
      CSVSobject$shrunken_oldCoefs
    } else{
      if (!shrunken & !standardized){
        CSVSobject$oldCoefs
      } else stop('standardized and non shrunken coefs not possible')
    }
  }

  outcomes <- unique(sapply(names(corrected),
                            function(x) str_split(x, ':')[[1]][2]))
  outcomesSpec_corr <- outcomesSpec_UNcorr <- lapply(outcomes, function(o) NA )
  #names(outcomesSpec) <- paste(':', outcomes, sep= '')
  count <- 1
  for (i in 1:(length(corrected))){
    whichOUT <- which(outcomes == str_split(names(corrected)[i], ':')[[1]][2])
    outcomesSpec_corr[[whichOUT]] <- c(outcomesSpec_corr[[whichOUT]], corrected[i])
    outcomesSpec_UNcorr[[whichOUT]] <- c(outcomesSpec_UNcorr[[whichOUT]],
                                         uncorrected[i])
    count <- count + 1
  }
  outcomesSpec_corr <- lapply(outcomesSpec_corr, function(x) x[-1])
  outcomesSpec_UNcorr <- lapply(outcomesSpec_UNcorr, function(x) x[-1])
  matrix_corr <- matrix(unlist(outcomesSpec_corr), nrow = length(outcomes), byrow = TRUE)
  matrix_UNcorr <- matrix(unlist(outcomesSpec_UNcorr),
                          nrow = length(outcomes), byrow = TRUE)
  rownames(matrix_corr) <- rownames(matrix_UNcorr) <- outcomes
  colnames(matrix_corr) <- colnames(matrix_UNcorr) <-
    unique(sapply(names(corrected), function(n) str_split(n, ':')[[1]][1]))

  par(las=2, mfrow= c(2,1))
  matrix_coefs_var <- matrix_UNcorr[ , -(1:numberIntercepts)]
  matrix_coefs_corr_var <- matrix_corr[ , -(1:numberIntercepts)]
  if (!is.null(namesVar)) {
    colnames(matrix_coefs_var) <- colnames(matrix_coefs_corr_var) <- namesVar
  } else namesVar <- colnames(matrix_coefs_var)
  color2D.matplot(abs(matrix_coefs_var), cs1 = c(1,1),cs2=c(1,0),
                  cs3=c(1,0), #show.values = TRUE,
                  show.legend = TRUE, ylab= '', xlab = '',
                  main= 'Before CSVS', axes=FALSE, ...)
  axis(side = 1, at = seq(0.5, ncol(matrix_coefs_var) - 0.5, by = 1),
       labels = namesVar, cex.axis = 0.8)
  axis(side = 2, at = seq(0.5, nrow(matrix_coefs_var) - 0.5, by = 1),
       labels = rev(rownames(matrix_coefs_var)), cex.axis = .8)

  color2D.matplot(abs(matrix_coefs_corr_var), cs1 = c(1,1),cs2=c(1,0),cs3=c(1,0),
                  #show.values = TRUE,
                  show.legend = TRUE,
                  ylab= '', xlab = '', main = 'After CSVS',
                  axes=FALSE, ...)
  axis(side = 1, at = seq(0.5, ncol(matrix_coefs_corr_var) - 0.5, by = 1),
       labels = namesVar, cex.axis = 0.8)
  axis(side = 2, at = seq(0.5, nrow(matrix_coefs_corr_var) - 0.5, by = 1),
       labels =  rev(rownames(matrix_coefs_var)), cex.axis = .8)
  print(list(before = matrix_coefs_var, after = matrix_coefs_corr_var))
}
