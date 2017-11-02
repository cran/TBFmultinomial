PIPs_by_landmarking <- function(fullModel, data, discreteSurv = TRUE,
                                numberCores = 1,
                                package = 'nnet', maxit = 150,
                                prior = 'flat',
                                method = 'LEB',
                                landmarkLength = 1, lastlandmark,
                                timeVariableName){
  #' Posterior inclusion probabilities (PIPs) by landmarking
  #'
  #'
  #' This function gives us the PIPs for each landmark.
  #' @param fullModel formula of the model including all potential variables
  #' @param data the data frame with all the information
  #' @param discreteSurv Boolean variable telling us whether a 'simple'
  #' multinomial regression is looked for or if the goal is a discrete
  #' survival-time model for multiple modes of failure is needed.
  #' @param numberCores How many cores should be used in parallel?
  #' @param package Which package should be used to fit the models; by default
  #' the \code{nnet} package is used; we could also specify to use the package
  #' 'VGAM'
  #' @param maxit Only needs to be specified with package \code{nnet}: maximal
  #' number of iterations
  #' @param method Method for the g definition
  #' @param prior Prior on the model space
  #' @param landmarkLength Length of the landmark, by default we use each day
  #' @param lastlandmark Where will be the last landmark?
  #' @param timeVariableName What is the name of the variable indicating time?
  #' @return a list with the PIPs for each landmark
  #' @examples
  #' # extract the data:
  #' data("VAP_data")
  #'
  #' # the definition of the full model with three potential predictors:
  #' FULL <- outcome ~ ns(day, df = 4) + gender + type + SOFA
  #' # here we define time as a spline with 3 knots
  #'
  #' PIPs_landmark <- PIPs_by_landmarking(fullModel = FULL, data = VAP_data,
  #'                                      discreteSurv = TRUE, numberCores = 1,
  #'                                      package = 'nnet', maxit = 150,
  #'                                      prior = 'flat',  method = 'LEB',
  #'                                      landmarkLength = 7, lastlandmark = 21,
  #'                                      timeVariableName = 'day')
  #' @export
  #'
  #' @author Rachel Heyard
  #'

  incl <- vector(mode = 'list', length = lastlandmark/landmarkLength)
  if (landmarkLength > 1) landmark <- 0 else landmark <- 1
  for (i in 1:((lastlandmark+landmarkLength)/landmarkLength)){
    landmark_data_subset <- data[which(data[, timeVariableName] >= landmark),]
    PMP <- PMP(fullModel= fullModel, data = landmark_data_subset,
               discreteSurv = discreteSurv, method = method, prior = prior,
               package = package, maxit = maxit,
               numberCores = numberCores)
    landmark <- landmark + landmarkLength
    incl[[i]] <- postInclusionProb(object = PMP)
  }
  return(incl)
}

