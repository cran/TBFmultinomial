nullModel_formula <-
function(fullModel, discreteSurv = FALSE){
  # Formulas of the null model
  #
  #
  # This function retrieves the formulas of the null model using the fullmodel
  # @param fullModel formula of the model including all potential variables
  # @param discreteSurv Boolean variable telling us whether a 'simple'
  # multinomial regression is looked for or if the goal is a discrete
  # survival-time model for multiple modes of failure is needed.
  # @return formula of the null model
  # @examples
  # data("VAP_data")
  # FULL <- outcome ~ ns(day, df = 4) + male + type + SOFA
  # nullModel <- TBFmultinomial:::nullModel_formula(fullModel = FULL,
  #                                                 discreteSurv = TRUE)
  # nullModel
  #
  # @author Rachel Heyard

  # Tests:
  if ('formula' %in% class(fullModel)) {
    fullModel <- deparse(fullModel)
  }

  y <- str_trim(strsplit(fullModel, '\\~')[[1]][1])
  if (discreteSurv){
    timedep_intercept <- str_trim(str_split(str_split(fullModel, '\\~')[[1]][2], '\\+')[[1]][1])
    nullModel <- paste(y, '~', timedep_intercept, collapse = '')
  }
  if (!discreteSurv){
    nullModel <- paste(y, '~', 1, collapse = '')
  }
  return(nullModel)
}
