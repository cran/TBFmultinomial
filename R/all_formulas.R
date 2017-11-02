all_formulas <-
function(fullModel, discreteSurv = TRUE){
  #' Formulas of all the candidate models
  #'
  #'
  #' This function retrieves the formulas of all the candidate models if the
  #' reference model is the null / baseline model.
  #' @param fullModel formula of the model including all potential variables
  #' @param discreteSurv Boolean variable telling us whether a `simple'
  #' multinomial regression is looked for or if the goal is a discrete
  #' survival-time model for multiple modes of failure is needed.
  #' @return character vector with all the formulas; the first one will be the
  #' reference model; the last element will be the full model.
  #' @import stringr
  #' @examples
  #' data("VAP_data")
  #' FULL <- outcome ~ ns(day, df = 4) + male + type + SOFA
  #' models <- TBFmultinomial:::all_formulas(fullModel = FULL,
  #' discreteSurv = TRUE)
  #' # models
  #'
  #' @author Rachel Heyard

  # Tests:
  if ('formula' %in% class(fullModel)) {
    fullModel <- deparse(fullModel)
  }

  # The null Model:
  nullModel <- nullModel_formula(fullModel = fullModel,
                                 discreteSurv = discreteSurv) # the null Model

  # getting out the variables:
  if (discreteSurv){
    vars <- str_trim(str_split(fullModel, '\\+')[[1]][-1])
  }
  if (!discreteSurv){
    vars <- str_trim(str_split(str_split(fullModel,
                                         '\\~')[[1]][-1], '\\+')[[1]])
  }

  p <- length(vars) # Number of variables
  # Initialisations:
  formulaAllModels <- rep("f", 2^p-1) # Initialisation of the vector with
                                      # all the formulas!
  # length of 2**p - 1 because the nullModel will only be added at the end

  # list containing possible combinations of the variables
  comb.list <- lapply(1:p, function(i){combn(vars, i) })

  # First, we create all those formulas where only ONE variable is added to
  # the null model
  formulaAllModels_1var <- sapply(1:p, function(i){
    str_trim(paste(nullModel, "+", comb.list[[1]][ , i]))
  })
  # Then, we create all those formulas where more then one variable is added to
  # the null model
  formulaAllModels_moreVar <- unlist(sapply(comb.list[-1], function(c){
    sapply(1:ncol(c), function(col) {
      paste(nullModel, paste(c[, col], collapse = ' + '), sep = ' + ')
    })
  }) )

  return(c(nullModel, formulaAllModels_1var, formulaAllModels_moreVar))
}
