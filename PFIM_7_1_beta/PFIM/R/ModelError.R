#' ModelError
#' @description The class \code{ModelError} is used to defined a model error.
#' @title ModelError
#' @param output A string giving the model error output.
#' @param equation A expression giving the model error equation.
#' @param derivatives A list giving the derivatives of the model error equation.
#' @param sigmaInter A double giving the sigma inter.
#' @param sigmaSlope A double giving the sigma slope
#' @param sigmaInterFixed A boolean giving if the  sigma inter is fixed or not. - not in the v7.0
#' @param sigmaSlopeFixed A boolean giving if the  sigma slope is fixed or not. - not in the v7.0
#' @param cError A integer giving the power parameter.
#' @export

ModelError = new_class("ModelError", package = "PFIM",
                       properties = list(
                         output = new_property(class_character, default = "output"),
                         equation = new_property(class_expression, default = expression()),
                         derivatives = new_property(class_list, default = list()),
                         sigmaInter = new_property(class_double, default = 0.1),
                         sigmaSlope = new_property(class_double, default = 0.0),
                         sigmaInterFixed = new_property(class_logical, default = FALSE),
                         sigmaSlopeFixed = new_property(class_logical, default = FALSE),
                         cError = new_property(class_double, default = 1.0)
                       ),
                       constructor = function(output = "output",
                                              equation = expression(),
                                              derivatives = list(),
                                              sigmaInter = 0.1,
                                              sigmaSlope = 0.0,
                                              sigmaInterFixed = FALSE,
                                              sigmaSlopeFixed = FALSE,
                                              cError = 1.0) {
                         new_object(
                           output = output,
                           equation = equation,
                           derivatives = derivatives,
                           sigmaInter = sigmaInter,
                           sigmaSlope = sigmaSlope,
                           sigmaInterFixed = sigmaInterFixed,
                           sigmaSlopeFixed = sigmaSlopeFixed,
                           cError = cError,
                           .parent = environment()
                         )
                       })

evaluateErrorModelDerivatives = new_generic( "evaluateErrorModelDerivatives", c( "modelError" ) )
getModelErrorData = new_generic( "getModelErrorData", c( "modelError" ) )

#' evaluateErrorModelDerivatives; evaluate the derivatives of the model error.
#' @name evaluateErrorModelDerivatives
#' @param modelError An object \code{ModelError} that defines the model error.
#' @param evaluationModel A dataframe giving the outputs for the model evaluation.
#' @return The matrices sigmaDerivatives and errorVariance.
#' @export

method( evaluateErrorModelDerivatives, ModelError ) = function( modelError, evaluationModel ) {

  sigmaInter = prop( modelError, "sigmaInter" )
  sigmaSlope = prop( modelError, "sigmaSlope" )

  equation = expression( ( sigmaInter + sigmaSlope * evaluationModel ) ** 2 )
  identityMatrix = diag( length( evaluationModel ) )

  modelErroParameterConditions = list(
    list( name = "sigmaInter", value = prop( modelError, "sigmaInter" ), fixed = prop( modelError, "sigmaInterFixed" ) ),
    list( name = "sigmaSlope", value = prop( modelError, "sigmaSlope" ), fixed = prop( modelError, "sigmaSlopeFixed" ) ) )

  sigmaDerivatives = modelErroParameterConditions %>%
    keep( ~ .x$value != 0 && ! .x$fixed ) %>%
    map(~ {
      sigmaDerivatives = eval( D( equation, .x$name ) ) * identityMatrix
      setNames( list(sigmaDerivatives ), .x$name )
    }) %>% flatten()

  errorVariance = ( sigmaInter + sigmaSlope * evaluationModel )**2 * identityMatrix

  return( list( sigmaDerivatives = sigmaDerivatives, errorVariance = errorVariance ) )
}

#' getModelErrorData: get the parameters sigma slope and sigma inter (used for the report).
#' @name getModelErrorData
#' @param modelError An object \code{ModelError} that defines the model error.
#' @return A list of dataframe with outcome, type of model error and sigma slope and inter.
#' @export

method( getModelErrorData, ModelError ) = function( modelError ) {
  modelErrorData = list(modelError) %>%
    map(function( model ) {
      list(
        outcome = prop( model, "output" ),
        type = str_remove(class(model)[1], "PFIM::"),
        sigmaSlope = as.character(prop(model, "sigmaSlope")),
        sigmaInter = as.character(prop(model, "sigmaInter"))
      )
    }) %>% map(~ as.data.frame(.x, stringsAsFactors = FALSE)) %>% list_rbind()
  return( modelErrorData )
}
