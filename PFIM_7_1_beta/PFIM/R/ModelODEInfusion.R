#' @description The class \code{ModelODEInfusion} is used to defined a model ModelODEInfusion.
#' @title ModelODEInfusion
#' @inheritParams ModelInfusion
#' @export

ModelODEInfusion = new_class( "ModelODEInfusion", package = "PFIM", parent = ModelInfusion )

# ==============================================================================
#' evaluateInitialConditions: evaluate the initial conditions.
#' @name evaluateInitialConditions
#' @param arm A object of class \code{Arm} giving the arm.
#' @param model A object of class \code{ModelODEInfusion} giving the model.
#' @param doseEvent A data frame giving the dose event for the ode solver.
#' @export
# ==============================================================================

method( evaluateInitialConditions, ModelODEInfusion ) = function( model, arm ) {

  initialConditions = prop( arm, "initialConditions")
  parameters = prop( model, "modelParameters")

  # assign mu values of the parameters
  mu = set_names(
    map( parameters, ~ {
      pluck(.x, "distribution", "mu")
    }),
    map( parameters, ~ prop( .x, "name") )
  )

  list2env( mu, envir = environment())

  # evaluate the initial conditions
  initialConditions = map( initialConditions, ~ {

    if ( is.numeric(.x) ) {
      return(.x)
    } else {
      eval( parse( text = .x ) )
    }
  })%>%unlist()

  return( initialConditions )
}





