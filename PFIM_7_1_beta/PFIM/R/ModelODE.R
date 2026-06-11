#' @description The class \code{ModelODE} is used to defined a ode model.
#' @title ModelODE
#' @inheritParams Model
#' @include Model.R
#' @export

ModelODE = new_class( "ModelODE", package = "PFIM", parent = Model )

# ==============================================================================
#' evaluateInitialConditions: evaluate the initial conditions.
#' @name evaluateInitialConditions
#' @param arm A object of class \code{Arm} giving the arm.
#' @param model A object of class \code{Model} giving the model.
#' @return A list giving the evaluated initial conditions.
#' @export
# ==============================================================================

method( evaluateInitialConditions, ModelODE ) = function( model, arm ) {

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
