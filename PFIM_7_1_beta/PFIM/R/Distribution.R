#' @description
#' The class \code{Distribution} represents and stores information for the parameter distribution.
#' @title Distribution
#' @param name A string giving the name of the distribution.
#' @param mu A double giving the mean mu.
#' @param omega A double giving omega.
#' @export
#'
Distribution = new_class(
  "Distribution",
  package = "PFIM",
  properties = list(
    name = new_property(class_character, default = character(0)),
    mu = new_property(class_double, default = 0.0),
    omega = new_property(class_double, default = 0.0)
  ))


adjustGradient = new_generic( "adjustGradient", c( "distribution" ),
                              function( distribution, gradient, thetaValue ) {
                                S7_dispatch()
                              })
