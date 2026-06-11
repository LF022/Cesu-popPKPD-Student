#' @description The class \code{Proportional} is used to defined a model error.
#' @title Proportional
#' @param output A string giving the model error output.
#' @param equation A expression giving the model error equation.
#' @param derivatives A list giving the derivatives of the model error equation.
#' @param sigmaInter A double giving the sigma inter.
#' @param sigmaSlope A double giving the sigma slope
#' @param sigmaInterFixed A Boolean giving if the  sigma inter is fixed or not. - not in the v7.0
#' @param sigmaSlopeFixed A Boolean giving if the  sigma slope is fixed or not. - not in the v7.0
#' @param cError A integer giving the power parameter.
#' @include ModelError.R
#' @export

Proportional = new_class("Proportional", package = "PFIM", parent = ModelError,
                         properties = list(
                           output = new_property(class_character, default = character(0)),
                           equation = new_property(class_expression, default = expression(sigmaSlope)),
                           derivatives = new_property(class_list, default = list()),
                           sigmaInter = new_property(class_double, default = 0.0),
                           sigmaSlope = new_property(class_double, default = 0.0),
                           sigmaInterFixed = new_property(class_logical, default = FALSE),
                           sigmaSlopeFixed = new_property(class_logical, default = FALSE),
                           cError = new_property(class_double, default = 1.0)
                         ),
                         constructor = function(output = character(0),
                                                equation = expression(sigmaSlope),
                                                derivatives = list(),
                                                sigmaInter = 0.0,
                                                sigmaSlope = 0.0,
                                                sigmaInterFixed = FALSE,
                                                sigmaSlopeFixed = FALSE,
                                                cError = 1.0) {
                           new_object(.parent = ModelError,
                                      output = output,
                                      equation = equation,
                                      derivatives = derivatives,
                                      sigmaInter = sigmaInter,
                                      sigmaSlope = sigmaSlope,
                                      sigmaInterFixed = sigmaInterFixed,
                                      sigmaSlopeFixed = sigmaSlopeFixed,
                                      cError = cError)
                         })
