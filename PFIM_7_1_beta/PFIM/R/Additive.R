#' Additive
#' @description The class \code{Additive} is used to define Additive model for covariates.
#' @title Additive
#' @param beta beta coefficient
#' @param combinedEffect combined effect
#' @param value value
#' @include CovariateModelEquation.R
#' @export

Additive = new_class("Additive", package = "PFIM", parent = CovariateModelEquation,
                         properties = list(
                           beta = new_property(class_double, default = 0),
                           combinedEffect = new_property(class_double, default = 0),
                           value = new_property(class_double, default = 0)
                         ),
                         constructor = function(beta = 0.0,
                                                combinedEffect = 0.0,
                                                value = 0.0) {
                           new_object(.parent = CovariateModelEquation,
                                      beta = beta,
                                      combinedEffect = combinedEffect,
                                      value = beta + combinedEffect)
                         })

method(computeCovariateValue, Additive) = function(equation, beta, combinedEffect) {
  Additive(beta = beta, combinedEffect = combinedEffect)
}
