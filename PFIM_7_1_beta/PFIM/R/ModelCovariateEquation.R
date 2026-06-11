#' ModelCovariateEquation
#' @description The class \code{ModelCovariateEquation} is an abstract class used to define models for covariates.
#' @title ModelCovariateEquation
#' @param beta beta
#' @param combinedEffect combinedEffect
#' @param value value
#' @export

ModelCovariateEquation = new_class("ModelCovariateEquation",
                                   package = "PFIM",
                                   properties = list(
                                     beta = new_property(class_double, default = 0),
                                     combinedEffect = new_property(class_double, default = 0),
                                     value = new_property(class_double, default = 0)
                                   ))
