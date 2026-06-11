#' CovariateModelEquation
#' @description The class \code{CovariateModelEquation} is used to defined a CovariateModelEquation.
#' @title CovariateModelEquation
#' @param beta beta coefficient
#' @param combinedEffect combined effect
#' @param value value
#' @export

CovariateModelEquation = new_class("CovariateModelEquation", package = "PFIM",
                       properties = list(
                         beta = new_property(class_double, default = 0),
                         combinedEffect = new_property(class_double, default = 0),
                         value = new_property(class_double, default = 0)
                       ),
                       constructor = function(beta = 0.0,
                                              combinedEffect = 0.0,
                                              value = 0.0) {
                         new_object(
                           CovariateModelEquation,
                           beta = beta,
                           combinedEffect = combinedEffect,
                           value = value#,
                           #.parent = environment()
                         )
                       })

computeCovariateValue = new_generic("computeCovariateValue", c("equation"))
