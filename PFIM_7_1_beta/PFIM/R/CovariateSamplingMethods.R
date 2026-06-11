#' @description The class \code{CovariateSamplingMethods} represents the CovariateSamplingMethods
#' @title CovariateSamplingMethods
#' @param data A data.frame
#' @param copula A Copula object
#' @param distributions A list of distributions
#' @param replace A logical value indicating whether to sample with replacement
#' @param MCSamples A numeric vector for Monte Carlo samples
# @include Copula.R
#' @noRd
#' @export

CovariateSamplingMethods = new_class("CovariateSamplingMethods", package = "PFIM",
                                     properties = list(
                                       data = new_property(class_data.frame, default = data.frame()),
                                       #copula = new_property(Copula, default = NULL),
                                       distributions = new_property(class_list, default = list()),
                                       replace = new_property(class_logical, default = FALSE),
                                       MCSamples = new_property(class_double, default = numeric(0))
                                     ))
