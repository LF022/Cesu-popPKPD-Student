#' @description The class \code{SamplingTimes} is used to defined SamplingTimes.
#' @title SamplingTimes
#' @param outcome A string giving the outcome.
#' @param samplings A vector of numeric giving the samplings.
#' @export

SamplingTimes = new_class("SamplingTimes", package = "PFIM",

                          properties = list(
                            outcome = new_property(class_character, default = character(0)),
                            samplings = new_property(class_double, default = numeric(0))
                          ))
