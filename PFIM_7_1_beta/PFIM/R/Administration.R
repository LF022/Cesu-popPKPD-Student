#' @description
#' The class \code{Administration} represents the administration and
#' stores information concerning the administration for the dosage regimen.
#' @title Administration
#' @param outcome A string giving the outcome for the administration.
#' @param timeDose A vector of double giving the time doses.
#' @param dose A vector of double giving the doses.
#' @param Tinf A vector of double giving the time for infusion Tinf.
#' @param tau An integer giving the tau value for repeated dose or steady state.
#' @export An object of class \code{Administration}.

Administration = new_class("Administration",
                           package = "PFIM",
                           properties = list(
                             outcome = new_property(class_character, default = character(0)),
                             timeDose = new_property(class_double, default = numeric(0)),
                             dose = new_property(class_double, default = numeric(0)),
                             Tinf = new_property(class_double, default = numeric(0)),
                             tau = new_property(class_double, default = 0.0)
                          ))
