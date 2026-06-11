#' @description The class \code{AdministrationConstraints} represents the constraint of an input to the system.
#' The class stores information concerning the constraints for the dosage regimen.
#' @title AdministrationConstraints
#' @param outcome A string giving the outcome for the administration.
#' @param doses A vector of double giving the doses.
#' @export An object of class \code{AdministrationConstraints}.

AdministrationConstraints = new_class( "AdministrationConstraints",
                                       package = "PFIM",
                                       properties = list(
                                         outcome = new_property(class_character, default = character(0)),
                                         doses = new_property(class_list, default = list())
                                       ))
