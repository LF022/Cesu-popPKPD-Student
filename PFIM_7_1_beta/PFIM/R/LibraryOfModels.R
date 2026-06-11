#' @description The class \code{LibraryOfModels} represents and stores information for the LibraryOfModels.
#' @title LibraryOfModels
#' @param models A list giving all the PK and PD models.
#' @export

LibraryOfModels = new_class("LibraryOfModels", package = "PFIM",

                       properties = list(
                         models = new_property( class_list, default = list())))

