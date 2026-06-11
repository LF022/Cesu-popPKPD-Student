#' ModelParameter
#' @description The class \code{ModelParameter} is used to defined the model parameters.
#' @title ModelParameter
#' @param name A string giving the name of the parameter.
#' @param gamma A numeric giving the value of the variance for the inter occasion variability (iov).
#' @param distribution A string giving the distribution of the parameter.
#' @param fixedMu A Boolean setting TRUE/FALSE if the mu is estimated or not.
#' @param fixedOmega A Boolean setting TRUE/FALSE if the omega is estimated or not.
#' @param value A numeric giving the value of the covariates effects.
#' @import Distribution.R
#' @export

ModelParameter = new_class( "ModelParameter",
                            package = "PFIM",
                            properties = list(
                              name = new_property(class_character, default = character(0)),
                              gamma = new_property(class_double, default = 0),
                              distribution =  new_property(Distribution, default = NULL),
                              fixedMu = new_property(class_logical, default = FALSE),
                              fixedOmega = new_property(class_logical, default = FALSE),
                              value = new_property(class_double, default = numeric(0))))

getModelParametersData = new_generic( "getModelParametersData", c( "modelParameter" ) )

# ==============================================================================
#' getModelParametersData: get model parameters data for report.
#' @name getModelParametersData
#' @param modelParameter An object if class \code{Model} giving the model.
#' @return A data frame with the data of all the parameters.
#' @export
# ==============================================================================

method( getModelParametersData, ModelParameter ) = function( modelParameter ) {

  modelParametersData = list(modelParameter) %>%
    map(function(parameter) {
      dist = prop(parameter, "distribution")
      list(
        Parameters = prop(parameter, "name"),
        mu = as.character(prop(dist, "mu")),
        omega2 = as.character(prop(dist, "omega")^2),
        Distribution = str_remove(class(dist)[1], "PFIM::"),
        mu_fixed = as.character(prop(parameter, "fixedMu")),
        omega2_fixed = as.character(prop(parameter, "fixedOmega"))
      )
    }) %>%
    map(~ as.data.frame(.x, stringsAsFactors = FALSE)) %>%
    list_rbind()

  return( modelParametersData )

}









