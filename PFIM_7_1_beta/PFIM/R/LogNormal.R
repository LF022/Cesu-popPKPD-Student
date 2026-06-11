#' @description The class \code{LogNormal} implements the LogNormal distribution.
#' @title LogNormal
#' @inheritParams Distribution
#' @include Distribution.R
#' @export

LogNormal = new_class( "LogNormal", package = "PFIM", parent = Distribution )

# ==============================================================================
#' adjustGradient: adjust the gradient for the log normal distribution.
#' @name adjustGradient
#' @param distribution An object \code{Distribution} giving the distribution.
#' @param gradient The gradient of the model responses.
#' @param thetaValue parameter value
#' @return The adjusted gradient of the model responses.
#' @export
# ==============================================================================

method( adjustGradient, LogNormal ) = function( distribution, gradient, thetaValue ) {
  return( gradient * thetaValue )
}
