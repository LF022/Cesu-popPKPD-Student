#' @description The class \code{Normal} implements the Normal distribution.
#' @title Normal
#' @inheritParams Distribution
#' @include Distribution.R
#' @export

Normal = new_class( "Normal", package = "PFIM", parent = Distribution )

method( adjustGradient, Normal ) = function( distribution, gradient, thetaValue ) {
  return( gradient )
}
