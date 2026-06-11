#' @description
#' The class \code{Fim} represents and stores information for the Fim.
#' @title Fim
#' @param fisherMatrix A matrix giving the numerical values of the Fim.
#' @param shrinkage A vector giving the shrinkage values.
#' @param fixedEffects A matrix giving the numerical values of the fixedEffects of the Fim.
#' @param varianceEffects A matrix giving the numerical values of varianceEffects of the Fim.
#' @param SEAndRSE A data frame giving the value of the SE and RSE.
#' @param condNumberFixedEffects The conditional number of the fixedEffects of the Fim.
#' @param condNumberVarianceEffects The conditional number of the varianceEffects of the Fim.
#' @export

Fim = new_class("Fim", package = "PFIM",
                properties = list(
                  fisherMatrix = new_property(class_double, default = numeric(0)),
                  fixedEffects = new_property(class_double, default = numeric(0)),
                  varianceEffects = new_property(class_double, default = numeric(0)),
                  SEAndRSE = new_property(class_list, default = list()),
                  condNumberFixedEffects = new_property(class_double, default = 0.0),
                  condNumberVarianceEffects = new_property(class_double, default = 0.0),
                  shrinkage = new_property(class_double, default = numeric(0))
                ))

evaluateFim = new_generic( "evaluateFim", c( "fim", "model", "arm" ) )
evaluateVarianceFIM = new_generic( "evaluateVarianceFIM", c( "fim", "model", "arm" ) )
setEvaluationFim = new_generic( "setEvaluationFim", c( "fim" ) )
setOptimalArms = new_generic( "setOptimalArms", c( "fim", "optimizationAlgorithm" ) )
Dcriterion = new_generic( "Dcriterion", c( "fim" ) )
showFIM = new_generic( "showFIM", c( "fim" ) )
plotSEFIM = new_generic( "plotSEFIM", c( "fim", "evaluation" ) )
plotRSEFIM = new_generic( "plotRSEFIM", c( "fim", "evaluation" ) )
plotShrinkage = new_generic( "plotShrinkage", c( "fim" ,"evaluation" ) )
tablesForReport = new_generic( "tablesForReport", c( "fim", "evaluation" ) )
generateReportEvaluation = new_generic( "generateReportEvaluation", c( "fim" ) )
generateReportOptimization = new_generic( "generateReportOptimization", c( "fim","optimizationAlgorithm" ) )

# ==============================================================================
#' Dcriterion: get the D-criterion of the Fim.
#' @name Dcriterion
#' @param Fim A object \code{Fim} giving the Fim.
#' @return A double giving the D-criterion of the Fim.
#' @export
# ==============================================================================

method( Dcriterion, Fim ) = function( fim )
{
  fisherMatrix = prop( fim, "fisherMatrix" )
  Dcriterion = det(fisherMatrix)**(1/dim(fisherMatrix)[1])
  return(Dcriterion)
}

