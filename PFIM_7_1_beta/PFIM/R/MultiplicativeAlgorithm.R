#' @description The class \code{MultiplicativeAlgorithm} implements the multiplicative algorithm.
#' @title MultiplicativeAlgorithm
#' @inheritParams Optimization
#' @param lambda A numeric giving the parameter lambda.
#' @param delta A numeric giving the parameter delta
#' @param numberOfIterations A numeric giving the number of iterations.
#' @param weightThreshold A numeric giving the weight threshold.
#' @param showProcess A Boolean for displaying the process or not.
#' @param multiplicativeAlgorithmOutputs A list giving the output of the optimization algorithm.
#' @include Optimization.R
#' @export

MultiplicativeAlgorithm = new_class("MultiplicativeAlgorithm",
                                    package = "PFIM",
                                    parent = Optimization,

                                    properties = list(
                                      lambda = new_property(class_numeric, default = 0.0),
                                      delta = new_property(class_numeric, default = 0.0),
                                      numberOfIterations = new_property(class_numeric, default = 0),
                                      weightThreshold = new_property(class_numeric, default = 0.0),
                                      showProcess = new_property(class_logical, default = FALSE),
                                      multiplicativeAlgorithmOutputs = new_property(class_list, default = list())
                                    ))

plotWeightsMultiplicativeAlgorithm = new_generic( "plotWeightsMultiplicativeAlgorithm", c( "optimization", "optimizationAlgorithm" ) )

# ==============================================================================
#' Function MultiplicativeAlgorithm_Rcpp
#' @name MultiplicativeAlgorithm_Rcpp
#' @description Run the MultiplicativeAlgorithm_Rcpp in Rcpp.
#' @param fisherMatrices_input The parameter fotfisherMatrices_input.
#' @param numberOfFisherMatrices_input The parameter numberOfFisherMatrices_input.
#' @param weights_input The parameter weights_input.
#' @param numberOfParameters_input The parameter numberOfParameters_input.
#' @param dim_input The parameter dim_input.
#' @param lambda_input The parameter lambda_input.
#' @param delta_input The parameter delta_input.
#' @param iterationInit_input The parameter iterationInit_input.
#' @return The list output with the outputs of the MultiplicativeAlgorithm_Rcpp.
#' @export
# ==============================================================================

MultiplicativeAlgorithm_Rcpp = function(fisherMatrices_input,
                                        numberOfFisherMatrices_input,
                                        weights_input,
                                        numberOfParameters_input,
                                        dim_input,
                                        lambda_input,
                                        delta_input,
                                        iterationInit_input){
  incltxt = '

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h> 	/* Mathematical functions */
#include <time.h>	/* Function time used to initialise the random number generator */
#include <float.h>	/* Implementation related constants */
#include <signal.h>	/* Signal handling used to detect arithmetic errors */

/******************************************************************************
 MultiplicativeAlgorithm_Rcpp
*******************************************************************************/

// [[Rcpp::export]]
List MultiplicativeAlgorithm_Rcpp( List fisherMatrices,
                    arma :: vec  numberOfFisherMatrices,
                    arma :: vec  weights,
                    arma :: vec  numberOfParameters,
                    arma :: vec  dim,
                    arma :: vec  lambda,
                    arma :: vec  delta,
                    arma :: vec  iterationInit)
 {

arma :: mat sum_weighted_fims(dim[0],dim[0]);
arma :: vec determinant;
arma :: vec Dcriteria;
arma :: mat derivative_phi(dim[0],dim[0]);
arma :: vec vector_of_multiplier(numberOfFisherMatrices[0]);
arma :: mat matmult(dim[0],dim[0]);

int iter;
for( iter=0 ; iter < iterationInit[0] ; iter ++){

//Rcout << "iteration = " << iter << std::endl;

// sum_weighted_fims =  weights[i]*fims[i]
int i;
for(i=0 ; i<numberOfFisherMatrices[0] ; i++){
sum_weighted_fims += Rcpp::as<arma::mat>( fisherMatrices[i])*weights[i];
}

// determinant of sum_weighted_fims
determinant = det( sum_weighted_fims );

// D-criteria
Dcriteria = pow(determinant,1/dim[0]);

// derivatives of function phi_D
derivative_phi = Dcriteria[0] * inv(sum_weighted_fims)/dim[0];

// vector of multiplier
for(i=0 ; i<numberOfFisherMatrices[0] ; i++){
matmult = derivative_phi * Rcpp::as<arma::mat>( fisherMatrices[i]);
vector_of_multiplier[i] = sum(matmult.diag());
}

// normalization of the weights
weights = weights % pow(vector_of_multiplier,lambda[0]) / sum(weights % pow(vector_of_multiplier,lambda[0]));

// stop criterion
 if (vector_of_multiplier.max()<(1+delta[0])*sum(weights%vector_of_multiplier))
        {
          break;
        }

} // end iteration

// output
return Rcpp::List::create( Rcpp::Named ("weights") = weights,
                           Rcpp::Named ("iterationEnd") = iter);

} // end MultiplicativeAlgorithm_Rcpp
'

MultiplicativeAlgorithm_Rcpp = inline::cxxfunction(

  signature( fisherMatrices_input = "list",
             numberOfFisherMatrices_input = "integer",
             weights_input = "numeric",
             numberOfParameters_input = "integer",
             dim_input = "integer",
             lambda_input = "numeric",
             delta_input = "numeric",
             iterationInit_input = "integer"),

  plugin = "RcppArmadillo",
  incl = incltxt,
  body = '
          List fisherMatrices = Rcpp::as<List>(fisherMatrices_input);
          arma::vec numberOfFisherMatrices  = Rcpp::as<arma::vec>(numberOfFisherMatrices_input);
          arma::vec weights  = Rcpp::as<arma::vec>(weights_input);
          arma::vec numberOfParameters = Rcpp::as<arma::vec>(numberOfParameters_input);
          arma::vec dim = Rcpp::as<arma::vec>(dim_input);
          arma::vec lambda = Rcpp::as<arma::vec>(lambda_input);
          arma::vec delta = Rcpp::as<arma::vec>(delta_input);
          arma::vec iterationInit = Rcpp::as<arma::vec>(iterationInit_input);

          return Rcpp::wrap( MultiplicativeAlgorithm_Rcpp(  fisherMatrices,
                                                            numberOfFisherMatrices,
                                                            weights,
                                                            numberOfParameters,
                                                            dim,
                                                            lambda,
                                                            delta,
                                                            iterationInit ) );')

output = MultiplicativeAlgorithm_Rcpp( fisherMatrices_input,
                                       numberOfFisherMatrices_input,
                                       weights_input,
                                       numberOfParameters_input,
                                       dim_input, lambda_input,
                                       delta_input,
                                       iterationInit_input )
return( output )

}

# ==============================================================================
#' Optimization MultiplicativeAlgorithm
#' @name optimizeDesign
#' @param optimizationObject A object \code{Optimization}.
#' @param optimizationAlgorithm A object \code{MultiplicativeAlgorithm}.
#' @return The object \code{optimizationObject} with the slots updated.
#' @export
# ==============================================================================

method( optimizeDesign, list( Optimization, MultiplicativeAlgorithm ) ) = function( optimizationObject, optimizationAlgorithm ) {

  # parameters of the optimization algorithm
  optimizerParameters = prop( optimizationObject, "optimizerParameters")
  lambda = optimizerParameters$lambda
  delta = optimizerParameters$delta
  numberOfIterations = optimizerParameters$numberOfIterations
  weightThreshold = optimizerParameters$weightThreshold

  # generate the Fims from administration and sampling times constraints
  fimsFromConstraints = generateFimsFromConstraints( optimizationObject )

  # run the multiplicative algorithm
  designs = prop( optimizationObject, "designs" )
  design = pluck( designs, 1 )
  optimalDesign = pluck( designs, 1 )

  # list for the evaluation of the optimal design
  evaluationOptimalDesignList = list()
  evaluationInitialDesignList = list()

  # design name
  designName = prop( design, "name" )

  # number of arms in the design
  numberOfArms = prop( design, "numberOfArms" )

  # set arms and fims from the evaluation of the constraints
  armFims = fimsFromConstraints$listArms[[designName]]
  fisherMatrices = fimsFromConstraints$listFimsAlgoMult[[designName]]

  # multiplicative algorithm parameters
  numberOfFisherMatrices = length( fisherMatrices )
  weights = rep( 1/numberOfFisherMatrices, numberOfFisherMatrices )
  dim = dim( pluck( fisherMatrices,1 ) )[1]
  numberOfParameters = length( prop( optimizationObject, "modelParameters" ) )

  # run the multiplicative algorithm
  multiplicativeAlgorithmOutput = MultiplicativeAlgorithm_Rcpp( fisherMatrices, numberOfFisherMatrices, weights, numberOfParameters, dim, lambda, delta, numberOfIterations )

  #get the optimal weights
  weights = multiplicativeAlgorithmOutput[["weights"]]
  weightsIndex = which( weights > weightThreshold )
  optimalWeights = weights[ weightsIndex ]

  # set the multiplicativeAlgorithmOutputs
  prop( optimizationAlgorithm, "multiplicativeAlgorithmOutputs" ) = list( armFims = armFims,
                                                                          multiplicativeAlgorithmOutput = multiplicativeAlgorithmOutput,
                                                                          numberOfArms = numberOfArms,
                                                                          weightThreshold = weightThreshold,
                                                                          weightsIndex = weightsIndex,
                                                                          optimalWeights = optimalWeights )

  # set the optimal arms to the optimal design
  fim =  prop( optimizationObject, "fim" )
  optimalArms = setOptimalArms( fim, optimizationAlgorithm )

  # set optimal arms
  prop( optimalDesign, "arms" ) = optimalArms

  # evaluate the optimal design
  evaluationOptimalDesign = Evaluation( name = "....",
                                        modelEquations = prop( optimizationObject, "modelEquations" ),
                                        modelParameters = prop( optimizationObject, "modelParameters" ),
                                        modelError = prop( optimizationObject, "modelError" ),

                                        modelCovariates = prop( optimizationObject, "modelCovariates" ),
                                        modelCovariatesEquation = prop( optimizationObject, "modelCovariatesEquation" ),
                                        numberOfOccasions = prop( optimizationObject, "numberOfOccasions" ),

                                        designs = list( optimalDesign ),
                                        fimType = prop( optimizationObject, "fimType" ),
                                        outputs = prop( optimizationObject, "outputs" ),
                                        odeSolverParameters = prop( optimizationObject, "odeSolverParameters" ) )

  evaluationOptimalDesign = run( evaluationOptimalDesign )

  # evaluate the initial design
  evaluationInitialDesign = Evaluation( name = "",
                                        modelEquations = prop( optimizationObject, "modelEquations" ),
                                        modelParameters = prop( optimizationObject, "modelParameters" ),
                                        modelError = prop( optimizationObject, "modelError" ),

                                        modelCovariates = prop( optimizationObject, "modelCovariates" ),
                                        modelCovariatesEquation = prop( optimizationObject, "modelCovariatesEquation" ),
                                        numberOfOccasions = prop( optimizationObject, "numberOfOccasions" ),

                                        designs = list( design ),
                                        fimType = prop( optimizationObject, "fimType" ),
                                        outputs = prop( optimizationObject, "outputs" ),
                                        odeSolverParameters = prop( optimizationObject, "odeSolverParameters" ) )

  evaluationInitialDesign = run( evaluationInitialDesign )

  # set the results in evaluation
  prop( optimizationObject, "optimisationDesign" ) = list( evaluationInitialDesign = evaluationInitialDesign,
                                                           evaluationOptimalDesign = evaluationOptimalDesign )

  prop( optimizationObject, "optimisationAlgorithmOutputs" ) = list( "optimizationAlgorithm" = optimizationAlgorithm,
                                                                     "optimalArms" = optimalArms, optimalWeights = optimalWeights )

  return( optimizationObject )
}

# ==============================================================================
#' plotWeightsMultiplicativeAlgorithm: plot the optimal weight.
#' @name plotWeightsMultiplicativeAlgorithm
#' @param optimization A object \code{Optimization}.
#' @param optimizationAlgorithm A object \code{MultiplicativeAlgorithm}.
#' @return The graph plotWeight.
#' @export
# ==============================================================================

method( plotWeightsMultiplicativeAlgorithm, list( Optimization, MultiplicativeAlgorithm ) ) = function( optimization, optimizationAlgorithm )
{
  optimisationAlgorithmOutputs = prop( optimization, "optimisationAlgorithmOutputs" )
  optimalArms = optimisationAlgorithmOutputs$optimalArms
  optimalArmsName = map( optimalArms, ~ prop(.x,"name" ) ) %>% unlist()
  optimalWeights = optimisationAlgorithmOutputs$optimalWeights
  optimalArms = data.frame( optimalArmsName, optimalWeights )

  weightPlot = ggplot(optimalArms, aes(x = reorder(optimalArmsName, optimalWeights), y = optimalWeights)) +
    geom_bar(stat = "identity", fill = "gray50") +
    scale_y_continuous(limits = c(0, 1),breaks = seq(0, 1, by = 0.1),minor_breaks = seq(0, 1, by = 0.05),expand = c(0, 0)  ) +
    scale_x_discrete(expand = c(0, 0)) +
    labs(  x = "Arms",  y = "Weights" ) +
    coord_flip() +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.title.x = element_text(color = "black", margin = margin(t = 10)),
      axis.title.y = element_text(color = "black", margin = margin(r = 10)),
      axis.text.x = element_text(color = "black", margin = margin(t = 5)),
      axis.text.y = element_text(color = "black", margin = margin(r = 5)),
      panel.grid.major.x = element_line(color = "gray90", size = 0.5),
      panel.grid.minor.x = element_line(color = "gray95", size = 0.3),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.border = element_rect(color = "gray80", fill = NA, size = 0.5),
      plot.margin = margin(10, 10, 10, 10) )
  return( weightPlot )
}

# ==============================================================================
#' constraintsTableForReport: table of the MultiplicativeAlgorithm constraints for the report.
#' @name constraintsTableForReport
#' @param optimizationAlgorithm A object \code{MultiplicativeAlgorithm}.
#' @param arms List of the arms.
#' @return The table for the constraints in the arms.
#' @export
# ==============================================================================

method( constraintsTableForReport, MultiplicativeAlgorithm ) = function( optimizationAlgorithm, arms  )
{
  armsConstraints = map( pluck( arms, 1 ) , ~ getArmConstraints( .x, optimizationAlgorithm ) )
  armsConstraints = map_df( pluck( armsConstraints, 1), ~ as.data.frame(.x, stringsAsFactors = FALSE ) )
  colnames( armsConstraints ) = c( "Arms name" , "Number of subjects", "Outcome", "Initial samplings", "Fixed times", "Number of samplings optimisable","Dose constraints" )
  armsConstraintsTable = kbl( armsConstraints, align = c( "l","c","c","c","c","c","c") ) %>%
    kable_styling( bootstrap_options = c(  "hover" ), full_width = FALSE, position = "center", font_size = 13 )
  return( armsConstraintsTable )
}
