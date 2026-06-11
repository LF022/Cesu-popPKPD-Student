#' @description The class \code{PSOAlgorithm} implements the PSO algorithm.
#' @title PSOAlgorithm
#' @inheritParams Optimization
#' @param maxIteration A numeric giving the maxIteration.
#' @param populationSize A numeric giving the populationSize.
#' @param seed A numeric giving the seed.
#' @param personalLearningCoefficient A numeric giving the personalLearningCoefficient.
#' @param globalLearningCoefficient A numeric giving the globalLearningCoefficient.
#' @param showProcess A Boolean giving the showProcess.
#' @include Optimization.R
#' @export

PSOAlgorithm = new_class( "PSOAlgorithm", package = "PFIM", parent = Optimization,

                           properties = list( maxIteration = new_property(class_double, default = numeric(0)),
                                              populationSize = new_property(class_double, default = numeric(0)),
                                              seed = new_property(class_double, default = numeric(0)),
                                              personalLearningCoefficient = new_property(class_double, default = numeric(0)),
                                              globalLearningCoefficient = new_property(class_double, default = numeric(0)),
                                              showProcess = new_property(class_logical, default = FALSE ) ) )

# ==============================================================================
#' Optimization PSOAlgorithm
#' @name optimizeDesign
#' @param optimizationObject A object \code{Optimization}.
#' @param optimizationAlgorithm A object \code{PSOAlgorithm}.
#' @return The object \code{optimizationObject} with the slots updated.
#' @export
# ==============================================================================

method( optimizeDesign, list( Optimization, PSOAlgorithm ) ) = function( optimizationObject, optimizationAlgorithm ) {

  # parameters of the optimization algorithm
  optimizerParameters = prop( optimizationObject, "optimizerParameters")
  populationSize = optimizerParameters$populationSize
  maxIteration = optimizerParameters$maxIteration
  personalLearningCoefficient = optimizerParameters$personalLearningCoefficient
  globalLearningCoefficient = optimizerParameters$globalLearningCoefficient
  showProcess = optimizerParameters$showProcess
  seed = optimizerParameters$seed
  set.seed( seed )

  # get the designs
  designs = prop( optimizationObject, "designs" )
  design = pluck( designs, 1 )
  optimalDesign = design

  # check validity of the samplingTimesConstraints
  checkValiditySamplingConstraint( optimalDesign )

  # in case of partial sampling constraints set the new arms with the sampling constraints
  design = setSamplingConstraintForOptimization( optimalDesign )

  # get the arms
  arms = prop( optimalDesign, "arms" )

  # generate the initial population
  cost = list()
  bestCost = rep( Inf, populationSize )
  globalBestCost = Inf

  globalBestPosition = list()
  bestPosition = list()
  velocity = list()
  position = list()

  # get arm outcomes
  outcomes = map( set_names( arms, map_chr( arms, ~ prop( .x, 'name' ) ) ), ~ {
    map_chr( prop( .x, "samplingTimesConstraints" ), ~ prop( .x, "outcome"))
  })

  # constrictionFactor
  phi1 = personalLearningCoefficient
  phi2 = globalLearningCoefficient
  phi = phi1 + phi2
  kappa = 1
  constrictionFactor = 2*kappa / abs( 2 - phi - sqrt( phi * ( phi - 4 ) ) )

  for ( iterPop in 1:populationSize )
  {
    armsList = list()
    for ( arm in arms )
    {
      # Extract values from the arm object
      armName = prop( arm, 'name' )
      samplingTimes = prop( arm, "samplingTimes" )
      samplingTimesConstraints = prop( arm, "samplingTimesConstraints" )

      for ( samplingTimesConstraint in samplingTimesConstraints )
      {
        # Generate samplings based on SamplingConstraints
        samplingsFromSamplingConstraint = generateSamplingsFromSamplingConstraints( samplingTimesConstraint )

        outcome = prop( samplingTimesConstraint, "outcome" )

        # Set initial position for the outcome
        position[[armName]][[outcome]][[iterPop]] = samplingsFromSamplingConstraint

        # Set the best position for the outcome
        bestPosition[[armName]][[outcome]][[iterPop]] = samplingsFromSamplingConstraint

        # Set initial velocity for the outcome
        velocity[[armName]][[outcome]][[iterPop]] = 0.0

        samplingTimes = samplingTimes %>%
          modify_at(
            .at = which(map_chr(., ~ prop( .x, "outcome" ) ) == outcome),
            .f = ~ {
              prop(.x, "samplings") = samplingsFromSamplingConstraint
              .x
            })
      }
      prop(arm, "samplingTimes") = samplingTimes
      armsList = append( armsList, arm )
    }

    # set new arms
    prop( optimalDesign, "arms" ) = armsList

    # evaluate the FIMs
    evaluationFIM = Evaluation( name = "",

                                modelEquations = prop( optimizationObject, "modelEquations" ),
                                modelParameters = prop( optimizationObject, "modelParameters" ),
                                modelError = prop( optimizationObject, "modelError" ),
                                fimType = prop( optimizationObject, "fimType" ),
                                outputs = prop( optimizationObject, "outputs" ),
                                designs = list( optimalDesign ),
                                odeSolverParameters = prop( optimizationObject, "odeSolverParameters" ) )

    evaluationFIM = run( evaluationFIM )

    # get D-criterion
    fim = prop( evaluationFIM, "fim" )
    cost[[ iterPop ]] = 1/Dcriterion( fim )

    # update bestPosition
    bestPosition = position

    # update bestCost
    bestCost[[ iterPop ]] = cost[[ iterPop ]]
    indexMinBestCost = which.min( unlist( bestCost ) )

    # update Global Best
    if ( bestCost[[ indexMinBestCost ]] < globalBestCost )
    {
      for ( arm in arms )
      {
        armName = prop( arm, "name")

        for ( outcome in outcomes[[ armName ]] )
        {
          globalBestPosition[[ armName ]][[ outcome ]] = bestPosition[[ armName ]][[ outcome ]][[ indexMinBestCost ]]
        }
      }
      globalBestCost = bestCost[[ indexMinBestCost ]]
    }
  } # end iterPop

  # Run the PSO
  for ( iteration in 1:maxIteration )
  {
    # show process
    if ( showProcess == TRUE )
    {
      print( paste0( "iter = ", iteration ) )
    }

    for ( iterPop in 1:populationSize )
    {
      # update Velocity

      for ( arm in arms )
      {
        armName = prop( arm, "name")

        for ( outcome in outcomes[[ armName ]] )
        {
          n = length( globalBestPosition[[ armName ]][[ outcome ]] )

          velocity[[ armName ]][[ outcome ]][[ iterPop ]] =
            constrictionFactor *
            ( velocity[[ armName ]][[ outcome ]][[ iterPop ]] +
                phi1 * runif(1,0,1) * ( bestPosition[[ armName ]][[ outcome ]][[ iterPop ]] - position[[ armName ]][[ outcome ]][[ iterPop ]] ) +
                phi2 * runif(1,0,1) * ( globalBestPosition[[ armName ]][[ outcome ]] - position[[ armName ]][[ outcome ]][[ iterPop ]] ) )
        }
      }

      # update Position
      for ( arm in arms )
      {
        armName = prop( arm, "name")

        for ( outcome in outcomes[[ armName ]] )
        {
          position[[ armName ]][[ outcome ]][[ iterPop ]] = position[[ armName ]][[ outcome ]][[ iterPop ]] + velocity[[ armName ]][[ outcome ]][[ iterPop ]]
          position[[ armName ]][[ outcome ]][[ iterPop ]] = sort( position[[ armName ]][[ outcome ]][[ iterPop ]] )
        }
      }

      # apply position limits
      for ( arm in arms )
      {
        armName = prop( arm, "name")
        samplingTimesConstraints = prop( arm, "samplingTimesConstraints" )

        for ( outcome in outcomes[[armName]] )
        {
          samplingTimesConstraint = keep( samplingTimesConstraints, ~ prop(.x,"outcome") == outcome ) %>% pluck(1)
          samplingsWindows = prop( samplingTimesConstraint, "samplingsWindows" )
          positionLength = length( position[[ armName ]][[ outcome ]][[ iterPop ]] )

          for( j in 1:positionLength )
          {
            distance = list()

            iter = 1
            for ( samplingsWindow in samplingsWindows )
            {
              tmp = map( samplingsWindow, ~ abs( .x - position[[armName]][[outcome]][[iterPop]][j] ) )

              distance[[iter]] = unlist( tmp )
              iter = iter+1
            }

            indDistanceMin = which.min( map_dbl( distance, ~ min(.x) ) )

            maxSamplings = max( samplingsWindows[[indDistanceMin]])
            minSamplings = min( samplingsWindows[[indDistanceMin]])

            position[[ armName ]][[ outcome ]][[ iterPop ]][ j ] = min( maxSamplings, position[[ armName ]][[ outcome ]][[ iterPop ]][ j ] )
            position[[ armName ]][[ outcome ]][[ iterPop ]][ j ] = max( minSamplings, position[[ armName ]][[ outcome ]][[ iterPop ]][ j ] )
          }
          position[[ armName ]][[ outcome ]][[ iterPop ]] = sort( position[[ armName ]][[ outcome ]][[ iterPop ]] )
        }
      }

      # constraints
      samplingConstraintsForMetaheuristic = list()

      for ( arm in arms )
      {
        armName = prop( arm, "name")
        samplingTimesConstraints = prop( arm, "samplingTimesConstraints" )

        # check the constraints on the samplings times
        for ( outcome in outcomes[[armName]] )
        {
          newSamplings = position[[ armName ]][[ outcome ]][[ iterPop ]]
          samplingTimesConstraint = keep( samplingTimesConstraints, ~ prop( .x,"outcome" ) == outcome ) %>% pluck(1)
          samplingConstraintsForMetaheuristic[[ armName ]][[ outcome ]] = checkSamplingTimeConstraintsForMetaheuristic( samplingTimesConstraint, arm, newSamplings, outcome )
        }
      }

      samplingConstraintsForMetaheuristic = unlist( samplingConstraintsForMetaheuristic )

      if( all( samplingConstraintsForMetaheuristic ) == TRUE )
      {
        arms = prop( optimalDesign, "arms")
        # set new sampling times
        armsList = list()
        for ( arm in arms )
        {
          armName = prop(arm, 'name')
          samplingTimes = prop(arm, "samplingTimes")

          for ( outcome in outcomes[[armName]] )
          {
            samplingTimes = samplingTimes %>%
              modify_at(
                .at = which( map_chr(., ~ prop( .x, "outcome" ) ) == outcome),
                .f = ~ {
                  prop(.x, "samplings") = position[[armName]][[outcome]][[iterPop]]
                  .x
                })
          }
          prop(arm, "samplingTimes") = samplingTimes
          armsList = append( armsList, arm )
        }

        # set new arms to the design
        prop( optimalDesign, "arms" ) = armsList

        # set and evaluate new design with the constraints
        evaluationFIM = Evaluation( name = "",
                                    modelEquations = prop( optimizationObject, "modelEquations" ),
                                    modelParameters = prop( optimizationObject, "modelParameters" ),
                                    modelError = prop( optimizationObject, "modelError" ),
                                    fimType = prop( optimizationObject, "fimType" ),
                                    outputs = prop( optimizationObject, "outputs" ),
                                    designs = list( optimalDesign ),
                                    odeSolverParameters = prop( optimizationObject, "odeSolverParameters" ) )

        evaluationFIM = run( evaluationFIM )

        # set the cost
        fim = prop( evaluationFIM, "fim" )

        cost[[ iterPop ]] = 1/Dcriterion( fim )

        if (is.nan( cost[[ iterPop ]]))
        {
          cost[[ iterPop ]] = Inf
        }

        if (is.nan( bestCost[[ iterPop ]]))
        {
          bestCost[[ iterPop ]] = Inf
        }

        # Update Personal Best
        if ( cost[[ iterPop ]] < bestCost[[ iterPop ]] )
        {
          for ( arm in arms )
          {
            armName = prop( arm, "name" )

            for ( outcome in outcomes[[armName]] )
            {
              bestPosition[[ armName ]][[ outcome ]][[ iterPop ]] = position[[ armName ]][[ outcome ]][[ iterPop ]]
            }
            bestCost[[ iterPop ]] = cost[[ iterPop ]]
          }

          indexMinBestCost = which.min( bestCost )

          if ( bestCost[[ indexMinBestCost ]] < globalBestCost )
          {
            for ( arm in arms )
            {
              armName = prop( arm, "name" )

              for ( outcome in outcomes[[armName]] )
              {
                globalBestPosition[[ armName ]][[ outcome ]] = bestPosition[[ armName ]][[ outcome ]][[ indexMinBestCost ]]
              }
            }
            globalBestCost = bestCost[[ indexMinBestCost ]]
          }
        }
      }
    } # end iter pop

    if ( showProcess == TRUE )
    {
      print( paste0( "globalBestCost = ", 1/globalBestCost ) )
    }
  } # end iteration

  # evaluate the optimal design
  evaluationOptimalDesign = Evaluation( name = "",
                                        modelEquations = prop( optimizationObject, "modelEquations" ),
                                        modelParameters = prop( optimizationObject, "modelParameters" ),
                                        modelError = prop( optimizationObject, "modelError" ),
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
                                        designs = list( design ),
                                        fimType = prop( optimizationObject, "fimType" ),
                                        outputs = prop( optimizationObject, "outputs" ),
                                        odeSolverParameters = prop( optimizationObject, "odeSolverParameters" ) )

  evaluationInitialDesign = run( evaluationInitialDesign )

  # set the results
  prop( optimizationObject, "optimisationDesign" ) = list( evaluationInitialDesign = evaluationInitialDesign, evaluationOptimalDesign = evaluationOptimalDesign )
  prop( optimizationObject, "optimisationAlgorithmOutputs" ) = list( "optimizationAlgorithm" = optimizationAlgorithm, "optimalArms" = armsList )
  return( optimizationObject )
}

# ==============================================================================
#' constraintsTableForReport: table of the PSOAlgorithm constraints for the report.
#' @name constraintsTableForReport
#' @param optimizationAlgorithm A object \code{PSOAlgorithm}.
#' @param arms List of the arms.
#' @return The table for the constraints in the arms.
#' @export
# ==============================================================================

method( constraintsTableForReport, PSOAlgorithm ) = function( optimizationAlgorithm, arms  )
{
  armsConstraints = map( pluck( arms, 1 ) , ~ getArmConstraints( .x, optimizationAlgorithm ) )
  armsConstraints = map_dfr( armsConstraints, ~ map_df(.x, ~ as.data.frame(.x, stringsAsFactors = FALSE)))
  colnames( armsConstraints ) = c( "Arms name" , "Number of subjects", "Outcome", "Initial samplings", "Samplings windows", "Number of times by windows","Min sampling" )
  armsConstraintsTable = kbl( armsConstraints, align = c( "l","c","c","c","c","c","c") ) %>% kable_styling( bootstrap_options = c(  "hover" ), full_width = FALSE, position = "center", font_size = 13 )
  return( armsConstraintsTable )
}







