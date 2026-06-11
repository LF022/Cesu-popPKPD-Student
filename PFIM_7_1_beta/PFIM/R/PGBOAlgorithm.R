#' @description The class \code{PGBOAlgorithm} implements the PGBO algorithm.
#' @title PGBOAlgorithm
#' @inheritParams Optimization
#' @param N A numeric giving the parameter N.
#' @param muteEffect A numeric giving the parameter muteEffect.
#' @param maxIteration A numeric giving the parameter maxIteration.
#' @param purgeIteration A numeric giving the parameter purgeIteration.
#' @param seed A numeric giving the parameter seed.
#' @param showProcess A Boolean giving showProcess.
#' @include Optimization.R
#' @export

PGBOAlgorithm = new_class( "PGBOAlgorithm", package = "PFIM", parent = Optimization,

                           properties = list( N = new_property(class_double, default = numeric(0)),
                                              muteEffect = new_property(class_double, default = numeric(0)),
                                              maxIteration = new_property(class_double, default = numeric(0)),
                                              purgeIteration = new_property(class_double, default = numeric(0)),
                                              seed = new_property(class_double, default = numeric(0)),
                                              showProcess = new_property(class_logical, default = FALSE )))

# ==============================================================================
#' Optimization PGBOAlgorithm
#' @name optimizeDesign
#' @param optimizationObject A object \code{Optimization}.
#' @param optimizationAlgorithm A object \code{PGBOAlgorithm}.
#' @return The object \code{optimizationObject} with the slots updated.
#' @export
# ==============================================================================

method( optimizeDesign, list( Optimization, PGBOAlgorithm ) ) = function( optimizationObject, optimizationAlgorithm ) {

  results = list()
  # designs

  designs = prop( optimizationObject, "designs" )
  design = pluck( designs, 1 )

  # check validity of the samplingTimesConstraints
  checkValiditySamplingConstraint( design )

  # in case of partial sampling constraints set the new arms with the sampling constraints
  design = setSamplingConstraintForOptimization( design )

  # get the arms
  arms = prop( design, "arms" )

  # get arm outcomes
  outcomes = map( set_names( arms, map_chr( arms, ~ prop( .x, 'name' ) ) ), ~ {
    map_chr( prop( .x, "samplingTimesConstraints" ), ~ prop( .x, "outcome") )
  })

  # set size for checking the constraints
  numberOfArms = length( arms )
  numberOutcomesArmsInConstraints = sum( lengths( outcomes ) )

  # initialize best design
  initialDesign = design
  designA = design
  best = design

  # get pgbo parameters
  optimizerParameters = prop( optimizationObject, "optimizerParameters")
  N = optimizerParameters$N
  muteEffect = optimizerParameters$muteEffect
  maxIteration = optimizerParameters$maxIteration
  purgeIteration = optimizerParameters$purgeIteration
  showProcess = optimizerParameters$showProcess
  seed = optimizerParameters$seed
  set.seed( seed )

  # Generate the initial solution
  samplingInitialSolution = list()

  # Generate the initial solution
  armsList = list()

  for ( arm in arms )
  {
    # Extract values from the arm object
    armName = prop( arm, 'name' )
    samplingTimes = prop( arm, "samplingTimes" )
    samplingTimesConstraints = prop( arm, "samplingTimesConstraints" )

    for ( samplingTimesConstraint in samplingTimesConstraints )
    {
      # Generate samplings based on samplingConstraints
      outcome = prop( samplingTimesConstraint, "outcome" )
      samplingInitialSolution[[ armName ]][[ outcome ]]  = generateSamplingsFromSamplingConstraints( samplingTimesConstraint )

      # Set initial position for the outcome
      samplingTimes = samplingTimes %>%
        modify_at(
          .at = which( map_chr(., ~ prop( .x, "outcome" ) ) == outcome ),
          .f = ~ {
            prop(.x, "samplings") = samplingInitialSolution[[ armName ]][[ outcome ]]
            .x
          })
    }
    prop( arm, "samplingTimes" ) = samplingTimes
    armsList = append( armsList, arm )
  }

  # set new arms
  prop( design, "arms" ) = armsList

  # evaluate the FIMs
  evaluationFIM = Evaluation( name = "",
                              modelEquations = prop( optimizationObject, "modelEquations" ),
                              modelParameters = prop( optimizationObject, "modelParameters" ),
                              modelError = prop( optimizationObject, "modelError" ),
                              fimType = prop( optimizationObject, "fimType" ),
                              outputs = prop( optimizationObject, "outputs" ),
                              designs = list( design ),
                              odeSolverParameters = prop( optimizationObject, "odeSolverParameters" ) )

  evaluationFIM = run( evaluationFIM )

  fim = prop( evaluationFIM, "fim" )

  d = 1/Dcriterion( fim )

  if ( is.finite(d) == FALSE )
  {
    d = 10e6
  }

  # set the PGBO parameters
  fitBase = 0.03
  theta = -log10( fitBase ) / d
  fitA = 10**( -theta * d )
  fitBest = fitA

  # Run PGBO
  arms = prop( designA, "arms" )

  samplingTimesArms = list()

  for ( iteration in 1:maxIteration )
  {
    # Boolean checking constraints in for while
    foundSamplingWithConstraints = FALSE

    while ( foundSamplingWithConstraints == FALSE )
    {
      # select arm
      indexArm = sample( numberOfArms, 1 )
      arm = arms[[indexArm]]
      armName = prop( arm, "name" )

      # sampling constraints
      samplingTimesConstraints = prop( arm, "samplingTimesConstraints" )

      # get samplings
      numberOfOutcome = length( outcomes[[armName]] )
      indexOutcome = sample( numberOfOutcome, 1 )
      outcome = outcomes[[armName]][indexOutcome]

      samplings = prop( arm, "samplingTimes" ) %>%
        keep( ~ prop( .x, "outcome" ) == outcome ) %>%
        pluck(1) %>%
        prop( "samplings" )

      # sampling time mutation
      indexSamplings = sample( length( samplings ), 1 )

      if ( runif( 1 ) < 0.8 )
      {
        samplings[indexSamplings] = samplings[indexSamplings] + rcauchy(1)*muteEffect
      } else{
        samplings = samplings + rnorm( length( samplings ) )*muteEffect
      }

      samplings = sort( samplings )

      # check sampling time constraints
      samplingTimesConstraint = keep( samplingTimesConstraints, ~ prop( .x,"outcome" ) == outcome ) %>% pluck(1)
      samplingConstraintsForMetaheuristic = checkSamplingTimeConstraintsForMetaheuristic( samplingTimesConstraint, arm, samplings, outcome )

      if( all( unlist( samplingConstraintsForMetaheuristic ) ) == TRUE )
      {
        samplingTime = prop( arm, "samplingTimes" ) %>% keep( ~ prop( .x, "outcome" ) == outcome ) %>% pluck(1)
        prop( samplingTime, "samplings" ) = samplings
        samplingTimesArms[[armName]][[outcome]] = samplingTime
      }

      # check if all constraints TRUE and number of constraints = nb Arm * nb Response
      if ( length( unlist( samplingTimesArms ) ) == numberOutcomesArmsInConstraints )
      {
        foundSamplingWithConstraints = TRUE
      }
    } # end while

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
              prop(.x, "samplings") = prop( samplingTimesArms[[armName]][[outcome]], "samplings")
              .x
            })
      }
      prop(arm, "samplingTimes") = samplingTimes
      armsList = append( armsList, arm )
    }

    prop( designA, "arms" ) = armsList

    # Evaluation
    designB = designA

    evaluationFIM = Evaluation( name = "",
                                modelEquations = prop( optimizationObject, "modelEquations" ),
                                modelParameters = prop( optimizationObject, "modelParameters" ),
                                modelError = prop( optimizationObject, "modelError" ),
                                fimType = prop( optimizationObject, "fimType" ),
                                outputs = prop( optimizationObject, "outputs" ),
                                designs = list( designB ),
                                odeSolverParameters = prop( optimizationObject, "odeSolverParameters" ) )

    evaluationFIM = run( evaluationFIM )

    # set the cost
    fim = prop( evaluationFIM, "fim" )

    d = 1/Dcriterion( fim )

    if ( is.finite(d) == FALSE )
    {
      d = 10e6
    }

    fitB = 10**( -theta * d )

    # Update
    if ( !is.nan(fitB) && fitB > 0 )
    {
      if ( fitA == fitB )
      {
        proba = 1/N
      }
      else
      {
        f = fitA/fitB
        proba = 1 - f**2
        proba = proba / (1 - f**( 2*N ) )

        if ( is.nan(proba) )
        {
          proba = 0
        }
      }

      if (runif(1) < proba)
      {
        fitA = fitB
        designA = designB

        if ( fitBest < 1/d )
        {
          fitBest = 1/d
          best = designA

          if ( showProcess == TRUE )
          {
            # iteration and Dcriteria
            print( paste0('Iteration = ',iteration))
            print( paste0('Criterion = ',1/d))
          }
        }
      }
    }

    # Purge
    if ( iteration%%purgeIteration == 0 )
    {
      d = - log10( fitA ) / theta
      theta = -log10( fitBase ) / d
      fitA = fitBase
    }

  } # end iteration

  # evaluate the optimal design
  evaluationOptimalDesign = Evaluation( name = "",
                                        modelEquations = prop( optimizationObject, "modelEquations" ),
                                        modelParameters = prop( optimizationObject, "modelParameters" ),
                                        modelError = prop( optimizationObject, "modelError" ),
                                        designs = list( best ),
                                        fimType = prop( optimizationObject, "fimType" ),
                                        outputs = prop( optimizationObject, "outputs" ),
                                        odeSolverParameters = prop( optimizationObject, "odeSolverParameters" ) )

  evaluationOptimalDesign = run( evaluationOptimalDesign )

  # evaluate the initial design
  evaluationInitialDesign = Evaluation( name = "",
                                        modelEquations = prop( optimizationObject, "modelEquations" ),
                                        modelParameters = prop( optimizationObject, "modelParameters" ),
                                        modelError = prop( optimizationObject, "modelError" ),
                                        designs = list( initialDesign ),
                                        fimType = prop( optimizationObject, "fimType" ),
                                        outputs = prop( optimizationObject, "outputs" ),
                                        odeSolverParameters = prop( optimizationObject, "odeSolverParameters" ) )

  evaluationInitialDesign = run( evaluationInitialDesign )

  # set the results in evaluation
  prop( optimizationObject, "optimisationDesign" ) = list( evaluationInitialDesign = evaluationInitialDesign, evaluationOptimalDesign = evaluationOptimalDesign )
  prop( optimizationObject, "optimisationAlgorithmOutputs" ) = list( "optimizationAlgorithm" = optimizationAlgorithm, "optimalArms" = armsList )
  return( optimizationObject )
}

# ==============================================================================
#' constraintsTableForReport: table of the PGBOAlgorithm constraints for the report.
#' @name constraintsTableForReport
#' @param optimizationAlgorithm A object \code{PGBOAlgorithm}.
#' @param arms List of the arms.
#' @return The table for the constraints in the arms.
#' @export
# ==============================================================================

method( constraintsTableForReport, PGBOAlgorithm ) = function( optimizationAlgorithm, arms  )
{
  armsConstraints = map( pluck( arms, 1 ) , ~ getArmConstraints( .x, optimizationAlgorithm ) )
  armsConstraints = map_dfr( armsConstraints, ~ map_df(.x, ~ as.data.frame(.x, stringsAsFactors = FALSE)))
  colnames( armsConstraints ) = c( "Arms name" , "Number of subjects", "Outcome", "Initial samplings", "Samplings windows", "Number of times by windows","Min sampling" )
  armsConstraintsTable = kbl( armsConstraints, align = c( "l","c","c","c","c","c","c") ) %>% kable_styling( bootstrap_options = c( "hover" ), full_width = FALSE, position = "center", font_size = 13 )
  return( armsConstraintsTable )
}













