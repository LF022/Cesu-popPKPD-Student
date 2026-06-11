#' Design
#' @description
#' The class \code{Design} represents and stores information for the Design.
#' @title Design
#' @param name A string giving the name of the design.
#' @param size A integer giving the size of the design.
#' @param arms A list giving the arms of the design.
#' @param numberOfArms A integer giving the number of arms.
#' @param evaluationArms A list giving the valuation of the arms of the design.
#' @param fim A object \code{Fim} giving the Fim of the design.
#' @include Fim.R
#' @export

Design = new_class("Design", package = "PFIM",
                   properties = list(
                     name = new_property(class_character, default = character(0)),
                     size = new_property(class_double, default = 0.0),
                     arms = new_property(class_list, default = list()),
                     evaluationArms = new_property(class_list, default = list()),
                     numberOfArms = new_property(class_double, default = 0.0),
                     fim = new_property(Fim, default = NULL)
                   ))


evaluateDesign = new_generic( "evaluateDesign", c( "design" ) )
generateDosesCombination = new_generic( "generateDoseCombination", c( "design" ) )
generateSamplingTimesCombination = new_generic( "generateSamplingTimesCombination", c( "design" ) )
checkValiditySamplingConstraint = new_generic( "checkValiditySamplingConstraint", c( "design" ) )
setSamplingConstraintForOptimization = new_generic( "setSamplingConstraintForOptimization", c( "design" ) )

# ==============================================================================
#' evaluateDesign: evaluation of a design.
#' @name evaluateDesign
#' @param design An object \code{Design} giving the design.
#' @param model An object \code{Model} giving the model.
#' @param fim An object \code{Fim} giving the Fim.
#' @return The object \code{Design} with its evaluation results.
#' @export
# ==============================================================================

method( evaluateDesign, Design ) = function( design, model, fim ) {

  # evaluate the arms of each design
  arms = prop( design, "arms")
  prop( design, "evaluationArms") = arms %>% map( ~ evaluateArm( .x , model, fim ) )

  evaluationDesign = pluck( prop( design, "evaluationArms"),1)
  evaluationFim = prop( evaluationDesign, "evaluationFim" )

  prop( fim, "fisherMatrix" ) = prop( design, "evaluationArms" ) %>% map( ~ prop( prop( .x, "evaluationFim" ), "fisherMatrix" ) ) %>% reduce(`+`)
  prop( fim, "shrinkage" ) = prop( evaluationFim, "shrinkage" )
  prop( design, "fim" ) = fim

  return( design )
}

# ==============================================================================
#' generateDosesCombination: generate the combination for the doses.
#' @name generateDosesCombination
#' @param design An object \code{Design} giving the design.
#' @return dosesForFIMs, numberOfDoses used in the design optimization.
#' @export
# ==============================================================================

method( generateDosesCombination, Design ) = function( design ) {

  arms = prop( design, "arms" )
  armNames = map_chr( arms, ~ prop( .x, "name" ) )
  outcomes = map( arms, ~ map_chr(prop( .x, "administrationsConstraints" ), ~ prop( .x,"outcome" ) ) )

  # combination of the doses
  dosesForFIMsTmp = map( arms, ~ {
    administrationsConstraints = prop( .x, "administrationsConstraints" )
    armName = prop( .x, "name" )
    doses = map( administrationsConstraints, ~ {
      doses = prop( .x, "doses" )
    })
  }) %>% flatten() %>% expand.grid()

  # assign arm and responses names to each combination
  dosesForFIMs = list()
  iter = 1
  for ( arm in arms ) {
    administrations = prop( arm, "administrations")
    armName = prop( arm, "name")
    for ( administration in administrations ) {
      outcome = prop( administration, "outcome")
      dosesForFIMs[[armName]][[outcome]] = unlist(dosesForFIMsTmp[,iter])
      iter = iter + 1 } }
  return( c( dosesForFIMs, numberOfDoses = dim(dosesForFIMsTmp)[1] ) )
}

# ==============================================================================
#' generateSamplingTimesCombination: generate the combination for the samplings.
#' @name generateSamplingTimesCombination
#' @param design An object \code{Design} giving the design.
#' @return samplingTimesCombinations used in the design optimization.
#' @export
# ==============================================================================

method( generateSamplingTimesCombination, Design ) = function( design ) {

  arms = prop( design, "arms" )
  armNames = map_chr( arms, ~ prop( .x, "name" ) )

  armResults = map( arms, function( arm ) {

    # Extract the relevant properties for each arm
    armName = prop( arm, "name" )
    samplingTimesConstraints = prop( arm, "samplingTimesConstraints" )
    samplingTimes = prop( arm, "samplingTimes")
    outcomeNames = map_chr( samplingTimesConstraints, ~ prop( .x, "outcome" ) )

    # Generate all combinations of sampling times for the arm
    samplingTimesCombinations = map( samplingTimesConstraints, function( samplingTimesConstraint ) {
      initialSamplings = prop( samplingTimesConstraint, "initialSamplings" )
      fixedTimes = prop( samplingTimesConstraint, "fixedTimes" )
      numberOfSamplingsOptimisable = prop( samplingTimesConstraint, "numberOfsamplingsOptimisable" )
      availableSamplings = setdiff( initialSamplings, fixedTimes )
      combinations = combn( availableSamplings, numberOfSamplingsOptimisable - length( fixedTimes ), simplify = FALSE )

      # Generate all combinations of sampling times
      map( combinations, ~ c( fixedTimes, .x ) )
    })

    # Flatten the list and convert to a data frame
    samplingTimesCombinations = expand.grid( samplingTimesCombinations )
    colnames( samplingTimesCombinations ) = outcomeNames

    # Convert to a list of named lists
    samplingTimesCombinations = pmap( samplingTimesCombinations, ~ list( ... ) )

    # Map the sampling times combinations to the sampling times for each outcome
    samplingsForFIM = map( samplingTimesCombinations, function( samplingTimeCombination ) {
      map2( samplingTimes, outcomeNames, function( samplingTime, outcomeName ) {
        prop( samplingTime, "samplings" ) = sort( samplingTimeCombination[[outcomeName]] )
        samplingTime
      })
    })
  })
  set_names( armResults, armNames )
}

# ==============================================================================
#' checkValiditySamplingConstraint: check if the constraints used for the design optimization are valid.
#' @name checkValiditySamplingConstraint
#' @param design An object \code{Design} giving the design.
#' @return A boolean TRUE / FALSE, if FALSE it also gives an error message.
#' @export
# ==============================================================================

method( checkValiditySamplingConstraint, Design ) = function( design ) {

  arms = prop( design, "arms" )

  walk ( arms, function( arm ) {
    armName = prop( arm, "name" )

    # get the outcomes
    samplingTimesConstraints = prop( arm, "samplingTimesConstraints" )
    outcomes = map( samplingTimesConstraints, ~ prop( .x, "outcome" ) ) %>% unlist()

    # get samplings window constraints
    samplingsWindow = map( samplingTimesConstraints, ~ prop( .x, "samplingsWindows" ) ) %>% setNames( outcomes )

    # get numberOfTimesByWindows constraints
    numberOfTimesByWindows = map( samplingTimesConstraints, ~ prop( .x, "numberOfTimesByWindows" ) ) %>% setNames( outcomes )

    # get minimal time step for each windows
    minSampling = map( samplingTimesConstraints, ~ prop( .x, "minSampling" ) ) %>% setNames( outcomes )

    inputRandomSpaced = list()
    samplingTimesArms = list()

    walk ( outcomes, function( outcome ) {
      intervalsConstraints = list()

      # get samplingTimes and samplings
      samplingTimes = prop( arm, "samplingTimes")
      samplings = map( samplingTimes, ~ if ( prop( .x, "outcome") == outcome ) prop(.x ,"samplings" ) ) %>% compact() %>% unlist()

      minSamplingAndNumberOfTimesByWindows = as.data.frame( list( minSampling[[outcome]], numberOfTimesByWindows[[outcome]] ) )
      tmp = t( as.data.frame(samplingsWindow[[outcome]] ) )
      inputRandomSpaced[[outcome]] = as.data.frame( do.call( "cbind", list( tmp, minSamplingAndNumberOfTimesByWindows ) ) )

      colnames( inputRandomSpaced[[outcome]] ) = c("min","max","delta","n")
      rownames( inputRandomSpaced[[outcome]] ) = NULL

      if ( sum( numberOfTimesByWindows[[outcome]] ) != length( samplings ) ) {
        print ( " ==================================================================================================== ")
        print( paste0( " The sampling times constraint is not possible for arm ", armName, " and outcome ", outcome ) )
        print ( " ==================================================================================================== ")
        stop() }

      walk( seq_len( length( inputRandomSpaced[[outcome]]$n)), function( iter ) {
        min = inputRandomSpaced[[outcome]]$min[iter]
        max = inputRandomSpaced[[outcome]]$max[iter]
        delta = inputRandomSpaced[[outcome]]$delta[iter]
        n = inputRandomSpaced[[outcome]]$n[iter]

        distance = max-min-(n-1)*delta

        if ( distance < 0 ) {
          print ( " ==================================================================================================== ")
          print( paste0( " The sampling times constraint is not possible for arm ", armName, " and outcome ", outcome ) )
          print ( " ==================================================================================================== ")
          stop() }
      })
    })
  })
}

# ==============================================================================
#' setSamplingConstraintForOptimization: set the sampling time constraints for an arm for the design optimization.
#' @name setSamplingConstraintForOptimization
#' @param design An object \code{Design} giving the design.
#' @return The arm with the sampling time constraint for the design optimization.
#' @export
# ==============================================================================

method( setSamplingConstraintForOptimization, Design ) = function( design ) {

  arms = prop( design, "arms" )

  arms = map ( arms, function( arm )
  {
    # get the outcomes in the sampling times
    samplingTimes = prop( arm, "samplingTimes" )
    outcomes = map_chr( samplingTimes, ~ prop( .x, "outcome" ) )

    # set the sampling time constraints for missing outcomes ie from its sampling times
    samplingTimesConstraints = prop( arm, "samplingTimesConstraints" )
    outcomesSamplingTimesConstraints = map_chr( samplingTimesConstraints, ~ prop( .x, "outcome" ) )
    outcomesSamplingNotInTimesConstraints = outcomes[!outcomes %in% outcomesSamplingTimesConstraints]

    if ( length( outcomesSamplingNotInTimesConstraints ) !=0 )
    {
      samplingTimesConstraints = map( outcomesSamplingNotInTimesConstraints, function( outcomeSamplingNotInTimesConstraints )
      {
        samplings = map( samplingTimes, ~ if ( prop( .x, "outcome") == outcomeSamplingNotInTimesConstraints ) prop(.x ,"samplings" ) ) %>% compact() %>% unlist()

        newSamplingTimeConstraints  = SamplingTimeConstraints( outcome = outcomeSamplingNotInTimesConstraints,
                                                               initialSamplings = samplings,
                                                               samplingsWindows = list( c( min( samplings ), max( samplings ) ) ),
                                                               numberOfTimesByWindows = length( samplings ),
                                                               minSampling = 0 )

        samplingTimesConstraints = append( samplingTimesConstraints, newSamplingTimeConstraints )
      }) %>% reduce(c(.))
    }
    prop( arm, "samplingTimesConstraints" ) = samplingTimesConstraints
    return(arm)
  })

  # set the design with arms
  prop( design, "arms" ) = arms
  return( design )
}






