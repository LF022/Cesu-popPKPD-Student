#' @description The class \code{ModelAnalytic} is used to defined an analytic model.
#' @title ModelAnalytic
#' @param wrapperModelAnalytic Wrapper for the ode solver.
#' @inheritParams Model
#' @param functionArgumentsModelAnalytic A list giving the functionArguments of the wrapper for the analytic model.
#' @param functionArgumentsSymbolModelAnalytic A list giving the functionArgumentsSymbol of the wrapper for the analytic model
#' @param solverInputs A list giving the solver inputs.
#' @include Model.R
#' @include ModelODE.R
#' @export

ModelAnalytic = new_class(
  "ModelAnalytic",
  package = "PFIM",
  parent = Model,

  properties = list(
    wrapperModelAnalytic = new_property(class_list, default = list()),
    functionArgumentsModelAnalytic = new_property(class_list, default = list()),
    functionArgumentsSymbolModelAnalytic = new_property(class_list, default = list()),
    solverInputs = new_property(class_list, default = list())
  ))

convertPKModelAnalyticToPKModelODE = new_generic( "convertPKModelAnalyticToPKModelODE", c( "pkModel" ) )

# ==============================================================================
#' defineModelWrapper: define the model wrapper for the ode solver
#' @name defineModelWrapper
#' @param model An object of class \code{ModelAnalytic} that defines the model.
#' @param evaluation An object of class Evaluation that defines the evaluation
#' @return The model with wrapperModelAnalytic, functionArgumentsModelAnalytic, functionArgumentsSymbolModelAnalytic, outputNames, outcomesWithAdministration
# ==============================================================================

method( defineModelWrapper, ModelAnalytic ) = function( model, evaluation ) {

  # outcomes with administration
  outcomesWithAdministration = evaluation %>%
    pluck( "designs" ) %>%
    map( ~ pluck( .x, "arms" ) ) %>%
    unlist() %>%
    map( ~ pluck( .x, "administrations" ) ) %>%
    unlist()%>%
    map( ~ pluck( .x, "outcome" ) ) %>%
    unlist() %>% unique()

  # arguments for the function
  parameters = prop( evaluation, "modelParameters" )
  parameterNames = map_chr( parameters, "name" )
  doseNames = paste( "dose_", outcomesWithAdministration, sep = "" )
  timeNames = paste( "t_", outcomesWithAdministration, sep = "" )

  # names of the equations with admin and no admin
  equations = prop( evaluation, "modelEquations" )
  equationsWithAdmin = equations[ names( equations ) %in% outcomesWithAdministration ]
  equationsWithNoAdmin = equations[ !( names( equations ) %in% outcomesWithAdministration ) ]

  # output
  outputs = names( equations )
  outputNames = unlist( outputs )

  # outputs with / without admin
  indexOutputNoAdmin = which( !( names( equations ) %in% outcomesWithAdministration ) )
  outputNoAdmin = outputNames[ indexOutputNoAdmin ] %>% unlist()

  # outputForEvaluation
  outputsForEvaluation = prop( evaluation, "outputs" )
  # pk model
  if ( length(outputsForEvaluation ) == 1 )
  {
    outputAdmin = unlist(outputsForEvaluation[1])
    outputNoAdmin = c()
    # pkpd model
  }else if ( length(outputsForEvaluation ) == 2 )
  {
    outputAdmin = unlist(outputsForEvaluation[1])
    outputNoAdmin = unlist(outputsForEvaluation[2])
  }

  # wrapper for function with outcome administration

  # args for function with admin
  functionArgumentsWithAdmin = unique( c( doseNames, parameterNames, timeNames ) )
  functionArgumentsSymbolWithAdmin = map( functionArgumentsWithAdmin, ~ as.symbol(.x) )

  # create function with admin
  equationsBodyWithAdmin = map_chr( names( equationsWithAdmin ), ~ sprintf( "%s = %s", .x, equationsWithAdmin[[.x]] ) )
  equationsBodyWithAdmin = map2_chr( equationsBodyWithAdmin, timeNames, ~ str_replace_all( .x, "\\bt\\b", .y ) )

  functionBodyWithAdmin = paste( equationsBodyWithAdmin, collapse = "\n" )
  functionBodyWithAdmin = sprintf( paste( "%s\nreturn(list(c(", paste( outputAdmin, collapse = ", ") , ")))", collapse = ", " ), functionBodyWithAdmin )
  functionDefinitionWithAdmin = sprintf( "function(%s) { %s }", paste( functionArgumentsWithAdmin, collapse = ", " ), functionBodyWithAdmin )
  functionDefinitionWithAdmin = eval( parse( text = functionDefinitionWithAdmin ) )

  # wrapper for function outcome without administration

  # args for function without admin
  functionArgumentsWithNoAdmin = unique( c( outcomesWithAdministration, parameterNames, timeNames ) )
  functionArgumentsSymbolWithNoAdmin = map( functionArgumentsWithNoAdmin, ~ as.symbol(.x) )

  # create function without admin
  equationsBodyWithNoAdmin = map_chr( names( equationsWithNoAdmin ), ~ sprintf( "%s = %s", .x, equationsWithNoAdmin[[.x]] ) )
  equationsBodyWithNoAdmin = map2_chr( equationsBodyWithNoAdmin, timeNames, ~ str_replace( .x, "\\bt\\b", .y ) )
  functionBodyWithNoAdmin = paste( equationsBodyWithNoAdmin, collapse = "\n" )
  functionBodyWithNoAdmin = sprintf( paste( "%s\nreturn(list(c(", paste( outputNoAdmin, collapse = ", "), ")))", collapse = ", " ), functionBodyWithNoAdmin )
  functionDefinitionWithNoAdmin = sprintf( "function(%s) { %s }", paste( functionArgumentsWithNoAdmin, collapse = ", " ), functionBodyWithNoAdmin )
  functionDefinitionWithNoAdmin = eval( parse( text = functionDefinitionWithNoAdmin ) )

  prop( model, "wrapperModelAnalytic" ) = list( functionDefinitionWithAdmin = functionDefinitionWithAdmin,
                                                functionDefinitionWithNoAdmin = functionDefinitionWithNoAdmin )

  prop( model, "functionArgumentsModelAnalytic" ) = list( functionArgumentsWithAdmin = functionArgumentsWithAdmin,
                                                          functionArgumentsWithNoAdmin = functionArgumentsWithNoAdmin )

  prop( model, "functionArgumentsSymbolModelAnalytic" ) = list( functionArgumentsSymbolWithAdmin = functionArgumentsSymbolWithAdmin,
                                                                functionArgumentsSymbolWithNoAdmin = functionArgumentsSymbolWithNoAdmin )

  # define the model
  prop( model, "outputNames") = unlist( outputs )
  prop( model, "outcomesWithAdministration") = outcomesWithAdministration

  return( model )
}

# ==============================================================================
#' defineModelAdministration: define the administration
#' @name defineModelAdministration
#' @param model An object of class \code{ModelAnalytic} that defines the model.
#' @param arm An object of class \code{Arm} that defines the arm.
#' @return The model with samplings, solverInputs
#' @export
# ==============================================================================

method( defineModelAdministration, ModelAnalytic ) = function( model, arm ) {

  # administrations and outcome
  administrations = prop( arm, "administrations" )
  outcomesWithAdministration =  prop( model, "outcomesWithAdministration" )
  # sampling times
  samplingTimes = prop( arm, "samplingTimes" )
  # define the samplings for all response
  samplings = map( samplingTimes, ~ prop( .x, "samplings" ) ) %>% unlist() %>% sort() %>% unique()
  # model outputs
  outputNames = prop( model, "outputNames" )
  # define solverInputs
  solverInputs = map( administrations, function(  administration ) {

    timeDose = prop( administration, "timeDose" )
    tau = prop( administration, "tau" )
    dose = prop( administration, "dose" )
    maxSampling = max( samplings )

    if ( tau != 0 ) {
      timeDose = seq( 0, maxSampling, tau )
      dose = rep( dose, length( timeDose ) )
    }

    # define the time doses
    timeDose = timeDose %>%
      map( ~ ifelse( samplings - .x > 0, samplings - .x, samplings ) ) %>%
      reduce( cbind )

    indicesDoses = if ( is.null( dim( timeDose ) ) ) {
      # dose unique
      indicesDoses = 1
    } else {
      # multi dose
      indicesDoses = map_int( seq_len( dim( timeDose )[1] ), ~{
        length( unique( timeDose[.x, ] ) )
      })
    }
    list( data = data.frame( timeDose, indicesDoses ), dose = dose )
  }) %>% setNames( outcomesWithAdministration )

  prop( model, "samplings" ) = samplings
  prop( model, "solverInputs" ) = solverInputs

  return( model )
}

# ==============================================================================
#' evaluateAnalyticCore: core function to evaluate analytic model
#' @name evaluateAnalyticCore
#' @param model An object of class \code{ModelAnalytic} that defines the model.
#' @param arm An object of class \code{Arm} that defines the arm.
#' @return A list of dataframes containing evaluation results
#' @keywords internal
# ==============================================================================

evaluateAnalyticCore = function( model, arm ) {

  # Extract model components
  parameters = prop( model, "modelParameters" )
  outcomesWithAdministration = prop( model, "outcomesWithAdministration" )
  outputNames = prop( model, "outputNames" )
  samplings = prop( model, "samplings" )
  solverInputs = prop( model, "solverInputs" )

  # Get wrapper functions
  wrapperModelAnalytic = prop( model, "wrapperModelAnalytic" )
  functionDefinitionWithAdmin = wrapperModelAnalytic$functionDefinitionWithAdmin
  functionDefinitionWithNoAdmin = wrapperModelAnalytic$functionDefinitionWithNoAdmin

  # Get function arguments
  functionArguments = prop( model, "functionArgumentsModelAnalytic" )
  functionArgumentsWithAdmin = functionArguments$functionArgumentsWithAdmin
  functionArgumentsWithNoAdmin = functionArguments$functionArgumentsWithNoAdmin

  functionArgumentsSymbols = prop( model, "functionArgumentsSymbolModelAnalytic" )
  functionArgumentsSymbolWithAdmin = functionArgumentsSymbols$functionArgumentsSymbolWithAdmin
  functionArgumentsSymbolWithNoAdmin = functionArgumentsSymbols$functionArgumentsSymbolWithNoAdmin

  # CREATE AN EXPLICIT ENVIRONMENT for parameter values
  evalEnv = new.env(parent = environment())

  # Assign parameter values to the explicit environment
  mu = set_names(
    map( parameters, ~ .x@distribution@mu ),
    map_chr( parameters, "name" )
  )
  list2env( mu, envir = evalEnv )

  # Evaluate model at each sampling time
  evaluationModelTmp = imap_dfr( samplings, function( time, iterTime ) {

    # Evaluate each outcome with administration
    outcomes = map( outcomesWithAdministration, function( outcomeWithAdministration ) {

      data = solverInputs[[outcomeWithAdministration]]$data
      dose = solverInputs[[outcomeWithAdministration]]$dose

      indicesDoses = data$indicesDoses[iterTime]
      times = data[iterTime, 1:indicesDoses]
      doses = dose[1:indicesDoses]

      # Sum contributions from all dose administrations
      evaluationOutcomeWithAdmin = sum( map_dbl( seq_len( indicesDoses ), function( indiceDose ) {

        # Assign dose and time variables to the explicit environment
        assign( paste0( "t_", outcomeWithAdministration ), times[indiceDose], envir = evalEnv )
        assign( paste0( "dose_", outcomeWithAdministration ), doses[indiceDose], envir = evalEnv )

        # Call function in the explicit environment
        do.call(
          functionDefinitionWithAdmin,
          setNames( functionArgumentsSymbolWithAdmin, functionArgumentsWithAdmin ),
          envir = evalEnv
        ) %>% unlist()
      }))

      # Assign PK response value to the explicit environment
      assign( outcomeWithAdministration, evaluationOutcomeWithAdmin, envir = evalEnv )

      # Evaluate PD response if exists
      evaluationOutcomeWithNoAdmin = do.call(
        functionDefinitionWithNoAdmin,
        setNames( functionArgumentsSymbolWithNoAdmin, functionArgumentsWithNoAdmin ),
        envir = evalEnv
      ) %>% unlist()

      # Return results
      if ( is.null( evaluationOutcomeWithNoAdmin ) || length( evaluationOutcomeWithNoAdmin ) == 0 ) {
        data.frame( Admin = evaluationOutcomeWithAdmin )
      } else {
        data.frame( Admin = evaluationOutcomeWithAdmin, NoAdmin = evaluationOutcomeWithNoAdmin )
      }
    })

    # Combine outcomes and add time column
    outcomes_combined = reduce( outcomes, cbind )
    data.frame( time = time, outcomes_combined )
  })
  # Rename columns
  colnames( evaluationModelTmp ) = c( "time", outputNames )

  # Filter by sampling times for each output
  samplingTimes = prop( arm, "samplingTimes" )
  samplingsByOutput = map( samplingTimes, ~ prop( .x, "samplings" ) ) %>%
    setNames( outputNames )

  map( outputNames, function( outputName ) {
    timeFilter = evaluationModelTmp$time %in% samplingsByOutput[[outputName]]
    evaluationModelTmp[timeFilter, c( "time", outputName )]
  }) %>% setNames( outputNames )
}

# ==============================================================================
#' evaluateModel: evaluate the model
#' @name evaluateModel
#' @param model An object of class \code{ModelAnalytic} that defines the model.
#' @param arm An object of class \code{Arm} that defines the arm.
#' @return A list of dataframes that contains the results for the evaluation of the model.
#' @export
# ==============================================================================

method( evaluateModel, ModelAnalytic ) = function( model, arm ) {

  # Analytic models follow the same occasion-aware dispatch as ODE models.
  if ( usesCovariateOccasionStructure( model ) ) {
    return( evaluateModelWithCovariates( model, arm, evaluateAnalyticCore ) )
  } else {
    return( evaluateAnalyticCore( model, arm ) )
  }
}

# ==============================================================================
#' convertPKModelAnalyticToPKModelODE: conversion from analytic to ode
#' @name convertPKModelAnalyticToPKModelODE
#' @param pkModel An object of class \code{ModelAnalytic} that defines the model.
#' @export
# ==============================================================================

method( convertPKModelAnalyticToPKModelODE, ModelAnalytic ) = function( pkModel  ) {

  pkModelEquations = prop( pkModel, "modelEquations")
  dtEquationPKsubstitute = D( parse( text = pkModelEquations ), "t" )
  dtEquationPKsubstitute = str_c( deparse( dtEquationPKsubstitute ), collapse = "" )
  pkModelEquations =  pluck( pkModelEquations, 1 )

  if ( str_detect( pkModelEquations, "Cl" ) )
  {
    pkModelEquations = str_c( dtEquationPKsubstitute, "+(Cl/V)*", pkModelEquations, "- (Cl/V)*RespPK" )
  } else {
    pkModelEquations = str_c( dtEquationPKsubstitute, "+k*", pkModelEquations, "- k*RespPK" )
  }
  pkModelEquations = str_replace_all( pkModelEquations, " ", "" )
  pkModelEquations = paste( Simplify( pkModelEquations ) )

  return( pkModelEquations )
}

# ==============================================================================
#' definePKModel: define a PK model from library of model
#' @name definePKModel
#' @param pkModel An object of class \code{ModelAnalytic} that defines the PK model.
#' @param pfimproject An object of class \code{PFIMProject} that defines the pfimproject.
#' @export
# ==============================================================================

method( definePKModel, list( ModelAnalytic, PFIMProject ) ) = function( pkModel, pfimproject ) {
  pkModelEquations = prop( pkModel, "modelEquations")
  return( pkModelEquations )
}

# ==============================================================================
#' definePKPDModel:  define a PKPD model from library of model
#' @name definePKPDModel
#' @param pkModel An object of class \code{ModelAnalytic} that defines the PK model.
#' @param pkModel An object of class \code{ModelAnalytic} that defines the PD model.
#' @param pfimproject An object of class \code{PFIMProject} that defines the pfimproject.
#' @export
# ==============================================================================

method( definePKPDModel, list( ModelAnalytic, ModelAnalytic, PFIMProject ) ) = function( pkModel, pdModel, pfimproject ) {
  pkModelEquations = prop( pkModel, "modelEquations")
  pdModelEquations = prop( pdModel, "modelEquations")
  equations = c( pkModelEquations, pdModelEquations )
  return( equations )
}

# ==============================================================================
#' definePKPDModel:  define a PKPD model from library of model
#' @name definePKPDModel
#' @param pkModel An object of class \code{ModelAnalytic} that defines the PK model.
#' @param pkModel An object of class \code{ModelODE} that defines the PD model.
#' @param pfimproject An object of class \code{PFIMProject} that defines the pfimproject.
#' @export
# ==============================================================================

method( definePKPDModel, list( ModelAnalytic, ModelODE, PFIMProject ) ) = function( pkModel, pdModel, pfimproject ) {

  # PKPD model equations
  pkModelEquations = convertPKModelAnalyticToPKModelODE( pkModel )
  pdModelEquations = prop( pdModel, "modelEquations")
  equations = c( pkModelEquations, pdModelEquations )

  # get the initial conditions to get variable names
  designs = prop( pfimproject, "designs" )
  variablesNames = designs %>% map(~ map( prop(.x,"arms"), ~ prop(.x,"initialConditions"))) %>% unlist() %>% names() %>% unique()
  variablesNamesToChange =  c("RespPK", "E")

  # modify variable names in the model equations
  equations = equations %>% imap( ~ reduce2( variablesNamesToChange, variablesNames, replaceVariablesLibraryOfModels, .init = .x ) ) %>% set_names( paste0( "Deriv_", variablesNames ) )

  return( equations )
}
