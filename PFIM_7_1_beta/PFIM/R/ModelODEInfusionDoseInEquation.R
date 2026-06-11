#' @description The class \code{ModelODEInfusionDoseInEquation} is used to defined a ModelODEInfusionDoseInEquation
#' @title ModelODEInfusionDoseInEquation
#' @inheritParams ModelODEInfusion
#' @param modelODE An object \code{modelODE}.
#' @param wrapperModelInfusion Wrapper for solver.
#' @param solverInputs A list giving the solver inputs.
#' @include ModelODEInfusion.R
#' @export

ModelODEInfusionDoseInEquation = new_class( "ModelODEInfusionDoseInEquation",
                                            package = "PFIM",
                                            parent = ModelODEInfusion,

                                            properties = list(
                                              modelODE = new_property(class_function, default = NULL ),
                                              wrapperModelInfusion = new_property(class_list, default = list()),
                                              solverInputs = new_property(class_list, default = list())
                                            ))

# ==============================================================================
#' defineModelWrapper: define the model wrapper for the ode solver
#' @name defineModelWrapper
#' @param model An object of class \code{ModelODEInfusionDoseInEquation} that defines the model.
#' @param evaluation An object of class Evaluation that defines the evaluation
#' @return The model with updated slots.
# ==============================================================================

method( defineModelWrapper, ModelODEInfusionDoseInEquation ) = function( model, evaluation ) {

  # outcomes with administration
  outcomesWithAdministration = evaluation %>%
    pluck( "designs" ) %>%
    map( ~ pluck( .x, "arms" ) ) %>%
    unlist() %>%
    map( ~ pluck( .x, "administrations" ) ) %>%
    unlist()%>%
    map( ~ pluck( .x, "outcome" ) ) %>%
    unlist()

  prop( model, "outcomesWithAdministration") = outcomesWithAdministration
  prop( model, "wrapperModelInfusion" ) = prop( evaluation, "modelEquations" )
  outputs = prop( evaluation, "outputs")
  prop( model, "outputFormula") = outputs
  prop( model, "outputNames") = names( outputs )

  return( model )
}

# ==============================================================================
#' defineModelAdministration: define the administration
#' @name defineModelAdministration
#' @param model An object of class \code{ModelODEInfusionDoseInEquation} that defines the model.
#' @param arm An object of class \code{Arm} that defines the arm.
#' @return The model with updated slots.
#' @export
# ==============================================================================

method( defineModelAdministration, ModelODEInfusionDoseInEquation ) = function( model, arm ) {

  # model wrapper
  wrapperModelInfusion = prop( model, "wrapperModelInfusion" )
  wrapperModelDuringInfusion = wrapperModelInfusion$duringInfusion
  wrapperModelAfterInfusion = wrapperModelInfusion$afterInfusion

  # variable derivative names
  variableDerivativeNames = names( wrapperModelDuringInfusion )

  #model parameters
  parameters = prop( model, "modelParameters" )
  parameterNames = map_chr( parameters, "name" )

  # administrations and outcome
  outcomesWithAdministration =  prop( model, "outcomesWithAdministration" )

  # number of equation for wrapper
  numberOfEquationsWithAdmin = length( outcomesWithAdministration )
  numberOfEquations = length( wrapperModelDuringInfusion  )

  # sampling times
  samplingTimes = prop( arm, "samplingTimes" )

  # define the sampling for all response
  samplings = map( samplingTimes, ~ prop( .x, "samplings" ) ) %>% unlist() %>% sort() %>% unique()
  samplings = unique( c( 0, samplings ) )

  # max values of the sampling times
  maxSampling = map_dbl( samplingTimes, ~ max( prop( .x ,"samplings" ) ) ) %>% max()
  # model outputs
  outputNames = prop( model, "outputNames" )
  outputFormula = prop( model, "outputFormula" )
  outputFormula = map( outputFormula, ~ parse( text=.x ) )
  outcomesWithAdministration =  prop( model, "outcomesWithAdministration" )

  # administration time
  administrationTime = list()
  administrations = prop( arm, "administrations" )
  solverInputs = map( administrations, ~ {

    outcome = prop( .x, "outcome" )
    timeDose = prop( .x, "timeDose" )
    tau = prop( .x, "tau" )
    dose = prop( .x, "dose" )
    Tinf = prop( .x, "Tinf")

    if ( tau != 0 ) {
      timeDose = seq( 0, maxSampling, tau )
      administrationTime = cbind( timeDose[ -length( timeDose ) ], timeDose[-1] )
      dose = rep( dose, length( timeDose ) )
      Tinf = rep( Tinf, length( timeDose ) )
    }

    administrationTime = cbind( timeDose, timeDose + Tinf ) %>% unname()

    setNames( list( list( administrationTime = administrationTime, dose = dose, Tinf = Tinf ) ), outcome )
  })  %>% flatten()


  # evaluate the initial conditions
  initialConditions = evaluateInitialConditions( model, arm )

  # Assign the values to variables in the current environment
  # Assign the values to the parameters in the current environment
  mu = set_names(
    map(parameters, ~ .x@distribution@mu),
    map(parameters, ~ .x@name)
  )

  list2env( mu, envir = environment() )

  # arguments for function evaluation model
  variableNames = names( initialConditions )
  doseNames = paste( "dose_", outcomesWithAdministration, sep = "" )
  tinfNames = paste( "Tinf_", outcomesWithAdministration, sep = "" )

  functionArguments = c( doseNames, tinfNames, parameterNames, variableNames )
  solverInputs$functionArguments = unique( functionArguments )
  solverInputs$functionArgumentsSymbols = map( functionArguments, ~ as.symbol(.x) )

  # outcome without administration
  outcomesWithoutAdministration = setdiff( variableNames, outcomesWithAdministration )

  # function to evaluate the model: create the wrapper with equations during/after infusion
  # define the equations without administration
  equationsWithAdministration = keep( wrapperModelDuringInfusion, ~ str_detect( ., "dose_") )
  equationsWithoutAdministration = keep( wrapperModelDuringInfusion, ~ !str_detect( ., "dose_") )

  modelODEInfusion = function( samplingTimes, initialConditions, solverInputs )
  {
    with( c( samplingTimes, initialConditions, solverInputs ),{

      # define the equations
      equation = list()
      equations = map2( outcomesWithAdministration, seq_along( outcomesWithAdministration ), function( outcomeWithAdministration,iter ) {

        equationArguments = list()

        # administration and index for infusion during/after
        administrationTime = solverInputs[[outcomeWithAdministration]]$administrationTime
        indexTime = which( samplingTimes >= administrationTime[, 1] & samplingTimes < administrationTime[, 2] )

        if ( length( indexTime ) != 0 )
        {
          # equations during infusion
          equation[[iter]] = wrapperModelDuringInfusion[[iter]]
          equationArguments[[iter]] = unique( c( doseNames, tinfNames, parameterNames, variableNames ) )
        }else
        {
          # equations after infusion
          equation[[iter]] = wrapperModelAfterInfusion[[iter]]
          equationArguments[[iter]] = c( parameterNames, variableNames )
        }
        return( list( equation = equation, equationArguments = equationArguments )  )
      })

      # model equation
      equation = map( equations, ~.x$equation ) %>% unlist()
      names( equation ) = names(equationsWithAdministration)
      equation = c( equation, equationsWithoutAdministration )

      # argument of model equation
      equationArguments = map( equations, ~.x$equationArguments ) %>% unlist() %>% unique()
      equationArgumentsSymbols = map( equationArguments, ~ as.symbol(.x) )

      # # create model wrapper
      equationsBody = map_chr( names( equation ), ~ sprintf( "%s = %s", .x, equation[[.x]] ) )
      functionBody = paste( equationsBody, collapse = "\n" )
      functionBody = sprintf( "%s\nreturn(list(c(%s)))", functionBody, paste( variableDerivativeNames, collapse = ", " ) )
      functionDefinition = sprintf( "function(%s) { %s }", paste( equationArguments, collapse = ", " ), functionBody )
      wrapper = eval( parse( text = functionDefinition ) )

      # evaluate wrapper
      # FIX: Initialize dose and Tinf with FIRST actual values from solverInputs
      # This avoids singularities and ensures variables exist for do.call()
      for ( outcome in outcomesWithAdministration ) {
        dose_val = solverInputs[[outcome]]$dose[1]
        Tinf_val = solverInputs[[outcome]]$Tinf[1]
        assign( paste0("dose_", outcome), dose_val )
        assign( paste0("Tinf_", outcome), Tinf_val )
      }

      for( outcomeWithAdministration in outcomesWithAdministration )
      {
        administrationTime = solverInputs[[outcomeWithAdministration]]$administrationTime
        indexTime = which( samplingTimes >= administrationTime[, 1] & samplingTimes < administrationTime[, 2] )

        if ( length( indexTime ) != 0 )
        {
          # equations during infusion
          doseNames = paste( "dose_", outcomeWithAdministration, sep = "" )
          tinfNames = paste( "Tinf_", outcomeWithAdministration, sep = "" )

          dose = solverInputs[[outcomeWithAdministration]]$dose
          Tinf = solverInputs[[outcomeWithAdministration]]$Tinf

          assign( doseNames, dose[indexTime] )
          assign( tinfNames, Tinf[indexTime] )
        }
      }
      evaluationModel = do.call( wrapper, setNames( equationArgumentsSymbols, equationArguments ) )
      evaluationOutputs  = map( outputFormula, ~ eval( .x ) )

      return( c( evaluationModel, evaluationOutputs ) )
    })
  }

  prop( model, "initialConditions" ) = initialConditions
  prop( model, "samplings" ) = samplings
  prop( model, "modelODE" ) = modelODEInfusion
  prop( model, "solverInputs" ) = solverInputs

  return( model )
}

# ==============================================================================
#' evaluateModel
#' @name evaluateModel
#' @param arm A object of class \code{Arm} giving the arm.
#' @param model A object of class \code{ModelODEInfusionDoseInEquation} giving the model.
#' @return A data frame giving the output of the model evaluation.
#' @export
# ==============================================================================

method( evaluateModel, ModelODEInfusionDoseInEquation ) = function( model, arm ) {

  initialConditions = prop( model, "initialConditions" )
  samplings = prop( model, "samplings" )
  modelODE = prop( model, "modelODE" )
  solverInputs = prop( model, "solverInputs" )
  odeSolverParameters = prop( model, "odeSolverParameters" )
  atol = odeSolverParameters$atol
  rtol = odeSolverParameters$rtol
  samplingTimes = prop( arm, "samplingTimes" )
  outputNames = prop( model, "outputNames" )

  # model evaluation
  evaluationModelTmp = ode( initialConditions, samplings, modelODE, solverInputs, hmax = 0.0, atol = atol, rtol = rtol )
  evaluationModelTmp = evaluationModelTmp %>% as.data.frame()

  # filter sampling time
  samplings = map( samplingTimes, ~ prop( .x, "samplings" ) ) %>% set_names( outputNames )

  evaluationModel = list()
  for ( outputName in outputNames )
  {
    time = evaluationModelTmp$time %in% samplings[[outputName]]
    evaluationModel[[outputName]] = evaluationModelTmp[ time , c( "time", outputName ) ]
  }

  return( evaluationModel )
}

# ==============================================================================
#' definePKModel: define PK model ode bolus
#' @name definePKModel
#' @param pkModel An object of class \code{ModelODEInfusionDoseInEquation} that defines the PK model.
#' @param pfimproject An object of class \code{PFIMProject} that defines the pfimproject.
#' @export
# ==============================================================================

method( definePKModel, list( ModelODEInfusionDoseInEquation, PFIMProject ) ) = function( pkModel, pfimproject ) {

  # get the initial conditions to get variable names
  designs = prop( pfimproject, "designs" )
  variablesNames = designs %>% map(~ map( prop(.x,"arms"), ~ prop(.x,"initialConditions"))) %>% unlist() %>% names() %>% unique()
  variablesNamesToChange =  c("C1", "C2")

  pkModelEquations = prop( pkModel, "modelEquations")
  pkModelEquations$duringInfusion = pkModelEquations$duringInfusion %>% imap(~reduce2(variablesNamesToChange, variablesNames, replaceVariablesLibraryOfModels, .init = .x))
  pkModelEquations$afterInfusion = pkModelEquations$afterInfusion %>% imap(~reduce2(variablesNamesToChange, variablesNames, replaceVariablesLibraryOfModels, .init = .x))

  pkModelEquations$duringInfusion =  pkModelEquations$duringInfusion %>% set_names( paste0( "Deriv_",variablesNames))
  pkModelEquations$afterInfusion =  pkModelEquations$afterInfusion %>% set_names( paste0( "Deriv_",variablesNames))

  return( pkModelEquations )
}

# ==============================================================================
#' definePKPDModel: define a PKPD model from library of model
#' @name definePKPDModel
#' @param pkModel An object of class \code{ModelODEInfusionDoseInEquation} that defines the PK model.
#' @param pkModel An object of class \code{ModelODE} that defines the PD model.
#' @param pfimproject An object of class \code{PFIMProject} that defines the pfimproject.
#' @export
# ==============================================================================

method( definePKPDModel, list( ModelODEInfusionDoseInEquation, ModelODE, PFIMProject ) ) = function( pkModel, pdModel, pfimproject ) {

  pkModelEquations = prop( pkModel, "modelEquations")
  pdModelEquations = prop( pdModel, "modelEquations")
  equations = c( pkModelEquations, pdModelEquations )
  pdModelEquations = str_replace_all( pdModelEquations, "E", paste0( "C", length( equations ) ) )
  pdModelEquations = str_replace_all( pdModelEquations, "RespPK", "C1" )
  names( pdModelEquations ) = paste0("Deriv_C",length( equations ) )

  equations = list( "duringInfusion" = c( pkModelEquations$duringInfusion, pdModelEquations ),
                    "afterInfusion" =  c( pkModelEquations$afterInfusion, pdModelEquations ) )

  return( equations )
}
