#' @description The class \code{ModelODEDoseNotInEquations} is used to defined a ModelODEDoseNotInEquations
#' @title ModelODEDoseNotInEquations
#' @inheritParams ModelODE
#' @param modelODE An object \code{modelODE}.
#' @param doseEvent A dataframge given the doseEvent for the ode solver.
#' @param solverInputs A list giving the solver inputs.
#' @include Model.R
#' @export

ModelODEDoseNotInEquations = new_class( "ModelODEDoseNotInEquations",
                                        package = "PFIM",
                                        parent = ModelODE,
                                        properties = list(
                                          modelODE = new_property(class_function, default = NULL ),
                                          doseEvent = new_property(class_list, default = list()),
                                          solverInputs = new_property(class_list, default = list())
                                        ))

# ==============================================================================
#' defineModelWrapper: define the model wrapper for the ode solver
#' @name defineModelWrapper
#' @param model An object of class \code{ModelODEDoseNotInEquations} that defines the model.
#' @param evaluation An object of class Evaluation that defines the evaluation
#' @return The model with the updated slots.
# ==============================================================================

method( defineModelWrapper, ModelODEDoseNotInEquations ) = function( model, evaluation ) {

  # names of the equations and the variables
  equations = prop( evaluation, "modelEquations" )
  variableNames = str_remove( names( equations ), "Deriv_" )
  variableNamesDerivatives = paste( names( equations ), collapse = ", " )

  # outcomes with administration
  outcomes = evaluation %>%
    pluck( "designs" ) %>%
    map( ~ pluck( .x, "arms" ) ) %>%
    unlist() %>%
    map( ~ pluck( .x, "administrations" ) ) %>%
    unlist()%>%
    map( ~ pluck( .x, "outcome" ) ) %>%
    unlist()

  # arguments for the function
  parameters = prop( evaluation, "modelParameters" )
  parameterNames = map_chr( parameters, "name" )

  functionArguments = c( parameterNames, variableNames, "t" )
  functionArguments = unique( functionArguments )
  functionArgumentsSymbol = map( functionArguments, ~ as.symbol(.x) )

  # create body function
  equationsBody = map_chr( names( equations ), ~ sprintf( "%s = %s", .x, equations[[.x]] ) )

  functionBody = paste( equationsBody, collapse = "\n" )
  functionBody = sprintf( paste( "%s\nreturn(list(c(", variableNamesDerivatives, ")))", collapse = ", " ), functionBody )
  functionDefinition = sprintf( "function(%s) { %s }", paste( functionArguments, collapse = ", " ), functionBody )

  outputs = prop( evaluation, "outputs")
  prop( model, "outputFormula") = outputs
  prop( model, "outputNames") = names( outputs )
  prop( model, "wrapper" ) = eval( parse( text = functionDefinition ) )
  prop( model, "functionArguments" ) = functionArguments
  prop( model, "functionArgumentsSymbol" ) = functionArgumentsSymbol

  return( model )
}

# ==============================================================================
#' defineModelAdministration: define the administration
#' @name defineModelAdministration
#' @param model An object of class \code{ModelODEDoseNotInEquations} that defines the model.
#' @param arm An object of class \code{Arm} that defines the arm.
#' @return The model with samplings, solverInputs
#' @export
# ==============================================================================

method( defineModelAdministration, ModelODEDoseNotInEquations ) = function( model, arm ) {

  # administrations and samplings
  administrations = prop( arm, "administrations" )
  samplingTimes = prop( arm, "samplingTimes" )
  samplings = map( samplingTimes, ~ prop( .x, "samplings" ) ) %>% unlist()
  samplings = unique( c( 0.0, samplings ) )
  # model wrapper function
  wrapper = prop( model, "wrapper")
  # model parameters
  parameters = prop( model, "modelParameters" )
  # args for model evaluation
  functionArguments = prop( model, "functionArguments" )
  functionArgumentsSymbols = prop( model, "functionArgumentsSymbol" )
  # model outputs
  outputFormula = prop( model, "outputFormula" )
  outputFormula = map( outputFormula, ~ parse( text=.x ) )

  # dose event: variable as compartment
  doseEvent = map( administrations, ~ {

    outcome = prop( .x, "outcome" )
    timeDose = prop( .x, "timeDose" )
    tau = prop( .x, "tau" )
    dose = prop( .x, "dose" )

    if ( tau !=0 )
    {
      timeDose = seq( 0, max( samplings ), tau )
      dose = rep( dose, length( timeDose ) )
    }
    data.frame( var = rep( outcome, length( timeDose ) ), time = timeDose, value = dose, method = ifelse(timeDose == 0, "replace", "add") )
  }) %>%
    reduce( rbind ) %>%
    .[ order( .$time ), ]

  # evaluate the initial conditions
  initialConditions = evaluateInitialConditions( model, arm )

  # Assign the values to variables in the current environment
  mu = set_names(
    map(parameters, ~ .x@distribution@mu),
    map(parameters, ~ .x@name)
  )

  list2env( mu, envir = environment() )

  # initial conditions with variable admin
  initialConditionsAdmin = setNames( rep(0, length(unique(doseEvent$var[doseEvent$time == 0]))),
                                     unique(doseEvent$var[doseEvent$time == 0]))

  initialConditions = c(initialConditionsAdmin, initialConditions)
  initialConditions = initialConditions[!duplicated(names(initialConditions)) ]

  # function evaluation model
  modelODEDoseAsCmpt = function( samplingTimes, initialConditions, parameters )
  {
    with( as.list( c(  samplingTimes, initialConditions, parameters ) ),{

      # evaluate wrapper and  model outputs
      evaluationModel = do.call( wrapper, setNames( functionArgumentsSymbols, functionArguments ) )
      evaluationOutputs  = map( outputFormula, ~ eval( .x ) )

      return( c( evaluationModel , evaluationOutputs ) )
    })}

  # set the model
  prop( model, "initialConditions" ) = initialConditions
  prop( model, "samplings" ) = samplings
  prop( model, "modelODE" ) = modelODEDoseAsCmpt
  prop( model, "doseEvent" ) = doseEvent

  return( model )
}

# ==============================================================================
#' evaluateModel: evaluate the model
#' @name evaluateModel
#' @param model An object of class \code{ModelODEDoseNotInEquations} that defines the model.
#' @param arm An object of class \code{Arm} that defines the arm.
#' @return A list of dataframes that contains the results for the evaluation of the model.
#' @export
# ==============================================================================

method( evaluateModel, ModelODEDoseNotInEquations ) = function( model, arm ) {

  evaluateODEDoseNotInEquationsCore = function( model, arm ) {

    initialConditions   = prop( model, "initialConditions" )
    samplings           = prop( model, "samplings" )
    modelODE            = prop( model, "modelODE" )
    parameters          = NULL
    odeSolverParameters = prop( model, "odeSolverParameters" )
    atol                = odeSolverParameters$atol
    rtol                = odeSolverParameters$rtol
    doseEvent           = prop( model, "doseEvent" )
    outputNames         = prop( model, "outputNames" )
    samplingTimes       = prop( arm, "samplingTimes" )

    # Évaluation du modèle
    evaluationModelTmp = ode(
      initialConditions,
      samplings,
      modelODE,
      parameters,
      events = list( data = doseEvent ),
      atol = atol,
      rtol = rtol
    ) %>% data.frame()

    # Filtrer les temps d’échantillonnage par output
    samplingsByOutput = map( samplingTimes, ~ prop( .x, "samplings" ) ) %>%
      set_names( outputNames )

    evaluationModel = list()
    for ( outputName in outputNames ) {
      timeSel = evaluationModelTmp$time %in% samplingsByOutput[[outputName]]
      evaluationModel[[outputName]] = evaluationModelTmp[ timeSel, c( "time", outputName ) ]
    }

    return( evaluationModel )
  }

  # Multi-occasion random IOV without covariates must also use the complex
  # combination/occasion path.
  if ( usesCovariateOccasionStructure( model ) ) {
    return( evaluateModelWithCovariates( model, arm, evaluateODEDoseNotInEquationsCore ) )
  } else {
    return( evaluateODEDoseNotInEquationsCore( model, arm ) )
  }
}

# ==============================================================================
#' definePKModel: define a PK model from library of model
#' @name definePKModel
#' @param pkModel An object of class \code{ModelODEDoseNotInEquations} that defines the PK model.
#' @param pfimproject An object of class \code{PFIMProject} that defines the pfimproject.
#' @export
# ==============================================================================

method( definePKModel, list( ModelODEDoseNotInEquations, PFIMProject ) ) = function( pkModel, pfimproject ) {
  pkModelEquations = prop( pkModel, "modelEquations")
  return( pkModelEquations )
}











