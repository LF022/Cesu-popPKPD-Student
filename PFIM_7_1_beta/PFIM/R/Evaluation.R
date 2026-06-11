#' @description The class \code{Evaluation} represents and stores information for the evaluation of a design
#' @title Evaluation
#' @inheritParams PFIMProject
#' @param evaluationDesign A list giving the evaluation of the design.
#' @param modelCovariatesEquation modelCovariatesEquation
#' @param numberOfOccasions numberOfOccasions
#' @include PFIMProject.R
#' @include Model.R
#' @export

Evaluation = new_class("Evaluation", package = "PFIM", parent = PFIMProject,
                       properties = list(
                         evaluationDesign = new_property(class_list, default = list()),
                         name = new_property(class_character, default = character(0)),
                         modelParameters = new_property(class_list, default = list()),
                         modelCovariates = new_property(class_list, default = list()),
                         modelCovariatesEquation = new_property(class_character, default = character(0)),
                         modelEquations = new_property(class_list, default = list()),
                         modelFromLibrary  = new_property(class_list, default = list()),
                         modelError = new_property(class_list, default = list()),
                         designs = new_property(class_list, default = list()),
                         outputs = new_property(class_list, default = list()),
                         fimType = new_property(class_character, default = character(0)),
                         odeSolverParameters =  new_property(class_list, default = list())
                       ),
                       constructor = function(evaluationDesign = list(),
                                              name = character(0),
                                              modelParameters = list(),
                                              modelCovariates = list(),
                                              modelCovariatesEquation = character(0),
                                              modelEquations = list(),
                                              modelFromLibrary = list(),
                                              modelError = list(),
                                              designs = list(),
                                              outputs = list(),
                                              fimType = character(0),
                                              numberOfOccasions = 1,
                                              odeSolverParameters = list() ) {
                         new_object(
                           .parent = PFIMProject(),
                           evaluationDesign = evaluationDesign,
                           name = name,
                           modelParameters = modelParameters,
                           modelCovariates = modelCovariates,
                           modelCovariatesEquation = modelCovariatesEquation,
                           modelEquations = modelEquations,
                           modelFromLibrary = modelFromLibrary,
                           modelError = modelError,
                           designs = designs,
                           outputs = outputs,
                           fimType = fimType,
                           numberOfOccasions = numberOfOccasions,
                           odeSolverParameters = odeSolverParameters
                         )
                       })

getFim = new_generic( "getFim", c( "evaluation" ) )

# ==============================================================================
#' getListLastName: routine to get the names of last element of a nested list.
#' @name getListLastName
#' @param list The list to be used.
#' @return The names of last element.
#' @export
# ==============================================================================

# names of last element nested list
getListLastName = function( list ) {
  if ( is.list( list ) ) {
    result = map( list, getListLastName )
    result = unlist( result, recursive = FALSE )
    if ( length( result ) == 0 ) {
      return( names( list ) )
    } else {
      return( result )
    }
  }
}

# ==============================================================================
#' run: run the evaluation of a design.
#' @name run
#' @param pfimproject A object \code{PFIMProject} giving the Evaluation.
#' @return The object \code{Evaluation} giving the design evaluation.
#'
#' @examples
#' \dontrun{
#' # Create or load a PFIM project (e.g., 'myPFIMproject')
#' # This object acts as an Evaluation instance containing all necessary
#' # arguments, e.g. model equations, parameters, and the design to be evaluated.
#'
#' # Execute the design evaluation
#' # This calculates the Fisher Information Matrix (FIM) and related metrics.
#' evaluationResults = run(myPFIMproject)
#'
#' # Access and display results from the Evaluation object
#' # Use the show method to see the FIM, RSE, and other design properties.
#' show(evaluationResults)
#'
#' }
#' @export
# ==============================================================================

# run the model evaluation
method( run, Evaluation ) = function( pfimproject )
{
  # define the model equations if model from library of model
  modelFromLibraryOfModel = prop( pfimproject, "modelFromLibrary" )

  if ( length( modelFromLibraryOfModel ) != 0 ) {
    prop( pfimproject, "modelEquations" ) = defineModelEquationsFromLibraryOfModel( pfimproject ) }

  model = pfimproject |>
    defineModelType() |>
    finiteDifferenceHessian() |>
    (\(m) defineModelWrapper(m, pfimproject))()

  # Build the combination/occasion structure as soon as it is needed:
  # either because covariates exist, or because random IOV implies
  # several occasions with no covariate at all.
  if (usesCovariateOccasionStructure(model)) {
    model = defineCovariatesData(model)
  }

  # define the type of the Fim
  fim = defineFim( pfimproject )

  # evaluate the designs
  designs = prop( pfimproject, "designs" )
  evaluationDesign = map( designs, ~ evaluateDesign( .x, model, fim ) )
  prop( pfimproject, "evaluationDesign" ) = evaluationDesign

  # results of the evaluation
  designNumber = 1
  evaluationDesign = pluck( evaluationDesign, designNumber )
  prop( pfimproject, "fim" ) = prop( evaluationDesign, "fim" )

  return( pfimproject )
}

# ==============================================================================
#' getFim: get the Fisher matrix.
#' @name getFim
#' @param evaluation An object \code{Evaluation} giving the evaluation to be run.
#' @return The matrices fisherMatrix, fixedEffects, varianceEffects.
#' @export
# ==============================================================================

method( getFim, Evaluation ) = function( evaluation )
{
  fim = prop( evaluation, "fim" )
  fisherMatrix = prop( fim, "fisherMatrix" )
  fixedEffects = prop( fim, "fixedEffects" )
  varianceEffects = prop( fim, "varianceEffects" )

  return( list( fisherMatrix = fisherMatrix, fixedEffects = fixedEffects, varianceEffects = varianceEffects ) )
}

# ==============================================================================
#' getFisherMatrix: Display the Fisher Information Matrix components
#'
#' @description
#' Extracts and returns the partitioned components of the Fisher Information Matrix (FIM),
#' including the full matrix, fixed effects, and variance effects.
#'
#' @name getFisherMatrix
#' @param evaluation An object of class \code{Evaluation} containing the results of a design evaluation.
#' @return A \code{list} containing:
#' \itemize{
#'   \item \code{fisherMatrix}: The full Fisher Information Matrix.
#'   \item \code{fixedEffects}: The sub-matrix corresponding to fixed effects.
#'   \item \code{varianceEffects}: The sub-matrix corresponding to variance components (random effects).
#' }
#'
#' @examples
#' \dontrun{
#' # Assuming 'evaluationResults' is an object returned by the run() function
#' fimComponents = getFisherMatrix(evaluationResults)
#'
#' # Access the specific matrices
#' fullFim = fimComponents$fisherMatrix
#' fixedEff = fimComponents$fixedEffects
#' }
#' @export
# ==============================================================================

method( getFisherMatrix, Evaluation ) = function( pfimproject )
{
  fim = prop( pfimproject, "fim" )
  fim = setEvaluationFim( fim, pfimproject )
  fisherMatrix = prop( fim, "fisherMatrix" )
  fixedEffects = prop( fim, "fixedEffects" )
  varianceEffects = prop( fim, "varianceEffects" )

  return( list( fisherMatrix = fisherMatrix, fixedEffects = fixedEffects, varianceEffects = varianceEffects ) )
}

# ==============================================================================
#' show: Display the evaluation results in the R console
#'
#' @description
#' Provides a formatted summary of the design evaluation, including the Fisher
#' Information Matrix (FIM), standard errors, and other relevant metrics.
#'
#' @name show
#' @param pfimproject An object of class \code{PFIMProject} representing the evaluation.
#' @return Displays the design evaluation results directly in the console.
#'
#' @examples
#' \dontrun{
#' # Assuming 'myPFIMproject' is your current project object
#' # Running 'run' usually returns an object that can be shown
#' evaluationResults = run(myPFIMproject)
#'
#' # Display the results
#' show(evaluationResults)
#' }
#' @export
# ==============================================================================

method( show, Evaluation ) = function( pfimproject )
{
  fim = prop( pfimproject, "fim" )
  fim = setEvaluationFim( fim, pfimproject )
  showFIM( fim )
}

# ==============================================================================
#' plotEvaluation: Plot model responses for design evaluation
#'
#' @description
#' Generates graphical representations of the model responses based on the
#' design and parameters defined in the PFIM project.
#'
#' @name plotEvaluation
#' @param pfimproject An object of class \code{PFIMProject} containing the design.
#' @param plotOptions A \code{list} specifying graphical parameters (e.g., \code{xlim}, \code{ylim}, \code{main}, \code{col}).
#' @return A series of plots representing the expected model responses and sampling schedules.
#'
#' @examples
#' \dontrun{
#' # 1. Define custom plot options, e.g. form the Vignette 01:
#' myPlotOptions = list( unitTime = c( "hour" ), unitOutcomes= c( "mcg/mL" , "DI%" ) )
#'
#' # 2. Generate the plots for the project
#' plotEvaluation(myPFIMproject, plotOptions = myPlotOptions)
#' }
#' @export
# ==============================================================================

method( plotEvaluation, Evaluation ) = function( pfimproject, plotOptions )
{
  designs = prop( pfimproject, "designs" )
  model = defineModelType( pfimproject ) |> finiteDifferenceHessian() |> defineModelWrapper( pfimproject )
  fim = defineFim( pfimproject )
  design = pluck( designs, 1 )
  designName = prop( design, "name" )
  arms = prop( design, "arms" )
  # generate and print all plots
  allPlots = map( arms, ~ processArmEvaluationResults( .x, model, fim, designName, plotOptions ) )
  allPlots = setNames( list( allPlots |> map( ~ .x[[designName]] ) |> list_flatten() ), designName )
  return( allPlots )
}

# ==============================================================================
#' plotSensitivityIndices: Plot gradients of the model responses
#'
#' @description
#' Generates plots for the sensitivity indices (partial derivatives) of the model
#' responses with respect to the population parameters. This is essential for
#' identifying the optimal sampling times where the model is most sensitive to
#' parameter changes.
#'
#' @name plotSensitivityIndices
#' @param pfimproject An object of class \code{PFIMProject} containing the design and model specifications.
#' @param plotOptions A \code{list} specifying graphical parameters for the sensitivity plots.
#' @return A series of plots showing the sensitivity of model responses over time for each parameter.
#'
#' @examples
#' \dontrun{
#' # 1. Define custom plot options, e.g. form the Vignette 01:
#' myPlotOptions = list( unitTime = c( "hour" ), unitOutcomes= c( "mcg/mL" , "DI%" ) )
#'
#' # 2. Generate sensitivity index plots for the current project
#' plotSensitivityIndices(myPFIMproject, plotOptions = myPlotOptions)
#' }
#' @export
# ==============================================================================

method( plotSensitivityIndices, Evaluation ) = function( pfimproject, plotOptions )
{
  designs = prop( pfimproject, "designs" )
  model = defineModelType( pfimproject ) |> finiteDifferenceHessian() |> defineModelWrapper( pfimproject )
  fim = defineFim( pfimproject )
  design = pluck( designs, 1 )
  designName = prop( design, "name" )
  arms = prop( design, "arms" )

  # generate and print all plots
  allPlots = map( arms, ~ processArmEvaluationSI( .x, model, fim, designName, plotOptions ) )
  allPlots = setNames( list( allPlots |> map( ~ .x[[designName]] ) |> list_flatten() ), designName )
  return( allPlots )
}

# ==============================================================================
#' plotSE: Bar plot of Standard Errors (SE)
#'
#' @description
#' Generates a bar plot showing the Standard Errors (SE) for the fixed effects
#' and variance components of the model. This visualization helps assess the
#' expected precision of the parameter estimates for the current design.
#'
#' @name plotSE
#' @param pfimproject An object of class \code{PFIMProject} containing the evaluation results.
#' @return A bar plot displaying the calculated SE for each model parameter.
#'
#' @examples
#' \dontrun{
#' # Assuming 'myPFIMproject' has been evaluated using run()
#'
#' # Generate the bar plot of Standard Errors
#' plotSE(myPFIMproject)
#' }
#' @export
# ==============================================================================

# plot SE  from evaluation
method( plotSE, Evaluation ) = function( pfimproject )
{
  # set the FIM and plot SE
  fim = prop( pfimproject, "fim" )
  plotSE = plotSEFIM( fim, pfimproject )
  return( plotSE )
}

# ==============================================================================
#' plotRSE: Bar plot of Relative Standard Errors (RSE)
#'
#' @description
#' Generates a bar plot showing the Relative Standard Errors (RSE, expressed as a
#' percentage) for the model parameters. This visualization is essential for
#' checking if the design meets the target precision criteria (e.g., RSE < 20%).
#'
#' @name plotRSE
#' @param pfimproject An object of class \code{PFIMProject} containing the evaluation results.
#' @return A bar plot displaying the RSE (%) for each model parameter.
#'
#' @examples
#' \dontrun{
#' # Assuming 'myPFIMproject' has been evaluated using the run() function
#'
#' # Generate the bar plot of Relative Standard Errors
#' plotRSE(myPFIMproject)
#' }
#' @export
# ==============================================================================

method( plotRSE, Evaluation ) = function( pfimproject )
{
  # set the FIM and plot RSE
  fim = prop( pfimproject, "fim" )
  plotRSE = plotRSEFIM( fim, pfimproject )
  return( plotRSE )
}

# ==============================================================================
#' getSE: get the SE
#' @name getSE
#' @param pfimproject A object \code{PFIMProject} giving the Evaluation.
#' @return The SE.
#' @export
# ==============================================================================

method( getSE, Evaluation ) = function( pfimproject )
{
  # set the FIM and plot SE
  fim = prop( pfimproject, "fim" )
  fim = setEvaluationFim( fim, pfimproject )
  SEAndRSE = prop( fim, "SEAndRSE" )
  SE = SEAndRSE$SE
  return( SE )
}

# ==============================================================================
#' getRSE: get the RSE
#' @name getRSE
#' @param pfimproject A object \code{PFIMProject} giving the Evaluation.
#' @return The RSE
#' @export
# ==============================================================================

method( getRSE, Evaluation ) = function( pfimproject )
{
  # set the FIM and plot SE
  fim = prop( pfimproject, "fim" )
  fim = setEvaluationFim( fim, pfimproject )
  SEAndRSE = prop( fim, "SEAndRSE" )
  RSE = SEAndRSE$RSE
  return( RSE )
}

# ==============================================================================
#' getShrinkage: get the shrinkage
#' @name getShrinkage
#' @param pfimproject A object \code{PFIMProject} giving the Evaluation.
#' @return The shrinkage
#' @export
# ==============================================================================

method( getShrinkage, Evaluation ) = function( pfimproject )
{
  # set the FIM and plot SE
  fim = prop( pfimproject, "fim" )
  fim = setEvaluationFim( fim, pfimproject )
  shrinkage = prop( fim, "shrinkage" )
  return( shrinkage )
}

# ==============================================================================
#' getDeterminant: get the determinant
#' @name getDeterminant
#' @param pfimproject A object \code{PFIMProject} giving the Evaluation.
#' @return The determinant
#' @export
# ==============================================================================

method( getDeterminant, Evaluation ) = function( pfimproject )
{
  fisherMatrix = getFisherMatrix( pfimproject )
  return( det( fisherMatrix$fisherMatrix ) )
}

# ==============================================================================
#' getDcriterion : get the Dcriterion
#' @name getDcriterion
#' @param pfimproject A object \code{PFIMProject} giving the Evaluation.
#' @return The Dcriterion
#' @export
# ==============================================================================

method( getDcriterion, Evaluation ) = function( pfimproject )
{
  fim = prop( pfimproject, "fim" )
  fim = setEvaluationFim( fim, pfimproject )
  return( Dcriterion( fim ) )
}

# ==============================================================================
#' getCorrelationMatrix : get the correlation matrix
#' @name getCorrelationMatrix
#' @param pfimproject A object \code{PFIMProject} giving the Evaluation.
#' @return The Dcriterion
#' @export
# ==============================================================================

method( getCorrelationMatrix, Evaluation ) = function( pfimproject )
{
  fisherMatrix = getFisherMatrix( pfimproject )
  fisherMatrix = fisherMatrix$fisherMatrix

  return( cor( fisherMatrix ) )
}

# ==============================================================================
#' Report: generate the report.
#' @name Report
#' @param pfimproject A object \code{PFIMProject} giving the Evaluation or Optimization.
#' @param outputPath A string giving the path where the output are saved.
#' @param outputFile A string giving the name of the output file.
#' @param plotOptions A list giving the plot options.
#' @return The html report of the design evaluation or optimization.
#' @export
# ==============================================================================

method( Report, Evaluation ) = function( pfimproject, outputPath, outputFile, plotOptions  )
{
  # projectName
  projectName = prop( pfimproject, "name" )

  # outputs
  evaluationOutputs = prop( pfimproject , "outputs" )

  # model
  model = defineModelType( pfimproject  )
  modelEquations = prop( model, "modelEquations" )

  # model error table
  modelError = prop( pfimproject, "modelError" )
  modelError = map( modelError, getModelErrorData )
  modelErrorData = map(modelError, ~ as.data.frame(.x, stringsAsFactors = FALSE )) |> list_rbind()
  colnames( modelErrorData ) = c( "Output", "Type", "$\\sigma_{slope}$", "$\\sigma_{inter}$" )

  modelErrorTable = kbl( modelErrorData, align = c( "c","c","c","c" ) ) |>
    kable_styling( bootstrap_options = c(  "hover" ), full_width = FALSE, position = "center", font_size = 13 )

  # model parameters table
  modelParameters = prop( pfimproject, "modelParameters" )
  modelParameters = map( modelParameters, getModelParametersData )
  modelParametersData = map(modelParameters, ~ as.data.frame(.x, stringsAsFactors = FALSE )) |> list_rbind()
  colnames( modelParametersData ) = c("Parameter","$\\mu$","$\\omega$","Distribution", paste0("$\\mu$"," fixed"), paste0("$\\omega$"," fixed"))

  modelParametersTable = kbl( modelParametersData, align = c( "l","l","l","c","c","c" ) ) |>
    kable_styling( bootstrap_options = c( "hover" ), full_width = FALSE, position = "center", font_size = 13 )

  # arms
  designs = prop( pfimproject, "designs" )
  designsNames = map_chr( designs, "name" )
  arms = map( designs, ~ prop( .x, "arms" ))
  armsData = list_flatten( map( pluck(arms,1), getArmData ) )

  # administration table
  administration = list_flatten( map( pluck(arms, 1), armAdministration ) )
  administrationData = map(administration, ~ as.data.frame(.x, stringsAsFactors = FALSE )) |> list_rbind()
  colnames( administrationData ) = c( "Design name","Arms name" , "Number of subject ", "Outcome", "Dose","Time dose", "$\\tau$", "$T_{inf}$" )

  administrationTable = kbl( administrationData, align = c( "l","l","l","c","c","c","c" ) ) |>
    kable_styling( bootstrap_options = c(  "hover" ), full_width = FALSE, position = "center", font_size = 13 )

  # initial design
  initialDesignData = map(armsData, ~ as.data.frame( .x, stringsAsFactors = FALSE )) |> list_rbind()
  colnames( initialDesignData ) = c( "Arms name" , "Number of subjects", "Outcome", "Dose","Sampling times" )

  initialDesignTable = kbl( initialDesignData, align = c( "l","c","c","c") ) |>
    kable_styling( bootstrap_options = c( "hover" ), full_width = FALSE, position = "center", font_size = 13 )

  # Fisher matrix and SE
  fim = prop( pfimproject, "fim" )
  fim = setEvaluationFim( fim, pfimproject )
  fimInitialDesignTable = tablesForReport( fim, pfimproject )

  # plotsEvaluation & plotSensitivityIndices
  plotsEvaluation = plotEvaluation( pfimproject, plotOptions )
  plotSensitivityIndices = plotSensitivityIndices( pfimproject, plotOptions )
  plotSE = plotSE( pfimproject )
  plotRSE = plotRSE( pfimproject )

  # tablesForReport
  tablesForReport = list(
    evaluationOutputs = evaluationOutputs,
    modelEquations = modelEquations,
    modelErrorTable = modelErrorTable,
    modelParametersTable = modelParametersTable,
    administrationTable = administrationTable,
    initialDesignTable = initialDesignTable,
    fimInitialDesignTable = fimInitialDesignTable,
    plotsEvaluation = plotsEvaluation,
    plotSensitivityIndices = plotSensitivityIndices,
    plotSE = plotSE,
    plotRSE = plotRSE,
    fim = fim,
    pfimproject = pfimproject,
    projectName = projectName )

  # generate the report
  generateReportEvaluation( fim, tablesForReport,
                            outputFile = outputFile,
                            outputPath = outputPath )
}
