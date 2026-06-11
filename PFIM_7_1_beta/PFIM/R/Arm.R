#' @description The class \code{Arm} represents an arm and stores information concerning an arm.
#' @title Arm
#' @param name A string giving the name of the arm.
#' @param size A integer giving the size of the arm.
#' @param administrations A list giving the objects of class \code{Administration} that define the administrations of the arm.
#' @param initialConditions A list giving the initial conditions for the ode model where the names are string that define the variable and their value are giving by double
#' @param samplingTimes A list giving the objects of class \code{SamplingTime} that define the sampling time of the arm.
#' @param administrationsConstraints A list giving the objects of class \code{AdministrationsConstraints} that define the administration constraints of the arm.
#' @param samplingTimesConstraints A list giving the objects of class \code{SamplingTimeConstraints} that define the sampling time constraints of the arm.
#' @param evaluationModel A list giving the evaluation of the responses of the arm.
#' @param evaluationGradients A list giving the evaluation of the responses gradient of the arm.
#' @param evaluationVariance A list giving the evaluation of the variance.
#' @param evaluationFim A object of class \code{Fim} giving the Fisher Information Matrix.
#' @param covariateSampling A list giving the covariate sampling method.
#' @include Model.R
#' @include Fim.R
#' @include MultiplicativeAlgorithm.R
#' @include FedorovWynnAlgorithm.R
#' @include SimplexAlgorithm.R
#' @include PSOAlgorithm.R
#' @include PGBOAlgorithm.R
#' @include CovariateSamplingMethods.R
#' @export

Arm = new_class("Arm", package = "PFIM",

                properties = list(
                  name = new_property(class_character, default = character(0)),
                  size = new_property(class_double, default = numeric(0)),
                  administrations = new_property(class_list, default = list()),
                  initialConditions = new_property(class_list, default = list()),
                  samplingTimes = new_property(class_list, default = list()),
                  administrationsConstraints = new_property(class_list, default = list()),
                  samplingTimesConstraints = new_property(class_list, default = list()),
                  evaluationModel = new_property(class_list, default = list()),
                  evaluationGradients = new_property(class_list, default = list()),
                  evaluationVariance = new_property(class_list, default = list()),
                  evaluationFim = new_property(Fim, default = NULL),
                  covariateSampling = new_property(CovariateSamplingMethods, default = NULL)
                ))



evaluateArm = new_generic( "evaluateArm", c( "arm" ) )
getSamplingData = new_generic( "getSamplingData", c( "arm" ) )
updateSamplingTimes = new_generic( "updateSamplingTimes", c( "arm" ) )
processArmEvaluationResults = new_generic( "processArmEvaluationResults", c( "arm", "model", "fim" ) )
processArmEvaluationSI = new_generic( "processArmEvaluationSI", c( "arm", "model", "fim" ) )
plotEvaluationResults = new_generic( "plotEvaluationResults", c( "arm" ) )
plotEvaluationSI = new_generic( "plotEvaluationSI", c( "arm" ) )

getArmData = new_generic( "getArmData", c( "arm" ) )
getArmConstraints = new_generic( "getArmConstraints", c( "arm", "optimizationAlgorithm" ) )
armAdministration = new_generic( "armAdministration", c( "arm" ) )

# ==============================================================================
#' evaluateArm: evaluation of the model with the arm parameters.
#' @name evaluateArm
#' @param arm A object of class \code{Arm} giving the arm.
#' @param model A object of class \code{Model} giving the model.
#' @param fim A object of class \code{Fim} giving the fim.
#' @return The object arm with the slots evaluationModel, evaluationGradients, evaluationVariance and evaluationFim.
#' @export
# ==============================================================================

method( evaluateArm, Arm ) = function( arm, model, fim ) {

  # define the data for the ode solver
  model = defineModelAdministration( model, arm )

  # evaluate the model
  prop( arm, "evaluationModel" ) = evaluateModel( model, arm )

  # evaluate the model gradients
  prop( arm, "evaluationGradients" ) = evaluateModelGradient( model, arm )

  # evaluate the variance of the model
  prop( arm, "evaluationVariance" ) = evaluateModelVariance( model, arm )

  # evaluate the Fim
  prop( arm, "evaluationFim" ) = evaluateFim( fim, model, arm )

  return( arm )
}

# ==============================================================================
#' getArmAdministration: get the administration parameters of an arm.
#' @name armAdministration
#' @param arm A object of class \code{Arm} giving the arm.
#' @return A list giving the administration parameters of an arm.
#' @export
# ==============================================================================

method( armAdministration, Arm ) = function( arm ) {

  armName = prop( arm, "name")
  armSize = round( prop( arm, "size" ), 2 )
  adminList = prop( arm, "administrations" )

  admin = map( adminList, function( adm ) {
    outcome = prop( adm, "outcome" )
    dose = as.character( prop( adm, "dose" ) )
    timeDose = as.character( prop( adm, "timeDose" ) )
    Tinf = prop( adm, "Tinf" )
    Tinf_str = if ( length( Tinf ) > 0 ) as.character( Tinf ) else "."

    list(
      "Design name"        = "Design optimized",
      "Arms name"          = armName,
      "Number of subjects" = as.character( armSize ),
      "Outcome"            = prop( adm, "outcome" ),
      "Dose"               = as.character( prop( adm, "dose" ) ),
      "Time of dose"       = if ( length( prop( adm, "timeDose" ) ) > 0 ) as.character( prop( adm, "timeDose" ) ) else ".",
      "tau"                = as.character( prop( adm, "tau" ) ),
      "Tinf"               = if ( length( prop( adm, "Tinf" ) ) > 0 ) as.character( prop( adm, "Tinf" ) ) else "."
    )
  })
  return( admin )
}

# ==============================================================================
#' getArmConstraints: get the administration and sampling time constraints for the MultiplicativeAlgorithm.
#' @name getArmConstraints
#' @param arm A object of class \code{Arm} giving the arm.
#' @param optimizationAlgorithm A object of class \code{Optimization} giving the optimization algorithm.
#' @return A list giving the administration and sampling time constraints for the MultiplicativeAlgorithm.
#' @export
# ==============================================================================

method( getArmConstraints, list( Arm, MultiplicativeAlgorithm ) ) = function( arm, optimizationAlgorithm ) {

  # Récupération des données du bras
  armName = prop( arm, "name" )
  armSize = prop( arm, "size" )

  # table de correspondance Outcome - Dose constraints
  admins = map( prop( arm, "administrationsConstraints" ), function( administrationConstraint ) {
    outcome = prop( administrationConstraint, "outcome" )
    doses = paste0( "(", paste( unlist( prop( administrationConstraint, "doses" ) ), collapse = ", " ), ")" )
    setNames( list( doses ), outcome )
  })

  admins = flatten( admins )

  # liste finale à partir des sampling constraints
  armConstraints = map( prop( arm, "samplingTimesConstraints" ), function( samplingConstraint ) {
    outcome = prop( samplingConstraint, "outcome" )
    initialSamplings = paste0( "(", paste( prop( samplingConstraint, "initialSamplings" ), collapse = ", " ), ")" )
    fixedTimes = paste0( "(", paste( prop( samplingConstraint, "fixedTimes" ), collapse = ", " ), ")" )
    numberOfsamplingsOptimisable = as.character( prop( samplingConstraint, "numberOfsamplingsOptimisable" ) )

    # Lookup dose constraints, or NA if not found
    doseConstraints = admins[[outcome]]
    if ( is.null( doseConstraints ) ) doseConstraints = '.'

    list( "Arms name" = armName,
          "Number of subjects" = armSize,
          "Outcome" = outcome,
          "Initial samplings" = initialSamplings,
          "Fixed times" = fixedTimes,
          "Number of samplings optimisable" = numberOfsamplingsOptimisable,
          "Dose constraints" = doseConstraints)
  })
  return( armConstraints )
}

# ==============================================================================
#' getArmConstraints: get the administration and sampling time constraints for the FedorovWynnAlgorithm.
#' @name getArmConstraints
#' @param arm A object of class \code{Arm} giving the arm.
#' @param optimizationAlgorithm A object of class \code{Optimization} giving the optimization algorithm.
#' @return A list giving the administration and sampling time constraints for the FedorovWynnAlgorithm.
#' @export
# ==============================================================================

method( getArmConstraints, list( Arm, FedorovWynnAlgorithm ) ) = function( arm, optimizationAlgorithm ) {

  # Récupération des métadonnées du bras
  armName = prop( arm, "name" )
  armSize = prop( arm, "size" )

  # table de correspondance Outcome - Dose constraints
  admins = map( prop( arm, "administrationsConstraints" ), function( administrationConstraint ) {
    outcome = prop( administrationConstraint, "outcome" )
    doses = paste0( "(", paste( unlist( prop( administrationConstraint, "doses" ) ), collapse = ", " ), ")" )
    setNames( list( doses ), outcome )
  })

  admins = flatten( admins )

  # liste finale à partir des sampling constraints
  armConstraints = map( prop( arm, "samplingTimesConstraints" ), function( samplingConstraint ) {
    outcome = prop( samplingConstraint, "outcome" )
    initialSamplings = paste0( "(", paste( prop( samplingConstraint, "initialSamplings" ), collapse = ", " ), ")" )
    fixedTimes = paste0( "(", paste( prop( samplingConstraint, "fixedTimes" ), collapse = ", " ), ")" )
    numberOfsamplingsOptimisable = as.character( prop( samplingConstraint, "numberOfsamplingsOptimisable" ) )

    # Lookup dose constraints, or NA if not found
    doseConstraints = admins[[outcome]]
    if ( is.null( doseConstraints ) ) doseConstraints = '.'

    list( "Arms name" = armName,
          "Number of subjects" = armSize,
          "Outcome" = outcome,
          "Initial samplings" = initialSamplings,
          "Fixed times" = fixedTimes,
          "Number of samplings optimisable" = numberOfsamplingsOptimisable,
          "Dose constraints" = doseConstraints )
  })

  return( armConstraints )
}

# ==============================================================================
#' getArmConstraints: get the administration and sampling time constraints for the SimplexAlgorithm.
#' @name getArmConstraints
#' @param arm A object of class \code{Arm} giving the arm.
#' @param optimizationAlgorithm A object of class \code{Optimization} giving the optimization algorithm.
#' @return A list giving the administration and sampling time constraints for the SimplexAlgorithm.
#' @export
# ==============================================================================

method( getArmConstraints, list( Arm, SimplexAlgorithm ) ) = function( arm, optimizationAlgorithm ) {

  armName = prop( arm, "name" )
  armSize = prop( arm, "size" )

  armConstraints = map(prop(arm, "samplingTimesConstraints"), function(samplingConstraint) {
    outcome = prop(samplingConstraint, "outcome")
    initialSamplings = paste0("(", paste(prop(samplingConstraint, "initialSamplings"), collapse = ", "), ")")
    samplingsWindows = paste( map_chr(prop(samplingConstraint, "samplingsWindows"), ~ paste0("(", paste(., collapse = ","), ")")), collapse = ", " )
    numberOfTimesByWindows = paste0("(", paste(prop(samplingConstraint, "numberOfTimesByWindows"), collapse = ", "), ")")
    minSampling = paste0("(", paste(prop(samplingConstraint, "minSampling"), collapse = ", "), ")")

    list( "Arms name" = armName,
          "Number of subjects" = armSize,
          "Outcome" = outcome,
          "Initial samplings" = initialSamplings,
          "Samplings windows" = samplingsWindows,
          "Number of tules by windows" = numberOfTimesByWindows,
          "Min sampling" = minSampling )
  })
  return(armConstraints)
}

# ==============================================================================
#' getArmConstraints: get the administration and sampling time constraints for the PSOAlgorithm.
#' @name getArmConstraints
#' @param arm A object of class \code{Arm} giving the arm.
#' @param optimizationAlgorithm A object of class \code{Optimization} giving the optimization algorithm.
#' @return A list giving the administration and sampling time constraints for the PSOAlgorithm.
#' @export
# ==============================================================================

method( getArmConstraints, list( Arm, PSOAlgorithm ) ) = function( arm, optimizationAlgorithm ) {
  armName = prop(arm, "name")
  armSize = prop(arm, "size")

  armConstraints = map(prop(arm, "samplingTimesConstraints"), function(samplingConstraint) {
    outcome = prop(samplingConstraint, "outcome")
    initialSamplings = paste0("(", paste(prop(samplingConstraint, "initialSamplings"), collapse = ", "), ")")
    samplingsWindows = paste( map_chr(prop(samplingConstraint, "samplingsWindows"), ~ paste0("(", paste(., collapse = ","), ")")), collapse = ", " )
    numberOfTimesByWindows = paste0("(", paste(prop(samplingConstraint, "numberOfTimesByWindows"), collapse = ", "), ")")
    minSampling = paste0("(", paste(prop(samplingConstraint, "minSampling"), collapse = ", "), ")")

    list( "Arms name" = armName,
          "Number of subjects" = armSize,
          "Outcome" = outcome,
          "Initial samplings" = initialSamplings,
          "Samplings windows" = samplingsWindows,
          "Number of tules by windows" = numberOfTimesByWindows,
          "Min sampling" = minSampling )
  })
  return(armConstraints)
}

# ==============================================================================
#' getArmConstraints: get the administration and sampling time constraints for the PGBOAlgorithm.
#' @name getArmConstraints
#' @param arm A object of class \code{Arm} giving the arm.
#' @param optimizationAlgorithm A object of class \code{Optimization} giving the optimization algorithm.
#' @return A list giving the administration and sampling time constraints for the PGBOAlgorithm.
#' @export
# ==============================================================================

method( getArmConstraints, list( Arm, PGBOAlgorithm ) ) = function( arm, optimizationAlgorithm ) {
  armName = prop(arm, "name")
  armSize = prop(arm, "size")

  armConstraints = map(prop(arm, "samplingTimesConstraints"), function(samplingConstraint) {
    outcome = prop(samplingConstraint, "outcome")
    initialSamplings = paste0("(", paste(prop(samplingConstraint, "initialSamplings"), collapse = ", "), ")")
    samplingsWindows = paste( map_chr(prop(samplingConstraint, "samplingsWindows"), ~ paste0("(", paste(., collapse = ","), ")")), collapse = ", " )
    numberOfTimesByWindows = paste0("(", paste(prop(samplingConstraint, "numberOfTimesByWindows"), collapse = ", "), ")")
    minSampling = paste0("(", paste(prop(samplingConstraint, "minSampling"), collapse = ", "), ")")

    list( "Arms name" = armName,
          "Number of subjects" = armSize,
          "Outcome" = outcome,
          "Initial samplings" = initialSamplings,
          "Samplings windows" = samplingsWindows,
          "Number of tules by windows" = numberOfTimesByWindows,
          "Min sampling" = minSampling )
  })
  return(armConstraints)
}

# ==============================================================================
#' getArmData: extract arm data for The Report
#' @name getArmData
#' @param arm A object of class \code{Arm} giving the arm.
#' @return A list giving the name, Number of subjects, Outcome, Dose and Sampling times of the arm.
#' @export
# ==============================================================================

method( getArmData, Arm ) = function( arm ) {
  armName = prop(arm, "name")
  armSize = round(prop(arm, "size"), 2)

  administrations = prop(arm, "administrations")
  doseList = map(administrations, function(adm) {
    list(
      outcome = prop(adm, "outcome"),
      dose = prop(adm, "dose") )
  })

  doseMap = flatten(doseList)
  doseDict = setNames(
    map(doseList, ~ paste(.x$dose, collapse = ", ")),
    map_chr(doseList, ~ .x$outcome) )

  samplingList = prop(arm, "samplingTimes")
  samplingOutcomes = map_chr(samplingList, ~ prop(.x, "outcome"))
  samplingTimes = map(samplingList, ~ prop(.x, "samplings"))

  armData = map2(samplingOutcomes, samplingTimes, function(outc, samps) {
    doseVal = if (outc %in% names(doseDict)) doseDict[[outc]] else "."
    list( "Arms name" = armName,
          "Number of subjects" = as.character(armSize),
          "Outcome" = outc,
          "Dose" = doseVal,
          "Sampling times" = paste0("(", paste( round(samps,2), collapse = ", "), ")") )
  })
  return( armData )
}

# ==============================================================================
#' getSamplingData: extract sampling times and max sampling time used for plot.
#' @name getSamplingData
#' @param arm A object of class \code{Arm} giving the arm.
#' @return A list giving the \code{samplingTimes} object, the vector samplings and the double samplingMax.
#' @export
# ==============================================================================

method( getSamplingData, Arm ) = function( arm ) {
  samplingTimes = prop( arm, "samplingTimes" )
  samplings = map(samplingTimes, ~ prop( .x, "samplings" ) ) %>%
    set_names(map_chr(samplingTimes, ~ prop(.x ,"outcome" ) ) )
  samplingMax = samplings %>% flatten_dbl() %>% max()
  list(samplingTimes = samplingTimes, samplings = samplings, samplingMax = samplingMax)
}

# ==============================================================================
#' updateSamplingTimes: update sampling times for plotting used for plot
#' @name updateSamplingTimes
#' @param arm A object of class \code{Arm} giving the arm.
#' @param samplingData The list giving as output in the method getSamplingData.
#' @return The updated sampling times.
#' @export
# ==============================================================================

method( updateSamplingTimes, Arm ) = function( arm, samplingData ) {
  prop( arm, "samplingTimes" ) = map(samplingData$samplingTimes, function(samplingsPlot) {
    prop( samplingsPlot, "samplings" ) = sort(unique(c(samplingData$samplings %>% flatten_dbl(), seq(0.0, samplingData$samplingMax, 0.1))))
    return(samplingsPlot)
  })
  arm
}

# ==============================================================================
#' processArmEvaluationResults: process for the evaluation of an arm.
#' @name processArmEvaluationResults
#' @param arm A object of class \code{Arm} giving the arm.
#' @param model A object of class \code{Model} giving the model.
#' @param fim A object of class \code{Fim} giving the fim.
#' @param designName A string giving the name of the design.
#' @param plotOptions A list giving the plot options.
#' @return A list of ggplot object giving the plot of the responses ans the gradient responses of the the model.
#' @export
# ==============================================================================

method( processArmEvaluationResults, list( Arm, Model, Fim ) ) = function( arm, model, fim, designName, plotOptions ) {

  outputNames = as.list(prop(model, "outputNames"))

  # Get sampling data
  samplingData = getSamplingData(arm)
  # Update sampling times for plotting
  arm = updateSamplingTimes(arm, samplingData)
  # Model evaluation with new sampling times
  model = defineModelAdministration(model, arm)
  evaluationModel = evaluateModel(model, arm)
  # Plot evaluation results
  plots = plotEvaluationResults( arm, evaluationModel, outputNames, samplingData, designName, plotOptions )
  return(plots)
}

# ==============================================================================
#' processArmEvaluationSI: process for the evaluation of the gradient of the responses.
#' @name processArmEvaluationSI
#' @param arm A object of class \code{Arm} giving the arm.
#' @param model A object of class \code{Model} giving the model.
#' @param fim A object of class \code{Fim} giving the fim.
#' @param designName A string giving the name of the design.
#' @return A list giving the ggplot object of the plots of the gradient.
#' @export
# ==============================================================================

method( processArmEvaluationSI, list( Arm, Model, Fim ) ) = function( arm, model, fim, designName, plotOptions ) {

  outputNames = as.list(prop(model, "outputNames"))

  # Get sampling data
  samplingData = getSamplingData(arm)

  # Update sampling times for plotting
  arm = updateSamplingTimes(arm, samplingData)

  # Model evaluation with new sampling times
  model = defineModelAdministration(model, arm)
  parametersNames = prop(model, "modelParameters") %>% map_chr(~ prop(.x, "name"))
  evaluationModelGradient = evaluateModelGradient(model, arm)
  timeSeq = seq(from = 0, by = 0.1, length.out = nrow( pluck( evaluationModelGradient, 1 ) ) )
  evaluationModelGradient = map2(outputNames, evaluationModelGradient, function(outputName, gradient) { data.frame(time = timeSeq, gradient) }) %>% set_names(outputNames)

  # Plot evaluation results
  plots = plotEvaluationSI( arm, evaluationModelGradient, parametersNames, outputNames, samplingData, designName, plotOptions )
  return(plots)
}

# ==============================================================================
#' plotEvaluationResults: process for the evaluation of the responses.
#' @name plotEvaluationResults
#' @param arm A object of class \code{Arm} giving the arm.
#' @param evaluationModel A list giving the evaluation of the model.
#' @param outputNames A list of string giving the output of the evaluation of the model.
#' @param samplingData A list giving the sampling data from the method getSamplingData.
#' @param unitXAxis A list giving the unit of the x-axis.
#' @param unitYAxis A list giving the unit of the y-axis.
#' @param designName A string giving the design name.
#' @return A list giving the plot of the evaluation of the model responses.
#' @export
# ==============================================================================

method( plotEvaluationResults, Arm ) = function( arm, evaluationModel, outputNames, samplingData, designName, plotOptions ) {

  # unit axis
  plotOptions = lapply(plotOptions, function(x) if (is.null(x)) " " else x)
  unitXAxis = plotOptions$unitTime
  unitYAxis = setNames( plotOptions$unitOutcomes, unlist( outputNames ) )
  plotList = list()
  armName = prop( arm, "name")
  plotList[[designName]] = list()
  plotList[[designName]][[armName]] = list()

  plots =  map2(outputNames, samplingData$samplings, function(outputName, sampling) {

    data = evaluationModel[[outputName]]

    samplingPoints = data[data$time %in% sampling, ]

    ggplot(data, aes(x = time, y = .data[[outputName]])) +
      geom_line() +
      geom_point(data = samplingPoints, aes(x = time, y = .data[[outputName]]), color = "red") +
      labs( x = paste0("Time (", unitXAxis, ")\n\nDesign: ", sub("_", " ", designName), "      Arm: ", armName ),
            y = paste0(outputName, " (", unitYAxis[[outputName]], ")\n" ) ) +
      scale_x_continuous(breaks = pretty_breaks(n = 10), sec.axis = sec_axis(~ . * 1, breaks = round(sampling, 2), name = "Sampling times")) +
      scale_y_continuous(breaks = pretty_breaks(n = 10)) +
      theme( legend.position = "none",
             axis.title.x.top = element_text(color = "red", vjust = 2.0),
             axis.text.x.top = element_text(angle = 90, hjust = 0, color = "red"),
             plot.title = element_text(size = 16, hjust = 0.5),
             axis.title.x = element_text(size = 16),
             axis.title.y = element_text(size = 16),
             axis.text.x = element_text(size = 16, angle = 90, vjust = 0.5),
             axis.text.y = element_text(size = 16, angle = 0, vjust = 0.5, hjust = 0.5),
             strip.text.x = element_text(size = 16) )
  })

  plotList[[designName]][[armName]] = set_names( plots, outputNames )
  return( plotList )
}

# ==============================================================================
#' plotEvaluationSI: process for the evaluation of the gradient of the responses.
#' @name plotEvaluationSI
#' @param arm A object of class \code{Arm} giving the arm.
#' @param evaluationModelGradient A list giving the evaluation of the gradient of the model responses.
#' @param parametersNames A vector of string giving the parameter names?
#' @param outputNames A list of string giving the name of the outputs.
#' @param samplingData A list giving the sampling data from the method getSamplingData.
#' @param unitXAxis A list giving the unit of the x-axis.
#' @param unitYAxis A list giving the unit of the y-axis.
#' @param designName A string giving the design name.
#' @return A list giving the plot of the evaluation of gradient of the model responses.
#' @export
# ==============================================================================

method( plotEvaluationSI, Arm ) = function( arm, evaluationModelGradient, parametersNames, outputNames, samplingData, designName, plotOptions ) {

  unitXAxis = plotOptions$unitTime
  armName = prop(arm, "name")
  plotList = list()
  plotList[[designName]] = list()
  plotList[[designName]][[armName]] = list()

  plots = map2(outputNames, samplingData$samplings, function(outputName, sampling) {

    gradientData = evaluationModelGradient[[outputName]]

    minYAxis = min(gradientData[, parametersNames], na.rm = TRUE)
    maxYAxis = max(gradientData[, parametersNames], na.rm = TRUE)

    map(parametersNames, function(parameterName) {

      data = as_tibble(gradientData[, c("time", parameterName)])
      names(data)[2] = "parameterValue"

      samplingPoints = data[data$time %in% sampling, ]

      ggplot(data, aes(x = time, y = parameterValue)) +
        geom_line() +
        geom_point(data = samplingPoints, color = "red") +
        labs( y = paste0("df/d", parameterName),
              x = paste0(
                "Time (", unitXAxis, ")\n\n",
                "Design: ", gsub("_", " ", designName), "   ",
                "Arm: ", armName, "   ",
                "Output: ", outputName, "   ",
                "Parameter: ", parameterName ) ) +
        scale_x_continuous( breaks = pretty_breaks(n = 10),
                            sec.axis = sec_axis(~., breaks = round(sampling, 2), name = "Sampling times") ) +
        scale_y_continuous( breaks = pretty_breaks(n = 10),
                            limits = c(minYAxis, maxYAxis) ) +
        theme( legend.position = "none",
               axis.title.x.top = element_text(color = "red", vjust = 2.0),
               axis.text.x.top = element_text(angle = 90, hjust = 0, color = "red"),
               plot.title = element_text(size = 16, hjust = 0.5),
               axis.title.x = element_text(size = 16),
               axis.title.y = element_text(size = 16),
               axis.text.x = element_text(size = 16, angle = 90, vjust = 0.5),
               axis.text.y = element_text(size = 16),
               strip.text.x = element_text(size = 16)
        )
    })
  })
  plotList[[designName]][[armName]] = map2(outputNames, plots, ~ setNames(.y, parametersNames)) %>% set_names(outputNames)
  return( plotList )
}














