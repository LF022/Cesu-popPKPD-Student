#' @description The class \code{Optimization} implements the Optimization.
#' @title Optimization
#' @inheritParams PFIMProject
#' @param optimisationDesign A list giving the evaluation of initial and optimal design.
#' @param optimisationAlgorithmOutputs A list giving the outputs of the optimization process.
#' @include PFIMProject.R
#' @export

Optimization = new_class( "Optimization",
                          package    = "PFIM",
                          parent     = PFIMProject,
                          properties = list(
                            optimisationDesign           = new_property( class_list, default = list() ),
                            optimisationAlgorithmOutputs = new_property( class_list, default = list() )
                          )
)

defineOptimizationAlgorithm  = new_generic( "defineOptimizationAlgorithm", c( "optimization" ) )
generateFimsFromConstraints  = new_generic( "generateFimsFromConstraints",  c( "optimization" ) )
plotWeights                  = new_generic( "plotWeights",                  c( "optimization" ) )
plotFrequencies              = new_generic( "plotFrequencies",              c( "optimization" ) )
optimizeDesign               = new_generic( "optimizeDesign",  c( "optimizationObject", "optimizationAlgorithm" ) )
constraintsTableForReport    = new_generic( "constraintsTableForReport",    c( "optimizationAlgorithm" ) )

# ── Private helpers ────────────────────────────────────────────────────────────

# Extract evaluationOptimalDesign from an Optimization object (used by all
# accessor methods to avoid repeating the same three lines everywhere).
.getOptimalEval = function( optimization ) {
  prop( optimization, "optimisationDesign" )$evaluationOptimalDesign
}

# Set FIM on the optimal-design evaluation and return the populated fim object.
.getOptimalFim = function( optimization ) {
  eval  = .getOptimalEval( optimization )
  fim   = prop( eval, "fim" )
  setEvaluationFim( fim, eval )
}

# ==============================================================================
#' generateFimsFromConstraints: enumerate and evaluate all FIMs from constraints
#' @name generateFimsFromConstraints
#' @param optimization An \code{Optimization} object.
#' @return A list with listArms, dimFim, listFimsAlgoFW, listFimsAlgoMult,
#'   samplingsForFedorovWynnAlgo.
#' @export
# ==============================================================================

method( generateFimsFromConstraints, Optimization ) = function( optimization ) {

  # ── Evaluation template ───────────────────────────────────────────────────────
  evaluation = Evaluation(
    name                    = "",
    modelEquations          = prop( optimization, "modelEquations"          ),
    modelParameters         = prop( optimization, "modelParameters"         ),
    modelError              = prop( optimization, "modelError"              ),
    modelCovariates         = prop( optimization, "modelCovariates"         ),
    modelCovariatesEquation = prop( optimization, "modelCovariatesEquation" ),
    numberOfOccasions       = prop( optimization, "numberOfOccasions"       ),
    designs                 = prop( optimization, "designs"                 ),
    fimType                 = prop( optimization, "fimType"                 ),
    outputs                 = prop( optimization, "outputs"                 ),
    odeSolverParameters     = prop( optimization, "odeSolverParameters"     )
  )

  designs     = prop( optimization, "designs" )
  designNames = map_chr( designs, ~ prop( .x, "name" ) )

  dosesForFIMs    = map( designs, ~ generateDosesCombination( .x ) ) |>
    set_names( designNames )
  samplingsForFIMs = map( designs, ~ generateSamplingTimesCombination( .x ) ) |>
    set_names( designNames )

  # ── Process every design → flat list of per-FIM results ──────────────────────
  allResults = imap( set_names( designs, designNames ), function( design, designName ) {

    arms           = prop( design, "arms" )
    dosesForDesign = dosesForFIMs[[ designName ]]
    numberOfDoses  = dosesForDesign$numberOfDoses

    combinationGrid = expand.grid(
      map( samplingsForFIMs[[ designName ]], seq_along )
    )
    nCombinations   = nrow( combinationGrid )
    totalIterations = numberOfDoses * nCombinations

    map( seq_len( numberOfDoses ), function( iterDose ) {

      # Assign new doses to all arms
      armsWithDoses = map( arms, function( arm ) {
        armName         = prop( arm, "name" )
        administrations = prop( arm, "administrations" )
        administrations = map( administrations, function( adm ) {
          prop( adm, "dose" ) = dosesForDesign[[ armName ]][[ prop( adm, "outcome" ) ]][ iterDose ]
          adm
        })
        prop( arm, "administrations" ) = administrations
        arm
      })

      map( seq_len( nCombinations ), function( iterComb ) {

        # Assign new sampling times; capture FW sampling vector simultaneously
        armsUpdated = map( armsWithDoses, function( arm ) {
          armName       = prop( arm, "name" )
          idx           = combinationGrid[ iterComb, armName ]
          samplingEntry = pluck( samplingsForFIMs, designName, armName, idx )
          prop( arm, "samplingTimes" ) = samplingEntry
          list(
            arm = arm,
            samplingsForFW = unlist( map( samplingEntry, ~ prop( .x, "samplings" ) ),
                                     use.names = FALSE )
          )
        })

        # Single-arm design assumed for constraint enumeration
        armResult  = armsUpdated[[ 1L ]]
        tempDesign = design
        prop( tempDesign, "arms"    ) = list( armResult$arm )
        tempEval   = evaluation
        prop( tempEval, "designs"   ) = list( tempDesign )

        fim          = getFim( run( tempEval ) )
        fisherMatrix = fim$fisherMatrix
        dimFim       = nrow( fisherMatrix )
        dimVec       = dimFim * ( dimFim + 1L ) / 2L

        # Lower-triangular + diagonal vectorisation for the FW C routine:
        # [(1,1), (2,1:2), (3,1:3), ...]  — reverse of upper-tri row-major
        fisherMatrixForAlgoFW = matrix(
          fisherMatrix[ rev( lower.tri( t( fisherMatrix ), diag = TRUE ) ) ],
          ncol  = dimVec,
          byrow = TRUE
        )

        fimIndex = ( iterDose - 1L ) * nCombinations + iterComb
        message( sprintf( "FIM evaluation: %d / %d", fimIndex, totalIterations ) )

        list(
          armResult             = armResult,
          samplingsForFW        = armResult$samplingsForFW,
          fisherMatrixForAlgoFW = fisherMatrixForAlgoFW,
          fisherMatrix          = fisherMatrix,
          dimFim                = dimFim
        )
      })
    }) |> list_flatten()
  })

  # ── Reconstruct named lists expected by downstream algorithms ─────────────────
  list(
    listArms                    = map( allResults, ~ map( .x, "armResult"             ) ),
    dimFim                      = allResults[[ 1L ]][[ 1L ]]$dimFim,
    listFimsAlgoFW              = map( allResults, ~ map( .x, "fisherMatrixForAlgoFW" ) ),
    listFimsAlgoMult            = map( allResults, ~ map( .x, "fisherMatrix"          ) ),
    samplingsForFedorovWynnAlgo = map( allResults, ~ map( .x, "samplingsForFW"        ) )
  )
}

# ==============================================================================
#' defineOptimizationAlgorithm: instantiate the algorithm from the optimizer slot
#' @name defineOptimizationAlgorithm
#' @param optimization An \code{Optimization} object.
#' @return An optimization algorithm object.
#' @export
# ==============================================================================

method( defineOptimizationAlgorithm, Optimization ) = function( optimization ) {
  optimizerParameters = prop( optimization, "optimizer" )

  switch( optimizerParameters,
          MultiplicativeAlgorithm = MultiplicativeAlgorithm(),
          FedorovWynnAlgorithm    = FedorovWynnAlgorithm(),
          PSOAlgorithm            = PSOAlgorithm(),
          PGBOAlgorithm           = PGBOAlgorithm(),
          SimplexAlgorithm        = SimplexAlgorithm(),
          stop( sprintf( "Unknown optimizer: '%s'", optimizerParameters ) )
  )
}

# ==============================================================================
#' run: run the optimization
#' @name run
#' @param optimization An \code{Optimization} object.
#' @return The optimization design results.
#' @export
# ==============================================================================

method( run, Optimization ) = function( pfimproject ) {
  optimizationAlgorithm = defineOptimizationAlgorithm( pfimproject )
  prop( pfimproject, "fim" ) = defineFim( pfimproject )

  if ( length( prop( pfimproject, "modelFromLibrary" ) ) != 0L )
    prop( pfimproject, "modelEquations" ) = defineModelEquationsFromLibraryOfModel( pfimproject )

  optimizeDesign( pfimproject, optimizationAlgorithm )
}

# ==============================================================================
#' show: display optimization results in the console
#' @name show
#' @param pfimproject An \code{Optimization} object.
#' @export
# ==============================================================================

method( show, Optimization ) = function( pfimproject ) {

  optimisationDesign      = prop( pfimproject, "optimisationDesign" )
  evaluationInitialDesign = optimisationDesign$evaluationInitialDesign
  evaluationOptimalDesign = optimisationDesign$evaluationOptimalDesign

  # Helper: extract arm data as a formatted data frame
  .armTable = function( evaluation ) {
    designs  = prop( evaluation, "designs" )
    armsData = list_flatten( map( pluck( map( designs, ~ prop( .x, "arms" ) ), 1L ),
                                  getArmData ) )
    df = map( armsData, ~ as.data.frame( .x, stringsAsFactors = FALSE ) ) |> list_rbind()
    colnames( df ) = c( "Arms name", "Number of subjects", "Outcome", "Dose", "Sampling times" )
    df
  }

  fimInitialDesign = setEvaluationFim( prop( evaluationInitialDesign, "fim" ), evaluationInitialDesign )
  fimOptimalDesign = setEvaluationFim( prop( evaluationOptimalDesign, "fim" ), evaluationOptimalDesign )

  cat( "\n===================================== \n" )
  cat( "  Initial design \n" )
  cat( "===================================== \n\n" )
  print( .armTable( evaluationInitialDesign ) )
  showFIM( fimInitialDesign )

  cat( "\n===================================== \n" )
  cat( "  Optimal design \n" )
  cat( "===================================== \n\n" )
  print( .armTable( evaluationOptimalDesign ) )
  showFIM( fimOptimalDesign )

  invisible( pfimproject )
}

# ==============================================================================
#' getFisherMatrix: return the three FIM components for the optimal design
#' @name getFisherMatrix
#' @param pfimproject An \code{Optimization} object.
#' @return List with fisherMatrix, fixedEffects, varianceEffects.
#' @export
# ==============================================================================

method( getFisherMatrix, Optimization ) = function( pfimproject ) {
  fim = .getOptimalFim( pfimproject )
  list(
    fisherMatrix    = prop( fim, "fisherMatrix"    ),
    fixedEffects    = prop( fim, "fixedEffects"    ),
    varianceEffects = prop( fim, "varianceEffects" )
  )
}

# ==============================================================================
#' getSE: return the SE data frame for the optimal design
#' @name getSE
#' @param pfimproject An \code{Optimization} object.
#' @export
# ==============================================================================

method( getSE, Optimization ) = function( pfimproject ) {
  prop( .getOptimalFim( pfimproject ), "SEAndRSE" )$SE
}

# ==============================================================================
#' getRSE: return the RSE data frame for the optimal design
#' @name getRSE
#' @param pfimproject An \code{Optimization} object.
#' @export
# ==============================================================================

method( getRSE, Optimization ) = function( pfimproject ) {
  prop( .getOptimalFim( pfimproject ), "SEAndRSE" )$RSE
}

# ==============================================================================
#' getShrinkage: return the shrinkage of the FIM
#' @name getShrinkage
#' @param pfimproject An \code{Optimization} object.
#' @export
# ==============================================================================

method( getShrinkage, Optimization ) = function( pfimproject ) {
  prop( .getOptimalFim( pfimproject ), "shrinkage" )
}

# ==============================================================================
#' getDeterminant: return the determinant of the FIM
#' @name getDeterminant
#' @param pfimproject An \code{Optimization} object.
#' @export
# ==============================================================================

method( getDeterminant, Optimization ) = function( pfimproject ) {
  det( prop( .getOptimalFim( pfimproject ), "fisherMatrix" ) )
}

# ==============================================================================
#' getCorrelationMatrix: return the correlation matrix of the FIM
#' @name getCorrelationMatrix
#' @param pfimproject An \code{Optimization} object.
#' @export
# ==============================================================================

method( getCorrelationMatrix, Optimization ) = function( pfimproject ) {
  cor( prop( .getOptimalFim( pfimproject ), "fisherMatrix" ) )
}

# ==============================================================================
#' getDcriterion: return the D-criterion of the FIM
#' @name getDcriterion
#' @param pfimproject An \code{Optimization} object.
#' @export
# ==============================================================================

method( getDcriterion, Optimization ) = function( pfimproject ) {
  Dcriterion( .getOptimalFim( pfimproject ) )
}

# ==============================================================================
#' plotSensitivityIndices: plot sensitivity indices for the optimal design
#' @name plotSensitivityIndices
#' @param pfimproject An \code{Optimization} object.
#' @export
# ==============================================================================

method( plotSensitivityIndices, Optimization ) = function( pfimproject ) {
  plotSensitivityIndices( .getOptimalEval( pfimproject ) )
}

# ==============================================================================
#' plotSE: barplot of SE for the optimal design
#' @name plotSE
#' @param pfimproject An \code{Optimization} object.
#' @export
# ==============================================================================

method( plotSE, Optimization ) = function( pfimproject ) {
  plotSEFIM( prop( pfimproject, "fim" ), .getOptimalEval( pfimproject ) )
}

# ==============================================================================
#' plotRSE: barplot of RSE for the optimal design
#' @name plotRSE
#' @param pfimproject An \code{Optimization} object.
#' @export
# ==============================================================================

method( plotRSE, Optimization ) = function( pfimproject ) {
  plotRSEFIM( prop( pfimproject, "fim" ), .getOptimalEval( pfimproject ) )
}

# ==============================================================================
#' plotWeights: plot algorithm weights (MultiplicativeAlgorithm)
#' @name plotWeights
#' @param optimization An \code{Optimization} object.
#' @export
# ==============================================================================

method( plotWeights, Optimization ) = function( optimization ) {
  plotWeightsMultiplicativeAlgorithm(
    optimization,
    prop( optimization, "optimisationAlgorithmOutputs" )$optimizationAlgorithm
  )
}

# ==============================================================================
#' plotFrequencies: plot optimal frequencies (FedorovWynnAlgorithm)
#' @name plotFrequencies
#' @param optimization An \code{Optimization} object.
#' @export
# ==============================================================================

method( plotFrequencies, Optimization ) = function( optimization ) {
  plotFrequenciesFedorovWynnAlgorithm(
    optimization,
    prop( optimization, "optimisationAlgorithmOutputs" )$optimizationAlgorithm
  )
}

# ==============================================================================
#' Report: generate the HTML optimization report
#' @name Report
#' @param pfimproject  An \code{Optimization} object.
#' @param outputPath   Directory for the rendered HTML file.
#' @param outputFile   File name for the rendered HTML file.
#' @param plotOptions  Plot options passed to plotEvaluation.
#' @export
# ==============================================================================

method( Report, Optimization ) = function( pfimproject, outputPath, outputFile, plotOptions ) {

  projectName = prop( pfimproject, "name" )

  optimisationDesign      = prop( pfimproject, "optimisationDesign" )
  evaluationInitialDesign = optimisationDesign$evaluationInitialDesign
  evaluationOptimalDesign = optimisationDesign$evaluationOptimalDesign
  evaluationOutputs       = prop( pfimproject, "outputs" )

  # ── Model ─────────────────────────────────────────────────────────────────────
  model          = defineModelType( evaluationInitialDesign )
  modelEquations = prop( model, "modelEquations" )

  # ── Model error table ─────────────────────────────────────────────────────────
  modelErrorData = prop( evaluationInitialDesign, "modelError" ) |>
    map( getModelErrorData ) |>
    map( ~ as.data.frame( .x, stringsAsFactors = FALSE ) ) |>
    list_rbind()
  colnames( modelErrorData ) = c( "Output", "Type", "$\\sigma_{slope}$", "$\\sigma_{inter}$" )
  modelErrorTable = kbl( modelErrorData, align = "c" ) |>
    kable_styling( bootstrap_options = "hover", full_width = FALSE,
                   position = "center", font_size = 13 )

  # ── Model parameters table ────────────────────────────────────────────────────
  modelParametersData = prop( evaluationInitialDesign, "modelParameters" ) |>
    map( getModelParametersData ) |>
    map( ~ as.data.frame( .x, stringsAsFactors = FALSE ) ) |>
    list_rbind()
  colnames( modelParametersData ) = c( "Parameter", "$\\mu$", "$\\omega$", "Distribution",
                                       paste0( "$\\mu$", " fixed" ),
                                       paste0( "$\\omega$", " fixed" ) )
  modelParametersTable = kbl( modelParametersData, align = c( "l","l","l","c","c","c" ) ) |>
    kable_styling( bootstrap_options = "hover", full_width = FALSE,
                   position = "center", font_size = 13 )

  # ── Arm helpers ───────────────────────────────────────────────────────────────
  .buildArmTable = function( evaluation, colnamesVec ) {
    designs  = prop( evaluation, "designs" )
    armsData = list_flatten( map( pluck( map( designs, ~ prop( .x, "arms" ) ), 1L ),
                                  getArmData ) )
    df = map( armsData, ~ as.data.frame( .x, stringsAsFactors = FALSE ) ) |> list_rbind()
    colnames( df ) = colnamesVec
    df
  }

  # ── Administration table ──────────────────────────────────────────────────────
  administrationData = prop( evaluationInitialDesign, "designs" ) |>
    map( ~ prop( .x, "arms" ) ) |>
    pluck( 1L ) |>
    map( armAdministration ) |>
    list_flatten() |>
    map( ~ as.data.frame( .x, stringsAsFactors = FALSE ) ) |>
    list_rbind()
  colnames( administrationData ) = c( "Design name", "Arms name", "Number of subject",
                                      "Outcome", "Dose", "Time dose",
                                      "$\\tau$", "$T_{inf}$" )
  administrationTable = kbl( administrationData, align = c( "l","l","l","c","c","c","c" ) ) |>
    kable_styling( bootstrap_options = "hover", full_width = FALSE,
                   position = "center", font_size = 13 )

  # ── Initial design table ──────────────────────────────────────────────────────
  armCols        = c( "Arms name", "Number of subjects", "Outcome", "Dose", "Sampling times" )
  initialDesignData = .buildArmTable( evaluationInitialDesign, armCols )
  initialDesignTable = kbl( initialDesignData, align = c( "l","c","c","c" ) ) |>
    kable_styling( bootstrap_options = "hover", full_width = FALSE,
                   position = "center", font_size = 13 )

  # ── FIM — initial design ──────────────────────────────────────────────────────
  fimInitial         = setEvaluationFim( prop( evaluationInitialDesign, "fim" ), evaluationInitialDesign )
  fimInitialDesignTable = tablesForReport( fimInitial, evaluationInitialDesign )

  # ── Design constraints table ──────────────────────────────────────────────────
  # FIX: renamed from constraintsTableForReport to constraintsData to avoid
  # shadowing the constraintsTableForReport() generic.
  optimAlgoOutputs  = prop( pfimproject, "optimisationAlgorithmOutputs" )
  optimizationAlgorithm = optimAlgoOutputs$optimizationAlgorithm
  constraintsData   = constraintsTableForReport(
    optimizationAlgorithm,
    map( prop( evaluationInitialDesign, "designs" ), ~ prop( .x, "arms" ) )
  )

  # ── Optimal design table ──────────────────────────────────────────────────────
  optimalDesignData  = .buildArmTable( evaluationOptimalDesign, armCols )
  optimalDesignTable = kbl( optimalDesignData, align = c( "l","c","c","c","c","c" ) ) |>
    kable_styling( bootstrap_options = "hover", full_width = FALSE,
                   position = "center", font_size = 13 )

  # ── FIM — optimal design ──────────────────────────────────────────────────────
  fimOptimal     = setEvaluationFim( prop( evaluationOptimalDesign, "fim" ), evaluationOptimalDesign )
  fimOptimalTable = tablesForReport( fimOptimal, evaluationOptimalDesign )

  # ── Plots ─────────────────────────────────────────────────────────────────────
  # FIX: renamed local variables to avoid shadowing plotSE() and plotRSE() generics.
  plotsEvaluationData       = plotEvaluation(          evaluationOptimalDesign, plotOptions )
  plotSensitivityIndicesData = plotSensitivityIndices( evaluationOptimalDesign, plotOptions )
  plotSEData                 = plotSEFIM(  fimOptimal, evaluationOptimalDesign )
  plotRSEData                = plotRSEFIM( fimOptimal, evaluationOptimalDesign )

  # ── Assemble report data ──────────────────────────────────────────────────────
  # FIX: renamed from tablesForReport to reportTables to avoid shadowing the
  # tablesForReport() generic; also added missing outputFile / outputPath.
  reportTables = list(
    evaluationOutputs         = evaluationOutputs,
    modelEquations            = modelEquations,
    modelErrorTable           = modelErrorTable,
    modelParametersTable      = modelParametersTable,
    administrationTable       = administrationTable,
    initialDesignTable        = initialDesignTable,
    constraintsTableForReport = constraintsData,
    fimInitialDesignTable     = fimInitialDesignTable,
    optimalDesignTable        = optimalDesignTable,
    fimOptimalTable           = fimOptimalTable,
    plotsEvaluation           = plotsEvaluationData,
    plotSensitivityIndices    = plotSensitivityIndicesData,
    plotSE                    = plotSEData,
    plotRSE                   = plotRSEData,
    fim                       = fimOptimal,
    pfimproject               = pfimproject,
    projectName               = projectName
  )

  generateReportOptimization( fimOptimal, optimizationAlgorithm,
                              reportTables, outputFile, outputPath )
}
