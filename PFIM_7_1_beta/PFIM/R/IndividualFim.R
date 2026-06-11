# ── Package-level gradient helper ─────────────────────────────────────────────
# Build a (nTimes x nParams) gradient matrix from prop(arm, "evaluationGradients").
# Structure for simple models: named list by output, each element a data.frame
# with one column per parameter.  rbind across outputs stacks the time points.
.gradientMatrix = function( arm, cols ) {
  raw = prop( arm, "evaluationGradients" )
  # Simple path: list(outputName = data.frame(param1, param2, ...))
  # Complex path with covariates: handled by aggregateGradientsWithCovariates in Model.R
  df = if ( is.data.frame( raw ) ) raw
  else do.call( rbind, raw )   # rbind across outputs
  df[ , cols, drop = FALSE ] |> as.matrix()
}

#' @description
#' The class \code{IndividualFim} represents and stores information for the Individual FIM.
#' @title IndividualFim
#' @inheritParams Fim
#' @include Fim.R
#' @export

IndividualFim = new_class( "IndividualFim", package = "PFIM", parent = Fim )

# ==============================================================================
#' evaluateVarianceFIM: compute V and MFVar for the individual FIM
#' @name evaluateVarianceFIM
#' @param fim   An object \code{IndividualFim}.
#' @param model An object \code{Model}.
#' @param arm   An object \code{Arm}.
#' @return List with \code{MFVar} and \code{V}.
#' @export
# ==============================================================================

method( evaluateVarianceFIM, list( IndividualFim, Model, Arm ) ) = function( fim, model, arm ) {

  parameterNames    = map_chr( prop( model, "modelParameters" ), ~ prop( .x, "name" ) )
  evaluationVariance = prop( arm, "evaluationVariance" )

  gradient = .gradientMatrix( arm, parameterNames )

  # FIX: was bdiag(evaluationVariance$errorVariance) then V = bdiag(V) — both calls
  # to bdiag() on a single Matrix object cause silent corruption (bdiag coerces via
  # as.list() producing element-wise scalars instead of one block).
  # as.matrix() safely converts any Matrix/matrix to a plain dense matrix.
  V                = as.matrix( evaluationVariance$errorVariance )
  sigmaDerivatives = evaluationVariance$sigmaDerivatives

  chol2invV = chol2inv( chol( V ) )
  n         = length( sigmaDerivatives )

  # Sigma-Fisher matrix: (i,j) entry = 0.5 * Tr(V^{-1} * dV/d(sigma_i) * V^{-1} * dV/d(sigma_j))
  # FIX: was matrix(map2_dbl(...)) %>% bdiag() — the trailing bdiag() on a single base-R
  # matrix corrupted the result (same as above); removed and return plain matrix.
  MFVar_mat = matrix(
    map_dbl( seq_len( n * n ), function( k ) {
      i = ( k - 1L ) %% n + 1L
      j = ( k - 1L ) %/% n + 1L
      0.5 * sum( diag(
        chol2invV %*% sigmaDerivatives[[ i ]] %*% chol2invV %*% sigmaDerivatives[[ j ]]
      ))
    }),
    nrow = n
  )

  list( MFVar = MFVar_mat, V = V )
}

# ==============================================================================
#' evaluateFim: compute the individual FIM for one arm
#' @name evaluateFim
#' @param fim   An object \code{IndividualFim}.
#' @param model An object \code{Model}.
#' @param arm   An object \code{Arm}.
#' @return The \code{IndividualFim} with fisherMatrix updated.
#' @export
# ==============================================================================

method( evaluateFim, list( IndividualFim, Model, Arm ) ) = function( fim, model, arm ) {

  # FIX: evaluateVarianceFIM is overridden as a plain function in PopulationFim.R
  # which masks the S7 generic and returns list(MFbeta, MFVar) without $V.
  # Call the S7 method explicitly via the method dispatch trick, or better:
  # recompute V directly here to guarantee correctness.
  evaluationVariance = prop( arm, "evaluationVariance" )
  V                  = as.matrix( evaluationVariance$errorVariance )
  sigmaDerivatives   = evaluationVariance$sigmaDerivatives
  chol2invV          = chol2inv( chol( V ) )
  n                  = length( sigmaDerivatives )
  MFVar              = matrix(
    map_dbl( seq_len( n * n ), function( k ) {
      i = ( k - 1L ) %% n + 1L
      j = ( k - 1L ) %/% n + 1L
      0.5 * sum( diag( chol2invV %*% sigmaDerivatives[[ i ]] %*% chol2invV %*% sigmaDerivatives[[ j ]] ) )
    }),
    nrow = n
  )

  parameters = prop( model, "modelParameters" )

  # Non-fixed mu parameters only (columns of the mu block of the FIM)
  paramNamesNonFixed = parameters |>
    keep( ~ !isTRUE( prop( .x, "fixedMu" ) ) ) |>
    map_chr( ~ prop( .x, "name" ) )

  gradients = .gradientMatrix( arm, paramNamesNonFixed )

  MFbeta = t( gradients ) %*% chol2invV %*% gradients

  prop( fim, "fisherMatrix" ) = as.matrix( bdiag( MFbeta, MFVar ) )

  fim
}

# ==============================================================================
#' setOptimalArms: MultiplicativeAlgorithm
#' @name setOptimalArms
#' @param fim                   An object \code{IndividualFim}.
#' @param optimizationAlgorithm An object \code{MultiplicativeAlgorithm}.
#' @return List of optimal arms.
#' @export
# ==============================================================================

method( setOptimalArms, list( IndividualFim, MultiplicativeAlgorithm ) ) =
  function( fim, optimizationAlgorithm ) {

    out          = prop( optimizationAlgorithm, "multiplicativeAlgorithmOutputs" )
    armFims      = out$armFims
    weights      = out$multiplicativeAlgorithmOutput[[ "weights" ]]
    weightsIndex = which( weights > out$weightThreshold )

    armList = map( weightsIndex, function( idx ) {
      arm = pluck( armFims[[ idx ]], 1L )
      prop( arm, "size" ) = 1.0
      prop( arm, "name" ) = paste0( "Arm", idx )
      arm
    })

    armList[ rev( order( map_dbl( armList, ~ prop( .x, "size" ) ) ) ) ]
  }

# ==============================================================================
#' setOptimalArms: FedorovWynnAlgorithm
#' @name setOptimalArms
#' @param fim                   An object \code{IndividualFim}.
#' @param optimizationAlgorithm An object \code{FedorovWynnAlgorithm}.
#' @return List of optimal arms.
#' @export
# ==============================================================================

method( setOptimalArms, list( IndividualFim, FedorovWynnAlgorithm ) ) =
  function( fim, optimizationAlgorithm ) {
    out = prop( optimizationAlgorithm, "FedorovWynnAlgorithmOutputs" )
    imap( out$listArms, ~ {
      prop( .x$arm, "name" ) = paste0( "Arm", .y )
      prop( .x$arm, "size" ) = 1.0
      .x
    })
  }

# ==============================================================================
#' setEvaluationFim: populate FIM result slots
#' @name setEvaluationFim
#' @param fim        An object \code{IndividualFim}.
#' @param evaluation An object \code{Evaluation}.
#' @return The updated \code{IndividualFim}.
#' @export
# ==============================================================================

method( setEvaluationFim, IndividualFim ) = function( fim, evaluation ) {

  parameters  = prop( evaluation, "modelParameters" )
  modelError  = prop( evaluation, "modelError"       )
  greek       = .GREEK_CONSOLE

  # ── Column names for mu block ──────────────────────────────────────────────
  # FIX: was .x@distribution@mu — S4 accessor on S7 object
  columnNamesMu = parameters |>
    keep( ~ !isTRUE( prop( .x, "fixedMu" ) ) ) |>
    keep( ~  prop( prop( .x, "distribution" ), "mu" ) != 0 ) |>
    map_chr( ~ prop( .x, "name" ) ) |>
    map_chr( ~ paste0( greek[ "mu" ], .x ) )

  # ── Column names for sigma block ───────────────────────────────────────────
  columnNamesSigma = modelError |>
    map( function( err ) {
      out = prop( err, "output" )
      c(
        if ( prop( err, "sigmaInter" ) != 0 && !prop( err, "sigmaInterFixed" ) )
          paste0( greek[ "sigma" ], "_inter_", out ),
        if ( prop( err, "sigmaSlope"  ) != 0 && !prop( err, "sigmaSlopeFixed"  ) )
          paste0( greek[ "sigma" ], "_slope_", out )
      )
    }) |>
    unlist( use.names = FALSE )

  # ── Parameter values ───────────────────────────────────────────────────────
  muValues = parameters |>
    keep( ~ !isTRUE( prop( .x, "fixedMu" ) ) ) |>
    map_dbl( ~ prop( prop( .x, "distribution" ), "mu" ) )

  sigmaValues = modelError |>
    map( function( err ) {
      vals = c()
      if ( prop( err, "sigmaInter" ) != 0 && !prop( err, "sigmaInterFixed" ) )
        vals = c( vals, sigmaInter = prop( err, "sigmaInter" ) )
      if ( prop( err, "sigmaSlope"  ) != 0 && !prop( err, "sigmaSlopeFixed"  ) )
        vals = c( vals, sigmaSlope = prop( err, "sigmaSlope" ) )
      vals
    }) |>
    unlist( use.names = FALSE )

  # ── FIM labelling and decomposition ───────────────────────────────────────
  fisherMatrix = prop( fim, "fisherMatrix" )
  allNames     = c( columnNamesMu, columnNamesSigma )
  colnames( fisherMatrix ) = allNames
  rownames( fisherMatrix ) = allNames

  fixedEffects    = fisherMatrix[ columnNamesMu,    columnNamesMu,    drop = FALSE ]
  varianceEffects = fisherMatrix[ columnNamesSigma, columnNamesSigma, drop = FALSE ]

  # ── SE / RSE ───────────────────────────────────────────────────────────────
  SE    = sqrt( diag( chol2inv( chol( fisherMatrix ) ) ) )
  pVals = c( muValues, sigmaValues )
  RSE   = SE / pVals * 100

  seDF = data.frame( parametersValues = pVals, SE = SE, RSE = RSE )
  rownames( seDF ) = allNames

  prop( fim, "fisherMatrix"              ) = fisherMatrix
  prop( fim, "fixedEffects"              ) = fixedEffects
  prop( fim, "varianceEffects"           ) = varianceEffects
  prop( fim, "condNumberFixedEffects"    ) = cond( fixedEffects    )
  prop( fim, "condNumberVarianceEffects" ) = cond( varianceEffects )
  prop( fim, "SEAndRSE"                  ) = list(
    SE       = data.frame( parametersValues = pVals, SE  = SE  ),
    RSE      = data.frame( parametersValues = pVals, RSE = RSE ),
    SEAndRSE = seDF
  )

  fim
}

# ==============================================================================
#' showFIM: print the individual FIM to the console
#' @name showFIM
#' @param fim An object \code{IndividualFim}.
#' @export
# ==============================================================================

method( showFIM, IndividualFim ) = function( fim ) {

  SEAndRSE                  = prop( fim, "SEAndRSE" )
  fisherMatrix              = prop( fim, "fisherMatrix"   )
  fixedEffects              = prop( fim, "fixedEffects"   )
  varianceEffects           = prop( fim, "varianceEffects")
  condNumberFixedEffects    = prop( fim, "condNumberFixedEffects"    )
  condNumberVarianceEffects = prop( fim, "condNumberVarianceEffects" )
  dcrit                     = Dcriterion( fim )

  cat( "\n*************************************** \n Individual Fisher Matrix \n*************************************** \n\n" )
  print( fisherMatrix )
  cat( "\n*************************************** \n Fixed effects (\u03bc) \n*************************************** \n\n" )
  print( fixedEffects )
  cat( "\n*************************************** \n Variance components (\u03c3) \n*************************************** \n\n" )
  print( varianceEffects )
  cat( "\n*********************************************** \n Determinant, condition numbers and D-criterion \n*********************************************** \n\n" )
  cat( c( "Determinant:",  as.numeric( det( fisherMatrix ) ) ), "\n" )
  cat( c( "D-criterion:",  as.numeric( dcrit              ) ), "\n" )
  cat( c( "Conditional number (fixed effects):",    as.numeric( condNumberFixedEffects    ), "\n" ) )
  cat( c( "Conditional number (variance effects):", as.numeric( condNumberVarianceEffects ), "\n" ) )
  cat( "\n*************************************** \n Parameters estimation \n*************************************** \n\n" )
  print( SEAndRSE$SEAndRSE )

  invisible( fim )
}

# ==============================================================================
#' @title Plot SE or RSE bar chart (Individual FIM — internal helper)
#' @description Internal helper shared by \code{plotSEFIM} and
#'   \code{plotRSEFIM}. Not exported.
#' @param fim        An object of class \code{IndividualFim}.
#' @param evaluation An object of class \code{PFIMProject}.
#' @param metric     Character string, either \code{"SE"} or \code{"RSE"}.
#' @return A \code{ggplot} bar chart.
#' @keywords internal
#' @noRd
.individualSEPlot = function( fim, evaluation, metric ) {
  parameters  = prop( evaluation, "modelParameters" )
  modelError  = prop( evaluation, "modelError"       )
  fim         = setEvaluationFim( prop( evaluation, "fim" ), evaluation )
  se          = prop( fim, "SEAndRSE" )
  greek       = .GREEK_PLOT

  paramsMu = parameters |>
    keep( ~ !isTRUE( prop( .x, "fixedMu" ) ) ) |>
    keep( ~ prop( prop( .x, "distribution" ), "mu" ) != 0 ) |>
    map_chr( ~ prop( .x, "name" ) )

  paramsSigma = modelError |>
    map( function( err ) {
      out = prop( err, "output" )
      c(
        if ( prop( err, "sigmaInter" ) != 0 && !prop( err, "sigmaInterFixed" ) )
          paste0( greek[ "sigma" ], "_inter_", out ),
        if ( prop( err, "sigmaSlope"  ) != 0 && !prop( err, "sigmaSlopeFixed"  ) )
          paste0( greek[ "sigma" ], "_slope_", out )
      )
    }) |> unlist( use.names = FALSE )

  yVals = if ( metric == "SE" ) se$SE$SE else se$RSE$RSE
  cats  = paste0( metric, " ",
                  c( rep( greek[ "mu" ], length( paramsMu ) ),
                     rep( greek[ "sigma" ], length( paramsSigma ) ) ) )

  data = data.frame(
    Parameter = c( paramsMu, paramsSigma ),
    y         = yVals,
    cat       = cats
  )
  names( data )[ 2L ] = metric

  ggplot( data, aes( x = Parameter, y = .data[[ metric ]] ) ) +
    geom_bar( stat = "identity", show.legend = FALSE ) +
    facet_wrap( ~factor( cat, levels = unique( cats ) ), scales = "free_x" ) +
    theme( legend.position = "none",
           plot.title  = element_text( size = 16, hjust = 0.5 ),
           axis.text.x = element_text( size = 16, angle = 90, vjust = 0.5 ) )
}

method( plotSEFIM,  list( IndividualFim, PFIMProject ) ) = function( fim, evaluation )
  .individualSEPlot( fim, evaluation, "SE"  )
method( plotRSEFIM, list( IndividualFim, PFIMProject ) ) = function( fim, evaluation )
  .individualSEPlot( fim, evaluation, "RSE" )

# ==============================================================================
#' tablesForReport: kableExtra tables for the individual FIM report
#' @name tablesForReport
#' @param fim        An object \code{IndividualFim}.
#' @param evaluation An object \code{Evaluation}.
#' @return List: fixedEffectsTable, varianceEffectsTable, FIMCriteriaTable, SEAndRSETable.
#' @export
# ==============================================================================

method( tablesForReport, list( IndividualFim, PFIMProject ) ) = function( fim, evaluation ) {

  parameters                = prop( evaluation, "modelParameters" )
  modelError                = prop( evaluation, "modelError"       )
  SEAndRSE                  = prop( fim, "SEAndRSE"               )$SEAndRSE
  fisherMatrix              = prop( fim, "fisherMatrix"            )
  fixedEffects              = as.matrix( prop( fim, "fixedEffects"    ) )
  varianceEffects           = as.matrix( prop( fim, "varianceEffects" ) )
  condNumberFixedEffects    = prop( fim, "condNumberFixedEffects"    )
  condNumberVarianceEffects = prop( fim, "condNumberVarianceEffects" )
  greek                     = .GREEK_LATEX

  columnNamesMu = parameters |>
    keep( ~ !isTRUE( prop( .x, "fixedMu" ) ) ) |>
    keep( ~  prop( prop( .x, "distribution" ), "mu" ) != 0 ) |>
    map_chr( ~ prop( .x, "name" ) ) |>
    map_chr( ~ paste0( greek[ "mu" ], .x, "}$" ) )

  columnNamesSigma = modelError |>
    map( function( err ) {
      out = prop( err, "output" )
      c(
        if ( prop( err, "sigmaInter" ) != 0 && !prop( err, "sigmaInterFixed" ) )
          paste0( greek[ "sigma" ], "{inter}}_{", out, "}$" ),
        if ( prop( err, "sigmaSlope"  ) != 0 && !prop( err, "sigmaSlopeFixed"  ) )
          paste0( greek[ "sigma" ], "{slope}}_{", out, "}$" )
      )
    }) |> unlist( use.names = FALSE )

  colnames( fixedEffects    ) = columnNamesMu
  rownames( fixedEffects    ) = columnNamesMu
  colnames( varianceEffects ) = columnNamesSigma
  rownames( varianceEffects ) = columnNamesSigma

  .kbl = function( df )
    kbl( df ) |>
    kable_styling( bootstrap_options = "hover", full_width = FALSE,
                   position = "center", font_size = 13 )

  FIMCriteriaTable = data.frame(
    determinant               = det( fisherMatrix ),
    dcriterion                = Dcriterion( fim ),
    condNumberFixedEffects    = condNumberFixedEffects,
    condNumberVarianceEffects = condNumberVarianceEffects
  ) |>
    kbl( col.names = c( "", "", "Fixed effects", "Variance effects" ),
         align = "c", format = "html" ) |>
    add_header_above( c( "Determinant" = 1, "D-criterion" = 1, "Condition number" = 2 ) ) |>
    kable_styling( bootstrap_options = "hover", full_width = FALSE,
                   position = "center", font_size = 13 )

  SEAndRSETable = data.frame(
    c( columnNamesMu, columnNamesSigma ), round( SEAndRSE, 3 )
  ) |>
    (\( df ) { row.names( df ) = NULL; df })() |>
    kbl( col.names = c( "Parameters", "Parameter values", "SE", "RSE (%)" ),
         align = "c" ) |>
    kable_styling( bootstrap_options = "hover", full_width = FALSE,
                   position = "center", font_size = 13 )

  list(
    fixedEffectsTable    = .kbl( fixedEffects    ),
    varianceEffectsTable = .kbl( varianceEffects ),
    FIMCriteriaTable     = FIMCriteriaTable,
    SEAndRSETable        = SEAndRSETable
  )
}

# ==============================================================================
# Report rendering methods
# ==============================================================================

.reportIndividualPath = function( filename )
  file.path( system.file( package = "PFIM" ),
             "rmarkdown", "templates", "skeleton", filename )

#' @name generateReportEvaluation
#' @export

method( generateReportEvaluation, IndividualFim ) =
  function( fim, tablesForReport, outputFile, outputPath ) {
    rmarkdown::render(
      input       = .reportIndividualPath( "EvaluationIndividualFIM.rmd" ),
      output_file = outputFile, output_dir = outputPath,
      params      = list( tablesForReport = tablesForReport )
    )
  }


#' @name generateReportOptimization
#' @export
#'
method( generateReportOptimization, list( IndividualFim, MultiplicativeAlgorithm ) ) =
  function( fim, optimizationAlgorithm, tablesForReport, outputFile, outputPath ) {
    rmarkdown::render(
      input       = .reportIndividualPath( "OptimizationMultiplicativeAlgorithmIndividualFim.rmd" ),
      output_file = outputFile, output_dir = outputPath,
      params      = list( tablesForReport = tablesForReport )
    )
  }
#' @name generateReportOptimization
#' @export
#'
method( generateReportOptimization, list( IndividualFim, FedorovWynnAlgorithm ) ) =
  function( fim, optimizationAlgorithm, tablesForReport, outputFile, outputPath ) {
    rmarkdown::render(
      input       = .reportIndividualPath( "OptimizationFedorovWynnAlgorithmIndividualFim.rmd" ),
      output_file = outputFile, output_dir = outputPath,
      params      = list( tablesForReport = tablesForReport )
    )
  }
#' @name generateReportOptimization
#' @export
#'
method( generateReportOptimization, list( IndividualFim, SimplexAlgorithm ) ) =
  function( fim, optimizationAlgorithm, tablesForReport, outputFile, outputPath ) {
    rmarkdown::render(
      input       = .reportIndividualPath( "OptimizationSimplexAlgorithmIndividualFim.rmd" ),
      output_file = outputFile, output_dir = outputPath,
      params      = list( tablesForReport = tablesForReport )
    )
  }
#' @name generateReportOptimization
#' @export
#'
method( generateReportOptimization, list( IndividualFim, PSOAlgorithm ) ) =
  function( fim, optimizationAlgorithm, tablesForReport, outputFile, outputPath ) {
    rmarkdown::render(
      input       = .reportIndividualPath( "OptimizationPSOAlgorithmIndividualFim.rmd" ),
      output_file = outputFile, output_dir = outputPath,
      params      = list( tablesForReport = tablesForReport )
    )
  }
#' @name generateReportOptimization
#' @export
#'
method( generateReportOptimization, list( IndividualFim, PGBOAlgorithm ) ) =
  function( fim, optimizationAlgorithm, tablesForReport, outputFile, outputPath ) {
    rmarkdown::render(
      input       = .reportIndividualPath( "OptimizationPGBOAlgorithmIndividualFim.rmd" ),
      output_file = outputFile, output_dir = outputPath,
      params      = list( tablesForReport = tablesForReport )
    )
  }
