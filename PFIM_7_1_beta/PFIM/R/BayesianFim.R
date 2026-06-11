# ── Package-level gradient helper ─────────────────────────────────────────────
# Build a (nTimes x nParams) gradient matrix from prop(arm, "evaluationGradients").
# Structure: named list by output, each element a data.frame with one column per
# parameter. do.call(rbind, raw) stacks all output data.frames by time points.
.gradientMatrix = function( arm, cols ) {
  raw = prop( arm, "evaluationGradients" )
  df  = if ( is.data.frame( raw ) ) raw else do.call( rbind, raw )
  df[ , cols, drop = FALSE ] |> as.matrix()
}

#' @description
#' The class \code{BayesianFim} represents and stores information for the Bayesian FIM.
#' @title BayesianFim
#' @param fisherMatrix             A matrix giving the numerical values of the FIM.
#' @param shrinkage                A vector giving the shrinkage values.
#' @param fixedEffects             A matrix giving the fixed-effects block of the FIM.
#' @param varianceEffects          A matrix giving the variance-effects block of the FIM.
#' @param SEAndRSE                 A data frame of SE and RSE values.
#' @param condNumberFixedEffects   Condition number of the fixed-effects block.
#' @param condNumberVarianceEffects Condition number of the variance-effects block.
#' @include Fim.R
#' @include MultiplicativeAlgorithm.R
#' @include FedorovWynnAlgorithm.R
#' @export

BayesianFim = new_class( "BayesianFim",
                         package    = "PFIM",
                         parent     = Fim,
                         properties = list(
                           fisherMatrix              = new_property( class_double, default = 0.0   ),
                           fixedEffects              = new_property( class_double, default = 0.0   ),
                           varianceEffects           = new_property( class_double, default = 0.0   ),
                           SEAndRSE                  = new_property( class_list,   default = list() ),
                           condNumberFixedEffects    = new_property( class_double, default = 0.0   ),
                           condNumberVarianceEffects = new_property( class_double, default = 0.0   ),
                           shrinkage                 = new_property( class_double, default = 0.0   )
                         )
)

# ── Private helper ─────────────────────────────────────────────────────────────
# Filter estimable (non-fixed, non-zero) parameters and return their names
# with the console Greek-mu prefix — shared by setEvaluationFim, plotSEFIM, etc.
.bayesianEstimableParamNames = function( parameters, greekPrefix ) {
  parameters |>
    keep( ~ !isTRUE( prop( .x, "fixedMu" ) ) ) |>
    keep( ~ prop( prop( .x, "distribution" ), "mu" ) != 0 ) |>
    keep( ~ !isTRUE( prop( .x, "fixedOmega" ) ) ) |>
    keep( ~ prop( prop( .x, "distribution" ), "omega" ) != 0 ) |>
    map_chr( ~ prop( .x, "name" ) ) |>
    map_chr( ~ paste0( greekPrefix, .x ) )
}

# ==============================================================================
#' evaluateVarianceFIM: compute the variance matrix V and MFbeta for the Bayesian FIM
#' @name evaluateVarianceFIM
#' @param fim   An object \code{BayesianFim}.
#' @param model An object \code{Model}.
#' @param arm   An object \code{Arm}.
#' @return List with \code{MFbeta} and \code{V}.
#' @export
# ==============================================================================

method( evaluateVarianceFIM, list( BayesianFim, Model, Arm ) ) = function( fim, model, arm ) {

  parameterNames = map_chr( prop( model, "modelParameters" ), ~ prop( .x, "name" ) )

  # Gradient matrix: (nTimes × nParam) — one column per parameter
  gradient = .gradientMatrix( arm, parameterNames )

  # FIX: was bdiag(evaluationVariance$errorVariance) which fails when errorVariance
  # is already a combined sparse Matrix — bdiag(single_matrix) calls as.list() on
  # it producing a list of scalars, not a list of matrices.
  # as.matrix() safely converts any Matrix/matrix to a plain dense matrix.
  V      = as.matrix( prop( arm, "evaluationVariance" )$errorVariance )
  MFbeta = crossprod( gradient, chol2inv( chol( V ) ) ) %*% gradient

  list( MFbeta = MFbeta, V = V )
}

# ==============================================================================
#' evaluateFim: compute the Bayesian FIM for one arm
#' @name evaluateFim
#' @param fim   An object \code{BayesianFim}.
#' @param model An object \code{Model}.
#' @param arm   An object \code{Arm}.
#' @return The \code{BayesianFim} with fisherMatrix and shrinkage updated.
#' @export
# ==============================================================================

method( evaluateFim, list( BayesianFim, Model, Arm ) ) = function( fim, model, arm ) {

  parameters  = prop( model, "modelParameters" )
  armSize     = prop( arm,   "size"            )

  # Full gradient matrix (all parameters, all outputs stacked)
  paramNamesAll = map_chr( parameters, ~ prop( .x, "name" ) )
  gradient = .gradientMatrix( arm, paramNamesAll )

  # FIX: renamed from evaluateVarianceFIM to varianceResult to avoid shadowing
  # the evaluateVarianceFIM() S7 generic, which caused silent dispatch failure.
  # FIX: evaluateVarianceFIM is overridden as a plain function in PopulationFim.R
  # which masks the S7 generic. Extract V directly from arm slot.
  V           = as.matrix( prop( arm, "evaluationVariance" )$errorVariance )

  MFbeta_full = crossprod( gradient, chol2inv( chol( V ) ) ) %*% gradient

  # ── mu diagonal matrix ───────────────────────────────────────────────────────
  # For log-normal parameters:  d(f)/d(mu) = d(f)/d(theta) * (theta/mu) ≈ 1 * mu_i
  # For normal parameters: mu_i = 1 (identity transformation)
  # FIX: was .x@distribution@mu — S4 accessor on an S7 object
  mu_vec = map_dbl( parameters, function( p ) {
    d = prop( p, "distribution" )
    if ( inherits( d, "PFIM::Normal" ) ) 1 else prop( d, "mu" )
  })
  mu = if ( length( mu_vec ) == 1L ) mu_vec[[ 1L ]] else diag( mu_vec )

  # ── omega diagonal matrix ────────────────────────────────────────────────────
  # FIX: was pluck(.x, "distribution", "omega") — now uses prop() for S7 safety
  omega_vec = map_dbl( parameters, ~ prop( prop( .x, "distribution" ), "omega" ) )
  omega     = if ( length( omega_vec ) == 1L ) omega_vec[[ 1L ]]^2
  else diag( omega_vec^2 )

  # ── Remove fixed / zero-mu / zero-omega parameters ───────────────────────────
  # FIX: was .x@distribution@mu == 0  — S4 accessor on S7
  isFixedMu    = map_lgl( parameters, ~ isTRUE( prop( .x, "fixedMu"    ) ) ||
                            prop( prop( .x, "distribution" ), "mu"    ) == 0 )
  isFixedOmega = map_lgl( parameters, ~ isTRUE( prop( .x, "fixedOmega" ) ) ||
                            prop( prop( .x, "distribution" ), "omega" ) == 0 )
  indexFixed   = unique( which( isFixedMu | isFixedOmega ) )

  if ( length( indexFixed ) > 0L ) {
    mu          = mu[          -indexFixed, -indexFixed, drop = FALSE ]
    omega       = omega[       -indexFixed, -indexFixed, drop = FALSE ]
    MFbeta_full = MFbeta_full[ -indexFixed, -indexFixed, drop = FALSE ]
  }

  # ── Bayesian FIM: M_F_Bayes = mu' * M_F_ind * mu + (mu * omega * mu)^{-1} ─────
  priorVariance = mu %*% omega %*% mu
  MFbeta        = t( mu ) %*% MFbeta_full %*% mu + solve( priorVariance )

  prop( fim, "fisherMatrix" ) = as.matrix( MFbeta )

  # ── Shrinkage ─────────────────────────────────────────────────────────────────
  # FIX: was `* 100 %>% as.vector()` which applied as.vector(100) due to %>%
  # precedence, not as.vector(shrinkage_vector * 100).
  prop( fim, "shrinkage" ) = as.vector(
    diag( chol2inv( chol( MFbeta ) ) %*% chol2inv( chol( priorVariance ) ) ) * 100
  )

  fim
}

# ==============================================================================
#' setOptimalArms: MultiplicativeAlgorithm
#' @name setOptimalArms
#' @param fim                    An object \code{BayesianFim}.
#' @param optimizationAlgorithm  An object \code{MultiplicativeAlgorithm}.
#' @return List of optimal arms.
#' @export
# ==============================================================================

method( setOptimalArms, list( BayesianFim, MultiplicativeAlgorithm ) ) =
  function( fim, optimizationAlgorithm ) {

    out             = prop( optimizationAlgorithm, "multiplicativeAlgorithmOutputs" )
    armFims         = out$armFims
    weights         = out$multiplicativeAlgorithmOutput[[ "weights" ]]
    weightsIndex    = which( weights > out$weightThreshold )

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
#' @param fim                    An object \code{BayesianFim}.
#' @param optimizationAlgorithm  An object \code{FedorovWynnAlgorithm}.
#' @return List of optimal arms.
#' @export
# ==============================================================================

method( setOptimalArms, list( BayesianFim, FedorovWynnAlgorithm ) ) =
  function( fim, optimizationAlgorithm ) {
    out = prop( optimizationAlgorithm, "FedorovWynnAlgorithmOutputs" )
    imap( out$listArms, ~ {
      prop( .x$arm, "name" ) = paste0( "Arm", .y )
      prop( .x$arm, "size" ) = 1.0  # FIX: use double, not integer
      .x
    })
  }

# ==============================================================================
#' setEvaluationFim: populate FIM result slots after optimisation/evaluation
#' @name setEvaluationFim
#' @param fim        An object \code{BayesianFim}.
#' @param evaluation An object \code{Evaluation}.
#' @return The updated \code{BayesianFim}.
#' @export
# ==============================================================================

method( setEvaluationFim, BayesianFim ) = function( fim, evaluation ) {

  parameters  = prop( evaluation, "modelParameters" )
  greek       = .GREEK_CONSOLE

  # Only include estimable (non-fixed, non-zero mu AND omega) parameters
  columnNamesMu = .bayesianEstimableParamNames( parameters, greek[ "mu" ] )

  muValues = parameters |>
    keep( ~ !isTRUE( prop( .x, "fixedMu"    ) ) ) |>
    keep( ~  prop( prop( .x, "distribution" ), "mu"    ) != 0 ) |>
    keep( ~ !isTRUE( prop( .x, "fixedOmega" ) ) ) |>
    keep( ~  prop( prop( .x, "distribution" ), "omega" ) != 0 ) |>
    map_dbl( ~ prop( prop( .x, "distribution" ), "mu" ) )

  fisherMatrix = prop( fim, "fisherMatrix" )
  colnames( fisherMatrix ) = columnNamesMu
  rownames( fisherMatrix ) = columnNamesMu
  fixedEffects = fisherMatrix[ columnNamesMu, columnNamesMu, drop = FALSE ]

  shrinkage = prop( fim, "shrinkage" )
  SE        = sqrt( diag( chol2inv( chol( fisherMatrix ) ) ) )
  RSE       = SE / muValues * 100
  seDF      = data.frame( parametersValues = muValues, SE = SE, RSE = RSE )
  rownames( seDF ) = rownames( fisherMatrix )

  prop( fim, "fisherMatrix"           ) = fisherMatrix
  prop( fim, "fixedEffects"           ) = fixedEffects
  prop( fim, "shrinkage"              ) = t( matrix( shrinkage, nrow = 1L,
                                                     dimnames = list( "Shrinkage", columnNamesMu ) ) )
  prop( fim, "condNumberFixedEffects" ) = cond( fixedEffects )
  prop( fim, "SEAndRSE"               ) = list(
    SE       = data.frame( parametersValues = muValues, SE  = SE  ),
    RSE      = data.frame( parametersValues = muValues, RSE = RSE ),
    SEAndRSE = seDF
  )

  fim
}

# ==============================================================================
#' showFIM: print the Bayesian FIM to the console
#' @name showFIM
#' @param fim An object \code{BayesianFim}.
#' @export
# ==============================================================================

method( showFIM, BayesianFim ) = function( fim ) {

  SEAndRSE               = prop( fim, "SEAndRSE" )
  fisherMatrix           = prop( fim, "fisherMatrix" )
  fixedEffects           = prop( fim, "fixedEffects" )
  shrinkage              = prop( fim, "shrinkage" )
  condNumberFixedEffects = prop( fim, "condNumberFixedEffects" )
  dcrit                  = Dcriterion( fim )

  cat( "\n*************************************** \n Bayesian Fisher Matrix \n*************************************** \n\n" )
  print( fisherMatrix )
  cat( "\n*************************************** \n Fixed effects \n*************************************** \n\n" )
  print( fixedEffects )
  cat( "\n*********************************************** \n Determinant, condition numbers and D-criterion \n*********************************************** \n\n" )
  cat( c( "Determinant:",  as.numeric( det( fisherMatrix ) ) ), "\n" )
  cat( c( "D-criterion:",  as.numeric( dcrit              ) ), "\n" )
  cat( c( "Conditional number of the fixed effects:", as.numeric( condNumberFixedEffects ), "\n" ) )
  cat( "\n*************************************** \n Shrinkage \n*************************************** \n\n" )
  print( shrinkage )
  cat( "\n*************************************** \n Parameters estimation \n*************************************** \n\n" )
  print( SEAndRSE$SEAndRSE )

  invisible( fim )
}

# ==============================================================================
#' plotSEFIM / plotRSEFIM / plotShrinkage: bar charts for SE, RSE, shrinkage
#' @name plotSEFIM
#' @export
# ==============================================================================

method( plotSEFIM, list( BayesianFim, PFIMProject ) ) = function( fim, evaluation ) {
  parameters = prop( evaluation, "modelParameters" )
  fim        = setEvaluationFim( prop( evaluation, "fim" ), evaluation )
  se         = prop( fim, "SEAndRSE" )
  greek      = .GREEK_PLOT

  params  = .bayesianEstimableParamNames( parameters, "" )  # raw names
  data    = data.frame(
    Parameter = params,
    SE        = se$SE$SE,
    cat       = paste0( "SE ", greek[ "mu" ] )
  )
  ggplot( data, aes( x = Parameter, y = SE ) ) +
    geom_bar( stat = "identity", show.legend = FALSE ) +
    facet_wrap( ~factor( cat, levels = paste0( "SE ", greek[ "mu" ] ) ), scales = "free_x" ) +
    theme( legend.position = "none",
           plot.title   = element_text( size = 16, hjust = 0.5 ),
           axis.text.x  = element_text( size = 16, angle = 90, vjust = 0.5 ) )
}

method( plotRSEFIM, list( BayesianFim, PFIMProject ) ) = function( fim, evaluation ) {
  parameters = prop( evaluation, "modelParameters" )
  fim        = setEvaluationFim( prop( evaluation, "fim" ), evaluation )
  se         = prop( fim, "SEAndRSE" )
  greek      = .GREEK_PLOT

  params  = .bayesianEstimableParamNames( parameters, "" )
  data    = data.frame(
    Parameter = params,
    RSE       = se$RSE$RSE,
    cat       = paste0( "RSE ", greek[ "mu" ] )
  )
  ggplot( data, aes( x = Parameter, y = RSE ) ) +
    geom_bar( stat = "identity", show.legend = FALSE ) +
    facet_wrap( ~factor( cat, levels = paste0( "RSE ", greek[ "mu" ] ) ), scales = "free_x" ) +
    theme( legend.position = "none",
           plot.title  = element_text( size = 16, hjust = 0.5 ),
           axis.text.x = element_text( size = 16, angle = 90, vjust = 0.5 ) )
}

method( plotShrinkage, list( BayesianFim, PFIMProject ) ) = function( fim, evaluation ) {
  parameters = prop( evaluation, "modelParameters" )
  shrinkage  = prop( fim, "shrinkage" )
  params     = .bayesianEstimableParamNames( parameters, "" )

  data = data.frame( Parameter = params, Shrinkage = as.vector( shrinkage ) )
  ggplot( data, aes( x = Parameter, y = Shrinkage ) ) +
    geom_bar( stat = "identity", show.legend = FALSE ) +
    theme( legend.position = "none",
           plot.title  = element_text( size = 16, hjust = 0.5 ),
           axis.text.x = element_text( size = 16, angle = 90, vjust = 0.5 ) )
}

# ==============================================================================
#' tablesForReport: kableExtra tables for the Bayesian FIM report
#' @name tablesForReport
#' @param fim        An object \code{BayesianFim}.
#' @param evaluation An object \code{Evaluation}.
#' @return List: fixedEffectsTable, FIMCriteriaTable, SEAndRSETable.
#' @export
# ==============================================================================

method( tablesForReport, list( BayesianFim, PFIMProject ) ) = function( fim, evaluation ) {

  parameters             = prop( evaluation, "modelParameters" )
  SEAndRSE               = prop( fim, "SEAndRSE" )$SEAndRSE
  fisherMatrix           = prop( fim, "fisherMatrix"   )
  fixedEffects           = as.matrix( prop( fim, "fixedEffects" ) )
  shrinkage              = prop( fim, "shrinkage" )
  condNumberFixedEffects = prop( fim, "condNumberFixedEffects" )
  greek                  = .GREEK_LATEX

  columnNamesMu = .bayesianEstimableParamNames( parameters, greek[ "mu" ] ) |>
    map_chr( ~ paste0( .x, "}$" ) )

  colnames( fixedEffects ) = columnNamesMu
  rownames( fixedEffects ) = columnNamesMu

  fixedEffectsTable = fixedEffects |>
    kbl() |>
    kable_styling( bootstrap_options = "hover", full_width = FALSE,
                   position = "center", font_size = 13 )

  FIMCriteriaTable = data.frame(
    determinant  = det( fisherMatrix ),
    dcriterion   = Dcriterion( fim ),
    FixedEffects = condNumberFixedEffects
  ) |>
    kbl( col.names = c( "", "", "Fixed effects" ), align = "c", format = "html" ) |>
    add_header_above( c( "Determinant" = 1, "D-criterion" = 1, "Condition number" = 1 ) ) |>
    kable_styling( bootstrap_options = "hover", full_width = FALSE,
                   position = "center", font_size = 13 )

  SEAndRSETable = data.frame(
    c( columnNamesMu ), round( SEAndRSE, 3 ), as.vector( shrinkage )
  ) |>
    (\( df ) { row.names( df ) = NULL; df })() |>
    kbl( col.names = c( "Parameters", "Parameter values", "SE", "RSE (%)", "Shrinkage" ),
         align = "c" ) |>
    kable_styling( bootstrap_options = "hover", full_width = FALSE,
                   position = "center", font_size = 13 )

  list( fixedEffectsTable = fixedEffectsTable,
        FIMCriteriaTable  = FIMCriteriaTable,
        SEAndRSETable     = SEAndRSETable )
}

# ==============================================================================
# Report rendering methods
# ==============================================================================

.reportBayesianPath = function( filename )
  file.path( system.file( package = "PFIM" ),
             "rmarkdown", "templates", "skeleton", filename )

#' @name generateReportEvaluation
#' @export
#'
method( generateReportEvaluation, BayesianFim ) =
  function( fim, tablesForReport, outputFile, outputPath ) {
    rmarkdown::render(
      input       = .reportBayesianPath( "EvaluationBayesianFim.rmd" ),
      output_file = outputFile,
      output_dir  = outputPath,
      params      = list( tablesForReport = tablesForReport )
    )
  }


#' @name generateReportOptimization
#' @export

method( generateReportOptimization, list( BayesianFim, MultiplicativeAlgorithm ) ) =
  function( fim, optimizationAlgorithm, tablesForReport, outputFile, outputPath ) {
    rmarkdown::render(
      input       = .reportBayesianPath( "OptimizationMultiplicativeAlgorithmBayesianFIM.rmd" ),
      output_file = outputFile, output_dir = outputPath,
      params      = list( tablesForReport = tablesForReport )
    )
  }


#' @name generateReportOptimization
#' @export
#'
method( generateReportOptimization, list( BayesianFim, FedorovWynnAlgorithm ) ) =
  function( fim, optimizationAlgorithm, tablesForReport, outputFile, outputPath ) {
    rmarkdown::render(
      input       = .reportBayesianPath( "OptimizationFedorovWynnAlgorithmBayesianFIM.rmd" ),
      output_file = outputFile, output_dir = outputPath,
      params      = list( tablesForReport = tablesForReport )
    )
  }


#' @name generateReportOptimization
#' @export
#'
method( generateReportOptimization, list( BayesianFim, SimplexAlgorithm ) ) =
  function( fim, optimizationAlgorithm, tablesForReport, outputFile, outputPath ) {
    rmarkdown::render(
      input       = .reportBayesianPath( "OptimizationSimplexAlgorithmBayesianFIM.rmd" ),
      output_file = outputFile, output_dir = outputPath,
      params      = list( tablesForReport = tablesForReport )
    )
  }


#' @name generateReportOptimization
#' @export
#'
method( generateReportOptimization, list( BayesianFim, PSOAlgorithm ) ) =
  function( fim, optimizationAlgorithm, tablesForReport, outputFile, outputPath ) {
    rmarkdown::render(
      input       = .reportBayesianPath( "OptimizationPSOAlgorithmBayesianFIM.rmd" ),
      output_file = outputFile, output_dir = outputPath,
      params      = list( tablesForReport = tablesForReport )
    )
  }


#' @name generateReportOptimization
#' @export
#'
method( generateReportOptimization, list( BayesianFim, PGBOAlgorithm ) ) =
  function( fim, optimizationAlgorithm, tablesForReport, outputFile, outputPath ) {
    rmarkdown::render(
      input       = .reportBayesianPath( "OptimizationPGBOAlgorithmBayesianFIM.rmd" ),
      output_file = outputFile, output_dir = outputPath,
      params      = list( tablesForReport = tablesForReport )
    )
  }
