#' @description The class \code{Model} represents and stores information for a model.
#' @title Model
#' @param name Character vector specifying the model name
#' @param modelParameters List of model parameters
#' @param modelCovariatesEquation modelCovariatesEquation
#' @param omegaWithIOV omegaWithIOV
#' @param modelCovariates List of model covariate.
#' @param covariatesCombination Dataframe giving the combination of the covariates
#' @param covariatesEffect List of the effects of the covariate
#' @param samplings Numeric vector of sampling times
#' @param modelEquations List containing the model equations
#' @param wrapper Function wrapper for the model (default: function () NULL)
#' @param outputFormula List of output formulas
#' @param outputNames Character vector of output names
#' @param variableNames Character vector of variable names
#' @param outcomesWithAdministration Character vector of outcomes with administration
#' @param outcomesWithNoAdministration Character vector of outcomes without administration
#' @param modelError List defining the error model
#' @param odeSolverParameters List of ODE solver parameters
#' @param parametersForComputingGradient List of parameters for gradient computation
#' @param initialConditions Numeric vector of initial conditions
#' @param functionArguments Character vector of function arguments
#' @param functionArgumentsSymbol List of function argument symbols
#' @include CovariateModelEquation.R
#' @export

Model = new_class( "Model", package = "PFIM",

                   properties = list(
                     name                           = new_property( class_character,        default = character(0) ),
                     modelParameters                = new_property( class_list,             default = list()       ),
                     modelCovariatesEquation        = new_property( CovariateModelEquation, default = NULL         ),
                     modelCovariates                = new_property( class_list,             default = list()       ),
                     covariatesEffect               = new_property( class_list,             default = list()       ),
                     covariatesCombination          = new_property( class_list,             default = list()       ),
                     modelParametersWithCovariates  = new_property( class_list,             default = list()       ),
                     numberOfOccasions              = new_property( class_double,           default = 1            ),
                     omegaWithIOV                   = new_property( class_double,           default = numeric(0)   ),
                     samplings                      = new_property( class_numeric,          default = numeric(0)   ),
                     modelEquations                 = new_property( class_list,             default = list()       ),
                     wrapper                        = new_property( class_function,         default = NULL         ),
                     outputFormula                  = new_property( class_list,             default = list()       ),
                     outputNames                    = new_property( class_character,        default = character(0) ),
                     variableNames                  = new_property( class_character,        default = character(0) ),
                     outcomesWithAdministration     = new_property( class_character,        default = character(0) ),
                     outcomesWithNoAdministration   = new_property( class_character,        default = character(0) ),
                     modelError                     = new_property( class_list,             default = list()       ),
                     odeSolverParameters            = new_property( class_list,             default = list()       ),
                     parametersForComputingGradient = new_property( class_list,             default = list()       ),
                     initialConditions              = new_property( class_double,           default = numeric(0)   ),
                     functionArguments              = new_property( class_character,        default = character(0) ),
                     functionArgumentsSymbol        = new_property( class_list,             default = list()       )
                   ) )

defineModelWrapper                = new_generic( "defineModelWrapper",                c( "model" ) )
defineModelAdministration         = new_generic( "defineModelAdministration",         c( "model" ) )
evaluateModel                     = new_generic( "evaluateModel",                     c( "model" ) )
evaluateModelGradient             = new_generic( "evaluateModelGradient",             c( "model" ) )
evaluateModelVariance             = new_generic( "evaluateModelVariance",             c( "model" ) )
evaluateInitialConditions         = new_generic( "evaluateInitialConditions",         c( "model" ) )
finiteDifferenceHessian           = new_generic( "finiteDifferenceHessian",           c( "model" ) )
evaluateCovariatesEffects         = new_generic( "evaluateCovariatesEffects",         c( "model" ) )
evaluateOmegaMatrixFromCovariates = new_generic( "evaluateOmegaMatrixFromCovariates", c( "model" ) )
generateCovariatesCombination     = new_generic( "generateCovariatesCombination",     c( "model" ) )
definePKModel                     = new_generic( "definePKModel",   c( "pkModel", "pfimproject" ) )
definePKPDModel                   = new_generic( "definePKPDModel", c( "pkModel", "pdModel", "pfimproject" ) )
modelParametersWithCovariates     = new_generic( "modelParametersWithCovariates",     c( "model" ) )
defineCovariatesData              = new_generic( "defineCovariatesData",              c( "model" ) )
hasCovariates                     = new_generic( "hasCovariates",                     c( "model" ) )
evaluateModelWithCovariates       = new_generic( "evaluateModelWithCovariates",       c( "model" ) )

# ── Package-level constants ────────────────────────────────────────────────────
# Integer codes for the covariate-equation type (avoids magic numbers throughout).
MODEL_COV_EQ_EXPONENTIAL = 1L   # theta = mu * exp(beta * cov)
MODEL_COV_EQ_ADDITIVE    = 2L   # theta = mu * (1 + beta * cov)

# Two-character separator used to build / parse reversible combination names.
# Must not appear in any covariate name or category label.
COMBINATION_SEP = "::"

# ── Occasion helpers ──────────────────────────────────────────────────────────

# Returns a safe positive integer occasion count from the model property.
# Values < 1, NA, and length-0 all normalise to 1L.
getRequestedNumberOfOccasions = function( model ) {
  n = prop( model, "numberOfOccasions" )
  if ( length( n ) == 0L || is.na( n[[1L]] ) || n[[1L]] < 1 ) return( 1L )
  as.integer( n[[1L]] )
}

# Infers the occasion count from IOV-covariate sequence lengths and validates
# that all such covariates agree on the same count.
getOccasionsFromIOVCovariates = function( modelCovariates ) {
  if ( length( modelCovariates ) == 0L ) return( 1L )

  covariateWithIov = .filterCovariatesByClass( modelCovariates, "CategoricalCovariateWithIOV" )
  if ( length( covariateWithIov ) == 0L ) return( 1L )

  sequenceLengths = covariateWithIov |>
    map( ~ map_int( prop( .x, "sequences" ), length ) ) |>
    list_c() |>
    unique()

  if ( length( sequenceLengths ) > 1L )
    stop( "All occasion-based covariate sequences must use the same number of occasions." )

  as.integer( sequenceLengths[[1L]] )
}

# Single source of truth: merges the explicit user setting and the value
# induced by occasion-based covariate sequences.
getNumberOfOccasionsForModel = function( model ) {
  requested = getRequestedNumberOfOccasions( model )
  inferred  = getOccasionsFromIOVCovariates( prop( model, "modelCovariates" ) )

  if ( requested > 1L && inferred > 1L && requested != inferred )
    stop( "`numberOfOccasions` must match the number of occasions defined by occasion-based covariate sequences." )

  max( requested, inferred, 1L )
}

# TRUE when the complex combination/occasion evaluation path is required.
usesCovariateOccasionStructure = function( model ) {
  length( prop( model, "modelCovariates" ) ) > 0L ||
    getNumberOfOccasionsForModel( model ) > 1L
}

# ── Private covariate-splitting helpers ───────────────────────────────────────

# Split a flat covariate list by short class name (strips the "PFIM::" prefix).
.splitCovariatesByClass = function( covariates ) {
  if ( length( covariates ) == 0L ) return( list() )
  split( covariates, map_chr( covariates, ~ str_remove( pluck( class( .x ), 1L ), "^PFIM::" ) ) )
}

# Return only covariates matching a given short class name.
.filterCovariatesByClass = function( covariates, shortClass ) {
  pluck( .splitCovariatesByClass( covariates ), shortClass, .default = list() )
}

# ── Combination-name encoding / decoding ──────────────────────────────────────

# Builds a fully reversible combination name from a named list of
# (covariateName -> categoryOrSequenceLabel) pairs.
# Uses COMBINATION_SEP ("::") and "=" so that underscores inside names are safe.
.encodeCombinationName = function( parts ) {
  if ( length( parts ) == 0L ) return( "Reference" )
  paste( names( parts ), unlist( parts ), sep = "=", collapse = COMBINATION_SEP )
}

# Inverse of .encodeCombinationName.
.decodeCombinationName = function( combinationName ) {
  if ( combinationName == "Reference" ) return( list() )

  tokens = str_split( combinationName, fixed( COMBINATION_SEP ) )[[1L]]
  pairs  = str_split( tokens, fixed( "=" ) )

  if ( any( map_int( pairs, length ) != 2L ) )
    stop( sprintf( "Malformed combination name: '%s'", combinationName ) )

  set_names(
    map( pairs, ~ .x[[2L]] ),
    map_chr( pairs, ~ .x[[1L]] )
  )
}

# Recovers the full covariate-value mapping for a combination name, injecting
# reference levels (first category / first sequence) for absent covariates.
parseCombinationName = function( combinationName, modelCovariates ) {
  decoded     = .decodeCombinationName( combinationName )
  covNames    = map_chr( modelCovariates, ~ prop( .x, "name" ) )
  missingMask = map_lgl( covNames, ~ is.null( decoded[[ .x ]] ) )

  refDefaults = map( modelCovariates[ missingMask ], function( cov ) {
    if ( inherits( cov, "CategoricalCovariate" ) ) {
      prop( cov, "categories" )[[1L]]
    } else {
      seqs = prop( cov, "sequences" )
      if ( is.null( names( seqs ) ) ) "sequence_1" else names( seqs )[[1L]]
    }
  }) |> set_names( covNames[ missingMask ] )

  c( decoded, refDefaults )
}

# ==============================================================================
# S7 Method Implementations
# ==============================================================================

method( evaluateModelWithCovariates, Model ) = function( model, arm, evaluateModelCore ) {

  covariatesCombinations = prop( model, "covariatesCombination" )$combinations
  modelParamsWithCov     = prop( model, "modelParametersWithCovariates" )
  baseModelParameters    = prop( model, "modelParameters" )

  map( seq_len( nrow( covariatesCombinations ) ), function( i ) {
    combinationName    = covariatesCombinations$name[[i]]
    proportion         = covariatesCombinations$proportion[[i]]
    parametersForCombo = modelParamsWithCov[[ combinationName ]]

    evaluationByOccasion = map( names( parametersForCombo ), function( occasion ) {
      occasionParams = parametersForCombo[[ occasion ]]

      updatedParameters = map2( baseModelParameters, occasionParams, function( param, newMu ) {
        if ( prop( param, "fixedMu" ) ) return( param )
        distribution = prop( param, "distribution" )
        prop( distribution, "mu" ) = newMu
        prop( param, "distribution" ) = distribution
        param
      })

      tempModel = model
      prop( tempModel, "modelParameters" ) = updatedParameters
      tempModel = defineModelAdministration( tempModel, arm )
      list( occasion = occasion, evaluation = evaluateModelCore( tempModel, arm ) )
    })

    list( combination = combinationName, proportion = proportion, evaluations = evaluationByOccasion )
  })
}

method( hasCovariates, Model ) = function( model ) {
  length( prop( model, "modelCovariates" ) ) > 0L
}

method( defineCovariatesData, Model ) = function( model ) {
  model |>
    evaluateCovariatesEffects()         |>
    generateCovariatesCombination()     |>
    modelParametersWithCovariates()     |>
    evaluateOmegaMatrixFromCovariates()
}

method( evaluateCovariatesEffects, Model ) = function( model ) {
  covariates = prop( model, "modelCovariates" )

  if ( length( covariates ) == 0L ) {
    prop( model, "covariatesEffect" ) = list()
    return( model )
  }

  zeroEffectVector = prop( model, "modelParameters" ) |>
    map_chr( ~ prop( .x, "name" ) ) |>
    set_names() |>
    map_dbl( ~ 0 )

  covariatesEffects = imap( .splitCovariatesByClass( covariates ), function( covList, className ) {
    effects = map( covList, ~ getCovariateEffects( .x, zeroEffectVector ) )
    set_names( effects, map_chr( covList, ~ prop( .x, "name" ) ) )
  })

  prop( model, "covariatesEffect" ) = covariatesEffects
  model
}

method( generateCovariatesCombination, Model ) = function( model ) {
  covariates = prop( model, "modelCovariates" )

  if ( length( covariates ) == 0L ) {
    prop( model, "covariatesCombination" ) = list(
      combinations = data.frame(
        combinationIndex               = 1L,
        covariateWithoutIovCombination = 1L,
        covariateWithIovCombination    = 1L,
        proportion                     = 1,
        name                           = "Reference",
        stringsAsFactors               = FALSE
      ),
      covariateWithoutIovGridValues = NULL,
      covariateWithIovGridValues    = NULL
    )
    return( model )
  }

  covariateWithoutIov = .filterCovariatesByClass( covariates, "CategoricalCovariate"        )
  covariateWithIov    = .filterCovariatesByClass( covariates, "CategoricalCovariateWithIOV" )

  # Build an index grid (data.frame) for one list of covariates.
  makeIndexGrid = function( covList, propName ) {
    if ( length( covList ) == 0L ) return( NULL )
    covList |>
      map( ~ seq_along( prop( .x, propName ) ) ) |>
      set_names( map_chr( covList, ~ prop( .x, "name" ) ) ) |>
      cross_df() |>
      as.data.frame()
  }

  covWithoutIovGrid = makeIndexGrid( covariateWithoutIov, "categories" )
  covWithIovGrid    = makeIndexGrid( covariateWithIov,    "sequences"  )

  hasWithoutIov = !is.null( covWithoutIovGrid )
  hasWithIov    = !is.null( covWithIovGrid    )

  # Unified cross-product: absent type contributes a single pseudo-index 1L.
  withoutIovIdx = if ( hasWithoutIov ) seq_len( nrow( covWithoutIovGrid ) ) else 1L
  withIovIdx    = if ( hasWithIov    ) seq_len( nrow( covWithIovGrid    ) ) else 1L

  fullCombinations = expand.grid(
    covariateWithoutIovIndex = withoutIovIdx,
    covariateWithIovIndex    = withIovIdx,
    stringsAsFactors         = FALSE
  )

  combinations = pmap_dfr( fullCombinations, function( covariateWithoutIovIndex,
                                                       covariateWithIovIndex ) {
    dataWithoutIov = if ( hasWithoutIov ) {
      imap( covariateWithoutIov, function( cov, i ) {
        catIdx = covWithoutIovGrid[ covariateWithoutIovIndex, i ]
        list(
          prop = prop( cov, "categoriesProportions" )[[ catIdx ]],
          key  = prop( cov, "name" ),
          val  = prop( cov, "categories"            )[[ catIdx ]]
        )
      })
    } else list()

    dataWithIov = if ( hasWithIov ) {
      imap( covariateWithIov, function( cov, i ) {
        seqIdx   = covWithIovGrid[ covariateWithIovIndex, i ]
        seqs     = prop( cov, "sequences"            )
        seqProps = prop( cov, "sequencesProportions" )
        seqNames = if ( is.null( names( seqs ) ) ) paste0( "sequence_", seq_along( seqs ) ) else names( seqs )
        list(
          prop = seqProps[[ seqIdx ]],
          key  = prop( cov, "name"    ),
          val  = seqNames[[ seqIdx ]]
        )
      })
    } else list()

    allData         = c( dataWithoutIov, dataWithIov )
    proportion      = reduce( map_dbl( allData, "prop" ), `*`, .init = 1 )
    nameParts       = set_names( map_chr( allData, "val" ), map_chr( allData, "key" ) )
    combinationName = .encodeCombinationName( as.list( nameParts ) )

    data.frame(
      combinationIndex               = NA_integer_,
      covariateWithoutIovCombination = covariateWithoutIovIndex,
      covariateWithIovCombination    = covariateWithIovIndex,
      proportion                     = proportion,
      name                           = combinationName,
      stringsAsFactors               = FALSE
    )
  })

  combinations$combinationIndex = seq_len( nrow( combinations ) )

  prop( model, "covariatesCombination" ) = list(
    combinations                  = combinations,
    covariateWithoutIovGridValues = covWithoutIovGrid,
    covariateWithIovGridValues    = covWithIovGrid
  )
  model
}

method( modelParametersWithCovariates, Model ) = function( model ) {

  modelCovariatesEquation = prop( model, "modelCovariatesEquation" )
  covariatesEffect        = prop( model, "covariatesEffect"        )
  modelParameters         = prop( model, "modelParameters"         )

  muValues = map_dbl( modelParameters, ~ prop( prop( .x, "distribution" ), "mu" ) ) |>
    set_names( map_chr( modelParameters, ~ prop( .x, "name" ) ) )

  combinations      = prop( model, "covariatesCombination" )$combinations
  covWithoutIovGrid = prop( model, "covariatesCombination" )$covariateWithoutIovGridValues
  covWithIovGrid    = prop( model, "covariatesCombination" )$covariateWithIovGridValues

  covariateWithoutIov = .filterCovariatesByClass( prop( model, "modelCovariates" ), "CategoricalCovariate"        )
  covariateWithIov    = .filterCovariatesByClass( prop( model, "modelCovariates" ), "CategoricalCovariateWithIOV" )

  covWithoutIovNames = map_chr( covariateWithoutIov, ~ prop( .x, "name" ) )
  covWithIovNames    = map_chr( covariateWithIov,    ~ prop( .x, "name" ) )

  maxOccasions = getNumberOfOccasionsForModel( model )
  zeroEffect   = set_names( rep( 0, length( muValues ) ), names( muValues ) )

  # Helper: cumulated covariate effect (without IOV) for one combination row.
  effectWithoutIov = function( rowIndex ) {
    if ( length( covWithoutIovNames ) == 0L ) return( zeroEffect )
    imap( covWithoutIovNames, function( covName, j ) {
      catIdx = as.integer( covWithoutIovGrid[ rowIndex, j ] )
      covariatesEffect$CategoricalCovariate[[ covName ]][[ catIdx ]]
    }) |> reduce( `+`, .init = zeroEffect )
  }

  # Helper: cumulated covariate effect (with IOV) for one combination row and occasion.
  effectWithIov = function( rowIndex, occIndex ) {
    if ( length( covWithIovNames ) == 0L ) return( zeroEffect )
    imap( covWithIovNames, function( covName, j ) {
      seqIdx    = as.integer( covWithIovGrid[ rowIndex, j ] )
      sequences = prop( covariateWithIov[[ j ]], "sequences" )
      seqName   = if ( is.null( names( sequences ) ) ) paste0( "sequence_", seqIdx ) else names( sequences )[[ seqIdx ]]
      effList   = covariatesEffect$CategoricalCovariateWithIOV[[ covName ]][[ seqName ]]
      if ( is.null( effList ) || occIndex > length( effList ) ) zeroEffect else effList[[ occIndex ]]
    }) |> reduce( `+`, .init = zeroEffect )
  }

  modelParamsForCov = pmap( combinations, function( name,
                                                    covariateWithoutIovCombination,
                                                    covariateWithIovCombination, ... ) {
    baseEffect = effectWithoutIov( covariateWithoutIovCombination )

    map( seq_len( maxOccasions ), function( occ ) {
      combinedEffect = baseEffect + effectWithIov( covariateWithIovCombination, occ )

      if ( is.null( modelCovariatesEquation ) ) {
        if ( !all( combinedEffect == 0 ) )
          stop( "`modelCovariatesEquation` must be defined when covariate effects are present." )
        return( muValues )
      }

      prop(
        computeCovariateValue( modelCovariatesEquation, beta = muValues, combinedEffect = combinedEffect ),
        "value"
      )
    }) |> set_names( paste0( "occ", seq_len( maxOccasions ) ) )

  }) |> set_names( combinations$name )

  prop( model, "modelParametersWithCovariates" ) = modelParamsForCov
  model
}

method( evaluateOmegaMatrixFromCovariates, Model ) = function( model ) {

  modelParameters = prop( model, "modelParameters" )
  paramNames      = map_chr( modelParameters, ~ prop( .x, "name" ) )

  # omega / gamma store standard deviations; square them for the variance matrix.
  omegaDiag = map_dbl( modelParameters, ~ prop( prop( .x, "distribution" ), "omega" ) )
  gammaDiag = map_dbl( modelParameters, ~ prop( .x, "gamma" ) )

  omegaMat = diag( omegaDiag^2 )
  gammaMat = diag( gammaDiag^2 )

  maxOccasions = getNumberOfOccasionsForModel( model )
  nIov         = if ( maxOccasions > 1L ) maxOccasions else 0L
  nIiv         = length( paramNames )
  nTotal       = nIiv + nIov * nIiv

  omegaWithIOV = matrix( 0, nrow = nTotal, ncol = nTotal )
  omegaWithIOV[ seq_len( nIiv ), seq_len( nIiv ) ] = omegaMat

  if ( nIov > 0L ) {
    iovIndices = map( seq_len( nIov ), ~ nIiv + ( .x - 1L ) * nIiv + seq_len( nIiv ) )
    for ( idx in iovIndices ) omegaWithIOV[ idx, idx ] = gammaMat
  }

  iovNames = if ( nIov > 0L ) {
    paste0( rep( paramNames, nIov ), "_occ", rep( seq_len( nIov ), each = nIiv ) )
  } else character( 0L )

  dimnames( omegaWithIOV ) = list( c( paramNames, iovNames ), c( paramNames, iovNames ) )

  prop( model, "omegaWithIOV" ) = omegaWithIOV
  model
}

method( finiteDifferenceHessian, Model ) = function( model ) {

  pars    = map_dbl( prop( model, "modelParameters" ), ~ prop( prop( .x, "distribution" ), "mu" ) )
  npar    = length( pars )
  relStep = .Machine$double.eps^( 1 / 3 )
  incr    = pmax( abs( pars ), 0 ) * relStep
  baseInd = diag( npar )

  # Build shift columns and quadratic-approximation fraction vector.
  extraCols = map( seq_len( npar - 1L ), ~ baseInd[ , .x ] + baseInd[ , -seq_len( .x ) ] )
  extraFrac = map( seq_len( npar - 1L ), ~ incr[ .x ] * incr[ -seq_len( .x ) ]            )

  cols    = c( list( 0, baseInd, -baseInd ), extraCols )
  frac    = c( 1, incr, incr^2, unlist( extraFrac ) )
  indMat  = do.call( cbind, cols )
  shifted = pars + incr * indMat

  indMatT = t( indMat )
  Xcols   = c(
    list( 1, indMatT, indMatT^2 ),
    map( seq_len( npar - 1L ), ~ indMatT[ , .x ] * indMatT[ , -seq_len( .x ) ] )
  )

  prop( model, "parametersForComputingGradient" ) = list(
    XcolsInv = solve( do.call( cbind, Xcols ) ),
    shifted  = shifted,
    frac     = frac
  )
  model
}

# ==============================================================================
# Gradient Computing Pipeline
# ==============================================================================

evaluateModelGradientCore = function( model, arm ) {

  parameters  = prop( model, "modelParameters" )
  paramNames  = map_chr( parameters, ~ prop( .x, "name" ) )
  outputNames = prop( model, "outputNames" )
  gradParams  = prop( model, "parametersForComputingGradient" )

  XcolsInv = gradParams$XcolsInv
  shifted  = gradParams$shifted
  frac     = gradParams$frac

  # Strip covariates / force single occasion so the core stays atomic.
  tempModel = model
  prop( tempModel, "modelCovariates"   ) = list()
  prop( tempModel, "numberOfOccasions" ) = 1

  evaluations = map( seq_len( ncol( shifted ) ), function( iter ) {
    shiftedParams = map2( parameters, shifted[ , iter ], function( param, newMu ) {
      distr = prop( param, "distribution" )
      prop( distr, "mu" ) = newMu
      prop( param, "distribution" ) = distr
      param
    })
    iterModel = tempModel
    prop( iterModel, "modelParameters" ) = shiftedParams
    iterModel = defineModelAdministration( iterModel, arm )
    evaluateModel( iterModel, arm )
  })

  map( outputNames, function( outName ) {
    outputMatrix = evaluations |>
      map( ~ .x[[ outName ]][ , 2L ] ) |>
      as.data.frame() |>
      t()

    raw   = XcolsInv %*% outputMatrix / frac
    grads = as.data.frame(
      matrix( t( raw )[ , 2L:( 1L + length( parameters ) ) ], ncol = length( parameters ) )
    )
    colnames( grads ) = paramNames
    grads
  }) |> set_names( outputNames )
}

method( evaluateModelGradient, Model ) = function( model, arm ) {
  if ( !usesCovariateOccasionStructure( model ) ) {
    evaluateModelGradientCore( model, arm )
  } else {
    evaluateModelGradientWithCovariates( model, arm, evaluateModelGradientCore )
  }
}

evaluateModelGradientWithCovariates = function( model, arm, evaluateModelGradientCore ) {

  covariatesCombinations = prop( model, "covariatesCombination" )$combinations
  modelParamsWithCov     = prop( model, "modelParametersWithCovariates" )
  covariatesEffect       = prop( model, "covariatesEffect" )
  modelCovariatesEq      = prop( model, "modelCovariatesEquation" )
  modelParameters        = prop( model, "modelParameters" )
  modelCovariates        = prop( model, "modelCovariates" )

  paramNames = map_chr( modelParameters, ~ prop( .x, "name" ) )
  muValues   = map_dbl( modelParameters, ~ prop( prop( .x, "distribution" ), "mu" ) ) |>
    set_names( paramNames )

  covEqType = if ( inherits( modelCovariatesEq, "Additive" ) ) MODEL_COV_EQ_ADDITIVE
  else                                              MODEL_COV_EQ_EXPONENTIAL

  covariateWithoutIov = .filterCovariatesByClass( modelCovariates, "CategoricalCovariate"        )
  covariateWithIov    = .filterCovariatesByClass( modelCovariates, "CategoricalCovariateWithIOV" )

  # Built once, reused across all combination x occasion iterations.
  betaList = .buildBetaList( covariateWithoutIov, covariateWithIov, covariatesEffect )

  # ── Predicates and correction helpers (closed over local env) ────────────────

  .isBetaActive = function( betaInfo, combinationInfo, occasionIndex ) {
    if ( !betaInfo$isIOV ) {
      identical( combinationInfo[[ betaInfo$covName ]], betaInfo$category )
    } else {
      seqLabel = combinationInfo[[ betaInfo$covName ]]
      if ( is.null( seqLabel ) ) return( FALSE )
      cov      = covariateWithIov[[ betaInfo$covIndex ]]
      seqs     = prop( cov, "sequences" )
      seqNames = if ( is.null( names( seqs ) ) ) paste0( "sequence_", seq_along( seqs ) ) else names( seqs )
      seqIdx   = match( seqLabel, seqNames )
      if ( is.na( seqIdx ) ) return( FALSE )
      seqValues = seqs[[ seqIdx ]]
      if ( occasionIndex > length( seqValues ) ) return( FALSE )
      identical( seqValues[[ occasionIndex ]], betaInfo$category )
    }
  }

  # Multiplicative correction per parameter for the additive model:
  # d(theta)/d(mu) = (1 + beta * cov).
  .additiveMuCorrection = function( combinationInfo, occIndex ) {
    correction = set_names( rep( 1, length( muValues ) ), names( muValues ) )

    applyCorrection = function( covList, effectClass ) {
      for ( cov in covList ) {
        covName = prop( cov, "name" )
        val     = combinationInfo[[ covName ]]
        if ( is.null( val ) ) next

        if ( inherits( cov, "CategoricalCovariateWithIOV" ) ) {
          seqs     = prop( cov, "sequences" )
          seqNames = if ( is.null( names( seqs ) ) ) paste0( "sequence_", seq_along( seqs ) ) else names( seqs )
          seqIdx   = match( val, seqNames )
          if ( is.na( seqIdx ) || occIndex > length( seqs[[ seqIdx ]] ) ) next
          activeCat = seqs[[ seqIdx ]][[ occIndex ]]
        } else {
          activeCat = val
        }

        if ( identical( activeCat, prop( cov, "categories" )[[1L]] ) ) next

        effects = covariatesEffect[[ effectClass ]][[ covName ]]
        catIdx  = match( activeCat, prop( cov, "categories" ) )
        if ( is.na( catIdx ) || is.null( effects[[ catIdx ]] ) ) next

        for ( pName in names( effects[[ catIdx ]] ) ) {
          b = effects[[ catIdx ]][[ pName ]]
          if ( b != 0 ) correction[[ pName ]] = correction[[ pName ]] * ( 1 + b )
        }
      }
      correction
    }

    correction = applyCorrection( covariateWithoutIov, "CategoricalCovariate"        )
    correction = applyCorrection( covariateWithIov,    "CategoricalCovariateWithIOV" )
    correction
  }

  # ── Main loop ─────────────────────────────────────────────────────────────────

  map( covariatesCombinations$name, function( combinationName ) {
    paramsForCombo  = modelParamsWithCov[[ combinationName ]]
    combinationInfo = parseCombinationName( combinationName, modelCovariates )
    proportion      = covariatesCombinations$proportion[[ match( combinationName,
                                                                 covariatesCombinations$name ) ]]

    evaluationByOccasion = imap( paramsForCombo, function( occasionParams, occasion ) {
      occIndex = as.integer( str_remove( occasion, "^occ" ) )

      updatedParams = map2( modelParameters, occasionParams, function( param, newMu ) {
        distr = prop( param, "distribution" )
        prop( distr, "mu" ) = newMu
        prop( param, "distribution" ) = distr
        param
      })

      tempModel = model
      prop( tempModel, "modelParameters"   ) = updatedParams
      prop( tempModel, "modelCovariates"   ) = list()
      prop( tempModel, "numberOfOccasions" ) = 1
      tempModel = finiteDifferenceHessian( tempModel )
      tempModel = defineModelAdministration( tempModel, arm )

      gradTheta = evaluateModelGradientCore( tempModel, arm )

      # ── Gradient w.r.t. mu ──────────────────────────────────────────────────
      gradMu = if ( covEqType == MODEL_COV_EQ_EXPONENTIAL ) {
        # d(f)/d(mu) = d(f)/d(theta) x (theta / mu)
        map( gradTheta, function( gradDf ) {
          imap_dfc( gradDf, function( gradVec, pName ) {
            mu_i    = muValues[[ pName ]]
            theta_i = occasionParams[[ pName ]]
            if ( is.null( mu_i ) || mu_i == 0 || is.na( mu_i ) ) gradVec
            else gradVec * ( theta_i / mu_i )
          })
        })
      } else {
        # d(f)/d(mu) = d(f)/d(theta) x (1 + beta x cov)
        correction = .additiveMuCorrection( combinationInfo, occIndex )
        map( gradTheta, function( gradDf ) {
          imap_dfc( gradDf, function( gradVec, pName ) gradVec * correction[[ pName ]] )
        })
      }

      # ── Gradient w.r.t. beta ────────────────────────────────────────────────
      outputNames = prop( model, "outputNames" )
      nRows       = nrow( gradMu[[ outputNames[[1L]] ]] )

      gradWithCovariates = map( outputNames, function( outName ) {
        gradMu_out    = gradMu[[ outName ]]
        gradTheta_out = gradTheta[[ outName ]]

        muCols = map( paramNames, ~ gradMu_out[[ .x ]] ) |>
          set_names( paste0( "mu_", paramNames ) )

        betaCols = imap( betaList, function( betaInfo, betaName ) {
          if ( !.isBetaActive( betaInfo, combinationInfo, occIndex ) )
            return( rep( 0, nRows ) )
          pName = betaInfo$param
          if ( covEqType == MODEL_COV_EQ_EXPONENTIAL )
            gradTheta_out[[ pName ]] * occasionParams[[ pName ]]   # d(f)/d(beta) = d(f)/d(theta) x theta
          else
            gradTheta_out[[ pName ]] * muValues[[ pName ]]         # d(f)/d(beta) = d(f)/d(theta) x mu
        })

        as.data.frame( c( muCols, betaCols ) )
      }) |> set_names( outputNames )

      list( occasion = occasion, gradient = gradWithCovariates )
    })

    list( combination = combinationName, proportion = proportion, gradients = evaluationByOccasion )
  })
}

# ── Private: beta metadata list ───────────────────────────────────────────────
# One entry per (covariate, non-reference category, affected parameter) triple.
# Built once before the combination loop; never re-created inside it.
#
# Index safety: seq_len(length(x) - 1L) + 1L produces integer(0) when
# length(x) == 1, so map() becomes a no-op with no descending-sequence risk.
.buildBetaList = function( covariateWithoutIov, covariateWithIov, covariatesEffect ) {

  betasWithoutIov = if ( length( covariateWithoutIov ) > 0L ) {
    covariateWithoutIov |>
      map( function( cov ) {
        covName    = prop( cov, "name"       )
        categories = prop( cov, "categories" )
        effects    = covariatesEffect$CategoricalCovariate[[ covName ]]

        if ( length( categories ) <= 1L ) return( list() )

        ( seq_len( length( categories ) - 1L ) + 1L ) |>
          map( function( icat ) {
            affectedParams = names( effects[[ icat ]] )[ effects[[ icat ]] != 0 ]
            map( affectedParams, ~ list(
              name  = paste0( "beta_", .x, "_", covName, "_", categories[[ icat ]] ),
              value = list( covName  = covName,
                            category = categories[[ icat ]],
                            param    = .x,
                            isIOV    = FALSE )
            ))
          }) |> list_flatten()
      }) |> list_flatten() |>
      (\( lst ) set_names( lst, map_chr( lst, function( x ) x[["name"]] ) ))() |>
      map( "value" )
  } else list()

  betasWithIov = if ( length( covariateWithIov ) > 0L ) {
    imap( covariateWithIov, function( cov, covIndex ) {
      covName    = prop( cov, "name"       )
      sequences  = prop( cov, "sequences"  )
      categories = prop( cov, "categories" )
      effects    = covariatesEffect$CategoricalCovariateWithIOV[[ covName ]]

      if ( length( categories ) <= 1L ) return( list() )

      ( seq_len( length( categories ) - 1L ) + 1L ) |>
        map( function( icat ) {
          category = categories[[ icat ]]

          affectedParams = seq_along( sequences ) |>
            map( function( iseq ) {
              effSeq = effects[[ iseq ]]
              if ( is.null( effSeq ) ) return( character( 0L ) )
              seqVals     = sequences[[ iseq ]]
              matchingOcc = detect_index( seq_along( seqVals ), ~ seqVals[[ .x ]] == category )
              if ( matchingOcc == 0L ) return( character( 0L ) )
              occEff = effSeq[[ matchingOcc ]]
              names( occEff )[ occEff != 0 ]
            }) |> list_c() |> unique()

          map( affectedParams, ~ list(
            covName    = covName,
            covIndex   = covIndex,
            category   = category,
            categories = categories,
            param      = .x,
            isIOV      = TRUE
          )) |> set_names(
            map_chr( affectedParams, ~ paste0( "beta_", .x, "_", covName, "_", category ) )
          )
        }) |> list_flatten()
    }) |> list_flatten()
  } else list()

  c( betasWithoutIov, betasWithIov )
}

aggregateGradientsWithCovariates = function( model, arm ) {

  allGradientsData = prop( arm, "evaluationGradients" )
  outputNames      = prop( model, "outputNames" )

  map( outputNames, function( outName ) {
    zeroGrad = allGradientsData[[1L]]$gradients[[1L]]$gradient[[ outName ]] * 0

    reduce( allGradientsData, function( acc, combinationData ) {
      nOcc     = length( combinationData$gradients )
      meanGrad = reduce(
        map( combinationData$gradients, ~ .x$gradient[[ outName ]] ),
        `+`
      ) / nOcc
      acc + combinationData$proportion * meanGrad
    }, .init = zeroGrad )
  }) |>
    reduce( rbind ) |>
    as.matrix()
}

# ==============================================================================
# Population Model Variance
# ==============================================================================

evaluateModelVarianceCore = function( model, evaluationModel ) {

  modelErrors = prop( model, "modelError"  )
  outputNames = prop( model, "outputNames" )

  # Error model derivatives, one entry per output (NULL entries dropped).
  errorDerivativesList = map( modelErrors, function( err ) {
    outcome = prop( err, "output" )
    if ( outcome %in% outputNames )
      evaluateErrorModelDerivatives( err, evaluationModel[[ outcome ]][ , outcome ] )
    else NULL
  }) |> compact() |> set_names( outputNames )

  errorVariance  = map( errorDerivativesList, ~ bdiag( .x$errorVariance ) ) |> bdiag()
  totalSamplings = map_int( evaluationModel, ~ length( .x[["time"]] ) ) |> sum()

  # Pre-compute start offset for each output block (no mutable iterator).
  samplingOffsets = accumulate(
    outputNames,
    function( offset, outName ) offset + length( evaluationModel[[ outName ]][["time"]] ),
    .init = 0L
  ) |> head( -1L ) |> set_names( outputNames )

  sigmaDerivatives = outputNames |>
    map( function( outName ) {
      n   = length( evaluationModel[[ outName ]][["time"]] )
      rng = ( samplingOffsets[[ outName ]] + 1L ):( samplingOffsets[[ outName ]] + n )

      map( errorDerivativesList[[ outName ]]$sigmaDerivatives, function( derivComp ) {
        mat             = matrix( 0, nrow = totalSamplings, ncol = totalSamplings )
        mat[ rng, rng ] = derivComp
        mat
      })
    }) |> list_flatten()

  list( errorVariance = errorVariance, sigmaDerivatives = sigmaDerivatives )
}

method( evaluateModelVariance, Model ) = function( model, arm ) {
  evaluationModel = prop( arm, "evaluationModel" )
  if ( !usesCovariateOccasionStructure( model ) ) {
    evaluateModelVarianceCore( model, evaluationModel )
  } else {
    evaluateModelVarianceWithCovariates( model, evaluationModel )
  }
}

evaluateModelVarianceWithCovariates = function( model, evaluationModelWithCovariates ) {
  map( evaluationModelWithCovariates, function( combinationData ) {
    list(
      combination = combinationData$combination,
      proportion  = combinationData$proportion,
      variances   = map( combinationData$evaluations, function( occasionData ) {
        list(
          occasion = occasionData$occasion,
          variance = evaluateModelVarianceCore( model, occasionData$evaluation )
        )
      })
    )
  })
}

# ==============================================================================
# Library-of-models variable substitution
# ==============================================================================

# Replace all word-boundary occurrences of `old` with `new` in `text`,
# skipping protected pharmacokinetic tokens.
#
# ICU lookbehind restriction: variable-width alternatives are unsupported, so
# protection is an early-return guard.  For "RespPK", a zero-width lookahead
# (?!\w) prevents matching inside longer identifiers.
replaceVariablesLibraryOfModels = function( text, old, new ) {
  PROTECTED_TERMS = c( "dose_", "Tinf_", "Emax" )

  if ( any( str_detect( old, fixed( PROTECTED_TERMS ) ) ) ) return( text )

  if ( old == "RespPK" ) {
    str_replace_all( text, regex( paste0( old, "(?!\\w)" ) ), new )
  } else {
    str_replace_all( text, regex( paste0( "\\b", old, "\\b" ) ), new )
  }
}
