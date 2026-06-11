#' @description
#' The class \code{PopulationFim} represents and stores information for the PopulationFim.
#' @title PopulationFim
#' @inheritParams Fim
#' @include Fim.R
#' @name PopulationFim
#' @export

PopulationFim = new_class( "PopulationFim", package = "PFIM", parent = Fim )

# ── Package-level constants ────────────────────────────────────────────────────
# Three Greek-letter dictionaries for the three rendering contexts.
# Defined once here to avoid repeated redefinition inside methods.

# Console / text display (Unicode with subscript separator)
.GREEK_CONSOLE = c(
  mu    = "\u03bc_",
  beta  = "\u03b2_",
  omega = "\u03c9\u00B2_",
  gamma = "\u03b3\u00B2_",
  sigma = "\u03c3"
)

# ggplot axis labels (bare Unicode, no separator)
.GREEK_PLOT = c(
  mu    = "\u03bc",
  omega = "\u03c9\u00B2",
  gamma = "\u03b3\u00B2",
  sigma = "\u03c3"
)

# LaTeX for kableExtra HTML tables
.GREEK_LATEX = c(
  mu    = "$\\mu_{",
  beta  = "$\\beta_{",
  omega = "$\\omega^2_{",
  gamma = "$\\gamma^2_{",
  sigma = "${\\sigma_"
)

# ==============================================================================
#' evaluateFim: evaluation of the Fim
#' @name evaluateFim
#' @param fim An object \code{PopulationFim} giving the Fim.
#' @param model An object \code{Model} giving the model.
#' @param arm An object \code{Arm} giving the arm.
#' @return The object \code{PopulationFim} with the fisherMatrix.
#' @export
# ==============================================================================

method( evaluateFim, list( PopulationFim, Model, Arm ) ) = function( fim, model, arm ) {

  parameters          = prop( model, "modelParameters" )
  armSize             = prop( arm, "size" )
  hasComplexStructure = usesCovariateOccasionStructure( model )

  # Logical masks — computed once, reused in evaluateVarianceFIM.
  # FIX: use prop() throughout instead of @-accessor (S7 API).
  isFixedMu    = map_lgl( parameters, ~ isTRUE( prop( .x, "fixedMu" ) ) ||
                            prop( prop( .x, "distribution" ), "mu" ) == 0 )
  isFixedOmega = map_lgl( parameters, ~ isTRUE( prop( .x, "fixedOmega" ) ) ||
                            prop( prop( .x, "distribution" ), "omega" ) == 0 )

  result = evaluateVarianceFIM( fim, model, arm,
                                isFixedMu    = isFixedMu,
                                isFixedOmega = isFixedOmega )

  if ( !hasComplexStructure ) {
    # MFbeta and MFVar already computed (and filtered by isFixedMu) in
    # evaluateVarianceFIM.  Remove fixed-omega rows/cols from MFVar.
    # MFVar structure: [omegas | sigmas]
    MFbeta      = result$MFbeta
    MFVar       = result$MFVar
    nSigma      = ncol( MFVar ) - length( isFixedOmega )
    keepLambda  = c( !isFixedOmega, rep( TRUE, nSigma ) )
    MFVar       = MFVar[ keepLambda, keepLambda, drop = FALSE ]
  } else {
    MFbeta = result$MFbeta
    MFVar  = result$MFVar
  }

  prop( fim, "fisherMatrix" ) = as.matrix( bdiag( MFbeta * armSize, MFVar * armSize ) )
  fim
}

# ==============================================================================
#' computeVMat: compute (1/2) * Tr(V^-1 A V^-1 B)
#' @name computeVMat
#' @param varParam1 Matrix A (derivative dV/dlambda_i)
#' @param varParam2 Matrix B (derivative dV/dlambda_j)
#' @param invCholV  \eqn{V^{-1}} computed via Cholesky
#' @return Scalar (1/2) * Tr(V^-1 A V^-1 B)
#' @export
# ==============================================================================

computeVMat = function( varParam1, varParam2, invCholV ) {
  0.5 * sum( diag( invCholV %*% varParam1 %*% invCholV %*% varParam2 ) )
}

# ==============================================================================
#' evaluateVarianceFIM: optimised computation of MFbeta and MFVar
#' @name evaluateVarianceFIM
#' @param fim   An object of class \code{PopulationFim}.
#' @param model An object of class \code{Model}.
#' @param arm          An object of class \code{Arm}.
#' @param isFixedMu    Optional logical vector; \code{TRUE} for parameters with fixed or zero mu.
#'   Computed internally if \code{NULL} (default).
#' @param isFixedOmega Optional logical vector; \code{TRUE} for parameters with fixed or zero omega.
#'   Computed internally if \code{NULL} (default).
#' @return List with MFbeta and MFVar.
# ==============================================================================

evaluateVarianceFIM = function( fim, model, arm,
                                isFixedMu    = NULL,
                                isFixedOmega = NULL ) {

  parameters     = prop( model, "modelParameters" )
  parameterNames = map_chr( parameters, ~ prop( .x, "name" ) )

  if ( is.null( isFixedMu ) ) {
    isFixedMu = map_lgl( parameters, ~ isTRUE( prop( .x, "fixedMu" ) ) ||
                           prop( prop( .x, "distribution" ), "mu" ) == 0 )
  }
  if ( is.null( isFixedOmega ) ) {
    isFixedOmega = map_lgl( parameters, ~ isTRUE( prop( .x, "fixedOmega" ) ) ||
                              prop( prop( .x, "distribution" ), "omega" ) == 0 )
  }

  # ── Helper: fill symmetric matrix upper-triangle then mirror ─────────────────
  .fillSymmetricFIM = function( dV_dlambda, temp_matrices ) {
    n   = length( dV_dlambda )
    mat = matrix( 0, nrow = n, ncol = n )
    for ( i in seq_len( n ) ) {
      for ( j in i:n ) {
        mat[ i, j ] = 0.5 * sum( temp_matrices[[ i ]] * dV_dlambda[[ j ]] )
        if ( i != j ) mat[ j, i ] = mat[ i, j ]
      }
    }
    mat
  }

  # ── Helper: outer product of a column vector (n×1 slice from a pre-transposed
  # matrix). tcrossprod(v) = v %*% t(v); BLAS DSYR path; result is n×n.
  # Pre-transposing matrices once (column extraction is cache-friendly) then
  # calling tcrossprod avoids per-iteration row-extraction + transposition.
  .outerCol = function( v ) tcrossprod( v )

  # ==========================================================================
  # one combination, one occasion
  # ==========================================================================

  if ( !usesCovariateOccasionStructure( model ) ) {

    evaluationModel  = prop( arm, "evaluationModel" )
    allGradientsData = prop( arm, "evaluationGradients" )
    varianceResults  = prop( arm, "evaluationVariance" )
    outputNames      = prop( model, "outputNames" )

    distributions = map( parameters, ~ prop( .x, "distribution" ) )
    omega_IIV     = map_dbl( distributions, ~ prop( .x, "omega" )^2 )
    muValues      = map_dbl( distributions, ~ prop( .x, "mu" ) )
    names( muValues ) = parameterNames

    nOmega  = length( parameterNames )
    nSigma  = length( varianceResults$sigmaDerivatives )
    nLambda = nOmega + nSigma

    errorVariance     = as.matrix( varianceResults$errorVariance )
    gradients         = do.call( cbind, map( outputNames, ~ t( allGradientsData[[ .x ]] ) ) )
    OMEGA             = diag(omega_IIV, nrow = length(omega_IIV))
    gradientsAdjusted = gradients * muValues

    tmp   = OMEGA %*% gradientsAdjusted
    V     = t( tmp ) %*% gradientsAdjusted + errorVariance
    V_inv = chol2inv( chol( V ) )

    MFbeta_full = ( gradients %*% V_inv ) %*% t( gradients )

    # Pre-transpose once: column extraction from gradT is cache-friendly.
    gradT      = t( gradientsAdjusted )   # nTime × nOmega
    dV_omega   = map( seq_len( nOmega ), ~ .outerCol( gradT[ , .x, drop = FALSE ] ) )
    dV_dlambda = c( dV_omega, varianceResults$sigmaDerivatives )

    temp_matrices = map( dV_dlambda, ~ V_inv %*% .x %*% V_inv )
    MFVar         = .fillSymmetricFIM( dV_dlambda, temp_matrices )

    mu_idx_keep = seq_len( nOmega )[ !isFixedMu ]
    MFbeta      = MFbeta_full[ mu_idx_keep, mu_idx_keep, drop = FALSE ]

    return( list( MFbeta = MFbeta, MFVar = MFVar ) )
  }

  # ==========================================================================
  # covariates without IOV, pure random IOV, fixed covariates + random
  # IOV, and occasion-based covariates.
  # ==========================================================================

  allGradientsData = prop( arm, "evaluationGradients" )
  varianceResults  = prop( arm, "evaluationVariance" )
  outputNames      = prop( model, "outputNames" )

  distributions = map( parameters, ~ prop( .x, "distribution" ) )
  omega_IIV     = map_dbl( distributions, ~ prop( .x, "omega" )^2 )
  gamma_values  = map_dbl( parameters, ~ pluck( .x, "gamma", .default = 0 ) )
  has_IOV       = any( gamma_values > 0 )
  muValues      = map_dbl( distributions, ~ prop( .x, "mu" ) )
  names( muValues ) = parameterNames

  nOmega = length( parameterNames )
  nGamma = sum( gamma_values > 0 )
  nSigma = length( varianceResults[[ 1L ]]$variances[[ 1L ]]$variance$sigmaDerivatives )

  numberOfOccasions = varianceResults[[ 1L ]]$variances |>
    map_chr( "occasion" ) |> unique() |> length()

  # nGamma_eff: number of estimable gamma entries in MFVar
  # - multi-occasion IOV : one gamma block per occasion → nGamma distinct entries
  # - single-occasion IOV: gamma merged with omega in OMEGA, but dV/dgamma^2
  #                        is still estimable → nGamma
  # - no IOV             : 0
  nGamma_eff = if ( has_IOV ) nGamma else 0L
  nLambda    = nOmega + nGamma_eff + nSigma

  compute_fisher_one_combination = function( iter ) {

    variance_by_occasion = map( seq_len( numberOfOccasions ),
                                ~ varianceResults[[ iter ]]$variances[[ .x ]]$variance$errorVariance )

    errorVariance = as.matrix( bdiag( variance_by_occasion ) )

    gradients_by_occasion = map( seq_len( numberOfOccasions ), function( occ ) {
      gradients_per_output = map( outputNames, function( output ) {
        t( allGradientsData[[ iter ]]$gradients[[ occ ]]$gradient[[ output ]] )
      })
      if ( length( gradients_per_output ) == 1L ) gradients_per_output[[ 1L ]]
      else do.call( cbind, gradients_per_output )
    })

    gradients  = if ( length( gradients_by_occasion ) == 1L ) gradients_by_occasion[[ 1L ]]
    else reduce( gradients_by_occasion, cbind )
    proportion = allGradientsData[[ iter ]]$proportion
    nTotalRows = nrow( gradients )

    gradients_mu = gradients[ seq_len( nOmega ), , drop = FALSE ]

    if ( nTotalRows > nOmega ) {
      gradients_beta = gradients[ ( nOmega + 1L ):nTotalRows, , drop = FALSE ]
      nBeta          = nrow( gradients_beta )
    } else {
      gradients_beta = matrix( 0, nrow = 0L, ncol = ncol( gradients ) )
      nBeta          = 0L
    }

    OMEGA = if ( has_IOV && numberOfOccasions > 1L ) {
      diag( c( omega_IIV, rep( gamma_values^2, numberOfOccasions ) ) )
    } else if ( has_IOV ) {
      diag( omega_IIV + gamma_values^2 )
    } else {
      diag( omega_IIV )
    }

    gradientsAdjusted_mu = gradients_mu * muValues
    gradientsAdjusted    = if ( nBeta > 0L ) rbind( gradientsAdjusted_mu, gradients_beta )
    else gradientsAdjusted_mu

    if ( has_IOV && numberOfOccasions > 1L ) {
      # FIX: use the actual column count of occasion 1 instead of assuming
      # symmetric sampling across occasions (asymmetric designs are valid).
      n_occ1 = ncol( gradients_by_occasion[[ 1L ]] )
      A_mu   = gradientsAdjusted_mu[ , seq_len( n_occ1 ), drop = FALSE ]
      B_mu   = gradientsAdjusted_mu[ , ( n_occ1 + 1L ):ncol( gradientsAdjusted_mu ), drop = FALSE ]

      # V_iiv: single BLAS DSYRK call via crossprod on scaled rows.
      # scaled[i,] = sqrt(omega_IIV[i]) * g_adj_i  →
      #   t(scaled) %*% scaled = sum_i omega_IIV[i] * g_adj_i^T g_adj_i
      scaled_mu = gradientsAdjusted_mu * sqrt( omega_IIV )
      V_iiv = as.matrix( crossprod( scaled_mu ) )

      # Pre-transpose A and B once: column access is cache-friendly.
      A_mu_T = t( A_mu )   # n_occ1 × nOmega
      B_mu_T = t( B_mu )   # n_occ2 × nOmega

      gamma_indices = which( gamma_values > 0 )
      V_iov = if ( length( gamma_indices ) > 0L ) {
        reduce(
          map( gamma_indices, function( i ) {
            gamma_values[ i ]^2 * as.matrix( bdiag(
              .outerCol( A_mu_T[ , i, drop = FALSE ] ),
              .outerCol( B_mu_T[ , i, drop = FALSE ] )
            ))
          }),
          `+`
        )
      } else {
        matrix( 0, nrow = nrow( errorVariance ), ncol = ncol( errorVariance ) )
      }

      V     = V_iiv + V_iov + as.matrix( errorVariance )
      V_inv = chol2inv( chol( V ) )

      MFbeta_full = ( gradients %*% V_inv ) %*% t( gradients )

      # Pre-transpose once for dV_omega and dV_gamma column extractions.
      gradMu_T = t( gradientsAdjusted_mu )   # nTime × nOmega

      # dV/domega^2_i = v_i %*% t(v_i)   (v_i = i-th column of gradMu_T)
      dV_omega = map( seq_len( nOmega ), ~ .outerCol( gradMu_T[ , .x, drop = FALSE ] ) )

      # dV/dgamma^2_j = bdiag( a_j %*% t(a_j), b_j %*% t(b_j) )
      dV_gamma = map( gamma_indices, function( i ) {
        as.matrix( bdiag(
          .outerCol( A_mu_T[ , i, drop = FALSE ] ),
          .outerCol( B_mu_T[ , i, drop = FALSE ] )
        ))
      })

    } else {
      # Single occasion or no IOV
      nOmegaRows = nrow( OMEGA )
      tmp        = OMEGA %*% gradientsAdjusted[ seq_len( nOmegaRows ), , drop = FALSE ]
      V          = t( tmp ) %*% gradientsAdjusted[ seq_len( nOmegaRows ), , drop = FALSE ] + errorVariance
      V_inv      = chol2inv( chol( V ) )

      MFbeta_full = ( gradients %*% V_inv ) %*% t( gradients )

      # Pre-transpose once for dV_omega and dV_gamma column extractions.
      gradMu_T = t( gradientsAdjusted_mu )   # nTime × nOmega

      # dV/domega^2_i = v_i %*% t(v_i)
      dV_omega = map( seq_len( nOmega ), ~ .outerCol( gradMu_T[ , .x, drop = FALSE ] ) )

      # dV/dgamma^2_i (single occasion) = dV/domega^2_i
      dV_gamma = if ( has_IOV ) {
        map( which( gamma_values > 0 ), ~ .outerCol( gradMu_T[ , .x, drop = FALSE ] ) )
      } else list()
    }

    dV_sigma = map( seq_len( nSigma ), function( i ) {
      sigma_deriv_by_occasion = map( seq_len( numberOfOccasions ),
                                     ~ varianceResults[[ iter ]]$variances[[ .x ]]$variance$sigmaDerivatives[[ i ]] )
      if ( length( sigma_deriv_by_occasion ) == 1L ) as.matrix( sigma_deriv_by_occasion[[ 1L ]] )
      else as.matrix( bdiag( sigma_deriv_by_occasion ) )
    })

    dV_dlambda    = c( dV_omega, dV_gamma, dV_sigma )
    temp_matrices = map( dV_dlambda, ~ V_inv %*% .x %*% V_inv )
    MFVar         = .fillSymmetricFIM( dV_dlambda, temp_matrices )

    as.matrix( bdiag( MFbeta_full, MFVar ) ) * proportion
  }

  fisherMatrix = reduce(
    map( seq_along( allGradientsData ), compute_fisher_one_combination ),
    `+`
  )

  # Dimension partitioning: [nOmega mu | nBeta beta | nLambda variance]
  nBeta_total = nrow( fisherMatrix ) - nLambda - nOmega
  nMuAndBeta  = nOmega + nBeta_total

  mu_idx_all  = seq_len( nOmega )
  beta_idx    = if ( nBeta_total > 0L ) ( nOmega + 1L ):nMuAndBeta else integer( 0L )
  var_idx     = ( nMuAndBeta + 1L ):nrow( fisherMatrix )

  mu_idx_keep = mu_idx_all[ !isFixedMu ]
  idx_beta_mu = c( mu_idx_keep, beta_idx )

  MFbeta = fisherMatrix[ idx_beta_mu, idx_beta_mu, drop = FALSE ]

  # Filter fixed omegas from MFVar. Structure: [nOmega omegas | nGamma_eff gammas | nSigma sigmas]
  MFVar_full = fisherMatrix[ var_idx, var_idx, drop = FALSE ]
  nRest      = nGamma_eff + nSigma          # gammas + sigmas (never filtered here)
  keepLambda = c( !isFixedOmega, rep( TRUE, nRest ) )
  MFVar      = MFVar_full[ keepLambda, keepLambda, drop = FALSE ]

  list( MFbeta = MFbeta, MFVar = MFVar )
}

# ==============================================================================
#' setOptimalArms: set the optimal arms of a MultiplicativeAlgorithm.
#' @name setOptimalArms
#' @param fim An object \code{PopulationFim} giving the Fim.
#' @param optimizationAlgorithm An object \code{MultiplicativeAlgorithm}.
#' @return The list optimalArms.
#' @export
# ==============================================================================

method( setOptimalArms, list( PopulationFim, MultiplicativeAlgorithm ) ) = function( fim, optimizationAlgorithm ) {

  outputs      = prop( optimizationAlgorithm, "multiplicativeAlgorithmOutputs" )
  armFims      = outputs$armFims
  weights      = outputs$optimalWeights
  weightsIndex = outputs$weightsIndex

  N_total                    = prop( armFims[[ weightsIndex[ 1L ] ]][[ 1L ]], "size" )
  numberOfIndividualPerGroup = round( N_total * weights / sum( weights ) )

  armList = map2( weightsIndex, numberOfIndividualPerGroup, function( idx, size ) {
    arm = armFims[[ idx ]][[ 1L ]]
    prop( arm, "size" ) = size
    prop( arm, "name" ) = paste0( "Arm", idx )
    arm
  })

  armList[ order( map_dbl( armList, ~ prop( .x, "size" ) ), decreasing = TRUE ) ]
}

# ==============================================================================
#' setOptimalArms: set the optimal arms of a FedorovWynnAlgorithm.
#' @name setOptimalArms
#' @param fim An object \code{PopulationFim} giving the Fim.
#' @param optimizationAlgorithm An object \code{FedorovWynnAlgorithm}.
#' @return The list optimalArms.
#' @export
# ==============================================================================

method( setOptimalArms, list( PopulationFim, FedorovWynnAlgorithm ) ) = function( fim, optimizationAlgorithm ) {

  outputs             = prop( optimizationAlgorithm, "FedorovWynnAlgorithmOutputs" )
  numberOfIndividuals = outputs$numberOfIndividuals
  listArms            = outputs$listArms

  map2( listArms, seq_along( listArms ), function( listArm, iter ) {
    prop( listArm$arm, "size" ) = as.double( numberOfIndividuals[[ iter ]] )
    prop( listArm$arm, "name" ) = paste0( "Arm", iter )
    listArm
  })
}

# ==============================================================================
# Private helpers shared by setEvaluationFim, plotSEFIM, plotRSEFIM, tablesForReport
# ==============================================================================

# Extract estimable parameter names filtered by fixedMu / fixedOmega / zero.
# Uses S7 prop() throughout — no @ accessor.
.muNames = function( parameters, greek ) {
  parameters |>
    keep( ~ !isTRUE( prop( .x, "fixedMu" ) ) &&
            prop( prop( .x, "distribution" ), "mu" ) != 0 ) |>
    map_chr( ~ paste0( greek, prop( .x, "name" ) ) )
}

.omegaNames = function( parameters, greek ) {
  parameters |>
    keep( ~ !isTRUE( prop( .x, "fixedOmega" ) ) &&
            prop( prop( .x, "distribution" ), "omega" ) != 0 ) |>
    map_chr( ~ paste0( greek, prop( .x, "name" ) ) )
}

.gammaNames = function( parameters, greek ) {
  parameters |>
    keep( ~ pluck( .x, "gamma", .default = 0 ) > 0 &&
            !isTRUE( prop( .x, "fixedOmega" ) ) ) |>
    map_chr( ~ paste0( greek, prop( .x, "name" ) ) )
}

.sigmaNames = function( modelError, greekSigma, sep = "_" ) {
  modelError |>
    map( function( error ) {
      output = prop( error, "output" )
      c(
        if ( prop( error, "sigmaInter" ) != 0 && !prop( error, "sigmaInterFixed" ) )
          paste0( greekSigma, sep, "inter_", output ),
        if ( prop( error, "sigmaSlope"  ) != 0 && !prop( error, "sigmaSlopeFixed"  ) )
          paste0( greekSigma, sep, "slope_", output )
      )
    }) |> unlist( use.names = FALSE )
}

# Extract estimable beta column names in the same order as evaluateModelGradientWithCovariates.
# Uses .splitCovariatesByClass from Model.R — no inline split().
.betaColumnNames = function( modelCovariates ) {
  if ( length( modelCovariates ) == 0L ) return( character( 0L ) )

  extractNames = function( covList ) {
    if ( length( covList ) == 0L ) return( character( 0L ) )
    covList |>
      map( function( cov ) {
        covName    = prop( cov, "name" )
        categories = prop( cov, "categories" )
        effects    = prop( cov, "effects" )
        # seq_len(n-1)+1 is safe: produces integer(0) if length(categories)==1
        ( seq_len( length( categories ) - 1L ) + 1L ) |>
          map( function( icat ) {
            cat = categories[[ icat ]]
            if ( !cat %in% names( effects ) ) return( character( 0L ) )
            affectedParams = names( effects[[ cat ]] )[ effects[[ cat ]] != 0 ]
            paste0( "beta_", affectedParams, "_", covName, "_", cat )
          }) |> unlist( use.names = FALSE )
      }) |> unlist( use.names = FALSE )
  }

  byClass     = .splitCovariatesByClass( modelCovariates )
  withoutIov  = pluck( byClass, "CategoricalCovariate",        .default = list() )
  withIov     = pluck( byClass, "CategoricalCovariateWithIOV", .default = list() )
  unique( c( extractNames( withoutIov ), extractNames( withIov ) ) )
}

# Build a named numeric vector (betaName -> numeric value) from modelCovariates.
# Early-return guard ensures type stability even when called on an empty list.
# as.numeric() inside the map guarantees a homogeneous atomic numeric output
# regardless of the storage type of the original effect coefficient.
.betaDict = function( modelCovariates ) {
  if ( length( modelCovariates ) == 0L )
    return( set_names( numeric( 0L ), character( 0L ) ) )

  lst = modelCovariates |>
    map( function( cov ) {
      covName    = prop( cov, "name" )
      categories = prop( cov, "categories" )
      effects    = prop( cov, "effects" )
      ( seq_len( length( categories ) - 1L ) + 1L ) |>
        map( function( icat ) {
          cat = categories[[ icat ]]
          if ( !cat %in% names( effects ) ) return( list() )
          eff = effects[[ cat ]]
          names( eff )[ eff != 0 ] |>
            map( function( param ) {
              list(
                key = paste0( "beta_", param, "_", covName, "_", cat ),
                val = as.numeric( eff[[ param ]] )
              )
            })
        }) |> list_flatten()
    }) |> list_flatten()

  if ( length( lst ) == 0L )
    return( set_names( numeric( 0L ), character( 0L ) ) )

  set_names( map_dbl( lst, "val" ), map_chr( lst, "key" ) )
}

# Shared barplot backbone for plotSEFIM and plotRSEFIM.
# metric = "SE" or "RSE"
.plotFimBars = function( fim, evaluation, metric ) {

  parameters = prop( evaluation, "modelParameters" )
  modelError = prop( evaluation, "modelError" )

  fim            = setEvaluationFim( prop( evaluation, "fim" ), evaluation )
  standardErrors = prop( fim, "SEAndRSE" )

  gamma_values = map_dbl( parameters, ~ pluck( .x, "gamma", .default = 0 ) )
  has_IOV      = any( gamma_values > 0 )
  greek        = .GREEK_PLOT

  paramsMu    = parameters |>
    keep( ~ !isTRUE( prop( .x, "fixedMu" ) ) &&
            prop( prop( .x, "distribution" ), "mu" ) != 0 ) |>
    map_chr( ~ prop( .x, "name" ) )

  paramsOmega = parameters |>
    keep( ~ !isTRUE( prop( .x, "fixedOmega" ) ) &&
            prop( prop( .x, "distribution" ), "omega" ) != 0 ) |>
    map_chr( ~ prop( .x, "name" ) )

  paramsGamma = if ( has_IOV ) {
    parameters |>
      keep( ~ pluck( .x, "gamma", .default = 0 ) > 0 &&
              !isTRUE( prop( .x, "fixedOmega" ) ) ) |>
      map_chr( ~ prop( .x, "name" ) )
  } else character( 0L )

  paramsSigma = .sigmaNames( modelError, greek[ "sigma" ] )

  yValues = if ( metric == "SE" ) standardErrors$SE$SE else standardErrors$RSE$RSE

  data = data.frame(
    Parameter        = c( paramsMu, paramsOmega, paramsGamma, paramsSigma ),
    parametersValues = standardErrors$SEAndRSE$parametersValues,
    y                = yValues,
    cat              = paste0( metric, " ", c(
      rep( greek[ "mu"    ], length( paramsMu    ) ),
      rep( greek[ "omega" ], length( paramsOmega ) ),
      rep( greek[ "gamma" ], length( paramsGamma ) ),
      rep( greek[ "sigma" ], length( paramsSigma ) )
    ))
  )
  names( data )[ names( data ) == "y" ] = metric

  facet_levels = paste0( metric, " ", c(
    greek[ "mu" ],
    greek[ "omega" ],
    if ( has_IOV ) greek[ "gamma" ],
    greek[ "sigma" ]
  ))

  ggplot( data, aes( x = Parameter, y = .data[[ metric ]] ) ) +
    geom_bar( stat = "identity", show.legend = FALSE ) +
    facet_wrap( ~ factor( cat, levels = facet_levels ), scales = "free_x" ) +
    theme(
      legend.position   = "none",
      plot.title        = element_text( size = 16, hjust = 0.5 ),
      axis.title        = element_text( size = 16 ),
      axis.text.x       = element_text( size = 16, angle = 90, vjust = 0.5 ),
      axis.text.y       = element_text( size = 16 ),
      strip.text.x      = element_text( size = 16 )
    )
}

# Helper: full path to an RMarkdown template in the PFIM package.
.reportTemplatePath = function( filename ) {
  file.path( system.file( package = "PFIM" ),
             "rmarkdown", "templates", "skeleton", filename )
}

# ==============================================================================
#' setEvaluationFim: set the Fim results.
#' @name setEvaluationFim
#' @param fim        An object \code{PopulationFim} giving the Fim.
#' @param evaluation An object \code{Evaluation} giving the evaluation of the model.
#' @return The object \code{PopulationFim} with fisherMatrix, fixedEffects,
#'   condNumberFixedEffects, and SEAndRSE populated.
#' @export
# ==============================================================================

method( setEvaluationFim, PopulationFim ) = function( fim, evaluation ) {

  parameters      = prop( evaluation, "modelParameters" )
  modelError      = prop( evaluation, "modelError" )
  modelCovariates = prop( evaluation, "modelCovariates" )
  greek           = .GREEK_CONSOLE

  has_IOV = any( map_dbl( parameters, ~ pluck( .x, "gamma", .default = 0 ) ) > 0 )
  has_cov = length( modelCovariates ) > 0

  # ── Column names ─────────────────────────────────────────────────────────────
  columnNamesMu    = .muNames(    parameters, greek[ "mu"    ] )
  columnNamesBeta  = if ( has_cov ) .betaColumnNames( modelCovariates ) else character( 0L )
  columnNamesOmega = .omegaNames( parameters, greek[ "omega" ] )
  columnNamesGamma = if ( has_IOV ) .gammaNames( parameters, greek[ "gamma" ] ) else character( 0L )
  columnNamesSigma = .sigmaNames( modelError, greek[ "sigma" ] )

  # ── Parameter values ─────────────────────────────────────────────────────────
  muValues = parameters |>
    keep( ~ !isTRUE( prop( .x, "fixedMu" ) ) &&
            prop( prop( .x, "distribution" ), "mu" ) != 0 ) |>
    map_dbl( ~ prop( prop( .x, "distribution" ), "mu" ) )

  betaValues = if ( has_cov && length( columnNamesBeta ) > 0L ) {
    dict = .betaDict( modelCovariates )
    map_dbl( columnNamesBeta, function( nm ) {
      if ( nm %in% names( dict ) ) dict[[ nm ]] else NA_real_
    })
  } else numeric( 0L )

  omegaValues = parameters |>
    keep( ~ !isTRUE( prop( .x, "fixedOmega" ) ) &&
            prop( prop( .x, "distribution" ), "omega" ) != 0 ) |>
    map_dbl( ~ prop( prop( .x, "distribution" ), "omega" )^2 )

  gammaValues = if ( has_IOV ) {
    parameters |>
      keep( ~ pluck( .x, "gamma", .default = 0 ) > 0 &&
              !isTRUE( prop( .x, "fixedOmega" ) ) ) |>
      map_dbl( ~ pluck( .x, "gamma" )^2 )
  } else numeric( 0L )

  sigmaValues = modelError |>
    map( function( error ) {
      vals = list()
      if ( prop( error, "sigmaInter" ) != 0 && !prop( error, "sigmaInterFixed" ) )
        vals$sigmaInter = prop( error, "sigmaInter" )
      if ( prop( error, "sigmaSlope"  ) != 0 && !prop( error, "sigmaSlopeFixed"  ) )
        vals$sigmaSlope = prop( error, "sigmaSlope"  )
      vals
    }) |> unlist( use.names = FALSE )

  # ── Fisher matrix ─────────────────────────────────────────────────────────────
  fisherMatrix   = prop( fim, "fisherMatrix" )
  allColumnNames = c( columnNamesMu, columnNamesBeta,
                      columnNamesOmega, columnNamesGamma, columnNamesSigma )

  if ( ncol( fisherMatrix ) != length( allColumnNames ) )
    stop( sprintf(
      "setEvaluationFim: FIM dimension (%d) != allColumnNames length (%d).\nNames: %s",
      ncol( fisherMatrix ), length( allColumnNames ),
      paste( allColumnNames, collapse = ", " )
    ))

  colnames( fisherMatrix ) = allColumnNames
  rownames( fisherMatrix ) = allColumnNames

  fixedEffectsNames    = c( columnNamesMu, columnNamesBeta )
  varianceEffectsNames = c( columnNamesOmega, columnNamesGamma, columnNamesSigma )

  fixedEffects    = fisherMatrix[ fixedEffectsNames,    fixedEffectsNames,    drop = FALSE ]
  varianceEffects = fisherMatrix[ varianceEffectsNames, varianceEffectsNames, drop = FALSE ]

  # ── SE and RSE ────────────────────────────────────────────────────────────────
  SE               = sqrt( diag( chol2inv( chol( fisherMatrix ) ) ) )
  parametersValues = c( muValues, betaValues, omegaValues, gammaValues, sigmaValues )
  RSE              = SE / abs( parametersValues ) * 100

  SEAndRSE_df = data.frame(
    parametersValues = parametersValues,
    SE               = SE,
    RSE              = RSE,
    row.names        = allColumnNames
  )

  # ── Store results ─────────────────────────────────────────────────────────────
  prop( fim, "fisherMatrix"              ) = fisherMatrix
  prop( fim, "fixedEffects"              ) = fixedEffects
  prop( fim, "varianceEffects"           ) = varianceEffects
  prop( fim, "condNumberFixedEffects"    ) = cond( fixedEffects )
  prop( fim, "condNumberVarianceEffects" ) = cond( varianceEffects )
  prop( fim, "SEAndRSE" ) = list(
    SE       = SEAndRSE_df[ , c( "parametersValues", "SE"  ) ],
    RSE      = SEAndRSE_df[ , c( "parametersValues", "RSE" ) ],
    SEAndRSE = SEAndRSE_df
  )

  fim
}

# ==============================================================================
#' showFIM: show the Fim in the R console.
#' @name showFIM
#' @param fim An object \code{PopulationFim} giving the Fim.
#' @return Printed console output (invisibly returns fim).
#' @export
# ==============================================================================

method( showFIM, PopulationFim ) = function( fim ) {

  SEAndRSE                  = prop( fim, "SEAndRSE" )
  fisherMatrix              = prop( fim, "fisherMatrix" )
  fixedEffects              = prop( fim, "fixedEffects" )
  varianceEffects           = prop( fim, "varianceEffects" )
  condNumberFixedEffects    = prop( fim, "condNumberFixedEffects" )
  condNumberVarianceEffects = prop( fim, "condNumberVarianceEffects" )
  dcriterionValue           = Dcriterion( fim )
  determinantValue          = det( fisherMatrix )

  cat( "\n*************************************** \n" )
  cat( " Population Fisher Matrix \n" )
  cat( "*************************************** \n\n" )
  print( fisherMatrix )

  cat( "\n*************************************** \n" )
  cat( " Fixed effects (\u03bc) \n" )
  cat( "*************************************** \n\n" )
  print( fixedEffects )

  cat( "\n*************************************** \n" )
  cat( " Variance components (\u03c9\u00B2, \u03b3\u00B2, \u03c3\u00B2) \n" )
  cat( "*************************************** \n\n" )
  print( varianceEffects )

  cat( "\n********************************************* \n" )
  cat( " determinant, condition numbers and d-criterion  \n" )
  cat( "*********************************************** \n\n" )
  cat( "determinant:", determinantValue, "\n" )
  cat( "D-criterion:", dcriterionValue,  "\n" )
  cat( "Conditional number (fixed effects):",     condNumberFixedEffects,    "\n" )
  cat( "Conditional number (variance components):", condNumberVarianceEffects, "\n" )

  cat( "\n*************************************** \n" )
  cat( " Parameters estimation \n" )
  cat( "*************************************** \n\n" )
  print( SEAndRSE$SEAndRSE )

  if ( any( grepl( "\u03b3\u00B2", rownames( fisherMatrix ) ) ) ) {
    cat( "\n*************************************** \n" )
    cat( " Legend: \n" )
    cat( " \u03bc  = fixed effects (population means) \n" )
    cat( " \u03c9\u00B2 = inter-individual variability (IIV) \n" )
    cat( " \u03b3\u00B2 = inter-occasion variability (IOV) \n" )
    cat( " \u03c3\u00B2 = residual error variance \n" )
    cat( "*************************************** \n\n" )
  }

  invisible( fim )
}

# ==============================================================================
#' plotSEFIM: barplot for the SE
#' @name plotSEFIM
#' @param fim        An object \code{PopulationFim} giving the Fim.
#' @param evaluation An object \code{Evaluation} giving the evaluation of the model.
#' @return A ggplot bar chart of SE by parameter type.
#' @export
# ==============================================================================

method( plotSEFIM, list( PopulationFim, PFIMProject ) ) = function( fim, evaluation ) {
  .plotFimBars( fim, evaluation, "SE" )
}

# ==============================================================================
#' plotRSEFIM: barplot for the RSE
#' @name plotRSEFIM
#' @param fim        An object \code{PopulationFim} giving the Fim.
#' @param evaluation An object \code{Evaluation} giving the evaluation of the model.
#' @return A ggplot bar chart of RSE by parameter type.
#' @export
# ==============================================================================

method( plotRSEFIM, list( PopulationFim, PFIMProject ) ) = function( fim, evaluation ) {
  .plotFimBars( fim, evaluation, "RSE" )
}

# ==============================================================================
#' tablesForReport: generate the tables for the HTML report.
#' @name tablesForReport
#' @param fim        An object \code{PopulationFim} giving the Fim.
#' @param evaluation An object \code{Evaluation} giving the evaluation of the model.
#' @return List: fixedEffectsTable, varianceEffectsTable, FIMCriteriaTable, SEAndRSETable.
#' @export
# ==============================================================================

method( tablesForReport, list( PopulationFim, PFIMProject ) ) = function( fim, evaluation ) {

  SEAndRSE                  = prop( fim, "SEAndRSE" )$SEAndRSE
  fisherMatrix              = prop( fim, "fisherMatrix" )
  fixedEffects              = as.matrix( prop( fim, "fixedEffects" ) )
  varianceEffects           = prop( fim, "varianceEffects" )
  condNumberFixedEffects    = prop( fim, "condNumberFixedEffects" )
  condNumberVarianceEffects = prop( fim, "condNumberVarianceEffects" )

  # FIX: renamed to dcriterionValue to avoid shadowing the Dcriterion() function.
  dcriterionValue  = Dcriterion( fim )
  determinantValue = det( fisherMatrix )

  parameters      = prop( evaluation, "modelParameters" )
  modelError      = prop( evaluation, "modelError" )
  modelCovariates = prop( evaluation, "modelCovariates" )
  greek           = .GREEK_LATEX

  has_IOV = any( map_dbl( parameters, ~ pluck( .x, "gamma", .default = 0 ) ) > 0 )
  has_cov = length( modelCovariates ) > 0

  # Column names (LaTeX)
  columnNamesMu    = paste0( .muNames(    parameters, greek[ "mu"    ] ), "}$" )
  # FIX: beta names included (dimension mismatch when covariates present).
  columnNamesBeta  = if ( has_cov )
    paste0( greek[ "beta" ], str_remove( .betaColumnNames( modelCovariates ), "^beta_" ), "}$" )
  else character( 0L )
  columnNamesOmega = paste0( .omegaNames( parameters, greek[ "omega" ] ), "}$" )
  columnNamesGamma = if ( has_IOV )
    paste0( .gammaNames( parameters, greek[ "gamma" ] ), "}$" )
  else character( 0L )
  columnNamesSigma = modelError |>
    map( function( error ) {
      output = prop( error, "output" )
      c(
        if ( prop( error, "sigmaInter" ) != 0 && !prop( error, "sigmaInterFixed" ) )
          paste0( greek[ "sigma" ], "{inter}}_{", output, "}$" ),
        if ( prop( error, "sigmaSlope"  ) != 0 && !prop( error, "sigmaSlopeFixed"  ) )
          paste0( greek[ "sigma" ], "{slope}}_{", output, "}$" )
      )
    }) |> unlist( use.names = FALSE )

  fixedEffectsNames              = c( columnNamesMu, columnNamesBeta )
  colnames( fixedEffects   )     = fixedEffectsNames
  rownames( fixedEffects   )     = fixedEffectsNames
  varianceNames                  = c( columnNamesOmega, columnNamesGamma, columnNamesSigma )
  colnames( varianceEffects )    = varianceNames
  rownames( varianceEffects )    = varianceNames

  fixedEffectsTable = fixedEffects |>
    kbl() |>
    kable_styling( bootstrap_options = "hover", full_width = FALSE,
                   position = "center", font_size = 13 )

  varianceEffectsTable = varianceEffects |>
    kbl() |>
    kable_styling( bootstrap_options = "hover", full_width = FALSE,
                   position = "center", font_size = 13 )

  FIMCriteriaTable = data.frame(
    determinant     = determinantValue,
    dcriterion      = dcriterionValue,
    FixedEffects    = condNumberFixedEffects,
    VarianceEffects = condNumberVarianceEffects
  ) |>
    kbl( col.names = c( "", "", "Fixed effects", "Variance effects" ), align = "c" ) |>
    add_header_above( c( "determinant" = 1, "d-criterion" = 1, "Condition number" = 2 ) ) |>
    kable_styling( bootstrap_options = "hover", full_width = FALSE,
                   position = "center", font_size = 13 )

  SEAndRSETable = data.frame(
    c( columnNamesMu, columnNamesBeta, columnNamesOmega, columnNamesGamma, columnNamesSigma ),
    round( SEAndRSE, 3 )
  ) |>
    kbl( col.names = c( "Parameters", "Parameter values", "SE", "RSE (%)" ),
         align = "c", row.names = FALSE ) |>
    kable_styling( bootstrap_options = "hover", full_width = FALSE,
                   position = "center", font_size = 13 )

  list(
    fixedEffectsTable    = fixedEffectsTable,
    varianceEffectsTable = varianceEffectsTable,
    FIMCriteriaTable     = FIMCriteriaTable,
    SEAndRSETable        = SEAndRSETable
  )
}

# ==============================================================================
# Report generation methods
# FIX: outputFile and outputPath added as explicit parameters (were undefined).
# FIX: params = list(tablesForReport = tablesForReport) — pass the object, not
#      the string "tablesForReport".
# ==============================================================================

#' @name generateReportEvaluation
#' @export
method( generateReportEvaluation, PopulationFim ) = function( fim, tablesForReport,
                                                              outputFile, outputPath ) {
  rmarkdown::render(
    input       = .reportTemplatePath( "EvaluationPopulationFIM.rmd" ),
    output_file = outputFile,
    output_dir  = outputPath,
    params      = list( tablesForReport = tablesForReport )
  )
}

#' @name generateReportOptimization
#' @export
method( generateReportOptimization, list( PopulationFim, MultiplicativeAlgorithm ) ) =
  function( fim, optimizationAlgorithm, tablesForReport, outputFile, outputPath ) {
    rmarkdown::render(
      input       = .reportTemplatePath( "OptimizationMultiplicativeAlgorithmPopulationFIM.rmd" ),
      output_file = outputFile,
      output_dir  = outputPath,
      params      = list( tablesForReport = tablesForReport )
    )
  }

#' @name generateReportOptimization
#' @export
method( generateReportOptimization, list( PopulationFim, FedorovWynnAlgorithm ) ) =
  function( fim, optimizationAlgorithm, tablesForReport, outputFile, outputPath ) {
    rmarkdown::render(
      input       = .reportTemplatePath( "OptimizationFedorovWynnAlgorithmPopulationFIM.rmd" ),
      output_file = outputFile,
      output_dir  = outputPath,
      params      = list( tablesForReport = tablesForReport )
    )
  }

#' @name generateReportOptimization
#' @export
method( generateReportOptimization, list( PopulationFim, SimplexAlgorithm ) ) =
  function( fim, optimizationAlgorithm, tablesForReport, outputFile, outputPath ) {
    rmarkdown::render(
      input       = .reportTemplatePath( "OptimizationSimplexAlgorithmPopulationFIM.rmd" ),
      output_file = outputFile,
      output_dir  = outputPath,
      params      = list( tablesForReport = tablesForReport )
    )
  }

#' @name generateReportOptimization
#' @export
method( generateReportOptimization, list( PopulationFim, PSOAlgorithm ) ) =
  function( fim, optimizationAlgorithm, tablesForReport, outputFile, outputPath ) {
    rmarkdown::render(
      input       = .reportTemplatePath( "OptimizationPSOAlgorithmPopulationFIM.rmd" ),
      output_file = outputFile,
      output_dir  = outputPath,
      params      = list( tablesForReport = tablesForReport )
    )
  }

#' @name generateReportOptimization
#' @export
method( generateReportOptimization, list( PopulationFim, PGBOAlgorithm ) ) =
  function( fim, optimizationAlgorithm, tablesForReport, outputFile, outputPath ) {
    rmarkdown::render(
      input       = .reportTemplatePath( "OptimizationPGBOAlgorithmPopulationFIM.rmd" ),
      output_file = outputFile,
      output_dir  = outputPath,
      params      = list( tablesForReport = tablesForReport )
    )
  }
