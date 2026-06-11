# ==============================================================================
#' CovariateTest: évaluation de la pertinence clinique des covariables
#'   dans les modèles NLME.
#'
#' Slots :
#'   @significance  — test de significativité (Wald bilatéral, H0: beta = 0)
#'                    appliqué à tous les mu et tous les beta
#'   @nonRelevance  — test de non-pertinence TOST (H1: ratio ∈ [0.80, 1.25])
#'                    appliqué aux beta dont exp(beta) ∈ [0.80, 1.25]
#'   @relevance     — test de pertinence (H1: ratio ∉ [0.80, 1.25])
#'                    appliqué aux beta dont exp(beta) < 0.80 ou > 1.25
#'
#' @include Evaluation.R
#' @export
# ==============================================================================

# ==============================================================================
#' covariateTest: generic function
#' @name covariateTest
#' @export
# ==============================================================================

covariateTest = new_generic( "covariateTest", "pfimproject" )

CovariateTest = new_class( "CovariateTest", package = "PFIM",
                               properties = list(
                                 significance = new_property( class_any, default = data.frame() ),
                                 nonRelevance = new_property( class_any, default = data.frame() ),
                                 relevance    = new_property( class_any, default = data.frame() )
                               )
)

# ==============================================================================
# Private pure functions — package level
#
# Extracted outside the method so they can be unit-tested independently via
# testthat without instantiating a full Evaluation object.
# All statistical constants (z_half, z_one, Binf, Bsup) are passed explicitly
# so every function is referentially transparent.
# ==============================================================================

# ── Significativité ────────────────────────────────────────────────────────────
# PS = 1 − Φ(z_half − β/SE) + Φ(−z_half − β/SE)
.powerSignificance = function( beta, SE, z_half ) {
  1 - pnorm( z_half - beta / SE ) + pnorm( -z_half - beta / SE )
}

# NSN significance
# SE_S = beta  / (z_half − Φ⁻¹(1−PS))   si beta > 0
#      = −beta / (z_half + Φ⁻¹(PS))     si beta < 0
# NSN  = sigma²_unit / SE_S²
.nRequiredSignificance = function( beta, sigma2_unit, z_half, PS ) {
  if ( beta == 0 ) return( NA_real_ )
  SE_S = if ( beta > 0 ) {
    beta  / ( z_half - qnorm( 1 - PS ) )
  } else {
    -beta / ( z_half + qnorm( PS ) )
  }
  if ( SE_S <= 0 ) return( NA_real_ )
  sigma2_unit / SE_S^2
}

# ── Non-pertinence / TOST ──────────────────────────────────────────────────────
# PNR = Φ(−z_one + (Bsup−β)/SE) − Φ(z_one + (Binf−β)/SE)
# Condition de validité : 2·z_one < (Bsup−Binf)/SE
.powerNonRelevance = function( beta, SE, Binf, Bsup, z_one ) {
  if ( 2 * z_one >= ( Bsup - Binf ) / SE ) return( 0 )
  pnorm( -z_one + ( Bsup - beta ) / SE ) -
    pnorm(  z_one + ( Binf - beta ) / SE )
}

# NSN non-pertinence — root-finder sur PNR(N) = PS.
# Retourne :
#   NA_real_  si beta hors zone d'équivalence (test inapplicable)
#   2         si la puissance est déjà atteinte à N = 2
#   Inf       si la puissance cible est asymptotique (p_max − PS < 0.005) :
#             une solution théorique existe mais nécessiterait un N non réaliste
#   racine    valeur réelle trouvée par uniroot
.nRequiredNonRelevance = function( beta, sigma2_unit, Binf, Bsup, z_one, PS ) {
  if ( beta <= Binf || beta >= Bsup ) return( NA_real_ )

  f     = function( N ) .powerNonRelevance( beta, sqrt( sigma2_unit / N ), Binf, Bsup, z_one ) - PS
  p_max = .powerNonRelevance( beta, 1e-12, Binf, Bsup, z_one )

  if ( f( 2 ) >= 0 ) return( 2 )

  # Puissance maximale insuffisante — solution inexistante.
  if ( p_max < PS ) return( NA_real_ )

  # Puissance cible asymptotique : uniroot convergerait vers un N astronomique
  # (ou échouerait). On signale Inf plutôt que de laisser la recherche diverger.
  if ( p_max - PS < 0.005 ) return( Inf )

  tryCatch(
    uniroot( f, interval = c( 2, 1e7 ), extendInt = "yes" )$root,
    error = function( e ) NA_real_
  )
}

# ── Pertinence ─────────────────────────────────────────────────────────────────
# PR = Φ(−z_one + (Binf−β)/SE) + 1 − Φ(z_one + (Bsup−β)/SE)
.powerRelevance = function( beta, SE, Binf, Bsup, z_one ) {
  pnorm( -z_one + ( Binf - beta ) / SE ) +
    1 - pnorm(  z_one + ( Bsup - beta ) / SE )
}

# NSN pertinence
# beta > Bsup : num = (Bsup−beta) < 0, den = (Φ⁻¹(1−PS) − z_one) < 0 → SE_R > 0
# beta < Binf : num = (Binf−beta) > 0, den = (Φ⁻¹(PS)   + z_one) > 0 → SE_R > 0
# NSN = sigma²_unit / SE_R²
.nRequiredRelevance = function( beta, sigma2_unit, Binf, Bsup, z_one, PS ) {
  if ( beta >= Binf && beta <= Bsup ) return( NA_real_ )
  SE_R = if ( beta > Bsup ) {
    ( Bsup - beta ) / ( qnorm( 1 - PS ) - z_one )
  } else {
    ( Binf - beta ) / ( qnorm( PS ) + z_one )
  }
  if ( is.na( SE_R ) || SE_R <= 0 ) return( NA_real_ )
  sigma2_unit / SE_R^2
}

# ==============================================================================
#' covariateTest: évalue la significativité et la pertinence clinique
#'   des effets de covariables à partir de la FIM.
#'
#' 1. Significativité — slot @significance — mu et beta
#'    H0 : beta = 0,  H1 : beta ≠ 0  (test de Wald bilatéral)
#'    PS  = 1 − Φ(q_{1−α/2} − β/SE) + Φ(−q_{1−α/2} − β/SE)
#'    NSN = N × (SE / SE_S)²
#'
#' 2. Non-pertinence / TOST — slot @nonRelevance — beta seulement
#'    H0 : ratio pertinent (hors [Rinf, Rsup]), H1 : ratio non-pertinent (dedans)
#'    PNR = Φ(−q_{1−α} + (Bsup−β)/SE) − Φ(q_{1−α} + (Binf−β)/SE)
#'    NSN via root-finder
#'
#' 3. Pertinence — slot @relevance — beta seulement
#'    H0 : ratio non-pertinent (dans [Rinf, Rsup]), H1 : ratio pertinent (dehors)
#'    PR  = Φ(−q_{1−α} + (Binf−β)/SE) + 1 − Φ(q_{1−α} + (Bsup−β)/SE)
#'    NSN = N × (SE / SE_R)²
#'
#' @name covariateTest
#' @param pfimproject Un objet \code{Evaluation}.
#' @return Un objet \code{CovariateTest} avec les trois slots.
#' @export
# ==============================================================================

method( covariateTest, Evaluation ) = function( pfimproject,
                                                    thetaL       = log( 0.80 ),
                                                    thetaU       = log( 1.25 ),
                                                    target_power = 0.90,
                                                    alpha        = 0.05 ) {

  # ── Paramètres globaux ───────────────────────────────────────────────────────
  # z_half = q_{1-alpha/2}  (test bilatéral — significance)
  # z_one  = q_{1-alpha}    (test unilatéral — TOST)
  z_half       = qnorm( 1 - alpha / 2 )   # 1.96
  z_one        = qnorm( 1 - alpha       )  # 1.645

  # ── Extraction des estimées et SE depuis la FIM ───────────────────────────────
  fim  = prop( pfimproject, "fim" )
  fim  = setEvaluationFim( fim, pfimproject )
  seDF = prop( fim, "SEAndRSE" )$SEAndRSE   # colonnes : parametersValues, SE, RSE

  allParametersNames  = rownames( seDF )
  betaHat             = set_names( seDF$parametersValues, allParametersNames )
  standardErrorsValue = set_names( seDF$SE,               allParametersNames )
  RSE                 = set_names( seDF$RSE,               allParametersNames )

  # N0 = taille totale du premier design
  N0 = sum( map_dbl( prop( prop( pfimproject, "designs" )[[1L]], "arms" ),
                     ~ prop( .x, "size" ) ) )

  # ── Indices par type de paramètre ────────────────────────────────────────────
  # idx_mu  : noms produits par .muNames()         via .GREEK_CONSOLE["mu"]    = "\u03bc_"
  # idx_beta: noms produits par .betaColumnNames() via paste0("beta_", ...)  = "beta_" (ASCII)
  idx_mu   = which( str_starts( allParametersNames, "\u03bc_" ) )
  idx_beta = which( str_starts( allParametersNames, "beta_"   ) )

  # Ratios exp(beta) pré-calculés une seule fois pour les deux filtres.
  ratioBeta   = exp( betaHat[ idx_beta ] )
  idx_beta_NR = idx_beta[ ratioBeta >= 0.80 & ratioBeta <= 1.25 ]
  idx_beta_R  = idx_beta[ ratioBeta <  0.80 | ratioBeta >  1.25 ]

  # IC 90% sur l'échelle naturelle (ratio) — utilisé dans les slots 2 et 3.
  lowerCI = exp( betaHat - z_one * standardErrorsValue )
  upperCI = exp( betaHat + z_one * standardErrorsValue )

  # Variance unitaire par sujet = SE² × N, pré-calculée pour tous les paramètres.
  sigma2_units = standardErrorsValue^2 * N0

  # ── Slot 1 — Significativité (mu et beta) ────────────────────────────────────
  df_signif = map( c( idx_mu, idx_beta ), function( i ) {
    b     = betaHat[i]
    SE    = standardErrorsValue[i]
    N_req = .nRequiredSignificance( b, sigma2_units[i], z_half, target_power )
    p_N0  = .powerSignificance( b, SE, z_half )

    status = if ( is.na( N_req ) ) {
      "Impossible (beta = 0)"
    } else if ( N_req <= 2 ) {
      sprintf( "OK (N = 2) — puissance a N = %d : %.1f%%", N0, p_N0 * 100 )
    } else {
      sprintf( "N = %d requis — puissance a N = %d : %.1f%%",
               ceiling( N_req ), N0, p_N0 * 100 )
    }

    data.frame(
      Parameter  = allParametersNames[i],
      Value      = round( b,      4 ),
      SE         = round( SE,     4 ),
      RSE        = round( RSE[i], 2 ),
      Power_N0   = round( p_N0 * 100, 1 ),
      N_Required = ceiling( N_req ),
      Status     = status,
      stringsAsFactors = FALSE,
      row.names        = NULL
    )
  }) |> list_rbind()

  # ── Slot 2 — Non-pertinence (beta dont ratio ∈ [0.80, 1.25]) ─────────────────
  df_NR = map( idx_beta_NR, function( i ) {
    b     = betaHat[i]
    SE    = standardErrorsValue[i]
    ratio = exp( b )
    p_max = .powerNonRelevance( b, 1e-12, thetaL, thetaU, z_one )

    N_req = .nRequiredNonRelevance( b, sigma2_units[i], thetaL, thetaU, z_one, target_power )
    p_N0  = .powerNonRelevance( b, SE, thetaL, thetaU, z_one )

    status = if ( is.na( N_req ) ) {
      sprintf( "Impossible — puissance max (N->inf) : %.1f%%", p_max * 100 )
    } else if ( is.infinite( N_req ) ) {
      # Puissance cible asymptotique : la solution theorique existe mais
      # exigerait un N non realisable (p_max - PS < 0.5%).
      sprintf( "Asymptotique — puissance max (N->inf) : %.1f%% (cible : %.0f%%)",
               p_max * 100, target_power * 100 )
    } else if ( N_req <= 2 ) {
      sprintf( "OK (N = 2) — puissance a N = %d : %.1f%%", N0, p_N0 * 100 )
    } else {
      sprintf( "N = %d requis — puissance a N = %d : %.1f%%",
               ceiling( N_req ), N0, p_N0 * 100 )
    }

    data.frame(
      Parameter  = allParametersNames[i],
      Value      = round( b,          4 ),
      SE         = round( SE,         4 ),
      RSE        = round( RSE[i],     2 ),
      Ratio      = round( ratio,      4 ),
      IC90_Inf   = round( lowerCI[i], 4 ),
      IC90_Sup   = round( upperCI[i], 4 ),
      Power_N0   = round( p_N0  * 100, 1 ),
      Power_max  = round( p_max * 100, 1 ),
      N_Required = ceiling( N_req ),
      Status     = status,
      stringsAsFactors = FALSE,
      row.names        = NULL
    )
  }) |> list_rbind()

  # ── Slot 3 — Pertinence (beta dont ratio ∉ [0.80, 1.25]) ─────────────────────
  df_R = map( idx_beta_R, function( i ) {
    b     = betaHat[i]
    SE    = standardErrorsValue[i]
    ratio = exp( b )

    N_req = .nRequiredRelevance( b, sigma2_units[i], thetaL, thetaU, z_one, target_power )
    p_N0  = .powerRelevance( b, SE, thetaL, thetaU, z_one )

    status = if ( is.na( N_req ) ) {
      "Impossible — puissance max insuffisante"
    } else if ( N_req <= 2 ) {
      sprintf( "OK (N = 2) — puissance a N = %d : %.1f%%", N0, p_N0 * 100 )
    } else {
      sprintf( "N = %d requis — puissance a N = %d : %.1f%%",
               ceiling( N_req ), N0, p_N0 * 100 )
    }

    data.frame(
      Parameter  = allParametersNames[i],
      Value      = round( b,          4 ),
      SE         = round( SE,         4 ),
      RSE        = round( RSE[i],     2 ),
      Ratio      = round( ratio,      4 ),
      IC90_Inf   = round( lowerCI[i], 4 ),
      IC90_Sup   = round( upperCI[i], 4 ),
      Power_N0   = round( p_N0 * 100, 1 ),
      N_Required = ceiling( N_req ),
      Status     = status,
      stringsAsFactors = FALSE,
      row.names        = NULL
    )
  }) |> list_rbind()

  CovariateTest(
    significance = df_signif,
    nonRelevance = df_NR,
    relevance    = df_R
  )
}

# ==============================================================================
#' show: affichage console d'un objet CovariateTest
#' @name show
#' @param pfimproject Un objet \code{CovariateTest}.
#' @export
# ==============================================================================

method( show, CovariateTest ) = function( pfimproject ) {

  print_slot = function( df, title ) {
    cat( "\n", strrep( "=", 60 ), "\n", sep = "" )
    cat( " ", title, "\n" )
    cat( strrep( "=", 60 ), "\n\n", sep = "" )
    if ( nrow( df ) == 0 ) {
      cat( "  (aucun parametre dans ce slot)\n" )
    } else {
      print( df, row.names = FALSE )
    }
  }

  print_slot( prop( pfimproject, "significance" ), "SLOT 1 — Significativite statistique" )
  print_slot( prop( pfimproject, "nonRelevance" ), "SLOT 2 — Non-pertinence clinique"      )
  print_slot( prop( pfimproject, "relevance"    ), "SLOT 3 — Pertinence clinique"           )

  invisible( pfimproject )
}

# ==============================================================================
#' saveCovariateTest: sauvegarde les trois slots dans un fichier texte
#' @name saveCovariateTest
#' @param relevanceResults Un objet \code{CovariateTest}.
#' @param file_name Nom du fichier de sortie.
#' @param folder Dossier de destination (optionnel).
#' @export
# ==============================================================================

saveCovariateTest = function( relevanceResults, file_name, folder = NULL ) {

  # Formatage sans dépendance externe : utils::capture.output + print natif.
  # Cela supprime la dépendance à knitr (paquet lourd) pour une simple mise en
  # forme tabulaire dans un fichier texte.
  format_slot = function( df, title ) {
    c(
      strrep( "=", 60 ),
      paste0( "  ", title ),
      strrep( "=", 60 ),
      "",
      if ( nrow( df ) == 0 ) {
        "  (aucun parametre dans ce slot)"
      } else {
        utils::capture.output( print( df, row.names = FALSE ) )
      },
      ""
    )
  }

  output_lines = c(
    strrep( "=", 60 ),
    "   PFIM — Clinical Relevance",
    strrep( "=", 60 ),
    "",
    format_slot( prop( relevanceResults, "significance" ), "SLOT 1 — Significativite statistique" ),
    format_slot( prop( relevanceResults, "nonRelevance" ), "SLOT 2 — Non-pertinence clinique"      ),
    format_slot( prop( relevanceResults, "relevance"    ), "SLOT 3 — Pertinence clinique"           ),
    strrep( "=", 60 ),
    "   END OF FILE",
    strrep( "=", 60 )
  )

  if ( !is.null( folder ) ) {
    if ( !dir.exists( folder ) ) {
      warning( sprintf(
        "saveCovariateTest: le dossier '%s' n'existait pas et a ete cree.",
        folder
      ) )
      dir.create( folder, recursive = TRUE )
    }
    file_name = file.path( folder, file_name )
  }

  writeLines( output_lines, file_name )
  message( "Sauvegarde : ", file_name )
}

# ==============================================================================
#' tost: alias de covariateTest pour compatibilite ascendante
#' @name tost
#' @param pfimproject Un objet \code{Evaluation}.
#' @return Un objet \code{CovariateTest}.
#' @export
# ==============================================================================

method( tost, Evaluation ) = function( pfimproject,
                                       thetaL       = log( 0.80 ),
                                       thetaU       = log( 1.25 ),
                                       target_power = 0.90,
                                       alpha        = 0.05 ) {
  covariateTest( pfimproject,
                     thetaL       = thetaL,
                     thetaU       = thetaU,
                     target_power = target_power,
                     alpha        = alpha )
}
