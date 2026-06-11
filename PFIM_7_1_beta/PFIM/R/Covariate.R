#' @description The class \code{Covariate} represents the Covariate.
#' @title Covariate
#' @param name Character vector specifying the covariate name.
#' @param effects List specifying the effects.
#' @export

Covariate = new_class("Covariate",
                      properties = list(
                        name = class_character,
                        effects = class_list
                      ),
                      constructor = function(name, effects = list(),
                                             categories = NULL, categoriesProportions = NULL,
                                             sequences = NULL, sequencesProportions = NULL) {

                        # Détection automatique du type de covariate
                        if (!is.null(categories)) {
                          if (!is.null(sequences) && !is.null(sequencesProportions)) {
                            # CategoricalCovariateWithIOV
                            return(CategoricalCovariateWithIOV(
                              name = name,
                              categories = categories,
                              sequences = sequences,
                              sequencesProportions = sequencesProportions,
                              effects = effects ))
                          } else if (!is.null(categoriesProportions)) {
                            # CategoricalCovariate simple
                            return(CategoricalCovariate(
                              name = name,
                              categories = categories,
                              categoriesProportions = categoriesProportions,
                              effects = effects
                            ))
                          }
                        }

                        # Si pas de catégories spécifiées, créer un Covariate de base
                        new_object(Covariate,
                                   name = name,
                                   effects = effects )
                      }
)

getCovariateEffects = new_generic( "getCovariateEffects", c( "covariate" ) )
createEffectVector = new_generic( "createEffectVector", c( "covariate" ) )
getCategoryOfReference = new_generic( "getCategoryOfReference", c( "covariate" ) )

# ==============================================================================
#' getCategoryOfReference: get the category of reference
#' @name getCategoryOfReference
#' @param categoricalCovariate A object of class \code{CategoricalCovariate} giving the categorical covariate.
#' @return A string giving the category of reference.
#' @export
# ==============================================================================

method( getCategoryOfReference, Covariate ) = function( covariate ) {
  categoryOfReference = pluck( prop( covariate, "categories" ), 1 )
  return( categoryOfReference )
}

# ==============================================================================
#' createEffectVector: create the vector for the effect
#' @name createEffectVector
#' @param covariate A object \code{Covariate} giving the covariate
#' @param category A string giving the category of the covariate
#' @param effectVector A vector of double named with the parameters names
#' @return The vector effectVector
#' @export
# ==============================================================================

method( createEffectVector, Covariate ) = function( covariate, category, effectVector ) {

  if ( category == getCategoryOfReference( covariate ) || !category %in% names( prop( covariate, "effects" ) ) ) {
    return( effectVector )
  }

  covariateEffects = prop( covariate, "effects" )
  effectVector[names( covariateEffects[[category]] )] = covariateEffects[[category]]

  return( effectVector )
}



