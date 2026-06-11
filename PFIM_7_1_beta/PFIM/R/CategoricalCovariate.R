#' @description The class \code{CategoricalCovariate} implements the Categorical distribution.
#' @title CategoricalCovariate
#' @param name name
#' @param categories categories
#' @param proportions proportions
#' @param effects effects
#' @include Covariate.R
#' @export

CategoricalCovariate = new_class("CategoricalCovariate",
                                 parent = Covariate,
                                 properties = list(
                                   categories = class_character,
                                   categoriesProportions = class_double ),
                                 constructor = function(name, categories, categoriesProportions, effects = list() ) {
                                   if ( sum(categoriesProportions) != 1) {
                                     stop("The proportions must sum to 1.")
                                   }
                                   new_object(CategoricalCovariate,
                                              name = name,
                                              effects = effects,
                                              categories = categories,
                                              categoriesProportions = categoriesProportions)
                                 } )


# ==============================================================================
#' Compute the effects of a covariate without iov.
#' @name getCovariateEffects
#' @param categoricalCovariate A object \code{CategoricalCovariate} giving the catecovariate
#' @return The effect of a covariate without iov.
#' @export
# ==============================================================================

method( getCovariateEffects, CategoricalCovariate ) = function( covariate, nullVector ) {

  covariateCategories = prop( covariate, "categories" )

  covariateEffects = covariateCategories %>%
    map(~ createEffectVector( covariate, .x, nullVector ) ) %>%
    set_names(covariateCategories)

  return( covariateEffects )
}


















