#' @description The class \code{CategoricalCovariateWithIOV} implements the CategoricalCovariateWithIOV
#' @title CategoricalCovariate
#' @param name name
#' @description The class \code{CategoricalCovariateWithIOV} implements the CategoricalCovariateWithIOV
#' @title CategoricalCovariateWithIOV
#' @param name name
#' @param categories categories
#' @param sequences sequences
#' @param seq_proportions sequence proportions
#' @param effects effects
#' @include Covariate.R
#' @export

CategoricalCovariateWithIOV = new_class("CategoricalCovariateWithIOV",
                                        parent = Covariate,
                                        properties = list(
                                          categories = class_character,
                                          sequences = class_list,
                                          sequencesProportions = class_double
                                        ),
                                        constructor = function(name, categories, sequences, sequencesProportions,
                                                               effects = list() ) {
                                          if ( sum(sequencesProportions) !=  1 ) {
                                            stop("The proportions must sum to 1.")
                                          }

                                          names(sequences) = paste0("sequence_", 1:length(sequences))

                                          new_object(CategoricalCovariateWithIOV,
                                                     name = name,
                                                     effects = effects,
                                                     categories = categories,
                                                     sequences = sequences,
                                                     sequencesProportions = sequencesProportions)
                                        } )

# ==============================================================================
#' Compute the effects of a covariate with iov.
#' @name getCovariateEffects
#' @param covariate A object \code{CategoricalCovariateWithIOV} giving the covariate with iov
#' @return The effect of a covariate with iov.
#' @export
# ==============================================================================

method( getCovariateEffects, CategoricalCovariateWithIOV ) = function( covariate, effectVector ) {

  covariateSequences = prop( covariate, "sequences" )
  numberOfOccasions = length( pluck( covariateSequences, 1 ) )

  covariateEffects = covariateSequences %>%
    imap(function( covariateSequence, sequenceName ) {
      seq_len( numberOfOccasions ) %>%
        map( function( occasionIndex ) {
          covariateCategory = covariateSequence[occasionIndex]
          createEffectVector( covariate, covariateCategory, effectVector )
        }) %>% set_names( paste0( "occasion_", seq_len( numberOfOccasions ) ) )
    }) %>% set_names( paste0( "sequence_", seq_len( length( covariateSequences ) ) ) )

  return( covariateEffects )
}


