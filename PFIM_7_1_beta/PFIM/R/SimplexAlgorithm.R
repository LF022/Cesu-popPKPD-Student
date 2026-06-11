#' @description The class \code{SimplexAlgorithm} implements the Simplex algorithm.
#' @title SimplexAlgorithm
#' @inheritParams Optimization
#' @param pctInitialSimplexBuilding A numeric giving the pctInitialSimplexBuilding.
#' @param maxIteration  A numeric giving the maxIteration.
#' @param tolerance  A numeric giving the tolerance.
#' @param seed  A numeric giving the seed.
#' @param showProcess A Boolean giving the showProcess.
#' @include Optimization.R
#' @export

SimplexAlgorithm = new_class( "SimplexAlgorithm", package = "PFIM", parent = Optimization,

                          properties = list( pctInitialSimplexBuilding = new_property(class_double, default = numeric(0)),
                                             maxIteration = new_property(class_double, default = numeric(0)),
                                             seed = new_property(class_double, default = numeric(0)),
                                             tolerance = new_property(class_double, default = numeric(0)),
                                             showProcess = new_property(class_logical, default = FALSE ) ) )

fisherSimplex = new_generic( "fisherSimplex", c( "optimizationObject" ) )

# ==============================================================================
#' Compute the fun.amoeba
#' @name fun.amoeba
#' @param p parameter p
#' @param y parameter y
#' @param ftol parameter ftol
#' @param itmax parameter itmax
#' @param funk parameter funk
#' @param outcomes The model outcomes.
#' @param data parameter data
#' @param showProcess Boolean.
#' @return fun.amoeba
#' @export
# ==============================================================================

fun.amoeba = function(p,y,ftol,itmax,funk,outcomes,data,showProcess){
  alpha = 1.0
  beta = 0.5
  gamma = 2.0

  eps = 1.e-10
  mpts = nrow(p)
  iter = 0
  contin = T
  converge = F
  results = data.frame()

  while (contin) {

    if ( showProcess == TRUE )
    {
      print( paste0('iter = ',iter))
      print( paste0('Criterion = ', 1/min(y) ) )
    }

    results = rbind( c( iter , 1/min(y) ) , results)

    # ============================================================
    ## First we must determine which point is the highest (worst),
    ## next highest, and lowest (best).
    # ============================================================

    ord=sort.list(y)
    ilo=ord[1]
    ihi=ord[mpts]
    inhi=ord[mpts-1]

    # ===================================================================
    ## Compute the fractional range from highest to lowest and return if
    ## satisfactory.
    # ===================================================================

    rtol=2.*abs(y[ihi]-y[ilo])/(abs(y[ihi])+abs(y[ilo])+eps)

    if ((rtol<ftol)||(iter==itmax))
    {
      contin=F
      converge=T
      if (iter==itmax) { converge=F }
    } else {
      if (iter==itmax) cat('Amoeba exceeding maximum iterations.\n')
      iter=iter+1

      # ======================================================================
      ## Begin a new iteration.  Compute the vector average of all points
      ## except the highest, i.e. the center of the face of the simplex across
      ## from the high point.  We will subsequently explore along the ray from
      ## the high point through that center.
      # ======================================================================

      pbar=matrix(p[-ihi,],(mpts-1),(mpts-1))
      pbar=apply(pbar,2,mean)

      # ======================================================================
      ## Extrapolate by a factor alpha through the face, i.e. reflect the
      ## simplex from the high point.  Evaluate the function at the reflected
      ## point.
      # ======================================================================

      pr=(1+alpha)*pbar-alpha*p[ihi,]
      ypr=funk(data,pr,outcomes)
      if (ypr<=y[ilo])
      {
        # ======================================================================
        ## Gives a result better than the best point, so try an additional
        ## extrapolation by a factor gamma, and check out the function there.
        # ======================================================================
        prr=gamma*pr+(1-gamma)*pbar
        yprr=funk(data,prr,outcomes)
        if (yprr<y[ilo])
        {
          # ========================================================================
          ## The additional extrapolation succeeded, and replaces the highest point.
          # ========================================================================
          p[ihi,]=prr
          y[ihi]=yprr
        } else {
          # ========================================================================
          ## The additional extrapolation failed, but we can still use the
          ## reflected point.
          # ========================================================================
          p[ihi,]=pr
          y[ihi]=ypr
        }
      } else {
        if (ypr>=y[inhi])
        {
          # ========================================================================
          ## The reflected point is worse than the second-highest.
          # ========================================================================
          if (ypr<y[ihi])
          {
            # ========================================================================
            ## If it's better than the highest, then replace the highest,
            # ========================================================================
            p[ihi,]=pr
            y[ihi]=ypr
          }
          # ========================================================================
          ## but look for an intermediate lower point, in other words, perform a
          ## contraction of the simplex along one dimension.  Then evaluate the
          ## function.
          # ========================================================================
          prr=beta*p[ihi,]+(1-beta)*pbar
          yprr=funk(data,prr,outcomes)
          if (yprr<y[ihi])
          {
            # ========================================================================
            ## Contraction gives an improvement, so accept it.
            # ========================================================================
            p[ihi,]=prr
            y[ihi]=yprr
          } else {
            # ========================================================================
            ## Can't seem to get rid of that high point.  Better contract around the
            ## lowest (best) point.
            # ========================================================================
            p=0.5*(p+matrix(p[ilo,],nrow=nrow(p),ncol=ncol(p),byrow=T))
            for (j in 1:length(y)) {
              if (j!=ilo) {
                y[j]=funk(data,p[j,],outcomes)
              }
            }
          }
        } else {
          # ========================================================================
          ## We arrive here if the original reflection gives a middling point.
          ## Replace the old high point and continue
          # ========================================================================
          p[ihi,]=pr
          y[ihi]=ypr
        }
      }
    }
  }
  return(list(p=p,y=y,iter=iter,converge=converge, results = results ))
} # end function amoeba

# ==============================================================================
#' Compute the fisher.simplex
#' @name fisherSimplex
#' @param simplex A list giving the parameters of the simplex.
#' @param optimizationObject An object \code{Optimization}.
#' @param outcomes A vector giving the outcomes of the arms.
#' @return A list giving the results of the optimization.
#' @export
# ==============================================================================

method( fisherSimplex, Optimization ) = function( optimizationObject, simplex, outcomes )
{
  samplingTimeConstraintsForContinuousOptimization = list( )

  # designs and arms
  designs = prop( optimizationObject, "designs" )

  for ( design in designs )
  {
    designName = prop( design, "name" )
    arms = prop( design, "arms" )

    # set sampling times in each arm
    armsList = list()

    for ( arm in arms )
    {
      armName = prop( arm, "name" )

      # get constraints
      samplingTimesArms = list()
      samplingTimesConstraints = prop( arm, "samplingTimesConstraints" )
      samplingTimes = prop( arm, "samplingTimes" )

      for ( outcome in outcomes[[armName]] )
      {
        # get the samplings times by design-arm-outcome
        namesSamplings = toString( c( designName, armName, outcome ) )
        indexSamplingTimes = which( names( simplex ) == namesSamplings )
        newSamplings = simplex[indexSamplingTimes]
        names( newSamplings ) = NULL

        # check the constraints on the samplings times
        samplingTimesConstraint = keep( samplingTimesConstraints, ~ prop( .x,"outcome" ) == outcome ) %>% pluck(1)

        samplingTimeConstraintsForContinuousOptimization[[ armName ]][[ outcome ]] =
          checkSamplingTimeConstraintsForMetaheuristic( samplingTimesConstraint, arm, newSamplings, outcome )

        samplingTimes = samplingTimes %>%
          modify_at(
            .at = which( map_chr( ., ~ prop( .x, "outcome" ) ) == outcome ),
            .f = ~ {
              prop( .x, "samplings" ) = newSamplings
              .x
            })
      }
      prop( arm, "samplingTimes" ) = samplingTimes
      armsList = append( armsList, arm )
    }

    # set new arms to the design
    prop( design, "arms" ) = armsList

    samplingTimeConstraintsForContinuousOptimization = unlist( samplingTimeConstraintsForContinuousOptimization )

    # if constraints are satisfied evaluate design
    if( all( samplingTimeConstraintsForContinuousOptimization ) == TRUE )
    {
      # evaluate the fim
      evaluationFIM = Evaluation( name = "",
                                  modelEquations = prop( optimizationObject, "modelEquations" ),
                                  modelParameters = prop( optimizationObject, "modelParameters" ),
                                  modelError = prop( optimizationObject, "modelError" ),

                                  modelCovariates = prop( optimizationObject, "modelCovariates" ),
                                  modelCovariatesEquation = prop( optimizationObject, "modelCovariatesEquation" ),
                                  numberOfOccasions = prop( optimizationObject, "numberOfOccasions" ),

                                  fimType = prop( optimizationObject, "fimType" ),
                                  outputs = prop( optimizationObject, "outputs" ),
                                  designs = list( design ),
                                  odeSolverParameters = prop( optimizationObject, "odeSolverParameters" ) )

      evaluationFIM =  run( evaluationFIM )

      # get D-criterion
      fim = prop( evaluationFIM, "fim" )
      Dcriterion = 1/Dcriterion( fim )

    }else{
      Dcriterion = 1
    }
  }

  # in case of Inf
  Dcriterion[is.infinite(Dcriterion)] = 1.0
  return( Dcriterion )
}

# ==============================================================================
#' Optimization SimplexAlgorithm
#' @name optimizeDesign
#' @param optimizationObject A object \code{Optimization}.
#' @param optimizationAlgorithm A object \code{SimplexAlgorithm}.
#' @return The object \code{optimizationObject} with the slots updated.
#' @export
# ==============================================================================

method( optimizeDesign, list( Optimization, SimplexAlgorithm ) ) = function( optimizationObject, optimizationAlgorithm )
{
  # get simplex parameters
  optimizerParameters = prop( optimizationObject, "optimizerParameters")
  showProcess = optimizerParameters$showProcess
  pctInitialSimplexBuilding = optimizerParameters$pctInitialSimplexBuilding
  tolerance = optimizerParameters$tolerance
  maxIteration = optimizerParameters$maxIteration

  # get the design to be optimized
  designs = prop( optimizationObject, "designs" )
  design = pluck( designs, 1 )
  optimalDesign =  pluck( designs, 1 )
  designName = prop( design, "name" )

  # check validity of the sampling Times Constraints
  checkValiditySamplingConstraint( design )

  # in case of partial sampling constraints set the new arms with the sampling constraints
  design = setSamplingConstraintForOptimization( design )

  # get the arms
  arms = prop( design, "arms" )

  # set the new designs in the optimizationObject
  prop( optimizationObject, "designs" ) = list( design )

  # get arm outcomes
  outcomes = map( set_names( arms, map_chr( arms, ~ prop( .x, 'name' ) ) ), ~ {
    map_chr( prop( .x, "samplingTimesConstraints" ), ~ prop( .x, "outcome"))
  })

  # create the initial simplex
  namesSamplingsSimplex = list()
  samplingsSimplex = list()
  initialSimplex = list()

  k=1
  for ( arm in arms )
  {
    armName = prop( arm, "name" )
    samplingTimesConstraints = prop( arm, "samplingTimesConstraints" )

    for ( outcome in outcomes[[armName]] )
    {
      samplingTimesConstraint = keep( samplingTimesConstraints, ~ prop(.x,"outcome") == outcome ) %>% pluck(1)
      samplingsSimplex[[k]] = prop( samplingTimesConstraint, "initialSamplings" )
      namesSamplingsSimplex[[k]] = rep( toString( c( designName, armName, outcome ) ), length( samplingsSimplex[[k]] ) )
      k=k+1
    }
  }

  samplingsSimplex = unlist( samplingsSimplex )
  samplingsSimplex = matrix( rep( samplingsSimplex,length( samplingsSimplex ) + 1 ), ncol = length( samplingsSimplex ), byrow=TRUE )
  colnames( samplingsSimplex ) = unlist( namesSamplingsSimplex )
  samplingsSimplex = t( apply( samplingsSimplex, 1, sort ) )

  # percentage initial simplex
  for ( i in 1:dim( samplingsSimplex )[2] )
  {
    samplingsSimplex[i+1,i] = samplingsSimplex[i+1,i]*( 1-pctInitialSimplexBuilding/100 )
  }

  # evaluate criteria
  y = c()
  for ( i in 1:dim( samplingsSimplex )[1] )
  {
    y[i] = fisherSimplex( optimizationObject, samplingsSimplex[i,], outcomes )
  }

  # Run the simplex
  opti = fun.amoeba( samplingsSimplex, y, tolerance, maxIteration, fisherSimplex, outcomes, data = optimizationObject, showProcess )

  # get the results of the simplex
  optimalDCriteria = opti$y
  samplingTimesSimplex = opti$p
  indexOptimalDCriteria = which( opti$y == min( opti$y ) )
  results = opti$results
  optimalsamplingTimes = samplingTimesSimplex[indexOptimalDCriteria,]

  # set optimal design
  namesDesignArmOutcome = unique( names( optimalsamplingTimes ) )
  arms = prop( optimalDesign, "arms" )
  armsList = list()

  for ( arm in arms )
  {
    armName = prop( arm, "name" )

    listOfSamplingTimes = list()

    for ( outcome in outcomes[[armName]] )
    {
      namesSamplings = toString(c( designName, armName, outcome ))
      indexSamplingTimes = which( names( optimalsamplingTimes ) == namesSamplings )
      samplings = optimalsamplingTimes[indexSamplingTimes]
      names(samplings) = NULL
      samplingTimes = SamplingTimes( outcome, samplings = samplings )
      listOfSamplingTimes = append( listOfSamplingTimes, samplingTimes )
    }
    prop( arm, "samplingTimes" ) = listOfSamplingTimes
    armsList = append( armsList, arm )
  }

  # set optimal arms
  prop( optimalDesign, "arms" ) = armsList

  # evaluate the optimal design
  evaluationOptimalDesign = Evaluation( name = "",
                                        modelEquations = prop( optimizationObject, "modelEquations" ),
                                        modelParameters = prop( optimizationObject, "modelParameters" ),
                                        modelError = prop( optimizationObject, "modelError" ),

                                        modelCovariates = prop( optimizationObject, "modelCovariates" ),
                                        modelCovariatesEquation = prop( optimizationObject, "modelCovariatesEquation" ),
                                        numberOfOccasions = prop( optimizationObject, "numberOfOccasions" ),

                                        designs = list( optimalDesign ),
                                        fimType = prop( optimizationObject, "fimType" ),
                                        outputs = prop( optimizationObject, "outputs" ),
                                        odeSolverParameters = prop( optimizationObject, "odeSolverParameters" ) )

  evaluationOptimalDesign = run( evaluationOptimalDesign )

  # evaluate the initial design
  evaluationInitialDesign = Evaluation( name = "",
                                        modelEquations = prop( optimizationObject, "modelEquations" ),
                                        modelParameters = prop( optimizationObject, "modelParameters" ),
                                        modelError = prop( optimizationObject, "modelError" ),

                                        modelCovariates = prop( optimizationObject, "modelCovariates" ),
                                        modelCovariatesEquation = prop( optimizationObject, "modelCovariatesEquation" ),
                                        numberOfOccasions = prop( optimizationObject, "numberOfOccasions" ),

                                        designs = list( design ),
                                        fimType = prop( optimizationObject, "fimType" ),
                                        outputs = prop( optimizationObject, "outputs" ),
                                        odeSolverParameters = prop( optimizationObject, "odeSolverParameters" ) )

  evaluationInitialDesign = run( evaluationInitialDesign )

  # set the results in evaluation
  prop( optimizationObject, "optimisationDesign" ) = list( evaluationInitialDesign = evaluationInitialDesign, evaluationOptimalDesign = evaluationOptimalDesign )
  prop( optimizationObject, "optimisationAlgorithmOutputs" ) = list( "optimizationAlgorithm" = optimizationAlgorithm, "optimalArms" = armsList )
  return( optimizationObject )
}

# ==============================================================================
#' constraintsTableForReport: table of the SimplexAlgorithm constraints for the report.
#' @name constraintsTableForReport
#' @param optimizationAlgorithm A object \code{SimplexAlgorithm}.
#' @param arms List of the arms.
#' @return The table for the constraints in the arms.
#' @export
# ==============================================================================

method( constraintsTableForReport, SimplexAlgorithm ) = function( optimizationAlgorithm, arms  )
{
  armsConstraints = map( pluck( arms, 1 ) , ~ getArmConstraints( .x, optimizationAlgorithm ) )
  armsConstraints = map_dfr( armsConstraints, ~ map_df(.x, ~ as.data.frame(.x, stringsAsFactors = FALSE)))
  colnames( armsConstraints ) = c( "Arms name" , "Number of subjects", "Outcome", "Initial samplings", "Samplings windows", "Number of times by windows","Min sampling" )
  armsConstraintsTable = kbl( armsConstraints, align = c( "l","c","c","c","c","c","c") ) %>% kable_styling( bootstrap_options = c(  "hover" ), full_width = FALSE, position = "center", font_size = 13 )
  return( armsConstraintsTable )
}




