###############################################################
# PFIM 7.1 SCRIPT PK-PD MODEL 
# -------------------------------------------------------------
# This script implements:
# - A 1-compartment IV bolus PK model
# - An Emax PD model driven by PK concentrations
# - Population design with predefined sampling times
#
###############################################################

library(ggplot2)
library(gridExtra)

devtools::load_all("D:/PMX/PROJECTS/LUCIE/CesuPopPKPD/PFIM_7_1_beta/PFIM")
cesuPath = "D:/PMX/PROJECTS/LUCIE/CesuPopPKPD"


###############################################################
# 1. MODEL DEFINITION
modelEquations = list(  "RespPK" = "dose_RespPK/V * exp(- Cl/V * t_RespPK) ",
                        "RespPD" = "E0 +  Emax * RespPK/( RespPK + C50 )" )

###############################################################
# 2. POPULATION PARAMETERS
# model parameters
modelParameters = list(
  ModelParameter( name = "V",    distribution = LogNormal( mu = 0.2, omega = sqrt(0.25) )),
  ModelParameter( name = "Cl",   distribution = LogNormal( mu = 0.05, omega = sqrt(0.25) ) ),
  ModelParameter( name = "E0",   distribution = LogNormal( mu = 1, omega = sqrt(0.09) ) ),
  ModelParameter( name = "C50",  distribution = LogNormal( mu = 1, omega = sqrt(0.09) ) ),
  ModelParameter( name = "Emax", distribution = LogNormal( mu = 4, omega = sqrt(0.09) ) ) )

###############################################################
# 3. ERROR MODELS
# error model
errorModelRespPK = Proportional( output = "RespPK", sigmaSlope = 0.3 )
errorModelRespPD = Constant( output = "RespPD", sigmaInter = 0.15 )
modelError = list( errorModelRespPK, errorModelRespPD )

###############################################################
# 4. DESIGN DEFINITION

## administration 
administrationRespPK = Administration( outcome = "RespPK", time = c( 0 ), dose = c( 1 ) )

## sampling times
samplingTimesRespPK = SamplingTimes( outcome = "RespPK", samplings = c( 0.167,  6, 12 ) )
samplingTimesRespPD = SamplingTimes( outcome = "RespPD", samplings = c( 0.167, 6, 12, 20 ) )

## arms definition
arm1 = Arm( name = "BrasTest",
            size = 100,
            administrations = list( administrationRespPK ) ,
            samplingTimes   = list( samplingTimesRespPK, samplingTimesRespPD ) )


design1 = Design( name = "design1", arms = list( arm1 ) )

###############################################################
# 5. FISHER INFORMATION MATRIX COMPUTATION

evaluationPopFim = Evaluation( name = "",
                               modelParameters = modelParameters,
                               modelEquations = modelEquations,
                               modelError =  modelError,
                               designs = list( design1 ),
                               outputs = list( "RespPK", "RespPD" ),
                               fimType = 'population' )

evaluationPopFim = run( evaluationPopFim )

show(evaluationPopFim)

getFisherMatrix(evaluationPopFim)$fisherMatrix
getDcriterion(evaluationPopFim)
getRSE(evaluationPopFim)

###############################################################
# 6. PLOTS AND REPORTS
plotOptions = list( unitTime = c("h"), unitOutcomes = c("Concentration","Effect")  )
plot_evaluation = plotEvaluation(evaluationPopFim, plotOptions)

plot_designs = grid.arrange(plot_evaluation$design1$BrasTest$RespPK + ylab("Concentration"),
                            plot_evaluation$design1$BrasTest$RespPD + ylab("Effect"),
                            ncol = 2)

save_my_plot(plot(plot_designs), 1, 0.75, "PFIMExample_plotdesigns", "cesuPopPK", cesuPath)

plot_sensitivity = plotSensitivityIndices(evaluationPopFim, plotOptions)



# Save results & Reports
outputPath = cesuPath
outputFile = "popFIM.html"
outputFileRDS = "popFIM.rds"
saveRDS(evaluationPopFim, file = file.path(outputPath, outputFileRDS))
Report( evaluationPopFim, outputPath, outputFile, plotOptions )

###############################################################
# 7. OPTIMISATION

# constraints to optimise only PD samplings
administrationConstraintsRespPK = AdministrationConstraints( outcome = "RespPK", 
                                                             doses = list( c( 1 )) )

samplingConstraintsRespPK = SamplingTimeConstraints( outcome = "RespPK",
                                                     initialSamplings = c( 0.167,  6, 12 ),
                                                     numberOfsamplingsOptimisable = 3,
                                                     fixedTimes = c(0.167, 6, 12))

samplingConstraintsRespPD = SamplingTimeConstraints( outcome = "RespPD",
                                                     initialSamplings = c( 0.167, 1, 2, 6, 8, 10, 12, 16, 20 ),
                                                     numberOfsamplingsOptimisable = 3 )

arm_opti_1= Arm( name = "BrasTest1",
                 size = 100,
                 administrations = list( administrationRespPK ),
                 samplingTimes   = list( samplingTimesRespPK, samplingTimesRespPD ),
                 administrationsConstraints = list( administrationConstraintsRespPK ),
                 samplingTimesConstraints = list(  samplingConstraintsRespPK, samplingConstraintsRespPD ) )

design_opti_1 = Design( name = "design1", arms = list( arm_opti_1), numberOfArms = 100 )

# optimizationPopFIM
optimizationPopFIMSettings = Optimization( name = "PKPD populationFIM",
                                           modelEquations = modelEquations,
                                           modelParameters = modelParameters,
                                           modelError = modelError,
                                           optimizer = "MultiplicativeAlgorithm",
                                           optimizerParameters = list( lambda = 0.99,
                                                                       numberOfIterations = 1000,
                                                                       weightThreshold = 0.01,
                                                                       delta = 1e-04, showProcess = T ),
                                           designs = list( design_opti_1 ),
                                           fimType = "population",
                                           outputs = list( "RespPK", "RespPD" ) )

resOptimizationPopFIM = run( optimizationPopFIMSettings )
show(resOptimizationPopFIM)
