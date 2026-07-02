# Makes a fake example for Paul's simulation study: criteria

rm(list=ls())
source('FunctionsCriteria.R')

# Dirs
parmDir <- '../parms/'

# Parms
n <- 1e4; d <- 10; B <- 1
dimName <- paste0('-d', d, '-n', n)
load(paste0(parmDir, 'ParmsLists', dimName, '.Rdata'))

# Scenarios
scenarioList <- list(Decentrage = list(mild = c(k=5, l=1, rho1=0.3), 
                                       strong = c(k=10, l=1, rho1=0.3),
                                       stronger = c(k=50, l=1, rho1=0.3)),
                     Concentration = list(mild = c(k=0, l=0.5, rho1=0.3), 
                                          strong = c(k=0, l=0.1, rho1=0.3),
                                          stronger = c(k=0, l=0.01, rho1=0.3)), 
                     Dilatation = list(mild = c(k=0, l=2, rho1=0.3), 
                                       strong = c(k=0, l=10, rho1=0.3),
                                       stronger = c(k=0, l=100, rho1=0.3)), 
                     Deformation = list(mild = c(k=0, l=1, rho1=0.5), 
                                        strong = c(k=0, l=1, rho1=0.7),
                                        stronger = c(k=0, l=1, rho1=0.95)),
                     DeformationConcentration = list(mild = c(k=0, l=0.5, rho1=0.5), 
                                                     strong = c(k=0, l=0.1, rho1=0.7), 
                                                     stronger = c(k=0, l=0.01, rho1=0.95)), 
                     DeformationDilatation = list(mild = c(k=0, l=2, rho1=0.5), 
                                                  strong = c(k=0, l=10, rho1=0.7), 
                                                  stronger = c(k=0, l=100, rho1=0.95)), 
                     DecentrageDilatation = list(mild = c(k=5, l=2, rho1=0.3), 
                                                 strong = c(k=10, l=10, rho1=0.3),
                                                 stronger = c(k=50, l=100, rho1=0.3))
                     )
save(scenarioList, file=paste0(parmDir, 'ScenarioList', dimName, '.Rdata'))
