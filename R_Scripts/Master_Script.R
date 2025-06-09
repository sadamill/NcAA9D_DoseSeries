library(reticulate)
library(tidyverse)
library(bio3d)
library(plyr)
library(geometry)
library(emmeans)
library(officedown)
library(officer)
library(RcppCNPy)

py_require('numpy')
py_require('mrcfile')
py_require('matplotlib')

options(scipen = 7)
mem.maxVSize(vsize = 30000) #Set max memory size to allow space for large distance matrices

# Raw Data Extraction/Preparation -----------------------------------------

#Extract PDB files
setwd('/Users/sm9/Documents/17_Meilleur_Lab/NcLPMO9D_Dose_Series_Study/Analysis/Files_From_Analysis_Cluster')
pseudohelixList <- lapply(1:36, function(x) {
  read.pdb(paste0('Pseudohelix_', x,'/Refine_1/Phenix_refine_001.pdb'), rm.alt = FALSE)
}) #Extract pseudohelix PDBs
wedgeList <- lapply(1:36, function(x) {
  read.pdb(paste0('Wedge_', x,'/Refine_1/Phenix_refine_001.pdb'), rm.alt = FALSE)
}) #Extract wedge PDBs

#Extract the atom element of the PDB files into new lists
pseudohelixAtoms <- lapply(pseudohelixList,function(pdb) pdb$atom) #Extract atoms from pseudohelices
wedgeAtoms <- lapply(wedgeList,function(pdb) pdb$atom) #Extract atoms from wedges

#Prep vectors containing all density-weighted dose values
pseudohelixDose <- c(0.95949477,2.309835566,3.099525168,3.614838606,3.998769111,4.305138794,4.559751096,4.774831928,4.972192543,5.16959389,5.362558116,5.529478983,5.690346875,5.85507347,6.006623963,6.144268815,6.277907374,6.424709533,6.582981515,6.723683491,6.848334088,6.979721043,7.09723829,7.227514626,7.375691563,7.522691199,7.680733524,7.875180309,8.110508433,8.344207821,8.606089844,8.943259379,9.363393076,9.883614072,10.67258175,11.97970325)
wedgeDose <- c(6.700685286,6.720294865,6.729773805,6.743861225,6.746327023,6.719767433,6.711379946,6.702595256,6.677506338,6.629260722,6.603718009,6.595284247,6.588494019,6.592423462,6.61832627,6.628786101,6.641042627,6.660685626,6.67941031,6.700307956,6.713330933,6.723102252,6.730564695,6.703977105,6.698935005,6.690346661,6.667656327,6.623201255,6.600886756,6.596310344,6.591970625,6.595671753,6.629001804,6.644251692,6.664629154,6.6804988)

# Global Functions and Objects --------------------------------------------

#Define functions
emtrends.coefficient <- function(model, regressor) {
  c(
    emtrends(model,'Molecule', regressor) %>% 
      test() %>% {.[,2]},
    emtrends(model,'Molecule', regressor) %>% 
      pairs() %>% summary() %>% {.[,2]}
  )
} #Extract coefficients from emtrends
emtrends.pvalue <- function(model, regressor) {
  c(
    emtrends(model,'Molecule', regressor) %>% 
      test() %>% {.[,6]},
    emtrends(model,'Molecule', regressor) %>% 
      pairs() %>% summary() %>% {.[,6]}
  )
} #Extract p-values from emtrends
emtrends.se <- function(model, regressor) {
  c(
    emtrends(model,'Molecule', regressor) %>% 
      summary() %>% {.[,3]},
    emtrends(model,'Molecule', regressor) %>% 
      pairs() %>% summary() %>% {.[,3]}
  )
} #Extract standard errors from emtrends
ggtheme <- function() {
  list(
    scale_color_manual(
      'Molecule',
      labels = c('A','B'),
      values = c('#004c0b','#90188d'),
    ),
      theme_bw(),
      theme(
        legend.position = 'inside',
        legend.position.inside = c(0.85,0.25),
        strip.background = element_rect(fill = 'white'),
        legend.background = element_rect(color = 'black'),
        plot.background = element_blank()
      ),
    coord_cartesian(expand = FALSE)
  )
} #Global ggplot theme

#Make blank lists to contain all visualization data
multipleRegressions <- list()
regressionSummaries <- list()
ggplots <- list()

# Data Analysis -----------------------------------------------------------

setwd('../R_Scripts')
source('OccupancyAnalysis.R')
source('BFactorAnalysis.R')
source('DistanceAnalysis.R')
source('AngleAnalysis.R')
source('WedgeCVAnalysis.R')
source('TrendVisualization.R')
source('ElectronDensityAnalysis.R')

setwd('../Python_Scripts')
py_run_file('PseuodohelixMapGeneration.py')
py_run_file('WedgeMapGeneration.py')

# Data write-out ----------------------------------------------------------

setwd('../Visualizations')

#Save all the plots
save.plots <- function(list) {
  for(i in 1:length(list)) {
    ggsave(
      filename = paste0(deparse(substitute(list)), '$', names(list)[i], '.svg'),
      plot = list[[i]],
      height = 5,
      width = 7
    )
  }
}

save.plots(ggplots$Occupancies)
save.plots(ggplots$BFactors)
save.plots(ggplots$Distances)
save.plots(ggplots$Angles)
save.plots(ggplots$CVs)

#Save all the PDBs
write.pdb(OccupancyColoredPDB,'OccupancyColoredPDB.pdb')
write.pdb(bFactorColoredPDB,'BFactorColoredPDB.pdb')

#Save all the tables
for (i in 1:length(regressionSummaries)) {
  write.csv(
    regressionSummaries[[i]],
    paste0(deparse(substitute(regressionSummaries)), '$', names(regressionSummaries)[i], '.csv')
  )
}

doc <- read_pptx(path = '/Users/sm9/Desktop/Template.pptx') %>% 
  layout_default('Title and Content') %>% 
  add_slide() %>% 
  ph_with(value = ggplots$Occupancies, 
          location = ph_location_fullsize()
  ) %>% 
  add_slide() %>% xgplots$OccupancyTrends, 
          location = ph_location_fullsize()
  ) %>% 
  add_slide() %>% 
  ph_with(value = ggplots$CopperDistances, 
          location = ph_location_fullsize()
  ) %>% 
  add_slide() %>% 
  ph_with(value = ggplots$DioxygenDistances, 
          location = ph_location_fullsize()
  ) %>%
  add_slide() %>% 
  ph_with(value = ggplots$DistanceTrends, 
          location = ph_location_fullsize()
  ) %>% 
  add_slide() %>% 
  ph_with(value = ggplots$Angles, 
          location = ph_location_fullsize()
  ) %>% 
  add_slide() %>% 
  ph_with(value = ggplots$AngleTrends, 
          location = ph_location_fullsize()
  ) %>% 
  add_slide() %>% 
  ph_with(value = ggplots$Densities, 
          location = ph_location_fullsize()
  ) %>% 
  add_slide() %>% 
  ph_with(value = ggplots$DensityTrends, 
          location = ph_location_fullsize()
  ) %>% 
  add_slide() %>% 
  ph_with(value = ggplots$WedgeCVs, 
          location = ph_location_fullsize()
  ) %>% 
  print(doc, target = '/Users/sm9/Desktop/Example Figures.pptx')
