# Packages ----------------------------------------------------------------

require(reticulate)
require(tidyverse)
require(bio3d)
require(plyr)
require(geometry)
require(emmeans)
require(officedown)
require(officer)
require(RcppCNPy)
require(ggh4x)
library(ggnewscale)

py_require("numpy")
py_require("mrcfile")
py_require("matplotlib")

options(scipen = 7) #I don"t like viewing things in scientific notation
mem.maxVSize(vsize = 30000) #Set max memory size to allow space for large distance matrices

# Raw Data Extraction/Preparation -----------------------------------------

#Extract PDB files
pseudohelixList <- lapply(1:36, function(x) {
  read.pdb(
    file = paste0("Input/Files_From_Analysis_Cluster/Pseudohelix_", x, "/Refine_1/Phenix_refine_001.pdb"), 
    rm.alt = FALSE
  )
}) #Extract pseudohelix PDBs
wedgeList <- lapply(1:36, function(x) {
  read.pdb(
    file = paste0("Input/Files_From_Analysis_Cluster/Wedge_", x, "/Refine_1/Phenix_refine_001.pdb"), 
    rm.alt = FALSE
  )
}) #Extract wedge PDBs

#Extract the atom element of the PDB files into new lists
pseudohelixAtoms <- lapply(pseudohelixList, function(pdb) pdb$atom) #Extract atoms from pseudohelices
wedgeAtoms <- lapply(wedgeList, function(pdb) pdb$atom) #Extract atoms from wedges

#Prep vectors containing all density-weighted dose values
pseudohelixDose <- c(0.959494770, 2.309835566, 3.099525168, 3.614838606, 3.998769111, 
                     4.305138794, 4.559751096, 4.774831928, 4.972192543, 5.169593890, 
                     5.362558116, 5.529478983, 5.690346875, 5.855073470, 6.006623963, 
                     6.144268815, 6.277907374, 6.424709533, 6.582981515, 6.723683491, 
                     6.848334088, 6.979721043, 7.097238290, 7.227514626, 7.375691563, 
                     7.522691199, 7.680733524, 7.875180309, 8.110508433, 8.344207821, 
                     8.606089844, 8.943259379, 9.363393076, 9.883614072, 10.67258175, 
                     11.97970325)
wedgeDose       <- c(6.700685286, 6.720294865, 6.729773805, 6.743861225, 6.746327023, 
                     6.719767433, 6.711379946, 6.702595256, 6.677506338, 6.629260722, 
                     6.603718009, 6.595284247, 6.588494019, 6.592423462, 6.618326270, 
                     6.628786101, 6.641042627, 6.660685626, 6.679410310, 6.700307956, 
                     6.713330933, 6.723102252, 6.730564695, 6.703977105, 6.698935005, 
                     6.690346661, 6.667656327, 6.623201255, 6.600886756, 6.596310344, 
                     6.591970625, 6.595671753, 6.629001804, 6.644251692, 6.664629154, 
                     6.680498800)

# Global Functions and Objects --------------------------------------------

emtrends.coefficient <- function(model, regressor) {
  c(
    emtrends(model, "Molecule", regressor) %>%
      test() %>% .[, 2], 
    emtrends(model, "Molecule", regressor) %>%
      pairs() %>% summary() %>% .[, 2]
  )
} #Extract coefficients from emtrends
emtrends.pvalue <- function(model, regressor) {
  c(
    emtrends(model, "Molecule", regressor) %>%
      test() %>% .[, 6], 
    emtrends(model, "Molecule", regressor) %>%
      pairs() %>%summary() %>% .[, 6]
  )
} #Extract p-values from emtrends
emtrends.se <- function(model, regressor) {
  c(
    emtrends(model, "Molecule", regressor) %>%
      summary() %>% .[, 3], 
    emtrends(model, "Molecule", regressor) %>%
      pairs() %>%summary() %>% .[, 3]
  )
} #Extract standard errors from emtrends

#Make blank lists to organize all visualization data
multipleRegressions <- list()
regressionSummaries <- list()
ggplots <- list()

# Data Analysis -----------------------------------------------------------

source("Scripts/R/OccupancyAnalysis.R")
source("Scripts/R/BFactorAnalysis.R")
source("Scripts/R/DistanceAnalysis.R")
source("Scripts/R/AngleAnalysis.R")
source("Scripts/R/ElectronDensityAnalysis.R")

py_run_file("Scripts/Python/PseuodohelixMapGeneration.py")
py_run_file("Scripts/Python/WedgeMapGeneration.py")

# Data Visualization ------------------------------------------------------

source("Scripts/R/TrendVisualization.R")
source("Scripts/R/PlotGeneration.R")
source("Scripts/R/WedgeCVAnalysis.R")

# Data write-out ----------------------------------------------------------

#Save all the plots
save.plots <- function(key) {
  for(theme in names(ggplots)) {
    for(i in names(ggplots[[theme]][[key]])) {
      if("ggplot" %in% class(ggplots[[theme]][[key]][[i]])) {
        ggsave(
          filename = paste0(
            "Output/Plots/",
            theme, "/",
            ifelse(
              i == "Pseudohelices",
              'Pseudohelix',
              'Wedge'
            ), "_",
            key, 
            ".svg"
          ), 
          plot = ggplots[[theme]][[key]][[i]], 
          height = 5, 
          width = 7
        )
      } else if("list" %in% class(ggplots[[theme]][[key]][[i]])) {
        for(j in names(ggplots[[theme]][[key]][[i]])) {
          ggsave(
            filename = paste0(
              "Output/Plots/",
              theme, "/",
              ifelse(
                j == "Pseudohelices",
                'Pseudohelix',
                'Wedge'
              ),
              key, "_",
              i,
              ".svg"
            ), 
            plot = ggplots[[theme]][[key]][[i]][[j]], 
            height = 5, 
            width = 7
          )
        }
      }
    }
  }
}

save.plots('Occupancies')
save.plots('BFactors')
save.plots('Distances')
save.plots('Angles')
save.plots('CVs')

#Save all associated PDBs
write.pdb(pdb = OccupancyColoredPDB, file = "Output/ColoredPDBs/OccupancyColoredPDB.pdb")
write.pdb(pdb = bFactorColoredPDB, file = "Output/ColoredPDBs/BFactorColoredPDB.pdb")
for(i in 1:length(pseudohelixList)) {
  write.pdb(
    pdb = pseudohelixList[[i]], 
    file = paste0("Output/PDBs/Pseudohelix", i, ".pdb")
  )
}
for(i in 1:length(wedgeList)) {
  write.pdb(
    pdb = wedgeList[[i]], 
    file = paste0("Output/PDBs/Wedge", i, ".pdb")
  )
}

#Save all the tables
for (i in 1:length(regressionSummaries)) {
  write.csv(
    regressionSummaries[[i]], 
    paste0(
      "Output/RegressionTables/",
      deparse(substitute(regressionSummaries)), 
      "$", names(regressionSummaries)[i], 
      ".csv"
    )
  )
}

#Make PowerPoint containing figures
new.figure.slide <- function(ppt, plot_list) {
  for (i in plot_list) {
    if ("ggplot" %in% class(i)) {
      ppt <- ppt %>%
        add_slide() %>%
        ph_with(value = i, location = ph_location_fullsize())
    } else if("list" %in% class(i)) {
      for(j in i) {
        ppt <- ppt %>%
          add_slide() %>%
          ph_with(value = j, location = ph_location_fullsize())
      }
    }
  }
  ppt
} #Function to insert slides with 

doc <- read_pptx(path = "/Users/sm9/Desktop/Template.pptx") %>%
  layout_default("Title and Content") %>%
  new.figure.slide(ppt = ., plot_list = ggplots$Dark$Occupancies) %>% 
  new.figure.slide(ppt = ., plot_list = ggplots$Dark$BFactors) %>% 
  new.figure.slide(ppt = ., plot_list = ggplots$Dark$Distances) %>%
  new.figure.slide(ppt = ., plot_list = ggplots$Dark$Angles) %>% 
  new.figure.slide(ppt = ., plot_list = ggplots$Dark$CVs)
print(doc, target = "/Users/sm9/Desktop/Example Figures.pptx")
