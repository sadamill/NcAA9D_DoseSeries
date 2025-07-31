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
require(ggnewscale)
require(plotly)
require(webshot2)
require(htmlwidgets)

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

#Run dose analysis to prepare dose state vectors
source("Scripts/R/RADDOSE.R")

#Prep vectors containing all density-weighted dose values
pseudohelixDose <- filter(dwds, datasetType == "Pseudohelices")$dwd_MGy
wedgeDose <- filter(dwds, datasetType == "Wedges")$dwd_MGy
write_csv(
  data.frame(pseudohelixDose, wedgeDose),
  file = "./Output/SlopeMaps/DWDs.csv"
)

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

ggsave("Output/Plots/Light/DWDs.svg", height = 5, width = 7, plot = ggplots$Light$Dose$DWDs)
ggsave("Output/Plots/Dark/DWDs.svg", height = 5, width = 7, plot = ggplots$Dark$Dose$DWDs)

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
  new.figure.slide(ppt = ., plot_list = ggplots$Dark$Dose) %>% 
  new.figure.slide(ppt = ., plot_list = ggplots$Dark$Occupancies) %>% 
  new.figure.slide(ppt = ., plot_list = ggplots$Dark$BFactors) %>% 
  new.figure.slide(ppt = ., plot_list = ggplots$Dark$Distances) %>%
  new.figure.slide(ppt = ., plot_list = ggplots$Dark$Angles) %>% 
  new.figure.slide(ppt = ., plot_list = ggplots$Dark$CVs)
print(doc, target = "/Users/sm9/Desktop/Example Figures.pptx")
