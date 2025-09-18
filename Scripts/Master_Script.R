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
require(janitor)

options(scipen = 7) #I don't like viewing things in scientific notation
mem.maxVSize(vsize = 30000) #Set max memory size to allow space for large distance matrices

# Raw Data Extraction/Preparation -----------------------------------------

#Extract PDB files
pseudohelixList <- lapply(1:36, function(x) {
  read.pdb(
    file = str_glue("Input/Files_From_Analysis_Cluster/Pseudohelix_{x}/Refine_1/Phenix_refine_001.pdb"), 
    rm.alt = FALSE
  )
}) #Extract pseudohelix PDBs
wedgeList <- lapply(1:36, function(x) {
  read.pdb(
    file = str_glue("Input/Files_From_Analysis_Cluster/Wedge_{x}/Refine_1/Phenix_refine_001.pdb"), 
    rm.alt = FALSE
  )
}) #Extract wedge PDBs

#Extract the atom element of the PDB files into new lists
pseudohelixAtoms <- lapply(pseudohelixList, function(pdb) pdb$atom) #Extract atoms from pseudohelices
wedgeAtoms <- lapply(wedgeList, function(pdb) pdb$atom) #Extract atoms from wedges

#Run dose analysis to prepare dose state vectors
source("Scripts/CheckingForHeterogeneity.R")
source("Scripts/RADDOSE.R")

#Prep vectors containing all density-weighted dose values
pseudohelixDose <- filter(dwds, datasetType == "Pseudohelices")$dwd_MGy
wedgeDose <- filter(dwds, datasetType == "Wedges")$dwd_MGy

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

source("Scripts/OccupancyAnalysis.R")
source("Scripts/DistanceAnalysis.R")
source("Scripts/AngleAnalysis.R")
source("Scripts/StructureStats.R")

# Data Visualization ------------------------------------------------------

source("Scripts/TrendVisualization.R")
source("Scripts/PlotGeneration.R")
source("Scripts/WedgeCVAnalysis.R")

# Data write-out ----------------------------------------------------------

#Save all the plots
save_plots <- function(parameter, height = 5, width = 7) {
  for(theme in names(ggplots)) {
    for(i in names(ggplots[[theme]][[parameter]])) {
      if("ggplot" %in% class(ggplots[[theme]][[parameter]][[i]])) {
        dataset_type = if(i == "Pseudohelices") {"Pseudohelix"} else if(i == "Wedges") {"Wedge"} else if(i == "DWDs") {"DWDs"} else if(i == "CrystalStats") {"CrystalStats"}
        this_parameter = if(parameter %in% c("Angles", "Distances", "Occupancies")) {parameter} else {""}
        ggsave(
          filename = str_glue("Output/Plots/{theme}/{dataset_type}_{parameter}.svg"), 
          plot = ggplots[[theme]][[parameter]][[i]], 
          height = height, 
          width = width
        )
      } else if("list" %in% class(ggplots[[theme]][[parameter]][[i]])) {
        for(j in names(ggplots[[theme]][[parameter]][[i]])) {
          dataset_type = if(j == "Pseudohelices") {"Pseudohelix"} else if(j == "Wedges") {"Wedge"}
          ggsave(
            filename = str_glue("Output/Plots/{theme}/{dataset_type}_{parameter}{i}.svg"), 
            plot = ggplots[[theme]][[parameter]][[i]][[j]], 
            height = height, 
            width = width
          )
        }
      }
    }
  }
}

save_plots('Occupancies')
save_plots('BFactors')
save_plots('Distances')
save_plots('Angles')
save_plots('CVs')
save_plots('Dose')
save_plots('Stats', height = 8)

#Save all associated PDBs
write.pdb(pdb = OccupancyColoredPDB, file = "Output/ColoredPDBs/OccupancyColoredPDB.pdb")

#Save all the tables
for (i in 1:length(regressionSummaries)) {
  write.csv(
    regressionSummaries[[i]], 
    str_glue("Output/RegressionTables/RegressionSummary_{names(regressionSummaries)[i]}.csv")
  )
}

#Make PowerPoint containing figures
new_figure_slide <- function(ppt, plot_list) {
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
  return(ppt)
} #Function to insert slides with 

doc <- read_pptx(path = "/Users/sm9/Desktop/Template.pptx") %>%
  layout_default("Title and Content") %>%
  new_figure_slide(ppt = ., plot_list = ggplots$Dark$Dose) %>% 
  new_figure_slide(ppt = ., plot_list = ggplots$Dark$Stats) %>% 
  new_figure_slide(ppt = ., plot_list = ggplots$Dark$Occupancies) %>% 
  new_figure_slide(ppt = ., plot_list = ggplots$Dark$Distances) %>%
  new_figure_slide(ppt = ., plot_list = ggplots$Dark$Angles) %>% 
  new_figure_slide(ppt = ., plot_list = ggplots$Dark$CVs)
print(doc, target = "/Users/sm9/Desktop/Example Figures.pptx")
