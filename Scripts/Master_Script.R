# Packages ----------------------------------------------------------------

library(cowplot)
library(tidyverse)
library(bio3d)
library(geometry)
library(emmeans)
library(ggh4x)
library(ggnewscale)
library(plotly)
library(htmlwidgets)
library(janitor)
library(systemfonts)

options(scipen = 7) #I don't like viewing things in scientific notation
mem.maxVSize(vsize = 100000) #Set max memory size to allow space for large distance matrices

# Raw Data Extraction/Preparation -----------------------------------------

#Extract PDB files
pseudohelixList <- lapply(1:36, function(x) {
  bio3d::read.pdb(
    file = stringr::str_glue("Input/structures/pseudohelix_{x}.pdb"), 
    rm.alt = FALSE
  )
}) #Extract pseudohelix PDBs
wedgeList <- lapply(1:36, function(x) {
  bio3d::read.pdb(
    file = stringr::str_glue("Input/structures/wedge_{x}.pdb"), 
    rm.alt = FALSE
  )
}) #Extract wedge PDBs

atoms <- list( # Make list of atom and xyz indices
  cu_a = bio3d::atom.select(pseudohelixList[[1]], resno = 230, chain = "A"),
  nterm_a = bio3d::atom.select(pseudohelixList[[1]], resno = 1, chain = "A", elety = "N"),
  his1nd_a = bio3d::atom.select(pseudohelixList[[1]], resno = 1, chain = "A", elety = "ND1"),
  his1ne_a =  bio3d::atom.select(pseudohelixList[[1]], resno = 1, chain = "A", elety = "NE2"),
  his1cg_a = bio3d::atom.select(pseudohelixList[[1]], resno = 1, chain = "A", elety = "CG"),
  his84nd_a = bio3d::atom.select(pseudohelixList[[1]], resno = 84, chain = "A", elety = "ND1"),
  his84ne_a = bio3d::atom.select(pseudohelixList[[1]], resno = 84, chain = "A", elety = "NE2"),
  his84cg_a = bio3d::atom.select(pseudohelixList[[1]], resno = 84, chain = "A", elety = "CG"),
  h2oax_a = list(atom = dplyr::filter(pseudohelixList[[1]]$atom, resno == 416, chain == "A", alt == "A")$eleno), # atom.select isn't sensitive to altcons
  h2oeq_a = list(atom = dplyr::filter(pseudohelixList[[1]]$atom, resno == 282, chain == "B", alt == "A")$eleno),
  tyr168oh_a = bio3d::atom.select(pseudohelixList[[1]], resno = 168, chain = "A", elety = "OH"),
  oxy_a = bio3d::atom.select(pseudohelixList[[1]], resno = 232, chain = "A", elety = "O2"),
  co2_a = bio3d::atom.select(pseudohelixList[[1]], resno = 228, chain = "B", elety = "C"),
  
  cu_b = bio3d::atom.select(pseudohelixList[[1]], resno = 229, chain = "B"),
  nterm_b = bio3d::atom.select(pseudohelixList[[1]], resno = 1, chain = "B", elety = "N"),
  his1nd_b = bio3d::atom.select(pseudohelixList[[1]], resno = 1, chain = "B", elety = "ND1"),
  his1ne_b =  bio3d::atom.select(pseudohelixList[[1]], resno = 1, chain = "B", elety = "NE2"),
  his1cg_b = bio3d::atom.select(pseudohelixList[[1]], resno = 1, chain = "B", elety = "CG"),
  his84nd_b = bio3d::atom.select(pseudohelixList[[1]], resno = 84, chain = "B", elety = "ND1"),
  his84ne_b = bio3d::atom.select(pseudohelixList[[1]], resno = 84, chain = "B", elety = "NE2"),
  his84cg_b = bio3d::atom.select(pseudohelixList[[1]], resno = 84, chain = "B", elety = "CG"),
  h2oax_b = list(atom = dplyr::filter(pseudohelixList[[1]]$atom, resno == 288, chain == "B", alt == "A")$eleno),
  h2oeq_b = list(atom = dplyr::filter(pseudohelixList[[1]]$atom, resno == 288, chain == "A", alt == "A")$eleno),
  tyr168oh_b = bio3d::atom.select(pseudohelixList[[1]], resno = 168, chain = "B", elety = "OH"),
  oxy_b = bio3d::atom.select(pseudohelixList[[1]], resno = 231, chain = "A", elety = "O2"),
  glu30_b = bio3d::atom.select(pseudohelixList[[1]], resno = 30, chain = "B", elety = "CD", alt = "A")
)

#Extract the atom element of the PDB files into new lists
pseudohelixAtoms <- lapply(pseudohelixList, function(pdb) pdb$atom) #Extract atoms from pseudohelices
wedgeAtoms <- lapply(wedgeList, function(pdb) pdb$atom) #Extract atoms from wedges

#Run dose analysis to prepare dose state vectors
source("Scripts/RADDOSE.R")

#Prep vectors containing all density-weighted dose values
pseudohelixDose <- dplyr::filter(dwds, datasetType == "Pseudohelices")$dwd_MGy
wedgeDose <- dplyr::filter(dwds, datasetType == "Wedges")$dwd_MGy

# Global Functions and Objects --------------------------------------------

emtrends.coefficient <- function(model, regressor) {
  c(
    emmeans::emtrends(model, "Molecule", regressor) %>%
      test() %>% .[, 2], 
    emmeans::emtrends(model, "Molecule", regressor) %>%
      pairs() %>% summary() %>% .[, 2]
  )
} #Extract coefficients from emtrends
emtrends.pvalue <- function(model, regressor) {
  c(
    emmeans::emtrends(model, "Molecule", regressor) %>%
      test() %>% .[, 6], 
    emmeans::emtrends(model, "Molecule", regressor) %>%
      pairs() %>% summary() %>% .[, 6]
  )
} #Extract p-values from emtrends
emtrends.se <- function(model, regressor) {
  c(
    emmeans::emtrends(model, "Molecule", regressor) %>%
      summary() %>% .[, 3], 
    emmeans::emtrends(model, "Molecule", regressor) %>%
      pairs() %>% summary() %>% .[, 3]
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

source("Scripts/TablePrep.R")
source("Scripts/CheckingForHeterogeneity.R")
source("Scripts/PlotGeneration.R")

# Data write-out ----------------------------------------------------------
 
#Save all the plots
save_plots <- function(parameter, outerheight = 10, outerwidth = 16, innerheight = 5, innerwidth = 8) {
  for(theme in names(ggplots)) {
    for(i in names(ggplots[[theme]][[parameter]])) {
      if("ggplot" %in% class(ggplots[[theme]][[parameter]][[i]])) {
        dataset_type = if(i == "Pseudohelices") {"Pseudohelix"} else if(i == "Wedges") {"Wedge"} else if(i == "DWDs") {"DWDs"} else if(i == "CrystalStats") {"CrystalStats"} else if(i == "RMSDs") {"RMSDs"}
        this_parameter = if(parameter %in% c("Angles", "Distances", "Occupancies")) {parameter} else {""}
        ggplot2::ggsave(
          filename = stringr::str_glue("Output/Plots/{theme}/{dataset_type}_{parameter}.svg"), 
          plot = ggplots[[theme]][[parameter]][[i]], 
          height = outerheight, 
          width = outerwidth,
          units = "cm",
          create.dir = TRUE
        )
      } else if("list" %in% class(ggplots[[theme]][[parameter]][[i]])) {
        for(j in names(ggplots[[theme]][[parameter]][[i]])) {
          dataset_type = if(j == "Pseudohelices") {"Pseudohelix"} else if(j == "Wedges") {"Wedge"}
          ggplot2::ggsave(
            filename = stringr::str_glue("Output/Plots/{theme}/{dataset_type}_{parameter}{i}.svg"), 
            plot = ggplots[[theme]][[parameter]][[i]][[j]], 
            height = innerheight, 
            width = innerwidth,
            units = "cm",
            create.dir = TRUE
          )
        }
      }
    }
  }
}

save_plots('Occupancies')
save_plots('Distances')
save_plots('Angles', outerheight = 12)
save_plots('Dose', outerwidth = 8.8, outerheight = 7)
save_plots('Stats', outerwidth = 16, outerheight = 16)
save_plots("Comparisons", outerheight = 16)

#Save all the tables
dir.create("Output/Tables", showWarnings = FALSE)
readr::write_csv(longData$Pseudohelices, "Output/Tables/RegressionSummary_Pseudohelices.csv")
readr::write_csv(longData$Wedges, "Output/Tables/RegressionSummary_Wedges.csv")
readr::write_csv(occ_table, "Output/Tables/all_occupancies.csv")
readr::write_csv(angle_table, "Output/Tables/all_angles.csv")
readr::write_csv(distance_table, "Output/Tables/all_distances.csv")

dir.create("Output/PlotlyHTMLs", showWarnings = FALSE)
htmlwidgets::saveWidget(plotly::as_widget(plotlys$dark$interactive$pseudohelixDoseState), "./Output/PlotlyHTMLs/Pseudohelix.html", selfcontained = FALSE)

