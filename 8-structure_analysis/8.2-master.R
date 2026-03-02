# Packages ----------------------------------------------------------------

library(factoextra)
library(nlme)
library(cowplot)
library(tidyverse)
library(bio3d)
library(geometry)
library(emmeans)
library(ggnewscale)
library(janitor)
library(systemfonts)
library(PCDimension)
library(ggforce)

# Raw Data Extraction/Preparation -----------------------------------------

#Extract PDB files
pseudohelixList <- lapply(1:36, function(x) {
  bio3d::read.pdb(
    file = stringr::str_glue("8-structure_analysis/input/structures/pseudohelix_{x}.pdb"), 
    rm.alt = FALSE
  )
}) #Extract pseudohelix PDBs
wedgeList <- lapply(1:36, function(x) {
  bio3d::read.pdb(
    file = stringr::str_glue("8-structure_analysis/input/structures/wedge_{x}.pdb"), 
    rm.alt = FALSE
  )
}) #Extract wedge PDBs

atoms <- list( # Make list of atom and xyz indices
  cu_a       = bio3d::atom.select(pseudohelixList[[1]], resno = 230, chain = "A"),
  nterm_a    = bio3d::atom.select(pseudohelixList[[1]], resno = 1, chain = "A", elety = "N"),
  his1nd_a   = bio3d::atom.select(pseudohelixList[[1]], resno = 1, chain = "A", elety = "ND1"),
  his1ne_a   =  bio3d::atom.select(pseudohelixList[[1]], resno = 1, chain = "A", elety = "NE2"),
  his1cg_a   = bio3d::atom.select(pseudohelixList[[1]], resno = 1, chain = "A", elety = "CG"),
  his84nd_a  = bio3d::atom.select(pseudohelixList[[1]], resno = 84, chain = "A", elety = "ND1"),
  his84ne_a  = bio3d::atom.select(pseudohelixList[[1]], resno = 84, chain = "A", elety = "NE2"),
  his84cg_a  = bio3d::atom.select(pseudohelixList[[1]], resno = 84, chain = "A", elety = "CG"),
  h2oax_a    = list(atom = dplyr::filter(pseudohelixList[[1]]$atom, resno == 416, chain == "A", alt == "A")$eleno), # atom.select isn't sensitive to altcons
  h2oeq_a.   = list(atom = dplyr::filter(pseudohelixList[[1]]$atom, resno == 282, chain == "B", alt == "A")$eleno),
  tyr168oh_a = bio3d::atom.select(pseudohelixList[[1]], resno = 168, chain = "A", elety = "OH"),
  oxy_a      = bio3d::atom.select(pseudohelixList[[1]], resno = 232, chain = "A", elety = "O2"),
  co2_a      = bio3d::atom.select(pseudohelixList[[1]], resno = 228, chain = "B", elety = "C"),
  
  cu_b       = bio3d::atom.select(pseudohelixList[[1]], resno = 229, chain = "B"),
  nterm_b    = bio3d::atom.select(pseudohelixList[[1]], resno = 1, chain = "B", elety = "N"),
  his1nd_b   = bio3d::atom.select(pseudohelixList[[1]], resno = 1, chain = "B", elety = "ND1"),
  his1ne_b   =  bio3d::atom.select(pseudohelixList[[1]], resno = 1, chain = "B", elety = "NE2"),
  his1cg_b   = bio3d::atom.select(pseudohelixList[[1]], resno = 1, chain = "B", elety = "CG"),
  his84nd_b  = bio3d::atom.select(pseudohelixList[[1]], resno = 84, chain = "B", elety = "ND1"),
  his84ne_b  = bio3d::atom.select(pseudohelixList[[1]], resno = 84, chain = "B", elety = "NE2"),
  his84cg_b  = bio3d::atom.select(pseudohelixList[[1]], resno = 84, chain = "B", elety = "CG"),
  h2oax_b    = list(atom = dplyr::filter(pseudohelixList[[1]]$atom, resno == 288, chain == "B", alt == "A")$eleno),
  h2oeq_b    = list(atom = dplyr::filter(pseudohelixList[[1]]$atom, resno == 288, chain == "A", alt == "A")$eleno),
  tyr168oh_b = bio3d::atom.select(pseudohelixList[[1]], resno = 168, chain = "B", elety = "OH"),
  oxy_b      = bio3d::atom.select(pseudohelixList[[1]], resno = 231, chain = "A", elety = "O2"),
  glu30_b    = bio3d::atom.select(pseudohelixList[[1]], resno = 30, chain = "B", elety = "CD", alt = "A")
)

#Extract the atom element of the PDB files into new lists
pseudohelixAtoms <- lapply(pseudohelixList, function(pdb) pdb$atom) #Extract atoms from pseudohelices
wedgeAtoms <- lapply(wedgeList, function(pdb) pdb$atom) #Extract atoms from wedges

#Prep vectors containing all ddwd values
pseudohelixDose <- readr::read_csv("8-structure_analysis/input/samples.csv")$ddwd
wedgeDose <- readr::read_csv("8-structure_analysis/input/ddwds.csv") |> 
  filter(dataset_type == "wedge", dose_type == "ddwd") |> 
  pull("dose")

# Global Functions and Objects --------------------------------------------

source("global_functions.R")

# make blank lists to organize visualizations
multipleRegressions <- list()
regressionSummaries <- list()
ggplots <- list()

# Data Analysis -----------------------------------------------------------

source("8-structure_analysis/occupancy_analysis.R")
source("8-structure_analysis/distance_analysis.R")
source("8-structure_analysis/angle_analysis.R")
source("8-structure_analysis/structure_stats.R")
source("8-structure_analysis/pca_clustering.R")

# Data Visualization ------------------------------------------------------

source("8-structure_analysis/table_cleanup.R")
source("8-structure_analysis/heterogeneity_analysis.R")
source("8-structure_analysis/plot_generation.R")

# Data write-out ----------------------------------------------------------

#Save all the plots
save_plots <- function(parameter, outerheight = 10, outerwidth = 16, innerheight = 5, innerwidth = 8) {
  for(theme in names(ggplots)) {
    for(i in names(ggplots[[theme]][[parameter]])) {
      if("ggplot" %in% class(ggplots[[theme]][[parameter]][[i]])) {
        this_parameter <- parameter
        dataset_type <- i
        if(dataset_type == "broken_stick") {outerheight <- 16} else {outerheight <- 10}
        if(dataset_type == parameter) {this_parameter <- ""}
        ggplot2::ggsave(
          filename = stringr::str_glue("8-structure_analysis/output/plots/{theme}/{this_parameter}_{dataset_type}.svg"), 
          plot = ggplots[[theme]][[parameter]][[i]], 
          height = outerheight, 
          width = outerwidth,
          units = "cm",
          create.dir = TRUE
        )
      }
    }
  }
}

save_plots("pca", outerwidth = 16)
save_plots("Occupancies")
save_plots("Distances")
save_plots("Angles", outerheight = 12)
save_plots("CrystalStats", outerwidth = 16, outerheight = 16)
save_plots("RMSDs", outerheight = 16)

#Save all the tables
dir.create("8-structure_analysis/output/tables", showWarnings = FALSE)
readr::write_csv(longData, "8-structure_analysis/output/tables/regression_summary_pseudohelices.csv")
readr::write_csv(occ_table, "8-structure_analysis/output/tables/all_occupancies.csv")
readr::write_csv(angle_table, "8-structure_analysis/output/tables/all_angles.csv")
readr::write_csv(distance_table, "8-structure_analysis/output/tables/all_distances.csv")
