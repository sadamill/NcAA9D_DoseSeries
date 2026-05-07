# Packages ----------------------------------------------------------------

# data import/cleaning/wrangling
library(bio3d)
library(janitor)
library(tidyverse)
library(fst)
library(vroom)

# linear models and estimated marginal trends
library(nlme)
library(emmeans)

# plotting
library(cowplot)
library(ggnewscale)
library(systemfonts)
library(ggforce)
library(plotly)
library(gganimate)

# pca/k-means clustering packages
library(factoextra)
library(PCDimension)
library(NbClust)

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
  nterm_a    = bio3d::atom.select(pseudohelixList[[1]], resno = 1,   chain = "A", elety = "N"),
  his1nd_a   = bio3d::atom.select(pseudohelixList[[1]], resno = 1,   chain = "A", elety = "ND1"),
  his1ne_a   = bio3d::atom.select(pseudohelixList[[1]], resno = 1,   chain = "A", elety = "NE2"),
  his1cg_a   = bio3d::atom.select(pseudohelixList[[1]], resno = 1,   chain = "A", elety = "CG"),
  his84nd_a  = bio3d::atom.select(pseudohelixList[[1]], resno = 84,  chain = "A", elety = "ND1"),
  his84ne_a  = bio3d::atom.select(pseudohelixList[[1]], resno = 84,  chain = "A", elety = "NE2"),
  his84cg_a  = bio3d::atom.select(pseudohelixList[[1]], resno = 84,  chain = "A", elety = "CG"),
  h2oax_a    = list(atom = dplyr::filter(pseudohelixList[[1]]$atom, resno == 416, chain == "A", alt == "A")$eleno), # atom.select isn't sensitive to altcons .-.
  h2oeq_a    = list(atom = dplyr::filter(pseudohelixList[[1]]$atom, resno == 282, chain == "B", alt == "A")$eleno),
  tyr168oh_a = bio3d::atom.select(pseudohelixList[[1]], resno = 168, chain = "A", elety = "OH"),
  oxy_a      = bio3d::atom.select(pseudohelixList[[1]], resno = 232, chain = "A", elety = "O2"),
  co2_a      = bio3d::atom.select(pseudohelixList[[1]], resno = 228, chain = "B", elety = "C"),
  
  cu_b       = bio3d::atom.select(pseudohelixList[[1]], resno = 229, chain = "B"),
  nterm_b    = bio3d::atom.select(pseudohelixList[[1]], resno = 1,   chain = "B", elety = "N"),
  his1nd_b   = bio3d::atom.select(pseudohelixList[[1]], resno = 1,   chain = "B", elety = "ND1"),
  his1ne_b   = bio3d::atom.select(pseudohelixList[[1]], resno = 1,   chain = "B", elety = "NE2"),
  his1cg_b   = bio3d::atom.select(pseudohelixList[[1]], resno = 1,   chain = "B", elety = "CG"),
  his84nd_b  = bio3d::atom.select(pseudohelixList[[1]], resno = 84,  chain = "B", elety = "ND1"),
  his84ne_b  = bio3d::atom.select(pseudohelixList[[1]], resno = 84,  chain = "B", elety = "NE2"),
  his84cg_b  = bio3d::atom.select(pseudohelixList[[1]], resno = 84,  chain = "B", elety = "CG"),
  h2oax_b    = list(atom = dplyr::filter(pseudohelixList[[1]]$atom, resno == 288, chain == "B", alt == "A")$eleno),
  h2oeq_b    = list(atom = dplyr::filter(pseudohelixList[[1]]$atom, resno == 288, chain == "A", alt == "A")$eleno),
  tyr168oh_b = bio3d::atom.select(pseudohelixList[[1]], resno = 168, chain = "B", elety = "OH"),
  oxy_b      = bio3d::atom.select(pseudohelixList[[1]], resno = 231, chain = "A", elety = "O2"),
  glu30_b    = bio3d::atom.select(pseudohelixList[[1]], resno = 30,  chain = "B", elety = "CD", alt = "A")
)

# Extract the atom element of the PDB files into new lists
pseudohelixAtoms <- lapply(pseudohelixList, function(pdb) pdb$atom) #Extract atoms from pseudohelices
wedgeAtoms <- lapply(wedgeList, function(pdb) pdb$atom) #Extract atoms from wedges

# Prep vectors containing all ddwd values
pseudohelixDose <- readr::read_csv("8-structure_analysis/input/samples.csv")$ddwd

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
source("8-structure_analysis/dose_distribution.R")
source("8-structure_analysis/publication_comparisons.R")

# Data Visualization ------------------------------------------------------

source("8-structure_analysis/table_cleanup.R")
source("8-structure_analysis/heterogeneity_analysis.R")
source("8-structure_analysis/plot_generation.R")

# Data write-out ----------------------------------------------------------

plot_dimensions <- tibble(
  measurement = c(rep("pca", 3),
                 rep("Occupancies", 2),
                 rep("Distances", 2),
                 rep("Angles", 2),
                 "CrystalStats",
                 "RMSDs",
                 "dose_distribution",
                 rep("publication_comparisons", 2)),
  plot = c("broken_stick", "wedge_clusters", "pseudohelix_clusters",
           "scatter", "trends",
           "scatter", "trends",
           "scatter", "trends",
           "CrystalStats",
           "RMSDs",
           "dose_distribution",
           "dose", "distance"),
  height = c(16, 5, 5,
             10, 5,
             10, 5,
             12, 5,
             16,
             16,
             16,
             8, 8),
  width = c(16, 8, 8,
            16, 8,
            16, 8,
            16, 8,
            16,
            16,
            12,
            12, 12)
)

#Save all the plots
purrr::imap(ggplots, \(theme_list, theme_name) {
  purrr::imap(theme_list, \(measurement_list, measurement_name) {
    purrr::imap(measurement_list, \(plot, plot_name) {
      
      height <- filter(plot_dimensions, measurement == measurement_name & plot == plot_name)$height
      width <- filter(plot_dimensions, measurement == measurement_name & plot == plot_name)$width
      
      ggsave(filename = stringr::str_glue("8-structure_analysis/output/plots/{theme_name}/{measurement_name}_{plot_name}.svg"),
             plot = plot,
             height = height,
             width = width,
             units = "cm",
             create.dir = TRUE)
      
    })
  })
})

#Save all the tables
readr::write_csv(longData, "8-structure_analysis/output/tables/regression_summary.csv")
readr::write_csv(occ_table, "8-structure_analysis/output/tables/all_occupancies.csv")
readr::write_csv(angle_table, "8-structure_analysis/output/tables/all_angles.csv")
readr::write_csv(distance_table, "8-structure_analysis/output/tables/all_distances.csv")
