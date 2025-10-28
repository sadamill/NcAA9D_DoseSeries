# Packages ----------------------------------------------------------------

require(cowplot)
require(tidyverse)
require(bio3d)
require(geometry)
require(emmeans)
require(officedown)
require(officer)
require(ggh4x)
require(ggnewscale)
require(plotly)
require(webshot2)
require(htmlwidgets)
require(janitor)
require(systemfonts)

options(scipen = 7) #I don't like viewing things in scientific notation
mem.maxVSize(vsize = 100000) #Set max memory size to allow space for large distance matrices

# Raw Data Extraction/Preparation -----------------------------------------

#Extract PDB files
pseudohelixList <- lapply(1:36, function(x) {
  read.pdb(
    file = str_glue("Input/Files_From_Analysis_Cluster/Pseudohelix_{x}/Refine_2/Phenix_refine_001.pdb"), 
    rm.alt = FALSE
  )
}) #Extract pseudohelix PDBs
wedgeList <- lapply(1:36, function(x) {
  read.pdb(
    file = str_glue("Input/Files_From_Analysis_Cluster/Wedge_{x}/Refine_2/Phenix_refine_001.pdb"), 
    rm.alt = FALSE
  )
}) #Extract wedge PDBs

atoms <- list( # Make list of atom and xyz indices
  cu_a = atom.select(pseudohelixList[[1]], resno = 1, chain = "C"),
  nterm_a = atom.select(pseudohelixList[[1]], resno = 1, chain = "A", elety = "N"),
  his1nd_a = atom.select(pseudohelixList[[1]], resno = 1, chain = "A", elety = "ND1"),
  his1ne_a =  atom.select(pseudohelixList[[1]], resno = 1, chain = "A", elety = "NE2"),
  his1cg_a = atom.select(pseudohelixList[[1]], resno = 1, chain = "A", elety = "CG"),
  his84nd_a = atom.select(pseudohelixList[[1]], resno = 84, chain = "A", elety = "ND1"),
  his84ne_a = atom.select(pseudohelixList[[1]], resno = 84, chain = "A", elety = "NE2"),
  his84cg_a = atom.select(pseudohelixList[[1]], resno = 84, chain = "A", elety = "CG"),
  h2oax_a = list(atom = filter(pseudohelixList[[1]]$atom, resno == 841, alt == "A")$eleno), # atom.select isn't sensitive to altcons
  h2oeq_a = list(atom = filter(pseudohelixList[[1]]$atom, resno == 639, alt == "A")$eleno),
  tyr168oh_a = atom.select(pseudohelixList[[1]], resno = 168, chain = "A", elety = "OH"),
  oxy_a = atom.select(pseudohelixList[[1]], resno = 1, chain = "D", elety = "O2"),
  glu30_b = list(atom = filter(pseudohelixList[[1]]$atom, chain == "B", resno == 30, elety == "CD", alt == "A")$eleno),
  
  cu_b = atom.select(pseudohelixList[[1]], resno = 2, chain = "C"),
  nterm_b = atom.select(pseudohelixList[[1]], resno = 1, chain = "B", elety = "N"),
  his1nd_b = atom.select(pseudohelixList[[1]], resno = 1, chain = "B", elety = "ND1"),
  his1ne_b =  atom.select(pseudohelixList[[1]], resno = 1, chain = "B", elety = "NE2"),
  his1cg_b = atom.select(pseudohelixList[[1]], resno = 1, chain = "B", elety = "CG"),
  his84nd_b = atom.select(pseudohelixList[[1]], resno = 84, chain = "B", elety = "ND1"),
  his84ne_b = atom.select(pseudohelixList[[1]], resno = 84, chain = "B", elety = "NE2"),
  his84cg_b = atom.select(pseudohelixList[[1]], resno = 84, chain = "B", elety = "CG"),
  h2oax_b = list(atom = filter(pseudohelixList[[1]]$atom, resno == 642, alt == "A")$eleno),
  h2oeq_b = list(atom = filter(pseudohelixList[[1]]$atom, resno == 512, alt == "A")$eleno),
  tyr168oh_b = atom.select(pseudohelixList[[1]], resno = 168, chain = "B", elety = "OH"),
  oxy_b = atom.select(pseudohelixList[[1]], resno = 2, chain = "D", elety = "O2"),
  glu30_a = atom.select(pseudohelixList[[1]], resno = 30, chain = "A", elety = "CD", alt = "A"),
  co2_b = atom.select(pseudohelixList[[1]], resno = 243, chain = "B", elety = "C")
)

#Extract the atom element of the PDB files into new lists
pseudohelixAtoms <- lapply(pseudohelixList, function(pdb) pdb$atom) #Extract atoms from pseudohelices
wedgeAtoms <- lapply(wedgeList, function(pdb) pdb$atom) #Extract atoms from wedges

#Run dose analysis to prepare dose state vectors
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
      pairs() %>% summary() %>% .[, 6]
  )
} #Extract p-values from emtrends
emtrends.se <- function(model, regressor) {
  c(
    emtrends(model, "Molecule", regressor) %>%
      summary() %>% .[, 3], 
    emtrends(model, "Molecule", regressor) %>%
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

source("Scripts/TrendVisualization.R")
source("Scripts/CheckingForHeterogeneity.R")
source("Scripts/PlotGeneration.R")
source("Scripts/WedgeCVAnalysis.R")

# Data write-out ----------------------------------------------------------

#Save all the plots
save_plots <- function(parameter, outerheight = 10, outerwidth = 16, innerheight = 5, innerwidth = 8) {
  for(theme in names(ggplots)) {
    for(i in names(ggplots[[theme]][[parameter]])) {
      if("ggplot" %in% class(ggplots[[theme]][[parameter]][[i]])) {
        dataset_type = if(i == "Pseudohelices") {"Pseudohelix"} else if(i == "Wedges") {"Wedge"} else if(i == "DWDs") {"DWDs"} else if(i == "CrystalStats") {"CrystalStats"} else if(i == "RMSDs") {"RMSDs"}
        this_parameter = if(parameter %in% c("Angles", "Distances", "Occupancies")) {parameter} else {""}
        ggsave(
          filename = str_glue("Output/Plots/{theme}/{dataset_type}_{parameter}.svg"), 
          plot = ggplots[[theme]][[parameter]][[i]], 
          height = outerheight, 
          width = outerwidth,
          units = "cm"
        )
      } else if("list" %in% class(ggplots[[theme]][[parameter]][[i]])) {
        for(j in names(ggplots[[theme]][[parameter]][[i]])) {
          dataset_type = if(j == "Pseudohelices") {"Pseudohelix"} else if(j == "Wedges") {"Wedge"}
          ggsave(
            filename = str_glue("Output/Plots/{theme}/{dataset_type}_{parameter}{i}.svg"), 
            plot = ggplots[[theme]][[parameter]][[i]][[j]], 
            height = innerheight, 
            width = innerwidth,
            units = "cm"
          )
        }
      }
    }
  }
}

save_plots('Occupancies')
save_plots('Distances')
save_plots('Angles', outerheight = 12)
save_plots('CVs')
save_plots('Dose', outerwidth = 8.8, outerheight = 7)
save_plots('Stats', outerwidth = 16, outerheight = 16)
save_plots("Comparisons", outerheight = 8)

#Save all associated PDBs
write.pdb(pdb = OccupancyColoredPDB, file = "Output/ColoredPDBs/OccupancyColoredPDB.pdb")

#Save all the tables
write_csv(longData$Pseudohelices, "Output/Tables/RegressionSummary_Pseudohelices.csv")
write_csv(longData$Wedges, "Output/Tables/RegressionSummary_Wedges.csv")
write_csv(occ_table, "Output/Tables/all_occupancies.csv")
write_csv(angle_table, "Output/Tables/all_angles.csv")
write_csv(distance_table, "Output/Tables/all_distances.csv")

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
