#This script is meant to be run from Master_Script.R

#This script analyzes the variance of occupancies as a function of X-ray dose
#It writes out a PDB file which replaces the B-factor column with the B-factor slope as a function of X-ray dose
#This allows the user to visualize potential impact sites in the structure by coloring the model by B-factor
#The script also generates a plot of dose vs. occupancy for important residues, allowing to highlight differences
#between different molecules in the asymmetric unit

# Occupancy analysis ------------------------------------------------------

#Make data frames containing the occupancy data for each atom across datasets
pseudohelixOccupancies <- lapply(pseudohelixAtoms, function(atoms) atoms$o) %>% #Pull occupancies from pseudohelices
  as.data.frame() %>% #Coerce the list into a data frame
  t() #Transpose the data frame to swap rows and columns
row.names(pseudohelixOccupancies) <- 1:36 #Set row names to be comprehensible
pseudohelixOccupancies <- data.frame(pseudohelixDose, pseudohelixOccupancies) #Add column for doses
colnames(pseudohelixOccupancies) <- c("dose_MGy", 1:5102) #Change column names

wedgeOccupancies <- lapply(wedgeAtoms, function(atoms) atoms$o) %>% #Pull occupancies from pseudohelices
  as.data.frame() %>% #Coerce the list into a data frame
  t() #Transpose the data frame to swap rows and columns
row.names(wedgeOccupancies) <- 1:36 #Set row names to be comprehensible
wedgeOccupancies <- data.frame(1:36, wedgeOccupancies) #Add column for doses
colnames(wedgeOccupancies) <- c("WedgeNumber", 1:5102) #Change column names

#Combine the dataframes to provide a singular data frame to pull from
trimmedOccupancies <- list(
  Pseudohelices = tibble(
    dose_MGy = rep(pseudohelixOccupancies$dose_MGy, 2), 
    CO2 = c(
      rep(NA, 36), 
      pseudohelixOccupancies$`3983`
    ), 
    Ax = c(
      pseudohelixOccupancies$`4903`, 
      pseudohelixOccupancies$`4715`
    ), 
    Eq = c(
      pseudohelixOccupancies$`4711`, 
      pseudohelixOccupancies$`4580`
    ), 
    Oxy = c(
      pseudohelixOccupancies$`5099`, 
      pseudohelixOccupancies$`5101`
    ), 
    Glu = c(
      pseudohelixOccupancies$`266`, 
      pseudohelixOccupancies$`2195`
    ), 
    Molecule = c(
      rep("A", 36), 
      rep("B", 36)
    )
  ), 
  Wedges = tibble(
    WedgeNumber = rep(1:36, 2), 
    CO2 = c(
      rep(NA, 36), 
      wedgeOccupancies$`3983`
    ), 
    Ax = c(
      wedgeOccupancies$`4903`, 
      wedgeOccupancies$`4715`
    ), 
    Eq = c(
      wedgeOccupancies$`4711`, 
      wedgeOccupancies$`4580`
    ), 
    Oxy = c(
      wedgeOccupancies$`5099`, 
      wedgeOccupancies$`5101`
    ), 
    Glu = c(
      wedgeOccupancies$`266`, 
      wedgeOccupancies$`2195`
    ), 
    Molecule = c(
      rep("A", 36), 
      rep("B", 36)
    )
  )
)

#Make a data frame containing occupancies, but stacked so they can be properly faceted in ggplot
stackedOccupancies <- list(
  Pseudohelices = tibble(
    Dose = rep(trimmedOccupancies$Pseudohelices$dose_MGy, 5), 
    Occupancy = c(
      trimmedOccupancies$Pseudohelices$CO2, 
      trimmedOccupancies$Pseudohelices$Ax, 
      trimmedOccupancies$Pseudohelices$Eq, 
      trimmedOccupancies$Pseudohelices$Oxy, 
      trimmedOccupancies$Pseudohelices$Glu
    ), 
    Residue = c(
      rep("CO2", 72), 
      rep("Ax", 72), 
      rep("Eq", 72), 
      rep("Oxy", 72), 
      rep("Glu", 72)
    ), 
    Molecule = rep(trimmedOccupancies$Pseudohelices$Molecule, 5), 
  ) %>% 
    mutate(Residue = factor(Residue, levels = c("Oxy", "Glu", "CO2", "Ax", "Eq"))), 
  Wedges =  tibble(
    WedgeNumber = rep(trimmedOccupancies[["Wedges"]]$WedgeNumber, 5), 
    Occupancy = c(
      trimmedOccupancies[["Wedges"]]$CO2, 
      trimmedOccupancies[["Wedges"]]$Ax, 
      trimmedOccupancies[["Wedges"]]$Eq, 
      trimmedOccupancies[["Wedges"]]$Oxy, 
      trimmedOccupancies[["Wedges"]]$Glu
    ), 
    Residue = c(
      rep("CO2", 72), 
      rep("Ax", 72), 
      rep("Eq", 72), 
      rep("Oxy", 72), 
      rep("Glu", 72)
    ), 
    Molecule = rep(trimmedOccupancies[["Wedges"]]$Molecule, 5), 
  ) %>% 
    mutate(Residue = factor(Residue, levels = c("Oxy", "Glu", "CO2", "Ax", "Eq")))
)

#Make vector containing occupancy change of each atom
pseudohelixOccSlopes <- apply(pseudohelixOccupancies[, 2:(ncol(pseudohelixOccupancies))], 2, function(x) { #Exclude the first column (dose) and apply across columns (margin = 2)
  lm(x ~ dose_MGy, data = pseudohelixOccupancies)$coefficients[2]  # Extract slope (2nd coefficient)
})

#Clean the slope vectors
pseudohelixOccSlopes[abs(pseudohelixOccSlopes) < 10e-15] <- 0 #Turn any values below 10e-15 to 0
pseudohelixOccSlopes <- pseudohelixOccSlopes * 1000 #Scale the numbers by 1000

# Linear regression analysis ----------------------------------------------

#Prepare a list of multiple linear regression models for each atom of interest
multipleRegressions$Occupancies <- list(
  Pseudohelices = list(
    Oxy = lm(Oxy ~ dose_MGy * Molecule, data = trimmedOccupancies$Pseudohelices), 
    Glu = lm(Glu ~ dose_MGy * Molecule, data = trimmedOccupancies$Pseudohelices), 
    CO2 = lm(CO2 ~ dose_MGy, data = subset(trimmedOccupancies$Pseudohelices, Molecule == "B")), 
    Ax = lm(Ax ~ dose_MGy * Molecule, data = trimmedOccupancies$Pseudohelices), 
    Eq = lm(Eq ~ dose_MGy * Molecule, data = trimmedOccupancies$Pseudohelices)
  ), 
  Wedges = list(
    Oxy = lm(Oxy ~ WedgeNumber * Molecule, data = trimmedOccupancies$Wedges), 
    Glu = lm(Glu ~ WedgeNumber * Molecule, data = trimmedOccupancies$Wedges), 
    CO2 = lm(CO2 ~ WedgeNumber, data = subset(trimmedOccupancies$Wedges, Molecule == "B")), 
    Ax = lm(Ax ~ WedgeNumber * Molecule, data = trimmedOccupancies$Wedges), 
    Eq = lm(Eq ~ WedgeNumber * Molecule, data = trimmedOccupancies$Wedges)
  )
)

#Prepare summary table of linear regression analysis
regressionSummaries$Occupancies <- list(
  Pseudohelices = tibble(
    Residue = c(
      rep("Oxy", 3), 
      rep("Glu", 3), 
      rep("CO2", 3), 
      rep("Ax", 3), 
      rep("Eq", 3)
    ), 
    Estimate = rep(c("TrendA", "TrendB", "Contrast"), 5), 
    Coefficient = c(
      emtrends.coefficient(multipleRegressions$Occupancies$Pseudohelices$Oxy, "dose_MGy"), 
      emtrends.coefficient(multipleRegressions$Occupancies$Pseudohelices$Glu, "dose_MGy"), 
      NA, 
      summary(multipleRegressions$Occupancies$Pseudohelices$CO2)$coefficients[2, 2], 
      NA, 
      emtrends.coefficient(multipleRegressions$Occupancies$Pseudohelices$Ax, "dose_MGy"), 
      emtrends.coefficient(multipleRegressions$Occupancies$Pseudohelices$Eq, "dose_MGy")
    ), 
    StandardError = c(
      emtrends.se(multipleRegressions$Occupancies$Pseudohelices$Oxy, "dose_MGy"), 
      emtrends.se(multipleRegressions$Occupancies$Pseudohelices$Glu, "dose_MGy"), 
      NA, 
      summary(multipleRegressions$Occupancies$Pseudohelices$CO2)$coefficients[2, 4], 
      NA, 
      emtrends.se(multipleRegressions$Occupancies$Pseudohelices$Ax, "dose_MGy"), 
      emtrends.se(multipleRegressions$Occupancies$Pseudohelices$Eq, "dose_MGy")
    ), 
    ModelRSquared = c(
      rep(summary(multipleRegressions$Occupancies$Pseudohelices$Oxy)$r.squared, 3), 
      rep(summary(multipleRegressions$Occupancies$Pseudohelices$Glu)$r.squared, 3), 
      NA, 
      summary(multipleRegressions$Occupancies$Pseudohelices$CO2)$r.squared, 
      NA, 
      rep(summary(multipleRegressions$Occupancies$Pseudohelices$Ax)$r.squared, 3), 
      rep(summary(multipleRegressions$Occupancies$Pseudohelices$Eq)$r.squared, 3)
    ), 
    PValue = c(
      emtrends.pvalue(multipleRegressions$Occupancies$Pseudohelices$Oxy, "dose_MGy"), 
      emtrends.pvalue(multipleRegressions$Occupancies$Pseudohelices$Glu, "dose_MGy"), 
      NA, 
      summary(multipleRegressions$Occupancies$Pseudohelices$CO2)$coefficients[2, 4], 
      NA, 
      emtrends.pvalue(multipleRegressions$Occupancies$Pseudohelices$Ax, "dose_MGy"), 
      emtrends.pvalue(multipleRegressions$Occupancies$Pseudohelices$Eq, "dose_MGy")
    )
  ) %>% mutate(Significance = ifelse(
    PValue <= 0.001, "***", 
    ifelse(
      PValue <= 0.01, "**", 
      ifelse(
        PValue <= 0.05, "*", " "
      )
    )
  )), 
  Wedges = tibble(
    Residue = c(
      rep("Oxy", 3), 
      rep("Glu", 3), 
      rep("CO2", 3), 
      rep("Ax", 3), 
      rep("Eq", 3)
    ), 
    Estimate = rep(c("TrendA", "TrendB", "Contrast"), 5), 
    Coefficient = c(
      emtrends.coefficient(multipleRegressions$Occupancies$Wedges$Oxy, "WedgeNumber"), 
      emtrends.coefficient(multipleRegressions$Occupancies$Wedges$Glu, "WedgeNumber"), 
      NA, 
      summary(multipleRegressions$Occupancies$Wedges$CO2)$coefficients[2, 2], 
      NA, 
      emtrends.coefficient(multipleRegressions$Occupancies$Wedges$Ax, "WedgeNumber"), 
      emtrends.coefficient(multipleRegressions$Occupancies$Wedges$Eq, "WedgeNumber")
    ), 
    StandardError = c(
      emtrends.se(multipleRegressions$Occupancies$Wedges$Oxy, "WedgeNumber"), 
      emtrends.se(multipleRegressions$Occupancies$Wedges$Glu, "WedgeNumber"), 
      NA, 
      summary(multipleRegressions$Occupancies$Wedges$CO2)$coefficients[2, 4], 
      NA, 
      emtrends.se(multipleRegressions$Occupancies$Wedges$Ax, "WedgeNumber"), 
      emtrends.se(multipleRegressions$Occupancies$Wedges$Eq, "WedgeNumber")
    ), 
    ModelRSquared = c(
      rep(summary(multipleRegressions$Occupancies$Wedges$Oxy)$r.squared, 3), 
      rep(summary(multipleRegressions$Occupancies$Wedges$Glu)$r.squared, 3), 
      NA, 
      summary(multipleRegressions$Occupancies$Wedges$CO2)$r.squared, 
      NA, 
      rep(summary(multipleRegressions$Occupancies$Wedges$Ax)$r.squared, 3), 
      rep(summary(multipleRegressions$Occupancies$Wedges$Eq)$r.squared, 3)
    ), 
    PValue = c(
      emtrends.pvalue(multipleRegressions$Occupancies$Wedges$Oxy, "WedgeNumber"), 
      emtrends.pvalue(multipleRegressions$Occupancies$Wedges$Glu, "WedgeNumber"), 
      NA, 
      summary(multipleRegressions$Occupancies$Wedges$CO2)$coefficients[2, 4], 
      NA, 
      emtrends.pvalue(multipleRegressions$Occupancies$Wedges$Ax, "WedgeNumber"), 
      emtrends.pvalue(multipleRegressions$Occupancies$Wedges$Eq, "WedgeNumber")
    )
  ) %>% mutate(Significance = ifelse(
    PValue <= 0.001, "***", 
    ifelse(
      PValue <= 0.01, "**", 
      ifelse(
        PValue <= 0.05, "*", " "
      )
    )
  ))
)

#Write out a PDB file for B-factor coloring in ChimeraX
OccupancyColoredPDB <- pseudohelixList[[1]]
OccupancyColoredPDB$atom$b <- pseudohelixOccSlopes