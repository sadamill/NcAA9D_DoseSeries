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
pseudohelixOccupancies <- data.frame(pseudohelixOccupancies, pseudohelixDose) #Add column for doses
colnames(pseudohelixOccupancies) <- c(stringr::str_glue("atom_{1:nrow(pseudohelixAtoms[[1]])}"), "dose_MGy") #Change column names

#Combine the dataframes to provide a singular data frame to pull from
trimmedOccupancies <- tibble::tibble(
  dose_MGy = rep(pseudohelixOccupancies$dose_MGy, 2), 
  CO2 = c(
    pseudohelixOccupancies[,atoms$co2_a$atom],
    rep(NA, 36)
  ), 
  Ax = c(
    pseudohelixOccupancies[,atoms$h2oax_a$atom], 
    pseudohelixOccupancies[,atoms$h2oax_b$atom]
  ), 
  Eq = c(
    pseudohelixOccupancies[,atoms$h2oeq_a$atom], 
    pseudohelixOccupancies[,atoms$h2oeq_b$atom]
  ), 
  Oxy = c(
    pseudohelixOccupancies[,atoms$oxy_a$atom], 
    pseudohelixOccupancies[,atoms$oxy_b$atom]
  ), 
  Glu = c(
    rep(NA, 36),
    pseudohelixOccupancies[,atoms$glu30_b$atom]
  ), 
  Molecule = c(
    rep("A", 36), 
    rep("B", 36)
  )
)

#Make a data frame containing occupancies, but stacked so they can be properly faceted in ggplot
stackedOccupancies <- tibble::tibble(
  Dose = rep(trimmedOccupancies$dose_MGy, 5), 
  Occupancy = c(
    trimmedOccupancies$CO2, 
    trimmedOccupancies$Ax, 
    trimmedOccupancies$Eq, 
    trimmedOccupancies$Oxy, 
    trimmedOccupancies$Glu
  ), 
  Residue = c(
    rep("CO2", 72), 
    rep("Ax", 72), 
    rep("Eq", 72), 
    rep("Oxy", 72), 
    rep("Glu", 72)
  ), 
  Molecule = rep(trimmedOccupancies$Molecule, 5), 
) %>% 
  dplyr::mutate(Residue = factor(Residue, levels = c("Oxy", "Glu", "CO2", "Ax", "Eq")))

# Linear regression analysis ----------------------------------------------

#Prepare a list of multiple linear regression models for each atom of interest
multipleRegressions$Occupancies <- list(
  Oxy = gls(Occupancy ~ Dose * Molecule, data = stackedOccupancies, subset = Residue == "Oxy", weights = varIdent(form = ~ 1 | Molecule)), 
  Glu = gls(Occupancy ~ Dose, data = stackedOccupancies, subset = Molecule == "B" & Residue == "Glu"), 
  CO2 = gls(Occupancy ~ Dose, data = stackedOccupancies, subset = Molecule == "A" & Residue == "CO2"), 
  Ax = gls(Occupancy ~ Dose * Molecule, data = stackedOccupancies, subset = Residue == "Ax", weights = varIdent(form = ~ 1 | Molecule)), 
  Eq = gls(Occupancy ~ Dose * Molecule, data = stackedOccupancies, subset = Residue == "Eq", weights = varIdent(form = ~ 1 | Molecule))
)

#Prepare summary table of linear regression analysis
regressionSummaries$Occupancies <- tibble::tibble(
  Residue = c(
    rep("Oxy", 3), 
    rep("Glu", 3), 
    rep("CO2", 3), 
    rep("Ax", 3), 
    rep("Eq", 3)
  ), 
  Measurement = rep("Occupancies", 15),
  Estimate = rep(c("TrendA", "TrendB", "Contrast"), 5), 
  Coefficient = c(
    emtrends.coefficient(multipleRegressions$Occupancies$Oxy, "Dose"), 
    NA,
    summary(multipleRegressions$Occupancies$Glu)$tTable[2, 1], 
    NA,
    summary(multipleRegressions$Occupancies$CO2)$tTable[2, 1], 
    NA, 
    NA, 
    emtrends.coefficient(multipleRegressions$Occupancies$Ax, "Dose"), 
    emtrends.coefficient(multipleRegressions$Occupancies$Eq, "Dose")
  ), 
  StandardError = c(
    emtrends.se(multipleRegressions$Occupancies$Oxy, "Dose"), 
    NA,
    summary(multipleRegressions$Occupancies$Glu)$tTable[2, 2], 
    NA,
    summary(multipleRegressions$Occupancies$CO2)$tTable[2, 2], 
    NA, 
    NA, 
    emtrends.se(multipleRegressions$Occupancies$Ax, "Dose"), 
    emtrends.se(multipleRegressions$Occupancies$Eq, "Dose")
  ), 
  PValue = c(
    emtrends.pvalue(multipleRegressions$Occupancies$Oxy, "Dose"), 
    NA, 
    summary(multipleRegressions$Occupancies$Glu)$tTable[2, 4], 
    NA, 
    summary(multipleRegressions$Occupancies$CO2)$tTable[2, 4], 
    NA, 
    NA, 
    emtrends.pvalue(multipleRegressions$Occupancies$Ax, "Dose"), 
    emtrends.pvalue(multipleRegressions$Occupancies$Eq, "Dose")
  )
)

occupancy_scatter <- multipleRegressions$Occupancies |> 
  purrr::map(\(model) {
    data <- getData(model)
    
    design_matrix <- model.matrix(model, data = data)
    pred <- predict(model)
    vcov_matrix <- model$varBeta
    
    se <- sqrt(rowSums((design_matrix %*% vcov_matrix) * design_matrix))
    
    alpha <- 0.05
    z_value <- qnorm(1 - alpha / 2)
    
    pred <- bind_cols(
      data,
      y_hat = pred,
      lower = pred - se * z_value,
      upper = pred + se * z_value,
    )
    
    return(pred)
    
  }) |> purrr::list_rbind()
