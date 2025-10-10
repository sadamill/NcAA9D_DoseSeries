#This script is meant to be run from Master_Script.R

#This script analyzes the variance of various distances as a function of X-ray dose
#It generates a plot of dose vs. distance for important atom pairs, allowing to highlight differences
#between different molecules in the asymmetric unit

# Distance Calculations ---------------------------------------------------

#Generate lists containing all-atom distance matrices for all datasets
distanceMatrices <- list(
  Pseudohelices = lapply(pseudohelixList, function(pdb) {
    pdb$xyz %>% dm()
  }), 
  Wedges = lapply(wedgeList, function(pdb) {
    pdb$xyz %>% dm()
  })
)

trimmedDistances <- list(
  Pseudohelices = lapply(distanceMatrices$Pseudohelices, function(dis) {
    tibble(
      CuTyr = c(dis[1449, 3986], dis[3435, 3987]), #Extract value from two cells in pseudohelixDMatList, the first containing distance for subunit A, and the other for subunit B
      CuNterm = c(dis[1, 3986], dis[1943, 3987]), #Repeat for all the desired values
      CuHis1ND1 = c(dis[7, 3986], dis[1949, 3987]), 
      CuHis84NE2 = c(dis[760, 3986], dis[2694, 3987]), 
      CuEq = c(dis[3986, 4711], dis[3987, 4580]), 
      CuAx = c(dis[3986, 4903], dis[3987, 4715]), 
      CuO1 = c(dis[3986, 5099], dis[3987, 5101]), 
      CuO2 = c(dis[3986, 5100], dis[3987, 5102]), 
      Molecule = c("A", "B")
    )
  }) %>% 
    bind_rows() %>% #Data is generated as a list; bind_rows turns it into a data frame
    .[order(.$Molecule), ] %>% #Data frame is ordered by alternating molecule, this will order the data frame by subunit
    data.frame(pseudohelixDose, .), #Amend a column containing doses to the data frame
  Wedges = lapply(distanceMatrices$Wedges, function(dis) {
    tibble(
      CuTyr = c(dis[1449, 3986], dis[3435, 3987]), #Extract value from two cells in wedgeDMatList, the first containing distance for subunit A, and the other for subunit B
      CuNterm = c(dis[1, 3986], dis[1943, 3987]), #Repeat for all the desired values
      CuHis1ND1 = c(dis[7, 3986], dis[1949, 3987]), 
      CuHis84NE2 = c(dis[760, 3986], dis[2694, 3987]), 
      CuEq = c(dis[3986, 4711], dis[3987, 4580]), 
      CuAx = c(dis[3986, 4903], dis[3987, 4715]), 
      CuO1 = c(dis[3986, 5099], dis[3987, 5101]), 
      CuO2 = c(dis[3986, 5100], dis[3987, 5102]), 
      Molecule = c("A", "B")
    )
  }) %>% 
    bind_rows() %>% #Data is generated as a list; bind_rows turns it into a data frame
    .[order(.$Molecule), ] %>% #Data frame is ordered by alternating molecule, this will order the data frame by subunit
    data.frame(1:36, .) #Amend a column containing doses to the data frame
)

colnames(trimmedDistances$Pseudohelices)[1] <- "dose_MGy"
colnames(trimmedDistances$Wedges)[1] <- "WedgeNumber"

stackedDistances <- list(
  Pseudohelices = tibble(
    Dose = rep(trimmedDistances$Pseudohelices$dose_MGy, 6), 
    AtomPair = c(
      rep("Cu-Tyr", 72), 
      rep("Cu-NTerm", 72), 
      rep("Cu-His1ND1", 72), 
      rep("Cu-His84NE2", 72), 
      rep("Cu-Eq", 72), 
      rep("Cu-Ax", 72)
    ), 
    Distance = c(
      trimmedDistances$Pseudohelices$CuTyr, 
      trimmedDistances$Pseudohelices$CuNterm, 
      trimmedDistances$Pseudohelices$CuHis1ND1, 
      trimmedDistances$Pseudohelices$CuHis84NE2, 
      trimmedDistances$Pseudohelices$CuEq, 
      trimmedDistances$Pseudohelices$CuAx
    ), 
    Molecule = rep(trimmedDistances$Pseudohelices$Molecule, 6)
  ),
  Wedges = tibble(
    WedgeNumber = rep(trimmedDistances$Wedges$WedgeNumber, 6), 
    AtomPair = c(
      rep("Cu-Tyr", 72), 
      rep("Cu-NTerm", 72), 
      rep("Cu-His1ND1", 72), 
      rep("Cu-His84NE2", 72), 
      rep("Cu-Eq", 72), 
      rep("Cu-Ax", 72)
    ), 
    Distance = c(
      trimmedDistances$Wedges$CuTyr, 
      trimmedDistances$Wedges$CuNterm, 
      trimmedDistances$Wedges$CuHis1ND1, 
      trimmedDistances$Wedges$CuHis84NE2, 
      trimmedDistances$Wedges$CuEq, 
      trimmedDistances$Wedges$CuAx
    ), 
    Molecule = rep(trimmedDistances$Wedges$Molecule, 6)
  )
)

# Linear Regression Analysis ----------------------------------------------

#Prepare a list of multiple linear regression models for each atom pair of interest
multipleRegressions$Distances <- list(
  Pseudohelices = list(
    CuTyr = lm(CuTyr ~ dose_MGy * Molecule, data = trimmedDistances$Pseudohelices), 
    CuNterm = lm(CuNterm ~ dose_MGy * Molecule, data = trimmedDistances$Pseudohelices), 
    CuHis1ND1 = lm(CuHis1ND1 ~ dose_MGy * Molecule, data = trimmedDistances$Pseudohelices), 
    CuHis84NE2 = lm(CuHis84NE2 ~ dose_MGy * Molecule, data = trimmedDistances$Pseudohelices), 
    CuEq = lm(CuEq ~ dose_MGy * Molecule, data = trimmedDistances$Pseudohelices), 
    CuAx = lm(CuAx ~ dose_MGy * Molecule, data = trimmedDistances$Pseudohelices)
  ), 
  Wedges = list(
    CuTyr = lm(CuTyr ~ WedgeNumber * Molecule, data = trimmedDistances$Wedges), 
    CuNterm = lm(CuNterm ~ WedgeNumber * Molecule, data = trimmedDistances$Wedges), 
    CuHis1ND1 = lm(CuHis1ND1 ~ WedgeNumber * Molecule, data = trimmedDistances$Wedges), 
    CuHis84NE2 = lm(CuHis84NE2 ~ WedgeNumber * Molecule, data = trimmedDistances$Wedges), 
    CuEq = lm(CuEq ~ WedgeNumber * Molecule, data = trimmedDistances$Wedges), 
    CuAx = lm(CuAx ~ WedgeNumber * Molecule, data = trimmedDistances$Wedges)
  )
)

regressionSummaries$Distances <- list(
  Pseudohelices = bind_rows(
    lapply(
      c("CuTyr", "CuNterm", "CuHis1ND1", "CuHis84NE2", "CuEq", "CuAx"), 
      function(atom) {
        tibble(
          Measurement = "Distances",
          Residue = atom, 
          Estimate = c("TrendA", "TrendB", "Contrast"), 
          Coefficient = emtrends.coefficient(multipleRegressions$Distances$Pseudohelices[[atom]], "dose_MGy"), 
          StandardError = emtrends.se(multipleRegressions$Distances$Pseudohelices[[atom]], "dose_MGy"), 
          ModelRSquared = summary(multipleRegressions$Distances$Pseudohelices[[atom]])$r.squared, 
          PValue = emtrends.pvalue(multipleRegressions$Distances$Pseudohelices[[atom]], "dose_MGy")
        )
      }
    )
  ), 
  Wedges = bind_rows(
    lapply(
      c("CuTyr", "CuNterm", "CuHis1ND1", "CuHis84NE2", "CuEq", "CuAx"), 
      function(atom) {
        tibble(
          Measurement = "Distances",
          Residue = atom, 
          Estimate = c("TrendA", "TrendB", "Contrast"), 
          Coefficient = emtrends.coefficient(multipleRegressions$Distances$Wedges[[atom]], "WedgeNumber"), 
          StandardError = emtrends.se(multipleRegressions$Distances$Wedges[[atom]], "WedgeNumber"), 
          ModelRSquared = summary(multipleRegressions$Distances$Wedges[[atom]])$r.squared, 
          PValue = emtrends.pvalue(multipleRegressions$Distances$Wedges[[atom]], "WedgeNumber")
        )
      }
    )
  )
)
