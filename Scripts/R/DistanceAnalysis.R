#This script is meant to be run from Master_Script.R

#This script analyzes the variance of various distances as a function of X-ray dose
#It  generates a plot of dose vs. distance for important atom pairs, allowing to highlight differences
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
    data.frame(
      CuTyr = c(dis[1449, 3986], dis[3435, 3987]), #Extract value from two cells in pseudohelixDMatList, the first containing distance for subunit A, and the other for subunit B
      CuNterm = c(dis[1, 3986], dis[1943, 3987]), #Repeat for all the desired values
      CuHis1ND1 = c(dis[7, 3986], dis[1949, 3987]), 
      CuHis84NE2 = c(dis[760, 3986], dis[2694, 3987]), 
      CuEq = c(dis[3986, 4711], dis[3987, 4580]), 
      CuAx = c(dis[3986, 4903], dis[3987, 4715]), 
      CuO1 = c(dis[3986, 5099], dis[3987, 5101]), 
      CuO2 = c(dis[3986, 5100], dis[3987, 5102]), 
      O1Eq = c(dis[4711, 5099], dis[4580, 5101]), 
      O2Eq = c(dis[4711, 5100], dis[4580, 5102]), 
      O1His157NE2 = c(dis[760, 5099], dis[3330, 5101]), 
      O2His157NE2 = c(dis[760, 5100], dis[3330, 5102]), 
      O1GluOE1 = c(dis[2201, 5099], dis[272, 5101]), 
      O2GluOE1 = c(dis[2201, 5100], dis[273, 5101]), 
      O1GluOE2 = c(dis[2202, 5099], dis[272, 5102]), 
      O2GluOE2 = c(dis[2202, 5100], dis[273, 5102]), 
      Molecule = c("A", "B")
    )
  }) %>% 
    ldply() %>% #Data is generated as a list; ldply turns it into a data frame
    .[order(.$Molecule), ] %>% #Data frame is ordered by alternating molecule, this will order the data frame by subunit
    data.frame(pseudohelixDose, .), #Amend a column containing doses to the data frame
  Wedges = lapply(distanceMatrices$Wedges, function(dis) {
    data.frame(
      CuTyr = c(dis[1449, 3986], dis[3435, 3987]), #Extract value from two cells in wedgeDMatList, the first containing distance for subunit A, and the other for subunit B
      CuNterm = c(dis[1, 3986], dis[1943, 3987]), #Repeat for all the desired values
      CuHis1ND1 = c(dis[7, 3986], dis[1949, 3987]), 
      CuHis84NE2 = c(dis[760, 3986], dis[2694, 3987]), 
      CuEq = c(dis[3986, 4711], dis[3987, 4580]), 
      CuAx = c(dis[3986, 4903], dis[3987, 4715]), 
      CuO1 = c(dis[3986, 5099], dis[3987, 5101]), 
      CuO2 = c(dis[3986, 5100], dis[3987, 5102]), 
      O1Eq = c(dis[4711, 5099], dis[4580, 5101]), 
      O2Eq = c(dis[4711, 5100], dis[4580, 5102]), 
      O1His157NE2 = c(dis[760, 5099], dis[3330, 5101]), 
      O2His157NE2 = c(dis[760, 5100], dis[3330, 5102]), 
      O1GluOE1 = c(dis[2201, 5099], dis[272, 5101]), 
      O2GluOE1 = c(dis[2201, 5100], dis[273, 5101]), 
      O1GluOE2 = c(dis[2202, 5099], dis[272, 5102]), 
      O2GluOE2 = c(dis[2202, 5100], dis[273, 5102]), 
      Molecule = c("A", "B")
    )
  }) %>% 
    ldply() %>% #Data is generated as a list; ldply turns it into a data frame
    .[order(.$Molecule), ] %>% #Data frame is ordered by alternating molecule, this will order the data frame by subunit
    data.frame(1:36, .) #Amend a column containing doses to the data frame
)

colnames(trimmedDistances$Pseudohelices)[1] <- "dose_MGy"
colnames(trimmedDistances$Wedges)[1] <- "WedgeNumber"

stackedDistances <- list(
  Pseudohelices = list(
    Cu = tibble(
      Dose = rep(trimmedDistances$Pseudohelices$dose_MGy, 6), 
      CuAtomPair = c(
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
    Oxy = tibble(
      Dose = rep(trimmedDistances$Pseudohelices$dose_MGy, 8), 
      OAtomPair = c(
        rep("O1-Eq", 72), 
        rep("O2-Eq", 72), 
        rep("O1-His157NE2", 72), 
        rep("O2-His157NE2", 72), 
        rep("O1-GluOE1", 72), 
        rep("O2-GluOE1", 72), 
        rep("O1-GluOE2", 72), 
        rep("O2-GluOE2", 72)
      ), 
      Distance = c(
        trimmedDistances$Pseudohelices$O1Eq, 
        trimmedDistances$Pseudohelices$O2Eq, 
        trimmedDistances$Pseudohelices$O1His157NE2, 
        trimmedDistances$Pseudohelices$O2His157NE2, 
        trimmedDistances$Pseudohelices$O1GluOE1, 
        trimmedDistances$Pseudohelices$O2GluOE1, 
        trimmedDistances$Pseudohelices$O1GluOE2, 
        trimmedDistances$Pseudohelices$O2GluOE2
      ), 
      Molecule = rep(trimmedDistances$Pseudohelices$Molecule, 8)
    )
  ), 
  Wedges = list(
    Cu = tibble(
      WedgeNumber = rep(trimmedDistances$Wedges$WedgeNumber, 6), 
      CuAtomPair = c(
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
    ), 
    Oxy = tibble(
      WedgeNumber = rep(trimmedDistances$Wedges$WedgeNumber, 8), 
      OAtomPair = c(
        rep("O1-Eq", 72), 
        rep("O2-Eq", 72), 
        rep("O1-His157NE2", 72), 
        rep("O2-His157NE2", 72), 
        rep("O1-GluOE1", 72), 
        rep("O2-GluOE1", 72), 
        rep("O1-GluOE2", 72), 
        rep("O2-GluOE2", 72)
      ), 
      Distance = c(
        trimmedDistances$Wedges$O1Eq, 
        trimmedDistances$Wedges$O2Eq, 
        trimmedDistances$Wedges$O1His157NE2, 
        trimmedDistances$Wedges$O2His157NE2, 
        trimmedDistances$Wedges$O1GluOE1, 
        trimmedDistances$Wedges$O2GluOE1, 
        trimmedDistances$Wedges$O1GluOE2, 
        trimmedDistances$Wedges$O2GluOE2
      ), 
      Molecule = rep(trimmedDistances$Wedges$Molecule, 8)
    )
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
    CuAx = lm(CuAx ~ dose_MGy * Molecule, data = trimmedDistances$Pseudohelices), 
    O1Eq = lm(O1Eq ~ dose_MGy * Molecule, data = trimmedDistances$Pseudohelices), 
    O2Eq = lm(O2Eq ~ dose_MGy * Molecule, data = trimmedDistances$Pseudohelices), 
    O1His157NE2 = lm(O1His157NE2 ~ dose_MGy * Molecule, data = trimmedDistances$Pseudohelices), 
    O2His157NE2 = lm(O2His157NE2 ~ dose_MGy * Molecule, data = trimmedDistances$Pseudohelices), 
    O1GluOE1 = lm(O1GluOE1 ~ dose_MGy * Molecule, data = trimmedDistances$Pseudohelices), 
    O1GluOE2 = lm(O1GluOE2 ~ dose_MGy * Molecule, data = trimmedDistances$Pseudohelices), 
    O2GluOE1 = lm(O2GluOE1 ~ dose_MGy * Molecule, data = trimmedDistances$Pseudohelices), 
    O2GluOE2 = lm(O2GluOE2 ~ dose_MGy * Molecule, data = trimmedDistances$Pseudohelices)
  ), 
  Wedges = list(
    CuTyr = lm(CuTyr ~ WedgeNumber * Molecule, data = trimmedDistances$Wedges), 
    CuNterm = lm(CuNterm ~ WedgeNumber * Molecule, data = trimmedDistances$Wedges), 
    CuHis1ND1 = lm(CuHis1ND1 ~ WedgeNumber * Molecule, data = trimmedDistances$Wedges), 
    CuHis84NE2 = lm(CuHis84NE2 ~ WedgeNumber * Molecule, data = trimmedDistances$Wedges), 
    CuEq = lm(CuEq ~ WedgeNumber * Molecule, data = trimmedDistances$Wedges), 
    CuAx = lm(CuAx ~ WedgeNumber * Molecule, data = trimmedDistances$Wedges), 
    O1Eq = lm(O1Eq ~ WedgeNumber * Molecule, data = trimmedDistances$Wedges), 
    O2Eq = lm(O2Eq ~ WedgeNumber * Molecule, data = trimmedDistances$Wedges), 
    O1His157NE2 = lm(O1His157NE2 ~ WedgeNumber * Molecule, data = trimmedDistances$Wedges), 
    O2His157NE2 = lm(O2His157NE2 ~ WedgeNumber * Molecule, data = trimmedDistances$Wedges), 
    O1GluOE1 = lm(O1GluOE1 ~ WedgeNumber * Molecule, data = trimmedDistances$Wedges), 
    O1GluOE2 = lm(O1GluOE2 ~ WedgeNumber * Molecule, data = trimmedDistances$Wedges), 
    O2GluOE1 = lm(O2GluOE1 ~ WedgeNumber * Molecule, data = trimmedDistances$Wedges), 
    O2GluOE2 = lm(O2GluOE2 ~ WedgeNumber * Molecule, data = trimmedDistances$Wedges)
  )
)

regressionSummaries$Distances <- list(
  Pseudohelices = bind_rows(
    lapply(
      c("CuTyr", "CuNterm", "CuHis1ND1", "CuHis84NE2", "CuEq", "CuAx", "O1Eq", "O2Eq", "O1His157NE2", "O2His157NE2", "O1GluOE1", "O2GluOE1", "O1GluOE2", "O2GluOE2"), 
      function(atom) {
        tibble(
          Residue = atom, 
          Estimate = c("TrendA", "TrendB", "Contrast"), 
          Coefficient = emtrends.coefficient(multipleRegressions$Distances$Pseudohelices[[atom]], "dose_MGy"), 
          StandardError = emtrends.se(multipleRegressions$Distances$Pseudohelices[[atom]], "dose_MGy"), 
          ModelRSquared = summary(multipleRegressions$Distances$Pseudohelices[[atom]])$r.squared, 
          PValue = emtrends.pvalue(multipleRegressions$Distances$Pseudohelices[[atom]], "dose_MGy")
        ) %>% mutate(Significance = ifelse(
          PValue <= 0.001, "***", 
          ifelse(
            PValue <= 0.01, "**", 
            ifelse(
              PValue <= 0.05, "*", " "
            )
          )
        ))
      }
    )
  ), 
  Wedges = bind_rows(
    lapply(
      c("CuTyr", "CuNterm", "CuHis1ND1", "CuHis84NE2", "CuEq", "CuAx", "O1Eq", "O2Eq", "O1His157NE2", "O2His157NE2", "O1GluOE1", "O2GluOE1", "O1GluOE2", "O2GluOE2"), 
      function(atom) {
        tibble(
          Residue = atom, 
          Estimate = c("TrendA", "TrendB", "Contrast"), 
          Coefficient = emtrends.coefficient(multipleRegressions$Distances$Wedges[[atom]], "WedgeNumber"), 
          StandardError = emtrends.se(multipleRegressions$Distances$Wedges[[atom]], "WedgeNumber"), 
          ModelRSquared = summary(multipleRegressions$Distances$Wedges[[atom]])$r.squared, 
          PValue = emtrends.pvalue(multipleRegressions$Distances$Wedges[[atom]], "WedgeNumber")
        ) %>% mutate(Significance = ifelse(
          PValue <= 0.001, "***", 
          ifelse(
            PValue <= 0.01, "**", 
            ifelse(
              PValue <= 0.05, "*", " "
            )
          )
        ))
      }
    )
  )
)
