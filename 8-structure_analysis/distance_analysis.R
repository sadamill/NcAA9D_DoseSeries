#This script is meant to be run from Master_Script.R

#This script analyzes the variance of various distances as a function of X-ray dose
#It generates a plot of dose vs. distance for important atom pairs, allowing to highlight differences
#between different molecules in the asymmetric unit

# Distance Calculations ---------------------------------------------------

stackedDistances <- list(
  Pseudohelices = purrr::map2(pseudohelixDose, pseudohelixList, \(dose, structure) {
    dmat <- structure$xyz |> bio3d::dm()
    
    # extract distances of interest
    cu_eq_a <- dmat[atoms$cu_a$atom, atoms$h2oeq_a$atom]
    cu_eq_b <- dmat[atoms$h2oeq_b$atom, atoms$cu_b$atom]
    cu_ax_a <- dmat[atoms$cu_a$atom, atoms$h2oax_a$atom]
    cu_ax_b <- dmat[atoms$cu_b$atom, atoms$h2oax_b$atom]
    cu_tyr_a <- dmat[atoms$tyr168oh_a$atom, atoms$cu_a$atom]
    cu_tyr_b <- dmat[atoms$tyr168oh_b$atom, atoms$cu_b$atom]
    cu_nterm_a <- dmat[atoms$nterm_a$atom, atoms$cu_a$atom]
    cu_nterm_b <- dmat[atoms$nterm_b$atom, atoms$cu_b$atom]
    cu_his1_nd1_a <- dmat[atoms$his1nd_a$atom, atoms$cu_a$atom]
    cu_his1_nd1_b <- dmat[atoms$his1nd_b$atom, atoms$cu_b$atom]
    cu_his84_ne2_a <- dmat[atoms$his84ne_a$atom, atoms$cu_a$atom]
    cu_his84_ne2_b <- dmat[atoms$his84ne_b$atom, atoms$cu_b$atom]
    
    # make a tibble to contain target values
    distances <- tibble::tibble(
      Dose     = rep(dose, 12),
      AtomPair = c(rep("Cu-Tyr", 2),
                   rep("Cu-NTerm", 2), 
                   rep("Cu-His1ND1", 2), 
                   rep("Cu-His84NE2", 2), 
                   rep("Cu-Eq", 2), 
                   rep("Cu-Ax", 2)),
      Molecule = factor(rep(c("A", "B"), 6),
                        levels = c("A", "B")),
      Distance = c(cu_eq_a, cu_eq_b,
                   cu_ax_a, cu_ax_b,
                   cu_tyr_a, cu_tyr_b,
                   cu_nterm_a, cu_nterm_b,
                   cu_his1_nd1_a, cu_his1_nd1_b,
                   cu_his84_ne2_a, cu_his84_ne2_b)
    )
  }) |> purrr::list_rbind(),
  Wedges = purrr::map2(1:36, wedgeList, \(wedge_number, structure) {
    dmat <- structure$xyz |> bio3d::dm()
    
    # extract distances of interest
    cu_eq_a <- dmat[atoms$cu_a$atom, atoms$h2oeq_a$atom]
    cu_eq_b <- dmat[atoms$h2oeq_b$atom, atoms$cu_b$atom]
    cu_ax_a <- dmat[atoms$cu_a$atom, atoms$h2oax_a$atom]
    cu_ax_b <- dmat[atoms$cu_b$atom, atoms$h2oax_b$atom]
    cu_tyr_a <- dmat[atoms$tyr168oh_a$atom, atoms$cu_a$atom]
    cu_tyr_b <- dmat[atoms$tyr168oh_b$atom, atoms$cu_b$atom]
    cu_nterm_a <- dmat[atoms$nterm_a$atom, atoms$cu_a$atom]
    cu_nterm_b <- dmat[atoms$nterm_b$atom, atoms$cu_b$atom]
    cu_his1_nd1_a <- dmat[atoms$his1nd_a$atom, atoms$cu_a$atom]
    cu_his1_nd1_b <- dmat[atoms$his1nd_b$atom, atoms$cu_b$atom]
    cu_his84_ne2_a <- dmat[atoms$his84ne_a$atom, atoms$cu_a$atom]
    cu_his84_ne2_b <- dmat[atoms$his84ne_b$atom, atoms$cu_b$atom]
    
    # make a tibble to contain target values
    distances <- tibble::tibble(
      WedgeNumber = rep(wedge_number, 12),
      AtomPair     = c(rep("Cu-Tyr", 2),
                       rep("Cu-NTerm", 2), 
                       rep("Cu-His1ND1", 2), 
                       rep("Cu-His84NE2", 2), 
                       rep("Cu-Eq", 2), 
                       rep("Cu-Ax", 2)),
      Molecule     = factor(rep(c("A", "B"), 6),
                            levels = c("A", "B")),
      Distance     = c(cu_eq_a, cu_eq_b,
                       cu_ax_a, cu_ax_b,
                       cu_tyr_a, cu_tyr_b,
                       cu_nterm_a, cu_nterm_b,
                       cu_his1_nd1_a, cu_his1_nd1_b,
                       cu_his84_ne2_a, cu_his84_ne2_b)
    )
  }) |> purrr::list_rbind()
)

# Linear Regression Analysis ----------------------------------------------

#Prepare a list of multiple linear regression models for each atom pair of interest
multipleRegressions$Distances <- list(
  Pseudohelices = list(
    CuTyr = gls(Distance ~ Dose * Molecule, data = stackedDistances$Pseudohelices, subset = AtomPair == "Cu-Tyr", weights = varIdent(form = ~ 1 | Molecule)), 
    CuNterm = gls(Distance ~ Dose * Molecule, data = stackedDistances$Pseudohelices, subset = AtomPair == "Cu-NTerm", weights = varIdent(form = ~ 1 | Molecule)), 
    CuHis1ND1 = gls(Distance ~ Dose * Molecule, data = stackedDistances$Pseudohelices, subset = AtomPair == "Cu-His1ND1", weights = varIdent(form = ~ 1 | Molecule)), 
    CuHis84NE2 = gls(Distance ~ Dose * Molecule, data = stackedDistances$Pseudohelices, subset = AtomPair == "Cu-His84NE2", weights = varIdent(form = ~ 1 | Molecule)), 
    CuEq = gls(Distance ~ Dose * Molecule, data = stackedDistances$Pseudohelices, subset = AtomPair == "Cu-Eq", weights = varIdent(form = ~ 1 | Molecule)), 
    CuAx = gls(Distance ~ Dose * Molecule, data = stackedDistances$Pseudohelices, subset = AtomPair == "Cu-Ax", weights = varIdent(form = ~ 1 | Molecule))
  ), 
  Wedges = list(
    CuTyr = gls(Distance ~ WedgeNumber * Molecule, data = stackedDistances$Wedges, subset = AtomPair == "Cu-Tyr", weights = varIdent(form = ~ 1 | Molecule)), 
    CuNterm = gls(Distance ~ WedgeNumber * Molecule, data = stackedDistances$Wedges, subset = AtomPair == "Cu-NTerm", weights = varIdent(form = ~ 1 | Molecule)), 
    CuHis1ND1 = gls(Distance ~ WedgeNumber * Molecule, data = stackedDistances$Wedges, subset = AtomPair == "Cu-His1ND1", weights = varIdent(form = ~ 1 | Molecule)), 
    CuHis84NE2 = gls(Distance ~ WedgeNumber * Molecule, data = stackedDistances$Wedges, subset = AtomPair == "Cu-His84NE2", weights = varIdent(form = ~ 1 | Molecule)), 
    CuEq = gls(Distance ~ WedgeNumber * Molecule, data = stackedDistances$Wedges, subset = AtomPair == "Cu-Eq", weights = varIdent(form = ~ 1 | Molecule)), 
    CuAx = gls(Distance ~ WedgeNumber * Molecule, data = stackedDistances$Wedges, subset = AtomPair == "Cu-Ax", weights = varIdent(form = ~ 1 | Molecule))
  )
)

regressionSummaries$Distances <- list(
  Pseudohelices = dplyr::bind_rows(
    lapply(
      c("CuTyr", "CuNterm", "CuHis1ND1", "CuHis84NE2", "CuEq", "CuAx"), 
      function(atom) {
        tibble::tibble(
          Measurement = "Distances",
          Residue = atom, 
          Estimate = c("TrendA", "TrendB", "Contrast"), 
          Coefficient = emtrends.coefficient(multipleRegressions$Distances$Pseudohelices[[atom]], "Dose"), 
          StandardError = emtrends.se(multipleRegressions$Distances$Pseudohelices[[atom]], "Dose"), 
          PValue = emtrends.pvalue(multipleRegressions$Distances$Pseudohelices[[atom]], "Dose")
        )
      }
    )
  ), 
  Wedges = dplyr::bind_rows(
    lapply(
      c("CuTyr", "CuNterm", "CuHis1ND1", "CuHis84NE2", "CuEq", "CuAx"), 
      function(atom) {
        tibble::tibble(
          Measurement = "Distances",
          Residue = atom, 
          Estimate = c("TrendA", "TrendB", "Contrast"), 
          Coefficient = emtrends.coefficient(multipleRegressions$Distances$Wedges[[atom]], "WedgeNumber"), 
          StandardError = emtrends.se(multipleRegressions$Distances$Wedges[[atom]], "WedgeNumber"), 
          PValue = emtrends.pvalue(multipleRegressions$Distances$Wedges[[atom]], "WedgeNumber")
        )
      }
    )
  )
)

