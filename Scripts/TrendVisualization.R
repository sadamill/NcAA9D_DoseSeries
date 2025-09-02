# Data preparation --------------------------------------------------------

make.wide <- function(trend, datasetType) {
  regressionSummaries[[trend]][[datasetType]] %>%
    filter(Estimate %in% c("TrendA", "TrendB")) %>%
    pivot_wider(
      id_cols = Residue, 
      names_from = Estimate, 
      values_from = Coefficient
    ) %>%
    left_join(
      regressionSummaries[[trend]][[datasetType]] %>%
        filter(Estimate == "Contrast") %>%
        select(Residue, PValue), 
      by = "Residue"
    ) %>%
    mutate(
      Significance = ifelse(PValue <= 0.05, "*", " "), 
      PValue = round(PValue, 3), 
      MinTrend = pmin(TrendA, TrendB)
    ) %>%
    arrange(MinTrend) %>%
    mutate(Residue = factor(Residue, levels = Residue))
}

wideData <- list(
  Occupancies = list(
    Pseudohelices = make.wide("Occupancies", "Pseudohelices"), 
    Wedges = make.wide("Occupancies", "Wedges")
  ), 
  BFactors = list(
    Pseudohelices = make.wide("BFactors", "Pseudohelices"), 
    Wedges = make.wide("BFactors", "Wedges")
  ), 
  Distances = list(
    Pseudohelices = make.wide("Distances", "Pseudohelices"), 
    Wedges = make.wide("Distances", "Wedges")
  ), 
  Angles = list(
    Pseudohelices = make.wide("Angles", "Pseudohelices"), 
    Wedges = make.wide("Angles", "Wedges")
  )
)
