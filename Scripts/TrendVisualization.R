# Data preparation --------------------------------------------------------
compile_regressions <- function(datasetType) {
  compiled <- rbind(
    regressionSummaries$Occupancies[[datasetType]],
    regressionSummaries$Distances[[datasetType]],
    regressionSummaries$Angles[[datasetType]]
  )
  df_wider <- compiled %>% 
    pivot_wider(
      id_cols = Residue, 
      names_from = Estimate, 
      values_from = c(Coefficient, PValue)
    ) %>% 
    left_join(
      compiled %>%
        filter(Estimate == "Contrast") %>%
        select(Residue, Measurement),
      by = "Residue"
    ) %>% 
    mutate(
      PValue_TrendA = p.adjust(PValue_TrendA, method = "BH"),
      PValue_TrendB = p.adjust(PValue_TrendB, method = "BH"),
      PValue_Contrast = p.adjust(PValue_Contrast, method = "BH")
    ) %>% arrange(
      pmax(Coefficient_TrendA, Coefficient_TrendB, na.rm = TRUE)
    ) %>% mutate(
      Residue = factor(Residue, levels = Residue)
    )
  return(df_wider)
}

wideData <- list(
  Pseudohelices = compile_widen("Pseudohelices"),
  Wedges = compile_widen("Wedges")
)

make_long <- function(datasetType) {
  pivot_longer(
    wideData[[datasetType]],
    cols = Coefficient_TrendA:PValue_Contrast,
    names_to = c("Type", "Trend"),
    names_sep = "_"
  ) %>% pivot_wider(
    id_cols = c(Trend, Residue, Measurement),
    names_from = Type,
    values_from = value
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

longData <- list(
  Pseudohelices = make_long("Pseudohelices"),
  Wedges = make_long("Wedges")
)
