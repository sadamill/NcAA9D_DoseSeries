# Data preparation --------------------------------------------------------
compile_widen <- function(datasetType) {
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
    select(!Coefficient_Contrast) %>% 
    mutate(
      PValue_TrendA = p.adjust(PValue_TrendA, method = "BH"),
      PValue_TrendB = p.adjust(PValue_TrendA, method = "BH"),
      PValue_Contrast = p.adjust(PValue_TrendA, method = "BH")
    ) %>% 
    mutate(
      PValue_TrendA = signif(PValue_TrendA, 2),
      PValue_TrendB = signif(PValue_TrendB, 2),
      PValue_Contrast = signif(PValue_Contrast, 2)
    )
  return(df_wider)
}

wideData <- list(
  Pseudohelices = compile_widen("Pseudohelices"),
  Wedges = compile_widen("Wedges")
)
