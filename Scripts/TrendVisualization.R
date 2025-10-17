# Data preparation --------------------------------------------------------
compile_widen <- function(datasetType) {
  compiled <<- rbind(
    regressionSummaries$Occupancies[[datasetType]],
    regressionSummaries$Distances[[datasetType]],
    regressionSummaries$Angles[[datasetType]]
  )
  df_wider <<- compiled %>% 
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
  compiled <<- rbind(
    regressionSummaries$Occupancies[[datasetType]],
    regressionSummaries$Distances[[datasetType]],
    regressionSummaries$Angles[[datasetType]]
  )
  
  full_table <- compiled %>% mutate(
    Significance = ifelse(
      PValue <= 0.001, "***", 
      ifelse(
        PValue <= 0.01, "**", 
        ifelse(
          PValue <= 0.05, "*", " "
        )
      )
    )
  )
  
  return(full_table)
}

longData <- list(
  Pseudohelices = make_long("Pseudohelices"),
  Wedges = make_long("Wedges")
)
