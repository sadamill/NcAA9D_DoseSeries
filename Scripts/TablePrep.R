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

occ_table <- bind_rows(
  pivot_wider(
    stackedOccupancies$Pseudohelices,
    names_from = Residue,
    values_from = Occupancy
  ) %>% cbind(
    Dataset = rep(c(str_glue("Pseudohelix {1:36}")), 2),
    .
  ),
  
  pivot_wider(
    stackedOccupancies$Wedges,
    names_from = Residue,
    values_from = Occupancy
  ) %>% mutate(
    Dataset = rep(c(str_glue("Wedge {1:36}")), 2),
    Dose = rep(wedgeDose, 2)
  ) %>% 
    select(!WedgeNumber)
)

angle_table <- bind_rows(
  pivot_wider(
    stackedAngles$Pseudohelices,
    names_from = AngleID,
    values_from = Angle
  ) %>% cbind(
    Dataset = rep(c(str_glue("Pseudohelix {1:36}")), 2),
    .
  ),
  
  pivot_wider(
    stackedAngles$Wedges,
    names_from = AngleID,
    values_from = Angle
  ) %>% mutate(
    Dataset = rep(c(str_glue("Wedge {1:36}")), 2),
    Dose = rep(wedgeDose, 2)
  ) %>% 
    select(!WedgeNumber)
)

distance_table <- bind_rows(
  pivot_wider(
    stackedDistances$Pseudohelices,
    names_from = AtomPair,
    values_from = Distance
  ) %>% cbind(
    Dataset = rep(c(str_glue("Pseudohelix {1:36}")), 2),
    .
  ),
  
  pivot_wider(
    stackedDistances$Wedges,
    names_from = AtomPair,
    values_from = Distance
  ) %>% mutate(
    Dataset = rep(c(str_glue("Wedge {1:36}")), 2),
    Dose = rep(wedgeDose, 2)
  ) %>% 
    select(!WedgeNumber)
)

