# Data preparation --------------------------------------------------------
compile_widen <- function(datasetType) {
  compiled <<- rbind(
    regressionSummaries$Occupancies[[datasetType]],
    regressionSummaries$Distances[[datasetType]],
    regressionSummaries$Angles[[datasetType]]
  )
  df_wider <<- compiled %>% 
    tidyr::pivot_wider(
      id_cols = Residue, 
      names_from = Estimate, 
      values_from = c(Coefficient, PValue)
    ) %>% 
    dplyr::left_join(
      compiled %>%
        dplyr::filter(Estimate == "Contrast") %>%
        dplyr::select(Residue, Measurement),
      by = "Residue"
    ) %>% dplyr::arrange(
      pmax(Coefficient_TrendA, Coefficient_TrendB, na.rm = TRUE)
    ) %>% dplyr::mutate(
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
  
  full_table <- compiled %>% dplyr::mutate(
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

occ_table <- dplyr::bind_rows(
  tidyr::pivot_wider(
    stackedOccupancies$Pseudohelices,
    names_from = Residue,
    values_from = Occupancy
  ) %>% cbind(
    Dataset = rep(c(stringr::str_glue("Pseudohelix {1:36}")), 2),
    .
  ),
  
  tidyr::pivot_wider(
    stackedOccupancies$Wedges,
    names_from = Residue,
    values_from = Occupancy
  ) %>% mutate(
    Dataset = rep(c(stringr::str_glue("Wedge {1:36}")), 2),
    Dose = rep(wedgeDose, 2)
  ) %>% 
    dplyr::select(!WedgeNumber)
)

angle_table <- dplyr::bind_rows(
  tidyr::pivot_wider(
    stackedAngles$Pseudohelices,
    names_from = AngleID,
    values_from = Angle
  ) %>% cbind(
    Dataset = rep(c(stringr::str_glue("Pseudohelix {1:36}")), 2),
    .
  ),
  
  tidyr::pivot_wider(
    stackedAngles$Wedges,
    names_from = AngleID,
    values_from = Angle
  ) %>% dplyr::mutate(
    Dataset = rep(c(stringr::str_glue("Wedge {1:36}")), 2),
    Dose = rep(wedgeDose, 2)
  ) %>% 
    dplyr::select(!WedgeNumber)
)

distance_table <- dplyr::bind_rows(
  tidyr::pivot_wider(
    stackedDistances$Pseudohelices,
    names_from = AtomPair,
    values_from = Distance
  ) %>% cbind(
    Dataset = rep(c(stringr::str_glue("Pseudohelix {1:36}")), 2),
    .
  ),
  
  pivot_wider(
    stackedDistances$Wedges,
    names_from = AtomPair,
    values_from = Distance
  ) %>% dplyr::mutate(
    Dataset = rep(c(stringr::str_glue("Wedge {1:36}")), 2),
    Dose = rep(wedgeDose, 2)
  ) %>% 
    dplyr::select(!WedgeNumber)
)

