# Data preparation --------------------------------------------------------
compile_widen <- function() {
  
  compiled <- rbind(
    regressionSummaries$Occupancies,
    regressionSummaries$Distances,
    regressionSummaries$Angles
  )
  
  df_wider <- compiled %>% 
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

wideData <- compile_widen()

make_long <- function() {
  
  full_table <- rbind(
    regressionSummaries$Occupancies,
    regressionSummaries$Distances,
    regressionSummaries$Angles
  ) %>% 
    dplyr::mutate(
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

longData <- make_long()

occ_table <- dplyr::bind_rows(
  tidyr::pivot_wider(
    bind_cols(
      rep(stringr::str_glue("Pseudohelix {1:36}"), 10),
      stackedOccupancies
    ),
    names_from = Residue,
    values_from = Occupancy
  )
) |> rename(Dataset = ...1)

distance_table <- dplyr::bind_rows(
  tidyr::pivot_wider(
    bind_cols(
      stringr::str_glue("Pseudohelix {sort(rep(1:36, 12))}"),
      stackedDistances
    ),
    names_from = AtomPair,
    values_from = Distance
  )
) |> rename(Dataset = ...1)

angle_table <- dplyr::bind_rows(
  tidyr::pivot_wider(
    bind_cols(
      rep(stringr::str_glue("Pseudohelix {1:36}"), 14),
      stackedAngles
    ),
    names_from = AngleID,
    values_from = Angle
  )
) |> rename(Dataset = ...1)
