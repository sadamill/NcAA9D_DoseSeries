read_crystal_stats <- function(filepath, dataset_type) {
  
  regressor_name <- deparse(substitute(regressor))
  
  start_angles <- if (dataset_type == "Pseudohelices") {
    readr::read_csv("8-structure_analysis/input/samples.csv") |> pull(start_angle)
  } else {
    readr::read_csv("8-structure_analysis/input/ddwds.csv") |> 
      filter(dataset_type == "wedge", dose_type == "ddwd") |> pull(start_angle)
  }
  
  data <- readr::read_csv(filepath) %>%
    {.[[1]] <- make.unique(as.character(.[[1]])); .} %>% # Make all parameter names unique
    tibble::column_to_rownames(var = colnames(.)[1]) %>% # Convert parameter name to row names to prep for transposition
    t() %>% as_tibble() |>  # Transpose and convert back to tibble
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ gsub("\\s*\\(.*\\)", "", .))) %>% # Remove parenthetical stats
    dplyr::mutate(dplyr::across(dplyr::everything(), trimws)) %>% # Clean up trailing white spaces
    dplyr::mutate(dplyr::across(dplyr::where(~ all(grepl("^[0-9 .-]*$", .x))), readr::parse_number)) %>% # Coerce all-number vectors to numeric
    dplyr::mutate(
      start_angle = start_angles,
      dataset_type = dataset_type
    ) %>% 
    dplyr::select(start_angle, dataset_type, `Completeness (%)`, `Multiplicity`, `Mean I/sigma(I)`, `Wilson B-factor`, `Average B-factor`, `CC1/2`, `R-work`, `R-free`, `RMS(bonds)`, `RMS(angles)`) %>% 
    janitor::clean_names() %>% # Convert column names to snake case
    tidyr::pivot_longer(
      cols = completeness_percent:rms_angles,
      names_to = "statistic",
    ) %>% 
    dplyr::mutate(statistic = factor(
      statistic,
      levels = c('multiplicity', 'completeness_percent', 'mean_i_sigma_i', 'wilson_b_factor', 'average_b_factor', 'cc1_2', 'r_work', 'r_free', 'rms_bonds', 'rms_angles')
    ))
  
  return(data)
}

crystal_stats <- list()
crystal_stats$pseudohelices <- read_crystal_stats(
  filepath = "8-structure_analysis/input/crystal_statistics/pseudohelix_statistics.csv",
  dataset_type = "Pseudohelices"
)
crystal_stats$wedges <- read_crystal_stats(
  filepath = "8-structure_analysis/input/crystal_statistics/wedge_statistics.csv",
  dataset_type = "Wedges"
)

crystal_stats$combined <- rbind(crystal_stats$pseudohelices, crystal_stats$wedges)