read_crystal_stats <- function(filepath, regressor) {
  regressor_name <- deparse(substitute(regressor))
  
  data <- read_csv(filepath) %>%
    {.[[1]] <- make.unique(as.character(.[[1]])); .} %>% # Make all parameter names unique
    column_to_rownames(var = colnames(.)[1]) %>% # Convert parameter name to row names to prep for transposition
    t() %>% as_tibble(rownames = "dataset_number") %>% # Transpose and convert back to tibble
    mutate(across(everything(), ~ gsub("\\s*\\(.*\\)", "", .))) %>% # Remove parenthetical stats
    mutate(across(everything(), trimws)) %>% # Clean up trailing white spaces
    mutate(across(where(~ all(grepl("^[0-9 .-]*$", .x))), readr::parse_number)) %>% # Coerce all-number vectors to numeric
    mutate(!!regressor_name := regressor) %>% 
    select(!!regressor_name, `Completeness (%)`, `Multiplicity`, `Mean I/sigma(I)`, `Wilson B-factor`, `CC1/2`, `R-work`, `R-free`, `RMS(bonds)`, `RMS(angles)`) %>% 
    clean_names() %>% # Convert column names to snake case
    pivot_longer(
      cols = completeness_percent:rms_angles,
      names_to = "statistic",
    )
  return(data)
}

crystal_stats <- list()
crystal_stats$pseudohelices <- read_crystal_stats(
  filepath = "Input/Crystal_Statistics/Pseudohelix_Statistics.csv",
  regressor = pseudohelixDose
)
crystal_stats$wedges <- read_crystal_stats(
  filepath = "Input/Crystal_Statistics/Wedge_Statistics.csv",
  regressor = wedgeNumber
)