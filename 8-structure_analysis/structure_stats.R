read_crystal_stats <- function(filepath, dataset_type) {
  
  regressor_name <- deparse(substitute(regressor))
  
  start_angles <- if (dataset_type == "Pseudohelices") {
    readr::read_csv("8-structure_analysis/input/samples.csv") |> pull(start_angle)
  } else {
    readr::read_csv("8-structure_analysis/input/ddwds.csv") |> 
      filter(dataset_type == "wedge", dose_type == "ddwd") |> pull(start_angle)
  }
  
  data <- readr::read_csv(filepath) |> 
    mutate(...1 = make.unique(...1)) |> 
    t() |> 
    as_tibble() |> 
    row_to_names(row_number = 1) |> 
    clean_names() |> 
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ gsub("\\s*\\(.*\\)", "", .)),
                  unit_cell = { str_split(unit_cell, pattern = " ") |> 
                      purrr::map(\(x) setNames(x, c("a", "b", "c", "alpha", "beta", "gamma"))) },
                  resolution_range = { str_split(resolution_range, pattern = "  - ") |> 
                      purrr::map(\(x) setNames(x, c("resolution_high", "resolution_low"))) }) |> 
    tidyr::unnest_wider(c(unit_cell, resolution_range)) |> 
    select(!where(~ all(is.na(.))) & !space_group) |> 
    purrr::map(\(x) as.numeric(x)) |> 
    bind_cols() |> 
    mutate(unit_cell_volume = { a * b * c * sqrt(1 + 2 * cos((alpha * pi) / 180) * cos((beta * pi) / 180) * cos((gamma * pi) / 180) - (cos((alpha * pi) / 180)^2) - (cos((beta * pi) / 180)^2) - (cos((gamma * pi) / 180)^2)) },
           start_angle = start_angles,
           dataset_type = dataset_type) |> 
    pivot_longer(cols = !c(dataset_type, start_angle),
                 names_to = "statistic")
  
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
