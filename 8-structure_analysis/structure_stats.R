read_crystal_stats <- function(filepath, dataset_type) {
  
  regressor_name <- deparse(substitute(regressor))
  
  start_angles <- if (dataset_type == "Pseudohelices") {
    readr::read_csv("8-structure_analysis/input/samples.csv") |> pull(start_angle)
  } else if (dataset_type == "Wedges") {
    readr::read_csv("8-structure_analysis/input/ddwds.csv") |> 
      filter(dataset_type == "wedge", dose_type == "ddwd") |> pull(start_angle)
  } else {stop("Improper dataset_type input")}
  
  data <- readr::read_csv(filepath) |> 
    mutate(...1 = make.unique(...1)) |> 
    t() |> 
    as_tibble() |> 
    row_to_names(row_number = 1) |> 
    clean_names() |> 
    select(!where(\(x) all(is.na(x)))) |> 
    mutate(across(where(\(x) {str_detect(x, "\\(") |> all()}),
                  ~ {str_split(., pattern = "\\s\\(") |> 
                      map(\(x) {str_remove(x, "\\)") |>
                          setNames(c(str_glue("{cur_column()}_overall"), str_glue("{cur_column()}_highest_shell")))})}),
           unit_cell = { str_split(unit_cell, pattern = " ") |> 
               purrr::map(\(x) setNames(x, c("a", "b", "c", "alpha", "beta", "gamma"))) }) |> 
    unnest_wider(where(\(x) is.list(x))) |> 
    mutate(across(!resolution_range_overall:space_group, ~ as.numeric(.)),
           unit_cell_volume = { a * b * c * sqrt(1 + 2 * cos((alpha * pi) / 180) * cos((beta * pi) / 180) * cos((gamma * pi) / 180) - (cos((alpha * pi) / 180)^2) - (cos((beta * pi) / 180)^2) - (cos((gamma * pi) / 180)^2)) },
           start_angle = start_angles,
           dataset_type = dataset_type) |> 
    relocate(unit_cell_volume, .after = gamma) |> 
    pivot_longer(cols = !c(dataset_type, start_angle, resolution_range_overall:space_group),
                 names_to = "statistic") |> 
    mutate(statistic = as_factor(statistic))
  
  return(data)
}

crystal_stats <- bind_rows(
  read_crystal_stats(
    filepath = "8-structure_analysis/input/crystal_statistics/pseudohelix_statistics.csv",
    dataset_type = "Pseudohelices"
  ),
  read_crystal_stats(
    filepath = "8-structure_analysis/input/crystal_statistics/wedge_statistics.csv",
    dataset_type = "Wedges"
  )
)
