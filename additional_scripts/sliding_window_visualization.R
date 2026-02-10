library(tidyverse)

calculate_dose <- function(type) {
  if (!type %in% c("fwd", "ddwd")) {stop("invalid dose type provided")}
  
  directory <- switch(type,
                      "fwd" = "1-fwd_calculation/input/r_input",
                      "ddwd" = "3-ddwd_calculation/input/r_input")
  
  raddose_output <- lapply(1:38, function(i) {
    read.csv(stringr::str_glue("{directory}/wedge{i}.csv"), header = TRUE) |> 
      dplyr::rename(angle = DWD.Angle)
  })
  
  # calculate the respective subwedge average doses and average these across wedges
  # 2-37 to calculate the average DWD for a pseudhelix
  subwedge_matrix <- sapply(5:184, function(start_angle) {
    subwedge_averages <- lapply(2:37, function(subwedge_number) {
      wedge <- raddose_output[[subwedge_number]]
      subwedge_start <- start_angle + ((subwedge_number - 2) * 5) # Start angle depends on the subwedge you are calculating
      subwedge_average <- dplyr::filter(wedge, angle >= subwedge_start & angle <= subwedge_start + 1)$DWD |>  # Extract the DWD column for angles within a target subwedge
        mean() # Average the extracted DWDs
      return(c(subwedge_average, subwedge_start))
    }) # Creates a vector of subwedge average DWDs
    return(subwedge_averages)
  })
  
  subwedges <- subwedge_matrix |> t() |> 
    as_tibble() |> 
    setNames(1:36) |> 
    pivot_longer(1:36, names_to = "wedge") |> 
    unnest_wider(value, names_sep = "_") |> 
    setNames(c("wedge_number", "ddwd", "angle")) |> 
    mutate(
      wedge_number = as.integer(wedge_number),
      angle = as.integer(angle)
    ) |> 
    arrange(wedge_number)
    
  return(subwedges)
}

ggplot2::ggplot(ddwds, aes(x = wedge_number, y = delta_angle, fill = ddwd)) +
  geom_tile() +
  scale_fill_viridis_c(option = "inferno", name = "DDWD") +
  coord_cartesian(expand = FALSE) + 
  theme_classic() +
  labs(x = "Wedge Number", y = "Δφ Angle (°)")

ggplot2::ggplot(ddwds, aes(x = wedge_number, y = angle, fill = ddwd)) +
  geom_tile() +
  scale_fill_viridis_c(option = "inferno", name = "DDWD") +
  coord_cartesian(expand = FALSE) + 
  theme_classic() +
  labs(x = "Wedge Number", y = "φ Angle (°)")
