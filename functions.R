# this file defines functions used in multiple other scripts

# function to calculate doses from raddose
calculate_dose <- function(type) {
  if (!type %in% c("fwd", "ddwd")) {stop("invalid dose type provided")}
  
  raddose_output <- lapply(1:38, function(i) {
    read.csv(stringr::str_glue("input/raddose/{type}/wedge_{i}.csv"), header = TRUE) |> 
      dplyr::rename(angle = DWD.Angle)
  })
  
  # average all fwd/ddwd values across a wedge
  wedge_doses <- tibble::tibble(
    dataset_type = rep("wedge", 36),
    dataset_number = 1:36,
    start_angle = seq(5, 180, 5),
    dose_type = rep(type, 36)
  )
  wedge_doses$dose <- sapply(2:37, function(wedge_number) {
    wedge <- raddose_output[[wedge_number]]
    return(mean(wedge$DWD))
  })
  
  # calculate the respective subwedge average doses and average these across wedges
  # 2-37 to calculate the average DWD for a pseudhelix
  pseudohelix_doses <- tibble::tibble(
    dataset_type = rep("pseudohelix", 176),
    dataset_number = 1:176,
    start_angle = 5:180,
    dose_type = rep(type, 176)
  )
  pseudohelix_doses$dose <- sapply(pseudohelix_doses$start_angle, function(start_angle) {
    subwedge_averages <- sapply(2:37, function(subwedge_number) {
      wedge <- raddose_output[[subwedge_number]]
      subwedge_start <- start_angle + ((subwedge_number - 2) * 5) # Start angle depends on the subwedge you are calculating
      subwedge_average <- dplyr::filter(wedge, angle >= subwedge_start & angle <= subwedge_start + 5)$DWD |>  # Extract the DWD column for angles within a target subwedge
        mean() # Average the extracted DWDs
      return(subwedge_average)
    }) # Creates a vector of subwedge average DWDs
    pseudohelix_average <- mean(subwedge_averages) # Average the subwedge average DWDs to generate a single pseudohelix average DWD
    return(pseudohelix_average)
  })
  
  doses <- dplyr::bind_rows(pseudohelix_doses, wedge_doses) |> 
    dplyr::mutate(dose = round(dose, 2))
  
  return(doses)
}

# emtrends helper functions
emtrends.coefficient <- function(model, regressor) {
  c(
    emmeans::emtrends(model, "Molecule", regressor) %>%
      test() %>% .[, 2], 
    emmeans::emtrends(model, "Molecule", regressor) %>%
      pairs() %>% summary() %>% .[, 2]
  )
}
emtrends.pvalue <- function(model, regressor) {
  c(
    emmeans::emtrends(model, "Molecule", regressor) %>%
      test() %>% .[, 6], 
    emmeans::emtrends(model, "Molecule", regressor) %>%
      pairs() %>% summary() %>% .[, 6]
  )
}
emtrends.se <- function(model, regressor) {
  c(
    emmeans::emtrends(model, "Molecule", regressor) %>%
      summary() %>% .[, 3], 
    emmeans::emtrends(model, "Molecule", regressor) %>%
      pairs() %>% summary() %>% .[, 3]
  )
}
