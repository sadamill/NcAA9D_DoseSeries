# this file defines functions used in multiple other scripts

# function to calculate doses from raddose
calculate_dose <- function(type) {
  if (!type %in% c("fwd", "ddwd")) {stop("invalid dose type provided")}
  
  directory <- switch(type,
                      "fwd" = "1-fwd_calculation/input/r_input",
                      "ddwd" = "3-ddwd_calculation/input/r_input")
  
  raddose_output <- lapply(1:38, function(i) {
    read.csv(stringr::str_glue("{directory}/wedge{i}.csv"), header = TRUE) |> 
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
    emmeans::emtrends(model, "Molecule", regressor, mode = "appx-satterthwaite") %>%
      test() %>% .[, 2], 
    emmeans::emtrends(model, "Molecule", regressor, mode = "appx-satterthwaite") %>%
      pairs() %>% summary() %>% .[, 2]
  )
}
emtrends.pvalue <- function(model, regressor) {
  c(
    emmeans::emtrends(model, "Molecule", regressor, mode = "appx-satterthwaite") %>%
      test() %>% .[, 6], 
    emmeans::emtrends(model, "Molecule", regressor, mode = "appx-satterthwaite") %>%
      pairs() %>% summary() %>% .[, 6]
  )
}
emtrends.se <- function(model, regressor) {
  c(
    emmeans::emtrends(model, "Molecule", regressor, mode = "appx-satterthwaite") %>%
      summary() %>% .[, 3], 
    emmeans::emtrends(model, "Molecule", regressor, mode = "appx-satterthwaite") %>%
      pairs() %>% summary() %>% .[, 3]
  )
}

# custom ggplot themes
ggtheme_dark <- function(...) {
  list(
    ggplot2::theme_dark(base_size = 8, base_family = "ArialMT"), 
    ggplot2::theme(
      #Overall elements
      rect = ggplot2::element_blank(), 
      text = ggplot2::element_text(color = "white"), 
      line = ggplot2::element_line(color = "black"), 
      
      #Legend positioning
      legend.position = "inside", 
      legend.position.inside = c(0.85, 0.25), 
      
      #Manual override of desired ggplot2::theme elements
      panel.background = ggplot2::element_rect(fill = "black", color = "white"), 
      legend.key = ggplot2::element_blank(), #panel.background automatically maps to legend.key; I want to override this
      legend.background = ggplot2::element_rect(fill = "black", color = "white"), 
      strip.background = ggplot2::element_rect(fill = "black", color = "white"),
      axis.text = ggplot2::element_text(color = 'gray'),
      legend.key.height = grid::unit(3, "mm")
    ),
    ggplot2::theme(...)
  )
}
ggtheme_light <- function(...) {
  list(
    ggplot2::theme_bw(base_size = 8, base_family = "ArialMT"), 
    ggplot2::theme(
      #Overall elements
      text = ggplot2::element_text(color = "black"), 
      
      #Legend positioning
      legend.position = "inside", 
      legend.position.inside = c(0.85, 0.25), 
      
      #Manual override of desired ggplot2::theme elements
      legend.background = ggplot2::element_rect(fill = "white", color = "black"), 
      strip.background = ggplot2::element_rect(fill = "white"), 
      plot.background = ggplot2::element_blank(),
      legend.key.height = grid::unit(3, "mm")
    ),
    ggplot2::theme(...)
  )
}
