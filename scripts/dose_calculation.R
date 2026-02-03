# Data import -------------------------------------------------------------

raddose_output <- lapply(1:38, function(i) {
  fwds <- read.csv(paste0('input/raddose/fwd/wedge_', i, '.csv'), header = TRUE) %>% 
    dplyr::rename(fwd = DWD, angle = DWD.Angle)
  ddwds <- read.csv(paste0('input/raddose/ddwd/wedge_', i, '.csv'), header = TRUE) %>% 
    dplyr::rename(ddwd = DWD, angle = DWD.Angle)
  all_doses <- bind_cols(select(fwds, angle, fwd), select(ddwds, ddwd))
  return(all_doses)
}) # import all wedge fwds and ddwds

# Dose Calculations -------------------------------------------------------

wedge_doses <- tibble::tibble(
  dataset_type = rep("wedge", 36),
  dataset_number = 1:36,
  start_angle = seq(0, 175, 5)
)

wedge_doses$fwd <- sapply(2:37, function(wedge_number) {
  wedge <- raddose_output[[wedge_number]]
  return(mean(wedge$fwd))
}) # Wedges 2-37 were used, so calculate the average DDWD for all these
wedge_doses$ddwd <- sapply(2:37, function(wedge_number) {
  wedge <- raddose_output[[wedge_number]]
  return(mean(wedge$ddwd))
}) # Wedges 2-37 were used, so calculate the average DDWD for all these

pseudohelix_doses <- tibble::tibble(
  dataset_type = rep("pseudohelix", 176),
  dataset_number = 1:176,
  start_angle = 0:175
)

# calculate the respective subwedge average doses and average these across wedges
# 2-37 to calculate the average DWD for a pseudhelix

pseudohelix_doses$fwd <- sapply(pseudohelix_doses$start_angle, function(start_angle) {
  subwedge_averages <- sapply(2:37, function(subwedge_number) {
    wedge <- raddose_output[[subwedge_number]]
    subwedge_start <- start_angle + ((subwedge_number - 1) * 5) # Start angle depends on the subwedge you are calculating
    subwedge_average <- dplyr::filter(wedge, angle >= subwedge_start & angle <= subwedge_start + 4)$fwd %>% # Extract the DWD column for angles within a target subwedge
      mean() # Average the extracted DWDs
    return(subwedge_average)
  }) # Creates a vector of subwedge average DWDs
  pseudohelix_average <- mean(subwedge_averages) # Average the subwedge average DWDs to generate a single pseudohelix average DWD
  return(pseudohelix_average)
})
pseudohelix_doses$ddwd <- sapply(pseudohelix_doses$start_angle, function(start_angle) {
  subwedge_averages <- sapply(2:37, function(subwedge_number) {
    wedge <- raddose_output[[subwedge_number]]
    subwedge_start <- start_angle + ((subwedge_number - 1) * 5) # Start angle depends on the subwedge you are calculating
    subwedge_average <- dplyr::filter(wedge, angle >= subwedge_start & angle <= subwedge_start + 4)$ddwd %>% # Extract the DWD column for angles within a target subwedge
      mean() # Average the extracted DWDs
    return(subwedge_average)
  }) # Creates a vector of subwedge average DWDs
  pseudohelix_average <- mean(subwedge_averages) # Average the subwedge average DWDs to generate a single pseudohelix average DWD
  return(pseudohelix_average)
})

doses <- dplyr::bind_rows(pseudohelix_doses, wedge_doses) %>% 
  dplyr::mutate(ddwd = round(ddwd, 2)
