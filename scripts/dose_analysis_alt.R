# Data import -------------------------------------------------------------

wedgeDWDs <- lapply(1:38, function(i) {
  read.csv(paste0('input/raddose/wedge_', i, '.csv'), header = TRUE) %>% 
    dplyr::mutate(number = paste(!!i), datatype = 'Wedge')
}) # Import all wedge DWD traces for the wedges (for wedge and pseudohelix DWD calculations)

# DWD Analysis ------------------------------------------------------------

dwds <- tibble::tibble(
  datasetNumber = rep(1:36)
) # Initialize a data frame to put all the doses in

pseudohelix_start_angles <- c(0,1,2,3,4,5,6,7,8,9,11,13,16,20,26,35,46,57,68,79,89,100,111,122,133,144,153,160,165,167,169,170,171,172,173,174)
wedge_start_angles <- seq(0,179,5)

dwds$wedges <- sapply(2:37, function(wedgeNumber) {
  wedge <- wedgeDWDs[[wedgeNumber]]
  mean(wedge$DWD)
}) # Wedges 2-37 were used, so calculate the average DWD for all these

# Pseudohelix average DWD calculations calculate the average DWD for 5 frames
# worth of each wedge (subwedge). Calculate the respective subwedge average DWDs
# and average these across wedges 2-37 to calculate the average DWD for a pseudhelix
dwds$pseudohelices <- sapply(pseudohelix_start_angles, function(start_angle) {
  subwedgeAverages <- sapply(2:37, function(subwedgeNumber) {
    wedge <- wedgeDWDs[[subwedgeNumber]]
    startAngle <- start_angle + ((subwedgeNumber - 1) * 5) # Start angle depends on the subwedge you are calculating
    subwedgeAverage <- dplyr::filter(wedge, DWD.Angle >= startAngle & DWD.Angle <= startAngle + 4)$DWD %>% # Extract the DWD column for angles within a target subwedge (5 frames = 4° angular range) 
      mean() # Average the extracted DWDs
    return(subwedgeAverage)
  }) # Creates a vector of subwedge average DWDs
  pseudohelixAverage <- mean(subwedgeAverages) # Average the subwedge average DWDs to generate a single pseudohelix average DWD
  return(pseudohelixAverage)
})

# Reshape the dwd tibble so the wedge and pseudohelix dwds are stacked and labeled. This allows for proper ggplotting
dwds <- tibble::tibble(
  datasetNumber = rep(1:36, 2),
  start_angle = c(wedge_start_angles, pseudohelix_start_angles),
  dwd_MGy = c(dwds$wedges, dwds$pseudohelices),
  datasetType = c(rep("Wedges", 36), rep("Pseudohelices", 36))
)

test <- filter(dwds, datasetType == "Pseudohelices")
