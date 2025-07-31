
wedgeNumber <- 1:36

# Extract electron densities from phenix output
electronDensities <- list(
  Pseudohelices = lapply(1:36, function(x) {
    lines <- readLines(paste0('Input/Electron_Density_Quantification/Outputs/Pseudohelices/Pseudohelix', x, '.txt')) # Put all lines into a list
    maplines <- grep("Map value:", lines, value = TRUE) # Extract only lines with a map value listed
    mapValues <- as.numeric(sub('.*Map value: ', '', maplines)) # Extract the map values and convert to a number
  }) %>% 
    as.data.frame() %>% # Add the map values to a dataframe
    t(),
  Wedges = lapply(1:36, function(x) {
    lines <- readLines(paste0('Input/Electron_Density_Quantification/Outputs/Wedges/Wedge', x, '.txt'))
    maplines <- grep("Map value:", lines, value = TRUE)
    mapValues <- as.numeric(sub('.*Map value: ', '', maplines))
  }) %>% 
    as.data.frame() %>% 
    t()
)

# Define grid dimensions
xdim <- seq(-43, 13, by = 2) %>% length()
ydim <- seq(-27, 30, by = 2) %>% length()
zdim <- seq(-53, 30, by = 2) %>% length()

electronDensities$Pseudohelices <- data.frame(electronDensities$Pseudohelices) #Convert matrix output from t() back to a df
rownames(electronDensities$Pseudohelices) <- 1:36
electronDensities$Wedges <- data.frame(electronDensities$Wedges)
rownames(electronDensities$Wedges) <- 1:36

# Prepare an empty list to contain the electron density slope arrays
slopeGrids <- list(
  Pseudohelices = array(NA, dim = c(xdim, ydim, zdim)),
  Wedges = array(NA, dim = c(xdim, ydim, zdim))
)

for (voxelIndex in 1:ncol(electronDensities$Pseudohelices)) {
  #Extract a vector of densities for a given voxel in the quantification grid
  voxelDensity <- electronDensities$Pseudohelices[,voxelIndex]
  
  #Generate a linear regression of density vs dose for the extracted voxel
  fit <- lm(voxelDensity ~ pseudohelixDose)
  
  #Get the slope from the linear regression
  slope <- coef(fit)[2]
  
  #Convert voxel index to relevant slopeGrid indices
  ix <- ((voxelIndex - 1) %% xdim) + 1 #X coordinate is modulo (remainder) of voxel index divided by X dimension of grid
  iy <- (((voxelIndex - 1) %/% xdim) %% ydim) + 1 #Y coordinate is modulo of voxel index divided by Y dimension of grid after dividing by xdim (to slice across X plane)
  iz <- ((voxelIndex - 1) %/% (xdim * ydim)) + 1 #Z coordinate is the quotient of voxel index divided by x and and y dimensions. This slices across x and y planes
  
  #Store slopes in slopeGrid at calculated index (ix, iy, iz)
  slopeGrids$Pseudohelices[ix,iy,iz] <- slope
}

for (voxelIndex in 1:ncol(electronDensities$Wedges)) {
  #Extract a vector of densities for a given voxel in the quantification grid
  voxelDensity <- electronDensities$Wedges[,voxelIndex]
  
  #Generate a linear regression of density vs wedge number for the extracted voxel
  fit <- lm(voxelDensity ~ wedgeNumber)
  
  #Get the slope from the linear regression
  slope <- coef(fit)[2]
  
  #Convert voxel index to relevant slopeGrid indices
  ix <- ((voxelIndex - 1) %% xdim) + 1 #X coordinate is modulo (remainder) of voxel index divided by X dimension of grid
  iy <- (((voxelIndex - 1) %/% xdim) %% ydim) + 1 #Y coordinate is modulo of voxel index divided by Y dimension of grid after dividing by xdim (to slice across X plane)
  iz <- ((voxelIndex - 1) %/% (xdim * ydim)) + 1 #Z coordinate is the quotient of voxel index divided by x and and y dimensions. This slices across x and y planes
  
  #Store slopes in slopeGrid at calculated index (ix, iy, iz)
  slopeGrids$Wedges[ix,iy,iz] <- slope
}

#Save numpy outputs
npySave('Scripts/Python/PseudohelixSlopeGrid.npy', slopeGrids$Pseudohelices)
npySave('Scripts/Python/WedgeSlopeGrid.npy', slopeGrids$Wedges)
