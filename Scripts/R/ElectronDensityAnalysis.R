wedgeNumber <- 1:36

electronDensities <- list(
  Pseudohelices = lapply(1:36, function(x) {
    lines <- readLines(paste0('../Electron_Density_Quantification/Outputs/Pseudohelices/Pseudohelix', x, '.txt'))
    maplines <- grep("Map value:", lines, value = TRUE)
    mapValues <- as.numeric(sub('.*Map value: ', '', maplines))
  }) %>% 
    as.data.frame() %>% 
    t(),
  Wedges = lapply(1:36, function(x) {
    lines <- readLines(paste0('../Electron_Density_Quantification/Outputs/Wedges/Wedge', x, '.txt'))
    maplines <- grep("Map value:", lines, value = TRUE)
    mapValues <- as.numeric(sub('.*Map value: ', '', maplines))
  }) %>% 
    as.data.frame() %>% 
    t()
)

electronDensities$Pseudohelices <- data.frame(electronDensities$Pseudohelices) #Convert matrix output from t() back to a df
rownames(electronDensities$Pseudohelices) <- 1:36

electronDensities$Wedges <- data.frame(electronDensities$Wedges) #Convert matrix output from t() back to a df
rownames(electronDensities$Wedges) <- 1:36

slopeGrids <- list(
  Pseudohelices = array(NA, dim = c(23,31,37)),
  Wedges = array(NA, dim = c(23,31,37))
)

for (voxelIndex in 1:ncol(electronDensities$Pseudohelices)) {
  #Extract a vector of densities for a given voxel in the quantification grid
  voxelDensity <- electronDensities$Pseudohelices[,voxelIndex]
  
  #Generate a linear regression of density vs dose for the extracted voxel
  fit <- lm(voxelDensity ~ pseudohelixDose)
  
  #Get the slope from the linear regression
  slope <- coef(fit)[2]
  
  #Convert voxel index to relevant slopeGrid indices
  ix <- ((voxelIndex - 1) %% 23) + 1 #X coordinate is modulo (remainder) of voxel index divided by 23 (X dimension of grid)
  iy <- (((voxelIndex - 1) %/% 23) %% 31) + 1 #Y coordinate is modulo of voxel index divided by 27 (Y dimension of grid) after dividing by 23 (to slice across X plane)
  iz <- ((voxelIndex - 1) %/% (23 * 31)) + 1 #Z coordinate is the quotient of voxel index divided by 23 and 27. This slices across x and y planes
  
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
  ix <- ((voxelIndex - 1) %% 23) + 1 #X coordinate is modulo (remainder) of voxel index divided by 23 (X dimension of grid)
  iy <- (((voxelIndex - 1) %/% 23) %% 31) + 1 #Y coordinate is modulo of voxel index divided by 27 (Y dimension of grid) after dividing by 23 (to slice across X plane)
  iz <- ((voxelIndex - 1) %/% (23 * 31)) + 1 #Z coordinate is the quotient of voxel index divided by 23 and 27. This slices across x and y planes
  
  #Store slopes in slopeGrid at calculated index (ix, iy, iz)
  slopeGrids$Wedges[ix,iy,iz] <- slope
}

#Save numpy outputs
npySave('Scripts/Python/PseudohelixSlopeGrid.npy', slopeGrids$Pseudohelices)
npySave('Scripts/Python/WedgeSlopeGrid.npy', slopeGrids$Wedges)
