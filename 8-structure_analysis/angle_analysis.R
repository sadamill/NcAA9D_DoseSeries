#This script is meant to be run from Master_Script.R

#This script analyzes the variance of various angles as a function of X-ray dose
#It  generates a plot of dose vs. angle for important atom sets, allowing to highlight differences
#between different molecules in the asymmetric unit

# Function definitions ----------------------------------------------------

# Calculates a vector orthogonal to a plane abc, where a, b, and c are three xyz points on a plane
orthogonal.vector <- function(a, b, c) {
  c(a-b, a-c) %>% # Calculate two vectors within the plane abc; point a is the vertex
    matrix(nrow = 2, byrow = TRUE) %>% #Make a 2x3 matrix containing ijk components of two vectors inside the plane abc
    {c( #Calculate the orthogonal vector by finding the cross product of the two vectors inside the matrix
      .[1, 2]*.[2, 3] - .[1, 3]*.[2, 2], #I component
      -(.[1, 1]*.[2, 3] - .[1, 3]*.[2, 1]), #J component
      .[1, 1]*.[2, 2] - .[1, 2]*.[2, 1] #K component
    )}
}

#Finds angle between 2 vectors (in degrees), where v1 and v2 are vectors <ijk> to measure the angle between
vector.angle <- function(v1, v2) {
  c(v1, v2) %>% 
    matrix(nrow = 2, byrow = TRUE) %>% #Make a 3x2 matrix containing the two vectors, one in each row
    {(geometry::dot(.[1, ], .[2, ]))/(sqrt(sum(.[1, 1]^2, .[1, 2]^2, .[1, 3]^2))*sqrt(sum(.[2, 1]^2, .[2, 2]^2, .[2, 3]^2)))} %>% 
    acos() %>% 
    {.*180/pi}
}

# Angle calculations ------------------------------------------------------

#Calculate angles for all datasets based off XYZ indices calculated earlier
allAngles <- lapply(pseudohelixList, function(pdb) {
    
    tibble::tibble(
      T1 = c(
        bio3d::angle.xyz(c(pdb$xyz[atoms$nterm_a$xyz], pdb$xyz[atoms$cu_a$xyz], pdb$xyz[atoms$his1nd_a$xyz])), 
        bio3d::angle.xyz(c(pdb$xyz[atoms$nterm_b$xyz], pdb$xyz[atoms$cu_b$xyz], pdb$xyz[atoms$his1nd_b$xyz]))
      ), 
      T2 = c(
        bio3d::angle.xyz(c(pdb$xyz[atoms$nterm_a$xyz], pdb$xyz[atoms$cu_a$xyz], pdb$xyz[atoms$his84ne_a$xyz])), 
        bio3d::angle.xyz(c(pdb$xyz[atoms$nterm_b$xyz], pdb$xyz[atoms$cu_b$xyz], pdb$xyz[atoms$his84ne_b$xyz]))
      ), 
      T3 = c(
        bio3d::angle.xyz(c(pdb$xyz[atoms$his1nd_a$xyz], pdb$xyz[atoms$cu_a$xyz], pdb$xyz[atoms$his84ne_a$xyz])), 
        bio3d::angle.xyz(c(pdb$xyz[atoms$his1nd_b$xyz], pdb$xyz[atoms$cu_b$xyz], pdb$xyz[atoms$his84ne_b$xyz]))
      ), 
      TT = c(
        orthogonal.vector(pdb$xyz[atoms$cu_a$xyz], pdb$xyz[atoms$nterm_a$xyz], pdb$xyz[atoms$his1nd_a$xyz]) %>% #Find the vector orthogonal to the Cu-Nterm-Nd1 plane
          vector.angle(pdb$xyz[atoms$his84ne_a$xyz]-pdb$xyz[atoms$cu_a$xyz]) %>% #Find the angle between the orthogonal vector and the Ne2-Cu vector
          {90-.}, 
        orthogonal.vector(pdb$xyz[atoms$cu_b$xyz], pdb$xyz[atoms$nterm_b$xyz], pdb$xyz[atoms$his1nd_b$xyz]) %>% #Repeat for molecule B
          vector.angle(pdb$xyz[atoms$his84ne_b$xyz]-pdb$xyz[atoms$cu_b$xyz]) %>%
          {90-.}
      ), 
      THH = c(
        orthogonal.vector(pdb$xyz[atoms$his1cg_a$xyz], pdb$xyz[atoms$his1nd_a$xyz], pdb$xyz[atoms$his1ne_a$xyz]) %>% #Find vector orthogonal to His1 imidazole ring plane
          vector.angle(
            orthogonal.vector(pdb$xyz[atoms$his84cg_a$xyz], pdb$xyz[atoms$his84nd_a$xyz], pdb$xyz[atoms$his84ne_a$xyz])#Find angle between the two vectors orthogonal to the His1-His84 imidazole planes
          ), 
        orthogonal.vector(pdb$xyz[atoms$his1cg_b$xyz], pdb$xyz[atoms$his1nd_b$xyz], pdb$xyz[atoms$his1ne_b$xyz]) %>% #Repeat for molecule B
          vector.angle(
            orthogonal.vector(pdb$xyz[atoms$his84cg_b$xyz], pdb$xyz[atoms$his84nd_b$xyz], pdb$xyz[atoms$his84ne_b$xyz])
          )
      ), 
      TH1 = c(
        orthogonal.vector(pdb$xyz[atoms$his1cg_a$xyz], pdb$xyz[atoms$his1nd_a$xyz], pdb$xyz[atoms$his1ne_a$xyz]) %>% #Find the vector orthogonal to the His1 imidazole ring plane
          vector.angle(pdb$xyz[atoms$his1nd_a$xyz]-pdb$xyz[atoms$cu_a$xyz]) %>% #Find the angle between the orthogonal vector and the Nd1-Cu vector
          {.-90}, 
        orthogonal.vector(pdb$xyz[atoms$his1cg_b$xyz], pdb$xyz[atoms$his1nd_b$xyz], pdb$xyz[atoms$his1ne_b$xyz]) %>% #Find the vector orthogonal to the His1 imidazole ring plane
          vector.angle(pdb$xyz[atoms$his1nd_b$xyz]-pdb$xyz[atoms$cu_b$xyz]) %>% #Find the angle between the orthogonal vector and the Nd1-Cu vector
          {.-90}
      ), 
      THN = c(
        orthogonal.vector(pdb$xyz[atoms$his84cg_a$xyz], pdb$xyz[atoms$his84nd_a$xyz], pdb$xyz[atoms$his84ne_a$xyz]) %>% #Find the vector orthogonal to the His1 imidazole ring plane
          vector.angle(pdb$xyz[atoms$his84nd_a$xyz]-pdb$xyz[atoms$cu_a$xyz]) %>% #Find the angle between the orthogonal vector and the Nd1-Cu vector
          {.-90}, 
        orthogonal.vector(pdb$xyz[atoms$his84cg_b$xyz], pdb$xyz[atoms$his84nd_b$xyz], pdb$xyz[atoms$his84ne_b$xyz]) %>% #Find the vector orthogonal to the His1 imidazole ring plane
          vector.angle(pdb$xyz[atoms$his84nd_b$xyz]-pdb$xyz[atoms$cu_b$xyz]) %>% #Find the angle between the orthogonal vector and the Nd1-Cu vector
          {.-90}
      ),
      Molecule = c("A", "B")
    )
  }) %>%
    dplyr::bind_rows() %>% #Data is generated as a list; dplyr::bind_rows turns it into a data frame
    .[order(.$Molecule), ] %>% #Data frame is ordered by alternating subunits, this will order the data frame by subunit
    tibble::tibble(rep(pseudohelixDose, 2), .) #Amend a column containing doses to the data frame

colnames(allAngles)[1] <- "dose_MGy"

stackedAngles <- tibble::tibble(
    Dose = rep(allAngles$dose_MGy, 7), 
    AngleID = c(
      rep("T1", 72), 
      rep("T2", 72), 
      rep("T3", 72), 
      rep("TT", 72), 
      rep("THH", 72), 
      rep("TH1", 72), 
      rep("THN", 72)
    ), 
    Angle = c(
      allAngles$T1, 
      allAngles$T2, 
      allAngles$T3, 
      allAngles$TT, 
      allAngles$THH, 
      allAngles$TH1, 
      allAngles$THN
    ), 
    Molecule = rep(allAngles$Molecule, 7)
  )

# Linear regression analysis ----------------------------------------------

multipleRegressions$Angles <- list(
  T1 = gls(Angle ~ Dose * Molecule, data = stackedAngles, subset = AngleID == "T1", weights = varIdent(form = ~ 1 | Molecule)), 
  T2 = gls(Angle ~ Dose * Molecule, data = stackedAngles, subset = AngleID == "T2", weights = varIdent(form = ~ 1 | Molecule)), 
  T3 = gls(Angle ~ Dose * Molecule, data = stackedAngles, subset = AngleID == "T3", weights = varIdent(form = ~ 1 | Molecule)), 
  TT = gls(Angle ~ Dose * Molecule, data = stackedAngles, subset = AngleID == "TT", weights = varIdent(form = ~ 1 | Molecule)), 
  THH = gls(Angle ~ Dose * Molecule, data = stackedAngles, subset = AngleID == "THH", weights = varIdent(form = ~ 1 | Molecule)), 
  TH1 = gls(Angle ~ Dose * Molecule, data = stackedAngles, subset = AngleID == "TH1", weights = varIdent(form = ~ 1 | Molecule)), 
  THN = gls(Angle ~ Dose * Molecule, data = stackedAngles, subset = AngleID == "THN", weights = varIdent(form = ~ 1 | Molecule))
)

regressionSummaries$Angles <- dplyr::bind_rows(
  lapply(
    c("T1", "T2", "T3", "TT", "THH", "TH1", "THN"), 
    function(atom) {tibble::tibble(
      Measurement = "Angles",
      Residue = atom, 
      Estimate = c("TrendA", "TrendB", "Contrast"), 
      Coefficient = emtrends.coefficient(multipleRegressions$Angles[[atom]], "Dose"), 
      StandardError = emtrends.se(multipleRegressions$Angles[[atom]], "Dose"), 
      PValue = emtrends.pvalue(multipleRegressions$Angles[[atom]], "Dose")
    )}
  )
)

angle_scatter <- multipleRegressions$Angles |> 
  purrr::map(\(model) {
    data <- getData(model)
    
    design_matrix <- model.matrix(model, data = data)
    pred <- predict(model)
    vcov_matrix <- model$varBeta
    
    se <- sqrt(rowSums((design_matrix %*% vcov_matrix) * design_matrix))
    
    alpha <- 0.05
    z_value <- qnorm(1 - alpha / 2)
    
    pred <- bind_cols(
      data,
      y_hat = pred,
      lower = pred - se * z_value,
      upper = pred + se * z_value,
    )
    
    return(pred)
    
  }) |> purrr::list_rbind()

