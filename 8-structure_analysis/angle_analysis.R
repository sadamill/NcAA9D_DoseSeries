#This script is meant to be run from Master_Script.R

#This script analyzes the variance of various angles as a function of X-ray dose
#It  generates a plot of dose vs. angle for important atom sets, allowing to highlight differences
#between different molecules in the asymmetric unit

# Function definitions ----------------------------------------------------

# Calculates a vector orthogonal to a plane abc, where a, b, and c are three xyz points on a plane
# point a is the point from which both vectors extend
orthogonal.vector <- function(a, b, c) {
  v1 <- b - a
  v2 <- c - a
  
  c(v1[2]*v2[3] - v1[3]*v2[2],
    v1[3]*v2[1] - v1[1]*v2[3],
    v1[1]*v2[2] - v1[2]*v2[1])
}



# finds angle (°) between vectors v1 and v2
vector.angle <- function(v1, v2) {
  v1_mag <- sqrt(sum(v1^2))
  v2_mag <- sqrt(sum(v2^2))
  
  acos((v1 %*% v2) / (v1_mag * v2_mag)) * 180/pi
}

# Angle calculations ------------------------------------------------------

#Calculate angles for all datasets based off XYZ indices calculated earlier
allAngles <- lapply(pseudohelixList, function(pdb) {
    
    tibble::tibble(
      T1 = c(
        bio3d::angle.xyz(
          c(pdb$xyz[atoms$nterm_a$xyz], 
            pdb$xyz[atoms$cu_a$xyz],
            pdb$xyz[atoms$his1nd_a$xyz])
        ), 
        bio3d::angle.xyz(
          c(pdb$xyz[atoms$nterm_b$xyz], 
            pdb$xyz[atoms$cu_b$xyz],
            pdb$xyz[atoms$his1nd_b$xyz])
        )
      ), 
      T2 = c(
        bio3d::angle.xyz(
          c(pdb$xyz[atoms$nterm_a$xyz],
            pdb$xyz[atoms$cu_a$xyz],
            pdb$xyz[atoms$his84ne_a$xyz])
        ), 
        bio3d::angle.xyz(
          c(pdb$xyz[atoms$nterm_b$xyz],
            pdb$xyz[atoms$cu_b$xyz],
            pdb$xyz[atoms$his84ne_b$xyz])
        )
      ), 
      T3 = c(
        bio3d::angle.xyz(
          c(pdb$xyz[atoms$his1nd_a$xyz],
            pdb$xyz[atoms$cu_a$xyz], 
            pdb$xyz[atoms$his84ne_a$xyz])
        ), 
        bio3d::angle.xyz(
          c(pdb$xyz[atoms$his1nd_b$xyz],
            pdb$xyz[atoms$cu_b$xyz],
            pdb$xyz[atoms$his84ne_b$xyz])
        )
      ), 
      TT = c(
        orthogonal.vector(
          pdb$xyz[atoms$cu_a$xyz], 
          pdb$xyz[atoms$nterm_a$xyz], 
          pdb$xyz[atoms$his1nd_a$xyz]
        ) %>% #Find the vector orthogonal to the Cu-Nterm-Nd1 plane
          vector.angle(
            pdb$xyz[atoms$his84ne_a$xyz]-pdb$xyz[atoms$cu_a$xyz]
          ) %>% #Find the angle between the orthogonal vector and the Ne2-Cu vector
          {90-.}, 
        orthogonal.vector(
          pdb$xyz[atoms$cu_b$xyz], 
          pdb$xyz[atoms$nterm_b$xyz], 
          pdb$xyz[atoms$his1nd_b$xyz]
        ) %>% #Repeat for molecule B
          vector.angle(
            pdb$xyz[atoms$his84ne_b$xyz]-pdb$xyz[atoms$cu_b$xyz]
          ) %>%
          {90-.}
      ), 
      THH = c(
        orthogonal.vector(pdb$xyz[atoms$his1cg_a$xyz], 
                          pdb$xyz[atoms$his1nd_a$xyz],
                          pdb$xyz[atoms$his1ne_a$xyz]) %>% #Find vector orthogonal to His1 imidazole ring plane
          vector.angle(
            orthogonal.vector(pdb$xyz[atoms$his84cg_a$xyz],
                              pdb$xyz[atoms$his84nd_a$xyz],
                              pdb$xyz[atoms$his84ne_a$xyz])#Find angle between the two vectors orthogonal to the His1-His84 imidazole planes
          ), 
        orthogonal.vector(pdb$xyz[atoms$his1cg_b$xyz],
                          pdb$xyz[atoms$his1nd_b$xyz],
                          pdb$xyz[atoms$his1ne_b$xyz]) %>% #Repeat for molecule B
          vector.angle(
            orthogonal.vector(pdb$xyz[atoms$his84cg_b$xyz],
                              pdb$xyz[atoms$his84nd_b$xyz],
                              pdb$xyz[atoms$his84ne_b$xyz])
          )
      ), 
      TH1 = c(
        orthogonal.vector(pdb$xyz[atoms$his1cg_a$xyz], 
                          pdb$xyz[atoms$his1nd_a$xyz],
                          pdb$xyz[atoms$his1ne_a$xyz]) %>% #Find the vector orthogonal to the His1 imidazole ring plane
          vector.angle(pdb$xyz[atoms$his1nd_a$xyz]-
                         pdb$xyz[atoms$cu_a$xyz]) %>% #Find the angle between the orthogonal vector and the Nd1-Cu vector
          {90-.}, 
        orthogonal.vector(pdb$xyz[atoms$his1cg_b$xyz],
                          pdb$xyz[atoms$his1nd_b$xyz],
                          pdb$xyz[atoms$his1ne_b$xyz]) %>% #Find the vector orthogonal to the His1 imidazole ring plane
          vector.angle(pdb$xyz[atoms$his1nd_b$xyz]-
                         pdb$xyz[atoms$cu_b$xyz]) %>% #Find the angle between the orthogonal vector and the Nd1-Cu vector
          {90-.}
      ), 
      THN = c(
        orthogonal.vector(pdb$xyz[atoms$his84cg_a$xyz],
                          pdb$xyz[atoms$his84nd_a$xyz],
                          pdb$xyz[atoms$his84ne_a$xyz]) %>% #Find the vector orthogonal to the His1 imidazole ring plane
          vector.angle(pdb$xyz[atoms$his84nd_a$xyz]
                       -pdb$xyz[atoms$cu_a$xyz]) %>% #Find the angle between the orthogonal vector and the Nd1-Cu vector
          {90-.}, 
        orthogonal.vector(pdb$xyz[atoms$his84cg_b$xyz],
                          pdb$xyz[atoms$his84nd_b$xyz],
                          pdb$xyz[atoms$his84ne_b$xyz]) %>% #Find the vector orthogonal to the His1 imidazole ring plane
          vector.angle(pdb$xyz[atoms$his84nd_b$xyz]
                       -pdb$xyz[atoms$cu_b$xyz]) %>% #Find the angle between the orthogonal vector and the Nd1-Cu vector
          {90-.}
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

