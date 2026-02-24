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
allAngles <- list(
  Pseudohelices = lapply(pseudohelixList, function(pdb) {
    
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
    tibble::tibble(rep(pseudohelixDose, 2), .), #Amend a column containing doses to the data frame
  Wedges = lapply(wedgeList, function(pdb) {
    
    data.frame(
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
      TH1 = c( #TH1/HN angle for His1
        orthogonal.vector(pdb$xyz[atoms$his1cg_a$xyz], pdb$xyz[atoms$his1nd_a$xyz], pdb$xyz[atoms$his1ne_a$xyz]) %>% #Find the vector orthogonal to the His1 imidazole ring plane
          vector.angle(pdb$xyz[atoms$his1nd_a$xyz]-pdb$xyz[atoms$cu_a$xyz]) %>% #Find the angle between the orthogonal vector and the Nd1-Cu vector
          {.-90}, 
        orthogonal.vector(pdb$xyz[atoms$his1cg_b$xyz], pdb$xyz[atoms$his1nd_b$xyz], pdb$xyz[atoms$his1ne_b$xyz]) %>% #Find the vector orthogonal to the His1 imidazole ring plane
          vector.angle(pdb$xyz[atoms$his1nd_b$xyz]-pdb$xyz[atoms$cu_b$xyz]) %>% #Find the angle between the orthogonal vector and the Nd1-Cu vector
          {.-90}
      ), 
      THN = c( #TH1/HN angle for His84
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
    .[order(.$Molecule), ] %>% #Data frame is ordered by alternating molecules, this will order the data frame by subunit
    tibble::tibble(rep(1:36, 2), .) #Amend a column containing doses to the data frame
)

colnames(allAngles$Pseudohelices)[1] <- "dose_MGy"
colnames(allAngles$Wedges)[1] <- "WedgeNumber"

stackedAngles <- list(
  Pseudohelices = tibble::tibble(
    Dose = rep(allAngles$Pseudohelices$dose_MGy, 7), 
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
      allAngles$Pseudohelices$T1, 
      allAngles$Pseudohelices$T2, 
      allAngles$Pseudohelices$T3, 
      allAngles$Pseudohelices$TT, 
      allAngles$Pseudohelices$THH, 
      allAngles$Pseudohelices$TH1, 
      allAngles$Pseudohelices$THN
    ), 
    Molecule = rep(allAngles$Pseudohelices$Molecule, 7)
  ), 
  Wedges = tibble::tibble(
    WedgeNumber = rep(allAngles$Wedges$WedgeNumber, 7), 
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
      allAngles$Wedges$T1, 
      allAngles$Wedges$T2, 
      allAngles$Wedges$T3, 
      allAngles$Wedges$TT, 
      allAngles$Wedges$THH, 
      allAngles$Wedges$TH1, 
      allAngles$Wedges$THN
    ), 
    Molecule = rep(allAngles$Wedges$Molecule, 7)
  )
)

# Linear regression analysis ----------------------------------------------

multipleRegressions$Angles <- list(
  Pseudohelices = list (
    T1 = gls(T1 ~ dose_MGy * Molecule, data = allAngles$Pseudohelices, weights = varIdent(form = ~ 1 | Molecule)), 
    T2 = gls(T2 ~ dose_MGy * Molecule, data = allAngles$Pseudohelices, weights = varIdent(form = ~ 1 | Molecule)), 
    T3 = gls(T3 ~ dose_MGy * Molecule, data = allAngles$Pseudohelices, weights = varIdent(form = ~ 1 | Molecule)), 
    TT = gls(TT ~ dose_MGy * Molecule, data = allAngles$Pseudohelices, weights = varIdent(form = ~ 1 | Molecule)), 
    THH = gls(THH ~ dose_MGy * Molecule, data = allAngles$Pseudohelices, weights = varIdent(form = ~ 1 | Molecule)), 
    TH1 = gls(TH1 ~ dose_MGy * Molecule, data = allAngles$Pseudohelices, weights = varIdent(form = ~ 1 | Molecule)), 
    THN = gls(THN ~ dose_MGy * Molecule, data = allAngles$Pseudohelices, weights = varIdent(form = ~ 1 | Molecule))
  ), 
  Wedges = list (
    T1 = gls(T1 ~ WedgeNumber * Molecule, data = allAngles$Wedges, weights = varIdent(form = ~ 1 | Molecule)), 
    T2 = gls(T2 ~ WedgeNumber * Molecule, data = allAngles$Wedges, weights = varIdent(form = ~ 1 | Molecule)), 
    T3 = gls(T3 ~ WedgeNumber * Molecule, data = allAngles$Wedges, weights = varIdent(form = ~ 1 | Molecule)), 
    TT = gls(TT ~ WedgeNumber * Molecule, data = allAngles$Wedges, weights = varIdent(form = ~ 1 | Molecule)), 
    THH = gls(THH ~ WedgeNumber * Molecule, data = allAngles$Wedges, weights = varIdent(form = ~ 1 | Molecule)), 
    TH1 = gls(TH1 ~ WedgeNumber * Molecule, data = allAngles$Wedges, weights = varIdent(form = ~ 1 | Molecule)), 
    THN = gls(THN ~ WedgeNumber * Molecule, data = allAngles$Wedges, weights = varIdent(form = ~ 1 | Molecule))
  )
)

regressionSummaries$Angles <- list(
  Pseudohelices = dplyr::bind_rows(
    lapply(
      c("T1", "T2", "T3", "TT", "THH", "TH1", "THN"), 
      function(atom) {
        tibble::tibble(
          Measurement = "Angles",
          Residue = atom, 
          Estimate = c("TrendA", "TrendB", "Contrast"), 
          Coefficient = emtrends.coefficient(multipleRegressions$Angles$Pseudohelices[[atom]], "dose_MGy"), 
          StandardError = emtrends.se(multipleRegressions$Angles$Pseudohelices[[atom]], "dose_MGy"), 
          PValue = emtrends.pvalue(multipleRegressions$Angles$Pseudohelices[[atom]], "dose_MGy")
        )
      }
    )
  ), 
  Wedges = dplyr::bind_rows(
    lapply(
      c("T1", "T2", "T3", "TT", "THH", "TH1", "THN"), 
      function(atom) {
        tibble::tibble(
          Measurement = "Angles",
          Residue = atom, 
          Estimate = c("TrendA", "TrendB", "Contrast"), 
          Coefficient = emtrends.coefficient(multipleRegressions$Angles$Wedges[[atom]], "WedgeNumber"), 
          StandardError = emtrends.se(multipleRegressions$Angles$Wedges[[atom]], "WedgeNumber"), 
          PValue = emtrends.pvalue(multipleRegressions$Angles$Wedges[[atom]], "WedgeNumber")
        )
      }
    )
  )
)
