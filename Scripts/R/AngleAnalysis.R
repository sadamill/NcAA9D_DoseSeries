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
    {(dot(.[1, ], .[2, ]))/(sqrt(sum(.[1, 1]^2, .[1, 2]^2, .[1, 3]^2))*sqrt(sum(.[2, 1]^2, .[2, 2]^2, .[2, 3]^2)))} %>% 
    acos() %>% 
    {.*180/pi}
}

# Angle calculations ------------------------------------------------------

#Create vectors containing coordinates for active site atoms
ACu <- atom.select(pseudohelixList[[1]], resno=1, chain="C")$xyz
A1N <- atom.select(pseudohelixList[[1]], resno=1, chain="A", elety="N")$xyz
A1ND1 <- atom.select(pseudohelixList[[1]], resno=1, chain="A", elety="ND1")$xyz
A1NE2 <- atom.select(pseudohelixList[[1]], resno=1, chain="A", elety="NE2")$xyz
A1CG <- atom.select(pseudohelixList[[1]], resno=1, chain="A", elety="CG")$xyz
A84ND1 <- atom.select(pseudohelixList[[1]], resno=84, chain="A", elety="ND1")$xyz
A84NE2 <- atom.select(pseudohelixList[[1]], resno=84, chain="A", elety="NE2")$xyz
A84CG <- atom.select(pseudohelixList[[1]], resno=84, chain="A", elety="CG")$xyz
AO1 <- atom.select(pseudohelixList[[1]], resno=1, chain="D", elety="O1")$xyz
AO2 <- atom.select(pseudohelixList[[1]], resno=1, chain="D", elety="O2")$xyz

BCu <- atom.select(pseudohelixList[[1]], resno=2, chain="C")$xyz
B1N <- atom.select(pseudohelixList[[1]], resno=1, chain="B", elety="N")$xyz
B1ND1 <- atom.select(pseudohelixList[[1]], resno=1, chain="B", elety="ND1")$xyz
B1NE2 <- atom.select(pseudohelixList[[1]], resno=1, chain="B", elety="NE2")$xyz
B1CG <- atom.select(pseudohelixList[[1]], resno=1, chain="B", elety="CG")$xyz
B84ND1 <- atom.select(pseudohelixList[[1]], resno=84, chain="B", elety="ND1")$xyz
B84NE2 <- atom.select(pseudohelixList[[1]], resno=84, chain="B", elety="NE2")$xyz
B84CG <- atom.select(pseudohelixList[[1]], resno=84, chain="B", elety="CG")$xyz
BO1 <- atom.select(pseudohelixList[[1]], resno=2, chain="D", elety="O1")$xyz
BO2 <- atom.select(pseudohelixList[[1]], resno=2, chain="D", elety="O2")$xyz

#Calculate angles for all datasets based off XYZ indices calculated earlier
allAngles <- list(
  Pseudohelices = lapply(pseudohelixList, function(pdb) {
    data.frame(
      T1 = c(
        angle.xyz(c(pdb$xyz[A1N], pdb$xyz[ACu], pdb$xyz[A1ND1])), 
        angle.xyz(c(pdb$xyz[B1N], pdb$xyz[BCu], pdb$xyz[B1ND1]))
      ), 
      T2 = c(
        angle.xyz(c(pdb$xyz[A1N], pdb$xyz[ACu], pdb$xyz[A84NE2])), 
        angle.xyz(c(pdb$xyz[B1N], pdb$xyz[BCu], pdb$xyz[B84NE2]))
      ), 
      T3 = c(
        angle.xyz(c(pdb$xyz[A1ND1], pdb$xyz[ACu], pdb$xyz[A84NE2])), 
        angle.xyz(c(pdb$xyz[B1ND1], pdb$xyz[BCu], pdb$xyz[B84NE2]))
      ), 
      TT = c(
        orthogonal.vector(pdb$xyz[ACu], pdb$xyz[A1N], pdb$xyz[A1ND1]) %>% #Find the vector orthogonal to the Cu-Nterm-Nd1 plane
          vector.angle(pdb$xyz[A84NE2]-pdb$xyz[ACu]) %>% #Find the angle between the orthogonal vector and the Ne2-Cu vector
          {90-.}, 
        orthogonal.vector(pdb$xyz[BCu], pdb$xyz[B1N], pdb$xyz[B1ND1]) %>% #Repeat for molecule B
          vector.angle(pdb$xyz[B84NE2]-pdb$xyz[BCu]) %>%
          {90-.}
      ), 
      THH = c(
        orthogonal.vector(pdb$xyz[A1CG], pdb$xyz[A1ND1], pdb$xyz[A1NE2]) %>% #Find vector orthogonal to His1 imidazole ring plane
          vector.angle(
            orthogonal.vector(pdb$xyz[A84CG], pdb$xyz[A84ND1], pdb$xyz[A84NE2])#Find angle between the two vectors orthogonal to the His1-His84 imidazole planes
          ), 
        orthogonal.vector(pdb$xyz[B1CG], pdb$xyz[B1ND1], pdb$xyz[B1NE2]) %>% #Repeat for molecule B
          vector.angle(
            orthogonal.vector(pdb$xyz[B84CG], pdb$xyz[B84ND1], pdb$xyz[B84NE2])
          )
      ), 
      TH1HN1 = c( #TH1/HN angle for His1
        orthogonal.vector(pdb$xyz[A1CG], pdb$xyz[A1ND1], pdb$xyz[A1NE2]) %>% #Find the vector orthogonal to the His1 imidazole ring plane
          vector.angle(pdb$xyz[A1ND1]-pdb$xyz[ACu]) %>% #Find the angle between the orthogonal vector and the Nd1-Cu vector
          {90+.}, 
        orthogonal.vector(pdb$xyz[B1CG], pdb$xyz[B1ND1], pdb$xyz[B1NE2]) %>% #Find the vector orthogonal to the His1 imidazole ring plane
          vector.angle(pdb$xyz[B1ND1]-pdb$xyz[BCu]) %>% #Find the angle between the orthogonal vector and the Nd1-Cu vector
          {90+.}
      ), 
      TH1HN84 = c( #TH1/HN angle for His84
        orthogonal.vector(pdb$xyz[A84CG], pdb$xyz[A84ND1], pdb$xyz[A84NE2]) %>% #Find the vector orthogonal to the His1 imidazole ring plane
          vector.angle(pdb$xyz[A84ND1]-pdb$xyz[ACu]) %>% #Find the angle between the orthogonal vector and the Nd1-Cu vector
          {90+.}, 
        orthogonal.vector(pdb$xyz[B84CG], pdb$xyz[B84ND1], pdb$xyz[B84NE2]) %>% #Find the vector orthogonal to the His1 imidazole ring plane
          vector.angle(pdb$xyz[B84ND1]-pdb$xyz[BCu]) %>% #Find the angle between the orthogonal vector and the Nd1-Cu vector
          {90+.}
      ), 
      TOxy = c(
        angle.xyz(c(pdb$xyz[ACu], pdb$xyz[AO1], pdb$xyz[AO2])), 
        angle.xyz(c(pdb$xyz[BCu], pdb$xyz[BO1], pdb$xyz[BO2]))
      ), 
      Molecule = c("A", "B")
    )
  }) %>%
    ldply() %>% #Data is generated as a list; ldply turns it into a data frame
    .[order(.$Molecule), ] %>% #Data frame is ordered by alternating subunits, this will order the data frame by subunit
    data.frame(pseudohelixDose, .), #Amend a column containing doses to the data frame
  Wedges = lapply(wedgeList, function(pdb) {
    data.frame(
      T1 = c(
        angle.xyz(c(pdb$xyz[A1N], pdb$xyz[ACu], pdb$xyz[A1ND1])), 
        angle.xyz(c(pdb$xyz[B1N], pdb$xyz[BCu], pdb$xyz[B1ND1]))
      ), 
      T2 = c(
        angle.xyz(c(pdb$xyz[A1N], pdb$xyz[ACu], pdb$xyz[A84NE2])), 
        angle.xyz(c(pdb$xyz[B1N], pdb$xyz[BCu], pdb$xyz[B84NE2]))
      ), 
      T3 = c(
        angle.xyz(c(pdb$xyz[A1ND1], pdb$xyz[ACu], pdb$xyz[A84NE2])), 
        angle.xyz(c(pdb$xyz[B1ND1], pdb$xyz[BCu], pdb$xyz[B84NE2]))
      ), 
      TT = c(
        orthogonal.vector(pdb$xyz[ACu], pdb$xyz[A1N], pdb$xyz[A1ND1]) %>% #Find the vector orthogonal to the Cu-Nterm-Nd1 plane
          vector.angle(pdb$xyz[A84NE2]-pdb$xyz[ACu]) %>% #Find the angle between the orthogonal vector and the Ne2-Cu vector
          {90-.}, 
        orthogonal.vector(pdb$xyz[BCu], pdb$xyz[B1N], pdb$xyz[B1ND1]) %>% #Repeat for molecule B
          vector.angle(pdb$xyz[B84NE2]-pdb$xyz[BCu]) %>%
          {90-.}
      ), 
      THH = c(
        orthogonal.vector(pdb$xyz[A1CG], pdb$xyz[A1ND1], pdb$xyz[A1NE2]) %>% #Find vector orthogonal to His1 imidazole ring plane
          vector.angle(
            orthogonal.vector(pdb$xyz[A84CG], pdb$xyz[A84ND1], pdb$xyz[A84NE2])#Find angle between the two vectors orthogonal to the His1-His84 imidazole planes
          ), 
        orthogonal.vector(pdb$xyz[B1CG], pdb$xyz[B1ND1], pdb$xyz[B1NE2]) %>% #Repeat for molecule B
          vector.angle(
            orthogonal.vector(pdb$xyz[B84CG], pdb$xyz[B84ND1], pdb$xyz[B84NE2])
          )
      ), 
      TH1HN1 = c( #TH1/HN angle for His1
        orthogonal.vector(pdb$xyz[A1CG], pdb$xyz[A1ND1], pdb$xyz[A1NE2]) %>% #Find the vector orthogonal to the His1 imidazole ring plane
          vector.angle(pdb$xyz[A1ND1]-pdb$xyz[ACu]) %>% #Find the angle between the orthogonal vector and the Nd1-Cu vector
          {90+.}, 
        orthogonal.vector(pdb$xyz[B1CG], pdb$xyz[B1ND1], pdb$xyz[B1NE2]) %>% #Find the vector orthogonal to the His1 imidazole ring plane
          vector.angle(pdb$xyz[B1ND1]-pdb$xyz[BCu]) %>% #Find the angle between the orthogonal vector and the Nd1-Cu vector
          {90+.}
      ), 
      TH1HN84 = c( #TH1/HN angle for His84
        orthogonal.vector(pdb$xyz[A84CG], pdb$xyz[A84ND1], pdb$xyz[A84NE2]) %>% #Find the vector orthogonal to the His1 imidazole ring plane
          vector.angle(pdb$xyz[A84ND1]-pdb$xyz[ACu]) %>% #Find the angle between the orthogonal vector and the Nd1-Cu vector
          {90+.}, 
        orthogonal.vector(pdb$xyz[B84CG], pdb$xyz[B84ND1], pdb$xyz[B84NE2]) %>% #Find the vector orthogonal to the His1 imidazole ring plane
          vector.angle(pdb$xyz[B84ND1]-pdb$xyz[BCu]) %>% #Find the angle between the orthogonal vector and the Nd1-Cu vector
          {90+.}
      ), 
      TOxy = c(
        angle.xyz(c(pdb$xyz[ACu], pdb$xyz[AO1], pdb$xyz[AO2])), 
        angle.xyz(c(pdb$xyz[BCu], pdb$xyz[BO1], pdb$xyz[BO2]))
      ), 
      Molecule = c("A", "B")
    )
  }) %>%
    ldply() %>% #Data is generated as a list; ldply turns it into a data frame
    .[order(.$Molecule), ] %>% #Data frame is ordered by alternating molecules, this will order the data frame by subunit
    data.frame(1:36, .) #Amend a column containing doses to the data frame
)

colnames(allAngles$Pseudohelices)[1] <- "dose_MGy"
colnames(allAngles$Wedges)[1] <- "WedgeNumber"

stackedAngles <- list(
  Pseudohelices = tibble(
    Dose = rep(allAngles$Pseudohelices$dose_MGy, 8), 
    AngleID = c(
      rep("T1", 72), 
      rep("T2", 72), 
      rep("T3", 72), 
      rep("TT", 72), 
      rep("THH", 72), 
      rep("TH1HN1", 72), 
      rep("TH1HN84", 72), 
      rep("TOxy", 72)
    ), 
    Angle = c(
      allAngles$Pseudohelices$T1, 
      allAngles$Pseudohelices$T2, 
      allAngles$Pseudohelices$T3, 
      allAngles$Pseudohelices$TT, 
      allAngles$Pseudohelices$THH, 
      allAngles$Pseudohelices$TH1HN1, 
      allAngles$Pseudohelices$TH1HN84, 
      allAngles$Pseudohelices$TOxy
    ), 
    Molecule = rep(allAngles$Pseudohelices$Molecule, 8)
  ), 
  Wedges = tibble(
    WedgeNumber = rep(allAngles$Wedges$WedgeNumber, 8), 
    AngleID = c(
      rep("T1", 72), 
      rep("T2", 72), 
      rep("T3", 72), 
      rep("TT", 72), 
      rep("THH", 72), 
      rep("TH1HN1", 72), 
      rep("TH1HN84", 72), 
      rep("TOxy", 72)
    ), 
    Angle = c(
      allAngles$Wedges$T1, 
      allAngles$Wedges$T2, 
      allAngles$Wedges$T3, 
      allAngles$Wedges$TT, 
      allAngles$Wedges$THH, 
      allAngles$Wedges$TH1HN1, 
      allAngles$Wedges$TH1HN84, 
      allAngles$Wedges$TOxy
    ), 
    Molecule = rep(allAngles$Wedges$Molecule, 8)
  )
)

# Linear regression analysis ----------------------------------------------

multipleRegressions$Angles <- list(
  Pseudohelices = list (
    T1 = lm(T1 ~ dose_MGy * Molecule, data = allAngles$Pseudohelices), 
    T2 = lm(T2 ~ dose_MGy * Molecule, data = allAngles$Pseudohelices), 
    T3 = lm(T3 ~ dose_MGy * Molecule, data = allAngles$Pseudohelices), 
    TT = lm(TT ~ dose_MGy * Molecule, data = allAngles$Pseudohelices), 
    THH = lm(THH ~ dose_MGy * Molecule, data = allAngles$Pseudohelices), 
    TH1HN1 = lm(TH1HN1 ~ dose_MGy * Molecule, data = allAngles$Pseudohelices), 
    TH1HN84 = lm(TH1HN84 ~ dose_MGy * Molecule, data = allAngles$Pseudohelices), 
    TOxy = lm(TOxy ~ dose_MGy * Molecule, data = allAngles$Pseudohelices)
  ), 
  Wedges = list (
    T1 = lm(T1 ~ WedgeNumber * Molecule, data = allAngles$Wedges), 
    T2 = lm(T2 ~ WedgeNumber * Molecule, data = allAngles$Wedges), 
    T3 = lm(T3 ~ WedgeNumber * Molecule, data = allAngles$Wedges), 
    TT = lm(TT ~ WedgeNumber * Molecule, data = allAngles$Wedges), 
    THH = lm(THH ~ WedgeNumber * Molecule, data = allAngles$Wedges), 
    TH1HN1 = lm(TH1HN1 ~ WedgeNumber * Molecule, data = allAngles$Wedges), 
    TH1HN84 = lm(TH1HN84 ~ WedgeNumber * Molecule, data = allAngles$Wedges), 
    TOxy = lm(TOxy ~ WedgeNumber * Molecule, data = allAngles$Wedges)
  )
)

regressionSummaries$Angles <- list(
  Pseudohelices = bind_rows(
    lapply(
      c("T1", "T2", "T3", "TT", "THH", "TH1HN1", "TH1HN84", "TOxy"), 
      function(atom) {
        tibble(
          Residue = atom, 
          Estimate = c("TrendA", "TrendB", "Contrast"), 
          Coefficient = emtrends.coefficient(multipleRegressions$Angles$Pseudohelices[[atom]], "dose_MGy"), 
          StandardError = emtrends.se(multipleRegressions$Angles$Pseudohelices[[atom]], "dose_MGy"), 
          ModelRSquared = summary(multipleRegressions$Angles$Pseudohelices[[atom]])$r.squared, 
          PValue = emtrends.pvalue(multipleRegressions$Angles$Pseudohelices[[atom]], "dose_MGy")
        ) %>% mutate(Significance = ifelse(
          PValue <= 0.001, "***", 
          ifelse(
            PValue <= 0.01, "**", 
            ifelse(
              PValue <= 0.05, "*", " "
            )
          )
        ))
      }
    )
  ), 
  Wedges = bind_rows(
    lapply(
      c("T1", "T2", "T3", "TT", "THH", "TH1HN1", "TH1HN84", "TOxy"), 
      function(atom) {
        tibble(
          Residue = atom, 
          Estimate = c("TrendA", "TrendB", "Contrast"), 
          Coefficient = emtrends.coefficient(multipleRegressions$Angles$Wedges[[atom]], "WedgeNumber"), 
          StandardError = emtrends.se(multipleRegressions$Angles$Wedges[[atom]], "WedgeNumber"), 
          ModelRSquared = summary(multipleRegressions$Angles$Wedges[[atom]])$r.squared, 
          PValue = emtrends.pvalue(multipleRegressions$Angles$Wedges[[atom]], "WedgeNumber")
        ) %>% mutate(Significance = ifelse(
          PValue <= 0.001, "***", 
          ifelse(
            PValue <= 0.01, "**", 
            ifelse(
              PValue <= 0.05, "*", " "
            )
          )
        ))
      }
    )
  )
)
