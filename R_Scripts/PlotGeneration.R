
# Occupancy plotting ------------------------------------------------------

#Write out a PDB file for B-factor coloring in ChimeraX
OccupancyColoredPDB <- pseudohelixList[[1]]
OccupancyColoredPDB$atom$b <- pseudohelixOccSlopes

ggplots$Occupancies$Pseudohelices <- ggplot(stackedOccupancies$Pseudohelices, aes(x = Dose, y = Occupancy, color = Molecule)) +
  scatter_plot() +
  ggtheme_light() +
  facet_wrap(
    ~Residue, 
    scales = "free", 
    labeller = labeller(
      Residue = c(
        "CO2" = "CO[2]", 
        "Ax"  = "H[2]*O[Ax]", 
        "Eq"  = "H[2]*O[Eq]", 
        "Glu" = "Glu[30]", 
        "Oxy" = ""Dioxygen""
      ), 
      .default = label_parsed
    )
  ) +
  labs(
    x = "Density-Weighted Dose (MGy)", 
    y = "Occupancy"
  )

ggplots$Occupancies$Wedges <- ggplot(stackedOccupancies$Wedges, aes(x = WedgeNumber, y = Occupancy, color = Molecule)) +
  scatter_plot() +
  ggtheme_light() +
  facet_wrap(
    ~Residue, 
    scales = "free", 
    labeller = labeller(
      Residue = c(
        "CO2" = "CO[2]", 
        "Ax"  = "H[2]*O[Ax]", 
        "Eq"  = "H[2]*O[Eq]", 
        "Glu" = "Glu[30]", 
        "Oxy" = ""Dioxygen""
      ), 
      .default = label_parsed
    )
  ) +
  labs(
    x = "Wedge Number", 
    y = "Occupancy"
  )



# B-factor plotting -------------------------------------------------------

#Write a model for B-factor coloring in ChimerAxBF
bFactorColoredPDB <- pseudohelixList[[1]]
bFactorColoredPDB$atom$b <- pseudohelixBFSlopes

ggplots$BFactors$Pseudohelices <- ggplot(stackedBFactors$Pseudohelices, aes(x = Dose, y = bFactor, color = Molecule)) +
  scatter_plot() +
  ggtheme_light() +
  facet_wrap(
    vars(Residue), 
    scales = "free", 
    labeller = labeller(
      Residue = c(
        "CO2BF" = "CO[2]", 
        "AxBF"  = "H[2]*O[Ax]", 
        "EqBF"  = "H[2]*O[Eq]", 
        "GluBF" = "Glu[30]", 
        "OxyBF" = ""Dioxygen""
      ), 
      .default = label_parsed
    )
  ) +
  labs(
    x = "Density-Weighted Dose (MGy)", 
    y = "Occupancy"
  )

ggplots$BFactors$Wedges <- ggplot(stackedBFactors$Wedges, aes(x = WedgeNumber, y = bFactor, color = Molecule)) +
  scatter_plot() +
  ggtheme_light() +
  facet_wrap(
    vars(Residue), 
    scales = "free", 
    labeller = labeller(
      Residue = c(
        "CO2BF" = "CO[2]", 
        "AxBF"  = "H[2]*O[Ax]", 
        "EqBF"  = "H[2]*O[Eq]", 
        "GluBF" = "Glu[30]", 
        "OxyBF" = ""Dioxygen""
      ), 
      .default = label_parsed
    )
  ) +
  labs(
    x = "Wedge Number", 
    y = "Occupancy"
  )

# Distance plotting -------------------------------------------------------

plot.distance.pseudohelix <- function(sphere) {
  ggplot(stackedDistances$Pseudohelices[[sphere]], aes(x = Dose, y = Distance, color = Molecule)) +
    scatter_plot() +
    ggtheme_light() +
    theme(
      legend.position = "right"
    ) +
    labs(
      x = "Density-Weighted Dose (MGy)", 
      y = "Distance (Å)"
    )
}

plot.distance.wedge <- function(sphere) {
  ggplot(stackedDistances$Wedges[[sphere]], aes(x = WedgeNumber, y = Distance, color = Molecule)) +
    scatter_plot() +
    ggtheme_light() +
    theme(
      legend.position = "right"
    ) +
    labs(
      x = "Wedge Number", 
      y = "Distance (Å)"
    )
}

ggplots$Distances$DioxygenPseudohelices <- plot.distance.pseudohelix("Oxy") +
  facet_wrap(
    ~ AtomPair, 
    scales = "free", 
    labeller = labeller(
      AtomPair = c(
        "O1-Eq" = "O[prox]*-H[2]*O[Eq]", 
        "O2-Eq" = "O[dist]*-H[2]*O[Eq]", 
        "O1-His157NE2" = "O[prox]*-His[157]*N[ε2]", 
        "O2-His157NE2" = "O[dist]*-His[157]*N[ε2]", 
        "O1-GluOE1" = "O[prox]*-Glu[30]*O[ε1]", 
        "O2-GluOE1" = "O[dist]*-Glu[30]*O[ε1]", 
        "O1-GluOE2" = "O[prox]*-Glu[30]*O[ε2]", 
        "O2-GluOE2" = "O[dist]*-Glu[30]*O[ε2]"
      ), 
      .default = label_parsed
    )
  ) +
  theme(
    legend.position = "inside", 
    legend.position.inside = c(0.85, 0.15)
  )

ggplots$Distances$CopperPseudohelices <- plot.distance.pseudohelix("Cu") +
  facet_wrap(
    ~ AtomPair, 
    scales = "free", 
    labeller = labeller(
      AtomPair = c(
        "Cu-Tyr" = "Cu-Tyr[168]", 
        "Cu-NTerm" = ""Cu-N Terminus"", 
        "Cu-His1ND1" = "Cu-His[1]*Nδ[1]", 
        "Cu-His84NE2" = "Cu-His[84]*Nε[2]", 
        "Cu-Eq" = "Cu-H[2]*O[Eq]", 
        "Cu-Ax" = "Cu-H[2]*O[Ax]"
      ), 
      .default = label_parsed
    )
  )

ggplots$Distances$DioxygenWedges <- plot.distance.wedge("Oxy") +
  facet_wrap(
    ~ AtomPair, 
    scales = "free", 
    labeller = labeller(
      AtomPair = c(
        "O1-Eq" = "O[prox]*-H[2]*O[Eq]", 
        "O2-Eq" = "O[dist]*-H[2]*O[Eq]", 
        "O1-His157NE2" = "O[prox]*-His[157]*N[ε2]", 
        "O2-His157NE2" = "O[dist]*-His[157]*N[ε2]", 
        "O1-GluOE1" = "O[prox]*-Glu[30]*O[ε1]", 
        "O2-GluOE1" = "O[dist]*-Glu[30]*O[ε1]", 
        "O1-GluOE2" = "O[prox]*-Glu[30]*O[ε2]", 
        "O2-GluOE2" = "O[dist]*-Glu[30]*O[ε2]"
      ), 
      .default = label_parsed
    )
  ) +
  theme(
    legend.position = "inside", 
    legend.position.inside = c(0.85, 0.15)
  )

ggplots$Distances$CopperWedges <- plot.distance.wedge("Cu") +
  facet_wrap(
    ~ AtomPair, 
    scales = "free", 
    labeller = labeller(
      AtomPair = c(
        "Cu-Tyr" = "Cu-Tyr[168]", 
        "Cu-NTerm" = ""Cu-N Terminus"", 
        "Cu-His1ND1" = "Cu-His[1]*Nδ[1]", 
        "Cu-His84NE2" = "Cu-His[84]*Nε[2]", 
        "Cu-Eq" = "Cu-H[2]*O[Eq]", 
        "Cu-Ax" = "Cu-H[2]*O[Ax]"
      ), 
      .default = label_parsed
    )
  )

# Angle Plotting ----------------------------------------------------------

ggplots$Angles$Pseudohelices <- ggplot(stackedAngles$Pseudohelices, aes(x = Dose, y = Angle, color = Molecule)) +
  scatter_plot() +
  ggtheme_light() +
  theme(legend.position.inside = c(0.85, 0.15)) +
  facet_wrap(
    ~AngleID, 
    scales = "free", 
    labeller = labeller(
      AngleID = c(
        "T1" = "θ[1]", 
        "T2" = "θ[2]", 
        "T3" = "θ[3]", 
        "TT" = "θ[T]", 
        "THH" = "θ[HH]", 
        "TH1HN1" = "θH[1]*HN[1]", 
        "TH1HN84" = "θH[1]*HN[84]", 
        "TOxy" = "θ[oxygen]"
      ), 
      .default = label_parsed
    )
  ) +
  labs(
    x = "Density-Weighted Dose (MGy)", 
    y = "Angle (°)"
  )

ggplots$Angles$Wedges <- ggplot(stackedAngles$Wedges, aes(x = WedgeNumber, y = Angle, color = Molecule)) +
  scatter_plot() +
  ggtheme_light() +
  theme(legend.position.inside = c(0.85, 0.15)) +
  facet_wrap(
    ~AngleID, 
    scales = "free", 
    labeller = labeller(
      AngleID = c(
        "T1" = "θ[1]", 
        "T2" = "θ[2]", 
        "T3" = "θ[3]", 
        "TT" = "θ[T]", 
        "THH" = "θ[HH]", 
        "TH1HN1" = "θH[1]*HN[1]", 
        "TH1HN84" = "θH[1]*HN[84]", 
        "TOxy" = "θ[oxygen]"
      ), 
      .default = label_parsed
    )
  ) +
  labs(
    x = "Wedge Number", 
    y = "Angle (°)"
  )
