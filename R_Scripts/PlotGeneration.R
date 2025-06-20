ggtheme_dark <- function() {
  list(
    scale_color_manual(
      "Molecule", 
      labels = c("A", "B"), 
      values = c("#5cb344", "#8100b6")
    ), 
    theme_dark(), 
    theme(
      #Overall elements
      rect = element_blank(), 
      text = element_text(color = "white"), 
      line = element_line(color = "black"), 
      
      #Legend positioning
      legend.position = "inside", 
      legend.position.inside = c(0.85, 0.25), 
      
      #Manual override of desired theme elements
      panel.background = element_rect(fill = "black", color = "white"), 
      legend.key = element_blank(), #panel.background automatically maps to legend.key; I want to override this
      legend.background = element_rect(fill = "black", color = "white"), 
      strip.background = element_rect(fill = "black", color = "white")
    ),
    coord_cartesian(expand = FALSE)
  )
} #Global ggplot dark theme
ggtheme_light <- function() {
  list(
    scale_color_manual(
      "Molecule", 
      labels = c("A", "B"), 
      values = c("#5cb344", "#8100b6"), 
    ), 
    theme_bw(), 
    theme(
      #Overall elements
      text = element_text(color = "black"), 
      
      #Legend positioning
      legend.position = "inside", 
      legend.position.inside = c(0.85, 0.25), 
      
      #Manual override of desired theme elements
      legend.background = element_rect(fill = "white", color = "black"), 
      strip.background = element_rect(fill = "white"), 
      plot.background = element_blank()
    ), 
    coord_cartesian(expand = FALSE)
  )
} #Global ggplot light theme
faceting <- function(facetVar, datasetType) {
    facet <- do.call(
      what = facet_wrap,
      args = list(
        facet = as.formula(paste("~", facetVar)), 
        scales = 'free',
        labeller = labeller(
          list(
            Residue = as_labeller(c(
              "CO2" = "CO[2]", 
              "Ax"  = "H[2]*O[Ax]", 
              "Eq"  = "H[2]*O[Eq]", 
              "Glu" = "Glu[30]", 
              "Oxy" = "'Dioxygen'"
            )),
            bResidue = as_labeller(c(
              "CO2BF" = "CO[2]", 
              "AxBF"  = "H[2]*O[Ax]", 
              "EqBF"  = "H[2]*O[Eq]", 
              "GluBF" = "Glu[30]", 
              "OxyBF" = "'Dioxygen'"
            )),
            AtomPair = as_labeller(c(
              "O1-Eq" = "O[prox]*-H[2]*O[Eq]", 
              "O2-Eq" = "O[dist]*-H[2]*O[Eq]", 
              "O1-His157NE2" = "O[prox]*-His[157]*N[ε2]", 
              "O2-His157NE2" = "O[dist]*-His[157]*N[ε2]", 
              "O1-GluOE1" = "O[prox]*-Glu[30]*O[ε1]", 
              "O2-GluOE1" = "O[dist]*-Glu[30]*O[ε1]", 
              "O1-GluOE2" = "O[prox]*-Glu[30]*O[ε2]", 
              "O2-GluOE2" = "O[dist]*-Glu[30]*O[ε2]"
            )),
            AngleID = as_labeller(c(
              "T1" = "θ[1]", 
              "T2" = "θ[2]", 
              "T3" = "θ[3]", 
              "TT" = "θ[T]", 
              "THH" = "θ[HH]", 
              "TH1HN1" = "θH[1]*HN[1]", 
              "TH1HN84" = "θH[1]*HN[84]", 
              "TOxy" = "θ[oxygen]"
            )
            )),
          .default = label_parsed
        )
      ))
    
    labs <- do.call(
      what = labs,
      args = list(
        x = if(datasetType == 'Pseudohelix') {
          "Dose (MGy)"
        } else if(datasetType == "Wedge") {
          "Wedge Number"
        },
        y = if(facetVar == "Residue") {
          "Occupancy"
        } else if(facetVar == "bResidue") {
          bquote("B-Factor (" * Å^2 * ")")
        } else if(facetVar == "AtomPair") {
          "Distance (Å)"
        } else if(facetVar == "AngleID") {
          "Angle (°)"
        }
      )
    )
    
    list(facet, labs)
}

scatter_dark <- function(data, mapping, facetVar, datasetType) {
  ggplot(data, mapping) +
    stat_smooth(
      method = "lm", 
      linewidth = 0, 
      fill = "gray", 
      show.legend = FALSE
    ) + #Standard error plotting
    stat_smooth(
      method = "lm", 
      linetype = 2, 
      se = FALSE
    ) + #Linear regression line
    geom_point(
      size = 0.3, 
      show.legend = FALSE
    ) + #Point for each occupancy value
    faceting(facetVar = facetVar, datasetType = datasetType) +
    ggtheme_dark()
} #Make scatter plot with fitted linear regression
scatter_light <- function(data, mapping, facetVar, datasetType) {
  ggplot(data, mapping) +
    stat_smooth(
      method = "lm", 
      linewidth = 0, 
      fill = "gray", 
      show.legend = FALSE
    ) + #Standard error plotting
    stat_smooth(
      method = "lm", 
      linetype = 2, 
      se = FALSE
    ) + #Linear regression line
    geom_point(
      size = 0.3, 
      show.legend = FALSE
    ) + #Point for each occupancy value
    faceting(facetVar = facetVar, datasetType = datasetType) +
    ggtheme_light()
} #Make scatter plot with fitted linear regression

scatter_light(
  data = stackedAngles$Pseudohelices,
  mapping = aes(x = Dose, y = Angle, color = Molecule),
  facetVar = 'AngleID',
  datasetType = "Pseudohelix"
)

# Occupancy plotting ------------------------------------------------------

ggplots$Occupancies$Pseudohelices <- scatter_plot() +
  plotcolors() +
  facetet() +
  labs(
    x = "Density-Weighted Dose (MGy)", 
    y = "Occupancy"
  )

ggplots$Occupancies$Pseudohelices <- ggplot(stackedOccupancies$Pseudohelices, aes(x = Dose, y = Occupancy, color = Molecule)) +
  scatter_plot() +
  plotcolors() +
  facetet() +
  labs(
    x = "Density-Weighted Dose (MGy)", 
    y = "Occupancy"
  )

ggplots$Occupancies$Wedges <- ggplot(stackedOccupancies$Wedges, aes(x = WedgeNumber, y = Occupancy, color = Molecule)) +
  scatter_plot() +
  plotcolors() +
  facet_wrap(
    ~Residue, 
    scales = "free", 
    labeller = labeller(
      Residue = c(
        "CO2" = "CO[2]", 
        "Ax"  = "H[2]*O[Ax]", 
        "Eq"  = "H[2]*O[Eq]", 
        "Glu" = "Glu[30]", 
        "Oxy" = "'Dioxygen'"
      ), 
      .default = label_parsed
    )
  ) +
  labs(
    x = "Wedge Number", 
    y = "Occupancy"
  )

# B-factor plotting -------------------------------------------------------

ggplots$BFactors$Pseudohelices <- scatter_plot(stackedBFactors$Pseudohelices, aes(x = Dose, y = bFactor, color = Molecule)) +
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
        "OxyBF" = "'Dioxygen'"
      ), 
      .default = label_parsed
    )
  ) +
  labs(
    x = "Density-Weighted Dose (MGy)", 
    y = bquote('B-Factor ('Å^2*')')
  )

ggplots$BFactors$Wedges <- scatter_plot(stackedBFactors$Wedges, aes(x = WedgeNumber, y = bFactor, color = Molecule)) +
  plotcolors() +
  facet_wrap(
    vars(Residue), 
    scales = "free", 
    labeller = labeller(
      Residue = c(
        "CO2BF" = "CO[2]", 
        "AxBF"  = "H[2]*O[Ax]", 
        "EqBF"  = "H[2]*O[Eq]", 
        "GluBF" = "Glu[30]", 
        "OxyBF" = "'Dioxygen'"
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
    plotcolors() +
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
    plotcolors() +
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
        "Cu-NTerm" = "'Cu-N Terminus'", 
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
        "Cu-NTerm" = "'Cu-N Terminus'", 
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
  plotcolors() +
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
  plotcolors() +
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
