# Function setup ----------------------------------------------------------

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
  mapping_list <- list(
    Residue = c(
      "CO2" = "CO[2]", 
      "Ax"  = "H[2]*O[Ax]", 
      "Eq"  = "H[2]*O[Eq]", 
      "Glu" = "Glu[30]", 
      "Oxy" = "'Dioxygen'"
    ),
    bResidue = c(
      "CO2BF" = "CO[2]", 
      "AxBF"  = "H[2]*O[Ax]", 
      "EqBF"  = "H[2]*O[Eq]", 
      "GluBF" = "Glu[30]", 
      "OxyBF" = "'Dioxygen'"
    ),
    CuAtomPair = c(
      "Cu-Tyr" = "Cu-Tyr[168]",
      "Cu-NTerm" = "Cu-N[term]",
      "Cu-His1ND1" = "Cu-His[1]*Nδ[1]",
      "Cu-His84NE2" = "Cu-His[84]*Nε[2]",
      "Cu-Eq" = "Cu-H[2]*O[Eq]",
      "Cu-Ax" = "Cu-H[2]*O[Ax]"
    ),
    OAtomPair = c(
      "O1-Eq" = "O[prox]*-H[2]*O[Eq]", 
      "O2-Eq" = "O[dist]*-H[2]*O[Eq]", 
      "O1-His157NE2" = "O[prox]*-His[157]*N[ε2]", 
      "O2-His157NE2" = "O[dist]*-His[157]*N[ε2]", 
      "O1-GluOE1" = "O[prox]*-Glu[30]*O[ε1]", 
      "O2-GluOE1" = "O[dist]*-Glu[30]*O[ε1]", 
      "O1-GluOE2" = "O[prox]*-Glu[30]*O[ε2]", 
      "O2-GluOE2" = "O[dist]*-Glu[30]*O[ε2]"
    ),
    AngleID = c(
      "T1" = "θ[1]", 
      "T2" = "θ[2]", 
      "T3" = "θ[3]", 
      "TT" = "θ[T]", 
      "THH" = "θ[HH]", 
      "TH1HN1" = "θH[1]*HN[1]", 
      "TH1HN84" = "θH[1]*HN[84]", 
      "TOxy" = "θ[oxygen]"
    )
  )
  
  list(
    facet_wrap(
      as.formula(paste("~", facetVar)), 
      scales = 'free',
      labeller = labeller(
        !!facetVar := as_labeller(mapping_list[[facetVar]], label_parsed),
        .default = label_parsed
      )),
    labs(
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
}
scales <- function(facetVar) {
  
  angle_window <- 7
  occupancy_window <- 0.4
  bfactor_window <- 8
  cudistance_window <- 0.25
  
  facetted_pos_scales(y = if(facetVar == "AngleID") {
    list(
      scale_y_continuous(limits = c(94 - angle_window / 2, 94 + angle_window / 2)),
      scale_y_continuous(limits = c(94 - angle_window / 2, 94 + angle_window / 2)),
      scale_y_continuous(limits = c(171 - angle_window / 2, 171 + angle_window / 2)),
      scale_y_continuous(limits = c(180 - angle_window / 2, 180 + angle_window / 2)),
      scale_y_continuous(limits = c(176 - angle_window / 2, 176 + angle_window / 2)),
      scale_y_continuous(limits = c(68 - angle_window / 2, 68 + angle_window / 2)),
      scale_y_continuous(limits = c(120, 150)),
      scale_y_continuous(limits = c(5 - angle_window / 2, 5 + angle_window / 2))
    )} else if(facetVar == "Residue") {
      list(
        scale_y_continuous(limits = c(0.6 - occupancy_window / 2, 0.6 + occupancy_window / 2)),
        scale_y_continuous(limits = c(0.55 - occupancy_window / 2, 0.55 + occupancy_window / 2)),
        scale_y_continuous(limits = c(0.35 - occupancy_window / 2, 0.35 + occupancy_window / 2)),
        scale_y_continuous(limits = c(0.55 - occupancy_window / 2, 0.55 + occupancy_window / 2)),
        scale_y_continuous(limits = c(0.65 - occupancy_window / 2, 0.65 + occupancy_window / 2))
      )} else if(facetVar == "bResidue") {
        list(
          scale_y_continuous(limits = c(21 - bfactor_window / 2, 21 + bfactor_window / 2)),
          scale_y_continuous(limits = c(10.5 - bfactor_window / 2, 10.5 + bfactor_window / 2)),
          scale_y_continuous(limits = c(22.6 - bfactor_window / 2, 22.6 + bfactor_window / 2)),
          scale_y_continuous(limits = c(14 - bfactor_window / 2, 14 + bfactor_window / 2)),
          scale_y_continuous(limits = c(12.5 - bfactor_window / 2, 12.5 + bfactor_window / 2))
        )
      } else if(facetVar == "CuAtomPair") {
        list(
          scale_y_continuous(limits = c(2.375 - cudistance_window / 2, 2.375 + cudistance_window / 2)),
          scale_y_continuous(limits = c(1.99 - cudistance_window / 2, 1.99 + cudistance_window / 2)),
          scale_y_continuous(limits = c(1.95 - cudistance_window / 2, 1.95 + cudistance_window / 2)),
          scale_y_continuous(limits = c(1.98 - cudistance_window / 2, 1.98 + cudistance_window / 2)),
          scale_y_continuous(limits = c(2.18 - cudistance_window / 2, 2.18 + cudistance_window / 2)),
          scale_y_continuous(limits = c(2.66 - cudistance_window / 2, 2.66 + cudistance_window / 2))
        )
      } else {NULL})
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
    ggtheme_dark() +
    scales(facetVar = facetVar) +
    if(facetVar == "AngleID" | facetVar == "OAtomPair") {
      theme(legend.position.inside = c(0.85, 0.15))
    } else if(facetVar == "CuAtomPair") {
      theme(legend.position = "right")
    }
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
    ggtheme_light() +
    scales(facetVar = facetVar) +
    if(facetVar == "AngleID" | facetVar == "OAtomPair") {
      theme(legend.position.inside = c(0.85, 0.15))
    } else if(facetVar == "CuAtomPair") {
      theme(legend.position = "right")
    }
} #Make scatter plot with fitted linear regression

ggdarklight <- function(data, key, key2 = NA) {
  yvar <- if(deparse(substitute(data)) == "stackedOccupancies") {
    'Occupancy'
  } else if(deparse(substitute(data)) == "stackedBFactors") {
    'bFactor'
  } else if(deparse(substitute(data)) %in% c("stackedDistances$Cu", "stackedDistances$Oxy")) {
    'Distance'
  } else if(deparse(substitute(data)) == "stackedAngles") {
    'Angle'
  } else {stop("Invalid data input")}
  
  facetStr <- if(deparse(substitute(data)) == "stackedOccupancies") {
    "Residue"
  } else if(deparse(substitute(data)) == "stackedBFactors") {
    "bResidue"
  } else if(deparse(substitute(data)) == "stackedDistances$Cu") {
    "CuAtomPair"
  } else if(deparse(substitute(data)) == "stackedDistances$Oxy") {
    "OAtomPair"
  } else if(deparse(substitute(data)) == "stackedAngles") {
    "AngleID"
  } else{stop("Invalid data input")}
  
  plots <- list(
    Light = list(
      Pseudohelices = scatter_light(
        data = data$Pseudohelices,
        mapping = aes(
          x = Dose,
          y = .data[[yvar]],
          color = Molecule
        ),
        facetVar = facetStr,
        datasetType = "Pseudohelix"
      ),
      Wedges = scatter_light(
        data = data$Wedges,
        mapping = aes(
          x = WedgeNumber,
          y = .data[[yvar]],
          color = Molecule
        ),
        facetVar = facetStr,
        datasetType = "Pseudohelix"
      )
    ),
    Dark = list(
      Pseudohelices = scatter_dark(
        data = data$Pseudohelices,
        mapping = aes(
          x = Dose,
          y = .data[[yvar]],
          color = Molecule
        ),
        facetVar = facetStr,
        datasetType = "Pseudohelix"
      ),
      Wedges = scatter_dark(
        data = data$Wedges,
        mapping = aes(
          x = WedgeNumber,
          y = .data[[yvar]],
          color = Molecule
        ),
        facetVar = facetStr,
        datasetType = "Pseudohelix"
      )
    )
  )
  
  if(is.na(key2)) {
    ggplots[[key]] <<- plots
  } else {
    ggplots[[key]][[key2]] <<- plots
  }

} #Combination of scatter_dark and scatter_light to ease the plotting of multiple datasets with multiple themes

# Occupancy plotting ------------------------------------------------------

ggdarklight(stackedOccupancies, "Occupancies")

# B-Factor plotting -------------------------------------------------------

ggdarklight(stackedBFactors, "BFactors")

# Distance plotting -------------------------------------------------------

ggdarklight(stackedDistances$Cu, "Distances", "Cu")
ggdarklight(stackedDistances$Oxy, "Distances", "Oxy")

# Angle plotting ----------------------------------------------------------

ggdarklight(stackedAngles, "Angles")
