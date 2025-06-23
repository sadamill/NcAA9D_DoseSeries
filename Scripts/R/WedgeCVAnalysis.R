#######################################################
#                    CV CALCULATION                   #
#######################################################

#Compile data frame containing all different wedge variables
allWedges <- tibble(
  cbind(
    trimmedOccupancies[["Wedges"]][2:(ncol(trimmedOccupancies[["Wedges"]])-1)], 
    trimmedBFactors[["Wedges"]][2:(ncol(trimmedBFactors[["Wedges"]])-1)], 
    trimmedDistances[["Wedges"]][2:(ncol(trimmedDistances[["Wedges"]])-1)], 
    allAngles[["Wedges"]][2:(ncol(allAngles[["Wedges"]]))]
  )
)

wedgeCVs <- tibble(
  Parameter = c(
    names(trimmedOccupancies[["Wedges"]][2:(ncol(trimmedOccupancies[["Wedges"]])-1)]), 
    names(trimmedBFactors[["Wedges"]][2:(ncol(trimmedBFactors[["Wedges"]])-1)]), 
    names(trimmedDistances[["Wedges"]][2:(ncol(trimmedDistances[["Wedges"]])-1)]), 
    names(allAngles[["Wedges"]][2:(ncol(allAngles[["Wedges"]])-1)])
  ), 
  ParameterType = c(
    rep("Occupancy", ncol(trimmedOccupancies[["Wedges"]])-2), 
    rep("B-Factor", ncol(trimmedBFactors[["Wedges"]])-2), 
    rep("Distance", ncol(trimmedDistances[["Wedges"]])-2), 
    rep("Angle", ncol(allAngles[["Wedges"]])-2)
  ), 
  MoleculeACV = sapply(subset(allWedges, Molecule == "A")[1:(ncol(allWedges)-1)], function(x) {
    sd(x, na.rm = TRUE)/mean(x, na.rm = TRUE)
  }), 
  MoleculeBCV = sapply(subset(allWedges, Molecule == "B")[1:(ncol(allWedges)-1)], function(x) {
    sd(x)/mean(x)
  })
) %>% 
  cbind(
    ., 
    AverageCV = rowMeans(.[3:4], na.rm = TRUE)
  ) %>% 
  mutate(ParameterType = factor(ParameterType, levels = c("Occupancy", "B-Factor", "Angle", "Distance")))

######################################################
#                    DATA PLOTTING                   #
######################################################

ggplots$Light$CVs$Wedges <- ggplot(wedgeCVs, aes(x = reorder(Parameter, AverageCV), y = AverageCV, fill = ParameterType)) +
  geom_bar(
    stat = "identity", 
    color = "black", 
    position = "dodge2"
  ) +
  ggtheme_light() +
  facet_grid(
    ~ParameterType, 
    scales = "free_x", 
    space = "free_x"
  ) +
  coord_cartesian(
    ylim = c(0, max(wedgeCVs$AverageCV, na.rm = TRUE)*1.05), 
    expand = FALSE
  ) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
    legend.position = "none"
  ) +
  labs(
    x = "Parameter ID", 
    y = "Coefficient of Variance"
  ) +
  scale_x_discrete(labels = c(
    "T1" = bquote(θ[1]), 
    "T2" = bquote(θ[2]), 
    "T3" = bquote(θ[3]), 
    "TH1HN1" = bquote(θH[1]*HN[1]), 
    "TH1HN84" = bquote(θH[1]*HN[84]), 
    "THH" = bquote(θ[HH]), 
    "TT" = bquote(θ[T]), 
    "CuAx" = bquote(Cu-H[2]*O[Ax]), 
    "CuEq" = bquote(Cu-H[2]*O[Eq]), 
    "CuHis1ND1" = bquote(Cu-His[1]*Nδ[1]), 
    "CuHis84NE2" = bquote(Cu-His[84]*Nε[2]), 
    "CuO1" = bquote(Cu-O[prox]), 
    "CuO2" = bquote(Cu-O[dist]), 
    "CuTyr" = bquote(Cu-Tyr[168]), 
    "O1Eq" = bquote(O[prox]*-H[2]*O[Eq]), 
    "O2Eq" = bquote(O[dist]*-H[2]*O[Eq]), 
    "O1GluOE1" = bquote(O[prox]*-Glu[30]*Oε[1]), 
    "O2GluOE1" = bquote(O[dist]*-Glu[30]*Oε[1]), 
    "O1GluOE2" = bquote(O[prox]*-Glu[30]*Oε[2]), 
    "O2GluOE2" = bquote(O[dist]*-Glu[30]*Oε[2]), 
    "O1His157NE2" = bquote(O[prox]*-His[157]*Nε[2]), 
    "O2His157NE2" = bquote(O[dist]*-His[157]*Nε[2]), 
    "CuNterm" = bquote(Cu-N[term]), 
    "Ax" = bquote(H[2]*O[Ax]), 
    "Eq" = bquote(H[2]*O[Eq]), 
    "CO2" = bquote(CO[2]), 
    "Glu" = bquote(Glu[30]), 
    "Oxy" = bquote(Dioxygen), 
    "AxBF" = bquote(H[2]*O[Ax]), 
    "EqBF" = bquote(H[2]*O[Eq]), 
    "CO2BF" = bquote(CO[2]), 
    "GluBF" = bquote(Glu[30]), 
    "OxyBF" = bquote(Dioxygen)
  )) +
  scale_fill_manual(
    values = c("gray95", "gray85", "gray75", "gray65")
  )

ggplots$Dark$CVs$Wedges <- ggplot(wedgeCVs, aes(x = reorder(Parameter, AverageCV), y = AverageCV, fill = ParameterType)) +
  geom_bar(
    stat = "identity", 
    color = "gray", 
    position = "dodge2"
  ) +
  ggtheme_dark() +
  facet_grid(
    ~ParameterType, 
    scales = "free_x", 
    space = "free_x"
  ) +
  coord_cartesian(
    ylim = c(0, max(wedgeCVs$AverageCV, na.rm = TRUE)*1.05), 
    expand = FALSE
  ) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
    legend.position = "none"
  ) +
  labs(
    x = "Parameter ID", 
    y = "Coefficient of Variance"
  ) +
  scale_x_discrete(labels = c(
    "T1" = bquote(θ[1]), 
    "T2" = bquote(θ[2]), 
    "T3" = bquote(θ[3]), 
    "TH1HN1" = bquote(θH[1]*HN[1]), 
    "TH1HN84" = bquote(θH[1]*HN[84]), 
    "THH" = bquote(θ[HH]), 
    "TT" = bquote(θ[T]), 
    "CuAx" = bquote(Cu-H[2]*O[Ax]), 
    "CuEq" = bquote(Cu-H[2]*O[Eq]), 
    "CuHis1ND1" = bquote(Cu-His[1]*Nδ[1]), 
    "CuHis84NE2" = bquote(Cu-His[84]*Nε[2]), 
    "CuO1" = bquote(Cu-O[prox]), 
    "CuO2" = bquote(Cu-O[dist]), 
    "CuTyr" = bquote(Cu-Tyr[168]), 
    "O1Eq" = bquote(O[prox]*-H[2]*O[Eq]), 
    "O2Eq" = bquote(O[dist]*-H[2]*O[Eq]), 
    "O1GluOE1" = bquote(O[prox]*-Glu[30]*Oε[1]), 
    "O2GluOE1" = bquote(O[dist]*-Glu[30]*Oε[1]), 
    "O1GluOE2" = bquote(O[prox]*-Glu[30]*Oε[2]), 
    "O2GluOE2" = bquote(O[dist]*-Glu[30]*Oε[2]), 
    "O1His157NE2" = bquote(O[prox]*-His[157]*Nε[2]), 
    "O2His157NE2" = bquote(O[dist]*-His[157]*Nε[2]), 
    "CuNterm" = bquote(Cu-N[term]), 
    "Ax" = bquote(H[2]*O[Ax]), 
    "Eq" = bquote(H[2]*O[Eq]), 
    "CO2" = bquote(CO[2]), 
    "Glu" = bquote(Glu[30]), 
    "Oxy" = bquote(Dioxygen), 
    "AxBF" = bquote(H[2]*O[Ax]), 
    "EqBF" = bquote(H[2]*O[Eq]), 
    "CO2BF" = bquote(CO[2]), 
    "GluBF" = bquote(Glu[30]), 
    "OxyBF" = bquote(Dioxygen)
  )) +
  scale_fill_manual(
    values = c("gray5", "gray15", "gray25", "gray35")
  )
