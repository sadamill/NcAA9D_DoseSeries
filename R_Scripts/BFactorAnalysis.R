#This script is meant to be run from Master_Script.R

#This script analyzes the variance of B-factors as a function of X-ray dose
#It writes out a PDB file which replaces the B-factor column with the B-factor slope as a function of X-ray dose
#This allows the user to visualize potential impact sites in the structure by coloring the model by B-factor

#####################################################
#            B-FACTOR VARIATION ANALYSIS            #
#####################################################

#Make data frames containing the B-Factor data for each atom across datasets
pseudohelixBFactors <- lapply(pseudohelixAtoms,function(atoms) atoms$b) %>% #Pull B-factors from pseudohelices
  as.data.frame() %>% #Coerce the list into a data frame
  t() #Transpose the data frame to swap rows and columns
row.names(pseudohelixBFactors) <- 1:36 #Set row names to be comprehensible
pseudohelixBFactors <- data.frame(pseudohelixDose,pseudohelixBFactors) #Add column for doses
colnames(pseudohelixBFactors) <- c('dose_MGy',1:5102) #Change column names

wedgeBFactors <- lapply(wedgeAtoms,function(atoms) atoms$b) %>% #Pull B-factors from pseudohelices
  as.data.frame() %>% #Coerce the list into a data frame
  t() #Transpose the data frame to swap rows and columns
row.names(wedgeBFactors) <- 1:36 #Set row names to be comprehensible
wedgeBFactors <- data.frame(wedgeDose,wedgeBFactors) #Add column for doses
colnames(wedgeBFactors) <- c('dose_MGy',1:5102) #Change column names

#Make linear regression vectors
pseudohelixBFSlopes <- apply(pseudohelixBFactors[,-1],2,function(x) { #Exclude the first column (dose) and apply across columns (margin = 2)
  lm(x ~ dose_MGy,data = pseudohelixBFactors)$coefficients[2]  # Extract slope (2nd coefficient)
})

#Clean the linear regression vectors and scale to be visible in PDB B-Factor column (i.e. |Scaled B-Factor Slope| < 100)
pseudohelixBFSlopes[abs(pseudohelixBFSlopes) < 10e-15] <- 0 #Turn any values below 10e-15 to 0
pseudohelixBFSlopes <- pseudohelixBFSlopes * 100 #Scale the numbers by 1000

#Combine the dataframes to provide a singular data frame to pull from
trimmedBFactors <- list(
  Pseudohelices = tibble(
    dose_MGy = rep(pseudohelixBFactors$dose_MGy,2),
    CO2BF = c(
      rep(NA,36),
      pseudohelixBFactors$`3983`
    ),
    AxBF = c(
      pseudohelixBFactors$`4903`,
      pseudohelixBFactors$`4715`
    ),
    EqBF = c(
      pseudohelixBFactors$`4711`,
      pseudohelixBFactors$`4580`
    ),
    OxyBF = c(
      pseudohelixBFactors$`5099`,
      pseudohelixBFactors$`5101`
    ),
    GluBF = c(
      pseudohelixBFactors$`266`,
      pseudohelixBFactors$`2195`
    ),
    Molecule = c(
      rep('A',36),
      rep('B',36)
    )
  ),
  Wedges = tibble(
    WedgeNumber = rep(1:36,2),
    CO2BF = c(
      rep(NA,36),
      wedgeBFactors$`3983`
    ),
    AxBF = c(
      wedgeBFactors$`4903`,
      wedgeBFactors$`4715`
    ),
    EqBF = c(
      wedgeBFactors$`4711`,
      wedgeBFactors$`4580`
    ),
    OxyBF = c(
      wedgeBFactors$`5099`,
      wedgeBFactors$`5101`
    ),
    GluBF = c(
      wedgeBFactors$`266`,
      wedgeBFactors$`2195`
    ),
    Molecule = c(
      rep('A',36),
      rep('B',36)
    )
  )
)

#Make a data frame containing occupancies, but stacked so they can be properly faceted in ggplot
stackedBFactors <- list(
  Pseudohelices = tibble(
    Dose = rep(trimmedBFactors$Pseudohelices$dose_MGy,5),
    bFactor = c(
      trimmedBFactors$Pseudohelices$CO2BF,
      trimmedBFactors$Pseudohelices$AxBF,
      trimmedBFactors$Pseudohelices$EqBF,
      trimmedBFactors$Pseudohelices$OxyBF,
      trimmedBFactors$Pseudohelices$GluBF
    ),
    Residue = c(
      rep('CO2BF',72),
      rep('AxBF',72),
      rep('EqBF',72),
      rep('OxyBF',72),
      rep('GluBF',72)
    ),
    Molecule = rep(trimmedBFactors$Pseudohelices$Molecule,5),
  ) %>% 
    mutate(Residue = factor(Residue, levels = c('OxyBF','GluBF','CO2BF','AxBF','EqBF'))),
  Wedges = tibble(
    WedgeNumber = rep(trimmedBFactors$Wedges$WedgeNumber,5),
    bFactor = c(
      trimmedBFactors$Wedges$CO2BF,
      trimmedBFactors$Wedges$AxBF,
      trimmedBFactors$Wedges$EqBF,
      trimmedBFactors$Wedges$OxyBF,
      trimmedBFactors$Wedges$GluBF
    ),
    Residue = c(
      rep('CO2BF',72),
      rep('AxBF',72),
      rep('EqBF',72),
      rep('OxyBF',72),
      rep('GluBF',72)
    ),
    Molecule = rep(trimmedBFactors$Wedges$Molecule,5),
  ) %>% 
    mutate(Residue = factor(Residue, levels = c('OxyBF','GluBF','CO2BF','AxBF','EqBF')))
)

#Make vector containing occupancy change of each atom
pseudohelixBFSlopes <- apply(pseudohelixBFactors[,2:(ncol(pseudohelixBFactors))],2,function(x) { #Exclude the first column (dose) and apply across columns (margin = 2)
  lm(x ~ dose_MGy,data = pseudohelixBFactors)$coefficients[2]  # Extract slope (2nd coefficient)
})

#Clean the slope vectors
pseudohelixBFSlopes[abs(pseudohelixBFSlopes) < 10e-15] <- 0 #Turn any values below 10e-15 to 0
pseudohelixBFSlopes <- pseudohelixBFSlopes * 1000 #Scale the numbers by 1000

##############################################################
#                 LINEAR REGRESSION ANALYSIS                 #
##############################################################

#Prepare a list of multiple linear regression models for each atom of interest
multipleRegressions$BFactors <- list(
  Pseudohelices = list(
    OxyBF = lm(OxyBF ~ dose_MGy * Molecule, data = trimmedBFactors$Pseudohelices),
    GluBF = lm(GluBF ~ dose_MGy * Molecule, data = trimmedBFactors$Pseudohelices),
    CO2BF = lm(CO2BF ~ dose_MGy, data = subset(trimmedBFactors$Pseudohelices, Molecule == 'B')),
    AxBF = lm(AxBF ~ dose_MGy * Molecule, data = trimmedBFactors$Pseudohelices),
    EqBF = lm(EqBF ~ dose_MGy * Molecule, data = trimmedBFactors$Pseudohelices)
  ),
  Wedges = list(
    OxyBF = lm(OxyBF ~ WedgeNumber * Molecule, data = trimmedBFactors$Wedges),
    GluBF = lm(GluBF ~ WedgeNumber * Molecule, data = trimmedBFactors$Wedges),
    CO2BF = lm(CO2BF ~ WedgeNumber, data = subset(trimmedBFactors$Wedges, Molecule == 'B')),
    AxBF = lm(AxBF ~ WedgeNumber * Molecule, data = trimmedBFactors$Wedges),
    EqBF = lm(EqBF ~ WedgeNumber * Molecule, data = trimmedBFactors$Wedges)
  )
)

#Prepare summary table of linear regression analysis
regressionSummaries$BFactors <- list(
  Pseudohelices = tibble(
    Residue = c(
      rep('OxyBF',3),
      rep('GluBF',3),
      rep('CO2BF',3),
      rep('AxBF',3),
      rep('EqBF',3)
    ),
    Estimate = rep(c('TrendA','TrendB','Contrast'),5),
    Coefficient = c(
      emtrends.coefficient(multipleRegressions$BFactors$Pseudohelices$OxyBF, 'dose_MGy'),
      emtrends.coefficient(multipleRegressions$BFactors$Pseudohelices$GluBF, 'dose_MGy'),
      NA,
      summary(multipleRegressions$BFactors$Pseudohelices$CO2BF)$coefficients[2,2],
      NA,
      emtrends.coefficient(multipleRegressions$BFactors$Pseudohelices$AxBF, 'dose_MGy'),
      emtrends.coefficient(multipleRegressions$BFactors$Pseudohelices$EqBF, 'dose_MGy')
    ),
    StandardError = c(
      emtrends.se(multipleRegressions$BFactors$Pseudohelices$OxyBF, 'dose_MGy'),
      emtrends.se(multipleRegressions$BFactors$Pseudohelices$GluBF, 'dose_MGy'),
      NA,
      summary(multipleRegressions$BFactors$Pseudohelices$CO2BF)$coefficients[2,4],
      NA,
      emtrends.se(multipleRegressions$BFactors$Pseudohelices$AxBF, 'dose_MGy'),
      emtrends.se(multipleRegressions$BFactors$Pseudohelices$EqBF, 'dose_MGy')
    ),
    ModelRSquared = c(
      rep(summary(multipleRegressions$BFactors$Pseudohelices$Oxy)$r.squared, 3),
      rep(summary(multipleRegressions$BFactors$Pseudohelices$Glu)$r.squared, 3),
      NA,
      summary(multipleRegressions$BFactors$Pseudohelices$CO2)$r.squared,
      NA,
      rep(summary(multipleRegressions$BFactors$Pseudohelices$Ax)$r.squared, 3),
      rep(summary(multipleRegressions$BFactors$Pseudohelices$Eq)$r.squared, 3)
    ),
    PValue = c(
      emtrends.pvalue(multipleRegressions$BFactors$Pseudohelices$OxyBF, 'dose_MGy'),
      emtrends.pvalue(multipleRegressions$BFactors$Pseudohelices$GluBF, 'dose_MGy'),
      NA,
      summary(multipleRegressions$BFactors$Pseudohelices$CO2BF)$coefficients[2,4],
      NA,
      emtrends.pvalue(multipleRegressions$BFactors$Pseudohelices$AxBF, 'dose_MGy'),
      emtrends.pvalue(multipleRegressions$BFactors$Pseudohelices$EqBF, 'dose_MGy')
    )
  ) %>% mutate(Significance = ifelse(
    PValue <= 0.001, '***',
    ifelse(
      PValue <= 0.01, '**',
      ifelse(
        PValue <= 0.05, '*', ' '
      )
    )
  )),
  Wedges = tibble(
    Residue = c(
      rep('OxyBF',3),
      rep('GluBF',3),
      rep('CO2BF',3),
      rep('AxBF',3),
      rep('EqBF',3)
    ),
    Estimate = rep(c('TrendA','TrendB','Contrast'),5),
    Coefficient = c(
      emtrends.coefficient(multipleRegressions$BFactors$Wedges$OxyBF, 'WedgeNumber'),
      emtrends.coefficient(multipleRegressions$BFactors$Wedges$GluBF, 'WedgeNumber'),
      NA,
      summary(multipleRegressions$BFactors$Wedges$CO2BF)$coefficients[2,2],
      NA,
      emtrends.coefficient(multipleRegressions$BFactors$Wedges$AxBF, 'WedgeNumber'),
      emtrends.coefficient(multipleRegressions$BFactors$Wedges$EqBF, 'WedgeNumber')
    ),
    StandardError = c(
      emtrends.se(multipleRegressions$BFactors$Wedges$OxyBF, 'WedgeNumber'),
      emtrends.se(multipleRegressions$BFactors$Wedges$GluBF, 'WedgeNumber'),
      NA,
      summary(multipleRegressions$BFactors$Wedges$CO2BF)$coefficients[2,4],
      NA,
      emtrends.se(multipleRegressions$BFactors$Wedges$AxBF, 'WedgeNumber'),
      emtrends.se(multipleRegressions$BFactors$Wedges$EqBF, 'WedgeNumber')
    ),
    ModelRSquared = c(
      rep(summary(multipleRegressions$BFactors$Wedges$Oxy)$r.squared, 3),
      rep(summary(multipleRegressions$BFactors$Wedges$Glu)$r.squared, 3),
      NA,
      summary(multipleRegressions$BFactors$Wedges$CO2)$r.squared,
      NA,
      rep(summary(multipleRegressions$BFactors$Wedges$Ax)$r.squared, 3),
      rep(summary(multipleRegressions$BFactors$Wedges$Eq)$r.squared, 3)
    ),
    PValue = c(
      emtrends.pvalue(multipleRegressions$BFactors$Wedges$OxyBF, 'WedgeNumber'),
      emtrends.pvalue(multipleRegressions$BFactors$Wedges$GluBF, 'WedgeNumber'),
      NA,
      summary(multipleRegressions$BFactors$Wedges$CO2BF)$coefficients[2,4],
      NA,
      emtrends.pvalue(multipleRegressions$BFactors$Wedges$AxBF, 'WedgeNumber'),
      emtrends.pvalue(multipleRegressions$BFactors$Wedges$EqBF, 'WedgeNumber')
    )
  ) %>% mutate(Significance = ifelse(
    PValue <= 0.001, '***',
    ifelse(
      PValue <= 0.01, '**',
      ifelse(
        PValue <= 0.05, '*', ' '
      )
    )
  ))
)

######################################################
#                 DATA VISUALIZATION                 #
######################################################

#Write a model for B-factor coloring in ChimerAxBF
bFactorColoredPDB <- pseudohelixList[[1]]
bFactorColoredPDB$atom$b <- pseudohelixBFSlopes

ggplots$BFactors$Pseudohelices <- ggplot(stackedBFactors$Pseudohelices,aes(x = Dose, y = bFactor, color = Molecule)) +
  stat_smooth( #Standard error plotting
    method = 'lm',
    linewidth = 0,
    fill = 'gray85',
    show.legend = FALSE
  ) +
  stat_smooth( #Linear regression line
    method = 'lm',
    linetype = 2,
    se = FALSE
  ) +
  geom_point(
    size = 0.3,
    show.legend = FALSE
  ) + #Point for each occupancy value
  ggtheme() +
  facet_wrap(
    vars(Residue),
    scales = 'free',
    labeller = labeller(
      Residue = c(
        'CO2BF' = 'CO[2]',
        'AxBF'  = 'H[2]*O[Ax]',
        'EqBF'  = 'H[2]*O[Eq]',
        'GluBF' = 'Glu[30]',
        'OxyBF' = '"Dioxygen"'
      ),
      .default = label_parsed
    )
  ) +
  labs(
    x = 'Density-Weighted Dose (MGy)',
    y = 'Occupancy'
  )

ggplots$BFactors$Wedges <- ggplot(stackedBFactors$Wedges,aes(x = WedgeNumber, y = bFactor, color = Molecule)) +
  stat_smooth( #Standard error plotting
    method = 'lm',
    linewidth = 0,
    fill = 'gray85',
    show.legend = FALSE
  ) +
  stat_smooth( #Linear regression line
    method = 'lm',
    linetype = 2,
    se = FALSE
  ) +
  geom_point(
    size = 0.3,
    show.legend = FALSE
  ) + #Point for each occupancy value
  ggtheme() +
  facet_wrap(
    vars(Residue),
    scales = 'free',
    labeller = labeller(
      Residue = c(
        'CO2BF' = 'CO[2]',
        'AxBF'  = 'H[2]*O[Ax]',
        'EqBF'  = 'H[2]*O[Eq]',
        'GluBF' = 'Glu[30]',
        'OxyBF' = '"Dioxygen"'
      ),
      .default = label_parsed
    )
  ) +
  labs(
    x = 'Wedge Number',
    y = 'Occupancy'
  )
