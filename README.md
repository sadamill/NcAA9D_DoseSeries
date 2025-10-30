# *Nc*AA9D Dose Series Analysis Scripts

The scripts in this repository were used to analyze *Nc*AA9D structures for
["PAPER TITLE HERE"](LINK HERE). A total of 72 protein crystal structures were
determined from a single crystal of *Nc*AA9D and used to determine the effects
of radiation damage on the structure.

The analysis pipeline built for this analysis imports and parses the .pdb files, 
storing atom parameters in a large data frame. Multiple linear regression models
are then calculated for parameters of interest versus dose and chain ID (see
the publication for more details). The script also generates plots which were
used in the paper.
