for i in $(seq 1 36);
do
echo "phenix.map_value_at_point /Users/sm9/Documents/17_Meilleur_Lab/NcLPMO9D_Dose_Series_Study/Analysis/Files_From_Analysis_Cluster/Pseudohelix_${i}/Refine_1/Phenix_refine_001.mtz /Users/sm9/Documents/17_Meilleur_Lab/NcLPMO9D_Dose_Series_Study/Analysis/Electron_Density_Quantification/QuantificationGrid.pdb miller_array.labels.name=\"2FOFCWT\" > ../Outputs/Pseudohelices/Pseudohelix${i}.txt" >> PseudohelixScript.sh
echo "phenix.map_value_at_point /Users/sm9/Documents/17_Meilleur_Lab/NcLPMO9D_Dose_Series_Study/Analysis/Files_From_Analysis_Cluster/Wedge_${i}/Refine_1/Phenix_refine_001.mtz /Users/sm9/Documents/17_Meilleur_Lab/NcLPMO9D_Dose_Series_Study/Analysis/Electron_Density_Quantification/QuantificationGrid.pdb miller_array.labels.name=\"2FOFCWT\" > ../Outputs/Wedges/Wedge${i}.txt" >> WedgeScript.sh
done
