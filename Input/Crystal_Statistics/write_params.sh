#!/usr/bin/env bash

cat > pseudohelix_params.eff <<'EOF'
debug = False
data_labels = None
table_one {
EOF

for i in $(seq 1 36); do
    cat >> pseudohelix_params.eff <<EOF
  structure {
    name = "Pseudohelix_${i}"
    pdb_file = "../../../Files_From_Analysis_Cluster/Pseudohelix_${i}/Refine_1/Phenix_refine_001.pdb"
    mtz_file = "../../../Files_From_Analysis_Cluster/Pseudohelix_${i}/Refine_1/Phenix_refine_001.mtz"
    cif_file = None
    data_labels = None
    r_free_flags_label = None
    source_reflections = *model_file recomputed
    source_r_factors = *model_file recomputed
    source_resolution_bins = *model_file recomputed manual
    wavelength = 1
    high_resolution = None
    n_bins = None
    data_type = *xray neutron electron
    unmerged_data = "./unmerged_data/pseudohelix_${i}_unmerged.mtz"
    unmerged_labels = "I,SIGI"
    use_internal_variance = False
    enable_twinning = True
    twin_law = None
    count_anomalous_pairs_separately = True
    ligand_selection = None
  }
EOF
done

cat >> pseudohelix_params.eff <<'EOF'
  debug = False
  output {
    directory = None
    show_missing_fields = True
    format = *txt *csv *rtf
    base_name = "pseudohelix_statistics"
    verbose = "True"
    text_field_separation = 2
  }
  gui {
    output_dir = None
  }
}
EOF

cat > wedge_params.eff <<'EOF'
debug = False
data_labels = None
table_one {
EOF

for i in $(seq 1 36); do
    cat >> wedge_params.eff <<EOF
  structure {
    name = "Wedge_${i}"
    pdb_file = "../../../Files_From_Analysis_Cluster/Wedge_${i}/Refine_1/Phenix_refine_001.pdb"
    mtz_file = "../../../Files_From_Analysis_Cluster/Wedge_${i}/Refine_1/Phenix_refine_001.mtz"
    cif_file = None
    data_labels = None
    r_free_flags_label = None
    source_reflections = *model_file recomputed
    source_r_factors = *model_file recomputed
    source_resolution_bins = *model_file recomputed manual
    wavelength = 1
    high_resolution = None
    n_bins = None
    data_type = *xray neutron electron
    unmerged_data = "./unmerged_data/HKLOUT_unmerged_Wedge${i}.mtz"
    unmerged_labels = "I,SIGI"
    use_internal_variance = False
    enable_twinning = True
    twin_law = None
    count_anomalous_pairs_separately = True
    ligand_selection = None
  }
EOF
done

cat >> wedge_params.eff <<'EOF'
  debug = False
  output {
    directory = None
    show_missing_fields = True
    format = *txt *csv *rtf
    base_name = "wedge_statistics"
    verbose = "True"
    text_field_separation = 2
  }
  gui {
    output_dir = None
  }
}
EOF