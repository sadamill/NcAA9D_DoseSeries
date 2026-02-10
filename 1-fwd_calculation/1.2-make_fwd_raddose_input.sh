#!/usr/bin/env bash

# start every file with a crystal and beam block
for wedge_number in {1..38}
do
	cat <<-EOF > 1-fwd_calculation/input/raddose_input/wedge${wedge_number}.txt
	##############################################################################
	#                                 Crystal Block                              #
	##############################################################################

	Crystal

	Type Cuboid
	    # Crystal shape can be Cuboid or Spherical
	
	WireframeType Obj
	    # Specify file format of crystal input model (obj is currently the only option)

	ModelFile ../../../input/raddose_input/CrystalWireframe.obj
	    # Specify filepath to crystal input model

	Dimensions 250 1187 250
	    # Dimensions of the crystal in X,Y,Z in µm.
	    # Z is the beam axis, Y the rotation axis and
	    # X completes the right handed set
	    # (vertical if starting face-on).

	PixelsPerMicron 0.10000000001
	    # This defines the coarseness of the simulation
	    # (i.e. how many voxels the crystal is divided into.)
	    # Preferably set as high as possible, however for a higher
	    # value the simulation will take longer to complete.
	    # Recommended to try increasing between 0.5 and 5 and ensure
	    # the reported dose value converges as PixelsPerMicron increases.
	    # As a rule of thumb, this needs to be at least 10x the beam
	    # FWHM for a Gaussian beam.
	    # e.g. 20 µm FWHM beam -> 2 µm voxels -> 0.5 voxels/µm

	# NOTE: Use AngleP/AngleL if your crystal is not face-on to the beam.
	# See RD3D user guide for more details

	AbsCoefCalc RD3D
	    # Absorption Coefficients calculated
	    # using RADDOSE-3D (Zeldin et al. 2013).

	UnitCell 68.12  42.23  70.29  90.0  98.3  90.0
	    # Unit cell size: a, b, c with α, β and γ angles default to 90°

	NumMonomers 4
	    # Number of monomers in unit cell

	NumResidues 223
	    # Number of residues per monomer

	NumCarb 2
	    # Number of carbohydrate residues per monomer

	ProteinHeavyAtoms Cu 1 S 9
	    # Heavy atoms added to protein part of the
	    # monomer, i.e. S, coordinated metals, Se in Se-Met

	SolventHeavyConc S 100
	    # Concentration of elements in the solvent
	    # in mmol/L. Oxygen and lighter elements
	    # should not be specified

	SolventFraction 0.37
	    # Fraction of the unit cell occupied by solvent

	##############################################################################
	#                                  Beam Block                                #
	##############################################################################

	Beam

	Type Gaussian
	    # Beam profile can be Gaussian or TopHat

	Flux 0.49e12
	    # in photons per second (2e12 = 2 * 10^12)

	FWHM 40 30
	    # in µm, horizontal by vertical for a Gaussian beam

	Energy 12.4
	    # Photon energy in keV

	Collimation Rectangular 20 20
	    # Horizontal/Vertical collimation of the beam.
	    # For 'uncollimated' Gaussians, 3xFWHM recommended
	EOF
done

# Append proper wedge blocks to their respective files

for subwedge in $(seq 1 38); do
  start=$(( (subwedge - 1) * 5 ))
  end=$(( start + 180 ))
  position=$(( subwedge * -30 + 570 ))
  angular_range=180
  frames=180
  angular_resolution=$(printf "%.10g" "$(echo "$angular_range / 200 - 0.0000000001" | bc -l)")

  block=$(printf "################################################################################
#                                  Wedge Block %-2d                              #
################################################################################

Wedge $start $end
    # Start and End rotational angle of the crystal with Start < End

ExposureTime $frames
    # Total time for entire angular range

StartOffset 0 $position 0
    # Offset x y z translation in μm relative to crystal origin
    # Origin is defined as intersection of beam and aligned goniometer axis

AngularResolution $angular_resolution
    # Only change from the defaults when using very
    # small wedges, e.g 5°.

# NOTE: To define more complex geometries (helical, de-centered, or offset),
# see the StartOffset, TranslatePerDegree, and RotAxBeamOffset keywords
# in the User Guide
" "$subwedge")

  for wedge_number in $(seq ${subwedge} 38); do
    echo "$block" >> 1-fwd_calculation/input/raddose_input/wedge${wedge_number}.txt
  done
done
