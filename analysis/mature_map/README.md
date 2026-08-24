# Mature-crown sharpened-map inspection

Input map (not copied into this repository):

`branch6_28_mature_smallest_spherule_8.5A.mrc`

## Map metadata

- Grid: 288 x 288 x 288, float32 (MRC mode 2)
- Voxel size: 2.5 A isotropic
- Box size: 720 A per axis
- Axis order: X, Y, Z (`MAPC/MAPR/MAPS = 1/2/3`)
- Origin: 0, 0, 0
- Actual voxel range: -0.5212 to 0.6088
- Actual mean and standard deviation: 0.00857 and 0.05633
- Header min/max/mean/RMS fields are zero rather than calculated

The crown is centered near X/Y = 143.3 voxels. A direct rotational comparison
of the central crown region gives correlation 0.901 at a 30-degree rotation,
supporting a strong C12 component, but agreement declines for more distant C12
multiples. This is consistent with heterogeneity/asymmetry and means C12 should
not automatically be imposed on every density feature.

## Coarse 8FM9 docking

`8FM9_coarse_fit.pdb` is the unconstrained CA-density optimum. It scores well
numerically but reverses the expected floor/Pol relationship relative to the
putative membrane-facing side. It must not be treated as a biological model.

`8FM9_membrane_consistent_coarse_fit.pdb` uses the alternate Z orientation so
that the floor faces the putative membrane side. Its substantially lower raw
density score shows that the orientation cannot be decided reliably from a
whole-proto-crown rigid-body score in the mature two-ring density. The next
step should fit floor, basal Pol, apical Pol, and legs as separate rigid bodies,
using an experimentally confirmed membrane-side orientation.

All fits here are provisional because only a sharpened full map is available.
There is no half-map cross-validation and no real-space refinement has been
performed.

