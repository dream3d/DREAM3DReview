# Export MASSIF Data (HDF5)  #

## Group (Subgroup) ##

DREAM3DReview (DREAM3DReview)

## Description ##

This **Filter** outputs an HDF5 file for use with the Micromechanical Analysis of Stress-strain Inhomogeneities with fast Fourier transforms (MASSIF) code. Optionally outputs an eigenstrain field file for incorporation of residual strain information.

## Parameters ##

| Name | Type | Description |
|------|------|-------------|
| Output File | string | Output file path of the main MASSIF input file |
| Write Eigenstrains | bool | Whether or not to output the eigenstrain input file |
| Eigenstrain Output File | string | Output file path of the MASSIF eigenstrain input file |

## Required Geometry ###

Image

## Required Objects ##

| Kind | Default Name | Type | Component Dimensions | Description |
|------|--------------|------|----------------------|-------------|
| **Cell Attribute Array** | FeatureIds | int32_t | (1) | Specifies to which **Feature** each **Cell** belongs |
| **Cell Attribute Array** | Euler Angles | float | (3) | Specifies the crystal orientation each **Cell** |
| **Cell Attribute Array** | Phases | int32_t | (1) | Specifies the phase of each **Cell** |
| **Cell Attribute Array** | Eigenstrains | float | (6) | Specifies the eigenstrain tensor of each **Cell** in Voigt notation |

## Created Objects ##

None

## Example Pipelines ##
 
ComputeFeatureEigenstrainsExample

## License & Copyright ##

Please see the description file distributed with this plugin.

## DREAM3D Mailing Lists ##

If you need more help with a filter, please consider asking your question on the DREAM3D Users mailing list:
