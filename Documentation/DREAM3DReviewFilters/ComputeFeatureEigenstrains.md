# Compute Eigenstrains by Feature (Grain/Inclusion) #

## Group (Subgroup) ##

DREAM3DReview (DREAM3DReview)

## Description ##

This **Filter** calculates the eigenstrain tensor of each **Feature** by using Poisson's ratio, **Feature** shape information, and elastic strain tensor information. This filter was developed alongside Ref. [1]. 

For each **Feature**, the Eshelby tensor is calculated using the solution for Eshelby's isotropic inclusion given by Mura [2]. If *Use Ellipsoidal Grains* is enabled, the elliptic integrals used to calculate the Eshelby tensor are numerically approximated using 32-point Gaussian quadrature. Otherwise, the tensor is exact when grains are assumed to be spherical. The eigenstrain tensor is then calculated using the known relation between Eshelby's tensor and the elastic strain tensor [2]. The correctional matrix can be used to emperically correct the eigenstrains using the approaches from Refs. [3] and [4].

## Parameters ##

| Name | Type | Description |
|------|------|------|
| Poisson's Ratio | float | Poisson's ratio (normally calculated using a self-consistent approximation) |
| Use Ellipsoidal Grains | bool | Whether or not to use the full ellipsoid Eshelby tensor solution (requires grain shape information) |
| Use Correctional Matrix | bool | Whether or not to use the correctional matrix for each strain component |
| Beta11 | float | Correctional value for the 11/xx component |
| Beta22 | float | Correctional value for the 22/yy component |
| Beta33 | float | Correctional value for the 33/zz component |
| Beta23 | float | Correctional value for the 23/yz component |
| Beta13 | float | Correctional value for the 13/xz component |
| Beta12 | float | Correctional value for the 12/xy component |

## Required Geometry ##

Image

## Required Objects ##

| Kind | Default Name | Type | Component Dimensions | Description |
|------|--------------|-------------|---------|-----|
| **Feature Attribute Array** | AxisLengths | float | (3) | Semi-axis lengths (a, b, c) for best-fit ellipsoid to **Feature** |
| **Feature Attribute Array** | AxisEulerAngles | float | (3) | Euler angles (in radians) necessary to rotate the sample reference frame to the reference frame of the **Feature**, where the prinicpal axes of the best-fit ellipsoid are (X, Y, Z) |
| **Feature Attribute Array** | ElasticStrains | float | (6) | Elastic strain tensor of the **Feature** in Voigt notation |

## Created Objects ##

| Kind | Default Name | Type | Component Dimensions | Description |
|------|--------------|-------------|---------|-----|
| **Feature Attribute Array** | Eigenstrains | float | (6) | Eigenstrain tensor of the **Feature** in Voigt notation |

## References ## 

[1] The AFRL AM Modeling Challenge: Predicting Micromechanical Fields in AM IN625 Using an FFT-Based Method with Direct Input from a 3D Microstructural Image, C.K. Cocke, A.D. Rollett, R.A. Lebensohn, A.D. Spear. In preperation.

[2] Micromechanics of Defects in Solids, Mura, T., Martinus-Nijhoff, Dordrecht, 1987

[3] Instantiation of crystal plasticity simulations for micromechanical modelling with direct input from microstructural data collected at light sources, R. Pokharel and R.A Lebensohn, Scripta Mater., 132 (2017), pp. 73-77

[4] Validation of micro-mechanical FFT-based simulations using High Energy Diffraction Microscopy on Ti-7Al, V. Tari, R.A. Lebensohn, R. Pokharel, T.J. Turner, P.A. Shade, J.V. Bernier, A.D. Rollet, Acta Mater., 154 (2018), pp. 273-283

## Example Pipelines ##

List the names of the example pipelines where this filter is used.

## License & Copyright ##

Please see the description file distributed with this plugin.

## DREAM3D Mailing Lists ##

If you need more help with a filter, please consider asking your question on the DREAM3D Users mailing list:
https://groups.google.com/forum/?hl=en#!forum/dream3d-users