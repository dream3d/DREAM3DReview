Find Mask Neighbors 
=============

## Group (Subgroup) ##

DREAM3D Review (Spatial)

## Description ##

This **Filter** determines whether a **Feature** touches a masked **Cell**, where the mask is supplied by the user. (FIXME) This is accomplished by simply querying the **Feature** owners of the **Cells** that sit at either. Any **Feature** that owns one of those **Cells** is said to touch an outer surface and all other **Features** are said to not touch an outer surface of the sample.

This **Filter** determines whether a **Feature** touches a masked **Cell**. A **Feature** is considered touching a masked **Cell** if the following condition is met:

+ Any cell has **Mask = True** as a neighbor.

The output of this filter is a **Feature** level array of booleans where 0=Not Touching Masked Data and 1=Touching Masked Data.

### 2D Image Geometry ###

If the structure/data is actually 2D, then the dimension that is planar is not considered and only the **Features** touching the edges are considered surface **Features**.


## Parameters ##

None

## Required Geometry ##

Image

## Required Objects ##

| Kind | Default Name | Type | Component Dimensions | Description |
|------|--------------|------|----------------------|-------------|
| **Cell Attribute Array** | FeatureIds | int32_t | (1) | Specifies to which **Feature** each **Cell** belongs |
| **Cell Attribute Array** | Mask | Bool | (1) | True for cells user wants to highlight and False for other cells | 

## Created Objects ##

| Kind | Default Name | Type | Component Dimensions | Description |
|------|--------------|------|----------------------|-------------|
| **Feature Attribute Array** | MaskNeighbors | bool | (1) | Flag of 1 if **Feature** touches a masked cell or of 0 if it does not |

## Example Pipelines ##



## License & Copyright ##

Please see the description file distributed with this **Plugin**

## DREAM.3D Mailing Lists ##

If you need more help with a **Filter**, please consider asking your question on the [DREAM.3D Users Google group!](https://groups.google.com/forum/?hl=en#!forum/dream3d-users)


