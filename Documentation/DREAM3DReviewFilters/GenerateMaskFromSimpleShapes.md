# GenerateMaskFromSimpleShapes #


## Group (Subgroup) ##

DREAM3DReview (DREAM3DReview)

## Description ##

This **Filter** takes an input array (which could be read in through the **Import ASCII Data**  **Filter**) where each tuple in array corresponds to a masking shape. The user can select between ellipoids, boxes or cylinders. Once a selection is made, every shape in the list is the same, just with different centroids and other relevant dimensions. For example, a user may choose to have a list of 10 boxes they want to use as a mask. Each box is a tuple in the array and would contain information about the center and the box dimensions. Then on the array to be masked, the algorithm checks each point (cell or vertex) to see if that point is within one of the 10 boxes specified in the array. If it is, the cell or vertex value is set to true. Otherwise it is false. 

## Parameters ##

| Name | Type | Description |
|------|------|-------------|
| Mask Shape | int32_t | Which shape to use for mask: ellipsoid, box or cylinder  |


## Required Geometry ##

Image or Vertex

## Required Objects ##

| Kind | Default Name | Type | Component Dimensions | Description |
|------|--------------|-------------|---------|-----|
| **Attribute Array** | Centers | float | 3 | An array with 3 components (x, y, z) where each tuple in the array is the center of the mask shape. Required for all simple shapes|
| **Attribute Array** | Cylinder Radii | float | 1 | An array with 1 component where each tuple in the array is the radius of the cylinder. Required for the cylinder simple shape|
| **Attribute Array** | Cylinder Heights | float | 1 | An array with 1 component where each tuple in the array is the height of the cylinder. Required for the cylinder simple shape |
| **Attribute Array** | Box Dimensions | float | 3 | An array with 3 components (x, y, z) where each tuple in the array is the dimensions of the masking box. Required for the box simple shape.|
| **Attribute Array** | Ellipsoid Axis Lengths | float | 3 | An array with 3 components (x, y, z) where each tuple in the array is the axis length of the shape. Required for the ellipsoid simple shape.|

## Created Objects ##

| Kind | Default Name | Type | Component Dimensions | Description |
|------|--------------|-------------|---------|-----|
| **Attribute Array** | Mask | bool | 1 | An array where if a point is located inside the shapes the value is true, else false|



## Example Pipelines ##

List the names of the example pipelines where this filter is used.

## License & Copyright ##

Please see the description file distributed with this plugin.

## DREAM3D Mailing Lists ##

If you need more help with a filter, please consider asking your question on the DREAM3D Users mailing list:
https://groups.google.com/forum/?hl=en#!forum/dream3d-users