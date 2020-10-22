# GenerateFeatureIDsbyBoundingBoxes #


## Group (Subgroup) ##

DREAM3DReview (DREAM3DReview)

## Description ##

This **Filter** takes an input array (which could be read in through the **Import ASCII Data**  **Filter**) where each tuple in the array corresponds to a bounding box, which is associated with a feature ID. The filter then checks every cell or vertex location in the array to see if it is within a bounding box in the list, and if so, assigns the correpsponding feature ID. 

## Parameters ##

None

## Required Geometry ##

Required Geometry Type -or- Not Applicable

## Required Objects ##

| Kind | Default Name | Type | Component Dimensions | Description |
|------|--------------|-------------|---------|-----|
| **Attribute Array** | Box Corner | float | 3 | An array with 3 components (x, y, z) where each tuple in the array is the origin coordinates a bounding box|
| **Attribute Array** | Box Dimensions | float | 3 | An array with 3 components (x, y, z) where each tuple in the array is the dimensions of the bounding box|
| **Attribute Array** | Feature IDs | int32_t | 1 | An array where each tuple in the array is the feature ID associated with the corresponding box corner and box dimension|

## Created Objects ##

| Kind | Default Name | Type | Component Dimensions | Description |
|------|--------------|-------------|---------|-----|
| **Attribute Array** | Feature IDs | int32_t | 1 | An array of feature IDs assigned based on the bounding box|
| **Attribute Matrix** | Feature Atribute Matrix | Feature Cell or Feature Vertex | N/A | The attribute matrix associated with the feature IDS created by the bounding box |


## Example Pipelines ##

List the names of the example pipelines where this filter is used.

## License & Copyright ##

Please see the description file distributed with this plugin.

## DREAM3D Mailing Lists ##

If you need more help with a filter, please consider asking your question on the DREAM3D Users mailing list:
https://groups.google.com/forum/?hl=en#!forum/dream3d-users