# Normalize Grayscale Values [Python] #

## Group (Subgroup) ##

Processing(Normalization)

## Description ##

This **Filter** normalizes brightness and contrast values in a 3D grayscale data array. The **Filter** gets the averages each layer in the x, y, or z direction.

## Parameters ##

| Name | Type | Description |
|------|------|------|
| Slice Axis | int | The direction for each layer |
| Create New Array | boolean | If checked, a new array is created with the corrected values. Else, the inputted array is overwritten. |


## Required Geometry ##

Image

## Required Objects ##

| Kind | Default Name | Type | Component Dimensions | Description |
|------|------|------|------|------|
| **Cell Data Array** | AttributeArray Name | Any | (3) | The selected data array to apply normalization |

## Created Objects ##

| Kind | Default Name | Type | Component Dimensions | Description |
|------|------|------|------|------|
 | **Cell Data Array** |AttributeArray Name | same as Required Data Array | (3) | Optional array created to store normalized data |


## Example Pipelines ##

List the names of the example pipelines where this filter is used.

## License & Copyright ##

Please see the description file distributed with this plugin.

## DREAM3D Mailing Lists ##

If you need more help with a filter, please consider asking your question on the DREAM3D Users mailing list:
https://groups.google.com/forum/?hl=en#!forum/dream3d-users
