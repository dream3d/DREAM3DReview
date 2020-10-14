# CreateArrayofIndices #


## Group (Subgroup) ##

DREAM3DReview (DREAM3DReview)

## Description ##

This **Filter** creates an array whose value is the index at every location. This works on any geometry in any type of attribute matrix. 

## Parameters ##

| Name | Type | Description |
|------|------|------|
| None | None | None |

## Required Geometry ##

All 

## Required Objects ##

| Kind | Default Name | Type | Component Dimensions | Description |
|------|--------------|-------------|---------|-----|
| **Attribute Matrix** | Attribute Matrix Name | Any | N/A | Attribute Matrix to Put Index Array In |


## Created Objects ##

| Kind | Default Name | Type | Component Dimensions | Description |
|------|--------------|-------------|---------|-----|
| **Data Array** | Index Array | uint64 | 1 | Array of indices in attribute matrix |


## Example Pipelines ##

List the names of the example pipelines where this filter is used.

## License & Copyright ##

Please see the description file distributed with this plugin.

## DREAM3D Mailing Lists ##

If you need more help with a filter, please consider asking your question on the DREAM3D Users mailing list:
https://groups.google.com/forum/?hl=en#!forum/dream3d-users