# SurfaceMeshToWaveFront #

## Group (Subgroup) ##

@FilterGroup@ (@FilterSubgroup@)

## Description ##

This **Filter** will export a **Triangle Geometry** into a WaveFront .obj file. The only items that are exported are the triangles. Data attached to the triangles cannot be exported at this time.

## Parameters ##

| Name | Type | Description |
|------|------|------|
| Output File| String  | Path to the output file. |

## Required Geometry ##

Triangle Geometry

## Required Objects ##

| Kind | Default Name | Type | Component Dimensions | Description |
|------|--------------|-------------|---------|-----|
| **Data Container** | Data Container | N/A | N/A | The DataContainer that holds the triangle geometry |


## Created Objects ##

No objects are creatd


## Example Pipelines ##

Prebuilt Pipelines->Examples->DREAM3DReview->WaveFront Export

## License & Copyright ##

Please see the description file distributed with this plugin.

## DREAM3D Mailing Lists ##

If you need more help with a filter, please consider asking your question on the DREAM3D Users mailing list:
https://groups.google.com/forum/?hl=en#!forum/dream3d-users