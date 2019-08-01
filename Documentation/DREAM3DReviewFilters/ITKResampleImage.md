ITK Resample Image
=============

## Group (Subgroup) ##

DREAM3D Review (Registration)

## Description ##

This **Filter** resamples a single data array, series of images, or series of CTF files according to a transformation file provided by the user. The transformation file is an .hdf5 file that contains information about the transformation. More information can be found [here](ITKPairwiseRegistration.md). The transformation file should contain the information about the type of transformation

If the user selects 'Single Data Array' as the Opeartion Mode, the user will be required to set the moving image as an **Attribute Array** that has somehow been already read in using another filer. 

If the user selects 'Series of Images' the user will have to point to a folder on their file system where the moving images are located and where the resampled images should be output to. The images must be somehow numerically ordered, but the user can select a subset of images by selecting the start and end index of the image numbering. The total number of images must equal the total number of transforms (number of HDF5 groups) in the transform HDF5 file specified by the user. 

If the user selects 'Series of CTFs' as the Operation mode, it will work similarly as above except instead of resampling images, the CTF files will be resampled and written out to the specified folder.

## Parameters ##

| Name | Type | Description |
|------|------|-------------|
| Operation Mode  | int | Whether to resample a Single Data Array, a Series of Images, or a Series of CTF files|
| Transform File Name | Path | Location of HDF5 Transform file that describes the transformation to apply | 
| Interpolation Type | int | Choose either linear or or nearest-neighbor for the type of interpolation when resampling the image | 
| Output Image File Name Prefix | string | Choose the prefix for the output file name if 'Series of Images' was selected | 
| Output EBSD File Name Prefix | string | Chose the prefix for the output EBSD file name if 'Series of CTFs' was selected | 
| Output Path for Resampled Images | Path | Choose the path for the output images if 'Series of Images' was selected | 
| Output Path for Resampled Orientation Data | Path | Choose the path for the output images if 'Series of CTFs' was selected | 


## Required Objects 
| Kind | Default Name | Type | Component Dimensions | Description |
|------|--------------|------|----------------------|-------------|
| Moving Image | |**Attribute Array** | Any | The **Attribute Array** to be resampled if 'Single Data Array' was selected | 


## Created Objects 
| Kind | Default Name | Type | Component Dimensions | Description |
|------|--------------|------|----------------------|-------------|
| Data Container | ResampledDataDC | **Data Container** | N/A | The newly created data container for the resampled image if 'Single Data Array' was selected | 
| Cell Attribute Matrix | ResampledDataAM | **Attribute Matrix** | N/A | The newly created cell attribute matrix for the resampled image if 'Single Data Array' was selected| 
| Image Data | ResampledData | **Attribute Array** | Same as Moving Image  | The resampled image data if 'Single Data Array' was selected| 

 

## Example Pipelines ##



## License & Copyright ##

Please see the description file distributed with this plugin.

## DREAM3D Mailing Lists ##

If you need more help with a filter, please consider asking your question on the DREAM3D Users mailing list:
