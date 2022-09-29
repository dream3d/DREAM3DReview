# Convert TDMS to H5 #

## Group (Subgroup) ##

Example(Sub Example)

## Required Conda Packages ##

+ h5py
+ nptdms

## Description ##

This **Filter** converts a series of TDMS files into individual part files, where each part is stored in an HDF5 file. Three offsets can be used for the camera, diode, and laser.

## Parameters ##

| Name | Type | Description |
|------|------|------|
| Area Offset | int | The offset to translate the area |
| Intensity Offset | int | The offset to change the intensity of the image’s pixels |
| Laser Offset | int | The offset to change the focus of the laser |
| Output Folder | File Path | The output path for the HDF5 files |
| Group | String | Specifies group name in TDMS file |
| Input Directory | File Path | The input folder where TDMS files are located |
| File Ordering | int | Ascending: Low to High; Descending: High to Low |
| File Prefix | String | Prefix of file before the number |
| File Suffix | String | Suffix of file after the number |
| File Extension | String | File extension without ‘.’ |
| Start Index | int | Start of index for series of TDMS files |
| End Index | int | End of index for series of TDMS files |
| Increment Index | int | Value to increment the index |
| Padding Digits | int | Number of padding digits used in the file names |

## Required Geometry ##

Not Applicable

## Required Objects ##

Not Applicable

## Created Objects ##

Not Applicable

## Example Pipelines ##

List the names of the example pipelines where this filter is used.

## License & Copyright ##

Please see the description file distributed with this plugin.

## DREAM3D Mailing Lists ##

If you need more help with a filter, please consider asking your question on the DREAM3D Users mailing list:
https://groups.google.com/forum/?hl=en#!forum/dream3d-users
