# 

# These are the simpl_py python modules
import os
import simpl
import simplpy
import simpl_helpers as sh
import simpl_test_dirs as sd
import dream3dreviewpy

MELTPOOL_DATACONTAINER_NAME = "MeltPool Data"
VERTEX_ATTRIBUTE_MATRIX_NAME = "Vertex Data"
LASERTTL_NAME = "LaserTTL"
MASK_ARRAY_NAME = "Mask"
VOXEL_INDICES_NAME = "VoxelIndices"
AREA_ARRAY_NAME = "Area"
INTENSITY_ARRAY_NAME = "Intensity"
UNIFORM_INTERPOLATION_TECHNIQUE = 0
IMAGE_DATA_CONTAINER_NAME = "ImageDataContainer"
INTERPOLATED_ATTRIBUTE_MATRIX_NAME = "InterpolatedAttributeMatrix"


def start_test():
# Create Data Container Array
    dca = simpl.DataContainerArray()

    # Import the QM Meltpool HDF5 file
    print(f'Importing QMMeltpool HDF5 File....')
    input_files = [sd.GetBuildDirectory() + "/Data/DREAM3DReview/QMMeltPool_Data.h5"]
    err = dream3dreviewpy.import_qm_meltpool_h5_file(data_container_array=dca, 
                                                      input_files = input_files, 
                                                      data_container_path = MELTPOOL_DATACONTAINER_NAME, 
                                                      vertex_attribute_matrix_name = VERTEX_ATTRIBUTE_MATRIX_NAME, 
                                                      slice_range = simpl.IntVec2([200,204]), 
                                                      power = 10)
    assert err == 0, f'import_qm_meltpool_h5_file: {err}'

    # Threshold Objects
    # Create the selected thresholds / comparison inputs for MultiThresholdObjects filter
    print(f'Thresholding Data....')
    selectedThresholds = simpl.ComparisonInputs()
    selectedThresholds.addInput(MELTPOOL_DATACONTAINER_NAME, VERTEX_ATTRIBUTE_MATRIX_NAME, LASERTTL_NAME,
                                simpl.ComparisonOperators.GreaterThan, 0.0)
    err = simplpy.multi_threshold_objects(dca, MASK_ARRAY_NAME, selectedThresholds)
    assert err == 0, f'MultiThresholdObjects ErrorCondition: {err}'

    # Map the vertices into a regular grid
    print(f'Mapping vertices onto a regular grid....')
    err = dream3dreviewpy.map_point_cloud_to_regular_grid(dca, 
                                    data_container_name = simpl.DataArrayPath(MELTPOOL_DATACONTAINER_NAME, "", ""), 
                                    created_image_data_container_name = IMAGE_DATA_CONTAINER_NAME, 
                                    # image_data_container_path = None, 
                                    voxel_indices_array_path = simpl.DataArrayPath(MELTPOOL_DATACONTAINER_NAME, VERTEX_ATTRIBUTE_MATRIX_NAME, VOXEL_INDICES_NAME),
                                    grid_dimensions = simpl.IntVec3([2048, 2048, 5]), 
                                    use_mask = True,
                                    sampling_grid_type = 0, 
                                    mask_array_path = simpl.DataArrayPath(MELTPOOL_DATACONTAINER_NAME, VERTEX_ATTRIBUTE_MATRIX_NAME, MASK_ARRAY_NAME) )
    assert err == 0, f'map_point_cloud_to_regular_grid ErrorCondition: {err}'

    # Sample the vertices onto the regular grid that was created in the last step
    print(f'Interpolating Vertex Data....')
    err = dream3dreviewpy.interpolate_point_cloud_to_regular_grid(dca, 
                                  data_container_name = MELTPOOL_DATACONTAINER_NAME,
                                  arrays_to_interpolate = [simpl.DataArrayPath(MELTPOOL_DATACONTAINER_NAME, VERTEX_ATTRIBUTE_MATRIX_NAME, AREA_ARRAY_NAME),
                                                           simpl.DataArrayPath(MELTPOOL_DATACONTAINER_NAME, VERTEX_ATTRIBUTE_MATRIX_NAME, INTENSITY_ARRAY_NAME)],
                                  arrays_to_copy = [simpl.DataArrayPath(MELTPOOL_DATACONTAINER_NAME, VERTEX_ATTRIBUTE_MATRIX_NAME, AREA_ARRAY_NAME),
                                                           simpl.DataArrayPath(MELTPOOL_DATACONTAINER_NAME, VERTEX_ATTRIBUTE_MATRIX_NAME, INTENSITY_ARRAY_NAME)],
                                  voxel_indices_array_path = simpl.DataArrayPath(MELTPOOL_DATACONTAINER_NAME, VERTEX_ATTRIBUTE_MATRIX_NAME, VOXEL_INDICES_NAME),
                                  interpolated_data_container_name = IMAGE_DATA_CONTAINER_NAME,
                                  interpolated_attribute_matrix_name = INTERPOLATED_ATTRIBUTE_MATRIX_NAME,
                                  interpolation_technique = UNIFORM_INTERPOLATION_TECHNIQUE,
                                  kernel_size = simpl.FloatVec3([75.0, 75.0, 1.0]),
                                  sigmas = None,
                                  use_mask = True,
                                  mask_array_path = simpl.DataArrayPath(MELTPOOL_DATACONTAINER_NAME, VERTEX_ATTRIBUTE_MATRIX_NAME, MASK_ARRAY_NAME),
                                  store_kernel_distances = False,
                                  kernel_distances_array_name = None,
                                  interpolated_suffix = " [Interpolated]",
                                  copy_suffix = " [Copied]")
    assert err == 0, f'interpolate_point_cloud_to_regular_grid ErrorCondition: {err}'

    # Find some scalar statistics on each voxel cell so that we can create an image to export
    print(f'Finding Statistics for each voxel....')
    dream3dreviewpy.find_neighbor_list_statistics(dca, 
                                      find_length = False,
                                      find_min = False,
                                      find_max = True, 
                                      find_mean = True, 
                                      find_median = True, 
                                      find_std_deviation = False, 
                                      find_summation = False, 
                                      destination_attribute_matrix = simpl.DataArrayPath(IMAGE_DATA_CONTAINER_NAME, INTERPOLATED_ATTRIBUTE_MATRIX_NAME, ""), 
                                      length_array_name = None, 
                                      minimum_array_name = None, 
                                      maximum_array_name = AREA_ARRAY_NAME + " [Maximum]", 
                                      mean_array_name = AREA_ARRAY_NAME + " [Mean]", 
                                      median_array_name = AREA_ARRAY_NAME + " [Median]", 
                                      std_deviation_array_name = None, 
                                      summation_array_name = None,
                                      selected_array_path = simpl.DataArrayPath(IMAGE_DATA_CONTAINER_NAME, INTERPOLATED_ATTRIBUTE_MATRIX_NAME, AREA_ARRAY_NAME + " [Interpolated]"))

    # Write to DREAM3D File
    print(f'Writing DREAM3D file ....')
    err = simplpy.data_container_writer(dca, sd.GetTestTempDirectory() +
                                        '/DREAM3DReview/' +
                                        'QMMeltPool_Data.dream3d',
                                        True, False)
    assert err == 0, f'DataContainerWriter ErrorCondition: {err}'




if __name__ == '__main__':
    print('Starting Test %s ' % os.path.basename(__file__))
    start_test()
    print('Ending Test %s ' % os.path.basename(__file__))