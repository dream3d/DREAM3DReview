# Based on CreateLambertSphereSurface pipeline example
# Tests the Find Attribute Array Statistics filter
# These are the simpl_py python modules
import os 
import simpl
import simplpy
import simpl_helpers as sh
import simpl_test_dirs as sd
import orientationanalysispy
import dream3dreviewpy

def start_test():
    # Create Data Container Array
    dca = simpl.DataContainerArray()

    # Create Data Container
    err = simplpy.create_data_container(dca, 'ImageDataContainer')
    assert err == 0, f'DataContainer ErrorCondition: {err}'

    # Create Image Geometry
    err = simplpy.create_image_geometry(dca, 'ImageDataContainer', simpl.IntVec3([101, 101, 1]),
                                        simpl.FloatVec3([0, 0, 0]),
                                        simpl.FloatVec3([1, 1, 1]))
    assert err == 0, f'ImageGeometry ErrorCondition: {err}'

    # Create Attribute Matrix #1 for Lambert Sphere
    new_row = simpl.VectorDouble()
    new_row.append(101)
    new_row.append(101)
    new_row.append(1)
    table_data = list()
    table_data.append(new_row)
    err = simplpy.create_attribute_matrix(dca, simpl.DataArrayPath('ImageDataContainer', 'CellAttributeMatrix', ''), 3,
                                          simpl.DynamicTableData(table_data, ['0', '1', '2'], ['0']))
    assert err == 0, f'AttributeMatrix #1 ErrorCondition: {err}'

    # Create Data Array
    err = simplpy.create_data_array(dca, simpl.ScalarTypes.UInt8, 1,
                                    simpl.DataArrayPath('ImageDataContainer', 'CellAttributeMatrix', 'ScalarValues'),
                                    simpl.InitializationType.Manual, '128', (0.0, 151.1))
    assert err == 0, f'DataArray ErrorCondition: {err}'

    # Create Lambert Sphere
    err = orientationanalysispy.create_lambert_sphere(dca, sh.Hemisphere.Northern,
                                                      simpl.DataArrayPath('ImageDataContainer', 'CellAttributeMatrix',
                                                                          'ScalarValues'),
                                                      'QuadDataContainer', 'TriangleDataContainer', 'EdgeDataContainer',
                                                      'VertexDataContainer', 'VertexAttributeMatrix',
                                                      'EdgeAttributeMatrix', 'FaceAttributeMatrix',
                                                      True, True, True, True)
    assert err == 0, f'LambertSphere ErrorCondition: {err}'

    # Create Attribute Matrix #2 for Statistics
    new_row = simpl.VectorDouble()
    new_row.append(1)
    table_data = list()
    table_data.append(new_row)
    err = simplpy.create_attribute_matrix(dca,
                                          simpl.DataArrayPath('QuadDataContainer', 'StatisticsAttributeMatrix', ''),
                                          simpl.AttributeMatrix.Type.Generic,
                                          simpl.DynamicTableData(table_data, ['0'], ['0']))
    assert err == 0, f'AttributeMatrix #2 ErrorCondition: {err}'

    # Find Attribute Array Statistics
    err = dream3dreviewpy.find_array_statistics(data_container_array=dca,
                                                find_histogram=True, 
                                                find_length=True, 
                                                find_min=True, 
                                                find_max=True, 
                                                find_mean=True, 
                                                find_median=True, 
                                                find_std_deviation=True,
                                                find_summation=True,
                                                use_mask=False,
                                                standardize_data=False,
                                                compute_by_index=False,
                                                destination_attribute_matrix=simpl.DataArrayPath('QuadDataContainer', 'StatisticsAttributeMatrix', ''),
                                                mask_array_path=simpl.DataArrayPath('', '', ''), 
                                                length_array_name='Length', 
                                                minimum_array_name='MinValue', 
                                                maximum_array_name='MaxValue',
                                                mean_array_name='Mean', 
                                                median_array_name='MedianValue', 
                                                std_deviation_array_name='StandardDeviation', 
                                                summation_array_name='Summation', 
                                                standardized_array_name='',
                                                selected_array_path=simpl.DataArrayPath('QuadDataContainer', 'FaceAttributeMatrix','Quad_ScalarValues'),
                                                feature_ids_array_path=simpl.DataArrayPath('', '', ''),
                                                histogram_array_name='Histogram',
                                                use_full_range=False,
                                                num_bins=100,
                                                min_range=0.0,
                                                max_range=1.0)
    assert err == 0, f'FindArrayStatistics ErrorCondition: {err}'

    # Write to DREAM3D File
    err = simplpy.data_container_writer(dca, sd.GetBuildDirectory() +
                                        '/Data/Output/UnitTest/Plugins/DREAM3DReview/Test/Temp/FindArrayStats.dream3d', True,
                                        False)
    assert err == 0, f'DataContainerWriter ErrorCondition: {err}'

    # Export ASCII Data
    selected_data_array_paths = [simpl.DataArrayPath('QuadDataContainer',
                                                     'StatisticsAttributeMatrix',
                                                     'Length'),
                                 simpl.DataArrayPath('QuadDataContainer',
                                                     'StatisticsAttributeMatrix',
                                                     'Mean'),
                                 simpl.DataArrayPath('QuadDataContainer',
                                                     'StatisticsAttributeMatrix',
                                                     'StandardDeviation'),
                                 simpl.DataArrayPath('QuadDataContainer',
                                                     'StatisticsAttributeMatrix',
                                                     'Summation')]
    outStyle = 0 
    err = simplpy.write_ascii_data(dca, selected_data_array_paths,
                                   sd.GetBuildDirectory() + '/Data/Output/UnitTest/Plugins/DREAM3DReview/Test/Temp/', # Only needed for Multi-File output
                                   sd.GetBuildDirectory() + '/Data/Output/UnitTest/Plugins/DREAM3DReview/Test/Temp/Statistics.csv', # Only needed for Single File Style
                                   simpl.DelimiterTypes.Comma, '.csv', 10, outStyle)
    assert err == 0, f'WriteAsciiData ErrorCondition: {err}'
    print(f'Output files Written to {sd.GetBuildDirectory()}/Data/Output/UnitTest/Plugins/DREAM3DReview/Test/Temp/')

if __name__ == '__main__':
    print('Starting Test %s ' % os.path.basename(__file__))
    start_test()
    print('Ending Test %s ' % os.path.basename(__file__))