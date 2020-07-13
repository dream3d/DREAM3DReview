# Based on CreateVertexGeometry example
# Tests the Remove Flagged Vertices filter
# These are the simpl_py python modules
import os
import simpl
import simplpy
import simpl_helpers as sc
import simpl_test_dirs as sd
import dream3dreviewpy

def start_test():
    # Create Data Container Array
    dca = simpl.DataContainerArray()

    # Create the Data Container
    err = simplpy.create_data_container(dca, 'DataContainer')
    assert err == 0, f'DataContainer ErrorCondition: {err}'

    # Import ASCII Data - #1 - Vertex Coordinates
    import_file = sd.GetBuildDirectory() + '/Data/SIMPL/VertexCoordinates.csv'
    wizard_data = {
        'inputFilePath': import_file,
        'beginIndex': 2,
        'numberOfLines': 145,
        'delimiters': [','],
        'consecutiveDelimiters': False,
        'automaticAM': True,
        'selectedPath': simpl.DataArrayPath('DataContainer', 'Bounds', ''),
        'headers': ['x', 'y', 'z'],
        'attrMatType': 3,
        'tupleDimensions': [144],
        'dataTypes': ['float', 'float', 'float']
    }
    err = simplpy.read_ascii_data(dca, wizard_data)
    assert err == 0, f'Import ASCII Data #1 -  ErrorCondition: {err}'

    # Combine Attribute Arrays # 1:
    selected_data_array_paths = [simpl.DataArrayPath('DataContainer', 'Bounds', 'x'),
                                 simpl.DataArrayPath('DataContainer', 'Bounds', 'y'),
                                 simpl.DataArrayPath('DataContainer', 'Bounds', 'z')]
    err = simplpy.combine_attribute_arrays(dca, selected_data_array_paths, 'Vertices', False)
    assert err == 0, f'Combined Attribute Arrays #1 -  ErrorCondition: {err}'

    # Delete Data # 1
    dcap = simpl.DataContainerArrayProxy()
    dcap.getDataContainerProxy('DataContainer').Flag = 0
    dcap.getDataContainerProxy('DataContainer').getAttributeMatrixProxy('Bounds').Flag = 0
    dcap.getDataContainerProxy('DataContainer').getAttributeMatrixProxy('Bounds').getDataArrayProxy('x').Flag = 2
    dcap.getDataContainerProxy('DataContainer').getAttributeMatrixProxy('Bounds').getDataArrayProxy('y').Flag = 2
    dcap.getDataContainerProxy('DataContainer').getAttributeMatrixProxy('Bounds').getDataArrayProxy('z').Flag = 2
    err = simplpy.remove_arrays(dca, dcap)
    assert err == 0, f'Remove Arrays #1 -  ErrorCondition: {err}'

    # Create Geometry
    err = sc.CreateGeometry(dca, 0, simpl.IGeometry.Type.Vertex, 'DataContainer', False,
                            shared_vertex_list_array_path=simpl.DataArrayPath('DataContainer', 'Bounds', 'Vertices'),
                            vertex_attribute_matrix_name='VertexData')
    assert err == 0, f'Create Geometry -  ErrorCondition: {err}'

    # Extract Attribute Arrays from Geometry
    empty_data_array_path = simpl.DataArrayPath('', '', '')
    vertex_list_path = simpl.DataArrayPath('DataContainer', 'VertexData', 'VertexCoordinates')
    err = simplpy.extract_attribute_arrays_from_geometry(dca, 'DataContainer', vertex_list_path,
                                                         empty_data_array_path, empty_data_array_path,
                                                         empty_data_array_path, empty_data_array_path,
                                                         empty_data_array_path, empty_data_array_path,
                                                         empty_data_array_path, empty_data_array_path,
                                                         empty_data_array_path, empty_data_array_path,
                                                         empty_data_array_path)
    assert err == 0, f'ExtractAttributeArraysFromGeometry ErrorCondition {err}'

    # Split Multicomponent Attribute Array
    err = simplpy.split_attribute_array(dca,
                                        simpl.DataArrayPath('DataContainer', 'VertexData',
                                                            'VertexCoordinates'),
                                        '_Component_')
    assert err == 0, f'SplitAttributeArray ErrorCondition {err}'

    # Threshold Objects
    err = sc.MultiThresholdObjects(dca, 'Mask', [('DataContainer', 'VertexData',
                                                  'VertexCoordinates_Component_0', '>', -0.5),
                                                 ('DataContainer', 'VertexData',
                                                  'VertexCoordinates_Component_1', '<', 0.5)])
    assert err == 0, f'MultiThresholdObjects ErrorCondition {err}'

    # Remove Flagged Vertices
    err = dream3dreviewpy.remove_flagged_vertices(dca, 'DataContainer',
                                                  simpl.DataArrayPath('DataContainer', 'VertexData',
                                                                      'Mask'),
                                                  'ReducedVertexDataContainer')
    assert err == 0, f'RemoveFlaggedVertices ErrorCondition {err}'

    # Delete Data
    err = sc.RemoveArrays(dca, [['DataContainer', '', '']])
    assert err, f'RemoveArrays ErrorCondition {err}'

    # Write to DREAM3D File
    err = simplpy.data_container_writer(dca, sd.GetBuildDirectory() +
                                        '/Data/Output/DREAM3DReview/' +
                                        'RemoveFlaggedVertices.dream3d', True, False)
    assert err == 0, f'DataContainerWriter ErrorCondition: {err}'

if __name__ == '__main__':
    print('Starting Test %s ' % os.path.basename(__file__))
    start_test()
    print('Ending Test %s ' % os.path.basename(__file__))