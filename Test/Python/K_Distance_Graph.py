# Based on CreateQuadGeometry example
# Tests the K Distance Graph filter

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

    # Import ASCII Data - #2 - Edge Connectivity
    import_file = sd.GetBuildDirectory() + '/Data/SIMPL/QuadConnectivity.csv'
    wizard_data = {
        'inputFilePath': import_file,
        'beginIndex': 2,
        'numberOfLines': 122,
        'delimiters': [','],
        'consecutiveDelimiters': False,
        'automaticAM': True,
        'selectedPath': simpl.DataArrayPath('DataContainer', 'QuadList', ''),
        'headers': ['V0', 'V1', 'V2', 'V3'],
        'attrMatType': 3,
        'tupleDimensions': [121],
        'dataTypes': ['uint64_t', 'uint64_t', 'uint64_t', 'uint64_t']
    }
    err = simplpy.read_ascii_data(dca, wizard_data)
    assert err == 0, f'Import ASCII Data #2 -  ErrorCondition: {err}'

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

    # Combine Attribute Arrays #2:
    selected_data_array_paths = [simpl.DataArrayPath('DataContainer', 'QuadList', 'V0'),
                                 simpl.DataArrayPath('DataContainer', 'QuadList', 'V1'),
                                 simpl.DataArrayPath('DataContainer', 'QuadList', 'V2'),
                                 simpl.DataArrayPath('DataContainer', 'QuadList', 'V3')]
    err = simplpy.combine_attribute_arrays(dca, selected_data_array_paths, 'Quads', False)
    assert err == 0, f'Combined Attribute Arrays #2 -  ErrorCondition: {err}'

    # Delete Data # 2
    dcap = simpl.DataContainerArrayProxy()
    dcap.getDataContainerProxy('DataContainer').Flag = 0
    dcap.getDataContainerProxy('DataContainer').getAttributeMatrixProxy('QuadList').Flag = 0
    dcap.getDataContainerProxy('DataContainer').getAttributeMatrixProxy('QuadList').getDataArrayProxy('V0').Flag = 2
    dcap.getDataContainerProxy('DataContainer').getAttributeMatrixProxy('QuadList').getDataArrayProxy('V1').Flag = 2
    dcap.getDataContainerProxy('DataContainer').getAttributeMatrixProxy('QuadList').getDataArrayProxy('V2').Flag = 2
    dcap.getDataContainerProxy('DataContainer').getAttributeMatrixProxy('QuadList').getDataArrayProxy('V3').Flag = 2
    err = simplpy.remove_arrays(dca, dcap)
    assert err == 0, f'Remove Arrays #2 -  ErrorCondition: {err}'

    # Create Geometry
    err = sc.CreateGeometry(dca, 0, simpl.IGeometry.Type.Quad, 'DataContainer', False,
                            shared_vertex_list_array_path=simpl.DataArrayPath('DataContainer', 'Bounds', 'Vertices'),
                            shared_quad_list_array_path=simpl.DataArrayPath('DataContainer', 'QuadList', 'Quads'),
                            vertex_attribute_matrix_name='VertexData',
                            face_attribute_matrix_name='FaceData')
    assert err == 0, f'Create Geometry -  ErrorCondition: {err}'

    # K Distance Graph
    err = dream3dreviewpy.k_distance_graph(dca, simpl.DataArrayPath('DataContainer', 'QuadList',
                                                                    'Quads'),
                                           False, simpl.DataArrayPath('', '', ''),
                                           simpl.DataArrayPath('DataContainer', 'QuadList',
                                                               'KDistance'), 7, 3)
    assert err == 0, f'KDistanceGraph ErrorCondition: {err}'

    # Write DREAM3D File
    err = simplpy.data_container_writer(dca, sd.GetBuildDirectory() +
                                        '/Data/Output/DREAM3DReview/' +
                                        'KDistanceGraph.dream3d', True, False)
    assert err == 0, f'DataContainerWriter ErrorCondition: {err}'

if __name__ == '__main__':
    print('Starting Test %s ' % os.path.basename(__file__))
    start_test()
    print('Ending Test %s ' % os.path.basename(__file__))