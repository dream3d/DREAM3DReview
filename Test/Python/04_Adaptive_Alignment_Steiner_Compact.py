'''
Pipeline example based on 04_Steiner Compact in Anisotropy examples
'''
import os
import simpl
import simplpy
import simpl_helpers as sh
import simpl_test_dirs as sd
import orientationanalysispy
import dream3dreviewpy
import reconstructionpy
import processingpy

def start_test():
    # Create Data Container Array
    dca = simpl.DataContainerArray()

    # Read H5EBSD File
    err = orientationanalysispy.read_h5_ebsd(dca, 'AlMgSc Data', 'Phase Data', 'EBSD SEM Scan Data',
                                            sd.GetBuildDirectory() + '/Data/Anisotropy/AlMgsc.h5ebsd',
                                            0, 9, True, sh.AngleRepresentation.Radians,
                                            simpl.StringSet({'Fit', 'Image Quality', 'EulerAngles',
                                                             'SEM Signal', 'Confidence Index', 'Phases',
                                                             'X Position', 'Y Position'}))
    assert err == 0, f'ReadH5Ebsd ErrorCondition {err}'

    # MultiThreshold Objects
    err = sh.MultiThresholdObjects(dca, 'Mask', [('AlMgSc Data', 'EBSD SEM Scan Data', 'Confidence Index', '>', 0.05),
                                                 ('AlMgSc Data', 'EBSD SEM Scan Data', 'Image Quality', '>', 800)])
    assert err == 0, f'MultiThresholdObjects ErrorCondition: {err}'

    # Convert Orientation Representation
    err = orientationanalysispy.convert_orientations(dca, 0, 2,
                                                    simpl.DataArrayPath('AlMgSc Data', 'EBSD SEM Scan Data', 'EulerAngles'),
                                                    'Quats')
    assert err == 0, f'Convert Orientations ErrorCondition {err}'

    # Neighbor Orientation Comparison (Bad Data)
    err = orientationanalysispy.bad_data_neighbor_orientation_check(dca, 2, 4,
                                              simpl.DataArrayPath('AlMgSc Data', 'EBSD SEM Scan Data', 'Mask'),
                                              simpl.DataArrayPath('AlMgSc Data', 'EBSD SEM Scan Data', 'Phases'),
                                              simpl.DataArrayPath('AlMgSc Data', 'Phase Data', 'CrystalStructures'),
                                              simpl.DataArrayPath('AlMgSc Data', 'EBSD SEM Scan Data', 'Quats'))
    assert err == 0, f'BadDataNeighborOrientationCheck ErrorCondition {err}'

    # Neighbor Orientation Correlation
    err = orientationanalysispy.neighbor_orientation_correlation(dca, 2, 0.05, 2,
                                           simpl.DataArrayPath('AlMgSc Data', 'EBSD SEM Scan Data', 'Confidence Index'),
                                           simpl.DataArrayPath('AlMgSc Data', 'EBSD SEM Scan Data', 'Phases'),
                                           simpl.DataArrayPath('AlMgSc Data', 'Phase Data', 'CrystalStructures'),
                                           simpl.DataArrayPath('AlMgSc Data', 'EBSD SEM Scan Data', 'Quats'))
    assert err == 0, f'NeighborOrientationCorrelation ErrorCondition {err}'

    # EBSD Segment Features (Misorientation)
    err = reconstructionpy.ebsd_segment_features(dca, 'CellFeatureData', 2,
                                                 True,
                                                 simpl.DataArrayPath('AlMgSc Data', 'EBSD SEM Scan Data', 'Mask'),
                                                 simpl.DataArrayPath('AlMgSc Data', 'EBSD SEM Scan Data', 'Phases'),
                                                 simpl.DataArrayPath('AlMgSc Data', 'Phase Data', 'CrystalStructures'),
                                                 simpl.DataArrayPath('AlMgSc Data', 'EBSD SEM Scan Data', 'Quats'),
                                                 'FeatureIds',
                                                 'Active')
    assert err == 0, f'EbsdSegmentFeatures ErrorCondition {err}'

    # Erode Bad Data
    err = processingpy.erode_dilate_bad_data(dca, sh.BadDataOperation.Erode, 3, True, True, True,
                                           simpl.DataArrayPath('AlMgSc Data', 'EBSD SEM Scan Data', 'FeatureIds'))
    assert err == 0, f'ErodeDilateBadData ErrorCondition {err}'

    # Steiner Compact

    # The following 'enums' need binded
    # Section Plane options --> [0 : XY] , [1 : XZ], [2 : YZ]
    # Number of Sites options --> [0 : 8] , [1 : 12], [2 : 16], [3 : 24], [4 : 36]
    err = dream3dreviewpy.steiner_compact(dca, True, sd.GetBuildDirectory() + '/Debug/SteinerCompactXY.vtk',
                          True, sd.GetBuildDirectory() + '/Debug/SteinerCompactXY.txt',
                          simpl.DataArrayPath('AlMgSc Data', 'EBSD SEM Scan Data', 'FeatureIds'),
                          simpl.DataArrayPath('AlMgSc Data', 'EBSD SEM Scan Data', 'Phases'),
                          0, 1)
    assert err == 0, f'SteinerCompact ErrorCondition {err}'

if __name__ == '__main__':
    print('Starting Test %s ' % os.path.basename(__file__))
    start_test()
    print('Ending Test %s ' % os.path.basename(__file__))
