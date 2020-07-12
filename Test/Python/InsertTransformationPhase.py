# Insert Transformation Phase
# Based on InsertTransformationPhase example pipeline
import os
import simpl
import simpl_helpers as sc
import simpl_test_dirs as sd
import simplpy
import orientationanalysispy
import genericpy
import statisticspy
import syntheticbuildingpy as syntheticbuilding
import dream3dreviewpy

def start_test():
    # Create Data Container Array
    dca = simpl.DataContainerArray()

     # Using the GenerateStatsData and CreateDynamicTableData functions
    euler_dynamic_table_data = sc.CreateDynamicTableData([[0, 0, 0, 0, 0]],
                                                      ['Euler 1', 'Euler 2', 'Euler 3', 'Weight', 'Sigma'], ['0'])
    axis_dynamic_table_data = sc.CreateDynamicTableData([[0, 0, 0, 0, 0]],
                                                     ['Angle(w)', 'Axis (h)', 'Axis (k)', 'Axis (l)', 'Weight (MRD)'],
                                                     ['0'])
    err = syntheticbuilding.generate_primary_stats_data(dca, 'Primary', 0, 1, 0,
                               1, 2.3, 0.4, 4, 3, 10.5,
                               True, 'StatsGeneratorDataContainer',
                               'CellEnsembleData', False,
                               simpl.DataArrayPath('', '', ''),
                               euler_dynamic_table_data,
                               axis_dynamic_table_data,
                               euler_dynamic_table_data)

    assert err == 0, f'Primary StatsData ErrorCondition: {err}'

     # Initialize Synthetic Volume
    err = syntheticbuilding.initialize_synthetic_volume(dca, 'SyntheticVolumeDataContainer', 'CellData',
                                                        'CellEnsembleMatrixName',
                                                        6,
                                                        simpl.IntVec3([128, 128, 128]),
                                                        simpl.FloatVec3([0.5, 0.5, 0.5]),
                                                        simpl.FloatVec3([0, 0, 0]),
                                                        simpl.DataArrayPath('StatsGeneratorDataContainer',
                                                                            'CellEnsembleData', 'Statistics'),
                                                        simpl.DataArrayPath('StatsGeneratorDataContainer',
                                                                            'CellEnsembleData', 'PhaseTypes'),
                                                        simpl.DataArrayPath('StatsGeneratorDataContainer',
                                                                            'CellEnsembleData', 'PhaseName'),
                                                        False, 0, 'NOT NEEDED')
    assert err == 0, f'InitializeSyntheticVolume ErrorCondition: {err}'

    # Establish Shape Types
    err = syntheticbuilding.establish_shape_types(dca,
                                                    simpl.DataArrayPath('StatsGeneratorDataContainer',
                                                                        'CellEnsembleData', 'PhaseTypes'),
                                                    'ShapeTypes', [simpl.ShapeType.Ellipsoid])
    assert err == 0, f'EstablishShapeTypes ErrorCondition: {err}'

    # Pack Primary Phases
    err = syntheticbuilding.pack_primary_phases(dca,
                                                  simpl.DataArrayPath('SyntheticVolumeDataContainer', 'CellData',
                                                                      ''),
                                                  'CellFeatureData', 'CellEnsembleData', 'FeatureIds', 'Phases',
                                                  'Phases', 'NumFeatures',
                                                  simpl.DataArrayPath('StatsGeneratorDataContainer',
                                                                      'CellEnsembleData', 'Statistics'),
                                                  simpl.DataArrayPath('StatsGeneratorDataContainer',
                                                                      'CellEnsembleData', 'PhaseTypes'),
                                                  simpl.DataArrayPath('StatsGeneratorDataContainer',
                                                                      'CellEnsembleData', 'PhaseName'),
                                                  simpl.DataArrayPath('StatsGeneratorDataContainer',
                                                                      'CellEnsembleData', 'ShapeTypes'),
                                                  simpl.DataArrayPath('', '', ''), False, 0, '', '', False,
                                                  False, False,
                                                  simpl.DataArrayPath('', '', ''),
                                                  simpl.DataArrayPath('', '', ''))
    assert err == 0, f'PackPrimaryPhases ErrorCondition: {err}'

    # Find Feature Neighbors
    err = statisticspy.find_neighbors(dca, simpl.DataArrayPath('SyntheticVolumeDataContainer', 'CellFeatureData', ''),
                                      'SharedSurfaceAreaList', 'NeighborList',
                                      simpl.DataArrayPath('SyntheticVolumeDataContainer', 'CellData',
                                                          'FeatureIds'),
                                      '', 'NumNeighbors', 'SurfaceFeatures', False, True)
    assert err == 0, f'FindNeighbors ErrorCondition: {err}'

    # Match Crystallography
    err = syntheticbuilding.match_crystallography(dca, simpl.DataArrayPath('StatsGeneratorDataContainer',
                                                                             'CellEnsembleData', 'Statistics'),
                                                    simpl.DataArrayPath('StatsGeneratorDataContainer',
                                                                        'CellEnsembleData', 'CrystalStructures'),
                                                    simpl.DataArrayPath('StatsGeneratorDataContainer',
                                                                        'CellEnsembleData', 'PhaseTypes'),
                                                    simpl.DataArrayPath('SyntheticVolumeDataContainer',
                                                                        'CellData', 'FeatureIds'),
                                                    simpl.DataArrayPath('SyntheticVolumeDataContainer',
                                                                        'CellFeatureData', 'Phases'),
                                                    simpl.DataArrayPath('SyntheticVolumeDataContainer',
                                                                        'CellFeatureData', 'SurfaceFeatures'),
                                                    simpl.DataArrayPath('SyntheticVolumeDataContainer',
                                                                        'CellFeatureData', 'NeighborList'),
                                                    simpl.DataArrayPath('SyntheticVolumeDataContainer',
                                                                        'CellFeatureData', 'SharedSurfaceAreaList'),
                                                    simpl.DataArrayPath('SyntheticVolumeDataContainer',
                                                                        'CellEnsembleData', 'NumFeatures'),
                                                    'EulerAngles', 'Volumes', 'EulerAngles', 'AvgQuats', 100000)
    assert err == 0, f'MatchCrystallography ErrorCondition: {err}'

    # Find Feature Centroids
    err = genericpy.find_feature_centroids(dca,
                                           simpl.DataArrayPath('SyntheticVolumeDataContainer',
                                                               'CellData', 'FeatureIds'),
                                           simpl.DataArrayPath('SyntheticVolumeDataContainer',
                                                               'CellFeatureData', 'Centroids'))
    assert err == 0, f'FindFeatureCentroids ErrorCondition: {err}'

    # Find Feature Sizes #1
    err = statisticspy.find_sizes(dca, simpl.DataArrayPath('SyntheticVolumeDataContainer', 'CellFeatureData', ''),
                                  simpl.DataArrayPath('SyntheticVolumeDataContainer',
                                                      'CellData', 'FeatureIds'),
                                  'Volumes2', 'EquivalentDiameters', 'NumElements', False)
    assert err == 0, f'FindSizes #1 ErrorCondition {err}'

    # Insert Transformation Phases
    # Possible candidate for a helper function
    # Possible enums: trans_crystal_struct
    # May want to use a 'dictionary' for the numerical portion in helper function
    err = dream3dreviewpy.insert_transformation_phases(dca,
                                                             1, 1, 60, simpl.FloatVec3([1, 1, 1]), True, True, 1, 0.2, 1,
                                                             0,
                                                             simpl.DataArrayPath('StatsGeneratorDataContainer',
                                                                                 'CellEnsembleData', ''),
                                                             simpl.DataArrayPath('SyntheticVolumeDataContainer',
                                                                                 'CellFeatureData', ''),
                                                             simpl.DataArrayPath('SyntheticVolumeDataContainer',
                                                                                 'CellData', 'FeatureIds'),
                                                             simpl.DataArrayPath('SyntheticVolumeDataContainer',
                                                                                 'CellData', 'EulerAngles'),
                                                             simpl.DataArrayPath('SyntheticVolumeDataContainer',
                                                                                 'CellData', 'Phases'),
                                                             simpl.DataArrayPath('SyntheticVolumeDataContainer',
                                                                                 'CellFeatureData', 'AvgQuats'),
                                                             simpl.DataArrayPath('SyntheticVolumeDataContainer',
                                                                                 'CellFeatureData', 'Centroids'),
                                                             simpl.DataArrayPath('SyntheticVolumeDataContainer',
                                                                                 'CellFeatureData',
                                                                                 'EquivalentDiameters'),
                                                             simpl.DataArrayPath('SyntheticVolumeDataContainer',
                                                                                 'CellFeatureData', 'EulerAngles'),
                                                             simpl.DataArrayPath('SyntheticVolumeDataContainer',
                                                                                 'CellFeatureData', 'Phases'),
                                                             'ParentIds',
                                                             simpl.DataArrayPath('SyntheticVolumeDataContainer',
                                                                                 'CellEnsembleData',
                                                                                 'NumFeaturesPerParent'),
                                                             simpl.DataArrayPath('StatsGeneratorDataContainer',
                                                                                 'CellEnsembleData',
                                                                                 'CrystalStructures'),
                                                             simpl.DataArrayPath('StatsGeneratorDataContainer',
                                                                                 'CellEnsembleData', 'PhaseTypes'),
                                                             simpl.DataArrayPath('StatsGeneratorDataContainer',
                                                                                 'CellEnsembleData', 'ShapeTypes'),
                                                             simpl.DataArrayPath('SyntheticVolumeDataContainer',
                                                                                 'CellEnsembleData', 'NumFeatures'))
    assert err == 0, f'InsertTransformationPhases ErrorCondition {err}'

    # Delete Data (using helper function)
    sc.RemoveArrays(dca, [('SyntheticVolumeDataContainer', 'CellFeatureData', 'AvgQuats'),
                          ('SyntheticVolumeDataContainer', 'CellFeatureData', 'Centroids'),
                          ('SyntheticVolumeDataContainer', 'CellFeatureData', 'EquivalentDiameters'),
                          ('SyntheticVolumeDataContainer', 'CellFeatureData', 'EulerAngles'),
                          ('SyntheticVolumeDataContainer', 'CellFeatureData', 'NeighborList'),
                          ('SyntheticVolumeDataContainer', 'CellFeatureData', 'NumElements'),
                          ('SyntheticVolumeDataContainer', 'CellFeatureData', 'NumFeaturesPerParent'),
                          ('SyntheticVolumeDataContainer', 'CellFeatureData', 'NumNeighbors'),
                          ('SyntheticVolumeDataContainer', 'CellFeatureData', 'SharedSurfaceAreaList'),
                          ('SyntheticVolumeDataContainer', 'CellFeatureData', 'SurfaceFeatures'),
                          ('SyntheticVolumeDataContainer', 'CellFeatureData', 'Volumes'),
                          ('SyntheticVolumeDataContainer', 'CellFeatureData', 'Volumes2')])

    # Generate IPF Colors
    err = orientationanalysispy.generate_ipf_colors(dca, simpl.FloatVec3([0, 0, 1]),
                                                    simpl.DataArrayPath('SyntheticVolumeDataContainer',
                                                                        'CellData', 'Phases'),
                                                    simpl.DataArrayPath('SyntheticVolumeDataContainer',
                                                                        'CellData', 'EulerAngles'),
                                                    simpl.DataArrayPath('StatsGeneratorDataContainer',
                                                                        'CellEnsembleData', 'CrystalStructures'),
                                                    False,
                                                    simpl.DataArrayPath('', '', ''), 'IPFColor')
    assert err == 0, f'GenerateIPFColors ErrorCondition: {err}'

    # Find Feature Sizes #2
    err = statisticspy.find_sizes(dca, simpl.DataArrayPath('SyntheticVolumeDataContainer', 'CellFeatureData', ''),
                                  simpl.DataArrayPath('SyntheticVolumeDataContainer',
                                                      'CellData', 'FeatureIds'),
                                  'Volumes2', 'EquivalentDiameters', 'NumElements', False)
    assert err == 0, f'FindSizes #2 ErrorCondition {err}'

    # Convert Orientation Representation
    err = orientationanalysispy.convert_orientations(dca, 0, 2,
                                                     simpl.DataArrayPath('SyntheticVolumeDataContainer', 'CellData',
                                                                         'EulerAngles'),
                                                     'Quats')
    assert err == 0, f'Convert Orientations ErrorCondition {err}'

    # Find Feature Average Orientations
    err = orientationanalysispy.find_avg_orientations(dca,
                                                      simpl.DataArrayPath('SyntheticVolumeDataContainer', 'CellData',
                                                                          'FeatureIds'),
                                                      simpl.DataArrayPath('SyntheticVolumeDataContainer', 'CellData',
                                                                          'Phases'),
                                                      simpl.DataArrayPath('SyntheticVolumeDataContainer', 'CellData',
                                                                          'Quats'),
                                                      simpl.DataArrayPath('SyntheticVolumeDataContainer',
                                                                          'CellEnsembleData', 'CrystalStructures'),
                                                      simpl.DataArrayPath('SyntheticVolumeDataContainer',
                                                                          'CellFeatureData', 'AvgQuats'),
                                                      simpl.DataArrayPath('SyntheticVolumeDataContainer',
                                                                          'CellFeatureData', 'AvgEulers'))
    assert err == 0, f'FindAvgOrientations ErrorCondition: {err}'

    # Find Feature Centroids
    err = genericpy.find_feature_centroids(dca, simpl.DataArrayPath('SyntheticVolumeDataContainer', 'CellData',
                                                                    'FeatureIds'),
                                           simpl.DataArrayPath('SyntheticVolumeDataContainer', 'CellFeatureData',
                                                               'Centroids'))
    assert err == 0, f'FindFeatureCentroids ErrorCondition: {err}'

    # Find Feature Neighborhoods
    err = statisticspy.find_neighborhoods(dca, 'NeighborhoodList', 1,
                                          simpl.DataArrayPath('SyntheticVolumeDataContainer', 'CellFeatureData',
                                                              'EquivalentDiameters'),
                                          simpl.DataArrayPath('SyntheticVolumeDataContainer', 'CellFeatureData',
                                                              'Phases'),
                                          simpl.DataArrayPath('SyntheticVolumeDataContainer', 'CellFeatureData',
                                                              'Centroids'),
                                          'Neighborhoods')
    assert err == 0, f'FindNeighborhoods ErrorCondition: {err}'

    # Find Feature Neighbors #2
    err = statisticspy.find_neighbors(dca, simpl.DataArrayPath('SyntheticVolumeDataContainer', 'CellFeatureData', ''),
                                      'SharedSurfaceAreaList', 'NeighborList',
                                      simpl.DataArrayPath('SyntheticVolumeDataContainer', 'CellData', 'FeatureIds'),
                                      '', 'NumNeighbors', 'SurfaceFeatures', False, False)
    assert err == 0, f'FindNeighbors #2 ErrorCondition: {err}'

    # Write to DREAM3D file
    err = sc.WriteDREAM3DFile(
        sd.GetBuildDirectory() + '/Data/Output/TransformationPhase/InsertTransformationPhases.dream3d',
        dca)
    assert err == 0, f'WriteDREAM3DFile ErrorCondition: {err}'

    # Export Feature Data as CSV File
    err = simplpy.feature_data_csv_writer(dca,
                                          simpl.DataArrayPath('SyntheticVolumeDataContainer', 'CellFeatureData', ''),
                                          sd.GetBuildDirectory() +
                                          '/Data/Output/TransformationPhase/InsertTransformationPhases.csv',
                                          False, simpl.DelimiterTypes.Comma, True)
    assert err == 0, f'FeatureDataCsvWriter ErrorCondition: {err}'

if __name__ == '__main__':
    print('Starting Test %s ' % os.path.basename(__file__))
    start_test()
    print('Ending Test %s ' % os.path.basename(__file__))