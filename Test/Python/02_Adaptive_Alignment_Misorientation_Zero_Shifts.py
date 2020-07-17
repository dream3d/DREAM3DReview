'''
Pipeline example based on 02_Adaptive Alignment - Misorientation - Zero Shifts in Anisotropy examples
'''
import os
import simpl
import simpl_helpers as sh
import simpl_test_dirs as sd
import orientationanalysispy
import dream3dreviewpy

def start_test():
    # Create Data Container Array
    dca = simpl.DataContainerArray()

    # Read H5EBSD File
    err = orientationanalysispy.read_h5_ebsd(dca, 'AlMgSc Data', 'Phase Data', 'EBSD SEM Scan Data',
                                            sd.GetBuildDirectory() + '/Data/Anisotropy/AlMgsh.h5ebsd',
                                             0, 9, True, sh.AngleRepresentation.Radians,
                                             simpl.StringSet({'Fit', 'Image Quality', 'EulerAngles',
                                                             'SEM Signal', 'Confidence Index', 'Phases',
                                                             'X Position', 'Y Position'}))
    assert err == 0, f'ReadH5Ebsd ErrorCondition {err}'

    # Convert Orientation Representation
    err = orientationanalysispy.convert_orientations(dca, 0, 2,
                                                     simpl.DataArrayPath('AlMgSc Data', 'EBSD SEM Scan Data',
                                                                         'EulerAngles'),
                                                     'Quats')
    assert err == 0, f'Convert Orientations ErrorCondition {err}'

    # Adaptive Alignment Misorientation
    err = dream3dreviewpy.adaptive_alignment_misorientation(dca, 5, False,
                                                         simpl.DataArrayPath('AlMgSc Data', 'EBSD SEM Scan Data',
                                                                             'Quats'),
                                                         simpl.DataArrayPath('AlMgSc Data', 'EBSD SEM Scan Data',
                                                                             'Phases'),
                                                         simpl.DataArrayPath('', '', ''),
                                                         simpl.DataArrayPath('AlMgSc Data', 'Phase Data',
                                                                             'CrystalStructures'))
    assert err == 0, f'AdaptiveAlignment Misorientation {err}'

if __name__ == '__main__':
    print('Starting Test %s ' % os.path.basename(__file__))
    start_test()
    print('Ending Test %s ' % os.path.basename(__file__))
