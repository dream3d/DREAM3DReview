'''
Pipeline example based on 01_Import Data in Anisotropy examples
'''

import simpl
import simplpy
import simpl_helpers as sh
import simpl_test_dirs as sd
import orientationanalysispy

def start_test():
    # Create Data Container Array
    dca = simpl.DataContainerArray()

    err = orientationanalysispy.ebsd_to_h5_ebsd(dca, sd.GetBuildDirectory() + '/Data/Anisotropy/ang',
                          sd.GetBuildDirectory() + '/Data/Output/Anisotropy/AlMgSc.h5ebsd',
                          'AlMgSc-TD_', '', 'ang', 3,
                          0, 9, simpl.AxisAngleInput(0, 0, 0, 0), simpl.AxisAngleInput(0, 0, 0, 0))
    assert err == 0, f'EBSD to H5EBSD ErrorCondition {err}'


if __name__ == '__main__':
    print('Starting Test %s ' % os.path.basename(__file__) )
    start_test()
    print('Ending Test %s ' % os.path.basename(__file__) )

