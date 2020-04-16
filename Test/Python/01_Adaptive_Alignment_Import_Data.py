'''
Pipeline example based on 01_Import Data in Anisotropy examples
'''

import simpl
import simplpy
import simpl_helpers as sc
import simpl_test_dirs as sd
import orientationanalysispy

def import_data():
    # Create Data Container Array
    dca = simpl.DataContainerArray()

    err = orientationanalysispy.ebsd_to_h5_ebsd(dca, sd.GetBuildDirectory() + '/Data/Anisotropy/ang',
                          sd.GetBuildDirectory() + '/Data/Output/Anisotropy/AlMgSc.h5ebsd',
                          'AlMgSc-TD_', '', 'ang', 3,
                          0, 9, simpl.AxisAngleInput(0, 0, 0, 0), simpl.AxisAngleInput(0, 0, 0, 0))
    if err < 0:
        print('EBSD to H5EBSD ErrorCondition %d' % err)


if __name__ == '__main__':
    import_data()
