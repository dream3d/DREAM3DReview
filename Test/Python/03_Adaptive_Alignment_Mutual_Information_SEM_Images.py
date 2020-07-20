'''
Pipeline example based on 03_Adaptive Alignment - Mutual Information - SEM Images in Anisotropy examples
'''
import os
import simpl
import simplpy
import simpl_helpers as sh
import simpl_test_dirs as sd
import orientationanalysispy
import dream3dreviewpy
import itkimageprocessing
import itkimageprocessingpy

def start_test():
    # Create Data Container Array
    dca = simpl.DataContainerArray()

    # Register the ImageIO Factories
    imageWriter = itkimageprocessing.ITKImageWriter.New()    
    imageWriter.registerImageIOFactories()

    # Read H5EBSD File
    print('Loading H5EBSD File')
    err = orientationanalysispy.read_h5_ebsd(dca, 'AlMgSc Data', 'Phase Data', 'EBSD SEM Scan Data',
                                            sd.GetBuildDirectory() + '/Data/Anisotropy/AlMgsc.h5ebsd',
                                            0, 9, True, sh.AngleRepresentation.Radians,
                                            simpl.StringSet({'Fit', 'Image Quality', 'EulerAngles',
                                                             'SEM Signal', 'Confidence Index', 'Phases'}))
    assert err == 0, f'ReadH5Ebsd ErrorCondition {err}'
    
    # Import Image Stack [ITK]
    print('Loading Images...')
    fileListInfo = simpl.FileListInfo(3, 0, 0, 9, 1, sd.GetBuildDirectory() + '/Data/Anisotropy/tif',
                                      'AlMgSc-TD_', '', 'tif')
    err = itkimageprocessingpy.itk_import_image_stack(dca, 'SEMAlMgSc Data', 'EBSD SEM Scan Data',
                                                      simpl.FloatVec3([0, 0, 0]), simpl.FloatVec3([1, 1, 1]),
                                                      fileListInfo, 10, 'ImageData')
    assert err == 0, f'ITK Import Image Stack ErrorCondition {err}'

    # Convert Orientation Representation
    print('Creating Quats....')
    err = orientationanalysispy.convert_orientations(dca, 0, 2,
                                                    simpl.DataArrayPath('AlMgSc Data', 'EBSD SEM Scan Data', 'EulerAngles'),
                                                    'Quats')
    assert err == 0, f'Convert Orientations ErrorCondition {err}'

    # Adaptive Alignment Mutual Information
    err = dream3dreviewpy.adaptive_alignment_mutual_information(dca, 5, False,
                                                             simpl.DataArrayPath('AlMgSc Data', 'EBSD SEM Scan Data', 'Quats'),
                                                             simpl.DataArrayPath('AlMgSc Data', 'EBSD SEM Scan Data', 'Phases'),
                                                             simpl.DataArrayPath('', '', ''),
                                                             simpl.DataArrayPath('AlMgSc Data', 'Phase Data', 'CrystalStructures'),
                                                             1,  # Global Correction
                                                             simpl.DataArrayPath('SEMAlMgSc Data', 'EBSD SEM Scan Data', 'ImageData'),
                                                             )
    assert err == 0, f'AdaptiveAlignment Mutual Information {err}'

if __name__ == '__main__':
    print('Starting Test %s ' % os.path.basename(__file__))
    start_test()
    print('Ending Test %s ' % os.path.basename(__file__))
