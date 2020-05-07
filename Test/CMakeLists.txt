set(PLUGIN_NAME "DREAM3DReview")

set(${PLUGIN_NAME}Test_SOURCE_DIR ${${PLUGIN_NAME}_SOURCE_DIR}/Test)
set(${PLUGIN_NAME}Test_BINARY_DIR ${${PLUGIN_NAME}_BINARY_DIR}/Test)

#------------------------------------------------------------------------------
# List all the source files here. They will NOT be compiled but instead
# be directly included in the main test source file. We list them here so that
# they will show up in IDEs
set(TEST_NAMES
   ImportQMMeltpoolTDMSFileTest
  ImportVolumeGraphicsFileTest
  GenerateMaskFromSimpleShapesTest
  GenerateFeatureIDsbyBoundingBoxesTest
  AnisotropyFilterTest
  ImportMASSIFDataTest
  FFTHDFWriterFilterTest
  EstablishFoamMorphologyTest
)

#------------------------------------------------------------------------------
# Include this file from the CMP Project
include(${CMP_SOURCE_DIR}/cmpCMakeMacros.cmake)
include(${SIMPLProj_SOURCE_DIR}/Source/SIMPLib/SIMPLibMacros.cmake)


get_filename_component(${PLUGIN_NAME}_PARENT_SOURCE_DIR "${${PLUGIN_NAME}_SOURCE_DIR}" DIRECTORY)
get_filename_component(${PLUGIN_NAME}_PARENT_BINARY_DIR "${${PLUGIN_NAME}_BINARY_DIR}" DIRECTORY)

set(Anisotropy_Data_DIR ${${PLUGIN_NAME}_SOURCE_DIR}/Data/Anisotropy)
list(APPEND ${PLUGIN_NAME}_LINK_LIBS Qt5::Core SIMPLib ${ITK_LIBRARIES})

SIMPL_GenerateUnitTestFile(PLUGIN_NAME ${PLUGIN_NAME}
                           TEST_DATA_DIR ${${PLUGIN_NAME}_SOURCE_DIR}/Test/Data
                           SOURCES ${TEST_NAMES}
                           LINK_LIBRARIES ${${PLUGIN_NAME}_LINK_LIBS} ${plug_target_name}
                           INCLUDE_DIRS ${${PLUGIN_NAME}_PARENT_SOURCE_DIR}
                                        ${${PLUGIN_NAME}Test_SOURCE_DIR}
                                        ${${PLUGIN_NAME}Test_BINARY_DIR}
                                        ${${PLUGIN_NAME}_PARENT_BINARY_DIR}
)

set(TEST_SCRIPT_FILE_EXT "sh")
set(EXE_EXT "")
if(WIN32)
  set(TEST_SCRIPT_FILE_EXT "bat")
  set(EXE_EXT ".exe")
endif()
# !!! NO Spaces in the following names. NO subdirectories either
# The names MUST match the names of the files on DISK including cApiTaliZation
set(PREBUILT_PIPELINE_NAMES
  "Open-Cell-Foam-Example"
  "ExtractInternalSurfacesFromTriangleGeometry"
  "WaveFront_Export"
)

set(PREBUILT_SUPPORT_DIR "${${PLUGIN_NAME}_SOURCE_DIR}/ExamplePipelines/${PLUGIN_NAME}")

foreach(test ${PREBUILT_PIPELINE_NAMES} )
  set(Example_PIPELINE_FILE "${${PLUGIN_NAME}_SOURCE_DIR}/ExamplePipelines/${PLUGIN_NAME}/${test}.json")
  set(SIMPL_CTEST_TEST_DRIVER "${${PLUGIN_NAME}Test_BINARY_DIR}/PrebuiltPipelineTests/${test}.${TEST_SCRIPT_FILE_EXT}")
  if(EXISTS "${DREAM3DProj_SOURCE_DIR}/Test/ctest_pipeline_driver.${TEST_SCRIPT_FILE_EXT}")
    configure_file(${DREAM3DProj_SOURCE_DIR}/Test/ctest_pipeline_driver.${TEST_SCRIPT_FILE_EXT} "${SIMPL_CTEST_TEST_DRIVER}" @ONLY)
    string(REPLACE "/" "_" test "${test}")
    add_test(NAME ${PLUGIN_NAME}_Example_${test} COMMAND "${SIMPL_CTEST_TEST_DRIVER}" )
  endif()
endforeach()

#------------------------------------------------------------------------------
# If Python is enabled, then enable the Python unit tests for this plugin
if(SIMPL_ENABLE_PYTHON)
  set(PLUGIN_PYTHON_TESTS
    "Apply_Transformation_To_Geometry"
    "Average_EdgeFaceCell_Array_To_Vertex_Array"
    "Average_Vertex_Array_To_EdgeFaceCell_Array"
    "Compute_Umeyama_Transform"
    "DBSCAN"
    "Extract_Internal_Surfaces_From_Triangle_Geometry"
    "Find_Attribute_Array_Statistics"
    "Find_Element_Centroids"
    "Find_Norm"
    "K_Medoids"
    "K_Distance_Graph"
    "K_Means_Test"
    "Normalize_Attribute_Arrays"
    "Point_Sample_Triangle_Geometry"
    "Potts_Model"
    "Principal_Component_Analysis"
    "Remove_Flagged_Vertices"
    "Robust_Automatic_Threshold"
    # "01_Adaptive_Alignment_Import_Data" # The conversion is already complete so this does not need to run
    "02_Adaptive_Alignment_Misorientation_Zero_Shifts"
    "03_Adaptive_Alignment_Mutual_Information_SEM_Images"
    "04_Adaptive_Alignment_Steiner_Compact"
    "InsertTransformationPhase"
  )

  include(${SIMPLProj_SOURCE_DIR}/Wrapping/Python/Binding/CreatePybind11Module.cmake)

  CreatePythonTests(PREFIX "PY_${PLUGIN_NAME}"
    INPUT_DIR ${${PLUGIN_NAME}Test_SOURCE_DIR}/Python
    TEST_NAMES ${PLUGIN_PYTHON_TESTS}
  )
endif()
