set(_filterGroupName DREAM3DReviewFilters)
set(${_filterGroupName}_FILTERS_HDRS "")

#--------
# This macro must come first before we start adding any filters
SIMPL_START_FILTER_GROUP(
  ALL_FILTERS_HEADERFILE ${AllFiltersHeaderFile}
  REGISTER_KNOWN_FILTERS_FILE ${RegisterKnownFiltersFile}
  FILTER_GROUP "${_filterGroupName}"
  BINARY_DIR ${${PLUGIN_NAME}_BINARY_DIR}
  )

#---------
# List your public filters here
set(_PublicFilters
  ApplyTransformationToGeometry
  AverageEdgeFaceCellArrayToVertexArray
  AverageVertexArrayToEdgeFaceCellArray
  ComputeUmeyamaTransform

  DBSCAN
  ExtractInternalSurfacesFromTriangleGeometry
  FindArrayStatistics
  FindElementCentroids
  FindNorm
  KDistanceGraph
  KMeans
  KMedoids
  NormalizeArrays
  PointSampleTriangleGeometry
  PottsModel
  PrincipalComponentAnalysis
  RemoveFlaggedVertices
  RobustAutomaticThreshold
  Silhouette

  # Registration 
  ITKPairwiseImageRegistration
  ITKResampleImage
  FindMaskNeighbors
)

list(LENGTH _PublicFilters PluginNumFilters)
set_property(GLOBAL PROPERTY PluginNumFilters ${PluginNumFilters})

#--------------
# Loop on all the filters adding each one. In this loop we default to making each filter exposed in the user
# interface in DREAM3D. If you want to have the filter compiled but NOT exposed to the user then use the next loop
foreach(f ${_PublicFilters} )
  ADD_SIMPL_FILTER(  "${PLUGIN_NAME}" "${PLUGIN_NAME}"
                        ${_filterGroupName} ${f}
                        ${${PLUGIN_NAME}_SOURCE_DIR}/Documentation/${_filterGroupName}/${f}.md 
                        TRUE
                        ${${PLUGIN_NAME}_BINARY_DIR})
endforeach()


#---------------
# This is the list of Private Filters. These filters are available from other filters but the user will not
# be able to use them from the DREAM3D user interface.
set(_PrivateFilters

)

#-----------------
# Loop on the Private Filters adding each one to the DREAM3DLib project so that it gets compiled.
foreach(f ${_PrivateFilters} )
  ADD_SIMPL_FILTER(  "${PLUGIN_NAME}" "${PLUGIN_NAME}"
                        ${_filterGroupName} ${f}
                        ${${PLUGIN_NAME}_SOURCE_DIR}/Documentation/${_filterGroupName}/${f}.md FALSE ${${PLUGIN_NAME}_BINARY_DIR})
endforeach()

#-------------
# These are files that need to be compiled but are NOT filters
ADD_SIMPL_SUPPORT_HEADER_SUBDIR(${${PLUGIN_NAME}_SOURCE_DIR} ${_filterGroupName} KMeansTemplate.hpp util/ClusteringAlgorithms)
ADD_SIMPL_SUPPORT_HEADER_SUBDIR(${${PLUGIN_NAME}_SOURCE_DIR} ${_filterGroupName} KMedoidsTemplate.hpp util/ClusteringAlgorithms)
ADD_SIMPL_SUPPORT_HEADER_SUBDIR(${${PLUGIN_NAME}_SOURCE_DIR} ${_filterGroupName} DBSCANTemplate.hpp util/ClusteringAlgorithms)
ADD_SIMPL_SUPPORT_HEADER_SUBDIR(${${PLUGIN_NAME}_SOURCE_DIR} ${_filterGroupName} SilhouetteTemplate.hpp util/EvaluationAlgorithms)
ADD_SIMPL_SUPPORT_HEADER_SUBDIR(${${PLUGIN_NAME}_SOURCE_DIR} ${_filterGroupName} KDistanceTemplate.hpp util/EvaluationAlgorithms)
ADD_SIMPL_SUPPORT_HEADER_SUBDIR(${${PLUGIN_NAME}_SOURCE_DIR} ${_filterGroupName} DistanceTemplate.hpp util)
ADD_SIMPL_SUPPORT_HEADER_SUBDIR(${${PLUGIN_NAME}_SOURCE_DIR} ${_filterGroupName} ITKTransformHelpers.h util)


ADD_SIMPL_SUPPORT_HEADER(${${PLUGIN_NAME}_SOURCE_DIR} ${_filterGroupName} util/EBSDWriterFactory.h)
ADD_SIMPL_SUPPORT_SOURCE(${${PLUGIN_NAME}_SOURCE_DIR} ${_filterGroupName} util/EBSDWriterFactory.cpp)

#---------------------
# This macro must come last after we are done adding all the filters and support files.
SIMPL_END_FILTER_GROUP(${${PLUGIN_NAME}_BINARY_DIR} "${_filterGroupName}" "${PLUGIN_NAME}")

