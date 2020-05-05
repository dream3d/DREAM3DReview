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

#---------------
# This is the list of Filters Superclasses. Other filters use these filters as a
# superclass. If the filters are NOT meant to be ever invoked from the user interface
# then they go here. This is also so that the python wrapping will work correctly.
set(_SuperclassFilters
  AdaptiveAlignment
)

#-----------------
# Loop on the Filter Superclasses adding each one to the DREAM3DLib project so that it gets compiled.
foreach(f ${_SuperclassFilters} )
  ADD_SIMPL_SUPERCLASS_FILTER("${PLUGIN_NAME}" "${PLUGIN_NAME}" ${_filterGroupName} ${f})
endforeach()

#---------
# List your public filters here
set(_PublicFilters
  ApplyTransformationToGeometry
  AverageEdgeFaceCellArrayToVertexArray
  AverageVertexArrayToEdgeFaceCellArray
  ComputeUmeyamaTransform
  EstablishFoamMorphology
  DBSCAN
  DownsampleVertexGeometry
  ExtractInternalSurfacesFromTriangleGeometry
  FindArrayStatistics
  FindElementCentroids
  FindNorm
  FindVertexToTriangleDistances
  IterativeClosestPoint
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
  WaveFrontObjectFileWriter

  # 2019-03-06 Filter AFRL Filter Release
  AlignGeometries
  ApproximatePointCloudHull
  CombineStlFiles
  ExportCLIFile
  ExtractTripleLinesFromTriangleGeometry
 # FindDynamicArrayStatistics
  FindLayerStatistics
  FindMinkowskiBouligandDimension
  FindSurfaceRoughness
  ImportCLIFile
  ImportVolumeGraphicsFile
 # InterpolateMeshToRegularGrid
  InterpolatePointCloudToRegularGrid
  LabelTriangleGeometry
  LaplacianSmoothPointCloud
  MapPointCloudToRegularGrid
  ReadBinaryCTNorthStar
  SliceTriangleGeometry

  #HEDM
  ReadMicData
  TesselateFarFieldGrains

  #Anisotropy
  AdaptiveAlignmentFeature
  AdaptiveAlignmentMisorientation
  AdaptiveAlignmentMutualInformation
  SteinerCompact

  #MASSIF
  FFTHDFWriterFilter
  ImportMASSIFData

  #DDDAnalysis
  LocalDislocationDensityCalculator
  IdentifyDislocationSegments
  DiscretizeDDDomain
  ParaDisReader
  
  # TransformationPhase
  FindCSLBoundaries
  InsertTransformationPhases
  TiDwellFatigueCrystallographicAnalysis

  # Imported from AFRLDistributionC
  ImportQMMeltpoolTDMSFile
  ImportPrintRiteHDF5File
  ImportPrintRiteTDMSFiles
  DelaunayTriangulation
  GenerateFeatureIDsbyBoundingBoxes
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
                        TRUE)
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
                        ${${PLUGIN_NAME}_SOURCE_DIR}/Documentation/${_filterGroupName}/${f}.md FALSE)
endforeach()

#-------------
# These are files that need to be compiled but are NOT filters
ADD_SIMPL_SUPPORT_HEADER_SUBDIR(${${PLUGIN_NAME}_SOURCE_DIR} ${_filterGroupName} KMeansTemplate.hpp util/ClusteringAlgorithms)
ADD_SIMPL_SUPPORT_HEADER_SUBDIR(${${PLUGIN_NAME}_SOURCE_DIR} ${_filterGroupName} KMedoidsTemplate.hpp util/ClusteringAlgorithms)
ADD_SIMPL_SUPPORT_HEADER_SUBDIR(${${PLUGIN_NAME}_SOURCE_DIR} ${_filterGroupName} DBSCANTemplate.hpp util/ClusteringAlgorithms)
ADD_SIMPL_SUPPORT_HEADER_SUBDIR(${${PLUGIN_NAME}_SOURCE_DIR} ${_filterGroupName} SilhouetteTemplate.hpp util/EvaluationAlgorithms)
ADD_SIMPL_SUPPORT_HEADER_SUBDIR(${${PLUGIN_NAME}_SOURCE_DIR} ${_filterGroupName} KDistanceTemplate.hpp util/EvaluationAlgorithms)
ADD_SIMPL_SUPPORT_HEADER_SUBDIR(${${PLUGIN_NAME}_SOURCE_DIR} ${_filterGroupName} DistanceTemplate.hpp util)
ADD_SIMPL_SUPPORT_HEADER_SUBDIR(${${PLUGIN_NAME}_SOURCE_DIR} ${_filterGroupName} nanoflann.hpp util) 


ADD_SIMPL_SUPPORT_HEADER(${${PLUGIN_NAME}_SOURCE_DIR} ${_filterGroupName} HEDM/H5MicImporter.h)
ADD_SIMPL_SUPPORT_SOURCE(${${PLUGIN_NAME}_SOURCE_DIR} ${_filterGroupName} HEDM/H5MicImporter.cpp)
ADD_SIMPL_SUPPORT_HEADER(${${PLUGIN_NAME}_SOURCE_DIR} ${_filterGroupName} HEDM/H5MicReader.h)
ADD_SIMPL_SUPPORT_SOURCE(${${PLUGIN_NAME}_SOURCE_DIR} ${_filterGroupName} HEDM/H5MicReader.cpp)
ADD_SIMPL_SUPPORT_HEADER(${${PLUGIN_NAME}_SOURCE_DIR} ${_filterGroupName} HEDM/H5MicVolumeReader.h)
ADD_SIMPL_SUPPORT_SOURCE(${${PLUGIN_NAME}_SOURCE_DIR} ${_filterGroupName} HEDM/H5MicVolumeReader.cpp)
ADD_SIMPL_SUPPORT_HEADER(${${PLUGIN_NAME}_SOURCE_DIR} ${_filterGroupName} HEDM/MicFields.h)
ADD_SIMPL_SUPPORT_SOURCE(${${PLUGIN_NAME}_SOURCE_DIR} ${_filterGroupName} HEDM/MicFields.cpp)
ADD_SIMPL_SUPPORT_HEADER(${${PLUGIN_NAME}_SOURCE_DIR} ${_filterGroupName} HEDM/MicPhase.h)
ADD_SIMPL_SUPPORT_SOURCE(${${PLUGIN_NAME}_SOURCE_DIR} ${_filterGroupName} HEDM/MicPhase.cpp)
ADD_SIMPL_SUPPORT_HEADER(${${PLUGIN_NAME}_SOURCE_DIR} ${_filterGroupName} HEDM/MicReader.h)
ADD_SIMPL_SUPPORT_SOURCE(${${PLUGIN_NAME}_SOURCE_DIR} ${_filterGroupName} HEDM/MicReader.cpp)
ADD_SIMPL_SUPPORT_HEADER(${${PLUGIN_NAME}_SOURCE_DIR} ${_filterGroupName} HEDM/MicHeaderEntry.h)
ADD_SIMPL_SUPPORT_HEADER(${${PLUGIN_NAME}_SOURCE_DIR} ${_filterGroupName} HEDM/MicConstants.h)

ADD_SIMPL_SUPPORT_HEADER(${${PLUGIN_NAME}_SOURCE_DIR} ${_filterGroupName} util/PrintRiteHelpers.h)
ADD_SIMPL_SUPPORT_HEADER(${${PLUGIN_NAME}_SOURCE_DIR} ${_filterGroupName} util/Delaunay2D.h)
ADD_SIMPL_SUPPORT_SOURCE(${${PLUGIN_NAME}_SOURCE_DIR} ${_filterGroupName} util/Delaunay2D.cpp)
ADD_SIMPL_SUPPORT_HEADER(${${PLUGIN_NAME}_SOURCE_DIR} ${_filterGroupName} util/TriMesh.h)
ADD_SIMPL_SUPPORT_SOURCE(${${PLUGIN_NAME}_SOURCE_DIR} ${_filterGroupName} util/TriMesh.cpp)
ADD_SIMPL_SUPPORT_HEADER(${${PLUGIN_NAME}_SOURCE_DIR} ${_filterGroupName} util/TriMeshPrimitives.hpp)

#---------------------
# This macro must come last after we are done adding all the filters and support files.
SIMPL_END_FILTER_GROUP(${${PLUGIN_NAME}_BINARY_DIR} "${_filterGroupName}" "${PLUGIN_NAME}")

