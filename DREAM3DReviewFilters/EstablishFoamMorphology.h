/*
 * This material is based upon work supported by the National Science Foundation
 * under Grant No. 1629660.  Any opinions, findings, and conclusions or
 * recommendations expressed in this material are those of the author(s) and do not
 * necessarily reflect the views of the National Science Foundation.
 *
 *
 * The code is a derivative work based on the filter PackPrimaryPhases,
 * InsertPrecipitatePhases and EuclideanDistanceMap which were developed under DoD
 * contract:   United States Air Force Prime Contract FA8650-07-D-5800. See those files for their respective copyright
 * and distribution license.
 */

#pragma once

#include "SIMPLib/Common/SIMPLibSetGetMacros.h"
#include "SIMPLib/DataArrays/NeighborList.hpp"
#include "SIMPLib/DataArrays/StatsDataArray.h"
#include "SIMPLib/DataArrays/StringDataArray.h"
#include "SIMPLib/DataContainers/DataContainer.h"
#include "SIMPLib/Filtering/AbstractFilter.h"
#include "SIMPLib/Geometry/ShapeOps/ShapeOps.h"
#include "SIMPLib/SIMPLib.h"

#include "OrientationLib/LaueOps/OrthoRhombicOps.h"

/**
 * @brief The EstablishFoamMorphology class. See [Filter documentation](@ref EstablishFoamMorphology) for details.
 */
class EstablishFoamMorphology : public AbstractFilter
{
  Q_OBJECT

  PYB11_CREATE_BINDINGS(EstablishFoamMorphology SUPERCLASS AbstractFilter)
  PYB11_PROPERTY(DataArrayPath DataContainerName READ getDataContainerName WRITE setDataContainerName)
  PYB11_PROPERTY(DataArrayPath OutputCellAttributeMatrixPath READ getOutputCellAttributeMatrixPath WRITE setOutputCellAttributeMatrixPath)
  PYB11_PROPERTY(DataArrayPath InputCellFeatureIdsArrayPath READ getInputCellFeatureIdsArrayPath WRITE setInputCellFeatureIdsArrayPath)
  PYB11_PROPERTY(QString OutputCellEnsembleAttributeMatrixName READ getOutputCellEnsembleAttributeMatrixName WRITE setOutputCellEnsembleAttributeMatrixName)
  PYB11_PROPERTY(QString OutputCellFeatureAttributeMatrixName READ getOutputCellFeatureAttributeMatrixName WRITE setOutputCellFeatureAttributeMatrixName)
  PYB11_PROPERTY(QString FeatureIdsArrayName READ getFeatureIdsArrayName WRITE setFeatureIdsArrayName)
  PYB11_PROPERTY(QString MaskArrayName READ getMaskArrayName WRITE setMaskArrayName)
  PYB11_PROPERTY(QString CellPhasesArrayName READ getCellPhasesArrayName WRITE setCellPhasesArrayName)
  PYB11_PROPERTY(QString NumFeaturesArrayName READ getNumFeaturesArrayName WRITE setNumFeaturesArrayName)
  PYB11_PROPERTY(QString QPEuclideanDistancesArrayName READ getQPEuclideanDistancesArrayName WRITE setQPEuclideanDistancesArrayName)
  PYB11_PROPERTY(QString TJEuclideanDistancesArrayName READ getTJEuclideanDistancesArrayName WRITE setTJEuclideanDistancesArrayName)
  PYB11_PROPERTY(QString GBEuclideanDistancesArrayName READ getGBEuclideanDistancesArrayName WRITE setGBEuclideanDistancesArrayName)
  PYB11_PROPERTY(DataArrayPath InputStatsArrayPath READ getInputStatsArrayPath WRITE setInputStatsArrayPath)
  PYB11_PROPERTY(DataArrayPath InputPhaseTypesArrayPath READ getInputPhaseTypesArrayPath WRITE setInputPhaseTypesArrayPath)
  PYB11_PROPERTY(DataArrayPath InputPhaseNamesArrayPath READ getInputPhaseNamesArrayPath WRITE setInputPhaseNamesArrayPath)
  PYB11_PROPERTY(DataArrayPath InputShapeTypesArrayPath READ getInputShapeTypesArrayPath WRITE setInputShapeTypesArrayPath)
  PYB11_PROPERTY(QString CsvOutputFile READ getCsvOutputFile WRITE setCsvOutputFile)
  PYB11_PROPERTY(bool PeriodicBoundaries READ getPeriodicBoundaries WRITE setPeriodicBoundaries)
  PYB11_PROPERTY(bool WriteGoalAttributes READ getWriteGoalAttributes WRITE setWriteGoalAttributes)
  PYB11_PROPERTY(int HaveFeatures READ getHaveFeatures WRITE setHaveFeatures)
  PYB11_PROPERTY(double MinStrutThickness READ getMinStrutThickness WRITE setMinStrutThickness)
  PYB11_PROPERTY(float StrutThicknessVariability READ getStrutThicknessVariability WRITE setStrutThicknessVariability)
  PYB11_PROPERTY(float StrutShapeVariability READ getStrutShapeVariability WRITE setStrutShapeVariability)
  PYB11_PROPERTY(bool SmoothStruts READ getSmoothStruts WRITE setSmoothStruts)

public:
  SIMPL_SHARED_POINTERS(EstablishFoamMorphology)
  SIMPL_STATIC_NEW_MACRO(EstablishFoamMorphology)
  SIMPL_TYPE_MACRO_SUPER_OVERRIDE(EstablishFoamMorphology, AbstractFilter)

  using Feature_t = struct
  {
    float m_Volumes;
    float m_EquivalentDiameters;
    float m_AxisLengths[3];
    float m_AxisEulerAngles[3];
    float m_Omega3s;
    int32_t m_FeaturePhases;
    int32_t m_Neighborhoods;
  };

  ~EstablishFoamMorphology() override;

  SIMPL_FILTER_PARAMETER(DataArrayPath, DataContainerName)
  Q_PROPERTY(DataArrayPath DataContainerName READ getDataContainerName WRITE setDataContainerName)

  SIMPL_FILTER_PARAMETER(DataArrayPath, OutputCellAttributeMatrixPath)
  Q_PROPERTY(DataArrayPath OutputCellAttributeMatrixPath READ getOutputCellAttributeMatrixPath WRITE setOutputCellAttributeMatrixPath)

  SIMPL_FILTER_PARAMETER(DataArrayPath, InputCellFeatureIdsArrayPath)
  Q_PROPERTY(DataArrayPath InputCellFeatureIdsArrayPath READ getInputCellFeatureIdsArrayPath WRITE setInputCellFeatureIdsArrayPath)

  SIMPL_INSTANCE_STRING_PROPERTY(OutputCellEnsembleAttributeMatrixName)
  Q_PROPERTY(QString OutputCellEnsembleAttributeMatrixName READ getOutputCellEnsembleAttributeMatrixName WRITE setOutputCellEnsembleAttributeMatrixName)

  SIMPL_INSTANCE_STRING_PROPERTY(OutputCellFeatureAttributeMatrixName)
  Q_PROPERTY(QString OutputCellFeatureAttributeMatrixName READ getOutputCellFeatureAttributeMatrixName WRITE setOutputCellFeatureAttributeMatrixName)

  SIMPL_FILTER_PARAMETER(QString, FeatureIdsArrayName)
  Q_PROPERTY(QString FeatureIdsArrayName READ getFeatureIdsArrayName WRITE setFeatureIdsArrayName)

  SIMPL_FILTER_PARAMETER(QString, MaskArrayName)
  Q_PROPERTY(QString MaskArrayName READ getMaskArrayName WRITE setMaskArrayName)

  SIMPL_FILTER_PARAMETER(QString, CellPhasesArrayName)
  Q_PROPERTY(QString CellPhasesArrayName READ getCellPhasesArrayName WRITE setCellPhasesArrayName)

  SIMPL_FILTER_PARAMETER(QString, NumFeaturesArrayName)
  Q_PROPERTY(QString NumFeaturesArrayName READ getNumFeaturesArrayName WRITE setNumFeaturesArrayName)

  SIMPL_FILTER_PARAMETER(QString, QPEuclideanDistancesArrayName)
  Q_PROPERTY(QString QPEuclideanDistancesArrayName READ getQPEuclideanDistancesArrayName WRITE setQPEuclideanDistancesArrayName)

  SIMPL_FILTER_PARAMETER(QString, TJEuclideanDistancesArrayName)
  Q_PROPERTY(QString TJEuclideanDistancesArrayName READ getTJEuclideanDistancesArrayName WRITE setTJEuclideanDistancesArrayName)

  SIMPL_FILTER_PARAMETER(QString, GBEuclideanDistancesArrayName)
  Q_PROPERTY(QString GBEuclideanDistancesArrayName READ getGBEuclideanDistancesArrayName WRITE setGBEuclideanDistancesArrayName)

  SIMPL_FILTER_PARAMETER(DataArrayPath, InputStatsArrayPath)
  Q_PROPERTY(DataArrayPath InputStatsArrayPath READ getInputStatsArrayPath WRITE setInputStatsArrayPath)

  SIMPL_FILTER_PARAMETER(DataArrayPath, InputPhaseTypesArrayPath)
  Q_PROPERTY(DataArrayPath InputPhaseTypesArrayPath READ getInputPhaseTypesArrayPath WRITE setInputPhaseTypesArrayPath)

  SIMPL_FILTER_PARAMETER(DataArrayPath, InputPhaseNamesArrayPath)
  Q_PROPERTY(DataArrayPath InputPhaseNamesArrayPath READ getInputPhaseNamesArrayPath WRITE setInputPhaseNamesArrayPath)

  SIMPL_FILTER_PARAMETER(DataArrayPath, InputShapeTypesArrayPath)
  Q_PROPERTY(DataArrayPath InputShapeTypesArrayPath READ getInputShapeTypesArrayPath WRITE setInputShapeTypesArrayPath)

  SIMPL_FILTER_PARAMETER(QString, CsvOutputFile)
  Q_PROPERTY(QString CsvOutputFile READ getCsvOutputFile WRITE setCsvOutputFile)

  SIMPL_FILTER_PARAMETER(bool, PeriodicBoundaries)
  Q_PROPERTY(bool PeriodicBoundaries READ getPeriodicBoundaries WRITE setPeriodicBoundaries)

  SIMPL_FILTER_PARAMETER(bool, WriteGoalAttributes)
  Q_PROPERTY(bool WriteGoalAttributes READ getWriteGoalAttributes WRITE setWriteGoalAttributes)

  SIMPL_FILTER_PARAMETER(int, HaveFeatures)
  Q_PROPERTY(int HaveFeatures READ getHaveFeatures WRITE setHaveFeatures)

  SIMPL_FILTER_PARAMETER(double, MinStrutThickness)
  Q_PROPERTY(double MinStrutThickness READ getMinStrutThickness WRITE setMinStrutThickness)

  SIMPL_FILTER_PARAMETER(float, StrutThicknessVariability)
  Q_PROPERTY(float StrutThicknessVariability READ getStrutThicknessVariability WRITE setStrutThicknessVariability)

  SIMPL_FILTER_PARAMETER(float, StrutShapeVariability)
  Q_PROPERTY(float StrutShapeVariability READ getStrutShapeVariability WRITE setStrutShapeVariability)

  SIMPL_FILTER_PARAMETER(bool, SmoothStruts)
  Q_PROPERTY(bool SmoothStruts READ getSmoothStruts WRITE setSmoothStruts)

  /**
   * @brief getCompiledLibraryName Reimplemented from @see AbstractFilter class
   */
  const QString getCompiledLibraryName() const override;

  /**
   * @brief getBrandingString Returns the branding string for the filter, which is a tag
   * used to denote the filter's association with specific plugins
   * @return Branding string
   */
  const QString getBrandingString() const override;

  /**
   * @brief getFilterVersion Returns a version string for this filter. Default
   * value is an empty string.
   * @return
   */
  const QString getFilterVersion() const override;

  /**
   * @brief newFilterInstance Reimplemented from @see AbstractFilter class
   */
  AbstractFilter::Pointer newFilterInstance(bool copyFilterParameters) const override;

  /**
   * @brief getGroupName Reimplemented from @see AbstractFilter class
   */
  const QString getGroupName() const override;

  /**
   * @brief getSubGroupName Reimplemented from @see AbstractFilter class
   */
  const QString getSubGroupName() const override;

  /**
   * @brief getUuid Return the unique identifier for this filter.
   * @return A QUuid object.
   */
  const QUuid getUuid() override;

  /**
   * @brief getHumanLabel Reimplemented from @see AbstractFilter class
   */
  const QString getHumanLabel() const override;
  /**
   * @brief This method will instantiate all the end user settable options/parameters
   * for this filter
   */
  void setupFilterParameters() override;

  /**
   * @brief Reimplemented from @see AbstractFilter class
   */
  void execute() override;

  /**
   * @brief This function runs some sanity checks on the DataContainer and inputs
   * in an attempt to ensure the filter can process the inputs.
   */
  void preflight() override;

signals:
  /**
   * @brief updateFilterParameters Emitted when the Filter requests all the latest Filter parameters
   * be pushed from a user-facing control (such as a widget)
   * @param filter Filter instance pointer
   */
  void updateFilterParameters(AbstractFilter* filter);

  /**
   * @brief parametersChanged Emitted when any Filter parameter is changed internally
   */
  void parametersChanged();

  /**
   * @brief preflightAboutToExecute Emitted just before calling dataCheck()
   */
  void preflightAboutToExecute();

  /**
   * @brief preflightExecuted Emitted just after calling dataCheck()
   */
  void preflightExecuted();

protected:
  EstablishFoamMorphology();
  /**
   * @brief dataCheck Checks for the appropriate parameter values and availability of arrays
   */
  void dataCheck();

  /**
   * @brief Initializes all the private instance variables.
   */
  void initialize();

  /**
   * @brief initialize_packinggrid Initializes internal variables for packing grid dimensions based on
   * the incoming grid set by the user
   * @return
   */
  Int32ArrayType::Pointer initialize_packinggrid();

  /**
   * @brief place_features Organizes the placement of Features into the packing volume while ensuring
   * proper size and morphological statistics are maintained
   * @param featureOwnersPtr Array of Feature Ids for each packing point
   */
  void place_features(const Int32ArrayType::Pointer& featureOwnersPtr);

  /**
   * @brief generate_feature Creates a Feature by sampling the size and morphological statistical distributions
   * @param phase Index of the Ensemble type for the Feature to be generated
   * @param feature Feature_t struct pointer to be intialized
   * @param shapeclass Type of Feature shape to be generated
   */
  void generate_feature(int32_t phase, Feature_t* feature, uint32_t shapeclass);

  /**
   * @brief transfer_attributes Moves variables held in the Feature_t struct into other arrays
   * @param gnum Id for the Feature to be copied
   * @param feature Feature_t struct pointer to be copied
   */
  void transfer_attributes(int32_t gnum, Feature_t* feature);

  /**
   * @brief insert_feature Performs the insertion of a Feature into the packing volume
   * @param gnum Id for the Feature to be inserted
   */
  void insert_feature(size_t gnum);

  /**
   * @brief move_feature Moves a Feature to the supplied (x,y,z) centroid coordinate
   * @param gnum Id for the Feature to be moved
   * @param xc x centroid coordinate
   * @param yc y centroid coordinate
   * @param zc z centroid coordinate
   */
  void move_feature(size_t gnum, float xc, float yc, float zc);

  /**
   * @brief check_sizedisterror Computes the error between the current Feature size distribution
   * and the goal Feature size distribution
   * @param feature Feature_t struct pointer used to determine the Ensemble type
   * @return Float error value between two distributions
   */
  float check_sizedisterror(Feature_t* feature);

  /**
   * @brief determine_neighbors Determines the neighbors for a given Feature
   * @param gnum Id for the Feature for which to find neighboring Features
   * @param add Value that determines whether to add or remove a Feature from the
   * list of neighbors
   */
  void determine_neighbors(size_t gnum, bool add);

  /**
   * @brief check_neighborhooderror Computes the error between the current Feature neighbor distribution
   * and the goal Feature neighbor distribution
   * @param gadd Value that determines whether to add a Feature for the neighbor list computation
   * @param gremove Value that determines whether to remove a Feature for the neighbor list computation
   * @return Float error value between two distributions
   */
  float check_neighborhooderror(int32_t gadd, int32_t gremove);

  /**
   * @brief check_fillingerror Computes the percentage of unassigned or multiple assigned packing points
   * @param gadd Value that determines whether to add point Ids to be filled
   * @param gremove Value that determines whether to add point Ids to be removed
   * @param featureOwnersPtr Array of Feature Ids for each packing point
   * @param exclusionOwnersPtr Array of exlusion Ids for each packing point
   * @return Float percentage value for the ratio of unassinged/"garbage" packing points
   */
  float check_fillingerror(int32_t gadd, int32_t gremove, const Int32ArrayType::Pointer& featureOwnersPtr, const Int32ArrayType::Pointer& exclusionOwnersPtr);

  /**
   * @brief update_availablepoints Updates the maps used to associate packing points with an "available" state
   * @param availablePoints Map between Feature owners and number of available points
   * @param availablePointsInv Inverse associations for the availablePoints map
   */
  void update_availablepoints(std::map<size_t, size_t>& availablePoints, std::map<size_t, size_t>& availablePointsInv);

  /**
   * @brief assign_voxels Assigns Feature Id values to voxels within the packing grid
   */
  void assign_voxels();

  /**
   * @brief assign_gaps_only Assigns Feature Id values to unassigned gaps within the packing grid
   */
  void assign_gaps_only();

  /**
   * @brief write_goal_attributes Outputs important algorithm data
   */
  void write_goal_attributes();

  /**
   * @brief form_struts sizes the struts and pores correctly
   */
  void form_struts();

  /**
   * @brief find_euclidean_dist_map finds the Euclidean distance map
   */
  void find_euclidean_dist_map();

  /**
   * @brief find_euclideandistmap Provides setup for Euclidean distance map caluclation; note that the
   * actual algorithm is performed in a threaded implementation
   */
  void find_euclideandistmap();

  /**
   * @brief compare_1Ddistributions Computes the 1D Bhattacharyya distance
   * @param sqrerror Float 1D Bhattacharyya distance
   */
  void compare_1Ddistributions(std::vector<float>, std::vector<float>, float& sqrerror);

  /**
   * @brief compare_2Ddistributions Computes the 2D Bhattacharyya distance
   * @param sqrerror Float 1D Bhattacharyya distance
   */
  void compare_2Ddistributions(std::vector<std::vector<float>>, std::vector<std::vector<float>>, float& sqrerror);

  /**
   * @brief compare_3Ddistributions Computes the 3D Bhattacharyya distance
   * @param sqrerror Float 1D Bhattacharyya distance
   */
  void compare_3Ddistributions(std::vector<std::vector<std::vector<float>>>, std::vector<std::vector<std::vector<float>>>, float& sqrerror);

  /**
   * @brief estimate_numfeatures Estimates the number of Features that will be generated based on the supplied statistics
   * @param xpoints Extent of x dimension
   * @param ypoints Extent of y dimension
   * @param zpoints Extent of z dimension
   * @param xres Resolution of x dimension
   * @param yres Resolution of y dimension
   * @param zres Resolution of z dimension
   * @return Integer number of esimtated Features
   */
  int32_t estimate_numfeatures(size_t xpoints, size_t ypoints, size_t zpoints, float xres, float yres, float zres);

private:
  // Names for the arrays used by the packing algorithm
  // These arrays are temporary and are removed from the Feature Attribute Matrix after completion
  QString m_FeaturePhasesArrayName = SIMPL::FeatureData::Phases;
  QString m_NeighborhoodsArrayName = SIMPL::FeatureData::Neighborhoods;
  QString m_NearestNeighborsArrayName = SIMPL::CellData::NearestNeighbors;
  QString m_CentroidsArrayName = SIMPL::FeatureData::Centroids;
  QString m_VolumesArrayName = SIMPL::FeatureData::Volumes;
  QString m_AxisLengthsArrayName = SIMPL::FeatureData::AxisLengths;
  QString m_AxisEulerAnglesArrayName = SIMPL::FeatureData::AxisEulerAngles;
  QString m_Omega3sArrayName = SIMPL::FeatureData::Omega3s;
  QString m_EquivalentDiametersArrayName = SIMPL::FeatureData::EquivalentDiameters;
  int64_t* m_Neighbors = nullptr;

  // Cell Data - make sure these are all initialized to nullptr in the constructor
  DEFINE_DATAARRAY_VARIABLE(int32_t, FeatureIds)
  DEFINE_DATAARRAY_VARIABLE(int32_t, CellPhases)
  DEFINE_DATAARRAY_VARIABLE(bool, Mask)
  DEFINE_DATAARRAY_VARIABLE(float, QPEuclideanDistances)
  DEFINE_DATAARRAY_VARIABLE(float, GBEuclideanDistances)
  DEFINE_DATAARRAY_VARIABLE(float, TJEuclideanDistances)
  DEFINE_DATAARRAY_VARIABLE(int32_t, NearestNeighbors)
  int8_t* m_BoundaryCells;

  // Feature Data - make sure these are all initialized to nullptr in the constructor
  DEFINE_DATAARRAY_VARIABLE(int32_t, FeaturePhases)
  DEFINE_DATAARRAY_VARIABLE(int32_t, Neighborhoods)
  DEFINE_DATAARRAY_VARIABLE(float, Centroids)
  DEFINE_DATAARRAY_VARIABLE(float, Volumes)
  DEFINE_DATAARRAY_VARIABLE(float, AxisLengths)
  DEFINE_DATAARRAY_VARIABLE(float, AxisEulerAngles)
  DEFINE_DATAARRAY_VARIABLE(float, Omega3s)
  DEFINE_DATAARRAY_VARIABLE(float, EquivalentDiameters)
  NeighborList<int>::WeakPointer m_NeighborList;

  // Ensemble Data - make sure these are all initialized to nullptr in the constructor
  DEFINE_DATAARRAY_VARIABLE(uint32_t, PhaseTypes)
  DEFINE_STRINGARRAY_VARIABLE(PhaseNames)
  DEFINE_DATAARRAY_VARIABLE(uint32_t, ShapeTypes)
  DEFINE_DATAARRAY_VARIABLE(int32_t, NumFeatures)
  StatsDataArray::WeakPointer m_StatsDataArray;

  // All other private variables
  QVector<ShapeOps::Pointer> m_ShapeOps;
  ShapeOps::Pointer m_UnknownShapeOps;
  ShapeOps::Pointer m_CubicOctohedronOps;
  ShapeOps::Pointer m_CylinderOps;
  ShapeOps::Pointer m_EllipsoidOps;
  ShapeOps::Pointer m_SuperEllipsoidOps;
  OrthoRhombicOps::Pointer m_OrthoOps;

  std::vector<std::vector<int64_t>> m_ColumnList;
  std::vector<std::vector<int64_t>> m_RowList;
  std::vector<std::vector<int64_t>> m_PlaneList;
  std::vector<std::vector<float>> m_EllipFuncList;

  std::vector<size_t> m_PointsToAdd;
  std::vector<size_t> m_PointsToRemove;

  uint64_t m_Seed;

  int32_t m_FirstFoamFeature;

  float m_SizeX = 0.0f;
  float m_SizeY = 0.0f;
  float m_SizeZ = 0.0f;
  float m_TotalVol = 0.0f;
  FloatVec3Type m_HalfPackingRes;
  FloatVec3Type m_OneOverPackingRes;
  FloatVec3Type m_OneOverHalfPackingRes;
  FloatVec3Type m_PackingRes;

  std::array<int64_t, 3> m_PackingPoints = {0, 0, 0};
  int64_t m_TotalPackingPoints = 0;

  std::vector<std::vector<float>> m_FeatureSizeDist;
  std::vector<std::vector<float>> m_SimFeatureSizeDist;
  std::vector<std::vector<std::vector<float>>> m_NeighborDist;
  std::vector<std::vector<std::vector<float>>> m_SimNeighborDist;

  std::vector<float> m_FeatureSizeDistStep;
  std::vector<float> m_NeighborDistStep;

  std::vector<int64_t> m_PackQualities;
  std::vector<int64_t> m_GSizes;

  std::vector<int32_t> m_PrecipitatePhases;
  std::vector<float> m_PrecipitatePhaseFractions;

  size_t m_AvailablePointsCount = 0;
  float m_FillingError = 0.0f;
  float m_OldFillingError = 0.0f;

  float m_CurrentNeighborhoodError = 0.0f;
  float m_OldNeighborhoodError = 0.0f;
  float m_CurrentSizeDistError = 0.0f;
  float m_OldSizeDistError = 0.0f;

  QString m_ErrorOutputFile = "";

  /**
   * @brief updateFeatureInstancePointers Resets the raw pointers that belong to a
   * Feature Attribute Matrix
   */
  void updateFeatureInstancePointers();

public:
  EstablishFoamMorphology(const EstablishFoamMorphology&) = delete;            // Copy Constructor Not Implemented
  EstablishFoamMorphology(EstablishFoamMorphology&&) = delete;                 // Move Constructor Not Implemented
  EstablishFoamMorphology& operator=(const EstablishFoamMorphology&) = delete; // Copy Assignment Not Implemented
  EstablishFoamMorphology& operator=(EstablishFoamMorphology&&) = delete;      // Move Assignment Not Implemented
};
