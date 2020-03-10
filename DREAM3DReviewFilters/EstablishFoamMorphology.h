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

#include "SIMPLib/DataArrays/NeighborList.hpp"
#include "SIMPLib/DataArrays/StatsDataArray.h"
#include "SIMPLib/DataArrays/StringDataArray.h"
#include "SIMPLib/DataContainers/DataContainer.h"
#include "SIMPLib/Filtering/AbstractFilter.h"
#include "SIMPLib/Geometry/ShapeOps/ShapeOps.h"
#include "SIMPLib/DataContainers/DataArrayPath.h"
#include "SIMPLib/SIMPLib.h"

#include "OrientationLib/LaueOps/OrthoRhombicOps.h"

#include "DREAM3DReview/DREAM3DReviewDLLExport.h"

/**
 * @brief The EstablishFoamMorphology class. See [Filter documentation](@ref EstablishFoamMorphology) for details.
 */
class DREAM3DReview_EXPORT EstablishFoamMorphology : public AbstractFilter
{
  Q_OBJECT
#ifdef SIMPL_ENABLE_PYTHON
  PYB11_CREATE_BINDINGS(EstablishFoamMorphology SUPERCLASS AbstractFilter)
  PYB11_SHARED_POINTERS(EstablishFoamMorphology)
  PYB11_FILTER_NEW_MACRO(EstablishFoamMorphology)
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
#endif

public:
  using Self = EstablishFoamMorphology;
  using Pointer = std::shared_ptr<Self>;
  using ConstPointer = std::shared_ptr<const Self>;
  using WeakPointer = std::weak_ptr<Self>;
  using ConstWeakPointer = std::weak_ptr<const Self>;

  static Pointer NullPointer();

  static Pointer New();

  /**
   * @brief Returns the name of the class for EstablishFoamMorphology
   */
  QString getNameOfClass() const override;
  /**
   * @brief Returns the name of the class for EstablishFoamMorphology
   */
  static QString ClassName();

  using Feature_t = struct
  {
    float m_Volumes;
    float m_EquivalentDiameters;
    std::array<float, 3> m_AxisLengths;
    std::array<float, 3> m_AxisEulerAngles;
    float m_Omega3s;
    int32_t m_FeaturePhases;
    int32_t m_Neighborhoods;
  };

  ~EstablishFoamMorphology() override;

  /**
   * @brief Setter property for DataContainerName
   */
  void setDataContainerName(const DataArrayPath& value);
  /**
   * @brief Getter property for DataContainerName
   * @return Value of DataContainerName
   */
  DataArrayPath getDataContainerName() const;

  Q_PROPERTY(DataArrayPath DataContainerName READ getDataContainerName WRITE setDataContainerName)

  /**
   * @brief Setter property for OutputCellAttributeMatrixPath
   */
  void setOutputCellAttributeMatrixPath(const DataArrayPath& value);
  /**
   * @brief Getter property for OutputCellAttributeMatrixPath
   * @return Value of OutputCellAttributeMatrixPath
   */
  DataArrayPath getOutputCellAttributeMatrixPath() const;

  Q_PROPERTY(DataArrayPath OutputCellAttributeMatrixPath READ getOutputCellAttributeMatrixPath WRITE setOutputCellAttributeMatrixPath)

  /**
   * @brief Setter property for InputCellFeatureIdsArrayPath
   */
  void setInputCellFeatureIdsArrayPath(const DataArrayPath& value);
  /**
   * @brief Getter property for InputCellFeatureIdsArrayPath
   * @return Value of InputCellFeatureIdsArrayPath
   */
  DataArrayPath getInputCellFeatureIdsArrayPath() const;

  Q_PROPERTY(DataArrayPath InputCellFeatureIdsArrayPath READ getInputCellFeatureIdsArrayPath WRITE setInputCellFeatureIdsArrayPath)

  /**
   * @brief Setter property for OutputCellEnsembleAttributeMatrixName
   */
  void setOutputCellEnsembleAttributeMatrixName(const QString& value);
  /**
   * @brief Getter property for OutputCellEnsembleAttributeMatrixName
   * @return Value of OutputCellEnsembleAttributeMatrixName
   */
  QString getOutputCellEnsembleAttributeMatrixName() const;

  Q_PROPERTY(QString OutputCellEnsembleAttributeMatrixName READ getOutputCellEnsembleAttributeMatrixName WRITE setOutputCellEnsembleAttributeMatrixName)

  /**
   * @brief Setter property for OutputCellFeatureAttributeMatrixName
   */
  void setOutputCellFeatureAttributeMatrixName(const QString& value);
  /**
   * @brief Getter property for OutputCellFeatureAttributeMatrixName
   * @return Value of OutputCellFeatureAttributeMatrixName
   */
  QString getOutputCellFeatureAttributeMatrixName() const;

  Q_PROPERTY(QString OutputCellFeatureAttributeMatrixName READ getOutputCellFeatureAttributeMatrixName WRITE setOutputCellFeatureAttributeMatrixName)

  /**
   * @brief Setter property for FeatureIdsArrayName
   */
  void setFeatureIdsArrayName(const QString& value);
  /**
   * @brief Getter property for FeatureIdsArrayName
   * @return Value of FeatureIdsArrayName
   */
  QString getFeatureIdsArrayName() const;

  Q_PROPERTY(QString FeatureIdsArrayName READ getFeatureIdsArrayName WRITE setFeatureIdsArrayName)

  /**
   * @brief Setter property for MaskArrayName
   */
  void setMaskArrayName(const QString& value);
  /**
   * @brief Getter property for MaskArrayName
   * @return Value of MaskArrayName
   */
  QString getMaskArrayName() const;

  Q_PROPERTY(QString MaskArrayName READ getMaskArrayName WRITE setMaskArrayName)

  /**
   * @brief Setter property for CellPhasesArrayName
   */
  void setCellPhasesArrayName(const QString& value);
  /**
   * @brief Getter property for CellPhasesArrayName
   * @return Value of CellPhasesArrayName
   */
  QString getCellPhasesArrayName() const;

  Q_PROPERTY(QString CellPhasesArrayName READ getCellPhasesArrayName WRITE setCellPhasesArrayName)

  /**
   * @brief Setter property for NumFeaturesArrayName
   */
  void setNumFeaturesArrayName(const QString& value);
  /**
   * @brief Getter property for NumFeaturesArrayName
   * @return Value of NumFeaturesArrayName
   */
  QString getNumFeaturesArrayName() const;

  Q_PROPERTY(QString NumFeaturesArrayName READ getNumFeaturesArrayName WRITE setNumFeaturesArrayName)

  /**
   * @brief Setter property for QPEuclideanDistancesArrayName
   */
  void setQPEuclideanDistancesArrayName(const QString& value);
  /**
   * @brief Getter property for QPEuclideanDistancesArrayName
   * @return Value of QPEuclideanDistancesArrayName
   */
  QString getQPEuclideanDistancesArrayName() const;

  Q_PROPERTY(QString QPEuclideanDistancesArrayName READ getQPEuclideanDistancesArrayName WRITE setQPEuclideanDistancesArrayName)

  /**
   * @brief Setter property for TJEuclideanDistancesArrayName
   */
  void setTJEuclideanDistancesArrayName(const QString& value);
  /**
   * @brief Getter property for TJEuclideanDistancesArrayName
   * @return Value of TJEuclideanDistancesArrayName
   */
  QString getTJEuclideanDistancesArrayName() const;

  Q_PROPERTY(QString TJEuclideanDistancesArrayName READ getTJEuclideanDistancesArrayName WRITE setTJEuclideanDistancesArrayName)

  /**
   * @brief Setter property for GBEuclideanDistancesArrayName
   */
  void setGBEuclideanDistancesArrayName(const QString& value);
  /**
   * @brief Getter property for GBEuclideanDistancesArrayName
   * @return Value of GBEuclideanDistancesArrayName
   */
  QString getGBEuclideanDistancesArrayName() const;

  Q_PROPERTY(QString GBEuclideanDistancesArrayName READ getGBEuclideanDistancesArrayName WRITE setGBEuclideanDistancesArrayName)

  /**
   * @brief Setter property for InputStatsArrayPath
   */
  void setInputStatsArrayPath(const DataArrayPath& value);
  /**
   * @brief Getter property for InputStatsArrayPath
   * @return Value of InputStatsArrayPath
   */
  DataArrayPath getInputStatsArrayPath() const;

  Q_PROPERTY(DataArrayPath InputStatsArrayPath READ getInputStatsArrayPath WRITE setInputStatsArrayPath)

  /**
   * @brief Setter property for InputPhaseTypesArrayPath
   */
  void setInputPhaseTypesArrayPath(const DataArrayPath& value);
  /**
   * @brief Getter property for InputPhaseTypesArrayPath
   * @return Value of InputPhaseTypesArrayPath
   */
  DataArrayPath getInputPhaseTypesArrayPath() const;

  Q_PROPERTY(DataArrayPath InputPhaseTypesArrayPath READ getInputPhaseTypesArrayPath WRITE setInputPhaseTypesArrayPath)

  /**
   * @brief Setter property for InputPhaseNamesArrayPath
   */
  void setInputPhaseNamesArrayPath(const DataArrayPath& value);
  /**
   * @brief Getter property for InputPhaseNamesArrayPath
   * @return Value of InputPhaseNamesArrayPath
   */
  DataArrayPath getInputPhaseNamesArrayPath() const;

  Q_PROPERTY(DataArrayPath InputPhaseNamesArrayPath READ getInputPhaseNamesArrayPath WRITE setInputPhaseNamesArrayPath)

  /**
   * @brief Setter property for InputShapeTypesArrayPath
   */
  void setInputShapeTypesArrayPath(const DataArrayPath& value);
  /**
   * @brief Getter property for InputShapeTypesArrayPath
   * @return Value of InputShapeTypesArrayPath
   */
  DataArrayPath getInputShapeTypesArrayPath() const;

  Q_PROPERTY(DataArrayPath InputShapeTypesArrayPath READ getInputShapeTypesArrayPath WRITE setInputShapeTypesArrayPath)

  /**
   * @brief Setter property for CsvOutputFile
   */
  void setCsvOutputFile(const QString& value);
  /**
   * @brief Getter property for CsvOutputFile
   * @return Value of CsvOutputFile
   */
  QString getCsvOutputFile() const;

  Q_PROPERTY(QString CsvOutputFile READ getCsvOutputFile WRITE setCsvOutputFile)

  /**
   * @brief Setter property for PeriodicBoundaries
   */
  void setPeriodicBoundaries(const bool& value);
  /**
   * @brief Getter property for PeriodicBoundaries
   * @return Value of PeriodicBoundaries
   */
  bool getPeriodicBoundaries() const;

  Q_PROPERTY(bool PeriodicBoundaries READ getPeriodicBoundaries WRITE setPeriodicBoundaries)

  /**
   * @brief Setter property for WriteGoalAttributes
   */
  void setWriteGoalAttributes(const bool& value);
  /**
   * @brief Getter property for WriteGoalAttributes
   * @return Value of WriteGoalAttributes
   */
  bool getWriteGoalAttributes() const;

  Q_PROPERTY(bool WriteGoalAttributes READ getWriteGoalAttributes WRITE setWriteGoalAttributes)

  /**
   * @brief Setter property for HaveFeatures
   */
  void setHaveFeatures(const int& value);
  /**
   * @brief Getter property for HaveFeatures
   * @return Value of HaveFeatures
   */
  int getHaveFeatures() const;

  Q_PROPERTY(int HaveFeatures READ getHaveFeatures WRITE setHaveFeatures)

  /**
   * @brief Setter property for MinStrutThickness
   */
  void setMinStrutThickness(const double& value);
  /**
   * @brief Getter property for MinStrutThickness
   * @return Value of MinStrutThickness
   */
  double getMinStrutThickness() const;

  Q_PROPERTY(double MinStrutThickness READ getMinStrutThickness WRITE setMinStrutThickness)

  /**
   * @brief Setter property for StrutThicknessVariability
   */
  void setStrutThicknessVariability(const float& value);
  /**
   * @brief Getter property for StrutThicknessVariability
   * @return Value of StrutThicknessVariability
   */
  float getStrutThicknessVariability() const;

  Q_PROPERTY(float StrutThicknessVariability READ getStrutThicknessVariability WRITE setStrutThicknessVariability)

  /**
   * @brief Setter property for StrutShapeVariability
   */
  void setStrutShapeVariability(const float& value);
  /**
   * @brief Getter property for StrutShapeVariability
   * @return Value of StrutShapeVariability
   */
  float getStrutShapeVariability() const;

  Q_PROPERTY(float StrutShapeVariability READ getStrutShapeVariability WRITE setStrutShapeVariability)

  /**
   * @brief Setter property for SmoothStruts
   */
  void setSmoothStruts(const bool& value);
  /**
   * @brief Getter property for SmoothStruts
   * @return Value of SmoothStruts
   */
  bool getSmoothStruts() const;

  Q_PROPERTY(bool SmoothStruts READ getSmoothStruts WRITE setSmoothStruts)

  /**
   * @brief getCompiledLibraryName Reimplemented from @see AbstractFilter class
   */
  QString getCompiledLibraryName() const override;

  /**
   * @brief getBrandingString Returns the branding string for the filter, which is a tag
   * used to denote the filter's association with specific plugins
   * @return Branding string
   */
  QString getBrandingString() const override;

  /**
   * @brief getFilterVersion Returns a version string for this filter. Default
   * value is an empty string.
   * @return
   */
  QString getFilterVersion() const override;

  /**
   * @brief newFilterInstance Reimplemented from @see AbstractFilter class
   */
  AbstractFilter::Pointer newFilterInstance(bool copyFilterParameters) const override;

  /**
   * @brief getGroupName Reimplemented from @see AbstractFilter class
   */
  QString getGroupName() const override;

  /**
   * @brief getSubGroupName Reimplemented from @see AbstractFilter class
   */
  QString getSubGroupName() const override;

  /**
   * @brief getUuid Return the unique identifier for this filter.
   * @return A QUuid object.
   */
  QUuid getUuid() const override;

  /**
   * @brief getHumanLabel Reimplemented from @see AbstractFilter class
   */
  QString getHumanLabel() const override;
  /**
   * @brief This method will instantiate all the end user settable options/parameters
   * for this filter
   */
  void setupFilterParameters() override;

  /**
   * @brief Reimplemented from @see AbstractFilter class
   */
  void execute() override;


protected:
  EstablishFoamMorphology();
  /**
   * @brief dataCheck Checks for the appropriate parameter values and availability of arrays
   */
  void dataCheck() override;

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
  DataArrayPath m_DataContainerName = {};
  DataArrayPath m_OutputCellAttributeMatrixPath = {};
  DataArrayPath m_InputCellFeatureIdsArrayPath = {};
  QString m_OutputCellEnsembleAttributeMatrixName = {};
  QString m_OutputCellFeatureAttributeMatrixName = {};
  QString m_FeatureIdsArrayName = {};
  QString m_MaskArrayName = {};
  QString m_CellPhasesArrayName = {};
  QString m_NumFeaturesArrayName = {};
  QString m_QPEuclideanDistancesArrayName = {};
  QString m_TJEuclideanDistancesArrayName = {};
  QString m_GBEuclideanDistancesArrayName = {};
  DataArrayPath m_InputStatsArrayPath = {};
  DataArrayPath m_InputPhaseTypesArrayPath = {};
  DataArrayPath m_InputPhaseNamesArrayPath = {};
  DataArrayPath m_InputShapeTypesArrayPath = {};
  QString m_CsvOutputFile = {};
  bool m_PeriodicBoundaries = {};
  bool m_WriteGoalAttributes = {};
  int m_HaveFeatures = {};
  double m_MinStrutThickness = {};
  float m_StrutThicknessVariability = {};
  float m_StrutShapeVariability = {};
  bool m_SmoothStruts = {};
  std::weak_ptr<DataArray<int32_t>> m_FeatureIdsPtr;
  int32_t* m_FeatureIds = nullptr;
  std::weak_ptr<DataArray<int32_t>> m_CellPhasesPtr;
  int32_t* m_CellPhases = nullptr;
  std::weak_ptr<DataArray<bool>> m_MaskPtr;
  bool* m_Mask = nullptr;
  std::weak_ptr<DataArray<float>> m_QPEuclideanDistancesPtr;
  float* m_QPEuclideanDistances = nullptr;
  std::weak_ptr<DataArray<float>> m_GBEuclideanDistancesPtr;
  float* m_GBEuclideanDistances = nullptr;
  std::weak_ptr<DataArray<float>> m_TJEuclideanDistancesPtr;
  float* m_TJEuclideanDistances = nullptr;
  std::weak_ptr<DataArray<int32_t>> m_NearestNeighborsPtr;
  int32_t* m_NearestNeighbors = nullptr;
  std::weak_ptr<DataArray<int32_t>> m_FeaturePhasesPtr;
  int32_t* m_FeaturePhases = nullptr;
  std::weak_ptr<DataArray<int32_t>> m_NeighborhoodsPtr;
  int32_t* m_Neighborhoods = nullptr;
  std::weak_ptr<DataArray<float>> m_CentroidsPtr;
  float* m_Centroids = nullptr;
  std::weak_ptr<DataArray<float>> m_VolumesPtr;
  float* m_Volumes = nullptr;
  std::weak_ptr<DataArray<float>> m_AxisLengthsPtr;
  float* m_AxisLengths = nullptr;
  std::weak_ptr<DataArray<float>> m_AxisEulerAnglesPtr;
  float* m_AxisEulerAngles = nullptr;
  std::weak_ptr<DataArray<float>> m_Omega3sPtr;
  float* m_Omega3s = nullptr;
  std::weak_ptr<DataArray<float>> m_EquivalentDiametersPtr;
  float* m_EquivalentDiameters = nullptr;
  std::weak_ptr<DataArray<uint32_t>> m_PhaseTypesPtr;
  uint32_t* m_PhaseTypes = nullptr;
  std::weak_ptr<DataArray<uint32_t>> m_ShapeTypesPtr;
  uint32_t* m_ShapeTypes = nullptr;
  std::weak_ptr<DataArray<int32_t>> m_NumFeaturesPtr;
  int32_t* m_NumFeatures = nullptr;

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
  int8_t* m_BoundaryCells = nullptr;

  // Feature Data - make sure these are all initialized to nullptr in the constructor
  NeighborList<int>::WeakPointer m_NeighborList;

  // Ensemble Data - make sure these are all initialized to nullptr in the constructor
  StringDataArray::WeakPointer m_PhaseNamesPtr;
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
