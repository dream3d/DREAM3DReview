#pragma once

#include <array>

#include "SIMPLib/SIMPLib.h"
#include "SIMPLib/DataArrays/DataArray.hpp"
#include "SIMPLib/DataContainers/DataContainerArray.h"
#include "SIMPLib/FilterParameters/FileListInfoFilterParameter.h"
#include "SIMPLib/FilterParameters/FloatVec2FilterParameter.h"
#include "SIMPLib/Filtering/AbstractFilter.h"

#include "DREAM3DReview/DREAM3DReviewFilters/util/PrintRiteHelpers.h"
#include "DREAM3DReview/DREAM3DReviewPlugin.h"

/**
 * @brief The ImportPrintRiteTDMSFiles class. See [Filter documentation](@ref importprintritetdmsfiles) for details.
 */
class DREAM3DReview_EXPORT ImportPrintRiteTDMSFiles : public AbstractFilter
{
  Q_OBJECT

public:
  using Self = ImportPrintRiteTDMSFiles;
  using Pointer = std::shared_ptr<Self>;
  using ConstPointer = std::shared_ptr<const Self>;
  using WeakPointer = std::weak_ptr<Self>;
  using ConstWeakPointer = std::weak_ptr<const Self>;
  static Pointer NullPointer();

  static std::shared_ptr<ImportPrintRiteTDMSFiles> New();

  /**
   * @brief Returns the name of the class for ImportPrintRiteTDMSFiles
   */
  QString getNameOfClass() const override;

  /**
   * @brief Returns the name of the class for ImportPrintRiteTDMSFiles
   */
  static QString ClassName();

  ~ImportPrintRiteTDMSFiles() override;

  /**
   * @brief Setter property for InputFilesList
   */
  void setInputFilesList(const StackFileListInfo& value);

  /**
   * @brief Getter property for InputFilesList
   * @return Value of InputFilesList
   */
  StackFileListInfo getInputFilesList() const;
  Q_PROPERTY(StackFileListInfo InputFilesList READ getInputFilesList WRITE setInputFilesList)

  /**
   * @brief Setter property for STLFilePath1
   */
  void setSTLFilePath1(const QString& value);

  /**
   * @brief Getter property for STLFilePath1
   * @return Value of STLFilePath1
   */
  QString getSTLFilePath1() const;
  Q_PROPERTY(QString STLFilePath1 READ getSTLFilePath1 WRITE setSTLFilePath1)

  /**
   * @brief Setter property for STLFilePath2
   */
  void setSTLFilePath2(const QString& value);

  /**
   * @brief Getter property for STLFilePath2
   * @return Value of STLFilePath2
   */
  QString getSTLFilePath2() const;
  Q_PROPERTY(QString STLFilePath2 READ getSTLFilePath2 WRITE setSTLFilePath2)

  /**
   * @brief Setter property for LayerThickness
   */
  void setLayerThickness(const float& value);

  /**
   * @brief Getter property for LayerThickness
   * @return Value of LayerThickness
   */
  float getLayerThickness() const;
  Q_PROPERTY(float LayerThickness READ getLayerThickness WRITE setLayerThickness)

  /**
   * @brief Setter property for LaserOnThreshold
   */
  void setLaserOnThreshold(const float& value);

  /**
   * @brief Getter property for LaserOnThreshold
   * @return Value of LaserOnThreshold
   */
  float getLaserOnThreshold() const;
  Q_PROPERTY(float LaserOnThreshold READ getLaserOnThreshold WRITE setLaserOnThreshold)

  /**
   * @brief Setter property for LaserOnArrayOption
   */
  void setLaserOnArrayOption(const int& value);

  /**
   * @brief Getter property for LaserOnArrayOption
   * @return Value of LaserOnArrayOption
   */
  int getLaserOnArrayOption() const;
  Q_PROPERTY(int LaserOnArrayOption READ getLaserOnArrayOption WRITE setLaserOnArrayOption)

  /**
   * @brief Setter property for OutputDirectory
   */
  void setOutputDirectory(const QString& value);

  /**
   * @brief Getter property for OutputDirectory
   * @return Value of OutputDirectory
   */
  QString getOutputDirectory() const;
  Q_PROPERTY(QString OutputDirectory READ getOutputDirectory WRITE setOutputDirectory)

  /**
   * @brief Setter property for OutputFilePrefix
   */
  void setOutputFilePrefix(const QString& value);

  /**
   * @brief Getter property for OutputFilePrefix
   * @return Value of OutputFilePrefix
   */
  QString getOutputFilePrefix() const;
  Q_PROPERTY(QString OutputFilePrefix READ getOutputFilePrefix WRITE setOutputFilePrefix)

  /**
   * @brief Setter property for SpatialTransformOption
   */
  void setSpatialTransformOption(const int& value);

  /**
   * @brief Getter property for SpatialTransformOption
   * @return Value of SpatialTransformOption
   */
  int getSpatialTransformOption() const;
  Q_PROPERTY(int SpatialTransformOption READ getSpatialTransformOption WRITE setSpatialTransformOption)

  /**
   * @brief Setter property for DowncastRawData
   */
  void setDowncastRawData(const bool& value);

  /**
   * @brief Getter property for DowncastRawData
   * @return Value of DowncastRawData
   */
  bool getDowncastRawData() const;
  Q_PROPERTY(bool DowncastRawData READ getDowncastRawData WRITE setDowncastRawData)

  /**
   * @brief Setter property for ScaleLaserPower
   */
  void setScaleLaserPower(const bool& value);

  /**
   * @brief Getter property for ScaleLaserPower
   * @return Value of ScaleLaserPower
   */
  bool getScaleLaserPower() const;
  Q_PROPERTY(bool ScaleLaserPower READ getScaleLaserPower WRITE setScaleLaserPower)

  /**
   * @brief Setter property for PowerScalingCoefficients
   */
  void setPowerScalingCoefficients(const FloatVec2Type& value);

  /**
   * @brief Getter property for PowerScalingCoefficients
   * @return Value of PowerScalingCoefficients
   */
  FloatVec2Type getPowerScalingCoefficients() const;
  Q_PROPERTY(FloatVec2Type PowerScalingCoefficients READ getPowerScalingCoefficients WRITE setPowerScalingCoefficients)

  /**
   * @brief Setter property for ScalePyrometerTemperature
   */
  void setScalePyrometerTemperature(const bool& value);

  /**
   * @brief Getter property for ScalePyrometerTemperature
   * @return Value of ScalePyrometerTemperature
   */
  bool getScalePyrometerTemperature() const;
  Q_PROPERTY(bool ScalePyrometerTemperature READ getScalePyrometerTemperature WRITE setScalePyrometerTemperature)

  /**
   * @brief Setter property for TemperatureScalingCoefficients
   */
  void setTemperatureScalingCoefficients(const FloatVec2Type& value);

  /**
   * @brief Getter property for TemperatureScalingCoefficients
   * @return Value of TemperatureScalingCoefficients
   */
  FloatVec2Type getTemperatureScalingCoefficients() const;
  Q_PROPERTY(FloatVec2Type TemperatureScalingCoefficients READ getTemperatureScalingCoefficients WRITE setTemperatureScalingCoefficients)

  /**
   * @brief Setter property for SplitRegions1
   */
  void setSplitRegions1(const bool& value);

  /**
   * @brief Getter property for SplitRegions1
   * @return Value of SplitRegions1
   */
  bool getSplitRegions1() const;
  Q_PROPERTY(bool SplitRegions1 READ getSplitRegions1 WRITE setSplitRegions1)

  /**
   * @brief Setter property for SplitRegions2
   */
  void setSplitRegions2(const bool& value);

  /**
   * @brief Getter property for SplitRegions2
   * @return Value of SplitRegions2
   */
  bool getSplitRegions2() const;
  Q_PROPERTY(bool SplitRegions2 READ getSplitRegions2 WRITE setSplitRegions2)

  /**
   * @brief Setter property for LayerForScaling
   */
  void setLayerForScaling(const int& value);

  /**
   * @brief Getter property for LayerForScaling
   * @return Value of LayerForScaling
   */
  int getLayerForScaling() const;
  Q_PROPERTY(int LayerForScaling READ getLayerForScaling WRITE setLayerForScaling)

  /**
   * @brief Setter property for InputSpatialTransformFilePath
   */
  void setInputSpatialTransformFilePath(const QString& value);

  /**
   * @brief Getter property for InputSpatialTransformFilePath
   * @return Value of InputSpatialTransformFilePath
   */
  QString getInputSpatialTransformFilePath() const;
  Q_PROPERTY(QString InputSpatialTransformFilePath READ getInputSpatialTransformFilePath WRITE setInputSpatialTransformFilePath)

  /**
   * @brief Setter property for SearchRadius
   */
  void setSearchRadius(const float& value);

  /**
   * @brief Getter property for SearchRadius
   * @return Value of SearchRadius
   */
  float getSearchRadius() const;
  Q_PROPERTY(float SearchRadius READ getSearchRadius WRITE setSearchRadius)

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
   * @brief getHumanLabel Reimplemented from @see AbstractFilter class
   */
  QString getHumanLabel() const override;

  /**
   * @brief setupFilterParameters Reimplemented from @see AbstractFilter class
   */
  void setupFilterParameters() override;

  /**
   * @brief execute Reimplemented from @see AbstractFilter class
   */
  void execute() override;

  /**
   * @brief getUuid Return the unique identifier for this filter.
   * @return A QUuid object.
   */
  QUuid getUuid() const override;

protected:
  ImportPrintRiteTDMSFiles();

  /**
   * @brief dataCheck Checks for the appropriate parameter values and availability of arrays
   */
  void dataCheck() override;

  /**
   * @brief Initializes all the private instance variables.
   */
  void initialize();

  void importLabelSliceSTL();

  void createHDF5Files();

  void closeHDF5Files();

  void writeMetaDataToHDF5Files();

  void writeSpatialTransformToFile();

  void computeSpatialTransformation(const QString& fname);

  void processLayers(const QVector<QString>& files);

  void determinePointsForLeastSquares(const PrintRiteHelpers::Polygons& polygons, FloatArrayType::Pointer tdms, std::pair<Int32ArrayType::Pointer, std::vector<bool>> pointsToPolys,
                                      BoolArrayType::Pointer mask);

  std::vector<std::vector<float>> findPolygonKeyPoints(const PrintRiteHelpers::Polygons& polygons);

  PrintRiteHelpers::Polygons extractLayerPolygons(int32_t layer);

  std::pair<Int32ArrayType::Pointer, std::vector<bool>> associatePointsWithPolygons(const PrintRiteHelpers::Polygons& polygons, FloatArrayType::Pointer tdms);

  std::pair<std::vector<float>, PrintRiteHelpers::BoundingBox> findCentralCluster(FloatArrayType::Pointer tdms, Int32ArrayType::Pointer clusterIds);

  std::pair<std::vector<float>, PrintRiteHelpers::BoundingBox> findCentralPolygon(const PrintRiteHelpers::Polygons& polygons, int32_t layer, const QString& fname);

private:
  StackFileListInfo m_InputFilesList = {};
  QString m_STLFilePath1 = {};
  QString m_STLFilePath2 = {};
  float m_LayerThickness = {0.02f};
  float m_LaserOnThreshold = {2.5f};
  int m_LaserOnArrayOption = {0};
  QString m_OutputDirectory = {};
  QString m_OutputFilePrefix = {"Build_"};
  int m_SpatialTransformOption = {0};
  bool m_DowncastRawData = {true};
  bool m_ScaleLaserPower = {false};
  FloatVec2Type m_PowerScalingCoefficients = {};
  bool m_ScalePyrometerTemperature = {false};
  FloatVec2Type m_TemperatureScalingCoefficients = {};
  bool m_SplitRegions1 = {false};
  bool m_SplitRegions2 = {false};
  int m_LayerForScaling = {0};
  QString m_InputSpatialTransformFilePath = {};
  float m_SearchRadius = {0.5f};

  int32_t m_NumParts = 1;
  int32_t m_NumLayers = 0;
  int32_t m_NumLayersToImport = 0;
  int32_t m_Offset = 0;
  QString m_STLFilePath;
  DataContainerArray::Pointer m_LocalStructure;
  std::map<int32_t, hid_t> m_RegionFileMap;
  std::vector<int32_t> m_PolygonsWithoutPoints;
  PrintRiteHelpers::Polynomial m_Polynomial;
  std::string m_LaserOnArrayName;

  ImportPrintRiteTDMSFiles(const ImportPrintRiteTDMSFiles&) = delete; // Copy Constructor Not Implemented
  ImportPrintRiteTDMSFiles(ImportPrintRiteTDMSFiles&&) = delete;      // Move Constructor Not Implemented
  void operator=(const ImportPrintRiteTDMSFiles&) = delete;           // Operator '=' Not Implemented
};
