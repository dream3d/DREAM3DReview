/* ============================================================================
* Software developed by US federal government employees (including military personnel) 
* as part of their official duties is not subject to copyright protection and is 
* considered "public domain" (see 17 USC Section 105). Public domain software can be used 
* by anyone for any purpose, and cannot be released under a copyright license 
* (including typical open source software licenses).
* 
* This source code file was originally written by United States DoD employees. The
* original source code files are released into the Public Domain.
* 
* Subsequent changes to the codes by others may elect to add a copyright and license
* for those changes.
* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
#ifndef _exportclifile_h_
#define _exportclifile_h_

#include <memory>

#include "SIMPLib/Filtering/AbstractFilter.h"
#include "SIMPLib/SIMPLib.h"
#include "SIMPLib/DataArrays/DataArray.hpp"

/**
 * @brief The ExportCLIFile class. See [Filter documentation](@ref exportclifile) for details.
 */
class ExportCLIFile : public AbstractFilter
{
  Q_OBJECT

public:
  using Self = ExportCLIFile;
  using Pointer = std::shared_ptr<Self>;
  using ConstPointer = std::shared_ptr<const Self>;
  using WeakPointer = std::weak_ptr<Self>;
  using ConstWeakPointer = std::weak_ptr<const Self>;

  static Pointer NullPointer();

  static std::shared_ptr<ExportCLIFile> New();

  /**
   * @brief Returns the name of the class for ExportCLIFile
   */
  QString getNameOfClass() const override;
  /**
   * @brief Returns the name of the class for ExportCLIFile
   */
  static QString ClassName();

  ~ExportCLIFile() override;

  /**
   * @brief Setter property for EdgeGeometry
   */
  void setEdgeGeometry(const DataArrayPath& value);
  /**
   * @brief Getter property for EdgeGeometry
   * @return Value of EdgeGeometry
   */
  DataArrayPath getEdgeGeometry() const;
  Q_PROPERTY(DataArrayPath EdgeGeometry READ getEdgeGeometry WRITE setEdgeGeometry)

  /**
   * @brief Setter property for LayerIdsArrayPath
   */
  void setLayerIdsArrayPath(const DataArrayPath& value);
  /**
   * @brief Getter property for LayerIdsArrayPath
   * @return Value of LayerIdsArrayPath
   */
  DataArrayPath getLayerIdsArrayPath() const;
  Q_PROPERTY(DataArrayPath LayerIdsArrayPath READ getLayerIdsArrayPath WRITE setLayerIdsArrayPath)

  /**
   * @brief Setter property for SplitByGroup
   */
  void setSplitByGroup(bool value);
  /**
   * @brief Getter property for SplitByGroup
   * @return Value of SplitByGroup
   */
  bool getSplitByGroup() const;
  Q_PROPERTY(bool SplitByGroup READ getSplitByGroup WRITE setSplitByGroup)

  /**
   * @brief Setter property for GroupIdsArrayPath
   */
  void setGroupIdsArrayPath(const DataArrayPath& value);
  /**
   * @brief Getter property for GroupIdsArrayPath
   * @return Value of GroupIdsArrayPath
   */
  DataArrayPath getGroupIdsArrayPath() const;
  Q_PROPERTY(DataArrayPath GroupIdsArrayPath READ getGroupIdsArrayPath WRITE setGroupIdsArrayPath)

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
   * @brief Setter property for UnitsScaleFactor
   */
  void setUnitsScaleFactor(double value);
  /**
   * @brief Getter property for UnitsScaleFactor
   * @return Value of UnitsScaleFactor
   */
  double getUnitsScaleFactor() const;
  Q_PROPERTY(double UnitsScaleFactor READ getUnitsScaleFactor WRITE setUnitsScaleFactor)

  /**
   * @brief Setter property for Precision
   */
  void setPrecision(int value);
  /**
   * @brief Getter property for Precision
   * @return Value of Precision
   */
  int getPrecision() const;
  Q_PROPERTY(int Precision READ getPrecision WRITE setPrecision)

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
  ExportCLIFile();

  /**
  * @brief dataCheck Checks for the appropriate parameter values and availability of arrays
  */
  void dataCheck() override;

  /**
  * @brief Initializes all the private instance variables.
  */
  void initialize();

private:
  std::weak_ptr<DataArray<int32_t>> m_LayerIdsPtr;
  int32_t* m_LayerIds = nullptr;
  std::weak_ptr<DataArray<int32_t>> m_GroupIdsPtr;
  int32_t* m_GroupIds = nullptr;

  DataArrayPath m_EdgeGeometry = {};
  DataArrayPath m_LayerIdsArrayPath = {};
  bool m_SplitByGroup = {};
  DataArrayPath m_GroupIdsArrayPath = {};
  QString m_OutputDirectory = {};
  QString m_OutputFilePrefix = {};
  double m_UnitsScaleFactor = {};
  int m_Precision = {};

public:
  ExportCLIFile(const ExportCLIFile&) = delete;            // Copy Constructor Not Implemented
  ExportCLIFile(ExportCLIFile&&) = delete;                 // Move Constructor Not Implemented
  ExportCLIFile& operator=(const ExportCLIFile&) = delete; // Copy Assignment Not Implemented
  ExportCLIFile& operator=(ExportCLIFile&&) = delete;      // Move Assignment Not Implemented
};

#endif /* _exportclifile_h_ */
