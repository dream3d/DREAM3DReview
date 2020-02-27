/* ============================================================================
* Software developed by US federal government employees (including military personnel) 
* as part of their official duties is not subject to copyright protection and is 
* considered “public domain” (see 17 USC Section 105). Public domain software can be used 
* by anyone for any purpose, and cannot be released under a copyright license 
* (including typical open source software licenses).
* 
* This source code file was originally written by United States DoD employees. The
* original source code files are released into the Public Domain.
* 
* Subsequent changes to the codes by others may elect to add a copyright and license
* for those changes.
* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
#ifndef _importclifile_h_
#define _importclifile_h_

#include <memory>

#include "SIMPLib/Filtering/AbstractFilter.h"
#include "SIMPLib/SIMPLib.h"

/**
 * @brief The ImportCLIFile class. See [Filter documentation](@ref importclifile) for details.
 */
class ImportCLIFile : public AbstractFilter
{
  Q_OBJECT

public:
  using Self = ImportCLIFile;
  using Pointer = std::shared_ptr<Self>;
  using ConstPointer = std::shared_ptr<const Self>;
  using WeakPointer = std::weak_ptr<Self>;
  using ConstWeakPointer = std::weak_ptr<const Self>;
  static Pointer NullPointer();

  static std::shared_ptr<ImportCLIFile> New();

  /**
   * @brief Returns the name of the class for ImportCLIFile
   */
  QString getNameOfClass() const override;
  /**
   * @brief Returns the name of the class for ImportCLIFile
   */
  static QString ClassName();

  ~ImportCLIFile() override;

  /**
   * @brief Setter property for CLIFile
   */
  void setCLIFile(const QString& value);
  /**
   * @brief Getter property for CLIFile
   * @return Value of CLIFile
   */
  QString getCLIFile() const;

  Q_PROPERTY(QString CLIFile READ getCLIFile WRITE setCLIFile)

  /**
   * @brief Setter property for EdgeDataContainerName
   */
  void setEdgeDataContainerName(const QString& value);
  /**
   * @brief Getter property for EdgeDataContainerName
   * @return Value of EdgeDataContainerName
   */
  QString getEdgeDataContainerName() const;

  Q_PROPERTY(QString EdgeDataContainerName READ getEdgeDataContainerName WRITE setEdgeDataContainerName)

  /**
   * @brief Setter property for VertexAttributeMatrixName
   */
  void setVertexAttributeMatrixName(const QString& value);
  /**
   * @brief Getter property for VertexAttributeMatrixName
   * @return Value of VertexAttributeMatrixName
   */
  QString getVertexAttributeMatrixName() const;

  Q_PROPERTY(QString VertexAttributeMatrixName READ getVertexAttributeMatrixName WRITE setVertexAttributeMatrixName)

  /**
   * @brief Setter property for EdgeAttributeMatrixName
   */
  void setEdgeAttributeMatrixName(const QString& value);
  /**
   * @brief Getter property for EdgeAttributeMatrixName
   * @return Value of EdgeAttributeMatrixName
   */
  QString getEdgeAttributeMatrixName() const;

  Q_PROPERTY(QString EdgeAttributeMatrixName READ getEdgeAttributeMatrixName WRITE setEdgeAttributeMatrixName)

  /**
   * @brief Setter property for LayerIdsArrayName
   */
  void setLayerIdsArrayName(const QString& value);
  /**
   * @brief Getter property for LayerIdsArrayName
   * @return Value of LayerIdsArrayName
   */
  QString getLayerIdsArrayName() const;

  Q_PROPERTY(QString LayerIdsArrayName READ getLayerIdsArrayName WRITE setLayerIdsArrayName)

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
  * @brief preflight Reimplemented from @see AbstractFilter class
  */
  void preflight() override;

  /**
  * @brief getUuid Return the unique identifier for this filter.
  * @return A QUuid object.
  */
  QUuid getUuid() const override;

protected:
  ImportCLIFile();

  /**
  * @brief dataCheck Checks for the appropriate parameter values and availability of arrays
  */
  void dataCheck();

  /**
  * @brief Initializes all the private instance variables.
  */
  void initialize();

private:
  QString m_CLIFile = {};
  QString m_EdgeDataContainerName = {};
  QString m_VertexAttributeMatrixName = {};
  QString m_EdgeAttributeMatrixName = {};
  QString m_LayerIdsArrayName = {};
  QString m_FeatureIdsArrayName = {};

  ImportCLIFile(const ImportCLIFile&) = delete;  // Copy Constructor Not Implemented
  ImportCLIFile(ImportCLIFile&&) = delete;       // Move Constructor Not Implemented
  void operator=(const ImportCLIFile&) = delete; // Operator '=' Not Implemented
};

#endif /* _ImportCLIFile_H_ */
