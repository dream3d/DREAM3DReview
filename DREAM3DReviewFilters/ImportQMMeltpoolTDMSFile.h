#pragma once

#include "SIMPLib/SIMPLib.h"
#include "SIMPLib/DataArrays/DataArray.hpp"
#include "SIMPLib/Filtering/AbstractFilter.h"

#include "DREAM3DReview/DREAM3DReviewPlugin.h"

/**
 * @brief The ImportQMMeltpoolTDMSFile class. See [Filter documentation](@ref importqmmeltpooltdmsfile) for details.
 */
class DREAM3DReview_EXPORT ImportQMMeltpoolTDMSFile : public AbstractFilter
{
  Q_OBJECT

public:
  using Self = ImportQMMeltpoolTDMSFile;
  using Pointer = std::shared_ptr<Self>;
  using ConstPointer = std::shared_ptr<const Self>;
  using WeakPointer = std::weak_ptr<Self>;
  using ConstWeakPointer = std::weak_ptr<const Self>;

  /**
   * @brief Returns a NullPointer wrapped by a shared_ptr<>
   * @return
   */
  static Pointer NullPointer();

  /**
   * @brief Creates a new object wrapped in a shared_ptr<>
   * @return
   */
  static Pointer New();

  /**
   * @brief Returns the name of the class for GenerateRodriguesColors
   */
  QString getNameOfClass() const override;
  /**
   * @brief Returns the name of the class for GenerateRodriguesColors
   */
  static QString ClassName();

  ~ImportQMMeltpoolTDMSFile() override;

  /**
   * @brief Setter property for InputFile
   */
  void setInputFile(const QString& value);

  /**
   * @brief Getter property for InputFile
   * @return Value of InputFile
   */
  QString getInputFile() const;
  Q_PROPERTY(QString InputFile READ getInputFile WRITE setInputFile)

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
   * @brief setupFilterParameters Reimplemented from @see AbstractFilter class
   */
  void setupFilterParameters() override;

  /**
   * @brief execute Reimplemented from @see AbstractFilter class
   */
  void execute() override;

protected:
  ImportQMMeltpoolTDMSFile();

  /**
   * @brief dataCheck Checks for the appropriate parameter values and availability of arrays
   */
  void dataCheck() override;

  /**
   * @brief Initializes all the private instance variables.
   */
  void initialize();

private:
  QString m_InputFile = {""};
  DataArrayPath m_DataContainerName = {"", "", ""};

public:
  /* Rule of 5: All special member functions should be defined if any are defined.
   * CppCoreGuidelines #c21 if you define or delete any default operation define or delete them all
   */
  ImportQMMeltpoolTDMSFile(const ImportQMMeltpoolTDMSFile&) = delete;            // Copy Constructor Not Implemented
  ImportQMMeltpoolTDMSFile& operator=(const ImportQMMeltpoolTDMSFile&) = delete; // Copy Assignment Not Implemented
  ImportQMMeltpoolTDMSFile(ImportQMMeltpoolTDMSFile&&) = delete;                 // Move Constructor Not Implemented
  ImportQMMeltpoolTDMSFile& operator=(ImportQMMeltpoolTDMSFile&&) = delete;      // Move Assignment Not Implemented
};
