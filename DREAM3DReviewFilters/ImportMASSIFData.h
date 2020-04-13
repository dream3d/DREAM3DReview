/*
 * Your License or Copyright can go here
 */

#pragma once

#include <memory>

#include "SIMPLib/SIMPLib.h"
#include "SIMPLib/Common/SIMPLArray.hpp"
#include "SIMPLib/Filtering/AbstractFilter.h"
#include "SIMPLib/DataArrays/DataArray.hpp"

class IDataArray;
using IDataArrayShPtrType = std::shared_ptr<IDataArray>;

#include "DREAM3DReview/DREAM3DReviewDLLExport.h"

/**
 * @brief The ImportMASSIFData class. See [Filter documentation](@ref importmassifdata) for details.
 */
class DREAM3DReview_EXPORT ImportMASSIFData : public AbstractFilter
{
  Q_OBJECT

  // Start Python bindings declarations
  PYB11_BEGIN_BINDINGS(ImportMASSIFData SUPERCLASS AbstractFilter)
  PYB11_FILTER()
  PYB11_SHARED_POINTERS(ImportMASSIFData)
  PYB11_FILTER_NEW_MACRO(ImportMASSIFData)
  PYB11_PROPERTY(QString MassifInputFilePath READ getMassifInputFilePath WRITE setMassifInputFilePath)
  PYB11_PROPERTY(QString FilePrefix READ getFilePrefix WRITE setFilePrefix)
  PYB11_PROPERTY(int StepNumber READ getStepNumber WRITE setStepNumber)
  PYB11_END_BINDINGS()
  // End Python bindings declarations

public:
  using Self = ImportMASSIFData;
  using Pointer = std::shared_ptr<Self>;
  using ConstPointer = std::shared_ptr<const Self>;
  using WeakPointer = std::weak_ptr<Self>;
  using ConstWeakPointer = std::weak_ptr<const Self>;
  static Pointer NullPointer();

  static std::shared_ptr<ImportMASSIFData> New();

  /**
   * @brief Returns the name of the class for ImportMASSIFData
   */
  QString getNameOfClass() const override;
  /**
   * @brief Returns the name of the class for ImportMASSIFData
   */
  static QString ClassName();

  ~ImportMASSIFData() override;

  /**
   * @brief Setter property for MassifInputFilePath
   */
  void setMassifInputFilePath(const QString& value);
  /**
   * @brief Getter property for MassifInputFilePath
   * @return Value of MassifInputFilePath
   */
  QString getMassifInputFilePath() const;
  Q_PROPERTY(QString MassifInputFilePath READ getMassifInputFilePath WRITE setMassifInputFilePath)

  /**
   * @brief Setter property for FilePrefix
   */
  void setFilePrefix(const QString& value);
  /**
   * @brief Getter property for FilePrefix
   * @return Value of FilePrefix
   */
  QString getFilePrefix() const;
  Q_PROPERTY(QString FilePrefix READ getFilePrefix WRITE setFilePrefix)

  /**
   * @brief Setter property for StepNumber
   */
  void setStepNumber(int value);
  /**
   * @brief Getter property for StepNumber
   * @return Value of StepNumber
   */
  int getStepNumber() const;
  Q_PROPERTY(int StepNumber READ getStepNumber WRITE setStepNumber)

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
  ImportMASSIFData();

  /**
   * @brief dataCheck Checks for the appropriate parameter values and availability of arrays
   */
  void dataCheck() override;

  /**
   * @brief Initializes all the private instance variables.
   */
  void initialize();

private:
  std::weak_ptr<DataArray<float>> m_DFieldPtr;
  std::weak_ptr<DataArray<float>> m_EFieldPtr;
  std::weak_ptr<DataArray<float>> m_SFieldPtr;

  QString m_MassifInputFilePath = {};
  QString m_FilePrefix = {};
  int m_StepNumber = {};

  QString m_PaddedStep = "";

  /**
   * @brief createDataArrayPaths
   * @return
   */
  QVector<DataArrayPath> createDataArrayPaths();

  /**
   * @brief createHDF5DatasetPaths
   * @return
   */
  QVector<QString> createHDF5DatasetPaths();

  /**
   * @brief getDataContainerGeometry
   * @param tDims
   * @param origin
   * @param res
   */
  void getDataContainerGeometry(std::vector<size_t>& tDims, FloatVec3Type& origin, FloatVec3Type& spacing);

  /**
   * @brief readIDataArray
   * @param gid
   * @param name
   * @param metaDataOnly
   * @return
   */
  IDataArray::Pointer readIDataArray(hid_t gid, const QString& name, std::vector<size_t> geoDims, bool metaDataOnly);

public:
  ImportMASSIFData(const ImportMASSIFData&) = delete;            // Copy Constructor Not Implemented
  ImportMASSIFData(ImportMASSIFData&&) = delete;                 // Move Constructor Not Implemented
  ImportMASSIFData& operator=(const ImportMASSIFData&) = delete; // Copy Assignment Not Implemented
  ImportMASSIFData& operator=(ImportMASSIFData&&) = delete;      // Move Assignment Not Implemented
};
