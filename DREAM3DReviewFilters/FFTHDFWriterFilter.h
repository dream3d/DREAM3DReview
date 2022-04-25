/*
 * Your License or Copyright can go here
 */

#pragma once

#include <memory>

#include <QtCore/QTextStream>

#include "SIMPLib/CoreFilters/FileWriter.h"
#include "SIMPLib/DataArrays/DataArray.hpp"
#include "SIMPLib/DataArrays/StringDataArray.h"
#include "SIMPLib/Filtering/AbstractFilter.h"
#include "SIMPLib/SIMPLib.h"

#include "DREAM3DReview/DREAM3DReviewDLLExport.h"

/**
 * @brief The FFTHDFWriterFilter class. See [Filter documentation](@ref ffthdfwriterfilter) for details.
 */
class DREAM3DReview_EXPORT FFTHDFWriterFilter : public AbstractFilter
{
  Q_OBJECT

  // Start Python bindings declarations
  PYB11_CREATE_BINDINGS(FFTHDFWriterFilter SUPERCLASS AbstractFilter)
  PYB11_PROPERTY(QString OutputFile READ getOutputFile WRITE setOutputFile)
  PYB11_PROPERTY(bool WriteEigenstrains READ getWriteEigenstrains WRITE setWriteEigenstrains)
  PYB11_PROPERTY(QString EigenstrainsOutputFile READ getEigenstrainsOutputFile WRITE setEigenstrainsOutputFile)
  PYB11_PROPERTY(DataArrayPath FeatureIdsArrayPath READ getFeatureIdsArrayPath WRITE setFeatureIdsArrayPath)
  PYB11_PROPERTY(DataArrayPath CellPhasesArrayPath READ getCellPhasesArrayPath WRITE setCellPhasesArrayPath)
  PYB11_PROPERTY(DataArrayPath CellEulerAnglesArrayPath READ getCellEulerAnglesArrayPath WRITE setCellEulerAnglesArrayPath)
  PYB11_PROPERTY(DataArrayPath CellEigenstrainsArrayPath READ getCellEigenstrainsArrayPath WRITE setCellEigenstrainsArrayPath)
  // End Python bindings declarations

public:
  using Self = FFTHDFWriterFilter;
  using Pointer = std::shared_ptr<Self>;
  using ConstPointer = std::shared_ptr<const Self>;
  using WeakPointer = std::weak_ptr<Self>;
  using ConstWeakPointer = std::weak_ptr<const Self>;
  static Pointer NullPointer();

  static std::shared_ptr<FFTHDFWriterFilter> New();

  /**
   * @brief Returns the name of the class for FFTHDFWriterFilter
   */
  const QString getNameOfClass() const override;
  /**
   * @brief Returns the name of the class for FFTHDFWriterFilter
   */
  static QString ClassName();

  ~FFTHDFWriterFilter() override;

  /**
   * @brief Setter property for OutputFile
   */
  void setOutputFile(const QString& value);
  /**
   * @brief Getter property for OutputFile
   * @return Value of OutputFile
   */
  QString getOutputFile() const;
  Q_PROPERTY(QString OutputFile READ getOutputFile WRITE setOutputFile)

  /**
   * @brief Setter property for EigenstrainsOutputFile
   */
  void setEigenstrainsOutputFile(const QString& value);
  /**
   * @brief Getter property for EigenstrainsOutputFile
   * @return Value of EigenstrainsOutputFile
   */
  QString getEigenstrainsOutputFile() const;
  Q_PROPERTY(QString EigenstrainsOutputFile READ getEigenstrainsOutputFile WRITE setEigenstrainsOutputFile)


  /**
   * @brief Setter property for WriteEigenstrains
   */
  void setWriteEigenstrains(bool value);
  /**
   * @brief Getter property for WriteEigenstrains
   * @return Value of WriteEigenstrains
   */
  bool getWriteEigenstrains() const;
  Q_PROPERTY(bool WriteEigenstrains READ getWriteEigenstrains WRITE setWriteEigenstrains)

  //-------------------------------------------------------------------

  /**
   * @brief Setter property for FeatureIdsArrayPath
   */
  void setFeatureIdsArrayPath(const DataArrayPath& value);
  /**
   * @brief Getter property for FeatureIdsArrayPath
   * @return Value of FeatureIdsArrayPath
   */
  DataArrayPath getFeatureIdsArrayPath() const;
  Q_PROPERTY(DataArrayPath FeatureIdsArrayPath READ getFeatureIdsArrayPath WRITE setFeatureIdsArrayPath)

  /**
   * @brief Setter property for CellPhasesArrayPath
   */
  void setCellPhasesArrayPath(const DataArrayPath& value);
  /**
   * @brief Getter property for CellPhasesArrayPath
   * @return Value of CellPhasesArrayPath
   */
  DataArrayPath getCellPhasesArrayPath() const;
  Q_PROPERTY(DataArrayPath CellPhasesArrayPath READ getCellPhasesArrayPath WRITE setCellPhasesArrayPath)

  /**
   * @brief Setter property for CellEulerAnglesArrayPath
   */
  void setCellEulerAnglesArrayPath(const DataArrayPath& value);
  /**
   * @brief Getter property for CellEulerAnglesArrayPath
   * @return Value of CellEulerAnglesArrayPath
   */
  DataArrayPath getCellEulerAnglesArrayPath() const;
  Q_PROPERTY(DataArrayPath CellEulerAnglesArrayPath READ getCellEulerAnglesArrayPath WRITE setCellEulerAnglesArrayPath)

  /**
   * @brief Setter property for CellEigenstrainsArrayPath
   */
  void setCellEigenstrainsArrayPath(const DataArrayPath& value);
  /**
   * @brief Getter property for CellEigenstrainsArrayPath
   * @return Value of CellEigenstrainsArrayPath
   */
  DataArrayPath getCellEigenstrainsArrayPath() const;
  Q_PROPERTY(DataArrayPath CellEigenstrainsArrayPath READ getCellEigenstrainsArrayPath WRITE setCellEigenstrainsArrayPath)

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
  const QUuid getUuid();

  /**
   * @brief getHumanLabel Reimplemented from @see AbstractFilter class
   */
  const QString getHumanLabel() const override;

  /**
   * @brief setupFilterParameters Reimplemented from @see AbstractFilter class
   */
  void setupFilterParameters() override;

  /**
   * @brief readFilterParameters Reimplemented from @see AbstractFilter class
   */
  void readFilterParameters(AbstractFilterParametersReader* reader, int index) override;

  /**
   * @brief execute Reimplemented from @see AbstractFilter class
   */
  void execute() override;

  /**
  * @brief preflight Reimplemented from @see AbstractFilter class
  */
  void preflight() override;

Q_SIGNALS:
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
  FFTHDFWriterFilter();

  /**
   * @brief dataCheck Checks for the appropriate parameter values and availability of arrays
   */
  void dataCheck();

  /**
   * @brief openFile Opens or creates an HDF5 file to write data into
   * @param append Should a new file be created or append data to a currently existing file
   * @return
   */
  void openFile(QString file, hid_t& fileId, bool append = false);

private:
  std::weak_ptr<Int32ArrayType> m_FeatureIdsPtr;
  int32_t* m_FeatureIds = nullptr;
  std::weak_ptr<Int32ArrayType> m_CellPhasesPtr;
  int32_t* m_CellPhases = nullptr;
  std::weak_ptr<FloatArrayType> m_CellEulerAnglesPtr;
  float* m_CellEulerAngles = nullptr;
  std::weak_ptr<FloatArrayType> m_CellEigenstrainsPtr;
  float* m_CellEigenstrains = nullptr;

  QString m_OutputFile = {""};
  QString m_EigenstrainsOutputFile = {""};
  bool m_WritePipeline = {true};
  bool m_WriteEigenstrains = {false};
  bool m_AppendToExisting = {false};
  DataArrayPath m_FeatureIdsArrayPath = {SIMPL::Defaults::ImageDataContainerName, SIMPL::Defaults::CellAttributeMatrixName, SIMPL::CellData::FeatureIds};
  DataArrayPath m_CellPhasesArrayPath = {SIMPL::Defaults::ImageDataContainerName, SIMPL::Defaults::CellAttributeMatrixName, SIMPL::CellData::Phases};
  DataArrayPath m_CellEulerAnglesArrayPath = {SIMPL::Defaults::ImageDataContainerName, SIMPL::Defaults::CellAttributeMatrixName, SIMPL::CellData::EulerAngles};
  DataArrayPath m_CellEigenstrainsArrayPath = {SIMPL::Defaults::ImageDataContainerName, SIMPL::Defaults::CellAttributeMatrixName, "Eigenstrains"};

  hid_t m_FileId = -1;
  hid_t m_FileIdEig = -2;

public:
  FFTHDFWriterFilter(const FFTHDFWriterFilter&) = delete;            // Copy Constructor Not Implemented
  FFTHDFWriterFilter(FFTHDFWriterFilter&&) = delete;                 // Move Constructor Not Implemented
  FFTHDFWriterFilter& operator=(const FFTHDFWriterFilter&) = delete; // Copy Assignment Not Implemented
  FFTHDFWriterFilter& operator=(FFTHDFWriterFilter&&) = delete;      // Move Assignment Not Implemented
};
