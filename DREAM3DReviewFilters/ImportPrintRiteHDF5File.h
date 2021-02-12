#pragma once

#include "SIMPLib/SIMPLib.h"
#include "SIMPLib/DataArrays/DataArray.hpp"
#include "SIMPLib/Filtering/AbstractFilter.h"

#include "DREAM3DReview/DREAM3DReviewPlugin.h"

/**
 * @brief The ImportPrintRiteHDF5File class. See [Filter documentation](@ref importprintritehdf5file) for details.
 */
class DREAM3DReview_EXPORT ImportPrintRiteHDF5File : public AbstractFilter
{
  Q_OBJECT

public:
  using Self = ImportPrintRiteHDF5File;
  using Pointer = std::shared_ptr<Self>;
  using ConstPointer = std::shared_ptr<const Self>;
  using WeakPointer = std::weak_ptr<Self>;
  using ConstWeakPointer = std::weak_ptr<const Self>;
  static Pointer NullPointer();

  static std::shared_ptr<ImportPrintRiteHDF5File> New();

  /**
   * @brief Returns the name of the class for ImportPrintRiteHDF5File
   */
  QString getNameOfClass() const override;

  /**
   * @brief Returns the name of the class for ImportPrintRiteHDF5File
   */
  static QString ClassName();

  ~ImportPrintRiteHDF5File() override;

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
   * @brief Setter property for HFDataContainerName
   */
  void setHFDataContainerName(const QString& value);

  /**
   * @brief Getter property for HFDataContainerName
   * @return Value of HFDataContainerName
   */
  QString getHFDataContainerName() const;
  Q_PROPERTY(QString HFDataContainerName READ getHFDataContainerName WRITE setHFDataContainerName)

  /**
   * @brief Setter property for HFDataName
   */
  void setHFDataName(const QString& value);

  /**
   * @brief Getter property for HFDataName
   * @return Value of HFDataName
   */
  QString getHFDataName() const;
  Q_PROPERTY(QString HFDataName READ getHFDataName WRITE setHFDataName)

  /**
   * @brief Setter property for HFSliceDataName
   */
  void setHFSliceDataName(const QString& value);

  /**
   * @brief Getter property for HFSliceDataName
   * @return Value of HFSliceDataName
   */
  QString getHFSliceDataName() const;
  Q_PROPERTY(QString HFSliceDataName READ getHFSliceDataName WRITE setHFSliceDataName)

  /**
   * @brief Setter property for HFSliceIdsArrayName
   */
  void setHFSliceIdsArrayName(const QString& value);

  /**
   * @brief Getter property for HFSliceIdsArrayName
   * @return Value of HFSliceIdsArrayName
   */
  QString getHFSliceIdsArrayName() const;
  Q_PROPERTY(QString HFSliceIdsArrayName READ getHFSliceIdsArrayName WRITE setHFSliceIdsArrayName)

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
   * @brief readFilterParameters Reimplemented from @see AbstractFilter class
   */
  void readFilterParameters(AbstractFilterParametersReader* reader, int index) override;

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
  void readArray(const DataArrayPath& path, const hid_t& gid, QString& name, size_t index);

  ImportPrintRiteHDF5File();

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
  QString m_HFDataContainerName = {"HighFrequencyDataContainer"};
  QString m_HFDataName = {"HighFrequencyData"};
  QString m_HFSliceDataName = {"HF_SliceAttributeMatrix"};
  QString m_HFSliceIdsArrayName = {"LayerIds"};
  std::weak_ptr<Int32ArrayType> m_HFSliceIdsPtr;
  int32_t* m_HFSliceIds = nullptr;

  QStringList m_SliceGroups;
  QStringList m_SortedSliceGroups;
  QStringList m_HFDataSetNames;
  QStringList m_HFDataSetTypes;
  QStringList m_LFDataSetNames;
  QStringList m_LFDataSetTypes;
  QStringList m_SliceDataSetNames;
  FloatArrayType::Pointer m_XCoords;
  FloatArrayType::Pointer m_YCoords;
  std::map<size_t, size_t> m_NumHFVertsPerSlice;
  std::map<size_t, size_t> m_NumLFVertsPerSlice;
  std::vector<size_t> m_CumulativeNumHFVerts;

  ImportPrintRiteHDF5File(const ImportPrintRiteHDF5File&) = delete; // Copy Constructor Not Implemented
  ImportPrintRiteHDF5File(ImportPrintRiteHDF5File&&) = delete;      // Move Constructor Not Implemented
  void operator=(const ImportPrintRiteHDF5File&) = delete;          // Operator '=' Not Implemented
};
