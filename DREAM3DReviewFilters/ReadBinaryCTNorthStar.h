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
#ifndef _readbinaryctnorthstar_h_
#define _readbinaryctnorthstar_h_

#include <memory>

#include <QtCore/QFile>

#include "SIMPLib/Filtering/AbstractFilter.h"
#include "SIMPLib/Geometry/ImageGeom.h"
#include "SIMPLib/SIMPLib.h"
#include "SIMPLib/DataArrays/DataArray.hpp"

/**
 * @brief The ReadBinaryCTNorthStar class. See [Filter documentation](@ref readbinaryctnorthstar) for details.
 */
class ReadBinaryCTNorthStar : public AbstractFilter
{
  Q_OBJECT

public:
  using Self = ReadBinaryCTNorthStar;
  using Pointer = std::shared_ptr<Self>;
  using ConstPointer = std::shared_ptr<const Self>;
  using WeakPointer = std::weak_ptr<Self>;
  using ConstWeakPointer = std::weak_ptr<const Self>;
  static Pointer NullPointer();

  static std::shared_ptr<ReadBinaryCTNorthStar> New();

  /**
   * @brief Returns the name of the class for ReadBinaryCTNorthStar
   */
  QString getNameOfClass() const override;
  /**
   * @brief Returns the name of the class for ReadBinaryCTNorthStar
   */
  static QString ClassName();

  ~ReadBinaryCTNorthStar() override;

  /**
   * @brief Setter property for InputFiles
   */
  void setInputFiles(const std::vector<QString>& value);
  /**
   * @brief Getter property for InputFiles
   * @return Value of InputFiles
   */
  std::vector<QString> getInputFiles() const;

  /**
   * @brief Setter property for SlicesPerFile
   */
  void setSlicesPerFile(const std::vector<int64_t>& value);
  /**
   * @brief Getter property for SlicesPerFile
   * @return Value of SlicesPerFile
   */
  std::vector<int64_t> getSlicesPerFile() const;

  /**
   * @brief Setter property for InputHeaderFile
   */
  void setInputHeaderFile(const QString& value);
  /**
   * @brief Getter property for InputHeaderFile
   * @return Value of InputHeaderFile
   */
  QString getInputHeaderFile() const;
  Q_PROPERTY(QString InputHeaderFile READ getInputHeaderFile WRITE setInputHeaderFile)

  /**
   * @brief Setter property for DataContainerName
   */
  void setDataContainerName(const QString& value);
  /**
   * @brief Getter property for DataContainerName
   * @return Value of DataContainerName
   */
  QString getDataContainerName() const;
  Q_PROPERTY(QString DataContainerName READ getDataContainerName WRITE setDataContainerName)

  /**
   * @brief Setter property for CellAttributeMatrixName
   */
  void setCellAttributeMatrixName(const QString& value);
  /**
   * @brief Getter property for CellAttributeMatrixName
   * @return Value of CellAttributeMatrixName
   */
  QString getCellAttributeMatrixName() const;
  Q_PROPERTY(QString CellAttributeMatrixName READ getCellAttributeMatrixName WRITE setCellAttributeMatrixName)

  /**
   * @brief Setter property for DensityArrayName
   */
  void setDensityArrayName(const QString& value);
  /**
   * @brief Getter property for DensityArrayName
   * @return Value of DensityArrayName
   */
  QString getDensityArrayName() const;
  Q_PROPERTY(QString DensityArrayName READ getDensityArrayName WRITE setDensityArrayName)

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
  ReadBinaryCTNorthStar();

  /**
   * @brief sanityCheckFileSizeVersusAllocatedSize Ensures the allocated array and the raw
   * binary file have the same number of bytes
   * @param allocatedBytes Number of bytes allocated
   * @param fileSize Size of the raw binary file
   * @return Integer error code
   */
  int32_t sanityCheckFileSizeVersusAllocatedSize(size_t allocatedBytes, size_t fileSize);

  /**
   * @brief readBinaryCTFile Reads the raw binary CT file
   * @return Integer error code
   */
  int32_t readBinaryCTFiles(size_t dims[3]);

  /**
   * @brief readHeaderMetaData Reads the number of voxels and voxel extents
   * from the NSI header file
   * @return Integer error code
   */
  int32_t readHeaderMetaData(ImageGeom::Pointer image);
  /**
   * @brief dataCheck Checks for the appropriate parameter values and availability of arrays
   */
  void dataCheck() override;

  /**
   * @brief Initializes all the private instance variables.
   */
  void initialize();

private:
  std::weak_ptr<DataArray<float>> m_DensityPtr;
  float* m_Density = nullptr;

  std::vector<QString> m_InputFiles = {};
  std::vector<int64_t> m_SlicesPerFile = {};
  QString m_InputHeaderFile = {};
  QString m_DataContainerName = {};
  QString m_CellAttributeMatrixName = {};
  QString m_DensityArrayName = {};

  QFile m_InHeaderStream;
  QFile m_InStream;

  ReadBinaryCTNorthStar(const ReadBinaryCTNorthStar&) = delete; // Copy Constructor Not Implemented
  ReadBinaryCTNorthStar(ReadBinaryCTNorthStar&&) = delete;      // Move Constructor Not Implemented
  void operator=(const ReadBinaryCTNorthStar&) = delete;        // Operator '=' Not Implemented
};

#endif /* _ReadBinaryCTNorthStar_H_ */
