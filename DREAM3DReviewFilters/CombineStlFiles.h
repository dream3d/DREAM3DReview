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
#pragma once

#include <memory>

#include <QtCore/QFileInfoList>
#include <QtCore/QString>

#include "SIMPLib/SIMPLib.h"
#include "SIMPLib/Filtering/AbstractFilter.h"
#include "SIMPLib/DataArrays/DataArray.hpp"
#include "SIMPLib/DataArrays/DataArray.hpp"

/**
 * @brief The CombineStlFiles class. See [Filter documentation](@ref combinestlfiles) for details.
 */
class CombineStlFiles : public AbstractFilter
{
  Q_OBJECT

public:
  using Self = CombineStlFiles;
  using Pointer = std::shared_ptr<Self>;
  using ConstPointer = std::shared_ptr<const Self>;
  using WeakPointer = std::weak_ptr<Self>;
  using ConstWeakPointer = std::weak_ptr<const Self>;
  static Pointer NullPointer();

  static std::shared_ptr<CombineStlFiles> New();

  /**
   * @brief Returns the name of the class for CombineStlFiles
   */
  QString getNameOfClass() const override;
  /**
   * @brief Returns the name of the class for CombineStlFiles
   */
  static QString ClassName();

  ~CombineStlFiles() override;

  /**
   * @brief Setter property for StlFilesPath
   */
  void setStlFilesPath(const QString& value);
  /**
   * @brief Getter property for StlFilesPath
   * @return Value of StlFilesPath
   */
  QString getStlFilesPath() const;

  Q_PROPERTY(QString StlFilesPath READ getStlFilesPath WRITE setStlFilesPath)

  /**
   * @brief Setter property for TriangleDataContainerName
   */
  void setTriangleDataContainerName(const QString& value);
  /**
   * @brief Getter property for TriangleDataContainerName
   * @return Value of TriangleDataContainerName
   */
  QString getTriangleDataContainerName() const;

  Q_PROPERTY(QString TriangleDataContainerName READ getTriangleDataContainerName WRITE setTriangleDataContainerName)

  /**
   * @brief Setter property for FaceAttributeMatrixName
   */
  void setFaceAttributeMatrixName(const QString& value);
  /**
   * @brief Getter property for FaceAttributeMatrixName
   * @return Value of FaceAttributeMatrixName
   */
  QString getFaceAttributeMatrixName() const;

  Q_PROPERTY(QString FaceAttributeMatrixName READ getFaceAttributeMatrixName WRITE setFaceAttributeMatrixName)

  /**
   * @brief Setter property for FaceNormalsArrayName
   */
  void setFaceNormalsArrayName(const QString& value);
  /**
   * @brief Getter property for FaceNormalsArrayName
   * @return Value of FaceNormalsArrayName
   */
  QString getFaceNormalsArrayName() const;

  Q_PROPERTY(QString FaceNormalsArrayName READ getFaceNormalsArrayName WRITE setFaceNormalsArrayName)

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
  CombineStlFiles();

  /**
  * @brief dataCheck Checks for the appropriate parameter values and availability of arrays
  */
  void dataCheck() override;

  /**
  * @brief Initializes all the private instance variables.
  */
  void initialize();

private:
  std::weak_ptr<DoubleArrayType> m_FaceNormalsPtr;
  double* m_FaceNormals = nullptr;

  QString m_StlFilesPath = {};
  QString m_TriangleDataContainerName = {};
  QString m_FaceAttributeMatrixName = {};
  QString m_FaceNormalsArrayName = {};

  QFileInfoList m_FileList;

public:
  CombineStlFiles(const CombineStlFiles&) = delete;            // Copy Constructor Not Implemented
  CombineStlFiles(CombineStlFiles&&) = delete;                 // Move Constructor Not Implemented
  CombineStlFiles& operator=(const CombineStlFiles&) = delete; // Copy Assignment Not Implemented
  CombineStlFiles& operator=(CombineStlFiles&&) = delete;      // Move Assignment Not Implemented
};
