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

#include "SIMPLib/Filtering/AbstractFilter.h"
#include "SIMPLib/SIMPLib.h"
#include "SIMPLib/DataArrays/DataArray.hpp"

/**
 * @brief The FindMinkowskiBouligandDimension class. See [Filter documentation](@ref findminkowskibouliganddimension) for details.
 */
class FindMinkowskiBouligandDimension : public AbstractFilter
{
  Q_OBJECT

public:
  using Self = FindMinkowskiBouligandDimension;
  using Pointer = std::shared_ptr<Self>;
  using ConstPointer = std::shared_ptr<const Self>;
  using WeakPointer = std::weak_ptr<Self>;
  using ConstWeakPointer = std::weak_ptr<const Self>;
  static Pointer NullPointer();

  static std::shared_ptr<FindMinkowskiBouligandDimension> New();

  /**
   * @brief Returns the name of the class for FindMinkowskiBouligandDimension
   */
  QString getNameOfClass() const override;
  /**
   * @brief Returns the name of the class for FindMinkowskiBouligandDimension
   */
  static QString ClassName();

  ~FindMinkowskiBouligandDimension() override;

  /**
   * @brief Setter property for MaskArrayPath
   */
  void setMaskArrayPath(const DataArrayPath& value);
  /**
   * @brief Getter property for MaskArrayPath
   * @return Value of MaskArrayPath
   */
  DataArrayPath getMaskArrayPath() const;

  Q_PROPERTY(DataArrayPath MaskArrayPath READ getMaskArrayPath WRITE setMaskArrayPath)

  /**
   * @brief Setter property for AttributeMatrixName
   */
  void setAttributeMatrixName(const QString& value);
  /**
   * @brief Getter property for AttributeMatrixName
   * @return Value of AttributeMatrixName
   */
  QString getAttributeMatrixName() const;

  Q_PROPERTY(QString AttributeMatrixName READ getAttributeMatrixName WRITE setAttributeMatrixName)

  /**
   * @brief Setter property for MinkowskiBouligandDimensionArrayName
   */
  void setMinkowskiBouligandDimensionArrayName(const QString& value);
  /**
   * @brief Getter property for MinkowskiBouligandDimensionArrayName
   * @return Value of MinkowskiBouligandDimensionArrayName
   */
  QString getMinkowskiBouligandDimensionArrayName() const;

  Q_PROPERTY(QString MinkowskiBouligandDimensionArrayName READ getMinkowskiBouligandDimensionArrayName WRITE setMinkowskiBouligandDimensionArrayName)

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
  FindMinkowskiBouligandDimension();

  /**
  * @brief dataCheck Checks for the appropriate parameter values and availability of arrays
  */
  void dataCheck() override;

  /**
  * @brief Initializes all the private instance variables.
  */
  void initialize();

private:
  std::weak_ptr<DataArray<bool>> m_MaskPtr;
  bool* m_Mask = nullptr;
  std::weak_ptr<DataArray<double>> m_MinkowskiBouligandDimensionPtr;
  double* m_MinkowskiBouligandDimension = nullptr;

  DataArrayPath m_MaskArrayPath = {};
  QString m_AttributeMatrixName = {};
  QString m_MinkowskiBouligandDimensionArrayName = {};

  void findSourceIndices(std::vector<size_t>& indices, size_t sourceDim, size_t index, size_t box[3], size_t imageDim);

public:
  FindMinkowskiBouligandDimension(const FindMinkowskiBouligandDimension&) = delete;            // Copy Constructor Not Implemented
  FindMinkowskiBouligandDimension(FindMinkowskiBouligandDimension&&) = delete;                 // Move Constructor Not Implemented
  FindMinkowskiBouligandDimension& operator=(const FindMinkowskiBouligandDimension&) = delete; // Copy Assignment Not Implemented
  FindMinkowskiBouligandDimension& operator=(FindMinkowskiBouligandDimension&&) = delete;      // Move Assignment Not Implemented
};

