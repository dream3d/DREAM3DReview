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
#ifndef _findsurfaceroughness_h_
#define _findsurfaceroughness_h_

#include <memory>

#include "SIMPLib/SIMPLib.h"
#include "SIMPLib/DataArrays/DataArray.hpp"
#include "SIMPLib/Filtering/AbstractFilter.h"

#include "DREAM3DReview/DREAM3DReviewDLLExport.h"

/**
 * @brief The FindSurfaceRoughness class. See [Filter documentation](@ref FindSurfaceRoughness) for details.
 */
class DREAM3DReview_EXPORT FindSurfaceRoughness : public AbstractFilter
{
  Q_OBJECT

public:
  using Self = FindSurfaceRoughness;
  using Pointer = std::shared_ptr<Self>;
  using ConstPointer = std::shared_ptr<const Self>;
  using WeakPointer = std::weak_ptr<Self>;
  using ConstWeakPointer = std::weak_ptr<const Self>;
  static Pointer NullPointer();

  static std::shared_ptr<FindSurfaceRoughness> New();

  /**
   * @brief Returns the name of the class for FindSurfaceRoughness
   */
  QString getNameOfClass() const override;
  /**
   * @brief Returns the name of the class for FindSurfaceRoughness
   */
  static QString ClassName();

  ~FindSurfaceRoughness() override;

  /**
   * @brief Setter property for BoundaryCellsArrayPath
   */
  void setBoundaryCellsArrayPath(const DataArrayPath& value);
  /**
   * @brief Getter property for BoundaryCellsArrayPath
   * @return Value of BoundaryCellsArrayPath
   */
  DataArrayPath getBoundaryCellsArrayPath() const;
  Q_PROPERTY(DataArrayPath BoundaryCellsArrayPath READ getBoundaryCellsArrayPath WRITE setBoundaryCellsArrayPath)

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
   * @brief Setter property for RoughnessParamsArrayName
   */
  void setRoughnessParamsArrayName(const QString& value);
  /**
   * @brief Getter property for RoughnessParamsArrayName
   * @return Value of RoughnessParamsArrayName
   */
  QString getRoughnessParamsArrayName() const;
  Q_PROPERTY(QString RoughnessParamsArrayName READ getRoughnessParamsArrayName WRITE setRoughnessParamsArrayName)

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
  FindSurfaceRoughness();

  /**
   * @brief dataCheck Checks for the appropriate parameter values and availability of arrays
   */
  void dataCheck() override;

  /**
   * @brief Initializes all the private instance variables.
   */
  void initialize();

private:
  std::weak_ptr<DataArray<int8_t>> m_BoundaryCellsPtr;
  int8_t* m_BoundaryCells = nullptr;
  std::weak_ptr<DataArray<double>> m_RoughnessParamsPtr;
  double* m_RoughnessParams = nullptr;

  DataArrayPath m_BoundaryCellsArrayPath = {SIMPL::Defaults::ImageDataContainerName, SIMPL::Defaults::CellAttributeMatrixName, SIMPL::CellData::BoundaryCells};
  QString m_AttributeMatrixName = {"RoughnessData"};
  QString m_RoughnessParamsArrayName = {"RougnessParameters"};

public:
  FindSurfaceRoughness(const FindSurfaceRoughness&) = delete;            // Copy Constructor Not Implemented
  FindSurfaceRoughness(FindSurfaceRoughness&&) = delete;                 // Move Constructor Not Implemented
  FindSurfaceRoughness& operator=(const FindSurfaceRoughness&) = delete; // Copy Assignment Not Implemented
  FindSurfaceRoughness& operator=(FindSurfaceRoughness&&) = delete;      // Move Assignment Not Implemented
};

#endif /* _FindSurfaceRoughness_H_ */
