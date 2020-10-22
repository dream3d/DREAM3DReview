/*
 * Your License or Copyright can go here
 */

#pragma once

#include "SIMPLib/SIMPLib.h"
#include "SIMPLib/DataArrays/DataArray.hpp"
#include "SIMPLib/Filtering/AbstractFilter.h"

#include "DREAM3DReview/DREAM3DReviewPlugin.h"

/**
 * @brief The CreateArrayofIndices class. See [Filter documentation](@ref createarrayofindices) for details.
 */
class DREAM3DReview_EXPORT CreateArrayofIndices : public AbstractFilter
{
  Q_OBJECT

  // clang-format off
  PYB11_BEGIN_BINDINGS(CreateArrayofIndices SUPERCLASS AbstractFilter)
  PYB11_FILTER()
  PYB11_SHARED_POINTERS(CreateArrayofIndices)
  PYB11_FILTER_NEW_MACRO(CreateArrayofIndices)
  PYB11_PROPERTY(DataArrayPath IndexArrayPath READ getIndexArrayPath WRITE setIndexArrayPath)
  PYB11_END_BINDINGS()
  // clang-format on

public:
  using Self = CreateArrayofIndices;
  using Pointer = std::shared_ptr<Self>;
  using ConstPointer = std::shared_ptr<const Self>;
  using WeakPointer = std::weak_ptr<Self>;
  using ConstWeakPointer = std::weak_ptr<const Self>;
  static Pointer NullPointer();

  static Pointer New();

  /**
   * @brief Returns the name of the class for CreateArrayofIndices
   */
  QString getNameOfClass() const override;

  /**
   * @brief Returns the name of the class for CreateArrayofIndices
   */
  static QString ClassName();

  ~CreateArrayofIndices() override;

  /**
   * @brief Sets the value for Filter Parameter for IndexArrayPath
   * @param value The new value to use.
   */
  void setIndexArrayPath(const DataArrayPath& value);
  /**
   * @brief Gets the Filter Parameter value for IndexArrayPath
   * @return The value for IndexArrayPath
   */
  DataArrayPath getIndexArrayPath() const;
  Q_PROPERTY(DataArrayPath IndexArrayPath READ getIndexArrayPath WRITE setIndexArrayPath)

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
  CreateArrayofIndices();

  /**
   * @brief dataCheck Checks for the appropriate parameter values and availability of arrays
   */
  void dataCheck() override;

  /**
   * @brief Initializes all the private instance variables.
   */
  void initialize();

private:
  DataArrayPath m_IndexArrayPath;

  std::weak_ptr<DataArray<size_t>> m_IndicesPtr;
  size_t* m_Indices = nullptr;

public:
  CreateArrayofIndices(const CreateArrayofIndices&) = delete;            // Copy Constructor Not Implemented
  CreateArrayofIndices& operator=(const CreateArrayofIndices&) = delete; // Copy Assignment Not Implemented
  CreateArrayofIndices(CreateArrayofIndices&&) = delete;                 // Move Constructor Not Implemented
  CreateArrayofIndices& operator=(CreateArrayofIndices&&) = delete;      // Move Assignment Not Implemented
};
