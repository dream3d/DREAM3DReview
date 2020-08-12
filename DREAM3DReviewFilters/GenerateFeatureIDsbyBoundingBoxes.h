/*
 * Your License or Copyright can go here
 */

#pragma once

#include "SIMPLib/SIMPLib.h"
#include "SIMPLib/Filtering/AbstractFilter.h"
#include "SIMPLib/DataArrays/DataArray.hpp"
#include "SIMPLib/DataContainers/DataContainer.h"

#include "DREAM3DReview/DREAM3DReviewPlugin.h"

/**
 * @brief The GenerateFeatureIDsbyBoundingBoxes class. See [Filter documentation](@ref generatefeatureidsbyboundingboxes) for details.
 */
class DREAM3DReview_EXPORT GenerateFeatureIDsbyBoundingBoxes : public AbstractFilter
{
  Q_OBJECT

  // Start Python bindings declarations
  PYB11_BEGIN_BINDINGS(GenerateFeatureIDsbyBoundingBoxes SUPERCLASS AbstractFilter)
  PYB11_FILTER()
  PYB11_SHARED_POINTERS(GenerateFeatureIDsbyBoundingBoxes)
  PYB11_FILTER_NEW_MACRO(GenerateFeatureIDsbyBoundingBoxes)
  PYB11_PROPERTY(DataArrayPath FeatureIDsArrayPath READ getFeatureIDsArrayPath WRITE setFeatureIDsArrayPath)
  PYB11_PROPERTY(DataArrayPath FeatureAttributeMatrixArrayPath READ getFeatureAttributeMatrixArrayPath WRITE setFeatureAttributeMatrixArrayPath)
  PYB11_PROPERTY(DataArrayPath BoxCenterArrayPath READ getBoxCenterArrayPath WRITE setBoxCenterArrayPath)
  PYB11_PROPERTY(DataArrayPath BoxDimensionsArrayPath READ getBoxDimensionsArrayPath WRITE setBoxDimensionsArrayPath)
  PYB11_PROPERTY(DataArrayPath BoxFeatureIDsArrayPath READ getBoxFeatureIDsArrayPath WRITE setBoxFeatureIDsArrayPath)
  PYB11_END_BINDINGS()
  // End Python bindings declarations

public:
  using Self = GenerateFeatureIDsbyBoundingBoxes;
  using Pointer = std::shared_ptr<Self>;
  using ConstPointer = std::shared_ptr<const Self>;
  using WeakPointer = std::weak_ptr<Self>;
  using ConstWeakPointer = std::weak_ptr<const Self>;
  static Pointer NullPointer();

  static std::shared_ptr<GenerateFeatureIDsbyBoundingBoxes> New();

  /**
   * @brief Returns the name of the class for GenerateFeatureIDsbyBoundingBoxes
   */
  QString getNameOfClass() const override;

  /**
   * @brief Returns the name of the class for GenerateFeatureIDsbyBoundingBoxes
   */
  static QString ClassName();

  ~GenerateFeatureIDsbyBoundingBoxes() override;

  /**
   * @brief Sets the value for Filter Parameter for FeatureIDsArrayPath
   * @param value The new value to use.
   */
  void setFeatureIDsArrayPath(const DataArrayPath& value);
  /**
   * @brief Gets the Filter Parameter value for FeatureIDsArrayPath
   * @return The value for FeatureIDsArrayPath
   */
  DataArrayPath getFeatureIDsArrayPath() const;
  Q_PROPERTY(DataArrayPath FeatureIDsArrayPath READ getFeatureIDsArrayPath WRITE setFeatureIDsArrayPath)

  /**
   * @brief Sets the value for Filter Parameter for FeatureAttributeMatrixArrayPath
   * @param value The new value to use.
   */
  void setFeatureAttributeMatrixArrayPath(const DataArrayPath& value);
  /**
   * @brief Gets the Filter Parameter value for FeatureAttributeMatrixArrayPath
   * @return The value for FeatureAttributeMatrixArrayPath
   */
  DataArrayPath getFeatureAttributeMatrixArrayPath() const;
  Q_PROPERTY(DataArrayPath FeatureAttributeMatrixArrayPath READ getFeatureAttributeMatrixArrayPath WRITE setFeatureAttributeMatrixArrayPath)

  /**
   * @brief Sets the value for Filter Parameter for BoxCenterArrayPath
   * @param value The new value to use.
   */
  void setBoxCenterArrayPath(const DataArrayPath& value);
  /**
   * @brief Gets the Filter Parameter value for BoxCenterArrayPath
   * @return The value for BoxCenterArrayPath
   */
  DataArrayPath getBoxCenterArrayPath() const;
  Q_PROPERTY(DataArrayPath BoxCenterArrayPath READ getBoxCenterArrayPath WRITE setBoxCenterArrayPath)

  /**
   * @brief Sets the value for Filter Parameter for BoxDimensionsArrayPath
   * @param value The new value to use.
   */
  void setBoxDimensionsArrayPath(const DataArrayPath& value);
  /**
   * @brief Gets the Filter Parameter value for BoxDimensionsArrayPath
   * @return The value for BoxDimensionsArrayPath
   */
  DataArrayPath getBoxDimensionsArrayPath() const;
  Q_PROPERTY(DataArrayPath BoxDimensionsArrayPath READ getBoxDimensionsArrayPath WRITE setBoxDimensionsArrayPath)

  /**
   * @brief Sets the value for Filter Parameter for BoxFeatureIDsArrayPath
   * @param value The new value to use.
   */
  void setBoxFeatureIDsArrayPath(const DataArrayPath& value);
  /**
   * @brief Gets the Filter Parameter value for BoxFeatureIDsArrayPath
   * @return The value for BoxFeatureIDsArrayPath
   */
  DataArrayPath getBoxFeatureIDsArrayPath() const;
  Q_PROPERTY(DataArrayPath BoxFeatureIDsArrayPath READ getBoxFeatureIDsArrayPath WRITE setBoxFeatureIDsArrayPath)

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
  GenerateFeatureIDsbyBoundingBoxes();

  /**
   * @brief dataCheck Checks for the appropriate parameter values and availability of arrays
   */
  void dataCheck() override;

  /**
   * @brief Initializes all the private instance variables.
   */
  void initialize();

  /**
   * @brief Initializes all the private instance variables.
   */
  void checkBoundingBoxImage();

  /**
   * @brief Initializes all the private instance variables.
   */
  void checkBoundingBoxEdge();

  /**
   * @brief Initializes all the private instance variables.
   */
  void checkBoundingBoxVertex();

private:
  std::weak_ptr<DataArray<int32_t>> m_FeatureIdsPtr;
  int32_t* m_FeatureIds = nullptr;
  std::weak_ptr<DataArray<int32_t>> m_BoxFeatureIdsPtr;
  int32_t* m_BoxFeatureIds = nullptr;
  std::weak_ptr<DataArray<float>> m_BoxCenterPtr;
  float* m_BoxCenter = nullptr;
  std::weak_ptr<DataArray<float>> m_BoxDimsPtr;
  float* m_BoxDims = nullptr;
  AttributeMatrix::Type m_DestAttributeMatrixType;

  DataArrayPath m_FeatureIDsArrayPath;
  DataArrayPath m_FeatureAttributeMatrixArrayPath;
  DataArrayPath m_BoxCenterArrayPath;
  DataArrayPath m_BoxDimensionsArrayPath;
  DataArrayPath m_BoxFeatureIDsArrayPath;

public:
  GenerateFeatureIDsbyBoundingBoxes(const GenerateFeatureIDsbyBoundingBoxes&) = delete;            // Copy Constructor Not Implemented
  GenerateFeatureIDsbyBoundingBoxes& operator=(const GenerateFeatureIDsbyBoundingBoxes&) = delete; // Copy Assignment Not Implemented
  GenerateFeatureIDsbyBoundingBoxes(GenerateFeatureIDsbyBoundingBoxes&&) = delete;                 // Move Constructor Not Implemented
  GenerateFeatureIDsbyBoundingBoxes& operator=(GenerateFeatureIDsbyBoundingBoxes&&) = delete;      // Move Assignment Not Implemented
};

