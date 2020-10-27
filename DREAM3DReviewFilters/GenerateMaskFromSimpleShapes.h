/*
 * Your License or Copyright can go here
 */

#pragma once

#include "SIMPLib/SIMPLib.h"
#include "SIMPLib/DataArrays/DataArray.hpp"
#include "SIMPLib/DataContainers/DataContainer.h"
#include "SIMPLib/Filtering/AbstractFilter.h"

#include "DREAM3DReview/DREAM3DReviewPlugin.h"

/**
 * @brief The GenerateMaskFromSimpleShapes class. See [Filter documentation](@ref generatemaskfromsimpleshapes) for details.
 */
class DREAM3DReview_EXPORT GenerateMaskFromSimpleShapes : public AbstractFilter
{
  Q_OBJECT

  // Start Python bindings declarations
  PYB11_BEGIN_BINDINGS(GenerateMaskFromSimpleShapes SUPERCLASS AbstractFilter)
  PYB11_FILTER()
  PYB11_SHARED_POINTERS(GenerateMaskFromSimpleShapes)
  PYB11_FILTER_NEW_MACRO(GenerateMaskFromSimpleShapes)
  PYB11_PROPERTY(DataArrayPath MaskArrayPath READ getMaskArrayPath WRITE setMaskArrayPath)
  PYB11_PROPERTY(DataArrayPath CentersArrayPath READ getCentersArrayPath WRITE setCentersArrayPath)
  PYB11_PROPERTY(DataArrayPath AxesLengthArrayPath READ getAxesLengthArrayPath WRITE setAxesLengthArrayPath)
  PYB11_PROPERTY(DataArrayPath BoxDimensionsArrayPath READ getBoxDimensionsArrayPath WRITE setBoxDimensionsArrayPath)
  PYB11_PROPERTY(DataArrayPath CylinderRadiusArrayPath READ getCylinderRadiusArrayPath WRITE setCylinderRadiusArrayPath)
  PYB11_PROPERTY(DataArrayPath CylinderHeightArrayPath READ getCylinderHeightArrayPath WRITE setCylinderHeightArrayPath)
  PYB11_PROPERTY(int MaskShape READ getMaskShape WRITE setMaskShape)
  PYB11_END_BINDINGS()
  // End Python bindings declarations

public:
  using Self = GenerateMaskFromSimpleShapes;
  using Pointer = std::shared_ptr<Self>;
  using ConstPointer = std::shared_ptr<const Self>;
  using WeakPointer = std::weak_ptr<Self>;
  using ConstWeakPointer = std::weak_ptr<const Self>;
  static Pointer NullPointer();

  static Pointer New();

  /**
   * @brief Returns the name of the class for GenerateMaskFromSimpleShapes
   */
  QString getNameOfClass() const override;

  /**
   * @brief Returns the name of the class for GenerateMaskFromSimpleShapes
   */
  static QString ClassName();

  ~GenerateMaskFromSimpleShapes() override;

  /**
   * @brief Sets the value for Filter Parameter for MaskArrayPath
   * @param value The new value to use.
   */
  void setMaskArrayPath(const DataArrayPath& value);
  /**
   * @brief Gets the Filter Parameter value for MaskArrayPath
   * @return The value for MaskArrayPath
   */
  DataArrayPath getMaskArrayPath() const;
  Q_PROPERTY(DataArrayPath MaskArrayPath READ getMaskArrayPath WRITE setMaskArrayPath)

  /**
   * @brief Sets the value for Filter Parameter for CentersArrayPath
   * @param value The new value to use.
   */
  void setCentersArrayPath(const DataArrayPath& value);
  /**
   * @brief Gets the Filter Parameter value for CentersArrayPath
   * @return The value for CentersArrayPath
   */
  DataArrayPath getCentersArrayPath() const;
  Q_PROPERTY(DataArrayPath CentersArrayPath READ getCentersArrayPath WRITE setCentersArrayPath)

  /**
   * @brief Sets the value for Filter Parameter for AxesLengthArrayPath
   * @param value The new value to use.
   */
  void setAxesLengthArrayPath(const DataArrayPath& value);
  /**
   * @brief Gets the Filter Parameter value for AxesLengthArrayPath
   * @return The value for AxesLengthArrayPath
   */
  DataArrayPath getAxesLengthArrayPath() const;
  Q_PROPERTY(DataArrayPath AxesLengthArrayPath READ getAxesLengthArrayPath WRITE setAxesLengthArrayPath)

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
   * @brief Sets the value for Filter Parameter for CylinderRadiusArrayPath
   * @param value The new value to use.
   */
  void setCylinderRadiusArrayPath(const DataArrayPath& value);
  /**
   * @brief Gets the Filter Parameter value for CylinderRadiusArrayPath
   * @return The value for CylinderRadiusArrayPath
   */
  DataArrayPath getCylinderRadiusArrayPath() const;
  Q_PROPERTY(DataArrayPath CylinderRadiusArrayPath READ getCylinderRadiusArrayPath WRITE setCylinderRadiusArrayPath)

  /**
   * @brief Sets the value for Filter Parameter for CylinderHeightArrayPath
   * @param value The new value to use.
   */
  void setCylinderHeightArrayPath(const DataArrayPath& value);
  /**
   * @brief Gets the Filter Parameter value for CylinderHeightArrayPath
   * @return The value for CylinderHeightArrayPath
   */
  DataArrayPath getCylinderHeightArrayPath() const;
  Q_PROPERTY(DataArrayPath CylinderHeightArrayPath READ getCylinderHeightArrayPath WRITE setCylinderHeightArrayPath)

  /**
   * @brief Setter property for ReferenceOrientation
   */
  void setMaskShape(int value);
  /**
   * @brief Getter property for ReferenceOrientation
   * @return Value of ReferenceOrientation
   */
  int getMaskShape() const;
  Q_PROPERTY(int MaskShape READ getMaskShape WRITE setMaskShape)

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
  GenerateMaskFromSimpleShapes();

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
  void createImageMask();

  /**
   * @brief Initializes all the private instance variables.
   */
  void createVertexMask();

private:
  std::weak_ptr<DataArray<bool>> m_MaskPtr;
  bool* m_Mask = nullptr;

  std::weak_ptr<DataArray<float>> m_CentersPtr;
  float* m_Centers = nullptr;

  std::weak_ptr<DataArray<float>> m_AxisLengthsPtr;
  float* m_AxisLengths = nullptr;

  std::weak_ptr<DataArray<float>> m_BoxDimsPtr;
  float* m_BoxDims = nullptr;

  std::weak_ptr<DataArray<float>> m_CylinderRadPtr;
  float* m_CylinderRad = nullptr;

  std::weak_ptr<DataArray<float>> m_CylinderHeightPtr;
  float* m_CylinderHeight = nullptr;

  DataArrayPath m_MaskArrayPath;
  DataArrayPath m_CentersArrayPath;
  DataArrayPath m_AxesLengthArrayPath;
  DataArrayPath m_BoxDimensionsArrayPath;
  DataArrayPath m_CylinderRadiusArrayPath;
  DataArrayPath m_CylinderHeightArrayPath;

  AttributeMatrix::Type m_DestAttributeMatrixType;

  int m_MaskShape = {};

public:
  GenerateMaskFromSimpleShapes(const GenerateMaskFromSimpleShapes&) = delete;            // Copy Constructor Not Implemented
  GenerateMaskFromSimpleShapes& operator=(const GenerateMaskFromSimpleShapes&) = delete; // Copy Assignment Not Implemented
  GenerateMaskFromSimpleShapes(GenerateMaskFromSimpleShapes&&) = delete;                 // Move Constructor Not Implemented
  GenerateMaskFromSimpleShapes& operator=(GenerateMaskFromSimpleShapes&&) = delete;      // Move Assignment Not Implemented
};
