#pragma once

#include "SIMPLib/SIMPLib.h"
#include "SIMPLib/Filtering/AbstractFilter.h"

#include "DREAM3DReview/DREAM3DReviewDLLExport.h"

/**
 * @brief The DecimatePointCloud class. See [Filter documentation](@ref decimatepointcloud) for details.
 */
class DREAM3DReview_EXPORT DecimatePointCloud : public AbstractFilter
{
  Q_OBJECT
  // clang-format off
  PYB11_BEGIN_BINDINGS(DecimatePointCloud SUPERCLASS AbstractFilter)
  PYB11_FILTER()
  PYB11_SHARED_POINTERS(DecimatePointCloud)
  PYB11_FILTER_NEW_MACRO(DecimatePointCloud)
  PYB11_END_BINDINGS()
  // clang-format on
public:
  using Self = DecimatePointCloud;
  using Pointer = std::shared_ptr<Self>;
  using ConstPointer = std::shared_ptr<const Self>;
  using WeakPointer = std::weak_ptr<Self>;
  using ConstWeakPointer = std::weak_ptr<const Self>;
  static Pointer NullPointer();

  static std::shared_ptr<DecimatePointCloud> New();

  /**
   * @brief Returns the name of the class for DecimatePointCloud
   */
  QString getNameOfClass() const override;

  /**
   * @brief Returns the name of the class for DecimatePointCloud
   */
  static QString ClassName();

  ~DecimatePointCloud() override;

  /**
   * @brief Setter property for VertexAttrMatPath
   */
  void setVertexAttrMatPath(const DataArrayPath& value);

  /**
   * @brief Getter property for VertexAttrMatPath
   * @return Value of VertexAttrMatPath
   */
  DataArrayPath getVertexAttrMatPath() const;
  Q_PROPERTY(DataArrayPath VertexAttrMatPath READ getVertexAttrMatPath WRITE setVertexAttrMatPath)

  /**
   * @brief Setter property for DecimationFreq
   */
  void setDecimationFreq(const float& value);

  /**
   * @brief Getter property for DecimationFreq
   * @return Value of DecimationFreq
   */
  float getDecimationFreq() const;
  Q_PROPERTY(float DecimationFreq READ getDecimationFreq WRITE setDecimationFreq)

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
  DecimatePointCloud();

  /**
   * @brief dataCheck Checks for the appropriate parameter values and availability of arrays
   */
  void dataCheck() override;

  /**
   * @brief Initializes all the private instance variables.
   */
  void initialize();

private:
  DataArrayPath m_VertexAttrMatPath = {"", "", ""};
  float m_DecimationFreq = {0.5};

  DecimatePointCloud(const DecimatePointCloud&) = delete; // Copy Constructor Not Implemented
  DecimatePointCloud(DecimatePointCloud&&) = delete;      // Move Constructor Not Implemented
  void operator=(const DecimatePointCloud&) = delete;     // Operator '=' Not Implemented
};
