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
#ifndef _laplaciasmoothpointcloud_h_
#define _laplaciasmoothpointcloud_h_

#include <memory>

#include "SIMPLib/Filtering/AbstractFilter.h"
#include "SIMPLib/SIMPLib.h"
#include "SIMPLib/DataArrays/DataArray.hpp"

/**
 * @brief The LaplacianSmoothPointCloud class. See [Filter documentation](@ref laplaciansmoothpointcloud) for details.
 */
class LaplacianSmoothPointCloud : public AbstractFilter
{
  Q_OBJECT

public:
  using Self = LaplacianSmoothPointCloud;
  using Pointer = std::shared_ptr<Self>;
  using ConstPointer = std::shared_ptr<const Self>;
  using WeakPointer = std::weak_ptr<Self>;
  using ConstWeakPointer = std::weak_ptr<const Self>;
  static Pointer NullPointer();

  static std::shared_ptr<LaplacianSmoothPointCloud> New();

  /**
   * @brief Returns the name of the class for LaplacianSmoothPointCloud
   */
  QString getNameOfClass() const override;
  /**
   * @brief Returns the name of the class for LaplacianSmoothPointCloud
   */
  static QString ClassName();

  ~LaplacianSmoothPointCloud() override;

  /**
   * @brief Setter property for DataContainerName
   */
  void setDataContainerName(const DataArrayPath& value);
  /**
   * @brief Getter property for DataContainerName
   * @return Value of DataContainerName
   */
  DataArrayPath getDataContainerName() const;

  Q_PROPERTY(DataArrayPath DataContainerName READ getDataContainerName WRITE setDataContainerName)

  /**
   * @brief Setter property for Lambda
   */
  void setLambda(float value);
  /**
   * @brief Getter property for Lambda
   * @return Value of Lambda
   */
  float getLambda() const;

  Q_PROPERTY(float Lambda READ getLambda WRITE setLambda)

  /**
   * @brief Setter property for NumIterations
   */
  void setNumIterations(int value);
  /**
   * @brief Getter property for NumIterations
   * @return Value of NumIterations
   */
  int getNumIterations() const;

  Q_PROPERTY(int NumIterations READ getNumIterations WRITE setNumIterations)

  /**
   * @brief Setter property for UseMask
   */
  void setUseMask(bool value);
  /**
   * @brief Getter property for UseMask
   * @return Value of UseMask
   */
  bool getUseMask() const;

  Q_PROPERTY(bool UseMask READ getUseMask WRITE setUseMask)

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
  LaplacianSmoothPointCloud();

  /**
   * @brief smoothPointCoordinates Smooths a point coordinate based on a window
   * @param vertex Vertex pointer
   * @param startCoords Start window coordinates
   * @param endCoords End window coordinates
   * @param vertIdx Vertex index
   */
  void determineNewCoords(float* vertex, int64_t vertId, float* newCoords);
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

  DataArrayPath m_DataContainerName = {};
  float m_Lambda = 0.1f;
  int m_NumIterations = 1;
  bool m_UseMask = false;
  DataArrayPath m_MaskArrayPath = {};

public:
  LaplacianSmoothPointCloud(const LaplacianSmoothPointCloud&) = delete;            // Copy Constructor Not Implemented
  LaplacianSmoothPointCloud(LaplacianSmoothPointCloud&&) = delete;                 // Move Constructor Not Implemented
  LaplacianSmoothPointCloud& operator=(const LaplacianSmoothPointCloud&) = delete; // Copy Assignment Not Implemented
  LaplacianSmoothPointCloud& operator=(LaplacianSmoothPointCloud&&) = delete;      // Move Assignment Not Implemented
};

#endif /* _laplaciansmoothpointcloud_H_ */
