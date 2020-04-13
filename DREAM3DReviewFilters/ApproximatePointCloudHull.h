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
#ifndef _approximatepointcloudhull_h_
#define _approximatepointcloudhull_h_

#include <memory>

#include "SIMPLib/FilterParameters/FloatVec3FilterParameter.h"
#include "SIMPLib/Filtering/AbstractFilter.h"
#include "SIMPLib/Geometry/ImageGeom.h"
#include "SIMPLib/SIMPLib.h"

/**
 * @brief The ApproximatePointCloudHull class. See [Filter documentation](@ref approximatepointcloudhull) for details.
 */
class ApproximatePointCloudHull : public AbstractFilter
{
  Q_OBJECT

public:
  using Self = ApproximatePointCloudHull;
  using Pointer = std::shared_ptr<Self>;
  using ConstPointer = std::shared_ptr<const Self>;
  using WeakPointer = std::weak_ptr<Self>;
  using ConstWeakPointer = std::weak_ptr<const Self>;
  static Pointer NullPointer();

  static std::shared_ptr<ApproximatePointCloudHull> New();

  /**
   * @brief Returns the name of the class for ApproximatePointCloudHull
   */
  QString getNameOfClass() const override;
  /**
   * @brief Returns the name of the class for ApproximatePointCloudHull
   */
  static QString ClassName();

  ~ApproximatePointCloudHull() override;

  /**
   * @brief Setter property for VertexDataContainerName
   */
  void setVertexDataContainerName(const DataArrayPath& value);
  /**
   * @brief Getter property for VertexDataContainerName
   * @return Value of VertexDataContainerName
   */
  DataArrayPath getVertexDataContainerName() const;
  Q_PROPERTY(DataArrayPath VertexDataContainerName READ getVertexDataContainerName WRITE setVertexDataContainerName)

  /**
   * @brief Setter property for HullDataContainerName
   */
  void setHullDataContainerName(const QString& value);
  /**
   * @brief Getter property for HullDataContainerName
   * @return Value of HullDataContainerName
   */
  QString getHullDataContainerName() const;
  Q_PROPERTY(QString HullDataContainerName READ getHullDataContainerName WRITE setHullDataContainerName)

  /**
   * @brief Setter property for GridResolution
   */
  void setGridResolution(const FloatVec3Type& value);
  /**
   * @brief Getter property for GridResolution
   * @return Value of GridResolution
   */
  FloatVec3Type getGridResolution() const;
  Q_PROPERTY(FloatVec3Type GridResolution READ getGridResolution WRITE setGridResolution)

  /**
   * @brief Setter property for NumberOfEmptyNeighbors
   */
  void setNumberOfEmptyNeighbors(int value);
  /**
   * @brief Getter property for NumberOfEmptyNeighbors
   * @return Value of NumberOfEmptyNeighbors
   */
  int getNumberOfEmptyNeighbors() const;
  Q_PROPERTY(int NumberOfEmptyNeighbors READ getNumberOfEmptyNeighbors WRITE setNumberOfEmptyNeighbors)

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
  ApproximatePointCloudHull();

  /**
  * @brief dataCheck Checks for the appropriate parameter values and availability of arrays
  */
  void dataCheck() override;

  /**
  * @brief Initializes all the private instance variables.
  */
  void initialize();

private:
  DataArrayPath m_VertexDataContainerName = {};
  QString m_HullDataContainerName = {};
  FloatVec3Type m_GridResolution = {};
  int m_NumberOfEmptyNeighbors = {};

  ImageGeom::Pointer m_SamplingGrid;
  std::vector<float> m_MeshMinExtents;
  std::vector<float> m_MeshMaxExtents;

  ApproximatePointCloudHull(const ApproximatePointCloudHull&) = delete; // Copy Constructor Not Implemented
  ApproximatePointCloudHull(ApproximatePointCloudHull&&) = delete;      // Move Constructor Not Implemented
  void operator=(const ApproximatePointCloudHull&) = delete;            // Operator '=' Not Implemented
};

#endif /* _ApproximatePointCloudHull_H_ */
