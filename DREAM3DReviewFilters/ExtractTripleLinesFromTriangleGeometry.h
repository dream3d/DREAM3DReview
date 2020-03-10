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
#ifndef _extracttriplelinesfromtrianglegeometry_h_
#define _extracttriplelinesfromtrianglegeometry_h_

#include <memory>

#include "SIMPLib/Filtering/AbstractFilter.h"
#include "SIMPLib/SIMPLib.h"
#include "SIMPLib/DataArrays/DataArray.hpp"

/**
 * @brief The ExtractTripleLinesFromTriangleGeometry class. See [Filter documentation](@ref extracttriplelinesfromtrianglegeometry) for details.
 */
class ExtractTripleLinesFromTriangleGeometry : public AbstractFilter
{
  Q_OBJECT

public:
  using Self = ExtractTripleLinesFromTriangleGeometry;
  using Pointer = std::shared_ptr<Self>;
  using ConstPointer = std::shared_ptr<const Self>;
  using WeakPointer = std::weak_ptr<Self>;
  using ConstWeakPointer = std::weak_ptr<const Self>;
  static Pointer NullPointer();

  static std::shared_ptr<ExtractTripleLinesFromTriangleGeometry> New();

  /**
   * @brief Returns the name of the class for ExtractTripleLinesFromTriangleGeometry
   */
  QString getNameOfClass() const override;
  /**
   * @brief Returns the name of the class for ExtractTripleLinesFromTriangleGeometry
   */
  static QString ClassName();

  ~ExtractTripleLinesFromTriangleGeometry() override;

  /**
   * @brief Setter property for NodeTypesArrayPath
   */
  void setNodeTypesArrayPath(const DataArrayPath& value);
  /**
   * @brief Getter property for NodeTypesArrayPath
   * @return Value of NodeTypesArrayPath
   */
  DataArrayPath getNodeTypesArrayPath() const;

  Q_PROPERTY(DataArrayPath NodeTypesArrayPath READ getNodeTypesArrayPath WRITE setNodeTypesArrayPath)

  /**
   * @brief Setter property for EdgeGeometry
   */
  void setEdgeGeometry(const QString& value);
  /**
   * @brief Getter property for EdgeGeometry
   * @return Value of EdgeGeometry
   */
  QString getEdgeGeometry() const;

  Q_PROPERTY(QString EdgeGeometry READ getEdgeGeometry WRITE setEdgeGeometry)

  /**
   * @brief Setter property for VertexAttributeMatrixName
   */
  void setVertexAttributeMatrixName(const QString& value);
  /**
   * @brief Getter property for VertexAttributeMatrixName
   * @return Value of VertexAttributeMatrixName
   */
  QString getVertexAttributeMatrixName() const;

  Q_PROPERTY(QString VertexAttributeMatrixName READ getVertexAttributeMatrixName WRITE setVertexAttributeMatrixName)

  /**
   * @brief Setter property for EdgeAttributeMatrixName
   */
  void setEdgeAttributeMatrixName(const QString& value);
  /**
   * @brief Getter property for EdgeAttributeMatrixName
   * @return Value of EdgeAttributeMatrixName
   */
  QString getEdgeAttributeMatrixName() const;

  Q_PROPERTY(QString EdgeAttributeMatrixName READ getEdgeAttributeMatrixName WRITE setEdgeAttributeMatrixName)

  /**
   * @brief Setter property for NodeTypesArrayName
   */
  void setNodeTypesArrayName(const QString& value);
  /**
   * @brief Getter property for NodeTypesArrayName
   * @return Value of NodeTypesArrayName
   */
  QString getNodeTypesArrayName() const;

  Q_PROPERTY(QString NodeTypesArrayName READ getNodeTypesArrayName WRITE setNodeTypesArrayName)

  /**
   * @brief Setter property for SmoothTripleLines
   */
  void setSmoothTripleLines(bool value);
  /**
   * @brief Getter property for SmoothTripleLines
   * @return Value of SmoothTripleLines
   */
  bool getSmoothTripleLines() const;

  Q_PROPERTY(bool SmoothTripleLines READ getSmoothTripleLines WRITE setSmoothTripleLines)

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
  ExtractTripleLinesFromTriangleGeometry();

  void extractTripleLines();

  void smoothTripleLines();

  std::vector<float> computeBsplineInterpolation(std::vector<std::vector<float>> points, float minVal, float maxVal);

  /**
  * @brief dataCheck Checks for the appropriate parameter values and availability of arrays
  */
  void dataCheck() override;

  /**
  * @brief Initializes all the private instance variables.
  */
  void initialize();

private:
  std::weak_ptr<DataArray<int8_t>> m_NodeTypesPtr;
  int8_t* m_NodeTypes = nullptr;
  std::weak_ptr<DataArray<int8_t>> m_TripleLineNodeTypesPtr;
  int8_t* m_TripleLineNodeTypes = nullptr;

  DataArrayPath m_NodeTypesArrayPath = {};
  QString m_EdgeGeometry = {};
  QString m_VertexAttributeMatrixName = {};
  QString m_EdgeAttributeMatrixName = {};
  QString m_NodeTypesArrayName = {};
  bool m_SmoothTripleLines = false;

public:
  ExtractTripleLinesFromTriangleGeometry(const ExtractTripleLinesFromTriangleGeometry&) = delete;            // Copy Constructor Not Implemented
  ExtractTripleLinesFromTriangleGeometry(ExtractTripleLinesFromTriangleGeometry&&) = delete;                 // Move Constructor Not Implemented
  ExtractTripleLinesFromTriangleGeometry& operator=(const ExtractTripleLinesFromTriangleGeometry&) = delete; // Copy Assignment Not Implemented
  ExtractTripleLinesFromTriangleGeometry& operator=(ExtractTripleLinesFromTriangleGeometry&&) = delete;      // Move Assignment Not Implemented
};

#endif /* _ExtractTripleLinesFromTriangleGeometry_H_ */
