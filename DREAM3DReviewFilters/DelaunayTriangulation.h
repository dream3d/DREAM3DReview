#pragma once

#include "SIMPLib/SIMPLib.h"
#include "SIMPLib/DataArrays/DataArray.hpp"
#include "SIMPLib/Filtering/AbstractFilter.h"

#include "DREAM3DReview/DREAM3DReviewPlugin.h"

/**
 * @brief The DelaunayTriangulation class. See [Filter documentation](@ref delaunaytriangulation) for details.
 */
class DREAM3DReview_EXPORT DelaunayTriangulation : public AbstractFilter
{
  Q_OBJECT

public:
  using Self = DelaunayTriangulation;
  using Pointer = std::shared_ptr<Self>;
  using ConstPointer = std::shared_ptr<const Self>;
  using WeakPointer = std::weak_ptr<Self>;
  using ConstWeakPointer = std::weak_ptr<const Self>;
  static Pointer NullPointer();

  static std::shared_ptr<DelaunayTriangulation> New();

  /**
   * @brief Returns the name of the class for DelaunayTriangulation
   */
  QString getNameOfClass() const override;

  /**
   * @brief Returns the name of the class for DelaunayTriangulation
   */
  static QString ClassName();

  ~DelaunayTriangulation() override;

  /**
   * @brief Setter property for InputGeometry
   */
  void setInputGeometry(const DataArrayPath& value);

  /**
   * @brief Getter property for InputGeometry
   * @return Value of InputGeometry
   */
  DataArrayPath getInputGeometry() const;
  Q_PROPERTY(DataArrayPath InputGeometry READ getInputGeometry WRITE setInputGeometry)

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
   * @brief Setter property for Offset
   */
  void setOffset(const double& value);

  /**
   * @brief Getter property for Offset
   * @return Value of Offset
   */
  double getOffset() const;
  Q_PROPERTY(double Offset READ getOffset WRITE setOffset)

  /**
   * @brief Setter property for Tolerance
   */
  void setTolerance(const double& value);

  /**
   * @brief Getter property for Tolerance
   * @return Value of Tolerance
   */
  double getTolerance() const;
  Q_PROPERTY(double Tolerance READ getTolerance WRITE setTolerance)

  /**
   * @brief Setter property for TriangulateByFeature
   */
  void setTriangulateByFeature(const bool& value);

  /**
   * @brief Getter property for TriangulateByFeature
   * @return Value of TriangulateByFeature
   */
  bool getTriangulateByFeature() const;
  Q_PROPERTY(bool TriangulateByFeature READ getTriangulateByFeature WRITE setTriangulateByFeature)

  /**
   * @brief Setter property for FeatureIdsArrayPath
   */
  void setFeatureIdsArrayPath(const DataArrayPath& value);

  /**
   * @brief Getter property for FeatureIdsArrayPath
   * @return Value of FeatureIdsArrayPath
   */
  DataArrayPath getFeatureIdsArrayPath() const;
  Q_PROPERTY(DataArrayPath FeatureIdsArrayPath READ getFeatureIdsArrayPath WRITE setFeatureIdsArrayPath)

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
  DelaunayTriangulation();

  /**
   * @brief dataCheck Checks for the appropriate parameter values and availability of arrays
   */
  void dataCheck() override;

  /**
   * @brief Initializes all the private instance variables.
   */
  void initialize();

private:
  DataArrayPath m_InputGeometry = {"", "", ""};
  QString m_TriangleDataContainerName = {SIMPL::Defaults::TriangleDataContainerName};
  QString m_VertexAttributeMatrixName = {SIMPL::Defaults::VertexAttributeMatrixName};
  QString m_FaceAttributeMatrixName = {SIMPL::Defaults::FaceAttributeMatrixName};
  double m_Offset = {1.0};
  double m_Tolerance = {0.00001};
  bool m_TriangulateByFeature = {false};
  DataArrayPath m_FeatureIdsArrayPath = {"", "", ""};
  std::weak_ptr<Int32ArrayType> m_FeatureIdsPtr;
  int32_t* m_FeatureIds = nullptr;

  DelaunayTriangulation(const DelaunayTriangulation&) = delete; // Copy Constructor Not Implemented
  DelaunayTriangulation(DelaunayTriangulation&&) = delete;      // Move Constructor Not Implemented
  void operator=(const DelaunayTriangulation&) = delete;        // Operator '=' Not Implemented
};
