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

#pragma once

#include <QMutex>

#include "SIMPLib/SIMPLib.h"
#include "SIMPLib/DataArrays/DataArray.hpp"
#include "SIMPLib/Filtering/AbstractFilter.h"

#include "DREAM3DReview/DREAM3DReviewDLLExport.h"

/**
 * @brief The FindVertexToTriangleDistances class. See [Filter documentation](@ref findvertextotriangledistances) for details.
 */
class DREAM3DReview_EXPORT FindVertexToTriangleDistances : public AbstractFilter
{
  Q_OBJECT

  PYB11_BEGIN_BINDINGS(FindVertexToTriangleDistances SUPERCLASS AbstractFilter)
  PYB11_FILTER()
  PYB11_SHARED_POINTERS(FindVertexToTriangleDistances)
  PYB11_FILTER_NEW_MACRO(FindVertexToTriangleDistances)
  PYB11_PROPERTY(DataArrayPath VertexDataContainer READ getVertexDataContainer WRITE setVertexDataContainer)
  PYB11_PROPERTY(DataArrayPath TriangleDataContainer READ getTriangleDataContainer WRITE setTriangleDataContainer)
  PYB11_PROPERTY(DataArrayPath TriangleNormalsArrayPath READ getTriangleNormalsArrayPath WRITE setTriangleNormalsArrayPath)
  PYB11_PROPERTY(DataArrayPath DistancesArrayPath READ getDistancesArrayPath WRITE setDistancesArrayPath)
  PYB11_PROPERTY(DataArrayPath ClosestTriangleIdArrayPath READ getClosestTriangleIdArrayPath WRITE setClosestTriangleIdArrayPath)
  PYB11_END_BINDINGS()

public:
  using Self = FindVertexToTriangleDistances;
  using Pointer = std::shared_ptr<Self>;
  using ConstPointer = std::shared_ptr<const Self>;
  using WeakPointer = std::weak_ptr<Self>;
  using ConstWeakPointer = std::weak_ptr<const Self>;
  static Pointer NullPointer();

  static std::shared_ptr<FindVertexToTriangleDistances> New();

  /**
   * @brief Returns the name of the class for FindVertexToTriangleDistances
   */
  QString getNameOfClass() const override;

  /**
   * @brief Returns the name of the class for FindVertexToTriangleDistances
   */
  static QString ClassName();

  ~FindVertexToTriangleDistances() override;

  /**
   * @brief Setter property for VertexDataContainer
   */
  void setVertexDataContainer(const DataArrayPath& value);

  /**
   * @brief Getter property for VertexDataContainer
   * @return Value of VertexDataContainer
   */
  DataArrayPath getVertexDataContainer() const;
  Q_PROPERTY(DataArrayPath VertexDataContainer READ getVertexDataContainer WRITE setVertexDataContainer)

  /**
   * @brief Setter property for TriangleDataContainer
   */
  void setTriangleDataContainer(const DataArrayPath& value);

  /**
   * @brief Getter property for TriangleDataContainer
   * @return Value of TriangleDataContainer
   */
  DataArrayPath getTriangleDataContainer() const;
  Q_PROPERTY(DataArrayPath TriangleDataContainer READ getTriangleDataContainer WRITE setTriangleDataContainer)

  /**
   * @brief Setter property for TriangleNormalsArrayPath
   */
  void setTriangleNormalsArrayPath(const DataArrayPath& value);

  /**
   * @brief Getter property for TriangleNormalsArrayPath
   * @return Value of TriangleNormalsArrayPath
   */
  DataArrayPath getTriangleNormalsArrayPath() const;
  Q_PROPERTY(DataArrayPath TriangleNormalsArrayPath READ getTriangleNormalsArrayPath WRITE setTriangleNormalsArrayPath)

  /**
   * @brief Setter property for DistancesArrayPath
   */
  void setDistancesArrayPath(const DataArrayPath& value);

  /**
   * @brief Getter property for DistancesArrayPath
   * @return Value of DistancesArrayPath
   */
  DataArrayPath getDistancesArrayPath() const;
  Q_PROPERTY(DataArrayPath DistancesArrayPath READ getDistancesArrayPath WRITE setDistancesArrayPath)

  /**
   * @brief Setter property for ClosestTriangleIdArrayPath
   */
  void setClosestTriangleIdArrayPath(const DataArrayPath& value);

  /**
   * @brief Getter property for ClosestTriangleIdArrayPath
   * @return Value of ClosestTriangleIdArrayPath
   */
  DataArrayPath getClosestTriangleIdArrayPath() const;
  Q_PROPERTY(DataArrayPath ClosestTriangleIdArrayPath READ getClosestTriangleIdArrayPath WRITE setClosestTriangleIdArrayPath)

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
  FindVertexToTriangleDistances();

  /**
   * @brief dataCheck Checks for the appropriate parameter values and availability of arrays
   */
  void dataCheck() override;

  /**
   * @brief Initializes all the private instance variables.
   */
  void initialize();

private:
  DataArrayPath m_VertexDataContainer = {"", "", ""};
  DataArrayPath m_TriangleDataContainer = {"", "", ""};
  DataArrayPath m_TriangleNormalsArrayPath = {"", "", ""};
  DataArrayPath m_DistancesArrayPath = {"", "", "Distances"};
  DataArrayPath m_ClosestTriangleIdArrayPath = {"", "", "ClosestTriangleId"};
  std::weak_ptr<DoubleArrayType> m_NormalsPtr;
  double* m_Normals = nullptr;
  std::weak_ptr<FloatArrayType> m_DistancesPtr;
  float* m_Distances = nullptr;
  std::weak_ptr<Int64ArrayType> m_ClosestTriangleIdsPtr;
  int64_t* m_ClosestTriangleIds = nullptr;

  /**
   * @brief distance
   * @param vec1
   * @param vec2
   * @return
   */
  float distance(const std::vector<float>& vec1, const std::vector<float>& vec2);

  /**
   * @brief point_segment_distance
   * @param x0
   * @param x1
   * @param x2
   * @return
   */
  float point_segment_distance(const std::vector<float>& x0, const std::vector<float>& x1, const std::vector<float>& x2);

  /**
   * @brief point_triangle_distance
   * @param x0
   * @param x1
   * @param x2
   * @param x3
   * @return
   */
  float point_triangle_distance(const std::vector<float>& x0, const std::vector<float>& x1, const std::vector<float>& x2, const std::vector<float>& x3, const int64_t& triangle);

  /**
   * @brief sendThreadSafeProgressMessage
   * @param counter
   * @param max
   */
  void sendThreadSafeProgressMessage(int64_t counter);

  QMutex m_Mutex;
  int64_t m_ProgressCounter = {0};
  int64_t m_TotalElements = {0};
  int64_t m_LastProgressInt = {0};

  friend class FindVertexToTriangleDistancesImpl;

public:
  FindVertexToTriangleDistances(const FindVertexToTriangleDistances&) = delete;            // Copy Constructor Not Implemented
  FindVertexToTriangleDistances(FindVertexToTriangleDistances&&) = delete;                 // Move Constructor Not Implemented
  FindVertexToTriangleDistances& operator=(const FindVertexToTriangleDistances&) = delete; // Copy Assignment Not Implemented
  FindVertexToTriangleDistances& operator=(FindVertexToTriangleDistances&&) = delete;      // Move Assignment Not Implemented
};
