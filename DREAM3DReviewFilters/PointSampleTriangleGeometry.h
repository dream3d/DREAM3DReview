/* ============================================================================
* Copyright (c) 2009-2016 BlueQuartz Software, LLC
*
* Redistribution and use in source and binary forms, with or without modification,
* are permitted provided that the following conditions are met:
*
* Redistributions of source code must retain the above copyright notice, this
* list of conditions and the following disclaimer.
*
* Redistributions in binary form must reproduce the above copyright notice, this
* list of conditions and the following disclaimer in the documentation and/or
* other materials provided with the distribution.
*
* Neither the name of BlueQuartz Software, the US Air Force, nor the names of its
* contributors may be used to endorse or promote products derived from this software
* without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
* FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
* DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
* SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
* CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
* OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
* USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*
* The code contained herein was partially funded by the followig contracts:
*    United States Air Force Prime Contract FA8650-07-D-5800
*    United States Air Force Prime Contract FA8650-10-D-5210
*    United States Prime Contract Navy N00173-07-C-2068
*
* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#pragma once

#include <memory>

#include <random>

#include "SIMPLib/SIMPLib.h"
#include "SIMPLib/Filtering/AbstractFilter.h"

#include "SIMPLib/Geometry/VertexGeom.h"
#include "SIMPLib/DataArrays/DataArray.hpp"

#include "DREAM3DReview/DREAM3DReviewDLLExport.h"

/**
 * @brief The PointSampleTriangleGeometry class. See [Filter documentation](@ref pointsampletrianglegeometry) for details.
 */
class DREAM3DReview_EXPORT PointSampleTriangleGeometry : public AbstractFilter
{
  Q_OBJECT

  // Start Python bindings declarations
  PYB11_BEGIN_BINDINGS(PointSampleTriangleGeometry SUPERCLASS AbstractFilter)
  PYB11_FILTER()
  PYB11_SHARED_POINTERS(PointSampleTriangleGeometry)
  PYB11_FILTER_NEW_MACRO(PointSampleTriangleGeometry)
  PYB11_PROPERTY(int SamplesNumberType READ getSamplesNumberType WRITE setSamplesNumberType)
  PYB11_PROPERTY(DataArrayPath TriangleGeometry READ getTriangleGeometry WRITE setTriangleGeometry)
  PYB11_PROPERTY(QString VertexGeometry READ getVertexGeometry WRITE setVertexGeometry)
  PYB11_PROPERTY(QString VertexAttributeMatrixName READ getVertexAttributeMatrixName WRITE setVertexAttributeMatrixName)
  PYB11_PROPERTY(int NumberOfSamples READ getNumberOfSamples WRITE setNumberOfSamples)
  PYB11_PROPERTY(DataArrayPath ParentGeometry READ getParentGeometry WRITE setParentGeometry)
  PYB11_PROPERTY(DataArrayPath TriangleAreasArrayPath READ getTriangleAreasArrayPath WRITE setTriangleAreasArrayPath)
  PYB11_PROPERTY(bool UseMask READ getUseMask WRITE setUseMask)
  PYB11_PROPERTY(DataArrayPath MaskArrayPath READ getMaskArrayPath WRITE setMaskArrayPath)
  PYB11_PROPERTY(QVector<DataArrayPath> SelectedDataArrayPaths READ getSelectedDataArrayPaths WRITE setSelectedDataArrayPaths)
  PYB11_END_BINDINGS()
  // End Python bindings declarations

public:
  using Self = PointSampleTriangleGeometry;
  using Pointer = std::shared_ptr<Self>;
  using ConstPointer = std::shared_ptr<const Self>;
  using WeakPointer = std::weak_ptr<Self>;
  using ConstWeakPointer = std::weak_ptr<const Self>;
  static Pointer NullPointer();

  static std::shared_ptr<PointSampleTriangleGeometry> New();

  /**
   * @brief Returns the name of the class for PointSampleTriangleGeometry
   */
  QString getNameOfClass() const override;
  /**
   * @brief Returns the name of the class for PointSampleTriangleGeometry
   */
  static QString ClassName();

  ~PointSampleTriangleGeometry() override;

  /**
   * @brief Setter property for SamplesNumberType
   */
  void setSamplesNumberType(int value);
  /**
   * @brief Getter property for SamplesNumberType
   * @return Value of SamplesNumberType
   */
  int getSamplesNumberType() const;
  Q_PROPERTY(int SamplesNumberType READ getSamplesNumberType WRITE setSamplesNumberType)

  /**
   * @brief Setter property for TriangleGeometry
   */
  void setTriangleGeometry(const DataArrayPath& value);
  /**
   * @brief Getter property for TriangleGeometry
   * @return Value of TriangleGeometry
   */
  DataArrayPath getTriangleGeometry() const;
  Q_PROPERTY(DataArrayPath TriangleGeometry READ getTriangleGeometry WRITE setTriangleGeometry)

  /**
   * @brief Setter property for VertexGeometry
   */
  void setVertexGeometry(const QString& value);
  /**
   * @brief Getter property for VertexGeometry
   * @return Value of VertexGeometry
   */
  QString getVertexGeometry() const;
  Q_PROPERTY(QString VertexGeometry READ getVertexGeometry WRITE setVertexGeometry)

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
   * @brief Setter property for NumberOfSamples
   */
  void setNumberOfSamples(int value);
  /**
   * @brief Getter property for NumberOfSamples
   * @return Value of NumberOfSamples
   */
  int getNumberOfSamples() const;
  Q_PROPERTY(int NumberOfSamples READ getNumberOfSamples WRITE setNumberOfSamples)

  /**
   * @brief Setter property for ParentGeometry
   */
  void setParentGeometry(const DataArrayPath& value);
  /**
   * @brief Getter property for ParentGeometry
   * @return Value of ParentGeometry
   */
  DataArrayPath getParentGeometry() const;
  Q_PROPERTY(DataArrayPath ParentGeometry READ getParentGeometry WRITE setParentGeometry)

  /**
   * @brief Setter property for TriangleAreasArrayPath
   */
  void setTriangleAreasArrayPath(const DataArrayPath& value);
  /**
   * @brief Getter property for TriangleAreasArrayPath
   * @return Value of TriangleAreasArrayPath
   */
  DataArrayPath getTriangleAreasArrayPath() const;
  Q_PROPERTY(DataArrayPath TriangleAreasArrayPath READ getTriangleAreasArrayPath WRITE setTriangleAreasArrayPath)

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
   * @brief Setter property for SelectedDataArrayPaths
   */
  void setSelectedDataArrayPaths(const QVector<DataArrayPath>& value);
  /**
   * @brief Getter property for SelectedDataArrayPaths
   * @return Value of SelectedDataArrayPaths
   */
  QVector<DataArrayPath> getSelectedDataArrayPaths() const;
  Q_PROPERTY(QVector<DataArrayPath> SelectedDataArrayPaths READ getSelectedDataArrayPaths WRITE setSelectedDataArrayPaths)

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
   * @brief readFilterParameters Reimplemented from @see AbstractFilter class
   */
  void readFilterParameters(AbstractFilterParametersReader* reader, int index) override;

  /**
   * @brief execute Reimplemented from @see AbstractFilter class
   */
  void execute() override;

protected:
  PointSampleTriangleGeometry();

  void sampleTriangle(float a[3], float b[3], float c[3], 
                      int64_t curVertex, VertexGeom::Pointer vertex, 
                      int64_t tri, std::mt19937_64& gen, 
                      std::uniform_real_distribution<>& dist);

  /**
   * @brief dataCheck Checks for the appropriate parameter values and availability of arrays
   */
  void dataCheck() override;

  /**
  * @brief Initializes all the private instance variables.
  */
  void initialize();

private:
  std::weak_ptr<DataArray<double>> m_TriangleAreasPtr;
  double* m_TriangleAreas = nullptr;
  std::weak_ptr<DataArray<bool>> m_MaskPtr;
  bool* m_Mask = nullptr;

  int m_SamplesNumberType = {};
  DataArrayPath m_TriangleGeometry = {};
  QString m_VertexGeometry = {};
  QString m_VertexAttributeMatrixName = {};
  int m_NumberOfSamples = {};
  DataArrayPath m_ParentGeometry = {};
  DataArrayPath m_TriangleAreasArrayPath = {};
  bool m_UseMask = {};
  DataArrayPath m_MaskArrayPath = {};
  QVector<DataArrayPath> m_SelectedDataArrayPaths = {};

  std::vector<IDataArray::WeakPointer> m_SelectedWeakPtrVector;
  std::vector<IDataArray::WeakPointer> m_CreatedWeakPtrVector;
  int32_t m_NumSamples;

public:
  PointSampleTriangleGeometry(const PointSampleTriangleGeometry&) = delete; // Copy Constructor Not Implemented
  PointSampleTriangleGeometry(PointSampleTriangleGeometry&&) = delete;      // Move Constructor Not Implemented
  PointSampleTriangleGeometry& operator=(const PointSampleTriangleGeometry&) = delete; // Copy Assignment Not Implemented
  PointSampleTriangleGeometry& operator=(PointSampleTriangleGeometry&&) = delete;      // Move Assignment Not Implemented
};

