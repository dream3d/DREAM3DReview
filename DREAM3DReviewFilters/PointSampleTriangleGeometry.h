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

#include <random>

#include "SIMPLib/Common/SIMPLibSetGetMacros.h"
#include "SIMPLib/Filtering/AbstractFilter.h"
#include "SIMPLib/SIMPLib.h"

#include "SIMPLib/Geometry/VertexGeom.h"

#include "DREAM3DReview/DREAM3DReviewDLLExport.h"

/**
 * @brief The PointSampleTriangleGeometry class. See [Filter documentation](@ref pointsampletrianglegeometry) for details.
 */
class DREAM3DReview_EXPORT PointSampleTriangleGeometry : public AbstractFilter
{
  Q_OBJECT
  PYB11_CREATE_BINDINGS(PointSampleTriangleGeometry SUPERCLASS AbstractFilter)
  PYB11_PROPERTY(int SamplesNumberType READ getSamplesNumberType WRITE setSamplesNumberType)
  PYB11_PROPERTY(QString TriangleGeometry READ getTriangleGeometry WRITE setTriangleGeometry)
  PYB11_PROPERTY(QString VertexGeometry READ getVertexGeometry WRITE setVertexGeometry)
  PYB11_PROPERTY(QString VertexAttributeMatrixName READ getVertexAttributeMatrixName WRITE setVertexAttributeMatrixName)
  PYB11_PROPERTY(int NumberOfSamples READ getNumberOfSamples WRITE setNumberOfSamples)
  PYB11_PROPERTY(QString ParentGeometry READ getParentGeometry WRITE setParentGeometry)
  PYB11_PROPERTY(DataArrayPath TriangleAreasArrayPath READ getTriangleAreasArrayPath WRITE setTriangleAreasArrayPath)
  PYB11_PROPERTY(bool UseMask READ getUseMask WRITE setUseMask)
  PYB11_PROPERTY(DataArrayPath MaskArrayPath READ getMaskArrayPath WRITE setMaskArrayPath)
  PYB11_PROPERTY(QVector<DataArrayPath> SelectedDataArrayPaths READ getSelectedDataArrayPaths WRITE setSelectedDataArrayPaths)

public:
  SIMPL_SHARED_POINTERS(PointSampleTriangleGeometry)
  SIMPL_FILTER_NEW_MACRO(PointSampleTriangleGeometry)
  SIMPL_TYPE_MACRO_SUPER_OVERRIDE(PointSampleTriangleGeometry, AbstractFilter)

  ~PointSampleTriangleGeometry() override;

  SIMPL_FILTER_PARAMETER(int, SamplesNumberType)
  Q_PROPERTY(int SamplesNumberType READ getSamplesNumberType WRITE setSamplesNumberType)

  SIMPL_FILTER_PARAMETER(QString, TriangleGeometry)
  Q_PROPERTY(QString TriangleGeometry READ getTriangleGeometry WRITE setTriangleGeometry)

  SIMPL_FILTER_PARAMETER(QString, VertexGeometry)
  Q_PROPERTY(QString VertexGeometry READ getVertexGeometry WRITE setVertexGeometry)

  SIMPL_FILTER_PARAMETER(QString, VertexAttributeMatrixName)
  Q_PROPERTY(QString VertexAttributeMatrixName READ getVertexAttributeMatrixName WRITE setVertexAttributeMatrixName)

  SIMPL_FILTER_PARAMETER(int, NumberOfSamples)
  Q_PROPERTY(int NumberOfSamples READ getNumberOfSamples WRITE setNumberOfSamples)

  SIMPL_FILTER_PARAMETER(QString, ParentGeometry)
  Q_PROPERTY(QString ParentGeometry READ getParentGeometry WRITE setParentGeometry)

  SIMPL_FILTER_PARAMETER(DataArrayPath, TriangleAreasArrayPath)
  Q_PROPERTY(DataArrayPath TriangleAreasArrayPath READ getTriangleAreasArrayPath WRITE setTriangleAreasArrayPath)

  SIMPL_FILTER_PARAMETER(bool, UseMask)
  Q_PROPERTY(bool UseMask READ getUseMask WRITE setUseMask)

  SIMPL_FILTER_PARAMETER(DataArrayPath, MaskArrayPath)
  Q_PROPERTY(DataArrayPath MaskArrayPath READ getMaskArrayPath WRITE setMaskArrayPath)

  SIMPL_FILTER_PARAMETER(QVector<DataArrayPath>, SelectedDataArrayPaths)
  Q_PROPERTY(QVector<DataArrayPath> SelectedDataArrayPaths READ getSelectedDataArrayPaths WRITE setSelectedDataArrayPaths)

  /**
   * @brief getCompiledLibraryName Reimplemented from @see AbstractFilter class
   */
  const QString getCompiledLibraryName() const override;

  /**
   * @brief getBrandingString Returns the branding string for the filter, which is a tag
   * used to denote the filter's association with specific plugins
   * @return Branding string
  */
  const QString getBrandingString() const override;

  /**
   * @brief getFilterVersion Returns a version string for this filter. Default
   * value is an empty string.
   * @return
   */
  const QString getFilterVersion() const override;

  /**
   * @brief newFilterInstance Reimplemented from @see AbstractFilter class
   */
  AbstractFilter::Pointer newFilterInstance(bool copyFilterParameters) const override;

  /**
   * @brief getGroupName Reimplemented from @see AbstractFilter class
   */
  const QString getGroupName() const override;

  /**
   * @brief getSubGroupName Reimplemented from @see AbstractFilter class
   */
  const QString getSubGroupName() const override;

  /**
   * @brief getUuid Return the unique identifier for this filter.
   * @return A QUuid object.
   */
  const QUuid getUuid() override;

  /**
   * @brief getHumanLabel Reimplemented from @see AbstractFilter class
   */
  const QString getHumanLabel() const override;

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
  * @brief preflight Reimplemented from @see AbstractFilter class
  */
  void preflight() override;

Q_SIGNALS:
  /**
   * @brief updateFilterParameters Emitted when the Filter requests all the latest Filter parameters
   * be pushed from a user-facing control (such as a widget)
   * @param filter Filter instance pointer
   */
  void updateFilterParameters(AbstractFilter* filter);

  /**
   * @brief parametersChanged Emitted when any Filter parameter is changed internally
   */
  void parametersChanged();

  /**
   * @brief preflightAboutToExecute Emitted just before calling dataCheck()
   */
  void preflightAboutToExecute();

  /**
   * @brief preflightExecuted Emitted just after calling dataCheck()
   */
  void preflightExecuted();

protected:
  PointSampleTriangleGeometry();

  void sampleTriangle(float a[3], float b[3], float c[3], 
                      int64_t curVertex, VertexGeom::Pointer vertex, 
                      int64_t tri, std::mt19937_64& gen, 
                      std::uniform_real_distribution<>& dist);

  /**
   * @brief dataCheck Checks for the appropriate parameter values and availability of arrays
   */
  void dataCheck();

  /**
  * @brief Initializes all the private instance variables.
  */
  void initialize();

private:
  DEFINE_DATAARRAY_VARIABLE(double, TriangleAreas)
  DEFINE_DATAARRAY_VARIABLE(bool, Mask)

  std::vector<IDataArray::WeakPointer> m_SelectedWeakPtrVector;
  std::vector<IDataArray::WeakPointer> m_CreatedWeakPtrVector;
  int32_t m_NumSamples;

public:
  PointSampleTriangleGeometry(const PointSampleTriangleGeometry&) = delete; // Copy Constructor Not Implemented
  PointSampleTriangleGeometry(PointSampleTriangleGeometry&&) = delete;      // Move Constructor Not Implemented
  PointSampleTriangleGeometry& operator=(const PointSampleTriangleGeometry&) = delete; // Copy Assignment Not Implemented
  PointSampleTriangleGeometry& operator=(PointSampleTriangleGeometry&&) = delete;      // Move Assignment Not Implemented
};

