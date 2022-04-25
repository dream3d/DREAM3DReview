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

#include "SIMPLib/Common/SIMPLibSetGetMacros.h"
#include "SIMPLib/Filtering/AbstractFilter.h"
#include "SIMPLib/SIMPLib.h"

#include "DREAM3DReview/DREAM3DReviewDLLExport.h"

/**
 * @brief The FindArrayStatistics class. See [Filter documentation](@ref findarraystatistics) for details.
 */
class DREAM3DReview_EXPORT FindArrayStatistics : public AbstractFilter
{
  Q_OBJECT
  PYB11_CREATE_BINDINGS(FindArrayStatistics SUPERCLASS AbstractFilter)
  PYB11_PROPERTY(bool FindLength READ getFindLength WRITE setFindLength)
  PYB11_PROPERTY(bool FindMin READ getFindMin WRITE setFindMin)
  PYB11_PROPERTY(bool FindMax READ getFindMax WRITE setFindMax)
  PYB11_PROPERTY(bool FindMean READ getFindMean WRITE setFindMean)
  PYB11_PROPERTY(bool FindMedian READ getFindMedian WRITE setFindMedian)
  PYB11_PROPERTY(bool FindStdDeviation READ getFindStdDeviation WRITE setFindStdDeviation)
  PYB11_PROPERTY(bool FindSummation READ getFindSummation WRITE setFindSummation)
  PYB11_PROPERTY(bool UseMask READ getUseMask WRITE setUseMask)
  PYB11_PROPERTY(bool StandardizeData READ getStandardizeData WRITE setStandardizeData)
  PYB11_PROPERTY(bool ComputeByIndex READ getComputeByIndex WRITE setComputeByIndex)
  PYB11_PROPERTY(DataArrayPath DestinationAttributeMatrix READ getDestinationAttributeMatrix WRITE setDestinationAttributeMatrix)
  PYB11_PROPERTY(DataArrayPath MaskArrayPath READ getMaskArrayPath WRITE setMaskArrayPath)
  PYB11_PROPERTY(QString LengthArrayName READ getLengthArrayName WRITE setLengthArrayName)
  PYB11_PROPERTY(QString MinimumArrayName READ getMinimumArrayName WRITE setMinimumArrayName)
  PYB11_PROPERTY(QString MaximumArrayName READ getMaximumArrayName WRITE setMaximumArrayName)
  PYB11_PROPERTY(QString MeanArrayName READ getMeanArrayName WRITE setMeanArrayName)
  PYB11_PROPERTY(QString MedianArrayName READ getMedianArrayName WRITE setMedianArrayName)
  PYB11_PROPERTY(QString StdDeviationArrayName READ getStdDeviationArrayName WRITE setStdDeviationArrayName)
  PYB11_PROPERTY(QString SummationArrayName READ getSummationArrayName WRITE setSummationArrayName)
  PYB11_PROPERTY(QString StandardizedArrayName READ getStandardizedArrayName WRITE setStandardizedArrayName)
  PYB11_PROPERTY(DataArrayPath SelectedArrayPath READ getSelectedArrayPath WRITE setSelectedArrayPath)
  PYB11_PROPERTY(DataArrayPath FeatureIdsArrayPath READ getFeatureIdsArrayPath WRITE setFeatureIdsArrayPath)

public:
  SIMPL_SHARED_POINTERS(FindArrayStatistics)
  SIMPL_FILTER_NEW_MACRO(FindArrayStatistics)
  SIMPL_TYPE_MACRO_SUPER_OVERRIDE(FindArrayStatistics, AbstractFilter)

  ~FindArrayStatistics() override;

  SIMPL_FILTER_PARAMETER(bool, FindLength)
  Q_PROPERTY(bool FindLength READ getFindLength WRITE setFindLength)

  SIMPL_FILTER_PARAMETER(bool, FindMin)
  Q_PROPERTY(bool FindMin READ getFindMin WRITE setFindMin)

  SIMPL_FILTER_PARAMETER(bool, FindMax)
  Q_PROPERTY(bool FindMax READ getFindMax WRITE setFindMax)

  SIMPL_FILTER_PARAMETER(bool, FindMean)
  Q_PROPERTY(bool FindMean READ getFindMean WRITE setFindMean)

  SIMPL_FILTER_PARAMETER(bool, FindMedian)
  Q_PROPERTY(bool FindMedian READ getFindMedian WRITE setFindMedian)

  SIMPL_FILTER_PARAMETER(bool, FindStdDeviation)
  Q_PROPERTY(bool FindStdDeviation READ getFindStdDeviation WRITE setFindStdDeviation)

  SIMPL_FILTER_PARAMETER(bool, FindSummation)
  Q_PROPERTY(bool FindSummation READ getFindSummation WRITE setFindSummation)

  SIMPL_FILTER_PARAMETER(bool, UseMask)
  Q_PROPERTY(bool UseMask READ getUseMask WRITE setUseMask)

  SIMPL_FILTER_PARAMETER(bool, StandardizeData)
  Q_PROPERTY(bool StandardizeData READ getStandardizeData WRITE setStandardizeData)

  SIMPL_FILTER_PARAMETER(bool, ComputeByIndex)
  Q_PROPERTY(bool ComputeByIndex READ getComputeByIndex WRITE setComputeByIndex)

  SIMPL_FILTER_PARAMETER(DataArrayPath, DestinationAttributeMatrix)
  Q_PROPERTY(DataArrayPath DestinationAttributeMatrix READ getDestinationAttributeMatrix WRITE setDestinationAttributeMatrix)

  SIMPL_FILTER_PARAMETER(DataArrayPath, MaskArrayPath)
  Q_PROPERTY(DataArrayPath MaskArrayPath READ getMaskArrayPath WRITE setMaskArrayPath)

  SIMPL_FILTER_PARAMETER(QString, LengthArrayName)
  Q_PROPERTY(QString LengthArrayName READ getLengthArrayName WRITE setLengthArrayName)

  SIMPL_FILTER_PARAMETER(QString, MinimumArrayName)
  Q_PROPERTY(QString MinimumArrayName READ getMinimumArrayName WRITE setMinimumArrayName)

  SIMPL_FILTER_PARAMETER(QString, MaximumArrayName)
  Q_PROPERTY(QString MaximumArrayName READ getMaximumArrayName WRITE setMaximumArrayName)

  SIMPL_FILTER_PARAMETER(QString, MeanArrayName)
  Q_PROPERTY(QString MeanArrayName READ getMeanArrayName WRITE setMeanArrayName)

  SIMPL_FILTER_PARAMETER(QString, MedianArrayName)
  Q_PROPERTY(QString MedianArrayName READ getMedianArrayName WRITE setMedianArrayName)

  SIMPL_FILTER_PARAMETER(QString, StdDeviationArrayName)
  Q_PROPERTY(QString StdDeviationArrayName READ getStdDeviationArrayName WRITE setStdDeviationArrayName)

  SIMPL_FILTER_PARAMETER(QString, SummationArrayName)
  Q_PROPERTY(QString SummationArrayName READ getSummationArrayName WRITE setSummationArrayName)

  SIMPL_FILTER_PARAMETER(QString, StandardizedArrayName)
  Q_PROPERTY(QString StandardizedArrayName READ getStandardizedArrayName WRITE setStandardizedArrayName)

  SIMPL_FILTER_PARAMETER(DataArrayPath, SelectedArrayPath)
  Q_PROPERTY(DataArrayPath SelectedArrayPath READ getSelectedArrayPath WRITE setSelectedArrayPath)

  SIMPL_FILTER_PARAMETER(DataArrayPath, FeatureIdsArrayPath)
  Q_PROPERTY(DataArrayPath FeatureIdsArrayPath READ getFeatureIdsArrayPath WRITE setFeatureIdsArrayPath)

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
  /**
   * @brief createCompatibleArrays Creates the output statistics arrays with compatible types based on the
   * user options and incoming DataArray type
   */
  template <typename T> void createCompatibleArrays()
  {
    QVector<size_t> cDims(1, 1);
    DataArrayPath path(getDestinationAttributeMatrix().getDataContainerName(), getDestinationAttributeMatrix().getAttributeMatrixName(), "");

    if(m_FindMin)
    {
      path.setDataArrayName(getMinimumArrayName());
      m_MinimumPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<T>, AbstractFilter, T>(this, path, 0, cDims);
    }

    if(m_FindMax)
    {
      path.setDataArrayName(getMaximumArrayName());
      m_MaximumPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<T>, AbstractFilter, T>(this, path, 0, cDims);
    }

    if(m_FindMean)
    {
      path.setDataArrayName(getMeanArrayName());
      m_MeanPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<float>, AbstractFilter, float>(this, path, 0, cDims);
    }

    if(m_FindMedian)
    {
      path.setDataArrayName(getMedianArrayName());
      m_MedianPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<float>, AbstractFilter, float>(this, path, 0, cDims);
    }

    if(m_FindStdDeviation)
    {
      path.setDataArrayName(getStdDeviationArrayName());
      m_StandardDeviationPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<float>, AbstractFilter, float>(this, path, 0, cDims);
    }

    if(m_FindSummation)
    {
      path.setDataArrayName(getSummationArrayName());
      m_SummationPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<float>, AbstractFilter, float>(this, path, 0, cDims);
    }
  }

  FindArrayStatistics();
  /**
   * @brief dataCheck Checks for the appropriate parameter values and availability of arrays
   */
  void dataCheck();

  /**
   * @brief Initializes all the private instance variables.
   */
  void initialize();

private:
  DEFINE_DATAARRAY_VARIABLE(int64_t, Length)
  DEFINE_IDATAARRAY_VARIABLE(Minimum)
  DEFINE_IDATAARRAY_VARIABLE(Maximum)
  DEFINE_DATAARRAY_VARIABLE(float, Mean)
  DEFINE_DATAARRAY_VARIABLE(float, Median)
  DEFINE_DATAARRAY_VARIABLE(float, StandardDeviation)
  DEFINE_DATAARRAY_VARIABLE(float, Summation)
  DEFINE_DATAARRAY_VARIABLE(float, Standardized)
  DEFINE_IDATAARRAY_VARIABLE(InputArray)
  DEFINE_DATAARRAY_VARIABLE(int32_t, FeatureIds)
  DEFINE_DATAARRAY_VARIABLE(bool, Mask)

public:
  FindArrayStatistics(const FindArrayStatistics&) = delete; // Copy Constructor Not Implemented
  FindArrayStatistics(FindArrayStatistics&&) = delete;      // Move Constructor Not Implemented
  FindArrayStatistics& operator=(const FindArrayStatistics&) = delete; // Copy Assignment Not Implemented
  FindArrayStatistics& operator=(FindArrayStatistics&&) = delete;      // Move Assignment Not Implemented
};

