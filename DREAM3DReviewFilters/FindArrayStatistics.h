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

#include "SIMPLib/SIMPLib.h"
#include "SIMPLib/Filtering/AbstractFilter.h"
#include "SIMPLib/DataArrays/DataArray.hpp"
#include "SIMPLib/DataContainers/DataContainerArray.h"

class IDataArray;
using IDataArrayWkPtrType = std::weak_ptr<IDataArray>;

#include "DREAM3DReview/DREAM3DReviewDLLExport.h"

/**
 * @brief The FindArrayStatistics class. See [Filter documentation](@ref findarraystatistics) for details.
 */
class DREAM3DReview_EXPORT FindArrayStatistics : public AbstractFilter
{
  Q_OBJECT

#ifdef SIMPL_ENABLE_PYTHON
  PYB11_CREATE_BINDINGS(FindArrayStatistics SUPERCLASS AbstractFilter)
  PYB11_SHARED_POINTERS(FindArrayStatistics)
  PYB11_FILTER_NEW_MACRO(FindArrayStatistics)
  PYB11_FILTER_PARAMETER(bool, FindLength)
  PYB11_FILTER_PARAMETER(bool, FindMin)
  PYB11_FILTER_PARAMETER(bool, FindMax)
  PYB11_FILTER_PARAMETER(bool, FindMean)
  PYB11_FILTER_PARAMETER(bool, FindMedian)
  PYB11_FILTER_PARAMETER(bool, FindStdDeviation)
  PYB11_FILTER_PARAMETER(bool, FindSummation)
  PYB11_FILTER_PARAMETER(bool, UseMask)
  PYB11_FILTER_PARAMETER(bool, StandardizeData)
  PYB11_FILTER_PARAMETER(bool, ComputeByIndex)
  PYB11_FILTER_PARAMETER(DataArrayPath, DestinationAttributeMatrix)
  PYB11_FILTER_PARAMETER(DataArrayPath, MaskArrayPath)
  PYB11_FILTER_PARAMETER(QString, LengthArrayName)
  PYB11_FILTER_PARAMETER(QString, MinimumArrayName)
  PYB11_FILTER_PARAMETER(QString, MaximumArrayName)
  PYB11_FILTER_PARAMETER(QString, MeanArrayName)
  PYB11_FILTER_PARAMETER(QString, MedianArrayName)
  PYB11_FILTER_PARAMETER(QString, StdDeviationArrayName)
  PYB11_FILTER_PARAMETER(QString, SummationArrayName)
  PYB11_FILTER_PARAMETER(QString, StandardizedArrayName)
  PYB11_FILTER_PARAMETER(DataArrayPath, SelectedArrayPath)
  PYB11_FILTER_PARAMETER(DataArrayPath, FeatureIdsArrayPath)
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
#endif

public:
  using Self = FindArrayStatistics;
  using Pointer = std::shared_ptr<Self>;
  using ConstPointer = std::shared_ptr<const Self>;
  using WeakPointer = std::weak_ptr<Self>;
  using ConstWeakPointer = std::weak_ptr<const Self>;
  static Pointer NullPointer();

  static std::shared_ptr<FindArrayStatistics> New();

  /**
   * @brief Returns the name of the class for FindArrayStatistics
   */
  QString getNameOfClass() const override;
  /**
   * @brief Returns the name of the class for FindArrayStatistics
   */
  static QString ClassName();

  ~FindArrayStatistics() override;

  /**
   * @brief Setter property for FindLength
   */
  void setFindLength(bool value);
  /**
   * @brief Getter property for FindLength
   * @return Value of FindLength
   */
  bool getFindLength() const;

  Q_PROPERTY(bool FindLength READ getFindLength WRITE setFindLength)

  /**
   * @brief Setter property for FindMin
   */
  void setFindMin(bool value);
  /**
   * @brief Getter property for FindMin
   * @return Value of FindMin
   */
  bool getFindMin() const;

  Q_PROPERTY(bool FindMin READ getFindMin WRITE setFindMin)

  /**
   * @brief Setter property for FindMax
   */
  void setFindMax(bool value);
  /**
   * @brief Getter property for FindMax
   * @return Value of FindMax
   */
  bool getFindMax() const;

  Q_PROPERTY(bool FindMax READ getFindMax WRITE setFindMax)

  /**
   * @brief Setter property for FindMean
   */
  void setFindMean(bool value);
  /**
   * @brief Getter property for FindMean
   * @return Value of FindMean
   */
  bool getFindMean() const;

  Q_PROPERTY(bool FindMean READ getFindMean WRITE setFindMean)

  /**
   * @brief Setter property for FindMedian
   */
  void setFindMedian(bool value);
  /**
   * @brief Getter property for FindMedian
   * @return Value of FindMedian
   */
  bool getFindMedian() const;

  Q_PROPERTY(bool FindMedian READ getFindMedian WRITE setFindMedian)

  /**
   * @brief Setter property for FindStdDeviation
   */
  void setFindStdDeviation(bool value);
  /**
   * @brief Getter property for FindStdDeviation
   * @return Value of FindStdDeviation
   */
  bool getFindStdDeviation() const;

  Q_PROPERTY(bool FindStdDeviation READ getFindStdDeviation WRITE setFindStdDeviation)

  /**
   * @brief Setter property for FindSummation
   */
  void setFindSummation(bool value);
  /**
   * @brief Getter property for FindSummation
   * @return Value of FindSummation
   */
  bool getFindSummation() const;

  Q_PROPERTY(bool FindSummation READ getFindSummation WRITE setFindSummation)

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
   * @brief Setter property for StandardizeData
   */
  void setStandardizeData(bool value);
  /**
   * @brief Getter property for StandardizeData
   * @return Value of StandardizeData
   */
  bool getStandardizeData() const;

  Q_PROPERTY(bool StandardizeData READ getStandardizeData WRITE setStandardizeData)

  /**
   * @brief Setter property for ComputeByIndex
   */
  void setComputeByIndex(bool value);
  /**
   * @brief Getter property for ComputeByIndex
   * @return Value of ComputeByIndex
   */
  bool getComputeByIndex() const;

  Q_PROPERTY(bool ComputeByIndex READ getComputeByIndex WRITE setComputeByIndex)

  /**
   * @brief Setter property for DestinationAttributeMatrix
   */
  void setDestinationAttributeMatrix(const DataArrayPath& value);
  /**
   * @brief Getter property for DestinationAttributeMatrix
   * @return Value of DestinationAttributeMatrix
   */
  DataArrayPath getDestinationAttributeMatrix() const;

  Q_PROPERTY(DataArrayPath DestinationAttributeMatrix READ getDestinationAttributeMatrix WRITE setDestinationAttributeMatrix)

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
   * @brief Setter property for LengthArrayName
   */
  void setLengthArrayName(const QString& value);
  /**
   * @brief Getter property for LengthArrayName
   * @return Value of LengthArrayName
   */
  QString getLengthArrayName() const;

  Q_PROPERTY(QString LengthArrayName READ getLengthArrayName WRITE setLengthArrayName)

  /**
   * @brief Setter property for MinimumArrayName
   */
  void setMinimumArrayName(const QString& value);
  /**
   * @brief Getter property for MinimumArrayName
   * @return Value of MinimumArrayName
   */
  QString getMinimumArrayName() const;

  Q_PROPERTY(QString MinimumArrayName READ getMinimumArrayName WRITE setMinimumArrayName)

  /**
   * @brief Setter property for MaximumArrayName
   */
  void setMaximumArrayName(const QString& value);
  /**
   * @brief Getter property for MaximumArrayName
   * @return Value of MaximumArrayName
   */
  QString getMaximumArrayName() const;

  Q_PROPERTY(QString MaximumArrayName READ getMaximumArrayName WRITE setMaximumArrayName)

  /**
   * @brief Setter property for MeanArrayName
   */
  void setMeanArrayName(const QString& value);
  /**
   * @brief Getter property for MeanArrayName
   * @return Value of MeanArrayName
   */
  QString getMeanArrayName() const;

  Q_PROPERTY(QString MeanArrayName READ getMeanArrayName WRITE setMeanArrayName)

  /**
   * @brief Setter property for MedianArrayName
   */
  void setMedianArrayName(const QString& value);
  /**
   * @brief Getter property for MedianArrayName
   * @return Value of MedianArrayName
   */
  QString getMedianArrayName() const;

  Q_PROPERTY(QString MedianArrayName READ getMedianArrayName WRITE setMedianArrayName)

  /**
   * @brief Setter property for StdDeviationArrayName
   */
  void setStdDeviationArrayName(const QString& value);
  /**
   * @brief Getter property for StdDeviationArrayName
   * @return Value of StdDeviationArrayName
   */
  QString getStdDeviationArrayName() const;

  Q_PROPERTY(QString StdDeviationArrayName READ getStdDeviationArrayName WRITE setStdDeviationArrayName)

  /**
   * @brief Setter property for SummationArrayName
   */
  void setSummationArrayName(const QString& value);
  /**
   * @brief Getter property for SummationArrayName
   * @return Value of SummationArrayName
   */
  QString getSummationArrayName() const;

  Q_PROPERTY(QString SummationArrayName READ getSummationArrayName WRITE setSummationArrayName)

  /**
   * @brief Setter property for StandardizedArrayName
   */
  void setStandardizedArrayName(const QString& value);
  /**
   * @brief Getter property for StandardizedArrayName
   * @return Value of StandardizedArrayName
   */
  QString getStandardizedArrayName() const;

  Q_PROPERTY(QString StandardizedArrayName READ getStandardizedArrayName WRITE setStandardizedArrayName)

  /**
   * @brief Setter property for SelectedArrayPath
   */
  void setSelectedArrayPath(const DataArrayPath& value);
  /**
   * @brief Getter property for SelectedArrayPath
   * @return Value of SelectedArrayPath
   */
  DataArrayPath getSelectedArrayPath() const;

  Q_PROPERTY(DataArrayPath SelectedArrayPath READ getSelectedArrayPath WRITE setSelectedArrayPath)

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

  /**
  * @brief preflight Reimplemented from @see AbstractFilter class
  */
  void preflight() override;

signals:
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
    std::vector<size_t> cDims(1, 1);
    DataArrayPath path(getDestinationAttributeMatrix().getDataContainerName(), getDestinationAttributeMatrix().getAttributeMatrixName(), "");

    if(m_FindMin)
    {
      path.setDataArrayName(getMinimumArrayName());
      m_MinimumPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<T>>(this, path, 0, cDims);
    }

    if(m_FindMax)
    {
      path.setDataArrayName(getMaximumArrayName());
      m_MaximumPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<T>>(this, path, 0, cDims);
    }

    if(m_FindMean)
    {
      path.setDataArrayName(getMeanArrayName());
      m_MeanPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<float>>(this, path, 0, cDims);
    }

    if(m_FindMedian)
    {
      path.setDataArrayName(getMedianArrayName());
      m_MedianPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<float>>(this, path, 0, cDims);
    }

    if(m_FindStdDeviation)
    {
      path.setDataArrayName(getStdDeviationArrayName());
      m_StandardDeviationPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<float>>(this, path, 0, cDims);
    }

    if(m_FindSummation)
    {
      path.setDataArrayName(getSummationArrayName());
      m_SummationPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<float>>(this, path, 0, cDims);
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
  IDataArrayWkPtrType m_MinimumPtr;
  void* m_Minimum = nullptr;
  IDataArrayWkPtrType m_MaximumPtr;
  void* m_Maximum = nullptr;
  IDataArrayWkPtrType m_InputArrayPtr;
  void* m_InputArray = nullptr;

  std::weak_ptr<DataArray<size_t>> m_LengthPtr;
  size_t* m_Length = nullptr;
  std::weak_ptr<DataArray<float>> m_MeanPtr;
  float* m_Mean = nullptr;
  std::weak_ptr<DataArray<float>> m_MedianPtr;
  float* m_Median = nullptr;
  std::weak_ptr<DataArray<float>> m_StandardDeviationPtr;
  float* m_StandardDeviation = nullptr;
  std::weak_ptr<DataArray<float>> m_SummationPtr;
  float* m_Summation = nullptr;
  std::weak_ptr<DataArray<float>> m_StandardizedPtr;
  float* m_Standardized = nullptr;
  std::weak_ptr<DataArray<int32_t>> m_FeatureIdsPtr;
  int32_t* m_FeatureIds = nullptr;
  std::weak_ptr<DataArray<bool>> m_MaskPtr;
  bool* m_Mask = nullptr;

  bool m_FindLength = {};
  bool m_FindMin = {};
  bool m_FindMax = {};
  bool m_FindMean = {};
  bool m_FindMedian = {};
  bool m_FindStdDeviation = {};
  bool m_FindSummation = {};
  bool m_UseMask = {};
  bool m_StandardizeData = {};
  bool m_ComputeByIndex = {};
  DataArrayPath m_DestinationAttributeMatrix = {};
  DataArrayPath m_MaskArrayPath = {};
  QString m_LengthArrayName = {};
  QString m_MinimumArrayName = {};
  QString m_MaximumArrayName = {};
  QString m_MeanArrayName = {};
  QString m_MedianArrayName = {};
  QString m_StdDeviationArrayName = {};
  QString m_SummationArrayName = {};
  QString m_StandardizedArrayName = {};
  DataArrayPath m_SelectedArrayPath = {};
  DataArrayPath m_FeatureIdsArrayPath = {};

public:
  FindArrayStatistics(const FindArrayStatistics&) = delete; // Copy Constructor Not Implemented
  FindArrayStatistics(FindArrayStatistics&&) = delete;      // Move Constructor Not Implemented
  FindArrayStatistics& operator=(const FindArrayStatistics&) = delete; // Copy Assignment Not Implemented
  FindArrayStatistics& operator=(FindArrayStatistics&&) = delete;      // Move Assignment Not Implemented
};

