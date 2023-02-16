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
 * The code contained herein was partially funded by the following contracts:
 *    United States Air Force Prime Contract FA8650-07-D-5800
 *    United States Air Force Prime Contract FA8650-10-D-5210
 *    United States Prime Contract Navy N00173-07-C-2068
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
#include "FindArrayStatistics.h"

#include <cmath>
#include <cstring>
#include <functional>
#include <unordered_map>

#include <QtCore/QTextStream>

#include "SIMPLib/Common/Constants.h"
#include "SIMPLib/Common/TemplateHelpers.h"
#include "SIMPLib/DataArrays/NeighborList.hpp"
#include "SIMPLib/FilterParameters/AbstractFilterParametersReader.h"
#include "SIMPLib/FilterParameters/AttributeMatrixSelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/BooleanFilterParameter.h"
#include "SIMPLib/FilterParameters/DataArraySelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/DoubleFilterParameter.h"
#include "SIMPLib/FilterParameters/IntFilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedBooleanFilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedPathCreationFilterParameter.h"
#include "SIMPLib/FilterParameters/SeparatorFilterParameter.h"

#include "DREAM3DReview/DREAM3DReviewConstants.h"
#include "DREAM3DReview/DREAM3DReviewVersion.h"

#define STATISTICS_FILTER_CLASS_NAME FindArrayStatistics
#include "util/StatisticsHelpers.hpp"

#ifdef SIMPL_USE_PARALLEL_ALGORITHMS
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/partitioner.h>
#endif

// -----------------------------------------------------------------------------
FindArrayStatistics::FindArrayStatistics() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
FindArrayStatistics::~FindArrayStatistics() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FindArrayStatistics::setupFilterParameters()
{
  FilterParameterVectorType parameters;
  parameters.push_back(SeparatorFilterParameter::Create("Statistics Options", FilterParameter::Category::Parameter));
  std::vector<QString> linkedProps;
  linkedProps.push_back("HistogramArrayName");
  linkedProps.push_back("UseFullRange");
  linkedProps.push_back("NumBins");
  linkedProps.push_back("MinRange");
  linkedProps.push_back("MaxRange");
  parameters.push_back(SIMPL_NEW_LINKED_BOOL_FP("Find Histogram", FindHistogram, FilterParameter::Category::Parameter, FindArrayStatistics, linkedProps));
  parameters.push_back(SIMPL_NEW_DOUBLE_FP("Histogram Min Value", MinRange, FilterParameter::Category::Parameter, FindArrayStatistics));
  parameters.push_back(SIMPL_NEW_DOUBLE_FP("Histogram Max Value", MaxRange, FilterParameter::Category::Parameter, FindArrayStatistics));
  parameters.push_back(SIMPL_NEW_BOOL_FP("Use Full Range for Histogram", UseFullRange, FilterParameter::Category::Parameter, FindArrayStatistics));
  parameters.push_back(SIMPL_NEW_INTEGER_FP("Number of Bins", NumBins, FilterParameter::Category::Parameter, FindArrayStatistics));

  linkedProps.clear();
  linkedProps.push_back("LengthArrayName");
  parameters.push_back(SIMPL_NEW_LINKED_BOOL_FP("Find Length", FindLength, FilterParameter::Category::Parameter, FindArrayStatistics, linkedProps));
  linkedProps.clear();
  linkedProps.push_back("MinimumArrayName");
  parameters.push_back(SIMPL_NEW_LINKED_BOOL_FP("Find Minimum", FindMin, FilterParameter::Category::Parameter, FindArrayStatistics, linkedProps));
  linkedProps.clear();
  linkedProps.push_back("MaximumArrayName");
  parameters.push_back(SIMPL_NEW_LINKED_BOOL_FP("Find Maximum", FindMax, FilterParameter::Category::Parameter, FindArrayStatistics, linkedProps));
  linkedProps.clear();
  linkedProps.push_back("MeanArrayName");
  parameters.push_back(SIMPL_NEW_LINKED_BOOL_FP("Find Mean", FindMean, FilterParameter::Category::Parameter, FindArrayStatistics, linkedProps));
  linkedProps.clear();
  linkedProps.push_back("MedianArrayName");
  parameters.push_back(SIMPL_NEW_LINKED_BOOL_FP("Find Median", FindMedian, FilterParameter::Category::Parameter, FindArrayStatistics, linkedProps));
  linkedProps.clear();
  linkedProps.push_back("StdDeviationArrayName");
  parameters.push_back(SIMPL_NEW_LINKED_BOOL_FP("Find Standard Deviation", FindStdDeviation, FilterParameter::Category::Parameter, FindArrayStatistics, linkedProps));
  linkedProps.clear();
  linkedProps.push_back("SummationArrayName");
  parameters.push_back(SIMPL_NEW_LINKED_BOOL_FP("Find Summation", FindSummation, FilterParameter::Category::Parameter, FindArrayStatistics, linkedProps));
  linkedProps.clear();

  parameters.push_back(SeparatorFilterParameter::Create("Algorithm Options", FilterParameter::Category::Parameter));
  linkedProps.push_back("MaskArrayPath");
  parameters.push_back(SIMPL_NEW_LINKED_BOOL_FP("Use Mask", UseMask, FilterParameter::Category::Parameter, FindArrayStatistics, linkedProps));
  linkedProps.clear();
  linkedProps.push_back("FeatureIdsArrayPath");
  parameters.push_back(SIMPL_NEW_LINKED_BOOL_FP("Compute Statistics Per Feature/Ensemble", ComputeByIndex, FilterParameter::Category::Parameter, FindArrayStatistics, linkedProps));
  linkedProps.clear();
  linkedProps.push_back("StandardizedArrayName");
  parameters.push_back(SIMPL_NEW_LINKED_BOOL_FP("Standardize Data", StandardizeData, FilterParameter::Category::Parameter, FindArrayStatistics, linkedProps));

  DataArraySelectionFilterParameter::RequirementType dasReq = DataArraySelectionFilterParameter::CreateRequirement(SIMPL::Defaults::AnyPrimitive, 1, AttributeMatrix::Type::Any, IGeometry::Type::Any);
  parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Attribute Array to Compute Statistics", SelectedArrayPath, FilterParameter::Category::RequiredArray, FindArrayStatistics, dasReq));

  dasReq = DataArraySelectionFilterParameter::CreateRequirement(SIMPL::TypeNames::Int32, 1, AttributeMatrix::Type::Any, IGeometry::Type::Any);
  parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Feature Ids", FeatureIdsArrayPath, FilterParameter::Category::RequiredArray, FindArrayStatistics, dasReq));
  DataArraySelectionFilterParameter::RequirementType req = DataArraySelectionFilterParameter::CreateRequirement(SIMPL::TypeNames::Bool, 1, AttributeMatrix::Type::Any, IGeometry::Type::Any);
  parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Mask", MaskArrayPath, FilterParameter::Category::RequiredArray, FindArrayStatistics, req));

  AttributeMatrixSelectionFilterParameter::RequirementType amReq = AttributeMatrixSelectionFilterParameter::CreateRequirement(AttributeMatrix::Type::Any, IGeometry::Type::Any);
  parameters.push_back(SIMPL_NEW_AM_SELECTION_FP("Destination Attribute Matrix", DestinationAttributeMatrix, FilterParameter::Category::RequiredArray, FindArrayStatistics, amReq));

  parameters.push_back(
      SIMPL_NEW_DA_WITH_LINKED_AM_FP("Histogram", HistogramArrayName, DestinationAttributeMatrix, DestinationAttributeMatrix, FilterParameter::Category::CreatedArray, FindArrayStatistics));
  parameters.push_back(SIMPL_NEW_DA_WITH_LINKED_AM_FP("Length", LengthArrayName, DestinationAttributeMatrix, DestinationAttributeMatrix, FilterParameter::Category::CreatedArray, FindArrayStatistics));
  parameters.push_back(
      SIMPL_NEW_DA_WITH_LINKED_AM_FP("Minimum", MinimumArrayName, DestinationAttributeMatrix, DestinationAttributeMatrix, FilterParameter::Category::CreatedArray, FindArrayStatistics));
  parameters.push_back(
      SIMPL_NEW_DA_WITH_LINKED_AM_FP("Maximum", MaximumArrayName, DestinationAttributeMatrix, DestinationAttributeMatrix, FilterParameter::Category::CreatedArray, FindArrayStatistics));
  parameters.push_back(SIMPL_NEW_DA_WITH_LINKED_AM_FP("Mean", MeanArrayName, DestinationAttributeMatrix, DestinationAttributeMatrix, FilterParameter::Category::CreatedArray, FindArrayStatistics));
  parameters.push_back(SIMPL_NEW_DA_WITH_LINKED_AM_FP("Median", MedianArrayName, DestinationAttributeMatrix, DestinationAttributeMatrix, FilterParameter::Category::CreatedArray, FindArrayStatistics));
  parameters.push_back(SIMPL_NEW_DA_WITH_LINKED_AM_FP("Standard Deviation", StdDeviationArrayName, DestinationAttributeMatrix, DestinationAttributeMatrix, FilterParameter::Category::CreatedArray,
                                                      FindArrayStatistics));
  parameters.push_back(
      SIMPL_NEW_DA_WITH_LINKED_AM_FP("Summation", SummationArrayName, DestinationAttributeMatrix, DestinationAttributeMatrix, FilterParameter::Category::CreatedArray, FindArrayStatistics));
  parameters.push_back(SIMPL_NEW_DA_WITH_LINKED_AM_FP("Standardized Data", StandardizedArrayName, SelectedArrayPath, SelectedArrayPath, FilterParameter::Category::CreatedArray, FindArrayStatistics));

  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FindArrayStatistics::initialize()
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FindArrayStatistics::dataCheck()
{
  clearErrorCode();
  clearWarningCode();

  if(!getFindHistogram() && !getFindMin() && !getFindMax() && !getFindMean() && !getFindMedian() && !getFindStdDeviation() && !getFindSummation() && !getFindLength())
  {
    QString ss = QObject::tr("No statistics have been selected, so this filter will perform no operations");
    setWarningCondition(-701, ss);
    return;
  }

  QVector<DataArrayPath> outputDataArrayPaths;
  QVector<DataArrayPath> inputDataArrayPaths;

  m_InputArrayPtr = getDataContainerArray()->getPrereqIDataArrayFromPath(this, getSelectedArrayPath());

  if(getErrorCode() < 0)
  {
    return;
  }

  inputDataArrayPaths.push_back(getSelectedArrayPath());

  if(m_InputArrayPtr.lock()->getNumberOfComponents() != 1)
  {
    QString ss = QObject::tr("Input Attribute Array must be a scalar array");
    setErrorCondition(-11002, ss);
  }

  if(!getComputeByIndex())
  {
    AttributeMatrix::Pointer destAttrMat = getDataContainerArray()->getPrereqAttributeMatrixFromPath(this, getDestinationAttributeMatrix(), -301);
    if(getErrorCode() < 0)
    {
      return;
    }
    std::vector<size_t> tDims = destAttrMat->getTupleDimensions();
    if(tDims.size() != 1)
    {
      QString ss =
          QObject::tr("Since option \"Compute Statistics Per Feature/Ensemble Id\" is not selected, a single value, representative of the whole array, will be computed for each scalar statistic. "
                      "The selected destination Attribute Matrix (%1) must then have exactly 1 dimension, but the current selection has dimensions %2. "
                      "Consider creating a new Generic Attribute Matrix with scalar tuple dimensions to store the statistics.")
              .arg(getDestinationAttributeMatrix().getAttributeMatrixName())
              .arg(tDims.size());
      setErrorCondition(-11002, ss);
      return;
    }
    if(tDims[0] != 1)
    {
      QString ss =
          QObject::tr("Since option \"Compute Statistics Per Feature/Ensemble Id\" is not selected, a single value, representative of the whole array, will be computed for each scalar statistic. "
                      "The selected destination Attribute Matrix (%1) must then have an extent of 1 in its single dimension , but the current extent is %2. "
                      "Consider creating a new Generic Attribute Matrix with scalar tuple dimensions to store the statistics.")
              .arg(getDestinationAttributeMatrix().getAttributeMatrixName())
              .arg(tDims[0]);
      setErrorCondition(-11002, ss);
      return;
    }
  }

  std::vector<size_t> cDims = {1};

  if(getComputeByIndex())
  {
    m_FeatureIdsPtr = getDataContainerArray()->getPrereqArrayFromPath<Int32ArrayType>(this, getFeatureIdsArrayPath(), cDims);
    if(m_FeatureIdsPtr.lock())
    {
      m_FeatureIds = m_FeatureIdsPtr.lock()->getPointer(0);
    }
    if(getErrorCode() >= 0)
    {
      inputDataArrayPaths.push_back(getFeatureIdsArrayPath());
    }
  }

  EXECUTE_FUNCTION_TEMPLATE(this, createCompatibleArrays, m_InputArrayPtr.lock(), outputDataArrayPaths)

  if(m_FindHistogram)
  {
    std::vector<size_t> cDims_List = {static_cast<size_t>(m_NumBins)};
    DataArrayPath path(getDestinationAttributeMatrix().getDataContainerName(), getDestinationAttributeMatrix().getAttributeMatrixName(), getHistogramArrayName());
    m_HistogramListPtr = getDataContainerArray()->createNonPrereqArrayFromPath<FloatArrayType>(this, path, 0, cDims_List, "", DataArrayID37);
    if(getErrorCode() < 0)
    {
      return;
    }
  }

  if(getUseMask())
  {
    m_MaskPtr = getDataContainerArray()->getPrereqArrayFromPath<BoolArrayType>(this, getMaskArrayPath(), cDims);
    if(m_MaskPtr.lock())
    {
      m_Mask = m_MaskPtr.lock()->getPointer(0);
    }
    if(getErrorCode() >= 0)
    {
      inputDataArrayPaths.push_back(getMaskArrayPath());
    }
  }

  if(getStandardizeData())
  {
    if(!m_FindMean || !m_FindStdDeviation)
    {
      QString ss = QObject::tr("To standardize data, the \"Find Mean\" and \"Find Standard Deviation\" options must also be checked");
      setErrorCondition(-11003, ss);
    }
    DataArrayPath path(getSelectedArrayPath().getDataContainerName(), getSelectedArrayPath().getAttributeMatrixName(), getStandardizedArrayName());
    m_StandardizedPtr = getDataContainerArray()->createNonPrereqArrayFromPath<FloatArrayType>(this, path, 0, cDims, "", DataArrayID39);
    if(getErrorCode() >= 0)
    {
      inputDataArrayPaths.push_back(path);
    }
  }

  getDataContainerArray()->validateNumberOfTuples(this, inputDataArrayPaths);
  getDataContainerArray()->validateNumberOfTuples(this, outputDataArrayPaths);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
template <template <typename, typename...> class C, typename T, typename... Ts>
std::vector<float> findHistogram(C<T, Ts...>& source, float histmin, float histmax, bool histfullrange, int32_t numBins)
{
  if(source.empty())
  {
    return std::vector<float>(numBins, 0);
  }

  float min = 0.0f;
  float max = 0.0f;

  if(histfullrange)
  {
    min = static_cast<float>(StatisticsHelpers::findMin(source));
    max = static_cast<float>(StatisticsHelpers::findMax(source));
  }
  else
  {
    min = histmin;
    max = histmax;
  }

  float increment = (max - min) / (numBins);
  if(std::abs(increment) < 1E-10)
  {
    numBins = 1;
  }

  std::vector<float> Histogram(numBins, 0);

  if(numBins == 1) // if one bin, just set the first element to total number of points
  {
    Histogram[0] = static_cast<float>(source.size());
  }
  else
  {
    for(const auto s : source)
    {
      float value = static_cast<float>(s);
      size_t bin = static_cast<size_t>((value - min) / increment); // find bin for this input array value
      if((bin >= 0) && (bin < numBins))                            // make certain bin is in range
      {
        Histogram[bin]++; // increment histogram element corresponding to this input array value
      }
      else if(value == max)
      {
        Histogram[numBins - 1]++;
      }
    }
  }

  return Histogram;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
template <typename T>
void standardizeDataByIndex(IDataArray::Pointer dataPtr, FloatArrayType::Pointer standardized, bool useMask, bool* mask, int32_t* featureIds, IDataArray::Pointer meanPtr, IDataArray::Pointer stdPtr)
{
  typename DataArray<T>::Pointer inDataPtr = std::dynamic_pointer_cast<DataArray<T>>(dataPtr);
  FloatArrayType::Pointer muPtr = std::dynamic_pointer_cast<FloatArrayType>(meanPtr);
  FloatArrayType::Pointer sigPtr = std::dynamic_pointer_cast<FloatArrayType>(stdPtr);
  FloatArrayType::Pointer standardizedPtr = std::dynamic_pointer_cast<FloatArrayType>(standardized);

  T* dPtr = inDataPtr->getPointer(0);
  float* mPtr = muPtr->getPointer(0);
  float* sPtr = sigPtr->getPointer(0);
  float* stPtr = standardizedPtr->getPointer(0);
  size_t numTuples = inDataPtr->getNumberOfTuples();

  for(size_t i = 0; i < numTuples; i++)
  {
    if(useMask)
    {
      if(mask[i])
      {
        stPtr[i] = (static_cast<float>(dPtr[i]) - mPtr[featureIds[i]]) / sPtr[featureIds[i]];
      }
    }
    else
    {
      stPtr[i] = (static_cast<float>(dPtr[i]) - mPtr[featureIds[i]]) / sPtr[featureIds[i]];
    }
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
template <>
void standardizeDataByIndex<bool>(IDataArray::Pointer dataPtr, FloatArrayType::Pointer standardized, bool useMask, bool* mask, int32_t* featureIds, IDataArray::Pointer meanPtr,
                                  IDataArray::Pointer stdPtr)
{
  // Standardization of a boolean array is a no-op
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
template <typename T>
void standardizeData(IDataArray::Pointer dataPtr, FloatArrayType::Pointer standardized, bool useMask, bool* mask, IDataArray::Pointer meanPtr, IDataArray::Pointer stdPtr)
{
  typename DataArray<T>::Pointer inDataPtr = std::dynamic_pointer_cast<DataArray<T>>(dataPtr);
  FloatArrayType::Pointer muPtr = std::dynamic_pointer_cast<FloatArrayType>(meanPtr);
  FloatArrayType::Pointer sigPtr = std::dynamic_pointer_cast<FloatArrayType>(stdPtr);
  FloatArrayType::Pointer standardizedPtr = std::dynamic_pointer_cast<FloatArrayType>(standardized);

  T* dPtr = inDataPtr->getPointer(0);
  float* mPtr = muPtr->getPointer(0);
  float* sPtr = sigPtr->getPointer(0);
  float* stPtr = standardizedPtr->getPointer(0);
  size_t numTuples = inDataPtr->getNumberOfTuples();

  for(size_t i = 0; i < numTuples; i++)
  {
    if(useMask)
    {
      if(mask[i])
      {
        stPtr[i] = (static_cast<float>(dPtr[i]) - mPtr[0]) / sPtr[0];
      }
    }
    else
    {
      stPtr[i] = (static_cast<float>(dPtr[i]) - mPtr[0]) / sPtr[0];
    }
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
template <>
void standardizeData<bool>(IDataArray::Pointer dataPtr, FloatArrayType::Pointer standardized, bool useMask, bool* mask, IDataArray::Pointer meanPtr, IDataArray::Pointer stdPtr)
{
  // Standardization of a boolean array is a no-op
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
template <typename T>
class FindArrayStatisticsByIndexImpl
{
public:
  FindArrayStatisticsByIndexImpl(std::unordered_map<int32_t, std::list<T>>& featureDataMap, bool length, bool min, bool max, bool mean, bool median, bool stdDeviation, bool summation,
                                 std::vector<IDataArray::Pointer>& arrays, bool hist, float histmin, float histmax, bool histfullrange, int32_t numBins)
  : m_FeatureDataMap(featureDataMap)
  , m_Length(length)
  , m_Min(min)
  , m_Max(max)
  , m_Mean(mean)
  , m_Median(median)
  , m_StdDeviation(stdDeviation)
  , m_Summation(summation)
  , m_Histogram(hist)
  , m_HistMin(histmin)
  , m_HistMax(histmax)
  , m_HistFullRange(histfullrange)
  , m_NumBins(numBins)
  , m_Arrays(arrays)
  {
  }

  virtual ~FindArrayStatisticsByIndexImpl() = default;

  void compute(size_t start, size_t end) const
  {
    for(size_t i = start; i < end; i++)
    {
      if(m_FeatureDataMap.find(i) == m_FeatureDataMap.end())
      {
        continue;
      }
      if(m_Length)
      {
        if(m_Arrays[0])
        {
          int64_t val = static_cast<int64_t>(m_FeatureDataMap[i].size());
          m_Arrays[0]->initializeTuple(i, &val);
        }
      }
      if(m_Min)
      {
        if(m_Arrays[1])
        {
          T val = StatisticsHelpers::findMin(m_FeatureDataMap[i]);
          m_Arrays[1]->initializeTuple(i, &val);
        }
      }
      if(m_Max)
      {
        if(m_Arrays[2])
        {
          T val = StatisticsHelpers::findMax(m_FeatureDataMap[i]);
          m_Arrays[2]->initializeTuple(i, &val);
        }
      }
      if(m_Mean)
      {
        if(m_Arrays[3])
        {
          float val = StatisticsHelpers::findMean(m_FeatureDataMap[i]);
          m_Arrays[3]->initializeTuple(i, &val);
        }
      }
      if(m_Median)
      {
        if(m_Arrays[4])
        {
          float val = StatisticsHelpers::findMedian(m_FeatureDataMap[i]);
          m_Arrays[4]->initializeTuple(i, &val);
        }
      }
      if(m_StdDeviation)
      {
        if(m_Arrays[5])
        {
          float val = StatisticsHelpers::findStdDeviation(m_FeatureDataMap[i]);
          m_Arrays[5]->initializeTuple(i, &val);
        }
      }
      if(m_Summation)
      {
        if(m_Arrays[6])
        {
          float val = StatisticsHelpers::findSummation(m_FeatureDataMap[i]);
          m_Arrays[6]->initializeTuple(i, &val);
        }
      }

      if(m_Histogram)
      {
        if(m_Arrays[7])
        {
          std::vector<float> vals = findHistogram(m_FeatureDataMap[i], m_HistMin, m_HistMax, m_HistFullRange, m_NumBins);
          std::shared_ptr<FloatArrayType> histArray = std::dynamic_pointer_cast<FloatArrayType>(m_Arrays[7]);
          histArray->setTuple(i, vals);
        }
      }
    }
  }

#ifdef SIMPL_USE_PARALLEL_ALGORITHMS
  void operator()(const tbb::blocked_range<size_t>& r) const
  {
    compute(r.begin(), r.end());
  }
#endif

private:
  std::unordered_map<int32_t, std::list<T>>& m_FeatureDataMap;
  bool m_Length;
  bool m_Min;
  bool m_Max;
  bool m_Mean;
  bool m_Median;
  bool m_StdDeviation;
  bool m_Summation;
  bool m_Histogram;
  float m_HistMin;
  float m_HistMax;
  bool m_HistFullRange;
  int32_t m_NumBins;
  std::vector<IDataArray::Pointer>& m_Arrays;
};

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
template <typename T>
void findStatisticsImpl(bool length, bool min, bool max, bool mean, bool median, bool stdDeviation, bool summation, std::vector<IDataArray::Pointer>& arrays, std::vector<T>& data, bool hist,
                        float histmin, float histmax, bool histfullrange, int32_t numBins)
{
  if(length)
  {
    if(arrays[0])
    {
      int64_t val = static_cast<int64_t>(data.size());
      arrays[0]->initializeTuple(0, &val);
    }
  }
  if(min)
  {
    if(arrays[1])
    {
      T val = StatisticsHelpers::findMin(data);
      arrays[1]->initializeTuple(0, &val);
    }
  }
  if(max)
  {
    if(arrays[2])
    {
      T val = StatisticsHelpers::findMax(data);
      arrays[2]->initializeTuple(0, &val);
    }
  }
  if(mean)
  {
    if(arrays[3])
    {
      float val = StatisticsHelpers::findMean(data);
      arrays[3]->initializeTuple(0, &val);
    }
  }
  if(median)
  {
    if(arrays[4])
    {
      float val = StatisticsHelpers::findMedian(data);
      arrays[4]->initializeTuple(0, &val);
    }
  }
  if(stdDeviation)
  {
    if(arrays[5])
    {
      float val = StatisticsHelpers::findStdDeviation(data);
      arrays[5]->initializeTuple(0, &val);
    }
  }
  if(summation)
  {
    if(arrays[6])
    {
      float val = StatisticsHelpers::findSummation(data);
      arrays[6]->initializeTuple(0, &val);
    }
  }

  if(hist)
  {
    if(arrays[7])
    {
      std::vector<float> vals = findHistogram(data, histmin, histmax, histfullrange, numBins);
      std::shared_ptr<FloatArrayType> histArray = std::dynamic_pointer_cast<FloatArrayType>(arrays[7]);
      histArray->setTuple(0, vals);
    }
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
template <typename T>
void findStatistics(IDataArray::Pointer source, Int32ArrayType::Pointer featureIds, bool useMask, bool* mask, bool length, bool min, bool max, bool mean, bool median, bool stdDeviation,
                    bool summation, std::vector<IDataArray::Pointer>& arrays, int32_t numFeatures, bool computeByIndex, bool hist, float histmin, float histmax, bool histfullrange, int32_t numBins)
{
  size_t numTuples = source->getNumberOfTuples();
  typename DataArray<T>::Pointer sourcePtr = std::dynamic_pointer_cast<DataArray<T>>(source);
  T* dataPtr = sourcePtr->getPointer(0);

  if(computeByIndex)
  {
    int32_t* featureIdsPtr = featureIds->getPointer(0);
    std::unordered_map<int32_t, std::list<T>> featureValueMap;

    for(size_t i = 0; i < numTuples; i++)
    {
      if(useMask)
      {
        if(mask[i])
        {
          featureValueMap[featureIdsPtr[i]].push_back(dataPtr[i]);
        }
      }
      else
      {
        featureValueMap[featureIdsPtr[i]].push_back(dataPtr[i]);
      }
    }

#ifdef SIMPL_USE_PARALLEL_ALGORITHMS
    bool doParallel = true;
#endif

#ifdef SIMPL_USE_PARALLEL_ALGORITHMS
    if(doParallel)
    {
      tbb::parallel_for(tbb::blocked_range<size_t>(0, numFeatures),
                        FindArrayStatisticsByIndexImpl<T>(featureValueMap, length, min, max, mean, median, stdDeviation, summation, arrays, hist, histmin, histmax, histfullrange, numBins),
                        tbb::auto_partitioner());
    }
    else
#endif
    {
      FindArrayStatisticsByIndexImpl<T> serial(featureValueMap, length, min, max, mean, median, stdDeviation, summation, arrays, hist, histmin, histmax, histfullrange, numBins);
      serial.compute(0, numFeatures);
    }
  }
  else
  {
    std::vector<T> data;
    data.reserve(numTuples);
    for(size_t i = 0; i < numTuples; i++)
    {
      if(useMask)
      {
        if(mask[i])
        {
          data.push_back(dataPtr[i]);
        }
      }
      else
      {
        data.push_back(dataPtr[i]);
      }
    }

    data.shrink_to_fit();

    findStatisticsImpl(length, min, max, mean, median, stdDeviation, summation, arrays, data, hist, histmin, histmax, histfullrange, numBins);
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FindArrayStatistics::execute()
{
  dataCheck();
  if(getErrorCode() < 0)
  {
    return;
  }

  if(!m_FindHistogram && !m_FindMin && !m_FindMax && !m_FindMean && !m_FindMedian && !m_FindStdDeviation && !m_FindSummation && !m_FindLength)
  {
    return;
  }

  int32_t numFeatures = 0;
  DataContainer::Pointer dc = getDataContainerArray()->getDataContainer(m_DestinationAttributeMatrix.getDataContainerName());


  if(m_ComputeByIndex)
  {
    AttributeMatrix::Pointer attrMat = dc->getAttributeMatrix(m_DestinationAttributeMatrix.getAttributeMatrixName());
    numFeatures = static_cast<int32_t>(attrMat->getNumberOfTuples());
    bool mismatchedFeatures = false;
    int32_t largestFeature = 0;
    size_t totalPoints = m_FeatureIdsPtr.lock()->getNumberOfTuples();
    for(size_t i = 0; i < totalPoints; i++)
    {
      if(m_FeatureIds[i] > largestFeature)
      {
        largestFeature = m_FeatureIds[i];
        if(largestFeature >= numFeatures)
        {
          mismatchedFeatures = true;
          break;
        }
      }
    }

    if(mismatchedFeatures)
    {
      QString ss = QObject::tr("The given FeatureIds Array %1 has a value that is larger than allowed by the given Feature Attribute Matrix %2.\n %3 >= %4")
                       .arg(m_FeatureIdsArrayPath.serialize("/"))
                       .arg(m_DestinationAttributeMatrix.serialize("/"))
                       .arg(largestFeature)
                       .arg(numFeatures);
      setErrorCondition(-5555, ss);
      return;
    }
  }

  std::vector<IDataArray::Pointer> arrays(8, nullptr);

  for(size_t i = 0; i < arrays.size(); i++)
  {
    if(m_FindLength)
    {
      arrays[0] = m_LengthPtr.lock();
      arrays[0]->initializeWithZeros();
    }
    if(m_FindMin)
    {
      arrays[1] = m_MinimumPtr.lock();
      arrays[1]->initializeWithZeros();
    }
    if(m_FindMax)
    {
      arrays[2] = m_MaximumPtr.lock();
      arrays[2]->initializeWithZeros();
    }
    if(m_FindMean)
    {
      arrays[3] = m_MeanPtr.lock();
      arrays[3]->initializeWithZeros();
    }
    if(m_FindMedian)
    {
      arrays[4] = m_MedianPtr.lock();
      arrays[4]->initializeWithZeros();
    }
    if(m_FindStdDeviation)
    {
      arrays[5] = m_StandardDeviationPtr.lock();
      arrays[5]->initializeWithZeros();
    }
    if(m_FindSummation)
    {
      arrays[6] = m_SummationPtr.lock();
      arrays[6]->initializeWithZeros();
    }
    if(m_FindHistogram)
    {
      arrays[7] = m_HistogramListPtr.lock();
      arrays[7]->initializeWithZeros();
    }
  }

  EXECUTE_FUNCTION_TEMPLATE(this, findStatistics, m_InputArrayPtr.lock(), m_InputArrayPtr.lock(), m_FeatureIdsPtr.lock(), m_UseMask, m_Mask, m_FindLength, m_FindMin, m_FindMax, m_FindMean,
                            m_FindMedian, m_FindStdDeviation, m_FindSummation, arrays, numFeatures, m_ComputeByIndex, m_FindHistogram, m_MinRange, m_MaxRange, m_UseFullRange, m_NumBins)

  if(m_StandardizeData)
  {
    IDataArray::Pointer meanPtr = getDataContainerArray()->getAttributeMatrix(m_DestinationAttributeMatrix)->getAttributeArray(getMeanArrayName());
    IDataArray::Pointer stdPtr = getDataContainerArray()->getAttributeMatrix(m_DestinationAttributeMatrix)->getAttributeArray(getStdDeviationArrayName());

    if(m_ComputeByIndex)
    {
      EXECUTE_FUNCTION_TEMPLATE(this, standardizeDataByIndex, m_InputArrayPtr.lock(), m_InputArrayPtr.lock(), m_StandardizedPtr.lock(), m_UseMask, m_Mask, m_FeatureIds, meanPtr, stdPtr)
    }
    else
    {
      EXECUTE_FUNCTION_TEMPLATE(this, standardizeData, m_InputArrayPtr.lock(), m_InputArrayPtr.lock(), m_StandardizedPtr.lock(), m_UseMask, m_Mask, meanPtr, stdPtr)
    }
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AbstractFilter::Pointer FindArrayStatistics::newFilterInstance(bool copyFilterParameters) const
{
  FindArrayStatistics::Pointer filter = FindArrayStatistics::New();
  if(copyFilterParameters)
  {
    copyFilterParameterInstanceVariables(filter.get());
  }
  return filter;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString FindArrayStatistics::getCompiledLibraryName() const
{
  return DREAM3DReviewConstants::DREAM3DReviewBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString FindArrayStatistics::getBrandingString() const
{
  return "DREAM3DReview";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString FindArrayStatistics::getFilterVersion() const
{
  QString version;
  QTextStream vStream(&version);
  vStream << DREAM3DReview::Version::Major() << "." << DREAM3DReview::Version::Minor() << "." << DREAM3DReview::Version::Patch();
  return version;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString FindArrayStatistics::getGroupName() const
{
  return DREAM3DReviewConstants::FilterGroups::DREAM3DReviewFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QUuid FindArrayStatistics::getUuid() const
{
  return QUuid("{bf35f515-294b-55ed-8c69-211b7e69cb56}");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString FindArrayStatistics::getSubGroupName() const
{
  return DREAM3DReviewConstants::FilterSubGroups::StatisticsFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString FindArrayStatistics::getHumanLabel() const
{
  return "Find Attribute Array Statistics";
}

// -----------------------------------------------------------------------------
FindArrayStatistics::Pointer FindArrayStatistics::NullPointer()
{
  return Pointer(static_cast<Self*>(nullptr));
}

// -----------------------------------------------------------------------------
std::shared_ptr<FindArrayStatistics> FindArrayStatistics::New()
{
  struct make_shared_enabler : public FindArrayStatistics
  {
  };
  std::shared_ptr<make_shared_enabler> val = std::make_shared<make_shared_enabler>();
  val->setupFilterParameters();
  return val;
}

// -----------------------------------------------------------------------------
QString FindArrayStatistics::getNameOfClass() const
{
  return QString("FindArrayStatistics");
}

// -----------------------------------------------------------------------------
QString FindArrayStatistics::ClassName()
{
  return QString("FindArrayStatistics");
}
// -----------------------------------------------------------------------------
void FindArrayStatistics::setFindHistogram(bool value)
{
  m_FindHistogram = value;
}
// -----------------------------------------------------------------------------
void FindArrayStatistics::setNumBins(int32_t value)
{
  m_NumBins = value;
}

// -----------------------------------------------------------------------------
int32_t FindArrayStatistics::getNumBins() const
{
  return m_NumBins;
}
// -----------------------------------------------------------------------------
void FindArrayStatistics::setMinRange(double value)
{
  m_MinRange = value;
}

// -----------------------------------------------------------------------------
double FindArrayStatistics::getMinRange() const
{
  return m_MinRange;
}

// -----------------------------------------------------------------------------
void FindArrayStatistics::setMaxRange(double value)
{
  m_MaxRange = value;
}

// -----------------------------------------------------------------------------
double FindArrayStatistics::getMaxRange() const
{
  return m_MaxRange;
}
// -----------------------------------------------------------------------------
void FindArrayStatistics::setUseFullRange(bool value)
{
  m_UseFullRange = value;
}

// -----------------------------------------------------------------------------
bool FindArrayStatistics::getUseFullRange() const
{
  return m_UseFullRange;
}
// -----------------------------------------------------------------------------
bool FindArrayStatistics::getFindHistogram() const
{
  return m_FindHistogram;
}
// -----------------------------------------------------------------------------
void FindArrayStatistics::setFindLength(bool value)
{
  m_FindLength = value;
}

// -----------------------------------------------------------------------------
bool FindArrayStatistics::getFindLength() const
{
  return m_FindLength;
}

// -----------------------------------------------------------------------------
void FindArrayStatistics::setFindMin(bool value)
{
  m_FindMin = value;
}

// -----------------------------------------------------------------------------
bool FindArrayStatistics::getFindMin() const
{
  return m_FindMin;
}

// -----------------------------------------------------------------------------
void FindArrayStatistics::setFindMax(bool value)
{
  m_FindMax = value;
}

// -----------------------------------------------------------------------------
bool FindArrayStatistics::getFindMax() const
{
  return m_FindMax;
}

// -----------------------------------------------------------------------------
void FindArrayStatistics::setFindMean(bool value)
{
  m_FindMean = value;
}

// -----------------------------------------------------------------------------
bool FindArrayStatistics::getFindMean() const
{
  return m_FindMean;
}

// -----------------------------------------------------------------------------
void FindArrayStatistics::setFindMedian(bool value)
{
  m_FindMedian = value;
}

// -----------------------------------------------------------------------------
bool FindArrayStatistics::getFindMedian() const
{
  return m_FindMedian;
}

// -----------------------------------------------------------------------------
void FindArrayStatistics::setFindStdDeviation(bool value)
{
  m_FindStdDeviation = value;
}

// -----------------------------------------------------------------------------
bool FindArrayStatistics::getFindStdDeviation() const
{
  return m_FindStdDeviation;
}

// -----------------------------------------------------------------------------
void FindArrayStatistics::setFindSummation(bool value)
{
  m_FindSummation = value;
}

// -----------------------------------------------------------------------------
bool FindArrayStatistics::getFindSummation() const
{
  return m_FindSummation;
}

// -----------------------------------------------------------------------------
void FindArrayStatistics::setUseMask(bool value)
{
  m_UseMask = value;
}

// -----------------------------------------------------------------------------
bool FindArrayStatistics::getUseMask() const
{
  return m_UseMask;
}

// -----------------------------------------------------------------------------
void FindArrayStatistics::setStandardizeData(bool value)
{
  m_StandardizeData = value;
}

// -----------------------------------------------------------------------------
bool FindArrayStatistics::getStandardizeData() const
{
  return m_StandardizeData;
}

// -----------------------------------------------------------------------------
void FindArrayStatistics::setComputeByIndex(bool value)
{
  m_ComputeByIndex = value;
}

// -----------------------------------------------------------------------------
bool FindArrayStatistics::getComputeByIndex() const
{
  return m_ComputeByIndex;
}

// -----------------------------------------------------------------------------
void FindArrayStatistics::setDestinationAttributeMatrix(const DataArrayPath& value)
{
  m_DestinationAttributeMatrix = value;
}

// -----------------------------------------------------------------------------
DataArrayPath FindArrayStatistics::getDestinationAttributeMatrix() const
{
  return m_DestinationAttributeMatrix;
}

// -----------------------------------------------------------------------------
void FindArrayStatistics::setMaskArrayPath(const DataArrayPath& value)
{
  m_MaskArrayPath = value;
}

// -----------------------------------------------------------------------------
DataArrayPath FindArrayStatistics::getMaskArrayPath() const
{
  return m_MaskArrayPath;
}

// -----------------------------------------------------------------------------
void FindArrayStatistics::setHistogramArrayName(const QString& value)
{
  m_HistogramArrayName = value;
}

// -----------------------------------------------------------------------------
QString FindArrayStatistics::getHistogramArrayName() const
{
  return m_HistogramArrayName;
}

// -----------------------------------------------------------------------------
void FindArrayStatistics::setLengthArrayName(const QString& value)
{
  m_LengthArrayName = value;
}

// -----------------------------------------------------------------------------
QString FindArrayStatistics::getLengthArrayName() const
{
  return m_LengthArrayName;
}

// -----------------------------------------------------------------------------
void FindArrayStatistics::setMinimumArrayName(const QString& value)
{
  m_MinimumArrayName = value;
}

// -----------------------------------------------------------------------------
QString FindArrayStatistics::getMinimumArrayName() const
{
  return m_MinimumArrayName;
}

// -----------------------------------------------------------------------------
void FindArrayStatistics::setMaximumArrayName(const QString& value)
{
  m_MaximumArrayName = value;
}

// -----------------------------------------------------------------------------
QString FindArrayStatistics::getMaximumArrayName() const
{
  return m_MaximumArrayName;
}

// -----------------------------------------------------------------------------
void FindArrayStatistics::setMeanArrayName(const QString& value)
{
  m_MeanArrayName = value;
}

// -----------------------------------------------------------------------------
QString FindArrayStatistics::getMeanArrayName() const
{
  return m_MeanArrayName;
}

// -----------------------------------------------------------------------------
void FindArrayStatistics::setMedianArrayName(const QString& value)
{
  m_MedianArrayName = value;
}

// -----------------------------------------------------------------------------
QString FindArrayStatistics::getMedianArrayName() const
{
  return m_MedianArrayName;
}

// -----------------------------------------------------------------------------
void FindArrayStatistics::setStdDeviationArrayName(const QString& value)
{
  m_StdDeviationArrayName = value;
}

// -----------------------------------------------------------------------------
QString FindArrayStatistics::getStdDeviationArrayName() const
{
  return m_StdDeviationArrayName;
}

// -----------------------------------------------------------------------------
void FindArrayStatistics::setSummationArrayName(const QString& value)
{
  m_SummationArrayName = value;
}

// -----------------------------------------------------------------------------
QString FindArrayStatistics::getSummationArrayName() const
{
  return m_SummationArrayName;
}

// -----------------------------------------------------------------------------
void FindArrayStatistics::setStandardizedArrayName(const QString& value)
{
  m_StandardizedArrayName = value;
}

// -----------------------------------------------------------------------------
QString FindArrayStatistics::getStandardizedArrayName() const
{
  return m_StandardizedArrayName;
}

// -----------------------------------------------------------------------------
void FindArrayStatistics::setSelectedArrayPath(const DataArrayPath& value)
{
  m_SelectedArrayPath = value;
}

// -----------------------------------------------------------------------------
DataArrayPath FindArrayStatistics::getSelectedArrayPath() const
{
  return m_SelectedArrayPath;
}

// -----------------------------------------------------------------------------
void FindArrayStatistics::setFeatureIdsArrayPath(const DataArrayPath& value)
{
  m_FeatureIdsArrayPath = value;
}

// -----------------------------------------------------------------------------
DataArrayPath FindArrayStatistics::getFeatureIdsArrayPath() const
{
  return m_FeatureIdsArrayPath;
}
