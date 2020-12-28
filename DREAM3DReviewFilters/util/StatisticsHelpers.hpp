

#include <vector>

#ifndef STATISTICS_FILTER_CLASS_NAME
#include "DREAM3DReviewFilters/FindArrayStatistics.h"
#define STATISTICS_FILTER_CLASS_NAME FindArrayStatistics

#endif

#ifndef STATISTICS_RENAME_DATAPATH_ENUM
/* Create Enumerations to allow the created Attribute Arrays to take part in renaming */
enum createdPathID : RenameDataPath::DataID_t
{
  DataArrayID30 = 30, // Length
  DataArrayID31 = 31, // Min
  DataArrayID32 = 32, // Max
  DataArrayID33 = 33, // Mean
  DataArrayID34 = 34, // Median
  DataArrayID35 = 35, // StdDev
  DataArrayID36 = 36, // Summation
  DataArrayID37 = 37, // Histogram
  DataArrayID38 = 38, //
  DataArrayID39 = 39, // StandardizedArray
  DataArrayID40 = 40, //
};
#define STATISTICS_RENAME_DATAPATH_ENUM 1
#endif

/**
 * @brief createCompatibleArrays Creates the output statistics arrays with compatible types based on the
 * user options and incoming DataArray type
 */
template <typename T>
void STATISTICS_FILTER_CLASS_NAME::createCompatibleArrays(QVector<DataArrayPath>& dataArrayPaths)
{
  std::vector<size_t> cDims = {1};
  using DataArrayType = DataArray<T>;

  if(m_FindLength)
  {
    DataArrayPath path(getDestinationAttributeMatrix().getDataContainerName(), getDestinationAttributeMatrix().getAttributeMatrixName(), getLengthArrayName());
    m_LengthPtr = getDataContainerArray()->createNonPrereqArrayFromPath<MeshIndexArrayType>(this, path, 0, cDims, "", DataArrayID30);
    if(getErrorCode() >= 0)
    {
      dataArrayPaths.push_back(path);
    }
  }

  if(m_FindMin)
  {
    DataArrayPath path(getDestinationAttributeMatrix().getDataContainerName(), getDestinationAttributeMatrix().getAttributeMatrixName(), getMinimumArrayName());
    m_MinimumPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArrayType>(this, path, 0, cDims, "", DataArrayID31);
    if(getErrorCode() >= 0)
    {
      dataArrayPaths.push_back(path);
    }
  }

  if(m_FindMax)
  {
    DataArrayPath path(getDestinationAttributeMatrix().getDataContainerName(), getDestinationAttributeMatrix().getAttributeMatrixName(), getMaximumArrayName());
    m_MaximumPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArrayType>(this, path, 0, cDims, "", DataArrayID32);
    if(getErrorCode() >= 0)
    {
      dataArrayPaths.push_back(path);
    }
  }

  if(m_FindMean)
  {
    DataArrayPath path(getDestinationAttributeMatrix().getDataContainerName(), getDestinationAttributeMatrix().getAttributeMatrixName(), getMeanArrayName());
    m_MeanPtr = getDataContainerArray()->createNonPrereqArrayFromPath<FloatArrayType>(this, path, 0, cDims, "", DataArrayID33);
    if(getErrorCode() >= 0)
    {
      dataArrayPaths.push_back(path);
    }
  }

  if(m_FindMedian)
  {
    DataArrayPath path(getDestinationAttributeMatrix().getDataContainerName(), getDestinationAttributeMatrix().getAttributeMatrixName(), getMedianArrayName());
    m_MedianPtr = getDataContainerArray()->createNonPrereqArrayFromPath<FloatArrayType>(this, path, 0, cDims, "", DataArrayID34);
    if(getErrorCode() >= 0)
    {
      dataArrayPaths.push_back(path);
    }
  }

  if(m_FindStdDeviation)
  {
    DataArrayPath path(getDestinationAttributeMatrix().getDataContainerName(), getDestinationAttributeMatrix().getAttributeMatrixName(), getStdDeviationArrayName());
    m_StandardDeviationPtr = getDataContainerArray()->createNonPrereqArrayFromPath<FloatArrayType>(this, path, 0, cDims, "", DataArrayID35);
    if(getErrorCode() >= 0)
    {
      dataArrayPaths.push_back(path);
    }
  }

  if(m_FindSummation)
  {
    DataArrayPath path(getDestinationAttributeMatrix().getDataContainerName(), getDestinationAttributeMatrix().getAttributeMatrixName(), getSummationArrayName());
    m_SummationPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<float>>(this, path, 0, cDims, "", DataArrayID36);
    if(getErrorCode() >= 0)
    {
      dataArrayPaths.push_back(path);
    }
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
template <template <typename, typename...> class C, typename T, typename... Ts>
int64_t findLength(C<T, Ts...>& source)
{
  return static_cast<int64_t>(source.size());
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
template <template <typename, typename...> class C, typename T, typename... Ts>
T findMin(C<T, Ts...>& source)
{
  if(source.empty())
  {
    return T(0);
  }
  return (*std::min_element(std::begin(source), std::end(source)));
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
template <template <typename, typename...> class C, typename T, typename... Ts>
T findMax(C<T, Ts...>& source)
{
  if(source.empty())
  {
    return T(0);
  }
  return (*std::max_element(std::begin(source), std::end(source)));
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
template <template <typename, typename...> class C, typename T, typename... Ts>
float findMean(C<T, Ts...>& source)
{
  if(source.empty())
  {
    return 0.0f;
  }
  float sum = std::accumulate(std::begin(source), std::end(source), 0.0f);
  return static_cast<float>(sum / source.size());
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
template <template <typename, typename...> class C, typename... Ts>
bool findMean(C<bool, Ts...>& source)
{
  if(source.empty())
  {
    return false;
  }
  size_t count = std::count(std::begin(source), std::end(source), true);
  return true ? count >= (source.size() - count) : false;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
template <template <typename, typename...> class C, typename T, typename... Ts>
float findMedian(C<T, Ts...>& source)
{
  // Need a copy, not a reference, since we will be messing with the vector order
  std::vector<T> tmpList{std::begin(source), std::end(source)};
  if(tmpList.empty())
  {
    return 0.0f;
  }
  auto halfElements = tmpList.size() / 2;
  std::nth_element(std::begin(tmpList), tmpList.begin() + halfElements, std::end(tmpList));
  T medVal = tmpList[halfElements];
  if(tmpList.size() % 2 == 1)
  {
    return static_cast<float>(medVal);
  }

  return static_cast<float>((0.5f * (medVal + tmpList[halfElements - 1])));
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
template <template <typename, typename...> class C, typename T, typename... Ts>
float findStdDeviation(C<T, Ts...>& source)
{
  if(source.empty())
  {
    return 0.0f;
  }
  std::vector<double> difference(source.size());
  float sum = std::accumulate(std::begin(source), std::end(source), 0.0f);
  float mean = static_cast<double>(sum / source.size());
  std::transform(std::begin(source), std::end(source), std::begin(difference), [mean](float a) { return a - mean; });
  float squaredSum = std::inner_product(std::begin(difference), std::end(difference), std::begin(difference), 0.0f);
  return std::sqrt(squaredSum / source.size());
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
template <template <typename, typename...> class C, typename... Ts>
bool findStdDeviation(C<bool, Ts...>& source)
{
  if(source.empty())
  {
    return false;
  }
  size_t count = std::count(std::begin(source), std::end(source), true);
  return true ? count >= (source.size() - count) : false;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
template <template <typename, typename...> class C, typename T, typename... Ts>
double findSummation(C<T, Ts...>& source)
{
  if(source.empty())
  {
    return 0.0f;
  }
  float sum = std::accumulate(std::begin(source), std::end(source), 0.0f);
  return sum;
}
