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
#include "FindDynamicArrayStatistics.h"

#include <functional>
#include <numeric>

#include "SIMPLib/Common/Constants.h"
#include "SIMPLib/Common/TemplateHelpers.h"
#include "SIMPLib/DataArrays/NeighborList.hpp"
#include "SIMPLib/FilterParameters/AbstractFilterParametersReader.h"
#include "SIMPLib/FilterParameters/DataArrayCreationFilterParameter.h"
#include "SIMPLib/FilterParameters/DataArraySelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedBooleanFilterParameter.h"
#include "SIMPLib/Math/SIMPLibMath.h"

#include "DREAM3DReview/DREAM3DReviewConstants.h"
#include "DREAM3DReview/DREAM3DReviewVersion.h"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
FindDynamicArrayStatistics::FindDynamicArrayStatistics()
: m_FindLength(false)
, m_FindMin(false)
, m_FindMax(false)
, m_FindMean(false)
, m_FindMedian(false)
, m_FindStdDeviation(false)
, m_FindSummation(false)
, m_LengthArrayPath("", "", "Length")
, m_MinimumArrayPath("", "", "Minimum")
, m_MaximumArrayPath("", "", "Maximum")
, m_MeanArrayPath("", "", "Mean")
, m_MedianArrayPath("", "", "Median")
, m_StdDeviationArrayPath("", "", "StandardDeviation")
, m_SummationArrayPath("", "", "Summation")
, m_SelectedArrayPath("", "", "")
, m_Length(nullptr)
, m_Minimum(nullptr)
, m_Maximum(nullptr)
, m_Mean(nullptr)
, m_Median(nullptr)
, m_StandardDeviation(nullptr)
, m_Summation(nullptr)
, m_InputArray(nullptr)
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
FindDynamicArrayStatistics::~FindDynamicArrayStatistics() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FindDynamicArrayStatistics::setupFilterParameters()
{
  FilterParameterVectorType parameters;
  QStringList linkedProps("LengthArrayPath");
  parameters.push_back(SIMPL_NEW_LINKED_BOOL_FP("Find Length", FindLength, FilterParameter::Parameter, FindDynamicArrayStatistics, linkedProps));
  linkedProps.clear();
  linkedProps << "MinimumArrayPath";
  parameters.push_back(SIMPL_NEW_LINKED_BOOL_FP("Find Minimum", FindMin, FilterParameter::Parameter, FindDynamicArrayStatistics, linkedProps));
  linkedProps.clear();
  linkedProps << "MaximumArrayPath";
  parameters.push_back(SIMPL_NEW_LINKED_BOOL_FP("Find Maximum", FindMax, FilterParameter::Parameter, FindDynamicArrayStatistics, linkedProps));
  linkedProps.clear();
  linkedProps << "MeanArrayPath";
  parameters.push_back(SIMPL_NEW_LINKED_BOOL_FP("Find Mean", FindMean, FilterParameter::Parameter, FindDynamicArrayStatistics, linkedProps));
  linkedProps.clear();
  linkedProps << "MedianArrayPath";
  parameters.push_back(SIMPL_NEW_LINKED_BOOL_FP("Find Median", FindMedian, FilterParameter::Parameter, FindDynamicArrayStatistics, linkedProps));
  linkedProps.clear();
  linkedProps << "StdDeviationArrayPath";
  parameters.push_back(SIMPL_NEW_LINKED_BOOL_FP("Find Standard Deviation", FindStdDeviation, FilterParameter::Parameter, FindDynamicArrayStatistics, linkedProps));
  linkedProps.clear();
  linkedProps << "SummationArrayPath";
  parameters.push_back(SIMPL_NEW_LINKED_BOOL_FP("Find Summation", FindSummation, FilterParameter::Parameter, FindDynamicArrayStatistics, linkedProps));
  DataArrayCreationFilterParameter::RequirementType dacReq = DataArrayCreationFilterParameter::CreateRequirement(AttributeMatrix::Type::Any, IGeometry::Type::Any);

  parameters.push_back(SIMPL_NEW_DA_CREATION_FP("Length", LengthArrayPath, FilterParameter::CreatedArray, FindDynamicArrayStatistics, dacReq));
  parameters.push_back(SIMPL_NEW_DA_CREATION_FP("Minimum", MinimumArrayPath, FilterParameter::CreatedArray, FindDynamicArrayStatistics, dacReq));
  parameters.push_back(SIMPL_NEW_DA_CREATION_FP("Maximum", MaximumArrayPath, FilterParameter::CreatedArray, FindDynamicArrayStatistics, dacReq));
  parameters.push_back(SIMPL_NEW_DA_CREATION_FP("Mean", MeanArrayPath, FilterParameter::CreatedArray, FindDynamicArrayStatistics, dacReq));
  parameters.push_back(SIMPL_NEW_DA_CREATION_FP("Median", MedianArrayPath, FilterParameter::CreatedArray, FindDynamicArrayStatistics, dacReq));
  parameters.push_back(SIMPL_NEW_DA_CREATION_FP("Standard Deviation", StdDeviationArrayPath, FilterParameter::CreatedArray, FindDynamicArrayStatistics, dacReq));
  parameters.push_back(SIMPL_NEW_DA_CREATION_FP("Summation", SummationArrayPath, FilterParameter::CreatedArray, FindDynamicArrayStatistics, dacReq));
  DataArraySelectionFilterParameter::RequirementType dasReq = DataArraySelectionFilterParameter::CreateRequirement(SIMPL::Defaults::AnyPrimitive, 1, AttributeMatrix::Type::Any, IGeometry::Type::Any);
  parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Attribute Array to Compute Statistics", SelectedArrayPath, FilterParameter::RequiredArray, FindDynamicArrayStatistics, dasReq));
  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FindDynamicArrayStatistics::readFilterParameters(AbstractFilterParametersReader* reader, int index)
{
  reader->openFilterGroup(this, index);
  setFindLength(reader->readValue("FindLength", getFindLength()));
  setFindMin(reader->readValue("FindMin", getFindMin()));
  setFindMax(reader->readValue("FindMax", getFindMax()));
  setFindMean(reader->readValue("FindMean", getFindMean()));
  setFindMedian(reader->readValue("FindMedian", getFindMedian()));
  setFindStdDeviation(reader->readValue("FindStdDeviation", getFindStdDeviation()));
  setFindSummation(reader->readValue("FindSummation", getFindSummation()));
  setLengthArrayPath(reader->readDataArrayPath("LengthArrayPath", getLengthArrayPath()));
  setMinimumArrayPath(reader->readDataArrayPath("MinimumArrayPath", getMinimumArrayPath()));
  setMaximumArrayPath(reader->readDataArrayPath("MaximumArrayPath", getMaximumArrayPath()));
  setMeanArrayPath(reader->readDataArrayPath("MeanArrayPath", getMeanArrayPath()));
  setMedianArrayPath(reader->readDataArrayPath("MedianArrayPath", getMedianArrayPath()));
  setStdDeviationArrayPath(reader->readDataArrayPath("StdDeviationArrayPath", getStdDeviationArrayPath()));
  setSummationArrayPath(reader->readDataArrayPath("SummationArrayPath", getSummationArrayPath()));
  setSelectedArrayPath(reader->readDataArrayPath("SelectedArrayPath", getSelectedArrayPath()));
  reader->closeFilterGroup();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FindDynamicArrayStatistics::initialize()
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FindDynamicArrayStatistics::dataCheck()
{
  clearErrorCode();
  clearWarningCode();

  QVector<DataArrayPath> dataArrayPaths;

  m_InputArrayPtr = getDataContainerArray()->getPrereqIDataArrayFromPath<IDataArray, AbstractFilter>(this, getSelectedArrayPath());

  if(getErrorCode() < 0)
  {
    return;
  }

  dataArrayPaths.push_back(getSelectedArrayPath());

  if(m_InputArrayPtr.lock()->getNumberOfComponents() != 1)
  {
    QString ss = QObject::tr("Input Attribute Array must be a scalar array");
    setErrorCondition(-11002, ss);
  }

  EXECUTE_FUNCTION_TEMPLATE_NEIGHBORLIST(this, createCompatibleArrays, m_InputArrayPtr.lock(), dataArrayPaths)

  QVector<size_t> cDims(1, 1);

  if(m_FindLength)
  {
    m_LengthPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<int32_t>, AbstractFilter, int32_t>(this, getLengthArrayPath(), 0,
                                                                                                                     cDims); /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
    if(nullptr != m_LengthPtr.lock().get()) /* Validate the Weak Pointer wraps a non-nullptr pointer to a DataArray<T> object */
    {
      m_Length = m_LengthPtr.lock()->getPointer(0);
    } /* Now assign the raw pointer to data from the DataArray<T> object */
    if(getErrorCode() >= 0)
    {
      dataArrayPaths.push_back(getLengthArrayPath());
    }
  }

  getDataContainerArray()->validateNumberOfTuples<AbstractFilter>(this, dataArrayPaths);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FindDynamicArrayStatistics::preflight()
{
  // These are the REQUIRED lines of CODE to make sure the filter behaves correctly
  setInPreflight(true);              // Set the fact that we are preflighting.
  emit preflightAboutToExecute();    // Emit this signal so that other widgets can do one file update
  emit updateFilterParameters(this); // Emit this signal to have the widgets push their values down to the filter
  dataCheck();                       // Run our DataCheck to make sure everthing is setup correctly
  emit preflightExecuted();          // We are done preflighting this filter
  setInPreflight(false);             // Inform the system this filter is NOT in preflight mode anymore.
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
template <typename T> void findLength(IDataArray::Pointer source, IDataArray::Pointer length)
{
  typename NeighborList<T>::Pointer inputDataPtr = std::dynamic_pointer_cast<NeighborList<T>>(source);
  typename DataArray<int32_t>::Pointer lengthDataPtr = std::dynamic_pointer_cast<DataArray<int32_t>>(length);
  int32_t* lengthPtr = lengthDataPtr->getPointer(0);

  size_t numTuples = inputDataPtr->getNumberOfTuples();

  for(size_t i = 0; i < numTuples; i++)
  {
    lengthPtr[i] = inputDataPtr->getListSize(i);
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
template <typename T> void findMin(IDataArray::Pointer source, IDataArray::Pointer min)
{
  typename NeighborList<T>::Pointer inputDataPtr = std::dynamic_pointer_cast<NeighborList<T>>(source);
  typename DataArray<T>::Pointer minDataPtr = std::dynamic_pointer_cast<DataArray<T>>(min);
  T* minPtr = minDataPtr->getPointer(0);

  size_t numTuples = inputDataPtr->getNumberOfTuples();

  for(size_t i = 0; i < numTuples; i++)
  {
    std::vector<T> tmpList = inputDataPtr->getListReference(i);
    if(tmpList.size() == 0)
    {
      continue;
    }
    minPtr[i] = *std::min_element(tmpList.begin(), tmpList.end());
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
template <typename T> void findMax(IDataArray::Pointer source, IDataArray::Pointer max)
{
  typename NeighborList<T>::Pointer inputDataPtr = std::dynamic_pointer_cast<NeighborList<T>>(source);
  typename DataArray<T>::Pointer maxDataPtr = std::dynamic_pointer_cast<DataArray<T>>(max);
  T* maxPtr = maxDataPtr->getPointer(0);

  size_t numTuples = inputDataPtr->getNumberOfTuples();

  for(size_t i = 0; i < numTuples; i++)
  {
    std::vector<T> tmpList = inputDataPtr->getListReference(i);
    if(tmpList.size() == 0)
    {
      continue;
    }
    maxPtr[i] = *std::max_element(tmpList.begin(), tmpList.end());
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
template <typename T> void findMean(IDataArray::Pointer source, IDataArray::Pointer mean)
{
  typename NeighborList<T>::Pointer inputDataPtr = std::dynamic_pointer_cast<NeighborList<T>>(source);
  typename DataArray<T>::Pointer meanDataPtr = std::dynamic_pointer_cast<DataArray<T>>(mean);
  T* meanPtr = meanDataPtr->getPointer(0);

  size_t numTuples = inputDataPtr->getNumberOfTuples();

  for(size_t i = 0; i < numTuples; i++)
  {
    std::vector<T> tmpList = inputDataPtr->getListReference(i);
    if(tmpList.size() == 0)
    {
      continue;
    }
    T sum = std::accumulate(tmpList.begin(), tmpList.end(), T(0));
    meanPtr[i] = sum / tmpList.size();
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
template <typename T> void findMedian(IDataArray::Pointer source, IDataArray::Pointer median)
{
  typename NeighborList<T>::Pointer inputDataPtr = std::dynamic_pointer_cast<NeighborList<T>>(source);
  typename DataArray<T>::Pointer medianDataPtr = std::dynamic_pointer_cast<DataArray<T>>(median);
  T* medianPtr = medianDataPtr->getPointer(0);

  size_t numTuples = inputDataPtr->getNumberOfTuples();

  for(size_t i = 0; i < numTuples; i++)
  {
    // Need a copy, not a reference, since we will be messing with the vector order
    std::vector<T> tmpList = inputDataPtr->copyOfList(i);
    if(tmpList.size() == 0)
    {
      continue;
    }
    auto halfElements = tmpList.size() / 2;
    std::nth_element(tmpList.begin(), tmpList.begin() + halfElements, tmpList.end());
    T medVal = tmpList[halfElements];
    if(tmpList.size() % 2 == 1)
    {
      medianPtr[i] = medVal;
    }
    else
    {
      medianPtr[i] = T(0.5) * (medVal + tmpList[halfElements - 1]);
    }
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
template <typename T> void findStdDeviation(IDataArray::Pointer source, IDataArray::Pointer stdDev)
{
  typename NeighborList<T>::Pointer inputDataPtr = std::dynamic_pointer_cast<NeighborList<T>>(source);
  typename DataArray<T>::Pointer stdDevDataPtr = std::dynamic_pointer_cast<DataArray<T>>(stdDev);
  T* stdDevPtr = stdDevDataPtr->getPointer(0);

  size_t numTuples = inputDataPtr->getNumberOfTuples();

  for(size_t i = 0; i < numTuples; i++)
  {
    std::vector<T> tmpList = inputDataPtr->getListReference(i);
    if(tmpList.size() == 0)
    {
      continue;
    }
    std::vector<T> difference(tmpList.size());
    T sum = std::accumulate(tmpList.begin(), tmpList.end(), T(0));
    T mean = sum / tmpList.size();
    std::transform(tmpList.begin(), tmpList.end(), difference.begin(), std::bind2nd(std::minus<T>(), mean));
    T squaredSum = std::inner_product(difference.begin(), difference.end(), difference.begin(), T(0));
    stdDevPtr[i] = sqrt(squaredSum / tmpList.size());
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
template <typename T> void findSummation(IDataArray::Pointer source, IDataArray::Pointer summation)
{
  typename NeighborList<T>::Pointer inputDataPtr = std::dynamic_pointer_cast<NeighborList<T>>(source);
  typename DataArray<T>::Pointer summationDataPtr = std::dynamic_pointer_cast<DataArray<T>>(summation);
  T* summationPtr = summationDataPtr->getPointer(0);

  size_t numTuples = inputDataPtr->getNumberOfTuples();

  for(size_t i = 0; i < numTuples; i++)
  {
    std::vector<T> tmpList = inputDataPtr->getListReference(i);
    if(tmpList.size() == 0)
    {
      continue;
    }
    T sum = std::accumulate(tmpList.begin(), tmpList.end(), T(0));
    summationPtr[i] = sum;
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FindDynamicArrayStatistics::execute()
{
  clearErrorCode();
  clearWarningCode();
  dataCheck();
  if(getErrorCode() < 0)
  {
    return;
  }

  if(m_FindLength)
  {
    EXECUTE_FUNCTION_TEMPLATE_NEIGHBORLIST(this, findLength, m_InputArrayPtr.lock(), m_InputArrayPtr.lock(), m_LengthPtr.lock())
  }
  if(m_FindMin)
  {
    EXECUTE_FUNCTION_TEMPLATE_NEIGHBORLIST(this, findMin, m_InputArrayPtr.lock(), m_InputArrayPtr.lock(), m_MinimumPtr.lock())
  }
  if(m_FindMax)
  {
    EXECUTE_FUNCTION_TEMPLATE_NEIGHBORLIST(this, findMax, m_InputArrayPtr.lock(), m_InputArrayPtr.lock(), m_MaximumPtr.lock())
  }
  if(m_FindMean)
  {
    EXECUTE_FUNCTION_TEMPLATE_NEIGHBORLIST(this, findMean, m_InputArrayPtr.lock(), m_InputArrayPtr.lock(), m_MeanPtr.lock())
  }
  if(m_FindMedian)
  {
    EXECUTE_FUNCTION_TEMPLATE_NEIGHBORLIST(this, findMedian, m_InputArrayPtr.lock(), m_InputArrayPtr.lock(), m_MedianPtr.lock())
  }
  if(m_FindStdDeviation)
  {
    EXECUTE_FUNCTION_TEMPLATE_NEIGHBORLIST(this, findStdDeviation, m_InputArrayPtr.lock(), m_InputArrayPtr.lock(), m_StandardDeviationPtr.lock())
  }
  if(m_FindSummation)
  {
    EXECUTE_FUNCTION_TEMPLATE_NEIGHBORLIST(this, findSummation, m_InputArrayPtr.lock(), m_InputArrayPtr.lock(), m_SummationPtr.lock())
  }

  notifyStatusMessage("Complete");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AbstractFilter::Pointer FindDynamicArrayStatistics::newFilterInstance(bool copyFilterParameters) const
{
  FindDynamicArrayStatistics::Pointer filter = FindDynamicArrayStatistics::New();
  if(copyFilterParameters)
  {
    copyFilterParameterInstanceVariables(filter.get());
  }
  return filter;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString FindDynamicArrayStatistics::getCompiledLibraryName() const
{
  return DREAM3DReviewConstants::DREAM3DReviewBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString FindDynamicArrayStatistics::getBrandingString() const
{
  return "DREAM3DReview";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString FindDynamicArrayStatistics::getFilterVersion() const
{
  QString version;
  QTextStream vStream(&version);
  vStream << DREAM3DReview::Version::Major() << "." << DREAM3DReview::Version::Minor() << "." << DREAM3DReview::Version::Patch();
  return version;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString FindDynamicArrayStatistics::getGroupName() const
{
  return SIMPL::FilterGroups::StatisticsFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString FindDynamicArrayStatistics::getSubGroupName() const
{
  return SIMPL::FilterSubGroups::MiscFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString FindDynamicArrayStatistics::getHumanLabel() const
{
  return "Find Dynamic Array Statistics";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QUuid FindDynamicArrayStatistics::getUuid()
{
  return QUuid("{73ee33b6-7622-5004-8b88-4d145514fb6a}");
}