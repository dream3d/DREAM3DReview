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

#include "FindNeighborListStatistics.h"

#include <functional>

#include <QtCore/QTextStream>

#include "SIMPLib/Common/Constants.h"
#include "SIMPLib/Common/TemplateHelpers.h"
#include "SIMPLib/DataArrays/NeighborList.hpp"
#include "SIMPLib/FilterParameters/AttributeMatrixSelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/DataArraySelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedBooleanFilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedPathCreationFilterParameter.h"
#include "SIMPLib/Utilities/ParallelDataAlgorithm.h"

#include "DREAM3DReview/DREAM3DReviewConstants.h"
#include "DREAM3DReview/DREAM3DReviewVersion.h"

#define STATISTICS_FILTER_CLASS_NAME FindNeighborListStatistics
#include "util/StatisticsHelpers.hpp"

// -----------------------------------------------------------------------------
FindNeighborListStatistics::FindNeighborListStatistics() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
FindNeighborListStatistics::~FindNeighborListStatistics() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FindNeighborListStatistics::setupFilterParameters()
{
  FilterParameterVectorType parameters;
  std::vector<QString> linkedProps;

  linkedProps.clear();
  linkedProps.push_back("LengthArrayName");
  parameters.push_back(SIMPL_NEW_LINKED_BOOL_FP("Find Length", FindLength, FilterParameter::Category::Parameter, FindNeighborListStatistics, linkedProps));
  linkedProps.clear();
  linkedProps.push_back("MinimumArrayName");
  parameters.push_back(SIMPL_NEW_LINKED_BOOL_FP("Find Minimum", FindMin, FilterParameter::Category::Parameter, FindNeighborListStatistics, linkedProps));
  linkedProps.clear();
  linkedProps.push_back("MaximumArrayName");
  parameters.push_back(SIMPL_NEW_LINKED_BOOL_FP("Find Maximum", FindMax, FilterParameter::Category::Parameter, FindNeighborListStatistics, linkedProps));
  linkedProps.clear();
  linkedProps.push_back("MeanArrayName");
  parameters.push_back(SIMPL_NEW_LINKED_BOOL_FP("Find Mean", FindMean, FilterParameter::Category::Parameter, FindNeighborListStatistics, linkedProps));
  linkedProps.clear();
  linkedProps.push_back("MedianArrayName");
  parameters.push_back(SIMPL_NEW_LINKED_BOOL_FP("Find Median", FindMedian, FilterParameter::Category::Parameter, FindNeighborListStatistics, linkedProps));
  linkedProps.clear();
  linkedProps.push_back("StdDeviationArrayName");
  parameters.push_back(SIMPL_NEW_LINKED_BOOL_FP("Find Standard Deviation", FindStdDeviation, FilterParameter::Category::Parameter, FindNeighborListStatistics, linkedProps));
  linkedProps.clear();
  linkedProps.push_back("SummationArrayName");
  parameters.push_back(SIMPL_NEW_LINKED_BOOL_FP("Find Summation", FindSummation, FilterParameter::Category::Parameter, FindNeighborListStatistics, linkedProps));
  linkedProps.clear();

  DataArraySelectionFilterParameter::RequirementType dasReq = DataArraySelectionFilterParameter::CreateRequirement(SIMPL::Defaults::AnyPrimitive, 1, AttributeMatrix::Type::Any, IGeometry::Type::Any);
  parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Attribute Array to Compute Statistics", SelectedArrayPath, FilterParameter::Category::RequiredArray, FindNeighborListStatistics, dasReq));

  AttributeMatrixSelectionFilterParameter::RequirementType amReq = AttributeMatrixSelectionFilterParameter::CreateRequirement(AttributeMatrix::Type::Any, IGeometry::Type::Any);
  parameters.push_back(SIMPL_NEW_AM_SELECTION_FP("Destination Attribute Matrix", DestinationAttributeMatrix, FilterParameter::Category::RequiredArray, FindNeighborListStatistics, amReq));

  parameters.push_back(
      SIMPL_NEW_DA_WITH_LINKED_AM_FP("Length", LengthArrayName, DestinationAttributeMatrix, DestinationAttributeMatrix, FilterParameter::Category::CreatedArray, FindNeighborListStatistics));
  parameters.push_back(
      SIMPL_NEW_DA_WITH_LINKED_AM_FP("Minimum", MinimumArrayName, DestinationAttributeMatrix, DestinationAttributeMatrix, FilterParameter::Category::CreatedArray, FindNeighborListStatistics));
  parameters.push_back(
      SIMPL_NEW_DA_WITH_LINKED_AM_FP("Maximum", MaximumArrayName, DestinationAttributeMatrix, DestinationAttributeMatrix, FilterParameter::Category::CreatedArray, FindNeighborListStatistics));
  parameters.push_back(
      SIMPL_NEW_DA_WITH_LINKED_AM_FP("Mean", MeanArrayName, DestinationAttributeMatrix, DestinationAttributeMatrix, FilterParameter::Category::CreatedArray, FindNeighborListStatistics));
  parameters.push_back(
      SIMPL_NEW_DA_WITH_LINKED_AM_FP("Median", MedianArrayName, DestinationAttributeMatrix, DestinationAttributeMatrix, FilterParameter::Category::CreatedArray, FindNeighborListStatistics));
  parameters.push_back(SIMPL_NEW_DA_WITH_LINKED_AM_FP("Standard Deviation", StdDeviationArrayName, DestinationAttributeMatrix, DestinationAttributeMatrix, FilterParameter::Category::CreatedArray,
                                                      FindNeighborListStatistics));
  parameters.push_back(
      SIMPL_NEW_DA_WITH_LINKED_AM_FP("Summation", SummationArrayName, DestinationAttributeMatrix, DestinationAttributeMatrix, FilterParameter::Category::CreatedArray, FindNeighborListStatistics));

  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FindNeighborListStatistics::initialize()
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FindNeighborListStatistics::dataCheck()
{

  clearErrorCode();
  clearWarningCode();

  if(!getFindMin() && !getFindMax() && !getFindMean() && !getFindMedian() && !getFindStdDeviation() && !getFindSummation() && !getFindLength())
  {
    QString ss = QObject::tr("No statistics have been selected, so this filter will perform no operations");
    setWarningCondition(-11005, ss);
    return;
  }

  QVector<DataArrayPath> dataArrayPaths;

  m_InputArrayPtr = getDataContainerArray()->getPrereqIDataArrayFromPath(this, getSelectedArrayPath());
  if(getErrorCode() < 0)
  {
    return;
  }
  // The input array must be a NeighborList.
  if(m_InputArrayPtr.lock()->getNameOfClass() != "NeighborList<T>")
  {
    QString ss = QObject::tr("Input Data must be a NeighborList Attribute Array");
    setErrorCondition(-11006, ss);
  }
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
  if(getErrorCode() < 0)
  {
    return;
  }
  EXECUTE_FUNCTION_TEMPLATE_NO_BOOL(NeighborList, this, createCompatibleArrays, m_InputArrayPtr.lock(), dataArrayPaths)

  getDataContainerArray()->validateNumberOfTuples(this, dataArrayPaths);
}

// -----------------------------------------------------------------------------
template <typename T>
class FindNeighborListStatisticsImpl
{
public:
  FindNeighborListStatisticsImpl(AbstractFilter* filter, IDataArray::Pointer& source, bool length, bool min, bool max, bool mean, bool median, bool stdDeviation, bool summation,
                                 std::vector<IDataArray::Pointer>& arrays)
  : m_Filter(filter)
  , m_Source(source)
  , m_Length(length)
  , m_Min(min)
  , m_Max(max)
  , m_Mean(mean)
  , m_Median(median)
  , m_StdDeviation(stdDeviation)
  , m_Summation(summation)
  , m_Arrays(arrays)
  {
  }

  virtual ~FindNeighborListStatisticsImpl() = default;

  void compute(size_t start, size_t end) const
  {
    using NeighborListType = NeighborList<T>;
    typename NeighborListType::Pointer inputDataPtr = std::dynamic_pointer_cast<NeighborListType>(m_Source);

    for(size_t i = start; i < end; i++)
    {
      if(m_Filter->getCancel())
      {
        break;
      }
      std::vector<T>& tmpList = (*inputDataPtr)[i];

      if(m_Length)
      {
        if(m_Arrays[0])
        {
          int64_t val = static_cast<int64_t>(tmpList.size());
          m_Arrays[0]->initializeTuple(i, &val);
        }
      }
      if(m_Min)
      {
        if(m_Arrays[1])
        {
          T val = StatisticsHelpers::findMin(tmpList);
          m_Arrays[1]->initializeTuple(i, &val);
        }
      }
      if(m_Max)
      {
        if(m_Arrays[2])
        {
          T val = StatisticsHelpers::findMax(tmpList);
          m_Arrays[2]->initializeTuple(i, &val);
        }
      }
      if(m_Mean)
      {
        if(m_Arrays[3])
        {
          float val = StatisticsHelpers::findMean(tmpList);
          m_Arrays[3]->initializeTuple(i, &val);
        }
      }
      if(m_Median)
      {
        if(m_Arrays[4])
        {
          float val = StatisticsHelpers::findMedian(tmpList);
          m_Arrays[4]->initializeTuple(i, &val);
        }
      }
      if(m_StdDeviation)
      {
        if(m_Arrays[5])
        {
          float val = StatisticsHelpers::findStdDeviation(tmpList);
          m_Arrays[5]->initializeTuple(i, &val);
        }
      }
      if(m_Summation)
      {
        if(m_Arrays[6])
        {
          float val = StatisticsHelpers::findSummation(tmpList);
          m_Arrays[6]->initializeTuple(i, &val);
        }
      }
    }
  }

  void operator()(const SIMPLRange& range) const
  {
    compute(range.min(), range.max());
  }

private:
  AbstractFilter* m_Filter = nullptr;
  IDataArray::Pointer m_Source;
  bool m_Length = false;
  bool m_Min = false;
  bool m_Max = false;
  bool m_Mean = false;
  bool m_Median = false;
  bool m_StdDeviation = false;
  bool m_Summation = false;

  std::vector<IDataArray::Pointer>& m_Arrays;
};

// -----------------------------------------------------------------------------
template <typename T>
void findStatistics(AbstractFilter* filter, IDataArray::Pointer source, bool length, bool min, bool max, bool mean, bool median, bool stdDeviation, bool summation,
                    std::vector<IDataArray::Pointer>& arrays)
{
  size_t numTuples = source->getNumberOfTuples();
  // Allow data-based parallelization
  ParallelDataAlgorithm dataAlg;
  dataAlg.setRange(0, numTuples);
  dataAlg.execute(FindNeighborListStatisticsImpl<T>(filter, source, length, min, max, mean, median, stdDeviation, summation, arrays));
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FindNeighborListStatistics::execute()
{
  dataCheck();
  if(getErrorCode() < 0)
  {
    return;
  }

  if(!m_FindMin && !m_FindMax && !m_FindMean && !m_FindMedian && !m_FindStdDeviation && !m_FindSummation && !m_FindLength)
  {
    return;
  }

  std::vector<IDataArray::Pointer> arrays(7, nullptr);

  if(m_FindLength)
  {
    arrays[0] = m_LengthPtr.lock();
  }
  if(m_FindMin)
  {
    arrays[1] = m_MinimumPtr.lock();
  }
  if(m_FindMax)
  {
    arrays[2] = m_MaximumPtr.lock();
  }
  if(m_FindMean)
  {
    arrays[3] = m_MeanPtr.lock();
  }
  if(m_FindMedian)
  {
    arrays[4] = m_MedianPtr.lock();
  }
  if(m_FindStdDeviation)
  {
    arrays[5] = m_StandardDeviationPtr.lock();
  }
  if(m_FindSummation)
  {
    arrays[6] = m_SummationPtr.lock();
  }

  EXECUTE_FUNCTION_TEMPLATE_NO_BOOL(NeighborList, this, findStatistics, m_InputArrayPtr.lock(), this, m_InputArrayPtr.lock(), m_FindLength, m_FindMin, m_FindMax, m_FindMean, m_FindMedian,
                                    m_FindStdDeviation, m_FindSummation, arrays)
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AbstractFilter::Pointer FindNeighborListStatistics::newFilterInstance(bool copyFilterParameters) const
{
  FindNeighborListStatistics::Pointer filter = FindNeighborListStatistics::New();
  if(copyFilterParameters)
  {
    copyFilterParameterInstanceVariables(filter.get());
  }
  return filter;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString FindNeighborListStatistics::getCompiledLibraryName() const
{
  return DREAM3DReviewConstants::DREAM3DReviewBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString FindNeighborListStatistics::getBrandingString() const
{
  return "DREAM3DReview";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString FindNeighborListStatistics::getFilterVersion() const
{
  QString version;
  QTextStream vStream(&version);
  vStream << DREAM3DReview::Version::Major() << "." << DREAM3DReview::Version::Minor() << "." << DREAM3DReview::Version::Patch();
  return version;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString FindNeighborListStatistics::getGroupName() const
{
  return DREAM3DReviewConstants::FilterGroups::DREAM3DReviewFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString FindNeighborListStatistics::getSubGroupName() const
{
  return DREAM3DReviewConstants::FilterSubGroups::StatisticsFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString FindNeighborListStatistics::getHumanLabel() const
{
  return "Find NeighborList Statistics";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QUuid FindNeighborListStatistics::getUuid() const
{
  return QUuid("{73ee33b6-7622-5004-8b88-4d145514fb6a}");
}
// -----------------------------------------------------------------------------
FindNeighborListStatistics::Pointer FindNeighborListStatistics::NullPointer()
{
  return Pointer(static_cast<Self*>(nullptr));
}

// -----------------------------------------------------------------------------
std::shared_ptr<FindNeighborListStatistics> FindNeighborListStatistics::New()
{
  struct make_shared_enabler : public FindNeighborListStatistics
  {
  };
  std::shared_ptr<make_shared_enabler> val = std::make_shared<make_shared_enabler>();
  val->setupFilterParameters();
  return val;
}

// -----------------------------------------------------------------------------
QString FindNeighborListStatistics::getNameOfClass() const
{
  return QString("FindNeighborListStatistics");
}

// -----------------------------------------------------------------------------
QString FindNeighborListStatistics::ClassName()
{
  return QString("FindNeighborListStatistics");
}

// -----------------------------------------------------------------------------
void FindNeighborListStatistics::setFindLength(bool value)
{
  m_FindLength = value;
}

// -----------------------------------------------------------------------------
bool FindNeighborListStatistics::getFindLength() const
{
  return m_FindLength;
}

// -----------------------------------------------------------------------------
void FindNeighborListStatistics::setFindMin(bool value)
{
  m_FindMin = value;
}

// -----------------------------------------------------------------------------
bool FindNeighborListStatistics::getFindMin() const
{
  return m_FindMin;
}

// -----------------------------------------------------------------------------
void FindNeighborListStatistics::setFindMax(bool value)
{
  m_FindMax = value;
}

// -----------------------------------------------------------------------------
bool FindNeighborListStatistics::getFindMax() const
{
  return m_FindMax;
}

// -----------------------------------------------------------------------------
void FindNeighborListStatistics::setFindMean(bool value)
{
  m_FindMean = value;
}

// -----------------------------------------------------------------------------
bool FindNeighborListStatistics::getFindMean() const
{
  return m_FindMean;
}

// -----------------------------------------------------------------------------
void FindNeighborListStatistics::setFindMedian(bool value)
{
  m_FindMedian = value;
}

// -----------------------------------------------------------------------------
bool FindNeighborListStatistics::getFindMedian() const
{
  return m_FindMedian;
}

// -----------------------------------------------------------------------------
void FindNeighborListStatistics::setFindStdDeviation(bool value)
{
  m_FindStdDeviation = value;
}

// -----------------------------------------------------------------------------
bool FindNeighborListStatistics::getFindStdDeviation() const
{
  return m_FindStdDeviation;
}

// -----------------------------------------------------------------------------
void FindNeighborListStatistics::setFindSummation(bool value)
{
  m_FindSummation = value;
}

// -----------------------------------------------------------------------------
bool FindNeighborListStatistics::getFindSummation() const
{
  return m_FindSummation;
}

// -----------------------------------------------------------------------------
void FindNeighborListStatistics::setDestinationAttributeMatrix(const DataArrayPath& value)
{
  m_DestinationAttributeMatrix = value;
}

// -----------------------------------------------------------------------------
DataArrayPath FindNeighborListStatistics::getDestinationAttributeMatrix() const
{
  return m_DestinationAttributeMatrix;
}

// -----------------------------------------------------------------------------
void FindNeighborListStatistics::setLengthArrayName(const QString& value)
{
  m_LengthArrayName = value;
}

// -----------------------------------------------------------------------------
QString FindNeighborListStatistics::getLengthArrayName() const
{
  return m_LengthArrayName;
}

// -----------------------------------------------------------------------------
void FindNeighborListStatistics::setMinimumArrayName(const QString& value)
{
  m_MinimumArrayName = value;
}

// -----------------------------------------------------------------------------
QString FindNeighborListStatistics::getMinimumArrayName() const
{
  return m_MinimumArrayName;
}

// -----------------------------------------------------------------------------
void FindNeighborListStatistics::setMaximumArrayName(const QString& value)
{
  m_MaximumArrayName = value;
}

// -----------------------------------------------------------------------------
QString FindNeighborListStatistics::getMaximumArrayName() const
{
  return m_MaximumArrayName;
}

// -----------------------------------------------------------------------------
void FindNeighborListStatistics::setMeanArrayName(const QString& value)
{
  m_MeanArrayName = value;
}

// -----------------------------------------------------------------------------
QString FindNeighborListStatistics::getMeanArrayName() const
{
  return m_MeanArrayName;
}

// -----------------------------------------------------------------------------
void FindNeighborListStatistics::setMedianArrayName(const QString& value)
{
  m_MedianArrayName = value;
}

// -----------------------------------------------------------------------------
QString FindNeighborListStatistics::getMedianArrayName() const
{
  return m_MedianArrayName;
}

// -----------------------------------------------------------------------------
void FindNeighborListStatistics::setStdDeviationArrayName(const QString& value)
{
  m_StdDeviationArrayName = value;
}

// -----------------------------------------------------------------------------
QString FindNeighborListStatistics::getStdDeviationArrayName() const
{
  return m_StdDeviationArrayName;
}

// -----------------------------------------------------------------------------
void FindNeighborListStatistics::setSummationArrayName(const QString& value)
{
  m_SummationArrayName = value;
}

// -----------------------------------------------------------------------------
QString FindNeighborListStatistics::getSummationArrayName() const
{
  return m_SummationArrayName;
}

// -----------------------------------------------------------------------------
void FindNeighborListStatistics::setSelectedArrayPath(const DataArrayPath& value)
{
  m_SelectedArrayPath = value;
}

// -----------------------------------------------------------------------------
DataArrayPath FindNeighborListStatistics::getSelectedArrayPath() const
{
  return m_SelectedArrayPath;
}
