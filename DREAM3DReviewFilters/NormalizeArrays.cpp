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

#include "NormalizeArrays.h"

#include <cassert>
#include <cmath>

#ifdef SIMPL_USE_PARALLEL_ALGORITHMS
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/partitioner.h>
#include <tbb/task_scheduler_init.h>
#endif

#include "SIMPLib/Common/Constants.h"
#include "SIMPLib/Common/TemplateHelpers.h"
#include "SIMPLib/FilterParameters/AbstractFilterParametersReader.h"
#include "SIMPLib/FilterParameters/DataArraySelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/DoubleFilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedBooleanFilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedChoicesFilterParameter.h"
#include "SIMPLib/FilterParameters/MultiDataArraySelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/StringFilterParameter.h"

#include "DREAM3DReview/DREAM3DReviewConstants.h"
#include "DREAM3DReview/DREAM3DReviewVersion.h"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
NormalizeArrays::NormalizeArrays()
: m_SelectedDataArrayPaths(QVector<DataArrayPath>())
, m_NormalizeType(0)
, m_RangeMin(0.0)
, m_RangeMax(1.0)
, m_Postfix("_Normalized")
, m_UseMask(false)
, m_MaskArrayPath("", "", "")
, m_DefaultValue(0.0)
, m_Mask(nullptr)
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
NormalizeArrays::~NormalizeArrays()
= default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void NormalizeArrays::setupFilterParameters()
{
  FilterParameterVector parameters;
  {
    LinkedChoicesFilterParameter::Pointer parameter = LinkedChoicesFilterParameter::New();
    parameter->setHumanLabel("Operation Type");
    parameter->setPropertyName("NormalizeType");
    parameter->setSetterCallback(SIMPL_BIND_SETTER(NormalizeArrays, this, NormalizeType));
    parameter->setGetterCallback(SIMPL_BIND_GETTER(NormalizeArrays, this, NormalizeType));
    QVector<QString> choices;
    choices.push_back("Rescale to Range");
    choices.push_back("Standardize");
    parameter->setChoices(choices);
    QStringList linkedProps = {"RangeMin", "RangeMax"};
    parameter->setLinkedProperties(linkedProps);
    parameter->setEditable(false);
    parameter->setCategory(FilterParameter::Parameter);
    parameters.push_back(parameter);
  }
  parameters.push_back(SIMPL_NEW_STRING_FP("Postfix", Postfix, FilterParameter::Parameter, NormalizeArrays));
  QStringList linkedProps = {"MaskArrayPath", "DefaultValue"};
  parameters.push_back(SIMPL_NEW_LINKED_BOOL_FP("Use Mask", UseMask, FilterParameter::Parameter, NormalizeArrays, linkedProps));
  parameters.push_back(SIMPL_NEW_DOUBLE_FP("Default Masked Value", DefaultValue, FilterParameter::Parameter, NormalizeArrays));
  parameters.push_back(SIMPL_NEW_DOUBLE_FP("Range Minimum", RangeMin, FilterParameter::Parameter, NormalizeArrays, 0));
  parameters.push_back(SIMPL_NEW_DOUBLE_FP("Range Maximum", RangeMax, FilterParameter::Parameter, NormalizeArrays, 0));
  MultiDataArraySelectionFilterParameter::RequirementType req =
      MultiDataArraySelectionFilterParameter::CreateRequirement(SIMPL::Defaults::AnyPrimitive, 1, AttributeMatrix::Type::Any, IGeometry::Type::Any);
  parameters.push_back(SIMPL_NEW_MDA_SELECTION_FP("Attribute Arrays to Normalize", SelectedDataArrayPaths, FilterParameter::RequiredArray, NormalizeArrays, req));
  DataArraySelectionFilterParameter::RequirementType dasReq = DataArraySelectionFilterParameter::CreateRequirement(SIMPL::TypeNames::Bool, 1, AttributeMatrix::Type::Any, IGeometry::Type::Any);
  parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Mask", MaskArrayPath, FilterParameter::RequiredArray, NormalizeArrays, dasReq));
  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void NormalizeArrays::readFilterParameters(AbstractFilterParametersReader* reader, int index)
{
  reader->openFilterGroup(this, index);
  setSelectedDataArrayPaths(reader->readDataArrayPathVector("SelectedDataArrayPaths", getSelectedDataArrayPaths()));
  setNormalizeType(reader->readValue("NormalizeType", getNormalizeType()));
  setRangeMin(reader->readValue("RangeMin", getRangeMin()));
  setRangeMax(reader->readValue("RangeMax", getRangeMax()));
  setPostfix(reader->readString("Postfix", getPostfix()));
  setUseMask(reader->readValue("UseMask", getUseMask()));
  setMaskArrayPath(reader->readDataArrayPath("MaskArrayPath", getMaskArrayPath()));
  setDefaultValue(reader->readValue("DefaultValue", getDefaultValue()));
  reader->closeFilterGroup();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void NormalizeArrays::initialize()
{
  m_SelectedWeakPtrVector.clear();
  m_NormalizedArraysPtrVector.clear();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void NormalizeArrays::dataCheck()
{
  setErrorCondition(0);
  setWarningCondition(0);
  initialize();

  if(getNormalizeType() != 0 && getNormalizeType() != 1)
  {
    QString ss = QObject::tr("Invalid selection for operation type");
    setErrorCondition(-701);
    notifyErrorMessage(getHumanLabel(), ss, getErrorCondition());
    return;
  }

  if(getSelectedDataArrayPaths().empty())
  {
    setErrorCondition(-11001);
    QString ss = QObject::tr("At least one Attribute Array must be selected");
    notifyErrorMessage(getHumanLabel(), ss, getErrorCondition());
    return;
  }

  if(getPostfix().isEmpty())
  {
    setErrorCondition(-11001);
    QString ss = QObject::tr("A postfix for the normalized Attribute Arrays must be entered");
    notifyErrorMessage(getHumanLabel(), ss, getErrorCondition());
  }

  QVector<DataArrayPath> paths = getSelectedDataArrayPaths();

  if(!DataArrayPath::ValidateVector(paths))
  {
    setErrorCondition(-11002);
    QString ss = QObject::tr("There are Attribute Arrays selected that are not contained in the same Attribute Matrix; all selected Attribute Arrays must belong to the same Attribute Matrix");
    notifyErrorMessage(getHumanLabel(), ss, getErrorCondition());
  }

  if(getErrorCondition() < 0)
  {
    return;
  }

  QVector<size_t> cDims(1, 1);
  QString dcName = paths[0].getDataContainerName();
  QString amName = paths[0].getAttributeMatrixName();

  for(auto&& path : paths)
  {
    IDataArray::WeakPointer ptr = getDataContainerArray()->getPrereqIDataArrayFromPath<IDataArray, AbstractFilter>(this, path);
    if(ptr.lock())
    {
      m_SelectedWeakPtrVector.push_back(ptr);
      int32_t numComps = ptr.lock()->getNumberOfComponents();
      if(numComps != 1)
      {
        setErrorCondition(-11003);
        QString ss = QObject::tr("All Attribute Arrays must be scalar arrays, but %1 has %2 total components").arg(ptr.lock()->getName()).arg(numComps);
        notifyErrorMessage(getHumanLabel(), ss, getErrorCondition());
      }
      else
      {
        QString arrayName = path.getDataArrayName() + getPostfix();
        DataArrayPath tempPath(dcName, amName, arrayName);
        double defaultValue = m_UseMask ? getDefaultValue() : 0.0;
        DoubleArrayType::WeakPointer ptr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<double>, AbstractFilter>(this, tempPath, defaultValue, cDims);
        if(getErrorCondition() >= 0)
        {
          m_NormalizedArraysPtrVector.push_back(ptr.lock());
        }
      }
    }
  }

  if(getUseMask())
  {
    m_MaskPtr = getDataContainerArray()->getPrereqArrayFromPath<DataArray<bool>, AbstractFilter>(this, getMaskArrayPath(), cDims);
    if(m_MaskPtr.lock())
    {
      m_Mask = m_MaskPtr.lock()->getPointer(0);
    }
    if(getErrorCondition() >= 0)
    {
      paths.push_back(getMaskArrayPath());
    }
  }

  getDataContainerArray()->validateNumberOfTuples<AbstractFilter>(this, paths);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void NormalizeArrays::preflight()
{
  setInPreflight(true);
  emit preflightAboutToExecute();
  emit updateFilterParameters(this);
  dataCheck();
  emit preflightExecuted();
  setInPreflight(false);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
template <typename T> void copyDataArrays(IDataArray::Pointer dataPtr, std::vector<double>& copy, bool useMask, bool* mask)
{
  typename DataArray<T>::Pointer inDataPtr = std::dynamic_pointer_cast<DataArray<T>>(dataPtr);
  T* dPtr = inDataPtr->getPointer(0);
  size_t numTuples = inDataPtr->getNumberOfTuples();

  for(size_t i = 0; i < numTuples; i++)
  {
    if(useMask)
    {
      if(mask[i])
      {
        copy.push_back(dPtr[i]);
      }
    }
    else
    {
      copy.push_back(dPtr[i]);
    }
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
class NormalizeArraysImpl
{
public:
  NormalizeArraysImpl(std::vector<std::vector<double>>& arrays, int32_t normalizeType, double rangeMin, double rangeMax)
  : m_Arrays(arrays)
  , m_NormalizeType(normalizeType)
  , m_RangeMin(rangeMin)
  , m_RangeMax(rangeMax)
  {
  }

  virtual ~NormalizeArraysImpl()
  = default;

  void compute(size_t start, size_t end) const
  {
    for(size_t i = start; i < end; i++)
    {
      if(m_NormalizeType == 0)
      {
        double min = *std::min_element(std::begin(m_Arrays[i]), std::end(m_Arrays[i]));
        double max = *std::max_element(std::begin(m_Arrays[i]), std::end(m_Arrays[i]));
        std::for_each(std::begin(m_Arrays[i]), std::end(m_Arrays[i]), [&](double& val) { val = m_RangeMin + (((val - min) * (m_RangeMax - m_RangeMin)) / (max - min)); });
      }
      else if(m_NormalizeType == 1)
      {
        double sum = std::accumulate(std::begin(m_Arrays[i]), std::end(m_Arrays[i]), 0.0);
        double mean = sum / m_Arrays[i].size();
        std::vector<double> difference(m_Arrays[i].size());
        std::transform(std::begin(m_Arrays[i]), std::end(m_Arrays[i]), std::begin(difference), [&](double val) { return val - mean; });
        double squaredSum = std::inner_product(std::begin(difference), std::end(difference), std::begin(difference), 0.0f);
        double stdDev = std::sqrt(squaredSum / m_Arrays[i].size());
        std::for_each(std::begin(m_Arrays[i]), std::end(m_Arrays[i]), [&](double& val) { val = (val - mean) / stdDev; });
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
  std::vector<std::vector<double>>& m_Arrays;
  int32_t m_NormalizeType;
  double m_RangeMin;
  double m_RangeMax;
};

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void NormalizeArrays::execute()
{
  setErrorCondition(0);
  setWarningCondition(0);
  dataCheck();
  if(getErrorCondition() < 0)
  {
    return;
  }

  if(m_SelectedDataArrayPaths.size() != m_SelectedWeakPtrVector.size())
  {
    QString ss = QObject::tr("The number of selected Attribute Arrays does not equal the number of internal weak pointers");
    setErrorCondition(-11008);
    notifyErrorMessage(getHumanLabel(), ss, getErrorCondition());
    return;
  }

  size_t numTuples = m_SelectedWeakPtrVector[0].lock()->getNumberOfTuples();

  std::vector<std::vector<double>> arrays(m_SelectedWeakPtrVector.size());
  for(auto&& vec : arrays)
  {
    vec.reserve(numTuples);
  }

  for(std::vector<double>::size_type i = 0; i < arrays.size(); i++)
  {
    EXECUTE_FUNCTION_TEMPLATE(this, copyDataArrays, m_SelectedWeakPtrVector[i].lock(), m_SelectedWeakPtrVector[i].lock(), arrays[i], m_UseMask, m_Mask);
  }

#ifdef SIMPL_USE_PARALLEL_ALGORITHMS
  tbb::task_scheduler_init init;
  bool doParallel = true;
#endif

#ifdef SIMPL_USE_PARALLEL_ALGORITHMS
  if(doParallel)
  {
    tbb::parallel_for(tbb::blocked_range<size_t>(0, arrays.size()), NormalizeArraysImpl(arrays, m_NormalizeType, m_RangeMin, m_RangeMax), tbb::auto_partitioner());
  }
  else
#endif
  {
    NormalizeArraysImpl serial(arrays, m_NormalizeType, m_RangeMin, m_RangeMax);
    serial.compute(0, arrays.size());
  }

  size_t index = 0;

  assert(m_NormalizedArraysPtrVector.size() == arrays.size());

  for(size_t i = 0; i < numTuples; i++)
  {
    for(std::vector<double>::size_type j = 0; j < arrays.size(); j++)
    {
      if(m_UseMask)
      {
        if(m_Mask[i])
        {
          m_NormalizedArraysPtrVector[j]->setValue(i, arrays[j][index]);
          index++;
        }
      }
      else
      {
        m_NormalizedArraysPtrVector[j]->setValue(i, arrays[j][i]);
      }
    }
  }

  notifyStatusMessage(getHumanLabel(), "Complete");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AbstractFilter::Pointer NormalizeArrays::newFilterInstance(bool copyFilterParameters) const
{
  NormalizeArrays::Pointer filter = NormalizeArrays::New();
  if(copyFilterParameters)
  {
    copyFilterParameterInstanceVariables(filter.get());
  }
  return filter;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString NormalizeArrays::getCompiledLibraryName() const
{
  return DREAM3DReviewConstants::DREAM3DReviewBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString NormalizeArrays::getBrandingString() const
{
  return "DREAM3DReview";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString NormalizeArrays::getFilterVersion() const
{
  QString version;
  QTextStream vStream(&version);
  vStream << DREAM3DReview::Version::Major() << "." << DREAM3DReview::Version::Minor() << "." << DREAM3DReview::Version::Patch();
  return version;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString NormalizeArrays::getGroupName() const
{
  return DREAM3DReviewConstants::FilterGroups::DREAM3DReviewFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QUuid NormalizeArrays::getUuid()
{
  return QUuid("{8c584519-15c3-5010-b5ed-a2ac626591a1}");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString NormalizeArrays::getSubGroupName() const
{
  return DREAM3DReviewConstants::FilterSubGroups::StatisticsFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString NormalizeArrays::getHumanLabel() const
{
  return "Normalize Attribute Arrays";
}
