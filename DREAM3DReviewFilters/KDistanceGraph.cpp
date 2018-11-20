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

#include "KDistanceGraph.h"

#include "SIMPLib/Common/TemplateHelpers.h"
#include "SIMPLib/FilterParameters/AbstractFilterParametersReader.h"
#include "SIMPLib/FilterParameters/ChoiceFilterParameter.h"
#include "SIMPLib/FilterParameters/DataArrayCreationFilterParameter.h"
#include "SIMPLib/FilterParameters/DataArraySelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/IntFilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedBooleanFilterParameter.h"

#include "util/EvaluationAlgorithms/KDistanceTemplate.hpp"

#include "DREAM3DReview/DREAM3DReviewConstants.h"
#include "DREAM3DReview/DREAM3DReviewVersion.h"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
KDistanceGraph::KDistanceGraph()
: m_SelectedArrayPath("", "", "")
, m_UseMask(false)
, m_MaskArrayPath("", "", "")
, m_KDistanceArrayPath("", "", "KDistance")
, m_MinDist(1)
, m_DistanceMetric(0)
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
KDistanceGraph::~KDistanceGraph() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void KDistanceGraph::setupFilterParameters()
{
  FilterParameterVector parameters;
  parameters.push_back(SIMPL_NEW_INTEGER_FP("K<sup>th</sup> Nearest Neighbor", MinDist, FilterParameter::Parameter, KDistanceGraph));
  {
    ChoiceFilterParameter::Pointer parameter = ChoiceFilterParameter::New();
    parameter->setHumanLabel("Distance Metric");
    parameter->setPropertyName("DistanceMetric");
    parameter->setSetterCallback(SIMPL_BIND_SETTER(KDistanceGraph, this, DistanceMetric));
    parameter->setGetterCallback(SIMPL_BIND_GETTER(KDistanceGraph, this, DistanceMetric));
    QVector<QString> choices = DistanceTemplate::GetDistanceMetricsOptions();
    parameter->setChoices(choices);
    parameter->setCategory(FilterParameter::Parameter);
    parameters.push_back(parameter);
  }
  QStringList linkedProps("MaskArrayPath");
  parameters.push_back(SIMPL_NEW_LINKED_BOOL_FP("Use Mask", UseMask, FilterParameter::Parameter, KDistanceGraph, linkedProps));
  DataArraySelectionFilterParameter::RequirementType dasReq =
      DataArraySelectionFilterParameter::CreateRequirement(SIMPL::Defaults::AnyPrimitive, SIMPL::Defaults::AnyComponentSize, AttributeMatrix::Type::Any, IGeometry::Type::Any);
  parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Input Attribute Array", SelectedArrayPath, FilterParameter::RequiredArray, KDistanceGraph, dasReq));
  dasReq = DataArraySelectionFilterParameter::CreateRequirement(SIMPL::TypeNames::Bool, 1, AttributeMatrix::Type::Any, IGeometry::Type::Any);
  DataArrayCreationFilterParameter::RequirementType dacReq = DataArrayCreationFilterParameter::CreateRequirement(AttributeMatrix::Category::Unknown);
  parameters.push_back(SIMPL_NEW_DA_CREATION_FP("K Distance", KDistanceArrayPath, FilterParameter::CreatedArray, KDistanceGraph, dacReq));
  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void KDistanceGraph::readFilterParameters(AbstractFilterParametersReader* reader, int index)
{
  reader->openFilterGroup(this, index);
  setSelectedArrayPath(reader->readDataArrayPath("SelectedArrayPath", getSelectedArrayPath()));
  setUseMask(reader->readValue("UseMask", getUseMask()));
  setMaskArrayPath(reader->readDataArrayPath("MaskArrayPath", getMaskArrayPath()));
  setMinDist(reader->readValue("MinDist", getMinDist()));
  setKDistanceArrayPath(reader->readDataArrayPath("KDistanceArrayPath", getKDistanceArrayPath()));
  setDistanceMetric(reader->readValue("DistanceMetric", getDistanceMetric()));
  reader->closeFilterGroup();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void KDistanceGraph::initialize()
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void KDistanceGraph::dataCheck()
{
  setErrorCondition(0);
  setWarningCondition(0);
  initialize();

  if(getMinDist() < 1)
  {
    setErrorCondition(-5555);
    notifyErrorMessage(getHumanLabel(), "Kth nearest neighbor must be greater than 0", getErrorCondition());
    return;
  }

  QVector<size_t> cDims(1, 1);
  QVector<DataArrayPath> dataArrayPaths;

  m_InDataPtr = getDataContainerArray()->getPrereqIDataArrayFromPath<IDataArray, AbstractFilter>(this, getSelectedArrayPath());
  if(getErrorCondition() >= 0)
  {
    dataArrayPaths.push_back(getSelectedArrayPath());
  }

  m_KDistanceArrayPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<double>, AbstractFilter, double>(this, getKDistanceArrayPath(), 0, cDims);
  if(m_KDistanceArrayPtr.lock())
  {
    m_KDistanceArray = m_KDistanceArrayPtr.lock()->getPointer(0);
  }
  if(getErrorCondition() >= 0)
  {
    dataArrayPaths.push_back(getKDistanceArrayPath());
  }

  if(getUseMask())
  {
    m_MaskPtr = getDataContainerArray()->getPrereqArrayFromPath<BoolArrayType, AbstractFilter>(this, getMaskArrayPath(), cDims);
    if(m_MaskPtr.lock())
    {
      m_Mask = m_MaskPtr.lock()->getPointer(0);
    }
    if(getErrorCondition() >= 0)
    {
      dataArrayPaths.push_back(getMaskArrayPath());
    }
  }

  getDataContainerArray()->validateNumberOfTuples<AbstractFilter>(this, dataArrayPaths);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void KDistanceGraph::preflight()
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
void KDistanceGraph::execute()
{
  setErrorCondition(0);
  setWarningCondition(0);
  dataCheck();
  if(getErrorCondition() < 0)
  {
    return;
  }

  if(m_UseMask)
  {
    EXECUTE_TEMPLATE(this, KDistanceTemplate, m_InDataPtr.lock(), this, m_InDataPtr.lock(), m_KDistanceArrayPtr.lock(), m_MaskPtr.lock(), m_MinDist, m_DistanceMetric)
  }
  else
  {
    size_t numTuples = m_InDataPtr.lock()->getNumberOfTuples();
    BoolArrayType::Pointer tmpMask = BoolArrayType::CreateArray(numTuples, "_INTERNAL_USE_ONLY_tmpMask");
    tmpMask->initializeWithValue(true);
    EXECUTE_TEMPLATE(this, KDistanceTemplate, m_InDataPtr.lock(), this, m_InDataPtr.lock(), m_KDistanceArrayPtr.lock(), tmpMask, m_MinDist, m_DistanceMetric)
  }

  notifyStatusMessage(getHumanLabel(), "Complete");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AbstractFilter::Pointer KDistanceGraph::newFilterInstance(bool copyFilterParameters) const
{
  KDistanceGraph::Pointer filter = KDistanceGraph::New();
  if(copyFilterParameters)
  {
    copyFilterParameterInstanceVariables(filter.get());
  }
  return filter;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString KDistanceGraph::getCompiledLibraryName() const
{
  return DREAM3DReviewConstants::DREAM3DReviewBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString KDistanceGraph::getBrandingString() const
{
  return "DREAM3DReview";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString KDistanceGraph::getFilterVersion() const
{
  QString version;
  QTextStream vStream(&version);
  vStream << DREAM3DReview::Version::Major() << "." << DREAM3DReview::Version::Minor() << "." << DREAM3DReview::Version::Patch();
  return version;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString KDistanceGraph::getGroupName() const
{
  return DREAM3DReviewConstants::FilterGroups::DREAM3DReviewFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QUuid KDistanceGraph::getUuid()
{
  return QUuid("{247ea39f-cb50-5dbb-bb9d-ecde12377fed}");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString KDistanceGraph::getSubGroupName() const
{
  return DREAM3DReviewConstants::FilterSubGroups::ClusteringFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString KDistanceGraph::getHumanLabel() const
{
  return "K Distance Graph";
}
