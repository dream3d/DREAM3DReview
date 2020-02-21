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

#include <memory>

#include "KMedoids.h"

#include <QtCore/QTextStream>

#include "SIMPLib/Common/Constants.h"

#include "SIMPLib/Common/TemplateHelpers.h"
#include "SIMPLib/FilterParameters/AbstractFilterParametersReader.h"
#include "SIMPLib/FilterParameters/ChoiceFilterParameter.h"
#include "SIMPLib/FilterParameters/DataArraySelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/IntFilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedBooleanFilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedPathCreationFilterParameter.h"
#include "SIMPLib/FilterParameters/StringFilterParameter.h"
#include "SIMPLib/DataContainers/DataContainerArray.h"
#include "SIMPLib/DataContainers/DataContainer.h"

#include "util/ClusteringAlgorithms/KMedoidsTemplate.hpp"

#include "DREAM3DReview/DREAM3DReviewConstants.h"
#include "DREAM3DReview/DREAM3DReviewVersion.h"

/* Create Enumerations to allow the created Attribute Arrays to take part in renaming */
enum createdPathID : RenameDataPath::DataID_t
{
  AttributeMatrixID21 = 21,

  DataArrayID30 = 30,
  DataArrayID31 = 31,
};

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
KMedoids::KMedoids()
: m_SelectedArrayPath("", "", "")
, m_UseMask(false)
, m_MaskArrayPath("", "", "")
, m_FeatureIdsArrayName("ClusterIds")
, m_MedoidsArrayName("ClusterMedoids")
, m_FeatureAttributeMatrixName("ClusterData")
, m_InitClusters(1)
, m_DistanceMetric(0)
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
KMedoids::~KMedoids()
= default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void KMedoids::setupFilterParameters()
{
  FilterParameterVectorType parameters;
  parameters.push_back(SIMPL_NEW_INTEGER_FP("Number of Clusters", InitClusters, FilterParameter::Parameter, KMedoids));
  {
    ChoiceFilterParameter::Pointer parameter = ChoiceFilterParameter::New();
    parameter->setHumanLabel("Distance Metric");
    parameter->setPropertyName("DistanceMetric");
    parameter->setSetterCallback(SIMPL_BIND_SETTER(KMedoids, this, DistanceMetric));
    parameter->setGetterCallback(SIMPL_BIND_GETTER(KMedoids, this, DistanceMetric));
    QVector<QString> choices = DistanceTemplate::GetDistanceMetricsOptions();
    parameter->setChoices(choices);
    parameter->setCategory(FilterParameter::Parameter);
    parameters.push_back(parameter);
  }
  QStringList linkedProps("MaskArrayPath");
  parameters.push_back(SIMPL_NEW_LINKED_BOOL_FP("Use Mask", UseMask, FilterParameter::Parameter, KMedoids, linkedProps));
  DataArraySelectionFilterParameter::RequirementType dasReq =
      DataArraySelectionFilterParameter::CreateRequirement(SIMPL::Defaults::AnyPrimitive, SIMPL::Defaults::AnyComponentSize, AttributeMatrix::Type::Any, IGeometry::Type::Any);
  parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Attribute Array to Cluster", SelectedArrayPath, FilterParameter::RequiredArray, KMedoids, dasReq));
  dasReq = DataArraySelectionFilterParameter::CreateRequirement(SIMPL::TypeNames::Bool, 1, AttributeMatrix::Type::Any, IGeometry::Type::Any);
  parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Mask", MaskArrayPath, FilterParameter::RequiredArray, KMedoids, dasReq));
  parameters.push_back(SIMPL_NEW_DA_WITH_LINKED_AM_FP("Cluster Ids", FeatureIdsArrayName, SelectedArrayPath, SelectedArrayPath, FilterParameter::CreatedArray, KMedoids));
  parameters.push_back(SIMPL_NEW_AM_WITH_LINKED_DC_FP("Cluster Attribute Matrix", FeatureAttributeMatrixName, SelectedArrayPath, FilterParameter::CreatedArray, KMedoids));
  parameters.push_back(SIMPL_NEW_DA_WITH_LINKED_AM_FP("Cluster Medoids", MedoidsArrayName, SelectedArrayPath, FeatureAttributeMatrixName, FilterParameter::CreatedArray, KMedoids));
  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void KMedoids::readFilterParameters(AbstractFilterParametersReader* reader, int index)
{
  reader->openFilterGroup(this, index);
  setSelectedArrayPath(reader->readDataArrayPath("SelectedArrayPath", getSelectedArrayPath()));
  setUseMask(reader->readValue("UseMask", getUseMask()));
  setMaskArrayPath(reader->readDataArrayPath("MaskArrayPath", getMaskArrayPath()));
  setFeatureIdsArrayName(reader->readString("FeatureIdsArrayName", getFeatureIdsArrayName()));
  setMedoidsArrayName(reader->readString("MedoidsArrayName", getMedoidsArrayName()));
  setFeatureAttributeMatrixName(reader->readString("FeatureAttributeMatrixName", getFeatureAttributeMatrixName()));
  setInitClusters(reader->readValue("InitClusters", getInitClusters()));
  setDistanceMetric(reader->readValue("DistanceMetric", getDistanceMetric()));
  reader->closeFilterGroup();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void KMedoids::initialize()
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void KMedoids::dataCheck()
{
  clearErrorCode();
  clearWarningCode();
  initialize();

  if(getInitClusters() < 1)
  {
    setErrorCondition(-5555, "Must have at least 1 cluster");
  }

  DataContainer::Pointer m = getDataContainerArray()->getPrereqDataContainer(this, getSelectedArrayPath().getDataContainerName(), false);
  AttributeMatrix::Pointer attrMat = getDataContainerArray()->getPrereqAttributeMatrixFromPath(this, getSelectedArrayPath(), -301);

  if(getErrorCode() < 0)
  {
    return;
  }

  AttributeMatrix::Type attrMatType = attrMat->getType();
  AttributeMatrix::Type destAttrMatType = AttributeMatrix::Type::Unknown;

  switch(attrMatType)
  {
  case AttributeMatrix::Type::Vertex:
  {
    destAttrMatType = AttributeMatrix::Type::VertexFeature;
    break;
  }
  case AttributeMatrix::Type::Edge:
  {
    destAttrMatType = AttributeMatrix::Type::EdgeFeature;
    break;
  }
  case AttributeMatrix::Type::Face:
  {
    destAttrMatType = AttributeMatrix::Type::FaceFeature;
    break;
  }
  case AttributeMatrix::Type::Cell:
  {
    destAttrMatType = AttributeMatrix::Type::CellFeature;
    break;
  }
  case AttributeMatrix::Type::VertexFeature:
  {
    destAttrMatType = AttributeMatrix::Type::VertexEnsemble;
    break;
  }
  case AttributeMatrix::Type::EdgeFeature:
  {
    destAttrMatType = AttributeMatrix::Type::EdgeEnsemble;
    break;
  }
  case AttributeMatrix::Type::FaceFeature:
  {
    destAttrMatType = AttributeMatrix::Type::FaceEnsemble;
    break;
  }
  case AttributeMatrix::Type::CellFeature:
  {
    destAttrMatType = AttributeMatrix::Type::CellEnsemble;
    break;
  }
  case AttributeMatrix::Type::VertexEnsemble:
  {
    destAttrMatType = AttributeMatrix::Type::VertexEnsemble;
    break;
  }
  case AttributeMatrix::Type::EdgeEnsemble:
  {
    destAttrMatType = AttributeMatrix::Type::EdgeEnsemble;
    break;
  }
  case AttributeMatrix::Type::FaceEnsemble:
  {
    destAttrMatType = AttributeMatrix::Type::FaceEnsemble;
    break;
  }
  case AttributeMatrix::Type::CellEnsemble:
  {
    destAttrMatType = AttributeMatrix::Type::CellEnsemble;
    break;
  }
  default:
  {
    destAttrMatType = AttributeMatrix::Type::Generic;
    break;
  }
  }

  std::vector<size_t> tDims(1, m_InitClusters + 1);
  m->createNonPrereqAttributeMatrix(this, getFeatureAttributeMatrixName(), tDims, destAttrMatType, AttributeMatrixID21);

  DataArrayPath tempPath;
  std::vector<size_t> cDims;
  QVector<DataArrayPath> dataArrayPaths;

  m_InDataPtr = getDataContainerArray()->getPrereqIDataArrayFromPath<IDataArray>(this, getSelectedArrayPath());
  if(getErrorCode() >= 0)
  {
    dataArrayPaths.push_back(getSelectedArrayPath());
  }
  if(m_InDataPtr.lock())
  {
    cDims = m_InDataPtr.lock()->getComponentDimensions();
    tempPath.update(getSelectedArrayPath().getDataContainerName(), getFeatureAttributeMatrixName(), getMedoidsArrayName());
    m_MedoidsArrayPtr = TemplateHelpers::CreateNonPrereqArrayFromArrayType()(this, tempPath, cDims, m_InDataPtr.lock());
  }

  cDims[0] = 1;
  tempPath.update(getSelectedArrayPath().getDataContainerName(), getSelectedArrayPath().getAttributeMatrixName(), getFeatureIdsArrayName());

  m_FeatureIdsPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<int32_t>>(this, tempPath, 0, cDims, "", DataArrayID31);
  if(m_FeatureIdsPtr.lock())
  {
    m_FeatureIds = m_FeatureIdsPtr.lock()->getPointer(0);
  }
  if(getErrorCode() >= 0)
  {
    dataArrayPaths.push_back(tempPath);
  }

  if(m_UseMask)
  {
    m_MaskPtr = getDataContainerArray()->getPrereqArrayFromPath<BoolArrayType>(this, getMaskArrayPath(), cDims);
    if(m_MaskPtr.lock())
    {
      m_Mask = m_MaskPtr.lock()->getPointer(0);
    }
    if(getErrorCode() >= 0)
    {
      dataArrayPaths.push_back(getMaskArrayPath());
    }
  }

  getDataContainerArray()->validateNumberOfTuples(this, dataArrayPaths);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void KMedoids::preflight()
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
void KMedoids::execute()
{
  clearErrorCode();
  clearWarningCode();
  dataCheck();
  if(getErrorCode() < 0)
  {
    return;
  }

  if(m_UseMask)
  {
    EXECUTE_TEMPLATE(this, KMedoidsTemplate, m_InDataPtr.lock(), this, m_InDataPtr.lock(), m_MedoidsArrayPtr.lock(), m_MaskPtr.lock(), m_InitClusters, m_FeatureIdsPtr.lock(), m_DistanceMetric)
  }
  else
  {
    size_t numTuples = m_InDataPtr.lock()->getNumberOfTuples();
    BoolArrayType::Pointer tmpMask = BoolArrayType::CreateArray(numTuples, "_INTERNAL_USE_ONLY_tmpMask", true);
    tmpMask->initializeWithValue(true);
    EXECUTE_TEMPLATE(this, KMedoidsTemplate, m_InDataPtr.lock(), this, m_InDataPtr.lock(), m_MedoidsArrayPtr.lock(), tmpMask, m_InitClusters, m_FeatureIdsPtr.lock(), m_DistanceMetric)
  }

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AbstractFilter::Pointer KMedoids::newFilterInstance(bool copyFilterParameters) const
{
  KMedoids::Pointer filter = KMedoids::New();
  if(copyFilterParameters)
  {
    copyFilterParameterInstanceVariables(filter.get());
  }
  return filter;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString KMedoids::getCompiledLibraryName() const
{
  return DREAM3DReviewConstants::DREAM3DReviewBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString KMedoids::getBrandingString() const
{
  return "DREAM3DReview";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString KMedoids::getFilterVersion() const
{
  QString version;
  QTextStream vStream(&version);
  vStream << DREAM3DReview::Version::Major() << "." << DREAM3DReview::Version::Minor() << "." << DREAM3DReview::Version::Patch();
  return version;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString KMedoids::getGroupName() const
{
  return DREAM3DReviewConstants::FilterGroups::DREAM3DReviewFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QUuid KMedoids::getUuid() const
{
  return QUuid("{f7486aa6-3049-5be7-8511-ae772b70c90b}");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString KMedoids::getSubGroupName() const
{
  return DREAM3DReviewConstants::FilterSubGroups::ClusteringFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString KMedoids::getHumanLabel() const
{
  return "K Medoids";
}

// -----------------------------------------------------------------------------
KMedoids::Pointer KMedoids::NullPointer()
{
  return Pointer(static_cast<Self*>(nullptr));
}

// -----------------------------------------------------------------------------
std::shared_ptr<KMedoids> KMedoids::New()
{
  struct make_shared_enabler : public KMedoids
  {
  };
  std::shared_ptr<make_shared_enabler> val = std::make_shared<make_shared_enabler>();
  val->setupFilterParameters();
  return val;
}

// -----------------------------------------------------------------------------
QString KMedoids::getNameOfClass() const
{
  return QString("KMedoids");
}

// -----------------------------------------------------------------------------
QString KMedoids::ClassName()
{
  return QString("KMedoids");
}

// -----------------------------------------------------------------------------
void KMedoids::setSelectedArrayPath(const DataArrayPath& value)
{
  m_SelectedArrayPath = value;
}

// -----------------------------------------------------------------------------
DataArrayPath KMedoids::getSelectedArrayPath() const
{
  return m_SelectedArrayPath;
}

// -----------------------------------------------------------------------------
void KMedoids::setUseMask(bool value)
{
  m_UseMask = value;
}

// -----------------------------------------------------------------------------
bool KMedoids::getUseMask() const
{
  return m_UseMask;
}

// -----------------------------------------------------------------------------
void KMedoids::setMaskArrayPath(const DataArrayPath& value)
{
  m_MaskArrayPath = value;
}

// -----------------------------------------------------------------------------
DataArrayPath KMedoids::getMaskArrayPath() const
{
  return m_MaskArrayPath;
}

// -----------------------------------------------------------------------------
void KMedoids::setFeatureIdsArrayName(const QString& value)
{
  m_FeatureIdsArrayName = value;
}

// -----------------------------------------------------------------------------
QString KMedoids::getFeatureIdsArrayName() const
{
  return m_FeatureIdsArrayName;
}

// -----------------------------------------------------------------------------
void KMedoids::setMedoidsArrayName(const QString& value)
{
  m_MedoidsArrayName = value;
}

// -----------------------------------------------------------------------------
QString KMedoids::getMedoidsArrayName() const
{
  return m_MedoidsArrayName;
}

// -----------------------------------------------------------------------------
void KMedoids::setFeatureAttributeMatrixName(const QString& value)
{
  m_FeatureAttributeMatrixName = value;
}

// -----------------------------------------------------------------------------
QString KMedoids::getFeatureAttributeMatrixName() const
{
  return m_FeatureAttributeMatrixName;
}

// -----------------------------------------------------------------------------
void KMedoids::setInitClusters(int value)
{
  m_InitClusters = value;
}

// -----------------------------------------------------------------------------
int KMedoids::getInitClusters() const
{
  return m_InitClusters;
}

// -----------------------------------------------------------------------------
void KMedoids::setDistanceMetric(int value)
{
  m_DistanceMetric = value;
}

// -----------------------------------------------------------------------------
int KMedoids::getDistanceMetric() const
{
  return m_DistanceMetric;
}
