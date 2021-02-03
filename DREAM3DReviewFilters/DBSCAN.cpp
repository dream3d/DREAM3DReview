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
#include "DBSCAN.h"

#include <QtCore/QTextStream>

#include "SIMPLib/Common/Constants.h"
#include "SIMPLib/Common/TemplateHelpers.h"
#include "SIMPLib/DataContainers/DataContainer.h"
#include "SIMPLib/DataContainers/DataContainerArray.h"
#include "SIMPLib/FilterParameters/AbstractFilterParametersReader.h"
#include "SIMPLib/FilterParameters/ChoiceFilterParameter.h"
#include "SIMPLib/FilterParameters/DataArraySelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/FloatFilterParameter.h"
#include "SIMPLib/FilterParameters/IntFilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedBooleanFilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedPathCreationFilterParameter.h"
#include "SIMPLib/FilterParameters/StringFilterParameter.h"

#include "util/ClusteringAlgorithms/DBSCANTemplate.hpp"

#include "DREAM3DReview/DREAM3DReviewConstants.h"
#include "DREAM3DReview/DREAM3DReviewVersion.h"

/* Create Enumerations to allow the created Attribute Arrays to take part in renaming */
enum createdPathID : RenameDataPath::DataID_t
{
  AttributeMatrixID20 = 20,
  AttributeMatrixID21 = 21,
};

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
DBSCAN::DBSCAN()
: m_SelectedArrayPath("", "", "")
, m_UseMask(false)
, m_MaskArrayPath("", "", "")
, m_FeatureIdsArrayName("ClusterIds")
, m_FeatureAttributeMatrixName("ClusterData")
, m_Epsilon(0.01f)
, m_MinPnts(50)
, m_DistanceMetric(0)
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
DBSCAN::~DBSCAN() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void DBSCAN::setupFilterParameters()
{
  FilterParameterVectorType parameters;
  parameters.push_back(SIMPL_NEW_FLOAT_FP("Epsilon", Epsilon, FilterParameter::Parameter, DBSCAN));
  parameters.push_back(SIMPL_NEW_INTEGER_FP("Minimum Number of Points", MinPnts, FilterParameter::Parameter, DBSCAN));
  {
    ChoiceFilterParameter::Pointer parameter = ChoiceFilterParameter::New();
    parameter->setHumanLabel("Distance Metric");
    parameter->setPropertyName("DistanceMetric");
    parameter->setSetterCallback(SIMPL_BIND_SETTER(DBSCAN, this, DistanceMetric));
    parameter->setGetterCallback(SIMPL_BIND_GETTER(DBSCAN, this, DistanceMetric));
    QVector<QString> choices = DistanceTemplate::GetDistanceMetricsOptions();
    parameter->setChoices(choices);
    parameter->setCategory(FilterParameter::Parameter);
    parameters.push_back(parameter);
  }
  QStringList linkedProps("MaskArrayPath");
  parameters.push_back(SIMPL_NEW_LINKED_BOOL_FP("Use Mask", UseMask, FilterParameter::Parameter, DBSCAN, linkedProps));
  DataArraySelectionFilterParameter::RequirementType dasReq =
      DataArraySelectionFilterParameter::CreateRequirement(SIMPL::Defaults::AnyPrimitive, SIMPL::Defaults::AnyComponentSize, AttributeMatrix::Type::Any, IGeometry::Type::Any);
  parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Attribute Array to Cluster", SelectedArrayPath, FilterParameter::RequiredArray, DBSCAN, dasReq));
  dasReq = DataArraySelectionFilterParameter::CreateRequirement(SIMPL::TypeNames::Bool, 1, AttributeMatrix::Type::Any, IGeometry::Type::Any);
  parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Mask", MaskArrayPath, FilterParameter::RequiredArray, DBSCAN, dasReq));
  parameters.push_back(SIMPL_NEW_DA_WITH_LINKED_AM_FP("Cluster Ids", FeatureIdsArrayName, SelectedArrayPath, SelectedArrayPath, FilterParameter::CreatedArray, DBSCAN));
  parameters.push_back(SIMPL_NEW_AM_WITH_LINKED_DC_FP("Cluster Attribute Matrix", FeatureAttributeMatrixName, SelectedArrayPath, FilterParameter::CreatedArray, DBSCAN));
  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void DBSCAN::readFilterParameters(AbstractFilterParametersReader* reader, int index)
{
  reader->openFilterGroup(this, index);
  setSelectedArrayPath(reader->readDataArrayPath("SelectedArrayPath", getSelectedArrayPath()));
  setUseMask(reader->readValue("UseMask", getUseMask()));
  setMaskArrayPath(reader->readDataArrayPath("MaskArrayPath", getMaskArrayPath()));
  setEpsilon(reader->readValue("Epsilon", getEpsilon()));
  setMinPnts(reader->readValue("MinPnts", getMinPnts()));
  setFeatureIdsArrayName(reader->readString("FeatureIdsArrayName", getFeatureIdsArrayName()));
  setFeatureAttributeMatrixName(reader->readString("FeatureAttributeMatrixName", getFeatureAttributeMatrixName()));
  setDistanceMetric(reader->readValue("DistanceMetric", getDistanceMetric()));
  reader->closeFilterGroup();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void DBSCAN::initialize()
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void DBSCAN::dataCheck()
{
  clearErrorCode();
  clearWarningCode();
  initialize();

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
  case AttributeMatrix::Type::Vertex: {
    destAttrMatType = AttributeMatrix::Type::VertexFeature;
    break;
  }
  case AttributeMatrix::Type::Edge: {
    destAttrMatType = AttributeMatrix::Type::EdgeFeature;
    break;
  }
  case AttributeMatrix::Type::Face: {
    destAttrMatType = AttributeMatrix::Type::FaceFeature;
    break;
  }
  case AttributeMatrix::Type::Cell: {
    destAttrMatType = AttributeMatrix::Type::CellFeature;
    break;
  }
  case AttributeMatrix::Type::VertexFeature: {
    destAttrMatType = AttributeMatrix::Type::VertexEnsemble;
    break;
  }
  case AttributeMatrix::Type::EdgeFeature: {
    destAttrMatType = AttributeMatrix::Type::EdgeEnsemble;
    break;
  }
  case AttributeMatrix::Type::FaceFeature: {
    destAttrMatType = AttributeMatrix::Type::FaceEnsemble;
    break;
  }
  case AttributeMatrix::Type::CellFeature: {
    destAttrMatType = AttributeMatrix::Type::CellEnsemble;
    break;
  }
  case AttributeMatrix::Type::VertexEnsemble: {
    destAttrMatType = AttributeMatrix::Type::VertexEnsemble;
    break;
  }
  case AttributeMatrix::Type::EdgeEnsemble: {
    destAttrMatType = AttributeMatrix::Type::EdgeEnsemble;
    break;
  }
  case AttributeMatrix::Type::FaceEnsemble: {
    destAttrMatType = AttributeMatrix::Type::FaceEnsemble;
    break;
  }
  case AttributeMatrix::Type::CellEnsemble: {
    destAttrMatType = AttributeMatrix::Type::CellEnsemble;
    break;
  }
  default: {
    destAttrMatType = AttributeMatrix::Type::Generic;
    break;
  }
  }

  std::vector<size_t> tDims(1, 0);
  m->createNonPrereqAttributeMatrix(this, getFeatureAttributeMatrixName(), tDims, destAttrMatType, AttributeMatrixID21);

  if(getEpsilon() <= 0)
  {
    setErrorCondition(-5555, "Epsilon must be positive");
    return;
  }
  if(getMinPnts() <= 1)
  {
    setErrorCondition(-5556, "Minimum number of points must be greater than 1");
    return;
  }

  std::vector<size_t> cDims(1, 1);
  QVector<DataArrayPath> dataArrayPaths;

  m_InDataPtr = getDataContainerArray()->getPrereqIDataArrayFromPath(this, getSelectedArrayPath());
  if(getErrorCode() >= 0)
  {
    dataArrayPaths.push_back(getSelectedArrayPath());
  }

  DataArrayPath path(getSelectedArrayPath().getDataContainerName(), getSelectedArrayPath().getAttributeMatrixName(), getFeatureIdsArrayName());
  m_FeatureIdsPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<int32_t>>(this, path, 0, cDims);
  if(nullptr != m_FeatureIdsPtr.lock())
  {
    m_FeatureIds = m_FeatureIdsPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */

  if(m_UseMask)
  {
    m_MaskPtr = getDataContainerArray()->getPrereqArrayFromPath<BoolArrayType>(this, getMaskArrayPath(), cDims);
    if(nullptr != m_MaskPtr.lock())
    {
      m_Mask = m_MaskPtr.lock()->getPointer(0);
    } /* Now assign the raw pointer to data from the DataArray<T> object */
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
void DBSCAN::execute()
{
  dataCheck();
  if(getErrorCode() < 0)
  {
    return;
  }

  DataContainer::Pointer m = getDataContainerArray()->getDataContainer(m_SelectedArrayPath.getDataContainerName());
  AttributeMatrix::Pointer featAttrMat = m->getAttributeMatrix(m_FeatureAttributeMatrixName);
  size_t numTuples = m_InDataPtr.lock()->getNumberOfTuples();

  if(m_UseMask)
  {
    EXECUTE_TEMPLATE(this, DBSCANTemplate, m_InDataPtr.lock(), this, m_InDataPtr.lock(), m_MaskPtr.lock(), m_FeatureIdsPtr.lock(), m_Epsilon, m_MinPnts, m_DistanceMetric);
  }
  else
  {
    BoolArrayType::Pointer tmpMask = BoolArrayType::CreateArray(numTuples, std::string("_INTERNAL_USE_ONLY_tmpMask"), true);
    tmpMask->initializeWithValue(true);
    EXECUTE_TEMPLATE(this, DBSCANTemplate, m_InDataPtr.lock(), this, m_InDataPtr.lock(), tmpMask, m_FeatureIdsPtr.lock(), m_Epsilon, m_MinPnts, m_DistanceMetric);
  }

  int32_t maxCluster = std::numeric_limits<int32_t>::min();

  for(size_t i = 0; i < numTuples; i++)
  {
    if(m_FeatureIds[i] > maxCluster)
    {
      maxCluster = m_FeatureIds[i];
    }
  }

  std::vector<size_t> tDims(1, maxCluster + 1);
  featAttrMat->resizeAttributeArrays(tDims);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AbstractFilter::Pointer DBSCAN::newFilterInstance(bool copyFilterParameters) const
{
  DBSCAN::Pointer filter = DBSCAN::New();
  if(copyFilterParameters)
  {
    copyFilterParameterInstanceVariables(filter.get());
  }
  return filter;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString DBSCAN::getCompiledLibraryName() const
{
  return DREAM3DReviewConstants::DREAM3DReviewBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString DBSCAN::getBrandingString() const
{
  return "DREAM3DReview";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString DBSCAN::getFilterVersion() const
{
  QString version;
  QTextStream vStream(&version);
  vStream << DREAM3DReview::Version::Major() << "." << DREAM3DReview::Version::Minor() << "." << DREAM3DReview::Version::Patch();
  return version;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString DBSCAN::getGroupName() const
{
  return DREAM3DReviewConstants::FilterGroups::DREAM3DReviewFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QUuid DBSCAN::getUuid() const
{
  return QUuid("{c2d4f1e8-2b04-5d82-b90f-2191e8f4262e}");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString DBSCAN::getSubGroupName() const
{
  return DREAM3DReviewConstants::FilterSubGroups::ClusteringFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString DBSCAN::getHumanLabel() const
{
  return "DBSCAN";
}

// -----------------------------------------------------------------------------
DBSCAN::Pointer DBSCAN::NullPointer()
{
  return Pointer(static_cast<Self*>(nullptr));
}

// -----------------------------------------------------------------------------
std::shared_ptr<DBSCAN> DBSCAN::New()
{
  struct make_shared_enabler : public DBSCAN
  {
  };
  std::shared_ptr<make_shared_enabler> val = std::make_shared<make_shared_enabler>();
  val->setupFilterParameters();
  return val;
}

// -----------------------------------------------------------------------------
QString DBSCAN::getNameOfClass() const
{
  return QString("DBSCAN");
}

// -----------------------------------------------------------------------------
QString DBSCAN::ClassName()
{
  return QString("DBSCAN");
}

// -----------------------------------------------------------------------------
void DBSCAN::setSelectedArrayPath(const DataArrayPath& value)
{
  m_SelectedArrayPath = value;
}

// -----------------------------------------------------------------------------
DataArrayPath DBSCAN::getSelectedArrayPath() const
{
  return m_SelectedArrayPath;
}

// -----------------------------------------------------------------------------
void DBSCAN::setUseMask(bool value)
{
  m_UseMask = value;
}

// -----------------------------------------------------------------------------
bool DBSCAN::getUseMask() const
{
  return m_UseMask;
}

// -----------------------------------------------------------------------------
void DBSCAN::setMaskArrayPath(const DataArrayPath& value)
{
  m_MaskArrayPath = value;
}

// -----------------------------------------------------------------------------
DataArrayPath DBSCAN::getMaskArrayPath() const
{
  return m_MaskArrayPath;
}

// -----------------------------------------------------------------------------
void DBSCAN::setFeatureIdsArrayName(const QString& value)
{
  m_FeatureIdsArrayName = value;
}

// -----------------------------------------------------------------------------
QString DBSCAN::getFeatureIdsArrayName() const
{
  return m_FeatureIdsArrayName;
}

// -----------------------------------------------------------------------------
void DBSCAN::setFeatureAttributeMatrixName(const QString& value)
{
  m_FeatureAttributeMatrixName = value;
}

// -----------------------------------------------------------------------------
QString DBSCAN::getFeatureAttributeMatrixName() const
{
  return m_FeatureAttributeMatrixName;
}

// -----------------------------------------------------------------------------
void DBSCAN::setEpsilon(float value)
{
  m_Epsilon = value;
}

// -----------------------------------------------------------------------------
float DBSCAN::getEpsilon() const
{
  return m_Epsilon;
}

// -----------------------------------------------------------------------------
void DBSCAN::setMinPnts(int value)
{
  m_MinPnts = value;
}

// -----------------------------------------------------------------------------
int DBSCAN::getMinPnts() const
{
  return m_MinPnts;
}

// -----------------------------------------------------------------------------
void DBSCAN::setDistanceMetric(int value)
{
  m_DistanceMetric = value;
}

// -----------------------------------------------------------------------------
int DBSCAN::getDistanceMetric() const
{
  return m_DistanceMetric;
}
