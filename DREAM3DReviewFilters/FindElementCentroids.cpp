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

#include "FindElementCentroids.h"

#include <cstring>

#include <QtCore/QTextStream>

#include "SIMPLib/Common/Constants.h"

#include "SIMPLib/FilterParameters/AbstractFilterParametersReader.h"
#include "SIMPLib/FilterParameters/DataArrayCreationFilterParameter.h"
#include "SIMPLib/FilterParameters/DataContainerCreationFilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedBooleanFilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedPathCreationFilterParameter.h"
#include "SIMPLib/FilterParameters/SeparatorFilterParameter.h"
#include "SIMPLib/FilterParameters/StringFilterParameter.h"
#include "SIMPLib/Geometry/VertexGeom.h"
#include "SIMPLib/DataContainers/DataContainerArray.h"
#include "SIMPLib/DataContainers/DataContainer.h"

#include "DREAM3DReview/DREAM3DReviewConstants.h"
#include "DREAM3DReview/DREAM3DReviewVersion.h"

enum createdPathID : RenameDataPath::DataID_t
{
  AttributeMatrixID21 = 21,

  DataArrayID31 = 31,

  DataContainerID = 1
};

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
FindElementCentroids::FindElementCentroids()
: m_CellCentroidsArrayPath("", "", "Centroids")
, m_CreateVertexDataContainer(false)
, m_NewDataContainerName(SIMPL::Defaults::VertexDataContainerName)
, m_VertexAttributeMatrixName(SIMPL::Defaults::VertexAttributeMatrixName)
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
FindElementCentroids::~FindElementCentroids()
= default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FindElementCentroids::setupFilterParameters()
{
  FilterParameterVectorType parameters;
  QStringList linkedProps;
  linkedProps << "NewDataContainerName"
              << "VertexAttributeMatrixName";
  parameters.push_back(SIMPL_NEW_LINKED_BOOL_FP("Create Vertex Geometry from Centroids", CreateVertexDataContainer, FilterParameter::Parameter, FindElementCentroids, linkedProps));
  parameters.push_back(SIMPL_NEW_DC_CREATION_FP("Vertex Data Container", NewDataContainerName, FilterParameter::CreatedArray, FindElementCentroids));
  parameters.push_back(SIMPL_NEW_AM_WITH_LINKED_DC_FP("Vertex Attribute Matrix", VertexAttributeMatrixName, NewDataContainerName, FilterParameter::CreatedArray, FindElementCentroids));
  parameters.push_back(SeparatorFilterParameter::New("Element Data", FilterParameter::CreatedArray));
  {
    DataArrayCreationFilterParameter::RequirementType req = DataArrayCreationFilterParameter::CreateRequirement(AttributeMatrix::Category::Element);
    parameters.push_back(SIMPL_NEW_DA_CREATION_FP("Element Centroids", CellCentroidsArrayPath, FilterParameter::CreatedArray, FindElementCentroids, req));
  }
  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FindElementCentroids::readFilterParameters(AbstractFilterParametersReader* reader, int index)
{
  reader->openFilterGroup(this, index);
  setCreateVertexDataContainer(reader->readValue("CreateVertexDataContainer", getCreateVertexDataContainer()));
  setCellCentroidsArrayPath(reader->readDataArrayPath("CellCentroidsArrayPath", getCellCentroidsArrayPath()));
  setNewDataContainerName(reader->readDataArrayPath("NewDataContainerName", getNewDataContainerName()));
  reader->closeFilterGroup();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FindElementCentroids::initialize()
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FindElementCentroids::dataCheck()
{
  clearErrorCode();
  clearWarningCode();

  IGeometry::Pointer geom = getDataContainerArray()->getPrereqGeometryFromDataContainer<IGeometry>(this, getCellCentroidsArrayPath().getDataContainerName());

  if(getErrorCode() < 0)
  {
    return;
  }

  IGeometry::Type geomType = geom->getGeometryType();
  size_t numElements = geom->getNumberOfElements();

  AttributeMatrix::Pointer attrMat = getDataContainerArray()->getPrereqAttributeMatrixFromPath(this, getCellCentroidsArrayPath(), -301);
  if(getErrorCode() < 0)
  {
    return;
  }

  AttributeMatrix::Type attrMatType = attrMat->getType();

  if(geomType == IGeometry::Type::Tetrahedral)
  {
    if(attrMatType != AttributeMatrix::Type::Cell)
    {
      QString ss = QObject::tr("The selected Data Container Geometry is %1, but the selected Attribute Matrix is not a Cell Attribute Matrix").arg(geom->getGeometryTypeAsString());
      setErrorCondition(-11000, ss);
    }
  }
  else if(geomType == IGeometry::Type::Triangle || geomType == IGeometry::Type::Quad)
  {
    if(attrMatType != AttributeMatrix::Type::Face)
    {
      QString ss = QObject::tr("The selected Data Container Geometry is %1, but the selected Attribute Matrix is not a Face Attribute Matrix").arg(geom->getGeometryTypeAsString());
      setErrorCondition(-11000, ss);
    }
  }
  else if(geomType == IGeometry::Type::Edge)
  {
    if(attrMatType != AttributeMatrix::Type::Edge)
    {
      QString ss = QObject::tr("The selected Data Container Geometry is %1, but the selected Attribute Matrix is not an Edge Attribute Matrix").arg(geom->getGeometryTypeAsString());
      setErrorCondition(-11000, ss);
    }
  }
  else
  {
    QString ss = QObject::tr("The selected Data Container Geometry is %1, but this type is not supported").arg(geom->getGeometryTypeAsString());
    setErrorCondition(-11000, ss);
  }

  if(getCreateVertexDataContainer())
  {
    DataContainer::Pointer vm = getDataContainerArray()->createNonPrereqDataContainer(this, getNewDataContainerName(), DataContainerID);
    if(getErrorCode() < 0)
    {
      return;
    }
    VertexGeom::Pointer vertices = VertexGeom::CreateGeometry(static_cast<int64_t>(numElements), SIMPL::Geometry::VertexGeometry, !getInPreflight());
    vm->setGeometry(vertices);

    std::vector<size_t> tDims(1, numElements);
    vm->createNonPrereqAttributeMatrix(this, getVertexAttributeMatrixName(), tDims, AttributeMatrix::Type::Vertex, AttributeMatrixID21);
  }

  std::vector<size_t> cDims(1, 3);

  m_CellCentroidsArrayPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<float>>(this, getCellCentroidsArrayPath(), 0, cDims, "", DataArrayID31);
  if(m_CellCentroidsArrayPtr.lock())
  {
    m_CellCentroidsArray = m_CellCentroidsArrayPtr.lock()->getPointer(0);
  }

  if(getErrorCode() < 0)
  {
    return;
  }

  size_t numTuples = m_CellCentroidsArrayPtr.lock()->getNumberOfTuples();

  if(numTuples != numElements)
  {
    QString ss = QObject::tr("The number of Elements in the selected Geometry is %1 and the number of tuples in the destination Attribute Matrix is %2; the Elements and tuples must match")
                     .arg(numElements)
                     .arg(numTuples);
    setErrorCondition(-11002, ss);
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FindElementCentroids::preflight()
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
void FindElementCentroids::execute()
{
  clearErrorCode();
  clearWarningCode();
  dataCheck();
  if(getErrorCode() < 0)
  {
    return;
  }

  int32_t err = 0;
  IGeometry::Pointer geom = getDataContainerArray()->getDataContainer(getCellCentroidsArrayPath().getDataContainerName())->getGeometry();

  if(geom->getElementCentroids() == nullptr)
  {
    err = geom->findElementCentroids();
  }
  if(err <= 0)
  {
    QString ss = QObject::tr("Error computing Element centroids for Geometry type %1").arg(geom->getGeometryTypeAsString());
    setErrorCondition(-11004, ss);
    return;
  }

  FloatArrayType::Pointer cellCentroidsPtr = geom->getElementCentroids();
  float* cellCentroids = cellCentroidsPtr->getPointer(0);
  size_t numElements = geom->getNumberOfElements();

  std::memcpy(m_CellCentroidsArray, cellCentroids, sizeof(float) * numElements * 3);

  if(getCreateVertexDataContainer())
  {
    DataContainer::Pointer vm = getDataContainerArray()->getDataContainer(getNewDataContainerName());
    VertexGeom::Pointer vertices = vm->getGeometryAs<VertexGeom>();
    vertices->resizeVertexList(numElements);
    float* vertex = vertices->getVertexPointer(0);
    std::memcpy(vertex, m_CellCentroidsArray, sizeof(float) * numElements * 3);
  }

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AbstractFilter::Pointer FindElementCentroids::newFilterInstance(bool copyFilterParameters) const
{
  FindElementCentroids::Pointer filter = FindElementCentroids::New();
  if(copyFilterParameters)
  {
    copyFilterParameterInstanceVariables(filter.get());
  }
  return filter;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString FindElementCentroids::getCompiledLibraryName() const
{
  return DREAM3DReviewConstants::DREAM3DReviewBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString FindElementCentroids::getBrandingString() const
{
  return "DREAM3DReview";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString FindElementCentroids::getFilterVersion() const
{
  QString version;
  QTextStream vStream(&version);
  vStream << DREAM3DReview::Version::Major() << "." << DREAM3DReview::Version::Minor() << "." << DREAM3DReview::Version::Patch();
  return version;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString FindElementCentroids::getGroupName() const
{
  return DREAM3DReviewConstants::FilterGroups::DREAM3DReviewFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QUuid FindElementCentroids::getUuid() const
{
  return QUuid("{ef37f78b-bc9a-5884-b372-d882df6abe74}");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString FindElementCentroids::getSubGroupName() const
{
  return DREAM3DReviewConstants::FilterSubGroups::GeometryFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString FindElementCentroids::getHumanLabel() const
{
  return "Find Element Centroids";
}

// -----------------------------------------------------------------------------
FindElementCentroids::Pointer FindElementCentroids::NullPointer()
{
  return Pointer(static_cast<Self*>(nullptr));
}

// -----------------------------------------------------------------------------
std::shared_ptr<FindElementCentroids> FindElementCentroids::New()
{
  struct make_shared_enabler : public FindElementCentroids
  {
  };
  std::shared_ptr<make_shared_enabler> val = std::make_shared<make_shared_enabler>();
  val->setupFilterParameters();
  return val;
}

// -----------------------------------------------------------------------------
QString FindElementCentroids::getNameOfClass() const
{
  return QString("FindElementCentroids");
}

// -----------------------------------------------------------------------------
QString FindElementCentroids::ClassName()
{
  return QString("FindElementCentroids");
}

// -----------------------------------------------------------------------------
void FindElementCentroids::setCellCentroidsArrayPath(const DataArrayPath& value)
{
  m_CellCentroidsArrayPath = value;
}

// -----------------------------------------------------------------------------
DataArrayPath FindElementCentroids::getCellCentroidsArrayPath() const
{
  return m_CellCentroidsArrayPath;
}

// -----------------------------------------------------------------------------
void FindElementCentroids::setCreateVertexDataContainer(bool value)
{
  m_CreateVertexDataContainer = value;
}

// -----------------------------------------------------------------------------
bool FindElementCentroids::getCreateVertexDataContainer() const
{
  return m_CreateVertexDataContainer;
}

// -----------------------------------------------------------------------------
void FindElementCentroids::setNewDataContainerName(const DataArrayPath& value)
{
  m_NewDataContainerName = value;
}

// -----------------------------------------------------------------------------
DataArrayPath FindElementCentroids::getNewDataContainerName() const
{
  return m_NewDataContainerName;
}

// -----------------------------------------------------------------------------
void FindElementCentroids::setVertexAttributeMatrixName(const QString& value)
{
  m_VertexAttributeMatrixName = value;
}

// -----------------------------------------------------------------------------
QString FindElementCentroids::getVertexAttributeMatrixName() const
{
  return m_VertexAttributeMatrixName;
}
