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

#include "AverageEdgeFaceCellArrayToVertexArray.h"

#include <QtCore/QTextStream>

#include "SIMPLib/Common/Constants.h"

#include "SIMPLib/Common/TemplateHelpers.h"
#include "SIMPLib/FilterParameters/AbstractFilterParametersReader.h"
#include "SIMPLib/FilterParameters/DataArrayCreationFilterParameter.h"
#include "SIMPLib/FilterParameters/DataArraySelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/SeparatorFilterParameter.h"
#include "SIMPLib/Geometry/EdgeGeom.h"
#include "SIMPLib/Geometry/GeometryHelpers.h"
#include "SIMPLib/Geometry/IGeometry2D.h"
#include "SIMPLib/Geometry/IGeometry3D.h"
#include "SIMPLib/DataContainers/DataContainerArray.h"
#include "SIMPLib/DataContainers/DataContainer.h"

#include "DREAM3DReview/DREAM3DReviewConstants.h"
#include "DREAM3DReview/DREAM3DReviewVersion.h"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AverageEdgeFaceCellArrayToVertexArray::AverageEdgeFaceCellArrayToVertexArray()
: m_SelectedArrayPath("", "", "")
, m_AverageVertexArrayPath(SIMPL::Defaults::VertexDataContainerName, SIMPL::Defaults::VertexAttributeMatrixName, "")
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AverageEdgeFaceCellArrayToVertexArray::~AverageEdgeFaceCellArrayToVertexArray()
= default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void AverageEdgeFaceCellArrayToVertexArray::setupFilterParameters()
{
  FilterParameterVectorType parameters;
  parameters.push_back(SeparatorFilterParameter::New("Edge/Face/Cell Data", FilterParameter::RequiredArray));
  {
    DataArraySelectionFilterParameter::RequirementType req =
        DataArraySelectionFilterParameter::CreateRequirement(SIMPL::Defaults::AnyPrimitive, SIMPL::Defaults::AnyComponentSize, AttributeMatrix::Type::Any, IGeometry::Type::Any);
    AttributeMatrix::Types amTypes = {AttributeMatrix::Type::Cell, AttributeMatrix::Type::Face, AttributeMatrix::Type::Edge};
    IGeometry::Types geomTypes = {IGeometry::Type::Edge, IGeometry::Type::Triangle, IGeometry::Type::Quad, IGeometry::Type::Tetrahedral};
    req.amTypes = amTypes;
    req.dcGeometryTypes = geomTypes;
    parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Edge/Face/Cell Array to Average", SelectedArrayPath, FilterParameter::RequiredArray, AverageEdgeFaceCellArrayToVertexArray, req));
  }
  parameters.push_back(SeparatorFilterParameter::New("Vertex Data", FilterParameter::CreatedArray));
  {
    DataArrayCreationFilterParameter::RequirementType req = DataArrayCreationFilterParameter::CreateRequirement(AttributeMatrix::Type::Vertex, IGeometry::Type::Any);
    IGeometry::Types geomTypes = {IGeometry::Type::Edge, IGeometry::Type::Triangle, IGeometry::Type::Quad, IGeometry::Type::Tetrahedral};
    req.dcGeometryTypes = geomTypes;
    parameters.push_back(SIMPL_NEW_DA_CREATION_FP("Average Vertex Array", AverageVertexArrayPath, FilterParameter::CreatedArray, AverageEdgeFaceCellArrayToVertexArray, req));
  }
  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void AverageEdgeFaceCellArrayToVertexArray::readFilterParameters(AbstractFilterParametersReader* reader, int index)
{
  reader->openFilterGroup(this, index);
  setSelectedArrayPath(reader->readDataArrayPath("SelectedArrayPath", getSelectedArrayPath()));
  setAverageVertexArrayPath(reader->readDataArrayPath("AverageVertexArrayPath", getAverageVertexArrayPath()));
  reader->closeFilterGroup();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
template <typename T>
void findCellAverage(AbstractFilter* filter, IDataArray::Pointer inDataPtr, DataArray<float>::Pointer outDataPtr, const DataContainer::Pointer& m)
{
  IGeometry::Pointer igeom = m->getGeometry();

  int32_t err = 0;

  if(!igeom->getElementsContainingVert())
  {
    err = igeom->findElementsContainingVert();
    if(err < 0)
    {
      QString ss = QObject::tr("Error computing Vertex to Element connectivity for Geometry type %1").arg(igeom->getGeometryTypeAsString());
      filter->setErrorCondition(err, ss);
      return;
    }
  }

  ElementDynamicList::Pointer elemsContainingVert = igeom->getElementsContainingVert();
  typename DataArray<T>::Pointer inputDataPtr = std::dynamic_pointer_cast<DataArray<T>>(inDataPtr);

  if(IGeometry2D::Pointer igeom2D = std::dynamic_pointer_cast<IGeometry2D>(igeom))
  {
    SharedVertexList::Pointer verts = igeom2D->getVertices();
    GeometryHelpers::Generic::AverageCellArrayValues<MeshIndexType, T, uint16_t, float>(elemsContainingVert, verts, inputDataPtr, outDataPtr);
  }
  else if(IGeometry3D::Pointer igeom3D = std::dynamic_pointer_cast<IGeometry3D>(igeom))
  {
    SharedVertexList::Pointer verts = igeom3D->getVertices();
    GeometryHelpers::Generic::AverageCellArrayValues<MeshIndexType, T, uint16_t, float>(elemsContainingVert, verts, inputDataPtr, outDataPtr);
  }
  else if(EdgeGeom::Pointer edge = std::dynamic_pointer_cast<EdgeGeom>(igeom))
  {
    SharedVertexList::Pointer verts = edge->getVertices();
    GeometryHelpers::Generic::AverageCellArrayValues<MeshIndexType, T, uint16_t, float>(elemsContainingVert, verts, inputDataPtr, outDataPtr);
  }
  else
  {
    return;
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void AverageEdgeFaceCellArrayToVertexArray::initialize()
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void AverageEdgeFaceCellArrayToVertexArray::dataCheck()
{
  clearErrorCode();
  clearWarningCode();

  IGeometry::Pointer igeom = getDataContainerArray()->getPrereqGeometryFromDataContainer<IGeometry>(this, getSelectedArrayPath().getDataContainerName());

  if(getErrorCode() < 0)
  {
    return;
  }

  IGeometry::Type geomType = igeom->getGeometryType();

  if(geomType != IGeometry::Type::Edge && geomType != IGeometry::Type::Triangle && geomType != IGeometry::Type::Quad && geomType != IGeometry::Type::Tetrahedral)
  {
    QString ss = QObject::tr("The Geometry type must be either Edge, Triangle, Quadrilateral or Tetrahedral, but the type is %1").arg(igeom->getGeometryTypeAsString());
    setErrorCondition(-11000, ss);
  }

  if(getSelectedArrayPath().getDataContainerName() != getAverageVertexArrayPath().getDataContainerName())
  {
    QString ss = QObject::tr("The input and output arrays must belong to the same Geometry (i.e., they must be part of the same Data Container)");
    setErrorCondition(-11000, ss);
  }

  AttributeMatrix::Pointer cellAttrMat = getDataContainerArray()->getPrereqAttributeMatrixFromPath(this, getSelectedArrayPath(), -301);
  AttributeMatrix::Pointer vertAttrMat = getDataContainerArray()->getPrereqAttributeMatrixFromPath(this, getAverageVertexArrayPath(), -301);

  if(getErrorCode() < 0)
  {
    return;
  }

  AttributeMatrix::Type vertAttrMatType = vertAttrMat->getType();
  AttributeMatrix::Type cellAttrMatType = cellAttrMat->getType();

  if(geomType == IGeometry::Type::Triangle || geomType == IGeometry::Type::Quad)
  {
    if(cellAttrMatType != AttributeMatrix::Type::Face)
    {
      QString ss = QObject::tr("The selected Geometry is %1, but the source Attribute Matrix is not a Face Attribute Matrix").arg(igeom->getGeometryTypeAsString());
      setErrorCondition(-11000, ss);
    }
  }
  else if(geomType == IGeometry::Type::Tetrahedral)
  {
    if(cellAttrMatType != AttributeMatrix::Type::Cell)
    {
      QString ss = QObject::tr("The selected Geometry is %1, but the source Attribute Matrix is not a Cell Attribute Matrix").arg(igeom->getGeometryTypeAsString());
      setErrorCondition(-11000, ss);
    }
  }
  else if(geomType == IGeometry::Type::Edge)
  {
    if(cellAttrMatType != AttributeMatrix::Type::Edge)
    {
      QString ss = QObject::tr("The selected Geometry is %1, but the source Attribute Matrix is not an Edge Attribute Matrix").arg(igeom->getGeometryTypeAsString());
      setErrorCondition(-11000, ss);
    }
  }

  if(vertAttrMatType != AttributeMatrix::Type::Vertex)
  {
    QString ss = QObject::tr("The destination Attribute Matrix must be a Vertex Attribute Matrix");
    setErrorCondition(-11001, ss);
  }

  m_InCellArrayPtr = getDataContainerArray()->getPrereqIDataArrayFromPath(this, getSelectedArrayPath());

  if(getErrorCode() < 0)
  {
    return;
  }

  size_t numVertices = 0;
  size_t numElements = 0;

  if(IGeometry2D::Pointer igeom2D = std::dynamic_pointer_cast<IGeometry2D>(igeom))
  {
    numVertices = igeom2D->getNumberOfVertices();
    numElements = igeom2D->getNumberOfElements();
  }
  else if(IGeometry3D::Pointer igeom3D = std::dynamic_pointer_cast<IGeometry3D>(igeom))
  {
    numVertices = igeom3D->getNumberOfVertices();
    numElements = igeom3D->getNumberOfElements();
  }
  else if(EdgeGeom::Pointer edge = std::dynamic_pointer_cast<EdgeGeom>(igeom))
  {
    numVertices = edge->getNumberOfVertices();
    numElements = edge->getNumberOfElements();
  }

  std::vector<size_t> cDims = m_InCellArrayPtr.lock()->getComponentDimensions();
  size_t numElemTuples = m_InCellArrayPtr.lock()->getNumberOfTuples();

  m_AverageVertexArrayPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<float>>(
      this, getAverageVertexArrayPath(), 0, cDims);   /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
  if(nullptr != m_AverageVertexArrayPtr.lock()) /* Validate the Weak Pointer wraps a non-nullptr pointer to a DataArray<T> object */
  {
    m_AverageVertexArray = m_AverageVertexArrayPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */

  if(getErrorCode() < 0)
  {
    return;
  }

  size_t numVertexTuples = m_AverageVertexArrayPtr.lock()->getNumberOfTuples();

  if(numVertexTuples != numVertices)
  {
    QString ss = QObject::tr("The number of Vertices in the selected Geometry is %1 and the number of tuples in the destination Attribute Matrix is %2; the Vertices and tuples must match")
                     .arg(numVertices)
                     .arg(numVertexTuples);
    setErrorCondition(-11002, ss);
  }

  if(numElemTuples != numElements)
  {
    QString ss = QObject::tr("The number of Elements in the selected Geometry is %1 and the number of tuples in the source Attribute Matrix is %2; the Elements and tuples must match")
                     .arg(numElements)
                     .arg(numElemTuples);
    setErrorCondition(-11003, ss);
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void AverageEdgeFaceCellArrayToVertexArray::preflight()
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
void AverageEdgeFaceCellArrayToVertexArray::execute()
{
  clearErrorCode();
  clearWarningCode();
  dataCheck();
  if(getErrorCode() < 0)
  {
    return;
  }

  DataContainer::Pointer m = getDataContainerArray()->getDataContainer(m_SelectedArrayPath.getDataContainerName());

  EXECUTE_FUNCTION_TEMPLATE(this, findCellAverage, m_InCellArrayPtr.lock(), this, m_InCellArrayPtr.lock(), m_AverageVertexArrayPtr.lock(), m)

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AbstractFilter::Pointer AverageEdgeFaceCellArrayToVertexArray::newFilterInstance(bool copyFilterParameters) const
{
  AverageEdgeFaceCellArrayToVertexArray::Pointer filter = AverageEdgeFaceCellArrayToVertexArray::New();
  if(copyFilterParameters)
  {
    copyFilterParameterInstanceVariables(filter.get());
  }
  return filter;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString AverageEdgeFaceCellArrayToVertexArray::getCompiledLibraryName() const
{
  return DREAM3DReviewConstants::DREAM3DReviewBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString AverageEdgeFaceCellArrayToVertexArray::getBrandingString() const
{
  return "DREAM3DReview";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString AverageEdgeFaceCellArrayToVertexArray::getFilterVersion() const
{
  QString version;
  QTextStream vStream(&version);
  vStream << DREAM3DReview::Version::Major() << "." << DREAM3DReview::Version::Minor() << "." << DREAM3DReview::Version::Patch();
  return version;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString AverageEdgeFaceCellArrayToVertexArray::getGroupName() const
{
  return DREAM3DReviewConstants::FilterGroups::DREAM3DReviewFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QUuid AverageEdgeFaceCellArrayToVertexArray::getUuid() const
{
  return QUuid("{9df4c18a-f51b-5698-8067-530d37f57c61}");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString AverageEdgeFaceCellArrayToVertexArray::getSubGroupName() const
{
  return DREAM3DReviewConstants::FilterSubGroups::StatisticsFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString AverageEdgeFaceCellArrayToVertexArray::getHumanLabel() const
{
  return "Average Edge/Face/Cell Array to Vertex Array";
}

// -----------------------------------------------------------------------------
AverageEdgeFaceCellArrayToVertexArray::Pointer AverageEdgeFaceCellArrayToVertexArray::NullPointer()
{
  return Pointer(static_cast<Self*>(nullptr));
}

// -----------------------------------------------------------------------------
std::shared_ptr<AverageEdgeFaceCellArrayToVertexArray> AverageEdgeFaceCellArrayToVertexArray::New()
{
  struct make_shared_enabler : public AverageEdgeFaceCellArrayToVertexArray
  {
  };
  std::shared_ptr<make_shared_enabler> val = std::make_shared<make_shared_enabler>();
  val->setupFilterParameters();
  return val;
}

// -----------------------------------------------------------------------------
QString AverageEdgeFaceCellArrayToVertexArray::getNameOfClass() const
{
  return QString("AverageEdgeFaceCellArrayToVertexArray");
}

// -----------------------------------------------------------------------------
QString AverageEdgeFaceCellArrayToVertexArray::ClassName()
{
  return QString("AverageEdgeFaceCellArrayToVertexArray");
}

// -----------------------------------------------------------------------------
void AverageEdgeFaceCellArrayToVertexArray::setSelectedArrayPath(const DataArrayPath& value)
{
  m_SelectedArrayPath = value;
}

// -----------------------------------------------------------------------------
DataArrayPath AverageEdgeFaceCellArrayToVertexArray::getSelectedArrayPath() const
{
  return m_SelectedArrayPath;
}

// -----------------------------------------------------------------------------
void AverageEdgeFaceCellArrayToVertexArray::setAverageVertexArrayPath(const DataArrayPath& value)
{
  m_AverageVertexArrayPath = value;
}

// -----------------------------------------------------------------------------
DataArrayPath AverageEdgeFaceCellArrayToVertexArray::getAverageVertexArrayPath() const
{
  return m_AverageVertexArrayPath;
}
