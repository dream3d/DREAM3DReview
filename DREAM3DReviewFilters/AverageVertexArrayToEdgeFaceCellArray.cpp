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

#include "AverageVertexArrayToEdgeFaceCellArray.h"

#include "SIMPLib/Common/Constants.h"
#include "SIMPLib/Common/TemplateHelpers.h"
#include "SIMPLib/FilterParameters/AbstractFilterParametersReader.h"
#include "SIMPLib/FilterParameters/BooleanFilterParameter.h"
#include "SIMPLib/FilterParameters/DataArrayCreationFilterParameter.h"
#include "SIMPLib/FilterParameters/DataArraySelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/SeparatorFilterParameter.h"
#include "SIMPLib/Geometry/EdgeGeom.h"
#include "SIMPLib/Geometry/GeometryHelpers.h"
#include "SIMPLib/Geometry/QuadGeom.h"
#include "SIMPLib/Geometry/TetrahedralGeom.h"
#include "SIMPLib/Geometry/TriangleGeom.h"

#include "DREAM3DReview/DREAM3DReviewConstants.h"
#include "DREAM3DReview/DREAM3DReviewVersion.h"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AverageVertexArrayToEdgeFaceCellArray::AverageVertexArrayToEdgeFaceCellArray()
: m_SelectedArrayPath("", "", "")
, m_AverageCellArrayPath("", "", "")
, m_WeightedAverage(false)
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AverageVertexArrayToEdgeFaceCellArray::~AverageVertexArrayToEdgeFaceCellArray()
= default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void AverageVertexArrayToEdgeFaceCellArray::setupFilterParameters()
{
  FilterParameterVector parameters;
  parameters.push_back(SIMPL_NEW_BOOL_FP("Perform Weighted Average", WeightedAverage, FilterParameter::Parameter, AverageVertexArrayToEdgeFaceCellArray));
  parameters.push_back(SeparatorFilterParameter::New("Vertex Data", FilterParameter::RequiredArray));
  {
    DataArraySelectionFilterParameter::RequirementType req =
        DataArraySelectionFilterParameter::CreateRequirement(SIMPL::Defaults::AnyPrimitive, SIMPL::Defaults::AnyComponentSize, AttributeMatrix::Type::Vertex, IGeometry::Type::Any);
    IGeometry::Types geomTypes = {IGeometry::Type::Edge, IGeometry::Type::Triangle, IGeometry::Type::Quad, IGeometry::Type::Tetrahedral};
    req.dcGeometryTypes = geomTypes;
    parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Vertex Array to Average", SelectedArrayPath, FilterParameter::RequiredArray, AverageVertexArrayToEdgeFaceCellArray, req));
  }
  parameters.push_back(SeparatorFilterParameter::New("Edge/Face/Cell Data", FilterParameter::CreatedArray));
  {
    DataArrayCreationFilterParameter::RequirementType req = DataArrayCreationFilterParameter::CreateRequirement(AttributeMatrix::Type::Any, IGeometry::Type::Any);
    AttributeMatrix::Types amTypes = {AttributeMatrix::Type::Cell, AttributeMatrix::Type::Face, AttributeMatrix::Type::Edge};
    IGeometry::Types geomTypes = {IGeometry::Type::Edge, IGeometry::Type::Triangle, IGeometry::Type::Quad, IGeometry::Type::Tetrahedral};
    req.amTypes = amTypes;
    req.dcGeometryTypes = geomTypes;
    parameters.push_back(SIMPL_NEW_DA_CREATION_FP("Average Edge/Face/Cell Array", AverageCellArrayPath, FilterParameter::CreatedArray, AverageVertexArrayToEdgeFaceCellArray, req));
  }
  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void AverageVertexArrayToEdgeFaceCellArray::readFilterParameters(AbstractFilterParametersReader* reader, int index)
{
  reader->openFilterGroup(this, index);
  setSelectedArrayPath(reader->readDataArrayPath("SelectedArrayPath", getSelectedArrayPath()));
  setAverageCellArrayPath(reader->readDataArrayPath("AverageCellArrayPath", getAverageCellArrayPath()));
  setWeightedAverage(reader->readValue("WeightedAverage", getWeightedAverage()));
  reader->closeFilterGroup();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
template <typename T> void findVertexAverage(AbstractFilter* filter, IDataArray::Pointer inDataPtr, DataArray<float>::Pointer outDataPtr, const DataContainer::Pointer& m, bool weightAverage)
{
  typename DataArray<T>::Pointer inputDataPtr = std::dynamic_pointer_cast<DataArray<T>>(inDataPtr);
  IGeometry::Pointer geom = m->getGeometry();
  IGeometry::Type geomType = geom->getGeometryType();

  int32_t err = 0;

  if(weightAverage)
  {
    if(!geom->getElementCentroids())
    {
      err = geom->findElementCentroids();
      if(err < 0)
      {
        filter->setErrorCondition(err);
        QString ss = QObject::tr("Error computing Element centroids for Geometry type %1").arg(geom->getGeometryTypeAsString());
        filter->notifyErrorMessage(filter->getHumanLabel(), ss, filter->getErrorCondition());
        return;
      }
    }
    switch(geomType)
    {
    case IGeometry::Type::Edge:
    {
      EdgeGeom::Pointer edgeGeom = m->getGeometryAs<EdgeGeom>();
      SharedVertexList::Pointer vertices = edgeGeom->getVertices();
      SharedEdgeList::Pointer edges = edgeGeom->getEdges();
      FloatArrayType::Pointer centroids = edgeGeom->getElementCentroids();
      GeometryHelpers::Generic::WeightedAverageVertexArrayValues<int64_t, T>(edges, vertices, centroids, inputDataPtr, outDataPtr);
      break;
    }
    case IGeometry::Type::Triangle:
    {
      TriangleGeom::Pointer triGeom = m->getGeometryAs<TriangleGeom>();
      SharedVertexList::Pointer vertices = triGeom->getVertices();
      SharedTriList::Pointer triangles = triGeom->getTriangles();
      FloatArrayType::Pointer centroids = triGeom->getElementCentroids();
      GeometryHelpers::Generic::WeightedAverageVertexArrayValues<int64_t, T>(triangles, vertices, centroids, inputDataPtr, outDataPtr);
      break;
    }
    case IGeometry::Type::Quad:
    {
      QuadGeom::Pointer quadGeom = m->getGeometryAs<QuadGeom>();
      SharedVertexList::Pointer vertices = quadGeom->getVertices();
      SharedQuadList::Pointer quads = quadGeom->getQuads();
      FloatArrayType::Pointer centroids = quadGeom->getElementCentroids();
      GeometryHelpers::Generic::WeightedAverageVertexArrayValues<int64_t, T>(quads, vertices, centroids, inputDataPtr, outDataPtr);
      break;
    }
    case IGeometry::Type::Tetrahedral:
    {
      TetrahedralGeom::Pointer tets = m->getGeometryAs<TetrahedralGeom>();
      SharedVertexList::Pointer vertices = tets->getVertices();
      SharedTetList::Pointer tet = tets->getTetrahedra();
      FloatArrayType::Pointer centroids = tets->getElementCentroids();
      GeometryHelpers::Generic::WeightedAverageVertexArrayValues<int64_t, T>(tet, vertices, centroids, inputDataPtr, outDataPtr);
      break;
    }
    default:
    {
      break;
    }
    }
  }
  else
  {
    switch(geomType)
    {
    case IGeometry::Type::Edge:
    {
      EdgeGeom::Pointer edgeGeom = m->getGeometryAs<EdgeGeom>();
      SharedEdgeList::Pointer edges = edgeGeom->getEdges();
      GeometryHelpers::Generic::AverageVertexArrayValues<int64_t, T>(edges, inputDataPtr, outDataPtr);
      break;
    }
    case IGeometry::Type::Triangle:
    {
      TriangleGeom::Pointer triGeom = m->getGeometryAs<TriangleGeom>();
      SharedTriList::Pointer triangles = triGeom->getTriangles();
      GeometryHelpers::Generic::AverageVertexArrayValues<int64_t, T>(triangles, inputDataPtr, outDataPtr);
      break;
    }
    case IGeometry::Type::Quad:
    {
      QuadGeom::Pointer quadGeom = m->getGeometryAs<QuadGeom>();
      SharedQuadList::Pointer quads = quadGeom->getQuads();
      GeometryHelpers::Generic::AverageVertexArrayValues<int64_t, T>(quads, inputDataPtr, outDataPtr);
      break;
    }
    case IGeometry::Type::Tetrahedral:
    {
      TetrahedralGeom::Pointer tets = m->getGeometryAs<TetrahedralGeom>();
      SharedTetList::Pointer tet = tets->getTetrahedra();
      GeometryHelpers::Generic::AverageVertexArrayValues<int64_t, T>(tet, inputDataPtr, outDataPtr);
      break;
    }
    default:
    {
      break;
    }
    }
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void AverageVertexArrayToEdgeFaceCellArray::initialize()
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void AverageVertexArrayToEdgeFaceCellArray::dataCheck()
{
  setErrorCondition(0);
  setWarningCondition(0);

  IGeometry::Pointer igeom = getDataContainerArray()->getPrereqGeometryFromDataContainer<IGeometry, AbstractFilter>(this, getSelectedArrayPath().getDataContainerName());

  if(getErrorCondition() < 0)
  {
    return;
  }

  IGeometry::Type geomType = igeom->getGeometryType();

  if(geomType != IGeometry::Type::Edge && geomType != IGeometry::Type::Triangle && geomType != IGeometry::Type::Quad && geomType != IGeometry::Type::Tetrahedral)
  {
    setErrorCondition(-11000);
    QString ss = QObject::tr("The Geometry type must be either Edge, Triangle, Quadrialteral or Tetrahedral, but the type is %1").arg(igeom->getGeometryTypeAsString());
    notifyErrorMessage(getHumanLabel(), ss, getErrorCondition());
    return;
  }

  if(getSelectedArrayPath().getDataContainerName() != getAverageCellArrayPath().getDataContainerName())
  {
    setErrorCondition(-11000);
    QString ss = QObject::tr("The input and output arrays must belong to the same Geometry (i.e., they must be part of the same Data Container)");
    notifyErrorMessage(getHumanLabel(), ss, getErrorCondition());
  }

  AttributeMatrix::Pointer vertAttrMat = getDataContainerArray()->getPrereqAttributeMatrixFromPath<AbstractFilter>(this, getSelectedArrayPath(), -301);
  AttributeMatrix::Pointer cellAttrMat = getDataContainerArray()->getPrereqAttributeMatrixFromPath<AbstractFilter>(this, getAverageCellArrayPath(), -301);

  if(getErrorCondition() < 0)
  {
    return;
  }

  AttributeMatrix::Type vertAttrMatType = vertAttrMat->getType();
  AttributeMatrix::Type cellAttrMatType = cellAttrMat->getType();

  if(geomType == IGeometry::Type::Triangle || geomType == IGeometry::Type::Quad)
  {
    if(cellAttrMatType != AttributeMatrix::Type::Face)
    {
      setErrorCondition(-11000);
      QString ss = QObject::tr("The selected Data Container Geometry is %1, but the destination Attribute Matrix is not a Face Attribute Matrix").arg(igeom->getGeometryTypeAsString());
      notifyErrorMessage(getHumanLabel(), ss, getErrorCondition());
    }
  }
  else if(geomType == IGeometry::Type::Tetrahedral)
  {
    if(cellAttrMatType != AttributeMatrix::Type::Cell)
    {
      setErrorCondition(-11000);
      QString ss = QObject::tr("The selected Data Container Geometry is %1, but the destination Attribute Matrix is not an Edge Attribute Matrix").arg(igeom->getGeometryTypeAsString());
      notifyErrorMessage(getHumanLabel(), ss, getErrorCondition());
    }
  }
  else if(geomType == IGeometry::Type::Edge)
  {
    if(cellAttrMatType != AttributeMatrix::Type::Edge)
    {
      setErrorCondition(-11000);
      QString ss = QObject::tr("The selected Data Container Geometry is %1, but the destination Attribute Matrix is not an Edge Attribute Matrix").arg(igeom->getGeometryTypeAsString());
      notifyErrorMessage(getHumanLabel(), ss, getErrorCondition());
    }
  }

  if(vertAttrMatType != AttributeMatrix::Type::Vertex)
  {
    setErrorCondition(-11001);
    QString ss = QObject::tr("The source Attribute Matrix must be a Vertex Attribute Matrix");
    notifyErrorMessage(getHumanLabel(), ss, getErrorCondition());
  }

  m_InVertexArrayPtr = getDataContainerArray()->getPrereqIDataArrayFromPath<IDataArray, AbstractFilter>(this, getSelectedArrayPath());

  if(getErrorCondition() < 0)
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

  QVector<size_t> cDims = m_InVertexArrayPtr.lock()->getComponentDimensions();
  size_t numVertexTuples = m_InVertexArrayPtr.lock()->getNumberOfTuples();

  m_AverageCellArrayPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<float>, AbstractFilter, float>(
      this, getAverageCellArrayPath(), 0, cDims);   /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
  if(nullptr != m_AverageCellArrayPtr.lock()) /* Validate the Weak Pointer wraps a non-nullptr pointer to a DataArray<T> object */
  {
    m_AverageCellArray = m_AverageCellArrayPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */

  if(getErrorCondition() < 0)
  {
    return;
  }

  size_t numElemTuples = m_AverageCellArrayPtr.lock()->getNumberOfTuples();

  if(numVertexTuples != numVertices)
  {
    setErrorCondition(-11002);
    QString ss = QObject::tr("The number of Vertices in the selected Geometry is %1 and the number of tuples in the source Attribute Matrix is %2; the Vertices and tuples must match")
                     .arg(numVertices)
                     .arg(numVertexTuples);
    notifyErrorMessage(getHumanLabel(), ss, getErrorCondition());
  }

  if(numElemTuples != numElements)
  {
    setErrorCondition(-11003);
    QString ss = QObject::tr("The number of Elements in the selected Geometry is %1 and the number of tuples in the destination Attribute Matrix is %2; the Elements and tuples must match")
                     .arg(numElements)
                     .arg(numElemTuples);
    notifyErrorMessage(getHumanLabel(), ss, getErrorCondition());
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void AverageVertexArrayToEdgeFaceCellArray::preflight()
{
  setInPreflight(true);
  Q_EMIT preflightAboutToExecute();
  Q_EMIT updateFilterParameters(this);
  dataCheck();
  Q_EMIT preflightExecuted();
  setInPreflight(false);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void AverageVertexArrayToEdgeFaceCellArray::execute()
{
  setErrorCondition(0);
  setWarningCondition(0);
  dataCheck();
  if(getErrorCondition() < 0)
  {
    return;
  }

  DataContainer::Pointer m = getDataContainerArray()->getDataContainer(m_SelectedArrayPath.getDataContainerName());

  EXECUTE_FUNCTION_TEMPLATE(this, findVertexAverage, m_InVertexArrayPtr.lock(), this, m_InVertexArrayPtr.lock(), m_AverageCellArrayPtr.lock(), m, getWeightedAverage())

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AbstractFilter::Pointer AverageVertexArrayToEdgeFaceCellArray::newFilterInstance(bool copyFilterParameters) const
{
  AverageVertexArrayToEdgeFaceCellArray::Pointer filter = AverageVertexArrayToEdgeFaceCellArray::New();
  if(copyFilterParameters)
  {
    copyFilterParameterInstanceVariables(filter.get());
  }
  return filter;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString AverageVertexArrayToEdgeFaceCellArray::getCompiledLibraryName() const
{
  return DREAM3DReviewConstants::DREAM3DReviewBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString AverageVertexArrayToEdgeFaceCellArray::getBrandingString() const
{
  return "DREAM3DReview";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString AverageVertexArrayToEdgeFaceCellArray::getFilterVersion() const
{
  QString version;
  QTextStream vStream(&version);
  vStream << DREAM3DReview::Version::Major() << "." << DREAM3DReview::Version::Minor() << "." << DREAM3DReview::Version::Patch();
  return version;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString AverageVertexArrayToEdgeFaceCellArray::getGroupName() const
{
  return DREAM3DReviewConstants::FilterGroups::DREAM3DReviewFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QUuid AverageVertexArrayToEdgeFaceCellArray::getUuid()
{
  return QUuid("{caed6d53-6846-523f-a837-0692b3a81570}");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString AverageVertexArrayToEdgeFaceCellArray::getSubGroupName() const
{
  return DREAM3DReviewConstants::FilterSubGroups::StatisticsFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString AverageVertexArrayToEdgeFaceCellArray::getHumanLabel() const
{
  return "Average Vertex Array to Edge/Face/Cell Array";
}
