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
#include "SliceTriangleGeometry.h"

#include <QtCore/QTextStream>

#include "SIMPLib/Common/Constants.h"
#include "SIMPLib/DataContainers/DataContainer.h"
#include "SIMPLib/DataContainers/DataContainerArray.h"
#include "SIMPLib/FilterParameters/AbstractFilterParametersReader.h"
#include "SIMPLib/FilterParameters/DataArraySelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/DataContainerSelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/FloatFilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedBooleanFilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedChoicesFilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedPathCreationFilterParameter.h"
#include "SIMPLib/FilterParameters/StringFilterParameter.h"
#include "SIMPLib/Geometry/EdgeGeom.h"
#include "SIMPLib/Geometry/TriangleGeom.h"
#include "SIMPLib/Math/GeometryMath.h"
#include "SIMPLib/Math/MatrixMath.h"

#include "EbsdLib/Core/Orientation.hpp"
#include "EbsdLib/Core/OrientationTransformation.hpp"

#include "DREAM3DReview/DREAM3DReviewConstants.h"
#include "DREAM3DReview/DREAM3DReviewVersion.h"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
SliceTriangleGeometry::SliceTriangleGeometry() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
SliceTriangleGeometry::~SliceTriangleGeometry() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void SliceTriangleGeometry::setupFilterParameters()
{
  FilterParameterVectorType parameters;

  parameters.push_back(SIMPL_NEW_FLOAT_VEC3_FP("Slice Direction (ijk)", SliceDirection, FilterParameter::Category::Parameter, SliceTriangleGeometry));
  {
    std::vector<QString> choices;
    choices.push_back("Full Range");
    choices.push_back("User Defined Range");

    std::vector<QString> linkedChoiceProps;
    linkedChoiceProps.push_back("Zstart");
    linkedChoiceProps.push_back("Zend");

    LinkedChoicesFilterParameter::Pointer parameter = LinkedChoicesFilterParameter::New();
    parameter->setHumanLabel("Slice Range");
    parameter->setPropertyName("SliceRange");
    parameter->setSetterCallback(SIMPL_BIND_SETTER(SliceTriangleGeometry, this, SliceRange));
    parameter->setGetterCallback(SIMPL_BIND_GETTER(SliceTriangleGeometry, this, SliceRange));
    parameter->setChoices(choices);
    parameter->setLinkedProperties(linkedChoiceProps);
    parameter->setCategory(FilterParameter::Category::Parameter);
    parameters.push_back(parameter);
  }
  parameters.push_back(SIMPL_NEW_FLOAT_FP("Slicing Start", Zstart, FilterParameter::Category::Parameter, SliceTriangleGeometry, 1));
  parameters.push_back(SIMPL_NEW_FLOAT_FP("Slicing End", Zend, FilterParameter::Category::Parameter, SliceTriangleGeometry, 1));

  parameters.push_back(SIMPL_NEW_FLOAT_FP("Slice Spacing", SliceResolution, FilterParameter::Category::Parameter, SliceTriangleGeometry));

  std::vector<QString> linkedProps = {"RegionIdArrayPath"};
  parameters.push_back(SIMPL_NEW_LINKED_BOOL_FP("Have Region Ids", HaveRegionIds, FilterParameter::Category::Parameter, SliceTriangleGeometry, linkedProps));
  linkedProps.clear();
  DataContainerSelectionFilterParameter::RequirementType dcsReq;
  IGeometry::Types geomTypes = {IGeometry::Type::Triangle};
  dcsReq.dcGeometryTypes = geomTypes;
  parameters.push_back(SIMPL_NEW_DC_SELECTION_FP("CAD Geometry", CADDataContainerName, FilterParameter::Category::RequiredArray, SliceTriangleGeometry, dcsReq));
  {
    DataArraySelectionFilterParameter::RequirementType req = DataArraySelectionFilterParameter::CreateRequirement(SIMPL::TypeNames::Int32, 1, AttributeMatrix::Type::Face, IGeometry::Type::Triangle);
    parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Region Ids", RegionIdArrayPath, FilterParameter::Category::RequiredArray, SliceTriangleGeometry, req));
  }
  parameters.push_back(SIMPL_NEW_STRING_FP("Slice Geometry", SliceDataContainerName, FilterParameter::Category::CreatedArray, SliceTriangleGeometry));
  parameters.push_back(SIMPL_NEW_AM_WITH_LINKED_DC_FP("Edge Attribute Matrix", EdgeAttributeMatrixName, SliceDataContainerName, FilterParameter::Category::CreatedArray, SliceTriangleGeometry));
  parameters.push_back(SIMPL_NEW_DA_WITH_LINKED_AM_FP("Slice Ids", SliceIdArrayName, SliceDataContainerName, EdgeAttributeMatrixName, FilterParameter::Category::CreatedArray, SliceTriangleGeometry));
  parameters.push_back(SIMPL_NEW_AM_WITH_LINKED_DC_FP("Slice Attribute Matrix", SliceAttributeMatrixName, SliceDataContainerName, FilterParameter::Category::CreatedArray, SliceTriangleGeometry));
  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void SliceTriangleGeometry::dataCheck()
{
  clearErrorCode();
  clearWarningCode();

  DataArrayPath tempPath;
  QVector<IDataArray::Pointer> dataArrays;

  TriangleGeom::Pointer triangle = getDataContainerArray()->getPrereqGeometryFromDataContainer<TriangleGeom>(this, getCADDataContainerName());

  if(m_SliceRange == 1)
  {
    if(m_Zstart >= m_Zend)
    {
      QString message = QObject::tr("Z end must be larger than Z start.");
      setErrorCondition(-62102, message);
    }
  }

  DataContainer::Pointer m = getDataContainerArray()->createNonPrereqDataContainer(this, getSliceDataContainerName());
  if(getErrorCode() < 0)
  {
    return;
  }
  SharedVertexList::Pointer vertices = EdgeGeom::CreateSharedVertexList(0);
  EdgeGeom::Pointer edge = EdgeGeom::CreateGeometry(0, vertices, SIMPL::Geometry::EdgeGeometry, !getInPreflight());
  m->setGeometry(edge);

  std::vector<size_t> tDims(1, 0);
  m->createNonPrereqAttributeMatrix(this, getEdgeAttributeMatrixName(), tDims, AttributeMatrix::Type::Edge);
  m->createNonPrereqAttributeMatrix(this, getSliceAttributeMatrixName(), tDims, AttributeMatrix::Type::EdgeFeature);
  if(getErrorCode() < 0)
  {
    return;
  }

  std::vector<size_t> cDims(1, 1);
  if(m_HaveRegionIds)
  {
    m_TriRegionIdPtr = getDataContainerArray()->getPrereqArrayFromPath<Int32ArrayType>(this, getRegionIdArrayPath(), cDims);
    if(getErrorCode() < 0)
    {
      return;
    }
    dataArrays.push_back(m_TriRegionIdPtr.lock());

    tempPath.update(getSliceDataContainerName(), getEdgeAttributeMatrixName(), getRegionIdArrayPath().getDataArrayName());
    m_RegionIdPtr = getDataContainerArray()->createNonPrereqArrayFromPath<Int32ArrayType>(this, tempPath, 0, cDims);
    if(getErrorCode() < 0)
    {
      return;
    }
    dataArrays.push_back(m_TriRegionIdPtr.lock());
  }

  tempPath.update(getSliceDataContainerName(), getEdgeAttributeMatrixName(), getSliceIdArrayName());
  m_SliceIdPtr = getDataContainerArray()->createNonPrereqArrayFromPath<Int32ArrayType>(this, tempPath, 0, cDims);
  // If more code is placed beyond this comment then you should check for an error and return if the error < 0
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void SliceTriangleGeometry::rotateVertices(unsigned int direction, float* n, int64_t numVerts, float* verts)
{
  // determine the rotation need to allow sectioning direction to be 001
  float sliceDirection[3] = {m_SliceDirection[0], m_SliceDirection[1], m_SliceDirection[2]};
  MatrixMath::Normalize3x1(sliceDirection);
  float angle = GeometryMath::AngleBetweenVectors(sliceDirection, n);
  float rotMat[3][3] = {{1.0f, 0.0f, 0.0f}, {0.0f, 1.0f, 0.0f}, {0.0f, 0.0f, 1.0f}};
  float invRotMat[3][3] = {{1.0f, 0.0f, 0.0f}, {0.0f, 1.0f, 0.0f}, {0.0f, 0.0f, 1.0f}};
  if(angle > (0.01f * SIMPLib::Constants::k_DegToRadD))
  {
    float crossDirection[3];
    crossDirection[0] = m_SliceDirection[1];
    crossDirection[1] = -m_SliceDirection[0];
    crossDirection[2] = 0;
    MatrixMath::Normalize3x1(crossDirection);
    OrientationTransformation::ax2om<OrientationF, OrientationF>(OrientationF(crossDirection[0], crossDirection[1], crossDirection[2], angle)).toGMatrix(rotMat);

    if(direction == rotBackward)
    {
      // invert the roation matrix and copy it back into itself
      MatrixMath::Invert3x3(rotMat, invRotMat);
      MatrixMath::Copy3x3(invRotMat, rotMat);
    }

    // TODO: PARALLELIZE THIS BIT: Not sure if it will really help though
    // rotate all vertices so sectioning direction will always be 001
    float coords[3] = {0.0f, 0.0f, 0.0f};
    float newcoords[3] = {0.0f, 0.0f, 0.0f};
    for(int64_t i = 0; i < numVerts; i++)
    {
      coords[0] = verts[3 * i];
      coords[1] = verts[3 * i + 1];
      coords[2] = verts[3 * i + 2];
      MatrixMath::Multiply3x3with3x1(rotMat, coords, newcoords);
      verts[3 * i] = newcoords[0];
      verts[3 * i + 1] = newcoords[1];
      verts[3 * i + 2] = newcoords[2];
    }
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
char SliceTriangleGeometry::rayIntersectsPlane(const float d, const float* q, const float* r, float* p)
{
  double rqDelZ;
  double dqDelZ;
  double t;

  rqDelZ = r[2] - q[2];
  dqDelZ = d - q[2];

  t = dqDelZ / rqDelZ;
  for(int i = 0; i < 3; i++)
  {
    p[i] = q[i] + (t * (r[i] - q[i]));
  }
  if(t > 0.0 && t < 1.0)
  {
    return '1';
  }
  if(t == 0.0)
  {
    return 'q';
  }
  if(t == 1.0)
  {
    return 'r';
  }

  return '0';
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void SliceTriangleGeometry::determineBoundsAndNumSlices(float& minDim, float& maxDim, MeshIndexType numTris, MeshIndexType* tris, float* triVerts)
{
  for(MeshIndexType i = 0; i < numTris; i++)
  {
    for(MeshIndexType j = 0; j < 3; j++)
    {
      MeshIndexType vert = tris[3 * i + j];
      if(minDim > triVerts[3 * vert + 2])
      {
        minDim = triVerts[3 * vert + 2];
      }
      if(maxDim < triVerts[3 * vert + 2])
      {
        maxDim = triVerts[3 * vert + 2];
      }
    }
  }

  // adjust sectioning range if user selected a specific range - check that user range is within actual range
  if(m_SliceRange == 1)
  {
    if(m_Zstart > minDim)
    {
      minDim = m_Zstart;
    }
    if(m_Zend < maxDim)
    {
      maxDim = m_Zend;
    }
  }
  // TODO: Is this still correct? Why have the subtraction at all?
  m_NumberOfSlices = static_cast<MeshIndexType>((maxDim - 0.0) / m_SliceResolution) + 1;
  //  m_NumberOfSlices = static_cast<MeshIndexType>((maxDim - minDim) / m_SliceResolution) + 1;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void SliceTriangleGeometry::execute()
{
  dataCheck();
  if(getErrorCode() < 0)
  {
    return;
  }

  // geometry will be rotated so that the sectioning direction is always 001 before rotating back
  float n[3] = {0.0f, 0.0f, 0.0f};
  n[0] = 0.0f;
  n[1] = 0.0f;
  n[2] = 1.0f;

  TriangleGeom::Pointer triangle = getDataContainerArray()->getDataContainer(getCADDataContainerName())->getGeometryAs<TriangleGeom>();
  triangle->findEdges();

  MeshIndexType* tris = triangle->getTriPointer(0);
  float* triVerts = triangle->getVertexPointer(0);
  MeshIndexType numTris = triangle->getNumberOfTris();
  MeshIndexType numTriVerts = triangle->getNumberOfVertices();

  // rotate CAD tiangles to get into sectioning orientation
  rotateVertices(rotForward, n, numTriVerts, triVerts);

  // determine bounds and number of slices needed for CAD geometry
  float minDim = std::numeric_limits<float>::max();
  float maxDim = -minDim;
  determineBoundsAndNumSlices(minDim, maxDim, numTris, tris, triVerts);
  int64_t minSlice = static_cast<int64_t>(minDim / m_SliceResolution);
  int64_t maxSlice = static_cast<int64_t>(maxDim / m_SliceResolution);

  float q[3] = {0.0f, 0.0f, 0.0f};
  float r[3] = {0.0f, 0.0f, 0.0f};
  float p[3] = {0.0f, 0.0f, 0.0f};
  float corner[3] = {0.0f, 0.0f, 0.0f};
  float d = 0;

  std::vector<float> slicedVerts;
  std::vector<int32_t> sliceIds;
  std::vector<int32_t> regionIds;

  // Get an object reference to the pointer
  Int32ArrayType& m_TriRegionId = *(m_TriRegionIdPtr.lock());

  int32_t edgeCounter = 0;
  for(MeshIndexType i = 0; i < numTris; i++)
  {
    int32_t regionId = 0;
    // get region Id of this triangle (if they are available)
    if(m_HaveRegionIds)
    {
      regionId = m_TriRegionId[i];
    }
    // determine which slices would hit the triangle
    float minTriDim = std::numeric_limits<float>::max();
    float maxTriDim = -minTriDim;
    for(size_t j = 0; j < 3; j++)
    {
      int64_t vert = tris[3 * i + j];
      if(minTriDim > triVerts[3 * vert + 2])
      {
        minTriDim = triVerts[3 * vert + 2];
      }
      if(maxTriDim < triVerts[3 * vert + 2])
      {
        maxTriDim = triVerts[3 * vert + 2];
      }
    }
    if(minTriDim > maxDim || maxTriDim < minDim)
    {
      continue;
    }
    if(minTriDim < minDim)
    {
      minTriDim = minDim;
    }
    if(maxTriDim > maxDim)
    {
      maxTriDim = maxDim;
    }
    int64_t firstSlice = static_cast<int64_t>(minTriDim / m_SliceResolution);
    int64_t lastSlice = static_cast<int64_t>(maxTriDim / m_SliceResolution);
    if(firstSlice < minSlice)
    {
      firstSlice = minSlice;
    }
    if(lastSlice > maxSlice)
    {
      lastSlice = maxSlice;
    }
    // get cross product of triangle vectors to get normals
    float vecAB[3];
    float vecAC[3];
    float triCross[3];
    char val;
    vecAB[0] = triVerts[3 * tris[3 * i + 1]] - triVerts[3 * tris[3 * i]];
    vecAB[1] = triVerts[3 * tris[3 * i + 1] + 1] - triVerts[3 * tris[3 * i] + 1];
    vecAB[2] = triVerts[3 * tris[3 * i + 1] + 2] - triVerts[3 * tris[3 * i] + 2];
    vecAC[0] = triVerts[3 * tris[3 * i + 2]] - triVerts[3 * tris[3 * i]];
    vecAC[1] = triVerts[3 * tris[3 * i + 2] + 1] - triVerts[3 * tris[3 * i] + 1];
    vecAC[2] = triVerts[3 * tris[3 * i + 2] + 2] - triVerts[3 * tris[3 * i] + 2];
    triCross[0] = vecAB[1] * vecAC[2] - vecAB[2] * vecAC[1];
    triCross[1] = vecAB[2] * vecAC[0] - vecAB[0] * vecAC[2];
    triCross[2] = vecAB[0] * vecAC[1] - vecAB[1] * vecAC[0];
    for(int64_t j = firstSlice; j <= lastSlice; j++)
    {
      int cut = 0;
      bool cornerHit = false;
      d = (m_SliceResolution * float(j));
      q[0] = triVerts[3 * tris[3 * i]];
      q[1] = triVerts[3 * tris[3 * i] + 1];
      q[2] = triVerts[3 * tris[3 * i] + 2];
      r[0] = triVerts[3 * tris[3 * i + 1]];
      r[1] = triVerts[3 * tris[3 * i + 1] + 1];
      r[2] = triVerts[3 * tris[3 * i + 1] + 2];
      if(q[2] > r[2])
      {
        val = rayIntersectsPlane(d, r, q, p);
      }
      else
      {
        val = rayIntersectsPlane(d, q, r, p);
      }
      if(val == '1')
      {
        slicedVerts.push_back(p[0]);
        slicedVerts.push_back(p[1]);
        slicedVerts.push_back(p[2]);
        cut++;
      }
      else if(val == 'q' || val == 'r')
      {
        cornerHit = true;
        corner[0] = p[0];
        corner[1] = p[1];
        corner[2] = p[2];
      }
      r[0] = triVerts[3 * tris[3 * i + 2]];
      r[1] = triVerts[3 * tris[3 * i + 2] + 1];
      r[2] = triVerts[3 * tris[3 * i + 2] + 2];
      if(q[2] > r[2])
      {
        val = rayIntersectsPlane(d, r, q, p);
      }
      else
      {
        val = rayIntersectsPlane(d, q, r, p);
      }
      if(val == '1')
      {
        slicedVerts.push_back(p[0]);
        slicedVerts.push_back(p[1]);
        slicedVerts.push_back(p[2]);
        cut++;
      }
      else if(val == 'q' || val == 'r')
      {
        cornerHit = true;
        corner[0] = p[0];
        corner[1] = p[1];
        corner[2] = p[2];
      }
      q[0] = triVerts[3 * tris[3 * i + 1]];
      q[1] = triVerts[3 * tris[3 * i + 1] + 1];
      q[2] = triVerts[3 * tris[3 * i + 1] + 2];
      if(q[2] > r[2])
      {
        val = rayIntersectsPlane(d, r, q, p);
      }
      else
      {
        val = rayIntersectsPlane(d, q, r, p);
      }
      if(val == '1')
      {
        slicedVerts.push_back(p[0]);
        slicedVerts.push_back(p[1]);
        slicedVerts.push_back(p[2]);
        cut++;
      }
      else if(val == 'q' || val == 'r')
      {
        cornerHit = true;
        corner[0] = p[0];
        corner[1] = p[1];
        corner[2] = p[2];
      }
      if(cut == 1 && !cornerHit)
      {
        for(int k = 0; k < 3; k++)
        {
          slicedVerts.pop_back();
        }
      }
      if(cut == 1 && cornerHit)
      {
        slicedVerts.push_back(corner[0]);
        slicedVerts.push_back(corner[1]);
        slicedVerts.push_back(corner[2]);
        cut++;
      }
      if(cut == 3)
      {
        for(int k = 0; k < 9; k++)
        {
          slicedVerts.pop_back();
        }
      }
      if(cut == 2)
      {
        size_t size = slicedVerts.size();
        // get delta x for the current ordering of the segment
        float delX = slicedVerts[size - 6] - slicedVerts[size - 3];
        // get cross product of vec with 001 slicing direction
        if((triCross[1] > 0 && delX < 0) || (triCross[1] < 0 && delX > 0))
        {
          float temp[3] = {slicedVerts[size - 3], slicedVerts[size - 2], slicedVerts[size - 1]};
          slicedVerts[size - 3] = slicedVerts[size - 6];
          slicedVerts[size - 2] = slicedVerts[size - 5];
          slicedVerts[size - 1] = slicedVerts[size - 4];
          slicedVerts[size - 6] = temp[0];
          slicedVerts[size - 5] = temp[1];
          slicedVerts[size - 4] = temp[2];
        }
        sliceIds.push_back(j);
        if(m_HaveRegionIds)
        {
          regionIds.push_back(regionId);
        }
        edgeCounter++;
      }
    }
  }

  size_t numVerts = slicedVerts.size() / 3;
  size_t numEdges = slicedVerts.size() / 6;

  if(numVerts != (2 * numEdges))
  {
    QString message = QObject::tr("Number of sectioned vertices and edges do not make sense.  Number of Vertices: %1 and Number of Edges: %2").arg(numVerts).arg(numEdges);
    setErrorCondition(-13003, message);
    return;
  }

  DataContainer::Pointer m = getDataContainerArray()->getDataContainer(getSliceDataContainerName());
  SharedVertexList::Pointer vertices = EdgeGeom::CreateSharedVertexList(numVerts);
  EdgeGeom::Pointer edge = EdgeGeom::CreateGeometry(numEdges, vertices, SIMPL::Geometry::EdgeGeometry, !getInPreflight());
  float* verts = edge->getVertexPointer(0);
  MeshIndexType* edges = edge->getEdgePointer(0);

  std::vector<size_t> tDims(1, numEdges);
  m->getAttributeMatrix(getEdgeAttributeMatrixName())->resizeAttributeArrays(tDims);

  tDims[0] = m_NumberOfSlices;
  m->getAttributeMatrix(getSliceAttributeMatrixName())->resizeAttributeArrays(tDims);

  // Weak pointers are still good because the resize operations are affecting the internal structure of the DataArray<T>
  // and not the actual pointer to the DataArray<T> object itself.
  Int32ArrayType& m_SliceId = *(m_SliceIdPtr.lock());
  Int32ArrayType& m_RegionId = *(m_RegionIdPtr.lock());

  for(size_t i = 0; i < numEdges; i++)
  {
    edges[2 * i] = 2 * i;
    edges[2 * i + 1] = 2 * i + 1;
    verts[3 * (2 * i)] = slicedVerts[3 * (2 * i)];
    verts[3 * (2 * i) + 1] = slicedVerts[3 * (2 * i) + 1];
    verts[3 * (2 * i) + 2] = slicedVerts[3 * (2 * i) + 2];
    verts[3 * (2 * i + 1)] = slicedVerts[3 * (2 * i + 1)];
    verts[3 * (2 * i + 1) + 1] = slicedVerts[3 * (2 * i + 1) + 1];
    verts[3 * (2 * i + 1) + 2] = slicedVerts[3 * (2 * i + 1) + 2];
    m_SliceId[i] = sliceIds[i];
    if(m_HaveRegionIds)
    {
      m_RegionId[i] = regionIds[i];
    }
  }

  // rotate all CAD triangles back to original orientation
  rotateVertices(rotBackward, n, numTriVerts, triVerts);

  // rotate all edges back to original orientation
  rotateVertices(rotBackward, n, numVerts, verts);

  m->setGeometry(edge);

  notifyStatusMessage("Complete");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AbstractFilter::Pointer SliceTriangleGeometry::newFilterInstance(bool copyFilterParameters) const
{
  SliceTriangleGeometry::Pointer filter = SliceTriangleGeometry::New();
  if(copyFilterParameters)
  {
    copyFilterParameterInstanceVariables(filter.get());
  }
  return filter;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString SliceTriangleGeometry::getCompiledLibraryName() const
{
  return DREAM3DReviewConstants::DREAM3DReviewBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString SliceTriangleGeometry::getBrandingString() const
{
  return "DREAM3DReview";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString SliceTriangleGeometry::getFilterVersion() const
{
  QString version;
  QTextStream vStream(&version);
  vStream << DREAM3DReview::Version::Major() << "." << DREAM3DReview::Version::Minor() << "." << DREAM3DReview::Version::Patch();
  return version;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString SliceTriangleGeometry::getGroupName() const
{
  return SIMPL::FilterGroups::SamplingFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString SliceTriangleGeometry::getSubGroupName() const
{
  return SIMPL::FilterSubGroups::GeometryFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString SliceTriangleGeometry::getHumanLabel() const
{
  return "Slice Triangle Geometry";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QUuid SliceTriangleGeometry::getUuid() const
{
  return QUuid("{222307a4-67fd-5cb5-a12e-d80f9fb970ae}");
}

// -----------------------------------------------------------------------------
SliceTriangleGeometry::Pointer SliceTriangleGeometry::NullPointer()
{
  return Pointer(static_cast<Self*>(nullptr));
}

// -----------------------------------------------------------------------------
std::shared_ptr<SliceTriangleGeometry> SliceTriangleGeometry::New()
{
  struct make_shared_enabler : public SliceTriangleGeometry
  {
  };
  std::shared_ptr<make_shared_enabler> val = std::make_shared<make_shared_enabler>();
  val->setupFilterParameters();
  return val;
}

// -----------------------------------------------------------------------------
QString SliceTriangleGeometry::getNameOfClass() const
{
  return QString("SliceTriangleGeometry");
}

// -----------------------------------------------------------------------------
QString SliceTriangleGeometry::ClassName()
{
  return QString("SliceTriangleGeometry");
}

// -----------------------------------------------------------------------------
void SliceTriangleGeometry::setCADDataContainerName(const DataArrayPath& value)
{
  m_CADDataContainerName = value;
}

// -----------------------------------------------------------------------------
DataArrayPath SliceTriangleGeometry::getCADDataContainerName() const
{
  return m_CADDataContainerName;
}

// -----------------------------------------------------------------------------
void SliceTriangleGeometry::setSliceDataContainerName(const QString& value)
{
  m_SliceDataContainerName = value;
}

// -----------------------------------------------------------------------------
QString SliceTriangleGeometry::getSliceDataContainerName() const
{
  return m_SliceDataContainerName;
}

// -----------------------------------------------------------------------------
void SliceTriangleGeometry::setEdgeAttributeMatrixName(const QString& value)
{
  m_EdgeAttributeMatrixName = value;
}

// -----------------------------------------------------------------------------
QString SliceTriangleGeometry::getEdgeAttributeMatrixName() const
{
  return m_EdgeAttributeMatrixName;
}

// -----------------------------------------------------------------------------
void SliceTriangleGeometry::setSliceAttributeMatrixName(const QString& value)
{
  m_SliceAttributeMatrixName = value;
}

// -----------------------------------------------------------------------------
QString SliceTriangleGeometry::getSliceAttributeMatrixName() const
{
  return m_SliceAttributeMatrixName;
}

// -----------------------------------------------------------------------------
void SliceTriangleGeometry::setSliceIdArrayName(const QString& value)
{
  m_SliceIdArrayName = value;
}

// -----------------------------------------------------------------------------
QString SliceTriangleGeometry::getSliceIdArrayName() const
{
  return m_SliceIdArrayName;
}

// -----------------------------------------------------------------------------
void SliceTriangleGeometry::setHaveRegionIds(bool value)
{
  m_HaveRegionIds = value;
}

// -----------------------------------------------------------------------------
bool SliceTriangleGeometry::getHaveRegionIds() const
{
  return m_HaveRegionIds;
}

// -----------------------------------------------------------------------------
void SliceTriangleGeometry::setRegionIdArrayPath(const DataArrayPath& value)
{
  m_RegionIdArrayPath = value;
}

// -----------------------------------------------------------------------------
DataArrayPath SliceTriangleGeometry::getRegionIdArrayPath() const
{
  return m_RegionIdArrayPath;
}

// -----------------------------------------------------------------------------
void SliceTriangleGeometry::setSliceDirection(const FloatVec3Type& value)
{
  m_SliceDirection = value;
}

// -----------------------------------------------------------------------------
FloatVec3Type SliceTriangleGeometry::getSliceDirection() const
{
  return m_SliceDirection;
}

// -----------------------------------------------------------------------------
void SliceTriangleGeometry::setSliceResolution(float value)
{
  m_SliceResolution = value;
}

// -----------------------------------------------------------------------------
float SliceTriangleGeometry::getSliceResolution() const
{
  return m_SliceResolution;
}

// -----------------------------------------------------------------------------
void SliceTriangleGeometry::setZstart(float value)
{
  m_Zstart = value;
}

// -----------------------------------------------------------------------------
float SliceTriangleGeometry::getZstart() const
{
  return m_Zstart;
}

// -----------------------------------------------------------------------------
void SliceTriangleGeometry::setZend(float value)
{
  m_Zend = value;
}

// -----------------------------------------------------------------------------
float SliceTriangleGeometry::getZend() const
{
  return m_Zend;
}

// -----------------------------------------------------------------------------
void SliceTriangleGeometry::setSliceRange(int value)
{
  m_SliceRange = value;
}

// -----------------------------------------------------------------------------
int SliceTriangleGeometry::getSliceRange() const
{
  return m_SliceRange;
}
