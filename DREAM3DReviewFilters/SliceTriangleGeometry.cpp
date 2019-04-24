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
#include "SliceTriangleGeometry.h"

#include "SIMPLib/Common/Constants.h"
#include "SIMPLib/FilterParameters/AbstractFilterParametersReader.h"
#include "SIMPLib/FilterParameters/DataArraySelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/DataContainerSelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/FloatFilterParameter.h"
#include "SIMPLib/FilterParameters/IntFilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedBooleanFilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedChoicesFilterParameter.h"
#include "SIMPLib/FilterParameters/StringFilterParameter.h"
#include "SIMPLib/Geometry/EdgeGeom.h"
#include "SIMPLib/Geometry/TriangleGeom.h"
#include "SIMPLib/Math/GeometryMath.h"

#include "OrientationLib/OrientationMath/OrientationTransforms.hpp"

#include "DREAM3DReview/DREAM3DReviewConstants.h"
#include "DREAM3DReview/DREAM3DReviewVersion.h"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
SliceTriangleGeometry::SliceTriangleGeometry()
: m_CADDataContainerName("TriangleDataContainer")
, m_SliceDataContainerName("SliceDataContainer")
, m_EdgeAttributeMatrixName("EdgeData")
, m_SliceAttributeMatrixName("SliceData")
, m_SliceIdArrayName("SliceIds")
, m_AreasArrayName("SliceAreas")
, m_PerimetersArrayName("SlicePerimeters")
, m_HaveRegionIds(false)
, m_RegionIdArrayPath("", "", "")
, m_SliceResolution(1.0f)
, m_Zstart(0.0)
, m_Zend(0.0)
, m_SliceRange(0)
{
  m_SliceDirection[0] = 0.0;
  m_SliceDirection[1] = 0.0;
  m_SliceDirection[2] = 1.0;
}

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

  parameters.push_back(SIMPL_NEW_FLOAT_VEC3_FP("Slice Direction (ijk)", SliceDirection, FilterParameter::Parameter, SliceTriangleGeometry));
  {
    QVector<QString> choices;
    choices.push_back("Full Range");
    choices.push_back("User Defined Range");

    QStringList linkedChoiceProps;
    linkedChoiceProps << "Zstart"
                      << "Zend";

    LinkedChoicesFilterParameter::Pointer parameter = LinkedChoicesFilterParameter::New();
    parameter->setHumanLabel("Slice Range");
    parameter->setPropertyName("SliceRange");
    parameter->setSetterCallback(SIMPL_BIND_SETTER(SliceTriangleGeometry, this, SliceRange));
    parameter->setGetterCallback(SIMPL_BIND_GETTER(SliceTriangleGeometry, this, SliceRange));
    parameter->setChoices(choices);
    parameter->setLinkedProperties(linkedChoiceProps);
    parameter->setCategory(FilterParameter::Parameter);
    parameters.push_back(parameter);
  }
  parameters.push_back(SIMPL_NEW_FLOAT_FP("Slicing Start", Zstart, FilterParameter::Parameter, SliceTriangleGeometry, 1));
  parameters.push_back(SIMPL_NEW_FLOAT_FP("Slicing End", Zend, FilterParameter::Parameter, SliceTriangleGeometry, 1));

  parameters.push_back(SIMPL_NEW_FLOAT_FP("Slice Spacing", SliceResolution, FilterParameter::Parameter, SliceTriangleGeometry));

  QStringList linkedProps("RegionIdArrayPath");
  parameters.push_back(SIMPL_NEW_LINKED_BOOL_FP("Have Region Ids", HaveRegionIds, FilterParameter::Parameter, SliceTriangleGeometry, linkedProps));
  linkedProps.clear();
  DataContainerSelectionFilterParameter::RequirementType dcsReq;
  IGeometry::Types geomTypes = {IGeometry::Type::Triangle};
  dcsReq.dcGeometryTypes = geomTypes;
  parameters.push_back(SIMPL_NEW_DC_SELECTION_FP("CAD Geometry", CADDataContainerName, FilterParameter::RequiredArray, SliceTriangleGeometry, dcsReq));
  {
    DataArraySelectionFilterParameter::RequirementType req = DataArraySelectionFilterParameter::CreateRequirement(SIMPL::TypeNames::Int32, 1, AttributeMatrix::Type::Face, IGeometry::Type::Triangle);
    parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Region Ids", RegionIdArrayPath, FilterParameter::RequiredArray, SliceTriangleGeometry, req));
  }
  parameters.push_back(SIMPL_NEW_STRING_FP("Slice Geometry", SliceDataContainerName, FilterParameter::CreatedArray, SliceTriangleGeometry));
  parameters.push_back(SIMPL_NEW_STRING_FP("Edge Attribute Matrix", EdgeAttributeMatrixName, FilterParameter::CreatedArray, SliceTriangleGeometry));
  parameters.push_back(SIMPL_NEW_STRING_FP("Slice Ids", SliceIdArrayName, FilterParameter::CreatedArray, SliceTriangleGeometry));
  parameters.push_back(SIMPL_NEW_STRING_FP("Slice Attribute Matrix", SliceAttributeMatrixName, FilterParameter::CreatedArray, SliceTriangleGeometry));
  parameters.push_back(SIMPL_NEW_STRING_FP("Areas", AreasArrayName, FilterParameter::CreatedArray, SliceTriangleGeometry));
  parameters.push_back(SIMPL_NEW_STRING_FP("Perimeters", PerimetersArrayName, FilterParameter::CreatedArray, SliceTriangleGeometry));
  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void SliceTriangleGeometry::readFilterParameters(AbstractFilterParametersReader* reader, int index)
{
  reader->openFilterGroup(this, index);
  setCADDataContainerName(reader->readDataArrayPath("CADDataContainerName", getCADDataContainerName()));
  setSliceDataContainerName(reader->readString("SliceDataContainerName", getSliceDataContainerName()));
  setEdgeAttributeMatrixName(reader->readString("EdgeAttributeMatrixName", getEdgeAttributeMatrixName()));
  setSliceAttributeMatrixName(reader->readString("SliceAttributeMatrixName", getSliceAttributeMatrixName()));
  setSliceIdArrayName(reader->readString("SliceIdArrayName", getSliceIdArrayName()));
  setHaveRegionIds(reader->readValue("HaveRegionIds", getHaveRegionIds()));
  setRegionIdArrayPath(reader->readDataArrayPath("RegionIdArrayPath", getRegionIdArrayPath()));
  setSliceDirection(reader->readFloatVec3("SliceDirection", getSliceDirection()));
  setSliceResolution(reader->readValue("SliceResolution", getSliceResolution()));
  setSliceRange(reader->readValue("SliceRange", getSliceRange()));
  setZstart(reader->readValue("Zstart", getZstart()));
  setZend(reader->readValue("Zend", getZend()));
  reader->closeFilterGroup();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void SliceTriangleGeometry::updateEdgeInstancePointers()
{
  clearErrorCode();
  clearWarningCode();

  if(nullptr != m_SliceIdPtr.lock().get()) /* Validate the Weak Pointer wraps a non-nullptr pointer to a DataArray<T> object */
  {
    m_SliceId = m_SliceIdPtr.lock()->getPointer(0);
  }                                         /* Now assign the raw pointer to data from the DataArray<T> object */
  if(nullptr != m_RegionIdPtr.lock().get()) /* Validate the Weak Pointer wraps a non-nullptr pointer to a DataArray<T> object */
  {
    m_RegionId = m_RegionIdPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void SliceTriangleGeometry::updateSliceInstancePointers()
{
  clearErrorCode();
  clearWarningCode();

  if(nullptr != m_AreaPtr.lock().get()) /* Validate the Weak Pointer wraps a non-nullptr pointer to a DataArray<T> object */
  {
    m_Area = m_AreaPtr.lock()->getPointer(0);
  }                                          /* Now assign the raw pointer to data from the DataArray<T> object */
  if(nullptr != m_PerimeterPtr.lock().get()) /* Validate the Weak Pointer wraps a non-nullptr pointer to a DataArray<T> object */
  {
    m_Perimeter = m_PerimeterPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */
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

  TriangleGeom::Pointer triangle = getDataContainerArray()->getPrereqGeometryFromDataContainer<TriangleGeom, AbstractFilter>(this, getCADDataContainerName());

  if(m_SliceRange == 1)
  {
    if(m_Zstart >= m_Zend)
    {
      QString message = QObject::tr("Z end must be larger than Z start.");
      setErrorCondition(-62102, message);
    }
  }

  DataContainer::Pointer m = getDataContainerArray()->createNonPrereqDataContainer<AbstractFilter>(this, getSliceDataContainerName());
  if(getErrorCode() < 0)
  {
    return;
  }
  SharedVertexList::Pointer vertices = EdgeGeom::CreateSharedVertexList(0);
  EdgeGeom::Pointer edge = EdgeGeom::CreateGeometry(0, vertices, SIMPL::Geometry::EdgeGeometry, !getInPreflight());
  m->setGeometry(edge);

  QVector<size_t> tDims(1, 0);
  m->createNonPrereqAttributeMatrix(this, getEdgeAttributeMatrixName(), tDims, AttributeMatrix::Type::Edge);
  m->createNonPrereqAttributeMatrix(this, getSliceAttributeMatrixName(), tDims, AttributeMatrix::Type::EdgeFeature);
  if(getErrorCode() < 0)
  {
    return;
  }

  QVector<size_t> cDims(1, 1);
  if(m_HaveRegionIds)
  {
    m_TriRegionIdPtr = getDataContainerArray()->getPrereqArrayFromPath<DataArray<int32_t>, AbstractFilter>(this, getRegionIdArrayPath(),
                                                                                                           cDims); /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
    if(nullptr != m_TriRegionIdPtr.lock().get())                                                                   /* Validate the Weak Pointer wraps a non-nullptr pointer to a DataArray<T> object */
    {
      m_TriRegionId = m_TriRegionIdPtr.lock()->getPointer(0);
    } /* Now assign the raw pointer to data from the DataArray<T> object */
    if(getErrorCode() >= 0)
    {
      dataArrays.push_back(m_TriRegionIdPtr.lock());
    }

    tempPath.update(getSliceDataContainerName(), getEdgeAttributeMatrixName(), getRegionIdArrayPath().getDataArrayName());
    m_RegionIdPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<int32_t>, AbstractFilter, int32_t>(
        this, tempPath, 0, cDims);            /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
    if(nullptr != m_RegionIdPtr.lock().get()) /* Validate the Weak Pointer wraps a non-nullptr pointer to a DataArray<T> object */
    {
      m_RegionId = m_RegionIdPtr.lock()->getPointer(0);
    } /* Now assign the raw pointer to data from the DataArray<T> object */
    if(getErrorCode() < 0)
    {
      return;
    }
  }

  tempPath.update(getSliceDataContainerName(), getEdgeAttributeMatrixName(), getSliceIdArrayName());
  m_SliceIdPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<int32_t>, AbstractFilter, int32_t>(this, tempPath, 0,
                                                                                                                    cDims); /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
  if(nullptr != m_SliceIdPtr.lock().get()) /* Validate the Weak Pointer wraps a non-nullptr pointer to a DataArray<T> object */
  {
    m_SliceId = m_SliceIdPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */
  if(getErrorCode() < 0)
  {
    return;
  }

  tempPath.update(getSliceDataContainerName(), getSliceAttributeMatrixName(), getAreasArrayName());
  m_AreaPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<float>, AbstractFilter, int32_t>(this, tempPath, 0,
                                                                                                               cDims); /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
  if(nullptr != m_AreaPtr.lock().get()) /* Validate the Weak Pointer wraps a non-nullptr pointer to a DataArray<T> object */
  {
    m_Area = m_AreaPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */
  if(getErrorCode() < 0)
  {
    return;
  }

  tempPath.update(getSliceDataContainerName(), getSliceAttributeMatrixName(), getPerimetersArrayName());
  m_PerimeterPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<float>, AbstractFilter, int32_t>(this, tempPath, 0,
                                                                                                                    cDims); /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
  if(nullptr != m_PerimeterPtr.lock().get()) /* Validate the Weak Pointer wraps a non-nullptr pointer to a DataArray<T> object */
  {
    m_Perimeter = m_PerimeterPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */
  if(getErrorCode() < 0)
  {
    return;
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void SliceTriangleGeometry::preflight()
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
void SliceTriangleGeometry::rotateVertices(unsigned int direction, float* n, int64_t numVerts, float* verts)
{
  // determine the rotation need to allow sectioning direction to be 001
  float sliceDirection[3] = {m_SliceDirection[0], m_SliceDirection[1], m_SliceDirection[2]};
  MatrixMath::Normalize3x1(sliceDirection);
  float angle = GeometryMath::AngleBetweenVectors(sliceDirection, n);
  float rotMat[3][3] = {{1.0f, 0.0f, 0.0f}, {0.0f, 1.0f, 0.0f}, {0.0f, 0.0f, 1.0f}};
  float invRotMat[3][3] = {{1.0f, 0.0f, 0.0f}, {0.0f, 1.0f, 0.0f}, {0.0f, 0.0f, 1.0f}};
  if(angle > (0.01f * SIMPLib::Constants::k_PiOver180))
  {
    float crossDirection[3];
    crossDirection[0] = m_SliceDirection[1];
    crossDirection[1] = -m_SliceDirection[0];
    crossDirection[2] = 0;
    MatrixMath::Normalize3x1(crossDirection);
    FOrientArrayType om(9);
    FOrientTransformsType::ax2om(FOrientArrayType(crossDirection[0], crossDirection[1], crossDirection[2], angle), om);
    om.toGMatrix(rotMat);

    if (direction == rotBackward)
    {
      //invert the roation matrix and copy it back into itself
      MatrixMath::Invert3x3(rotMat, invRotMat);
      MatrixMath::Copy3x3(invRotMat, rotMat);
    }

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
void SliceTriangleGeometry::determineBoundsAndNumSlices(float& minDim, float& maxDim, int64_t numTris, int64_t* tris, float* triVerts)
{
  for(int64_t i = 0; i < numTris; i++)
  {
    for(size_t j = 0; j < 3; j++)
    {
      int64_t vert = tris[3 * i + j];
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

  m_NumberOfSlices = static_cast<size_t>((maxDim - minDim) / m_SliceResolution) + 1;
}

    // -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void SliceTriangleGeometry::execute()
{
  clearErrorCode();
  clearWarningCode();
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

  int64_t* tris = triangle->getTriPointer(0);
  float* triVerts = triangle->getVertexPointer(0);
  int64_t numTris = triangle->getNumberOfTris();
  int64_t numTriVerts = triangle->getNumberOfVertices();

  //rotate CAD tiangles to get into sectioning orientation
  rotateVertices(rotForward, n, numTriVerts, triVerts);

  //determine bounds and number of slices needed for CAD geometry
  float minDim = std::numeric_limits<float>::max();
  float maxDim = -minDim;
  determineBoundsAndNumSlices(minDim, maxDim, numTris, tris, triVerts);

  float q[3] = {0.0f, 0.0f, 0.0f};
  float r[3] = {0.0f, 0.0f, 0.0f};
  float p[3] = {0.0f, 0.0f, 0.0f};
  float d = 0;

  std::vector<float> slicedVerts;
  std::vector<int32_t> sliceIds;
  int32_t regionId = 0;
  std::vector<int32_t> regionIds;
  float min_shift = m_SliceResolution / 1000.0f;

  for(int64_t i = 0; i < numTris; i++)
  {
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
    size_t firstSlice = (minTriDim - minDim) / m_SliceResolution;
    size_t lastSlice = ((maxTriDim - minDim) / m_SliceResolution) + 1;
    // get cross product of triangle vectors to get normals
    float vecAB[3];
    float vecAC[3];
    float triCross[3];
    vecAB[0] = triVerts[3 * tris[3 * i + 1]] - triVerts[3 * tris[3 * i]];
    vecAB[1] = triVerts[3 * tris[3 * i + 1] + 1] - triVerts[3 * tris[3 * i] + 1];
    vecAB[2] = triVerts[3 * tris[3 * i + 1] + 2] - triVerts[3 * tris[3 * i] + 2];
    vecAC[0] = triVerts[3 * tris[3 * i + 2]] - triVerts[3 * tris[3 * i]];
    vecAC[1] = triVerts[3 * tris[3 * i + 2] + 1] - triVerts[3 * tris[3 * i] + 1];
    vecAC[2] = triVerts[3 * tris[3 * i + 2] + 2] - triVerts[3 * tris[3 * i] + 2];
    triCross[0] = vecAB[1] * vecAC[2] - vecAB[2] * vecAC[1];
    triCross[1] = vecAB[2] * vecAC[0] - vecAB[0] * vecAC[2];
    triCross[2] = vecAB[0] * vecAC[1] - vecAB[1] * vecAC[0];
    for(size_t j = firstSlice; j <= lastSlice; j++)
    {
      int cut = 0;
      d = minDim + (m_SliceResolution * float(j)) + min_shift;
      //	  d = minDim + (m_SliceResolution * float(j));
      q[0] = triVerts[3 * tris[3 * i]];
      q[1] = triVerts[3 * tris[3 * i] + 1];
      q[2] = triVerts[3 * tris[3 * i] + 2];
      r[0] = triVerts[3 * tris[3 * i + 1]];
      r[1] = triVerts[3 * tris[3 * i + 1] + 1];
      r[2] = triVerts[3 * tris[3 * i + 1] + 2];
      char val = GeometryMath::RayIntersectsPlane(n, d, q, r, p);
      if(val == '1')
      {
        slicedVerts.push_back(p[0]);
        slicedVerts.push_back(p[1]);
        slicedVerts.push_back(p[2]);
        cut++;
      }
      r[0] = triVerts[3 * tris[3 * i + 2]];
      r[1] = triVerts[3 * tris[3 * i + 2] + 1];
      r[2] = triVerts[3 * tris[3 * i + 2] + 2];
      val = GeometryMath::RayIntersectsPlane(n, d, q, r, p);
      if(val == '1')
      {
        slicedVerts.push_back(p[0]);
        slicedVerts.push_back(p[1]);
        slicedVerts.push_back(p[2]);
        cut++;
      }
      q[0] = triVerts[3 * tris[3 * i + 1]];
      q[1] = triVerts[3 * tris[3 * i + 1] + 1];
      q[2] = triVerts[3 * tris[3 * i + 1] + 2];
      val = GeometryMath::RayIntersectsPlane(n, d, q, r, p);
      if(val == '1')
      {
        slicedVerts.push_back(p[0]);
        slicedVerts.push_back(p[1]);
        slicedVerts.push_back(p[2]);
        cut++;
      }
      if(cut == 1)
      {
        for(int k = 0; k < 3; k++)
        {
          slicedVerts.pop_back();
        }
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
        slicedVerts[size - 1] -= min_shift;
        slicedVerts[size - 4] -= min_shift;
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
      }
    }
  }

  size_t numVerts = slicedVerts.size() / 3;
  size_t numEdges = slicedVerts.size() / 6;

  DataContainer::Pointer m = getDataContainerArray()->getDataContainer(getSliceDataContainerName());
  SharedVertexList::Pointer vertices = EdgeGeom::CreateSharedVertexList(numVerts);
  EdgeGeom::Pointer edge = EdgeGeom::CreateGeometry(numEdges, vertices, SIMPL::Geometry::EdgeGeometry, !getInPreflight());
  float* verts = edge->getVertexPointer(0);
  int64_t* edges = edge->getEdgePointer(0);

  QVector<size_t> tDims(1, numEdges);
  m->getAttributeMatrix(getEdgeAttributeMatrixName())->resizeAttributeArrays(tDims);
  updateEdgeInstancePointers();
  tDims[0] = m_NumberOfSlices;
  m->getAttributeMatrix(getSliceAttributeMatrixName())->resizeAttributeArrays(tDims);
  updateSliceInstancePointers();

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

  // calculate slice areas before rotating back so we can use simple area calculation
  for(size_t i = 0; i < numEdges; i++)
  {
    int32_t sliceId = m_SliceId[i];
    float height = 0.5 * (verts[3 * (2 * i) + 1] + verts[3 * (2 * i + 1) + 1]);
    float width = (verts[3 * (2 * i + 1)] - verts[3 * (2 * i)]);
    float length = ((verts[3 * (2 * i + 1)] - verts[3 * (2 * i)]) * (verts[3 * (2 * i + 1)] - verts[3 * (2 * i)])) +
                   ((verts[3 * (2 * i + 1) + 1] - verts[3 * (2 * i) + 1]) * (verts[3 * (2 * i + 1) + 1] - verts[3 * (2 * i) + 1])) +
                   ((verts[3 * (2 * i + 1) + 2] - verts[3 * (2 * i) + 2]) * (verts[3 * (2 * i + 1) + 2] - verts[3 * (2 * i) + 2]));
    length = sqrt(length);
    float area = height * width;
    m_Area[sliceId] += area;
    m_Perimeter[sliceId] += length;
  }

  // take absolute value to ensure areas are positive (in case winding is such that areas come out negative)
  for(size_t i = 0; i < m_NumberOfSlices; i++)
  {
    m_Area[i] = fabsf(m_Area[i]);
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
const QString SliceTriangleGeometry::getCompiledLibraryName() const
{
  return DREAM3DReviewConstants::DREAM3DReviewBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString SliceTriangleGeometry::getBrandingString() const
{
  return "DREAM3DReview";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString SliceTriangleGeometry::getFilterVersion() const
{
  QString version;
  QTextStream vStream(&version);
  vStream << DREAM3DReview::Version::Major() << "." << DREAM3DReview::Version::Minor() << "." << DREAM3DReview::Version::Patch();
  return version;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString SliceTriangleGeometry::getGroupName() const
{
  return SIMPL::FilterGroups::SamplingFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString SliceTriangleGeometry::getSubGroupName() const
{
  return SIMPL::FilterSubGroups::GeometryFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString SliceTriangleGeometry::getHumanLabel() const
{
  return "Slice CAD Geometry";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QUuid SliceTriangleGeometry::getUuid()
{
  return QUuid("{222307a4-67fd-5cb5-a12e-d80f9fb970ae}");
}
