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
#include "InterpolateMeshToRegularGrid.h"

#include <cmath>

#include "SIMPLib/Common/TemplateHelpers.h"
#include "SIMPLib/FilterParameters/AbstractFilterParametersReader.h"
#include "SIMPLib/FilterParameters/DataContainerSelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/IntFilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedChoicesFilterParameter.h"
#include "SIMPLib/FilterParameters/SeparatorFilterParameter.h"
#include "SIMPLib/FilterParameters/StringFilterParameter.h"
#include "SIMPLib/Geometry/GeometryHelpers.h"
#include "SIMPLib/Geometry/ImageGeom.h"
#include "SIMPLib/Geometry/QuadGeom.h"
#include "SIMPLib/Geometry/TriangleGeom.h"

#include "DREAM3DReview/DREAM3DReviewConstants.h"
#include "DREAM3DReview/DREAM3DReviewFilters/util/DifferentialGeometryOps.hpp"
#include "DREAM3DReview/DREAM3DReviewFilters/util/DistanceTemplate.hpp"
#include "DREAM3DReview/DREAM3DReviewVersion.h"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
InterpolateMeshToRegularGrid::InterpolateMeshToRegularGrid()
: m_SelectedDataContainerName("")
, m_InterpolatedDataContainerName("")
, m_InterpolatedAttributeMatrixName("")
, m_ScaleOrSpecifyNumCells(0)
, m_SetXDimension(10)
, m_SetYDimension(10)
, m_ScaleFactorNumCells(1)
, m_OutsideMeshIdentifier(999)
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
InterpolateMeshToRegularGrid::~InterpolateMeshToRegularGrid() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void InterpolateMeshToRegularGrid::setupFilterParameters()
{
  FilterParameterVectorType parameters;
  {
    QVector<QString> choices;
    choices.push_back("Specify Dimensions");
    choices.push_back("Automatically Scale");
    QStringList linkedProps;
    linkedProps << "SetXDimension"
                << "SetYDimension"
                << "ScaleFactorNumCells";
    LinkedChoicesFilterParameter::Pointer parameter = LinkedChoicesFilterParameter::New();
    parameter->setHumanLabel("Grid Generation Method");
    parameter->setPropertyName("ScaleOrSpecifyNumCells");
    parameter->setSetterCallback(SIMPL_BIND_SETTER(InterpolateMeshToRegularGrid, this, ScaleOrSpecifyNumCells));
    parameter->setGetterCallback(SIMPL_BIND_GETTER(InterpolateMeshToRegularGrid, this, ScaleOrSpecifyNumCells));
    parameter->setChoices(choices);
    parameter->setLinkedProperties(linkedProps);
    parameter->setCategory(FilterParameter::Parameter);
    parameters.push_back(parameter);
  }
  parameters.push_back(SIMPL_NEW_INTEGER_FP("X Dimension", SetXDimension, FilterParameter::Parameter, InterpolateMeshToRegularGrid, 0));
  parameters.push_back(SIMPL_NEW_INTEGER_FP("Y Dimension", SetYDimension, FilterParameter::Parameter, InterpolateMeshToRegularGrid, 0));
  parameters.push_back(SIMPL_NEW_INTEGER_FP("Scale Factor for Number of Elements", ScaleFactorNumCells, FilterParameter::Parameter, InterpolateMeshToRegularGrid, 1));
  parameters.push_back(SIMPL_NEW_INTEGER_FP("Outside Mesh Identifier", OutsideMeshIdentifier, FilterParameter::Parameter, InterpolateMeshToRegularGrid));
  {
    DataContainerSelectionFilterParameter::RequirementType req;
    IGeometry::Types geomTypes;
    geomTypes.push_back(IGeometry::Type::Triangle);
    geomTypes.push_back(IGeometry::Type::Quad);
    req.dcGeometryTypes = geomTypes;
    parameters.push_back(SIMPL_NEW_DC_SELECTION_FP("Data Container to Interpolate", SelectedDataContainerName, FilterParameter::RequiredArray, InterpolateMeshToRegularGrid, req));
  }
  parameters.push_back(SIMPL_NEW_STRING_FP("Interpolated Data Container", InterpolatedDataContainerName, FilterParameter::CreatedArray, InterpolateMeshToRegularGrid));
  parameters.push_back(SIMPL_NEW_STRING_FP("Interpolated Attribute Matrix", InterpolatedAttributeMatrixName, FilterParameter::CreatedArray, InterpolateMeshToRegularGrid));
  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void InterpolateMeshToRegularGrid::readFilterParameters(AbstractFilterParametersReader* reader, int index)
{
  reader->openFilterGroup(this, index);
  setSelectedDataContainerName(reader->readString("SelectedArrayPath", getSelectedDataContainerName()));
  setScaleOrSpecifyNumCells(reader->readValue("ScaleOrSpecifyNumCells", getScaleOrSpecifyNumCells()));
  setSetXDimension(reader->readValue("SetXDimension", getSetXDimension()));
  setSetYDimension(reader->readValue("SetYDimension", getSetYDimension()));
  setScaleFactorNumCells(reader->readValue("ScaleFactorNumCells", getScaleFactorNumCells()));
  setOutsideMeshIdentifier(reader->readValue("OutsideMeshIdentifier", getOutsideMeshIdentifier()));
  setInterpolatedDataContainerName(reader->readString("InterpolatedDataContainerName", getInterpolatedDataContainerName()));
  setInterpolatedAttributeMatrixName(reader->readString("InterpolatedAttributeMatrixName", getInterpolatedAttributeMatrixName()));
  reader->closeFilterGroup();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void InterpolateMeshToRegularGrid::createRegularGrid()
{
  DataContainer::Pointer m = getDataContainerArray()->getDataContainer(getSelectedDataContainerName());
  DataContainer::Pointer interpolatedDC = getDataContainerArray()->getDataContainer(getInterpolatedDataContainerName());
  IGeometry2D::Pointer geom2D = m->getGeometryAs<IGeometry2D>();
  ImageGeom::Pointer image = interpolatedDC->getGeometryAs<ImageGeom>();
  AttributeMatrix::Pointer interpolatedAttrMat = interpolatedDC->getAttributeMatrix(getInterpolatedAttributeMatrixName());

  size_t numElements = geom2D->getNumberOfElements();
  int64_t numVerts = geom2D->getNumberOfVertices();
  float* vertex = geom2D->getVertexPointer(0);

  // Currently supporting exactly 3 dimensions (Euclidean space)...
  size_t nDims = 3;

  // Find the largest/smallest (x,y,z) dimensions of the incoming data to be used to define the maximum dimensions for the regular grid
  for(size_t i = 0; i < nDims; i++)
  {
    m_MeshMaxExtents.push_back(-1.0f * std::numeric_limits<float>::max());
    m_MeshMinExtents.push_back(std::numeric_limits<float>::max());
  }

  for(int64_t i = 0; i < numVerts; i++)
  {
    if(vertex[3 * i] > m_MeshMaxExtents[0])
    {
      m_MeshMaxExtents[0] = vertex[3 * i];
    }
    if(vertex[3 * i + 1] > m_MeshMaxExtents[1])
    {
      m_MeshMaxExtents[1] = vertex[3 * i + 1];
    }
    if(vertex[3 * i + 2] > m_MeshMaxExtents[2])
    {
      m_MeshMaxExtents[2] = vertex[3 * i + 2];
    }
    if(vertex[3 * i] < m_MeshMinExtents[0])
    {
      m_MeshMinExtents[0] = vertex[3 * i];
    }
    if(vertex[3 * i + 1] < m_MeshMinExtents[1])
    {
      m_MeshMinExtents[1] = vertex[3 * i + 1];
    }
    if(vertex[3 * i + 2] < m_MeshMinExtents[2])
    {
      m_MeshMinExtents[2] = vertex[3 * i + 2];
    }
  }

  QVector<size_t> iDims(3, 0);

  if(m_ScaleOrSpecifyNumCells == 1)
  {
    // Initialize the interpolated grid dimensions using the original mesh extent aspect ratios
    float ratioYX = m_MeshMaxExtents[1] / m_MeshMaxExtents[0];
    float ratioZX = m_MeshMaxExtents[2] / m_MeshMaxExtents[1];

    iDims[0] = 10;
    iDims[1] = static_cast<size_t>(iDims[0] * ratioYX);
    iDims[2] = static_cast<size_t>(iDims[0] * ratioZX);

    for(size_t i = 0; i < nDims; i++)
    {
      if(iDims[i] == 0)
      {
        iDims[i] = 1;
      }
    }

    // Scale the number of cells in the interpolated grid to near the same number as in the original geometry
    float scaleFactor = 0.0f;
    float rootVal = 3.0f;

    if((iDims[0] * iDims[1] * iDims[2]) < numElements)
    {
      scaleFactor = static_cast<float>(numElements) / static_cast<float>(iDims[0] * iDims[1] * iDims[2]);
      for(size_t i = 0; i < nDims; i++)
      {
        if(iDims[i] == 1)
        {
          rootVal--;
        }
      }

      scaleFactor = std::pow(scaleFactor, 1 / rootVal);

      for(size_t i = 0; i < nDims; i++)
      {
        if(iDims[i] != 1)
        {
          iDims[i] *= static_cast<size_t>(scaleFactor) * m_ScaleFactorNumCells;
        }
      }
    }
  }
  else
  {
    iDims[0] = m_SetXDimension;
    iDims[1] = m_SetYDimension;
    iDims[2] = 1;
  }

  image->setDimensions(iDims[0], iDims[1], iDims[2]);

  QVector<float> iRes(3, 0.0f);
  QVector<float> iOrigin(3, 0.0f);
  float minRes = -1.0f * std::numeric_limits<float>::max();

  iRes[0] = m_MeshMaxExtents[0] / (static_cast<float>(iDims[0]));
  iRes[1] = m_MeshMaxExtents[1] / (static_cast<float>(iDims[1]));
  iRes[2] = m_MeshMaxExtents[2] / (static_cast<float>(iDims[2]));

  for(size_t i = 0; i < nDims; i++)
  {
    if(iRes[i] > minRes)
    {
      minRes = iRes[i];
    }
  }

  iOrigin[0] = m_MeshMinExtents[0];
  iOrigin[1] = m_MeshMinExtents[1];
  iOrigin[2] = m_MeshMinExtents[2];

  image->setSpacing(minRes, minRes, minRes);
  image->setOrigin(iOrigin[0], iOrigin[1], iOrigin[2]);

  interpolatedAttrMat->resizeAttributeArrays(iDims);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void InterpolateMeshToRegularGrid::interpolateMeshPointsToGrid(DataContainer::Pointer inputMesh, DataContainer::Pointer interpolatedGrid)
{
  // Get the cell centroids needed for interpolation
  IGeometry2D::Pointer geom2D = inputMesh->getGeometryAs<IGeometry2D>();
  ImageGeom::Pointer image = interpolatedGrid->getGeometryAs<ImageGeom>();

  int32_t err = 0;
  if((geom2D->getElementCentroids()).get() == nullptr)
  {
    err = geom2D->findElementCentroids();
  }
  if(err < 0)
  {
    QString ss = QObject::tr("Error computing Element centroids for Geometry type %1").arg(geom2D->getGeometryTypeAsString());
    setErrorCondition(err, ss);
    return;
  }

  FloatArrayType::Pointer cellCentroidsPtr = geom2D->getElementCentroids();
  float* cellCentroids = cellCentroidsPtr->getPointer(0);
  float* vertex = geom2D->getVertexPointer(0);
  size_t numElements = geom2D->getNumberOfElements();

  // Currently supporting exactly 3 dimensions (Euclidean space)...
  size_t nDims = 3;

  // Find the wound vertices of the original mesh bounding polygon
  if((geom2D->getUnsharedEdges()).get() == nullptr)
  {
    err = geom2D->findUnsharedEdges();
  }
  if(err < 0)
  {
    QString ss = QObject::tr("Error computing unshared Edges for Geometry type %1").arg(geom2D->getGeometryTypeAsString());
    setErrorCondition(err, ss);
    return;
  }

  SharedEdgeList::Pointer bEdges = geom2D->getUnsharedEdges();
  int64_t* edge = bEdges->getPointer(0);
  size_t numBoundaryEdges = bEdges->getNumberOfTuples();

  int64_t startVertex = 0;
  QVector<int64_t> woundVertices;

  startVertex = edge[0];
  woundVertices.push_back(startVertex);
  ;
  QSet<int32_t> visitedEdges;
  visitedEdges.insert(0);

  while(visitedEdges.size() < numBoundaryEdges)
  {
    for(size_t i = 1; i < numBoundaryEdges; i++)
    {
      if(!visitedEdges.contains(i))
      {
        if(edge[2 * i] == startVertex)
        {
          woundVertices.push_back(edge[2 * i + 1]);
          startVertex = edge[2 * i + 1];
          visitedEdges.insert(i);
          break;
        }
        else if(edge[2 * i + 1] == startVertex)
        {
          woundVertices.push_back(edge[2 * i]);
          startVertex = edge[2 * i];
          visitedEdges.insert(i);
          break;
        }
      }
    }
  }

  // Determine if cells lie within the original mesh geometry
  size_t index = 0;
  FloatVec3Type iPos;
  FloatVec3Type iRes;
  image->getSpacing(iRes);
  size_t iDims[3] = {0, 0, 0};
  std::tie(iDims[0], iDims[1], iDims[2]) = image->getDimensions();
  float iOrigin[3] = {0, 0, 0};
  image->getOrigin(iOrigin);

  bool insidePolygon = false;
  float xEpsilon = iRes[0] / 2.0f;
  float yEpsilon = iRes[1] / 2.0f;

  DifferentialGeometryOps diffGeomPimpl;

  for(size_t z = 0; z < iDims[2]; z++)
  {
    for(size_t y = 0; y < iDims[1]; y++)
    {
      for(size_t x = 0; x < iDims[0]; x++)
      {
        index = (z * iDims[0] * iDims[1]) + (iDims[0] * y) + x;

        // Find the (x,y,z) position of the current cell in the interpolated grid
        iPos[0] = float(x) * iRes[0] + iOrigin[0];
        iPos[1] = float(y) * iRes[1] + iOrigin[1];
        iPos[2] = float(z) * iRes[2] + iOrigin[2];

        // If the cell is on the boundary, move its position by epsilon
        if(iPos[0] == 0)
        {
          iPos[0] += xEpsilon;
        }
        else if(iPos[1] == 0)
        {
          iPos[1] += yEpsilon;
        }

        // Apply point-in-poygon
        insidePolygon = diffGeomPimpl.pointInPolygon(iPos[0], iPos[1], woundVertices, vertex);

        // If the point passes point-in-polygon, then we are inside the mesh ; else, we are outside the mesh
        if(insidePolygon)
        {
          m_InsideMesh[index] = true;
        }
      }
    }
  }

  // Set the interpolated data point to the point closest to the original mesh
  DistanceTemplate distancePimpl;

  for(size_t z = 0; z < iDims[2]; z++)
  {
    for(size_t y = 0; y < iDims[1]; y++)
    {
      for(size_t x = 0; x < iDims[0]; x++)
      {
        index = (z * iDims[0] * iDims[1]) + (iDims[0] * y) + x;

        // Find the (x,y,z) position of the current cell in the interpolated grid
        iPos[0] = float(x) * iRes[0] + iOrigin[0];
        iPos[1] = float(y) * iRes[1] + iOrigin[1];
        iPos[2] = float(z) * iRes[2] + iOrigin[2];

        // Find the index of the closest input cell centroid to the interpolated grid position
        float minDist = std::numeric_limits<float>::max();
        float dist = 0.0f;
        size_t ptrIndex = 0;

        for(size_t i = 0; i < numElements; i++)
        {
          dist = distancePimpl.getDistance<float, float, float>(cellCentroids + (nDims * i), iPos, nDims, 0);
          if(dist < minDist)
          {
            minDist = dist;
            ptrIndex = i;
          }
        }

        m_InterpolatedIndex[index] = ptrIndex;
      }
    }
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
template <typename inDataType>
void copyDataToInterpolatedGrid(IDataArray::Pointer inDataPtr, IDataArray::Pointer outDataPtr, DataContainer::Pointer interpolatedGrid, std::vector<bool> insideMesh,
                                const std::vector<size_t> interpolatedIndex, int outsideMeshVal)
{
  // Cast the IDataArray pointers to the correct types
  typename DataArray<inDataType>::Pointer inputDataPtr = std::dynamic_pointer_cast<DataArray<inDataType>>(inDataPtr);
  inDataType* inputData = static_cast<inDataType*>(inputDataPtr->getPointer(0));
  typename DataArray<inDataType>::Pointer interpolatedDataPtr = std::dynamic_pointer_cast<DataArray<inDataType>>(outDataPtr);
  inDataType* interpolatedData = static_cast<inDataType*>(interpolatedDataPtr->getPointer(0));

  size_t iDims[3] = {0, 0, 0};
  ImageGeom::Pointer image = interpolatedGrid->getGeometryAs<ImageGeom>();
  std::tie(iDims[0], iDims[1], iDims[2]) = image->getDimensions();

  size_t nComps = inDataPtr->getNumberOfComponents();
  size_t index = 0;
  size_t tmpIndex = 0;
  size_t ptrIndex = 0;

  // Copy the data from the original mesh to the regular grid using the pre-calculated interpolated indices
  for(size_t z = 0; z < iDims[2]; z++)
  {
    for(size_t y = 0; y < iDims[1]; y++)
    {
      for(size_t x = 0; x < iDims[0]; x++)
      {
        for(size_t d = 0; d < nComps; d++)
        {
          index = (z * iDims[0] * iDims[1]) + (iDims[0] * y) + x;
          tmpIndex = nComps * ((z * iDims[0] * iDims[1]) + (iDims[0] * y) + x) + d;
          ptrIndex = nComps * interpolatedIndex[index] + d;

          if(insideMesh[index])
          {
            interpolatedData[tmpIndex] = static_cast<inDataType>(inputData[ptrIndex]);
          }
          else
          {
            interpolatedData[tmpIndex] = static_cast<inDataType>(outsideMeshVal);
          }
        }
      }
    }
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
template <typename DataType> void findVertexAverage(AbstractFilter* filter, IDataArray::Pointer inDataPtr, DataArray<float>::Pointer outDataPtr, DataContainer::Pointer m, bool weightAverage)
{
  typename DataArray<DataType>::Pointer inputDataPtr = std::dynamic_pointer_cast<DataArray<DataType>>(inDataPtr);
  IGeometry2D::Pointer geom2D = m->getGeometryAs<IGeometry2D>();
  int32_t err = 0;
  IGeometry::Type geomType = geom2D->getGeometryType();

  if(weightAverage)
  {
    if(geom2D->getElementCentroids().get() == nullptr)
    {
      err = geom2D->findElementCentroids();
    }
    if(err < 0)
    {
      QString ss = QObject::tr("Error computing Element centroids for Geometry type %1").arg(geom2D->getGeometryTypeAsString());
      filter->setErrorCondition(err, ss);
      return;
    }
    switch(geomType)
    {
    case IGeometry::Type::Triangle:
    {
      TriangleGeom::Pointer tris = std::dynamic_pointer_cast<TriangleGeom>(geom2D);
      GeometryHelpers::Generic::WeightedAverageVertexArrayValues<int64_t, DataType>(tris->getTriangles(), tris->getVertices(), tris->getElementCentroids(), inputDataPtr, outDataPtr);
      break;
    }
    case IGeometry::Type::Quad:
    {
      QuadGeom::Pointer quads = std::dynamic_pointer_cast<QuadGeom>(geom2D);
      GeometryHelpers::Generic::WeightedAverageVertexArrayValues<int64_t, DataType>(quads->getQuads(), quads->getVertices(), quads->getElementCentroids(), inputDataPtr, outDataPtr);
      break;
    }
    case IGeometry::Type::Image:
      break;
    case IGeometry::Type::Vertex:
      break;
    case IGeometry::Type::Edge:
      break;
    case IGeometry::Type::Tetrahedral:
      break;
    case IGeometry::Type::Unknown:
      break;
    case IGeometry::Type::Any:
      break;
    default:
      break;
    }
  }
  else
  {
    switch(geomType)
    {
    case IGeometry::Type::Triangle:
    {
      TriangleGeom::Pointer tris = std::dynamic_pointer_cast<TriangleGeom>(geom2D);
      GeometryHelpers::Generic::AverageVertexArrayValues<int64_t, DataType>(tris->getTriangles(), inputDataPtr, outDataPtr);
      break;
    }
    case IGeometry::Type::Quad:
    {
      QuadGeom::Pointer quads = std::dynamic_pointer_cast<QuadGeom>(geom2D);
      GeometryHelpers::Generic::AverageVertexArrayValues<int64_t, DataType>(quads->getQuads(), inputDataPtr, outDataPtr);
      break;
    }
    case IGeometry::Type::Image:
      break;
    case IGeometry::Type::Vertex:
      break;
    case IGeometry::Type::Edge:
      break;
    case IGeometry::Type::Tetrahedral:
      break;
    case IGeometry::Type::Unknown:
      break;
    case IGeometry::Type::Any:
      break;
    default:
      break;
    }
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void InterpolateMeshToRegularGrid::initialize()
{
  m_AttrMatList.clear();
  m_AttrArrayMap.clear();
  m_MeshMinExtents.clear();
  m_MeshMaxExtents.clear();
  m_InsideMesh.clear();
  m_InterpolatedIndex.clear();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void InterpolateMeshToRegularGrid::dataCheck()
{
  clearErrorCode();
  clearWarningCode();

  IGeometry2D::Pointer geom2D = getDataContainerArray()->getPrereqGeometryFromDataContainer<IGeometry2D, AbstractFilter>(this, getSelectedDataContainerName());

  if(getErrorCode() < 0)
  {
    return;
  }

  uint32_t spatialDim = geom2D->getSpatialDimensionality();

  if(spatialDim != 2)
  {
    QString ss = QObject::tr("The spatial dimensionality of the Geometry is %1, but it must be 2").arg(spatialDim);
    setErrorCondition(-11000, ss);
  }

  if(m_ScaleFactorNumCells < 1)
  {
    setErrorCondition(-11004, "The scale factor for the total number of Elements must be greater than zero");
    return;
  }

  DataContainer::Pointer interpolatedDC = getDataContainerArray()->createNonPrereqDataContainer<AbstractFilter>(this, getInterpolatedDataContainerName());
  if(getErrorCode() < 0)
  {
    return;
  }
  ImageGeom::Pointer image = ImageGeom::CreateGeometry(SIMPL::Geometry::ImageGeometry);
  interpolatedDC->setGeometry(image);

  // Get the list of all attribute matrices in the input data container and add them to the regular grid data container
  DataContainer::Pointer m = getDataContainerArray()->getDataContainer(getSelectedDataContainerName());
  m_AttrMatList = m->getAttributeMatrixNames();
  QVector<size_t> tDims(1, 0);
  QList<QString> tempDataArrayList;
  DataArrayPath tempPath;
  AttributeMatrix::Type tempAttrMatType = AttributeMatrix::Type::Vertex;

  // All vertex/edge/face/cell data in the original mesh will be interpolated to the regular grid's cells
  // Feature/ensemble attribute matrices will remain unchanged and are deep copied to the new data container below
  // Create the attribute matrix where all the interpolated data will be stored
  interpolatedDC->createAndAddAttributeMatrix(tDims, getInterpolatedAttributeMatrixName(), AttributeMatrix::Type::Cell);

  // Loop through all the attribute matrices in the original data container
  // If we are in a vertex/edge/face/cell attribute matrix, create data arrays for all in the new interpolated data attribute matrix
  // Else, we are in a feature/ensemble attribute matrix, and just deep copy it into the new data container
  if(!m_SelectedDataContainerName.isEmpty())
  {
    for(QList<QString>::Iterator it = m_AttrMatList.begin(); it != m_AttrMatList.end(); ++it)
    {
      AttributeMatrix::Pointer tmpAttrMat = m->getPrereqAttributeMatrix(this, *it, -301);
      if(getErrorCode() >= 0)
      {
        tempAttrMatType = tmpAttrMat->getType();
        if(tempAttrMatType > AttributeMatrix::Type::Cell)
        {
          AttributeMatrix::Pointer attrMat = tmpAttrMat->deepCopy(getInPreflight());
          interpolatedDC->addOrReplaceAttributeMatrix(attrMat);
        }
        else
        {
          tempDataArrayList = tmpAttrMat->getAttributeArrayNames();
          m_AttrArrayMap.insert(*it, tempDataArrayList);
          for(QList<QString>::iterator jt = tempDataArrayList.begin(); jt != tempDataArrayList.end(); ++jt)
          {
            tempPath.update(getInterpolatedDataContainerName(), getInterpolatedAttributeMatrixName(), *jt);
            IDataArray::Pointer tmpDataArray = tmpAttrMat->getPrereqIDataArray<IDataArray, AbstractFilter>(this, *jt, -90002);
            if(getErrorCode() >= 0)
            {
              QVector<size_t> cDims = tmpDataArray->getComponentDimensions();
              if(tempAttrMatType == AttributeMatrix::Type::Vertex)
              {
                getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<float>, AbstractFilter, float>(this, tempPath, 0, cDims);
              }
              else
              {
                TemplateHelpers::CreateNonPrereqArrayFromArrayType()(this, tempPath, cDims, tmpDataArray);
              }
            }
          }
        }
      }
    }
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void InterpolateMeshToRegularGrid::preflight()
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
void InterpolateMeshToRegularGrid::execute()
{
  clearErrorCode();
  clearWarningCode();
  dataCheck();
  if(getErrorCode() < 0)
  {
    return;
  }

  DataContainer::Pointer m = getDataContainerArray()->getDataContainer(getSelectedDataContainerName());
  DataContainer::Pointer interpolatedDC = getDataContainerArray()->getDataContainer(getInterpolatedDataContainerName());
  IGeometry2D::Pointer geom2D = m->getGeometryAs<IGeometry2D>();
  ImageGeom::Pointer image = interpolatedDC->getGeometryAs<ImageGeom>();
  size_t numElements = geom2D->getNumberOfElements();

  // Create and properly scale the structured rectilinear grid
  createRegularGrid();

  // Set up internal variables
  size_t iDims[3] = {0, 0, 0};
  std::tie(iDims[0], iDims[1], iDims[2]) = image->getDimensions();
  m_InsideMesh.resize(iDims[0] * iDims[1] * iDims[2], false);
  m_InterpolatedIndex.resize(iDims[0] * iDims[1] * iDims[2], 0);

  // Determine interpolation
  interpolateMeshPointsToGrid(m, interpolatedDC);

  // Move the data on the mesh to the regular grid
  AttributeMatrix::Type tempAttrMatType = AttributeMatrix::Type::Vertex;

  for(QMap<QString, QList<QString>>::iterator it = m_AttrArrayMap.begin(); it != m_AttrArrayMap.end(); ++it)
  {
    tempAttrMatType = m->getAttributeMatrix(it.key())->getType();
    for(size_t i = 0; i < it.value().size(); i++)
    {
      IDataArray::Pointer tmpInPtr = m->getAttributeMatrix(it.key())->getAttributeArray(it.value()[i]);
      IDataArray::Pointer tmpOutPtr = interpolatedDC->getAttributeMatrix(getInterpolatedAttributeMatrixName())->getAttributeArray(it.value()[i]);

      // If we are in a vertex attribute matrix, we know the return data type will be float due to the cell averaging
      if(tempAttrMatType == AttributeMatrix::Type::Vertex)
      {
        QVector<size_t> cDims = tmpInPtr->getComponentDimensions();
        DataArray<float>::Pointer tmpVertInPtr = DataArray<float>::CreateArray(numElements, cDims, "_INTERNAL_USE_ONLY_tmpVertArray", true);

        EXECUTE_FUNCTION_TEMPLATE(this, findVertexAverage, tmpInPtr, this, tmpInPtr, tmpVertInPtr, m, true);

        copyDataToInterpolatedGrid<float>(tmpVertInPtr, tmpOutPtr, interpolatedDC, m_InsideMesh, m_InterpolatedIndex, getOutsideMeshIdentifier());
      }
      // We are not in a vertex attribute matrix, so execute the function template as normal
      else
      {
        EXECUTE_FUNCTION_TEMPLATE(this, copyDataToInterpolatedGrid, tmpInPtr, tmpInPtr, tmpOutPtr, interpolatedDC, m_InsideMesh, m_InterpolatedIndex, getOutsideMeshIdentifier());
      }
    }
  }

  /* Let the GUI know we are done with this filter */
  notifyStatusMessage("Complete");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AbstractFilter::Pointer InterpolateMeshToRegularGrid::newFilterInstance(bool copyFilterParameters) const
{
  InterpolateMeshToRegularGrid::Pointer filter = InterpolateMeshToRegularGrid::New();
  if(copyFilterParameters)
  {
    copyFilterParameterInstanceVariables(filter.get());
  }
  return filter;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString InterpolateMeshToRegularGrid::getCompiledLibraryName() const
{
  return DREAM3DReviewConstants::DREAM3DReviewBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString InterpolateMeshToRegularGrid::getBrandingString() const
{
  return "DREAM3DReview Plugin";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString InterpolateMeshToRegularGrid::getFilterVersion() const
{
  QString version;
  QTextStream vStream(&version);
  vStream << DREAM3DReview::Version::Major() << "." << DREAM3DReview::Version::Minor() << "." << DREAM3DReview::Version::Patch();
  return version;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString InterpolateMeshToRegularGrid::getGroupName() const
{
  return SIMPL::FilterGroups::SamplingFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString InterpolateMeshToRegularGrid::getSubGroupName() const
{
  return DREAM3DReviewConstants::FilterSubGroups::InterpolationFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString InterpolateMeshToRegularGrid::getHumanLabel() const
{
  return "Interpolate Mesh to Regular Grid";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QUuid InterpolateMeshToRegularGrid::getUuid()
{
  return QUuid("{2c1961d5-fbc5-5c90-9758-57c8852bd06e}");
}