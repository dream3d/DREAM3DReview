#include "DelaunayTriangulation.h"

#include <array>
#include <unordered_map>

#include "SIMPLib/Common/Constants.h"
#include "SIMPLib/DataContainers/DataContainer.h"
#include "SIMPLib/DataContainers/DataContainerArray.h"
#include "SIMPLib/FilterParameters/DataArraySelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/DataContainerSelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/DoubleFilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedBooleanFilterParameter.h"
#include "SIMPLib/FilterParameters/SeparatorFilterParameter.h"
#include "SIMPLib/FilterParameters/StringFilterParameter.h"
#include "SIMPLib/Geometry/EdgeGeom.h"
#include "SIMPLib/Geometry/IGeometry2D.h"
#include "SIMPLib/Geometry/IGeometry3D.h"
#include "SIMPLib/Geometry/IGeometryGrid.h"
#include "SIMPLib/Geometry/TriangleGeom.h"
#include "SIMPLib/Geometry/VertexGeom.h"

#include "DREAM3DReview/DREAM3DReviewConstants.h"
#include "DREAM3DReview/DREAM3DReviewFilters/util/Delaunay2D.h"
#include "DREAM3DReview/DREAM3DReviewVersion.h"

namespace
{
template <class T>
inline void hashCombine(size_t& seed, const T& obj)
{
  std::hash<T> hasher;
  seed ^= hasher(obj) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

using Vertex = std::array<float, 3>;

struct VertexHasher
{
  size_t operator()(const Vertex& vert) const
  {
    size_t hash = std::hash<float>()(vert[0]);
    hashCombine(hash, vert[1]);
    hashCombine(hash, vert[2]);
    return hash;
  }
};

using VertexMap = std::unordered_map<Vertex, int64_t, VertexHasher>;

void mergeTriangleGeometries(std::vector<TriangleGeom::Pointer>& triangles)
{
  int64_t numVerts = 0;
  int64_t numTris = 0;
  for(auto&& triangle : triangles)
  {
    numVerts += triangle->getNumberOfVertices();
    numTris += triangle->getNumberOfTris();
  }

  VertexMap vertexMap;
  std::vector<Vertex> mergedVerts;
  mergedVerts.reserve(numVerts);
  int64_t vertCounter = 0;
  for(auto&& triangle : triangles)
  {
    float* vertPtr = triangle->getVertexPointer(0);
    int64_t localNumVerts = triangle->getNumberOfVertices();
    for(int64_t i = 0; i < localNumVerts; i++)
    {
      Vertex vert = {vertPtr[3 * i + 0], vertPtr[3 * i + 1], vertPtr[3 * i + 2]};
      auto iter = vertexMap.find(vert);
      if(iter == vertexMap.end())
      {
        mergedVerts.push_back(vert);
        vertexMap[vert] = vertCounter;
        vertCounter++;
      }
    }
  }

  SharedVertexList::Pointer mergedVertices = TriangleGeom::CreateSharedVertexList(mergedVerts.size());
  float* mergedVertsPtr = mergedVertices->getPointer(0);
  for(size_t i = 0; i < mergedVertices->getNumberOfTuples(); i++)
  {
    mergedVertsPtr[3 * i + 0] = mergedVerts[i][0];
    mergedVertsPtr[3 * i + 1] = mergedVerts[i][1];
    mergedVertsPtr[3 * i + 2] = mergedVerts[i][2];
  }

  mergedVerts.clear();
  mergedVerts.shrink_to_fit();

  TriangleGeom::Pointer mergedTriangle = TriangleGeom::CreateGeometry(numTris, mergedVertices, SIMPL::Geometry::TriangleGeometry);
  size_t* mergedTrisPtr = mergedTriangle->getTriPointer(0);
  size_t triCounter = 0;
  for(auto&& triangle : triangles)
  {
    size_t localNumTris = triangle->getNumberOfTris();
    float vertCoords[9];
    for(size_t i = 0; i < localNumTris; i++)
    {
      triangle->getVertCoordsAtTri(i, vertCoords, vertCoords + 3, vertCoords + 6);
      for(size_t j = 0; j < 3; j++)
      {
        Vertex vert = {vertCoords[3 * j + 0], vertCoords[3 * j + 1], vertCoords[3 * j + 2]};
        mergedTrisPtr[3 * triCounter + j] = vertexMap[vert];
      }
      triCounter++;
    }
  }

  triangles.resize(1);
  triangles[0] = mergedTriangle;
}
}; // namespace

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
DelaunayTriangulation::DelaunayTriangulation()
: m_InputGeometry("")
, m_TriangleDataContainerName(SIMPL::Defaults::TriangleDataContainerName)
, m_VertexAttributeMatrixName(SIMPL::Defaults::VertexAttributeMatrixName)
, m_FaceAttributeMatrixName(SIMPL::Defaults::FaceAttributeMatrixName)
, m_Offset(1.0)
, m_Tolerance(0.00001)
, m_TriangulateByFeature(false)
, m_FeatureIdsArrayPath("", "", "")
, m_FeatureIds(nullptr)
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
DelaunayTriangulation::~DelaunayTriangulation() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void DelaunayTriangulation::setupFilterParameters()
{
  FilterParameterVectorType parameters;
  parameters.push_back(SIMPL_NEW_DOUBLE_FP("Offset", Offset, FilterParameter::Parameter, DelaunayTriangulation));
  parameters.push_back(SIMPL_NEW_DOUBLE_FP("Tolerance", Tolerance, FilterParameter::Parameter, DelaunayTriangulation));
  QStringList linkedProps = {"FeatureIdsArrayPath"};
  parameters.push_back(SIMPL_NEW_LINKED_BOOL_FP("Triangulate by Feature", TriangulateByFeature, FilterParameter::Parameter, DelaunayTriangulation, linkedProps));
  DataContainerSelectionFilterParameter::RequirementType dcReq;
  IGeometry::Types geomTypes = {IGeometry::Type::Vertex, IGeometry::Type::Edge, IGeometry::Type::Triangle, IGeometry::Type::Quad, IGeometry::Type::Tetrahedral};
  parameters.push_back(SIMPL_NEW_DC_SELECTION_FP("Input Vertices", InputGeometry, FilterParameter::RequiredArray, DelaunayTriangulation, dcReq));
  DataArraySelectionFilterParameter::RequirementType daReq = DataArraySelectionFilterParameter::CreateRequirement(SIMPL::TypeNames::Int32, 1, AttributeMatrix::Type::Any, IGeometry::Type::Any);
  daReq.dcGeometryTypes = geomTypes;
  parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Feature Ids", FeatureIdsArrayPath, FilterParameter::RequiredArray, DelaunayTriangulation, daReq));
  parameters.push_back(SIMPL_NEW_STRING_FP("Triangle Data Container", TriangleDataContainerName, FilterParameter::CreatedArray, DelaunayTriangulation));
  parameters.push_back(SeparatorFilterParameter::New("Vertex Data", FilterParameter::CreatedArray));
  parameters.push_back(SIMPL_NEW_STRING_FP("Vertex Attribute Matrix", VertexAttributeMatrixName, FilterParameter::CreatedArray, DelaunayTriangulation));
  parameters.push_back(SeparatorFilterParameter::New("Face Data", FilterParameter::CreatedArray));
  parameters.push_back(SIMPL_NEW_STRING_FP("Face Attribute Matrix", FaceAttributeMatrixName, FilterParameter::CreatedArray, DelaunayTriangulation));
  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void DelaunayTriangulation::initialize()
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void DelaunayTriangulation::dataCheck()
{
  clearErrorCode();
  clearWarningCode();

  if(getOffset() < 0)
  {
    QString ss = QObject::tr("Offset must be 0 or greater");
    setErrorCondition(-1, ss);
  }

  if(getTolerance() < 0)
  {
    QString ss = QObject::tr("Tolerance must be 0 or greater");
    setErrorCondition(-1, ss);
  }

  IGeometry::Pointer geometry = getDataContainerArray()->getPrereqGeometryFromDataContainer<IGeometry>(this, getInputGeometry());

  if(getErrorCode() < 0)
  {
    return;
  }

  if(std::dynamic_pointer_cast<IGeometryGrid>(geometry))
  {
    QString ss = QObject::tr("Selected Geometry type (%1) does not contain explicit Vertices for triangulation\n"
                             "Geometry must be of type Vertex, Edge, Triangle, Quadrilateral, or Tetrahedral")
                     .arg(geometry->getGeometryTypeAsString());
    setErrorCondition(-1, ss);
  }

  if(getTriangulateByFeature())
  {
    std::vector<size_t> cDims(1, 1);
    m_FeatureIdsPtr = getDataContainerArray()->getPrereqArrayFromPath<Int32ArrayType>(this, getFeatureIdsArrayPath(), cDims);
    if(m_FeatureIdsPtr.lock())
    {
      m_FeatureIds = m_FeatureIdsPtr.lock()->getPointer(0);
    }

    SharedVertexList::Pointer vertices;

    if(VertexGeom::Pointer vertex = std::dynamic_pointer_cast<VertexGeom>(geometry))
    {
      vertices = vertex->getVertices();
    }
    else if(EdgeGeom::Pointer edge = std::dynamic_pointer_cast<EdgeGeom>(geometry))
    {
      vertices = edge->getVertices();
    }
    else if(IGeometry2D::Pointer geometry2d = std::dynamic_pointer_cast<IGeometry2D>(geometry))
    {
      vertices = geometry2d->getVertices();
    }
    else if(IGeometry3D::Pointer geometry3d = std::dynamic_pointer_cast<IGeometry3D>(geometry))
    {
      vertices = geometry3d->getVertices();
    }
    else
    {
      QString ss = QObject::tr("Selected Geometry type (%1) does not contain explicit Vertices for triangulation\n"
                               "Geometry must be of type Vertex, Edge, Triangle, Quadrilateral, or Tetrahedral")
                       .arg(geometry->getGeometryTypeAsString());
      setErrorCondition(-1, ss);
    }

    if(getErrorCode() < 0)
    {
      return;
    }

    QVector<IDataArray::Pointer> dataArrays = {m_FeatureIdsPtr.lock(), vertices};
    getDataContainerArray()->validateNumberOfTuples(this, dataArrays);
  }

  SharedVertexList::Pointer vertices = TriangleGeom::CreateSharedVertexList(0, !getInPreflight());
  TriangleGeom::Pointer triangle = TriangleGeom::CreateGeometry(0, vertices, SIMPL::Geometry::TriangleGeometry, !getInPreflight());

  DataContainer::Pointer dc = getDataContainerArray()->createNonPrereqDataContainer(this, getTriangleDataContainerName());

  if(getErrorCode() < 0)
  {
    return;
  }

  dc->setGeometry(triangle);

  std::vector<size_t> tDims(1, 0);
  dc->createNonPrereqAttributeMatrix(this, getVertexAttributeMatrixName(), tDims, AttributeMatrix::Type::Vertex);
  dc->createNonPrereqAttributeMatrix(this, getFaceAttributeMatrixName(), tDims, AttributeMatrix::Type::Face);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void DelaunayTriangulation::execute()
{
  dataCheck();
  if(getErrorCode() < 0)
  {
    return;
  }

  IGeometry::Pointer geometry = getDataContainerArray()->getDataContainer(m_InputGeometry)->getGeometry();
  SharedVertexList::Pointer vertices;
  int64_t numVerts = 0;

  if(VertexGeom::Pointer vertex = std::dynamic_pointer_cast<VertexGeom>(geometry))
  {
    vertices = vertex->getVertices();
    numVerts = vertex->getNumberOfVertices();
  }
  else if(EdgeGeom::Pointer edge = std::dynamic_pointer_cast<EdgeGeom>(geometry))
  {
    vertices = edge->getVertices();
    numVerts = edge->getNumberOfVertices();
  }
  else if(IGeometry2D::Pointer geometry2d = std::dynamic_pointer_cast<IGeometry2D>(geometry))
  {
    vertices = geometry2d->getVertices();
    numVerts = geometry2d->getNumberOfVertices();
  }
  else if(IGeometry3D::Pointer geometry3d = std::dynamic_pointer_cast<IGeometry3D>(geometry))
  {
    vertices = geometry3d->getVertices();
    numVerts = geometry3d->getNumberOfVertices();
  }
  else
  {
    QString ss = QObject::tr("Selected Geometry type (%1) does not contain explicit Vertices for triangulation\n"
                             "Geometry must be of type Vertex, Edge, Triangle, Quadrilateral, or Tetrahedral")
                     .arg(geometry->getGeometryTypeAsString());
    setErrorCondition(-1, ss);
    return;
  }

  float* vertPtr = vertices->getPointer(0);
  std::vector<TriMesh::VertexCoordList> vertexLists;
  int32_t numFeatures = std::numeric_limits<int32_t>::min();

  if(m_TriangulateByFeature)
  {
    for(size_t i = 0; i < m_FeatureIdsPtr.lock()->getNumberOfTuples(); i++)
    {
      if(m_FeatureIds[i] > numFeatures)
      {
        numFeatures = m_FeatureIds[i];
      }
    }

    if(numFeatures <= 0)
    {
      QString ss = QObject::tr("No Features were found in the supplied Feature Ids array");
      setErrorCondition(-1, ss);
      return;
    }

    vertexLists.resize(numFeatures + 1);
    for(size_t i = 0; i < m_FeatureIdsPtr.lock()->getNumberOfTuples(); i++)
    {
      vertexLists[m_FeatureIds[i]].push_back({vertPtr[3 * i + 0], vertPtr[3 * i + 1], vertPtr[3 * i + 2]});
    }
  }
  else
  {
    vertexLists.resize(1);
    for(size_t i = 0; i < vertices->getNumberOfTuples(); i++)
    {
      vertexLists[0].push_back({vertPtr[3 * i + 0], vertPtr[3 * i + 1], vertPtr[3 * i + 2]});
    }
  }

  size_t counter = 1;
  std::vector<TriangleGeom::Pointer> triangles(vertexLists.size());

  for(auto&& vertexList : vertexLists)
  {
    Delaunay2D::Pointer delaunay = Delaunay2D::New(vertexList, m_Offset, m_Tolerance, 0.0, this);
    delaunay->setMessagePrefix(getHumanLabel());
    QString ss = QObject::tr("Performing Delaunay Triangulation");
    if(m_TriangulateByFeature)
    {
      QString featureProgress = QObject::tr(" || Feature %1 of %2").arg(counter).arg(numFeatures);
      ss += featureProgress;
    }
    delaunay->setMessageTitle(ss);
    triangles[counter - 1] = delaunay->triangulate();
    counter++;
  }

  DataContainer::Pointer dc = getDataContainerArray()->getDataContainer(m_TriangleDataContainerName);

  if(m_TriangulateByFeature)
  {
    mergeTriangleGeometries(triangles);
  }

  dc->setGeometry(triangles[0]);
  AttributeMatrix::Pointer attrMat = getDataContainerArray()->getDataContainer(m_TriangleDataContainerName)->getAttributeMatrix(m_VertexAttributeMatrixName);
  std::vector<size_t> tDims(1, triangles[0]->getNumberOfVertices());
  attrMat->resizeAttributeArrays(tDims);
  attrMat = getDataContainerArray()->getDataContainer(m_TriangleDataContainerName)->getAttributeMatrix(m_FaceAttributeMatrixName);
  tDims[0] = triangles[0]->getNumberOfTris();
  attrMat->resizeAttributeArrays(tDims);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AbstractFilter::Pointer DelaunayTriangulation::newFilterInstance(bool copyFilterParameters) const
{
  DelaunayTriangulation::Pointer filter = DelaunayTriangulation::New();
  if(copyFilterParameters)
  {
    copyFilterParameterInstanceVariables(filter.get());
  }
  return filter;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString DelaunayTriangulation::getCompiledLibraryName() const
{
  return DREAM3DReviewConstants::DREAM3DReviewBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString DelaunayTriangulation::getBrandingString() const
{
  return "DREAM3DReview";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString DelaunayTriangulation::getFilterVersion() const
{
  QString version;
  QTextStream vStream(&version);
  vStream << DREAM3DReview::Version::Major() << "." << DREAM3DReview::Version::Minor() << "." << DREAM3DReview::Version::Patch();
  return version;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString DelaunayTriangulation::getGroupName() const
{
  return SIMPL::FilterGroups::SurfaceMeshingFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString DelaunayTriangulation::getSubGroupName() const
{
  return SIMPL::FilterSubGroups::GenerationFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString DelaunayTriangulation::getHumanLabel() const
{
  return "Delaunay Triangulation";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QUuid DelaunayTriangulation::getUuid() const
{
  return QUuid("{10b319ca-6b2f-5118-9f4e-d19ed481cd1f}");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
DelaunayTriangulation::Pointer DelaunayTriangulation::NullPointer()
{
  return Pointer(static_cast<Self*>(nullptr));
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
std::shared_ptr<DelaunayTriangulation> DelaunayTriangulation::New()
{
  struct make_shared_enabler : public DelaunayTriangulation
  {
  };
  std::shared_ptr<make_shared_enabler> val = std::make_shared<make_shared_enabler>();
  val->setupFilterParameters();
  return val;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString DelaunayTriangulation::getNameOfClass() const
{
  return QString("DelaunayTriangulation");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString DelaunayTriangulation::ClassName()
{
  return QString("DelaunayTriangulation");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void DelaunayTriangulation::setInputGeometry(const DataArrayPath& value)
{
  m_InputGeometry = value;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
DataArrayPath DelaunayTriangulation::getInputGeometry() const
{
  return m_InputGeometry;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void DelaunayTriangulation::setTriangleDataContainerName(const QString& value)
{
  m_TriangleDataContainerName = value;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString DelaunayTriangulation::getTriangleDataContainerName() const
{
  return m_TriangleDataContainerName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void DelaunayTriangulation::setVertexAttributeMatrixName(const QString& value)
{
  m_VertexAttributeMatrixName = value;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString DelaunayTriangulation::getVertexAttributeMatrixName() const
{
  return m_VertexAttributeMatrixName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void DelaunayTriangulation::setFaceAttributeMatrixName(const QString& value)
{
  m_FaceAttributeMatrixName = value;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString DelaunayTriangulation::getFaceAttributeMatrixName() const
{
  return m_FaceAttributeMatrixName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void DelaunayTriangulation::setOffset(const double& value)
{
  m_Offset = value;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
double DelaunayTriangulation::getOffset() const
{
  return m_Offset;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void DelaunayTriangulation::setTolerance(const double& value)
{
  m_Tolerance = value;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
double DelaunayTriangulation::getTolerance() const
{
  return m_Tolerance;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void DelaunayTriangulation::setTriangulateByFeature(const bool& value)
{
  m_TriangulateByFeature = value;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
bool DelaunayTriangulation::getTriangulateByFeature() const
{
  return m_TriangulateByFeature;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void DelaunayTriangulation::setFeatureIdsArrayPath(const DataArrayPath& value)
{
  m_FeatureIdsArrayPath = value;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
DataArrayPath DelaunayTriangulation::getFeatureIdsArrayPath() const
{
  return m_FeatureIdsArrayPath;
}
