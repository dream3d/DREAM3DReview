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

#include "FindVertexToTriangleDistances.h"

#include <cassert>
#include <cmath>

#include "SIMPLib/Common/Constants.h"
#include "SIMPLib/DataContainers/DataContainer.h"
#include "SIMPLib/DataContainers/DataContainerArray.h"
#include "SIMPLib/FilterParameters/DataArrayCreationFilterParameter.h"
#include "SIMPLib/FilterParameters/DataArraySelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/DataContainerSelectionFilterParameter.h"
#include "SIMPLib/Geometry/TriangleGeom.h"
#include "SIMPLib/Geometry/VertexGeom.h"
#include "SIMPLib/Math/GeometryMath.h"
#include "SIMPLib/Utilities/ParallelDataAlgorithm.h"

#include "DREAM3DReview/DREAM3DReviewConstants.h"
#include "DREAM3DReview/DREAM3DReviewVersion.h"

#ifdef SIMPL_USE_PARALLEL_ALGORITHMS
#include <tbb/blocked_range.h>
#endif

class FindVertexToTriangleDistancesImpl
{
public:
  FindVertexToTriangleDistancesImpl(FindVertexToTriangleDistances* filter, std::vector<std::vector<size_t>>& tris, std::vector<std::vector<float>>& verts, float* sourceVerts, float* distances,
                                    int32_t* closestTri, size_t numTris)
  : m_Filter(filter)
  , m_Tris(tris)
  , m_Verts(verts)
  , m_SourceVerts(sourceVerts)
  , m_Distances(distances)
  , m_ClosestTri(closestTri)
  , m_NumTris(numTris)
  {
  }
  virtual ~FindVertexToTriangleDistancesImpl() = default;

  void compute(int64_t start, int64_t end) const
  {
    int64_t counter = 0;
    int64_t totalElements = end - start;
    int64_t progIncrement = static_cast<int64_t>(totalElements / 100);

    for(int64_t v = start; v < end; v++)
    {
      for(int64_t t = 0; t < m_NumTris; t++)
      {
        if(m_Filter->getCancel())
        {
          return;
        }
        int64_t p = m_Tris[t][0];
        int64_t q = m_Tris[t][1];
        int64_t r = m_Tris[t][2];
        std::vector<float> gx = {m_SourceVerts[3 * v + 0], m_SourceVerts[3 * v + 1], m_SourceVerts[3 * v + 2]};
        float d = m_Filter->point_triangle_distance(gx, m_Verts[p], m_Verts[q], m_Verts[r], t);
        if(std::abs(d) < std::abs(m_Distances[v]))
        {
          m_Distances[v] = d;
          m_ClosestTri[v] = t;
        }
      }

      if(counter > progIncrement)
      {
        m_Filter->sendThreadSafeProgressMessage(counter);
        counter = 0;
      }
      counter++;
    }
  }

  void operator()(const SIMPLRange& range) const
  {
    compute(range.min(), range.max());
  }

private:
  FindVertexToTriangleDistances* m_Filter;
  std::vector<std::vector<size_t>>& m_Tris;
  std::vector<std::vector<float>>& m_Verts;
  float* m_SourceVerts;
  float* m_Distances;
  int32_t* m_ClosestTri;
  int64_t m_NumTris;
};

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
FindVertexToTriangleDistances::FindVertexToTriangleDistances()
{
  initialize();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
FindVertexToTriangleDistances::~FindVertexToTriangleDistances() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FindVertexToTriangleDistances::initialize()
{
  clearErrorCode();
  setCancel(false);
  m_ProgressCounter = 0;
  m_TotalElements = 0;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FindVertexToTriangleDistances::setupFilterParameters()
{
  FilterParameterVectorType parameters;
  DataContainerSelectionFilterParameter::RequirementType dcsReq;
  IGeometry::Types geomTypes = {IGeometry::Type::Vertex};
  dcsReq.dcGeometryTypes = geomTypes;
  parameters.push_back(SIMPL_NEW_DC_SELECTION_FP("Source Vertex Geometry", VertexDataContainer, FilterParameter::Category::RequiredArray, FindVertexToTriangleDistances, dcsReq));
  geomTypes = {IGeometry::Type::Triangle};
  dcsReq.dcGeometryTypes = geomTypes;
  parameters.push_back(SIMPL_NEW_DC_SELECTION_FP("Target Triangle Geometry", TriangleDataContainer, FilterParameter::Category::RequiredArray, FindVertexToTriangleDistances, dcsReq));
  DataArraySelectionFilterParameter::RequirementType dasReq = DataArraySelectionFilterParameter::CreateRequirement(SIMPL::TypeNames::Double, 3, AttributeMatrix::Type::Face, IGeometry::Type::Triangle);
  parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Triangle Normals", TriangleNormalsArrayPath, FilterParameter::Category::RequiredArray, FindVertexToTriangleDistances, dasReq));
  DataArrayCreationFilterParameter::RequirementType dacReq;
  parameters.push_back(SIMPL_NEW_DA_CREATION_FP("Distances", DistancesArrayPath, FilterParameter::Category::CreatedArray, FindVertexToTriangleDistances, dacReq));
  parameters.push_back(SIMPL_NEW_DA_CREATION_FP("Closest Triangle Ids", ClosestTriangleIdArrayPath, FilterParameter::Category::CreatedArray, FindVertexToTriangleDistances, dacReq));
  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FindVertexToTriangleDistances::dataCheck()
{
  clearErrorCode();
  clearWarningCode();
  getDataContainerArray()->getPrereqGeometryFromDataContainer<VertexGeom>(this, getVertexDataContainer());
  getDataContainerArray()->getPrereqGeometryFromDataContainer<TriangleGeom>(this, getTriangleDataContainer());

  std::vector<size_t> cDims(1, 3);

  m_NormalsPtr = getDataContainerArray()->getPrereqArrayFromPath<DoubleArrayType>(this, getTriangleNormalsArrayPath(), cDims);
  if(m_NormalsPtr.lock())
  {
    m_Normals = m_NormalsPtr.lock()->getPointer(0);
  }

  cDims[0] = 1;

  m_DistancesPtr = getDataContainerArray()->createNonPrereqArrayFromPath<FloatArrayType>(this, getDistancesArrayPath(), 0, cDims);
  if(m_DistancesPtr.lock())
  {
    m_Distances = m_DistancesPtr.lock()->getPointer(0);
  }

  m_ClosestTriangleIdsPtr = getDataContainerArray()->createNonPrereqArrayFromPath<Int32ArrayType>(this, getClosestTriangleIdArrayPath(), 0, cDims);
  if(m_ClosestTriangleIdsPtr.lock())
  {
    m_ClosestTriangleIds = m_ClosestTriangleIdsPtr.lock()->getPointer(0);
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
float FindVertexToTriangleDistances::distance(const std::vector<float>& vec1, const std::vector<float>& vec2)
{
  assert(vec1.size() == vec2.size());

  float dist = 0.0f;

  for(size_t i = 0; i < vec1.size(); i++)
  {
    dist += (vec1[i] - vec2[i]) * (vec1[i] - vec2[i]);
  }

  return std::sqrt(dist);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
float FindVertexToTriangleDistances::point_segment_distance(const std::vector<float>& x0, const std::vector<float>& x1, const std::vector<float>& x2)
{
  std::vector<float> dx = {0.0f, 0.0f, 0.0f};
  std::transform(x2.begin(), x2.end(), x1.begin(), dx.begin(), std::minus<float>());
  double m2 = 0.0;

  for(size_t i = 0; i < dx.size(); i++)
  {
    m2 += static_cast<double>(dx[i] * dx[i]);
  }

  std::vector<float> x2minx0 = {0.0f, 0.0f, 0.0f};
  std::transform(x2.begin(), x2.end(), x0.begin(), x2minx0.begin(), std::minus<float>());
  double dotProduct = 0.0;
  for(size_t i = 0; i < x2minx0.size(); i++)
  {
    dotProduct += static_cast<double>(x2minx0[i] * dx[i]);
  }

  float s12 = static_cast<float>(dotProduct / m2);
  if(s12 < 0.0f)
  {
    s12 = 0.0f;
  }
  else if(s12 > 1.0f)
  {
    s12 = 1.0f;
  }

  std::vector<float> multVal1 = {x1[0] * s12, x1[1] * s12, x1[2] * s12};
  std::vector<float> mutlVal2 = {x2[0] * (1 - s12), x2[1] * (1 - s12), x2[2] * (1 - s12)};

  std::vector<float> vecSum = {0.0f, 0.0f, 0.0f};
  for(size_t i = 0; i < x1.size(); i++)
  {
    vecSum[i] = multVal1[i] + mutlVal2[i];
  }

  return distance(x0, vecSum);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
float FindVertexToTriangleDistances::point_triangle_distance(const std::vector<float>& x0, const std::vector<float>& x1, const std::vector<float>& x2, const std::vector<float>& x3,
                                                             const int64_t& triangle)
{
  float dist = 0.0f;

  std::vector<float> x13 = {0.0f, 0.0f, 0.0f};
  std::vector<float> x23 = {0.0f, 0.0f, 0.0f};
  std::vector<float> x03 = {0.0f, 0.0f, 0.0f};
  std::transform(x1.begin(), x1.end(), x3.begin(), x13.begin(), std::minus<float>());
  std::transform(x2.begin(), x2.end(), x3.begin(), x23.begin(), std::minus<float>());
  std::transform(x0.begin(), x0.end(), x3.begin(), x03.begin(), std::minus<float>());

  float m13 = 0.0f;
  float m23 = 0.0f;
  for(size_t i = 0; i < x13.size(); i++)
  {
    m13 += (x13[i] * x13[i]);
  }
  for(size_t i = 0; i < x23.size(); i++)
  {
    m23 += (x23[i] * x23[i]);
  }
  float d = 0.0;
  for(size_t i = 0; i < x13.size(); i++)
  {
    d += (x13[i] * x23[i]);
  }
  float invdet = 1.0f / std::max(m13 * m23 - d * d, 1e-30f);
  float a = 0.0f;
  float b = 0.0f;
  for(size_t i = 0; i < x13.size(); i++)
  {
    a += (x13[i] * x03[i]);
    b += (x23[i] * x03[i]);
  }

  float w23 = invdet * (m23 * a - d * b);
  float w31 = invdet * (m13 * b - d * a);
  float w12 = 1 - w23 - w31;

  if(w23 >= 0.0f && w31 >= 0.0f && w12 >= 0.0f)
  {
    std::vector<float> tmpVec = {0.0f, 0.0f, 0.0f};
    for(size_t i = 0; i < tmpVec.size(); i++)
    {
      tmpVec[i] = (w23 * x1[i]) + (w31 * x2[i]) + (w12 * x3[i]);
    }

    dist = distance(x0, tmpVec);
  }
  else
  {
    if(w23 > 0)
    {
      dist = std::min(point_segment_distance(x0, x1, x2), point_segment_distance(x0, x1, x3));
    }
    else if(w31 > 0)
    {
      dist = std::min(point_segment_distance(x0, x1, x2), point_segment_distance(x0, x2, x3));
    }
    else
    {
      dist = std::min(point_segment_distance(x0, x1, x3), point_segment_distance(x0, x2, x3));
    }
  }

  float normal[3] = {static_cast<float>(m_Normals[3 * triangle + 0]), static_cast<float>(m_Normals[3 * triangle + 1]), static_cast<float>(m_Normals[3 * triangle + 2])};

  float cosTheta = GeometryMath::CosThetaBetweenVectors(normal, x0.data());

  if(cosTheta < 0.0f)
  {
    dist *= -1.0f;
  }

  return dist;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FindVertexToTriangleDistances::sendThreadSafeProgressMessage(int64_t counter)
{
  m_Mutex.lock();

  m_ProgressCounter += counter;
  int64_t progressInt = static_cast<int64_t>((static_cast<float>(m_ProgressCounter) / m_TotalElements) * 100.0f);

  int64_t progIncrement = m_TotalElements / 100;
  int64_t prog = 1;

  if(m_ProgressCounter > prog && m_LastProgressInt != progressInt)
  {
    QString ss = QObject::tr("Working on Vertex to Triangle Distances || %1% Completed").arg(progressInt);
    notifyStatusMessage(ss);
    prog += progIncrement;
  }

  m_LastProgressInt = progressInt;

  m_Mutex.unlock();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FindVertexToTriangleDistances::execute()
{
  initialize();
  dataCheck();
  if(getErrorCode() < 0)
  {
    return;
  }

  VertexGeom::Pointer sourceGeom = getDataContainerArray()->getDataContainer(m_VertexDataContainer)->getGeometryAs<VertexGeom>();
  TriangleGeom::Pointer targetGeom = getDataContainerArray()->getDataContainer(m_TriangleDataContainer)->getGeometryAs<TriangleGeom>();
  size_t numSourceVerts = sourceGeom->getNumberOfVertices();
  size_t numTris = targetGeom->getNumberOfTris();
  size_t numVerts = targetGeom->getNumberOfVertices();
  float* sourceVerts = sourceGeom->getVertexPointer(0);
  size_t* triangles = targetGeom->getTriPointer(0);
  float* vertices = targetGeom->getVertexPointer(0);

  m_TotalElements = numSourceVerts;

  std::vector<std::vector<size_t>> tmpTris;
  std::vector<std::vector<float>> tmpVerts;

  for(size_t i = 0; i < numTris; i++)
  {
    std::vector<size_t> tmpTri = {triangles[3 * i + 0], triangles[3 * i + 1], triangles[3 * i + 2]};
    tmpTris.push_back(tmpTri);
  }

  for(size_t i = 0; i < numVerts; i++)
  {
    std::vector<float> tmpVert = {vertices[3 * i + 0], vertices[3 * i + 1], vertices[3 * i + 2]};
    tmpVerts.push_back(tmpVert);
  }

  m_DistancesPtr.lock()->initializeWithValue(std::numeric_limits<float>::max());
  m_Distances = m_DistancesPtr.lock()->getPointer(0);
  m_ClosestTriangleIdsPtr.lock()->initializeWithValue(-1);
  m_ClosestTriangleIds = m_ClosestTriangleIdsPtr.lock()->getPointer(0);

  // Allow data-based parallelization
  ParallelDataAlgorithm dataAlg;
  dataAlg.setRange(0, numSourceVerts);
  dataAlg.execute(FindVertexToTriangleDistancesImpl(this, tmpTris, tmpVerts, sourceVerts, m_Distances, m_ClosestTriangleIds, numTris));
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AbstractFilter::Pointer FindVertexToTriangleDistances::newFilterInstance(bool copyFilterParameters) const
{
  FindVertexToTriangleDistances::Pointer filter = FindVertexToTriangleDistances::New();
  if(copyFilterParameters)
  {
    copyFilterParameterInstanceVariables(filter.get());
  }
  return filter;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString FindVertexToTriangleDistances::getCompiledLibraryName() const
{
  return DREAM3DReviewConstants::DREAM3DReviewBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString FindVertexToTriangleDistances::getBrandingString() const
{
  return "DREAM3DReview";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString FindVertexToTriangleDistances::getFilterVersion() const
{
  QString version;
  QTextStream vStream(&version);
  vStream << DREAM3DReview::Version::Major() << "." << DREAM3DReview::Version::Minor() << "." << DREAM3DReview::Version::Patch();
  return version;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString FindVertexToTriangleDistances::getGroupName() const
{
  return SIMPL::FilterGroups::SamplingFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString FindVertexToTriangleDistances::getSubGroupName() const
{
  return SIMPL::FilterSubGroups::SpatialFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString FindVertexToTriangleDistances::getHumanLabel() const
{
  return "Find Vertex to Triangle Distances";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QUuid FindVertexToTriangleDistances::getUuid() const
{
  return QUuid("{fcdde553-36b4-5731-bc88-fc499806cb4e}");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
FindVertexToTriangleDistances::Pointer FindVertexToTriangleDistances::NullPointer()
{
  return Pointer(static_cast<Self*>(nullptr));
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
std::shared_ptr<FindVertexToTriangleDistances> FindVertexToTriangleDistances::New()
{
  struct make_shared_enabler : public FindVertexToTriangleDistances
  {
  };
  std::shared_ptr<make_shared_enabler> val = std::make_shared<make_shared_enabler>();
  val->setupFilterParameters();
  return val;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString FindVertexToTriangleDistances::getNameOfClass() const
{
  return QString("FindVertexToTriangleDistances");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString FindVertexToTriangleDistances::ClassName()
{
  return QString("FindVertexToTriangleDistances");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FindVertexToTriangleDistances::setVertexDataContainer(const DataArrayPath& value)
{
  m_VertexDataContainer = value;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
DataArrayPath FindVertexToTriangleDistances::getVertexDataContainer() const
{
  return m_VertexDataContainer;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FindVertexToTriangleDistances::setTriangleDataContainer(const DataArrayPath& value)
{
  m_TriangleDataContainer = value;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
DataArrayPath FindVertexToTriangleDistances::getTriangleDataContainer() const
{
  return m_TriangleDataContainer;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FindVertexToTriangleDistances::setTriangleNormalsArrayPath(const DataArrayPath& value)
{
  m_TriangleNormalsArrayPath = value;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
DataArrayPath FindVertexToTriangleDistances::getTriangleNormalsArrayPath() const
{
  return m_TriangleNormalsArrayPath;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FindVertexToTriangleDistances::setDistancesArrayPath(const DataArrayPath& value)
{
  m_DistancesArrayPath = value;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
DataArrayPath FindVertexToTriangleDistances::getDistancesArrayPath() const
{
  return m_DistancesArrayPath;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FindVertexToTriangleDistances::setClosestTriangleIdArrayPath(const DataArrayPath& value)
{
  m_ClosestTriangleIdArrayPath = value;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
DataArrayPath FindVertexToTriangleDistances::getClosestTriangleIdArrayPath() const
{
  return m_ClosestTriangleIdArrayPath;
}
