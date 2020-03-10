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
#include <memory>

#include "ApproximatePointCloudHull.h"

#include <QtCore/QTextStream>

#include "SIMPLib/Common/Constants.h"

#include "SIMPLib/FilterParameters/DataContainerSelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/IntFilterParameter.h"
#include "SIMPLib/FilterParameters/StringFilterParameter.h"
#include "SIMPLib/Geometry/VertexGeom.h"
#include "SIMPLib/DataContainers/DataContainerArray.h"
#include "SIMPLib/DataContainers/DataContainer.h"

#include "DREAM3DReview/DREAM3DReviewConstants.h"
#include "DREAM3DReview/DREAM3DReviewVersion.h"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
ApproximatePointCloudHull::ApproximatePointCloudHull()
: m_VertexDataContainerName(SIMPL::Defaults::VertexDataContainerName)
, m_HullDataContainerName("HullDataContainer")
, m_NumberOfEmptyNeighbors(1)
{
  m_GridResolution[0] = 1.0f;
  m_GridResolution[1] = 1.0f;
  m_GridResolution[2] = 1.0f;

  initialize();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
ApproximatePointCloudHull::~ApproximatePointCloudHull() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ApproximatePointCloudHull::initialize()
{
  clearErrorCode();
  clearWarningCode();
  setCancel(false);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ApproximatePointCloudHull::setupFilterParameters()
{
  FilterParameterVectorType parameters;
  parameters.push_back(SIMPL_NEW_FLOAT_VEC3_FP("Grid Resolution", GridResolution, FilterParameter::Parameter, ApproximatePointCloudHull));
  parameters.push_back(SIMPL_NEW_INTEGER_FP("Minimum Number of Emtpy Neighbors", NumberOfEmptyNeighbors, FilterParameter::Parameter, ApproximatePointCloudHull));
  DataContainerSelectionFilterParameter::RequirementType dcsReq;
  IGeometry::Types req = {IGeometry::Type::Vertex};
  dcsReq.dcGeometryTypes = req;
  parameters.push_back(SIMPL_NEW_DC_SELECTION_FP("Vertex Geometry", VertexDataContainerName, FilterParameter::RequiredArray, ApproximatePointCloudHull, dcsReq));
  parameters.push_back(SIMPL_NEW_STRING_FP("Hull Vertex Geoemetry", HullDataContainerName, FilterParameter::CreatedArray, ApproximatePointCloudHull));
  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ApproximatePointCloudHull::dataCheck()
{
  clearErrorCode();
  clearWarningCode();

  if(getGridResolution()[0] <= 0.0f || getGridResolution()[1] <= 0.0f || getGridResolution()[2] <= 0.0f)
  {
    QString ss = QObject::tr("Grid resolutions must be greater than zero");
    setErrorCondition(-11001, ss);
  }

  if(getNumberOfEmptyNeighbors() < 0)
  {
    QString ss = QObject::tr("Minimum number of empty neighbors must be positive");
    setErrorCondition(-11001, ss);
  }

  getDataContainerArray()->getPrereqGeometryFromDataContainer<VertexGeom>(this, getVertexDataContainerName());

  DataContainer::Pointer dc = getDataContainerArray()->createNonPrereqDataContainer(this, getHullDataContainerName());

  if(getErrorCode() < 0)
  {
    return;
  }

  VertexGeom::Pointer vertex = VertexGeom::CreateGeometry(0, SIMPL::Geometry::VertexGeometry, !getInPreflight());
  dc->setGeometry(vertex);
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
bool validNeighbor(int64_t dims[3], int64_t neighborhood[78], size_t index, int64_t x, int64_t y, int64_t z)
{
  int64_t modX = x + neighborhood[3 * index + 0];
  int64_t modY = y + neighborhood[3 * index + 1];
  int64_t modZ = z + neighborhood[3 * index + 2];

  return (modX >= 0 && modX < dims[0]) && (modY >= 0 && modY < dims[1]) && (modZ >= 0 && modZ < dims[2]);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ApproximatePointCloudHull::execute()
{
  initialize();
  dataCheck();
  if(getErrorCode() < 0)
  {
    return;
  }

  float inverseResolution[3] = {1.0f / m_GridResolution[0], 1.0f / m_GridResolution[1], 1.0f / m_GridResolution[2]};

  VertexGeom::Pointer source = getDataContainerArray()->getDataContainer(m_VertexDataContainerName)->getGeometryAs<VertexGeom>();
  float* verts = source->getVertexPointer(0);

  m_SamplingGrid = ImageGeom::CreateGeometry(SIMPL::Geometry::ImageGeometry);
  m_SamplingGrid->setSpacing(m_GridResolution[0], m_GridResolution[1], m_GridResolution[2]);

  int64_t numVerts = source->getNumberOfVertices();
  float* vertex = source->getVertexPointer(0);

  for(size_t i = 0; i < 3; i++)
  {
    m_MeshMaxExtents.push_back(std::numeric_limits<float>::lowest());
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

  for(auto i = 0; i < 3; i++)
  {
    m_MeshMinExtents[i] -= (inverseResolution[i] / 2.0f);
    m_MeshMaxExtents[i] += (inverseResolution[i] / 2.0f);
  }

  int64_t bboxMin[3] = {0, 0, 0};
  int64_t bboxMax[3] = {0, 0, 0};

  bboxMin[0] = static_cast<int64_t>(std::floor(m_MeshMinExtents[0] * inverseResolution[0]));
  bboxMin[1] = static_cast<int64_t>(std::floor(m_MeshMinExtents[1] * inverseResolution[1]));
  bboxMin[2] = static_cast<int64_t>(std::floor(m_MeshMinExtents[2] * inverseResolution[2]));

  bboxMax[0] = static_cast<int64_t>(std::floor(m_MeshMaxExtents[0] * inverseResolution[0]));
  bboxMax[1] = static_cast<int64_t>(std::floor(m_MeshMaxExtents[1] * inverseResolution[1]));
  bboxMax[2] = static_cast<int64_t>(std::floor(m_MeshMaxExtents[2] * inverseResolution[2]));

  int64_t dims[3] = {bboxMax[0] - bboxMin[0] + 1, bboxMax[1] - bboxMin[1] + 1, bboxMax[2] - bboxMin[2] + 1};
  m_SamplingGrid->setDimensions(dims[0], dims[1], dims[2]);

  int64_t multiplier[3] = {1, 
    static_cast<int64_t>(m_SamplingGrid->getXPoints()), 
    static_cast<int64_t>(m_SamplingGrid->getXPoints() * m_SamplingGrid->getYPoints())};
  std::vector<std::vector<int64_t>> vertsInVoxels(m_SamplingGrid->getNumberOfElements());

  int64_t progIncrement = numVerts / 100;
  int64_t prog = 1;
  int64_t progressInt = 0;

  for(int64_t v = 0; v < numVerts; v++)
  {
    int64_t i = static_cast<int64_t>(std::floor(verts[3 * v + 0] * inverseResolution[0]) - static_cast<float>(bboxMin[0]));
    int64_t j = static_cast<int64_t>(std::floor(verts[3 * v + 1] * inverseResolution[1]) - static_cast<float>(bboxMin[1]));
    int64_t k = static_cast<int64_t>(std::floor(verts[3 * v + 2] * inverseResolution[2]) - static_cast<float>(bboxMin[2]));
    int64_t index = i * multiplier[0] + j * multiplier[1] + k * multiplier[2];
    vertsInVoxels[index].push_back(v);

    if(v > prog)
    {
      progressInt = static_cast<int64_t>((static_cast<float>(v) / numVerts) * 100.0f);
      QString ss = QObject::tr("Mapping Vertices to Voxels || %1% Complete").arg(progressInt);
      notifyStatusMessage(ss);
      prog = prog + progIncrement;
    }
  }

  std::vector<float> tmpVerts;
  int64_t neighborhood[78] = {1,  0, 0,  -1, 0, 0, 0, 1, 0,  0, -1, 0, 0, 0,  1,  0, 0, -1, 1, 1, 0,  -1, 1,  0, 1, -1, 0,  -1, -1, 0, 1,  0, 1,  1,  0,  -1, -1, 0,  1,
                              -1, 0, -1, 0,  1, 1, 0, 1, -1, 0, -1, 1, 0, -1, -1, 1, 1, 1,  1, 1, -1, 1,  -1, 1, 1, -1, -1, -1, 1,  1, -1, 1, -1, -1, -1, 1,  -1, -1, -1};

  progIncrement = (dims[0] * dims[1] * dims[2]) / 100;
  prog = 1;
  progressInt = 0;
  int64_t counter = 0;
  int64_t vertCounter = 0;
  float xAvg = 0.0f;
  float yAvg = 0.0f;
  float zAvg = 0.0f;

  for(int64_t z = 0; z < dims[2]; z++)
  {
    for(int64_t y = 0; y < dims[1]; y++)
    {
      for(int64_t x = 0; x < dims[0]; x++)
      {
        size_t index = (z * dims[1] * dims[0]) + (y * dims[0]) + x;
        if(vertsInVoxels[index].empty())
        {
          counter++;
          continue;
        }

        size_t emtpyNeighbors = 0;

        for(size_t n = 0; n < 26; n++)
        {
          if(validNeighbor(dims, neighborhood, n, x, y, z))
          {
            size_t neighborIndex = ((z + neighborhood[3 * n + 2]) * dims[1] * dims[0]) + ((y + neighborhood[3 * n + 1]) * dims[0]) + (x + neighborhood[3 * n + 0]);
            if(vertsInVoxels[neighborIndex].empty())
            {
              emtpyNeighbors++;
            }
          }
        }

        if(emtpyNeighbors > m_NumberOfEmptyNeighbors)
        {
          for(auto vert : vertsInVoxels[index])
          {
            vertCounter++;
            xAvg += verts[3 * vert + 0];
            yAvg += verts[3 * vert + 1];
            zAvg += verts[3 * vert + 2];
          }
          xAvg /= static_cast<float>(vertCounter);
          yAvg /= static_cast<float>(vertCounter);
          zAvg /= static_cast<float>(vertCounter);
          tmpVerts.push_back(xAvg);
          tmpVerts.push_back(yAvg);
          tmpVerts.push_back(zAvg);
          vertCounter = 0;
          xAvg = 0.0f;
          yAvg = 0.0f;
          zAvg = 0.0f;
        }

        if(counter > prog)
        {
          progressInt = static_cast<int64_t>((static_cast<float>(counter) / (dims[0] * dims[1] * dims[2])) * 100.0f);
          QString ss = QObject::tr("Trimming Interior Voxels || %1% Complete").arg(progressInt);
          notifyStatusMessage(ss);
          prog = prog + progIncrement;
        }
        counter++;
      }
    }
  }

  VertexGeom::Pointer hull = getDataContainerArray()->getDataContainer(m_HullDataContainerName)->getGeometryAs<VertexGeom>();
  hull->resizeVertexList(tmpVerts.size() / 3);
  float* hullVerts = hull->getVertexPointer(0);
  std::memcpy(hullVerts, tmpVerts.data(), hull->getNumberOfVertices() * 3 * sizeof(float));

  notifyStatusMessage("Complete");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AbstractFilter::Pointer ApproximatePointCloudHull::newFilterInstance(bool copyFilterParameters) const
{
  ApproximatePointCloudHull::Pointer filter = ApproximatePointCloudHull::New();
  if(copyFilterParameters)
  {
    copyFilterParameterInstanceVariables(filter.get());
  }
  return filter;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ApproximatePointCloudHull::getCompiledLibraryName() const
{
  return DREAM3DReviewConstants::DREAM3DReviewBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ApproximatePointCloudHull::getBrandingString() const
{
  return "DREAM3DReview";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ApproximatePointCloudHull::getFilterVersion() const
{
  QString version;
  QTextStream vStream(&version);
  vStream << DREAM3DReview::Version::Major() << "." << DREAM3DReview::Version::Minor() << "." << DREAM3DReview::Version::Patch();
  return version;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ApproximatePointCloudHull::getGroupName() const
{
  return DREAM3DReviewConstants::FilterSubGroups::PointCloudFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ApproximatePointCloudHull::getSubGroupName() const
{
  return SIMPL::FilterSubGroups::GeometryFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ApproximatePointCloudHull::getHumanLabel() const
{
  return "Approximate Point Cloud Hull";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QUuid ApproximatePointCloudHull::getUuid() const
{
  return QUuid("{fab669ad-66c6-5a39-bdb7-fc47b94311ed}");
}

// -----------------------------------------------------------------------------
ApproximatePointCloudHull::Pointer ApproximatePointCloudHull::NullPointer()
{
  return Pointer(static_cast<Self*>(nullptr));
}

// -----------------------------------------------------------------------------
std::shared_ptr<ApproximatePointCloudHull> ApproximatePointCloudHull::New()
{
  struct make_shared_enabler : public ApproximatePointCloudHull
  {
  };
  std::shared_ptr<make_shared_enabler> val = std::make_shared<make_shared_enabler>();
  val->setupFilterParameters();
  return val;
}

// -----------------------------------------------------------------------------
QString ApproximatePointCloudHull::getNameOfClass() const
{
  return QString("ApproximatePointCloudHull");
}

// -----------------------------------------------------------------------------
QString ApproximatePointCloudHull::ClassName()
{
  return QString("ApproximatePointCloudHull");
}

// -----------------------------------------------------------------------------
void ApproximatePointCloudHull::setVertexDataContainerName(const DataArrayPath& value)
{
  m_VertexDataContainerName = value;
}

// -----------------------------------------------------------------------------
DataArrayPath ApproximatePointCloudHull::getVertexDataContainerName() const
{
  return m_VertexDataContainerName;
}

// -----------------------------------------------------------------------------
void ApproximatePointCloudHull::setHullDataContainerName(const QString& value)
{
  m_HullDataContainerName = value;
}

// -----------------------------------------------------------------------------
QString ApproximatePointCloudHull::getHullDataContainerName() const
{
  return m_HullDataContainerName;
}

// -----------------------------------------------------------------------------
void ApproximatePointCloudHull::setGridResolution(const FloatVec3Type& value)
{
  m_GridResolution = value;
}

// -----------------------------------------------------------------------------
FloatVec3Type ApproximatePointCloudHull::getGridResolution() const
{
  return m_GridResolution;
}

// -----------------------------------------------------------------------------
void ApproximatePointCloudHull::setNumberOfEmptyNeighbors(int value)
{
  m_NumberOfEmptyNeighbors = value;
}

// -----------------------------------------------------------------------------
int ApproximatePointCloudHull::getNumberOfEmptyNeighbors() const
{
  return m_NumberOfEmptyNeighbors;
}
