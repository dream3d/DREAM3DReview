#include "IterativeClosestPoint.h"

#include <Eigen/Geometry>

#include "SIMPLib/Common/Constants.h"
#include "SIMPLib/DataContainers/DataContainer.h"
#include "SIMPLib/DataContainers/DataContainerArray.h"
#include "SIMPLib/FilterParameters/BooleanFilterParameter.h"
#include "SIMPLib/FilterParameters/DataContainerSelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/IntFilterParameter.h"
#include "SIMPLib/FilterParameters/StringFilterParameter.h"
#include "SIMPLib/Geometry/VertexGeom.h"

#include "DREAM3DReview/DREAM3DReviewConstants.h"
#include "DREAM3DReview/DREAM3DReviewVersion.h"

#include "DREAM3DReview/DREAM3DReviewFilters/util/nanoflann.hpp"

namespace
{
template <typename Derived>
struct VertexGeomAdaptor
{
  const Derived& obj;

  VertexGeomAdaptor(const Derived& obj_)
  : obj(obj_)
  {
  }

  inline const Derived& derived() const
  {
    return obj;
  }

  inline size_t kdtree_get_point_count() const
  {
    return derived()->getNumberOfVertices();
  }

  inline float kdtree_get_pt(const size_t idx, const size_t dim) const
  {
    if(dim == 0)
    {
      return derived()->getVertexPointer(idx)[0];
    }
    if(dim == 1)
    {
      return derived()->getVertexPointer(idx)[1];
    }

    return derived()->getVertexPointer(idx)[2];
  }

  template <class BBOX>
  bool kdtree_get_bbox(BBOX& /*bb*/) const
  {
    return false;
  }
};
} // namespace

enum createdPathID : RenameDataPath::DataID_t
{
  AttributeMatrixID20 = 20,
  ArrayID21 = 21,
};

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
IterativeClosestPoint::IterativeClosestPoint()
{
  initialize();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
IterativeClosestPoint::~IterativeClosestPoint() = default;
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void IterativeClosestPoint::initialize()
{
  clearErrorCode();
  clearWarningCode();
  setCancel(false);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void IterativeClosestPoint::setupFilterParameters()
{
  FilterParameterVectorType parameters;
  DataContainerSelectionFilterParameter::RequirementType dcsReq;
  dcsReq.dcGeometryTypes = {IGeometry::Type::Vertex};
  parameters.push_back(SIMPL_NEW_DC_SELECTION_FP("Moving Vertex Geometry", MovingVertexGeometry, FilterParameter::Category::RequiredArray, IterativeClosestPoint, dcsReq));
  parameters.push_back(SIMPL_NEW_DC_SELECTION_FP("Target Vertex Geometry", TargetVertexGeometry, FilterParameter::Category::RequiredArray, IterativeClosestPoint, dcsReq));
  parameters.push_back(SIMPL_NEW_INTEGER_FP("Number of Iterations", Iterations, FilterParameter::Category::Parameter, IterativeClosestPoint));
  parameters.push_back(SIMPL_NEW_BOOL_FP("Apply Transform to Moving Geometry", ApplyTransform, FilterParameter::Category::Parameter, IterativeClosestPoint));
  parameters.push_back(SIMPL_NEW_STRING_FP("Transform Attribute Matrix Name", TransformAttributeMatrixName, FilterParameter::Category::CreatedArray, IterativeClosestPoint));
  parameters.push_back(SIMPL_NEW_STRING_FP("Transform Array Name", TransformArrayName, FilterParameter::Category::CreatedArray, IterativeClosestPoint));
  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void IterativeClosestPoint::dataCheck()
{
  clearErrorCode();
  clearWarningCode();

  getDataContainerArray()->getPrereqGeometryFromDataContainer<VertexGeom>(this, m_MovingVertexGeometry);
  getDataContainerArray()->getPrereqGeometryFromDataContainer<VertexGeom>(this, m_TargetVertexGeometry);

  if(getIterations() < 1)
  {
    setErrorCondition(-1, "Number if iterations must be at least 1");
  }

  DataContainer::Pointer dc = getDataContainerArray()->getPrereqDataContainer(this, m_MovingVertexGeometry);

  if(getErrorCode() < 0)
  {
    return;
  }

  AttributeMatrix::Pointer am = dc->createNonPrereqAttributeMatrix(this, DataArrayPath(m_MovingVertexGeometry.getDataContainerName(), getTransformAttributeMatrixName(), ""), {1},
                                                                   AttributeMatrix::Type::Generic, AttributeMatrixID20);

  if(getErrorCode() < 0)
  {
    return;
  }

  am->createNonPrereqArray<DataArray<float>>(this, getTransformArrayName(), 0, {4, 4}, ArrayID21);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void IterativeClosestPoint::execute()
{
  initialize();
  dataCheck();
  if(getErrorCode() < 0)
  {
    return;
  }

  VertexGeom::Pointer moving = getDataContainerArray()->getDataContainer(m_MovingVertexGeometry.getDataContainerName())->getGeometryAs<VertexGeom>();
  VertexGeom::Pointer movingCopy = std::dynamic_pointer_cast<VertexGeom>(moving->deepCopy());
  VertexGeom::Pointer target = getDataContainerArray()->getDataContainer(m_TargetVertexGeometry.getDataContainerName())->getGeometryAs<VertexGeom>();

  float* movingPtr = moving->getVertexPointer(0);
  float* movingCopyPtr = movingCopy->getVertexPointer(0);
  float* targetPtr = target->getVertexPointer(0);

  size_t numMovingVerts = moving->getNumberOfVertices();
  std::vector<size_t> cDims(1, 3);
  FloatArrayType::Pointer dynTarget = FloatArrayType::CreateArray(numMovingVerts, cDims, "tmp", true);
  dynTarget->initializeWithZeros();
  float* dynTargetPtr = dynTarget->getPointer(0);

  using Adaptor = VertexGeomAdaptor<VertexGeom::Pointer>;
  const Adaptor adaptor(target);

  notifyStatusMessage("Building kd-tree index...");

  using KDtree = nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Adaptor<float, Adaptor>, Adaptor, 3>;
  KDtree index(3, adaptor, nanoflann::KDTreeSingleIndexAdaptorParams(30));
  index.buildIndex();

  size_t iters = m_Iterations;
  const size_t nn = 1;

  typedef Eigen::Matrix<float, 3, Eigen::Dynamic, Eigen::ColMajor> PointCloud;
  typedef Eigen::Matrix<float, 4, 4, Eigen::ColMajor> UmeyamaTransform;

  UmeyamaTransform globalTransform;
  globalTransform << 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1;

  int64_t progIncrement = iters / 100;
  int64_t prog = 1;
  int64_t progressInt = 0;
  int64_t counter = 0;

  for(size_t i = 0; i < iters; i++)
  {
    if(getCancel())
    {
      return;
    }

    for(size_t j = 0; j < numMovingVerts; j++)
    {
      size_t id;
      float dist;
      nanoflann::KNNResultSet<float> results(nn);
      results.init(&id, &dist);
      index.findNeighbors(results, movingCopyPtr + (3 * j), nanoflann::SearchParams());
      dynTargetPtr[3 * j + 0] = targetPtr[3 * id + 0];
      dynTargetPtr[3 * j + 1] = targetPtr[3 * id + 1];
      dynTargetPtr[3 * j + 2] = targetPtr[3 * id + 2];
    }

    Eigen::Map<PointCloud> moving_(movingCopyPtr, 3, numMovingVerts);
    Eigen::Map<PointCloud> target_(dynTargetPtr, 3, numMovingVerts);

    UmeyamaTransform transform = Eigen::umeyama(moving_, target_, false);

    for(size_t j = 0; j < numMovingVerts; j++)
    {
      Eigen::Vector4f position(movingCopyPtr[3 * j + 0], movingCopyPtr[3 * j + 1], movingCopyPtr[3 * j + 2], 1);
      Eigen::Vector4f transformedPosition = transform * position;
      std::memcpy(movingCopyPtr + (3 * j), transformedPosition.data(), sizeof(float) * 3);
    }

    globalTransform = transform * globalTransform;

    if(counter > prog)
    {
      progressInt = static_cast<int64_t>((static_cast<float>(counter) / iters) * 100.0f);
      QString ss = QObject::tr("Performing Registration Iterations || %1% Completed").arg(progressInt);
      notifyStatusMessage(ss);
      prog = prog + progIncrement;
    }
    counter++;
  }

  float* transformPtr = getDataContainerArray()
                            ->getDataContainer(m_MovingVertexGeometry.getDataContainerName())
                            ->getAttributeMatrix(m_TransformAttributeMatrixName)
                            ->getAttributeArrayAs<DataArray<float>>(m_TransformArrayName)
                            ->getPointer(0);

  if(m_ApplyTransform)
  {
    for(size_t j = 0; j < numMovingVerts; j++)
    {
      Eigen::Vector4f position(movingPtr[3 * j + 0], movingPtr[3 * j + 1], movingPtr[3 * j + 2], 1);
      Eigen::Vector4f transformedPosition = globalTransform * position;
      std::memcpy(movingPtr + (3 * j), transformedPosition.data(), sizeof(float) * 3);
    }
  }

  globalTransform.transposeInPlace();
  std::memcpy(transformPtr, globalTransform.data(), sizeof(float) * 16);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AbstractFilter::Pointer IterativeClosestPoint::newFilterInstance(bool copyFilterParameters) const
{
  IterativeClosestPoint::Pointer filter = IterativeClosestPoint::New();
  if(copyFilterParameters)
  {
    copyFilterParameterInstanceVariables(filter.get());
  }
  return filter;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString IterativeClosestPoint::getCompiledLibraryName() const
{
  return DREAM3DReviewConstants::DREAM3DReviewBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString IterativeClosestPoint::getBrandingString() const
{
  return "DREAM3DReview";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString IterativeClosestPoint::getFilterVersion() const
{
  QString version;
  QTextStream vStream(&version);
  vStream << DREAM3DReview::Version::Major() << "." << DREAM3DReview::Version::Minor() << "." << DREAM3DReview::Version::Patch();
  return version;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString IterativeClosestPoint::getGroupName() const
{
  return SIMPL::FilterGroups::ReconstructionFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString IterativeClosestPoint::getSubGroupName() const
{
  return SIMPL::FilterSubGroups::AlignmentFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString IterativeClosestPoint::getHumanLabel() const
{
  return "Iterative Closest Point";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QUuid IterativeClosestPoint::getUuid() const
{
  return QUuid("{6c8fb24b-5b12-551c-ba6d-ae2fa7724764}");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
IterativeClosestPoint::Pointer IterativeClosestPoint::NullPointer()
{
  return Pointer(static_cast<Self*>(nullptr));
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
std::shared_ptr<IterativeClosestPoint> IterativeClosestPoint::New()
{
  struct make_shared_enabler : public IterativeClosestPoint
  {
  };
  std::shared_ptr<make_shared_enabler> val = std::make_shared<make_shared_enabler>();
  val->setupFilterParameters();
  return val;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString IterativeClosestPoint::getNameOfClass() const
{
  return QString("IterativeClosestPoint");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString IterativeClosestPoint::ClassName()
{
  return QString("IterativeClosestPoint");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void IterativeClosestPoint::setMovingVertexGeometry(const DataArrayPath& value)
{
  m_MovingVertexGeometry = value;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
DataArrayPath IterativeClosestPoint::getMovingVertexGeometry() const
{
  return m_MovingVertexGeometry;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void IterativeClosestPoint::setTargetVertexGeometry(const DataArrayPath& value)
{
  m_TargetVertexGeometry = value;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
DataArrayPath IterativeClosestPoint::getTargetVertexGeometry() const
{
  return m_TargetVertexGeometry;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void IterativeClosestPoint::setIterations(const int& value)
{
  m_Iterations = value;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int IterativeClosestPoint::getIterations() const
{
  return m_Iterations;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void IterativeClosestPoint::setApplyTransform(const bool& value)
{
  m_ApplyTransform = value;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
bool IterativeClosestPoint::getApplyTransform() const
{
  return m_ApplyTransform;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void IterativeClosestPoint::setTransformAttributeMatrixName(const QString& value)
{
  m_TransformAttributeMatrixName = value;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString IterativeClosestPoint::getTransformAttributeMatrixName() const
{
  return m_TransformAttributeMatrixName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void IterativeClosestPoint::setTransformArrayName(const QString& value)
{
  m_TransformArrayName = value;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString IterativeClosestPoint::getTransformArrayName() const
{
  return m_TransformArrayName;
}
