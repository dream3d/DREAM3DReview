#include "DecimatePointCloud.h"

#include "SIMPLib/Common/Constants.h"
#include "SIMPLib/DataContainers/DataContainer.h"
#include "SIMPLib/DataContainers/DataContainerArray.h"
#include "SIMPLib/FilterParameters/AbstractFilterParametersReader.h"
#include "SIMPLib/FilterParameters/AttributeMatrixSelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/FloatFilterParameter.h"
#include "SIMPLib/FilterParameters/SeparatorFilterParameter.h"
#include "SIMPLib/Geometry/VertexGeom.h"

#include "DREAM3DReview/DREAM3DReviewConstants.h"
#include "DREAM3DReview/DREAM3DReviewVersion.h"


#include <chrono>
#include <random>

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
DecimatePointCloud::DecimatePointCloud() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
DecimatePointCloud::~DecimatePointCloud() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void DecimatePointCloud::setupFilterParameters()
{
  FilterParameterVectorType parameters;
  parameters.push_back(SIMPL_NEW_FLOAT_FP("Percentage to Keep (0 < N < 1.0)", DecimationFreq, FilterParameter::Category::Parameter, DecimatePointCloud));
  parameters.push_back(SeparatorFilterParameter::Create("Vertex Data", FilterParameter::Category::RequiredArray));
  {
    AttributeMatrixSelectionFilterParameter::RequirementType req;
    IGeometry::Types geomTypes = {IGeometry::Type::Vertex};
    AttributeMatrix::Types amTypes = {AttributeMatrix::Type::Vertex};
    req.dcGeometryTypes = geomTypes;
    req.amTypes = amTypes;
    parameters.push_back(SIMPL_NEW_AM_SELECTION_FP("Vertex Attribute Matrix", VertexAttrMatPath, FilterParameter::Category::RequiredArray, DecimatePointCloud, req));
  }
  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void DecimatePointCloud::readFilterParameters(AbstractFilterParametersReader* reader, int index)
{
  reader->openFilterGroup(this, index);
  setVertexAttrMatPath(reader->readDataArrayPath("VertexAttrMatPath", getVertexAttrMatPath()));
  setDecimationFreq(reader->readValue("DecimationFreq", getDecimationFreq()));
  reader->closeFilterGroup();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void DecimatePointCloud::initialize()
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void DecimatePointCloud::dataCheck()
{
  clearErrorCode();
  clearWarningCode();

  VertexGeom::Pointer vertices = getDataContainerArray()->getPrereqGeometryFromDataContainer<VertexGeom>(this, getVertexAttrMatPath().getDataContainerName());
  if(getErrorCode() < 0)
  {
    return;
  }

  if(getDecimationFreq() >= 1.0F || getDecimationFreq() <= 0.0F)
  {
    QString ss = QObject::tr("Decimation frequency must be between 0 and 1");
    setErrorCondition(-11000, ss);
  }

  int64_t numVerts = vertices->getNumberOfVertices();
  
  AttributeMatrix::Pointer attrMat = getDataContainerArray()->getPrereqAttributeMatrixFromPath(this, getVertexAttrMatPath(), -301);
  if(getErrorCode() < 0)
  {
    return;
  }

  int64_t totalTuples = static_cast<int64_t>(attrMat->getNumberOfTuples());
  if(totalTuples != numVerts)
  {
    QString ss = QObject::tr("The selected Vertex Attribute Matrix does not have the same number of tuples as the number of vertices in the Vertex Geometry");
    setErrorCondition(-11002, ss);
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void DecimatePointCloud::execute()
{
  dataCheck();
  if(getErrorCode() < 0)
  {
    return;
  }

  VertexGeom::Pointer vertices = getDataContainerArray()->getDataContainer(getVertexAttrMatPath().getDataContainerName())->getGeometryAs<VertexGeom>();
  AttributeMatrix::Pointer attrMat = getDataContainerArray()->getAttributeMatrix(getVertexAttrMatPath());

  int64_t numVerts = vertices->getNumberOfVertices();
  float* vertex = vertices->getVertexPointer(0);

  int64_t progIncrement = numVerts / 100;
  int64_t prog = 1;
  int64_t progressInt = 0;
  std::vector<size_t> removeList;
  size_t currentSlot = 0;

    std::mt19937_64::result_type seed = static_cast<std::mt19937_64::result_type>(std::chrono::steady_clock::now().time_since_epoch().count());
    std::mt19937_64 generator(seed); // Standard mersenne_twister_engine seeded with milliseconds
    std::uniform_real_distribution<> distribution(0.0f, 1.0);

  for(int64_t i = 0; i < numVerts; i++)
  {
    float randNum = distribution(generator);
    if(randNum >= m_DecimationFreq)
    {
      removeList.push_back(i);
    }
    else
    {
      vertex[3 * currentSlot] = vertex[3 * i];
      vertex[3 * currentSlot + 1] = vertex[3 * i + 1];
      vertex[3 * currentSlot + 2] = vertex[3 * i + 2];
      currentSlot++;
    }

    if(i > prog)
    {
      progressInt = static_cast<int64_t>((static_cast<float>(i) / numVerts) * 100.0f);
      QString ss = QObject::tr("Decimating Point Cloud || %1% Complete").arg(progressInt);
      notifyStatusMessage(ss);
      prog = prog + progIncrement;
    }
    if (getCancel())
    {
      break;
    }
  }

  QString ss = QObject::tr("Resizing Vertex Geometry...");
  notifyStatusMessage(ss);

  vertices->resizeVertexList(currentSlot);

  if(!removeList.empty())
  {
    QList<QString> headers = attrMat->getAttributeArrayNames();
    for(QList<QString>::iterator iter = headers.begin(); iter != headers.end(); ++iter)
    {
      if(getCancel())
      {
        return;
      }
      IDataArray::Pointer p = attrMat->getAttributeArray(*iter);
      QString type = p->getTypeAsString();
      if(type.compare("NeighborList<T>") == 0)
      {
        attrMat->removeAttributeArray(*iter);
      }
      else
      {
        p->eraseTuples(removeList);
      }
    }
    std::vector<size_t> tDims(1, (numVerts - removeList.size()));
    attrMat->setTupleDimensions(tDims);
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AbstractFilter::Pointer DecimatePointCloud::newFilterInstance(bool copyFilterParameters) const
{
  DecimatePointCloud::Pointer filter = DecimatePointCloud::New();
  if(copyFilterParameters)
  {
    copyFilterParameterInstanceVariables(filter.get());
  }
  return filter;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString DecimatePointCloud::getCompiledLibraryName() const
{
  return DREAM3DReviewConstants::DREAM3DReviewBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString DecimatePointCloud::getBrandingString() const
{
  return "DREAM3DReview";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString DecimatePointCloud::getFilterVersion() const
{
  QString version;
  QTextStream vStream(&version);
  vStream << DREAM3DReview::Version::Major() << "." << DREAM3DReview::Version::Minor() << "." << DREAM3DReview::Version::Patch();
  return version;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString DecimatePointCloud::getGroupName() const
{
  return SIMPL::FilterGroups::SamplingFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString DecimatePointCloud::getSubGroupName() const
{
  return DREAM3DReviewConstants::FilterSubGroups::PointCloudFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString DecimatePointCloud::getHumanLabel() const
{
  return "Decimate Point Cloud";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QUuid DecimatePointCloud::getUuid() const
{
  return QUuid("{5cbd9d8e-e2eb-59e7-be63-6ab9deeed8d2}");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
DecimatePointCloud::Pointer DecimatePointCloud::NullPointer()
{
  return Pointer(static_cast<Self*>(nullptr));
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
std::shared_ptr<DecimatePointCloud> DecimatePointCloud::New()
{
  struct make_shared_enabler : public DecimatePointCloud
  {
  };
  std::shared_ptr<make_shared_enabler> val = std::make_shared<make_shared_enabler>();
  val->setupFilterParameters();
  return val;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString DecimatePointCloud::getNameOfClass() const
{
  return QString("DecimatePointCloud");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString DecimatePointCloud::ClassName()
{
  return QString("DecimatePointCloud");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void DecimatePointCloud::setVertexAttrMatPath(const DataArrayPath& value)
{
  m_VertexAttrMatPath = value;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
DataArrayPath DecimatePointCloud::getVertexAttrMatPath() const
{
  return m_VertexAttrMatPath;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void DecimatePointCloud::setDecimationFreq(const float& value)
{
  m_DecimationFreq = value;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
float DecimatePointCloud::getDecimationFreq() const
{
  return m_DecimationFreq;
}
