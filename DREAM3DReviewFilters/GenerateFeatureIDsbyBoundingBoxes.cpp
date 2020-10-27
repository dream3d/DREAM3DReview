/*
 * Your License or Copyright can go here
 */

#include "GenerateFeatureIDsbyBoundingBoxes.h"

#include <array>

#include <QtCore/QTextStream>

#include "SIMPLib/Common/Constants.h"

#include "SIMPLib/DataContainers/DataContainerArray.h"
#include "SIMPLib/FilterParameters/AttributeMatrixCreationFilterParameter.h"
#include "SIMPLib/FilterParameters/DataArrayCreationFilterParameter.h"
#include "SIMPLib/FilterParameters/DataArraySelectionFilterParameter.h"
#include "SIMPLib/Geometry/ImageGeom.h"
#include "SIMPLib/Geometry/VertexGeom.h"

#include "DREAM3DReview/DREAM3DReviewConstants.h"
#include "DREAM3DReview/DREAM3DReviewVersion.h"

namespace
{
constexpr int32_t k_AttributeMatrixTypeSelectionError = -5555;

}

/* Create Enumerations to allow the created Attribute Arrays to take part in renaming */
enum createdPathID : RenameDataPath::DataID_t
{
  AttributeMatrixID20 = 20,
};

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
GenerateFeatureIDsbyBoundingBoxes::GenerateFeatureIDsbyBoundingBoxes()
: m_FeatureIDsArrayPath("", "", "")
, m_FeatureAttributeMatrixArrayPath("", "", "")
, m_BoxCenterArrayPath("", "", "")
, m_BoxDimensionsArrayPath("", "", "")
{
  initialize();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
GenerateFeatureIDsbyBoundingBoxes::~GenerateFeatureIDsbyBoundingBoxes() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void GenerateFeatureIDsbyBoundingBoxes::initialize()
{
  clearErrorCode();
  clearWarningCode();
  setCancel(false);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void GenerateFeatureIDsbyBoundingBoxes::setupFilterParameters()
{
  FilterParameterVectorType parameters;
  DataArrayCreationFilterParameter::RequirementType dacReq;
  dacReq.amTypes = {AttributeMatrix::Type::Vertex, AttributeMatrix::Type::Cell};
  parameters.push_back(SIMPL_NEW_DA_CREATION_FP("Feature IDs", FeatureIDsArrayPath, FilterParameter::CreatedArray, GenerateFeatureIDsbyBoundingBoxes, dacReq));
  AttributeMatrixCreationFilterParameter::RequirementType amcReq;
  parameters.push_back(SIMPL_NEW_AM_CREATION_FP("Feature Attribute Matrix", FeatureAttributeMatrixArrayPath, FilterParameter::CreatedArray, GenerateFeatureIDsbyBoundingBoxes, amcReq));
  DataArraySelectionFilterParameter::RequirementType dasReq;
  parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Box Corner", BoxCenterArrayPath, FilterParameter::RequiredArray, GenerateFeatureIDsbyBoundingBoxes, dasReq));
  parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Box Dimensions", BoxDimensionsArrayPath, FilterParameter::RequiredArray, GenerateFeatureIDsbyBoundingBoxes, dasReq));
  parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Box Feature IDs Array", BoxFeatureIDsArrayPath, FilterParameter::RequiredArray, GenerateFeatureIDsbyBoundingBoxes, dasReq));
  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void GenerateFeatureIDsbyBoundingBoxes::dataCheck()
{
  clearErrorCode();
  clearWarningCode();

  DataContainer::Pointer m = getDataContainerArray()->getPrereqDataContainer(this, getFeatureAttributeMatrixArrayPath().getDataContainerName(), false);
  if(getErrorCode() < 0)
  {
    return;
  }
  DataArrayPath path(getFeatureIDsArrayPath().getDataContainerName(), getFeatureIDsArrayPath().getAttributeMatrixName(), "");
  AttributeMatrix::Pointer attrMat = getDataContainerArray()->getPrereqAttributeMatrixFromPath(this, path, -301);

  if(getErrorCode() < 0)
  {
    return;
  }

  m_DestAttributeMatrixType = AttributeMatrix::Type::Unknown;
  AttributeMatrix::Type attrMatType = attrMat->getType();

  if(attrMatType == AttributeMatrix::Type::Vertex)
  {
    m_DestAttributeMatrixType = AttributeMatrix::Type::VertexFeature;
  }
  else if(attrMatType == AttributeMatrix::Type::Cell)
  {
    m_DestAttributeMatrixType = AttributeMatrix::Type::CellFeature;
  }
  // else if (attrMatType == AttributeMatrix::Type::Edge)
  //{
  //  m_DestAttributeMatrixType = AttributeMatrix::Type::EdgeFeature;
  //}
  else
  {
    m_DestAttributeMatrixType = AttributeMatrix::Type::Unknown;
    QString ss = QObject::tr("The Attribute Matrix must have a cell, vertex geometry.");
    setErrorCondition(::k_AttributeMatrixTypeSelectionError, ss);
    return;
  }

  std::vector<size_t> cDims = {1};
  m_BoxFeatureIdsPtr = getDataContainerArray()->getPrereqArrayFromPath<DataArray<int32_t>>(this, getBoxFeatureIDsArrayPath(), cDims);
  if(nullptr != m_BoxFeatureIdsPtr.lock())
  {
    m_BoxFeatureIds = m_BoxFeatureIdsPtr.lock()->getPointer(0);
  }

  cDims[0] = 3;
  m_BoxCenterPtr = getDataContainerArray()->getPrereqArrayFromPath<DataArray<float>>(this, getBoxCenterArrayPath(), cDims);
  if(nullptr != m_BoxCenterPtr.lock())
  {
    m_BoxCenter = m_BoxCenterPtr.lock()->getPointer(0);
  }

  m_BoxDimsPtr = getDataContainerArray()->getPrereqArrayFromPath<DataArray<float>>(this, getBoxDimensionsArrayPath(), cDims);
  if(nullptr != m_BoxDimsPtr.lock())
  {
    m_BoxDims = m_BoxDimsPtr.lock()->getPointer(0);
  }

  cDims[0] = 1;
  m_FeatureIdsPtr =
      getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<int32_t>>(this, getFeatureIDsArrayPath(), 0, cDims); /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
  if(nullptr != m_FeatureIdsPtr.lock())
  {
    m_FeatureIds = m_FeatureIdsPtr.lock()->getPointer(0);
  }

  DataArrayPath path_bounding(getBoxFeatureIDsArrayPath().getDataContainerName(), getBoxFeatureIDsArrayPath().getAttributeMatrixName(), "");
  AttributeMatrix::Pointer attrMat_bounding = getDataContainerArray()->getPrereqAttributeMatrixFromPath(this, path_bounding, -301);

  if(getErrorCode() < 0)
  {
    return;
  }

  std::vector<size_t> tDims = {attrMat_bounding->getNumberOfTuples() + 1};
  m->createNonPrereqAttributeMatrix(this, getFeatureAttributeMatrixArrayPath(), tDims, m_DestAttributeMatrixType, AttributeMatrixID20);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
bool IsPointInBounds(float xmax, float xmin, float ymax, float ymin, float zmax, float zmin, float x, float y, float z)
{
  return (x < xmax) && (x > xmin) && (y < ymax) && (y > ymin) && (z < zmax) && (z > zmin);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void GenerateFeatureIDsbyBoundingBoxes::checkBoundingBoxImage()
{
  size_t totalNumFIDs = m_BoxFeatureIdsPtr.lock()->getNumberOfTuples();
  size_t totalNumElementsDest = m_FeatureIdsPtr.lock()->getNumberOfTuples();
  DataContainer::Pointer dc = getDataContainerArray()->getDataContainer(getFeatureIDsArrayPath().getDataContainerName());
  ImageGeom::Pointer image = dc->getGeometryAs<ImageGeom>();
  SizeVec3Type udims = image->getDimensions();
  FloatVec3Type uorigin = image->getOrigin();
  FloatVec3Type uspacing = image->getSpacing();

  std::array<int64_t, 3> dims = {
      static_cast<int64_t>(udims[0]),
      static_cast<int64_t>(udims[1]),
      static_cast<int64_t>(udims[2]),
  };

  std::array<float, 3> origin = {
      static_cast<float>(uorigin[0]),
      static_cast<float>(uorigin[1]),
      static_cast<float>(uorigin[2]),
  };

  std::array<float, 3> spacing = {
      static_cast<float>(uspacing[0]),
      static_cast<float>(uspacing[1]),
      static_cast<float>(uspacing[2]),
  };

  for(size_t i = 0; i < totalNumElementsDest; i++)
  {
    int64_t zindex = i / (dims[0] * dims[1]);
    int64_t tempMod = i % (dims[0] * dims[1]);
    int64_t yindex = tempMod / (dims[0]);
    int64_t xindex = tempMod % (dims[0]);

    float currentx = origin[0] + spacing[0] * xindex;
    float currenty = origin[1] + spacing[1] * yindex;
    float currentz = origin[2] + spacing[2] * zindex;

    bool inBounds = false;
    size_t k = 0;

    while(!inBounds && k < totalNumFIDs)
    {
      float xmin = m_BoxCenter[3 * k] - m_BoxDims[3 * k] / 2.0;
      float xmax = m_BoxCenter[3 * k] + m_BoxDims[3 * k] / 2.0;
      float ymin = m_BoxCenter[3 * k + 1] - m_BoxDims[3 * k + 1] / 2.0;
      float ymax = m_BoxCenter[3 * k + 1] + m_BoxDims[3 * k + 1] / 2.0;
      float zmin = m_BoxCenter[3 * k + 2] - m_BoxDims[3 * k + 2] / 2.0;
      float zmax = m_BoxCenter[3 * k + 2] + m_BoxDims[3 * k + 2] / 2.0;

      inBounds = IsPointInBounds(xmax, xmin, ymax, ymin, zmax, zmin, currentx, currenty, currentz);

      if(inBounds)
      {
        m_FeatureIds[i] = m_BoxFeatureIds[k];
      }

      k++;
    }
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void GenerateFeatureIDsbyBoundingBoxes::checkBoundingBoxEdge()
{
  size_t totalNumFIDs = m_BoxFeatureIdsPtr.lock()->getNumberOfTuples();
  size_t totalNumElementsDest = m_FeatureIdsPtr.lock()->getNumberOfTuples();

  for(size_t i = 0; i < totalNumElementsDest; i++)
  {
    for(size_t k = 0; k < totalNumFIDs; k++)
    {
    }
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void GenerateFeatureIDsbyBoundingBoxes::checkBoundingBoxVertex()
{
  size_t totalNumFIDs = m_BoxFeatureIdsPtr.lock()->getNumberOfTuples();
  size_t totalNumElementsDest = m_FeatureIdsPtr.lock()->getNumberOfTuples();
  DataContainer::Pointer dc = getDataContainerArray()->getDataContainer(getFeatureIDsArrayPath().getDataContainerName());
  VertexGeom::Pointer vertex = dc->getGeometryAs<VertexGeom>();

  float* vertices = vertex->getVertexPointer(0);
  for(size_t i = 0; i < totalNumElementsDest; i++)
  {
    bool inBounds = false;
    size_t k = 0;
    float currentx = vertices[3 * i];
    float currenty = vertices[3 * i + 1];
    float currentz = vertices[3 * i + 2];

    while(!inBounds && k < totalNumFIDs)
    {
      float xmin = m_BoxCenter[3 * k] - m_BoxDims[3 * k] / 2.0;
      float xmax = m_BoxCenter[3 * k] + m_BoxDims[3 * k] / 2.0;
      float ymin = m_BoxCenter[3 * k + 1] - m_BoxDims[3 * k + 1] / 2.0;
      float ymax = m_BoxCenter[3 * k + 1] + m_BoxDims[3 * k + 1] / 2.0;
      float zmin = m_BoxCenter[3 * k + 2] - m_BoxDims[3 * k + 2] / 2.0;
      float zmax = m_BoxCenter[3 * k + 2] + m_BoxDims[3 * k + 2] / 2.0;

      inBounds = IsPointInBounds(xmax, xmin, ymax, ymin, zmax, zmin, currentx, currenty, currentz);

      if(inBounds)
      {
        m_FeatureIds[i] = m_BoxFeatureIds[k];
      }

      k++;
    }
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void GenerateFeatureIDsbyBoundingBoxes::execute()
{
  initialize();
  dataCheck();
  if(getErrorCode() < 0)
  {
    return;
  }

  if(getCancel())
  {
    return;
  }

  if(m_DestAttributeMatrixType == AttributeMatrix::Type::CellFeature)
  {
    checkBoundingBoxImage();
  }
  else if(m_DestAttributeMatrixType == AttributeMatrix::Type::VertexFeature)
  {
    checkBoundingBoxVertex();
  }
  // else if(m_DestAttributeMatrixType == AttributeMatrix::Type::EdgeFeature)
  // {
  //   checkBoundingBoxEdge();
  // }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AbstractFilter::Pointer GenerateFeatureIDsbyBoundingBoxes::newFilterInstance(bool copyFilterParameters) const
{
  GenerateFeatureIDsbyBoundingBoxes::Pointer filter = GenerateFeatureIDsbyBoundingBoxes::New();
  if(copyFilterParameters)
  {
    copyFilterParameterInstanceVariables(filter.get());
  }
  return filter;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString GenerateFeatureIDsbyBoundingBoxes::getCompiledLibraryName() const
{
  return DREAM3DReviewConstants::DREAM3DReviewBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString GenerateFeatureIDsbyBoundingBoxes::getBrandingString() const
{
  return "DREAM3DReview";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString GenerateFeatureIDsbyBoundingBoxes::getFilterVersion() const
{
  QString version;
  QTextStream vStream(&version);
  vStream << DREAM3DReview::Version::Major() << "." << DREAM3DReview::Version::Minor() << "." << DREAM3DReview::Version::Patch();
  return version;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString GenerateFeatureIDsbyBoundingBoxes::getGroupName() const
{
  return SIMPL::FilterGroups::Unsupported;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString GenerateFeatureIDsbyBoundingBoxes::getSubGroupName() const
{
  return "DREAM3DReview";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString GenerateFeatureIDsbyBoundingBoxes::getHumanLabel() const
{
  return "Generate FeatureIDs by Bounding Boxes";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QUuid GenerateFeatureIDsbyBoundingBoxes::getUuid() const
{
  return QUuid("{8b2fa51e-3bad-58ec-9fbf-b03b065e67fc}");
}

// -----------------------------------------------------------------------------
GenerateFeatureIDsbyBoundingBoxes::Pointer GenerateFeatureIDsbyBoundingBoxes::NullPointer()
{
  return Pointer(static_cast<Self*>(nullptr));
}

// -----------------------------------------------------------------------------
std::shared_ptr<GenerateFeatureIDsbyBoundingBoxes> GenerateFeatureIDsbyBoundingBoxes::New()
{
  struct make_shared_enabler : public GenerateFeatureIDsbyBoundingBoxes
  {
  };
  std::shared_ptr<make_shared_enabler> val = std::make_shared<make_shared_enabler>();
  val->setupFilterParameters();
  return val;
}

// -----------------------------------------------------------------------------
QString GenerateFeatureIDsbyBoundingBoxes::getNameOfClass() const
{
  return QString("GenerateFeatureIDsbyBoundingBoxes");
}

// -----------------------------------------------------------------------------
QString GenerateFeatureIDsbyBoundingBoxes::ClassName()
{
  return QString("GenerateFeatureIDsbyBoundingBoxes");
}

// -----------------------------------------------------------------------------
void GenerateFeatureIDsbyBoundingBoxes::setFeatureIDsArrayPath(const DataArrayPath& value)
{
  m_FeatureIDsArrayPath = value;
}
// -----------------------------------------------------------------------------
DataArrayPath GenerateFeatureIDsbyBoundingBoxes::getFeatureIDsArrayPath() const
{
  return m_FeatureIDsArrayPath;
}
// -----------------------------------------------------------------------------
void GenerateFeatureIDsbyBoundingBoxes::setFeatureAttributeMatrixArrayPath(const DataArrayPath& value)
{
  m_FeatureAttributeMatrixArrayPath = value;
}
// -----------------------------------------------------------------------------
DataArrayPath GenerateFeatureIDsbyBoundingBoxes::getFeatureAttributeMatrixArrayPath() const
{
  return m_FeatureAttributeMatrixArrayPath;
}
// -----------------------------------------------------------------------------
void GenerateFeatureIDsbyBoundingBoxes::setBoxCenterArrayPath(const DataArrayPath& value)
{
  m_BoxCenterArrayPath = value;
}
// -----------------------------------------------------------------------------
DataArrayPath GenerateFeatureIDsbyBoundingBoxes::getBoxCenterArrayPath() const
{
  return m_BoxCenterArrayPath;
}
// -----------------------------------------------------------------------------
void GenerateFeatureIDsbyBoundingBoxes::setBoxDimensionsArrayPath(const DataArrayPath& value)
{
  m_BoxDimensionsArrayPath = value;
}
// -----------------------------------------------------------------------------
DataArrayPath GenerateFeatureIDsbyBoundingBoxes::getBoxDimensionsArrayPath() const
{
  return m_BoxDimensionsArrayPath;
}
// -----------------------------------------------------------------------------
void GenerateFeatureIDsbyBoundingBoxes::setBoxFeatureIDsArrayPath(const DataArrayPath& value)
{
  m_BoxFeatureIDsArrayPath = value;
}
// -----------------------------------------------------------------------------
DataArrayPath GenerateFeatureIDsbyBoundingBoxes::getBoxFeatureIDsArrayPath() const
{
  return m_BoxFeatureIDsArrayPath;
}
