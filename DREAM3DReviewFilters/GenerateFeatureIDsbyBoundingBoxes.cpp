/*
 * Your License or Copyright can go here
 */

#include "GenerateFeatureIDsbyBoundingBoxes.h"

#include <QtCore/QTextStream>

#include "SIMPLib/Common/Constants.h"

#include "SIMPLib/FilterParameters/DataArrayCreationFilterParameter.h"
#include "SIMPLib/FilterParameters/DataArraySelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/AttributeMatrixCreationFilterParameter.h"
#include "SIMPLib/DataContainers/DataContainerArray.h"
#include "SIMPLib/Geometry/ImageGeom.h"
#include "SIMPLib/Geometry/VertexGeom.h"


#include "DREAM3DReview/DREAM3DReviewConstants.h"
#include "DREAM3DReview/DREAM3DReviewVersion.h"

 /* Create Enumerations to allow the created Attribute Arrays to take part in renaming */
enum createdPathID : RenameDataPath::DataID_t
{
  AttributeMatrixID20 = 20,
};


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
GenerateFeatureIDsbyBoundingBoxes::GenerateFeatureIDsbyBoundingBoxes() 
:  m_FeatureIDsArrayPath("", "", "")
,  m_FeatureAttributeMatrixArrayPath("", "", "")
,  m_BoxCenterArrayPath("", "", "")
,  m_BoxDimensionsArrayPath("", "", "")
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
  DataArrayPath path(getFeatureIDsArrayPath().getDataContainerName(), getFeatureIDsArrayPath().getAttributeMatrixName(), "");
  AttributeMatrix::Pointer attrMat = getDataContainerArray()->getPrereqAttributeMatrixFromPath(this, path, -301);

  if (getErrorCode() < 0)
  {
    return;
  }

  AttributeMatrix::Type attrMatType = attrMat->getType();
  m_DestAttributeMatrixType = AttributeMatrix::Type::Unknown;

  switch (attrMatType)
  {
  case AttributeMatrix::Type::Vertex:
  {
    m_DestAttributeMatrixType = AttributeMatrix::Type::VertexFeature;
    break;
  }
  case AttributeMatrix::Type::Edge:
  {
    m_DestAttributeMatrixType = AttributeMatrix::Type::EdgeFeature;
    break;
  }
  case AttributeMatrix::Type::Face:
  {
    m_DestAttributeMatrixType = AttributeMatrix::Type::FaceFeature;
    break;
  }
  case AttributeMatrix::Type::Cell:
  {
    m_DestAttributeMatrixType = AttributeMatrix::Type::CellFeature;
    break;
  }
  case AttributeMatrix::Type::VertexFeature:
  {
    m_DestAttributeMatrixType = AttributeMatrix::Type::VertexEnsemble;
    break;
  }
  case AttributeMatrix::Type::EdgeFeature:
  {
    m_DestAttributeMatrixType = AttributeMatrix::Type::EdgeEnsemble;
    break;
  }
  case AttributeMatrix::Type::FaceFeature:
  {
    m_DestAttributeMatrixType = AttributeMatrix::Type::FaceEnsemble;
    break;
  }
  case AttributeMatrix::Type::CellFeature:
  {
    m_DestAttributeMatrixType = AttributeMatrix::Type::CellEnsemble;
    break;
  }
  case AttributeMatrix::Type::VertexEnsemble:
  {
    m_DestAttributeMatrixType = AttributeMatrix::Type::VertexEnsemble;
    break;
  }
  case AttributeMatrix::Type::EdgeEnsemble:
  {
    m_DestAttributeMatrixType = AttributeMatrix::Type::EdgeEnsemble;
    break;
  }
  case AttributeMatrix::Type::FaceEnsemble:
  {
    m_DestAttributeMatrixType = AttributeMatrix::Type::FaceEnsemble;
    break;
  }
  case AttributeMatrix::Type::CellEnsemble:
  {
    m_DestAttributeMatrixType = AttributeMatrix::Type::CellEnsemble;
    break;
  }
  default:
  {
    m_DestAttributeMatrixType = AttributeMatrix::Type::Generic;
    break;
  }
  }


  std::vector<size_t> cDims(1, 1);
  m_BoxFeatureIdsPtr = getDataContainerArray()->getPrereqArrayFromPath<DataArray<int32_t>>(this, getBoxFeatureIDsArrayPath(), cDims);
  if (nullptr != m_BoxFeatureIdsPtr.lock())
  {
    m_BoxFeatureIds = m_BoxFeatureIdsPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */

  cDims[0] = 3;
  m_BoxCenterPtr = getDataContainerArray()->getPrereqArrayFromPath<DataArray<float>>(this, getBoxCenterArrayPath(), cDims);
  if (nullptr != m_BoxCenterPtr.lock()) /* Validate the Weak Pointer wraps a non-nullptr pointer to a DataArray<T> object */
  {
    m_BoxCenter = m_BoxCenterPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */

  m_BoxDimsPtr = getDataContainerArray()->getPrereqArrayFromPath<DataArray<float>>(this, getBoxDimensionsArrayPath(), cDims);
  if (nullptr != m_BoxDimsPtr.lock()) /* Validate the Weak Pointer wraps a non-nullptr pointer to a DataArray<T> object */
  {
    m_BoxDims = m_BoxDimsPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */

  cDims[0] = 1;
  m_FeatureIdsPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<int32_t>>(this, getFeatureIDsArrayPath(), 0, cDims);                  /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
  if (nullptr != m_FeatureIdsPtr.lock()) /* Validate the Weak Pointer wraps a non-nullptr pointer to a DataArray<T> object */
  {
    m_FeatureIds = m_FeatureIdsPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */


  
  DataArrayPath path_bounding(getBoxFeatureIDsArrayPath().getDataContainerName(), getBoxFeatureIDsArrayPath().getAttributeMatrixName(), "");
  AttributeMatrix::Pointer attrMat_bounding = getDataContainerArray()->getPrereqAttributeMatrixFromPath(this, path_bounding, -301);

  if (getErrorCode() < 0)
  {
    return;
  }
 
  int32_t numFeatureIDs = attrMat_bounding->getNumberOfTuples();
  std::vector<size_t> tDims(1, numFeatureIDs+1);
  m->createNonPrereqAttributeMatrix(this, getFeatureAttributeMatrixArrayPath(), tDims, m_DestAttributeMatrixType, AttributeMatrixID20);

  

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------

bool IsPointInBounds(float xmax, float xmin, float ymax, float ymin, float zmax, float zmin, float x, float y, float z)
{
  bool inBounds = false;

  if ((x < xmax) && (x > xmin) && (y < ymax) && (y > ymin) && (z < zmax) && (z > zmin))
  {
    inBounds = true;
  }

  return inBounds;
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

  int64_t dims[3] = 
  {
    static_cast<int64_t>(udims[0]), static_cast<int64_t>(udims[1]), static_cast<int64_t>(udims[2]),
  };

  float origin[3] =
  {
    static_cast<float>(uorigin[0]), static_cast<float>(uorigin[1]), static_cast<float>(uorigin[2]),
  };

  float spacing[3] =
  {
    static_cast<float>(uspacing[0]), static_cast<float>(uspacing[1]), static_cast<float>(uspacing[2]),
  };

  bool inBounds = false;
  float xmin = 0;
  float xmax = 0;
  float ymin = 0;
  float ymax = 0;
  float zmin = 0;
  float zmax = 0;

  float currentx = 0;
  float currenty = 0;
  float currentz = 0;

  int64_t xindex = 0;
  int64_t yindex = 0; 
  int64_t zindex = 0;
  int64_t tempMod = 0;

  size_t k = 0;

  for (size_t i = 0; i < totalNumElementsDest; i++)
  {
    
    zindex = i / (dims[0] * dims[1]);
    tempMod = i % (dims[0] * dims[1]);
    yindex = tempMod / (dims[0]); 
    xindex = tempMod % (dims[0]);

    currentx = origin[0] + spacing[0] * xindex;
    currenty = origin[1] + spacing[1] * yindex;
    currentz = origin[2] + spacing[2] * zindex;

    inBounds = false;
    k = 0;
    

    while (!inBounds && k < totalNumFIDs)
    {
      xmin = m_BoxCenter[3 * k] - m_BoxDims[3 * k] / 2.0;
      xmax = m_BoxCenter[3 * k] + m_BoxDims[3 * k] / 2.0;
      ymin = m_BoxCenter[3 * k + 1] - m_BoxDims[3 * k + 1] / 2.0;
      ymax = m_BoxCenter[3 * k + 1] + m_BoxDims[3 * k + 1] / 2.0;
      zmin = m_BoxCenter[3 * k + 2] - m_BoxDims[3 * k + 2] / 2.0;
      zmax = m_BoxCenter[3 * k + 2] + m_BoxDims[3 * k + 2] / 2.0;

      inBounds = IsPointInBounds(xmax, xmin, ymax, ymin, zmax, zmin, currentx, currenty, currentz);

      if (inBounds)
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

  for (size_t i; i < totalNumElementsDest; i++)
  {
    for (size_t k; k < totalNumFIDs; k++)
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

  bool inBounds = false;
  float xmin = 0;
  float xmax = 0;
  float ymin = 0; 
  float ymax = 0;
  float zmin = 0; 
  float zmax = 0;

  float currentx = 0;
  float currenty = 0;
  float currentz = 0;

  float* vertices = vertex->getVertexPointer(0);
  size_t k = 0;
  for (size_t i = 0; i < totalNumElementsDest; i++)
  {
    inBounds = false;
    k = 0;
    currentx = vertices[3 * i];
    currenty = vertices[3 * i + 1];
    currentz = vertices[3 * i + 2];

    while (!inBounds && k < totalNumFIDs)
    {
      xmin = m_BoxCenter[3 * k] - m_BoxDims[3 * k] / 2.0;
      xmax = m_BoxCenter[3 * k] + m_BoxDims[3 * k] / 2.0;
      ymin = m_BoxCenter[3 * k + 1] - m_BoxDims[3 * k + 1] / 2.0;
      ymax = m_BoxCenter[3 * k + 1] + m_BoxDims[3 * k + 1] / 2.0;
      zmin = m_BoxCenter[3 * k + 2] - m_BoxDims[3 * k + 2] / 2.0;
      zmax = m_BoxCenter[3 * k + 2] + m_BoxDims[3 * k + 2] / 2.0;

      inBounds = IsPointInBounds(xmax, xmin, ymax, ymin, zmax, zmin, currentx, currenty, currentz);

      if (inBounds)
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
  if(getErrorCode() < 0) { return; }

  if (getCancel()) { return; }

  if (getWarningCode() < 0)
  {
    QString ss = QObject::tr("Some warning message");
    setWarningCondition(-88888888, ss);
  }

  if (getErrorCode() < 0)
  {
    QString ss = QObject::tr("Some error message");
    setErrorCondition(-99999999, ss);
    return;
  }

  if (m_DestAttributeMatrixType == AttributeMatrix::Type::CellFeature)
  {
    checkBoundingBoxImage();
  }
  else if (m_DestAttributeMatrixType == AttributeMatrix::Type::VertexFeature)
  {
    checkBoundingBoxVertex();
  }
  else if (m_DestAttributeMatrixType == AttributeMatrix::Type::EdgeFeature)
  {
    checkBoundingBoxEdge();
  }

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
  vStream <<  DREAM3DReview::Version::Major() << "." << DREAM3DReview::Version::Minor() << "." << DREAM3DReview::Version::Patch();
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


