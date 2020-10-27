/*
 * Your License or Copyright can go here
 */

#include "GenerateMaskFromSimpleShapes.h"

#include <QtCore/QTextStream>

#include "SIMPLib/Common/Constants.h"
#include "SIMPLib/DataContainers/DataContainerArray.h"
#include "SIMPLib/FilterParameters/DataArrayCreationFilterParameter.h"
#include "SIMPLib/FilterParameters/DataArraySelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedChoicesFilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedPathCreationFilterParameter.h"
#include "SIMPLib/Geometry/ImageGeom.h"
#include "SIMPLib/Geometry/VertexGeom.h"

#include "DREAM3DReview/DREAM3DReviewConstants.h"
#include "DREAM3DReview/DREAM3DReviewVersion.h"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
GenerateMaskFromSimpleShapes::GenerateMaskFromSimpleShapes()
: m_MaskArrayPath("", "", "Mask")
, m_CentersArrayPath("", "", "Centroids")
, m_AxesLengthArrayPath("", "", "AxisLengths")
, m_BoxDimensionsArrayPath("", "", "")
, m_CylinderRadiusArrayPath("", "", "Radii")
, m_CylinderHeightArrayPath("", "", "Heights")
, m_MaskShape(0)
{
  initialize();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
GenerateMaskFromSimpleShapes::~GenerateMaskFromSimpleShapes() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void GenerateMaskFromSimpleShapes::initialize()
{
  clearErrorCode();
  clearWarningCode();
  setCancel(false);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void GenerateMaskFromSimpleShapes::setupFilterParameters()
{
  FilterParameterVectorType parameters;
  {
    LinkedChoicesFilterParameter::Pointer parameter = LinkedChoicesFilterParameter::New();
    parameter->setHumanLabel("Mask Shape");
    parameter->setPropertyName("MaskShape");
    parameter->setSetterCallback(SIMPL_BIND_SETTER(GenerateMaskFromSimpleShapes, this, MaskShape));
    parameter->setGetterCallback(SIMPL_BIND_GETTER(GenerateMaskFromSimpleShapes, this, MaskShape));

    QVector<QString> choices;
    choices.push_back("Ellipsoid");
    choices.push_back("Box");
    choices.push_back("Cylinder");
    parameter->setChoices(choices);
    QStringList linkedProps;
    linkedProps << "AxesLengthArrayPath"
                << "BoxDimensionsArrayPath"
                << "CylinderHeightArrayPath"
                << "CylinderRadiusArrayPath";
    parameter->setLinkedProperties(linkedProps);
    parameter->setEditable(false);
    parameter->setCategory(FilterParameter::Parameter);
    parameters.push_back(parameter);
  }

  DataArrayCreationFilterParameter::RequirementType dacReq;
  parameters.push_back(SIMPL_NEW_DA_CREATION_FP("Mask", MaskArrayPath, FilterParameter::CreatedArray, GenerateMaskFromSimpleShapes, dacReq));
  DataArraySelectionFilterParameter::RequirementType dasReq;
  parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Centers", CentersArrayPath, FilterParameter::RequiredArray, GenerateMaskFromSimpleShapes, dasReq));
  parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Ellipsoid Axes Lengths", AxesLengthArrayPath, FilterParameter::RequiredArray, GenerateMaskFromSimpleShapes, dasReq, 0));
  parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Box Dimensions", BoxDimensionsArrayPath, FilterParameter::RequiredArray, GenerateMaskFromSimpleShapes, dasReq, 1));
  parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Cylinder Radii", CylinderRadiusArrayPath, FilterParameter::RequiredArray, GenerateMaskFromSimpleShapes, dasReq, 2));
  parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Cylinder Heights", CylinderHeightArrayPath, FilterParameter::RequiredArray, GenerateMaskFromSimpleShapes, dasReq, 2));
  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void GenerateMaskFromSimpleShapes::dataCheck()
{
  clearErrorCode();
  clearWarningCode();

  AttributeMatrix::Pointer attrMat = getDataContainerArray()->getPrereqAttributeMatrixFromPath(this, getMaskArrayPath(), -301);
  if(getErrorCode() < 0)
  {
    return;
  }

  AttributeMatrix::Type attrMatType = attrMat->getType();

  if(attrMatType == AttributeMatrix::Type::Vertex || attrMatType == AttributeMatrix::Type::Cell)
  {
    QString ss = QObject::tr("The Attribute Matrix must have a cell or vertex geometry.");
    setErrorCondition(-5555, ss);
    return;
  }

  std::vector<size_t> cDims = {1};
  m_MaskPtr =
      getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<bool>>(this, getMaskArrayPath(), false, cDims); /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
  if(nullptr != m_MaskPtr.lock()) /* Validate the Weak Pointer wraps a non-nullptr pointer to a DataArray<T> object */
  {
    m_Mask = m_MaskPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */

  cDims[0] = 3;
  m_CentersPtr = getDataContainerArray()->getPrereqArrayFromPath<DataArray<float>>(this, getCentersArrayPath(), cDims);
  if(nullptr != m_CentersPtr.lock()) /* Validate the Weak Pointer wraps a non-nullptr pointer to a DataArray<T> object */
  {
    m_Centers = m_CentersPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */

  if(m_MaskShape == 0)
  {
    cDims[0] = 3;
    m_AxisLengthsPtr = getDataContainerArray()->getPrereqArrayFromPath<DataArray<float>>(this, getAxesLengthArrayPath(), cDims);
    if(nullptr != m_AxisLengthsPtr.lock()) /* Validate the Weak Pointer wraps a non-nullptr pointer to a DataArray<T> object */
    {
      m_AxisLengths = m_AxisLengthsPtr.lock()->getPointer(0);
    } /* Now assign the raw pointer to data from the DataArray<T> object */
  }

  if(m_MaskShape == 1)
  {
    cDims[0] = 3;
    m_BoxDimsPtr = getDataContainerArray()->getPrereqArrayFromPath<DataArray<float>>(this, getBoxDimensionsArrayPath(), cDims);
    if(nullptr != m_BoxDimsPtr.lock()) /* Validate the Weak Pointer wraps a non-nullptr pointer to a DataArray<T> object */
    {
      m_BoxDims = m_BoxDimsPtr.lock()->getPointer(0);
    } /* Now assign the raw pointer to data from the DataArray<T> object */
  }

  if(m_MaskShape == 2)
  {
    cDims[0] = 1;
    m_CylinderRadPtr = getDataContainerArray()->getPrereqArrayFromPath<DataArray<float>>(this, getCylinderRadiusArrayPath(), cDims);
    if(nullptr != m_CylinderRadPtr.lock()) /* Validate the Weak Pointer wraps a non-nullptr pointer to a DataArray<T> object */
    {
      m_CylinderRad = m_CylinderRadPtr.lock()->getPointer(0);
    } /* Now assign the raw pointer to data from the DataArray<T> object */

    cDims[0] = 1;
    m_CylinderHeightPtr = getDataContainerArray()->getPrereqArrayFromPath<DataArray<float>>(this, getCylinderHeightArrayPath(), cDims);
    if(nullptr != m_CylinderHeightPtr.lock()) /* Validate the Weak Pointer wraps a non-nullptr pointer to a DataArray<T> object */
    {
      m_CylinderHeight = m_CylinderHeightPtr.lock()->getPointer(0);
    } /* Now assign the raw pointer to data from the DataArray<T> object */
  }
}
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
bool IsPointInBoxBounds(float xcenter, float ycenter, float zcenter, float xwidth, float yheight, float zdepth, float x, float y, float z)
{
  bool inBounds = false;

  float xmin = xcenter - xwidth / 2.0;
  float xmax = xcenter + xwidth / 2.0;
  float ymin = ycenter - yheight / 2.0;
  float ymax = ycenter + yheight / 2.0;
  float zmin = zcenter - zdepth / 2.0;
  float zmax = zcenter + zdepth / 2.0;

  if((x < xmax) && (x > xmin) && (y < ymax) && (y > ymin) && (z < zmax) && (z > zmin))
  {
    inBounds = true;
  }

  return inBounds;
}
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
bool IsPointInCylinderBounds(float xcenter, float ycenter, float zcenter, float radius, float zdepth, float x, float y, float z)
{
  bool inBounds = false;

  float zmin = zcenter - zdepth / 2.0;
  float zmax = zcenter + zdepth / 2.0;

  if((sqrt((xcenter - x) * (xcenter - x) + (ycenter - y) * (ycenter - y)) < radius) && (z < zmax) && (z > zmin))
  {
    inBounds = true;
  }

  return inBounds;
}
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
bool IsPointInEllipsoidBounds(float xcenter, float ycenter, float zcenter, float a, float b, float c, float x, float y, float z)
{
  bool inBounds = false;

  if(((xcenter - x) * (xcenter - x) / (a * a) + (ycenter - y) * (ycenter - y) / (b * b) + (zcenter - z) * (zcenter - z) / (c * c)) < 1)
  {
    inBounds = true;
  }

  return inBounds;
}
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void GenerateMaskFromSimpleShapes::createImageMask()
{
  size_t totalNumFIDs = m_CentersPtr.lock()->getNumberOfTuples();
  size_t totalNumElementsDest = m_MaskPtr.lock()->getNumberOfTuples();
  DataContainer::Pointer dc = getDataContainerArray()->getDataContainer(getMaskArrayPath().getDataContainerName());
  ImageGeom::Pointer image = dc->getGeometryAs<ImageGeom>();
  SizeVec3Type udims = image->getDimensions();
  FloatVec3Type uorigin = image->getOrigin();
  FloatVec3Type uspacing = image->getSpacing();

  int64_t dims[3] = {
      static_cast<int64_t>(udims[0]),
      static_cast<int64_t>(udims[1]),
      static_cast<int64_t>(udims[2]),
  };

  float origin[3] = {
      static_cast<float>(uorigin[0]),
      static_cast<float>(uorigin[1]),
      static_cast<float>(uorigin[2]),
  };

  std::array<float, 3> spacing = {
      static_cast<float>(uspacing[0]),
      static_cast<float>(uspacing[1]),
      static_cast<float>(uspacing[2]),
  };

  bool inBounds = false;

  float currentx = 0;
  float currenty = 0;
  float currentz = 0;

  int64_t xindex = 0;
  int64_t yindex = 0;
  int64_t zindex = 0;
  int64_t tempMod = 0;

  size_t k = 0;

  if(m_MaskShape == 0)
  {
    for(size_t i = 0; i < totalNumElementsDest; i++)
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

      while(!inBounds && k < totalNumFIDs)
      {
        inBounds = IsPointInEllipsoidBounds(m_Centers[3 * k], m_Centers[3 * k + 1], m_Centers[3 * k + 2], m_AxisLengths[3 * k], m_AxisLengths[3 * k + 1], m_AxisLengths[3 * k + 2], currentx, currenty,
                                            currentz);

        if(inBounds)
        {
          m_Mask[i] = true;
        }

        k++;
      }
    }
  }

  if(m_MaskShape == 1)
  {
    for(size_t i = 0; i < totalNumElementsDest; i++)
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

      while(!inBounds && k < totalNumFIDs)
      {
        inBounds = IsPointInBoxBounds(m_Centers[3 * k], m_Centers[3 * k + 1], m_Centers[3 * k + 2], m_BoxDims[0], m_BoxDims[1], m_BoxDims[2], currentx, currenty, currentz);

        if(inBounds)
        {
          m_Mask[i] = true;
        }

        k++;
      }
    }
  }

  if(m_MaskShape == 2)
  {
    for(size_t i = 0; i < totalNumElementsDest; i++)
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

      while(!inBounds && k < totalNumFIDs)
      {
        inBounds = IsPointInCylinderBounds(m_Centers[3 * k], m_Centers[3 * k + 1], m_Centers[3 * k + 2], m_CylinderRad[k], m_CylinderHeight[k], currentx, currenty, currentz);

        if(inBounds)
        {
          m_Mask[i] = true;
        }

        k++;
      }
    }
  }
}
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void GenerateMaskFromSimpleShapes::createVertexMask()
{
  size_t totalNumCenters = m_CentersPtr.lock()->getNumberOfTuples();
  size_t totalNumElementsDest = m_MaskPtr.lock()->getNumberOfTuples();
  DataContainer::Pointer dc = getDataContainerArray()->getDataContainer(getMaskArrayPath().getDataContainerName());
  VertexGeom::Pointer vertex = dc->getGeometryAs<VertexGeom>();

  bool inBounds = false;

  float currentx = 0;
  float currenty = 0;
  float currentz = 0;

  float* vertices = vertex->getVertexPointer(0);
  size_t k = 0;

  if(m_MaskShape == 0)
  {
    for(size_t i = 0; i < totalNumElementsDest; i++)
    {
      inBounds = false;
      k = 0;
      currentx = vertices[3 * i];
      currenty = vertices[3 * i + 1];
      currentz = vertices[3 * i + 2];

      while(!inBounds && k < totalNumCenters)
      {
        inBounds = IsPointInEllipsoidBounds(m_Centers[3 * k], m_Centers[3 * k + 1], m_Centers[3 * k + 2], m_AxisLengths[3 * k], m_AxisLengths[3 * k + 1], m_AxisLengths[3 * k + 2], currentx, currenty,
                                            currentz);

        if(inBounds)
        {
          m_Mask[i] = true;
        }

        k++;
      }
    }
  }

  if(m_MaskShape == 1)
  {
    for(size_t i = 0; i < totalNumElementsDest; i++)
    {
      inBounds = false;
      k = 0;
      currentx = vertices[3 * i];
      currenty = vertices[3 * i + 1];
      currentz = vertices[3 * i + 2];

      while(!inBounds && k < totalNumCenters)
      {
        inBounds = IsPointInBoxBounds(m_Centers[3 * k], m_Centers[3 * k + 1], m_Centers[3 * k + 2], m_BoxDims[3 * k], m_BoxDims[3 * k + 1], m_BoxDims[3 * k + 2], currentx, currenty, currentz);

        if(inBounds)
        {
          m_Mask[i] = true;
        }

        k++;
      }
    }
  }

  if(m_MaskShape == 2)
  {
    for(size_t i = 0; i < totalNumElementsDest; i++)
    {
      inBounds = false;
      k = 0;
      currentx = vertices[3 * i];
      currenty = vertices[3 * i + 1];
      currentz = vertices[3 * i + 2];

      while(!inBounds && k < totalNumCenters)
      {
        inBounds = IsPointInCylinderBounds(m_Centers[3 * k], m_Centers[3 * k + 1], m_Centers[3 * k + 2], m_CylinderRad[k], m_CylinderHeight[k], currentx, currenty, currentz);

        if(inBounds)
        {
          m_Mask[i] = true;
        }

        k++;
      }
    }
  }
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void GenerateMaskFromSimpleShapes::execute()
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

  DataContainer::Pointer m = getDataContainerArray()->getPrereqDataContainer(this, getMaskArrayPath().getDataContainerName(), false);
  DataArrayPath path(getMaskArrayPath().getDataContainerName(), getMaskArrayPath().getAttributeMatrixName(), "");
  AttributeMatrix::Pointer attrMat = getDataContainerArray()->getPrereqAttributeMatrixFromPath(this, path, -301);

  if(getErrorCode() < 0)
  {
    return;
  }

  AttributeMatrix::Type attrMatType = attrMat->getType();

  if(attrMatType == AttributeMatrix::Type::Cell)
  {
    createImageMask();
  }
  else if(attrMatType == AttributeMatrix::Type::Vertex)
  {
    createVertexMask();
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AbstractFilter::Pointer GenerateMaskFromSimpleShapes::newFilterInstance(bool copyFilterParameters) const
{
  GenerateMaskFromSimpleShapes::Pointer filter = GenerateMaskFromSimpleShapes::New();
  if(copyFilterParameters)
  {
    copyFilterParameterInstanceVariables(filter.get());
  }
  return filter;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString GenerateMaskFromSimpleShapes::getCompiledLibraryName() const
{
  return DREAM3DReviewConstants::DREAM3DReviewBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString GenerateMaskFromSimpleShapes::getBrandingString() const
{
  return "DREAM3DReview";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString GenerateMaskFromSimpleShapes::getFilterVersion() const
{
  QString version;
  QTextStream vStream(&version);
  vStream << DREAM3DReview::Version::Major() << "." << DREAM3DReview::Version::Minor() << "." << DREAM3DReview::Version::Patch();
  return version;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString GenerateMaskFromSimpleShapes::getGroupName() const
{
  return SIMPL::FilterGroups::Unsupported;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString GenerateMaskFromSimpleShapes::getSubGroupName() const
{
  return "DREAM3DReview";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString GenerateMaskFromSimpleShapes::getHumanLabel() const
{
  return "Generate Mask From Simple Shapes";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QUuid GenerateMaskFromSimpleShapes::getUuid() const
{
  return QUuid("{e4d6fda0-1143-56cc-9d00-c09a94e2f501}");
}

// -----------------------------------------------------------------------------
GenerateMaskFromSimpleShapes::Pointer GenerateMaskFromSimpleShapes::NullPointer()
{
  return Pointer(static_cast<Self*>(nullptr));
}

// -----------------------------------------------------------------------------
std::shared_ptr<GenerateMaskFromSimpleShapes> GenerateMaskFromSimpleShapes::New()
{
  struct make_shared_enabler : public GenerateMaskFromSimpleShapes
  {
  };
  std::shared_ptr<make_shared_enabler> val = std::make_shared<make_shared_enabler>();
  val->setupFilterParameters();
  return val;
}

// -----------------------------------------------------------------------------
QString GenerateMaskFromSimpleShapes::getNameOfClass() const
{
  return QString("GenerateMaskFromSimpleShapes");
}

// -----------------------------------------------------------------------------
QString GenerateMaskFromSimpleShapes::ClassName()
{
  return QString("GenerateMaskFromSimpleShapes");
}

// -----------------------------------------------------------------------------
void GenerateMaskFromSimpleShapes::setMaskArrayPath(const DataArrayPath& value)
{
  m_MaskArrayPath = value;
}
// -----------------------------------------------------------------------------
DataArrayPath GenerateMaskFromSimpleShapes::getMaskArrayPath() const
{
  return m_MaskArrayPath;
}
// -----------------------------------------------------------------------------
void GenerateMaskFromSimpleShapes::setCentersArrayPath(const DataArrayPath& value)
{
  m_CentersArrayPath = value;
}
// -----------------------------------------------------------------------------
DataArrayPath GenerateMaskFromSimpleShapes::getCentersArrayPath() const
{
  return m_CentersArrayPath;
}
// -----------------------------------------------------------------------------
void GenerateMaskFromSimpleShapes::setAxesLengthArrayPath(const DataArrayPath& value)
{
  m_AxesLengthArrayPath = value;
}
// -----------------------------------------------------------------------------
DataArrayPath GenerateMaskFromSimpleShapes::getAxesLengthArrayPath() const
{
  return m_AxesLengthArrayPath;
}
// -----------------------------------------------------------------------------
void GenerateMaskFromSimpleShapes::setBoxDimensionsArrayPath(const DataArrayPath& value)
{
  m_BoxDimensionsArrayPath = value;
}
// -----------------------------------------------------------------------------
DataArrayPath GenerateMaskFromSimpleShapes::getBoxDimensionsArrayPath() const
{
  return m_BoxDimensionsArrayPath;
}
// -----------------------------------------------------------------------------
void GenerateMaskFromSimpleShapes::setCylinderRadiusArrayPath(const DataArrayPath& value)
{
  m_CylinderRadiusArrayPath = value;
}
// -----------------------------------------------------------------------------
DataArrayPath GenerateMaskFromSimpleShapes::getCylinderRadiusArrayPath() const
{
  return m_CylinderRadiusArrayPath;
}
// -----------------------------------------------------------------------------
void GenerateMaskFromSimpleShapes::setCylinderHeightArrayPath(const DataArrayPath& value)
{
  m_CylinderHeightArrayPath = value;
}
// -----------------------------------------------------------------------------
DataArrayPath GenerateMaskFromSimpleShapes::getCylinderHeightArrayPath() const
{
  return m_CylinderHeightArrayPath;
}

// -----------------------------------------------------------------------------
void GenerateMaskFromSimpleShapes::setMaskShape(int value)
{
  m_MaskShape = value;
}

// -----------------------------------------------------------------------------
int GenerateMaskFromSimpleShapes::getMaskShape() const
{
  return m_MaskShape;
}
