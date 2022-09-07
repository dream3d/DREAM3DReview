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

#include "AlignGeometries.h"

#include <QtCore/QTextStream>

#include "SIMPLib/Common/Constants.h"
#include "SIMPLib/DataContainers/DataContainerArray.h"
#include "SIMPLib/FilterParameters/ChoiceFilterParameter.h"
#include "SIMPLib/FilterParameters/DataContainerSelectionFilterParameter.h"
#include "SIMPLib/Geometry/EdgeGeom.h"
#include "SIMPLib/Geometry/IGeometry2D.h"
#include "SIMPLib/Geometry/IGeometry3D.h"
#include "SIMPLib/Geometry/ImageGeom.h"
#include "SIMPLib/Geometry/RectGridGeom.h"
#include "SIMPLib/Geometry/VertexGeom.h"
#include "SIMPLib/FilterParameters/LinkedChoicesFilterParameter.h"

#include "DREAM3DReview/DREAM3DReviewConstants.h"
#include "DREAM3DReview/DREAM3DReviewVersion.h"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AlignGeometries::AlignGeometries()
{
  initialize();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AlignGeometries::~AlignGeometries() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void AlignGeometries::initialize()
{
  clearErrorCode();
  clearWarningCode();
  setCancel(false);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void AlignGeometries::setupFilterParameters()
{
  FilterParameterVectorType parameters;

  {
    ChoiceFilterParameter::Pointer parameter = ChoiceFilterParameter::New();
    parameter->setHumanLabel("Alignment Type");
    parameter->setPropertyName("AlignmentType");
    parameter->setSetterCallback(SIMPL_BIND_SETTER(AlignGeometries, this, AlignmentType));
    parameter->setGetterCallback(SIMPL_BIND_GETTER(AlignGeometries, this, AlignmentType));
    std::vector<QString> choices;
    choices.push_back("Origin");
    choices.push_back("Centroid");
    choices.push_back("Plane (XY Min)");
    choices.push_back("Plane (XY Max)");
    choices.push_back("Plane (XZ Min)");
    choices.push_back("Plane (XZ Max)");
    choices.push_back("Plane (YZ Min)");
    choices.push_back("Plane (YZ Max)");
    parameter->setChoices(choices);
    parameter->setCategory(FilterParameter::Category::Parameter);
    parameter->setEditable(false);
    parameters.push_back(parameter);
  }

  DataContainerSelectionFilterParameter::RequirementType dcsReq;
  parameters.push_back(SIMPL_NEW_DC_SELECTION_FP("Moving Geometry", MovingGeometry, FilterParameter::Category::RequiredArray, AlignGeometries, dcsReq));
  parameters.push_back(SIMPL_NEW_DC_SELECTION_FP("Target Geometry", TargetGeometry, FilterParameter::Category::RequiredArray, AlignGeometries, dcsReq));



  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void AlignGeometries::dataCheck()
{
  clearErrorCode();
  clearWarningCode();

  getDataContainerArray()->getPrereqGeometryFromDataContainer<IGeometry>(this, getMovingGeometry());
  getDataContainerArray()->getPrereqGeometryFromDataContainer<IGeometry>(this, getTargetGeometry());

  if(getAlignmentType() < 0 || getAlignmentType() > 7)
  {
    QString ss = QObject::tr("Invalid selection for alignment type. Alignments types are: Origin=0, Centroid=1, Plane=...");
    setErrorCondition(-1, ss);
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
FloatVec3Type extractOrigin(const IGeometry::Pointer& geometry)
{
  if(ImageGeom::Pointer image = std::dynamic_pointer_cast<ImageGeom>(geometry))
  {
    return image->getOrigin();
  }
  if(RectGridGeom::Pointer rectGrid = std::dynamic_pointer_cast<RectGridGeom>(geometry))
  {
    FloatArrayType::Pointer xBounds = rectGrid->getXBounds();
    FloatArrayType::Pointer yBounds = rectGrid->getYBounds();
    FloatArrayType::Pointer zBounds = rectGrid->getZBounds();
    FloatVec3Type origin(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max());
    for(size_t i = 0; i < xBounds->getNumberOfTuples(); i++)
    {
      if(xBounds->getValue(i) < origin[0])
      {
        origin[0] = xBounds->getValue(i);
      }
    }
    for(size_t i = 0; i < yBounds->getNumberOfTuples(); i++)
    {
      if(yBounds->getValue(i) < origin[1])
      {
        origin[1] = yBounds->getValue(i);
      }
    }
    for(size_t i = 0; i < zBounds->getNumberOfTuples(); i++)
    {
      if(zBounds->getValue(i) < origin[2])
      {
        origin[2] = zBounds->getValue(i);
      }
    }
    return origin;
  }
  if(VertexGeom::Pointer vertex = std::dynamic_pointer_cast<VertexGeom>(geometry))
  {
    float* vertices = vertex->getVertexPointer(0);
    FloatVec3Type origin(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max());

    for(MeshIndexType i = 0; i < vertex->getNumberOfVertices(); i++)
    {
      for(size_t j = 0; j < 3; j++)
      {
        if(vertices[3 * i + j] < origin[j])
        {
          origin[j] = vertices[3 * i + j];
        }
      }
    }
    return origin;
  }
  if(EdgeGeom::Pointer edge = std::dynamic_pointer_cast<EdgeGeom>(geometry))
  {
    float* vertices = edge->getVertexPointer(0);
    FloatVec3Type origin(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max());

    for(MeshIndexType i = 0; i < edge->getNumberOfVertices(); i++)
    {
      for(size_t j = 0; j < 3; j++)
      {
        if(vertices[3 * i + j] < origin[j])
        {
          origin[j] = vertices[3 * i + j];
        }
      }
    }
    return origin;
  }
  if(IGeometry2D::Pointer geometry2d = std::dynamic_pointer_cast<IGeometry2D>(geometry))
  {
    float* vertices = geometry2d->getVertexPointer(0);
    FloatVec3Type origin(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max());

    for(MeshIndexType i = 0; i < geometry2d->getNumberOfVertices(); i++)
    {
      for(size_t j = 0; j < 3; j++)
      {
        if(vertices[3 * i + j] < origin[j])
        {
          origin[j] = vertices[3 * i + j];
        }
      }
    }
    return origin;
  }
  if(IGeometry3D::Pointer geometry3d = std::dynamic_pointer_cast<IGeometry3D>(geometry))
  {
    float* vertices = geometry3d->getVertexPointer(0);
    FloatVec3Type origin(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max());

    for(MeshIndexType i = 0; i < geometry3d->getNumberOfVertices(); i++)
    {
      for(size_t j = 0; j < 3; j++)
      {
        if(vertices[3 * i + j] < origin[j])
        {
          origin[j] = vertices[3 * i + j];
        }
      }
    }
    return origin;
  }
  FloatVec3Type origin(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max());
  return origin;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
FloatVec3Type extractCentroid(const IGeometry::Pointer& geometry)
{
  FloatVec3Type centroid(0.0f, 0.0f, 0.0f);
  if(ImageGeom::Pointer image = std::dynamic_pointer_cast<ImageGeom>(geometry))
  {
    SizeVec3Type dims = image->getDimensions();
    FloatVec3Type origin = image->getOrigin();
    FloatVec3Type res = image->getSpacing();

    centroid[0] = (static_cast<float>(dims[0]) * res[0] / 2.0f) + origin[0];
    centroid[1] = (static_cast<float>(dims[1]) * res[1] / 2.0f) + origin[1];
    centroid[2] = (static_cast<float>(dims[2]) * res[2] / 2.0f) + origin[2];
    return centroid;
  }
  if(RectGridGeom::Pointer rectGrid = std::dynamic_pointer_cast<RectGridGeom>(geometry))
  {
    FloatArrayType::Pointer xBounds = rectGrid->getXBounds();
    FloatArrayType::Pointer yBounds = rectGrid->getYBounds();
    FloatArrayType::Pointer zBounds = rectGrid->getZBounds();
    float min[3] = {std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max()};
    float max[3] = {std::numeric_limits<float>::min(), std::numeric_limits<float>::min(), std::numeric_limits<float>::min()};
    for(size_t i = 0; i < xBounds->getNumberOfTuples(); i++)
    {
      if(xBounds->getValue(i) < min[0])
      {
        min[0] = xBounds->getValue(i);
      }
      if(xBounds->getValue(i) > max[0])
      {
        max[0] = xBounds->getValue(i);
      }
    }
    for(size_t i = 0; i < yBounds->getNumberOfTuples(); i++)
    {
      if(yBounds->getValue(i) < min[1])
      {
        min[1] = yBounds->getValue(i);
      }
      if(yBounds->getValue(i) > max[1])
      {
        max[1] = yBounds->getValue(i);
      }
    }
    for(size_t i = 0; i < zBounds->getNumberOfTuples(); i++)
    {
      if(zBounds->getValue(i) < min[2])
      {
        min[2] = zBounds->getValue(i);
      }
      if(zBounds->getValue(i) > max[2])
      {
        max[2] = zBounds->getValue(i);
      }
    }
    centroid[0] = (max[0] - min[0]) / 2.0f;
    centroid[1] = (max[1] - min[1]) / 2.0f;
    centroid[2] = (max[2] - min[2]) / 2.0f;
    return centroid;
  }
  if(VertexGeom::Pointer vertex = std::dynamic_pointer_cast<VertexGeom>(geometry))
  {
    float* vertices = vertex->getVertexPointer(0);
    centroid[0] = 0.0f;
    centroid[1] = 0.0f;
    centroid[2] = 0.0f;
    for(MeshIndexType i = 0; i < vertex->getNumberOfVertices(); i++)
    {
      centroid[0] += vertices[3 * i + 0];
      centroid[1] += vertices[3 * i + 1];
      centroid[2] += vertices[3 * i + 2];
    }
    centroid[0] /= static_cast<float>(vertex->getNumberOfVertices());
    centroid[1] /= static_cast<float>(vertex->getNumberOfVertices());
    centroid[2] /= static_cast<float>(vertex->getNumberOfVertices());
    return centroid;
  }
  if(EdgeGeom::Pointer edge = std::dynamic_pointer_cast<EdgeGeom>(geometry))
  {
    float* vertices = edge->getVertexPointer(0);
    centroid[0] = 0.0f;
    centroid[1] = 0.0f;
    centroid[2] = 0.0f;
    for(MeshIndexType i = 0; i < edge->getNumberOfVertices(); i++)
    {
      centroid[0] += vertices[3 * i + 0];
      centroid[1] += vertices[3 * i + 1];
      centroid[2] += vertices[3 * i + 2];
    }
    centroid[0] /= static_cast<float>(edge->getNumberOfVertices());
    centroid[1] /= static_cast<float>(edge->getNumberOfVertices());
    centroid[2] /= static_cast<float>(edge->getNumberOfVertices());
    return centroid;
  }
  if(IGeometry2D::Pointer geometry2d = std::dynamic_pointer_cast<IGeometry2D>(geometry))
  {
    float* vertices = geometry2d->getVertexPointer(0);
    centroid[0] = 0.0f;
    centroid[1] = 0.0f;
    centroid[2] = 0.0f;
    for(MeshIndexType i = 0; i < geometry2d->getNumberOfVertices(); i++)
    {
      centroid[0] += vertices[3 * i + 0];
      centroid[1] += vertices[3 * i + 1];
      centroid[2] += vertices[3 * i + 2];
    }
    centroid[0] /= static_cast<float>(geometry2d->getNumberOfVertices());
    centroid[1] /= static_cast<float>(geometry2d->getNumberOfVertices());
    centroid[2] /= static_cast<float>(geometry2d->getNumberOfVertices());
    return centroid;
  }
  if(IGeometry3D::Pointer geometry3d = std::dynamic_pointer_cast<IGeometry3D>(geometry))
  {
    float* vertices = geometry3d->getVertexPointer(0);
    centroid[0] = 0.0f;
    centroid[1] = 0.0f;
    centroid[2] = 0.0f;
    for(MeshIndexType i = 0; i < geometry3d->getNumberOfVertices(); i++)
    {
      centroid[0] += vertices[3 * i + 0];
      centroid[1] += vertices[3 * i + 1];
      centroid[2] += vertices[3 * i + 2];
    }
    centroid[0] /= static_cast<float>(geometry3d->getNumberOfVertices());
    centroid[1] /= static_cast<float>(geometry3d->getNumberOfVertices());
    centroid[2] /= static_cast<float>(geometry3d->getNumberOfVertices());
    return centroid;
  }
  return centroid;
}

// -----------------------------------------------------------------------------
FloatVec6Type FindBoundingBox(const SharedVertexList & vertices)
{
  FloatVec6Type minMax;
 // float* vertices = vertex->getVertexPointer(0);
  minMax[0] = std::numeric_limits<float>::max();
  minMax[1] = std::numeric_limits<float>::max();
  minMax[2] = std::numeric_limits<float>::max();

  minMax[3] = std::numeric_limits<float>::min();
  minMax[4] = std::numeric_limits<float>::min();
  minMax[5] = std::numeric_limits<float>::min();

  for(MeshIndexType i = 0; i < vertices.getNumberOfTuples(); i++)
  {
    if( vertices[3 * i + 0] < minMax[0])
    {
      minMax[0] = vertices[3 * i + 0];
    }
    if( vertices[3 * i + 0] > minMax[3])
    {
      minMax[3] = vertices[3 * i + 0];
    }

    if( vertices[3 * i + 1] < minMax[1])
    {
      minMax[1] = vertices[3 * i + 1];
    }
    if( vertices[3 * i + 1] > minMax[4])
    {
      minMax[4] = vertices[3 * i + 1];
    }

    if( vertices[3 * i + 2] < minMax[2])
    {
      minMax[2] = vertices[3 * i + 2];
    }
    if( vertices[3 * i + 2] > minMax[5])
    {
      minMax[5] = vertices[3 * i + 2];
    }
  }

  return minMax;
}

// -----------------------------------------------------------------------------
FloatVec6Type extractBoundingBox(const IGeometry::Pointer& geometry)
{
  FloatVec6Type minMax(0.0f, 0.0f, 0.0f, 1.0f, 1.0f, 1.0f);

  if(ImageGeom::Pointer image = std::dynamic_pointer_cast<ImageGeom>(geometry))
  {
    SizeVec3Type dims = image->getDimensions();
    FloatVec3Type origin = image->getOrigin();
    FloatVec3Type res = image->getSpacing();

    // Min Coordinate
    minMax[0] = origin[0];
    minMax[1] = origin[1];
    minMax[2] = origin[2];
    // Max Coordinate
    minMax[3] = (static_cast<float>(dims[0]) * res[0] ) + origin[0];
    minMax[4] = (static_cast<float>(dims[1]) * res[1] ) + origin[1];
    minMax[5] = (static_cast<float>(dims[2]) * res[2] ) + origin[2];
    return minMax;
  }

  if(RectGridGeom::Pointer rectGrid = std::dynamic_pointer_cast<RectGridGeom>(geometry))
  {
    FloatArrayType& xBounds = *rectGrid->getXBounds();
    FloatArrayType& yBounds = *rectGrid->getYBounds();
    FloatArrayType& zBounds = *rectGrid->getZBounds();

    // Min Coordinate
    minMax[0] = xBounds[0];
    minMax[1] = yBounds[1];
    minMax[2] = zBounds[2];

    // Max Coordinate
    minMax[3] = xBounds[xBounds.size() - 1];
    minMax[4] = yBounds[xBounds.size() - 1];
    minMax[5] = zBounds[xBounds.size() - 1];

    return minMax;
  }
  if(VertexGeom::Pointer vertex = std::dynamic_pointer_cast<VertexGeom>(geometry))
  {
    return FindBoundingBox(*vertex->getVertices());
  }
  if(EdgeGeom::Pointer edge = std::dynamic_pointer_cast<EdgeGeom>(geometry))
  {
    return FindBoundingBox(*edge->getVertices());
  }
  if(IGeometry2D::Pointer geometry2d = std::dynamic_pointer_cast<IGeometry2D>(geometry))
  {
    return FindBoundingBox(*geometry2d->getVertices());
  }
  if(IGeometry3D::Pointer geometry3d = std::dynamic_pointer_cast<IGeometry3D>(geometry))
  {
    return FindBoundingBox(*geometry3d->getVertices());
  }
  return minMax;
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void translateGeometry(const IGeometry::Pointer& geometry, const FloatVec3Type& translation)
{
  if(ImageGeom::Pointer image = std::dynamic_pointer_cast<ImageGeom>(geometry))
  {
    FloatVec3Type origin = image->getOrigin();
    origin[0] += translation[0];
    origin[1] += translation[1];
    origin[2] += translation[2];
    image->setOrigin(origin);
    return;
  }
  if(RectGridGeom::Pointer rectGrid = std::dynamic_pointer_cast<RectGridGeom>(geometry))
  {
    float* xBounds = rectGrid->getXBounds()->getPointer(0);
    float* yBounds = rectGrid->getYBounds()->getPointer(0);
    float* zBounds = rectGrid->getZBounds()->getPointer(0);
    for(size_t i = 0; i < rectGrid->getXPoints(); i++)
    {
      xBounds[i] += translation[0];
    }
    for(size_t i = 0; i < rectGrid->getYPoints(); i++)
    {
      yBounds[i] += translation[1];
    }
    for(size_t i = 0; i < rectGrid->getZPoints(); i++)
    {
      zBounds[i] += translation[2];
    }
    return;
  }
  if(VertexGeom::Pointer vertex = std::dynamic_pointer_cast<VertexGeom>(geometry))
  {
    float* vertices = vertex->getVertexPointer(0);
    for(MeshIndexType i = 0; i < vertex->getNumberOfVertices(); i++)
    {
      vertices[3 * i + 0] += translation[0];
      vertices[3 * i + 1] += translation[1];
      vertices[3 * i + 2] += translation[2];
    }
    return;
  }
  if(EdgeGeom::Pointer edge = std::dynamic_pointer_cast<EdgeGeom>(geometry))
  {
    float* vertices = edge->getVertexPointer(0);
    for(MeshIndexType i = 0; i < edge->getNumberOfVertices(); i++)
    {
      vertices[3 * i + 0] += translation[0];
      vertices[3 * i + 1] += translation[1];
      vertices[3 * i + 2] += translation[2];
    }
    return;
  }
  if(IGeometry2D::Pointer geometry2d = std::dynamic_pointer_cast<IGeometry2D>(geometry))
  {
    float* vertices = geometry2d->getVertexPointer(0);
    for(MeshIndexType i = 0; i < geometry2d->getNumberOfVertices(); i++)
    {
      vertices[3 * i + 0] += translation[0];
      vertices[3 * i + 1] += translation[1];
      vertices[3 * i + 2] += translation[2];
    }
    return;
  }
  if(IGeometry3D::Pointer geometry3d = std::dynamic_pointer_cast<IGeometry3D>(geometry))
  {
    float* vertices = geometry3d->getVertexPointer(0);
    for(MeshIndexType i = 0; i < geometry3d->getNumberOfVertices(); i++)
    {
      vertices[3 * i + 0] += translation[0];
      vertices[3 * i + 1] += translation[1];
      vertices[3 * i + 2] += translation[2];
    }
    return;
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void AlignGeometries::execute()
{
  initialize();
  dataCheck();
  if(getErrorCode() < 0)
  {
    return;
  }

  IGeometry::Pointer moving = getDataContainerArray()->getDataContainer(m_MovingGeometry)->getGeometry();
  IGeometry::Pointer target = getDataContainerArray()->getDataContainer(m_TargetGeometry)->getGeometry();

  if(m_AlignmentType == 0)
  {
    FloatVec3Type movingOrigin = extractOrigin(moving);
    FloatVec3Type targetOrigin = extractOrigin(target);

    float translation[3] = {targetOrigin[0] - movingOrigin[0], targetOrigin[1] - movingOrigin[1], targetOrigin[2] - movingOrigin[2]};
    translateGeometry(moving, translation);
  }
  else if(m_AlignmentType == 1)
  {
    FloatVec3Type movingCentroid = extractCentroid(moving);
    FloatVec3Type targetCentroid = extractCentroid(target);

    float translation[3] = {targetCentroid[0] - movingCentroid[0], targetCentroid[0] - movingCentroid[0], targetCentroid[0] - movingCentroid[0]};
    translateGeometry(moving, translation);
  }
  else if(m_AlignmentType > 1 && m_AlignmentType < 8)
  {
    FloatVec6Type movingBoundingBox = extractBoundingBox(moving);
    FloatVec6Type targetBoundingBox = extractBoundingBox(target);

    if(m_AlignmentType == 2) // XY Min
    {
      float translation[3] = {0.0F, 0.0F, targetBoundingBox[2] - movingBoundingBox[2]}; // move in -Z direction
      translateGeometry(moving, translation);
    }
    else if(m_AlignmentType == 3) // XY Max
    {
      float translation[3] = {0.0F, 0.0F, targetBoundingBox[5] - movingBoundingBox[5]}; // move in +Z direction
      translateGeometry(moving, translation);
    }
    else if(m_AlignmentType == 4) // XZ Min
    {
      float translation[3] = {0.0F, targetBoundingBox[1] - movingBoundingBox[1], 0.0F}; // move in -Y direction
      translateGeometry(moving, translation);
    }
    else if(m_AlignmentType == 5) // XZ Max
    {
      float translation[3] = {0.0F, targetBoundingBox[4] - movingBoundingBox[4], 0.0F}; // move in +Y direction
      translateGeometry(moving, translation);
    }
    else if(m_AlignmentType == 6) // YZ Min
    {
      float translation[3] = { targetBoundingBox[0] - movingBoundingBox[0], 0.0F, 0.0F}; // move in -X direction
      translateGeometry(moving, translation);
    }
    else if(m_AlignmentType == 7) // YZ Max
    {
      float translation[3] = {targetBoundingBox[3] - movingBoundingBox[3], 0.0F, 0.0F}; // move in +X direction
      translateGeometry(moving, translation);
    }

  }
  else
  {
    QString ss = QObject::tr("Invalid selection for alignment type");
    setErrorCondition(-1, ss);
  }

  notifyStatusMessage("Complete");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AbstractFilter::Pointer AlignGeometries::newFilterInstance(bool copyFilterParameters) const
{
  AlignGeometries::Pointer filter = AlignGeometries::New();
  if(copyFilterParameters)
  {
    copyFilterParameterInstanceVariables(filter.get());
  }
  return filter;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString AlignGeometries::getCompiledLibraryName() const
{
  return DREAM3DReviewConstants::DREAM3DReviewBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString AlignGeometries::getBrandingString() const
{
  return "DREAM3DReview";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString AlignGeometries::getFilterVersion() const
{
  QString version;
  QTextStream vStream(&version);
  vStream << DREAM3DReview::Version::Major() << "." << DREAM3DReview::Version::Minor() << "." << DREAM3DReview::Version::Patch();
  return version;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString AlignGeometries::getGroupName() const
{
  return SIMPL::FilterGroups::ReconstructionFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString AlignGeometries::getSubGroupName() const
{
  return SIMPL::FilterSubGroups::AlignmentFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString AlignGeometries::getHumanLabel() const
{
  return "Align Geometries";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QUuid AlignGeometries::getUuid() const
{
  return QUuid("{ce1ee404-0336-536c-8aad-f9641c9458be}");
}

// -----------------------------------------------------------------------------
AlignGeometries::Pointer AlignGeometries::NullPointer()
{
  return Pointer(static_cast<Self*>(nullptr));
}

// -----------------------------------------------------------------------------
std::shared_ptr<AlignGeometries> AlignGeometries::New()
{
  struct make_shared_enabler : public AlignGeometries
  {
  };
  std::shared_ptr<make_shared_enabler> val = std::make_shared<make_shared_enabler>();
  val->setupFilterParameters();
  return val;
}

// -----------------------------------------------------------------------------
QString AlignGeometries::getNameOfClass() const
{
  return QString("AlignGeometries");
}

// -----------------------------------------------------------------------------
QString AlignGeometries::ClassName()
{
  return QString("AlignGeometries");
}

// -----------------------------------------------------------------------------
void AlignGeometries::setMovingGeometry(const DataArrayPath& value)
{
  m_MovingGeometry = value;
}

// -----------------------------------------------------------------------------
DataArrayPath AlignGeometries::getMovingGeometry() const
{
  return m_MovingGeometry;
}

// -----------------------------------------------------------------------------
void AlignGeometries::setTargetGeometry(const DataArrayPath& value)
{
  m_TargetGeometry = value;
}

// -----------------------------------------------------------------------------
DataArrayPath AlignGeometries::getTargetGeometry() const
{
  return m_TargetGeometry;
}

// -----------------------------------------------------------------------------
void AlignGeometries::setAlignmentType(int value)
{
  m_AlignmentType = value;
}

// -----------------------------------------------------------------------------
int AlignGeometries::getAlignmentType() const
{
  return m_AlignmentType;
}
