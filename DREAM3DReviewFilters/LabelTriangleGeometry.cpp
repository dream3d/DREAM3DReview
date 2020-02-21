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

#include "LabelTriangleGeometry.h"

#include <QtCore/QTextStream>

#include "SIMPLib/Common/Constants.h"

#include "SIMPLib/FilterParameters/AbstractFilterParametersReader.h"
#include "SIMPLib/FilterParameters/DataArraySelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/DataContainerSelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/StringFilterParameter.h"
#include "SIMPLib/Geometry/TriangleGeom.h"
#include "SIMPLib/Math/GeometryMath.h"
#include "SIMPLib/DataContainers/DataContainerArray.h"
#include "SIMPLib/DataContainers/DataContainer.h"

#include "DREAM3DReview/DREAM3DReviewConstants.h"
#include "DREAM3DReview/DREAM3DReviewVersion.h"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
LabelTriangleGeometry::LabelTriangleGeometry()
: m_CADDataContainerName("")
, m_TriangleAttributeMatrixName("")
, m_RegionIdArrayName("")
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
LabelTriangleGeometry::~LabelTriangleGeometry() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void LabelTriangleGeometry::setupFilterParameters()
{
  FilterParameterVectorType parameters;
  DataContainerSelectionFilterParameter::RequirementType dcsReq;
  IGeometry::Types geomTypes = {IGeometry::Type::Triangle};
  dcsReq.dcGeometryTypes = geomTypes;
  parameters.push_back(SIMPL_NEW_DC_SELECTION_FP("CAD Geometry", CADDataContainerName, FilterParameter::RequiredArray, LabelTriangleGeometry, dcsReq));
  parameters.push_back(SIMPL_NEW_STRING_FP("Triangle Attribute Matrix", TriangleAttributeMatrixName, FilterParameter::CreatedArray, LabelTriangleGeometry));
  parameters.push_back(SIMPL_NEW_STRING_FP("Region Ids", RegionIdArrayName, FilterParameter::CreatedArray, LabelTriangleGeometry));
  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void LabelTriangleGeometry::readFilterParameters(AbstractFilterParametersReader* reader, int index)
{
  reader->openFilterGroup(this, index);
  setCADDataContainerName(reader->readDataArrayPath("CADDataContainerName", getCADDataContainerName()));
  setTriangleAttributeMatrixName(reader->readString("TriangleAttributeMatrixName", getTriangleAttributeMatrixName()));
  setRegionIdArrayName(reader->readString("RegionIdArrayName", getRegionIdArrayName()));
  reader->closeFilterGroup();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void LabelTriangleGeometry::updateTriangleInstancePointers()
{
  clearErrorCode();
  clearWarningCode();

  if(nullptr != m_RegionIdPtr.lock().get()) /* Validate the Weak Pointer wraps a non-nullptr pointer to a DataArray<T> object */
  {
    m_RegionId = m_RegionIdPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void LabelTriangleGeometry::dataCheck()
{
  clearErrorCode();
  clearWarningCode();

  DataArrayPath tempPath;

  TriangleGeom::Pointer triangle = getDataContainerArray()->getPrereqGeometryFromDataContainer<TriangleGeom>(this, getCADDataContainerName());

  std::vector<size_t> tDims(1, 0);
  getDataContainerArray()->getDataContainer(getCADDataContainerName())->createNonPrereqAttributeMatrix(this, getTriangleAttributeMatrixName(), tDims, AttributeMatrix::Type::Face);

  if(getErrorCode() < 0)
  {
    return;
  }

  std::vector<size_t> cDims(1, 1);
  tempPath.update(getCADDataContainerName().getDataContainerName(), getTriangleAttributeMatrixName(), getRegionIdArrayName());
  m_RegionIdPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<int32_t>>(this, tempPath, 0,
                                                                                                                     cDims); /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
  if(nullptr != m_RegionIdPtr.lock().get()) /* Validate the Weak Pointer wraps a non-nullptr pointer to a DataArray<T> object */
  {
    m_RegionId = m_RegionIdPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */
  if(getErrorCode() < 0)
  {
    return;
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void LabelTriangleGeometry::preflight()
{
  // These are the REQUIRED lines of CODE to make sure the filter behaves correctly
  setInPreflight(true);              // Set the fact that we are preflighting.
  emit preflightAboutToExecute();    // Emit this signal so that other widgets can do one file update
  emit updateFilterParameters(this); // Emit this signal to have the widgets push their values down to the filter
  dataCheck();                       // Run our DataCheck to make sure everthing is setup correctly
  emit preflightExecuted();          // We are done preflighting this filter
  setInPreflight(false);             // Inform the system this filter is NOT in preflight mode anymore.
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void LabelTriangleGeometry::execute()
{
  clearErrorCode();
  clearWarningCode();
  dataCheck();
  if(getErrorCode() < 0)
  {
    return;
  }

  TriangleGeom::Pointer triangle = getDataContainerArray()->getDataContainer(getCADDataContainerName())->getGeometryAs<TriangleGeom>();

  MeshIndexType* tris = triangle->getTriPointer(0);
  float* triVerts = triangle->getVertexPointer(0);
  size_t numTris = triangle->getNumberOfTris();

  DataContainer::Pointer m = getDataContainerArray()->getDataContainer(getCADDataContainerName());

  std::vector<size_t> tDims(1, numTris);
  m->getAttributeMatrix(getTriangleAttributeMatrixName())->resizeAttributeArrays(tDims);
  updateTriangleInstancePointers();

  int check = triangle->findElementNeighbors();
  if(check < 0)
  {
    // FIXME: Should check this error code.....
  }
  ElementDynamicList::Pointer m_TriangleNeighbors = triangle->getElementNeighbors();

  size_t chunkSize = 1000;
  std::vector<int64_t> triList(chunkSize, -1);

  // first identify connected triangle sets as features
  size_t size = 0;
  int32_t regionCount = 1;
  for(size_t i = 0; i < numTris; i++)
  {
    if(m_RegionId[i] == 0)
    {
      m_RegionId[i] = regionCount;
      size = 0;
      triList[size] = i;
      size++;
      while(size > 0)
      {
        MeshIndexType tri = triList[size - 1];
        size -= 1;
        uint16_t tCount = m_TriangleNeighbors->getNumberOfElements(tri);
        MeshIndexType* data = m_TriangleNeighbors->getElementListPointer(tri);
        for(int j = 0; j < tCount; j++)
        {
          MeshIndexType neighTri = data[j];
          if(m_RegionId[neighTri] == 0)
          {
            m_RegionId[neighTri] = regionCount;
            triList[size] = neighTri;
            size++;
            if(size >= triList.size())
            {
              size = triList.size();
              triList.resize(size + chunkSize);
              for(std::vector<int64_t>::size_type k = size; k < triList.size(); ++k)
              {
                triList[k] = -1;
              }
            }
          }
        }
      }
      regionCount++;
    }
  }

  // next determine bounding boxes so we can see if any regions are within other regions
  std::vector<float> xMinList(regionCount, std::numeric_limits<float>::max());
  std::vector<float> yMinList(regionCount, std::numeric_limits<float>::max());
  std::vector<float> zMinList(regionCount, std::numeric_limits<float>::max());
  std::vector<float> xMaxList(regionCount, -std::numeric_limits<float>::max());
  std::vector<float> yMaxList(regionCount, -std::numeric_limits<float>::max());
  std::vector<float> zMaxList(regionCount, -std::numeric_limits<float>::max());
  for(MeshIndexType i = 0; i < numTris; i++)
  {
    int32_t regionId = m_RegionId[i];
    for(int j = 0; j < 3; j++)
    {
      MeshIndexType vert = tris[3 * i + j];
      if(xMinList[regionId] > triVerts[3 * vert + 0])
      {
        xMinList[regionId] = triVerts[3 * vert + 0];
      }
      if(yMinList[regionId] > triVerts[3 * vert + 1])
      {
        yMinList[regionId] = triVerts[3 * vert + 1];
      }
      if(zMinList[regionId] > triVerts[3 * vert + 2])
      {
        zMinList[regionId] = triVerts[3 * vert + 2];
      }
      if(xMaxList[regionId] < triVerts[3 * vert + 0])
      {
        xMaxList[regionId] = triVerts[3 * vert + 0];
      }
      if(yMaxList[regionId] < triVerts[3 * vert + 1])
      {
        yMaxList[regionId] = triVerts[3 * vert + 1];
      }
      if(zMaxList[regionId] < triVerts[3 * vert + 2])
      {
        zMaxList[regionId] = triVerts[3 * vert + 2];
      }
    }
  }

  std::vector<int32_t> newRegionIds(regionCount);
  std::vector<int32_t> contiguousRegionIds(regionCount);
  for(size_t i = 1; i < regionCount; i++)
  {
    newRegionIds[i] = i;
    for(size_t j = 1; j < regionCount; j++)
    {
      if(xMinList[i] > xMinList[j] && xMinList[i] < xMaxList[j] && xMaxList[i] > xMinList[j] && xMaxList[i] < xMaxList[j] && yMinList[i] > yMinList[j] && yMinList[i] < yMaxList[j] &&
         yMaxList[i] > yMinList[j] && yMaxList[i] < yMaxList[j] && zMinList[i] > zMinList[j] && zMinList[i] < zMaxList[j] && zMaxList[i] > zMinList[j] && zMaxList[i] < zMaxList[j])
      {
        newRegionIds[i] = j;
      }
    }
  }

  int32_t newRegionCount = 1;
  for(size_t i = 1; i < regionCount; i++)
  {
    if(newRegionIds[i] == i)
    {
      contiguousRegionIds[i] = newRegionCount;
      newRegionCount++;
    }
  }

  int32_t newRegionId = 0;
  for(size_t i = 0; i < numTris; i++)
  {
    int32_t regionId = m_RegionId[i];
    bool keepGoing = true;
    while(keepGoing)
    {
      if(newRegionIds[regionId] != regionId)
      {
        regionId = newRegionIds[regionId];
      }
      else
      {
        newRegionId = regionId;
        keepGoing = false;
      }
    }
    m_RegionId[i] = contiguousRegionIds[newRegionId];
  }

  notifyStatusMessage("Complete");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AbstractFilter::Pointer LabelTriangleGeometry::newFilterInstance(bool copyFilterParameters) const
{
  LabelTriangleGeometry::Pointer filter = LabelTriangleGeometry::New();
  if(copyFilterParameters)
  {
    copyFilterParameterInstanceVariables(filter.get());
  }
  return filter;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString LabelTriangleGeometry::getCompiledLibraryName() const
{
  return DREAM3DReviewConstants::DREAM3DReviewBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString LabelTriangleGeometry::getBrandingString() const
{
  return "DREAM3DReview";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString LabelTriangleGeometry::getFilterVersion() const
{
  QString version;
  QTextStream vStream(&version);
  vStream << DREAM3DReview::Version::Major() << "." << DREAM3DReview::Version::Minor() << "." << DREAM3DReview::Version::Patch();
  return version;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString LabelTriangleGeometry::getGroupName() const
{
  return SIMPL::FilterGroups::ReconstructionFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString LabelTriangleGeometry::getSubGroupName() const
{
  return SIMPL::FilterSubGroups::GroupingFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString LabelTriangleGeometry::getHumanLabel() const
{
  return "Label CAD Geometry";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QUuid LabelTriangleGeometry::getUuid() const
{
  return QUuid("{a250a228-3b6b-5b37-a6e4-8687484f04c4}");
}

// -----------------------------------------------------------------------------
LabelTriangleGeometry::Pointer LabelTriangleGeometry::NullPointer()
{
  return Pointer(static_cast<Self*>(nullptr));
}

// -----------------------------------------------------------------------------
std::shared_ptr<LabelTriangleGeometry> LabelTriangleGeometry::New()
{
  struct make_shared_enabler : public LabelTriangleGeometry
  {
  };
  std::shared_ptr<make_shared_enabler> val = std::make_shared<make_shared_enabler>();
  val->setupFilterParameters();
  return val;
}

// -----------------------------------------------------------------------------
QString LabelTriangleGeometry::getNameOfClass() const
{
  return QString("LabelTriangleGeometry");
}

// -----------------------------------------------------------------------------
QString LabelTriangleGeometry::ClassName()
{
  return QString("LabelTriangleGeometry");
}

// -----------------------------------------------------------------------------
void LabelTriangleGeometry::setCADDataContainerName(const DataArrayPath& value)
{
  m_CADDataContainerName = value;
}

// -----------------------------------------------------------------------------
DataArrayPath LabelTriangleGeometry::getCADDataContainerName() const
{
  return m_CADDataContainerName;
}

// -----------------------------------------------------------------------------
void LabelTriangleGeometry::setTriangleAttributeMatrixName(const QString& value)
{
  m_TriangleAttributeMatrixName = value;
}

// -----------------------------------------------------------------------------
QString LabelTriangleGeometry::getTriangleAttributeMatrixName() const
{
  return m_TriangleAttributeMatrixName;
}

// -----------------------------------------------------------------------------
void LabelTriangleGeometry::setRegionIdArrayName(const QString& value)
{
  m_RegionIdArrayName = value;
}

// -----------------------------------------------------------------------------
QString LabelTriangleGeometry::getRegionIdArrayName() const
{
  return m_RegionIdArrayName;
}
