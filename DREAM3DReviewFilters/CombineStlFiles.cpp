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
#include "CombineStlFiles.h"

#include <cstring>

#include "SIMPLib/Common/Constants.h"
#include "SIMPLib/FilterParameters/InputPathFilterParameter.h"
#include "SIMPLib/FilterParameters/SeparatorFilterParameter.h"
#include "SIMPLib/FilterParameters/StringFilterParameter.h"
#include "SIMPLib/Filtering/FilterFactory.hpp"
#include "SIMPLib/Filtering/FilterManager.h"
#include "SIMPLib/Geometry/TriangleGeom.h"

#include "DREAM3DReview/DREAM3DReviewConstants.h"
#include "DREAM3DReview/DREAM3DReviewVersion.h"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
CombineStlFiles::CombineStlFiles()
: m_StlFilesPath("")
, m_TriangleDataContainerName(SIMPL::Defaults::TriangleDataContainerName)
, m_FaceAttributeMatrixName(SIMPL::Defaults::FaceAttributeMatrixName)
, m_FaceNormalsArrayName(SIMPL::FaceData::SurfaceMeshFaceNormals)
, m_FaceNormals(nullptr)
, m_FileList(QFileInfoList())
{
  initialize();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
CombineStlFiles::~CombineStlFiles() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void CombineStlFiles::initialize()
{
  clearErrorCode();
  clearWarningCode();
  setCancel(false);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void CombineStlFiles::setupFilterParameters()
{
  FilterParameterVectorType parameters;
  parameters.push_back(SIMPL_NEW_INPUT_PATH_FP("Path to STL Files", StlFilesPath, FilterParameter::Parameter, CombineStlFiles, "", ""));
  parameters.push_back(SIMPL_NEW_STRING_FP("Data Container", TriangleDataContainerName, FilterParameter::CreatedArray, CombineStlFiles));
  parameters.push_back(SeparatorFilterParameter::New("Face Data", FilterParameter::CreatedArray));
  parameters.push_back(SIMPL_NEW_STRING_FP("Face Attribute Matrix", FaceAttributeMatrixName, FilterParameter::CreatedArray, CombineStlFiles));
  parameters.push_back(SIMPL_NEW_STRING_FP("Face Normals", FaceNormalsArrayName, FilterParameter::CreatedArray, CombineStlFiles));
  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void CombineStlFiles::dataCheck()
{
  clearErrorCode();
  clearWarningCode();

  FilterManager* fm = FilterManager::Instance();
  IFilterFactory::Pointer factory = fm->getFactoryFromClassName("ReadStlFile");
  if(!factory)
  {
    QString ss = QObject::tr("Combine STL Files requires the Read STL File filter to be loaded");
    setErrorCondition(-1, ss);
    return;
  }

  QFileInfo fi(getStlFilesPath());

  if(getStlFilesPath().isEmpty())
  {
    QString ss = QObject::tr("The input STL directory must be set");
    setErrorCondition(-387, ss);
  }
  else if(!fi.exists())
  {
    QString ss = QObject::tr("The input STL directory does not exist");
    setErrorCondition(-388, ss);
  }

  if(getErrorCode() < 0)
  {
    return;
  }

  m_FileList.clear();

  QDir directory(getStlFilesPath());
  m_FileList = directory.entryInfoList(QDir::Files | QDir::NoDotAndDotDot | QDir::NoSymLinks);

  size_t s = m_FileList.size();

  m_FileList.erase(std::remove_if(std::begin(m_FileList), std::end(m_FileList), [](decltype(*std::begin(m_FileList))& val) -> bool { return (val.suffix().toLower() != "stl"); }),
                   std::end(m_FileList));

  s = m_FileList.size();

  if(m_FileList.isEmpty())
  {
    QString ss = QObject::tr("No STL files were found in the selected directory");
    setErrorCondition(-1, ss);
    return;
  }

  DataContainer::Pointer sm = getDataContainerArray()->createNonPrereqDataContainer<AbstractFilter>(this, getTriangleDataContainerName());
  if(getErrorCode() < 0)
  {
    return;
  }

  SharedVertexList::Pointer sharedVertList = TriangleGeom::CreateSharedVertexList(0);
  TriangleGeom::Pointer triangleGeom = TriangleGeom::CreateGeometry(0, sharedVertList, SIMPL::Geometry::TriangleGeometry, !getInPreflight());
  sm->setGeometry(triangleGeom);

  QVector<size_t> tDims(1, 0);
  sm->createNonPrereqAttributeMatrix(this, getFaceAttributeMatrixName(), tDims, AttributeMatrix::Type::Face);

  QVector<size_t> cDims(1, 3);
  DataArrayPath path(getTriangleDataContainerName(), getFaceAttributeMatrixName(), getFaceNormalsArrayName());
  m_FaceNormalsPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<double>, AbstractFilter, double>(
      this, path, 0, cDims);                   /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
  if(nullptr != m_FaceNormalsPtr.lock().get()) /* Validate the Weak Pointer wraps a non-nullptr pointer to a DataArray<T> object */
  {
    m_FaceNormals = m_FaceNormalsPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void CombineStlFiles::preflight()
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
void CombineStlFiles::execute()
{
  initialize();
  dataCheck();
  if(getErrorCode() < 0)
  {
    return;
  }

  DataContainerArray::Pointer dca = DataContainerArray::New();

  FilterManager* fm = FilterManager::Instance();
  IFilterFactory::Pointer factory = fm->getFactoryFromClassName("ReadStlFile");
  AbstractFilter::Pointer reader = factory->create();
  reader->setDataContainerArray(dca);

  for(auto&& file : m_FileList)
  {
    QVariant var;
    var.setValue(file.canonicalFilePath());
    reader->setProperty("StlFilePath", var);
    var.setValue(file.baseName());
    reader->setProperty("SurfaceMeshDataContainerName", var);
    var.setValue(SIMPL::Defaults::FaceAttributeMatrixName);
    reader->setProperty("FaceAttributeMatrixName", var);
    var.setValue(SIMPL::FaceData::SurfaceMeshFaceNormals);
    reader->setProperty("FaceNormalsArrayName", var);
    reader->execute();
    if(reader->getErrorCode() < 0)
    {
      QString ss = QObject::tr("Error reading STL file: %1").arg(file.fileName());
      setErrorCondition(reader->getErrorCode(), ss);
      return;
    }
  }

  std::list<TriangleGeom::Pointer> stlGeoms;
  std::list<DoubleArrayType::Pointer> faceNormals;
  auto containers = dca->getDataContainers();
  DataArrayPath path;

  for(auto&& container : containers)
  {
    stlGeoms.push_back(container->getGeometryAs<TriangleGeom>());
    path.update(container->getName(), SIMPL::Defaults::FaceAttributeMatrixName, "");
    faceNormals.push_back(container->getAttributeMatrix(path)->getAttributeArrayAs<DoubleArrayType>(SIMPL::FaceData::SurfaceMeshFaceNormals));
  }

  int64_t totalTriangles = 0;
  int64_t totalVertices = 0;

  for(auto&& geom : stlGeoms)
  {
    totalTriangles += geom->getNumberOfTris();
    totalVertices += geom->getNumberOfVertices();
  }

  TriangleGeom::Pointer combined = getDataContainerArray()->getDataContainer(m_TriangleDataContainerName)->getGeometryAs<TriangleGeom>();
  AttributeMatrix::Pointer faceAttrmat = getDataContainerArray()->getAttributeMatrix(DataArrayPath(m_TriangleDataContainerName, m_FaceAttributeMatrixName, ""));
  QVector<size_t> tDims(1, totalTriangles);
  combined->resizeTriList(totalTriangles);
  combined->resizeVertexList(totalVertices);
  faceAttrmat->resizeAttributeArrays(tDims);
  int64_t triOffset = 0;
  int64_t vertOffset = 0;
  int64_t triCounter = 0;
  int64_t* tris = combined->getTriPointer(0);
  float* verts = combined->getVertexPointer(0);

  for(auto&& geom : stlGeoms)
  {
    int64_t* curTris = geom->getTriPointer(0);
    for(int64_t t = 0; t < geom->getNumberOfTris(); t++)
    {
      curTris[3 * t + 0] += triCounter;
      curTris[3 * t + 1] += triCounter;
      curTris[3 * t + 2] += triCounter;
    }
    triCounter += geom->getNumberOfVertices();
    float* curVerts = geom->getVertexPointer(0);
    std::memcpy(tris + triOffset, curTris, geom->getNumberOfTris() * 3 * sizeof(int64_t));
    std::memcpy(verts + vertOffset, curVerts, geom->getNumberOfVertices() * 3 * sizeof(float));
    triOffset += geom->getNumberOfTris() * 3;
    vertOffset += geom->getNumberOfVertices() * 3;
  }

  size_t faceNormalsOffset = 0;
  m_FaceNormals = faceAttrmat->getAttributeArrayAs<DoubleArrayType>(m_FaceNormalsArrayName)->getPointer(0);

  for(auto&& dataArray : faceNormals)
  {
    std::memcpy(m_FaceNormals + faceNormalsOffset, dataArray->getPointer(0), dataArray->getSize() * sizeof(double));
    faceNormalsOffset += dataArray->getSize();
  }

  notifyStatusMessage("Complete");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AbstractFilter::Pointer CombineStlFiles::newFilterInstance(bool copyFilterParameters) const
{
  CombineStlFiles::Pointer filter = CombineStlFiles::New();
  if(copyFilterParameters)
  {
    copyFilterParameterInstanceVariables(filter.get());
  }
  return filter;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString CombineStlFiles::getCompiledLibraryName() const
{
  return DREAM3DReviewConstants::DREAM3DReviewBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString CombineStlFiles::getBrandingString() const
{
  return "DREAM3DReview";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString CombineStlFiles::getFilterVersion() const
{
  QString version;
  QTextStream vStream(&version);
  vStream << DREAM3DReview::Version::Major() << "." << DREAM3DReview::Version::Minor() << "." << DREAM3DReview::Version::Patch();
  return version;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString CombineStlFiles::getGroupName() const
{
  return SIMPL::FilterGroups::IOFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString CombineStlFiles::getSubGroupName() const
{
  return SIMPL::FilterSubGroups::InputFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString CombineStlFiles::getHumanLabel() const
{
  return "Combine STL Files";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QUuid CombineStlFiles::getUuid()
{
  return QUuid("{71d46128-1d2d-58fd-9924-1714695768c3}");
}