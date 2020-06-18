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

#include "DownsampleVertexGeometry.h"

#include <random>

#include "SIMPLib/Common/Constants.h"
#include "SIMPLib/Common/TemplateHelpers.h"
#include "SIMPLib/DataContainers/DataContainer.h"
#include "SIMPLib/DataContainers/DataContainerArray.h"
#include "SIMPLib/FilterParameters/AbstractFilterParametersReader.h"
#include "SIMPLib/FilterParameters/AttributeMatrixSelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/FloatFilterParameter.h"
#include "SIMPLib/FilterParameters/IntFilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedChoicesFilterParameter.h"
#include "SIMPLib/FilterParameters/SeparatorFilterParameter.h"
#include "SIMPLib/Geometry/ImageGeom.h"
#include "SIMPLib/Geometry/VertexGeom.h"

#include "DREAM3DReview/DREAM3DReviewConstants.h"
#include "DREAM3DReview/DREAM3DReviewVersion.h"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
DownsampleVertexGeometry::DownsampleVertexGeometry()
{
  initialize();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
DownsampleVertexGeometry::~DownsampleVertexGeometry() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void DownsampleVertexGeometry::setupFilterParameters()
{
  FilterParameterVectorType parameters;
  {
    LinkedChoicesFilterParameter::Pointer parameter = LinkedChoicesFilterParameter::New();
    parameter->setHumanLabel("Downsample Type");
    parameter->setPropertyName("DownsampleType");
    parameter->setSetterCallback(SIMPL_BIND_SETTER(DownsampleVertexGeometry, this, DownsampleType));
    parameter->setGetterCallback(SIMPL_BIND_GETTER(DownsampleVertexGeometry, this, DownsampleType));
    QVector<QString> choices;
    choices.push_back("Remove Every Nth Point");
    choices.push_back("Remove a Fixed Random Fraction of Points");
    choices.push_back("Downsample the Geometry on a Grid");
    parameter->setChoices(choices);
    QStringList linkedProps = {"DecimationFreq", "DecimationFraction", "GridResolution"};
    parameter->setLinkedProperties(linkedProps);
    parameter->setEditable(false);
    parameter->setCategory(FilterParameter::Parameter);
    parameters.push_back(parameter);
  }
  parameters.push_back(SIMPL_NEW_INTEGER_FP("Decimation Frequency", DecimationFreq, FilterParameter::Parameter, DownsampleVertexGeometry, 0));
  parameters.push_back(SIMPL_NEW_FLOAT_FP("Fraction to Remove", DecimationFraction, FilterParameter::Parameter, DownsampleVertexGeometry, 1));
  parameters.push_back(SIMPL_NEW_FLOAT_VEC3_FP("Grid Resolution", GridResolution, FilterParameter::Parameter, DownsampleVertexGeometry, 2));
  parameters.push_back(SeparatorFilterParameter::New("Vertex Data", FilterParameter::RequiredArray));
  {
    AttributeMatrixSelectionFilterParameter::RequirementType req;
    IGeometry::Types geomTypes = {IGeometry::Type::Vertex};
    AttributeMatrix::Types amTypes = {AttributeMatrix::Type::Vertex};
    req.dcGeometryTypes = geomTypes;
    req.amTypes = amTypes;
    parameters.push_back(SIMPL_NEW_AM_SELECTION_FP("Vertex Attribute Matrix", VertexAttrMatPath, FilterParameter::RequiredArray, DownsampleVertexGeometry, req));
  }
  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void DownsampleVertexGeometry::readFilterParameters(AbstractFilterParametersReader* reader, int index)
{
  reader->openFilterGroup(this, index);
  setVertexAttrMatPath(reader->readDataArrayPath("VertexAttrMatPath", getVertexAttrMatPath()));
  setDecimationFreq(reader->readValue("DecimationFreq", getDecimationFreq()));
  reader->closeFilterGroup();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void DownsampleVertexGeometry::initialize()
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void DownsampleVertexGeometry::dataCheck()
{
  clearErrorCode();
  clearWarningCode();

  VertexGeom::Pointer vertices = getDataContainerArray()->getPrereqGeometryFromDataContainer<VertexGeom>(this, getVertexAttrMatPath().getDataContainerName());
  
  if(getErrorCode() < 0)
  {
    return;
  }

  int64_t numVerts = vertices->getNumberOfVertices();

  if(getDownsampleType() == 0)
  {
    if(getDecimationFreq() <= 1)
    {
      QString ss = QObject::tr("Decimation frequency must be greater than 1 or else all points would be removed");
      setErrorCondition(-11000, ss);
    }

    if(getDecimationFreq() > numVerts)
    {
      QString ss = QObject::tr("Decimation frequency is larger than the total number of verteices in the Vertex Geometry, so no points would be removed");
      setErrorCondition(-11001, ss);
    }
  }

  if(getDownsampleType() == 1)
  {
    if(getDecimationFraction() <= 0.0f || getDecimationFraction() >= 1.0f)
    {
      QString ss = QObject::tr("Decimation fraction must be in the interval (0, 1)");
      setErrorCondition(-11002, ss);
    }
  }

  if(getDownsampleType() == 2)
  {
    if(getGridResolution()[0] <= 0.0f || getGridResolution()[1] <= 0.0f || getGridResolution()[2] <= 0.0f)
    {
      QString ss = QObject::tr("Grid resolutions must be greater than zero");
      setErrorCondition(-11001, ss);
    }
  }

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
void DownsampleVertexGeometry::removeNthPoint()
{
  VertexGeom::Pointer vertices = getDataContainerArray()->getDataContainer(getVertexAttrMatPath().getDataContainerName())->getGeometryAs<VertexGeom>();
  AttributeMatrix::Pointer attrMat = getDataContainerArray()->getAttributeMatrix(getVertexAttrMatPath());

  int64_t numVerts = vertices->getNumberOfVertices();
  float* vertex = vertices->getVertexPointer(0);

  int64_t progIncrement = numVerts / 100;
  int64_t prog = 1;
  int64_t progressInt = 0;
  std::vector<size_t> removeList;
  int32_t freqCounter = 1;
  size_t currentSlot = 0;

  for(int64_t i = 0; i < numVerts; i++)
  {
    if(freqCounter == m_DecimationFreq)
    {
      removeList.push_back(i);
      freqCounter = 1;
    }
    else
    {
      freqCounter += 1;
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
void DownsampleVertexGeometry::removeFractionPoints()
{
  VertexGeom::Pointer vertices = getDataContainerArray()->getDataContainer(getVertexAttrMatPath().getDataContainerName())->getGeometryAs<VertexGeom>();
  AttributeMatrix::Pointer attrMat = getDataContainerArray()->getAttributeMatrix(getVertexAttrMatPath());

  int64_t numVerts = vertices->getNumberOfVertices();
  float* vertex = vertices->getVertexPointer(0);

  size_t rangeMin = 0;
  size_t rangeMax = numVerts - 1;
  std::mt19937_64::result_type seed = static_cast<std::mt19937_64::result_type>(std::chrono::steady_clock::now().time_since_epoch().count());
  std::mt19937_64 gen(seed);
  std::uniform_int_distribution<size_t> dist(rangeMin, rangeMax);

  int64_t numDecimatedPoints = static_cast<int64_t>(m_DecimationFraction * numVerts);
  std::vector<size_t> removeList;
  BoolArrayType::Pointer mask = BoolArrayType::CreateArray(numVerts, std::string("_INTERNAL_USE_ONLY_Mask"), true);
  mask->initializeWithValue(true);
  bool* maskPtr = mask->getPointer(0);
  int64_t counter = 0;

  int64_t progIncrement = numDecimatedPoints / 100;
  int64_t prog = 1;
  int64_t progressInt = 0;

  while(counter < numDecimatedPoints)
  {
    size_t index = dist(gen);
    if(maskPtr[index])
    {
      maskPtr[index] = false;
      removeList.push_back(index);
      counter++;
    }

    if(counter > prog)
    {
      progressInt = static_cast<int64_t>((static_cast<float>(counter) / numDecimatedPoints) * 100.0f);
      QString ss = QObject::tr("Decimating Point Cloud || %1% Complete").arg(progressInt);
      notifyStatusMessage(ss);
      prog = prog + progIncrement;
    }
  }

  std::vector<float> tmpVerts;
  for(int64_t i = 0; i < numVerts; i++)
  {
    if(maskPtr[i])
    {
      tmpVerts.push_back(vertex[3 * i + 0]);
      tmpVerts.push_back(vertex[3 * i + 1]);
      tmpVerts.push_back(vertex[3 * i + 2]);
    }
  }

  QString ss = QObject::tr("Resizing Vertex Geometry...");
  notifyStatusMessage(ss);

  vertices->resizeVertexList(numVerts - removeList.size());
  vertex = vertices->getVertexPointer(0);
  std::memcpy(vertex, tmpVerts.data(), sizeof(float) * tmpVerts.size());
  std::sort(std::begin(removeList), std::end(removeList));

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
template <typename T>
void downsampleDataByAveraging(IDataArray::Pointer source, IDataArray::Pointer dest, std::vector<int64_t>& indices, size_t destIdx)
{
  typename DataArray<T>::Pointer sourcePtr = std::dynamic_pointer_cast<DataArray<T>>(source);
  T* sourceData = sourcePtr->getPointer(0);
  typename DataArray<T>::Pointer destPtr = std::dynamic_pointer_cast<DataArray<T>>(dest);
  T* destData = destPtr->getPointer(0);

  double accumulator = 0.0;
  for(int64_t idx : indices)
  {
    accumulator += sourceData[idx];
  }
  accumulator /= static_cast<double>(indices.size());

  destData[destIdx] = static_cast<T>(accumulator);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
template <>
void downsampleDataByAveraging<bool>(IDataArray::Pointer source, IDataArray::Pointer dest, std::vector<int64_t>& indices, size_t destIdx)
{
  typename DataArray<bool>::Pointer sourcePtr = std::dynamic_pointer_cast<DataArray<bool>>(source);
  bool* sourceData = sourcePtr->getPointer(0);
  typename DataArray<bool>::Pointer destPtr = std::dynamic_pointer_cast<DataArray<bool>>(dest);
  bool* destData = destPtr->getPointer(0);
  
  size_t counter = 0;
  bool value = false;
  for(int64_t idx : indices)
  {
    if(sourceData[idx])
    {
      counter++;
    }
  }
  if(static_cast<float>(counter) >= 0.5f * static_cast<float>(indices.size()))
  {
    value = true;
  }

  destData[destIdx] = value;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void DownsampleVertexGeometry::gridDownsample()
{
  float inverseResolution[3] = {1.0f / m_GridResolution[0], 1.0f / m_GridResolution[1], 1.0f / m_GridResolution[2]};

  VertexGeom::Pointer vertices = getDataContainerArray()->getDataContainer(getVertexAttrMatPath().getDataContainerName())->getGeometryAs<VertexGeom>();
  AttributeMatrix::Pointer attrMat = getDataContainerArray()->getAttributeMatrix(getVertexAttrMatPath());
  float* verts = vertices->getVertexPointer(0);

  ImageGeom::Pointer grid = ImageGeom::CreateGeometry(SIMPL::Geometry::ImageGeometry);
  grid->setSpacing(m_GridResolution[0], m_GridResolution[1], m_GridResolution[2]);

  int64_t numVerts = vertices->getNumberOfVertices();
  float* vertex = vertices->getVertexPointer(0);

  std::vector<float> meshMinExtents;
  std::vector<float> meshMaxExtents;

  for(size_t i = 0; i < 3; i++)
  {
    meshMaxExtents.push_back(std::numeric_limits<float>::lowest());
    meshMinExtents.push_back(std::numeric_limits<float>::max());
  }

  for(int64_t i = 0; i < numVerts; i++)
  {
    if(vertex[3 * i] > meshMaxExtents[0])
    {
      meshMaxExtents[0] = vertex[3 * i];
    }
    if(vertex[3 * i + 1] > meshMaxExtents[1])
    {
      meshMaxExtents[1] = vertex[3 * i + 1];
    }
    if(vertex[3 * i + 2] > meshMaxExtents[2])
    {
      meshMaxExtents[2] = vertex[3 * i + 2];
    }
    if(vertex[3 * i] < meshMinExtents[0])
    {
      meshMinExtents[0] = vertex[3 * i];
    }
    if(vertex[3 * i + 1] < meshMinExtents[1])
    {
      meshMinExtents[1] = vertex[3 * i + 1];
    }
    if(vertex[3 * i + 2] < meshMinExtents[2])
    {
      meshMinExtents[2] = vertex[3 * i + 2];
    }
  }

  for(auto i = 0; i < 3; i++)
  {
    meshMinExtents[i] -= (inverseResolution[i] / 2.0f);
    meshMaxExtents[i] += (inverseResolution[i] / 2.0f);
  }

  int64_t bboxMin[3] = {0, 0, 0};
  int64_t bboxMax[3] = {0, 0, 0};

  bboxMin[0] = static_cast<int64_t>(std::floor(meshMinExtents[0] * inverseResolution[0]));
  bboxMin[1] = static_cast<int64_t>(std::floor(meshMinExtents[1] * inverseResolution[1]));
  bboxMin[2] = static_cast<int64_t>(std::floor(meshMinExtents[2] * inverseResolution[2]));

  bboxMax[0] = static_cast<int64_t>(std::floor(meshMaxExtents[0] * inverseResolution[0]));
  bboxMax[1] = static_cast<int64_t>(std::floor(meshMaxExtents[1] * inverseResolution[1]));
  bboxMax[2] = static_cast<int64_t>(std::floor(meshMaxExtents[2] * inverseResolution[2]));

  int64_t dims[3] = {bboxMax[0] - bboxMin[0] + 1, bboxMax[1] - bboxMin[1] + 1, bboxMax[2] - bboxMin[2] + 1};
  grid->setDimensions(dims[0], dims[1], dims[2]);

  if(static_cast<int64_t>(grid->getNumberOfElements()) >= numVerts)
  {
    QString ss = QObject::tr("The total number of voxels in the sampling grid exceeds the number of vertices in the Vertex Geometry, so no points will be removed");
    setErrorCondition(-1, ss);
    return;
  }

  int64_t multiplier[3] = {1, static_cast<int64_t>(grid->getXPoints()), static_cast<int64_t>(grid->getXPoints() * grid->getYPoints())};
  std::vector<std::vector<int64_t>> vertsInVoxels(grid->getNumberOfElements());

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
  progIncrement = (dims[0] * dims[1] * dims[2]) / 100;
  prog = 1;
  progressInt = 0;
  int64_t counter = 0;
  int64_t vertCounter = 0;
  float xAvg = 0.0f;
  float yAvg = 0.0f;
  float zAvg = 0.0f;

  AttributeMatrix::Pointer downsampledData = attrMat->deepCopy();
  QList<QString> headers = downsampledData->getAttributeArrayNames();
  for(QList<QString>::iterator iter = headers.begin(); iter != headers.end(); ++iter)
  {
    IDataArray::Pointer p = downsampledData->getAttributeArray(*iter);
    QString type = p->getTypeAsString();
    if(type.compare("NeighborList<T>") == 0)
    {
      attrMat->removeAttributeArray(*iter);
    }
  }

  size_t globalIndex = 0;

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

        QList<QString> headers = downsampledData->getAttributeArrayNames();
        for(QList<QString>::iterator iter = headers.begin(); iter != headers.end(); ++iter)
        {
          if(getCancel())
          {
            return;
          }
          IDataArray::Pointer source = attrMat->getAttributeArray(*iter);
          IDataArray::Pointer dest = downsampledData->getAttributeArray(*iter);
          EXECUTE_FUNCTION_TEMPLATE(this, downsampleDataByAveraging, source, source, dest, vertsInVoxels[index], globalIndex)
        }
        globalIndex++;

        if(counter > prog)
        {
          progressInt = static_cast<int64_t>((static_cast<float>(counter) / (dims[0] * dims[1] * dims[2])) * 100.0f);
          QString ss = QObject::tr("Performing Grid Downsampling || %1% Complete").arg(progressInt);
          notifyStatusMessage(ss);
          prog = prog + progIncrement;
        }
        counter++;
      }
    }
  }

  std::vector<size_t> tDims = {globalIndex};
  downsampledData->setTupleDimensions(tDims);

  DataContainer::Pointer dc = getDataContainerArray()->getDataContainer(getVertexAttrMatPath().getDataContainerName());
  dc->removeAttributeMatrix(attrMat->getName());
  dc->addOrReplaceAttributeMatrix(downsampledData);

  vertices->resizeVertexList(tmpVerts.size() / 3);
  float* gridVerts = vertices->getVertexPointer(0);
  std::memcpy(gridVerts, tmpVerts.data(), vertices->getNumberOfVertices() * 3 * sizeof(float));
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void DownsampleVertexGeometry::execute()
{
  clearErrorCode();
  clearWarningCode();
  dataCheck();
  if(getErrorCode() < 0)
  {
    return;
  }

  switch(m_DownsampleType)
  {
  case 0:
    removeNthPoint();
    break;
  case 1:
    removeFractionPoints();
    break;
  case 2: 
    gridDownsample();
    break;
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AbstractFilter::Pointer DownsampleVertexGeometry::newFilterInstance(bool copyFilterParameters) const
{
  DownsampleVertexGeometry::Pointer filter = DownsampleVertexGeometry::New();
  if(copyFilterParameters)
  {
    copyFilterParameterInstanceVariables(filter.get());
  }
  return filter;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString DownsampleVertexGeometry::getCompiledLibraryName() const
{
  return DREAM3DReviewConstants::DREAM3DReviewBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString DownsampleVertexGeometry::getBrandingString() const
{
  return "DREAM3DReview";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString DownsampleVertexGeometry::getFilterVersion() const
{
  QString version;
  QTextStream vStream(&version);
  vStream << DREAM3DReview::Version::Major() << "." << DREAM3DReview::Version::Minor() << "." << DREAM3DReview::Version::Patch();
  return version;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString DownsampleVertexGeometry::getGroupName() const
{
  return SIMPL::FilterGroups::SamplingFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString DownsampleVertexGeometry::getSubGroupName() const
{
  return SIMPL::FilterSubGroups::GeometryFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString DownsampleVertexGeometry::getHumanLabel() const
{
  return "Downsample Vertex Geometry";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QUuid DownsampleVertexGeometry::getUuid() const
{
  return QUuid("{5cbd9d8e-e2eb-59e7-be63-6ab9deeed8d2}");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
DownsampleVertexGeometry::Pointer DownsampleVertexGeometry::NullPointer()
{
  return Pointer(static_cast<Self*>(nullptr));
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
std::shared_ptr<DownsampleVertexGeometry> DownsampleVertexGeometry::New()
{
  struct make_shared_enabler : public DownsampleVertexGeometry
  {
  };
  std::shared_ptr<make_shared_enabler> val = std::make_shared<make_shared_enabler>();
  val->setupFilterParameters();
  return val;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString DownsampleVertexGeometry::getNameOfClass() const
{
  return QString("DownsampleVertexGeometry");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString DownsampleVertexGeometry::ClassName()
{
  return QString("DownsampleVertexGeometry");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void DownsampleVertexGeometry::setVertexAttrMatPath(const DataArrayPath& value)
{
  m_VertexAttrMatPath = value;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
DataArrayPath DownsampleVertexGeometry::getVertexAttrMatPath() const
{
  return m_VertexAttrMatPath;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void DownsampleVertexGeometry::setDownsampleType(const int& value)
{
  m_DownsampleType = value;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int DownsampleVertexGeometry::getDownsampleType() const
{
  return m_DownsampleType;
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void DownsampleVertexGeometry::setDecimationFreq(const int& value)
{
  m_DecimationFreq = value;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int DownsampleVertexGeometry::getDecimationFreq() const
{
  return m_DecimationFreq;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void DownsampleVertexGeometry::setDecimationFraction(const float& value)
{
  m_DecimationFraction = value;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
float DownsampleVertexGeometry::getDecimationFraction() const
{
  return m_DecimationFraction;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void DownsampleVertexGeometry::setGridResolution(const FloatVec3Type& value)
{
  m_GridResolution = value;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
FloatVec3Type DownsampleVertexGeometry::getGridResolution() const
{
  return m_GridResolution;
}
