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
#include "ImportCLIFile.h"

#include <QtCore/QFileInfo>

#include "SIMPLib/Common/Constants.h"

#include "SIMPLib/FilterParameters/InputFileFilterParameter.h"
#include "SIMPLib/FilterParameters/SeparatorFilterParameter.h"
#include "SIMPLib/FilterParameters/StringFilterParameter.h"
#include "SIMPLib/Geometry/EdgeGeom.h"

#include "DREAM3DReview/DREAM3DReviewConstants.h"
#include "DREAM3DReview/DREAM3DReviewVersion.h"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
ImportCLIFile::ImportCLIFile()
: m_CLIFile("")
, m_EdgeDataContainerName(SIMPL::Defaults::EdgeDataContainerName)
, m_VertexAttributeMatrixName(SIMPL::Defaults::VertexAttributeMatrixName)
, m_EdgeAttributeMatrixName(SIMPL::Defaults::EdgeAttributeMatrixName)
, m_LayerIdsArrayName("LayerIds")
, m_FeatureIdsArrayName(SIMPL::CellData::FeatureIds)
{
  initialize();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
ImportCLIFile::~ImportCLIFile() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ImportCLIFile::initialize()
{
  clearErrorCode();
  clearWarningCode();
  setCancel(false);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ImportCLIFile::setupFilterParameters()
{
  FilterParameterVectorType parameters;
  parameters.push_back(SIMPL_NEW_INPUT_FILE_FP("CLI File", CLIFile, FilterParameter::Parameter, ImportCLIFile, "*.cli"));
  parameters.push_back(SIMPL_NEW_STRING_FP("Edge Data Container", EdgeDataContainerName, FilterParameter::CreatedArray, ImportCLIFile));
  parameters.push_back(SeparatorFilterParameter::New("Vertex Data", FilterParameter::CreatedArray));
  parameters.push_back(SIMPL_NEW_STRING_FP("Vertex Attribute Matrix", VertexAttributeMatrixName, FilterParameter::CreatedArray, ImportCLIFile));
  parameters.push_back(SeparatorFilterParameter::New("Edge Data", FilterParameter::CreatedArray));
  parameters.push_back(SIMPL_NEW_STRING_FP("Edge Attribute Matrix", EdgeAttributeMatrixName, FilterParameter::CreatedArray, ImportCLIFile));
  parameters.push_back(SIMPL_NEW_STRING_FP("Layer Ids", LayerIdsArrayName, FilterParameter::CreatedArray, ImportCLIFile));
  parameters.push_back(SIMPL_NEW_STRING_FP("Feature Ids", FeatureIdsArrayName, FilterParameter::CreatedArray, ImportCLIFile));
  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ImportCLIFile::dataCheck()
{
  clearErrorCode();
  clearWarningCode();

  QFileInfo fi(getCLIFile());

  if(getCLIFile().isEmpty())
  {
    QString ss = QObject::tr("The input CLI file must be set");
    setErrorCondition(-387, ss);
  }
  else if(!fi.exists() || !fi.isFile())
  {
    QString ss = QObject::tr("The input CLI file does not exist");
    setErrorCondition(-388, ss);
  }

  if(getErrorCode() < 0)
  {
    return;
  }

  DataContainer::Pointer dc = getDataContainerArray()->createNonPrereqDataContainer(this, getEdgeDataContainerName());
  SharedVertexList::Pointer vertices = EdgeGeom::CreateSharedVertexList(0, !getInPreflight());
  EdgeGeom::Pointer edge = EdgeGeom::CreateGeometry(0, vertices, SIMPL::Geometry::EdgeGeometry, !getInPreflight());

  if(getErrorCode() < 0)
  {
    return;
  }

  dc->setGeometry(edge);

  QVector<size_t> tDims(1, 0);
  dc->createNonPrereqAttributeMatrix(this, getVertexAttributeMatrixName(), tDims, AttributeMatrix::Type::Vertex);
  dc->createNonPrereqAttributeMatrix(this, getEdgeAttributeMatrixName(), tDims, AttributeMatrix::Type::Edge);

  if(getErrorCode() < 0)
  {
    return;
  }

  DataArrayPath path(getEdgeDataContainerName(), getEdgeAttributeMatrixName(), getLayerIdsArrayName());
  QVector<size_t> cDims(1, 1);

  getDataContainerArray()->createNonPrereqArrayFromPath<Int32ArrayType, AbstractFilter, int32_t>(this, path, 0, cDims);

  path.setDataArrayName(getFeatureIdsArrayName());

  getDataContainerArray()->createNonPrereqArrayFromPath<Int32ArrayType, AbstractFilter, int32_t>(this, path, 0, cDims);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ImportCLIFile::preflight()
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
void ImportCLIFile::execute()
{
  initialize();
  dataCheck();
  if(getErrorCode() < 0)
  {
    return;
  }

  QFile stream;
  stream.setFileName(getCLIFile());
  if(!stream.open(QIODevice::ReadOnly | QIODevice::Text))
  {
    QString ss = QObject::tr("Input CLI file could not be opened: %1").arg(getCLIFile());
    setErrorCondition(-100, ss);
    return;
  }

  QByteArray buf;
  QList<QByteArray> tokens;
  bool ok = false;
  float units = 1.0f;
  std::vector<float> tmpVertices;
  std::vector<int64_t> tmpEdges;
  std::vector<int32_t> tmpLayerIds;
  std::vector<int32_t> tmpFeatureIds;
  float layerHeight = 0.0f;
  int32_t lines = 1;
  int32_t layer = 0;
  int32_t features = 0;
  int64_t vertexCounter = 0;

  while(!stream.atEnd())
  {
    buf = stream.readLine();
    buf = buf.trimmed();
    buf = buf.simplified();
    if(buf.startsWith("$$UNITS"))
    {
      tokens = buf.split('/');
      if(tokens.size() != 2)
      {
        QString ss = QObject::tr("Unable to parse units from CLI file line %1: %2").arg(lines).arg(QString::fromStdString(buf.toStdString()));
        setErrorCondition(-1, ss);
        return;
      }
      units = tokens.at(1).toFloat(&ok);
      if(!ok)
      {
        QString ss = QObject::tr("Unable to parse units from CLI file line %1: %2").arg(lines).arg(QString::fromStdString(buf.toStdString()));
        setErrorCondition(-1, ss);
        return;
      }
    }
    else if(buf.startsWith("$$LAYER"))
    {
      layer++;
      tokens = buf.split('/');
      if(tokens.size() != 2)
      {
        QString ss = QObject::tr("Unable to parse layer height from CLI file line %1: %2").arg(lines).arg(QString::fromStdString(buf.toStdString()));
        setErrorCondition(-1, ss);
        return;
      }
      layerHeight = tokens.at(1).toFloat(&ok);
      if(!ok)
      {
        QString ss = QObject::tr("Unable to parse layer height from CLI file line %1: %2").arg(lines).arg(QString::fromStdString(buf.toStdString()));
        setErrorCondition(-1, ss);
        return;
      }
    }
    else if(buf.startsWith("$$POLYLINE") || buf.startsWith("$POLYLINE"))
    {
      tokens = buf.split('/');
      if(tokens.size() != 2)
      {
        QString ss = QObject::tr("Unable to parse polyline from CLI file line %1: %2").arg(lines).arg(QString::fromStdString(buf.toStdString()));
        setErrorCondition(-1, ss);
        return;
      }
      tokens = tokens.at(1).split(',');
      if((tokens.size() - 3) > 0)
      {
        features++;
        if((tokens.size() - 3) % 2 != 0)
        {
          QString ss = QObject::tr("Polyline at line %1 does not contain an even number of elements: %2").arg(lines).arg(QString::fromStdString(buf.toStdString()));
          setErrorCondition(-1, ss);
          return;
        }
        for(auto i = 3; i < (tokens.size() - 2); i++)
        {
          float val = tokens.at(i).toFloat(&ok);
          if(!ok)
          {
            QString ss = QObject::tr("Unable to parse polyline coordinate from CLI file line %1: %2").arg(lines).arg(QString::fromStdString(buf.toStdString()));
            setErrorCondition(-1, ss);
            return;
          }
          tmpVertices.push_back(val);
          if(i % 2 == 0)
          {
            tmpVertices.push_back(layerHeight);
          }
        }
        int64_t numVerts = ((tokens.size() - 3) / 2) - 1;
        int32_t edgeCounter = 0;
        for(auto i = vertexCounter; i < (vertexCounter + numVerts); i++)
        {
          tmpEdges.push_back(i);
          tmpEdges.push_back(i + 1);
          edgeCounter++;
        }
        tmpEdges.pop_back();
        tmpEdges.push_back(vertexCounter);
        vertexCounter += numVerts;
        for(auto i = 0; i < edgeCounter; i++)
        {
          tmpLayerIds.push_back(layer);
          tmpFeatureIds.push_back(features);
        }
      }
    }
    else if(buf.startsWith("$$HATCHES") || buf.startsWith("$HATCHES"))
    {
      tokens = buf.split('/');
      if(tokens.size() != 2)
      {
        QString ss = QObject::tr("Unable to parse hatch from CLI file line %1: %2").arg(lines).arg(QString::fromStdString(buf.toStdString()));
        setErrorCondition(-1, ss);
        return;
      }
      tokens = tokens.at(1).split(',');
      if((tokens.size() - 2) > 0)
      {
        features++;
        if((tokens.size() - 2) % 4 != 0)
        {
          QString ss = QObject::tr("Hatch at line %1 does not contain an even number of elements: %2").arg(lines).arg(QString::fromStdString(buf.toStdString()));
          setErrorCondition(-1, ss);
          return;
        }
        for(auto i = 2; i < tokens.size(); i += 4)
        {
          for(auto j = 0; j < 4; j++)
          {
            float val = tokens.at(i + j).toFloat(&ok);
            if(!ok)
            {
              QString ss = QObject::tr("Unable to parse hatch coordinate from CLI file line %1: %2").arg(lines).arg(QString::fromStdString(buf.toStdString()));
              setErrorCondition(-1, ss);
              return;
            }
            tmpVertices.push_back(val);
            if(j % 2 != 0)
            {
              tmpVertices.push_back(layerHeight);
            }
          }
        }
        int64_t numVerts = ((tokens.size() - 2) / 2);
        int32_t edgeCounter = 0;
        for(auto i = vertexCounter; i < (vertexCounter + numVerts); i += 2)
        {
          tmpEdges.push_back(i);
          tmpEdges.push_back(i + 1);
          edgeCounter++;
        }
        vertexCounter += numVerts;
        for(auto i = 0; i < edgeCounter; i++)
        {
          tmpLayerIds.push_back(layer);
          tmpFeatureIds.push_back(features);
        }
      }
    }
    lines++;
  }

  tmpVertices.shrink_to_fit();
  tmpEdges.shrink_to_fit();
  tmpLayerIds.shrink_to_fit();
  tmpFeatureIds.shrink_to_fit();

  std::transform(std::begin(tmpVertices), std::end(tmpVertices), std::begin(tmpVertices), [&](float val) { return val * units; });

  EdgeGeom::Pointer edge = getDataContainerArray()->getDataContainer(m_EdgeDataContainerName)->getGeometryAs<EdgeGeom>();
  edge->resizeVertexList(tmpVertices.size() / 3);
  edge->resizeEdgeList(tmpEdges.size() / 2);
  float* verts = edge->getVertexPointer(0);
  int64_t* edges = edge->getEdgePointer(0);
  std::memcpy(verts, tmpVertices.data(), 3 * edge->getNumberOfVertices() * sizeof(float));
  std::memcpy(edges, tmpEdges.data(), 2 * edge->getNumberOfEdges() * sizeof(int64_t));

  AttributeMatrix::Pointer edgeAttrMat = getDataContainerArray()->getDataContainer(m_EdgeDataContainerName)->getAttributeMatrix(m_EdgeAttributeMatrixName);
  AttributeMatrix::Pointer vertAttrMat = getDataContainerArray()->getDataContainer(m_EdgeDataContainerName)->getAttributeMatrix(m_VertexAttributeMatrixName);
  QVector<size_t> tDims(1, edge->getNumberOfEdges());
  edgeAttrMat->resizeAttributeArrays(tDims);
  tDims[0] = edge->getNumberOfVertices();
  vertAttrMat->resizeAttributeArrays(tDims);
  int32_t* layerIds = edgeAttrMat->getAttributeArrayAs<Int32ArrayType>(m_LayerIdsArrayName)->getPointer(0);
  int32_t* featureIds = edgeAttrMat->getAttributeArrayAs<Int32ArrayType>(m_FeatureIdsArrayName)->getPointer(0);
  std::memcpy(layerIds, tmpLayerIds.data(), edge->getNumberOfEdges() * sizeof(int32_t));
  std::memcpy(featureIds, tmpFeatureIds.data(), edge->getNumberOfEdges() * sizeof(int32_t));

  notifyStatusMessage("Complete");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AbstractFilter::Pointer ImportCLIFile::newFilterInstance(bool copyFilterParameters) const
{
  ImportCLIFile::Pointer filter = ImportCLIFile::New();
  if(copyFilterParameters)
  {
    copyFilterParameterInstanceVariables(filter.get());
  }
  return filter;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString ImportCLIFile::getCompiledLibraryName() const
{
  return DREAM3DReviewConstants::DREAM3DReviewBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString ImportCLIFile::getBrandingString() const
{
  return "DREAM3DReview";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString ImportCLIFile::getFilterVersion() const
{
  QString version;
  QTextStream vStream(&version);
  vStream << DREAM3DReview::Version::Major() << "." << DREAM3DReview::Version::Minor() << "." << DREAM3DReview::Version::Patch();
  return version;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString ImportCLIFile::getGroupName() const
{
  return SIMPL::FilterGroups::IOFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString ImportCLIFile::getSubGroupName() const
{
  return SIMPL::FilterSubGroups::InputFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString ImportCLIFile::getHumanLabel() const
{
  return "Import CLI File";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QUuid ImportCLIFile::getUuid()
{
  return QUuid("{76ec2a58-a0cd-5548-b248-5a5eb08a1581}");
}