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

#include "ExportCLIFile.h"

#include <cassert>
#include <map>

#include <QtCore/QDir>
#include <QtCore/QTextStream>

#include "SIMPLib/Common/Constants.h"

#include "SIMPLib/FilterParameters/DataArraySelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/DataContainerSelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/DoubleFilterParameter.h"
#include "SIMPLib/FilterParameters/IntFilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedBooleanFilterParameter.h"
#include "SIMPLib/FilterParameters/OutputPathFilterParameter.h"
#include "SIMPLib/FilterParameters/SeparatorFilterParameter.h"
#include "SIMPLib/FilterParameters/StringFilterParameter.h"
#include "SIMPLib/Geometry/EdgeGeom.h"
#include "SIMPLib/Math/SIMPLibMath.h"
#include "SIMPLib/DataContainers/DataContainerArray.h"

#include "DREAM3DReview/DREAM3DReviewConstants.h"
#include "DREAM3DReview/DREAM3DReviewVersion.h"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
ExportCLIFile::ExportCLIFile()
: m_EdgeGeometry("")
, m_LayerIdsArrayPath("", "", "")
, m_GroupIdsArrayPath("", "", "")
, m_OutputDirectory("")
, m_OutputFilePrefix("")
, m_UnitsScaleFactor(1.0)
, m_Precision(5)
{
  initialize();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
ExportCLIFile::~ExportCLIFile() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ExportCLIFile::initialize()
{
  clearErrorCode();
  clearWarningCode();
  setCancel(false);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ExportCLIFile::setupFilterParameters()
{
  FilterParameterVectorType parameters;
  parameters.push_back(SIMPL_NEW_DOUBLE_FP("Units Scale Factor", UnitsScaleFactor, FilterParameter::Parameter, ExportCLIFile));
  parameters.push_back(SIMPL_NEW_INTEGER_FP("Precision (places after decimal)", Precision, FilterParameter::Parameter, ExportCLIFile));
  parameters.push_back(SIMPL_NEW_OUTPUT_PATH_FP("Output File Directory", OutputDirectory, FilterParameter::Parameter, ExportCLIFile));
  parameters.push_back(SIMPL_NEW_STRING_FP("Output File Prefix", OutputFilePrefix, FilterParameter::Parameter, ExportCLIFile));
  QStringList linkedProps = {"GroupIdsArrayPath"};
  parameters.push_back(SIMPL_NEW_LINKED_BOOL_FP("Split CLI Files by Group", SplitByGroup, FilterParameter::Parameter, ExportCLIFile, linkedProps));
  DataContainerSelectionFilterParameter::RequirementType dcsReq;
  IGeometry::Types geomTypes = {IGeometry::Type::Edge};
  dcsReq.dcGeometryTypes = geomTypes;
  parameters.push_back(SIMPL_NEW_DC_SELECTION_FP("Edge Geometry", EdgeGeometry, FilterParameter::RequiredArray, ExportCLIFile, dcsReq));
  parameters.push_back(SeparatorFilterParameter::New("Edge Data", FilterParameter::RequiredArray));
  DataArraySelectionFilterParameter::RequirementType dasReq = DataArraySelectionFilterParameter::CreateRequirement(SIMPL::TypeNames::Int32, 1, AttributeMatrix::Type::Edge, IGeometry::Type::Edge);
  parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Layer Ids", LayerIdsArrayPath, FilterParameter::RequiredArray, ExportCLIFile, dasReq));
  parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Group Ids", GroupIdsArrayPath, FilterParameter::RequiredArray, ExportCLIFile, dasReq));
  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ExportCLIFile::dataCheck()
{
  clearErrorCode();
  clearWarningCode();

  if(getUnitsScaleFactor() <= 0)
  {
    QString ss = QObject::tr("Units scale factor must be greater than 0");
    setErrorCondition(-1, ss);
  }

  if(getPrecision() <= 0)
  {
    QString ss = QObject::tr("Precision must be greater than 0");
    setErrorCondition(-1, ss);
  }

  if(getOutputDirectory().isEmpty())
  {
    QString ss = QObject::tr("The output directory must be set");
    setErrorCondition(-391, ss);
  }

  QDir dir(getOutputDirectory());

  if(getOutputDirectory().isEmpty())
  {
    QString ss = QObject::tr("The output directory must be set");
    setErrorCondition(-391, ss);
  }
  else if(!dir.exists())
  {
    QString ss = QObject::tr("The output directory does not exist; DREAM.3D will attempt to create this path during execution");
    setWarningCondition(-1, ss);
  }

  if(getOutputFilePrefix().isEmpty())
  {
    QString ss = QObject::tr("The output file prefix must be set");
    setErrorCondition(-392, ss);
  }

  EdgeGeom::Pointer edge = getDataContainerArray()->getPrereqGeometryFromDataContainer<EdgeGeom>(this, getEdgeGeometry());

  QVector<IDataArray::Pointer> dataArrays;

  if(getErrorCode() >= 0)
  {
    dataArrays.push_back(edge->getEdges());
  }

  std::vector<size_t> cDims(1, 1);
  m_LayerIdsPtr = getDataContainerArray()->getPrereqArrayFromPath<Int32ArrayType>(this, getLayerIdsArrayPath(), cDims);
  if(m_LayerIdsPtr.lock())
  {
    m_LayerIds = m_LayerIdsPtr.lock()->getPointer(0);
  }

  if(getErrorCode() >= 0)
  {
    dataArrays.push_back(m_LayerIdsPtr.lock());
  }

  if(getSplitByGroup())
  {
    m_GroupIdsPtr = getDataContainerArray()->getPrereqArrayFromPath<Int32ArrayType>(this, getGroupIdsArrayPath(), cDims);
    if(m_GroupIdsPtr.lock())
    {
      m_GroupIds = m_GroupIdsPtr.lock()->getPointer(0);
    }

    if(getErrorCode() >= 0)
    {
      dataArrays.push_back(m_GroupIdsPtr.lock());
    }
  }

  getDataContainerArray()->validateNumberOfTuples(this, dataArrays);
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ExportCLIFile::execute()
{
  initialize();
  dataCheck();
  if(getErrorCode() < 0)
  {
    return;
  }

  EdgeGeom::Pointer edge = getDataContainerArray()->getDataContainer(m_EdgeGeometry)->getGeometryAs<EdgeGeom>();
  MeshIndexType* edges = edge->getEdgePointer(0);
  float* vertices = edge->getVertexPointer(0);
  MeshIndexType numEdges = edge->getNumberOfEdges();
  // int64_t numVerts = edge->getNumberOfVertices();
  int32_t numGroups = m_SplitByGroup ? 0 : 1;
  int32_t numLayers = 0;

  for(MeshIndexType i = 0; i < numEdges; i++)
  {
    if(m_SplitByGroup)
    {
      if(m_GroupIds[i] > numGroups)
      {
        numGroups = m_GroupIds[i];
      }
    }
    if(m_LayerIds[i] > numLayers)
    {
      numLayers = m_LayerIds[i];
    }
  }

  if(numLayers == 0)
  {
    QString ss = QObject::tr("No valid layers found in Layer Ids");
    setErrorCondition(-1, ss);
    return;
  }

  QDir dir;
  if(!dir.mkpath(m_OutputDirectory))
  {
    QString ss = QObject::tr("Error creating output path '%1'").arg(m_OutputDirectory);
    setErrorCondition(-11110, ss);
    return;
  }

  numGroups++;
  numLayers++;

  QVector<std::shared_ptr<QFile>> files;
  for(auto i = 0; i < numGroups; i++)
  {
    QString fname = m_OutputDirectory + "/" + m_OutputFilePrefix + "Group" + QString::number(i + 1) + ".cli";
    std::shared_ptr<QFile> filePtr = std::make_shared<QFile>(fname);
    filePtr->open(QIODevice::WriteOnly | QIODevice::Text);
    files.push_back(filePtr);
  }

  for(auto&& f : files)
  {
    QString header = "$$HEADERSTART\n"
                     "$$ASCII\n"
                     "$$UNITS/" +
                     QString::number(m_UnitsScaleFactor, 'f', m_Precision) + "\n"
                                                                             "$$HEADEREND\n"
                                                                             "$$GEOMETRYSTART\n"
                                                                             "\n"
                                                                             "$$LAYER/0.00000\n"
                                                                             "\n";

    QTextStream out(f.get());
    out << header;
  }

  using GroupLayerTable = std::vector<std::vector<std::vector<int64_t>>>;
  using LayerTable = std::vector<std::vector<int64_t>>;
  GroupLayerTable table(numGroups, LayerTable(numLayers));
  Int32ArrayType::Pointer groups = Int32ArrayType::CreateArray(numEdges, "_INTERNAL_USE_ONLY_GroupIds", false);
  int32_t* groupIds = nullptr;
  if(m_SplitByGroup)
  {
    groupIds = m_GroupIds;
  }
  else
  {
    groups->allocate();
    groups->initializeWithValue(1);
    groupIds = groups->getPointer(0);
  }
  for(MeshIndexType i = 0; i < numEdges; i++)
  {
    table[groupIds[i]][m_LayerIds[i]].push_back(i);
  }

  for(auto i = 0; i < numGroups; i++)
  {
    QFile* file = files[i].get();
    QTextStream out(file);
    for(auto j = 0; j < numLayers; j++)
    {
      if(!table[i][j].empty())
      {
        double layerHeight = static_cast<double>(vertices[3 * edges[2 * table[i][j].front() + 0] + 2]);
        QString layerStr = "$$LAYER/" + QString::number(layerHeight / m_UnitsScaleFactor, 'f', m_Precision) + "\n";
        out << layerStr;
        QString hatchStr = "$$HATCHES/1," + QString::number(table[i][j].size()) + ",";
        for(auto&& hatch : table[i][j])
        {
          double zStart = static_cast<double>(vertices[3 * edges[2 * hatch + 0] + 2] / m_UnitsScaleFactor);
          double zEnd = static_cast<double>(vertices[3 * edges[2 * hatch + 1] + 2] / m_UnitsScaleFactor);
          if(!SIMPLibMath::closeEnough(zStart, layerHeight) || !SIMPLibMath::closeEnough(zEnd, layerHeight))
          {
            QString ss = QObject::tr("Found Edge (%1) that spans multipe layers").arg(hatch);
            setErrorCondition(-1, ss);
            return;
          }

          double xStart = static_cast<double>(vertices[3 * edges[2 * hatch + 0] + 0] / m_UnitsScaleFactor);
          double yStart = static_cast<double>(vertices[3 * edges[2 * hatch + 0] + 1] / m_UnitsScaleFactor);
          double xEnd = static_cast<double>(vertices[3 * edges[2 * hatch + 1] + 0] / m_UnitsScaleFactor);
          double yEnd = static_cast<double>(vertices[3 * edges[2 * hatch + 1] + 1] / m_UnitsScaleFactor);
          hatchStr += (QString::number(xStart, 'f', m_Precision) + "," + QString::number(yStart, 'f', m_Precision) + "," + QString::number(xEnd, 'f', m_Precision) + "," +
                       QString::number(yEnd, 'f', m_Precision) + ",");
        }
        hatchStr.chop(1);
        hatchStr += "\n\n";
        out << hatchStr;
      }
    }
    QString footer = "$$GEOMETRYEND\n";
    out << footer;
  }

  notifyStatusMessage("Complete");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AbstractFilter::Pointer ExportCLIFile::newFilterInstance(bool copyFilterParameters) const
{
  ExportCLIFile::Pointer filter = ExportCLIFile::New();
  if(copyFilterParameters)
  {
    copyFilterParameterInstanceVariables(filter.get());
  }
  return filter;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ExportCLIFile::getCompiledLibraryName() const
{
  return DREAM3DReviewConstants::DREAM3DReviewBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ExportCLIFile::getBrandingString() const
{
  return "DREAM3DReview";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ExportCLIFile::getFilterVersion() const
{
  QString version;
  QTextStream vStream(&version);
  vStream << DREAM3DReview::Version::Major() << "." << DREAM3DReview::Version::Minor() << "." << DREAM3DReview::Version::Patch();
  return version;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ExportCLIFile::getGroupName() const
{
  return SIMPL::FilterGroups::IOFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ExportCLIFile::getSubGroupName() const
{
  return SIMPL::FilterSubGroups::OutputFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ExportCLIFile::getHumanLabel() const
{
  return "Export CLI File(s)";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QUuid ExportCLIFile::getUuid() const
{
  return QUuid("{351c38c9-330d-599c-8632-19d74868be86}");
}

// -----------------------------------------------------------------------------
ExportCLIFile::Pointer ExportCLIFile::NullPointer()
{
  return Pointer(static_cast<Self*>(nullptr));
}

// -----------------------------------------------------------------------------
std::shared_ptr<ExportCLIFile> ExportCLIFile::New()
{
  struct make_shared_enabler : public ExportCLIFile
  {
  };
  std::shared_ptr<make_shared_enabler> val = std::make_shared<make_shared_enabler>();
  val->setupFilterParameters();
  return val;
}

// -----------------------------------------------------------------------------
QString ExportCLIFile::getNameOfClass() const
{
  return QString("ExportCLIFile");
}

// -----------------------------------------------------------------------------
QString ExportCLIFile::ClassName()
{
  return QString("ExportCLIFile");
}

// -----------------------------------------------------------------------------
void ExportCLIFile::setEdgeGeometry(const DataArrayPath& value)
{
  m_EdgeGeometry = value;
}

// -----------------------------------------------------------------------------
DataArrayPath ExportCLIFile::getEdgeGeometry() const
{
  return m_EdgeGeometry;
}

// -----------------------------------------------------------------------------
void ExportCLIFile::setLayerIdsArrayPath(const DataArrayPath& value)
{
  m_LayerIdsArrayPath = value;
}

// -----------------------------------------------------------------------------
DataArrayPath ExportCLIFile::getLayerIdsArrayPath() const
{
  return m_LayerIdsArrayPath;
}

// -----------------------------------------------------------------------------
void ExportCLIFile::setSplitByGroup(bool value)
{
  m_SplitByGroup = value;
}

// -----------------------------------------------------------------------------
bool ExportCLIFile::getSplitByGroup() const
{
  return m_SplitByGroup;
}

// -----------------------------------------------------------------------------
void ExportCLIFile::setGroupIdsArrayPath(const DataArrayPath& value)
{
  m_GroupIdsArrayPath = value;
}

// -----------------------------------------------------------------------------
DataArrayPath ExportCLIFile::getGroupIdsArrayPath() const
{
  return m_GroupIdsArrayPath;
}

// -----------------------------------------------------------------------------
void ExportCLIFile::setOutputDirectory(const QString& value)
{
  m_OutputDirectory = value;
}

// -----------------------------------------------------------------------------
QString ExportCLIFile::getOutputDirectory() const
{
  return m_OutputDirectory;
}

// -----------------------------------------------------------------------------
void ExportCLIFile::setOutputFilePrefix(const QString& value)
{
  m_OutputFilePrefix = value;
}

// -----------------------------------------------------------------------------
QString ExportCLIFile::getOutputFilePrefix() const
{
  return m_OutputFilePrefix;
}

// -----------------------------------------------------------------------------
void ExportCLIFile::setUnitsScaleFactor(double value)
{
  m_UnitsScaleFactor = value;
}

// -----------------------------------------------------------------------------
double ExportCLIFile::getUnitsScaleFactor() const
{
  return m_UnitsScaleFactor;
}

// -----------------------------------------------------------------------------
void ExportCLIFile::setPrecision(int value)
{
  m_Precision = value;
}

// -----------------------------------------------------------------------------
int ExportCLIFile::getPrecision() const
{
  return m_Precision;
}
