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

#include "ImportVolumeGraphicsFile.h"

#include <tuple>

#include <QtCore/QDir>
#include <QtCore/QFileInfo>
#include <QtCore/QTextStream>

#include "SIMPLib/Common/Constants.h"
#include "SIMPLib/Common/ScopedFileMonitor.hpp"
#include "SIMPLib/DataContainers/DataContainer.h"
#include "SIMPLib/DataContainers/DataContainerArray.h"
#include "SIMPLib/FilterParameters/AbstractFilterParametersReader.h"
#include "SIMPLib/FilterParameters/InputFileFilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedPathCreationFilterParameter.h"
#include "SIMPLib/FilterParameters/SeparatorFilterParameter.h"
#include "SIMPLib/FilterParameters/StringFilterParameter.h"

#include "DREAM3DReview/DREAM3DReviewConstants.h"
#include "DREAM3DReview/DREAM3DReviewVersion.h"

#define CREATE_BLOCK_CONST(VAR) static inline const QString k_##VAR##Block("[" #VAR "]");

#define CREATE_ELEMENT_CONST(VAR) static inline const QString k_##VAR(#VAR);

namespace ImportVolumeGraphicsFileConstants
{
static inline constexpr int32_t k_EmptyFileInput = -91500;
static inline constexpr int32_t k_FileDoesNotExist = -91501;
static inline constexpr int32_t k_VgiOpenError = -91502;
static inline constexpr int32_t k_VgiParseError = -91503;
static inline constexpr int32_t k_FileIsDirectory = -91507;

static inline constexpr int32_t k_VolBinaryAllocateMismatch = -91504;
static inline constexpr int32_t k_VolOpenError = -91505;
static inline constexpr int32_t k_VolReadError = -91506;

static inline const QString k_Millimeter("mm");

CREATE_BLOCK_CONST(representation)

// File Block
CREATE_BLOCK_CONST(file1)

CREATE_ELEMENT_CONST(RegionOfInterestStart)
CREATE_ELEMENT_CONST(RegionOfInterestEnd)
CREATE_ELEMENT_CONST(FileFormat)
CREATE_ELEMENT_CONST(Size)
CREATE_ELEMENT_CONST(Name)
CREATE_ELEMENT_CONST(Datatype)
CREATE_ELEMENT_CONST(datarange)
CREATE_ELEMENT_CONST(BitsPerElement)

// Geometry Block
CREATE_BLOCK_CONST(geometry)

CREATE_ELEMENT_CONST(status)
CREATE_ELEMENT_CONST(relativeposition)
CREATE_ELEMENT_CONST(position)
CREATE_ELEMENT_CONST(resolution)
CREATE_ELEMENT_CONST(scale)
CREATE_ELEMENT_CONST(center)
CREATE_ELEMENT_CONST(rotate)
CREATE_ELEMENT_CONST(unit)

} // namespace ImportVolumeGraphicsFileConstants

using namespace ImportVolumeGraphicsFileConstants;
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
ImportVolumeGraphicsFile::ImportVolumeGraphicsFile() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
ImportVolumeGraphicsFile::~ImportVolumeGraphicsFile() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ImportVolumeGraphicsFile::setupFilterParameters()
{
  FilterParameterVectorType parameters;
  parameters.push_back(SIMPL_NEW_INPUT_FILE_FP("VolumeGraphics .vgi File", VGHeaderFile, FilterParameter::Category::Parameter, ImportVolumeGraphicsFile, "*.vgi", "Volume Graphics Header"));
  // parameters.push_back(SIMPL_NEW_INPUT_FILE_FP("VolumeGraphics .vol File", VGDataFile, FilterParameter::Category::Parameter, ImportVolumeGraphicsFile, "*.vol", "Volume Graphics Data"));
  parameters.push_back(SIMPL_NEW_STRING_FP("Data Container", DataContainerName, FilterParameter::Category::CreatedArray, ImportVolumeGraphicsFile));
  parameters.push_back(SeparatorFilterParameter::Create("Cell Data", FilterParameter::Category::CreatedArray));
  parameters.push_back(SIMPL_NEW_AM_WITH_LINKED_DC_FP("Cell Attribute Matrix", CellAttributeMatrixName, DataContainerName, FilterParameter::Category::CreatedArray, ImportVolumeGraphicsFile));
  parameters.push_back(SIMPL_NEW_DA_WITH_LINKED_AM_FP("Density", DensityArrayName, DataContainerName, CellAttributeMatrixName, FilterParameter::Category::CreatedArray, ImportVolumeGraphicsFile));
  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ImportVolumeGraphicsFile::initialize()
{
  if(m_InHeaderStream.isOpen())
  {
    m_InHeaderStream.close();
  }
  if(m_InStream.isOpen())
  {
    m_InStream.close();
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ImportVolumeGraphicsFile::dataCheck()
{
  clearErrorCode();
  clearWarningCode();
  initialize();
  int32_t err = 0;

  QFileInfo fiHdr(getVGHeaderFile());

  if(getVGHeaderFile().isEmpty())
  {
    QString ss = QObject::tr("The Volume Graphics header file (.vgi) must be set");
    setErrorCondition(k_EmptyFileInput, ss);
  }
  else if(fiHdr.isDir())
  {
    QString ss = QObject::tr("The file path is a directory. Please select a file with the .vgi extension. '%1'").arg(getVGHeaderFile());
    setErrorCondition(k_FileIsDirectory, ss);
  }
  else if(!fiHdr.exists())
  {
    QString ss = QObject::tr("The Volume Graphics header file (.vgi) file does not exist. '%1'").arg(getVGHeaderFile());
    setErrorCondition(k_FileDoesNotExist, ss);
  }

  if(getErrorCode() < 0)
  {
    return;
  }

  if(m_InHeaderStream.isOpen())
  {
    m_InHeaderStream.close();
  }

  m_InHeaderStream.setFileName(getVGHeaderFile());

  if(!m_InHeaderStream.open(QIODevice::ReadOnly | QIODevice::Text))
  {
    QString ss = QObject::tr("Error opening input header file: %1").arg(getVGHeaderFile());
    setErrorCondition(k_VgiOpenError, ss);
    return;
  }

  ImageGeom::Pointer image = ImageGeom::CreateGeometry(SIMPL::Geometry::ImageGeometry);

  readHeaderMetaData(image);
  if(getErrorCode() < 0)
  {
    return;
  }
  QFileInfo fi(getVGDataFile());

  if(getVGDataFile().isEmpty())
  {
    QString ss = QObject::tr("The Volume Graphics voxel file (.vol) must be set");
    setErrorCondition(k_EmptyFileInput, ss);
  }
  else if(fi.isDir())
  {
    QString ss = QObject::tr("The file path is a directory. Please select a file with the .vgi extension. '%1'").arg(getVGDataFile());
    setErrorCondition(k_FileIsDirectory, ss);
  }
  else if(!fi.exists())
  {
    QString ss = QObject::tr("The Volume Graphics voxel file (.vol) file does not exist. '%1'").arg(getVGDataFile());
    setErrorCondition(k_FileDoesNotExist, ss);
  }

  if(err < 0)
  {
    return;
  }

  DataContainer::Pointer m = getDataContainerArray()->createNonPrereqDataContainer(this, getDataContainerName());

  if(getErrorCode() < 0)
  {
    return;
  }

  m->setGeometry(image);

  SizeVec3Type dims = image->getDimensions();

  std::vector<size_t> tDims = {dims[0], dims[1], dims[2]};
  std::vector<size_t> cDims(1, 1);

  m->createNonPrereqAttributeMatrix(this, getCellAttributeMatrixName(), tDims, AttributeMatrix::Type::Cell);

  DataArrayPath path(getDataContainerName(), getCellAttributeMatrixName(), getDensityArrayName());

  m_DensityPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<float>>(this, path, 0, cDims);
  if(nullptr != m_DensityPtr.lock())
  {
    m_Density = m_DensityPtr.lock()->getPointer(0);
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ImportVolumeGraphicsFile::execute()
{
  dataCheck();
  if(getErrorCode() < 0)
  {
    return;
  }

  int32_t err = readVolFile();
  if(err < 0)
  {
    QString ss = QObject::tr("Error reading binary input file");
    setErrorCondition(k_VolReadError, ss);
  }

  notifyStatusMessage("Complete");
}

// -----------------------------------------------------------------------------
int32_t ImportVolumeGraphicsFile::readVolFile()
{
  ImageGeom::Pointer image = getDataContainerArray()->getDataContainer(getDataContainerName())->getGeometryAs<ImageGeom>();
  int32_t error = 0;

  QFileInfo fi(getVGDataFile());

  size_t filesize = static_cast<size_t>(fi.size());
  FloatArrayType& density = *(m_DensityPtr.lock().get());
  size_t allocatedBytes = density.getSize() * sizeof(float);

  if(filesize < allocatedBytes)
  {
    QString ss = QObject::tr("Binary file size is smaller than the number of allocated bytes");
    setErrorCondition(k_VolBinaryAllocateMismatch, ss);
    return getErrorCode();
  }

  FILE* f = fopen(getVGDataFile().toLatin1().data(), "rb");
  if(nullptr == f)
  {
    QString ss = QObject::tr("Error opening binary input file: %1").arg(getVGDataFile());
    setErrorCondition(k_VolOpenError, ss);
    return getErrorCode();
  }
  ScopedFileMonitor monitor(f);

  QString ss = QObject::tr("Reading Data from .vol File.....");
  notifyStatusMessage(ss);
  if(fread(density.getTuplePointer(0), 1, filesize, f) != filesize)
  {
    ss = QObject::tr("Error Reading .vol file. Not enough bytes read....");
    setErrorCondition(k_VolReadError, ss);
    return getErrorCode();
  }

  return error;
}

// -----------------------------------------------------------------------------
SizeVec3Type ParseSizeVec3(ImportVolumeGraphicsFile* filter, const QList<QByteArray>& tokens)
{
  SizeVec3Type vec = {0, 0, 0};
  bool ok = false;

  vec[0] = tokens[2].toULongLong(&ok);
  if(!ok)
  {
    filter->setErrorCondition(k_VgiParseError, "Error parsing size_t values from vgi file.");
    return {0, 0, 0};
  }
  vec[1] = tokens[3].toULongLong(&ok);
  if(!ok)
  {
    filter->setErrorCondition(k_VgiParseError, "Error parsing size_t values from vgi file.");
    return {0, 0, 0};
  }
  vec[2] = tokens[4].toULongLong(&ok);
  if(!ok)
  {
    filter->setErrorCondition(k_VgiParseError, "Error parsing size_t values from vgi file.");
    return {0, 0, 0};
  }

  return vec;
}

// -----------------------------------------------------------------------------
FloatVec2Type ParseFloatVec2(ImportVolumeGraphicsFile* filter, const QList<QByteArray>& tokens)
{
  FloatVec2Type vec = {0.0F, 0.0F};
  bool ok = false;

  vec[0] = tokens[2].toFloat(&ok);
  if(!ok)
  {
    filter->setErrorCondition(k_VgiParseError, "Error parsing float values from vgi file.");
    return {0.0F, 0.0F};
  }
  vec[1] = tokens[3].toFloat(&ok);
  if(!ok)
  {
    filter->setErrorCondition(k_VgiParseError, "Error parsing float values from vgi file.");
    return {0.0F, 0.0F};
  }

  return vec;
}

// -----------------------------------------------------------------------------
FloatVec3Type ParseFloatVec3(ImportVolumeGraphicsFile* filter, const QList<QByteArray>& tokens)
{
  FloatVec3Type vec = {0.0F, 0.0F, 0.0F};
  bool ok = false;

  vec[0] = tokens[2].toFloat(&ok);
  if(!ok)
  {
    filter->setErrorCondition(k_VgiParseError, "Error parsing float values from vgi file.");
    return {0.0F, 0.0F, 0.0F};
  }
  vec[1] = tokens[3].toFloat(&ok);
  if(!ok)
  {
    filter->setErrorCondition(k_VgiParseError, "Error parsing float values from vgi file.");
    return {0.0F, 0.0F, 0.0F};
  }
  vec[2] = tokens[4].toFloat(&ok);
  if(!ok)
  {
    filter->setErrorCondition(k_VgiParseError, "Error parsing float values from vgi file.");
    return {0.0F, 0.0F, 0.0F};
  }

  return vec;
}

// -----------------------------------------------------------------------------
template <typename T, size_t Dimension = 9>
std::array<T, Dimension> ParseFloatArray(ImportVolumeGraphicsFile* filter, const QList<QByteArray>& tokens)
{
  std::array<T, Dimension> vec;
  bool ok = false;

  for(size_t i = 0; i < Dimension; i++)
  {
    vec[i] = tokens[i + 2].toFloat(&ok);
    if(!ok)
    {
      filter->setErrorCondition(k_VgiParseError, "Error parsing float values from vgi file.");
      return {};
    }
  }

  return vec;
}

// -----------------------------------------------------------------------------
template <typename T, size_t Dimension = 9>
std::array<T, Dimension> ParseIntArray(ImportVolumeGraphicsFile* filter, const QList<QByteArray>& tokens)
{
  std::array<T, Dimension> vec;
  bool ok = false;

  for(size_t i = 0; i < Dimension; i++)
  {
    vec[i] = tokens[i + 2].toFloat(&ok);
    if(!ok)
    {
      filter->setErrorCondition(k_VgiParseError, "Error parsing float values from vgi file.");
      return {};
    }
  }

  return vec;
}

// -----------------------------------------------------------------------------
ImportVolumeGraphicsFile::FileBlock ParseFileBlock(ImportVolumeGraphicsFile* filter, QFile& in)
{
  ImportVolumeGraphicsFile::FileBlock fileBlock;

  int64_t currentPos = in.pos();
  QByteArray buf = in.readLine();
  buf = buf.trimmed();
  QList<QByteArray> tokens;
  bool ok = false;

  while(buf[0] != '[' && !in.atEnd())
  {
    tokens = buf.split(' ');
    if(tokens[0] == k_RegionOfInterestStart)
    {
      fileBlock.RegionOfInterestStart = ParseSizeVec3(filter, tokens);
    }
    else if(tokens[0] == k_RegionOfInterestEnd)
    {
      fileBlock.RegionOfInterestEnd = ParseSizeVec3(filter, tokens);
    }
    else if(tokens[0] == k_FileFormat)
    {
      fileBlock.FileFormat = tokens[2].toStdString();
    }
    else if(tokens[0] == k_Size)
    {
      fileBlock.Size = ParseSizeVec3(filter, tokens);
    }
    else if(tokens[0] == k_Name)
    {
      fileBlock.Name = tokens[2].toStdString();
    }
    else if(tokens[0] == k_Datatype)
    {
      fileBlock.Datatype = tokens[2].toStdString();
    }
    else if(tokens[0] == k_datarange)
    {
      fileBlock.datarange = ParseFloatVec2(filter, tokens);
    }
    else if(tokens[0] == k_BitsPerElement)
    {
      fileBlock.BitsPerElement = tokens[2].toInt(&ok);
    }

    currentPos = in.pos();
    buf = in.readLine(); // Read the next line
    buf = buf.trimmed();
  }
  in.seek(currentPos); // Put the file position back before the last read.
  return fileBlock;
}

// -----------------------------------------------------------------------------
ImportVolumeGraphicsFile::GeometryBlock ParseGeometryBlock(ImportVolumeGraphicsFile* filter, QFile& in)
{
  ImportVolumeGraphicsFile::GeometryBlock block;

  int64_t currentPos = in.pos();
  QByteArray buf = in.readLine();
  buf = buf.trimmed();
  QList<QByteArray> tokens;

  while(buf[0] != '[' && !in.atEnd())
  {
    tokens = buf.split(' ');
    if(tokens[0] == k_status)
    {
      block.status = tokens[2].toStdString();
    }
    else if(tokens[0] == k_relativeposition)
    {
      block.relativeposition = ParseFloatVec3(filter, tokens);
    }
    else if(tokens[0] == k_position)
    {
      block.position = ParseFloatVec3(filter, tokens);
    }
    else if(tokens[0] == k_resolution)
    {
      block.resolution = ParseFloatVec3(filter, tokens);
    }
    else if(tokens[0] == k_scale)
    {
      block.scale = ParseFloatVec3(filter, tokens);
    }
    else if(tokens[0] == k_center)
    {
      block.center = ParseFloatVec3(filter, tokens);
    }
    else if(tokens[0] == k_rotate)
    {
      block.rotate = ParseIntArray<int32_t, 9>(filter, tokens);
    }
    else if(tokens[0] == k_unit)
    {
      block.unit = tokens[2].toStdString();
    }
    currentPos = in.pos();
    buf = in.readLine(); // Read the next line
    buf = buf.trimmed();
  }
  in.seek(currentPos); // Put the file position back before the last read.
  return block;
}

// -----------------------------------------------------------------------------
void ImportVolumeGraphicsFile::readHeaderMetaData(ImageGeom::Pointer image)
{
  QByteArray buf;
  SizeVec3Type dims = {0, 0, 0};
  FloatVec3Type res = {0.0, 0.0, 0.0};
  FileBlock fileBlock;
  GeometryBlock geomBlock;

  while(!m_InHeaderStream.atEnd())
  {
    buf = m_InHeaderStream.readLine();
    buf = buf.trimmed();

    if(buf == k_file1Block)
    {
      fileBlock = ParseFileBlock(this, m_InHeaderStream);
      dims = fileBlock.Size;
    }

    else if(buf == k_geometryBlock)
    {
      geomBlock = ParseGeometryBlock(this, m_InHeaderStream);
      res = geomBlock.resolution;
    }
  }

  m_InHeaderStream.reset();

  image->setDimensions(dims);
  image->setSpacing(res);
  if(geomBlock.unit == k_Millimeter.toStdString())
  {
    image->setUnits(IGeometry::LengthUnit::Millimeter);
  }

  // Set the input .vol file from the header part.
  QFileInfo fiHdr(getVGHeaderFile());
  setVGDataFile(fiHdr.absolutePath() + QDir::separator() + QString::fromStdString(fileBlock.Name));
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AbstractFilter::Pointer ImportVolumeGraphicsFile::newFilterInstance(bool copyFilterParameters) const
{
  ImportVolumeGraphicsFile::Pointer filter = ImportVolumeGraphicsFile::New();
  if(copyFilterParameters)
  {
    copyFilterParameterInstanceVariables(filter.get());
  }
  return filter;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ImportVolumeGraphicsFile::getCompiledLibraryName() const
{
  return DREAM3DReviewConstants::DREAM3DReviewBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ImportVolumeGraphicsFile::getBrandingString() const
{
  return "DREAM3DReview";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ImportVolumeGraphicsFile::getFilterVersion() const
{
  QString version;
  QTextStream vStream(&version);
  vStream << DREAM3DReview::Version::Major() << "." << DREAM3DReview::Version::Minor() << "." << DREAM3DReview::Version::Patch();
  return version;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ImportVolumeGraphicsFile::getGroupName() const
{
  return SIMPL::FilterGroups::IOFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ImportVolumeGraphicsFile::getSubGroupName() const
{
  return SIMPL::FilterSubGroups::InputFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ImportVolumeGraphicsFile::getHumanLabel() const
{
  return "Import Volume Graphics File (.vgi/.vol)";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QUuid ImportVolumeGraphicsFile::getUuid() const
{
  return QUuid("{5fa10d81-94b4-582b-833f-8eabe659069e}");
}

// -----------------------------------------------------------------------------
ImportVolumeGraphicsFile::Pointer ImportVolumeGraphicsFile::NullPointer()
{
  return Pointer(static_cast<Self*>(nullptr));
}

// -----------------------------------------------------------------------------
std::shared_ptr<ImportVolumeGraphicsFile> ImportVolumeGraphicsFile::New()
{
  struct make_shared_enabler : public ImportVolumeGraphicsFile
  {
  };
  std::shared_ptr<make_shared_enabler> val = std::make_shared<make_shared_enabler>();
  val->setupFilterParameters();
  return val;
}

// -----------------------------------------------------------------------------
QString ImportVolumeGraphicsFile::getNameOfClass() const
{
  return QString("ImportVolumeGraphicsFile");
}

// -----------------------------------------------------------------------------
QString ImportVolumeGraphicsFile::ClassName()
{
  return QString("ImportVolumeGraphicsFile");
}

// -----------------------------------------------------------------------------
void ImportVolumeGraphicsFile::setVGDataFile(const QString& value)
{
  m_VGDataFile = value;
}

// -----------------------------------------------------------------------------
QString ImportVolumeGraphicsFile::getVGDataFile() const
{
  return m_VGDataFile;
}

// -----------------------------------------------------------------------------
void ImportVolumeGraphicsFile::setVGHeaderFile(const QString& value)
{
  m_VGHeaderFile = value;
}

// -----------------------------------------------------------------------------
QString ImportVolumeGraphicsFile::getVGHeaderFile() const
{
  return m_VGHeaderFile;
}

// -----------------------------------------------------------------------------
void ImportVolumeGraphicsFile::setDataContainerName(const QString& value)
{
  m_DataContainerName = value;
}

// -----------------------------------------------------------------------------
QString ImportVolumeGraphicsFile::getDataContainerName() const
{
  return m_DataContainerName;
}

// -----------------------------------------------------------------------------
void ImportVolumeGraphicsFile::setCellAttributeMatrixName(const QString& value)
{
  m_CellAttributeMatrixName = value;
}

// -----------------------------------------------------------------------------
QString ImportVolumeGraphicsFile::getCellAttributeMatrixName() const
{
  return m_CellAttributeMatrixName;
}

// -----------------------------------------------------------------------------
void ImportVolumeGraphicsFile::setDensityArrayName(const QString& value)
{
  m_DensityArrayName = value;
}

// -----------------------------------------------------------------------------
QString ImportVolumeGraphicsFile::getDensityArrayName() const
{
  return m_DensityArrayName;
}
