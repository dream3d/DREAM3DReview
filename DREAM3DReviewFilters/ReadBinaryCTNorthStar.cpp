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

#include "ReadBinaryCTNorthStar.h"

#include <tuple>

#include <QtCore/QDir>
#include <QtCore/QFileInfo>
#include <QtCore/QTextStream>

#include "SIMPLib/Common/Constants.h"
#include "SIMPLib/Common/ScopedFileMonitor.hpp"
#include "SIMPLib/DataContainers/DataContainer.h"
#include "SIMPLib/DataContainers/DataContainerArray.h"
#include "SIMPLib/FilterParameters/AbstractFilterParametersReader.h"
#include "SIMPLib/FilterParameters/ChoiceFilterParameter.h"
#include "SIMPLib/FilterParameters/InputFileFilterParameter.h"
#include "SIMPLib/FilterParameters/IntVec3FilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedBooleanFilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedPathCreationFilterParameter.h"
#include "SIMPLib/FilterParameters/PreflightUpdatedValueFilterParameter.h"
#include "SIMPLib/FilterParameters/SeparatorFilterParameter.h"
#include "SIMPLib/FilterParameters/StringFilterParameter.h"

#include "DREAM3DReview/DREAM3DReviewConstants.h"
#include "DREAM3DReview/DREAM3DReviewVersion.h"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
ReadBinaryCTNorthStar::ReadBinaryCTNorthStar()
{
  initialize();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
ReadBinaryCTNorthStar::~ReadBinaryCTNorthStar() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ReadBinaryCTNorthStar::setupFilterParameters()
{
  FilterParameterVectorType parameters;
  parameters.push_back(SIMPL_NEW_INPUT_FILE_FP("Input Header File", InputHeaderFile, FilterParameter::Parameter, ReadBinaryCTNorthStar, "*.nsihdr", "NSI header"));
  parameters.push_back(SIMPL_NEW_STRING_FP("Data Container", DataContainerName, FilterParameter::CreatedArray, ReadBinaryCTNorthStar));
  parameters.push_back(SeparatorFilterParameter::New("Cell Data", FilterParameter::CreatedArray));
  parameters.push_back(SIMPL_NEW_AM_WITH_LINKED_DC_FP("Cell Attribute Matrix", CellAttributeMatrixName, DataContainerName, FilterParameter::CreatedArray, ReadBinaryCTNorthStar));
  parameters.push_back(SIMPL_NEW_DA_WITH_LINKED_AM_FP("Density", DensityArrayName, DataContainerName, CellAttributeMatrixName, FilterParameter::CreatedArray, ReadBinaryCTNorthStar));

  QVector<QString> choices = IGeometry::GetAllLengthUnitStrings();
  parameters.push_back(SIMPL_NEW_CHOICE_FP("Length Unit", LengthUnit, FilterParameter::Parameter, ReadBinaryCTNorthStar, choices, false));

  PreflightUpdatedValueFilterParameter::Pointer param = SIMPL_NEW_PREFLIGHTUPDATEDVALUE_FP("Original Volume", VolumeDescription, FilterParameter::Parameter, ReadBinaryCTNorthStar);
  param->setReadOnly(true);
  parameters.push_back(param);
  param = SIMPL_NEW_PREFLIGHTUPDATEDVALUE_FP("Dat Files", DataFileInfo, FilterParameter::Parameter, ReadBinaryCTNorthStar);
  param->setReadOnly(true);
  parameters.push_back(param);

  QStringList linkedProps = {"StartVoxelCoord", "EndVoxelCoord", "ImportedVolumeDescription"};
  parameters.push_back(SIMPL_NEW_LINKED_BOOL_FP("Import Subvolume", ImportSubvolume, FilterParameter::Parameter, ReadBinaryCTNorthStar, linkedProps));

  parameters.push_back(SIMPL_NEW_INT_VEC3_FP("Starting XYZ Voxel", StartVoxelCoord, FilterParameter::Parameter, ReadBinaryCTNorthStar));
  parameters.push_back(SIMPL_NEW_INT_VEC3_FP("Ending XYZ Voxel", EndVoxelCoord, FilterParameter::Parameter, ReadBinaryCTNorthStar));

  param = SIMPL_NEW_PREFLIGHTUPDATEDVALUE_FP("Imported Volume", ImportedVolumeDescription, FilterParameter::Parameter, ReadBinaryCTNorthStar);
  param->setReadOnly(true);
  parameters.push_back(param);

  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ReadBinaryCTNorthStar::initialize()
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
void ReadBinaryCTNorthStar::dataCheck()
{
  clearErrorCode();
  clearWarningCode();
  initialize();

  QFileInfo fiHdr(getInputHeaderFile());

  if(getInputHeaderFile().isEmpty())
  {
    QString ss = QObject::tr("The input header file must be set");
    setErrorCondition(-38700, ss);
    return;
  }
  if(!fiHdr.exists())
  {
    QString ss = QObject::tr("The input header file does not exist");
    setErrorCondition(-38701, ss);
    return;
  }

  if(m_InHeaderStream.isOpen())
  {
    m_InHeaderStream.close();
  }

  m_InHeaderStream.setFileName(getInputHeaderFile());

  if(!m_InHeaderStream.open(QIODevice::ReadOnly | QIODevice::Text))
  {
    QString ss = QObject::tr("Error opening input file: %1").arg(getInputHeaderFile());
    setErrorCondition(-38702, ss);
    return;
  }

  if(readHeaderMetaData() < 0)
  {
    QString ss = QObject::tr("Error reading input header file: %1").arg(getInputHeaderFile());
    setErrorCondition(-38703, ss);
    return;
  }

  // Sanity check the Start/End Voxels
  if(getImportSubvolume())
  {
    if(m_StartVoxelCoord[0] > m_EndVoxelCoord[0])
    {
      QString ss = QObject::tr("Starting X Voxel > Ending X Voxel (%1 > %2)").arg(m_StartVoxelCoord[0]).arg(m_EndVoxelCoord[0]);
      setErrorCondition(-38712, ss);
    }
    if(m_StartVoxelCoord[1] > m_EndVoxelCoord[1])
    {
      QString ss = QObject::tr("Starting Y Voxel > Ending Y Voxel (%1 > %2)").arg(m_StartVoxelCoord[1]).arg(m_EndVoxelCoord[1]);
      setErrorCondition(-38713, ss);
    }
    if(m_StartVoxelCoord[2] > m_EndVoxelCoord[2])
    {
      QString ss = QObject::tr("Starting Z Voxel > Ending Z Voxel (%1 > %2)").arg(m_StartVoxelCoord[2]).arg(m_EndVoxelCoord[2]);
      setErrorCondition(-38714, ss);
    }

    if(m_StartVoxelCoord[0] < 0 || m_StartVoxelCoord[1] < 0 || m_StartVoxelCoord[1] < 0)
    {
      QString ss = QObject::tr("Start Voxel < ZERO.(%1, %2, %3)").arg(m_StartVoxelCoord[0]).arg(m_StartVoxelCoord[1]).arg(m_StartVoxelCoord[2]);
      setErrorCondition(-38715, ss);
    }

    SizeVec3Type origDims = m_OriginalVolume->getDimensions();
    if(m_EndVoxelCoord[0] >= origDims[0])
    {
      QString ss = QObject::tr("Ending X Voxel > Original Volume Dimension (%1 >= %2)").arg(m_EndVoxelCoord[0]).arg(origDims[0]);
      setErrorCondition(-38716, ss);
    }
    if(m_EndVoxelCoord[1] >= origDims[1])
    {
      QString ss = QObject::tr("Ending Y Voxel > Original Volume Dimension (%1 >= %2)").arg(m_EndVoxelCoord[1]).arg(origDims[1]);
      setErrorCondition(-38717, ss);
    }
    if(m_EndVoxelCoord[2] >= origDims[2])
    {
      QString ss = QObject::tr("Ending Z Voxel > Original Volume Dimension (%1 >= %2)").arg(m_EndVoxelCoord[2]).arg(origDims[2]);
      setErrorCondition(-38718, ss);
    }
  }
  if(getErrorCode() < 0)
  {
    return;
  }

  QFileInfo headerPath(getInputHeaderFile());
  QDir headerDir = headerPath.absoluteDir();

  for(auto& p : m_DataFiles)
  {
    QFileInfo fi(headerDir.absolutePath() + "/" + p.first);
    p.first = headerDir.absolutePath() + "/" + p.first;

    if(!fi.exists())
    {
      QString ss = QObject::tr("The input binary CT file: %1 could not be found in the path: %2. The input binary CT file name was obtained from the provided header.  Please ensure the file is in "
                               "the same path as the header and the name has not been altered.")
                       .arg(p.first)
                       .arg(headerDir.absolutePath());
      setErrorCondition(-38704, ss);
    }
  }

  DataContainer::Pointer m = getDataContainerArray()->createNonPrereqDataContainer(this, getDataContainerName());

  if(getErrorCode() < 0)
  {
    return;
  }

  m->setGeometry(m_ImportedVolume);

  SizeVec3Type dims = m_ImportedVolume->getDimensions();

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
int32_t ReadBinaryCTNorthStar::sanityCheckFileSizeVersusAllocatedSize(size_t allocatedBytes, size_t fileSize)
{
  if(fileSize < allocatedBytes)
  {
    return -1;
  }
  if(fileSize > allocatedBytes)
  {
    return 1;
  }
  // File Size and Allocated Size are equal so we  are good to go
  return 0;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int32_t ReadBinaryCTNorthStar::readBinaryCTFiles()
{

  SizeVec3Type origDims = m_OriginalVolume->getDimensions();
  SizeVec3Type importedDims = m_ImportedVolume->getDimensions();

  size_t deltaX = m_EndVoxelCoord[0] - m_StartVoxelCoord[0] + 1;

  int32_t error = 0;
  size_t zShift = 0;
  int32_t fileIndex = 1;

  auto densityPtr = m_DensityPtr.lock();
  FloatArrayType& density = *densityPtr;
  density.initializeWithValue(0xABCDEF);
  float* ptr = density.getPointer(0);

  for(const auto& dataFileInput : m_DataFiles)
  {
    QFileInfo fi(dataFileInput.first);
    size_t filesize = static_cast<size_t>(fi.size());
    // allocated bytes should be the x * y dims * number of slices in the current data file....not necessarily the size of the whole density array
    size_t allocatedBytes = origDims[0] * origDims[1] * dataFileInput.second * sizeof(float);

    error = sanityCheckFileSizeVersusAllocatedSize(allocatedBytes, filesize);

    if(error < 0)
    {
      QString msg;
      QTextStream ss(&msg);
      ss << "Binary file size " << filesize << " is smaller than the number of allocated bytes: " << allocatedBytes;
      setErrorCondition(-38705, msg);
      return getErrorCode();
    }

    FILE* f = fopen(dataFileInput.first.toLatin1().data(), "rb");
    if(nullptr == f)
    {
      QString ss = QObject::tr("Error opening binary input file: %1").arg(dataFileInput.first);
      setErrorCondition(-38706, ss);
      return getErrorCode();
    }

    ScopedFileMonitor monitor(f);

    size_t fileZSlice = 0;
    size_t index = 0;

    for(size_t z = zShift; z < (zShift + dataFileInput.second); z++)
    {
      if(m_ImportSubvolume && (z < m_StartVoxelCoord[2] || z > m_EndVoxelCoord[2]))
      {
        fileZSlice++;
        continue;
      }
      QString ss = QObject::tr("Importing Data || Data File: %1 || Importing Slice %2").arg(dataFileInput.first).arg(z);
      notifyStatusMessage(ss);
      for(size_t y = 0; y < origDims[1]; y++)
      {
        if(m_ImportSubvolume && (y < m_StartVoxelCoord[1] || y > m_EndVoxelCoord[1]))
        {
          continue;
        }

        size_t fpOffset = ((origDims[1] * origDims[0] * fileZSlice) + (origDims[0] * y) + m_StartVoxelCoord[0]) * sizeof(float);
        if(fseek(f, fpOffset, SEEK_SET) != 0)
        {
          QString ss = QObject::tr("Could not seek to postion %1 in file %2").arg(fpOffset).arg(dataFileInput.first);
          setErrorCondition(-38707, ss);
          fclose(f);
          return getErrorCode();
        }

        index = (importedDims[0] * importedDims[1] * (z - m_StartVoxelCoord[2])) + (importedDims[0] * (y - m_StartVoxelCoord[1])) + (0);
        ptr = density.getPointer(index);
        if(fread(ptr, sizeof(float), deltaX, f) != deltaX)
        {
          QString ss = QObject::tr("Error reading file at position %1 in file %2").arg(fpOffset).arg(dataFileInput.first);
          setErrorCondition(-387008, ss);
          fclose(f);
          return getErrorCode();
        }
      }
      fileZSlice++;
    }
    zShift += dataFileInput.second;
    fileIndex++;
    fclose(f);
  }

  return error;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int32_t ReadBinaryCTNorthStar::readHeaderMetaData()
{
  int32_t error = 0;
  m_OriginalVolume = ImageGeom::NullPointer();
  m_ImportedVolume = ImageGeom::NullPointer();
  m_DataFiles.clear();

  QByteArray buf;
  QList<QByteArray> tokens;
  bool ok = false;

  SizeVec3Type voxels = {0, 0, 0};
  FloatVec3Type maxLocation = {0.0f, 0.0f, 0.0f};
  FloatVec3Type minLocation = {0.0f, 0.0f, 0.0f};
  FloatVec3Type spacing = {0.0f, 0.0f, 0.0f};
  bool done = false;

  while(!m_InHeaderStream.atEnd() && !done)
  {
    buf = m_InHeaderStream.readLine();
    buf = buf.trimmed();
    tokens = buf.split(' ');
    for(qint32 i = 0; i < tokens.size(); i++)
    {
      if(tokens[i] == "<Voxels>")
      {
        voxels[1] = tokens[i + 1].toULongLong(&ok);
        if(!ok)
        {
          return -1;
        }
        voxels[2] = tokens[i + 2].toULongLong(&ok);
        if(!ok)
        {
          return -1;
        }
        voxels[0] = tokens[i + 3].toULongLong(&ok);
        if(!ok)
        {
          return -1;
        }
        done = true;
        break;
      }
    }
  }

  m_InHeaderStream.reset();
  done = false;

  while(!m_InHeaderStream.atEnd() && !done)
  {
    buf = m_InHeaderStream.readLine();
    buf = buf.trimmed();
    tokens = buf.split(' ');
    for(qint32 i = 0; i < tokens.size(); i++)
    {
      if(tokens[i] == "<Location>")
      {
        buf = m_InHeaderStream.readLine();
        buf = buf.trimmed();
        tokens = buf.split(' ');
        for(qint32 j = 0; j < tokens.size(); j++)
        {
          if(tokens[j] == "<Min>")
          {
            minLocation[1] = tokens[j + 1].toFloat(&ok);
            if(!ok)
            {
              return -1;
            }
            minLocation[2] = tokens[j + 2].toFloat(&ok);
            if(!ok)
            {
              return -1;
            }
            minLocation[0] = tokens[j + 3].toFloat(&ok);
            if(!ok)
            {
              return -1;
            }
          }
        }

        buf = m_InHeaderStream.readLine();
        buf = buf.trimmed();
        tokens = buf.split(' ');
        for(qint32 j = 0; j < tokens.size(); j++)
        {
          if(tokens[j] == "<Max>")
          {
            maxLocation[1] = tokens[j + 1].toFloat(&ok);
            if(!ok)
            {
              return -1;
            }
            maxLocation[2] = tokens[j + 2].toFloat(&ok);
            if(!ok)
            {
              return -1;
            }
            maxLocation[0] = tokens[j + 3].toFloat(&ok);
            if(!ok)
            {
              return -1;
            }
          }
        }
        done = true;
        break;
      }
    }
  }

  m_InHeaderStream.reset();
  done = false;
  // Look for all the input files
  QString contents = m_InHeaderStream.readAll();
  QStringList list = contents.split(QRegExp("\\n"));
  QStringListIterator sourceLines(list);
  bool filesSectionActive = false;
  QString inputDataFile;
  QString nbSlices;
  while(sourceLines.hasNext())
  {
    QString line = sourceLines.next();
    std::pair<std::string, int64_t> input;

    if(line.trimmed().startsWith("<Files>"))
    {
      filesSectionActive = true;
    }
    else if(line.trimmed().startsWith("<Name>") && filesSectionActive)
    {
      inputDataFile = line;
      inputDataFile = inputDataFile.replace("<Name>", "").trimmed();
    }
    else if(line.trimmed().startsWith("<NbSlices>") && filesSectionActive)
    {
      nbSlices = line;
      nbSlices = nbSlices.replace("<NbSlices>", "").trimmed();
    }
    else if(line.trimmed().startsWith("</Files>"))
    {
      filesSectionActive = false;
      break; // We are done getting the file list
    }

    if(!inputDataFile.isEmpty() && !nbSlices.isEmpty())
    {
      std::pair<QString, int64_t> p(inputDataFile, nbSlices.trimmed().toULongLong(&ok));
      m_DataFiles.push_back(p);
      inputDataFile.clear();
      nbSlices.clear();
    }
  }

  for(size_t i = 0; i < 3; i++)
  {
    spacing[i] = (maxLocation[i] - minLocation[i]) / voxels[i];
  }

  m_OriginalVolume = ImageGeom::CreateGeometry(SIMPL::Geometry::ImageGeometry);

  m_OriginalVolume->setDimensions(voxels);
  m_OriginalVolume->setSpacing(spacing);
  m_OriginalVolume->setOrigin(minLocation);
  m_OriginalVolume->setUnits(static_cast<IGeometry::LengthUnit>(getLengthUnit()));

  if(m_ImportSubvolume)
  {
    // Now figure out the Imported Volume Information
    m_ImportedVolume = ImageGeom::CreateGeometry(SIMPL::Geometry::ImageGeometry);
    // Origin
    minLocation[0] = minLocation[0] + (m_StartVoxelCoord[0] * spacing[0]);
    minLocation[1] = minLocation[1] + (m_StartVoxelCoord[1] * spacing[1]);
    minLocation[2] = minLocation[2] + (m_StartVoxelCoord[2] * spacing[2]);
    m_ImportedVolume->setOrigin(minLocation);
    // Spacing
    m_ImportedVolume->setSpacing(spacing);
    // Dimensions
    SizeVec3Type dims = {static_cast<size_t>(m_EndVoxelCoord[0] - m_StartVoxelCoord[0] + 1), static_cast<size_t>(m_EndVoxelCoord[1] - m_StartVoxelCoord[1] + 1),
                         static_cast<size_t>(m_EndVoxelCoord[2] - m_StartVoxelCoord[2] + 1)};
    m_ImportedVolume->setDimensions(dims);
    m_ImportedVolume->setUnits(static_cast<IGeometry::LengthUnit>(getLengthUnit()));
  }
  else
  {
    m_ImportedVolume = m_OriginalVolume;
    m_StartVoxelCoord = {0, 0, 0};
    m_EndVoxelCoord = {static_cast<int32_t>(voxels[0] - 1), static_cast<int32_t>(voxels[1] - 1), static_cast<int32_t>(voxels[2] - 1)};
  }

  return error;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ReadBinaryCTNorthStar::execute()
{
  dataCheck();
  if(getErrorCode() < 0)
  {
    return;
  }

  if(readBinaryCTFiles() < 0)
  {
    notifyStatusMessage("Error during import");
    return;
  }

  notifyStatusMessage("Complete");
}

// -----------------------------------------------------------------------------
QString GenerateGeometryInfoString(ImageGeom* geom)
{
  QString desc;
  QTextStream ss(&desc);
  SizeVec3Type dims = geom->getDimensions();
  FloatVec3Type origin = geom->getOrigin();
  FloatVec3Type spacing = geom->getSpacing();
  QString lengthUnit = IGeometry::LengthUnitToString(static_cast<IGeometry::LengthUnit>(geom->getUnits()));

  ss << "X Range: " << origin[0] << " to " << (origin[0] + (dims[0] * spacing[0])) << " (Delta: " << (dims[0] * spacing[0]) << " " << lengthUnit << ") " << 0 << "-" << dims[0] - 1 << " Voxels\n";
  ss << "Y Range: " << origin[1] << " to " << (origin[1] + (dims[1] * spacing[1])) << " (Delta: " << (dims[1] * spacing[1]) << " " << lengthUnit << ") " << 0 << "-" << dims[1] - 1 << " Voxels\n";
  ss << "Z Range: " << origin[2] << " to " << (origin[2] + (dims[2] * spacing[2])) << " (Delta: " << (dims[2] * spacing[2]) << " " << lengthUnit << ") " << 0 << "-" << dims[2] - 1 << " Voxels\n";
  return desc;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ReadBinaryCTNorthStar::getVolumeDescription()
{
  QString desc = QString("Geometry Not initialized.");
  if(m_OriginalVolume.get() != nullptr)
  {
    desc = GenerateGeometryInfoString(m_OriginalVolume.get());
  }
  return desc;
}

// -----------------------------------------------------------------------------
QString ReadBinaryCTNorthStar::getDataFileInfo()
{
  QString desc;
  QTextStream out(&desc);
  for(const auto& dataFileInput : m_DataFiles)
  {
    QFileInfo fi(dataFileInput.first);
    out << fi.fileName();
    if(fi.exists())
    {
      out << " [ok]";
    }
    else
    {
      out << " [!]";
    }
    out << "\n";
  }
  return desc;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ReadBinaryCTNorthStar::getImportedVolumeDescription()
{
  QString desc = QString("Geometry Not initialized.");
  if(m_ImportedVolume.get() != nullptr)
  {
    desc = GenerateGeometryInfoString(m_ImportedVolume.get());
  }
  return desc;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AbstractFilter::Pointer ReadBinaryCTNorthStar::newFilterInstance(bool copyFilterParameters) const
{
  ReadBinaryCTNorthStar::Pointer filter = ReadBinaryCTNorthStar::New();
  if(copyFilterParameters)
  {
    copyFilterParameterInstanceVariables(filter.get());
  }
  return filter;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ReadBinaryCTNorthStar::getCompiledLibraryName() const
{
  return DREAM3DReviewConstants::DREAM3DReviewBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ReadBinaryCTNorthStar::getBrandingString() const
{
  return "DREAM3DReview";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ReadBinaryCTNorthStar::getFilterVersion() const
{
  QString version;
  QTextStream vStream(&version);
  vStream << DREAM3DReview::Version::Major() << "." << DREAM3DReview::Version::Minor() << "." << DREAM3DReview::Version::Patch();
  return version;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ReadBinaryCTNorthStar::getGroupName() const
{
  return SIMPL::FilterGroups::IOFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ReadBinaryCTNorthStar::getSubGroupName() const
{
  return SIMPL::FilterSubGroups::InputFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ReadBinaryCTNorthStar::getHumanLabel() const
{
  return "Import North Star Imaging CT (.nsihdr/.nsidat)";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QUuid ReadBinaryCTNorthStar::getUuid() const
{
  return QUuid("{f2259481-5011-5f22-9fcb-c92fb6f8be10}");
}

// -----------------------------------------------------------------------------
ReadBinaryCTNorthStar::Pointer ReadBinaryCTNorthStar::NullPointer()
{
  return Pointer(static_cast<Self*>(nullptr));
}

// -----------------------------------------------------------------------------
std::shared_ptr<ReadBinaryCTNorthStar> ReadBinaryCTNorthStar::New()
{
  struct make_shared_enabler : public ReadBinaryCTNorthStar
  {
  };
  std::shared_ptr<make_shared_enabler> val = std::make_shared<make_shared_enabler>();
  val->setupFilterParameters();
  return val;
}

// -----------------------------------------------------------------------------
QString ReadBinaryCTNorthStar::getNameOfClass() const
{
  return QString("ReadBinaryCTNorthStar");
}

// -----------------------------------------------------------------------------
QString ReadBinaryCTNorthStar::ClassName()
{
  return QString("ReadBinaryCTNorthStar");
}

// -----------------------------------------------------------------------------
void ReadBinaryCTNorthStar::setInputHeaderFile(const QString& value)
{
  m_InputHeaderFile = value;
}

// -----------------------------------------------------------------------------
QString ReadBinaryCTNorthStar::getInputHeaderFile() const
{
  return m_InputHeaderFile;
}

// -----------------------------------------------------------------------------
void ReadBinaryCTNorthStar::setDataContainerName(const QString& value)
{
  m_DataContainerName = value;
}

// -----------------------------------------------------------------------------
QString ReadBinaryCTNorthStar::getDataContainerName() const
{
  return m_DataContainerName;
}

// -----------------------------------------------------------------------------
void ReadBinaryCTNorthStar::setCellAttributeMatrixName(const QString& value)
{
  m_CellAttributeMatrixName = value;
}

// -----------------------------------------------------------------------------
QString ReadBinaryCTNorthStar::getCellAttributeMatrixName() const
{
  return m_CellAttributeMatrixName;
}

// -----------------------------------------------------------------------------
void ReadBinaryCTNorthStar::setDensityArrayName(const QString& value)
{
  m_DensityArrayName = value;
}

// -----------------------------------------------------------------------------
QString ReadBinaryCTNorthStar::getDensityArrayName() const
{
  return m_DensityArrayName;
}

// -----------------------------------------------------------------------------
void ReadBinaryCTNorthStar::setLengthUnit(int32_t value)
{
  m_LengthUnit = value;
}

// -----------------------------------------------------------------------------
int32_t ReadBinaryCTNorthStar::getLengthUnit() const
{
  return m_LengthUnit;
}

// -----------------------------------------------------------------------------
void ReadBinaryCTNorthStar::setStartVoxelCoord(const IntVec3Type& value)
{
  m_StartVoxelCoord = value;
}

// -----------------------------------------------------------------------------
IntVec3Type ReadBinaryCTNorthStar::getStartVoxelCoord() const
{
  return m_StartVoxelCoord;
}

// -----------------------------------------------------------------------------
void ReadBinaryCTNorthStar::setEndVoxelCoord(const IntVec3Type& value)
{
  m_EndVoxelCoord = value;
}

// -----------------------------------------------------------------------------
IntVec3Type ReadBinaryCTNorthStar::getEndVoxelCoord() const
{
  return m_EndVoxelCoord;
}

// -----------------------------------------------------------------------------
void ReadBinaryCTNorthStar::setImportSubvolume(bool value)
{
  m_ImportSubvolume = value;
}

// -----------------------------------------------------------------------------
bool ReadBinaryCTNorthStar::getImportSubvolume() const
{
  return m_ImportSubvolume;
}
