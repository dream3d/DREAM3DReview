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
#include <memory>

#include "ImportVolumeGraphicsFile.h"

#include <QtCore/QFileInfo>

#include <QtCore/QTextStream>

#include "SIMPLib/Common/Constants.h"

#include "SIMPLib/Common/ScopedFileMonitor.hpp"
#include "SIMPLib/FilterParameters/AbstractFilterParametersReader.h"
#include "SIMPLib/FilterParameters/InputFileFilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedPathCreationFilterParameter.h"
#include "SIMPLib/FilterParameters/SeparatorFilterParameter.h"
#include "SIMPLib/FilterParameters/StringFilterParameter.h"
#include "SIMPLib/DataContainers/DataContainerArray.h"
#include "SIMPLib/DataContainers/DataContainer.h"

#include "DREAM3DReview/DREAM3DReviewConstants.h"
#include "DREAM3DReview/DREAM3DReviewVersion.h"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
ImportVolumeGraphicsFile::ImportVolumeGraphicsFile()
: m_InputFile("")
, m_InputHeaderFile("")
, m_DataContainerName("ImageDataContainer")
, m_CellAttributeMatrixName("CellData")
, m_DensityArrayName("Density")
{
}

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
  parameters.push_back(SIMPL_NEW_INPUT_FILE_FP("Input CT File", InputFile, FilterParameter::Parameter, ImportVolumeGraphicsFile, "*.vol", "Voxel data"));
  parameters.push_back(SIMPL_NEW_INPUT_FILE_FP("Input Header File", InputHeaderFile, FilterParameter::Parameter, ImportVolumeGraphicsFile, "*.vgi", "NSI header"));
  parameters.push_back(SIMPL_NEW_STRING_FP("Data Container", DataContainerName, FilterParameter::CreatedArray, ImportVolumeGraphicsFile));
  parameters.push_back(SeparatorFilterParameter::New("Cell Data", FilterParameter::CreatedArray));
  parameters.push_back(SIMPL_NEW_AM_WITH_LINKED_DC_FP("Cell Attribute Matrix", CellAttributeMatrixName, DataContainerName, FilterParameter::CreatedArray, ImportVolumeGraphicsFile));
  parameters.push_back(SIMPL_NEW_DA_WITH_LINKED_AM_FP("Density", DensityArrayName, DataContainerName, CellAttributeMatrixName, FilterParameter::CreatedArray, ImportVolumeGraphicsFile));
  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ImportVolumeGraphicsFile::readFilterParameters(AbstractFilterParametersReader* reader, int index)
{
  reader->openFilterGroup(this, index);
  setInputFile(reader->readString("InputFile", getInputFile()));
  setInputHeaderFile(reader->readString("InputHeaderFile", getInputHeaderFile()));
  setDataContainerName(reader->readString("DataContainerName", getDataContainerName()));
  setCellAttributeMatrixName(reader->readString("CellAttributeMatrixName", getCellAttributeMatrixName()));
  setDensityArrayName(reader->readString("DensityArrayName", getDensityArrayName()));
  reader->closeFilterGroup();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ImportVolumeGraphicsFile::initialize()
{
  if(m_InHeaderStream.isOpen())
    m_InHeaderStream.close();
  if(m_InStream.isOpen())
    m_InStream.close();
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

  QFileInfo fi(getInputFile());

  if(getInputFile().isEmpty() == true)
  {
    QString ss = QObject::tr("The input binary CT voxel file must be set");
    setErrorCondition(-387, ss);
  }
  else if(fi.exists() == false)
  {
    QString ss = QObject::tr("The input binary CT voxel file does not exist");
    setErrorCondition(-388, ss);
  }

  QFileInfo fiHdr(getInputHeaderFile());

  if(getInputHeaderFile().isEmpty() == true)
  {
    QString ss = QObject::tr("The input header file must be set");
    setErrorCondition(-387, ss);
  }
  else if(fiHdr.exists() == false)
  {
    QString ss = QObject::tr("The input header file does not exist");
    setErrorCondition(-388, ss);
  }

  if(getErrorCode() < 0)
  {
    return;
  }

  if(m_InHeaderStream.isOpen() == true)
  {
    m_InHeaderStream.close();
  }

  m_InHeaderStream.setFileName(getInputHeaderFile());

  if(!m_InHeaderStream.open(QIODevice::ReadOnly | QIODevice::Text))
  {
    QString ss = QObject::tr("Error opening input header file: %1").arg(getInputHeaderFile());
    setErrorCondition(-100, ss);
    return;
  }

  ImageGeom::Pointer image = ImageGeom::CreateGeometry(SIMPL::Geometry::ImageGeometry);

  err = readHeaderMetaData(image);

  if(err < 0)
  {
    QString ss = QObject::tr("Error reading input header file: %1").arg(getInputHeaderFile());
    setErrorCondition(-100, ss);
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

  m_DensityPtr =
      getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<float>>(this, path, 0, cDims); /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
  if(nullptr != m_DensityPtr.lock().get()) /* Validate the Weak Pointer wraps a non-nullptr pointer to a DataArray<T> object */
  {
    m_Density = m_DensityPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int32_t ImportVolumeGraphicsFile::sanityCheckFileSizeVersusAllocatedSize(size_t allocatedBytes, size_t fileSize)
{
  if(fileSize < allocatedBytes)
  {
    return -1;
  }
  else if(fileSize > allocatedBytes)
  {
    return 1;
  }
  // File Size and Allocated Size are equal so we  are good to go
  return 0;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int32_t ImportVolumeGraphicsFile::readBinaryCTFile(size_t dims[3])
{
  int32_t error = 0;

  QFileInfo fi(getInputFile());
  size_t filesize = static_cast<size_t>(fi.size());
  size_t allocatedBytes = m_DensityPtr.lock()->getSize() * sizeof(float);

  error = sanityCheckFileSizeVersusAllocatedSize(allocatedBytes, filesize);

  if(error < 0)
  {
    QString ss = QObject::tr("Binary file size is smaller than the number of allocated bytes");
    setErrorCondition(-100, ss);
    return getErrorCode();
  }

  FILE* f = fopen(getInputFile().toLatin1().data(), "rb");
  if(nullptr == f)
  {
    QString ss = QObject::tr("Error opening binary input file: %1").arg(getInputFile());
    setErrorCondition(-100, ss);
    return getErrorCode();
  }

  ScopedFileMonitor monitor(f);

  float* ptr = m_DensityPtr.lock()->getPointer(0);

  size_t index = 0;
  size_t totalPoints = dims[0] * dims[1] * dims[2];
  int64_t progIncrement = static_cast<int64_t>(totalPoints / 100);
  int64_t prog = 1;
  int64_t progressInt = 0;
  int64_t counter = 0;

  for(size_t z = 0; z < dims[2]; z++)
  {
    for(size_t y = 0; y < dims[1]; y++)
    {
      for(size_t x = 0; x < dims[0]; x++)
      {
        // Compute the proper index into the final storage array
        index = (dims[0] * dims[1] * z) + (dims[0] * (dims[1] - y - 1)) + (dims[0] - x - 1);
        fread(ptr + index, sizeof(float), 1, f);

        // Check for progress
        if(counter > prog)
        {
          progressInt = static_cast<int64_t>((static_cast<float>(counter) / totalPoints) * 100.0f);
          QString ss = QObject::tr("Reading Data || %1% Completed").arg(progressInt);
          notifyStatusMessage(ss);
          prog = prog + progIncrement;
        }
        counter++;
      }
    }
  }

  return error;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int32_t ImportVolumeGraphicsFile::readHeaderMetaData(ImageGeom::Pointer image)
{
  int32_t error = 0;

  QByteArray buf;
  QList<QByteArray> tokens;
  bool ok = false;

  size_t dims[3] = {0, 0, 0};
  float res[3] = {0.0, 0.0, 0.0};
  bool done = false;

  while(m_InHeaderStream.atEnd() == false && done == false)
  {
    buf = m_InHeaderStream.readLine();
    buf = buf.trimmed();
    tokens = buf.split(' ');
    for(qint32 i = 0; i < tokens.size(); i++)
    {
      if(tokens[i] == "[representation]")
      {
        buf = m_InHeaderStream.readLine();
        buf = buf.trimmed();
        tokens = buf.split(' ');
        for(qint32 j = 0; j < tokens.size(); j++)
        {
          if(tokens[j] == "size")
          {
            dims[0] = tokens[j + 2].toULongLong(&ok);
            if(!ok)
            {
              return -1;
            }
            dims[1] = tokens[j + 3].toULongLong(&ok);
            if(!ok)
            {
              return -1;
            }
            dims[2] = tokens[j + 4].toULongLong(&ok);
            if(!ok)
            {
              return -1;
            }
            done = true;
            break;
          }
        }
      }
    }
  }

  m_InHeaderStream.reset();
  done = false;

  while(m_InHeaderStream.atEnd() == false && done == false)
  {
    buf = m_InHeaderStream.readLine();
    buf = buf.trimmed();
    tokens = buf.split(' ');
    for(qint32 i = 0; i < tokens.size(); i++)
    {
      if(tokens[i] == "[geometry]")
      {
        bool done2 = false;
        while(m_InHeaderStream.atEnd() == false && done2 == false)
        {
          buf = m_InHeaderStream.readLine();
          buf = buf.trimmed();
          tokens = buf.split(' ');
          for(qint32 j = 0; j < tokens.size(); j++)
          {
            if(tokens[j] == "resolution")
            {
              res[0] = tokens[j + 2].toFloat(&ok);
              if(!ok)
              {
                return -1;
              }
              res[1] = tokens[j + 3].toFloat(&ok);
              if(!ok)
              {
                return -1;
              }
              res[2] = tokens[j + 4].toFloat(&ok);
              if(!ok)
              {
                return -1;
              }
              done2 = true;
              break;
            }
          }
        }
        if(done2 == true)
        {
          done = true;
          break;
        }
      }
    }
  }

  image->setDimensions(dims);
  image->setSpacing(res);

  return error;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ImportVolumeGraphicsFile::execute()
{
  clearErrorCode();
  clearWarningCode();
  dataCheck();
  if(getErrorCode() < 0)
  {
    return;
  }

  int32_t err = 0;

  ImageGeom::Pointer image = getDataContainerArray()->getDataContainer(getDataContainerName())->getGeometryAs<ImageGeom>();

  SizeVec3Type dims = image->getDimensions();

  err = readBinaryCTFile(dims.data());

  if(err < 0)
  {
    QString ss = QObject::tr("Error reading binary input file");
    setErrorCondition(-100, ss);
  }

  notifyStatusMessage("Complete");
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
  return "Import Binary CT File (XRadia)";
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
void ImportVolumeGraphicsFile::setInputFile(const QString& value)
{
  m_InputFile = value;
}

// -----------------------------------------------------------------------------
QString ImportVolumeGraphicsFile::getInputFile() const
{
  return m_InputFile;
}

// -----------------------------------------------------------------------------
void ImportVolumeGraphicsFile::setInputHeaderFile(const QString& value)
{
  m_InputHeaderFile = value;
}

// -----------------------------------------------------------------------------
QString ImportVolumeGraphicsFile::getInputHeaderFile() const
{
  return m_InputHeaderFile;
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
