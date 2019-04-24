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
#include "ReadBinaryCTNorthStar.h"

#include <QtCore/QDir>
#include <QtCore/QFileInfo>

#include "SIMPLib/Common/Constants.h"
#include "SIMPLib/Common/ScopedFileMonitor.hpp"
#include "SIMPLib/FilterParameters/AbstractFilterParametersReader.h"
#include "SIMPLib/FilterParameters/InputFileFilterParameter.h"
#include "SIMPLib/FilterParameters/SeparatorFilterParameter.h"
#include "SIMPLib/FilterParameters/StringFilterParameter.h"

#include "DREAM3DReview/DREAM3DReviewConstants.h"
#include "DREAM3DReview/DREAM3DReviewVersion.h"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
ReadBinaryCTNorthStar::ReadBinaryCTNorthStar()
: m_InputFiles(0, "")
, m_SlicesPerFile(0, 0)
, m_InputHeaderFile("")
, m_DataContainerName("ImageDataContainer")
, m_CellAttributeMatrixName("CellData")
, m_DensityArrayName("Density")
{
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
  parameters.push_back(SIMPL_NEW_STRING_FP("Cell Attribute Matrix", CellAttributeMatrixName, FilterParameter::CreatedArray, ReadBinaryCTNorthStar));
  parameters.push_back(SIMPL_NEW_STRING_FP("Density", DensityArrayName, FilterParameter::CreatedArray, ReadBinaryCTNorthStar));
  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ReadBinaryCTNorthStar::readFilterParameters(AbstractFilterParametersReader* reader, int index)
{
  reader->openFilterGroup(this, index);
  setInputHeaderFile(reader->readString("InputHeaderFile", getInputHeaderFile()));
  setDataContainerName(reader->readString("DataContainerName", getDataContainerName()));
  setCellAttributeMatrixName(reader->readString("CellAttributeMatrixName", getCellAttributeMatrixName()));
  setDensityArrayName(reader->readString("DensityArrayName", getDensityArrayName()));
  reader->closeFilterGroup();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ReadBinaryCTNorthStar::initialize()
{
  if(m_InHeaderStream.isOpen())
    m_InHeaderStream.close();
  if(m_InStream.isOpen())
    m_InStream.close();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ReadBinaryCTNorthStar::dataCheck()
{
  clearErrorCode();
  clearWarningCode();
  initialize();
  int32_t err = 0;

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
    QString ss = QObject::tr("Error opening input file: %1").arg(getInputHeaderFile());
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

  QFileInfo headerPath(getInputHeaderFile());
  QDir headerDir = headerPath.absoluteDir();

  for(size_t iter = 0; iter < m_InputFiles.size(); iter++)
  {
    QFileInfo fi(headerDir.absolutePath() + "/" + m_InputFiles[iter]);
    m_InputFiles[iter] = headerDir.absolutePath() + "/" + m_InputFiles[iter];

    if(fi.exists() == false)
    {
      QString ss = QObject::tr("The input binary CT file: %1 could not be found in the path: %2. The input binary CT file name was obtained from the provided header.  Please ensure the file is in "
                               "the same path as the header and the name has not been altered.")
                       .arg(m_InputFiles[iter])
                       .arg(headerDir.absolutePath());
      setErrorCondition(-388, ss);
    }
  }

  DataContainer::Pointer m = getDataContainerArray()->createNonPrereqDataContainer<AbstractFilter>(this, getDataContainerName());

  if(getErrorCode() < 0)
  {
    return;
  }

  m->setGeometry(image);

  SizeVec3Type dims = image->getDimensions();

  QVector<size_t> tDims = {dims[0], dims[1], dims[2]};
  QVector<size_t> cDims(1, 1);

  m->createNonPrereqAttributeMatrix(this, getCellAttributeMatrixName(), tDims, AttributeMatrix::Type::Cell);

  DataArrayPath path(getDataContainerName(), getCellAttributeMatrixName(), getDensityArrayName());

  m_DensityPtr =
      getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<float>, AbstractFilter, float>(this, path, 0, cDims); /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
  if(nullptr != m_DensityPtr.lock().get()) /* Validate the Weak Pointer wraps a non-nullptr pointer to a DataArray<T> object */
  {
    m_Density = m_DensityPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ReadBinaryCTNorthStar::preflight()
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
int32_t ReadBinaryCTNorthStar::sanityCheckFileSizeVersusAllocatedSize(size_t allocatedBytes, size_t fileSize)
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
int32_t ReadBinaryCTNorthStar::readBinaryCTFiles(size_t dims[3])
{
  int32_t error = 0;
  size_t zShift = 0;

  for(size_t iter = 0; iter < m_InputFiles.size(); iter++)
  {
    QFileInfo fi(m_InputFiles[iter]);
    size_t filesize = static_cast<size_t>(fi.size());
    // allocated bytes should be the x * y dims * number of slices in the current data file....not necessarily the size of the whole density array
    size_t allocatedBytes = dims[0] * dims[1] * m_SlicesPerFile[iter] * sizeof(float);

    error = sanityCheckFileSizeVersusAllocatedSize(allocatedBytes, filesize);

    if(error < 0)
    {
      QString ss = QObject::tr("Binary file size is smaller than the number of allocated bytes");
      setErrorCondition(-100, ss);
      return getErrorCode();
    }

    FILE* f = fopen(m_InputFiles[iter].toLatin1().data(), "rb");
    if(nullptr == f)
    {
      QString ss = QObject::tr("Error opening binary input file: %1").arg(m_InputFiles[iter]);
      setErrorCondition(-100, ss);
      return getErrorCode();
    }

    ScopedFileMonitor monitor(f);

    float* ptr = m_DensityPtr.lock()->getPointer(0);

    size_t index = 0;
    size_t totalPoints = dims[0] * dims[1] * m_SlicesPerFile[iter];
    int64_t progIncrement = static_cast<int64_t>(totalPoints / 100);
    int64_t prog = 1;
    int64_t progressInt = 0;
    int64_t counter = 0;

    for(size_t z = zShift; z < (zShift + m_SlicesPerFile[iter]); z++)
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
            QString ss = QObject::tr("Reading Data || Data File: %1 of %2 || %3% Completed").arg(iter + 1).arg(m_InputFiles.size()).arg(progressInt);
            notifyStatusMessage(ss);
            prog = prog + progIncrement;
          }
          counter++;
        }
      }
    }
    zShift += m_SlicesPerFile[iter];
  }

  return error;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int32_t ReadBinaryCTNorthStar::readHeaderMetaData(ImageGeom::Pointer image)
{
  int32_t error = 0;

  m_InputFiles.clear();
  m_SlicesPerFile.clear();

  QByteArray buf;
  QList<QByteArray> tokens;
  bool ok = false;

  size_t dims[3] = {0, 0, 0};
  float maxRes[3] = {0.0, 0.0, 0.0};
  float minRes[3] = {0.0, 0.0, 0.0};
  float res[3] = {0.0, 0.0, 0.0};
  bool done = false;

  while(m_InHeaderStream.atEnd() == false && done == false)
  {
    buf = m_InHeaderStream.readLine();
    buf = buf.trimmed();
    tokens = buf.split(' ');
    for(qint32 i = 0; i < tokens.size(); i++)
    {
      if(tokens[i] == "<Voxels>")
      {
        dims[1] = tokens[i + 1].toULongLong(&ok);
        if(!ok)
        {
          return -1;
        }
        dims[2] = tokens[i + 2].toULongLong(&ok);
        if(!ok)
        {
          return -1;
        }
        dims[0] = tokens[i + 3].toULongLong(&ok);
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

  while(m_InHeaderStream.atEnd() == false && done == false)
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
            minRes[1] = tokens[j + 1].toFloat(&ok);
            if(!ok)
            {
              return -1;
            }
            minRes[2] = tokens[j + 2].toFloat(&ok);
            if(!ok)
            {
              return -1;
            }
            minRes[0] = tokens[j + 3].toFloat(&ok);
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
            maxRes[1] = tokens[j + 1].toFloat(&ok);
            if(!ok)
            {
              return -1;
            }
            maxRes[2] = tokens[j + 2].toFloat(&ok);
            if(!ok)
            {
              return -1;
            }
            maxRes[0] = tokens[j + 3].toFloat(&ok);
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

  while(m_InHeaderStream.atEnd() == false && done == false)
  {
    buf = m_InHeaderStream.readLine();
    buf = buf.trimmed();
    tokens = buf.split(' ');
    for(qint32 i = 0; i < tokens.size(); i++)
    {
      if(tokens[i] == "<Files>")
      {
        while(done == false)
        {
          buf = m_InHeaderStream.readLine();
          if(m_InHeaderStream.atEnd())
          {
            done = true;
            break;
          }
          buf = buf.trimmed();
          // splitting by the right chevron instead of space to keep spaces in file names and account for when file name is butted up against the chevron with no space
          tokens = buf.split('>');
          for(qint32 j = 0; j < tokens.size(); j++)
          {
            if(tokens[j] == "</Files")
            {
              done = true;
              break;
            }
            if(tokens[j] == "<Name")
            {
              m_InputFiles.push_back(tokens[j + 1]);
              break;
            }
            if(tokens[j] == "<NbSlices")
            {
              m_SlicesPerFile.push_back(tokens[j + 1].trimmed().toULongLong(&ok));
              if(!ok)
              {
                return -1;
              }
              break;
            }
          }
        }
        break;
      }
    }
  }

  for(size_t i = 0; i < 3; i++)
  {
    res[i] = (maxRes[i] - minRes[i]) / dims[i];
  }

  image->setDimensions(dims);
  image->setSpacing(res);

  return error;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ReadBinaryCTNorthStar::execute()
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

  err = readBinaryCTFiles(dims.data());

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
const QString ReadBinaryCTNorthStar::getCompiledLibraryName() const
{
  return DREAM3DReviewConstants::DREAM3DReviewBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString ReadBinaryCTNorthStar::getBrandingString() const
{
  return "DREAM3DReview";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString ReadBinaryCTNorthStar::getFilterVersion() const
{
  QString version;
  QTextStream vStream(&version);
  vStream << DREAM3DReview::Version::Major() << "." << DREAM3DReview::Version::Minor() << "." << DREAM3DReview::Version::Patch();
  return version;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString ReadBinaryCTNorthStar::getGroupName() const
{
  return SIMPL::FilterGroups::IOFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString ReadBinaryCTNorthStar::getSubGroupName() const
{
  return SIMPL::FilterSubGroups::InputFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString ReadBinaryCTNorthStar::getHumanLabel() const
{
  return "Import Binary CT File (North Star Imaging)";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QUuid ReadBinaryCTNorthStar::getUuid()
{
  return QUuid("{f2259481-5011-5f22-9fcb-c92fb6f8be10}");
}
