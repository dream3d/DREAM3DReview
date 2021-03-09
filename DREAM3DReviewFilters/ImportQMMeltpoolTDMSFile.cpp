#include "ImportQMMeltpoolTDMSFile.h"

#include <QtCore/QFileInfo>

#include "SIMPLib/Common/Constants.h"
#include "SIMPLib/DataContainers/DataContainer.h"
#include "SIMPLib/DataContainers/DataContainerArray.h"
#include "SIMPLib/FilterParameters/DataContainerCreationFilterParameter.h"
#include "SIMPLib/FilterParameters/InputFileFilterParameter.h"

#include "DREAM3DReview/DREAM3DReviewConstants.h"
#include "DREAM3DReview/DREAM3DReviewVersion.h"

#include "DREAM3DReview/TDMSSupport/TDMSExceptionHandler.h"
#include "DREAM3DReview/TDMSSupport/TDMSFileProxy.h"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
ImportQMMeltpoolTDMSFile::ImportQMMeltpoolTDMSFile()
{
  initialize();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
ImportQMMeltpoolTDMSFile::~ImportQMMeltpoolTDMSFile() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ImportQMMeltpoolTDMSFile::initialize()
{
  clearErrorCode();
  clearWarningCode();
  setCancel(false);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ImportQMMeltpoolTDMSFile::setupFilterParameters()
{
  FilterParameterVectorType parameters;
  parameters.push_back(SIMPL_NEW_INPUT_FILE_FP("Input TDMS File", InputFile, FilterParameter::Category::Parameter, ImportQMMeltpoolTDMSFile, "*.raw *.bin"));
  parameters.push_back(SIMPL_NEW_DC_CREATION_FP("Data Container", DataContainerName, FilterParameter::Category::CreatedArray, ImportQMMeltpoolTDMSFile));
  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ImportQMMeltpoolTDMSFile::dataCheck()
{
  clearErrorCode();
  clearWarningCode();

  QFileInfo fi(getInputFile());

  if(getInputFile().isEmpty())
  {
    QString ss = QObject::tr("The input binary CT file must be set");
    setErrorCondition(-387, ss);
  }
  else if(!fi.exists())
  {
    QString ss = QObject::tr("The input binary CT file does not exist");
    setErrorCondition(-388, ss);
  }

  TDMSFileProxy::Pointer proxy = nullptr;
  try
  {
    proxy = TDMSFileProxy::New(getInputFile().toStdString());
    proxy->readMetaData();
    // proxy->allocateObjects();
    // proxy->readRawData();
  } catch(const FatalTDMSException& exc)
  {
    std::string tmp = exc.getMessage();
    QString msg = QString::fromStdString(exc.getMessage());
    setErrorCondition(-1, msg);
    return;
  }

  std::unordered_map<std::string, TDMSObject::Pointer> channels = proxy->channelObjects();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ImportQMMeltpoolTDMSFile::execute()
{
  initialize();
  dataCheck();
  if(getErrorCode() < 0)
  {
    return;
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AbstractFilter::Pointer ImportQMMeltpoolTDMSFile::newFilterInstance(bool copyFilterParameters) const
{
  ImportQMMeltpoolTDMSFile::Pointer filter = ImportQMMeltpoolTDMSFile::New();
  if(copyFilterParameters)
  {
    copyFilterParameterInstanceVariables(filter.get());
  }
  return filter;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ImportQMMeltpoolTDMSFile::getCompiledLibraryName() const
{
  return DREAM3DReviewConstants::DREAM3DReviewBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ImportQMMeltpoolTDMSFile::getBrandingString() const
{
  return "DREAM3DReview";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ImportQMMeltpoolTDMSFile::getFilterVersion() const
{
  QString version;
  QTextStream vStream(&version);
  vStream << DREAM3DReview::Version::Major() << "." << DREAM3DReview::Version::Minor() << "." << DREAM3DReview::Version::Patch();
  return version;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ImportQMMeltpoolTDMSFile::getGroupName() const
{
  return SIMPL::FilterGroups::IOFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ImportQMMeltpoolTDMSFile::getSubGroupName() const
{
  return SIMPL::FilterSubGroups::InputFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ImportQMMeltpoolTDMSFile::getHumanLabel() const
{
  return "Import QM Meltpool TDMS File";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QUuid ImportQMMeltpoolTDMSFile::getUuid() const
{
  return QUuid("{60b75e1c-4c65-5d86-8cb0-8b8086193d23}");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
ImportQMMeltpoolTDMSFile::Pointer ImportQMMeltpoolTDMSFile::NullPointer()
{
  return Pointer(static_cast<Self*>(nullptr));
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
std::shared_ptr<ImportQMMeltpoolTDMSFile> ImportQMMeltpoolTDMSFile::New()
{
  struct make_shared_enabler : public ImportQMMeltpoolTDMSFile
  {
  };
  std::shared_ptr<make_shared_enabler> val = std::make_shared<make_shared_enabler>();
  val->setupFilterParameters();
  return val;
}

// -----------------------------------------------------------------------------
QString ImportQMMeltpoolTDMSFile::getNameOfClass() const
{
  return QString("ImportQMMeltpoolTDMSFile");
}

// -----------------------------------------------------------------------------
QString ImportQMMeltpoolTDMSFile::ClassName()
{
  return QString("ImportQMMeltpoolTDMSFile");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ImportQMMeltpoolTDMSFile::setInputFile(const QString& value)
{
  m_InputFile = value;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ImportQMMeltpoolTDMSFile::getInputFile() const
{
  return m_InputFile;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ImportQMMeltpoolTDMSFile::setDataContainerName(const DataArrayPath& value)
{
  m_DataContainerName = value;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
DataArrayPath ImportQMMeltpoolTDMSFile::getDataContainerName() const
{
  return m_DataContainerName;
}
