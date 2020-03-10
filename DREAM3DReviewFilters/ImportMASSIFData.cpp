/*
 * Your License or Copyright can go here
 */

#include <memory>

#include "ImportMASSIFData.h"

#include <QtCore/QFileInfo>

#include <QtCore/QTextStream>

#include <QtCore/QDebug>

#include "SIMPLib/Common/Constants.h"

#include "SIMPLib/FilterParameters/AttributeMatrixSelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/InputFileFilterParameter.h"
#include "SIMPLib/FilterParameters/IntFilterParameter.h"
#include "SIMPLib/FilterParameters/StringFilterParameter.h"
#include "SIMPLib/Geometry/ImageGeom.h"
#include "SIMPLib/Utilities/StringOperations.h"
#include "SIMPLib/DataContainers/DataContainerArray.h"
#include "SIMPLib/DataContainers/DataContainer.h"

#include "H5Support/H5ScopedSentinel.h"
#include "H5Support/H5Utilities.h"
#include "H5Support/QH5Utilities.h"
#include "H5Support/QH5Lite.h"

#include "DREAM3DReview/DREAM3DReviewConstants.h"
#include "DREAM3DReview/DREAM3DReviewVersion.h"

/* Create Enumerations to allow the created Attribute Arrays to take part in renaming */
enum createdPathID : RenameDataPath::DataID_t
{
  AttributeMatrixID20 = 20,
  AttributeMatrixID21 = 21,
};

namespace Detail
{
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
template <typename T>
IDataArray::Pointer readH5Dataset(hid_t locId, const QString& datasetPath, const std::vector<size_t>& tDims, const std::vector<size_t>& cDims)
{
  herr_t err = -1;
  IDataArray::Pointer ptr;

  ptr = DataArray<T>::CreateArray(tDims, cDims, datasetPath, true);

  T* data = (T*)(ptr->getVoidPointer(0));
  err = QH5Lite::readPointerDataset(locId, datasetPath, data);
  if(err < 0)
  {
    qDebug() << "readH5Data read error: " << __FILE__ << "(" << __LINE__ << ")";
    ptr = IDataArray::NullPointer();
  }
  return ptr;
}
} // namespace Detail

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
ImportMASSIFData::ImportMASSIFData()
: m_FilePrefix("Step-")
, m_StepNumber(2)
{
  initialize();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
ImportMASSIFData::~ImportMASSIFData() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ImportMASSIFData::initialize()
{
  clearErrorCode();
  clearWarningCode();
  setCancel(false);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ImportMASSIFData::setupFilterParameters()
{
  FilterParameterVectorType parameters;

  {
    parameters.push_back(SIMPL_NEW_INPUT_FILE_FP("Input File", MassifInputFilePath, FilterParameter::Parameter, ImportMASSIFData, "*.h5;*.hdf5;*.dream3d"));
  }

  {
    parameters.push_back(SIMPL_NEW_STRING_FP("File Prefix", FilePrefix, FilterParameter::Parameter, ImportMASSIFData));
  }

  {
    parameters.push_back(SIMPL_NEW_INTEGER_FP("Step Value", StepNumber, FilterParameter::Parameter, ImportMASSIFData));
  }

  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ImportMASSIFData::dataCheck()
{
  clearErrorCode();
  clearWarningCode();

  if(m_MassifInputFilePath.isEmpty())
  {
    QString ss = "No HDF5 file has been selected";
    setErrorCondition(-3001, ss);
    return;
  }

  QFileInfo hdf5FileInfo(m_MassifInputFilePath);
  QString ext = hdf5FileInfo.suffix();
  if(hdf5FileInfo.exists() && hdf5FileInfo.isFile() && ext != "h5" && ext != "hdf5" && ext != "dream3d")
  {
    QString ss = tr("The selected file '%1' is not an HDF5 file.").arg(hdf5FileInfo.fileName());
    setErrorCondition(-3002, ss);
    return;
  }

  if(m_StepNumber < 0 || m_StepNumber > MASSIFUtilitiesConstants::ImportMassifData::MaxStepNumber)
  {
    QString ss = tr("The selected step number is greater than %1, the maximum allowable step number.").arg(MASSIFUtilitiesConstants::ImportMassifData::MaxStepNumber);
    setErrorCondition(-3003, ss);
    return;
  }

  if(m_FilePrefix.isEmpty())
  {
    QString ss = tr("The file prefix is empty.  Are you sure you meant to do this?");
    setWarningCondition(-3004, ss);
  }

  // Get the geometry info of the arrays
  std::vector<size_t> geoDims;
  FloatVec3Type origin;
  FloatVec3Type res;
  getDataContainerGeometry(geoDims, origin, res);
  if(getErrorCode() < 0)
  {
    return;
  }
  if(geoDims.size() != 3)
  {
    QString ss = tr("The dimensions are not of size 3");
    setErrorCondition(-3003, ss);
    return;
  }
  if(origin.size() != 3)
  {
    QString ss = tr("The origin is not of size 3");
    setErrorCondition(-3003, ss);
    return;
  }
  if(res.size() != 3)
  {
    QString ss = tr("The resolution does not have a size of 3");
    setErrorCondition(-3003, ss);
    return;
  }

  DataContainer::Pointer dc = getDataContainerArray()->createNonPrereqDataContainer(this, MASSIFUtilitiesConstants::ImportMassifData::MassifDC);
  ImageGeom::Pointer image = ImageGeom::CreateGeometry(SIMPL::Geometry::ImageGeometry);
  image->setSpacing(FloatVec3Type(res[0], res[1], res[2]));
  image->setOrigin(FloatVec3Type(origin[0], origin[1], origin[2]));
  image->setDimensions(SizeVec3Type(geoDims[0], geoDims[1], geoDims[2]));
  dc->setGeometry(image);
  if(getErrorCode() < 0)
  {
    return;
  }

  int err = 0;
  AttributeMatrix::Pointer am = dc->createNonPrereqAttributeMatrix(this, MASSIFUtilitiesConstants::ImportMassifData::MassifAM, geoDims, AttributeMatrix::Type::Cell, AttributeMatrixID21);
  if(getErrorCode() < 0 || err < 0)
  {
    return;
  }

  hid_t fileId = H5Utilities::openFile(m_MassifInputFilePath.toStdString(), true);
  if(fileId < 0)
  {
    QString ss = tr("Error Reading HDF5 file: %1\n").arg(m_MassifInputFilePath);
    setErrorCondition(-3003, ss);
    return;
  }

  QVector<QString> hdf5ArrayPaths = createHDF5DatasetPaths();

  for(const auto& hdf5ArrayPath : hdf5ArrayPaths)
  {
    QString parentPath = QH5Utilities::getParentPath(hdf5ArrayPath);

    hid_t parentId = QH5Utilities::openHDF5Object(fileId, parentPath);
    H5ScopedGroupSentinel sentinel(&parentId, false);
    // Read dataset into DREAM.3D structure
    QString objectName = QH5Utilities::getObjectNameFromPath(hdf5ArrayPath);
    std::vector<size_t> geometryDims = {geoDims[0], geoDims[1], geoDims[2]};
    IDataArray::Pointer dPtr = readIDataArray(parentId, objectName, geometryDims, getInPreflight());
    if(dPtr == IDataArray::NullPointer())
    {
      QString ss = tr("Could not read dataset '%1' at path '%2'").arg(objectName).arg(parentPath);
      setErrorCondition(-3003, ss);
      return;
    }
    am->insertOrAssign(dPtr);
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QVector<DataArrayPath> ImportMASSIFData::createDataArrayPaths()
{
  QVector<DataArrayPath> dataArrayPaths;
  dataArrayPaths.push_back(
      DataArrayPath(MASSIFUtilitiesConstants::ImportMassifData::MassifDC, MASSIFUtilitiesConstants::ImportMassifData::MassifAM, MASSIFUtilitiesConstants::ImportMassifData::DField));
  dataArrayPaths.push_back(
      DataArrayPath(MASSIFUtilitiesConstants::ImportMassifData::MassifDC, MASSIFUtilitiesConstants::ImportMassifData::MassifAM, MASSIFUtilitiesConstants::ImportMassifData::EField));
  dataArrayPaths.push_back(
      DataArrayPath(MASSIFUtilitiesConstants::ImportMassifData::MassifDC, MASSIFUtilitiesConstants::ImportMassifData::MassifAM, MASSIFUtilitiesConstants::ImportMassifData::SField));

  dataArrayPaths.push_back(DataArrayPath(MASSIFUtilitiesConstants::ImportMassifData::MassifDC, MASSIFUtilitiesConstants::ImportMassifData::MassifAM, MASSIFUtilitiesConstants::ImportMassifData::EVM));
  dataArrayPaths.push_back(DataArrayPath(MASSIFUtilitiesConstants::ImportMassifData::MassifDC, MASSIFUtilitiesConstants::ImportMassifData::MassifAM, MASSIFUtilitiesConstants::ImportMassifData::SVM));
  dataArrayPaths.push_back(
      DataArrayPath(MASSIFUtilitiesConstants::ImportMassifData::MassifDC, MASSIFUtilitiesConstants::ImportMassifData::MassifAM, MASSIFUtilitiesConstants::ImportMassifData::GrainID));
  dataArrayPaths.push_back(
      DataArrayPath(MASSIFUtilitiesConstants::ImportMassifData::MassifDC, MASSIFUtilitiesConstants::ImportMassifData::MassifAM, MASSIFUtilitiesConstants::ImportMassifData::Phase));

  return dataArrayPaths;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QVector<QString> ImportMASSIFData::createHDF5DatasetPaths()
{
  QVector<QString> arrayPaths;
  QString parentPath = "/" + MASSIFUtilitiesConstants::ImportMassifData::DCGrpName + "/" + m_PaddedStep + "/" + MASSIFUtilitiesConstants::ImportMassifData::Datapoint;

  arrayPaths.push_back(parentPath + "/" + MASSIFUtilitiesConstants::ImportMassifData::DFieldsGrpName + "/" + MASSIFUtilitiesConstants::ImportMassifData::DField);

  arrayPaths.push_back(parentPath + "/" + MASSIFUtilitiesConstants::ImportMassifData::EFieldsGrpName + "/" + MASSIFUtilitiesConstants::ImportMassifData::EField);

  arrayPaths.push_back(parentPath + "/" + MASSIFUtilitiesConstants::ImportMassifData::SFieldsGrpName + "/" + MASSIFUtilitiesConstants::ImportMassifData::SField);

  arrayPaths.push_back(parentPath + "/" + MASSIFUtilitiesConstants::ImportMassifData::DFieldsGrpName + "/" + MASSIFUtilitiesConstants::ImportMassifData::EVM);

  arrayPaths.push_back(parentPath + "/" + MASSIFUtilitiesConstants::ImportMassifData::SFieldsGrpName + "/" + MASSIFUtilitiesConstants::ImportMassifData::SVM);

  arrayPaths.push_back(parentPath + "/" + MASSIFUtilitiesConstants::ImportMassifData::EulerAngleGrpName + "/" + MASSIFUtilitiesConstants::ImportMassifData::Phi1);

  arrayPaths.push_back(parentPath + "/" + MASSIFUtilitiesConstants::ImportMassifData::EulerAngleGrpName + "/" + MASSIFUtilitiesConstants::ImportMassifData::Phi);

  arrayPaths.push_back(parentPath + "/" + MASSIFUtilitiesConstants::ImportMassifData::EulerAngleGrpName + "/" + MASSIFUtilitiesConstants::ImportMassifData::Phi2);

  arrayPaths.push_back(parentPath + "/" + MASSIFUtilitiesConstants::ImportMassifData::GrainID);

  arrayPaths.push_back(parentPath + "/" + MASSIFUtilitiesConstants::ImportMassifData::Phase);

  return arrayPaths;
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ImportMASSIFData::execute()
{
  initialize();
  dataCheck();
  if(getErrorCode() < 0)
  {
    return;
  }

  if(getCancel())
  {
    return;
  }

  if(getErrorCode() < 0)
  {
    QString ss = QObject::tr("Some error message");
    setErrorCondition(-99999999, ss);
    return;
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ImportMASSIFData::getDataContainerGeometry(std::vector<size_t>& tDims, FloatVec3Type& origin, FloatVec3Type& spacing)
{
  hid_t fileId = QH5Utilities::openFile(m_MassifInputFilePath, true);
  if(fileId < 0)
  {
    QString ss = tr("Error Reading HDF5 file: %1\n").arg(m_MassifInputFilePath);
    setErrorCondition(-3005, ss);
    return;
  }

  H5ScopedFileSentinel sentinel(&fileId, true);

  QString totalPath = MASSIFUtilitiesConstants::ImportMassifData::DCGrpName;
  hid_t gid = QH5Utilities::openHDF5Object(fileId, MASSIFUtilitiesConstants::ImportMassifData::DCGrpName);
  if(gid < 0)
  {
    QString ss = tr("Could not open path: %1\n").arg(totalPath);
    setErrorCondition(-3006, ss);
    return;
  }
  sentinel.addGroupId(&gid);

  // Build up the step
  m_PaddedStep = m_FilePrefix;
  m_PaddedStep.append(StringOperations::GenerateIndexString(m_StepNumber, MASSIFUtilitiesConstants::ImportMassifData::MaxStepNumber));
  totalPath.append("/" + m_PaddedStep);

  hid_t stepGid = QH5Utilities::openHDF5Object(gid, m_PaddedStep);
  if(stepGid < 0)
  {
    QString ss = tr("Could not open path: %1\n").arg(totalPath);
    setErrorCondition(-3007, ss);
    return;
  }
  sentinel.addGroupId(&stepGid);

  totalPath.append("/" + MASSIFUtilitiesConstants::ImportMassifData::GeometryGrpName);
  hid_t geoGid = QH5Utilities::openHDF5Object(stepGid, MASSIFUtilitiesConstants::ImportMassifData::GeometryGrpName);
  if(geoGid < 0)
  {
    QString ss = tr("Could not open path: %1\n").arg(totalPath);
    setErrorCondition(-3007, ss);
    return;
  }
  sentinel.addGroupId(&geoGid);

  herr_t err = QH5Lite::readPointerDataset(geoGid, MASSIFUtilitiesConstants::ImportMassifData::DimGrpName, tDims.data());
  if(err < 0)
  {
    QString ss = tr("Could not read dataset %1 at path '%2'\n").arg(MASSIFUtilitiesConstants::ImportMassifData::DimGrpName).arg(totalPath);
    setErrorCondition(-3008, ss);
    return;
  }

  err = QH5Lite::readPointerDataset(geoGid, MASSIFUtilitiesConstants::ImportMassifData::OriginGrpName, origin.data());
  if(err < 0)
  {
    QString ss = tr("Could not read dataset %1 at path '%2'\n").arg(MASSIFUtilitiesConstants::ImportMassifData::OriginGrpName).arg(totalPath);
    setErrorCondition(-3009, ss);
    return;
  }

  err = QH5Lite::readPointerDataset(geoGid, MASSIFUtilitiesConstants::ImportMassifData::SpacingGrpName, spacing.data());
  if(err < 0)
  {
    QString ss = tr("Could not read dataset %1 at path '%2'\n").arg(MASSIFUtilitiesConstants::ImportMassifData::SpacingGrpName).arg(totalPath);
    setErrorCondition(-3010, ss);
    return;
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AbstractFilter::Pointer ImportMASSIFData::newFilterInstance(bool copyFilterParameters) const
{
  ImportMASSIFData::Pointer filter = ImportMASSIFData::New();
  if(copyFilterParameters)
  {
    copyFilterParameterInstanceVariables(filter.get());
  }
  return filter;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ImportMASSIFData::getCompiledLibraryName() const
{
  return MASSIFUtilitiesConstants::MASSIFUtilitiesBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ImportMASSIFData::getBrandingString() const
{
  return "MASSIFUtilities Plugin";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ImportMASSIFData::getFilterVersion() const
{
  QString version;
  QTextStream vStream(&version);
  vStream << DREAM3DReview::Version::Major() << "." << DREAM3DReview::Version::Minor() << "." << DREAM3DReview::Version::Patch();
  return version;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ImportMASSIFData::getGroupName() const
{
  return SIMPL::FilterGroups::IOFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QUuid ImportMASSIFData::getUuid() const
{
  return QUuid("{85654f78-04a8-50ed-9ae1-25dfbd0539b3}");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ImportMASSIFData::getSubGroupName() const
{
  return SIMPL::FilterSubGroups::InputFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ImportMASSIFData::getHumanLabel() const
{
  return "Import MASSIF Data (HDF5)";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
IDataArray::Pointer ImportMASSIFData::readIDataArray(hid_t gid, const QString& name, std::vector<size_t> geoDims, bool metaDataOnly)
{

  herr_t err = -1;
  // herr_t retErr = 1;
  hid_t typeId = -1;
  H5T_class_t attr_type;
  size_t attr_size;
  QString res;

  QVector<hsize_t> dims; // Reusable for the loop
  IDataArray::Pointer ptr = IDataArray::NullPointer();
  // qDebug() << "Reading Attribute " << *iter ;
  typeId = QH5Lite::getDatasetType(gid, name);
  if(typeId < 0)
  {
    return ptr;
  }
  // Get the HDF5 DataSet information. the dimensions will be the combined Tuple Dims and the Data Array Componenet dimes
  err = QH5Lite::getDatasetInfo(gid, name, dims, attr_type, attr_size);
  if(err < 0)
  {
    qDebug() << "Error in getAttributeInfo method in readUserMetaData.";
  }
  else
  {
    QString classType;
    // int version = 0;
    std::vector<size_t>& tDims = geoDims;
    std::vector<size_t> cDims;

    if(geoDims.size() == dims.size())
    {
      cDims.push_back(1);
    }
    else
    {
      for(int i = geoDims.size(); i < dims.size(); i++)
      {
        cDims.push_back(dims[i]);
      }
    }

    //    err = ReadRequiredAttributes(gid, name, classType, version, tDims, cDims);
    //    if(err < 0)
    //    {
    //      return ptr;
    //    }

    // Check to see if we are reading a bool array and if so read it and return
    if(classType.compare("DataArray<bool>") == 0)
    {
      if(!metaDataOnly)
      {
        ptr = Detail::readH5Dataset<bool>(gid, name, tDims, cDims);
      }
      else
      {
        ptr = DataArray<bool>::CreateArray(tDims, cDims, name, false);
      }
      err = H5Tclose(typeId);
      return ptr; // <== Note early return here.
    }
    switch(attr_type)
    {
    case H5T_STRING:
      res.clear(); // Clear the string out first
      err = QH5Lite::readStringDataset(gid, name, res);
      break;
    case H5T_INTEGER:
      // qDebug() << "User Meta Data Type is Integer" ;
      if((H5Tequal(typeId, H5T_STD_U8BE) != 0) || (H5Tequal(typeId, H5T_STD_U8LE) != 0))
      {
        if(!metaDataOnly)
        {
          ptr = Detail::readH5Dataset<uint8_t>(gid, name, tDims, cDims);
        }
        else
        {
          ptr = DataArray<uint8_t>::CreateArray(tDims, cDims, name, false);
        }
      }
      else if((H5Tequal(typeId, H5T_STD_U16BE) != 0) || (H5Tequal(typeId, H5T_STD_U16LE) != 0))
      {
        if(!metaDataOnly)
        {
          ptr = Detail::readH5Dataset<uint16_t>(gid, name, tDims, cDims);
        }
        else
        {
          ptr = DataArray<uint16_t>::CreateArray(tDims, cDims, name, false);
        }
      }
      else if((H5Tequal(typeId, H5T_STD_U32BE) != 0) || (H5Tequal(typeId, H5T_STD_U32LE) != 0))
      {
        if(!metaDataOnly)
        {
          ptr = Detail::readH5Dataset<uint32_t>(gid, name, tDims, cDims);
        }
        else
        {
          ptr = DataArray<uint32_t>::CreateArray(tDims, cDims, name, false);
        }
      }
      else if((H5Tequal(typeId, H5T_STD_U64BE) != 0) || (H5Tequal(typeId, H5T_STD_U64LE) != 0))
      {
        if(!metaDataOnly)
        {
          ptr = Detail::readH5Dataset<uint64_t>(gid, name, tDims, cDims);
        }
        else
        {
          ptr = DataArray<uint64_t>::CreateArray(tDims, cDims, name, false);
        }
      }
      else if((H5Tequal(typeId, H5T_STD_I8BE) != 0) || (H5Tequal(typeId, H5T_STD_I8LE) != 0))
      {
        if(!metaDataOnly)
        {
          ptr = Detail::readH5Dataset<int8_t>(gid, name, tDims, cDims);
        }
        else
        {
          ptr = DataArray<int8_t>::CreateArray(tDims, cDims, name, false);
        }
      }
      else if((H5Tequal(typeId, H5T_STD_I16BE) != 0) || (H5Tequal(typeId, H5T_STD_I16LE) != 0))
      {
        if(!metaDataOnly)
        {
          ptr = Detail::readH5Dataset<int16_t>(gid, name, tDims, cDims);
        }
        else
        {
          ptr = DataArray<int16_t>::CreateArray(tDims, cDims, name, false);
        }
      }
      else if((H5Tequal(typeId, H5T_STD_I32BE) != 0) || (H5Tequal(typeId, H5T_STD_I32LE) != 0))
      {
        if(!metaDataOnly)
        {
          ptr = Detail::readH5Dataset<int32_t>(gid, name, tDims, cDims);
        }
        else
        {
          ptr = DataArray<int32_t>::CreateArray(tDims, cDims, name, false);
        }
      }
      else if((H5Tequal(typeId, H5T_STD_I64BE) != 0) || (H5Tequal(typeId, H5T_STD_I64LE) != 0))
      {
        if(!metaDataOnly)
        {
          ptr = Detail::readH5Dataset<int64_t>(gid, name, tDims, cDims);
        }
        else
        {
          ptr = DataArray<int64_t>::CreateArray(tDims, cDims, name, false);
        }
      }
      else
      {
        qDebug() << "Unknown Type: " << typeId << " at " << name;
        err = -1;
      }
      break;
    case H5T_FLOAT:
      if(attr_size == 4)
      {
        if(!metaDataOnly)
        {
          ptr = Detail::readH5Dataset<float>(gid, name, tDims, cDims);
        }
        else
        {
          ptr = DataArray<float>::CreateArray(tDims, cDims, name, false);
        }
      }
      else if(attr_size == 8)
      {
        if(!metaDataOnly)
        {
          ptr = Detail::readH5Dataset<double>(gid, name, tDims, cDims);
        }
        else
        {
          ptr = DataArray<double>::CreateArray(tDims, cDims, name, false);
        }
      }
      else
      {
        qDebug() << "Unknown Floating point type";
        err = -1;
      }
      break;
    default:
      qDebug() << "Error: readUserMetaData() Unknown attribute type: " << attr_type;
      QH5Utilities::printHDFClassType(attr_type);
    }

    err = H5Tclose(typeId);
    // Close the H5A type Id that was retrieved during the loop
  }

  return ptr;
}

// -----------------------------------------------------------------------------
ImportMASSIFData::Pointer ImportMASSIFData::NullPointer()
{
  return Pointer(static_cast<Self*>(nullptr));
}

// -----------------------------------------------------------------------------
std::shared_ptr<ImportMASSIFData> ImportMASSIFData::New()
{
  struct make_shared_enabler : public ImportMASSIFData
  {
  };
  std::shared_ptr<make_shared_enabler> val = std::make_shared<make_shared_enabler>();
  val->setupFilterParameters();
  return val;
}

// -----------------------------------------------------------------------------
QString ImportMASSIFData::getNameOfClass() const
{
  return QString("ImportMASSIFData");
}

// -----------------------------------------------------------------------------
QString ImportMASSIFData::ClassName()
{
  return QString("ImportMASSIFData");
}

// -----------------------------------------------------------------------------
void ImportMASSIFData::setMassifInputFilePath(const QString& value)
{
  m_MassifInputFilePath = value;
}

// -----------------------------------------------------------------------------
QString ImportMASSIFData::getMassifInputFilePath() const
{
  return m_MassifInputFilePath;
}

// -----------------------------------------------------------------------------
void ImportMASSIFData::setFilePrefix(const QString& value)
{
  m_FilePrefix = value;
}

// -----------------------------------------------------------------------------
QString ImportMASSIFData::getFilePrefix() const
{
  return m_FilePrefix;
}

// -----------------------------------------------------------------------------
void ImportMASSIFData::setStepNumber(int value)
{
  m_StepNumber = value;
}

// -----------------------------------------------------------------------------
int ImportMASSIFData::getStepNumber() const
{
  return m_StepNumber;
}
