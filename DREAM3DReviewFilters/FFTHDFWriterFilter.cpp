/*
 * Your License or Copyright can go here
 */

#include "FFTHDFWriterFilter.h"

#include <QtCore/QDir>
#include <QtCore/QTextStream>

#include "SIMPLib/Common/Constants.h"
#include "SIMPLib/DataContainers/DataContainer.h"
#include "SIMPLib/DataContainers/DataContainerArray.h"
#include "SIMPLib/FilterParameters/AbstractFilterParametersReader.h"
#include "SIMPLib/FilterParameters/BooleanFilterParameter.h"
#include "SIMPLib/FilterParameters/DataArraySelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/OutputFileFilterParameter.h"
#include "SIMPLib/FilterParameters/SeparatorFilterParameter.h"
#include "SIMPLib/Geometry/ImageGeom.h"
#include "SIMPLib/Utilities/FileSystemPathHelper.h"

#include "H5Support/H5ScopedSentinel.h"
#include "H5Support/QH5Utilities.h"

#include "DREAM3DReview/DREAM3DReviewConstants.h"
#include "DREAM3DReview/DREAM3DReviewVersion.h"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
FFTHDFWriterFilter::FFTHDFWriterFilter() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
FFTHDFWriterFilter::~FFTHDFWriterFilter()
{
  closeFile(m_FileId);
  closeFile(m_FileIdEig);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FFTHDFWriterFilter::setupFilterParameters()
{
  FilterParameterVectorType parameters;

  parameters.push_back(OutputFileFilterParameter::Create("Output File", "OutputFile", getOutputFile(), FilterParameter::Category::Parameter, SIMPL_BIND_SETTER(FFTHDFWriterFilter, this, OutputFile),
                                                         SIMPL_BIND_GETTER(FFTHDFWriterFilter, this, OutputFile), "*.dream3d", ""));
  //  parameters.push_back(BooleanFilterParameter::New("Write Xdmf File", "WriteXdmfFile", getWriteXdmfFile(), FilterParameter::Parameter, "ParaView Compatible File"));
  QStringList linkedProps;
  linkedProps << "EigenstrainsOutputFile"
              << "CellEigenstrainsArrayPath";
  parameters.push_back(SIMPL_NEW_LINKED_BOOL_FP("Write Eigenstrains", WriteEigenstrains, FilterParameter::Category::Parameter, FFTHDFWriterFilter, linkedProps));

  parameters.push_back(OutputFileFilterParameter::Create("Eigenstrain Output File", "EigenstrainsOutputFile", getEigenstrainsOutputFile(), FilterParameter::Category::Parameter,
                                                         SIMPL_BIND_SETTER(FFTHDFWriterFilter, this, EigenstrainsOutputFile), SIMPL_BIND_GETTER(FFTHDFWriterFilter, this, EigenstrainsOutputFile),
                                                         "*.dream3d", ""));

  //--------------
  parameters.push_back(SeparatorFilterParameter::Create("Cell Data", FilterParameter::Category::RequiredArray));
  {
    DataArraySelectionFilterParameter::RequirementType req = DataArraySelectionFilterParameter::CreateRequirement(SIMPL::TypeNames::Int32, 1, AttributeMatrix::Type::Cell, IGeometry::Type::Image);
    parameters.push_back(DataArraySelectionFilterParameter::Create("Feature Ids", "FeatureIdsArrayPath", getFeatureIdsArrayPath(), FilterParameter::Category::RequiredArray,
                                                                   SIMPL_BIND_SETTER(FFTHDFWriterFilter, this, FeatureIdsArrayPath), SIMPL_BIND_GETTER(FFTHDFWriterFilter, this, FeatureIdsArrayPath),
                                                                   req));
  }
  {
    DataArraySelectionFilterParameter::RequirementType req = DataArraySelectionFilterParameter::CreateRequirement(SIMPL::TypeNames::Float, 3, AttributeMatrix::Type::Cell, IGeometry::Type::Image);
    parameters.push_back(DataArraySelectionFilterParameter::Create("Euler Angles", "CellEulerAnglesArrayPath", getCellEulerAnglesArrayPath(), FilterParameter::Category::RequiredArray,
                                                                   SIMPL_BIND_SETTER(FFTHDFWriterFilter, this, CellEulerAnglesArrayPath),
                                                                   SIMPL_BIND_GETTER(FFTHDFWriterFilter, this, CellEulerAnglesArrayPath), req));
  }
  {
    DataArraySelectionFilterParameter::RequirementType req = DataArraySelectionFilterParameter::CreateRequirement(SIMPL::TypeNames::Int32, 1, AttributeMatrix::Type::Cell, IGeometry::Type::Image);
    parameters.push_back(DataArraySelectionFilterParameter::Create("Phases", "CellPhasesArrayPath", getCellPhasesArrayPath(), FilterParameter::Category::RequiredArray,
                                                                   SIMPL_BIND_SETTER(FFTHDFWriterFilter, this, CellPhasesArrayPath), SIMPL_BIND_GETTER(FFTHDFWriterFilter, this, CellPhasesArrayPath),
                                                                   req));
  }
  {
    DataArraySelectionFilterParameter::RequirementType req = DataArraySelectionFilterParameter::CreateRequirement(SIMPL::TypeNames::Float, 6, AttributeMatrix::Type::Cell, IGeometry::Type::Image);
    parameters.push_back(DataArraySelectionFilterParameter::Create("Eigenstrains", "CellEigenstrainsArrayPath", getCellEigenstrainsArrayPath(), FilterParameter::Category::RequiredArray,
                                                                   SIMPL_BIND_SETTER(FFTHDFWriterFilter, this, CellEigenstrainsArrayPath),
                                                                   SIMPL_BIND_GETTER(FFTHDFWriterFilter, this, CellEigenstrainsArrayPath), req));
  }

  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FFTHDFWriterFilter::readFilterParameters(AbstractFilterParametersReader* reader, int index)
{
  reader->openFilterGroup(this, index);
  setOutputFile(reader->readString("OutputFile", getOutputFile()));
  setWriteEigenstrains(reader->readValue("WriteEigenstrains", getWriteEigenstrains()));
  setEigenstrainsOutputFile(reader->readString("EigenstrainsOutputFile", getEigenstrainsOutputFile()));
  //----------------------------
  setCellEulerAnglesArrayPath(reader->readDataArrayPath("CellEulerAnglesArrayPath", getCellEulerAnglesArrayPath()));
  setCellPhasesArrayPath(reader->readDataArrayPath("CellPhasesArrayPath", getCellPhasesArrayPath()));
  setFeatureIdsArrayPath(reader->readDataArrayPath("FeatureIdsArrayPath", getFeatureIdsArrayPath()));
  setCellEigenstrainsArrayPath(reader->readDataArrayPath("CellEigenstrainsArrayPath", getCellEigenstrainsArrayPath()));

  reader->closeFilterGroup();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FFTHDFWriterFilter::dataCheck()
{
  clearErrorCode();
  clearWarningCode();
  QString ss;

  getDataContainerArray()->getPrereqGeometryFromDataContainer<ImageGeom>(this, getFeatureIdsArrayPath().getDataContainerName());

  QFileInfo fi(getOutputFile());
  if(fi.suffix().compare("") == 0)
  {
    m_OutputFile.append(".dream3d");
  }
  FileSystemPathHelper::CheckOutputFile(this, "Output File Name", getEigenstrainsOutputFile(), true);

  if(m_WriteEigenstrains)
  {
    QFileInfo fiEig(getOutputFile());
    if(fiEig.suffix().compare("") == 0)
    {
      m_EigenstrainsOutputFile.append(".dream3d");
    }
    FileSystemPathHelper::CheckOutputFile(this, "Eigenstrains Output File Name", getEigenstrainsOutputFile(), true);
  }

  QVector<DataArrayPath> dataArrayPaths;

  std::vector<size_t> cDims(1, 1);
  m_FeatureIdsPtr = getDataContainerArray()->getPrereqArrayFromPath<DataArray<int32_t>>(this, getFeatureIdsArrayPath(), cDims);
  if(nullptr != m_FeatureIdsPtr.lock())
  {
    m_FeatureIds = m_FeatureIdsPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */
  if(getErrorCode() >= 0)
  {
    dataArrayPaths.push_back(getFeatureIdsArrayPath());
  }

  m_CellPhasesPtr = getDataContainerArray()->getPrereqArrayFromPath<DataArray<int32_t>>(this, getCellPhasesArrayPath(), cDims);
  if(nullptr != m_CellPhasesPtr.lock())
  {
    m_CellPhases = m_CellPhasesPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */
  if(getErrorCode() >= 0)
  {
    dataArrayPaths.push_back(getCellPhasesArrayPath());
  }

  cDims[0] = 3;
  m_CellEulerAnglesPtr = getDataContainerArray()->getPrereqArrayFromPath<DataArray<float>>(this, getCellEulerAnglesArrayPath(), cDims);
  if(nullptr != m_CellEulerAnglesPtr.lock())
  {
    m_CellEulerAngles = m_CellEulerAnglesPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */
  if(getErrorCode() >= 0)
  {
    dataArrayPaths.push_back(getCellEulerAnglesArrayPath());
  }

  if(m_WriteEigenstrains)
  {
    cDims[0] = 6;
    m_CellEigenstrainsPtr = getDataContainerArray()->getPrereqArrayFromPath<DataArray<float>>(this, getCellEigenstrainsArrayPath(), cDims);
    if(nullptr != m_CellEigenstrainsPtr.lock())
    {
      m_CellEigenstrains = m_CellEigenstrainsPtr.lock()->getPointer(0);
    } /* Now assign the raw pointer to data from the DataArray<T> object */
    if(getErrorCode() >= 0)
    {
      dataArrayPaths.push_back(getCellEigenstrainsArrayPath());
    }
  }

  getDataContainerArray()->validateNumberOfTuples(this, dataArrayPaths);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FFTHDFWriterFilter::execute()
{
  dataCheck();
  if(getErrorCode() < 0)
  {
    return;
  }

  int err = 0;

  // Make sure any directory path is also available as the user may have just typed
  // in a path without actually creating the full path
  QFileInfo fi(m_OutputFile);
  QString parentPath = fi.path();
  QDir dir;
  if(!dir.mkpath(parentPath))
  {
    QString ss = QObject::tr("Error creating parent path '%1'").arg(parentPath);
    setErrorCondition(-11110, ss);
    return;
  }

  openFile(m_OutputFile, m_FileId, m_AppendToExisting); // Do NOT append to any existing file
  if(m_FileId < 0)
  {
    QString ss = QObject::tr("The HDF5 file could not be opened or created.\n The given filename was:\n\t[%1]").arg(m_OutputFile);
    setErrorCondition(-11112, ss);
    return;
  }

  // This will make sure if we return early from this method that the HDF5 File is properly closed.
  H5ScopedFileSentinel scopedFileSentinel(m_FileId, true);

  // Create DataContainer group!
  err = H5Utilities::createGroupsFromPath(SIMPL::StringConstants::DataContainerGroupName.toLatin1().data(), m_FileId);
  if(err < 0)
  {
    QString ss = QObject::tr("Error creating HDF5 Group '%1'").arg(SIMPL::StringConstants::DataContainerGroupName);
    setErrorCondition(-60, ss);
    return;
  }

  hid_t dcaGid = H5Gopen(m_FileId, SIMPL::StringConstants::DataContainerGroupName.toLatin1().data(), H5P_DEFAULT);
  scopedFileSentinel.addGroupId(dcaGid);

  AttributeMatrix::Pointer attrMat = getDataContainerArray()->getAttributeMatrix(m_FeatureIdsArrayPath);
  std::vector<size_t> tDims = attrMat->getTupleDimensions();
  m_FeatureIdsPtr.lock()->writeH5Data(dcaGid, tDims);
  //     H5Lite::writePointerDataset;

  attrMat = getDataContainerArray()->getAttributeMatrix(m_CellPhasesArrayPath);
  tDims = attrMat->getTupleDimensions();
  m_CellPhasesPtr.lock()->writeH5Data(dcaGid, tDims);

  attrMat = getDataContainerArray()->getAttributeMatrix(m_CellEulerAnglesArrayPath);
  tDims = attrMat->getTupleDimensions();
  m_CellEulerAnglesPtr.lock()->writeH5Data(dcaGid, tDims);

  // MASSIF Eigenstrains file
  if(m_WriteEigenstrains)
  {
    QFileInfo fiEig(m_EigenstrainsOutputFile);
    QString parentPathEig = fiEig.path();
    QDir dirEig;
    if(!dirEig.mkpath(parentPathEig))
    {
      QString ss = QObject::tr("Error creating parent path '%1'").arg(parentPathEig);
      setErrorCondition(-11110, ss);
      return;
    }

    openFile(m_EigenstrainsOutputFile, m_FileIdEig, m_AppendToExisting); // Do NOT append to any existing file
    if(m_FileIdEig < 0)
    {
      QString ss = QObject::tr("The HDF5 file could not be opened or created.\n The given filename was:\n\t[%1]").arg(m_EigenstrainsOutputFile);
      setErrorCondition(-11113, ss);
      return;
    }

    H5ScopedFileSentinel scopedFileSentinelEig(m_FileIdEig, true);

    err = H5Utilities::createGroupsFromPath(SIMPL::StringConstants::DataContainerGroupName.toLatin1().data(), m_FileIdEig);
    if(err < 0)
    {
      QString ss = QObject::tr("Error creating HDF5 Group '%1'").arg(SIMPL::StringConstants::DataContainerGroupName);
      setErrorCondition(-60, ss);
      return;
    }

    hid_t dcaGidEig = H5Gopen(m_FileIdEig, SIMPL::StringConstants::DataContainerGroupName.toLatin1().data(), H5P_DEFAULT);
    scopedFileSentinelEig.addGroupId(dcaGidEig);

    attrMat = getDataContainerArray()->getAttributeMatrix(m_FeatureIdsArrayPath);
    tDims = attrMat->getTupleDimensions();
    m_FeatureIdsPtr.lock()->writeH5Data(dcaGidEig, tDims);

    attrMat = getDataContainerArray()->getAttributeMatrix(m_CellEigenstrainsArrayPath);
    tDims = attrMat->getTupleDimensions();
    m_CellEigenstrainsPtr.lock()->writeH5Data(dcaGidEig, tDims);
  }

  // DataContainer::Pointer dc = getDataContainerArray()->getDataContainer(getFeatureIdsArrayPath().getDataContainerName());

  //      err = dc->writeAttributeMatricesToHDF5(dcaGid);
  //    	if (err < 0)
  //        {
  //        setErrorCondition(-803, "Error writing DataContainer AttributeMatrices");
  //     	return;
  //        }

  //        H5Gclose(dcaGid);

  //       dcaGid = -1;
}

//--------------------------------------------------------------

void FFTHDFWriterFilter::writeXdmfHeader(QTextStream& xdmf)
{
  xdmf << "<?xml version=\"1.0\"?>"
       << "\n";
  xdmf << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\"[]>"
       << "\n";
  xdmf << "<Xdmf xmlns:xi=\"http://www.w3.org/2003/XInclude\" Version=\"2.2\">"
       << "\n";
  xdmf << " <Domain>"
       << "\n";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FFTHDFWriterFilter::writeXdmfFooter(QTextStream& xdmf)
{
  xdmf << " </Domain>"
       << "\n";
  xdmf << "</Xdmf>"
       << "\n";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FFTHDFWriterFilter::openFile(QString file, hid_t& fileId, bool appendData)
{
  // Try to open a file to append data into
  if(appendData)
  {
    fileId = QH5Utilities::openFile(file, false);
  }
  // No file was found or we are writing new data only to a clean file
  if(!appendData || fileId < 0)
  {
    fileId = QH5Utilities::createFile(file);
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
herr_t FFTHDFWriterFilter::closeFile(hid_t fileId)
{
  // Close the file when we are finished with it
  if(fileId > 0)
  {
    return QH5Utilities::closeFile(fileId);
  }
  return 0;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AbstractFilter::Pointer FFTHDFWriterFilter::newFilterInstance(bool copyFilterParameters) const
{
  FFTHDFWriterFilter::Pointer filter = FFTHDFWriterFilter::New();
  if(copyFilterParameters)
  {
    copyFilterParameterInstanceVariables(filter.get());
  }
  return filter;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString FFTHDFWriterFilter::getCompiledLibraryName() const
{
  return MASSIFUtilitiesConstants::MASSIFUtilitiesBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString FFTHDFWriterFilter::getBrandingString() const
{
  return "MASSIFUtilities Plugin";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString FFTHDFWriterFilter::getFilterVersion() const
{
  QString version;
  QTextStream vStream(&version);
  vStream << DREAM3DReview::Version::Major() << "." << DREAM3DReview::Version::Minor() << "." << DREAM3DReview::Version::Patch();
  return version;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString FFTHDFWriterFilter::getGroupName() const
{
  return SIMPL::FilterGroups::IOFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QUuid FFTHDFWriterFilter::getUuid() const
{
  return QUuid("{b6b1ba7c-14aa-5c6f-9436-8c46face6053}");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString FFTHDFWriterFilter::getSubGroupName() const
{
  return SIMPL::FilterSubGroups::OutputFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString FFTHDFWriterFilter::getHumanLabel() const
{
  return "Export MASSIF Data (HDF5)";
}

// -----------------------------------------------------------------------------
FFTHDFWriterFilter::Pointer FFTHDFWriterFilter::NullPointer()
{
  return Pointer(static_cast<Self*>(nullptr));
}

// -----------------------------------------------------------------------------
std::shared_ptr<FFTHDFWriterFilter> FFTHDFWriterFilter::New()
{
  struct make_shared_enabler : public FFTHDFWriterFilter
  {
  };
  std::shared_ptr<make_shared_enabler> val = std::make_shared<make_shared_enabler>();
  val->setupFilterParameters();
  return val;
}

// -----------------------------------------------------------------------------
QString FFTHDFWriterFilter::getNameOfClass() const
{
  return QString("FFTHDFWriterFilter");
}

// -----------------------------------------------------------------------------
QString FFTHDFWriterFilter::ClassName()
{
  return QString("FFTHDFWriterFilter");
}

// -----------------------------------------------------------------------------
void FFTHDFWriterFilter::setOutputFile(const QString& value)
{
  m_OutputFile = value;
}

// -----------------------------------------------------------------------------
QString FFTHDFWriterFilter::getOutputFile() const
{
  return m_OutputFile;
}

// -----------------------------------------------------------------------------
void FFTHDFWriterFilter::setWriteEigenstrains(bool value)
{
  m_WriteEigenstrains = value;
}

// -----------------------------------------------------------------------------
bool FFTHDFWriterFilter::getWriteEigenstrains() const
{
  return m_WriteEigenstrains;
}

// -----------------------------------------------------------------------------
void FFTHDFWriterFilter::setEigenstrainsOutputFile(const QString& value)
{
  m_EigenstrainsOutputFile = value;
}

// -----------------------------------------------------------------------------
QString FFTHDFWriterFilter::getEigenstrainsOutputFile() const
{
  return m_EigenstrainsOutputFile;
}

// -----------------------------------------------------------------------------
void FFTHDFWriterFilter::setWritePipeline(bool value)
{
  m_WritePipeline = value;
}

// -----------------------------------------------------------------------------
bool FFTHDFWriterFilter::getWritePipeline() const
{
  return m_WritePipeline;
}

// -----------------------------------------------------------------------------
void FFTHDFWriterFilter::setAppendToExisting(bool value)
{
  m_AppendToExisting = value;
}

// -----------------------------------------------------------------------------
bool FFTHDFWriterFilter::getAppendToExisting() const
{
  return m_AppendToExisting;
}

// -----------------------------------------------------------------------------
void FFTHDFWriterFilter::setFeatureIdsArrayPath(const DataArrayPath& value)
{
  m_FeatureIdsArrayPath = value;
}

// -----------------------------------------------------------------------------
DataArrayPath FFTHDFWriterFilter::getFeatureIdsArrayPath() const
{
  return m_FeatureIdsArrayPath;
}

// -----------------------------------------------------------------------------
void FFTHDFWriterFilter::setCellPhasesArrayPath(const DataArrayPath& value)
{
  m_CellPhasesArrayPath = value;
}

// -----------------------------------------------------------------------------
DataArrayPath FFTHDFWriterFilter::getCellPhasesArrayPath() const
{
  return m_CellPhasesArrayPath;
}

// -----------------------------------------------------------------------------
void FFTHDFWriterFilter::setCellEulerAnglesArrayPath(const DataArrayPath& value)
{
  m_CellEulerAnglesArrayPath = value;
}

// -----------------------------------------------------------------------------
DataArrayPath FFTHDFWriterFilter::getCellEulerAnglesArrayPath() const
{
  return m_CellEulerAnglesArrayPath;
}

// -----------------------------------------------------------------------------
void FFTHDFWriterFilter::setCellEigenstrainsArrayPath(const DataArrayPath& value)
{
  m_CellEigenstrainsArrayPath = value;
}

// -----------------------------------------------------------------------------
DataArrayPath FFTHDFWriterFilter::getCellEigenstrainsArrayPath() const
{
  return m_CellEigenstrainsArrayPath;
}
