/* ============================================================================
 * Copyright (c) 2012 Michael A. Jackson (BlueQuartz Software)
 * Copyright (c) 2012 Dr. Michael A. Groeber (US Air Force Research Laboratories)
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 * Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 *
 * Redistributions in binary form must reproduce the above copyright notice, this
 * list of conditions and the following disclaimer in the documentation and/or
 * other materials provided with the distribution.
 *
 * Neither the name of Michael A. Groeber, Michael A. Jackson, the US Air Force,
 * BlueQuartz Software nor the names of its contributors may be used to endorse
 * or promote products derived from this software without specific prior written
 * permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
 * USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *  This code was written under United States Air Force Contract number
 *                           FA8650-07-D-5800
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
#include <memory>

#include "ReadMicData.h"

#include <limits>
#include <sstream>
#include <vector>

#include <QtCore/QFileInfo>

#include "EbsdLib/EbsdLib.h"
#include "DREAM3DReview/DREAM3DReviewFilters/HEDM/MicFields.h"
#include "DREAM3DReview/DREAM3DReviewFilters/HEDM/MicReader.h"

#include "SIMPLib/FilterParameters/AbstractFilterParametersReader.h"
#include "SIMPLib/FilterParameters/DataContainerCreationFilterParameter.h"
#include "SIMPLib/FilterParameters/InputFileFilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedPathCreationFilterParameter.h"
#include "SIMPLib/FilterParameters/SeparatorFilterParameter.h"
#include "SIMPLib/FilterParameters/StringFilterParameter.h"
#include "SIMPLib/Geometry/ImageGeom.h"
#include "SIMPLib/DataContainers/DataContainerArray.h"
#include "SIMPLib/DataContainers/DataContainer.h"

#include "DREAM3DReview/DREAM3DReviewConstants.h"

enum createdPathID : RenameDataPath::DataID_t
{
  AttributeMatrixID21 = 21,
  AttributeMatrixID22 = 22,

  DataArrayID31 = 31,
  DataArrayID32 = 32,
  DataArrayID33 = 33,

  DataContainerID = 1
};

/* ############## Start Private Implementation ############################### */
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
class ReadMicDataPrivate
{
  Q_DISABLE_COPY(ReadMicDataPrivate)
  Q_DECLARE_PUBLIC(ReadMicData)
  ReadMicData* const q_ptr;
  ReadMicDataPrivate(ReadMicData* ptr);

  Mic_Private_Data m_Data;

  QString m_InputFile_Cache;
  QDateTime m_TimeStamp_Cache;
};

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
ReadMicDataPrivate::ReadMicDataPrivate(ReadMicData* ptr)
: q_ptr(ptr)
, m_InputFile_Cache("")
, m_TimeStamp_Cache(QDateTime())
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
ReadMicData::ReadMicData()
: m_DataContainerName(SIMPL::Defaults::ImageDataContainerName)
, m_CellEnsembleAttributeMatrixName(SIMPL::Defaults::CellEnsembleAttributeMatrixName)
, m_CellAttributeMatrixName(SIMPL::Defaults::CellAttributeMatrixName)
, m_PhaseNameArrayName("")
, m_MaterialNameArrayName(SIMPL::EnsembleData::PhaseName)
, m_InputFile("")
, m_CellEulerAnglesArrayName(SIMPL::CellData::EulerAngles)
, m_CellPhasesArrayName(SIMPL::CellData::Phases)
, m_CrystalStructuresArrayName(SIMPL::EnsembleData::CrystalStructures)
, m_LatticeConstantsArrayName(SIMPL::EnsembleData::LatticeConstants)
, m_RefFrameZDir(SIMPL::RefFrameZDir::UnknownRefFrameZDirection)
, m_Manufacturer(Ebsd::OEM::Unknown)
, d_ptr(new ReadMicDataPrivate(this))
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
ReadMicData::~ReadMicData() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void ReadMicData::setData(const Mic_Private_Data& value)
{
  Q_D(ReadMicData);
  d->m_Data = value;
}

// -----------------------------------------------------------------------------
Mic_Private_Data ReadMicData::getData() const
{
  Q_D(const ReadMicData);
  return d->m_Data;
}

// -----------------------------------------------------------------------------
void ReadMicData::setInputFile_Cache(const QString& value)
{
  Q_D(ReadMicData);
  d->m_InputFile_Cache = value;
}

// -----------------------------------------------------------------------------
QString ReadMicData::getInputFile_Cache() const
{
  Q_D(const ReadMicData);
  return d->m_InputFile_Cache;
}

// -----------------------------------------------------------------------------
void ReadMicData::setTimeStamp_Cache(const QDateTime& value)
{
  Q_D(ReadMicData);
  d->m_TimeStamp_Cache = value;
}

// -----------------------------------------------------------------------------
QDateTime ReadMicData::getTimeStamp_Cache() const
{
  Q_D(const ReadMicData);
  return d->m_TimeStamp_Cache;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ReadMicData::setupFilterParameters()
{
  FilterParameterVectorType parameters;
  parameters.push_back(SIMPL_NEW_INPUT_FILE_FP("Input File", InputFile, FilterParameter::Parameter, ReadMicData, "*.mic"));
  parameters.push_back(SIMPL_NEW_DC_CREATION_FP("Data Container", DataContainerName, FilterParameter::CreatedArray, ReadMicData));
  parameters.push_back(SeparatorFilterParameter::New("Cell Data", FilterParameter::CreatedArray));
  parameters.push_back(SIMPL_NEW_AM_WITH_LINKED_DC_FP("Cell Attribute Matrix", CellAttributeMatrixName, DataContainerName, FilterParameter::CreatedArray, ReadMicData));
  parameters.push_back(SeparatorFilterParameter::New("Cell Ensemble Data", FilterParameter::CreatedArray));
  parameters.push_back(SIMPL_NEW_AM_WITH_LINKED_DC_FP("Cell Ensemble Attribute Matrix", CellEnsembleAttributeMatrixName, DataContainerName, FilterParameter::CreatedArray, ReadMicData));
  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ReadMicData::readFilterParameters(AbstractFilterParametersReader* reader, int index)
{
  reader->openFilterGroup(this, index);
  setDataContainerName(reader->readDataArrayPath("DataContainerName", getDataContainerName()));
  setCellAttributeMatrixName(reader->readString("CellAttributeMatrixName", getCellAttributeMatrixName()));
  setCellEnsembleAttributeMatrixName(reader->readString("CellEnsembleAttributeMatrixName", getCellEnsembleAttributeMatrixName()));
  setInputFile(reader->readString("InputFile", getInputFile()));
  reader->closeFilterGroup();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ReadMicData::flushCache()
{
  setInputFile_Cache("");
  setTimeStamp_Cache(QDateTime());
  setData(Mic_Private_Data());
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ReadMicData::populateMicData(MicReader* reader, DataContainer::Pointer m, std::vector<size_t> dims, MIC_READ_FLAG flag)
{
  QFileInfo fi(m_InputFile);
  QDateTime timeStamp(fi.lastModified());
  if(flag == MIC_FULL_FILE)
  {
    setInputFile_Cache(""); // We need something to trigger the file read below
  }
  // Drop into this if statement if we need to read from a file
  if(m_InputFile != getInputFile_Cache() || !getTimeStamp_Cache().isValid() || getTimeStamp_Cache() < timeStamp)
  {
    int zDim = 1;
    float zStep = 1.0, xOrigin = 0.0f, yOrigin = 0.0f, zOrigin = 0.0f;

    reader->setFileName(m_InputFile);

    if(flag == MIC_HEADER_ONLY)
    {
      int err = reader->readHeaderOnly();
      if(err < 0)
      {
        setErrorCondition(err, reader->getErrorMessage());
        setErrorCondition(getErrorCode(), "MicReader could not read the .mic file header.");
        m_FileWasRead = false;
        return;
      }

      m_FileWasRead = true;
    }
    else
    {
      int err = reader->readFile();
      if(err < 0)
      {
        setErrorCondition(err, reader->getErrorMessage());
        setErrorCondition(getErrorCode(), "MicReader could not read the .mic file.");
        return;
      }
    }

    dims[0] = reader->getXDimension();
    dims[1] = reader->getYDimension();
    dims[2] = zDim; // We are reading a single slice

    // Set cache with values from file
    {
      Mic_Private_Data data;
      data.dims = SizeVec3Type(dims[0], dims[1], dims[2]);
      data.resolution[0] = reader->getXStep();
      data.resolution[1] = reader->getYStep();
      data.resolution[2] = zStep;
      data.origin[0] = xOrigin;
      data.origin[1] = yOrigin;
      data.origin[2] = zOrigin;
      data.phases = reader->getPhaseVector();
      setData(data);

      setInputFile_Cache(m_InputFile);

      QFileInfo newFi(m_InputFile);
      QDateTime timeStamp(newFi.lastModified());
      setTimeStamp_Cache(timeStamp);
    }
  }
  else
  {
    m_FileWasRead = false;
  }

  // Read from cache
  {
    dims[0] = getData().dims[0];
    dims[1] = getData().dims[1];
    dims[2] = getData().dims[2];
    m->getGeometryAs<ImageGeom>()->setDimensions(dims[0], dims[1], dims[2]);
    m->getGeometryAs<ImageGeom>()->setSpacing(getData().resolution);
    m->getGeometryAs<ImageGeom>()->setOrigin(getData().origin);
  }

  if(flag == MIC_FULL_FILE)
  {
    loadMaterialInfo(reader);
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ReadMicData::initialize()
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ReadMicData::dataCheck()
{
  // Reset FileWasRead flag
  m_FileWasRead = false;

  DataArrayPath tempPath;
  clearErrorCode();
  clearWarningCode();

  DataContainer::Pointer m = getDataContainerArray()->createNonPrereqDataContainer(this, getDataContainerName(), DataContainerID);
  if(getErrorCode() < 0)
  {
    return;
  }

  // Create the Image Geometry
  ImageGeom::Pointer image = ImageGeom::CreateGeometry(SIMPL::Geometry::ImageGeometry);
  m->setGeometry(image);

  std::vector<size_t> tDims(3, 0);
  AttributeMatrix::Pointer cellAttrMat = m->createNonPrereqAttributeMatrix(this, getCellAttributeMatrixName(), tDims, AttributeMatrix::Type::Cell, AttributeMatrixID21);
  if(getErrorCode() < 0)
  {
    return;
  }
  tDims.resize(1);
  tDims[0] = 0;
  AttributeMatrix::Pointer cellEnsembleAttrMat = m->createNonPrereqAttributeMatrix(this, getCellEnsembleAttributeMatrixName(), tDims, AttributeMatrix::Type::CellEnsemble, AttributeMatrixID22);
  if(getErrorCode() < 0)
  {
    return;
  }

  QFileInfo fi(m_InputFile);
  if(!fi.exists())
  {
    QString ss = QObject::tr("The input file does not exist: '%1'").arg(getInputFile());
    setErrorCondition(-388, ss);
  }

  if(m_InputFile.isEmpty() && m_Manufacturer == Ebsd::OEM::Unknown)
  {
    QString ss = QObject::tr("The input file must be set for property %1").arg("InputFile");
    setErrorCondition(-1, ss);
  }

  if(!m_InputFile.isEmpty()) // User set a filename, so lets check it
  {
    std::vector<size_t> dims(3, 0);

    QString ext = fi.suffix();
    QVector<QString> names;
    if(ext.compare(Mic::FileExt) == 0)
    {
      MicReader* reader = new MicReader();

      populateMicData(reader, m, dims, MIC_HEADER_ONLY);

      // Update the size of the Cell Attribute Matrix now that the dimensions of the volume are known
      cellAttrMat->resizeAttributeArrays(dims);
      MicFields micfeatures;
      names = micfeatures.getFilterFeatures<QVector<QString>>();
      dims.resize(1);
      dims[0] = 1;
      for(qint32 i = 0; i < names.size(); ++i)
      {
        if(reader->getPointerType(names[i]) == Ebsd::Int32)
        {
          cellAttrMat->createAndAddAttributeArray<DataArray<int32_t>>(this, names[i], 0, dims);
        }
        else if(reader->getPointerType(names[i]) == Ebsd::Float)
        {
          cellAttrMat->createAndAddAttributeArray<DataArray<float>>(this, names[i], 0, dims);
        }
      }

      delete reader;
    }
    else
    {
      QString ss = QObject::tr("The File extension '%1' was not recognized. The reader only recognizes the .mic file extension").arg(ext);
      setErrorCondition(-997, ss);
      return;
    }

    std::vector<size_t> dim(1, 3);
    tempPath.update(getDataContainerName().getDataContainerName(), getCellAttributeMatrixName(), getCellEulerAnglesArrayName());
    m_CellEulerAnglesPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<float>>(
        this, tempPath, 0, dim);               /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
    if(nullptr != m_CellEulerAnglesPtr.lock()) /* Validate the Weak Pointer wraps a non-nullptr pointer to a DataArray<T> object */
    {
      m_CellEulerAngles = m_CellEulerAnglesPtr.lock()->getPointer(0);
    } /* Now assign the raw pointer to data from the DataArray<T> object */
    dim[0] = 1;
    tempPath.update(getDataContainerName().getDataContainerName(), getCellAttributeMatrixName(), getCellPhasesArrayName());
    m_CellPhasesPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<int32_t>>(
        this, tempPath, 0, dim);          /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
    if(nullptr != m_CellPhasesPtr.lock()) /* Validate the Weak Pointer wraps a non-nullptr pointer to a DataArray<T> object */
    {
      m_CellPhases = m_CellPhasesPtr.lock()->getPointer(0);
    } /* Now assign the raw pointer to data from the DataArray<T> object */

    // typedef DataArray<unsigned int> XTalStructArrayType;
    tempPath.update(getDataContainerName().getDataContainerName(), getCellEnsembleAttributeMatrixName(), getCrystalStructuresArrayName());
    m_CrystalStructuresPtr =
        getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<uint32_t>>(this, tempPath, Ebsd::CrystalStructure::UnknownCrystalStructure, dim, "", DataArrayID33);
    if(nullptr != m_CrystalStructuresPtr.lock()) /* Validate the Weak Pointer wraps a non-nullptr pointer to a DataArray<T> object */
    {
      m_CrystalStructures = m_CrystalStructuresPtr.lock()->getPointer(0);
    } /* Now assign the raw pointer to data from the DataArray<T> object */
    dim[0] = 6;
    tempPath.update(getDataContainerName().getDataContainerName(), getCellEnsembleAttributeMatrixName(), getLatticeConstantsArrayName());
    m_LatticeConstantsPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<float>>(
        this, tempPath, 0.0, dim);              /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
    if(nullptr != m_LatticeConstantsPtr.lock()) /* Validate the Weak Pointer wraps a non-nullptr pointer to a DataArray<T> object */
    {
      m_LatticeConstants = m_LatticeConstantsPtr.lock()->getPointer(0);
    } /* Now assign the raw pointer to data from the DataArray<T> object */

    StringDataArray::Pointer materialNames = StringDataArray::CreateArray(cellEnsembleAttrMat->getNumberOfTuples(), SIMPL::EnsembleData::PhaseName, true);
    cellEnsembleAttrMat->insertOrAssign(materialNames);
  }
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ReadMicData::execute()
{
  clearErrorCode();
  clearWarningCode();

  dataCheck();
  if(getErrorCode() < 0)
  {
    return;
  }

  QFileInfo fi(getInputFile());
  QString ext = fi.suffix();

  readMicFile();

  // Set the file name and time stamp into the cache, if we are reading from the file and after all the reading has been done
  {
    QFileInfo newFi(m_InputFile);
    QDateTime timeStamp(newFi.lastModified());

    if(m_InputFile == getInputFile_Cache() && getTimeStamp_Cache().isValid() && getTimeStamp_Cache() >= timeStamp)
    {
      setTimeStamp_Cache(timeStamp);
      setInputFile_Cache(m_InputFile);
    }
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ReadMicData::readMicFile()
{
  int err = 0;
  std::shared_ptr<MicReader> reader(new MicReader());
  reader->setFileName(m_InputFile);
  err = reader->readFile();
  if(err < 0)
  {
    setErrorCondition(err, reader->getErrorMessage());
    return;
  }
  DataContainer::Pointer m = getDataContainerArray()->getDataContainer(getDataContainerName());
  AttributeMatrix::Pointer cellAttrMat = m->getAttributeMatrix(getCellAttributeMatrixName());

  int64_t dims[3];
  dims[0] = reader->getXDimension();
  dims[1] = reader->getYDimension();
  dims[2] = 1; // We are reading a single slice
  m->getGeometryAs<ImageGeom>()->setDimensions(dims[0], dims[1], dims[2]);
  m->getGeometryAs<ImageGeom>()->setSpacing(FloatVec3Type(reader->getXStep(), reader->getYStep(), 1.0));
  m->getGeometryAs<ImageGeom>()->setOrigin(FloatVec3Type(0.0f, 0.0f, 0.0f));

  err = loadMaterialInfo(reader.get());

  float* f1 = nullptr;
  float* f2 = nullptr;
  float* f3 = nullptr;
  int* phasePtr = nullptr;

  FloatArrayType::Pointer fArray = FloatArrayType::NullPointer();
  Int32ArrayType::Pointer iArray = Int32ArrayType::NullPointer();
  size_t totalPoints = m->getGeometryAs<ImageGeom>()->getNumberOfElements();
  // Prepare the Cell Attribute Matrix with the correct number of tuples based on the total points being read from the file.
  std::vector<size_t> tDims(3, 0);
  tDims[0] = m->getGeometryAs<ImageGeom>()->getXPoints();
  tDims[1] = m->getGeometryAs<ImageGeom>()->getYPoints();
  tDims[2] = m->getGeometryAs<ImageGeom>()->getZPoints();
  cellAttrMat->resizeAttributeArrays(tDims);

  float x, y;
  float xMin = 10000000;
  float yMin = 10000000;
  f1 = reinterpret_cast<float*>(reader->getPointerByName(Mic::X));
  f2 = reinterpret_cast<float*>(reader->getPointerByName(Mic::Y));
  for(size_t i = 0; i < totalPoints; i++)
  {
    x = f1[i];
    y = f2[i];
    if(x < xMin)
    {
      xMin = x;
    }
    if(y < yMin)
    {
      yMin = y;
    }
  }
  m->getGeometryAs<ImageGeom>()->setOrigin(FloatVec3Type(xMin, yMin, 0.0));

  {
    phasePtr = reinterpret_cast<int*>(reader->getPointerByName(Mic::Phase));
    for(size_t i = 0; i < totalPoints; i++)
    {
      if(phasePtr[i] < 1)
      {
        phasePtr[i] = 1;
      }
    }
    iArray = Int32ArrayType::CreateArray(totalPoints, SIMPL::CellData::Phases, true);
    ::memcpy(iArray->getPointer(0), phasePtr, sizeof(int32_t) * totalPoints);
    cellAttrMat->insertOrAssign(iArray);
  }

  std::vector<size_t> compDims(1, 3); // Initially set this up for the Euler Angle 1x3
  {
    //  radianconversion = M_PI / 180.0;
    f1 = reinterpret_cast<float*>(reader->getPointerByName(Mic::Euler1));
    f2 = reinterpret_cast<float*>(reader->getPointerByName(Mic::Euler2));
    f3 = reinterpret_cast<float*>(reader->getPointerByName(Mic::Euler3));
    fArray = FloatArrayType::CreateArray(totalPoints, compDims, SIMPL::CellData::EulerAngles, true);
    float* cellEulerAngles = fArray->getPointer(0);
    for(size_t i = 0; i < totalPoints; i++)
    {
      cellEulerAngles[3 * i] = f1[i];
      cellEulerAngles[3 * i + 1] = f2[i];
      cellEulerAngles[3 * i + 2] = f3[i];
    }
    cellAttrMat->insertOrAssign(fArray);
  }

  compDims[0] = 1; // Now reset the size of the first dimension to 1
  {
    phasePtr = reinterpret_cast<int*>(reader->getPointerByName(Mic::Phase));
    iArray = Int32ArrayType::CreateArray(totalPoints, compDims, SIMPL::CellData::Phases, true);
    ::memcpy(iArray->getPointer(0), phasePtr, sizeof(int32_t) * totalPoints);
    cellAttrMat->insertOrAssign(iArray);
  }

  {
    f1 = reinterpret_cast<float*>(reader->getPointerByName(Mic::Confidence));
    fArray = FloatArrayType::CreateArray(totalPoints, compDims, Mic::Confidence, true);
    ::memcpy(fArray->getPointer(0), f1, sizeof(float) * totalPoints);
    cellAttrMat->insertOrAssign(fArray);
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AbstractFilter::Pointer ReadMicData::newFilterInstance(bool copyFilterParameters) const
{
  ReadMicData::Pointer filter = ReadMicData::New();
  if(copyFilterParameters)
  {
    filter->setFilterParameters(getFilterParameters());
    copyFilterParameterInstanceVariables(filter.get());
  }
  return filter;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ReadMicData::getCompiledLibraryName() const
{
  return HEDMAnalysisConstants::HEDMAnalysisBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ReadMicData::getGroupName() const
{
  return SIMPL::FilterGroups::Unsupported;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QUuid ReadMicData::getUuid() const
{
  return QUuid("{5a2b714e-bae9-5213-8171-d1e190609e2d}");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ReadMicData::getSubGroupName() const
{
  return SIMPL::FilterSubGroups::InputFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ReadMicData::getHumanLabel() const
{
  return "Import HEDM Data (.mic)";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int ReadMicData::loadMaterialInfo(MicReader* reader)
{
  QVector<MicPhase::Pointer> phases = getData().phases;
  if(phases.empty())
  {
    setErrorCondition(reader->getErrorCode(), reader->getErrorMessage());
    return getErrorCode();
  }

  DataArray<unsigned int>::Pointer crystalStructures = DataArray<unsigned int>::CreateArray(phases.size() + 1, getCrystalStructuresArrayName(), true);
  StringDataArray::Pointer materialNames = StringDataArray::CreateArray(phases.size() + 1, getMaterialNameArrayName());
  std::vector<size_t> dims(1, 6);
  FloatArrayType::Pointer latticeConstants = FloatArrayType::CreateArray(phases.size() + 1, dims, getLatticeConstantsArrayName(), true);

  // Initialize the zero'th element to unknowns. The other elements will
  // be filled in based on values from the data file
  crystalStructures->setValue(0, Ebsd::CrystalStructure::UnknownCrystalStructure);
  materialNames->setValue(0, "Invalid Phase");
  latticeConstants->setComponent(0, 0, 0.0f);
  latticeConstants->setComponent(0, 1, 0.0f);
  latticeConstants->setComponent(0, 2, 0.0f);
  latticeConstants->setComponent(0, 3, 0.0f);
  latticeConstants->setComponent(0, 4, 0.0f);
  latticeConstants->setComponent(0, 5, 0.0f);

  for(size_t i = 0; i < phases.size(); i++)
  {
    int phaseID = phases[i]->getPhaseIndex();
    crystalStructures->setValue(phaseID, phases[i]->determineLaueGroup());
    materialNames->setValue(phaseID, phases[i]->getMaterialName());
    QVector<float> lc = phases[i]->getLatticeConstants();

    latticeConstants->setComponent(phaseID, 0, lc[0]);
    latticeConstants->setComponent(phaseID, 1, lc[1]);
    latticeConstants->setComponent(phaseID, 2, lc[2]);
    latticeConstants->setComponent(phaseID, 3, lc[3]);
    latticeConstants->setComponent(phaseID, 4, lc[4]);
    latticeConstants->setComponent(phaseID, 5, lc[5]);
  }
  DataContainer::Pointer vdc = getDataContainerArray()->getDataContainer(getDataContainerName());
  if(nullptr == vdc)
  {
    return -1;
  }
  AttributeMatrix::Pointer attrMatrix = vdc->getAttributeMatrix(getCellEnsembleAttributeMatrixName());
  if(nullptr == attrMatrix.get())
  {
    return -2;
  }

  // Resize the AttributeMatrix based on the size of the crystal structures array
  std::vector<size_t> tDims(1, crystalStructures->getNumberOfTuples());
  attrMatrix->resizeAttributeArrays(tDims);
  // Now add the attributeArray to the AttributeMatrix
  attrMatrix->insertOrAssign(crystalStructures);
  attrMatrix->insertOrAssign(materialNames);
  attrMatrix->insertOrAssign(latticeConstants);

  // Now reset the internal ensemble array references to these new arrays
  m_CrystalStructuresPtr = crystalStructures;
  if(nullptr != m_CrystalStructuresPtr.lock()) /* Validate the Weak Pointer wraps a non-nullptr pointer to a DataArray<T> object */
  {
    m_CrystalStructures = m_CrystalStructuresPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */

  m_LatticeConstantsPtr = latticeConstants;
  if(nullptr != m_LatticeConstantsPtr.lock()) /* Validate the Weak Pointer wraps a non-nullptr pointer to a DataArray<T> object */
  {
    m_LatticeConstants = m_LatticeConstantsPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */
  return 0;
}

// -----------------------------------------------------------------------------
ReadMicData::Pointer ReadMicData::NullPointer()
{
  return Pointer(static_cast<Self*>(nullptr));
}

// -----------------------------------------------------------------------------
std::shared_ptr<ReadMicData> ReadMicData::New()
{
  struct make_shared_enabler : public ReadMicData
  {
  };
  std::shared_ptr<make_shared_enabler> val = std::make_shared<make_shared_enabler>();
  val->setupFilterParameters();
  return val;
}

// -----------------------------------------------------------------------------
QString ReadMicData::getNameOfClass() const
{
  return QString("ReadMicData");
}

// -----------------------------------------------------------------------------
QString ReadMicData::ClassName()
{
  return QString("ReadMicData");
}

// -----------------------------------------------------------------------------
void ReadMicData::setDataContainerName(const DataArrayPath& value)
{
  m_DataContainerName = value;
}

// -----------------------------------------------------------------------------
DataArrayPath ReadMicData::getDataContainerName() const
{
  return m_DataContainerName;
}

// -----------------------------------------------------------------------------
void ReadMicData::setCellEnsembleAttributeMatrixName(const QString& value)
{
  m_CellEnsembleAttributeMatrixName = value;
}

// -----------------------------------------------------------------------------
QString ReadMicData::getCellEnsembleAttributeMatrixName() const
{
  return m_CellEnsembleAttributeMatrixName;
}

// -----------------------------------------------------------------------------
void ReadMicData::setCellAttributeMatrixName(const QString& value)
{
  m_CellAttributeMatrixName = value;
}

// -----------------------------------------------------------------------------
QString ReadMicData::getCellAttributeMatrixName() const
{
  return m_CellAttributeMatrixName;
}

// -----------------------------------------------------------------------------
void ReadMicData::setFileWasRead(bool value)
{
  m_FileWasRead = value;
}

// -----------------------------------------------------------------------------
bool ReadMicData::getFileWasRead() const
{
  return m_FileWasRead;
}

// -----------------------------------------------------------------------------
void ReadMicData::setPhaseNameArrayName(const QString& value)
{
  m_PhaseNameArrayName = value;
}

// -----------------------------------------------------------------------------
QString ReadMicData::getPhaseNameArrayName() const
{
  return m_PhaseNameArrayName;
}

// -----------------------------------------------------------------------------
void ReadMicData::setMaterialNameArrayName(const QString& value)
{
  m_MaterialNameArrayName = value;
}

// -----------------------------------------------------------------------------
QString ReadMicData::getMaterialNameArrayName() const
{
  return m_MaterialNameArrayName;
}

// -----------------------------------------------------------------------------
void ReadMicData::setInputFile(const QString& value)
{
  m_InputFile = value;
}

// -----------------------------------------------------------------------------
QString ReadMicData::getInputFile() const
{
  return m_InputFile;
}

// -----------------------------------------------------------------------------
void ReadMicData::setCellEulerAnglesArrayName(const QString& value)
{
  m_CellEulerAnglesArrayName = value;
}

// -----------------------------------------------------------------------------
QString ReadMicData::getCellEulerAnglesArrayName() const
{
  return m_CellEulerAnglesArrayName;
}

// -----------------------------------------------------------------------------
void ReadMicData::setCellPhasesArrayName(const QString& value)
{
  m_CellPhasesArrayName = value;
}

// -----------------------------------------------------------------------------
QString ReadMicData::getCellPhasesArrayName() const
{
  return m_CellPhasesArrayName;
}

// -----------------------------------------------------------------------------
void ReadMicData::setCrystalStructuresArrayName(const QString& value)
{
  m_CrystalStructuresArrayName = value;
}

// -----------------------------------------------------------------------------
QString ReadMicData::getCrystalStructuresArrayName() const
{
  return m_CrystalStructuresArrayName;
}

// -----------------------------------------------------------------------------
void ReadMicData::setLatticeConstantsArrayName(const QString& value)
{
  m_LatticeConstantsArrayName = value;
}

// -----------------------------------------------------------------------------
QString ReadMicData::getLatticeConstantsArrayName() const
{
  return m_LatticeConstantsArrayName;
}

// -----------------------------------------------------------------------------
void ReadMicData::setRefFrameZDir(uint32_t value)
{
  m_RefFrameZDir = value;
}

// -----------------------------------------------------------------------------
uint32_t ReadMicData::getRefFrameZDir() const
{
  return m_RefFrameZDir;
}

// -----------------------------------------------------------------------------
void ReadMicData::setManufacturer(const Ebsd::OEM& value)
{
  m_Manufacturer = value;
}

// -----------------------------------------------------------------------------
Ebsd::OEM ReadMicData::getManufacturer() const
{
  return m_Manufacturer;
}
