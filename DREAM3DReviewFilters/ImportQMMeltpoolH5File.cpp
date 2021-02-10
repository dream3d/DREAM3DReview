/* ============================================================================
 * Copyright (c) 2020 BlueQuartz Software, LLC
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
 * Neither the names of any of the BlueQuartz Software contributors
 * may be used to endorse or promote products derived from this software without
 * specific prior written permission.
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
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include "ImportQMMeltpoolH5File.h"

#include <algorithm>

#include <QtCore/QDateTime>
#include <QtCore/QFileInfo>

#include "H5Support/H5Lite.h"
#include "H5Support/H5ScopedSentinel.h"
#include "H5Support/H5Utilities.h"

#include "EbsdLib/EbsdLib.h"

#include "SIMPLib/Common/Constants.h"
#include "SIMPLib/DataContainers/AttributeMatrix.h"
#include "SIMPLib/DataContainers/DataContainer.h"
#include "SIMPLib/DataContainers/DataContainerArray.h"
#include "SIMPLib/FilterParameters/DataContainerCreationFilterParameter.h"
#include "SIMPLib/FilterParameters/FileListInfoFilterParameter.h"
#include "SIMPLib/FilterParameters/FloatFilterParameter.h"
#include "SIMPLib/FilterParameters/InputFileFilterParameter.h"
#include "SIMPLib/FilterParameters/IntVec2FilterParameter.h"
#include "SIMPLib/FilterParameters/MultiInputFileFilterParameter.h"
#include "SIMPLib/FilterParameters/PreflightUpdatedValueFilterParameter.h"
#include "SIMPLib/FilterParameters/StringFilterParameter.h"
#include "SIMPLib/Geometry/VertexGeom.h"

#include "DREAM3DReview/DREAM3DReviewConstants.h"
#include "DREAM3DReview/DREAM3DReviewVersion.h"

namespace
{
constexpr int32_t k_FileDoesNotExist = -23700;
constexpr int32_t k_FileTypeError = -23701;
constexpr int32_t k_FilePathEmptyError = -23702;
constexpr int32_t k_HDF5FileOpenError = -23703;
constexpr int32_t k_HDF5GroupOpenError = -23704;
constexpr int32_t k_HDF5DatasetError = -23705;
constexpr int32_t k_HDF5AttributeError = -23706;
constexpr int32_t k_StartEndError = -23707;
constexpr int32_t k_SliceRangeError = -23708;
constexpr int32_t k_MissingSlicesError = -23709;
constexpr int32_t k_DataStructureError = -23710;
constexpr int32_t k_NumElementsError = -23711;
constexpr int32_t k_InvalidOffsetError = -23712;
// constexpr int32_t k_InvalidSlicesError = -23713;
constexpr int32_t k_DataContainerError = -23714;
constexpr int32_t k_FileFormatVersionMissingError = -23715;
constexpr int32_t k_FileFormatVersionError = -23716;

constexpr int64_t k_FileFormatVersion = 3;

constexpr size_t k_IndexNumCols = 3;
// constexpr size_t k_IndexCol = 0;
constexpr size_t k_IndexLayerThicknessCol = 1;
constexpr size_t k_IndexVerticesCol = 2;

const std::string k_Index = "Index";

const std::string k_TDMSData = "TDMSData";

const std::string k_Slice = "Slice";

const std::string k_Area = "Area";
const std::string k_Intensity = "Intensity";
const std::string k_LaserTTL = "LaserTTL";
const std::string k_XAxis = "X-Axis";
const std::string k_YAxis = "Y-Axis";

const std::string k_BitgainOS1 = "Bitgain OS 1";
const std::string k_BitgainOS2 = "Bitgain OS 2";

const std::string k_PartStartTime = "PartStartTime";
const std::string k_PartEndTime = "PartEndTime";

const std::string k_LayerThickness = "layerThickness";

const QString k_Power = "Power";
const QString k_Time = "Time";

template <class Container>
std::string make_list_string(const Container& items)
{
  std::stringstream ss;
  ss << '[' << items.front();
  for(size_t i = 1; i < items.size(); i++)
  {
    ss << ", " << items[i];
  }
  ss << ']';

  return ss.str();
}

/**
 * @brief Returns Error Code and Error Message if the checks do not pass
 * @param filePath The file path to check
 * @return Error Code (<0 is bad) and Error Message if something didn't pass
 */
std::pair<int32_t, std::string> checkFile(const std::string& filePath)
{
  if(filePath.empty())
  {
    return {k_FilePathEmptyError, "The HDF5 file path is empty.  Please select an HDF5 file."};
  }

  fs::path ifPath = fs::path(filePath);
  // Make sure the file exists on disk
  if(!fs::exists(ifPath))
  {
    std::stringstream ss;
    ss << "Input File does not exist at '" << ifPath.string() << "'\n";
    return {k_FileDoesNotExist, ss.str()};
  }
  fs::path ext = ifPath.extension();

  if(ext != ".h5" && ext != ".hdf5")
  {
    std::stringstream ss;
    ss << "The selected file '%1' is not an HDF5 file." << filePath;
    return {k_FileTypeError, ss.str()};
  }

  return {0, ""};
}

} // namespace

struct ImportQMMeltpoolH5File::Cache
{
  fs::file_time_type lastModified;
  std::string filePath;
  IntVec2Type sliceRange;

  std::vector<int64_t> missingIndices;
  std::vector<int64_t> layerThicknesses;
  std::vector<int64_t> numElements;

  void flush()
  {
    filePath = "";
    sliceRange = {};

    missingIndices.clear();

    layerThicknesses.clear();
    numElements.clear();
  }
};

// -----------------------------------------------------------------------------
ImportQMMeltpoolH5File::ImportQMMeltpoolH5File() = default;

// -----------------------------------------------------------------------------
ImportQMMeltpoolH5File::~ImportQMMeltpoolH5File() = default;

// -----------------------------------------------------------------------------
void ImportQMMeltpoolH5File::createUpdateCacheEntries()
{
  // Check all of the files
  if(m_Caches.size() != m_InputFiles.size())
  {
    m_Caches.resize(m_InputFiles.size());
  }

  for(size_t i = 0; i < m_InputFiles.size(); i++)
  {
    std::string& inputFile = m_InputFiles.at(i);
    Cache& cache = m_Caches.at(i);

    fs::path ifPath = fs::path(inputFile);
    auto timeStamp = fs::last_write_time(ifPath);

    if(!(cache.filePath == inputFile && timeStamp == cache.lastModified && cache.sliceRange == m_SliceRange))
    {
      // Read from the file and cache all of the slices and slice data
      // Will need to read the number of tuples in the X-Axis and store for each slice
      // User might be able to control which of the 3 data sets that they want to import
      // Building up the Data Structure with this information will allow the Geometry to be fully specified
      cache.flush();
      cache.filePath = inputFile;
      generateCache(cache);
      if(cache.missingIndices.empty())
      {
        cache.filePath = inputFile;
        cache.lastModified = timeStamp;
        cache.sliceRange = m_SliceRange;
      }
    }

    if(!cache.missingIndices.empty())
    {
      std::string index_list = make_list_string(cache.missingIndices);
      setErrorCondition(k_MissingSlicesError, QString("Slices %1 in the given range are missing from the file").arg(QString::fromStdString(index_list)));
      return;
    }
  }
}

// -----------------------------------------------------------------------------
void ImportQMMeltpoolH5File::setupFilterParameters()
{
  FilterParameterVectorType parameters;
  parameters.push_back(SIMPL_NEW_MULTI_INPUT_FILE_FP("Input File(s)", InputFiles, FilterParameter::Parameter, ImportQMMeltpoolH5File, "*.h5 *.hdf5"));
  parameters.push_back(SIMPL_NEW_PREFLIGHTUPDATEDVALUE_FP("Possible Slice Indices", PossibleIndices, FilterParameter::Parameter, ImportQMMeltpoolH5File));
  parameters.push_back(SIMPL_NEW_INT_VEC2_FP("Slice Index Start/End [Inclusive]", SliceRange, FilterParameter::Parameter, ImportQMMeltpoolH5File));
  parameters.push_back(SIMPL_NEW_DC_CREATION_FP("Data Container Name", DataContainerPath, FilterParameter::Parameter, ImportQMMeltpoolH5File));
  parameters.push_back(SIMPL_NEW_STRING_FP("Vertex Attribute Matrix Name", VertexAttributeMatrixName, FilterParameter::Parameter, ImportQMMeltpoolH5File));
  parameters.push_back(SIMPL_NEW_FLOAT_FP("Power", Power, FilterParameter::Parameter, ImportQMMeltpoolH5File));
  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
void ImportQMMeltpoolH5File::dataCheck()
{
  clearErrorCode();
  clearWarningCode();
  setCancel(false);

  auto dca = getDataContainerArray();

  if(dca == nullptr)
  {
    setErrorCondition(k_DataStructureError, "Cannot acquire DataContainerArray");
    return;
  }

  // Check all of the files
  for(const auto& inputFile : m_InputFiles)
  {
    std::pair<int32_t, std::string> result = ::checkFile(inputFile);
    if(result.first < 0)
    {
      setErrorCondition(result.first, S2Q(result.second));
    }
  }

  // Check the slice range
  if(m_SliceRange[0] < 0 || m_SliceRange[1] < 0)
  {
    setErrorCondition(k_StartEndError, "Start/end must be 0 or positive");
    return;
  }

  if(m_SliceRange[1] - m_SliceRange[0] < 0)
  {
    setErrorCondition(k_SliceRangeError, "Slice cannot be negative");
    return;
  }

  if(dca->doesDataContainerExist(m_DataContainerPath))
  {
    setErrorCondition(k_DataContainerError, QString("DataContainer \"%1\" already exists").arg(m_DataContainerPath.getDataContainerName()));
    return;
  }
  createUpdateCacheEntries();
  generateDataStructure();
}

// -----------------------------------------------------------------------------
void ImportQMMeltpoolH5File::execute()
{
  dataCheck();
  if(getErrorCode() < 0)
  {
    return;
  }

  if(getCancel())
  {
    return;
  }

  readDataFromFile();
}

// -----------------------------------------------------------------------------
void ImportQMMeltpoolH5File::generateCache(ImportQMMeltpoolH5File::Cache& cache)
{
  const std::string& filePath = cache.filePath;

  hid_t fileId = H5Utilities::openFile(filePath, true);
  if(fileId < 0)
  {
    setErrorCondition(k_HDF5FileOpenError, "Error opening HDF5 file.");
    return;
  }
  H5ScopedFileSentinel sentinel(fileId, true);

  int64_t fileFormatVersion = 0;
  int32_t err = H5Lite::readScalarAttribute(fileId, "Version", fileFormatVersion);
  if(err < 0)
  {
    setErrorCondition(k_FileFormatVersionMissingError, "The HDF5 Attribute 'Version' for the '/' group is missing.");
    return;
  }
  if(fileFormatVersion != k_FileFormatVersion)
  {
    setErrorCondition(k_FileFormatVersionError, "The file format version is not correct. This filter only accepts version 3");
    return;
  }

  std::vector<int64_t> indexDataset;

  if(H5Lite::readVectorDataset(fileId, k_Index, indexDataset) < 0)
  {
    setErrorCondition(k_HDF5GroupOpenError, "Error opening HDF5 group.");
    return;
  }

  std::vector<int64_t> indices(indexDataset.size() / k_IndexNumCols, 0);

  for(size_t i = 0; i < indices.size(); i++)
  {
    indices[i] = indexDataset[i * k_IndexNumCols];
  }

  m_PossibleIndices = QString("Min = %1\nMax = %2").arg(QString::number(indices.front()), QString::number(indices.back()));

  std::vector<int64_t> layerThicknesses;
  std::vector<int64_t> numElements;

  std::vector<int64_t> missingIndices;

  size_t numSlices = static_cast<size_t>(m_SliceRange[1]) - m_SliceRange[0] + 1;

  layerThicknesses.reserve(numSlices);
  numElements.reserve(numSlices);

  for(int32_t slice = m_SliceRange[0]; slice <= m_SliceRange[1]; slice++)
  {
    auto match = std::find(indices.cbegin(), indices.cend(), slice);
    if(match != indices.cend())
    {
      size_t row = std::distance(indices.cbegin(), match);
      int64_t sliceVertices = indexDataset[row * k_IndexNumCols + k_IndexVerticesCol];
      numElements.push_back(sliceVertices);
      layerThicknesses.push_back(indexDataset[row * k_IndexNumCols + k_IndexLayerThicknessCol]);
    }
    else
    {
      missingIndices.push_back(slice);
      continue;
    }
  }

  cache.missingIndices = std::move(missingIndices);
  cache.layerThicknesses = std::move(layerThicknesses);
  cache.numElements = std::move(numElements);
}

// -----------------------------------------------------------------------------
void ImportQMMeltpoolH5File::generateDataStructure()
{
  DataContainerArray::Pointer dca = getDataContainerArray();
  // Create the Data Container
  DataContainer::Pointer dc = DataContainer::New(getDataContainerPath());
  dca->addOrReplaceDataContainer(dc);

  // Compute the total number of vertices
  size_t numVerts = 0;
  for(const auto& cache : m_Caches)
  {
    numVerts += std::accumulate(cache.numElements.cbegin(), cache.numElements.cend(), 0ULL);
  }

  if(!getInPreflight())
  {
    notifyStatusMessage("Allocating Vertex Geometry");
  }
  // Create the Geometry
  VertexGeom::Pointer vertGeom = VertexGeom::CreateGeometry(numVerts, "Vertex Geometry", !getInPreflight());
  vertGeom->setUnits(IGeometry::LengthUnit::Micrometer);
  dc->setGeometry(vertGeom);

  // Create an Attribute Matrix
  AttributeMatrix::Pointer vertAm = AttributeMatrix::New({numVerts}, m_VertexAttributeMatrixName, AttributeMatrix::Type::Vertex);
  dc->addOrReplaceAttributeMatrix(vertAm);

  // Create the various Vertex Data Arrays
  // The list is hard coded for now
  if(!getInPreflight())
  {
    notifyStatusMessage("Allocating " + QString::fromStdString(k_Area) + " Array");
  }
  Int16ArrayType::Pointer areaData = Int16ArrayType::CreateArray(numVerts, k_Area, !getInPreflight());
  vertAm->addOrReplaceAttributeArray(areaData);

  if(!getInPreflight())
  {
    notifyStatusMessage("Allocating " + QString::fromStdString(k_Intensity) + " Array");
  }
  Int16ArrayType::Pointer intensityData = Int16ArrayType::CreateArray(numVerts, k_Intensity, !getInPreflight());
  vertAm->addOrReplaceAttributeArray(intensityData);

  if(!getInPreflight())
  {
    notifyStatusMessage("Allocating " + QString::fromStdString(k_LaserTTL) + " Array");
  }
  UInt8ArrayType::Pointer laserTtlData = UInt8ArrayType::CreateArray(numVerts, k_LaserTTL, !getInPreflight());
  vertAm->addOrReplaceAttributeArray(laserTtlData);

  if(!getInPreflight())
  {
    notifyStatusMessage("Allocating " + QString::fromStdString(k_Slice) + " Array");
  }
  // This is an extra data set that we are going to add to aid in visualization
  Int16ArrayType::Pointer sliceData = Int16ArrayType::CreateArray(numVerts, k_Slice, !getInPreflight());
  vertAm->addOrReplaceAttributeArray(sliceData);

  if(!getInPreflight())
  {
    notifyStatusMessage("Allocating " + k_Time + " Array");
  }
  DoubleArrayType::Pointer timeData = DoubleArrayType::CreateArray(numVerts, k_Time, !getInPreflight());
  vertAm->addOrReplaceAttributeArray(timeData);

  if(!getInPreflight())
  {
    notifyStatusMessage("Allocating " + k_Power + " Array");
  }
  FloatArrayType::Pointer powerData = FloatArrayType::CreateArray(numVerts, k_Power, !getInPreflight());
  vertAm->addOrReplaceAttributeArray(powerData);

  if(!getInPreflight())
  {
    std::fill(powerData->begin(), powerData->end(), m_Power);
  }
}

// -----------------------------------------------------------------------------
void ImportQMMeltpoolH5File::readDataFromFile()
{

  DataContainer::Pointer dc = getDataContainerArray()->getDataContainer(getDataContainerPath());
  if(dc == nullptr)
  {
    setErrorCondition(k_DataStructureError, "Error acquiring DataContainer");
    return;
  }
  AttributeMatrix::Pointer vertAM = dc->getAttributeMatrix(m_VertexAttributeMatrixName);
  if(vertAM == nullptr)
  {
    setErrorCondition(k_DataStructureError, "Error acquiring AttributeMatrix");
    return;
  }

  size_t numTuples = vertAM->getNumberOfTuples();

  // Now Create the Vertex Geometry
  VertexGeom::Pointer vertGeom = dc->getGeometryAs<VertexGeom>();
  if(vertGeom == nullptr)
  {
    setErrorCondition(k_DataStructureError, "Error acquiring DataContainer");
    return;
  }

  size_t offset = 0;
  int32_t fileNum = 0;
  for(const auto& cache : m_Caches)
  {
    if(getCancel())
    {
      return;
    }
    QString msg;
    QTextStream out(&msg);
    out << fileNum++ << "/" << m_Caches.size() << ": Reading File: " + QString::fromStdString(cache.filePath);
    notifyStatusMessage(msg.toLatin1().data());

    hid_t fileId = H5Utilities::openFile(cache.filePath, true);
    if(fileId < 0)
    {
      setErrorCondition(k_HDF5FileOpenError, "Error opening HDF5 file.");
      return;
    }
    H5ScopedFileSentinel sentinel(fileId, true);

    // Find the min/max tuple count of all the slices
    auto max = std::max_element(cache.numElements.cbegin(), cache.numElements.cend());

    if(max == cache.numElements.cend())
    {
      setErrorCondition(k_NumElementsError, "Couldn't get number of elements");
      return;
    }

    // Create vectors large enough to hold the largest slice data. Cuts down on memory allocations
    std::vector<float> xCoord(*max);
    std::vector<float> yCoord(*max);

    auto areaData = vertAM->getAttributeArrayAs<Int16ArrayType>(QString::fromStdString(k_Area));
    if(areaData == nullptr)
    {
      setErrorCondition(k_DataStructureError, "Failed to acquire Area DataArray");
      return;
    }
    auto intensityData = vertAM->getAttributeArrayAs<Int16ArrayType>(QString::fromStdString(k_Intensity));
    if(intensityData == nullptr)
    {
      setErrorCondition(k_DataStructureError, "Failed to acquire Intensity DataArray");
      return;
    }
    auto laserTtlData = vertAM->getAttributeArrayAs<UInt8ArrayType>(QString::fromStdString(k_LaserTTL));
    if(laserTtlData == nullptr)
    {
      setErrorCondition(k_DataStructureError, "Failed to acquire LaserTTL DataArray");
      return;
    }
    auto sliceData = vertAM->getAttributeArrayAs<Int16ArrayType>(QString::fromStdString(k_Slice));
    if(sliceData == nullptr)
    {
      setErrorCondition(k_DataStructureError, "Failed to acquire Slice DataArray");
      return;
    }
    auto timeData = vertAM->getAttributeArrayAs<DoubleArrayType>(k_Time);
    if(timeData == nullptr)
    {
      setErrorCondition(k_DataStructureError, "Failed to acquire Time DataArray");
      return;
    }

    float cummulativeLayerThickness = 0.0F; // Assumes microns? Maybe?
    double currentTime = 0.0;
    QDateTime initialTime;

    hid_t dataGroup = H5Utilities::openHDF5Object(fileId, k_TDMSData);
    if(dataGroup < 0)
    {
      setErrorCondition(k_HDF5FileOpenError, "Failed to open HDF5 file");
      return;
    }
    H5ScopedGroupSentinel dataGroupSentinel(dataGroup, true);

    {
      hid_t sliceGroup = H5Utilities::openHDF5Object(dataGroup, std::to_string(m_SliceRange[0]));
      if(sliceGroup < 0)
      {
        setErrorCondition(k_HDF5GroupOpenError, "Failed to open HDF5 slice group");
        return;
      }
      H5ScopedGroupSentinel sliceGroupSentinel(sliceGroup, true);

      std::string partStartTimeString;
      herr_t err = H5Lite::readStringAttribute(sliceGroup, k_PartStartTime, partStartTimeString);
      if(err < 0)
      {
        setErrorCondition(k_HDF5AttributeError, "Failed to open HDF5 attribute PartStartTime");
        return;
      }

      initialTime = QDateTime::fromString(QString::fromStdString(partStartTimeString), Qt::DateFormat::ISODateWithMs);
    }

    // Loop over each slice group
    for(auto slice = static_cast<size_t>(m_SliceRange[0]); slice <= static_cast<size_t>(m_SliceRange[1]); slice++)
    {
      size_t sliceIndex = slice - m_SliceRange[0];

      int64_t layerThickness = cache.layerThicknesses[sliceIndex];
      cummulativeLayerThickness += static_cast<float>(layerThickness);

      int64_t sliceElements = cache.numElements[sliceIndex];

      if(sliceElements == 0)
      {
        continue;
      }

      if(sliceElements + offset > numTuples)
      {
        setErrorCondition(k_InvalidOffsetError, "Invalid offset");
        return;
      }

      hid_t sliceGroup = H5Utilities::openHDF5Object(dataGroup, std::to_string(slice));
      if(sliceGroup < 0)
      {
        setErrorCondition(k_HDF5GroupOpenError, "Failed to open HDF5 slice group");
        return;
      }
      H5ScopedGroupSentinel sliceGroupSentinel(sliceGroup, true);

      herr_t err = H5Lite::readPointerDataset(sliceGroup, k_Area, areaData->getTuplePointer(offset));
      if(err < 0)
      {
        setErrorCondition(k_HDF5DatasetError, "Failed to open HDF5 dataset Area");
        return;
      }

      err = H5Lite::readPointerDataset(sliceGroup, k_Intensity, intensityData->getTuplePointer(offset));
      if(err < 0)
      {
        setErrorCondition(k_HDF5DatasetError, "Failed to open HDF5 dataset Intensity");
        return;
      }

      err = H5Lite::readPointerDataset(sliceGroup, k_LaserTTL, laserTtlData->getTuplePointer(offset));
      if(err < 0)
      {
        setErrorCondition(k_HDF5DatasetError, "Failed to open HDF5 dataset LaserTTL");
        return;
      }

      // Read the XY coordinates
      err = H5Lite::readPointerDataset(sliceGroup, k_XAxis, xCoord.data());
      if(err < 0)
      {
        setErrorCondition(k_HDF5DatasetError, "Failed to open HDF5 dataset X-Axis");
        return;
      }

      err = H5Lite::readPointerDataset(sliceGroup, k_YAxis, yCoord.data());
      if(err < 0)
      {
        setErrorCondition(k_HDF5DatasetError, "Failed to open HDF5 dataset Y-Axis");
        return;
      }

      std::string partStartTimeString;
      err = H5Lite::readStringAttribute(sliceGroup, k_PartStartTime, partStartTimeString);
      if(err < 0)
      {
        setErrorCondition(k_HDF5AttributeError, "Failed to open HDF5 attribute PartStartTime");
        return;
      }

      std::string partEndTimeString;
      err = H5Lite::readStringAttribute(sliceGroup, k_PartEndTime, partEndTimeString);
      if(err < 0)
      {
        setErrorCondition(k_HDF5AttributeError, "Failed to open HDF5 attribute PartEndTime");
        return;
      }

      QDateTime partStartTime = QDateTime::fromString(QString::fromStdString(partStartTimeString), Qt::DateFormat::ISODateWithMs);
      QDateTime partEndTime = QDateTime::fromString(QString::fromStdString(partEndTimeString), Qt::DateFormat::ISODateWithMs);

      currentTime = static_cast<double>(initialTime.msecsTo(partStartTime)) / 1000.0;

      int64_t partDeltaTime = partStartTime.msecsTo(partEndTime);

      double deltaTime = std::nearbyint(static_cast<double>(partDeltaTime) / static_cast<double>(sliceElements) * 1000.0) / 1e6;

      // Now fill in the appropriate parts of the sliceData array and vertex array
      for(size_t i = offset; i < offset + sliceElements; i++)
      {
        (*sliceData)[i] = static_cast<int16_t>(slice);
        (*timeData)[i] = currentTime;

        std::array<float, 3> coord = {static_cast<float>(xCoord[i - offset]), static_cast<float>(yCoord[i - offset]), cummulativeLayerThickness};
        vertGeom->setCoords(i, coord.data());

        currentTime += deltaTime;
      }

      offset += sliceElements;
    }
  }
}

// -----------------------------------------------------------------------------
AbstractFilter::Pointer ImportQMMeltpoolH5File::newFilterInstance(bool copyFilterParameters) const
{
  ImportQMMeltpoolH5File::Pointer filter = ImportQMMeltpoolH5File::New();
  if(copyFilterParameters)
  {
    copyFilterParameterInstanceVariables(filter.get());
  }
  return filter;
}

// -----------------------------------------------------------------------------
QString ImportQMMeltpoolH5File::getCompiledLibraryName() const
{
  return DREAM3DReviewConstants::DREAM3DReviewBaseName;
}

// -----------------------------------------------------------------------------
QString ImportQMMeltpoolH5File::getBrandingString() const
{
  return DREAM3DReviewConstants::DREAM3DReviewBaseName;
}

// -----------------------------------------------------------------------------
QString ImportQMMeltpoolH5File::getFilterVersion() const
{
  QString version;
  QTextStream vStream(&version);
  vStream << DREAM3DReview::Version::Major() << "." << DREAM3DReview::Version::Minor() << "." << DREAM3DReview::Version::Patch();
  return version;
}

// -----------------------------------------------------------------------------
QString ImportQMMeltpoolH5File::getGroupName() const
{
  return SIMPL::FilterGroups::IOFilters;
}

// -----------------------------------------------------------------------------
QString ImportQMMeltpoolH5File::getSubGroupName() const
{
  return SIMPL::FilterSubGroups::InputFilters;
}

// -----------------------------------------------------------------------------
QString ImportQMMeltpoolH5File::getHumanLabel() const
{
  return "Import QM Meltpool HDF5 File";
}

// -----------------------------------------------------------------------------
QUuid ImportQMMeltpoolH5File::getUuid() const
{
  return "{14f85d76-2400-57b8-9650-438563a8b8eb}";
}

// -----------------------------------------------------------------------------
ImportQMMeltpoolH5File::Pointer ImportQMMeltpoolH5File::NullPointer()
{
  return nullptr;
}

// -----------------------------------------------------------------------------
std::shared_ptr<ImportQMMeltpoolH5File> ImportQMMeltpoolH5File::New()
{
  struct make_shared_enabler : public ImportQMMeltpoolH5File
  {
  };
  std::shared_ptr<make_shared_enabler> val = std::make_shared<make_shared_enabler>();
  val->setupFilterParameters();
  return val;
}

// -----------------------------------------------------------------------------
QString ImportQMMeltpoolH5File::getNameOfClass() const
{
  return ClassName();
}

// -----------------------------------------------------------------------------
QString ImportQMMeltpoolH5File::ClassName()
{
  return "ImportQMMeltpoolH5File";
}

// -----------------------------------------------------------------------------
void ImportQMMeltpoolH5File::setInputFiles(const VectString& value)
{
  m_InputFiles = value;
}

// -----------------------------------------------------------------------------
ImportQMMeltpoolH5File::VectString ImportQMMeltpoolH5File::getInputFiles() const
{
  return m_InputFiles;
}

//// -----------------------------------------------------------------------------
// void ImportQMMeltpoolH5File::setHDF5FilePath(const QString& value)
//{
//  m_HDF5FilePath = value;
//}

//// -----------------------------------------------------------------------------
// QString ImportQMMeltpoolH5File::getHDF5FilePath() const
//{
//  return m_HDF5FilePath;
//}

// -----------------------------------------------------------------------------
void ImportQMMeltpoolH5File::setDataContainerPath(const DataArrayPath& value)
{
  m_DataContainerPath = value;
}

// -----------------------------------------------------------------------------
DataArrayPath ImportQMMeltpoolH5File::getDataContainerPath() const
{
  return m_DataContainerPath;
}

// -----------------------------------------------------------------------------
void ImportQMMeltpoolH5File::setVertexAttributeMatrixName(const QString& value)
{
  m_VertexAttributeMatrixName = value;
}

// -----------------------------------------------------------------------------
QString ImportQMMeltpoolH5File::getVertexAttributeMatrixName() const
{
  return m_VertexAttributeMatrixName;
}

// -----------------------------------------------------------------------------
void ImportQMMeltpoolH5File::setSliceRange(const IntVec2Type& value)
{
  m_SliceRange = value;
}

// -----------------------------------------------------------------------------
IntVec2Type ImportQMMeltpoolH5File::getSliceRange() const
{
  return m_SliceRange;
}

// -----------------------------------------------------------------------------
void ImportQMMeltpoolH5File::setPower(float value)
{
  m_Power = value;
}

// -----------------------------------------------------------------------------
float ImportQMMeltpoolH5File::getPower() const
{
  return m_Power;
}

// -----------------------------------------------------------------------------
QString ImportQMMeltpoolH5File::getPossibleIndices() const
{
  return m_PossibleIndices;
}
