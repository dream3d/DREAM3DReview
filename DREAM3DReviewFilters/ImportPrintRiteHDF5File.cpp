#include "ImportPrintRiteHDF5File.h"

#include <QtCore/QFileInfo>

#include "H5Support/QH5Lite.h"

#include "H5Support/H5ScopedSentinel.h"
#include "H5Support/QH5Utilities.h"
#include "SIMPLib/Common/Constants.h"
#include "SIMPLib/Common/TemplateHelpers.h"
#include "SIMPLib/DataArrays/NeighborList.hpp"
#include "SIMPLib/DataContainers/DataContainer.h"
#include "SIMPLib/DataContainers/DataContainerArray.h"
#include "SIMPLib/FilterParameters/AbstractFilterParametersReader.h"
#include "SIMPLib/FilterParameters/InputFileFilterParameter.h"
#include "SIMPLib/FilterParameters/SeparatorFilterParameter.h"
#include "SIMPLib/FilterParameters/StringFilterParameter.h"
#include "SIMPLib/Geometry/VertexGeom.h"
#include "SIMPLib/HDF5/H5DataArrayReader.h"
#include "SIMPLib/Math/SIMPLibMath.h"

#include "DREAM3DReview/DREAM3DReviewConstants.h"
#include "DREAM3DReview/DREAM3DReviewVersion.h"
#include "DREAM3DReviewFilters/util/PrintRiteHelpers.h"

namespace PRH = PrintRiteHelpers;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
ImportPrintRiteHDF5File::ImportPrintRiteHDF5File() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
ImportPrintRiteHDF5File::~ImportPrintRiteHDF5File() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ImportPrintRiteHDF5File::initialize()
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ImportPrintRiteHDF5File::setupFilterParameters()
{
  FilterParameterVectorType parameters;
  parameters.push_back(SIMPL_NEW_INPUT_FILE_FP("Input TDMS-HDF5 File", InputFile, FilterParameter::Category::Parameter, ImportPrintRiteHDF5File, "*.hdf5", "TDMS-HDF5"));
  parameters.push_back(SIMPL_NEW_STRING_FP("High Frequency Data Container", HFDataContainerName, FilterParameter::Category::CreatedArray, ImportPrintRiteHDF5File));
  parameters.push_back(SeparatorFilterParameter::Create("Vertex Data", FilterParameter::Category::CreatedArray));
  parameters.push_back(SIMPL_NEW_STRING_FP("High Frequency Vertex Attribute Matrix", HFDataName, FilterParameter::Category::CreatedArray, ImportPrintRiteHDF5File));
  parameters.push_back(SIMPL_NEW_STRING_FP("High Frequency Slice Ids", HFSliceIdsArrayName, FilterParameter::Category::CreatedArray, ImportPrintRiteHDF5File));
  parameters.push_back(SeparatorFilterParameter::Create("Vertex Feature Data", FilterParameter::Category::CreatedArray));
  parameters.push_back(SIMPL_NEW_STRING_FP("High Frequency Slice Attribute Matrix", HFSliceDataName, FilterParameter::Category::CreatedArray, ImportPrintRiteHDF5File));
  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ImportPrintRiteHDF5File::readFilterParameters(AbstractFilterParametersReader* reader, int index)
{
  reader->openFilterGroup(this, index);
  setInputFile(reader->readString("InputFile", getInputFile()));
  setHFDataContainerName(reader->readString("HFDataContainerName", getHFDataContainerName()));
  setHFDataName(reader->readString("HFDataName", getHFDataName()));
  setHFSliceDataName(reader->readString("HFSliceDataName", getHFSliceDataName()));
  setHFSliceIdsArrayName(reader->readString("HFSliceIdsArrayName", getHFSliceIdsArrayName()));
  reader->closeFilterGroup();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ImportPrintRiteHDF5File::readArray(const DataArrayPath& path, const hid_t& gid, QString& name, size_t index)
{
  herr_t err = 0;
  IDataArray::Pointer tmpIPtr = getDataContainerArray()->getPrereqIDataArrayFromPath(this, path);
  if(TemplateHelpers::CanDynamicCast<BoolArrayType>()(tmpIPtr))
  {
    DataArray<bool>::Pointer tmpPtr = std::dynamic_pointer_cast<DataArray<bool>>(tmpIPtr);
    bool* tmpRawPtr = tmpPtr->getPointer(index);
    size_t size = QH5Lite::getNumberOfElements(gid, name);
    if(index + size <= tmpPtr->getNumberOfTuples())
    {
      err = QH5Lite::readPointerDataset<bool>(gid, name, tmpRawPtr);
    }
    else
      err = -33;
  }
  else if(TemplateHelpers::CanDynamicCast<UInt8ArrayType>()(tmpIPtr))
  {
    DataArray<uint8_t>::Pointer tmpPtr = std::dynamic_pointer_cast<DataArray<uint8_t>>(tmpIPtr);
    uint8_t* tmpRawPtr = tmpPtr->getPointer(index);
    size_t size = QH5Lite::getNumberOfElements(gid, name);
    if(index + size <= tmpPtr->getNumberOfTuples())
    {
      err = QH5Lite::readPointerDataset<uint8_t>(gid, name, tmpRawPtr);
    }
    else
      err = -33;
  }
  else if(TemplateHelpers::CanDynamicCast<Int8ArrayType>()(tmpIPtr))
  {
    DataArray<int8_t>::Pointer tmpPtr = std::dynamic_pointer_cast<DataArray<int8_t>>(tmpIPtr);
    int8_t* tmpRawPtr = tmpPtr->getPointer(index);
    size_t size = QH5Lite::getNumberOfElements(gid, name);
    if(index + size <= tmpPtr->getNumberOfTuples())
    {
      err = QH5Lite::readPointerDataset<int8_t>(gid, name, tmpRawPtr);
    }
    else
      err = -33;
  }
  else if(TemplateHelpers::CanDynamicCast<UInt16ArrayType>()(tmpIPtr))
  {
    DataArray<uint16_t>::Pointer tmpPtr = std::dynamic_pointer_cast<DataArray<uint16_t>>(tmpIPtr);
    uint16_t* tmpRawPtr = tmpPtr->getPointer(index);
    size_t size = QH5Lite::getNumberOfElements(gid, name);
    if(index + size <= tmpPtr->getNumberOfTuples())
    {
      err = QH5Lite::readPointerDataset<uint16_t>(gid, name, tmpRawPtr);
    }
    else
      err = -33;
  }
  else if(TemplateHelpers::CanDynamicCast<Int16ArrayType>()(tmpIPtr))
  {
    DataArray<int16_t>::Pointer tmpPtr = std::dynamic_pointer_cast<DataArray<int16_t>>(tmpIPtr);
    int16_t* tmpRawPtr = tmpPtr->getPointer(index);
    size_t size = QH5Lite::getNumberOfElements(gid, name);
    if(index + size <= tmpPtr->getNumberOfTuples())
    {
      err = QH5Lite::readPointerDataset<int16_t>(gid, name, tmpRawPtr);
    }
    else
      err = -33;
  }
  else if(TemplateHelpers::CanDynamicCast<UInt32ArrayType>()(tmpIPtr))
  {
    DataArray<uint32_t>::Pointer tmpPtr = std::dynamic_pointer_cast<DataArray<uint32_t>>(tmpIPtr);
    uint32_t* tmpRawPtr = tmpPtr->getPointer(index);
    size_t size = QH5Lite::getNumberOfElements(gid, name);
    if(index + size <= tmpPtr->getNumberOfTuples())
    {
      err = QH5Lite::readPointerDataset<uint32_t>(gid, name, tmpRawPtr);
    }
    else
      err = -33;
  }
  else if(TemplateHelpers::CanDynamicCast<Int32ArrayType>()(tmpIPtr))
  {
    Int32ArrayType::Pointer tmpPtr = std::dynamic_pointer_cast<Int32ArrayType>(tmpIPtr);
    int32_t* tmpRawPtr = tmpPtr->getPointer(index);
    size_t size = QH5Lite::getNumberOfElements(gid, name);
    if(index + size <= tmpPtr->getNumberOfTuples())
    {
      err = QH5Lite::readPointerDataset<int32_t>(gid, name, tmpRawPtr);
    }
    else
      err = -33;
  }
  else if(TemplateHelpers::CanDynamicCast<UInt64ArrayType>()(tmpIPtr))
  {
    DataArray<uint64_t>::Pointer tmpPtr = std::dynamic_pointer_cast<DataArray<uint64_t>>(tmpIPtr);
    uint64_t* tmpRawPtr = tmpPtr->getPointer(index);
    size_t size = QH5Lite::getNumberOfElements(gid, name);
    if(index + size <= tmpPtr->getNumberOfTuples())
    {
      err = QH5Lite::readPointerDataset<uint64_t>(gid, name, tmpRawPtr);
    }
    else
      err = -33;
  }
  else if(TemplateHelpers::CanDynamicCast<Int64ArrayType>()(tmpIPtr))
  {
    DataArray<int64_t>::Pointer tmpPtr = std::dynamic_pointer_cast<DataArray<int64_t>>(tmpIPtr);
    int64_t* tmpRawPtr = tmpPtr->getPointer(index);
    size_t size = QH5Lite::getNumberOfElements(gid, name);
    if(index + size <= tmpPtr->getNumberOfTuples())
    {
      err = QH5Lite::readPointerDataset<int64_t>(gid, name, tmpRawPtr);
    }
    else
      err = -33;
  }
  else if(TemplateHelpers::CanDynamicCast<FloatArrayType>()(tmpIPtr))
  {
    FloatArrayType::Pointer tmpPtr = std::dynamic_pointer_cast<FloatArrayType>(tmpIPtr);
    float* tmpRawPtr = tmpPtr->getPointer(index);
    size_t size = QH5Lite::getNumberOfElements(gid, name);
    if(index + size <= tmpPtr->getNumberOfTuples())
    {
      err = QH5Lite::readPointerDataset<float>(gid, name, tmpRawPtr);
    }
    else
      err = -33;
  }
  else if(TemplateHelpers::CanDynamicCast<DoubleArrayType>()(tmpIPtr))
  {
    DoubleArrayType::Pointer tmpPtr = std::dynamic_pointer_cast<DoubleArrayType>(tmpIPtr);
    double* tmpRawPtr = tmpPtr->getPointer(index);
    size_t size = QH5Lite::getNumberOfElements(gid, name);
    if(index + size <= tmpPtr->getNumberOfTuples())
    {
      err = QH5Lite::readPointerDataset<double>(gid, name, tmpRawPtr);
    }
    else
      err = -33;
  }
  if(err == -33)
  {
    QString ss = QObject::tr("Error reading data set with name %1; data set would not fit in supplied array").arg(name);
    setErrorCondition(-388, ss);
    return;
  }
  else if(err < 0)
  {
    QString ss = QObject::tr("Error reading data set with name %1.").arg(name);
    setErrorCondition(-389, ss);
    return;
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ImportPrintRiteHDF5File::dataCheck()
{
  clearErrorCode();
  clearWarningCode();

  m_SliceGroups.clear();
  m_SortedSliceGroups.clear();
  m_HFDataSetNames.clear();
  m_HFDataSetTypes.clear();
  m_LFDataSetNames.clear();
  m_LFDataSetTypes.clear();
  m_SliceDataSetNames.clear();
  m_NumHFVertsPerSlice.clear();
  m_NumLFVertsPerSlice.clear();
  m_CumulativeNumHFVerts.clear();

  int32_t err = -1;

  QFileInfo fi(getInputFile());

  if(getInputFile().isEmpty() == true)
  {
    QString ss = QObject::tr("The input file must be set");
    setErrorCondition(-387, ss);
  }
  else if(fi.exists() == false)
  {
    QString ss = QObject::tr("The input file does not exist");
    setErrorCondition(-388, ss);
  }

  if(getErrorCode() < 0)
  {
    return;
  }

  hid_t fileId = QH5Utilities::openFile(getInputFile(), true);
  if(fileId < 0)
  {
    QString ss = QObject::tr("Error opening input file");
    setErrorCondition(-388, ss);
    return;
  }

  H5ScopedFileSentinel sentinel(fileId, true);

  size_t totalHFVertices = 0;
  size_t tmpHFVerts = 0;

  hid_t layerGroupId = H5Gopen(fileId, "Layer Data", H5P_DEFAULT);
  if(layerGroupId < 0)
  {
    QString ss = QObject::tr("No group with name 'Layer Data' found in supplied HDF5 file");
    setErrorCondition(-388, ss);
    return;
  }

  err = QH5Utilities::getGroupObjects(layerGroupId, H5Utilities::CustomHDFDataTypes::Group, m_SliceGroups);
  if(m_SliceGroups.size() == 0)
  {
    QString ss = QObject::tr("No layer groups found in supplied HDF5 file");
    setErrorCondition(-388, ss);
    return;
  }

  std::set<int32_t> sortedSliceGroups;
  for(QStringList::Iterator iter = std::begin(m_SliceGroups); iter != std::end(m_SliceGroups); ++iter)
  {
    int32_t layerNum = (*iter).toInt();
    sortedSliceGroups.insert(layerNum);
  }

  // size_t sliceIndex = 0;
  for(auto&& layerNum : sortedSliceGroups)
  {
    QString layerNumAsString = QString::number(layerNum);
    m_SortedSliceGroups.append(layerNumAsString);
  }

  hid_t gid = H5Gopen(layerGroupId, (m_SortedSliceGroups[0]).toStdString().c_str(), H5P_DEFAULT);
  if(gid > 0)
  {
    int32_t layerNum = (m_SortedSliceGroups[0]).toInt();
    m_HFDataSetNames.clear();
    m_LFDataSetNames.clear();
    QString hfData = "High Frequency Data";
    QString lfData = "Low Frequency Data";
    QString unassociatedData = "Unassociated Data";
    hid_t hfGroupId = H5Gopen(gid, hfData.toStdString().c_str(), H5P_DEFAULT);
    if(hfGroupId > 0)
    {
      err = QH5Utilities::getGroupObjects(hfGroupId, H5Utilities::CustomHDFDataTypes::Dataset, m_HFDataSetNames);
      tmpHFVerts = QH5Lite::getNumberOfElements(hfGroupId, m_HFDataSetNames[0]);
      m_NumHFVertsPerSlice[layerNum] = tmpHFVerts;
      totalHFVertices += tmpHFVerts;
    }
    // TODO: once implemented, read lf data
    // hid_t lfGroupId = H5Gopen(gid, lfData.toStdString().c_str(), H5P_DEFAULT);
    // if(lfGroupId > 0)
    //{
    //  err = QH5Utilities::getGroupObjects(lfGroupId, H5Utilities::CustomHDFDataTypes::Dataset, m_LFDataSetNames);
    //  if(m_LFDataSetNames.size() > 0)
    //  {
    //    m_NumLFVertsPerSlice[layerNum] = QH5Lite::getNumberOfElements(lfGroupId, m_LFDataSetNames[0]);
    //  }
    //  else
    //  {
    //    m_NumLFVertsPerSlice[layerNum] = 0;
    //  }
    //}
    QH5Utilities::closeHDF5Object(hfGroupId);
    // QH5Utilities::closeHDF5Object(lfGroupId);
  }
  QH5Utilities::closeHDF5Object(gid);

  for(QStringList::Iterator iter = std::next(m_SortedSliceGroups.begin()); iter != m_SortedSliceGroups.end(); ++iter)
  {
    hid_t gid = H5Gopen(layerGroupId, (*iter).toStdString().c_str(), H5P_DEFAULT);
    if(gid > 0)
    {
      int32_t layerNum = (*iter).toInt();
      QString hfData = "High Frequency Data";
      QString lfData = "Low Frequency Data";
      QString unassociatedData = "Unassociated Data";
      hid_t hfGroupId = H5Gopen(gid, hfData.toStdString().c_str(), H5P_DEFAULT);
      if(hfGroupId > 0)
      {
        tmpHFVerts = QH5Lite::getNumberOfElements(hfGroupId, m_HFDataSetNames[0]);
        m_NumHFVertsPerSlice[layerNum] = tmpHFVerts;
        totalHFVertices += tmpHFVerts;
      }
      // TODO: once implemented, read lf data
      // hid_t lfGroupId = H5Gopen(gid, lfData.toStdString().c_str(), H5P_DEFAULT);
      // if(lfGroupId > 0)
      //{
      //  if(m_LFDataSetNames.size() > 0)
      //  {
      //    m_NumLFVertsPerSlice[layerNum] = QH5Lite::getNumberOfElements(lfGroupId, m_LFDataSetNames[0]);
      //  }
      //  else
      //  {
      //    m_NumLFVertsPerSlice[layerNum] = 0;
      //  }
      //}
      QH5Utilities::closeHDF5Object(hfGroupId);
      // QH5Utilities::closeHDF5Object(lfGroupId);
    }
    QH5Utilities::closeHDF5Object(gid);
  }

  size_t tmpvertcounter = 0;
  for(auto&& pair : m_NumHFVertsPerSlice)
  {
    tmpvertcounter += pair.second;
  }

  m_CumulativeNumHFVerts.resize(m_NumHFVertsPerSlice.size());
  m_CumulativeNumHFVerts[0] = 0;
  size_t index = 1;
  for(auto it = std::next(std::begin(m_NumHFVertsPerSlice)); it != std::end(m_NumHFVertsPerSlice); ++it)
  {
    m_CumulativeNumHFVerts[index] = m_CumulativeNumHFVerts[index - 1] + (*std::prev(it)).second;
    index++;
  }

  DataContainer::Pointer hfDC = getDataContainerArray()->createNonPrereqDataContainer(this, getHFDataContainerName());
  VertexGeom::Pointer hfVertexGeom = VertexGeom::CreateGeometry(static_cast<int64_t>(totalHFVertices), SIMPL::Geometry::VertexGeometry, !getInPreflight());

  if(getErrorCode() < 0)
  {
    return;
  }

  hfDC->setGeometry(hfVertexGeom);

  // TODO: check for build/spatial meta data

  std::vector<size_t> hf_tDims(1, totalHFVertices);
  std::vector<size_t> slice_tDims(1, static_cast<size_t>(m_SortedSliceGroups.size() + 1));
  AttributeMatrix::Pointer hfAttrmat = hfDC->createNonPrereqAttributeMatrix(this, getHFDataName(), hf_tDims, AttributeMatrix::Type::Vertex);
  hfDC->createNonPrereqAttributeMatrix(this, getHFSliceDataName(), slice_tDims, AttributeMatrix::Type::VertexFeature);

  if(getErrorCode() < 0)
  {
    return;
  }

  DataArrayPath hf_path(getHFDataContainerName(), getHFDataName(), "");
  DataArrayPath hfSlice_path(getHFDataContainerName(), getHFSliceDataName(), "");

  // TODO: correct validation for position array existence

  bool haveXCoords = false;
  bool haveYCoords = false;
  // bool haveLayerThickness = false;
  // size_t numElements = getDataContainerArray()->getDataContainer(hf_path.getDataContainerName())->getAttributeMatrix(hf_path.getAttributeMatrixName())->getNumberOfTuples();
  gid = H5Gopen(layerGroupId, m_SortedSliceGroups.front().toStdString().c_str(), H5P_DEFAULT);
  for(QStringList::Iterator iter = m_HFDataSetNames.begin(); iter != m_HFDataSetNames.end(); ++iter)
  {
    QString hfData = "High Frequency Data";
    hid_t hfGroupId = H5Gopen(gid, hfData.toStdString().c_str(), H5P_DEFAULT);
    IDataArray::Pointer ptr = IDataArray::NullPointer();
    ptr = H5DataArrayReader::ReadIDataArray(hfGroupId, *iter, true);
    if(ptr)
    {
      hf_path.setDataArrayName((*iter));
      TemplateHelpers::CreateNonPrereqArrayFromArrayType()(this, hf_path, ptr->getComponentDimensions(), ptr);
    }
    else
    {
      QString ss = QObject::tr("Error reading data set with name '%1'").arg((*iter).toStdString().c_str());
      setErrorCondition(-389, ss);
      return;
    }
    if((*iter).compare("X Position") == 0)
    {
      m_XCoords = FloatArrayType::CreateArray(totalHFVertices, std::string("_INTERNAL_USE_ONLY_XPosition"), !getInPreflight());
      haveXCoords = true;
    }
    else if((*iter).compare("Y Position") == 0)
    {
      m_YCoords = FloatArrayType::CreateArray(totalHFVertices, std::string("_INTERNAL_USE_ONLY_YPosition"), !getInPreflight());
      haveYCoords = true;
    }
  }
  QH5Utilities::closeHDF5Object(gid);

  // TODO: import for lf/unassociated data

  hf_path.setDataArrayName(getHFSliceIdsArrayName());
  std::vector<size_t> cDims(1, 1);
  m_HFSliceIdsPtr = getDataContainerArray()->createNonPrereqArrayFromPath<Int32ArrayType>(this, hf_path, 0, cDims);
  if(nullptr != m_HFSliceIdsPtr.lock().get())
  {
    m_HFSliceIds = m_HFSliceIdsPtr.lock()->getPointer(0);
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ImportPrintRiteHDF5File::execute()
{
  dataCheck();
  if(getErrorCode() < 0)
  {
    return;
  }

  // int32_t err = -1;
  hid_t fileId = QH5Utilities::openFile(getInputFile(), true);
  if(fileId < 0)
  {
    QString ss = QObject::tr("Error opening input file");
    setErrorCondition(-388, ss);
    return;
  }

  H5ScopedFileSentinel sentinel(fileId, true);

  // TODO: the following needs to be checked up in data check

  hid_t buildMetaDataGroupId = H5Gopen(fileId, "Build Meta Data", H5P_DEFAULT);
  IDataArray::Pointer layerThicknessArray = H5DataArrayReader::ReadIDataArray(buildMetaDataGroupId, "Layer Thickness");
  float layerThickness = std::dynamic_pointer_cast<FloatArrayType>(layerThicknessArray)->getValue(0);
  hid_t scalingMetaDataGroupId = H5Gopen(fileId, "Scaling Meta Data", H5P_DEFAULT);
  IDataArray::Pointer spatialCoefficientsArray = H5DataArrayReader::ReadIDataArray(scalingMetaDataGroupId, "Spatial Scaling Coefficients");
  DoubleArrayType::Pointer spatialCoefficients = std::dynamic_pointer_cast<DoubleArrayType>(spatialCoefficientsArray);
  PRH::Polynomial polynomial;
  polynomial.setOrder(3);
  polynomial.setCoefficients(spatialCoefficients);

  VertexGeom::Pointer hf_vertex = getDataContainerArray()->getDataContainer(getHFDataContainerName())->getGeometryAs<VertexGeom>();
  float* hf_vertices = hf_vertex->getVertexPointer(0);
  size_t hf_index = 0;
  DataArrayPath hf_path(getHFDataContainerName(), getHFDataName(), "");
  DataArrayPath hfSlice_path(getHFDataContainerName(), getHFSliceDataName(), "");

  // TODO: import lf/unassociated/meta data
  hid_t layerGroupId = H5Gopen(fileId, "Layer Data", H5P_DEFAULT);
  size_t index = 0;
  for(QStringList::Iterator it = m_SortedSliceGroups.begin(); it != m_SortedSliceGroups.end(); ++it)
  {
    hf_index = m_CumulativeNumHFVerts[index];
    hid_t gid = H5Gopen(layerGroupId, (*it).toStdString().c_str(), H5P_DEFAULT);
    if(gid > 0)
    {
      QString hfData = "High Frequency Data";
      hid_t hfGroupId = H5Gopen(gid, hfData.toStdString().c_str(), H5P_DEFAULT);
      if(hfGroupId > 0 && hf_index < hf_vertex->getNumberOfVertices())
      {
        for(QStringList::Iterator iter = m_HFDataSetNames.begin(); iter != m_HFDataSetNames.end(); ++iter)
        {
          hf_path.setDataArrayName((*iter));

          // TESTING
          std::string tmpit = (*it).toStdString();

          std::string tmpdcname = hf_path.getDataContainerName().toStdString();
          std::string tmpamname = hf_path.getAttributeMatrixName().toStdString();
          std::string tmpdaname = hf_path.getAttributeMatrixName().toStdString();

          std::string tmpiter = (*iter).toStdString();
          // TESTING

          readArray(hf_path, hfGroupId, (*iter), hf_index);
        }
      }
      QH5Utilities::closeHDF5Object(hfGroupId);
    }
    QH5Utilities::closeHDF5Object(gid);
    index++;
  }

  float* xPos = getDataContainerArray()->getDataContainer(getHFDataContainerName())->getAttributeMatrix(getHFDataName())->getAttributeArrayAs<FloatArrayType>("X Position")->getPointer(0);
  float* yPos = getDataContainerArray()->getDataContainer(getHFDataContainerName())->getAttributeMatrix(getHFDataName())->getAttributeArrayAs<FloatArrayType>("Y Position")->getPointer(0);
  int32_t groupIndex = (*m_SortedSliceGroups.begin()).toInt() - 1;
  size_t hf_counter = 0;
  for(size_t i = 0; i < hf_vertex->getNumberOfVertices(); i++)
  {
    float pt[2] = {static_cast<float>(xPos[i]), static_cast<float>(-yPos[i])};
    if(polynomial.nullCoefficients())
    {
      hf_vertices[3 * i + 0] = pt[0];
      hf_vertices[3 * i + 1] = pt[1];
    }
    else
    {
      hf_vertices[3 * i + 0] = polynomial.transformPoint(pt, 0);
      hf_vertices[3 * i + 1] = polynomial.transformPoint(pt, 1);
    }
    hf_vertices[3 * i + 2] = layerThickness * groupIndex;
    m_HFSliceIds[i] = groupIndex + 1;
    hf_counter++;
    if(hf_counter == m_NumHFVertsPerSlice[groupIndex + 1])
    {
      groupIndex++;
      hf_counter = 0;
    }
  }

  // TODO: fix hf time by layer offsets
  // UInt32ArrayType::Pointer hfTimePtr = getDataContainerArray()->getDataContainer(getHFDataContainerName())->getAttributeMatrix(getHFDataName())->getAttributeArrayAs<UInt32ArrayType>("HF Time");
  // uint32_t* hfTime = hfTimePtr->getPointer(0);
  // DoubleArrayType::Pointer layerTimeResPtr =
  // getDataContainerArray()->getDataContainer(getHFDataContainerName())->getAttributeMatrix(getHFSliceDataName())->getAttributeArrayAs<DoubleArrayType>("hf_time_res");
  // double* layerTimeRes = layerTimeResPtr->getPointer(0);
  ////create a float array to replace the hf time array that was an int32
  // DoubleArrayType::Pointer newHFTimePtr = DoubleArrayType::CreateArray(hf_vertex->getNumberOfVertices(), "HF Time", true);
  // double* newHFTime = newHFTimePtr->getPointer(0);
  // for (int64_t i = 0; i < hf_vertex->getNumberOfVertices(); i++)
  //{
  //  newHFTime[i] = (hfTime[i] * layerTimeRes[m_HFSliceIds[i]]) + layerTimeOffset[m_HFSliceIds[i]];
  //}
  // getDataContainerArray()->getDataContainer(getHFDataContainerName())->getAttributeMatrix(getHFDataName())->addAttributeArray("HF Time", newHFTimePtr);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AbstractFilter::Pointer ImportPrintRiteHDF5File::newFilterInstance(bool copyFilterParameters) const
{
  ImportPrintRiteHDF5File::Pointer filter = ImportPrintRiteHDF5File::New();
  if(copyFilterParameters)
  {
    copyFilterParameterInstanceVariables(filter.get());
  }
  return filter;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ImportPrintRiteHDF5File::getCompiledLibraryName() const
{
  return DREAM3DReviewConstants::DREAM3DReviewBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ImportPrintRiteHDF5File::getBrandingString() const
{
  return "DREAM3DReview";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ImportPrintRiteHDF5File::getFilterVersion() const
{
  QString version;
  QTextStream vStream(&version);
  vStream << DREAM3DReview::Version::Major() << "." << DREAM3DReview::Version::Minor() << "." << DREAM3DReview::Version::Patch();
  return version;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ImportPrintRiteHDF5File::getGroupName() const
{
  return SIMPL::FilterGroups::IOFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ImportPrintRiteHDF5File::getSubGroupName() const
{
  return SIMPL::FilterSubGroups::InputFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ImportPrintRiteHDF5File::getHumanLabel() const
{
  return "Import PrintRite HDF5 File";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QUuid ImportPrintRiteHDF5File::getUuid() const
{
  return QUuid("{ab8f6892-2b57-5ec6-88b7-01610d80c32c}");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
ImportPrintRiteHDF5File::Pointer ImportPrintRiteHDF5File::NullPointer()
{
  return Pointer(static_cast<Self*>(nullptr));
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
std::shared_ptr<ImportPrintRiteHDF5File> ImportPrintRiteHDF5File::New()
{
  struct make_shared_enabler : public ImportPrintRiteHDF5File
  {
  };
  std::shared_ptr<make_shared_enabler> val = std::make_shared<make_shared_enabler>();
  val->setupFilterParameters();
  return val;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ImportPrintRiteHDF5File::getNameOfClass() const
{
  return QString("ImportPrintRiteHDF5File");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ImportPrintRiteHDF5File::ClassName()
{
  return QString("ImportPrintRiteHDF5File");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ImportPrintRiteHDF5File::setInputFile(const QString& value)
{
  m_InputFile = value;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ImportPrintRiteHDF5File::getInputFile() const
{
  return m_InputFile;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ImportPrintRiteHDF5File::setHFDataContainerName(const QString& value)
{
  m_HFDataContainerName = value;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ImportPrintRiteHDF5File::getHFDataContainerName() const
{
  return m_HFDataContainerName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ImportPrintRiteHDF5File::setHFDataName(const QString& value)
{
  m_HFDataName = value;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ImportPrintRiteHDF5File::getHFDataName() const
{
  return m_HFDataName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ImportPrintRiteHDF5File::setHFSliceDataName(const QString& value)
{
  m_HFSliceDataName = value;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ImportPrintRiteHDF5File::getHFSliceDataName() const
{
  return m_HFSliceDataName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ImportPrintRiteHDF5File::setHFSliceIdsArrayName(const QString& value)
{
  m_HFSliceIdsArrayName = value;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ImportPrintRiteHDF5File::getHFSliceIdsArrayName() const
{
  return m_HFSliceIdsArrayName;
}
