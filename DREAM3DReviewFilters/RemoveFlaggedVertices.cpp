/* ============================================================================
* Copyright (c) 2009-2016 BlueQuartz Software, LLC
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
* Neither the name of BlueQuartz Software, the US Air Force, nor the names of its
* contributors may be used to endorse or promote products derived from this software
* without specific prior written permission.
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
* The code contained herein was partially funded by the followig contracts:
*    United States Air Force Prime Contract FA8650-07-D-5800
*    United States Air Force Prime Contract FA8650-10-D-5210
*    United States Prime Contract Navy N00173-07-C-2068
*
* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include <memory>

#include "RemoveFlaggedVertices.h"

#include <QtCore/QTextStream>

#include "SIMPLib/Common/Constants.h"

#include "SIMPLib/Common/TemplateHelpers.h"
#include "SIMPLib/FilterParameters/AbstractFilterParametersReader.h"
#include "SIMPLib/FilterParameters/AttributeMatrixSelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/DataArraySelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/DataContainerSelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/StringFilterParameter.h"
#include "SIMPLib/Geometry/VertexGeom.h"
#include "SIMPLib/DataContainers/DataContainerArray.h"
#include "SIMPLib/DataContainers/DataContainer.h"

#include "DREAM3DReview/DREAM3DReviewConstants.h"
#include "DREAM3DReview/DREAM3DReviewVersion.h"

enum createdPathID : RenameDataPath::DataID_t
{
  AttributeMatrixID21 = 21,

  DataContainerID = 1
};

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
RemoveFlaggedVertices::RemoveFlaggedVertices()
: m_VertexGeometry(SIMPL::Defaults::VertexDataContainerName)
, m_MaskArrayPath(SIMPL::Defaults::VertexDataContainerName, SIMPL::Defaults::VertexAttributeMatrixName, SIMPL::CellData::Mask)
, m_ReducedVertexGeometry("ReducedVertexDataContainer")
{
  initialize();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
RemoveFlaggedVertices::~RemoveFlaggedVertices()
= default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void RemoveFlaggedVertices::initialize()
{
  clearErrorCode();
  clearWarningCode();
  setCancel(false);
  m_AttrMatList.clear();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void RemoveFlaggedVertices::setupFilterParameters()
{
  FilterParameterVectorType parameters;
  DataContainerSelectionFilterParameter::RequirementType dcsReq;
  IGeometry::Types geomTypes = {IGeometry::Type::Vertex};
  dcsReq.dcGeometryTypes = geomTypes;
  parameters.push_back(SIMPL_NEW_DC_SELECTION_FP("Vertex Geometry", VertexGeometry, FilterParameter::RequiredArray, RemoveFlaggedVertices, dcsReq));
  DataArraySelectionFilterParameter::RequirementType dasReq = DataArraySelectionFilterParameter::CreateRequirement(SIMPL::TypeNames::Bool, 1, AttributeMatrix::Type::Vertex, IGeometry::Type::Vertex);
  parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Mask", MaskArrayPath, FilterParameter::RequiredArray, RemoveFlaggedVertices, dasReq));
  parameters.push_back(SIMPL_NEW_STRING_FP("Reduced Vertex Data Container", ReducedVertexGeometry, FilterParameter::CreatedArray, RemoveFlaggedVertices));
  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void RemoveFlaggedVertices::readFilterParameters(AbstractFilterParametersReader* reader, int index)
{
  reader->openFilterGroup(this, index);
  setVertexGeometry(reader->readDataArrayPath("VertexGeometry", getVertexGeometry()));
  setMaskArrayPath(reader->readDataArrayPath("MaskArrayPath", getMaskArrayPath()));
  setReducedVertexGeometry(reader->readString("ReducedVertexGeometry", getReducedVertexGeometry()));
  reader->closeFilterGroup();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void RemoveFlaggedVertices::dataCheck()
{
  clearErrorCode();
  clearWarningCode();
  initialize();

  QVector<IDataArray::Pointer> dataArrays;

  VertexGeom::Pointer vertex = getDataContainerArray()->getPrereqGeometryFromDataContainer<VertexGeom>(this, getVertexGeometry());

  if(getErrorCode() < 0)
  {
    return;
  }

  dataArrays.push_back(vertex->getVertices());

  std::vector<size_t> cDims(1, 1);

  m_MaskPtr = getDataContainerArray()->getPrereqArrayFromPath<DataArray<bool>>(this, getMaskArrayPath(), cDims);
  if(m_MaskPtr.lock())
  {
    m_Mask = m_MaskPtr.lock()->getPointer(0);
  }
  if(getErrorCode() >= 0)
  {
    dataArrays.push_back(m_MaskPtr.lock());
  }

  DataContainer::Pointer dc = getDataContainerArray()->createNonPrereqDataContainer(this, getReducedVertexGeometry(), DataContainerID);

  if(getErrorCode() < 0)
  {
    return;
  }

  VertexGeom::Pointer reduced = VertexGeom::CreateGeometry(0, SIMPL::Geometry::VertexGeometry, !getInPreflight());
  dc->setGeometry(reduced);

  getDataContainerArray()->validateNumberOfTuples(this, dataArrays);

  DataContainer::Pointer m = getDataContainerArray()->getDataContainer(getVertexGeometry());
  m_AttrMatList = m->getAttributeMatrixNames();
  std::vector<size_t> tDims(1, 0);
  QList<QString> tempDataArrayList;
  DataArrayPath tempPath;
  AttributeMatrix::Type tempAttrMatType = AttributeMatrix::Type::Vertex;

  if(getErrorCode() < 0)
  {
    return;
  }

  for(auto&& attr_mat : m_AttrMatList)
  {
    AttributeMatrix::Pointer tmpAttrMat = m->getPrereqAttributeMatrix(this, attr_mat, -301);
    if(getErrorCode() >= 0)
    {
      tempAttrMatType = tmpAttrMat->getType();
      if(tempAttrMatType != AttributeMatrix::Type::Vertex)
      {
        AttributeMatrix::Pointer attrMat = tmpAttrMat->deepCopy(getInPreflight());
        dc->addOrReplaceAttributeMatrix(attrMat);
      }
      else
      {
        dc->createNonPrereqAttributeMatrix(this, tmpAttrMat->getName(), tDims, AttributeMatrix::Type::Vertex, AttributeMatrixID21);
        tempDataArrayList = tmpAttrMat->getAttributeArrayNames();
        for(auto&& data_array : tempDataArrayList)
        {
          tempPath.update(getReducedVertexGeometry(), tmpAttrMat->getName(), data_array);
          IDataArray::Pointer tmpDataArray = tmpAttrMat->getPrereqIDataArray(this, data_array, -90002);
          if(getErrorCode() >= 0)
          {
            std::vector<size_t> cDims = tmpDataArray->getComponentDimensions();
            TemplateHelpers::CreateNonPrereqArrayFromArrayType()(this, tempPath, cDims, tmpDataArray);
          }
        }
      }
    }
  }
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
template <typename T> void copyDataToMaskedGeometry(IDataArray::Pointer inDataPtr, IDataArray::Pointer outDataPtr, std::vector<int64_t> maskPoints)
{
  typename DataArray<T>::Pointer inputDataPtr = std::dynamic_pointer_cast<DataArray<T>>(inDataPtr);
  typename DataArray<T>::Pointer maskedDataPtr = std::dynamic_pointer_cast<DataArray<T>>(outDataPtr);
  T* inputData = static_cast<T*>(inputDataPtr->getPointer(0));
  T* maskedData = static_cast<T*>(maskedDataPtr->getPointer(0));

  size_t nComps = inDataPtr->getNumberOfComponents();
  size_t tmpIndex = 0;
  size_t ptrIndex = 0;

  for(std::vector<int64_t>::size_type i = 0; i < maskPoints.size(); i++)
  {
    for(size_t d = 0; d < nComps; d++)
    {
      tmpIndex = nComps * i + d;
      ptrIndex = nComps * maskPoints[i] + d;
      maskedData[tmpIndex] = inputData[ptrIndex];
    }
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void RemoveFlaggedVertices::execute()
{
  initialize();
  dataCheck();
  if(getErrorCode() < 0)
  {
    return;
  }

  VertexGeom::Pointer vertex = getDataContainerArray()->getDataContainer(getVertexGeometry())->getGeometryAs<VertexGeom>();
  int64_t numVerts = vertex->getNumberOfVertices();
  std::vector<int64_t> maskPoints;
  maskPoints.reserve(m_MaskPtr.lock()->getNumberOfTuples());

  for(int64_t i = 0; i < numVerts; i++)
  {
    if(!m_Mask[i])
    {
      maskPoints.push_back(i);
    }
  }

  maskPoints.shrink_to_fit();

  DataContainer::Pointer reduced = getDataContainerArray()->getDataContainer(getReducedVertexGeometry());
  VertexGeom::Pointer reducedVertex = reduced->getGeometryAs<VertexGeom>();
  reducedVertex->resizeVertexList(maskPoints.size());
  float coords[3] = {0.0f, 0.0f, 0.0f};

  for(std::vector<int64_t>::size_type i = 0; i < maskPoints.size(); i++)
  {
    vertex->getCoords(maskPoints[i], coords);
    reducedVertex->setCoords(i, coords);
  }

  DataContainer::Pointer m = getDataContainerArray()->getDataContainer(getVertexGeometry());
  AttributeMatrix::Type tempAttrMatType = AttributeMatrix::Type::Vertex;
  std::vector<size_t> tDims(1, maskPoints.size());

  for(auto&& attr_mat : m_AttrMatList)
  {
    AttributeMatrix::Pointer tmpAttrMat = reduced->getPrereqAttributeMatrix(this, attr_mat, -301);
    if(getErrorCode() >= 0)
    {
      tempAttrMatType = tmpAttrMat->getType();
      if(tempAttrMatType == AttributeMatrix::Type::Vertex)
      {
        tmpAttrMat->resizeAttributeArrays(tDims);
        QList<QString> srcDataArrays = tmpAttrMat->getAttributeArrayNames();
        AttributeMatrix::Pointer srcAttrMat = m->getAttributeMatrix(tmpAttrMat->getName());
        for(auto&& data_array : srcDataArrays)
        {
          IDataArray::Pointer src = srcAttrMat->getAttributeArray(data_array);
          IDataArray::Pointer dest = tmpAttrMat->getAttributeArray(data_array);

          EXECUTE_FUNCTION_TEMPLATE(this, copyDataToMaskedGeometry, src, src, dest, maskPoints)
        }
      }
    }
  }

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AbstractFilter::Pointer RemoveFlaggedVertices::newFilterInstance(bool copyFilterParameters) const
{
  RemoveFlaggedVertices::Pointer filter = RemoveFlaggedVertices::New();
  if(copyFilterParameters)
  {
    copyFilterParameterInstanceVariables(filter.get());
  }
  return filter;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString RemoveFlaggedVertices::getCompiledLibraryName() const
{
  return DREAM3DReviewConstants::DREAM3DReviewBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString RemoveFlaggedVertices::getBrandingString() const
{
  return "DREAM3DReview";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString RemoveFlaggedVertices::getFilterVersion() const
{
  QString version;
  QTextStream vStream(&version);
  vStream << DREAM3DReview::Version::Major() << "." << DREAM3DReview::Version::Minor() << "." << DREAM3DReview::Version::Patch();
  return version;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString RemoveFlaggedVertices::getGroupName() const
{
  return DREAM3DReviewConstants::FilterGroups::DREAM3DReviewFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QUuid RemoveFlaggedVertices::getUuid() const
{
  return QUuid("{379ccc67-16dd-530a-984f-177db2314bac}");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString RemoveFlaggedVertices::getSubGroupName() const
{
  return DREAM3DReviewConstants::FilterSubGroups::GeometryFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString RemoveFlaggedVertices::getHumanLabel() const
{
  return "Remove Flagged Vertices";
}

// -----------------------------------------------------------------------------
RemoveFlaggedVertices::Pointer RemoveFlaggedVertices::NullPointer()
{
  return Pointer(static_cast<Self*>(nullptr));
}

// -----------------------------------------------------------------------------
std::shared_ptr<RemoveFlaggedVertices> RemoveFlaggedVertices::New()
{
  struct make_shared_enabler : public RemoveFlaggedVertices
  {
  };
  std::shared_ptr<make_shared_enabler> val = std::make_shared<make_shared_enabler>();
  val->setupFilterParameters();
  return val;
}

// -----------------------------------------------------------------------------
QString RemoveFlaggedVertices::getNameOfClass() const
{
  return QString("RemoveFlaggedVertices");
}

// -----------------------------------------------------------------------------
QString RemoveFlaggedVertices::ClassName()
{
  return QString("RemoveFlaggedVertices");
}

// -----------------------------------------------------------------------------
void RemoveFlaggedVertices::setVertexGeometry(const DataArrayPath& value)
{
  m_VertexGeometry = value;
}

// -----------------------------------------------------------------------------
DataArrayPath RemoveFlaggedVertices::getVertexGeometry() const
{
  return m_VertexGeometry;
}

// -----------------------------------------------------------------------------
void RemoveFlaggedVertices::setMaskArrayPath(const DataArrayPath& value)
{
  m_MaskArrayPath = value;
}

// -----------------------------------------------------------------------------
DataArrayPath RemoveFlaggedVertices::getMaskArrayPath() const
{
  return m_MaskArrayPath;
}

// -----------------------------------------------------------------------------
void RemoveFlaggedVertices::setReducedVertexGeometry(const QString& value)
{
  m_ReducedVertexGeometry = value;
}

// -----------------------------------------------------------------------------
QString RemoveFlaggedVertices::getReducedVertexGeometry() const
{
  return m_ReducedVertexGeometry;
}
