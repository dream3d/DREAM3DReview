/* ============================================================================
 * Copyright (c) 2009-2020 BlueQuartz Software, LLC
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
#include "ApplyTransformationToGeometry.h"

#include <Eigen/Dense>

#include <QtCore/QTextStream>

#include "SIMPLib/Common/Constants.h"
#include "SIMPLib/Common/SIMPLRange.h"
#include "SIMPLib/DataContainers/DataContainerArray.h"
#include "SIMPLib/FilterParameters/AbstractFilterParametersReader.h"
#include "SIMPLib/FilterParameters/DataArraySelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/DataContainerSelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/DynamicTableFilterParameter.h"
#include "SIMPLib/FilterParameters/FloatFilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedChoicesFilterParameter.h"
#include "SIMPLib/Geometry/EdgeGeom.h"
#include "SIMPLib/Geometry/IGeometry2D.h"
#include "SIMPLib/Geometry/IGeometry3D.h"
#include "SIMPLib/Geometry/VertexGeom.h"
#include "SIMPLib/Math/SIMPLibMath.h"
#include "SIMPLib/Utilities/ParallelDataAlgorithm.h"

#include "EbsdLib/Core/Orientation.hpp"
#include "EbsdLib/Core/OrientationTransformation.hpp"

#include "DREAM3DReview/DREAM3DReviewConstants.h"
#include "DREAM3DReview/DREAM3DReviewVersion.h"

namespace ApplyTransformationProgress
{
static size_t s_InstanceIndex = 0;
static std::map<size_t, int64_t> s_ProgressValues;
static std::map<size_t, int64_t> s_LastProgressInt;
} // namespace ApplyTransformationProgress

class ApplyTransformationToGeometryImpl
{

public:
  ApplyTransformationToGeometryImpl(ApplyTransformationToGeometry& filter, float* transformationMatrix, const SharedVertexList::Pointer verticesPtr)
  : m_Filter(filter)
  , m_TransformationMatrix(transformationMatrix)
  , m_Vertices(verticesPtr)
  {
  }
  ~ApplyTransformationToGeometryImpl() = default;

  void convert(size_t start, size_t end) const
  {
    using ProjectiveMatrix = Eigen::Matrix<float, 4, 4, Eigen::RowMajor>;
    Eigen::Map<ProjectiveMatrix> transformation(m_TransformationMatrix);

    int64_t progCounter = 0;
    int64_t totalElements = (end - start);
    int64_t progIncrement = static_cast<int64_t>(totalElements / 100);

    SharedVertexList& vertices = *(m_Vertices.get());
    for(size_t i = start; i < end; i++)
    {
      if(m_Filter.getCancel())
      {
        return;
      }
      Eigen::Vector4f position(vertices[3 * i + 0], vertices[3 * i + 1], vertices[3 * i + 2], 1);
      Eigen::Vector4f transformedPosition = transformation * position;
      vertices.setTuple(i, transformedPosition.data());

      if(progCounter > progIncrement)
      {
        m_Filter.sendThreadSafeProgressMessage(progCounter);
        progCounter = 0;
      }
      progCounter++;
    }
  }

  void operator()(const SIMPLRange& range) const
  {
    convert(range.min(), range.max());
  }

private:
  ApplyTransformationToGeometry& m_Filter;
  float* m_TransformationMatrix = nullptr;
  SharedVertexList::Pointer m_Vertices;
};

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
ApplyTransformationToGeometry::ApplyTransformationToGeometry()
: m_ComputedTransformationMatrix("", "", "TransformationMatrix")
, m_GeometryToTransform("")
, m_TransformationMatrixType(1)
{
  m_RotationAngle = 0.0f;
  m_RotationAxis[0] = 0.0f;
  m_RotationAxis[1] = 0.0f;
  m_RotationAxis[2] = 1.0f;

  m_Translation[0] = 0.0f;
  m_Translation[1] = 0.0f;
  m_Translation[2] = 0.0f;

  m_Scale[0] = 0.0f;
  m_Scale[1] = 0.0f;
  m_Scale[2] = 0.0f;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
ApplyTransformationToGeometry::~ApplyTransformationToGeometry() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ApplyTransformationToGeometry::setupFilterParameters()
{
  FilterParameterVectorType parameters;
  {
    LinkedChoicesFilterParameter::Pointer parameter = LinkedChoicesFilterParameter::New();
    parameter->setHumanLabel("Transformation Type");
    parameter->setPropertyName("TransformationMatrixType");
    parameter->setSetterCallback(SIMPL_BIND_SETTER(ApplyTransformationToGeometry, this, TransformationMatrixType));
    parameter->setGetterCallback(SIMPL_BIND_GETTER(ApplyTransformationToGeometry, this, TransformationMatrixType));
    QVector<QString> choices;
    choices.push_back("No Transformation");
    choices.push_back("Pre-Computed Transformation Matrix");
    choices.push_back("Manual Transformation Matrix");
    choices.push_back("Rotation");
    choices.push_back("Translation");
    choices.push_back("Scale");
    parameter->setChoices(choices);
    QStringList linkedProps = {"ComputedTransformationMatrix", "ManualTransformationMatrix", "RotationAngle", "RotationAxis", "Translation", "Scale"};
    parameter->setLinkedProperties(linkedProps);
    parameter->setEditable(false);
    parameter->setCategory(FilterParameter::Parameter);
    parameters.push_back(parameter);
  }
  {
    QStringList rHeaders, cHeaders;
    std::vector<std::vector<double>> defaultTable;
    for(size_t i = 0; i < 4; i++)
    {
      std::vector<double> row(4, 0);
      row[i] = 1.0;
      defaultTable.push_back(row);
    }
    m_ManualTransformationMatrix.setTableData(defaultTable);
    parameters.push_back(SIMPL_NEW_DYN_TABLE_FP("Transformation Matrix", ManualTransformationMatrix, FilterParameter::Parameter, ApplyTransformationToGeometry, 2));
  }
  parameters.push_back(SIMPL_NEW_FLOAT_FP("Rotation Angle (Degrees)", RotationAngle, FilterParameter::Parameter, ApplyTransformationToGeometry, 3));
  parameters.push_back(SIMPL_NEW_FLOAT_VEC3_FP("Rotation Axis (ijk)", RotationAxis, FilterParameter::Parameter, ApplyTransformationToGeometry, 3));
  parameters.push_back(SIMPL_NEW_FLOAT_VEC3_FP("Translation", Translation, FilterParameter::Parameter, ApplyTransformationToGeometry, 4));
  parameters.push_back(SIMPL_NEW_FLOAT_VEC3_FP("Scale", Scale, FilterParameter::Parameter, ApplyTransformationToGeometry, 5));
  DataContainerSelectionFilterParameter::RequirementType dcReq;
  IGeometry::Types geomTypes = {IGeometry::Type::Vertex, IGeometry::Type::Edge, IGeometry::Type::Triangle, IGeometry::Type::Quad, IGeometry::Type::Tetrahedral};
  dcReq.dcGeometryTypes = geomTypes;
  parameters.push_back(SIMPL_NEW_DC_SELECTION_FP("Geometry to Transform", GeometryToTransform, FilterParameter::RequiredArray, ApplyTransformationToGeometry, dcReq));
  {
    DataArraySelectionFilterParameter::RequirementType dasReq =
        DataArraySelectionFilterParameter::CreateRequirement(SIMPL::TypeNames::Float, SIMPL::Defaults::AnyComponentSize, AttributeMatrix::Type::Generic, IGeometry::Type::Any);
    parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Transformation Matrix", ComputedTransformationMatrix, FilterParameter::RequiredArray, ApplyTransformationToGeometry, dasReq, 1));
  }
  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ApplyTransformationToGeometry::readFilterParameters(AbstractFilterParametersReader* reader, int index)
{
  reader->openFilterGroup(this, index);
  setManualTransformationMatrix(reader->readDynamicTableData("ManualTransformationMatrix", getManualTransformationMatrix()));
  setComputedTransformationMatrix(reader->readDataArrayPath("ComputedTransformationMatrix", getComputedTransformationMatrix()));
  setGeometryToTransform(reader->readDataArrayPath("GeometryToTransform", getGeometryToTransform()));
  setTransformationMatrixType(reader->readValue("TransformationMatrixType", getTransformationMatrixType()));
  setRotationAxis(reader->readFloatVec3("RotationAxis", getRotationAxis()));
  setRotationAngle(reader->readValue("RotationAngle", getRotationAngle()));
  setTranslation(reader->readFloatVec3("Translation", getTranslation()));
  setScale(reader->readFloatVec3("Scale", getScale()));
  reader->closeFilterGroup();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ApplyTransformationToGeometry::dataCheck()
{
  clearErrorCode();
  clearWarningCode();

  IGeometry::Pointer igeom = getDataContainerArray()->getPrereqGeometryFromDataContainer<IGeometry>(this, getGeometryToTransform());

  if(getErrorCode() < 0)
  {
    return;
  }

  if(!std::dynamic_pointer_cast<IGeometry2D>(igeom) && !std::dynamic_pointer_cast<IGeometry3D>(igeom) && !std::dynamic_pointer_cast<VertexGeom>(igeom) && !std::dynamic_pointer_cast<EdgeGeom>(igeom))
  {
    QString ss =
        QObject::tr("Geometry to transform must be an unstructured geometry (Vertex, Edge, Triangle, Quadrilateral, or Tetrahedral), but the type is %1").arg(igeom->getGeometryTypeAsString());
    setErrorCondition(-702, ss);
  }

  std::vector<size_t> cDims = {4, 4};

  switch(getTransformationMatrixType())
  {
  case 0: // No-Op
  {
    QString ss = QObject::tr("No transformation has been selected, so this filter will perform no operations");
    setWarningCondition(-701, ss);
  }
  case 1: // Transformation matrix from array
  {
    m_TransformationMatrixPtr = getDataContainerArray()->getPrereqArrayFromPath<DataArray<float>>(this, getComputedTransformationMatrix(), cDims);
    if(m_TransformationMatrixPtr.lock())
    {
      m_TransformationMatrix = m_TransformationMatrixPtr.lock()->getPointer(0);
    }
    break;
  }
  case 2: // Manual transformation matrix
  {
    if(getManualTransformationMatrix().getNumRows() != 4)
    {
      QString ss = QObject::tr("Manually entered transformation matrix must have exactly 4 rows");
      setErrorCondition(-702, ss);
      return;
    }
    if(getManualTransformationMatrix().getNumCols() != 4)
    {
      QString ss = QObject::tr("Manually entered transformation matrix must have exactly 4 columns");
      setErrorCondition(-703, ss);
      return;
    }
    std::vector<std::vector<double>> tableData = getManualTransformationMatrix().getTableData();
    m_TransformationReference = FloatArrayType::CreateArray(1, cDims, "_INTERNAL_USE_ONLY_ManualTransformationMatrix", true);
    m_TransformationReference->initializeWithZeros();
    m_TransformationMatrixPtr = m_TransformationReference;
    if(m_TransformationMatrixPtr.lock())
    {
      m_TransformationMatrix = m_TransformationMatrixPtr.lock()->getPointer(0);
      for(size_t i = 0; i < tableData.size(); i++)
      {
        std::vector<double> row = tableData[i];
        for(size_t j = 0; j < row.size(); j++)
        {
          m_TransformationMatrix[4 * i + j] = static_cast<float>(row[j]);
        }
      }
    }
    break;
  }
  case 3: // Rotation via axis-angle
  {
    float rotAngle = m_RotationAngle * SIMPLib::Constants::k_PiOver180D;
    OrientationF om = OrientationTransformation::ax2om<OrientationF, OrientationF>(OrientationF(m_RotationAxis[0], m_RotationAxis[1], m_RotationAxis[2], rotAngle));

    m_TransformationReference = FloatArrayType::CreateArray(1, cDims, "_INTERNAL_USE_ONLY_ManualTransformationMatrix", true);
    m_TransformationReference->initializeWithZeros();
    m_TransformationMatrixPtr = m_TransformationReference;
    if(m_TransformationMatrixPtr.lock())
    {
      m_TransformationMatrix = m_TransformationMatrixPtr.lock()->getPointer(0);
      for(size_t i = 0; i < 3; i++)
      {
        m_TransformationMatrix[4 * i + 0] = om[3 * i + 0];
        m_TransformationMatrix[4 * i + 1] = om[3 * i + 1];
        m_TransformationMatrix[4 * i + 2] = om[3 * i + 2];
        m_TransformationMatrix[4 * i + 3] = 0.0f;
      }
      m_TransformationMatrix[4 * 3 + 3] = 1.0f;
    }
    break;
  }
  case 4: // Translation
  {
    m_TransformationReference = FloatArrayType::CreateArray(1, cDims, "_INTERNAL_USE_ONLY_ManualTransformationMatrix", true);
    m_TransformationReference->initializeWithZeros();
    m_TransformationMatrixPtr = m_TransformationReference;
    if(m_TransformationMatrixPtr.lock())
    {
      m_TransformationMatrix = m_TransformationMatrixPtr.lock()->getPointer(0);
      m_TransformationMatrix[4 * 0 + 0] = 1.0f;
      m_TransformationMatrix[4 * 1 + 1] = 1.0f;
      m_TransformationMatrix[4 * 2 + 2] = 1.0f;
      m_TransformationMatrix[4 * 0 + 3] = m_Translation[0];
      m_TransformationMatrix[4 * 1 + 3] = m_Translation[1];
      m_TransformationMatrix[4 * 2 + 3] = m_Translation[2];
      m_TransformationMatrix[4 * 3 + 3] = 1.0f;
    }
    break;
  }
  case 5: // Scale
  {
    m_TransformationReference = FloatArrayType::CreateArray(1, cDims, "_INTERNAL_USE_ONLY_ManualTransformationMatrix", true);
    m_TransformationReference->initializeWithZeros();
    m_TransformationMatrixPtr = m_TransformationReference;
    if(m_TransformationMatrixPtr.lock())
    {
      m_TransformationMatrix = m_TransformationMatrixPtr.lock()->getPointer(0);
      m_TransformationMatrix[4 * 0 + 0] = m_Scale[0];
      m_TransformationMatrix[4 * 1 + 1] = m_Scale[1];
      m_TransformationMatrix[4 * 2 + 2] = m_Scale[2];
      m_TransformationMatrix[4 * 3 + 3] = 1.0f;
    }
    break;
  }
  default: {
    QString ss = QObject::tr("Invalid selection for transformation type");
    setErrorCondition(-701, ss);
    break;
  }
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ApplyTransformationToGeometry::applyTransformation()
{

  IGeometry::Pointer igeom = getDataContainerArray()->getDataContainer(m_GeometryToTransform)->getGeometry();

  SharedVertexList::Pointer vertexList;

  if(IGeometry2D::Pointer igeom2D = std::dynamic_pointer_cast<IGeometry2D>(igeom))
  {
    vertexList = igeom2D->getVertices();
  }
  else if(IGeometry3D::Pointer igeom3D = std::dynamic_pointer_cast<IGeometry3D>(igeom))
  {
    vertexList = igeom3D->getVertices();
  }
  else if(VertexGeom::Pointer vertex = std::dynamic_pointer_cast<VertexGeom>(igeom))
  {
    vertexList = vertex->getVertices();
  }
  else if(EdgeGeom::Pointer edge = std::dynamic_pointer_cast<EdgeGeom>(igeom))
  {
    vertexList = edge->getVertices();
  }
  else
  {
    return;
  }

  using ProjectiveMatrix = Eigen::Matrix<float, 4, 4, Eigen::RowMajor>;
  Eigen::Map<ProjectiveMatrix> transformation(m_TransformationMatrix);
  m_TotalElements = vertexList->getNumberOfTuples();
  // Allow data-based parallelization
#if 1
  ParallelDataAlgorithm dataAlg;
  dataAlg.setRange(0, m_TotalElements);
  dataAlg.execute(ApplyTransformationToGeometryImpl(*this, m_TransformationMatrix, vertexList));
#else
  // THis chunk of code will FORCE single threaded mode. Don't do this unless you really mean it.
  ApplyTransformationToGeometryImpl doThis(*this, m_TransformationMatrix, vertexList);
  doThis({0, vertexList->getNumberOfTuples()});
#endif
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ApplyTransformationToGeometry::execute()
{
  dataCheck();
  if(getErrorCode() < 0)
  {
    return;
  }
  if(getWarningCode() < 0)
  {
    return;
  }

  if(m_TransformationMatrixType == 0)
  {
    return;
  }
  // Needed for Threaded Progress Messages
  m_InstanceIndex = ++ApplyTransformationProgress::s_InstanceIndex;
  ApplyTransformationProgress::s_ProgressValues[m_InstanceIndex] = 0;
  ApplyTransformationProgress::s_LastProgressInt[m_InstanceIndex] = 0;

  applyTransformation();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ApplyTransformationToGeometry::sendThreadSafeProgressMessage(int64_t counter)
{
  std::lock_guard<std::mutex> guard(m_ProgressMessage_Mutex);

  int64_t& progCounter = ApplyTransformationProgress::s_ProgressValues[m_InstanceIndex];
  progCounter += counter;
  int64_t progressInt = static_cast<int64_t>((static_cast<float>(progCounter) / m_TotalElements) * 100.0f);

  int64_t progIncrement = m_TotalElements / 100;
  int64_t prog = 1;

  int64_t& lastProgressInt = ApplyTransformationProgress::s_LastProgressInt[m_InstanceIndex];

  if(progCounter > prog && lastProgressInt != progressInt)
  {
    QString ss = QObject::tr("Transforming || %1% Completed").arg(progressInt);
    notifyStatusMessage(ss);
    prog += progIncrement;
  }

  lastProgressInt = progressInt;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AbstractFilter::Pointer ApplyTransformationToGeometry::newFilterInstance(bool copyFilterParameters) const
{
  ApplyTransformationToGeometry::Pointer filter = ApplyTransformationToGeometry::New();
  if(copyFilterParameters)
  {
    copyFilterParameterInstanceVariables(filter.get());
  }
  return filter;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ApplyTransformationToGeometry::getCompiledLibraryName() const
{
  return DREAM3DReviewConstants::DREAM3DReviewBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ApplyTransformationToGeometry::getBrandingString() const
{
  return "DREAM3DReview";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ApplyTransformationToGeometry::getFilterVersion() const
{
  QString version;
  QTextStream vStream(&version);
  vStream << DREAM3DReview::Version::Major() << "." << DREAM3DReview::Version::Minor() << "." << DREAM3DReview::Version::Patch();
  return version;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ApplyTransformationToGeometry::getGroupName() const
{
  return DREAM3DReviewConstants::FilterGroups::DREAM3DReviewFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QUuid ApplyTransformationToGeometry::getUuid() const
{
  return QUuid("{c681caf4-22f2-5885-bbc9-a0476abc72eb}");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ApplyTransformationToGeometry::getSubGroupName() const
{
  return DREAM3DReviewConstants::FilterSubGroups::RotationTransformationFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ApplyTransformationToGeometry::getHumanLabel() const
{
  return "Apply Transformation to Geometry";
}

// -----------------------------------------------------------------------------
ApplyTransformationToGeometry::Pointer ApplyTransformationToGeometry::NullPointer()
{
  return Pointer(static_cast<Self*>(nullptr));
}

// -----------------------------------------------------------------------------
std::shared_ptr<ApplyTransformationToGeometry> ApplyTransformationToGeometry::New()
{
  struct make_shared_enabler : public ApplyTransformationToGeometry
  {
  };
  std::shared_ptr<make_shared_enabler> val = std::make_shared<make_shared_enabler>();
  val->setupFilterParameters();
  return val;
}

// -----------------------------------------------------------------------------
QString ApplyTransformationToGeometry::getNameOfClass() const
{
  return QString("ApplyTransformationToGeometry");
}

// -----------------------------------------------------------------------------
QString ApplyTransformationToGeometry::ClassName()
{
  return QString("ApplyTransformationToGeometry");
}

// -----------------------------------------------------------------------------
void ApplyTransformationToGeometry::setManualTransformationMatrix(const DynamicTableData& value)
{
  m_ManualTransformationMatrix = value;
}

// -----------------------------------------------------------------------------
DynamicTableData ApplyTransformationToGeometry::getManualTransformationMatrix() const
{
  return m_ManualTransformationMatrix;
}

// -----------------------------------------------------------------------------
void ApplyTransformationToGeometry::setComputedTransformationMatrix(const DataArrayPath& value)
{
  m_ComputedTransformationMatrix = value;
}

// -----------------------------------------------------------------------------
DataArrayPath ApplyTransformationToGeometry::getComputedTransformationMatrix() const
{
  return m_ComputedTransformationMatrix;
}

// -----------------------------------------------------------------------------
void ApplyTransformationToGeometry::setGeometryToTransform(const DataArrayPath& value)
{
  m_GeometryToTransform = value;
}

// -----------------------------------------------------------------------------
DataArrayPath ApplyTransformationToGeometry::getGeometryToTransform() const
{
  return m_GeometryToTransform;
}

// -----------------------------------------------------------------------------
void ApplyTransformationToGeometry::setTransformationMatrixType(int value)
{
  m_TransformationMatrixType = value;
}

// -----------------------------------------------------------------------------
int ApplyTransformationToGeometry::getTransformationMatrixType() const
{
  return m_TransformationMatrixType;
}

// -----------------------------------------------------------------------------
void ApplyTransformationToGeometry::setRotationAxis(const FloatVec3Type& value)
{
  m_RotationAxis = value;
}

// -----------------------------------------------------------------------------
FloatVec3Type ApplyTransformationToGeometry::getRotationAxis() const
{
  return m_RotationAxis;
}

// -----------------------------------------------------------------------------
void ApplyTransformationToGeometry::setRotationAngle(float value)
{
  m_RotationAngle = value;
}

// -----------------------------------------------------------------------------
float ApplyTransformationToGeometry::getRotationAngle() const
{
  return m_RotationAngle;
}

// -----------------------------------------------------------------------------
void ApplyTransformationToGeometry::setTranslation(const FloatVec3Type& value)
{
  m_Translation = value;
}

// -----------------------------------------------------------------------------
FloatVec3Type ApplyTransformationToGeometry::getTranslation() const
{
  return m_Translation;
}

// -----------------------------------------------------------------------------
void ApplyTransformationToGeometry::setScale(const FloatVec3Type& value)
{
  m_Scale = value;
}

// -----------------------------------------------------------------------------
FloatVec3Type ApplyTransformationToGeometry::getScale() const
{
  return m_Scale;
}
