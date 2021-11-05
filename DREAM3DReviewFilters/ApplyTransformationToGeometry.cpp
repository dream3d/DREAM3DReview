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
 * The code contained herein was partially funded by the following contracts:
 *    United States Air Force Prime Contract FA8650-07-D-5800
 *    United States Air Force Prime Contract FA8650-10-D-5210
 *    United States Prime Contract Navy N00173-07-C-2068
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
#include "ApplyTransformationToGeometry.h"

#ifdef SIMPL_USE_PARALLEL_ALGORITHMS
#include <tbb/blocked_range3d.h>
#include <tbb/parallel_for.h>
#include <tbb/partitioner.h>
#endif

#include <Eigen/Dense>

#include <QtCore/QTextStream>

#include "SIMPLib/Common/Constants.h"
#include "SIMPLib/Common/SIMPLRange.h"
#include "SIMPLib/DataContainers/DataContainerArray.h"
#include "SIMPLib/FilterParameters/AbstractFilterParametersReader.h"
#include "SIMPLib/FilterParameters/AttributeMatrixSelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/DataArraySelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/DataContainerSelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/DynamicTableFilterParameter.h"
#include "SIMPLib/FilterParameters/FloatFilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedChoicesFilterParameter.h"
#include "SIMPLib/FilterParameters/SeparatorFilterParameter.h"
#include "SIMPLib/Geometry/EdgeGeom.h"
#include "SIMPLib/Geometry/IGeometry2D.h"
#include "SIMPLib/Geometry/IGeometry3D.h"
#include "SIMPLib/Geometry/ImageGeom.h"
#include "SIMPLib/Geometry/VertexGeom.h"
#include "SIMPLib/Math/SIMPLibMath.h"
#include "SIMPLib/Utilities/ParallelDataAlgorithm.h"
#include "SIMPLib/Math/MatrixMath.h"

#include "EbsdLib/Core/Orientation.hpp"
#include "EbsdLib/Core/OrientationTransformation.hpp"

#include "DREAM3DReview/DREAM3DReviewConstants.h"
#include "DREAM3DReview/DREAM3DReviewVersion.h"



namespace ApplyTransformationProgress
{

using Matrix3fR = Eigen::Matrix<float, 3, 3, Eigen::RowMajor>;
using Transform3f = Eigen::Transform<float, 3, Eigen::Affine>;
using MatrixTranslation = Eigen::Matrix<float, 1, 3, Eigen::RowMajor>;

struct RotateArgs
    {
  int64_t xp = 0;
  int64_t yp = 0;
  int64_t zp = 0;
  float xRes = 0.0f;
  float yRes = 0.0f;
  float zRes = 0.0f;
  int64_t xpNew = 0;
  int64_t ypNew = 0;
  int64_t zpNew = 0;
  float xResNew = 0.0f;
  float yResNew = 0.0f;
  float zResNew = 0.0f;
  float xMinNew = 0.0f;
  float yMinNew = 0.0f;
  float zMinNew = 0.0f;
    };

const Eigen::Vector3f k_XAxis = Eigen::Vector3f::UnitX();
const Eigen::Vector3f k_YAxis = Eigen::Vector3f::UnitY();
const Eigen::Vector3f k_ZAxis = Eigen::Vector3f::UnitZ();

// Function for determining new ImageGeom dimensions after transformation
void determineMinMax(const Matrix3fR& rotationMatrix, const FloatVec3Type& spacing, size_t col, size_t row, size_t plane, float& xMin, float& xMax, float& yMin, float& yMax, float& zMin, float& zMax)
{
  Eigen::Vector3f coords(static_cast<float>(col) * spacing[0], static_cast<float>(row) * spacing[1], static_cast<float>(plane) * spacing[2]);

  Eigen::Vector3f newCoords = rotationMatrix * coords;



  xMin = std::min(newCoords[0], xMin);
  xMax = std::max(newCoords[0], xMax);

  yMin = std::min(newCoords[1], yMin);
  yMax = std::max(newCoords[1], yMax);

  zMin = std::min(newCoords[2], zMin);
  zMax = std::max(newCoords[2], zMax);
}

float cosBetweenVectors(const Eigen::Vector3f& a, const Eigen::Vector3f& b)
{
  float normA = a.norm();
  float normB = b.norm();

  if(normA == 0.0f || normB == 0.0f)
  {
    return 1.0f;
  }

  return a.dot(b) / (normA * normB);
}

// Function for determining new ImageGeom Spacing between points for scaling
float determineSpacing(const FloatVec3Type& spacing, const Eigen::Vector3f& axisNew)
{
  float xAngle = std::abs(cosBetweenVectors(k_XAxis, axisNew));
  float yAngle = std::abs(cosBetweenVectors(k_YAxis, axisNew));
  float zAngle = std::abs(cosBetweenVectors(k_ZAxis, axisNew));

  std::array<float, 3> axes = {xAngle, yAngle, zAngle};

  auto iter = std::max_element(axes.cbegin(), axes.cend());

  size_t index = std::distance(axes.cbegin(), iter);

  return spacing[index];
}

RotateArgs createRotateParams(const ImageGeom& imageGeom, const Transform3f transformationMatrix)
{
  const SizeVec3Type origDims = imageGeom.getDimensions();
  const FloatVec3Type spacing = imageGeom.getSpacing();
  // const FloatVec3Type origin = imageGeom.getOrigin();

  ApplyTransformationProgress::Matrix3fR rotationMatrix = ApplyTransformationProgress::Matrix3fR::Zero();
  ApplyTransformationProgress::Matrix3fR scaleMatrix = ApplyTransformationProgress::Matrix3fR::Zero();

  transformationMatrix.computeRotationScaling(&rotationMatrix, &scaleMatrix);

  //  for(size_t i=0; i<3; i++){
  //    for(size_t j=0; j<3; j++){
  //      rotationMatrix(i,j) = transformationMatrix.data()[i + 4*j];
  //    }
  //  }

  float xMin = std::numeric_limits<float>::max();
  float xMax = std::numeric_limits<float>::min();
  float yMin = std::numeric_limits<float>::max();
  float yMax = std::numeric_limits<float>::min();
  float zMin = std::numeric_limits<float>::max();
  float zMax = std::numeric_limits<float>::min();

  const std::vector<std::vector<size_t>> coords{{0, 0, 0},
                                                {origDims[0] - 1, 0, 0},
                                                {0, origDims[1] - 1, 0},
                                                {origDims[0] - 1, origDims[1] - 1, 0},
                                                {0, 0, origDims[2] - 1},
                                                {origDims[0] - 1, 0, origDims[2] - 1},
                                                {0, origDims[1] - 1, origDims[2] - 1},
                                                {origDims[0] - 1, origDims[1] - 1, origDims[2] - 1}};

  for(const auto& item : coords)
  {
    determineMinMax(rotationMatrix, spacing, item[0], item[1], item[2], xMin, xMax, yMin, yMax, zMin, zMax);
  }

  Eigen::Vector3f xAxisNew = rotationMatrix * k_XAxis;
  Eigen::Vector3f yAxisNew = rotationMatrix * k_YAxis;
  Eigen::Vector3f zAxisNew = rotationMatrix * k_ZAxis;

  float xResNew = determineSpacing(spacing, xAxisNew) * scaleMatrix(0,0);
  float yResNew = determineSpacing(spacing, yAxisNew) * scaleMatrix(1,1);
  float zResNew = determineSpacing(spacing, zAxisNew) * scaleMatrix(2,2);

  MeshIndexType xpNew = static_cast<int64_t>(std::nearbyint((xMax - xMin) / xResNew) + 1);
  MeshIndexType ypNew = static_cast<int64_t>(std::nearbyint((yMax - yMin) / yResNew) + 1);
  MeshIndexType zpNew = static_cast<int64_t>(std::nearbyint((zMax - zMin) / zResNew) + 1);

  RotateArgs params;

  params.xp = origDims[0];
  params.xRes = spacing[0];
  params.yp = origDims[1];
  params.yRes = spacing[1];
  params.zp = origDims[2];
  params.zRes = spacing[2];

  params.xpNew = xpNew;
  params.xResNew = xResNew;
  params.xMinNew = xMin;
  params.ypNew = ypNew;
  params.yResNew = yResNew;
  params.yMinNew = yMin;
  params.zpNew = zpNew;
  params.zResNew = zResNew;
  params.zMinNew = zMin;

  return params;
}

void updateGeometry(ImageGeom& imageGeom, const RotateArgs& params)
{
  FloatVec3Type origin = imageGeom.getOrigin();

  imageGeom.setSpacing(params.xResNew, params.yResNew, params.zResNew);
  imageGeom.setDimensions(params.xpNew, params.ypNew, params.zpNew);
  origin[0] += params.xMinNew;
  origin[1] += params.yMinNew;
  origin[2] += params.zMinNew;
  imageGeom.setOrigin(origin);
}

/**
 * @brief The RotateSampleRefFrameImpl class implements a threaded algorithm to do the
 * actual computation of the rotation by applying the rotation to each Euler angle
 */
class SampleRefFrameRotator
    {
  DataArray<int64_t>::Pointer m_NewIndicesPtr;
  float m_RotMatrixInv[3][3] = {{0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f}};
  float m_ScalingMatrix[3][3] = {{0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f}};
  float m_TranslationMatrix[3][1] = {0.0f, 0.0f, 0.0f};
  RotateArgs m_Params;

    public:
      SampleRefFrameRotator(DataArray<int64_t>::Pointer newindices, const RotateArgs& args, const Matrix3fR& rotationMatrix, const Matrix3fR& scalingMatrix, const MatrixTranslation translationMatrix)
      : m_NewIndicesPtr(newindices)
      , m_Params(args)
      {
        // We have to inline the 3x3 Maxtrix transpose here because of the "const" nature of the 'convert' function
        Matrix3fR transpose = rotationMatrix.transpose();
        // Need to use row based Eigen matrix so that the values get mapped to the right place in the raw array
        // Raw array is faster than Eigen
        Eigen::Map<Matrix3fR>(&m_RotMatrixInv[0][0], transpose.rows(), transpose.cols()) = transpose;
        Eigen::Map<Matrix3fR>(&m_ScalingMatrix[0][0], scalingMatrix.rows(), scalingMatrix.cols()) = scalingMatrix;
        Eigen::Map<MatrixTranslation>(&m_TranslationMatrix[0][0], translationMatrix.rows(), translationMatrix.cols()) = translationMatrix;
      }

      ~SampleRefFrameRotator() = default;

      void convert(int64_t zStart, int64_t zEnd, int64_t yStart, int64_t yEnd, int64_t xStart, int64_t xEnd) const
      {
        int64_t* newindicies = m_NewIndicesPtr->getPointer(0);

        for(int64_t k = zStart; k < zEnd; k++)
        {
          int64_t ktot = (m_Params.xpNew * m_Params.ypNew) * k;
          for(int64_t j = yStart; j < yEnd; j++)
          {
            int64_t jtot = (m_Params.xpNew) * j;
            for(int64_t i = xStart; i < xEnd; i++)
            {
              int64_t index = ktot + jtot + i;
              newindicies[index] = -1;

              float coords[3] = {0.0f, 0.0f, 0.0f};
              float coordsNew[3] = {0.0f, 0.0f, 0.0f};


              //TODO: Add in translation data here?
              coords[0] = (static_cast<float>(i) * m_Params.xResNew) + m_Params.xMinNew + m_TranslationMatrix[0][0];
              coords[1] = (static_cast<float>(j) * m_Params.yResNew) + m_Params.yMinNew + m_TranslationMatrix[1][0];
              coords[2] = (static_cast<float>(k) * m_Params.zResNew) + m_Params.zMinNew + m_TranslationMatrix[2][0];

              MatrixMath::Multiply3x3with3x1(m_RotMatrixInv, coords, coordsNew);

              //TODO: Linear Interpolation Implementation also scale this val after calculation
              int64_t colOld = static_cast<int64_t>((std::floor(coordsNew[0]) / m_Params.xRes)) * m_ScalingMatrix[3][0];
              int64_t rowOld = static_cast<int64_t>(std::nearbyint(coordsNew[1] / m_Params.yRes)) * m_ScalingMatrix[3][1];
              int64_t planeOld = static_cast<int64_t>(std::nearbyint(coordsNew[2] / m_Params.zRes)) * m_ScalingMatrix[3][2];

              if(colOld >= 0 && colOld < m_Params.xp && rowOld >= 0 && rowOld < m_Params.yp && planeOld >= 0 && planeOld < m_Params.zp)
              {
                newindicies[index] = (m_Params.xp * m_Params.yp * planeOld) + (m_Params.xp * rowOld) + colOld;
              }
            }
          }
        }
      }

#ifdef SIMPL_USE_PARALLEL_ALGORITHMS
void operator()(const tbb::blocked_range3d<int64_t, int64_t, int64_t>& r) const
{
        convert(r.pages().begin(), r.pages().end(), r.rows().begin(), r.rows().end(), r.cols().begin(), r.cols().end());
}
#endif
    };



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

struct ApplyTransformationToGeometry::Impl
    {
  ApplyTransformationProgress::Matrix3fR m_RotationMatrix = ApplyTransformationProgress::Matrix3fR::Zero();
  ApplyTransformationProgress::Matrix3fR m_ScalingMatrix = ApplyTransformationProgress::Matrix3fR::Zero();
  ApplyTransformationProgress::MatrixTranslation m_TranslationMatrix = ApplyTransformationProgress::MatrixTranslation::Zero();
  ApplyTransformationProgress::RotateArgs m_Params;


  void reset()
  {
    m_RotationMatrix.setZero();

    m_Params = ApplyTransformationProgress::RotateArgs();
  }
    };


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
ApplyTransformationToGeometry::ApplyTransformationToGeometry()
: p_Impl(std::make_unique<Impl>())
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
    std::vector<QString> choices;
    choices.push_back("No Transformation");
    choices.push_back("Pre-Computed Transformation Matrix");
    choices.push_back("Manual Transformation Matrix");
    choices.push_back("Rotation");
    choices.push_back("Translation");
    choices.push_back("Scale");
    parameter->setChoices(choices);
    std::vector<QString> linkedProps = {"ComputedTransformationMatrix", "ManualTransformationMatrix", "RotationAngle", "RotationAxis", "Translation", "Scale"};
    parameter->setLinkedProperties(linkedProps);
    parameter->setEditable(false);
    parameter->setCategory(FilterParameter::Category::Parameter);
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
    parameters.push_back(SIMPL_NEW_DYN_TABLE_FP("Transformation Matrix", ManualTransformationMatrix, FilterParameter::Category::Parameter, ApplyTransformationToGeometry, 2));
  }
  parameters.push_back(SIMPL_NEW_FLOAT_FP("Rotation Angle (Degrees)", RotationAngle, FilterParameter::Category::Parameter, ApplyTransformationToGeometry, 3));
  parameters.push_back(SIMPL_NEW_FLOAT_VEC3_FP("Rotation Axis (ijk)", RotationAxis, FilterParameter::Category::Parameter, ApplyTransformationToGeometry, 3));
  parameters.push_back(SIMPL_NEW_FLOAT_VEC3_FP("Translation", Translation, FilterParameter::Category::Parameter, ApplyTransformationToGeometry, 4));
  parameters.push_back(SIMPL_NEW_FLOAT_VEC3_FP("Scale", Scale, FilterParameter::Category::Parameter, ApplyTransformationToGeometry, 5));
  //  DataContainerSelectionFilterParameter::RequirementType dcReq;
  //  IGeometry::Types geomTypes = {IGeometry::Type::Vertex, IGeometry::Type::Edge, IGeometry::Type::Triangle, IGeometry::Type::Quad, IGeometry::Type::Tetrahedral, IGeometry::Type::Image};
  //  dcReq.dcGeometryTypes = geomTypes;
  //  parameters.push_back(SIMPL_NEW_DC_SELECTION_FP("Geometry to Transform", GeometryToTransform, FilterParameter::Category::RequiredArray, ApplyTransformationToGeometry, dcReq));
  //  {
  //    DataArraySelectionFilterParameter::RequirementType dasReq =
  //        DataArraySelectionFilterParameter::CreateRequirement(SIMPL::TypeNames::Float, SIMPL::Defaults::AnyComponentSize, AttributeMatrix::Type::Generic, IGeometry::Type::Any);
  //    parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Transformation Matrix", ComputedTransformationMatrix, FilterParameter::Category::RequiredArray, ApplyTransformationToGeometry, dasReq, 1));
  //  }

  parameters.push_back(SeparatorFilterParameter::Create("Cell Data", FilterParameter::Category::RequiredArray));
  {
    AttributeMatrixSelectionFilterParameter::RequirementType req = AttributeMatrixSelectionFilterParameter::CreateRequirement(AttributeMatrix::Type::Generic, IGeometry::Type::Any);
    parameters.push_back(SIMPL_NEW_AM_SELECTION_FP("Cell Attribute Matrix", CellAttributeMatrixPath, FilterParameter::Category::RequiredArray, ApplyTransformationToGeometry, req));
  }

  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ApplyTransformationToGeometry::readFilterParameters(AbstractFilterParametersReader* reader, int index)
{
  reader->openFilterGroup(this, index);
  setCellAttributeMatrixPath(reader->readDataArrayPath("CellAttributeMatrixPath", getCellAttributeMatrixPath()));
  setManualTransformationMatrix(reader->readDynamicTableData("ManualTransformationMatrix", getManualTransformationMatrix()));
  setComputedTransformationMatrix(reader->readDataArrayPath("ComputedTransformationMatrix", getComputedTransformationMatrix()));
  setTransformationMatrixType(reader->readValue("TransformationMatrixType", getTransformationMatrixType()));
  //  setGeometryToTransform(getCellAttributeMatrixPath());
  setRotationAxis(reader->readFloatVec3("RotationAxis", getRotationAxis()));
  setRotationAngle(reader->readValue("RotationAngle", getRotationAngle()));
  setTranslation(reader->readFloatVec3("Translation", getTranslation()));
  setScale(reader->readFloatVec3("Scale", getScale()));
  reader->closeFilterGroup();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------

//Need to add code in to create new image geom, in process, modifying updateGeometry, createRotateParams, determineMinMax
void ApplyTransformationToGeometry::dataCheck()
{
  using ProjectiveMatrix = Eigen::Matrix<float, 4, 4, Eigen::RowMajor>;
  using RotateMatrix = Eigen::Matrix<float, 3, 3, Eigen::RowMajor>;
  //  using TransformationObj = Eigen::Transform<float, 3, Eigen::Affine>;
  //  TransformationObj transform;
  clearErrorCode();
  clearWarningCode();

  p_Impl->reset();

  IGeometry::Pointer igeom = getDataContainerArray()->getPrereqGeometryFromDataContainer<IGeometry>(this, getCellAttributeMatrixPath().getDataContainerName());

  if(getErrorCode() < 0)
  {
    return;
  }

  if(!std::dynamic_pointer_cast<IGeometry2D>(igeom) && !std::dynamic_pointer_cast<IGeometry3D>(igeom) && !std::dynamic_pointer_cast<VertexGeom>(igeom) && !std::dynamic_pointer_cast<EdgeGeom>(igeom) && !std::dynamic_pointer_cast<ImageGeom>(igeom))
  {
    QString ss =
        QObject::tr("Geometry to transform must be an unstructured geometry (Vertex, Edge, Triangle, Quadrilateral, Tetrahedral, or Image), but the type is %1").arg(igeom->getGeometryTypeAsString());
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

  //if ImageGeom found:
  if(std::dynamic_pointer_cast<ImageGeom>(igeom))
  {
    //    DataContainer::Pointer m = getDataContainerArray()->getDataContainer(m_GeometryToTransform);
    ImageGeom::Pointer imageGeom =  std::dynamic_pointer_cast<ImageGeom>(igeom);
    Eigen::Map<ProjectiveMatrix> transformation(m_TransformationMatrix);
    Eigen::Transform<float, 3, Eigen::Affine> transform = Eigen::Transform<float, 3, Eigen::Affine>::Transform(transformation);
    ApplyTransformationProgress::Matrix3fR rotationMatrix = ApplyTransformationProgress::Matrix3fR::Zero();
    ApplyTransformationProgress::Matrix3fR scaleMatrix = ApplyTransformationProgress::Matrix3fR::Zero();
    ApplyTransformationProgress::MatrixTranslation translationMatrix = ApplyTransformationProgress::MatrixTranslation::Zero();

    transform.computeRotationScaling(&rotationMatrix, &scaleMatrix);
    translationMatrix(0,0) = transform.data()[3];
    translationMatrix(1,0) = transform.data()[7];
    translationMatrix(2,0) = transform.data()[11];

    p_Impl->m_RotationMatrix = rotationMatrix;
    p_Impl->m_ScalingMatrix = scaleMatrix;
    p_Impl->m_TranslationMatrix = translationMatrix;

    p_Impl->m_Params = ApplyTransformationProgress::createRotateParams(*imageGeom, transform);
    updateGeometry(*imageGeom, p_Impl->m_Params);

    // Resize attribute matrix

    std::vector<size_t> tDims(3);
    tDims[0] = p_Impl->m_Params.xpNew;
    tDims[1] = p_Impl->m_Params.ypNew;
    tDims[2] = p_Impl->m_Params.zpNew;
    //    QString attrMatName = getCellAttributeMatrixPath().getAttributeMatrixName();
    DataContainer::Pointer m = getDataContainerArray()->getDataContainer(getCellAttributeMatrixPath().getDataContainerName());
    QString attrMatName = getCellAttributeMatrixPath().getAttributeMatrixName();
    m->getAttributeMatrix(attrMatName)->resizeAttributeArrays(tDims);
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------

void ApplyTransformationToGeometry::ApplyImageTransformation()
{
  DataContainer::Pointer m = getDataContainerArray()->getDataContainer(getCellAttributeMatrixPath().getDataContainerName());
  int64_t newNumCellTuples = p_Impl->m_Params.xpNew * p_Impl->m_Params.ypNew * p_Impl->m_Params.zpNew;

  DataArray<int64_t>::Pointer newIndiciesPtr = DataArray<int64_t>::CreateArray(newNumCellTuples, std::string("_INTERNAL_USE_ONLY_RotateSampleRef_NewIndicies"), true);
  newIndiciesPtr->initializeWithValue(-1);
  int64_t* newindicies = newIndiciesPtr->getPointer(0);

#ifdef SIMPL_USE_PARALLEL_ALGORITHMS
  tbb::parallel_for(tbb::blocked_range3d<int64_t, int64_t, int64_t>(0, p_Impl->m_Params.zpNew, 0, p_Impl->m_Params.ypNew, 0, p_Impl->m_Params.xpNew),
                    ApplyTransformationProgress::SampleRefFrameRotator(newIndiciesPtr, p_Impl->m_Params, p_Impl->m_RotationMatrix, p_Impl->m_ScalingMatrix, p_Impl->m_TranslationMatrix), tbb::auto_partitioner());
#else
  {
    SampleRefFrameRotator serial(newIndiciesPtr, p_Impl->m_Params, p_Impl->m_RotationMatrix, p_Impl->m_ScalingMatrix, p_Impl->m_TranslationMatrix);
    serial.convert(0, p_Impl->m_Params.zpNew, 0, p_Impl->m_Params.ypNew, 0, p_Impl->m_Params.xpNew);
  }
#endif

QString attrMatName = getCellAttributeMatrixPath().getAttributeMatrixName();
  QList<QString> voxelArrayNames = m->getAttributeMatrix(attrMatName)->getAttributeArrayNames();

  for(const auto& attrArrayName : voxelArrayNames)
  {
    IDataArray::Pointer p = m->getAttributeMatrix(attrMatName)->getAttributeArray(attrArrayName);

    // Make a copy of the 'p' array that has the same name. When placed into
    // the data container this will over write the current array with
    // the same name.

    IDataArray::Pointer data = p->createNewArray(newNumCellTuples, p->getComponentDimensions(), p->getName());
    int64_t newIndicies_I = 0;
    for(size_t i = 0; i < static_cast<size_t>(newNumCellTuples); i++)
    {
      newIndicies_I = newindicies[i];
      if(newIndicies_I >= 0)
      {
        if(!data->copyFromArray(i, p, newIndicies_I, 1))
        {
          QString ss = QObject::tr("copyFromArray Failed: ");
          QTextStream out(&ss);
          out << "Source Array Name: " << p->getName() << " Source Tuple Index: " << newIndicies_I << "\n";
          out << "Dest Array Name: " << data->getName() << "  Dest. Tuple Index: " << i << "\n";
          setErrorCondition(-45102, ss);
          return;
        }
      }
      else
      {
        int var = 0;
        data->initializeTuple(i, &var);
      }
    }
    m->getAttributeMatrix(attrMatName)->insertOrAssign(data);
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------


void ApplyTransformationToGeometry::applyTransformation()
{

  IGeometry::Pointer igeom = getDataContainerArray()->getDataContainer(getCellAttributeMatrixPath().getDataContainerName())->getGeometry();
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
  else if(ImageGeom::Pointer image = std::dynamic_pointer_cast<ImageGeom>(igeom)){
    ApplyImageTransformation();//Function for applying Image Transformation
    return;
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
void ApplyTransformationToGeometry::setCellAttributeMatrixPath(const DataArrayPath& value)
{
  m_CellAttributeMatrixPath = value;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
DataArrayPath ApplyTransformationToGeometry::getCellAttributeMatrixPath() const
{
  return m_CellAttributeMatrixPath;
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
//void ApplyTransformationToGeometry::setGeometryToTransform(const DataArrayPath& value)
//{
//  m_GeometryToTransform = value;
//}

// -----------------------------------------------------------------------------
//DataArrayPath ApplyTransformationToGeometry::getGeometryToTransform() const
//{
//  return m_GeometryToTransform;
//}

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
