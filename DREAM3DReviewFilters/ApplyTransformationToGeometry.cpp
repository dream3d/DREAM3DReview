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

#include <QtCore/QTextStream>

#include "SIMPLib/Common/Constants.h"
#include "SIMPLib/Common/SIMPLRange.h"
#include "SIMPLib/DataContainers/DataContainerArray.h"
#include "SIMPLib/FilterParameters/AbstractFilterParametersReader.h"
#include "SIMPLib/FilterParameters/AttributeMatrixSelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/ChoiceFilterParameter.h"
#include "SIMPLib/FilterParameters/DataArraySelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/DynamicTableFilterParameter.h"
#include "SIMPLib/FilterParameters/FloatFilterParameter.h"
#include "SIMPLib/FilterParameters/FloatVec3FilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedBooleanFilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedChoicesFilterParameter.h"
#include "SIMPLib/FilterParameters/MultiDataArraySelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/SeparatorFilterParameter.h"
#include "SIMPLib/Geometry/EdgeGeom.h"
#include "SIMPLib/Geometry/IGeometry2D.h"
#include "SIMPLib/Geometry/IGeometry3D.h"
#include "SIMPLib/Geometry/ImageGeom.h"
#include "SIMPLib/Geometry/VertexGeom.h"
#include "SIMPLib/Math/MatrixMath.h"
#include "SIMPLib/Math/SIMPLibMath.h"
#include "SIMPLib/Utilities/ParallelDataAlgorithm.h"
#include "SIMPLib/Utilities/ParallelTaskAlgorithm.h"

#include "DREAM3DReview/DREAM3DReviewConstants.h"
#include "DREAM3DReview/DREAM3DReviewVersion.h"

namespace OrientationTransformation
{

/**
 * The Orientation codes are written in such a way that the value of -1 indicates
 * an Active Rotation and +1 indicates a passive rotation.
 *
 * DO NOT UNDER ANY CIRCUMSTANCE CHANGE THESE VARIABLES. THERE WILL BE BAD
 * CONSEQUENCES IF THESE ARE CHANGED. EVERY PIECE OF CODE THAT RELIES ON THESE
 * FUNCTIONS WILL BREAK. IN ADDITION, THE QUATERNION ARITHMETIC WILL NO LONGER
 * BE CONSISTENT WITH ROTATION ARITHMETIC.
 *
 * YOU HAVE BEEN WARNED.
 *
 * Adam  Morawiec's book uses Passive rotations.
 *
 * This code was taken from "EbsdLib/Core/OrientationTransformation.hpp"
 **/
namespace Rotations::Constants
{

#define DREAM3D_PASSIVE_ROTATION

#ifdef DREAM3D_PASSIVE_ROTATION
constexpr float epsijk = 1.0f;
// constexpr double epsijkd = 1.0;
#elif DREAM3D_ACTIVE_ROTATION
static const float epsijk = -1.0f;
// static const double epsijkd = -1.0;
#endif
} // namespace Rotations::Constants
template <typename InputType, typename OutputType>
OutputType ax2om(const InputType& a)
{
  OutputType res(9);
  using value_type = typename OutputType::value_type;

  value_type c = cos(a[3]);
  value_type s = sin(a[3]);

  value_type omc = static_cast<value_type>(1.0 - c);

  res[0] = a[0] * a[0] * omc + c;
  res[4] = a[1] * a[1] * omc + c;
  res[8] = a[2] * a[2] * omc + c;
  int _01 = 1;
  int _10 = 3;
  int _12 = 5;
  int _21 = 7;
  int _02 = 2;
  int _20 = 6;
  // Check to see if we need to transpose
  if(Rotations::Constants::epsijk == 1.0F)
  {
    _01 = 3;
    _10 = 1;
    _12 = 7;
    _21 = 5;
    _02 = 6;
    _20 = 2;
  }

  value_type q = omc * a[0] * a[1];
  res[_01] = q + s * a[2];
  res[_10] = q - s * a[2];
  q = omc * a[1] * a[2];
  res[_12] = q + s * a[0];
  res[_21] = q - s * a[0];
  q = omc * a[2] * a[0];
  res[_02] = q - s * a[1];
  res[_20] = q + s * a[1];

  return res;
}
} // namespace OrientationTransformation

namespace
{
constexpr int32_t k_NearestNeighborInterpolation = 0;
constexpr int32_t k_LinearInterpolation = 1;

constexpr int32_t k_NoTransform = 0;
constexpr int32_t k_PrecomputedTransform = 1;
constexpr int32_t k_ManualTransform = 2;
constexpr int32_t k_RotationTransform = 3;
constexpr int32_t k_TranslationTransform = 4;
constexpr int32_t k_ScaleTransform = 5;

const Eigen::Vector3f k_XAxis = Eigen::Vector3f::UnitX();
const Eigen::Vector3f k_YAxis = Eigen::Vector3f::UnitY();
const Eigen::Vector3f k_ZAxis = Eigen::Vector3f::UnitZ();

/**
 * @brief Function to determine the min and max coordinates (bounding box) of the transformed Image Geometry.
 * @param imageGeometry
 * @param transformationMatrix
 * @return
 */
FloatVec6Type DetermineMinMaxCoords(const ImageGeom& imageGeometry, const ApplyTransformationToGeometry::Matrix4fR& transformationMatrix)
{
  auto origImageGeomBox = imageGeometry.getBoundingBox();
  // clang-format off
  std::vector<FloatVec3Type> imageGeomCornerCoords = {{origImageGeomBox[0], origImageGeomBox[2], origImageGeomBox[4]},
                                                      {origImageGeomBox[1], origImageGeomBox[2], origImageGeomBox[4]},
                                                      {origImageGeomBox[1], origImageGeomBox[3], origImageGeomBox[4]},
                                                      {origImageGeomBox[0], origImageGeomBox[3], origImageGeomBox[4]},
                                                      {origImageGeomBox[0], origImageGeomBox[2], origImageGeomBox[5]},
                                                      {origImageGeomBox[1], origImageGeomBox[2], origImageGeomBox[5]},
                                                      {origImageGeomBox[1], origImageGeomBox[3], origImageGeomBox[5]},
                                                      {origImageGeomBox[0], origImageGeomBox[3], origImageGeomBox[5]}};
  // clang-format on
  FloatVec6Type minMaxValues = {std::numeric_limits<float>::max(),  -std::numeric_limits<float>::max(), std::numeric_limits<float>::max(),
                                -std::numeric_limits<float>::max(), std::numeric_limits<float>::max(),  -std::numeric_limits<float>::max()};

  for(size_t i = 0; i < 8; i++)
  {
    Eigen::Vector4f coords(imageGeomCornerCoords[i][0], imageGeomCornerCoords[i][1], imageGeomCornerCoords[i][2], 1.0F);

    Eigen::Vector4f newCoords = transformationMatrix * coords;

    minMaxValues[0] = std::min(newCoords[0], minMaxValues[0]);
    minMaxValues[1] = std::max(newCoords[0], minMaxValues[1]);

    minMaxValues[2] = std::min(newCoords[1], minMaxValues[2]);
    minMaxValues[3] = std::max(newCoords[1], minMaxValues[3]);

    minMaxValues[4] = std::min(newCoords[2], minMaxValues[4]);
    minMaxValues[5] = std::max(newCoords[2], minMaxValues[5]);
  }
  return minMaxValues;
}

/**
 * @brief Finds the Cosine of the angle between 2 vectors
 * @param a
 * @param b
 * @return
 */
float CosBetweenVectors(const Eigen::Vector3f& a, const Eigen::Vector3f& b)
{
  float normA = a.norm();
  float normB = b.norm();

  if(normA == 0.0f || normB == 0.0f)
  {
    return 1.0f;
  }

  return a.dot(b) / (normA * normB);
}

/**
 * @brief Function for determining new ImageGeom Spacing between points for scaling
 * @param spacing
 * @param axisNew
 * @return spacing for a given axis.
 */
float DetermineSpacing(const FloatVec3Type& spacing, const Eigen::Vector3f& axisNew)
{
  float xAngle = std::abs(CosBetweenVectors(k_XAxis, axisNew));
  float yAngle = std::abs(CosBetweenVectors(k_YAxis, axisNew));
  float zAngle = std::abs(CosBetweenVectors(k_ZAxis, axisNew));

  std::array<float, 3> axes = {xAngle, yAngle, zAngle};

  auto iter = std::max_element(axes.cbegin(), axes.cend());

  size_t index = std::distance(axes.cbegin(), iter);

  return spacing[index];
}

/**
 * @brief Determines paramaters for image rotation
 * @param imageGeom
 * @param transformationMatrix
 * @return New RotateArgs object
 */
ImageRotationUtilities::RotateArgs CreateRotationArgs(const ImageGeom& imageGeom, const ApplyTransformationToGeometry::Matrix4fR& transformationMatrix)
{
  const SizeVec3Type origDims = imageGeom.getDimensions();
  const FloatVec3Type spacing = imageGeom.getSpacing();

  ApplyTransformationToGeometry::Matrix3fR rotationMatrix = transformationMatrix.block(0, 0, 3, 3);

  FloatVec6Type minMaxCoords = DetermineMinMaxCoords(imageGeom, transformationMatrix);

  Eigen::Vector3f xAxisNew = rotationMatrix * k_XAxis;
  Eigen::Vector3f yAxisNew = rotationMatrix * k_YAxis;
  Eigen::Vector3f zAxisNew = rotationMatrix * k_ZAxis;

  float xResNew = DetermineSpacing(spacing, xAxisNew);
  float yResNew = DetermineSpacing(spacing, yAxisNew);
  float zResNew = DetermineSpacing(spacing, zAxisNew);

  MeshIndexType xpNew = static_cast<int64_t>(std::nearbyint((minMaxCoords[1] - minMaxCoords[0]) / xResNew));
  MeshIndexType ypNew = static_cast<int64_t>(std::nearbyint((minMaxCoords[3] - minMaxCoords[2]) / yResNew));
  MeshIndexType zpNew = static_cast<int64_t>(std::nearbyint((minMaxCoords[5] - minMaxCoords[4]) / zResNew));

  ImageRotationUtilities::RotateArgs params;

  params.origImageGeom = ImageGeom::New();
  params.origImageGeom->setSpacing(imageGeom.getSpacing());
  params.origImageGeom->setDimensions(imageGeom.getDimensions());
  params.origImageGeom->setOrigin(imageGeom.getOrigin());

  params.xp = origDims[0];
  params.xRes = spacing[0];
  params.yp = origDims[1];
  params.yRes = spacing[1];
  params.zp = origDims[2];
  params.zRes = spacing[2];

  params.transformedImageGeom = ImageGeom::New();
  params.transformedImageGeom->setSpacing(xResNew, yResNew, zResNew);
  params.transformedImageGeom->setDimensions(xpNew, ypNew, zpNew);
  params.transformedImageGeom->setOrigin(minMaxCoords[0], minMaxCoords[2], minMaxCoords[4]);

  params.xpNew = xpNew;
  params.xResNew = xResNew;
  params.xMinNew = minMaxCoords[0];
  params.ypNew = ypNew;
  params.yResNew = yResNew;
  params.yMinNew = minMaxCoords[2];
  params.zpNew = zpNew;
  params.zResNew = zResNew;
  params.zMinNew = minMaxCoords[4];

  return params;
}

} // namespace

/**
 * @brief Used to transform node based geometries (Vertex, Edge, Triangle... )
 */
static size_t s_InstanceIndex = 0;
static std::map<size_t, int64_t> s_ProgressValues;
static std::map<size_t, int64_t> s_LastProgressInt;

/**
 * @brief The ApplyTransformationToGeometryImpl class will apply a transformation to a node based geometry.
 */
class ApplyTransformationToGeometryImpl
{

public:
  ApplyTransformationToGeometryImpl(ApplyTransformationToGeometry& filter, const Matrix4fR& transformationMatrix, const SharedVertexList::Pointer verticesPtr)
  : m_Filter(filter)
  , m_TransformationMatrix(transformationMatrix)
  , m_Vertices(verticesPtr)
  {
  }
  ~ApplyTransformationToGeometryImpl() = default;

  ApplyTransformationToGeometryImpl(const ApplyTransformationToGeometryImpl&) = default;           // Copy Constructor Not Implemented
  ApplyTransformationToGeometryImpl(ApplyTransformationToGeometryImpl&&) = delete;                 // Move Constructor Not Implemented
  ApplyTransformationToGeometryImpl& operator=(const ApplyTransformationToGeometryImpl&) = delete; // Copy Assignment Not Implemented
  ApplyTransformationToGeometryImpl& operator=(ApplyTransformationToGeometryImpl&&) = delete;      // Move Assignment Not Implemented

  void convert(size_t start, size_t end) const
  {
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
      Eigen::Vector4f transformedPosition = m_TransformationMatrix * position;
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
  const Matrix4fR m_TransformationMatrix;
  SharedVertexList::Pointer m_Vertices;
};

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
ApplyTransformationToGeometry::ApplyTransformationToGeometry()
{
  m_RotationAngle = 0.0f;
  m_RotationAxis[0] = 0.0f;
  m_RotationAxis[1] = 0.0f;
  m_RotationAxis[2] = 1.0f;

  m_Translation[0] = 0.0f;
  m_Translation[1] = 0.0f;
  m_Translation[2] = 0.0f;

  m_Scale[0] = 1.0f;
  m_Scale[1] = 1.0f;
  m_Scale[2] = 1.0f;
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
    ChoiceFilterParameter::Pointer parameter2 = ChoiceFilterParameter::New();
    parameter2->setHumanLabel("Interpolation Type");
    parameter2->setPropertyName("InterpolationType");
    parameter2->setSetterCallback(SIMPL_BIND_SETTER(ApplyTransformationToGeometry, this, InterpolationType));
    parameter2->setGetterCallback(SIMPL_BIND_GETTER(ApplyTransformationToGeometry, this, InterpolationType));
    std::vector<QString> choices = {"Nearest Neighbor", "Linear"};
    parameter2->setChoices(choices);
    parameter2->setEditable(false);
    parameter2->setCategory(FilterParameter::Category::Parameter);
    parameters.push_back(parameter2);
  }
  {
    LinkedBooleanFilterParameter::Pointer parameter3 = LinkedBooleanFilterParameter::New();
    parameter3->setHumanLabel("Select Data Arrays");
    parameter3->setPropertyName("UseDataArraySelection");
    parameter3->setSetterCallback(SIMPL_BIND_SETTER(ApplyTransformationToGeometry, this, UseDataArraySelection));
    parameter3->setGetterCallback(SIMPL_BIND_GETTER(ApplyTransformationToGeometry, this, UseDataArraySelection));
    std::vector<QString> linkedProps2 = {"DataArraySelection"};
    parameter3->setConditionalProperties(linkedProps2);
    parameter3->setCategory(FilterParameter::Category::Parameter);
    parameters.push_back(parameter3);
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
    parameters.push_back(SIMPL_NEW_DYN_TABLE_FP("Transformation Matrix", ManualTransformationMatrix, FilterParameter::Category::Parameter, ApplyTransformationToGeometry, {2}));
  }
  parameters.push_back(SIMPL_NEW_FLOAT_FP("Rotation Angle (Degrees)", RotationAngle, FilterParameter::Category::Parameter, ApplyTransformationToGeometry, {3}));
  parameters.push_back(SIMPL_NEW_FLOAT_VEC3_FP("Rotation Axis (ijk)", RotationAxis, FilterParameter::Category::Parameter, ApplyTransformationToGeometry, {3}));
  parameters.push_back(SIMPL_NEW_FLOAT_VEC3_FP("Translation", Translation, FilterParameter::Category::Parameter, ApplyTransformationToGeometry, {4}));
  parameters.push_back(SIMPL_NEW_FLOAT_VEC3_FP("Scale", Scale, FilterParameter::Category::Parameter, ApplyTransformationToGeometry, {5}));

  {
    DataArraySelectionFilterParameter::RequirementType datReq =
        DataArraySelectionFilterParameter::CreateRequirement(SIMPL::TypeNames::Float, SIMPL::Defaults::AnyComponentSize, AttributeMatrix::Type::Generic, IGeometry::Type::Any);
    parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Transformation Matrix", ComputedTransformationMatrix, FilterParameter::Category::RequiredArray, ApplyTransformationToGeometry, datReq, {1}));
  }

  parameters.push_back(SeparatorFilterParameter::Create("Cell Data", FilterParameter::Category::RequiredArray));
  {
    AttributeMatrixSelectionFilterParameter::RequirementType amReq;
    IGeometry::Types geomTypes = {IGeometry::Type::Vertex, IGeometry::Type::Edge, IGeometry::Type::Triangle, IGeometry::Type::Quad, IGeometry::Type::Tetrahedral, IGeometry::Type::Image};
    amReq.dcGeometryTypes = geomTypes;
    AttributeMatrix::Types amType = {AttributeMatrix::Type::Any};
    amReq.amTypes = amType;

    parameters.push_back(SIMPL_NEW_AM_SELECTION_FP("Cell Attribute Matrix", CellAttributeMatrixPath, FilterParameter::Category::RequiredArray, ApplyTransformationToGeometry, amReq));
  }

  {
    MultiDataArraySelectionFilterParameter::RequirementType dasReq;
    dasReq.amTypes = AttributeMatrix::Types(1, AttributeMatrix::Type::Cell);
    parameters.push_back(SIMPL_NEW_MDA_SELECTION_FP("Data Array Selection", DataArraySelection, FilterParameter::Category::RequiredArray, ApplyTransformationToGeometry, dasReq));
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
  setInterpolationType(reader->readValue("InterpolationType", getInterpolationType()));
  setRotationAxis(reader->readFloatVec3("RotationAxis", getRotationAxis()));
  setRotationAngle(reader->readValue("RotationAngle", getRotationAngle()));
  setTranslation(reader->readFloatVec3("Translation", getTranslation()));
  setScale(reader->readFloatVec3("Scale", getScale()));
  setUseDataArraySelection(reader->readValue("UseDataArraySelection", getUseDataArraySelection()));
  // setDataArraySelection(reader->readDataArrayPathVector("DataArraySelection", getDataArraySelection()));
  reader->closeFilterGroup();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------

void ApplyTransformationToGeometry::dataCheck()
{
  // using ProjectiveMatrix = Eigen::Matrix<float, 4, 4, Eigen::RowMajor>;

  clearErrorCode();
  clearWarningCode();

  reset();

  IGeometry::Pointer igeom = getDataContainerArray()->getPrereqGeometryFromDataContainer<IGeometry>(this, getCellAttributeMatrixPath().getDataContainerName());

  if(getErrorCode() < 0)
  {
    return;
  }

  int err = 0;
  AttributeMatrixShPtr attriPtr = this->getDataContainerArray()->getPrereqAttributeMatrixFromPath(this, getCellAttributeMatrixPath(), err);

  if(getErrorCode() < 0)
  {
    return;
  }

  if(!std::dynamic_pointer_cast<IGeometry2D>(igeom) && !std::dynamic_pointer_cast<IGeometry3D>(igeom) && !std::dynamic_pointer_cast<VertexGeom>(igeom) && !std::dynamic_pointer_cast<EdgeGeom>(igeom) &&
     !std::dynamic_pointer_cast<ImageGeom>(igeom))
  {
    QString ss =
        QObject::tr("Geometry to transform must be an unstructured geometry (Vertex, Edge, Triangle, Quadrilateral, Tetrahedral, or Image), but the type is %1").arg(igeom->getGeometryTypeAsString());
    setErrorCondition(-702, ss);
  }

  std::vector<size_t> cDims = {4, 4};

  // Reset the final Transformation Matrix to all Zeros before we fill it with what the user has entered.
  m_TransformationMatrix.fill(0.0F);

  switch(getTransformationMatrixType())
  {
  case k_NoTransform: // No-Op
  {
    QString ss = QObject::tr("No transformation has been selected, so this filter will perform no operations");
    setWarningCondition(-701, ss);
  }
  case k_PrecomputedTransform: // Transformation matrix from array
  {
    auto precomputedTransformationMatrix = getDataContainerArray()->getPrereqArrayFromPath<FloatArrayType>(this, getComputedTransformationMatrix(), cDims);
    if(getErrorCode() < 0)
    {
      return;
    }
    if(!getInPreflight())
    {
      std::vector<size_t> cDims = {4, 4};
      auto precomputedTransformationMatrix = getDataContainerArray()->getPrereqArrayFromPath<FloatArrayType>(this, getComputedTransformationMatrix(), cDims);
      for(size_t r = 0; r < 4; r++)
      {
        for(size_t c = 0; c < 4; c++)
        {
          m_TransformationMatrix(r, c) = (*precomputedTransformationMatrix)[r * 4 + c];
        }
      }
    }
    break;
  }
  case k_ManualTransform: // Manual transformation matrix
  {
    int numTableRows = getManualTransformationMatrix().getNumRows();
    int numTableCols = getManualTransformationMatrix().getNumCols();
    if(numTableRows != 4)
    {
      QString ss = QObject::tr("Manually entered transformation matrix must have exactly 4 rows");
      setErrorCondition(-702, ss);
      return;
    }
    if(numTableCols != 4)
    {
      QString ss = QObject::tr("Manually entered transformation matrix must have exactly 4 columns");
      setErrorCondition(-703, ss);
      return;
    }
    std::vector<std::vector<double>> tableData = getManualTransformationMatrix().getTableData();

    for(size_t rowIndex = 0; rowIndex < numTableRows; rowIndex++)
    {
      std::vector<double> row = tableData[rowIndex];
      for(size_t colIndex = 0; colIndex < numTableCols; colIndex++)
      {
        m_TransformationMatrix(rowIndex, colIndex) = static_cast<float>(row[colIndex]);
      }
    }

    break;
  }
  case k_RotationTransform: // Rotation via axis-angle
  {
    // Convert Degrees to Radians for the last element
    float rotAngle = m_RotationAngle * SIMPLib::Constants::k_PiOver180D;

    // Ensure the axis part is normalized
    FloatVec3Type normalizedAxis = m_RotationAxis;
    MatrixMath::Normalize3x1(normalizedAxis.data());

    float cosTheta = cos(rotAngle);
    float oneMinusCosTheta = 1 - cosTheta;
    float sinTheta = sin(rotAngle);
    float l = normalizedAxis[0];
    float m = normalizedAxis[1];
    float n = normalizedAxis[2];

    // First Row:
    m_TransformationMatrix(0) = l * l * (oneMinusCosTheta) + cosTheta;
    m_TransformationMatrix(1) = m * l * (oneMinusCosTheta) - (n * sinTheta);
    m_TransformationMatrix(2) = n * l * (oneMinusCosTheta) + (m * sinTheta);
    m_TransformationMatrix(3) = 0.0F;

    // Second Row:
    m_TransformationMatrix(4) = l * m * (oneMinusCosTheta) + (n * sinTheta);
    m_TransformationMatrix(5) = m * m * (oneMinusCosTheta) + cosTheta;
    m_TransformationMatrix(6) = n * m * (oneMinusCosTheta) - (l * sinTheta);
    m_TransformationMatrix(7) = 0.0F;

    // Third Row:
    m_TransformationMatrix(8) = l * n * (oneMinusCosTheta) - (m * sinTheta);
    m_TransformationMatrix(9) = m * n * (oneMinusCosTheta) + (l * sinTheta);
    m_TransformationMatrix(10) = n * n * (oneMinusCosTheta) + cosTheta;
    m_TransformationMatrix(11) = 0.0F;

    // Fourth Row:
    m_TransformationMatrix(12) = 0.0F;
    m_TransformationMatrix(13) = 0.0F;
    m_TransformationMatrix(14) = 0.0F;
    m_TransformationMatrix(15) = 1.0F;

    break;
  }
  case k_TranslationTransform: // Translation
  {
    m_TransformationMatrix(0, 0) = 1.0f;
    m_TransformationMatrix(1, 1) = 1.0f;
    m_TransformationMatrix(2, 2) = 1.0f;
    m_TransformationMatrix(3, 3) = 1.0f;
    m_TransformationMatrix(0, 3) = m_Translation[0];
    m_TransformationMatrix(1, 3) = m_Translation[1];
    m_TransformationMatrix(2, 3) = m_Translation[2];
    break;
  }
  case k_ScaleTransform: // Scale
  {
    m_TransformationMatrix(0, 0) = m_Scale[0];
    m_TransformationMatrix(1, 1) = m_Scale[1];
    m_TransformationMatrix(2, 2) = m_Scale[2];
    m_TransformationMatrix(3, 3) = 1.0f;
    break;
  }
  default: {
    QString ss = QObject::tr("Invalid selection for transformation type");
    setErrorCondition(-701, ss);
    return;
  }
  }

  // if ImageGeom was selected to be transformed:
  if(std::dynamic_pointer_cast<ImageGeom>(igeom))
  {
    // DataContainer::Pointer m = getDataContainerArray()->getDataContainer(getCellAttributeMatrixPath().getDataContainerName());
    // For nearest neighbor we can use any data array as we are only going to copy from a hit cell, no interpolation
    QList<QString> selectedCellArrayNames = attriPtr->getAttributeArrayNames();

    if(getInterpolationType() == k_LinearInterpolation)
    {
      if(m_UseDataArraySelection)
      {
        selectedCellArrayNames.clear();
        for(const auto& dataArrayPath : m_DataArraySelection)
        {
          selectedCellArrayNames.append(dataArrayPath.getDataArrayName());
        }
      }
      else
      {
        setWarningCondition(-705, "Interpolation has been selected but no Cell arrays are selected to be interpolated. No data will be transferred to the new geometry");
      }

      for(const auto& attrArrayName : selectedCellArrayNames)
      {
        IDataArray::Pointer p = attriPtr->getAttributeArray(attrArrayName);
        if(p->getTypeAsString() == "bool" || p->getTypeAsString() == "StringDataArray" || p->getTypeAsString() == "StatsDataArray")
        {
          QString ss = QObject::tr("Input Type Error, cannot run linear interpolation on a boolean, stirngs or StatsData arrays");
          setErrorCondition(-706, ss);
        }
      }
    }
    ImageGeom::Pointer imageGeom = std::dynamic_pointer_cast<ImageGeom>(igeom);

    m_Params = ::CreateRotationArgs(*imageGeom, m_TransformationMatrix);

    // If the user is purely doing a translation then just adjust the origin and be done.
    if(getTransformationMatrixType() == k_TranslationTransform)
    {
      imageGeom->setOrigin(m_Params.origImageGeom->getOrigin() + m_Translation);
    }
    else
    {
      imageGeom->setDimensions(m_Params.transformedImageGeom->getDimensions());
      imageGeom->setOrigin(m_Params.transformedImageGeom->getOrigin());
      imageGeom->setSpacing(m_Params.transformedImageGeom->getSpacing());

      DataContainer::Pointer m = getDataContainerArray()->getDataContainer(getCellAttributeMatrixPath().getDataContainerName());

      // Resize attribute matrix
      std::vector<size_t> tDims = {static_cast<size_t>(m_Params.xpNew), static_cast<size_t>(m_Params.ypNew), static_cast<size_t>(m_Params.zpNew)};
      // size_t numTuples = tDims[0] * tDims[1] * tDims[2];

      QString attrMatName = getCellAttributeMatrixPath().getAttributeMatrixName();
      // Get the List of Array Names FIRST
      // QList<QString> voxelArrayNames = m->getAttributeMatrix(attrMatName)->getAttributeArrayNames();
      // Now remove the current Cell Attribute Matrix and store it in the instance variable
      m_SourceAttributeMatrix = m->removeAttributeMatrix(attrMatName);
      // Now create a new Attribute Matrix that has the correct Tuple Dims.
      AttributeMatrix::Pointer targetAttributeMatrix = m->createNonPrereqAttributeMatrix(this, attrMatName, tDims, AttributeMatrix::Type::Cell);
      // Loop over all of the original cell data arrays and create new ones and insert that into the new Attribute Matrix.
      // DO NOT ALLOCATE the arrays, even during execute as this could potentially be a LARGE memory hog. Wait until
      // execute to allocate the arrays one at a time, do the copy, then deallocate the old array. This will keep the memory
      // consumption to a minimum.
      for(const auto& attrArrayName : selectedCellArrayNames)
      {
        IDataArray::Pointer p = m_SourceAttributeMatrix->getAttributeArray(attrArrayName);
        auto compDims = p->getComponentDimensions();
        IDataArray::Pointer targetArray = p->createNewArray(tDims[0] * tDims[1] * tDims[2], compDims, p->getName(), false);
        targetAttributeMatrix->addOrReplaceAttributeArray(targetArray);
      }
    }
  }
}

/**
 * @brief ApplyTransformationToGeometry::applyImageGeometryTransformation
 */
void ApplyTransformationToGeometry::applyImageGeometryTransformation()
{
  DataContainer::Pointer m = getDataContainerArray()->getDataContainer(getCellAttributeMatrixPath().getDataContainerName());

  if(m == nullptr)
  {
    QString ss = QObject::tr("Failed to get DataContainer '%1'").arg(getCellAttributeMatrixPath().getDataContainerName());
    setErrorCondition(-45101, ss);
    return;
  }

  ImageGeom::Pointer imageGeom = m->getGeometryAs<ImageGeom>();
  // If we are doing a pure translation then adjust the origin and be done.
  if(getTransformationMatrixType() == k_TranslationTransform)
  {
    imageGeom->setOrigin(m_Params.origImageGeom->getOrigin() + m_Translation);
    return;
  }

  int64_t newNumCellTuples = m_Params.xpNew * m_Params.ypNew * m_Params.zpNew;
  m_TotalElements = newNumCellTuples;

  QString attrMatName = getCellAttributeMatrixPath().getAttributeMatrixName();
  // AttributeMatrix::Pointer targetAttributeMatrix = m->getAttributeMatrix(attrMatName);

  QList<QString> selectedCellArrayNames = m->getAttributeMatrix(attrMatName)->getAttributeArrayNames();

  ParallelTaskAlgorithm taskRunner;
  taskRunner.setParallelizationEnabled(true);

  for(const auto& attrArrayName : selectedCellArrayNames)
  {
    // Get the source array from our "cached" Cell AttributeMatrix
    IDataArray::Pointer sourceArray = m_SourceAttributeMatrix->getAttributeArray(attrArrayName);
    // Get the target attribute array from our destination Cell AttributeMatrix which has *NOT* been allocated yet.
    IDataArray::Pointer targetArray = m->getAttributeMatrix(attrMatName)->getAttributeArray(attrArrayName);

    // So this little work-around is because if we just try to resize the DataArray<T> will think the sizes are the same
    // and never actually allocate the data. So we just resize to 1 tuple, and then to the real size.
    targetArray->resizeTuples(1);                // Allocate the memory for this data array
    targetArray->resizeTuples(newNumCellTuples); // Allocate the memory for this data array
    if(m_InterpolationType == k_NearestNeighborInterpolation)
    {
    }
    else if(m_InterpolationType == k_LinearInterpolation)
    {
      ImageRotationUtilities::ExecuteParallelFunction<ImageRotationUtilities::RotateImageGeometryWithTrilinearInterpolation>(sourceArray, taskRunner, this, this, sourceArray, targetArray, m_Params,
                                                                                                                             m_TransformationMatrix);
    }
    if(getCancel() || getErrorCode() < 0)
    {
      break;
    }
  }
  taskRunner.wait();
}

/**
 * @brief ApplyTransformationToGeometry::applyNodeGeometryTransformation
 */
void ApplyTransformationToGeometry::applyNodeGeometryTransformation()
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
  else
  {
    return;
  }
  m_TotalElements = vertexList->getNumberOfTuples();

  // Allow data-based parallelization
  ParallelDataAlgorithm dataAlg;
  dataAlg.setRange(0, m_TotalElements);
  dataAlg.execute(ApplyTransformationToGeometryImpl(*this, m_TransformationMatrix, vertexList));
}

// -----------------------------------------------------------------------------
void ApplyTransformationToGeometry::reset()
{
  m_Params = ImageRotationUtilities::RotateArgs();
}

// -----------------------------------------------------------------------------
void ApplyTransformationToGeometry::execute()
{
  dataCheck();
  if(getErrorCode() < 0)
  {
    return;
  }

  if(m_TransformationMatrixType == k_NoTransform)
  {
    return;
  }

  // Needed for Threaded Progress Messages
  IGeometry::Pointer iGeometryPtr = getDataContainerArray()->getDataContainer(getCellAttributeMatrixPath().getDataContainerName())->getGeometry();
  if(ImageGeom::Pointer image = std::dynamic_pointer_cast<ImageGeom>(iGeometryPtr))
  {
    // Function for applying Image Transformation
    applyImageGeometryTransformation();
  }
  else
  {
    applyNodeGeometryTransformation();
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ApplyTransformationToGeometry::sendThreadSafeProgressMessage(int64_t counter)
{
  std::lock_guard<std::mutex> guard(m_ProgressMessage_Mutex);

  int64_t& progCounter = ::s_ProgressValues[m_InstanceIndex];
  progCounter += counter;
  int64_t progressInt = static_cast<int64_t>((static_cast<float>(progCounter) / m_TotalElements) * 100.0f);

  // int64_t progIncrement = m_TotalElements / 100;
  int64_t& lastProgressInt = ::s_LastProgressInt[m_InstanceIndex];

  if(lastProgressInt != progressInt)
  {
    QString ss = QObject::tr("Transforming || %1% Completed").arg(progressInt);
    notifyStatusMessage(ss);
  }

  lastProgressInt = progressInt;
}

void ApplyTransformationToGeometry::sendThreadSafeProgressMessage(const QString& message)
{
  static std::mutex mutex;
  std::lock_guard<std::mutex> lock(mutex);
  qint64 currentMillis = QDateTime::currentMSecsSinceEpoch();
  if(currentMillis - m_Millis > 1000)
  {
    notifyStatusMessage(message);
    m_Millis = QDateTime::currentMSecsSinceEpoch();
  }
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

void ApplyTransformationToGeometry::setTransformationMatrixType(int value)
{
  m_TransformationMatrixType = value;
}

// -----------------------------------------------------------------------------
int ApplyTransformationToGeometry::getTransformationMatrixType() const
{
  return m_TransformationMatrixType;
}
//------------------------------------------------------------------------------
void ApplyTransformationToGeometry::setInterpolationType(int value)
{
  m_InterpolationType = value;
}

// -----------------------------------------------------------------------------
int ApplyTransformationToGeometry::getInterpolationType() const
{
  return m_InterpolationType;
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

//------------------------------------------------------------------------------
void ApplyTransformationToGeometry::setUseDataArraySelection(bool value)
{
  m_UseDataArraySelection = value;
}

//------------------------------------------------------------------------------
bool ApplyTransformationToGeometry::getUseDataArraySelection() const
{
  return m_UseDataArraySelection;
}

// -----------------------------------------------------------------------------
void ApplyTransformationToGeometry::setDataArraySelection(const std::vector<DataArrayPath>& value)
{
  m_DataArraySelection = value;
}

// -----------------------------------------------------------------------------
std::vector<DataArrayPath> ApplyTransformationToGeometry::getDataArraySelection() const
{
  return m_DataArraySelection;
}
