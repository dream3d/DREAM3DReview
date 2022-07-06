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

#include <QtCore/QTextStream>

#include "SIMPLib/Common/Constants.h"
#include "SIMPLib/Common/SIMPLRange.h"
#include "SIMPLib/DataContainers/AttributeMatrixProxy.h"
#include "SIMPLib/DataContainers/DataContainerArray.h"
#include "SIMPLib/FilterParameters/AbstractFilterParametersReader.h"
#include "SIMPLib/FilterParameters/AttributeMatrixSelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/ChoiceFilterParameter.h"
#include "SIMPLib/FilterParameters/DataArraySelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/DataContainerSelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/DynamicTableFilterParameter.h"
#include "SIMPLib/FilterParameters/FloatFilterParameter.h"
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
  value_type q = 0.0L;
  value_type c = 0.0L;
  value_type s = 0.0L;
  value_type omc = 0.0L;

  c = cos(a[3]);
  s = sin(a[3]);

  omc = static_cast<value_type>(1.0 - c);

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

  q = omc * a[0] * a[1];
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

template <class T>
void linearEquivalent(T& linEquivalent, IDataArray::Pointer linIData, int64_t linIntIndexes, double xt, double yt, double zt, const ApplyTransformationToGeometry::RotateArgs& m_Params)
{
  typename DataArray<T>::Pointer lin = std::dynamic_pointer_cast<DataArray<T>>(linIData);
  int index0 = linIntIndexes;
  int index1 = linIntIndexes + 1;
  int index2 = linIntIndexes + m_Params.xp;
  int index3 = linIntIndexes + 1 + m_Params.xp;
  int index4 = linIntIndexes + m_Params.xp * m_Params.yp;
  int index5 = linIntIndexes + 1 + m_Params.xp * m_Params.yp;
  int index6 = linIntIndexes + m_Params.xp + m_Params.xp * m_Params.yp;
  int index7 = linIntIndexes + 1 + m_Params.xp + m_Params.xp * m_Params.yp;
  if(index0 >= 0 && index0 <= m_Params.xp * m_Params.yp * m_Params.zp)
  {
    linEquivalent += lin->getPointer(0)[index0];
  }
  if(index1 >= 0 && index1 <= m_Params.xp * m_Params.yp * m_Params.zp)
  {
    linEquivalent += ((lin->getPointer(0)[index1] - lin->getPointer(0)[index0]) * xt);
  }
  if(index2 >= 0 && index2 <= m_Params.xp * m_Params.yp * m_Params.zp)
  {
    linEquivalent += ((lin->getPointer(0)[index2] - lin->getPointer(0)[index0]) * yt);
  }
  if(index3 >= 0 && index3 <= m_Params.xp * m_Params.yp * m_Params.zp)
  {
    linEquivalent += ((lin->getPointer(0)[index3] - lin->getPointer(0)[index2] - lin->getPointer(0)[index1] + lin->getPointer(0)[index0]) * xt * yt);
  }
  if(index4 >= 0 && index4 <= m_Params.xp * m_Params.yp * m_Params.zp)
  {
    linEquivalent += ((lin->getPointer(0)[index4] - lin->getPointer(0)[index0]) * zt);
  }
  if(index5 >= 0 && index5 <= m_Params.xp * m_Params.yp * m_Params.zp)
  {
    linEquivalent += ((lin->getPointer(0)[index5] - lin->getPointer(0)[index4] - lin->getPointer(0)[index1] + lin->getPointer(0)[index0]) * xt * zt);
  }
  if(index6 >= 0 && index6 <= m_Params.xp * m_Params.yp * m_Params.zp)
  {
    linEquivalent += ((lin->getPointer(0)[index6] - lin->getPointer(0)[index4] - lin->getPointer(0)[index2] + lin->getPointer(0)[index0]) * yt * zt);
  }
  if(index7 >= 0 && index7 <= m_Params.xp * m_Params.yp * m_Params.zp)
  {
    linEquivalent += ((lin->getPointer(0)[index7] - lin->getPointer(0)[index6] - lin->getPointer(0)[index5] - lin->getPointer(0)[index3] + lin->getPointer(0)[index1] + lin->getPointer(0)[index4] +
                       lin->getPointer(0)[index2] - lin->getPointer(0)[index0]) *
                      xt * yt * zt);
  }
}

template <class T>
void linearEquivalentRGB(T linEquivalent[3], IDataArray::Pointer linIData, int64_t linIntIndexes, double xt, double yt, double zt, const ApplyTransformationToGeometry::RotateArgs& m_Params)
{
  typename DataArray<T>::Pointer lin = std::dynamic_pointer_cast<DataArray<T>>(linIData);
  int index0 = linIntIndexes;
  int index1 = linIntIndexes + 1;
  int index2 = linIntIndexes + m_Params.xp;
  int index3 = linIntIndexes + 1 + m_Params.xp;
  int index4 = linIntIndexes + m_Params.xp * m_Params.yp;
  int index5 = linIntIndexes + 1 + m_Params.xp * m_Params.yp;
  int index6 = linIntIndexes + m_Params.xp + m_Params.xp * m_Params.yp;
  int index7 = linIntIndexes + 1 + m_Params.xp + m_Params.xp * m_Params.yp;
  if(index0 >= 0 && index0 <= m_Params.xp * m_Params.yp * m_Params.zp)
  {
    for(int i = 0; i < 3; i++)
    {
      linEquivalent[i] += lin->getComponent(index0, i);
    }
  }
  if(index1 >= 0 && index1 <= m_Params.xp * m_Params.yp * m_Params.zp)
  {
    for(int i = 0; i < 3; i++)
    {
      linEquivalent[i] += ((lin->getComponent(index1, i) - lin->getComponent(index0, i)) * xt);
    }
  }
  if(index2 >= 0 && index2 <= m_Params.xp * m_Params.yp * m_Params.zp)
  {
    for(int i = 0; i < 3; i++)
    {
      linEquivalent[i] += ((lin->getComponent(index2, i) - lin->getComponent(index0, i)) * yt);
    }
  }
  if(index3 >= 0 && index3 <= m_Params.xp * m_Params.yp * m_Params.zp)
  {
    for(int i = 0; i < 3; i++)
    {
      linEquivalent[i] += ((lin->getComponent(index3, i) - lin->getComponent(index2, i) - lin->getComponent(index1, i) + lin->getComponent(index0, i)) * xt * yt);
    }
  }
  if(index4 >= 0 && index4 <= m_Params.xp * m_Params.yp * m_Params.zp)
  {
    for(int i = 0; i < 3; i++)
    {
      linEquivalent[i] += ((lin->getComponent(index4, i) - lin->getComponent(index0, i)) * zt);
    }
  }
  if(index5 >= 0 && index5 <= m_Params.xp * m_Params.yp * m_Params.zp)
  {
    for(int i = 0; i < 3; i++)
    {
      linEquivalent[i] += ((lin->getComponent(index5, i) - lin->getComponent(index4, i) - lin->getComponent(index1, i) + lin->getComponent(index0, i)) * xt * zt);
    }
  }
  if(index6 >= 0 && index6 <= m_Params.xp * m_Params.yp * m_Params.zp)
  {
    for(int i = 0; i < 3; i++)
    {
      linEquivalent[i] += ((lin->getComponent(index6, i) - lin->getComponent(index4, i) - lin->getComponent(index2, i) + lin->getComponent(index0, i)) * yt * zt);
    }
  }
  if(index7 >= 0 && index7 <= m_Params.xp * m_Params.yp * m_Params.zp)
  {
    for(int i = 0; i < 3; i++)
    {
      linEquivalent[i] += ((lin->getComponent(index7, i) - lin->getComponent(index6, i) - lin->getComponent(index5, i) - lin->getComponent(index3, i) + lin->getComponent(index1, i) +
                            lin->getComponent(index4, i) + lin->getComponent(index2, i) - lin->getComponent(index0, i)) *
                           xt * yt * zt);
    }
  }
}

template <class T>
bool linearIndexes(double* LinearInterpolationData, int64_t tupleIndex, T& linEquivalent, IDataArray::Pointer linIData, const ApplyTransformationToGeometry::RotateArgs& m_Params)
{
  bool write = false;
  double xt = LinearInterpolationData[tupleIndex];
  double yt = LinearInterpolationData[tupleIndex + 1];
  double zt = LinearInterpolationData[tupleIndex + 2];
  double colOld = LinearInterpolationData[tupleIndex + 3];
  double rowOld = LinearInterpolationData[tupleIndex + 4];
  double planeOld = LinearInterpolationData[tupleIndex + 5];

  if(colOld >= 0 && colOld < m_Params.xp && colOld >= 0 && colOld < m_Params.xp && rowOld >= 0 && rowOld < m_Params.yp && planeOld >= 0 && planeOld < m_Params.zp)
  {
    int planeFloor = std::floor(planeOld);
    int rowFloor = std::floor(rowOld);
    int colFloor = std::floor(colOld);

    int64_t linIntIndexes = std::nearbyint((m_Params.xp * m_Params.yp * planeFloor) + (m_Params.xp * rowFloor) + colFloor);
    linearEquivalent<T>(linEquivalent, linIData, linIntIndexes, xt, yt, zt, m_Params);
    write = true;
  }
  return write;
}

template <class T>
bool linearIndexesRGB(double* LinearInterpolationData, int64_t tupleIndex, T linEquivalent[3], IDataArray::Pointer linIData, const ApplyTransformationToGeometry::RotateArgs& m_Params)
{
  bool write = false;

  double xt = LinearInterpolationData[tupleIndex];
  double yt = LinearInterpolationData[tupleIndex + 1];
  double zt = LinearInterpolationData[tupleIndex + 2];
  double colOld = LinearInterpolationData[tupleIndex + 3];
  double rowOld = LinearInterpolationData[tupleIndex + 4];
  double planeOld = LinearInterpolationData[tupleIndex + 5];

  if(colOld >= 0 && colOld < m_Params.xp && colOld >= 0 && colOld < m_Params.xp && rowOld >= 0 && rowOld < m_Params.yp && planeOld >= 0 && planeOld < m_Params.zp)
  {
    int planeFloor = std::floor(planeOld);
    int rowFloor = std::floor(rowOld);
    int colFloor = std::floor(colOld);

    int64_t linIntIndexes = std::nearbyint((m_Params.xp * m_Params.yp * planeFloor) + (m_Params.xp * rowFloor) + colFloor);
    linearEquivalentRGB<T>(linEquivalent, linIData, linIntIndexes, xt, yt, zt, m_Params);
    write = true;
  }
  return write;
}

template <typename T>
void wrapLinearIndexes(double* LinearInterpolationData, int64_t tupleIndex, typename DataArray<T>::Pointer lin, typename IDataArray::Pointer linData,
                       const ApplyTransformationToGeometry::RotateArgs& m_Params)
{
  bool wrote = false;
  int index = tupleIndex / 6;

  T linEquivalent = 0;
  typename IDataArray::Pointer linIData = std::dynamic_pointer_cast<IDataArray>(lin);
  wrote = linearIndexes<T>(LinearInterpolationData, tupleIndex, linEquivalent, linIData, m_Params);
  if(wrote)
  {
    linData->initializeTuple(index, &linEquivalent);
  }
  else
  {
    int var = 0;
    linData->initializeTuple(index, &var);
  }
}

template <typename T>
void wrapLinearIndexesRGB(double* LinearInterpolationData, int64_t tupleIndex, typename DataArray<T>::Pointer lin, typename IDataArray::Pointer linData,
                          const ApplyTransformationToGeometry::RotateArgs& m_Params)
{
  bool wrote = false;
  int index = tupleIndex / 6;
  T linEquivalent[3] = {0, 0, 0};
  typename IDataArray::Pointer linIData = std::dynamic_pointer_cast<IDataArray>(lin);
  wrote = linearIndexesRGB<T>(LinearInterpolationData, tupleIndex, linEquivalent, linIData, m_Params);
  if(wrote)
  {
    linData->initializeTuple(index, &linEquivalent);
  }
  else
  {
    int var = 0;
    linData->initializeTuple(index, &var);
  }
}

template <typename T>
bool applyLinearInterpolation(typename DataArray<T>::Pointer lin, int64_t index, int64_t tupleIndex, double* LinearInterpolationData, typename IDataArray::Pointer linData, bool RGB,
                              const ApplyTransformationToGeometry::RotateArgs& m_Params)
{
  if(!lin)
  {
    return false;
  }

  if(RGB)
  {
    wrapLinearIndexesRGB<T>(LinearInterpolationData, tupleIndex, lin, linData, m_Params);
  }
  else
  {
    wrapLinearIndexes<T>(LinearInterpolationData, tupleIndex, lin, linData, m_Params);
  }
  return true;
}

} // end anonymous namespace

namespace ApplyTransformationProgress
{

} // namespace ApplyTransformationProgress

const Eigen::Vector3f k_XAxis = Eigen::Vector3f::UnitX();
const Eigen::Vector3f k_YAxis = Eigen::Vector3f::UnitY();
const Eigen::Vector3f k_ZAxis = Eigen::Vector3f::UnitZ();

// Function for determining new ImageGeom dimensions after transformation
void determineMinMax(const ApplyTransformationToGeometry::Matrix3fR& rotationMatrix, const FloatVec3Type& spacing, size_t col, size_t row, size_t plane, float& xMin, float& xMax, float& yMin,
                     float& yMax, float& zMin, float& zMax)
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

// Determines paramaters for image rotation
ApplyTransformationToGeometry::RotateArgs createRotateParams(const ImageGeom& imageGeom, const ApplyTransformationToGeometry::Transform3f transformationMatrix)
{
  const SizeVec3Type origDims = imageGeom.getDimensions();
  const FloatVec3Type spacing = imageGeom.getSpacing();

  ApplyTransformationToGeometry::Matrix3fR rotationMatrix = ApplyTransformationToGeometry::Matrix3fR::Zero();
  ApplyTransformationToGeometry::Matrix3fR scaleMatrix = ApplyTransformationToGeometry::Matrix3fR::Zero();

  transformationMatrix.computeRotationScaling(&rotationMatrix, &scaleMatrix);

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

  float xResNew = determineSpacing(spacing, xAxisNew);
  float yResNew = determineSpacing(spacing, yAxisNew);
  float zResNew = determineSpacing(spacing, zAxisNew);

  MeshIndexType xpNew = static_cast<int64_t>(std::nearbyint((xMax - xMin) / xResNew) + 1);
  MeshIndexType ypNew = static_cast<int64_t>(std::nearbyint((yMax - yMin) / yResNew) + 1);
  MeshIndexType zpNew = static_cast<int64_t>(std::nearbyint((zMax - zMin) / zResNew) + 1);

  ApplyTransformationToGeometry::RotateArgs params;

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

// Alters image parameters, scales and translates
void updateGeometry(ImageGeom& imageGeom, const ApplyTransformationToGeometry::RotateArgs& params, const ApplyTransformationToGeometry::Matrix3fR& scalingMatrix,
                    const ApplyTransformationToGeometry::Matrix3fR& rotationMatrix, const ApplyTransformationToGeometry::MatrixTranslation translationMatrix)
{
  float m_ScalingMatrix[3][3] = {{0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f}};
  float m_TranslationMatrix[3] = {0.0f, 0.0f, 0.0f};
  Eigen::Map<ApplyTransformationToGeometry::Matrix3fR>(&m_ScalingMatrix[0][0], scalingMatrix.rows(), scalingMatrix.cols()) = scalingMatrix;

  Eigen::Map<ApplyTransformationToGeometry::MatrixTranslation>(m_TranslationMatrix, translationMatrix.rows(), translationMatrix.cols()) = translationMatrix;
  FloatVec3Type origin = imageGeom.getOrigin();

  Eigen::Vector3f original_translation(m_TranslationMatrix[0], m_TranslationMatrix[1], m_TranslationMatrix[2]);

  Eigen::Vector3f original_origin(origin[0], origin[1], origin[2]);
  Eigen::Vector3f original_origin_rot = rotationMatrix * original_origin;

  // Applies Scaling to Image
  imageGeom.setSpacing(params.xResNew * m_ScalingMatrix[0][0], params.yResNew * m_ScalingMatrix[1][1], params.zResNew * m_ScalingMatrix[2][2]);

  imageGeom.setDimensions(params.xpNew, params.ypNew, params.zpNew);

  // Applies Translation to Image
  origin[0] = params.xMinNew * m_ScalingMatrix[0][0] + original_translation[0] + original_origin_rot[0] * params.xResNew * m_ScalingMatrix[0][0] / params.xRes;
  origin[1] = params.yMinNew * m_ScalingMatrix[1][1] + original_translation[1] + original_origin_rot[1] * params.yResNew * m_ScalingMatrix[1][1] / params.yRes;
  origin[2] = params.zMinNew * m_ScalingMatrix[2][2] + original_translation[2] + original_origin_rot[2] * params.zResNew * m_ScalingMatrix[2][2] / params.zRes;
  imageGeom.setOrigin(origin);
}

/**
 * @brief The RotateSampleRefFrameImpl class implements a threaded algorithm to do the
 * actual computation of the rotation by applying the rotation to each Euler angle
 */
class SampleRefFrameRotator
{
  DataArray<int64_t>::Pointer m_NewIndicesPtr;
  DataArray<double>::Pointer m_LinearInterpolationDataPtr;
  float m_RotMatrixInv[3][3] = {{0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f}};
  ApplyTransformationToGeometry::RotateArgs m_Params;
  int interpolationType;

public:
  SampleRefFrameRotator(DataArray<int64_t>::Pointer newindices, DataArray<double>::Pointer linearInterpolationDataPtr, const ApplyTransformationToGeometry::RotateArgs& args,
                        const ApplyTransformationToGeometry::Matrix3fR& rotationMatrix, int interpType)
  : m_NewIndicesPtr(newindices)
  , m_LinearInterpolationDataPtr(linearInterpolationDataPtr)
  , m_Params(args)
  , interpolationType(interpType)
  {
    // We have to inline the 3x3 Maxtrix transpose here because of the "const" nature of the 'convert' function
    ApplyTransformationToGeometry::Matrix3fR transpose = rotationMatrix.transpose();
    // Need to use row based Eigen matrix so that the values get mapped to the right place in the raw array
    // Raw array is faster than Eigen
    Eigen::Map<ApplyTransformationToGeometry::Matrix3fR>(&m_RotMatrixInv[0][0], transpose.rows(), transpose.cols()) = transpose;
  }

  ~SampleRefFrameRotator() = default;

  void convert(int64_t zStart, int64_t zEnd, int64_t yStart, int64_t yEnd, int64_t xStart, int64_t xEnd) const
  {
    int64_t* newindicies = m_NewIndicesPtr->getPointer(0);
    double* linearInterpolationDataPtr = m_LinearInterpolationDataPtr->getPointer(0);

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
          linearInterpolationDataPtr[index * 6] = -1;
          linearInterpolationDataPtr[index * 6 + 1] = -1;
          linearInterpolationDataPtr[index * 6 + 2] = -1;
          linearInterpolationDataPtr[index * 6 + 3] = -1;
          linearInterpolationDataPtr[index * 6 + 4] = -1;
          linearInterpolationDataPtr[index * 6 + 5] = -1;

          float coords[3] = {0.0f, 0.0f, 0.0f};
          float coordsNew[3] = {0.0f, 0.0f, 0.0f};

          coords[0] = (static_cast<float>(i) * m_Params.xResNew) + m_Params.xMinNew;
          coords[1] = (static_cast<float>(j) * m_Params.yResNew) + m_Params.yMinNew;
          coords[2] = (static_cast<float>(k) * m_Params.zResNew) + m_Params.zMinNew;

          MatrixMath::Multiply3x3with3x1(m_RotMatrixInv, coords, coordsNew);

          double x0 = static_cast<float>(std::floor(coordsNew[0] / m_Params.xRes));
          double x1 = static_cast<float>(std::ceil(coordsNew[0] / m_Params.xRes));
          double y0 = static_cast<float>(std::floor(coordsNew[1] / m_Params.yRes));
          double y1 = static_cast<float>(std::ceil(coordsNew[1] / m_Params.yRes));
          double z0 = static_cast<float>(std::floor(coordsNew[2] / m_Params.zRes));
          double z1 = static_cast<float>(std::ceil(coordsNew[2] / m_Params.zRes));

          int64_t colOld = static_cast<int64_t>(std::nearbyint(coordsNew[0] / m_Params.xRes));
          int64_t rowOld = static_cast<int64_t>(std::nearbyint(coordsNew[1] / m_Params.yRes));
          int64_t planeOld = static_cast<int64_t>(std::nearbyint(coordsNew[2] / m_Params.zRes));

          if(colOld >= 0 && colOld < m_Params.xp && colOld >= 0 && colOld < m_Params.xp && rowOld >= 0 && rowOld < m_Params.yp && planeOld >= 0 && planeOld < m_Params.zp)
          {
            newindicies[index] = ((m_Params.xp * m_Params.yp * planeOld) + (m_Params.xp * rowOld) + colOld);
          }

          if(interpolationType == 1)
          {

            double colOld = static_cast<double>(coordsNew[0] / m_Params.xRes);
            double rowOld = static_cast<double>(coordsNew[1] / m_Params.yRes);
            double planeOld = static_cast<double>(coordsNew[2] / m_Params.zRes);

            double xt = 0.5;
            double yt = 0.5;
            double zt = 0.5;

            if(x0 == x1)
            {
              xt = 0;
            }
            else
            {
              xt = (colOld - x0) / (x1 - x0);
            }

            if(y0 == y1)
            {
              yt = 0;
            }
            else
            {
              yt = (rowOld - y0) / (y1 - y0);
            }
            if(z0 == z1)
            {
              zt = 0;
            }
            else
            {
              zt = (planeOld - z0) / (z1 - z0);
            }

            if(colOld >= 0 && colOld < m_Params.xp && colOld >= 0 && colOld < m_Params.xp && rowOld >= 0 && rowOld < m_Params.yp && planeOld >= 0 && planeOld < m_Params.zp)
            {
              linearInterpolationDataPtr[index * 6] = xt;
              linearInterpolationDataPtr[index * 6 + 1] = yt;
              linearInterpolationDataPtr[index * 6 + 2] = zt;
              linearInterpolationDataPtr[index * 6 + 3] = colOld;
              linearInterpolationDataPtr[index * 6 + 4] = rowOld;
              linearInterpolationDataPtr[index * 6 + 5] = planeOld;
            }
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
    parameters.push_back(SIMPL_NEW_DYN_TABLE_FP("Transformation Matrix", ManualTransformationMatrix, FilterParameter::Category::Parameter, ApplyTransformationToGeometry, 2));
  }
  parameters.push_back(SIMPL_NEW_FLOAT_FP("Rotation Angle (Degrees)", RotationAngle, FilterParameter::Category::Parameter, ApplyTransformationToGeometry, 3));
  parameters.push_back(SIMPL_NEW_FLOAT_VEC3_FP("Rotation Axis (ijk)", RotationAxis, FilterParameter::Category::Parameter, ApplyTransformationToGeometry, 3));
  parameters.push_back(SIMPL_NEW_FLOAT_VEC3_FP("Translation", Translation, FilterParameter::Category::Parameter, ApplyTransformationToGeometry, 4));
  parameters.push_back(SIMPL_NEW_FLOAT_VEC3_FP("Scale", Scale, FilterParameter::Category::Parameter, ApplyTransformationToGeometry, 5));

  {
    DataArraySelectionFilterParameter::RequirementType datReq =
        DataArraySelectionFilterParameter::CreateRequirement(SIMPL::TypeNames::Float, SIMPL::Defaults::AnyComponentSize, AttributeMatrix::Type::Generic, IGeometry::Type::Any);
    parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Transformation Matrix", ComputedTransformationMatrix, FilterParameter::Category::RequiredArray, ApplyTransformationToGeometry, datReq, 1));
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
  using ProjectiveMatrix = Eigen::Matrix<float, 4, 4, Eigen::RowMajor>;
  // using RotateMatrix = Eigen::Matrix<float, 3, 3, Eigen::RowMajor>;

  clearErrorCode();
  clearWarningCode();

  reset();

  IGeometry::Pointer igeom = getDataContainerArray()->getPrereqGeometryFromDataContainer<IGeometry>(this, getCellAttributeMatrixPath().getDataContainerName());

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
    using OrientationF = std::vector<float>;
    OrientationF om = OrientationTransformation::ax2om<OrientationF, OrientationF>(OrientationF({m_RotationAxis[0], m_RotationAxis[1], m_RotationAxis[2], rotAngle}));

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

  // if ImageGeom found:
  if(std::dynamic_pointer_cast<ImageGeom>(igeom) && m_TransformationMatrix != nullptr)
  {
    DataContainer::Pointer m = getDataContainerArray()->getDataContainer(getCellAttributeMatrixPath().getDataContainerName());
    QString attrMatName = getCellAttributeMatrixPath().getAttributeMatrixName();
    QList<QString> voxelArrayNames = m->getAttributeMatrix(attrMatName)->getAttributeArrayNames();

    if(getInterpolationType() == 1)
    {
      if(m_UseDataArraySelection)
      {
        voxelArrayNames.clear();
        for(const auto& dataArrayPath : m_DataArraySelection)
        {
          voxelArrayNames.append(dataArrayPath.getDataArrayName());
        }
      }

      for(const auto& attrArrayName : voxelArrayNames)
      {
        IDataArray::Pointer p = m->getAttributeMatrix(attrMatName)->getAttributeArray(attrArrayName);
        if(p->getTypeAsString().compare("bool") == 0)
        {
          QString ss = QObject::tr("Input Type Error, cannot run linear interpolation on a boolean array");
          setErrorCondition(-704, ss);
        }
      }
    }
    ImageGeom::Pointer imageGeom = std::dynamic_pointer_cast<ImageGeom>(igeom);
    Eigen::Map<ProjectiveMatrix> transformation(m_TransformationMatrix);

    // Column Major can convert if needed
    Eigen::Transform<float, 3, Eigen::Affine> transform = Eigen::Transform<float, 3, Eigen::Affine>(transformation);
    ApplyTransformationToGeometry::Matrix3fR rotationMatrix = ApplyTransformationToGeometry::Matrix3fR::Zero();
    ApplyTransformationToGeometry::Matrix3fR scaleMatrix = ApplyTransformationToGeometry::Matrix3fR::Zero();
    ApplyTransformationToGeometry::MatrixTranslation translationMatrix = ApplyTransformationToGeometry::MatrixTranslation::Zero();

    transform.computeRotationScaling(&rotationMatrix, &scaleMatrix);

    translationMatrix(0, 0) = transform.data()[12];
    translationMatrix(0, 1) = transform.data()[13];
    translationMatrix(0, 2) = transform.data()[14];

    m_RotationMatrix = rotationMatrix;
    m_ScalingMatrix = scaleMatrix;
    m_TranslationMatrix = translationMatrix;

    m_Params = ::createRotateParams(*imageGeom, transform);
    ::updateGeometry(*imageGeom, m_Params, scaleMatrix, rotationMatrix, translationMatrix);

    // Resize attribute matrix
    std::vector<size_t> tDims(3);
    tDims[0] = m_Params.xpNew;
    tDims[1] = m_Params.ypNew;
    tDims[2] = m_Params.zpNew;

    m->getAttributeMatrix(attrMatName)->resizeAttributeArrays(tDims);
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------

void ApplyTransformationToGeometry::ApplyImageTransformation()
{
  DataContainer::Pointer m = getDataContainerArray()->getDataContainer(getCellAttributeMatrixPath().getDataContainerName());
  int64_t newNumCellTuples = m_Params.xpNew * m_Params.ypNew * m_Params.zpNew;
  int64_t newNumCellTuplesLinData = newNumCellTuples * 6;

  QString name = "_INTERNAL_USE_ONLY_RotateSampleRef_LinearInterpolationData";
  DataArray<int64_t>::Pointer newIndiciesPtr = DataArray<int64_t>::CreateArray(newNumCellTuples, std::string("_INTERNAL_USE_ONLY_RotateSampleRef_NewIndicies"), true);
  DataArray<double>::Pointer LinearInterpolationDataPtr = DataArray<double>::CreateArray(newNumCellTuplesLinData, name, true);
  newIndiciesPtr->initializeWithValue(-1);
  LinearInterpolationDataPtr->initializeWithValue(-1);
  int64_t* newindicies = newIndiciesPtr->getPointer(0);
  double* LinearInterpolationData = LinearInterpolationDataPtr->getPointer(0);

#ifdef SIMPL_USE_PARALLEL_ALGORITHMS
  tbb::parallel_for(tbb::blocked_range3d<int64_t, int64_t, int64_t>(0, m_Params.zpNew, 0, m_Params.ypNew, 0, m_Params.xpNew),
                    ::SampleRefFrameRotator(newIndiciesPtr, LinearInterpolationDataPtr, m_Params, m_RotationMatrix, m_InterpolationType), tbb::auto_partitioner());
#else
  {
    SampleRefFrameRotator serial(newIndiciesPtr, LinearInterpolationDataPtr, m_Params, m_RotationMatrix, m_InterpolationType);
    serial.convert(0, m_Params.zpNew, 0, m_Params.ypNew, 0, m_Params.xpNew);
  }
#endif

  QString attrMatName = getCellAttributeMatrixPath().getAttributeMatrixName();
  QList<QString> voxelArrayNames = m->getAttributeMatrix(attrMatName)->getAttributeArrayNames();

  if(m_UseDataArraySelection)
  {
    voxelArrayNames.clear();
    for(const auto& dataArrayPath : m_DataArraySelection)
    {
      voxelArrayNames.append(dataArrayPath.getDataArrayName());
    }
  }

  for(const auto& attrArrayName : voxelArrayNames)
  {
    IDataArray::Pointer p = m->getAttributeMatrix(attrMatName)->getAttributeArray(attrArrayName);

    // Make a copy of the 'p' array that has the same name. When placed into
    // the data container this will over write the current array with
    // the same name.
    IDataArray::Pointer data = p->createNewArray(newNumCellTuples, p->getComponentDimensions(), p->getName());
    IDataArray::Pointer linData = p->createNewArray(newNumCellTuples, p->getComponentDimensions(), p->getName());
    IDataArray::Pointer linPtr = p->createNewArray(p->getNumberOfTuples(), p->getComponentDimensions(), p->getName());

    bool RGB = false;

    if(p->getNumberOfComponents() == 3)
    {
      RGB = true;
    }

    if(!linPtr->copyFromArray(0, p, 0, p->getNumberOfTuples()))
    {
      QString ss = QObject::tr("copyFromArray Failed: linPtr");
      QTextStream out(&ss);
      setErrorCondition(-45102, ss);
      return;
    }

    int64_t newIndicies_I = 0;
    if(m_InterpolationType == 0)
    {
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
    else if(m_InterpolationType == 1)
    {
      for(size_t i = 0; i < static_cast<size_t>(newNumCellTuples); i++)
      {
        if(i >= 0)
        {
          int64_t tupleIndex = i * 6;

          if(DataArray<int8_t>::Pointer lin = std::dynamic_pointer_cast<DataArray<int8_t>>(linPtr))
          {
            applyLinearInterpolation<int8_t>(lin, i, tupleIndex, LinearInterpolationData, linData, RGB, m_Params);
          }
          else if(DataArray<uint8_t>::Pointer lin = std::dynamic_pointer_cast<DataArray<uint8_t>>(linPtr))
          {
            applyLinearInterpolation<uint8_t>(lin, i, tupleIndex, LinearInterpolationData, linData, RGB, m_Params);
          }
          else if(DataArray<int16_t>::Pointer lin = std::dynamic_pointer_cast<DataArray<int16_t>>(linPtr))
          {
            applyLinearInterpolation<int16_t>(lin, i, tupleIndex, LinearInterpolationData, linData, RGB, m_Params);
          }
          else if(DataArray<uint16_t>::Pointer lin = std::dynamic_pointer_cast<DataArray<uint16_t>>(linPtr))
          {
            applyLinearInterpolation<uint16_t>(lin, i, tupleIndex, LinearInterpolationData, linData, RGB, m_Params);
          }
          else if(DataArray<int32_t>::Pointer lin = std::dynamic_pointer_cast<DataArray<int32_t>>(linPtr))
          {
            applyLinearInterpolation<int32_t>(lin, i, tupleIndex, LinearInterpolationData, linData, RGB, m_Params);
          }
          else if(DataArray<uint32_t>::Pointer lin = std::dynamic_pointer_cast<DataArray<uint32_t>>(linPtr))
          {
            applyLinearInterpolation<uint32_t>(lin, i, tupleIndex, LinearInterpolationData, linData, RGB, m_Params);
          }
          else if(DataArray<int64_t>::Pointer lin = std::dynamic_pointer_cast<DataArray<int64_t>>(linPtr))
          {
            applyLinearInterpolation<int64_t>(lin, i, tupleIndex, LinearInterpolationData, linData, RGB, m_Params);
          }
          else if(DataArray<uint64_t>::Pointer lin = std::dynamic_pointer_cast<DataArray<uint64_t>>(linPtr))
          {
            applyLinearInterpolation<uint64_t>(lin, i, tupleIndex, LinearInterpolationData, linData, RGB, m_Params);
          }
          else if(DataArray<float>::Pointer lin = std::dynamic_pointer_cast<DataArray<float>>(linPtr))
          {
            applyLinearInterpolation<float>(lin, i, tupleIndex, LinearInterpolationData, linData, RGB, m_Params);
          }
          else if(DataArray<double>::Pointer lin = std::dynamic_pointer_cast<DataArray<double>>(linPtr))
          {
            applyLinearInterpolation<double>(lin, i, tupleIndex, LinearInterpolationData, linData, RGB, m_Params);
          }
          else
          {
            QString ss = QObject::tr("Casted Array Linear Interpolation Failed");
            QTextStream out(&ss);
            setErrorCondition(-45102, ss);
            return;
          }
        }
        else
        {
          int var = 0;
          linData->initializeTuple(i, &var);
        }
      }
      m->getAttributeMatrix(attrMatName)->insertOrAssign(linData);
    }
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
  else if(ImageGeom::Pointer image = std::dynamic_pointer_cast<ImageGeom>(igeom))
  {
    // Function for applying Image Transformation
    ApplyImageTransformation();
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
void ApplyTransformationToGeometry::reset()
{
  m_RotationMatrix.setZero();

  m_Params = ApplyTransformationToGeometry::RotateArgs();
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
  m_InstanceIndex = ++::s_InstanceIndex;
  ::s_ProgressValues[m_InstanceIndex] = 0;
  ::s_LastProgressInt[m_InstanceIndex] = 0;

  applyTransformation();
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

  int64_t progIncrement = m_TotalElements / 100;
  int64_t prog = 1;

  int64_t& lastProgressInt = ::s_LastProgressInt[m_InstanceIndex];

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
