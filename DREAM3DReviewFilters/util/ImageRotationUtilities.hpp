#pragma once

#include "SIMPLib/Common/SIMPLArray.hpp"
#include "SIMPLib/DataArrays/DataArray.hpp"
#include "SIMPLib/Filtering/AbstractFilter.h"
#include "SIMPLib/Geometry/ImageGeom.h"

#include <Eigen/Dense>

#include <fstream>
#include <iostream>

using Matrix3fR = Eigen::Matrix<float, 3, 3, Eigen::RowMajor>;
using Matrix4fR = Eigen::Matrix<float, 4, 4, Eigen::RowMajor>;

// typedef Array<Scalar, RowsAtCompileTime, ColsAtCompileTime, Options> MyArrayType;

using Vector3s = Eigen::Array<size_t, 1, 3>;
using Vector3i64 = Eigen::Array<int64_t, 1, 3>;

using Int64Vec3Type = IVec3<int64_t>;

namespace ImageRotationUtilities
{

struct RotateArgs
{
  ImageGeom::Pointer origImageGeom;

  int64_t xp = 0;
  int64_t yp = 0;
  int64_t zp = 0;
  float xRes = 0.0f;
  float yRes = 0.0f;
  float zRes = 0.0f;

  ImageGeom::Pointer transformedImageGeom;
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

template <typename T>
T inline GetSourceArrayValue(const RotateArgs& params, Vector3i64 xyzIndex, DataArray<T>& sourceArray, size_t compIndex)
{
  if(xyzIndex[0] < 0)
  {
    xyzIndex[0] = 0;
  }
  if(xyzIndex[0] >= params.xp)
  {
    xyzIndex[0] = params.xp - 1;
  }

  if(xyzIndex[1] < 0)
  {
    xyzIndex[1] = 0;
  }
  if(xyzIndex[1] >= params.yp)
  {
    xyzIndex[1] = params.yp - 1;
  }

  if(xyzIndex[2] < 0)
  {
    xyzIndex[2] = 0;
  }
  if(xyzIndex[2] >= params.zp)
  {
    xyzIndex[2] = params.zp - 1;
  }

  // Now just compute the proper index
  size_t index = (xyzIndex[2] * params.xp * params.yp) + (xyzIndex[1] * params.xp) + xyzIndex[0];
  return sourceArray[index * sourceArray.getNumberOfComponents() + compIndex];
}

/**
 * @brief FindOctant
 * @param params
 * @param index
 * @param coord
 * @return
 */
inline size_t FindOctant(const RotateArgs& params, size_t index, FloatVec3Type coord)
{

  float xResHalf = params.xRes * 0.5;
  float yResHalf = params.yRes * 0.5;
  float zResHalf = params.zRes * 0.5;

  // Get the center coord of the original source voxel
  Eigen::Vector3f centerPoint;
  params.origImageGeom->getCoords(index, centerPoint.data());

  // Form the 8 corner coords for the voxel
  // clang-format off
  std::array<FloatVec3Type, 8> unitSquareCoords = {
  /* P1 */ FloatVec3Type(centerPoint[0]-xResHalf, centerPoint[1]-yResHalf, centerPoint[2]-zResHalf),
  /* P2 */ FloatVec3Type(centerPoint[0]+xResHalf, centerPoint[1]-yResHalf, centerPoint[2]-zResHalf),
  /* P3 */ FloatVec3Type(centerPoint[0]+xResHalf, centerPoint[1]+yResHalf, centerPoint[2]-zResHalf),
  /* P4 */ FloatVec3Type(centerPoint[0]-xResHalf, centerPoint[1]+yResHalf, centerPoint[2]-zResHalf),
  /* P5 */ FloatVec3Type(centerPoint[0]-xResHalf, centerPoint[1]-yResHalf, centerPoint[2]+zResHalf),
  /* P6 */ FloatVec3Type(centerPoint[0]+xResHalf, centerPoint[1]-yResHalf, centerPoint[2]+zResHalf),
  /* P7 */ FloatVec3Type(centerPoint[0]+xResHalf, centerPoint[1]+yResHalf, centerPoint[2]+zResHalf),
  /* P8 */ FloatVec3Type(centerPoint[0]-xResHalf, centerPoint[1]+yResHalf, centerPoint[2]+zResHalf),
  };
  // clang-format on

  FloatVec3Type temp;
  // Now figure out which corner the inverse transformed point is closest to
  // this will give us which octant the point lies.
  float minDistance = std::numeric_limits<float>::max();
  size_t minIndex = 0;
  for(size_t i = 0; i < 8; i++)
  {
    float distance = (unitSquareCoords[i] - coord).norm();
    if(distance < minDistance)
    {
      minDistance = distance;
      minIndex = i;
    }
  }

  return minIndex;
}

using OctantOffsetArrayType = std::array<Vector3i64, 8>;

static const OctantOffsetArrayType k_IndexOffset0 = {Vector3i64{-1, -1, -1}, Vector3i64{0, -1, -1}, Vector3i64{0, 0, -1}, Vector3i64{-1, 0, -1},
                                                     Vector3i64{-1, -1, 0},  Vector3i64{0, -1, 0},  Vector3i64{0, 0, 0},  Vector3i64{-1, 0, 0}};
static const OctantOffsetArrayType k_IndexOffset1 = {Vector3i64{0, -1, -1}, Vector3i64{1, -1, -1}, Vector3i64{1, 0, -1}, Vector3i64{0, 0, -1},
                                                     Vector3i64{0, -1, 0},  Vector3i64{1, -1, 0},  Vector3i64{1, 0, 0},  Vector3i64{0, 0, 0}};
static const OctantOffsetArrayType k_IndexOffset2 = {Vector3i64{0, 0, -1}, Vector3i64{1, 0, -1}, Vector3i64{1, 1, -1}, Vector3i64{0, 1, -1},
                                                     Vector3i64{0, 0, 0},  Vector3i64{1, 0, 0},  Vector3i64{1, 1, 0},  Vector3i64{0, 1, 0}};
static const OctantOffsetArrayType k_IndexOffset3 = {Vector3i64{-1, 0, -1}, Vector3i64{0, 0, -1}, Vector3i64{0, 1, -1}, Vector3i64{-1, 1, -1},
                                                     Vector3i64{-1, 0, 0},  Vector3i64{0, 0, 0},  Vector3i64{0, 1, 0},  Vector3i64{-1, 1, 0}};
static const OctantOffsetArrayType k_IndexOffset4 = {Vector3i64{-1, -1, 0}, Vector3i64{0, -1, 0}, Vector3i64{0, 0, 0}, Vector3i64{-1, 0, 0},
                                                     Vector3i64{-1, -1, 1}, Vector3i64{0, -1, 1}, Vector3i64{0, 0, 1}, Vector3i64{-1, 0, 1}};
static const OctantOffsetArrayType k_IndexOffset5 = {Vector3i64{0, -1, 0}, Vector3i64{1, -1, 0}, Vector3i64{1, 0, 0}, Vector3i64{0, 0, 0},
                                                     Vector3i64{0, -1, 1}, Vector3i64{1, -1, 1}, Vector3i64{1, 0, 1}, Vector3i64{0, 0, 1}};
static const OctantOffsetArrayType k_IndexOffset6 = {Vector3i64{0, 0, 0}, Vector3i64{1, 0, 0}, Vector3i64{1, 1, 0}, Vector3i64{0, 1, 0},
                                                     Vector3i64{0, 0, 1}, Vector3i64{1, 0, 1}, Vector3i64{1, 1, 1}, Vector3i64{0, 1, 1}};
static const OctantOffsetArrayType k_IndexOffset7 = {Vector3i64{-1, 0, 0}, Vector3i64{0, 0, 0}, Vector3i64{0, 1, 0}, Vector3i64{-1, -1, 0},
                                                     Vector3i64{-1, 0, 1}, Vector3i64{0, 0, 1}, Vector3i64{0, 1, 1}, Vector3i64{-1, -1, 1}};
static const std::array<OctantOffsetArrayType, 8> k_AllOctantOffsets{k_IndexOffset0, k_IndexOffset1, k_IndexOffset2, k_IndexOffset3, k_IndexOffset4, k_IndexOffset5, k_IndexOffset6, k_IndexOffset7};

/**
 * @brief FindInterpolationValues
 * @param params
 * @param index
 * @param octant
 * @param oldIndicesU
 * @param oldCoords
 * @param sourceArray
 * @param pValues
 * @param uvw
 */
template <typename T>
inline void FindInterpolationValues(const RotateArgs& params, size_t index, size_t octant, Vector3s oldIndicesU, Eigen::Array4f& oldCoords, DataArray<T>& sourceArray, std::vector<T>& pValues,
                                    Eigen::Vector3f& uvw)
{
  const std::array<Vector3i64, 8>& indexOffset = k_AllOctantOffsets[octant];

  Vector3i64 oldIndices(static_cast<int64_t>(oldIndicesU[0]), static_cast<int64_t>(oldIndicesU[1]), static_cast<int64_t>(oldIndicesU[2]));
  size_t numComps = sourceArray.getNumberOfComponents();

  Vector3i64 pIndices;
  Eigen::Vector3f p1Coord;

  auto origin = params.origImageGeom->getOrigin();

  for(size_t i = 0; i < 8; i++)
  {
    auto pIndices = oldIndices + indexOffset[i];
    for(size_t compIndex = 0; compIndex < numComps; compIndex++)
    {
      T value = GetSourceArrayValue<T>(params, pIndices, sourceArray, compIndex);
      pValues[i * numComps + compIndex] = value;
    }
    if(i == 0)
    {
      p1Coord = {pIndices[0] * params.xRes + (0.5F * params.xRes) + origin[0], pIndices[1] * params.yRes + (0.5F * params.yRes) + origin[1],
                 pIndices[2] * params.zRes + (0.5F * params.zRes) + origin[2]};
    }
  }
  uvw[0] = oldCoords[0] - p1Coord[0];
  uvw[1] = oldCoords[1] - p1Coord[1];
  uvw[2] = oldCoords[2] - p1Coord[2];
}

/**
 * @brief The RotateImageGeometryWithTrilinearInterpolation class
 */
template <typename T, class FilterType>
class RotateImageGeometryWithTrilinearInterpolation
{
private:
  FilterType* m_Filter = nullptr;
  IDataArray::Pointer m_SourceArray;
  IDataArray::Pointer m_TargetArray;
  ImageRotationUtilities::RotateArgs m_Params;
  Matrix4fR m_TransformationMatrix;
public:
  RotateImageGeometryWithTrilinearInterpolation(FilterType* filter, IDataArray::Pointer& sourceArray, IDataArray::Pointer& targetArray, RotateArgs& rotateArgs, Matrix4fR& transformationMatrix)
  : m_Filter(filter)
  , m_SourceArray(sourceArray)
  , m_TargetArray(targetArray)
  , m_Params(rotateArgs)
  , m_TransformationMatrix(transformationMatrix)
  {
  }

  /**
   * @brief CalculateInterpolatedValue
   *
   * This comes from https://www.cs.purdue.edu/homes/cs530/slides/04.DataStructure.pdf, page 36.
   * @param sourceArray
   * @param oldIndex
   * @param indices
   * @return
   */
  T CalculateInterpolatedValue(std::vector<T>& pValues, Eigen::Vector3f& uvw, size_t numComps, size_t compIndex) const
  {
    constexpr size_t P1 = 0;
    constexpr size_t P2 = 1;
    constexpr size_t P3 = 2;
    constexpr size_t P4 = 3;
    constexpr size_t P5 = 4;
    constexpr size_t P6 = 5;
    constexpr size_t P7 = 6;
    constexpr size_t P8 = 7;

    const float u = uvw[0];
    const float v = uvw[1];
    const float w = uvw[2];

    T value = pValues[0];
    // clang-format off
    value += u * (pValues[P2 * numComps + compIndex] - pValues[P1 * numComps + compIndex]);
    value += v * (pValues[P4 * numComps + compIndex] - pValues[P1 * numComps + compIndex]);
    value += w * (pValues[P5 * numComps + compIndex] - pValues[P1 * numComps + compIndex]);
    value += u * v * (pValues[P1 * numComps + compIndex] - pValues[P2 * numComps + compIndex] + pValues[P3 * numComps + compIndex] - pValues[P4 * numComps + compIndex]);
    value += u * w * (pValues[P1 * numComps + compIndex] - pValues[P2 * numComps + compIndex] - pValues[P5 * numComps + compIndex] + pValues[P6 * numComps + compIndex]);
    value += v * w * (pValues[P1 * numComps + compIndex] - pValues[P4 * numComps + compIndex] - pValues[P5 * numComps + compIndex] + pValues[P8 * numComps + compIndex]);
    value += u * v * w *
             (-pValues[P1 * numComps + compIndex] + pValues[P2 * numComps + compIndex] - pValues[P3 * numComps + compIndex] + pValues[P4 * numComps + compIndex] -
              pValues[P5 * numComps + compIndex] + pValues[P6 * numComps + compIndex] - pValues[P7 * numComps + compIndex] + pValues[P8 * numComps + compIndex]);

    // clang-format on
    return value;
  }

  /**
   * @brief This is the main algorithm to perform the interpolation and get a final value that is placed into the transformed
   * voxel. This uses Trilinear interpolation which will devolve into Bilinear and Linear interpolation depending on the
   * values of U, V and W.
   */
  void operator()() const
  {
    using DataArrayType = DataArray<T>;
    using DataArrayPointerType = typename DataArrayType::Pointer;

    DataArrayPointerType sourceArrayPtr = std::dynamic_pointer_cast<DataArrayType>(m_SourceArray);
    DataArrayType& sourceArray = *(sourceArrayPtr.get());
    const size_t numComps = sourceArray.getNumberOfComponents();
    if(numComps == 0)
    {
      m_Filter->setErrorCondition(-1000, "Invalid DataArray with ZERO components");
      m_Filter->sendThreadSafeProgressMessage(QString("%1: Number of Components was Zero for array. Exiting Transform.").arg(sourceArray.getName()));
      return;
    }

    m_Filter->sendThreadSafeProgressMessage(QString("%1: Transform Starting").arg(sourceArray.getName()));

    DataArrayPointerType targetArrayPtr = std::dynamic_pointer_cast<DataArrayType>(m_TargetArray);

    //    std::ofstream originalPointsFile("/tmp/original_point_centers.csv", std::ios_base::binary);
    //    originalPointsFile << "x,y,z" << std::endl;
    //    std::ofstream transformedPointsFile("/tmp/transformed_point_centers.csv", std::ios_base::binary);
    //    transformedPointsFile << "x,y,z,index,error" << std::endl;

    Vector3s origImageGeomDims(m_Params.origImageGeom->getDimensions().data());
    Eigen::Vector3f coordsOld(0.0F, 0.0F, 0.0F);

    // TrilinearInterpolationData<T> interpolationValues;
    Vector3s oldGeomIndices = {std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max()};

    Eigen::Vector4f coordsNew;
    std::vector<T> pValues(8 * numComps);
    Eigen::Vector3f uvw;

    Matrix4fR inverseTransform = m_TransformationMatrix.inverse();

    for(int64_t k = 0; k < m_Params.zpNew; k++)
    {
      if(m_Filter->getCancel() || m_Filter->getErrorCode() < 0)
      {
        break;
      }
      int64_t ktot = (m_Params.xpNew * m_Params.ypNew) * k;

      m_Filter->sendThreadSafeProgressMessage(QString("%1: Interpolating values for slice '%2/%3'").arg(m_SourceArray->getName()).arg(k).arg(m_Params.zpNew));
      for(int64_t j = 0; j < m_Params.ypNew; j++)
      {
        int64_t jtot = (m_Params.xpNew) * j;
        for(int64_t i = 0; i < m_Params.xpNew; i++)
        {
          int64_t newIndex = ktot + jtot + i;

          coordsNew[0] = (static_cast<float>(i) * m_Params.xResNew) + m_Params.xMinNew + 0.5F * m_Params.xResNew;
          coordsNew[1] = (static_cast<float>(j) * m_Params.yResNew) + m_Params.yMinNew + 0.5F * m_Params.yResNew;
          coordsNew[2] = (static_cast<float>(k) * m_Params.zResNew) + m_Params.zMinNew + 0.5F * m_Params.zResNew;
          coordsNew[3] = 1.0F; // We take translation into account

          Eigen::Array4f coordsOld = inverseTransform * coordsNew;

          auto errorResult = m_Params.origImageGeom->computeCellIndex(coordsOld.data(), oldGeomIndices.data());
          //          transformedPointsFile << coordsNew[0] << "," << coordsNew[1] << "," << coordsNew[2] << "," << newIndex << "," << static_cast<int32_t>(errorResult) << std::endl;
          //          originalPointsFile << coordsOld[0] << "," << coordsOld[1] << "," << coordsOld[2]
          //          << ","
          //                             << oldIndex
          /* << "," << value */
          //                             << std::endl;
          // Now we know what voxel the new cell center maps back to in the original geometry.
          if(errorResult == ImageGeom::ErrorType::NoError)
          {
            size_t oldIndex = (origImageGeomDims[0] * origImageGeomDims[1] * oldGeomIndices[2]) + (origImageGeomDims[0] * oldGeomIndices[1]) + oldGeomIndices[0];
            int octant = FindOctant(m_Params, oldIndex, {coordsOld.data()});

            FindInterpolationValues(m_Params, oldIndex, octant, oldGeomIndices, coordsOld, sourceArray, pValues, uvw);
            for(size_t compIndex = 0; compIndex < numComps; compIndex++)
            {
              T value = CalculateInterpolatedValue(pValues, uvw, numComps, compIndex);
              targetArrayPtr->setComponent(newIndex, compIndex, value);
            }
          }
          else
          {
            T value = static_cast<T>(0);
            targetArrayPtr->fillTuple(newIndex, value);
          }
        }
      }
    }

    m_SourceArray->resizeTuples(0);

    m_Filter->sendThreadSafeProgressMessage(QString("%1: Transform Ending").arg(sourceArray.getName()));
  }
};

/**
 * @brief ExecuteParallelFunction
 * @param sourceArray
 * @param runner
 * @param args
 */
template <template <class, class> class ClassT, class ParallelRunnerT, class FilterType, class... ArgsT>
auto ExecuteParallelFunction(IDataArray::Pointer& sourceArray, ParallelRunnerT&& runner, FilterType* filter, ArgsT&&... args)
{

  if(DataArray<int8_t>::Pointer array = std::dynamic_pointer_cast<DataArray<int8_t>>(sourceArray))
  {
    return runner.template execute<>(ClassT<int8_t, FilterType>(std::forward<ArgsT>(args)...));
  }
  if(DataArray<int16_t>::Pointer array = std::dynamic_pointer_cast<DataArray<int16_t>>(sourceArray))
  {
    return runner.template execute<>(ClassT<int16_t, FilterType>(std::forward<ArgsT>(args)...));
  }
  if(DataArray<int32_t>::Pointer array = std::dynamic_pointer_cast<DataArray<int32_t>>(sourceArray))
  {
    return runner.template execute<>(ClassT<int32_t, FilterType>(std::forward<ArgsT>(args)...));
  }
  if(DataArray<int64_t>::Pointer array = std::dynamic_pointer_cast<DataArray<int64_t>>(sourceArray))
  {
    return runner.template execute<>(ClassT<int64_t, FilterType>(std::forward<ArgsT>(args)...));
  }
  if(DataArray<uint8_t>::Pointer array = std::dynamic_pointer_cast<DataArray<uint8_t>>(sourceArray))
  {
    return runner.template execute<>(ClassT<uint8_t, FilterType>(std::forward<ArgsT>(args)...));
  }
  if(DataArray<uint16_t>::Pointer array = std::dynamic_pointer_cast<DataArray<uint16_t>>(sourceArray))
  {
    return runner.template execute<>(ClassT<uint16_t, FilterType>(std::forward<ArgsT>(args)...));
  }
  if(DataArray<uint32_t>::Pointer array = std::dynamic_pointer_cast<DataArray<uint32_t>>(sourceArray))
  {
    return runner.template execute<>(ClassT<uint32_t, FilterType>(std::forward<ArgsT>(args)...));
  }
  if(DataArray<uint64_t>::Pointer array = std::dynamic_pointer_cast<DataArray<uint64_t>>(sourceArray))
  {
    return runner.template execute<>(ClassT<uint64_t, FilterType>(std::forward<ArgsT>(args)...));
  }
  if(DataArray<float>::Pointer array = std::dynamic_pointer_cast<DataArray<float>>(sourceArray))
  {
    return runner.template execute<>(ClassT<float, FilterType>(std::forward<ArgsT>(args)...));
  }
  if(DataArray<double>::Pointer array = std::dynamic_pointer_cast<DataArray<double>>(sourceArray))
  {
    return runner.template execute<>(ClassT<double, FilterType>(std::forward<ArgsT>(args)...));
  }

  throw std::runtime_error("Can not interpolate data type. Bool, String, StatsData?");
}

} // namespace ImageRotationUtilities
