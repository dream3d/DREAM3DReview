#pragma once

#include "SIMPLib/Common/SIMPLArray.hpp"
#include "SIMPLib/DataArrays/DataArray.hpp"
#include "SIMPLib/Filtering/AbstractFilter.h"
#include "SIMPLib/Geometry/ImageGeom.h"
#include "SIMPLib/Math/MatrixMath.h"

#include <Eigen/Dense>

#include <iostream>

using Matrix3fR = Eigen::Matrix<float, 3, 3, Eigen::RowMajor>;
using Transform3f = Eigen::Transform<float, 3, Eigen::Affine>;
using MatrixTranslation = Eigen::Matrix<float, 1, 3, Eigen::RowMajor>;

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
struct TrilinearInterpolationData
{
  std::vector<T> pValues;
  FloatVec3Type uvw;
};

template <typename T>
T inline GetValue(const RotateArgs& params, Int64Vec3Type xyzIndex, DataArray<T>& sourceArray, size_t exemplarIndex, size_t compIndex)
{
  // if any of the X, Y, or Z are out of bounds, then just return the target voxel data value.
  if(xyzIndex[0] < 0 || xyzIndex[0] >= params.xp || xyzIndex[1] < 0 || xyzIndex[1] >= params.yp || xyzIndex[2] < 0 || xyzIndex[2] >= params.zp)
  {
    return sourceArray.getComponent(exemplarIndex, compIndex);
  }

  // Now just compute the proper index
  size_t index = (xyzIndex[2] * params.xp * params.yp) + (xyzIndex[1] * params.xp) + xyzIndex[0];
  return sourceArray.getComponent(index, compIndex);
}

void inline GetPlaneCoords(const RotateArgs& params, size_t idx, float coords[3])
{
  size_t column = idx % params.xp;
  size_t row = (idx / params.xp) % params.yp;
  size_t plane = idx / (params.xp * params.yp);

  coords[0] = column * params.xRes;
  coords[1] = row * params.yRes;
  coords[2] = plane * params.zRes;
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
  FloatVec3Type centerPoint;
  GetPlaneCoords(params, index, centerPoint.data());

  // Form the 8 corner coords for the voxel
  // clang-format off
  std::array<FloatVec3Type, 8> unitSquareCoords = {
  /* P1 */ centerPoint + FloatVec3Type(-xResHalf, -yResHalf, -zResHalf),
  /* P2 */ centerPoint + FloatVec3Type(xResHalf, -yResHalf, -zResHalf),
  /* P3 */ centerPoint + FloatVec3Type(xResHalf, yResHalf, -zResHalf),
  /* P4 */ centerPoint + FloatVec3Type(-xResHalf, yResHalf, -zResHalf),
  /* P5 */ centerPoint + FloatVec3Type(-xResHalf, -yResHalf, zResHalf),
  /* P6 */ centerPoint + FloatVec3Type(xResHalf, -yResHalf, zResHalf),
  /* P7 */ centerPoint + FloatVec3Type(xResHalf, yResHalf, zResHalf),
  /* P8 */ centerPoint + FloatVec3Type(-xResHalf, yResHalf, zResHalf),
  };
  // clang-format on

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

template <typename T>
inline TrilinearInterpolationData<T> FindNeighborIndices(const RotateArgs& params, size_t index, size_t octant, std::array<size_t, 3> oldIndicesU, FloatVec3Type& oldCoords, DataArray<T>& sourceArray)
{
  //  int64_t idx = static_cast<int64_t>(index);
  //  int64_t pointsInPlane = static_cast<int64_t>(params.xp * params.yp);
  std::array<int64_t, 3> oldIndices = {static_cast<int64_t>(oldIndicesU[0]), static_cast<int64_t>(oldIndicesU[1]), static_cast<int64_t>(oldIndicesU[2])};

  std::array<float, 3> res = {params.xRes, params.yRes, params.zRes};

  std::array<Int64Vec3Type, 8> pIndices;

  FloatVec3Type p1Coord;

  FloatVec3Type uvw;
  size_t numComps = sourceArray.getNumberOfComponents();
  std::vector<T> pValues(8 * numComps);

  switch(octant)
  {
  case 0: {
    pIndices[0] = {(oldIndices[0] - 1), (oldIndices[1] - 1), (oldIndices[2] - 1)};
    pIndices[1] = {(oldIndices[0]), (oldIndices[1] - 1), (oldIndices[2] - 1)};
    pIndices[2] = {(oldIndices[0]), (oldIndices[1]), (oldIndices[2] - 1)};
    pIndices[3] = {(oldIndices[0] - 1), (oldIndices[1]), (oldIndices[2] - 1)};
    pIndices[4] = {(oldIndices[0] - 1), (oldIndices[1] - 1), (oldIndices[2])};
    pIndices[5] = {(oldIndices[0]), (oldIndices[1] - 1), (oldIndices[2])};
    pIndices[6] = {(oldIndices[0]), (oldIndices[1]), (oldIndices[2])};
    pIndices[7] = {(oldIndices[0] - 1), (oldIndices[1]), (oldIndices[2])};

    p1Coord = {pIndices[0][0] * res[0], pIndices[0][1] * res[1], pIndices[0][2] * res[2]};
    uvw = oldCoords - p1Coord;

    for(size_t i = 0; i < 8; i++)
    {
      for(size_t compIndex = 0; compIndex < numComps; compIndex++)
      {
        pValues[i * numComps + compIndex] = GetValue<T>(params, pIndices[i], sourceArray, index, compIndex);
      }
    }

    return {std::move(pValues), std::move(uvw)};
  }
  case 1:
  //  return {idx - (pointsInPlane - params.xp), idx - (pointsInPlane - params.xp + 1), idx - (pointsInPlane + 1), idx - (pointsInPlane), idx - (params.xp), idx - (params.xp + 1), idx + (1), idx};
  case 2:
  //  return {idx - (pointsInPlane), idx - (pointsInPlane + 1), idx - (pointsInPlane + params.xp + 1), idx - (pointsInPlane + params.xp), idx, idx + 1, idx - (params.xp + 1), idx - (params.xp)};
  case 3: {

    pIndices[0] = {(oldIndices[0] - 1), (oldIndices[1]), (oldIndices[2] - 1)};
    pIndices[1] = {(oldIndices[0]), (oldIndices[1]), (oldIndices[2] - 1)};
    pIndices[2] = {(oldIndices[0]), (oldIndices[1] + 1), (oldIndices[2] - 1)};
    pIndices[3] = {(oldIndices[0] + 1), (oldIndices[1] + 1), (oldIndices[2] - 1)};
    pIndices[4] = {(oldIndices[0] - 1), (oldIndices[1]), (oldIndices[2] + 1)};
    pIndices[5] = {(oldIndices[0]), (oldIndices[1]), (oldIndices[2] + 1)};
    pIndices[6] = {(oldIndices[0]), (oldIndices[1] + 1), (oldIndices[2] + 1)};
    pIndices[7] = {(oldIndices[0] + 1), (oldIndices[1] + 1), (oldIndices[2] + 1)};

    p1Coord = {pIndices[0][0] * res[0], pIndices[0][1] * res[1], pIndices[0][2] * res[2]};
    uvw = oldCoords - p1Coord;

    for(size_t i = 0; i < 8; i++)
    {
      for(size_t compIndex = 0; compIndex < numComps; compIndex++)
      {
        pValues[i * numComps + compIndex] = GetValue<T>(params, pIndices[i], sourceArray, index, compIndex);
      }
    }

    return {std::move(pValues), std::move(uvw)};
  }
  case 4:
    // return {idx - (params.xp - 1), idx - (params.xp), idx, idx - 1, idx + (pointsInPlane - params.xp - 1), idx + (pointsInPlane - params.xp), idx + (pointsInPlane), idx + (pointsInPlane - 1)};
  case 5:
  //  return {idx - (params.xp), idx - (params.xp + 1), idx + (1), idx, idx + (pointsInPlane - params.xp), idx + (pointsInPlane - params.xp + 1), idx + (pointsInPlane + 1), idx + (pointsInPlane)};
  case 6:
  //  return {idx, idx + 1, idx + (params.xp + 1), idx + (params.xp), idx + (pointsInPlane), idx + (pointsInPlane + 1), idx + (pointsInPlane + params.xp + 1), idx + (pointsInPlane + params.xp)};
  case 7:
    //  return {idx - 1, idx, idx + params.xp, idx + (params.xp - 1), idx + (pointsInPlane - 1), idx + (pointsInPlane), idx + (pointsInPlane + params.xp), idx + (params.xp - 1)};
    break;
  }
  return {};
}

#define INDEX_VALID(P) (P >= 0 && P < numTuples)

template <typename T>
class RotateWithInterpolation
{
private:
  AbstractFilter* m_Filter = nullptr;
  IDataArray::Pointer m_SourceArray;
  IDataArray::Pointer m_TargetArray;
  float m_RotMatrixInv[3][3] = {{0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f}};
  ImageRotationUtilities::RotateArgs m_Params;

public:
  RotateWithInterpolation(AbstractFilter* filter, IDataArray::Pointer& sourceArray, IDataArray::Pointer& targetArray, RotateArgs& rotateArgs, Matrix3fR& rotationMatrix)
  : m_Filter(filter)
  , m_SourceArray(sourceArray)
  , m_TargetArray(targetArray)
  , m_Params(rotateArgs)
  {
    // We have to inline the 3x3 Maxtrix transpose here because of the "const" nature of the 'convert' function
    Matrix3fR transpose = rotationMatrix.transpose();
    // Need to use row based Eigen matrix so that the values get mapped to the right place in the raw array
    // Raw array is faster than Eigen
    Eigen::Map<Matrix3fR>(&m_RotMatrixInv[0][0], transpose.rows(), transpose.cols()) = transpose;
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
  T CalculateInterpolatedValue(const RotateArgs& params, DataArray<T>& sourceArray, FloatVec3Type sourceCoord, int64_t oldIndex, std::array<int64_t, 8>& indices, size_t compIndex) const
  {
    const size_t numTuples = params.xp * params.yp * params.zp;
    const size_t numComps = sourceArray.getNumberOfComponents();
    const size_t P1 = static_cast<size_t>(indices[0]);
    const size_t P2 = static_cast<size_t>(indices[1]);
    const size_t P3 = static_cast<size_t>(indices[2]);
    const size_t P4 = static_cast<size_t>(indices[3]);
    const size_t P5 = static_cast<size_t>(indices[4]);
    const size_t P6 = static_cast<size_t>(indices[5]);
    const size_t P7 = static_cast<size_t>(indices[6]);
    const size_t P8 = static_cast<size_t>(indices[7]);

    FloatVec3Type p1Coords;
    GetPlaneCoords(params, P1, p1Coords.data());
    auto uvw = sourceCoord - p1Coords;

    const float u = uvw[0];
    const float v = uvw[1];
    const float w = uvw[2];

    T value = sourceArray[oldIndex];

    if(INDEX_VALID(indices[0]) && INDEX_VALID(indices[1]))
    {
      value += u * (sourceArray[P2 * numComps + compIndex] - sourceArray[P1 * numComps + compIndex]);
    }
    if(INDEX_VALID(indices[3]) && INDEX_VALID(indices[0]))
    {
      value += v * (sourceArray[P4 * numComps + compIndex] - sourceArray[P1 * numComps + compIndex]);
    }
    if(INDEX_VALID(indices[4]) && INDEX_VALID(indices[0]))
    {
      value += w * (sourceArray[P5 * numComps + compIndex] - sourceArray[P1 * numComps + compIndex]);
    }
    if(INDEX_VALID(indices[0]) && INDEX_VALID(indices[1]) && INDEX_VALID(indices[2]) && INDEX_VALID(indices[3]))
    {
      value += u * v * (sourceArray[P1 * numComps + compIndex] - sourceArray[P2 * numComps + compIndex] + sourceArray[P3 * numComps + compIndex] - sourceArray[P4 * numComps + compIndex]);
    }
    if(INDEX_VALID(indices[0]) && INDEX_VALID(indices[1]) && INDEX_VALID(indices[4]) && INDEX_VALID(indices[5]))
    {
      value += u * w * (sourceArray[P1 * numComps + compIndex] - sourceArray[P2 * numComps + compIndex] - sourceArray[P5 * numComps + compIndex] + sourceArray[P6 * numComps + compIndex]);
    }
    if(INDEX_VALID(indices[0]) && INDEX_VALID(indices[3]) && INDEX_VALID(indices[4]) && INDEX_VALID(indices[7]))
    {
      value += v * w * (sourceArray[P1 * numComps + compIndex] - sourceArray[P4 * numComps + compIndex] - sourceArray[P5 * numComps + compIndex] + sourceArray[P8 * numComps + compIndex]);
    }

    if(INDEX_VALID(indices[0]) && INDEX_VALID(indices[1]) && INDEX_VALID(indices[2]) && INDEX_VALID(indices[3]) && INDEX_VALID(indices[4]) && INDEX_VALID(indices[5]) && INDEX_VALID(indices[6]) &&
       INDEX_VALID(indices[7]))
    {
      value += u * v * w *
               (-sourceArray[P1 * numComps + compIndex] + sourceArray[P2 * numComps + compIndex] - sourceArray[P3 * numComps + compIndex] + sourceArray[P4 * numComps + compIndex] -
                sourceArray[P5 * numComps + compIndex] + sourceArray[P6 * numComps + compIndex] - sourceArray[P7 * numComps + compIndex] + sourceArray[P8 * numComps + compIndex]);
    }
    // clang-format on

    return value;
  }

  void operator()() const
  {
    using DataArrayType = DataArray<T>;
    using DataArrayPointerType = typename DataArrayType::Pointer;

    DataArrayPointerType sourceArrayPtr = std::dynamic_pointer_cast<DataArrayType>(m_SourceArray);
    DataArrayType& sourceArray = *(sourceArrayPtr.get());
    size_t numComps = sourceArray.getNumberOfComponents();
    DataArrayPointerType targetArrayPtr = std::dynamic_pointer_cast<DataArrayType>(m_TargetArray);
    // DataArrayType& targetArray = *(targetArrayPtr.get());

    FloatVec3Type coordsOld = {0.0f, 0.0f, 0.0f};

    for(int64_t k = 0; k < m_Params.zpNew; k++)
    {
      int64_t ktot = (m_Params.xpNew * m_Params.ypNew) * k;
      for(int64_t j = 0; j < m_Params.ypNew; j++)
      {
        int64_t jtot = (m_Params.xpNew) * j;
        for(int64_t i = 0; i < m_Params.xpNew; i++)
        {

          int64_t newIndex = ktot + jtot + i;
          size_t oldIndex = std::numeric_limits<size_t>::max();

          std::array<float, 3> coordsNew = {(static_cast<float>(i) * m_Params.xResNew) + m_Params.xMinNew + 0.5F * m_Params.xResNew,
                                            (static_cast<float>(j) * m_Params.yResNew) + m_Params.yMinNew + 0.5F * m_Params.yResNew,
                                            (static_cast<float>(k) * m_Params.zResNew) + m_Params.zMinNew + 0.5F * m_Params.zResNew};

          MatrixMath::Multiply3x3with3x1(m_RotMatrixInv, coordsNew.data(), coordsOld.data());

          std::array<size_t, 3> oldGeomIndices = {std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max()};
          auto errorResult = m_Params.origImageGeom->computeCellIndex(coordsOld.data(), oldGeomIndices.data());
          errorResult = m_Params.origImageGeom->computeCellIndex(coordsOld.data(), oldIndex);

          //          auto colOld = static_cast<int64_t>(std::nearbyint(coordsOld[0] / m_Params.xRes));
          //          auto rowOld = static_cast<int64_t>(std::nearbyint(coordsOld[1] / m_Params.yRes));
          //          auto planeOld = static_cast<int64_t>(std::nearbyint(coordsOld[2] / m_Params.zRes));

          std::cout << coordsOld[0] << ", " << coordsOld[1] << ", " << coordsOld[2] << "," << coordsNew[0] << ", " << coordsNew[1] << ", " << coordsNew[2] << "," << oldGeomIndices[0] << ", "
                    << oldGeomIndices[1] << ", " << oldGeomIndices[2] << "," << newIndex << "," << oldIndex << ", " << static_cast<int32_t>(errorResult) << std::endl;

          // Now we know what voxel the new cell center maps back to in the original geometry.
          if(errorResult == ImageGeom::ErrorType::NoError)
          {

            int octant = FindOctant(m_Params, oldIndex, {coordsOld.data()});
            auto neighborIndices = FindNeighborIndices(m_Params, oldIndex, octant, oldGeomIndices, coordsOld, sourceArray);
            for(size_t compIndex = 0; compIndex < numComps; compIndex++)
            {
              // T value = CalculateInterpolatedValue(m_Params, sourceArray, coordsOld, oldIndex, neighborIndices, compIndex);
              // targetArrayPtr->setComponent(newIndex, compIndex, value);
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
  }
};

/**
 * @brief ExecuteParallelFunction
 * @param sourceArray
 * @param runner
 * @param args
 */
template <template <class> class ClassT, class ParallelRunnerT, class... ArgsT>
auto ExecuteParallelFunction(IDataArray::Pointer& sourceArray, ParallelRunnerT&& runner, ArgsT&&... args)
{

  if(DataArray<int8_t>::Pointer array = std::dynamic_pointer_cast<DataArray<int8_t>>(sourceArray))
  {
    return runner.template execute<>(ClassT<int8_t>(std::forward<ArgsT>(args)...));
  }
  if(DataArray<int16_t>::Pointer array = std::dynamic_pointer_cast<DataArray<int16_t>>(sourceArray))
  {
    return runner.template execute<>(ClassT<int16_t>(std::forward<ArgsT>(args)...));
  }
  if(DataArray<int32_t>::Pointer array = std::dynamic_pointer_cast<DataArray<int32_t>>(sourceArray))
  {
    return runner.template execute<>(ClassT<int32_t>(std::forward<ArgsT>(args)...));
  }
  if(DataArray<int64_t>::Pointer array = std::dynamic_pointer_cast<DataArray<int64_t>>(sourceArray))
  {
    return runner.template execute<>(ClassT<int64_t>(std::forward<ArgsT>(args)...));
  }
  if(DataArray<uint8_t>::Pointer array = std::dynamic_pointer_cast<DataArray<uint8_t>>(sourceArray))
  {
    return runner.template execute<>(ClassT<uint8_t>(std::forward<ArgsT>(args)...));
  }
  if(DataArray<uint16_t>::Pointer array = std::dynamic_pointer_cast<DataArray<uint16_t>>(sourceArray))
  {
    return runner.template execute<>(ClassT<uint16_t>(std::forward<ArgsT>(args)...));
  }
  if(DataArray<uint32_t>::Pointer array = std::dynamic_pointer_cast<DataArray<uint32_t>>(sourceArray))
  {
    return runner.template execute<>(ClassT<uint32_t>(std::forward<ArgsT>(args)...));
  }
  if(DataArray<uint64_t>::Pointer array = std::dynamic_pointer_cast<DataArray<uint64_t>>(sourceArray))
  {
    return runner.template execute<>(ClassT<uint64_t>(std::forward<ArgsT>(args)...));
  }
  if(DataArray<float>::Pointer array = std::dynamic_pointer_cast<DataArray<float>>(sourceArray))
  {
    return runner.template execute<>(ClassT<float>(std::forward<ArgsT>(args)...));
  }
  if(DataArray<double>::Pointer array = std::dynamic_pointer_cast<DataArray<double>>(sourceArray))
  {
    return runner.template execute<>(ClassT<double>(std::forward<ArgsT>(args)...));
  }

  throw std::runtime_error("Can not interpolate data type. Bool, String, StatsData?");
}

} // namespace ImageRotationUtilities
