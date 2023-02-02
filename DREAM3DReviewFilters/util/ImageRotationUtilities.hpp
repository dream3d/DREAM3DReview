#pragma once

#include "SIMPLib/Common/SIMPLArray.hpp"
#include "SIMPLib/DataArrays/DataArray.hpp"
#include "SIMPLib/Filtering/AbstractFilter.h"
#include "SIMPLib/Geometry/ImageGeom.h"
#include "SIMPLib/Math/MatrixMath.h"

#include <Eigen/Dense>

#include <fstream>
#include <iostream>

using Matrix3fR = Eigen::Matrix<float, 3, 3, Eigen::RowMajor>;
using Matrix4fR = Eigen::Matrix<float, 4, 4, Eigen::RowMajor>;
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


template<typename T>
inline void Vec3Add(const IVec3<T>& a, const IVec3<T>& b, IVec3<T>& c)
{
  c[0] = a[0] + b[0];
  c[1] = a[1] + b[1];
  c[2] = a[2] + b[2];
}

template<typename T>
inline void Vec3Subtract(const IVec3<T>& a, const IVec3<T>& b, IVec3<T>& c)
{
  c[0] = a[0] - b[0];
  c[1] = a[1] - b[1];
  c[2] = a[2] - b[2];
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
  params.origImageGeom->getCoords(index, centerPoint.data());
  // GetPlaneCoords(params, index, centerPoint.data());

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
    Vec3Subtract<float>(unitSquareCoords[i], coord, temp);
    float distance = temp.norm();
    if(distance < minDistance)
    {
      minDistance = distance;
      minIndex = i;
    }
  }

  return minIndex;
}

using OctantOffsetArrayType = std::array<Int64Vec3Type, 8>;

static const OctantOffsetArrayType k_IndexOffset0 = {Int64Vec3Type{-1, -1, -1}, Int64Vec3Type{0, -1, -1}, Int64Vec3Type{0, 0, -1}, Int64Vec3Type{-1, 0, -1},
                                                     Int64Vec3Type{-1, -1, 0},  Int64Vec3Type{0, -1, 0},  Int64Vec3Type{0, 0, 0},  Int64Vec3Type{-1, 0, 0}};
static const OctantOffsetArrayType k_IndexOffset1 = {Int64Vec3Type{0, -1, -1}, Int64Vec3Type{1, -1, -1}, Int64Vec3Type{1, 0, -1}, Int64Vec3Type{0, 0, -1},
                                                     Int64Vec3Type{0, -1, 0},  Int64Vec3Type{1, -1, 0},  Int64Vec3Type{1, 0, 0},  Int64Vec3Type{0, 0, 0}};
static const OctantOffsetArrayType k_IndexOffset2 = {Int64Vec3Type{0, 0, -1}, Int64Vec3Type{1, 0, -1}, Int64Vec3Type{1, 1, -1}, Int64Vec3Type{0, 1, -1},
                                                     Int64Vec3Type{0, 0, 0},  Int64Vec3Type{1, 0, 0},  Int64Vec3Type{1, 1, 0},  Int64Vec3Type{0, 1, 0}};
static const OctantOffsetArrayType k_IndexOffset3 = {Int64Vec3Type{-1, 0, -1}, Int64Vec3Type{0, 0, -1}, Int64Vec3Type{0, 1, -1}, Int64Vec3Type{-1, 1, -1},
                                                     Int64Vec3Type{-1, 0, 0},  Int64Vec3Type{0, 0, 0},  Int64Vec3Type{0, 1, 0},  Int64Vec3Type{-1, 1, 0}};
static const OctantOffsetArrayType k_IndexOffset4 = {Int64Vec3Type{-1, -1, 0}, Int64Vec3Type{0, -1, 0}, Int64Vec3Type{0, 0, 0}, Int64Vec3Type{-1, 0, 0},
                                                     Int64Vec3Type{-1, -1, 1}, Int64Vec3Type{0, -1, 1}, Int64Vec3Type{0, 0, 1}, Int64Vec3Type{-1, 0, 1}};
static const OctantOffsetArrayType k_IndexOffset5 = {Int64Vec3Type{0, -1, 0}, Int64Vec3Type{1, -1, 0}, Int64Vec3Type{1, 0, 0}, Int64Vec3Type{0, 0, 0},
                                                     Int64Vec3Type{0, -1, 1}, Int64Vec3Type{1, -1, 1}, Int64Vec3Type{1, 0, 1}, Int64Vec3Type{0, 0, 1}};
static const OctantOffsetArrayType k_IndexOffset6 = {Int64Vec3Type{0, 0, 0}, Int64Vec3Type{1, 0, 0}, Int64Vec3Type{1, 1, 0}, Int64Vec3Type{0, 1, 0},
                                                     Int64Vec3Type{0, 0, 1}, Int64Vec3Type{1, 0, 1}, Int64Vec3Type{1, 1, 1}, Int64Vec3Type{0, 1, 1}};
static const OctantOffsetArrayType k_IndexOffset7 = {Int64Vec3Type{-1, 0, 0}, Int64Vec3Type{0, 0, 0}, Int64Vec3Type{0, 1, 0}, Int64Vec3Type{-1, -1, 0},
                                                     Int64Vec3Type{-1, 0, 1}, Int64Vec3Type{0, 0, 1}, Int64Vec3Type{0, 1, 1}, Int64Vec3Type{-1, -1, 1}};

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
inline void FindInterpolationValues(const RotateArgs& params, size_t index, size_t octant, std::array<size_t, 3> oldIndicesU, FloatVec3Type& oldCoords, DataArray<T>& sourceArray,
                                    std::vector<T>& pValues, FloatVec3Type& uvw)
{
  const std::array<Int64Vec3Type, 8>& indexOffset = k_AllOctantOffsets[octant];

  Int64Vec3Type oldIndices(static_cast<int64_t>(oldIndicesU[0]), static_cast<int64_t>(oldIndicesU[1]), static_cast<int64_t>(oldIndicesU[2]));
  size_t numComps = sourceArray.getNumberOfComponents();

  Int64Vec3Type pIndices;
  FloatVec3Type p1Coord;

  auto origin = params.origImageGeom->getOrigin();

  for(size_t i = 0; i < 8; i++)
  {
    Vec3Add<int64_t>(oldIndices, indexOffset[i], pIndices);
    for(size_t compIndex = 0; compIndex < numComps; compIndex++)
    {
      T value = GetValue<T>(params, pIndices, sourceArray, index, compIndex);
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
template <typename T>
class RotateImageGeometryWithTrilinearInterpolation
{
private:
  AbstractFilter* m_Filter = nullptr;
  IDataArray::Pointer m_SourceArray;
  IDataArray::Pointer m_TargetArray;
  float m_RotMatrixInv[3][3] = {{0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f}};
  ImageRotationUtilities::RotateArgs m_Params;

public:
  RotateImageGeometryWithTrilinearInterpolation(AbstractFilter* filter, IDataArray::Pointer& sourceArray, IDataArray::Pointer& targetArray, RotateArgs& rotateArgs, Matrix4fR& transformationMatrix)
  : m_Filter(filter)
  , m_SourceArray(sourceArray)
  , m_TargetArray(targetArray)
  , m_Params(rotateArgs)
  {
    //    // We have to inline the 3x3 Maxtrix transpose here because of the "const" nature of the 'convert' function
    //    Matrix3fR transpose = rotationMatrix.transpose();
    //    // Need to use row based Eigen matrix so that the values get mapped to the right place in the raw array
    //    // Raw array is faster than Eigen
    //    Eigen::Map<Matrix3fR>(&m_RotMatrixInv[0][0], transpose.rows(), transpose.cols()) = transpose;
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
  T CalculateInterpolatedValue(std::vector<T>& pValues, FloatVec3Type& uvw, size_t numComps, size_t compIndex) const
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

    //  if(INDEX_VALID(indices[0]) && INDEX_VALID(indices[1]))
    {
      value += u * (pValues[P2 * numComps + compIndex] - pValues[P1 * numComps + compIndex]);
    }
    //   if(INDEX_VALID(indices[3]) && INDEX_VALID(indices[0]))
    {
      value += v * (pValues[P4 * numComps + compIndex] - pValues[P1 * numComps + compIndex]);
    }
    //  if(INDEX_VALID(indices[4]) && INDEX_VALID(indices[0]))
    {
      value += w * (pValues[P5 * numComps + compIndex] - pValues[P1 * numComps + compIndex]);
    }
    //  if(INDEX_VALID(indices[0]) && INDEX_VALID(indices[1]) && INDEX_VALID(indices[2]) && INDEX_VALID(indices[3]))
    {
      value += u * v * (pValues[P1 * numComps + compIndex] - pValues[P2 * numComps + compIndex] + pValues[P3 * numComps + compIndex] - pValues[P4 * numComps + compIndex]);
    }
    //   if(INDEX_VALID(indices[0]) && INDEX_VALID(indices[1]) && INDEX_VALID(indices[4]) && INDEX_VALID(indices[5]))
    {
      value += u * w * (pValues[P1 * numComps + compIndex] - pValues[P2 * numComps + compIndex] - pValues[P5 * numComps + compIndex] + pValues[P6 * numComps + compIndex]);
    }
    //   if(INDEX_VALID(indices[0]) && INDEX_VALID(indices[3]) && INDEX_VALID(indices[4]) && INDEX_VALID(indices[7]))
    {
      value += v * w * (pValues[P1 * numComps + compIndex] - pValues[P4 * numComps + compIndex] - pValues[P5 * numComps + compIndex] + pValues[P8 * numComps + compIndex]);
    }

    //  if(INDEX_VALID(indices[0]) && INDEX_VALID(indices[1]) && INDEX_VALID(indices[2]) && INDEX_VALID(indices[3]) && INDEX_VALID(indices[4]) && INDEX_VALID(indices[5]) && INDEX_VALID(indices[6]) &&
    //     INDEX_VALID(indices[7]))
    {
      value += u * v * w *
               (-pValues[P1 * numComps + compIndex] + pValues[P2 * numComps + compIndex] - pValues[P3 * numComps + compIndex] + pValues[P4 * numComps + compIndex] -
                pValues[P5 * numComps + compIndex] + pValues[P6 * numComps + compIndex] - pValues[P7 * numComps + compIndex] + pValues[P8 * numComps + compIndex]);
    }
    // clang-format on
    return value;
  }

/**
 * 
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
      exit(1);
    }
    DataArrayPointerType targetArrayPtr = std::dynamic_pointer_cast<DataArrayType>(m_TargetArray);

    std::ofstream originalPointsFile("/tmp/original_point_centers.csv", std::ios_base::binary);
    originalPointsFile << "x,y,z,i,j,k,index,octant" << std::endl;
    std::ofstream transformedPointsFile("/tmp/transformed_point_centers.csv", std::ios_base::binary);
    transformedPointsFile << "x,y,z,index,error" << std::endl;

    FloatVec3Type coordsOld = {0.0f, 0.0f, 0.0f};
    SizeVec3Type origImageGeomDims = m_Params.origImageGeom->getDimensions();
    // TrilinearInterpolationData<T> interpolationValues;
    std::array<size_t, 3> oldGeomIndices = {std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max()};
    std::array<float, 3> coordsNew;
    std::vector<T> pValues(8 * numComps);
    FloatVec3Type uvw;
    for(int64_t k = 0; k < m_Params.zpNew; k++)
    {
      int64_t ktot = (m_Params.xpNew * m_Params.ypNew) * k;


      m_Filter->notifyStatusMessage(QString("Transforming Slice '%1'").arg(k));
      for(int64_t j = 0; j < m_Params.ypNew; j++)
      {
        int64_t jtot = (m_Params.xpNew) * j;
        for(int64_t i = 0; i < m_Params.xpNew; i++)
        {

          int64_t newIndex = ktot + jtot + i;
          size_t oldIndex = std::numeric_limits<size_t>::max();

          coordsNew[0] = (static_cast<float>(i) * m_Params.xResNew) + m_Params.xMinNew + 0.5F * m_Params.xResNew;
          coordsNew[1] = (static_cast<float>(j) * m_Params.yResNew) + m_Params.yMinNew + 0.5F * m_Params.yResNew;
          coordsNew[2] = (static_cast<float>(k) * m_Params.zResNew) + m_Params.zMinNew + 0.5F * m_Params.zResNew;

          MatrixMath::Multiply3x3with3x1(m_RotMatrixInv, coordsNew.data(), coordsOld.data());

          auto errorResult = m_Params.origImageGeom->computeCellIndex(coordsOld.data(), oldGeomIndices.data());
          transformedPointsFile << coordsNew[0] << "," << coordsNew[1] << "," << coordsNew[2] << "," << newIndex << "," << static_cast<int32_t>(errorResult) << std::endl;

          // Now we know what voxel the new cell center maps back to in the original geometry.
          if(errorResult == ImageGeom::ErrorType::NoError)
          {
            oldIndex = (origImageGeomDims[0] * origImageGeomDims[1] * oldGeomIndices[2]) + (origImageGeomDims[0] * oldGeomIndices[1]) + oldGeomIndices[0];
            int octant = FindOctant(m_Params, oldIndex, {coordsOld.data()});

            FindInterpolationValues(m_Params, oldIndex, octant, oldGeomIndices, coordsOld, sourceArray, pValues, uvw);
            for(size_t compIndex = 0; compIndex < numComps; compIndex++)
            {
              T value = CalculateInterpolatedValue(pValues, uvw, numComps, compIndex);
              targetArrayPtr->setComponent(newIndex, compIndex, value);
              originalPointsFile << coordsOld[0] << "," << coordsOld[1] << "," << coordsOld[2] << "," << oldGeomIndices[0] << "," << oldGeomIndices[1] << "," << oldGeomIndices[2] << "," << oldIndex
                                 << "," << value << std::endl;
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
