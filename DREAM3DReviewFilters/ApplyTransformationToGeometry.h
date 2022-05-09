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
 * The code contained herein was partially funded by the following contracts:
 *    United States Air Force Prime Contract FA8650-07-D-5800
 *    United States Air Force Prime Contract FA8650-10-D-5210
 *    United States Prime Contract Navy N00173-07-C-2068
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#pragma once

#include <memory>
#include <mutex>

#include <Eigen/Dense>

#include "DREAM3DReview/DREAM3DReviewDLLExport.h"

#include "SIMPLib/SIMPLib.h"
#include "SIMPLib/DataArrays/DataArray.hpp"
#include "SIMPLib/FilterParameters/DynamicTableData.h"
#include "SIMPLib/FilterParameters/FloatVec3FilterParameter.h"
#include "SIMPLib/Filtering/AbstractFilter.h"

/**
 * @brief The ApplyTransformationToGeometry class. See [Filter documentation](@ref applytransformationtogeometry) for details.
 */
class DREAM3DReview_EXPORT ApplyTransformationToGeometry : public AbstractFilter
{
  Q_OBJECT

  // Start Python bindings declarations
  PYB11_BEGIN_BINDINGS(ApplyTransformationToGeometry SUPERCLASS AbstractFilter)
  PYB11_FILTER()
  PYB11_SHARED_POINTERS(ApplyTransformationToGeometry)
  PYB11_FILTER_NEW_MACRO(ApplyTransformationToGeometry)
  PYB11_PROPERTY(DynamicTableData ManualTransformationMatrix READ getManualTransformationMatrix WRITE setManualTransformationMatrix)
  PYB11_PROPERTY(DataArrayPath ComputedTransformationMatrix READ getComputedTransformationMatrix WRITE setComputedTransformationMatrix)
  PYB11_PROPERTY(DataArrayPath CellAttributeMatrixPath READ getCellAttributeMatrixPath WRITE setCellAttributeMatrixPath)
  //  PYB11_PROPERTY(DataArrayPath GeometryToTransform READ getGeometryToTransform WRITE setGeometryToTransform)
  PYB11_PROPERTY(int TransformationMatrixType READ getTransformationMatrixType WRITE setTransformationMatrixType)
  PYB11_PROPERTY(int InterpolationType READ getInterpolationType WRITE setInterpolationType)
  PYB11_PROPERTY(bool UseDataArraySelection READ getUseDataArraySelection WRITE setUseDataArraySelection)
  PYB11_PROPERTY(std::vector<DataArrayPath> DataArraySelection READ getDataArraySelection WRITE setDataArraySelection)
  PYB11_PROPERTY(FloatVec3Type RotationAxis READ getRotationAxis WRITE setRotationAxis)
  PYB11_PROPERTY(float RotationAngle READ getRotationAngle WRITE setRotationAngle)
  PYB11_PROPERTY(FloatVec3Type Translation READ getTranslation WRITE setTranslation)
  PYB11_PROPERTY(FloatVec3Type Scale READ getScale WRITE setScale)
  PYB11_END_BINDINGS()
  // End Python bindings declarations

public:
  using Self = ApplyTransformationToGeometry;
  using Pointer = std::shared_ptr<Self>;
  using ConstPointer = std::shared_ptr<const Self>;
  using WeakPointer = std::weak_ptr<Self>;
  using ConstWeakPointer = std::weak_ptr<const Self>;
  static Pointer NullPointer();

  static std::shared_ptr<ApplyTransformationToGeometry> New();

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

  using Matrix3fR = Eigen::Matrix<float, 3, 3, Eigen::RowMajor>;
  using Transform3f = Eigen::Transform<float, 3, Eigen::Affine>;
  using MatrixTranslation = Eigen::Matrix<float, 1, 3, Eigen::RowMajor>;

  /**
   * @brief Returns the name of the class for ApplyTransformationToGeometry
   */
  QString getNameOfClass() const override;
  /**
   * @brief Returns the name of the class for ApplyTransformationToGeometry
   */
  static QString ClassName();

  ~ApplyTransformationToGeometry() override;

  /**
   * @brief Setter property for ManualTransformationMatrix
   */
  void setManualTransformationMatrix(const DynamicTableData& value);
  /**
   * @brief Getter property for ManualTransformationMatrix
   * @return Value of ManualTransformationMatrix
   */
  DynamicTableData getManualTransformationMatrix() const;
  Q_PROPERTY(DynamicTableData ManualTransformationMatrix READ getManualTransformationMatrix WRITE setManualTransformationMatrix)

  /**
   * @brief Setter property for ComputedTransformationMatrix
   */
  void setComputedTransformationMatrix(const DataArrayPath& value);
  /**
   * @brief Getter property for ComputedTransformationMatrix
   * @return Value of ComputedTransformationMatrix
   */
  DataArrayPath getComputedTransformationMatrix() const;
  Q_PROPERTY(DataArrayPath ComputedTransformationMatrix READ getComputedTransformationMatrix WRITE setComputedTransformationMatrix)

  //  /**
  //   * @brief Setter property for GeometryToTransform
  //   */
  //  void setGeometryToTransform(const DataArrayPath& value);
  //  /**
  //   * @brief Getter property for GeometryToTransform
  //   * @return Value of GeometryToTransform
  //   */
  //  DataArrayPath getGeometryToTransform() const;
  //  Q_PROPERTY(DataArrayPath GeometryToTransform READ getGeometryToTransform WRITE setGeometryToTransform)

  /**
   * @brief Setter property for TransformationMatrixType
   */
  void setTransformationMatrixType(int value);
  /**
   * @brief Getter property for TransformationMatrixType
   * @return Value of TransformationMatrixType
   */
  int getTransformationMatrixType() const;
  Q_PROPERTY(int TransformationMatrixType READ getTransformationMatrixType WRITE setTransformationMatrixType)

  /**
   * @brief Setter property for InterpolationType
   */
  void setInterpolationType(int value);
  /**
   * @brief Getter property for InterpolationType
   * @return Value of InterpolationType
   */
  int getInterpolationType() const;
  Q_PROPERTY(int InterpolationType READ getInterpolationType WRITE setInterpolationType)

  /**
   * @brief Setter property for RotationAxis
   */
  void setRotationAxis(const FloatVec3Type& value);
  /**
   * @brief Getter property for RotationAxis
   * @return Value of RotationAxis
   */
  FloatVec3Type getRotationAxis() const;
  Q_PROPERTY(FloatVec3Type RotationAxis READ getRotationAxis WRITE setRotationAxis)

  /**
   * @brief Setter property for RotationAngle
   */
  void setRotationAngle(float value);
  /**
   * @brief Getter property for RotationAngle
   * @return Value of RotationAngle
   */
  float getRotationAngle() const;
  Q_PROPERTY(float RotationAngle READ getRotationAngle WRITE setRotationAngle)

  /**
   * @brief Setter property for Translation
   */
  void setTranslation(const FloatVec3Type& value);
  /**
   * @brief Getter property for Translation
   * @return Value of Translation
   */
  FloatVec3Type getTranslation() const;
  Q_PROPERTY(FloatVec3Type Translation READ getTranslation WRITE setTranslation)

  /**
   * @brief Setter property for Scale
   */
  void setScale(const FloatVec3Type& value);
  /**
   * @brief Getter property for Scale
   * @return Value of Scale
   */
  FloatVec3Type getScale() const;
  Q_PROPERTY(FloatVec3Type Scale READ getScale WRITE setScale)

  /**
   * @brief Setter property for UseDataArraySelection
   */
  void setUseDataArraySelection(bool value);

  /**
   * @brief Getter property for UseDataArraySelection
   * @return Value of UseDataArraySelection
   */
  bool getUseDataArraySelection() const;
  Q_PROPERTY(bool UseDataArraySelection READ getUseDataArraySelection WRITE setUseDataArraySelection)

  /**
   * @brief Setter property for DataArraySelection
   */
  void setDataArraySelection(const std::vector<DataArrayPath>& value);

  /**
   * @brief Getter property for DataArraySelection
   * @return Value of DataArraySelection
   */
  std::vector<DataArrayPath> getDataArraySelection() const;
  Q_PROPERTY(std::vector<DataArrayPath> DataArraySelection READ getDataArraySelection WRITE setDataArraySelection)

  /**
   * @brief getCompiledLibraryName Reimplemented from @see AbstractFilter class
   */
  QString getCompiledLibraryName() const override;

  /**
   * @brief getBrandingString Returns the branding string for the filter, which is a tag
   * used to denote the filter's association with specific plugins
   * @return Branding string
   */
  QString getBrandingString() const override;

  /**
   * @brief getFilterVersion Returns a version string for this filter. Default
   * value is an empty string.
   * @return
   */
  QString getFilterVersion() const override;

  /**
   * @brief newFilterInstance Reimplemented from @see AbstractFilter class
   */
  AbstractFilter::Pointer newFilterInstance(bool copyFilterParameters) const override;

  /**
   * @brief getGroupName Reimplemented from @see AbstractFilter class
   */
  QString getGroupName() const override;

  /**
   * @brief getSubGroupName Reimplemented from @see AbstractFilter class
   */
  QString getSubGroupName() const override;

  /**
   * @brief getUuid Return the unique identifier for this filter.
   * @return A QUuid object.
   */
  QUuid getUuid() const override;

  /**
   * @brief getHumanLabel Reimplemented from @see AbstractFilter class
   */
  QString getHumanLabel() const override;

  /**
   * @brief setupFilterParameters Reimplemented from @see AbstractFilter class
   */
  void setupFilterParameters() override;

  /**
   * @brief readFilterParameters Reimplemented from @see AbstractFilter class
   */
  void readFilterParameters(AbstractFilterParametersReader* reader, int index) override;

  /**
   * @brief execute Reimplemented from @see AbstractFilter class
   */
  void execute() override;

  /**
   * @brief sendThreadSafeProgressMessage
   * @param counter
   */
  void sendThreadSafeProgressMessage(int64_t counter);

protected:
  ApplyTransformationToGeometry();

  /**
   * @brief applyTransformation
   */
  void applyTransformation();

  /**
   * @brief dataCheck Checks for the appropriate parameter values and availability of arrays
   */
  void dataCheck() override;

  template <class T>
  void linearEquivalent(T& linEquivalent, IDataArray::Pointer linIData, int64_t linIntIndexes, double xt, double yt, double zt)
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
  void linearEquivalentRGB(T linEquivalent[3], IDataArray::Pointer linIData, int64_t linIntIndexes, double xt, double yt, double zt)
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
  bool linearIndexes(double* LinearInterpolationData, int64_t tupleIndex, T& linEquivalent, IDataArray::Pointer linIData)
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
      linearEquivalent<T>(linEquivalent, linIData, linIntIndexes, xt, yt, zt);
      write = true;
    }
    return write;
  }

  template <class T>
  bool linearIndexesRGB(double* LinearInterpolationData, int64_t tupleIndex, T linEquivalent[3], IDataArray::Pointer linIData)
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
      linearEquivalentRGB<T>(linEquivalent, linIData, linIntIndexes, xt, yt, zt);
      write = true;
    }
    return write;
  }

  template <typename T>
  void wrapLinearIndexes(double* LinearInterpolationData, int64_t tupleIndex, typename DataArray<T>::Pointer lin, typename IDataArray::Pointer linData)
  {
    bool wrote = false;
    int index = tupleIndex / 6;

    T linEquivalent = 0;
    typename IDataArray::Pointer linIData = std::dynamic_pointer_cast<IDataArray>(lin);
    wrote = linearIndexes<T>(LinearInterpolationData, tupleIndex, linEquivalent, linIData);
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
  void wrapLinearIndexesRGB(double* LinearInterpolationData, int64_t tupleIndex, typename DataArray<T>::Pointer lin, typename IDataArray::Pointer linData)
  {
    bool wrote = false;
    int index = tupleIndex / 6;
    T linEquivalent[3] = {0, 0, 0};
    typename IDataArray::Pointer linIData = std::dynamic_pointer_cast<IDataArray>(lin);
    wrote = linearIndexesRGB<T>(LinearInterpolationData, tupleIndex, linEquivalent, linIData);
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
  bool applyLinearInterpolation(typename DataArray<T>::Pointer lin, int64_t index, int64_t tupleIndex, double* LinearInterpolationData, typename IDataArray::Pointer linData, bool RGB)
  {
    if(!lin)
    {
      return false;
    }

    if(RGB)
    {
      wrapLinearIndexesRGB<T>(LinearInterpolationData, tupleIndex, lin, linData);
    }
    else
    {
      wrapLinearIndexes<T>(LinearInterpolationData, tupleIndex, lin, linData);
    }
    return true;
  }

  /**
   * @brief setCellAttributeMatrixPath
   * @param value
   */
  void setCellAttributeMatrixPath(const DataArrayPath& value);

  /**
   * @brief getCellAttributeMatrixPath
   */
  DataArrayPath getCellAttributeMatrixPath() const;

  Q_PROPERTY(DataArrayPath CellAttributeMatrixPath READ getCellAttributeMatrixPath WRITE setCellAttributeMatrixPath)

  /**
   * @brief ApplyImageTransformation
   */

  void ApplyImageTransformation();

  void reset()
  {
    m_RotationMatrix.setZero();

    m_Params = ApplyTransformationToGeometry::RotateArgs();
  }

private:
  DataArrayPath m_CellAttributeMatrixPath = {SIMPL::Defaults::DataContainerName, SIMPL::Defaults::CellAttributeMatrixName, ""};

  std::weak_ptr<DataArray<float>> m_TransformationMatrixPtr;
  float* m_TransformationMatrix = nullptr;

  DynamicTableData m_ManualTransformationMatrix = {};
  DataArrayPath m_ComputedTransformationMatrix = {"", "", "TransformationMatrix"};
  //  DataArrayPath m_GeometryToTransform = {"", "", ""};
  int m_TransformationMatrixType = {1};
  int m_InterpolationType = {1};
  FloatVec3Type m_RotationAxis = {};
  float m_RotationAngle = {};
  FloatVec3Type m_Translation = {};
  FloatVec3Type m_Scale = {};
  bool m_UseDataArraySelection = false;
  std::vector<DataArrayPath> m_DataArraySelection = {};
  bool m_SliceBySlice = false;

  FloatArrayType::Pointer m_TransformationReference;

  // Threadsafe Progress Message
  mutable std::mutex m_ProgressMessage_Mutex;
  size_t m_InstanceIndex = {0};
  int64_t m_TotalElements = {};

  Matrix3fR m_RotationMatrix = Matrix3fR::Zero();
  Matrix3fR m_ScalingMatrix = Matrix3fR::Zero();
  MatrixTranslation m_TranslationMatrix = MatrixTranslation::Zero();
  ApplyTransformationToGeometry::RotateArgs m_Params;

public:
  ApplyTransformationToGeometry(const ApplyTransformationToGeometry&) = delete;            // Copy Constructor Not Implemented
  ApplyTransformationToGeometry(ApplyTransformationToGeometry&&) = delete;                 // Move Constructor Not Implemented
  ApplyTransformationToGeometry& operator=(const ApplyTransformationToGeometry&) = delete; // Copy Assignment Not Implemented
  ApplyTransformationToGeometry& operator=(ApplyTransformationToGeometry&&) = delete;
  // Move Assignment Not Implemented
};
