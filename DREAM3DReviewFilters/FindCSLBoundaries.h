/* ============================================================================
 * Copyright (c) 2011 Michael A. Jackson (BlueQuartz Software)
 * Copyright (c) 2011 Dr. Michael A. Groeber (US Air Force Research Laboratories)
 * Copyright (c) 2014 Dr. Joseph C. Tucker (UES, Inc.)
 * All rights reserved.
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
 * Neither the name of Joseph C. Tucker, Michael A. Groeber, Michael A. Jackson,
 * UES, Inc., the US Air Force, BlueQuartz Software nor the names of its contributors
 * may be used to endorse or promote products derived from this software without
 * specific prior written permission.
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
 *  This code was written under United States Air Force Contract number
 *                   FA8650-07-D-5800 and FA8650-10-D-5226
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#pragma once

#include <memory>
#include <vector>

#include <QtCore/QString>

#include "SIMPLib/SIMPLib.h"
#include "SIMPLib/DataArrays/IDataArray.h"
#include "SIMPLib/DataContainers/DataContainer.h"
#include "SIMPLib/Filtering/AbstractFilter.h"
#include "SIMPLib/DataArrays/DataArray.hpp"

#include "OrientationLib/LaueOps/CubicOps.h"
#include "OrientationLib/LaueOps/HexagonalOps.h"
#include "OrientationLib/LaueOps/LaueOps.h"
#include "OrientationLib/LaueOps/OrthoRhombicOps.h"

#include "DREAM3DReview/DREAM3DReviewDLLExport.h"


class LaueOps;
using LaueOpsShPtrType = std::shared_ptr<LaueOps>;
using LaueOpsContainer = std::vector<LaueOpsShPtrType>;

/**
 * @class FindCSLBoundaries FindCSLBoundaries.h /TransformationPhase/FindCSLBoundaries.h
 * @brief
 * @author Joseph C. Tucker
 * @date July 28, 2014
 * @version 5.1
 */
class DREAM3DReview_EXPORT FindCSLBoundaries : public AbstractFilter
{
  Q_OBJECT

  // Start Python bindings declarations
  PYB11_BEGIN_BINDINGS(FindCSLBoundaries SUPERCLASS AbstractFilter)
  PYB11_SHARED_POINTERS(FindCSLBoundaries)
  PYB11_FILTER_NEW_MACRO(FindCSLBoundaries)
  PYB11_PROPERTY(float CSL READ getCSL WRITE setCSL)
  PYB11_PROPERTY(float AxisTolerance READ getAxisTolerance WRITE setAxisTolerance)
  PYB11_PROPERTY(float AngleTolerance READ getAngleTolerance WRITE setAngleTolerance)
  PYB11_PROPERTY(DataArrayPath AvgQuatsArrayPath READ getAvgQuatsArrayPath WRITE setAvgQuatsArrayPath)
  PYB11_PROPERTY(DataArrayPath FeaturePhasesArrayPath READ getFeaturePhasesArrayPath WRITE setFeaturePhasesArrayPath)
  PYB11_PROPERTY(DataArrayPath CrystalStructuresArrayPath READ getCrystalStructuresArrayPath WRITE setCrystalStructuresArrayPath)
  PYB11_PROPERTY(DataArrayPath SurfaceMeshFaceLabelsArrayPath READ getSurfaceMeshFaceLabelsArrayPath WRITE setSurfaceMeshFaceLabelsArrayPath)
  PYB11_PROPERTY(DataArrayPath SurfaceMeshFaceNormalsArrayPath READ getSurfaceMeshFaceNormalsArrayPath WRITE setSurfaceMeshFaceNormalsArrayPath)
  PYB11_PROPERTY(QString SurfaceMeshCSLBoundaryArrayName READ getSurfaceMeshCSLBoundaryArrayName WRITE setSurfaceMeshCSLBoundaryArrayName)
  PYB11_PROPERTY(QString SurfaceMeshCSLBoundaryIncoherenceArrayName READ getSurfaceMeshCSLBoundaryIncoherenceArrayName WRITE setSurfaceMeshCSLBoundaryIncoherenceArrayName)
  PYB11_END_BINDINGS()
  // End Python bindings declarations

public:
  using Self = FindCSLBoundaries;
  using Pointer = std::shared_ptr<Self>;
  using ConstPointer = std::shared_ptr<const Self>;
  using WeakPointer = std::weak_ptr<Self>;
  using ConstWeakPointer = std::weak_ptr<const Self>;
  static Pointer NullPointer();

  static std::shared_ptr<FindCSLBoundaries> New();

  /**
   * @brief Returns the name of the class for FindCSLBoundaries
   */
  QString getNameOfClass() const override;
  /**
   * @brief Returns the name of the class for FindCSLBoundaries
   */
  static QString ClassName();

  ~FindCSLBoundaries() override;

  /**
   * @brief Setter property for CSL
   */
  void setCSL(float value);
  /**
   * @brief Getter property for CSL
   * @return Value of CSL
   */
  float getCSL() const;

  Q_PROPERTY(float CSL READ getCSL WRITE setCSL)
  /**
   * @brief Setter property for AxisTolerance
   */
  void setAxisTolerance(float value);
  /**
   * @brief Getter property for AxisTolerance
   * @return Value of AxisTolerance
   */
  float getAxisTolerance() const;

  Q_PROPERTY(float AxisTolerance READ getAxisTolerance WRITE setAxisTolerance)
  /**
   * @brief Setter property for AngleTolerance
   */
  void setAngleTolerance(float value);
  /**
   * @brief Getter property for AngleTolerance
   * @return Value of AngleTolerance
   */
  float getAngleTolerance() const;

  Q_PROPERTY(float AngleTolerance READ getAngleTolerance WRITE setAngleTolerance)

  /**
   * @brief Setter property for AvgQuatsArrayPath
   */
  void setAvgQuatsArrayPath(const DataArrayPath& value);
  /**
   * @brief Getter property for AvgQuatsArrayPath
   * @return Value of AvgQuatsArrayPath
   */
  DataArrayPath getAvgQuatsArrayPath() const;

  Q_PROPERTY(DataArrayPath AvgQuatsArrayPath READ getAvgQuatsArrayPath WRITE setAvgQuatsArrayPath)

  /**
   * @brief Setter property for FeaturePhasesArrayPath
   */
  void setFeaturePhasesArrayPath(const DataArrayPath& value);
  /**
   * @brief Getter property for FeaturePhasesArrayPath
   * @return Value of FeaturePhasesArrayPath
   */
  DataArrayPath getFeaturePhasesArrayPath() const;

  Q_PROPERTY(DataArrayPath FeaturePhasesArrayPath READ getFeaturePhasesArrayPath WRITE setFeaturePhasesArrayPath)

  /**
   * @brief Setter property for CrystalStructuresArrayPath
   */
  void setCrystalStructuresArrayPath(const DataArrayPath& value);
  /**
   * @brief Getter property for CrystalStructuresArrayPath
   * @return Value of CrystalStructuresArrayPath
   */
  DataArrayPath getCrystalStructuresArrayPath() const;

  Q_PROPERTY(DataArrayPath CrystalStructuresArrayPath READ getCrystalStructuresArrayPath WRITE setCrystalStructuresArrayPath)

  /**
   * @brief Setter property for SurfaceMeshFaceLabelsArrayPath
   */
  void setSurfaceMeshFaceLabelsArrayPath(const DataArrayPath& value);
  /**
   * @brief Getter property for SurfaceMeshFaceLabelsArrayPath
   * @return Value of SurfaceMeshFaceLabelsArrayPath
   */
  DataArrayPath getSurfaceMeshFaceLabelsArrayPath() const;

  Q_PROPERTY(DataArrayPath SurfaceMeshFaceLabelsArrayPath READ getSurfaceMeshFaceLabelsArrayPath WRITE setSurfaceMeshFaceLabelsArrayPath)

  /**
   * @brief Setter property for SurfaceMeshFaceNormalsArrayPath
   */
  void setSurfaceMeshFaceNormalsArrayPath(const DataArrayPath& value);
  /**
   * @brief Getter property for SurfaceMeshFaceNormalsArrayPath
   * @return Value of SurfaceMeshFaceNormalsArrayPath
   */
  DataArrayPath getSurfaceMeshFaceNormalsArrayPath() const;

  Q_PROPERTY(DataArrayPath SurfaceMeshFaceNormalsArrayPath READ getSurfaceMeshFaceNormalsArrayPath WRITE setSurfaceMeshFaceNormalsArrayPath)

  /**
   * @brief Setter property for SurfaceMeshCSLBoundaryArrayName
   */
  void setSurfaceMeshCSLBoundaryArrayName(const QString& value);
  /**
   * @brief Getter property for SurfaceMeshCSLBoundaryArrayName
   * @return Value of SurfaceMeshCSLBoundaryArrayName
   */
  QString getSurfaceMeshCSLBoundaryArrayName() const;

  Q_PROPERTY(QString SurfaceMeshCSLBoundaryArrayName READ getSurfaceMeshCSLBoundaryArrayName WRITE setSurfaceMeshCSLBoundaryArrayName)

  /**
   * @brief Setter property for SurfaceMeshCSLBoundaryIncoherenceArrayName
   */
  void setSurfaceMeshCSLBoundaryIncoherenceArrayName(const QString& value);
  /**
   * @brief Getter property for SurfaceMeshCSLBoundaryIncoherenceArrayName
   * @return Value of SurfaceMeshCSLBoundaryIncoherenceArrayName
   */
  QString getSurfaceMeshCSLBoundaryIncoherenceArrayName() const;

  Q_PROPERTY(QString SurfaceMeshCSLBoundaryIncoherenceArrayName READ getSurfaceMeshCSLBoundaryIncoherenceArrayName WRITE setSurfaceMeshCSLBoundaryIncoherenceArrayName)

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

protected:
  FindCSLBoundaries();

  /**
   * @brief dataCheck Checks for the appropriate parameter values and availability of arrays
   */
  void dataCheck() override;

private:
  std::weak_ptr<DataArray<float>> m_AvgQuatsPtr;
  float* m_AvgQuats = nullptr;
  std::weak_ptr<DataArray<int32_t>> m_FeaturePhasesPtr;
  int32_t* m_FeaturePhases = nullptr;
  std::weak_ptr<DataArray<unsigned int>> m_CrystalStructuresPtr;
  unsigned int* m_CrystalStructures = nullptr;
  std::weak_ptr<DataArray<int32_t>> m_SurfaceMeshFaceLabelsPtr;
  int32_t* m_SurfaceMeshFaceLabels = nullptr;
  std::weak_ptr<DataArray<double>> m_SurfaceMeshFaceNormalsPtr;
  double* m_SurfaceMeshFaceNormals = nullptr;
  std::weak_ptr<DataArray<bool>> m_SurfaceMeshCSLBoundaryPtr;
  bool* m_SurfaceMeshCSLBoundary = nullptr;
  std::weak_ptr<DataArray<float>> m_SurfaceMeshCSLBoundaryIncoherencePtr;
  float* m_SurfaceMeshCSLBoundaryIncoherence = nullptr;

  float m_CSL = {};
  float m_AxisTolerance = {};
  float m_AngleTolerance = {};
  DataArrayPath m_AvgQuatsArrayPath = {};
  DataArrayPath m_FeaturePhasesArrayPath = {};
  DataArrayPath m_CrystalStructuresArrayPath = {};
  DataArrayPath m_SurfaceMeshFaceLabelsArrayPath = {};
  DataArrayPath m_SurfaceMeshFaceNormalsArrayPath = {};
  QString m_SurfaceMeshCSLBoundaryArrayName = {};
  QString m_SurfaceMeshCSLBoundaryIncoherenceArrayName = {};

  LaueOpsContainer m_OrientationOps;  
  CubicOps::Pointer m_CubicOps;
  HexagonalOps::Pointer m_HexOps;
  OrthoRhombicOps::Pointer m_OrthoOps;

  void dataCheckVoxel();
  void dataCheckSurfaceMesh();

public:
  FindCSLBoundaries(const FindCSLBoundaries&) = delete;            // Copy Constructor Not Implemented
  FindCSLBoundaries(FindCSLBoundaries&&) = delete;                 // Move Constructor Not Implemented
  FindCSLBoundaries& operator=(const FindCSLBoundaries&) = delete; // Copy Assignment Not Implemented
  FindCSLBoundaries& operator=(FindCSLBoundaries&&) = delete;      // Move Assignment Not Implemented
};
