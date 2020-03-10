/* Copyright (c) 2009-2016 BlueQuartz Software, LLC
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

#pragma once

#include <memory>

#include "SIMPLib/SIMPLib.h"
#include "SIMPLib/Filtering/AbstractFilter.h"
#include "SIMPLib/Math/SIMPLibMath.h"
#include "SIMPLib/DataArrays/DataArray.hpp"

class IDataArray;
using IDataArrayWkPtrType = std::weak_ptr<IDataArray>;

#include "DREAM3DReview/DREAM3DReviewDLLExport.h"

/**
 * @brief The DBSCAN class. See [Filter documentation](@ref DBSCAN) for details.
 */
class DREAM3DReview_EXPORT DBSCAN : public AbstractFilter
{
  Q_OBJECT

#ifdef SIMPL_ENABLE_PYTHON
  PYB11_CREATE_BINDINGS(DBSCAN SUPERCLASS AbstractFilter)
  PYB11_SHARED_POINTERS(DBSCAN)
  PYB11_FILTER_NEW_MACRO(DBSCAN)
  PYB11_FILTER_PARAMETER(DataArrayPath, SelectedArrayPath)
  PYB11_FILTER_PARAMETER(bool, UseMask)
  PYB11_FILTER_PARAMETER(DataArrayPath, MaskArrayPath)
  PYB11_FILTER_PARAMETER(QString, FeatureIdsArrayName)
  PYB11_FILTER_PARAMETER(QString, FeatureAttributeMatrixName)
  PYB11_FILTER_PARAMETER(float, Epsilon)
  PYB11_FILTER_PARAMETER(int, MinPnts)
  PYB11_FILTER_PARAMETER(int, DistanceMetric)
  PYB11_PROPERTY(DataArrayPath SelectedArrayPath READ getSelectedArrayPath WRITE setSelectedArrayPath)
  PYB11_PROPERTY(bool UseMask READ getUseMask WRITE setUseMask)
  PYB11_PROPERTY(DataArrayPath MaskArrayPath READ getMaskArrayPath WRITE setMaskArrayPath)
  PYB11_PROPERTY(QString FeatureIdsArrayName READ getFeatureIdsArrayName WRITE setFeatureIdsArrayName)
  PYB11_PROPERTY(QString FeatureAttributeMatrixName READ getFeatureAttributeMatrixName WRITE setFeatureAttributeMatrixName)
  PYB11_PROPERTY(float Epsilon READ getEpsilon WRITE setEpsilon)
  PYB11_PROPERTY(int MinPnts READ getMinPnts WRITE setMinPnts)
  PYB11_PROPERTY(int DistanceMetric READ getDistanceMetric WRITE setDistanceMetric)
#endif

public:
  using Self = DBSCAN;
  using Pointer = std::shared_ptr<Self>;
  using ConstPointer = std::shared_ptr<const Self>;
  using WeakPointer = std::weak_ptr<Self>;
  using ConstWeakPointer = std::weak_ptr<const Self>;
  static Pointer NullPointer();

  static std::shared_ptr<DBSCAN> New();

  /**
   * @brief Returns the name of the class for DBSCAN
   */
  QString getNameOfClass() const override;
  /**
   * @brief Returns the name of the class for DBSCAN
   */
  static QString ClassName();

  ~DBSCAN() override;

  /**
   * @brief Setter property for SelectedArrayPath
   */
  void setSelectedArrayPath(const DataArrayPath& value);
  /**
   * @brief Getter property for SelectedArrayPath
   * @return Value of SelectedArrayPath
   */
  DataArrayPath getSelectedArrayPath() const;

  Q_PROPERTY(DataArrayPath SelectedArrayPath READ getSelectedArrayPath WRITE setSelectedArrayPath)

  /**
   * @brief Setter property for UseMask
   */
  void setUseMask(bool value);
  /**
   * @brief Getter property for UseMask
   * @return Value of UseMask
   */
  bool getUseMask() const;

  Q_PROPERTY(bool UseMask READ getUseMask WRITE setUseMask)

  /**
   * @brief Setter property for MaskArrayPath
   */
  void setMaskArrayPath(const DataArrayPath& value);
  /**
   * @brief Getter property for MaskArrayPath
   * @return Value of MaskArrayPath
   */
  DataArrayPath getMaskArrayPath() const;

  Q_PROPERTY(DataArrayPath MaskArrayPath READ getMaskArrayPath WRITE setMaskArrayPath)

  /**
   * @brief Setter property for FeatureIdsArrayName
   */
  void setFeatureIdsArrayName(const QString& value);
  /**
   * @brief Getter property for FeatureIdsArrayName
   * @return Value of FeatureIdsArrayName
   */
  QString getFeatureIdsArrayName() const;

  Q_PROPERTY(QString FeatureIdsArrayName READ getFeatureIdsArrayName WRITE setFeatureIdsArrayName)

  /**
   * @brief Setter property for FeatureAttributeMatrixName
   */
  void setFeatureAttributeMatrixName(const QString& value);
  /**
   * @brief Getter property for FeatureAttributeMatrixName
   * @return Value of FeatureAttributeMatrixName
   */
  QString getFeatureAttributeMatrixName() const;

  Q_PROPERTY(QString FeatureAttributeMatrixName READ getFeatureAttributeMatrixName WRITE setFeatureAttributeMatrixName)

  /**
   * @brief Setter property for Epsilon
   */
  void setEpsilon(float value);
  /**
   * @brief Getter property for Epsilon
   * @return Value of Epsilon
   */
  float getEpsilon() const;

  Q_PROPERTY(float Epsilon READ getEpsilon WRITE setEpsilon)

  /**
   * @brief Setter property for MinPnts
   */
  void setMinPnts(int value);
  /**
   * @brief Getter property for MinPnts
   * @return Value of MinPnts
   */
  int getMinPnts() const;

  Q_PROPERTY(float MinPnts READ getMinPnts WRITE setMinPnts)

  /**
   * @brief Setter property for DistanceMetric
   */
  void setDistanceMetric(int value);
  /**
   * @brief Getter property for DistanceMetric
   * @return Value of DistanceMetric
   */
  int getDistanceMetric() const;

  Q_PROPERTY(int DistanceMetric READ getDistanceMetric WRITE setDistanceMetric)

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
  DBSCAN();

  /**
   * @brief dataCheck Checks for the appropriate parameter values and availability of arrays
   */
  void dataCheck() override;

  /**
   * @brief Initializes all the private instance variables.
   */
  void initialize();

private:
  IDataArrayWkPtrType m_InDataPtr;
  void* m_InData = nullptr;

  std::weak_ptr<DataArray<bool>> m_MaskPtr;
  bool* m_Mask = nullptr;
  std::weak_ptr<DataArray<int32_t>> m_FeatureIdsPtr;
  int32_t* m_FeatureIds = nullptr;

  DataArrayPath m_SelectedArrayPath = {};
  bool m_UseMask = {};
  DataArrayPath m_MaskArrayPath = {};
  QString m_FeatureIdsArrayName = {};
  QString m_FeatureAttributeMatrixName = {};
  float m_Epsilon = {};
  int m_MinPnts = {};
  int m_DistanceMetric = {};

public:
  DBSCAN(const DBSCAN&) = delete;            // Copy Constructor Not Implemented
  DBSCAN(DBSCAN&&) = delete;                 // Move Constructor Not Implemented
  DBSCAN& operator=(const DBSCAN&) = delete; // Copy Assignment Not Implemented
  DBSCAN& operator=(DBSCAN&&) = delete;      // Move Assignment Not Implemented
};

