/* ============================================================================
 * Copyright (c) 2011 Michael A. Jackson (BlueQuartz Software)
 * Copyright (c) 2011 Dr. Michael A. Groeber (US Air Force Research Laboratories)
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
 * Neither the name of Michael A. Groeber, Michael A. Jackson,
 * the US Air Force, BlueQuartz Software nor the names of its contributors
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

#include <QtCore/QString>
#include <set>

#include "SIMPLib/SIMPLib.h"
#include "SIMPLib/Common/Constants.h"
#include "SIMPLib/DataArrays/IDataArray.h"
#include "SIMPLib/DataContainers/DataContainer.h"
#include "SIMPLib/FilterParameters/FloatVec3FilterParameter.h"
#include "SIMPLib/Filtering/AbstractFilter.h"
#include "SIMPLib/Geometry/MeshStructs.h"
#include "SIMPLib/DataArrays/DataArray.hpp"

#include "DREAM3DReview/DREAM3DReviewConstants.h"

#include "DREAM3DReview/DREAM3DReviewDLLExport.h"

/**
 * @class LocalDislocationDensityCalculator LocalDislocationDensityCalculator.h /FilterCategoryFilters/LocalDislocationDensityCalculator.h
 * @brief
 * @author
 * @date
 * @version 1.0
 */
class DREAM3DReview_EXPORT LocalDislocationDensityCalculator : public AbstractFilter
{
  Q_OBJECT

#ifdef SIMPL_ENABLE_PYTHON
  PYB11_CREATE_BINDINGS(LocalDislocationDensityCalculator SUPERCLASS AbstractFilter)
  PYB11_SHARED_POINTERS(LocalDislocationDensityCalculator)
  PYB11_FILTER_NEW_MACRO(LocalDislocationDensityCalculator)
  PYB11_FILTER_PARAMETER(DataArrayPath, EdgeDataContainerName)
  PYB11_FILTER_PARAMETER(DataArrayPath, BurgersVectorsArrayPath)
  PYB11_FILTER_PARAMETER(DataArrayPath, SlipPlaneNormalsArrayPath)
  PYB11_FILTER_PARAMETER(FloatVec3Type, CellSize)
  PYB11_FILTER_PARAMETER(DataArrayPath, OutputDataContainerName)
  PYB11_FILTER_PARAMETER(QString, OutputAttributeMatrixName)
  PYB11_FILTER_PARAMETER(QString, OutputArrayName)
  PYB11_FILTER_PARAMETER(QString, DominantSystemArrayName)
  PYB11_PROPERTY(DataArrayPath EdgeDataContainerName READ getEdgeDataContainerName WRITE setEdgeDataContainerName)
  PYB11_PROPERTY(DataArrayPath BurgersVectorsArrayPath READ getBurgersVectorsArrayPath WRITE setBurgersVectorsArrayPath)
  PYB11_PROPERTY(DataArrayPath SlipPlaneNormalsArrayPath READ getSlipPlaneNormalsArrayPath WRITE setSlipPlaneNormalsArrayPath)
  PYB11_PROPERTY(FloatVec3Type CellSize READ getCellSize WRITE setCellSize)
  PYB11_PROPERTY(DataArrayPath OutputDataContainerName READ getOutputDataContainerName WRITE setOutputDataContainerName)
  PYB11_PROPERTY(QString OutputAttributeMatrixName READ getOutputAttributeMatrixName WRITE setOutputAttributeMatrixName)
  PYB11_PROPERTY(QString OutputArrayName READ getOutputArrayName WRITE setOutputArrayName)
  PYB11_PROPERTY(QString DominantSystemArrayName READ getDominantSystemArrayName WRITE setDominantSystemArrayName)
#endif

public:
  using Self = LocalDislocationDensityCalculator;
  using Pointer = std::shared_ptr<Self>;
  using ConstPointer = std::shared_ptr<const Self>;
  using WeakPointer = std::weak_ptr<Self>;
  using ConstWeakPointer = std::weak_ptr<const Self>;
  static Pointer NullPointer();

  static std::shared_ptr<LocalDislocationDensityCalculator> New();

  /**
   * @brief Returns the name of the class for LocalDislocationDensityCalculator
   */
  QString getNameOfClass() const override;
  /**
   * @brief Returns the name of the class for LocalDislocationDensityCalculator
   */
  static QString ClassName();

  ~LocalDislocationDensityCalculator() override;
  /**
   * @brief Setter property for EdgeDataContainerName
   */
  void setEdgeDataContainerName(const DataArrayPath& value);
  /**
   * @brief Getter property for EdgeDataContainerName
   * @return Value of EdgeDataContainerName
   */
  DataArrayPath getEdgeDataContainerName() const;

  Q_PROPERTY(DataArrayPath EdgeDataContainerName READ getEdgeDataContainerName WRITE setEdgeDataContainerName)

  /**
   * @brief Setter property for BurgersVectorsArrayPath
   */
  void setBurgersVectorsArrayPath(const DataArrayPath& value);
  /**
   * @brief Getter property for BurgersVectorsArrayPath
   * @return Value of BurgersVectorsArrayPath
   */
  DataArrayPath getBurgersVectorsArrayPath() const;

  Q_PROPERTY(DataArrayPath BurgersVectorsArrayPath READ getBurgersVectorsArrayPath WRITE setBurgersVectorsArrayPath)

  /**
   * @brief Setter property for SlipPlaneNormalsArrayPath
   */
  void setSlipPlaneNormalsArrayPath(const DataArrayPath& value);
  /**
   * @brief Getter property for SlipPlaneNormalsArrayPath
   * @return Value of SlipPlaneNormalsArrayPath
   */
  DataArrayPath getSlipPlaneNormalsArrayPath() const;

  Q_PROPERTY(DataArrayPath SlipPlaneNormalsArrayPath READ getSlipPlaneNormalsArrayPath WRITE setSlipPlaneNormalsArrayPath)

  /**
   * @brief Setter property for CellSize
   */
  void setCellSize(const FloatVec3Type& value);
  /**
   * @brief Getter property for CellSize
   * @return Value of CellSize
   */
  FloatVec3Type getCellSize() const;

  Q_PROPERTY(FloatVec3Type CellSize READ getCellSize WRITE setCellSize)

  // The user selects a new DataContainerName
  /**
   * @brief Setter property for OutputDataContainerName
   */
  void setOutputDataContainerName(const DataArrayPath& value);
  /**
   * @brief Getter property for OutputDataContainerName
   * @return Value of OutputDataContainerName
   */
  DataArrayPath getOutputDataContainerName() const;

  Q_PROPERTY(DataArrayPath OutputDataContainerName READ getOutputDataContainerName WRITE setOutputDataContainerName)
  // Name the new AttributeMatrix that will get created
  /**
   * @brief Setter property for OutputAttributeMatrixName
   */
  void setOutputAttributeMatrixName(const QString& value);
  /**
   * @brief Getter property for OutputAttributeMatrixName
   * @return Value of OutputAttributeMatrixName
   */
  QString getOutputAttributeMatrixName() const;

  Q_PROPERTY(QString OutputAttributeMatrixName READ getOutputAttributeMatrixName WRITE setOutputAttributeMatrixName)

  // Give the created data array a name
  /**
   * @brief Setter property for OutputArrayName
   */
  void setOutputArrayName(const QString& value);
  /**
   * @brief Getter property for OutputArrayName
   * @return Value of OutputArrayName
   */
  QString getOutputArrayName() const;

  Q_PROPERTY(QString OutputArrayName READ getOutputArrayName WRITE setOutputArrayName)

  /**
   * @brief Setter property for DominantSystemArrayName
   */
  void setDominantSystemArrayName(const QString& value);
  /**
   * @brief Getter property for DominantSystemArrayName
   * @return Value of DominantSystemArrayName
   */
  QString getDominantSystemArrayName() const;

  Q_PROPERTY(QString DominantSystemArrayName READ getDominantSystemArrayName WRITE setDominantSystemArrayName)

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
   * @brief This method will instantiate all the end user settable options/parameters
   * for this filter
   */
  void setupFilterParameters() override;

  /**
   * @brief This method will read the options from a file
   * @param reader The reader that is used to read the options from a file
   */
  void readFilterParameters(AbstractFilterParametersReader* reader, int index) override;

  /**
   * @brief Reimplemented from @see AbstractFilter class
   */
  void execute() override;

  /**
   * @brief This function runs some sanity checks on the DataContainer and inputs
   * in an attempt to ensure the filter can process the inputs.
   */
  void preflight() override;

protected:
  LocalDislocationDensityCalculator();

  /**
   * @brief dataCheck Checks for the appropriate parameter values and availability of arrays
   */
  void dataCheck();

  /**
   * @brief Initializes all the private instance variables.
   */
  void initialize();

  void updateCellInstancePointers();

  int determine_slip_system(int edgeNum);

private:
  std::weak_ptr<DataArray<float>> m_OutputArrayPtr;
  float* m_OutputArray = nullptr;
  std::weak_ptr<DataArray<float>> m_DominantSystemArrayPtr;
  float* m_DominantSystemArray = nullptr;
  std::weak_ptr<DataArray<float>> m_DomainBoundsPtr;
  float* m_DomainBounds = nullptr;
  std::weak_ptr<DataArray<float>> m_BurgersVectorsPtr;
  float* m_BurgersVectors = nullptr;
  std::weak_ptr<DataArray<float>> m_SlipPlaneNormalsPtr;
  float* m_SlipPlaneNormals = nullptr;

  DataArrayPath m_EdgeDataContainerName = {};
  DataArrayPath m_BurgersVectorsArrayPath = {};
  DataArrayPath m_SlipPlaneNormalsArrayPath = {};
  FloatVec3Type m_CellSize = {};
  DataArrayPath m_OutputDataContainerName = {};
  QString m_OutputAttributeMatrixName = {};
  QString m_OutputArrayName = {};
  QString m_DominantSystemArrayName = {};

public:
  LocalDislocationDensityCalculator(const LocalDislocationDensityCalculator&) = delete;            // Copy Constructor Not Implemented
  LocalDislocationDensityCalculator(LocalDislocationDensityCalculator&&) = delete;                 // Move Constructor Not Implemented
  LocalDislocationDensityCalculator& operator=(const LocalDislocationDensityCalculator&) = delete; // Copy Assignment Not Implemented
  LocalDislocationDensityCalculator& operator=(LocalDislocationDensityCalculator&&) = delete;      // Move Assignment Not Implemented
};
