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
#include "SIMPLib/DataArrays/DataArray.hpp"

#include "DREAM3DReview/DREAM3DReviewDLLExport.h"

/**
 * @brief The FindElementCentroids class. See [Filter documentation](@ref findelementcentroids) for details.
 */
class DREAM3DReview_EXPORT FindElementCentroids : public AbstractFilter
{
  Q_OBJECT

  // Start Python bindings declarations
  PYB11_BEGIN_BINDINGS(FindElementCentroids SUPERCLASS AbstractFilter)
  PYB11_FILTER()
  PYB11_SHARED_POINTERS(FindElementCentroids)
  PYB11_FILTER_NEW_MACRO(FindElementCentroids)
  PYB11_PROPERTY(DataArrayPath CellCentroidsArrayPath READ getCellCentroidsArrayPath WRITE setCellCentroidsArrayPath)
  PYB11_PROPERTY(bool CreateVertexDataContainer READ getCreateVertexDataContainer WRITE setCreateVertexDataContainer)
  PYB11_PROPERTY(DataArrayPath NewDataContainerName READ getNewDataContainerName WRITE setNewDataContainerName)
  PYB11_PROPERTY(QString VertexAttributeMatrixName READ getVertexAttributeMatrixName WRITE setVertexAttributeMatrixName)
  PYB11_END_BINDINGS()
  // End Python bindings declarations

public:
  using Self = FindElementCentroids;
  using Pointer = std::shared_ptr<Self>;
  using ConstPointer = std::shared_ptr<const Self>;
  using WeakPointer = std::weak_ptr<Self>;
  using ConstWeakPointer = std::weak_ptr<const Self>;
  static Pointer NullPointer();

  static std::shared_ptr<FindElementCentroids> New();

  /**
   * @brief Returns the name of the class for FindElementCentroids
   */
  QString getNameOfClass() const override;
  /**
   * @brief Returns the name of the class for FindElementCentroids
   */
  static QString ClassName();

  ~FindElementCentroids() override;

  /**
   * @brief Setter property for CellCentroidsArrayPath
   */
  void setCellCentroidsArrayPath(const DataArrayPath& value);
  /**
   * @brief Getter property for CellCentroidsArrayPath
   * @return Value of CellCentroidsArrayPath
   */
  DataArrayPath getCellCentroidsArrayPath() const;
  Q_PROPERTY(DataArrayPath CellCentroidsArrayPath READ getCellCentroidsArrayPath WRITE setCellCentroidsArrayPath)

  /**
   * @brief Setter property for CreateVertexDataContainer
   */
  void setCreateVertexDataContainer(bool value);
  /**
   * @brief Getter property for CreateVertexDataContainer
   * @return Value of CreateVertexDataContainer
   */
  bool getCreateVertexDataContainer() const;
  Q_PROPERTY(bool CreateVertexDataContainer READ getCreateVertexDataContainer WRITE setCreateVertexDataContainer)

  /**
   * @brief Setter property for NewDataContainerName
   */
  void setNewDataContainerName(const DataArrayPath& value);
  /**
   * @brief Getter property for NewDataContainerName
   * @return Value of NewDataContainerName
   */
  DataArrayPath getNewDataContainerName() const;
  Q_PROPERTY(DataArrayPath NewDataContainerName READ getNewDataContainerName WRITE setNewDataContainerName)

  /**
   * @brief Setter property for VertexAttributeMatrixName
   */
  void setVertexAttributeMatrixName(const QString& value);
  /**
   * @brief Getter property for VertexAttributeMatrixName
   * @return Value of VertexAttributeMatrixName
   */
  QString getVertexAttributeMatrixName() const;
  Q_PROPERTY(QString VertexAttributeMatrixName READ getVertexAttributeMatrixName WRITE setVertexAttributeMatrixName)

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
  FindElementCentroids();
  /**
   * @brief dataCheck Checks for the appropriate parameter values and availability of arrays
   */
  void dataCheck() override;

  /**
   * @brief Initializes all the private instance variables.
   */
  void initialize();

private:
  std::weak_ptr<DataArray<float>> m_CellCentroidsArrayPtr;
  float* m_CellCentroidsArray = nullptr;

  DataArrayPath m_CellCentroidsArrayPath = {};
  bool m_CreateVertexDataContainer = {};
  DataArrayPath m_NewDataContainerName = {};
  QString m_VertexAttributeMatrixName = {};

public:
  FindElementCentroids(const FindElementCentroids&) = delete; // Copy Constructor Not Implemented
  FindElementCentroids(FindElementCentroids&&) = delete;      // Move Constructor Not Implemented
  FindElementCentroids& operator=(const FindElementCentroids&) = delete; // Copy Assignment Not Implemented
  FindElementCentroids& operator=(FindElementCentroids&&) = delete;      // Move Assignment Not Implemented
};

