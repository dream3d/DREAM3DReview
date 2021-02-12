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
#include "SIMPLib/DataArrays/DataArray.hpp"
#include "SIMPLib/DataArrays/IDataArray.h"
#include "SIMPLib/DataContainers/DataContainer.h"
#include "SIMPLib/FilterParameters/FloatVec3FilterParameter.h"
#include "SIMPLib/Filtering/AbstractFilter.h"
#include "SIMPLib/Geometry/MeshStructs.h"

#include "DREAM3DReview/DREAM3DReviewConstants.h"
#include "DREAM3DReview/DREAM3DReviewDLLExport.h"

/**
 * @class DiscretizeDDDomain DiscretizeDDDomain.h /FilterCategoryFilters/DiscretizeDDDomain.h
 * @brief
 * @author
 * @date
 * @version 1.0
 */
class DREAM3DReview_EXPORT DiscretizeDDDomain : public AbstractFilter
{
  Q_OBJECT

  // Start Python bindings declarations
  PYB11_BEGIN_BINDINGS(DiscretizeDDDomain SUPERCLASS AbstractFilter)
  PYB11_FILTER()
  PYB11_SHARED_POINTERS(DiscretizeDDDomain)
  PYB11_FILTER_NEW_MACRO(DiscretizeDDDomain)
  PYB11_PROPERTY(DataArrayPath EdgeDataContainerName READ getEdgeDataContainerName WRITE setEdgeDataContainerName)
  PYB11_PROPERTY(FloatVec3Type CellSize READ getCellSize WRITE setCellSize)
  PYB11_PROPERTY(DataArrayPath OutputDataContainerName READ getOutputDataContainerName WRITE setOutputDataContainerName)
  PYB11_PROPERTY(QString OutputAttributeMatrixName READ getOutputAttributeMatrixName WRITE setOutputAttributeMatrixName)
  PYB11_PROPERTY(QString OutputArrayName READ getOutputArrayName WRITE setOutputArrayName)
  PYB11_END_BINDINGS()
  // End Python bindings declarations

public:
  using Self = DiscretizeDDDomain;
  using Pointer = std::shared_ptr<Self>;
  using ConstPointer = std::shared_ptr<const Self>;
  using WeakPointer = std::weak_ptr<Self>;
  using ConstWeakPointer = std::weak_ptr<const Self>;
  static Pointer NullPointer();

  static std::shared_ptr<DiscretizeDDDomain> New();

  /**
   * @brief Returns the name of the class for DiscretizeDDDomain
   */
  QString getNameOfClass() const override;
  /**
   * @brief Returns the name of the class for DiscretizeDDDomain
   */
  static QString ClassName();

  ~DiscretizeDDDomain() override;
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

protected:
  DiscretizeDDDomain();

  /**
   * @brief dataCheck Checks for the appropriate parameter values and availability of arrays
   */
  void dataCheck() override;

  /**
   * @brief Initializes all the private instance variables.
   */
  void initialize();

private:
  std::weak_ptr<DataArray<int32_t>> m_OutputArrayPtr;
  int32_t* m_OutputArray = nullptr;

  DataArrayPath m_EdgeDataContainerName = {SIMPL::Defaults::DataContainerName, "", ""};
  FloatVec3Type m_CellSize = {};
  DataArrayPath m_OutputDataContainerName = {SIMPL::Defaults::NewDataContainerName, "", ""};
  QString m_OutputAttributeMatrixName = {SIMPL::Defaults::CellAttributeMatrixName};
  QString m_OutputArrayName = {"DislocationLineDensity"};

public:
  DiscretizeDDDomain(const DiscretizeDDDomain&) = delete;            // Copy Constructor Not Implemented
  DiscretizeDDDomain(DiscretizeDDDomain&&) = delete;                 // Move Constructor Not Implemented
  DiscretizeDDDomain& operator=(const DiscretizeDDDomain&) = delete; // Copy Assignment Not Implemented
  DiscretizeDDDomain& operator=(DiscretizeDDDomain&&) = delete;      // Move Assignment Not Implemented
};
