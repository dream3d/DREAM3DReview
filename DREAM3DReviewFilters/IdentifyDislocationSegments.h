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
#include "SIMPLib/Filtering/AbstractFilter.h"
#include "SIMPLib/Geometry/MeshStructs.h"
#include "SIMPLib/DataArrays/DataArray.hpp"

#include "DREAM3DReview/DREAM3DReviewDLLExport.h"

#include "DREAM3DReview/DREAM3DReviewConstants.h"
/**
 * @class IdentifyDislocationSegments IdentifyDislocationSegments.h /FilterCategoryFilters/IdentifyDislocationSegments.h
 * @brief
 * @author
 * @date
 * @version 1.0
 */
class DREAM3DReview_EXPORT IdentifyDislocationSegments : public AbstractFilter
{
  Q_OBJECT

#ifdef SIMPL_ENABLE_PYTHON
  PYB11_CREATE_BINDINGS(IdentifyDislocationSegments SUPERCLASS AbstractFilter)
  PYB11_SHARED_POINTERS(IdentifyDislocationSegments)
  PYB11_FILTER_NEW_MACRO(IdentifyDislocationSegments)
  PYB11_FILTER_PARAMETER(DataArrayPath, BurgersVectorsArrayPath)
  PYB11_FILTER_PARAMETER(DataArrayPath, SlipPlaneNormalsArrayPath)
  PYB11_FILTER_PARAMETER(QString, DislocationIdsArrayName)
  PYB11_FILTER_PARAMETER(QString, ActiveArrayName)
  PYB11_PROPERTY(DataArrayPath BurgersVectorsArrayPath READ getBurgersVectorsArrayPath WRITE setBurgersVectorsArrayPath)
  PYB11_PROPERTY(DataArrayPath SlipPlaneNormalsArrayPath READ getSlipPlaneNormalsArrayPath WRITE setSlipPlaneNormalsArrayPath)
  PYB11_PROPERTY(QString DislocationIdsArrayName READ getDislocationIdsArrayName WRITE setDislocationIdsArrayName)
  PYB11_PROPERTY(QString ActiveArrayName READ getActiveArrayName WRITE setActiveArrayName)
#endif

public:
  using Self = IdentifyDislocationSegments;
  using Pointer = std::shared_ptr<Self>;
  using ConstPointer = std::shared_ptr<const Self>;
  using WeakPointer = std::weak_ptr<Self>;
  using ConstWeakPointer = std::weak_ptr<const Self>;
  static Pointer NullPointer();

  static std::shared_ptr<IdentifyDislocationSegments> New();

  /**
   * @brief Returns the name of the class for IdentifyDislocationSegments
   */
  QString getNameOfClass() const override;
  /**
   * @brief Returns the name of the class for IdentifyDislocationSegments
   */
  static QString ClassName();

  ~IdentifyDislocationSegments() override;
  /**
   * @brief Setter property for EdgeFeatureAttributeMatrixName
   */
  void setEdgeFeatureAttributeMatrixName(const QString& value);
  /**
   * @brief Getter property for EdgeFeatureAttributeMatrixName
   * @return Value of EdgeFeatureAttributeMatrixName
   */
  QString getEdgeFeatureAttributeMatrixName() const;

  Q_PROPERTY(QString EdgeFeatureAttributeMatrixName READ getEdgeFeatureAttributeMatrixName WRITE setEdgeFeatureAttributeMatrixName)

  /**
   * @brief This returns the group that the filter belonds to. You can select
   * a different group if you want. The string returned here will be displayed
   * in the GUI for the filter
   */
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
   * @brief Setter property for DislocationIdsArrayName
   */
  void setDislocationIdsArrayName(const QString& value);
  /**
   * @brief Getter property for DislocationIdsArrayName
   * @return Value of DislocationIdsArrayName
   */
  QString getDislocationIdsArrayName() const;

  Q_PROPERTY(QString DislocationIdsArrayName READ getDislocationIdsArrayName WRITE setDislocationIdsArrayName)

  /**
   * @brief Setter property for ActiveArrayName
   */
  void setActiveArrayName(const QString& value);
  /**
   * @brief Getter property for ActiveArrayName
   * @return Value of ActiveArrayName
   */
  QString getActiveArrayName() const;

  Q_PROPERTY(QString ActiveArrayName READ getActiveArrayName WRITE setActiveArrayName)

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
  IdentifyDislocationSegments();

  /**
   * @brief dataCheck Checks for the appropriate parameter values and availability of arrays
   */
  void dataCheck() override;

  /**
   * @brief Initializes all the private instance variables.
   */
  void initialize();

  void updateEdgeFeatureInstancePointers();

private:
  std::weak_ptr<DataArray<float>> m_BurgersVectorsPtr;
  float* m_BurgersVectors = nullptr;
  std::weak_ptr<DataArray<float>> m_SlipPlaneNormalsPtr;
  float* m_SlipPlaneNormals = nullptr;
  std::weak_ptr<DataArray<int32_t>> m_DislocationIdsPtr;
  int32_t* m_DislocationIds = nullptr;
  std::weak_ptr<DataArray<bool>> m_ActivePtr;
  bool* m_Active = nullptr;

  QString m_EdgeFeatureAttributeMatrixName = {};
  DataArrayPath m_BurgersVectorsArrayPath = {};
  DataArrayPath m_SlipPlaneNormalsArrayPath = {};
  QString m_DislocationIdsArrayName = {};
  QString m_ActiveArrayName = {};

public:
  IdentifyDislocationSegments(const IdentifyDislocationSegments&) = delete;            // Copy Constructor Not Implemented
  IdentifyDislocationSegments(IdentifyDislocationSegments&&) = delete;                 // Move Constructor Not Implemented
  IdentifyDislocationSegments& operator=(const IdentifyDislocationSegments&) = delete; // Copy Assignment Not Implemented
  IdentifyDislocationSegments& operator=(IdentifyDislocationSegments&&) = delete;      // Move Assignment Not Implemented
};
