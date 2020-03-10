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
 * @brief The ComputeUmeyamaTransform class. See [Filter documentation](@ref computeumeyamatransform) for details.
 */
class DREAM3DReview_EXPORT ComputeUmeyamaTransform : public AbstractFilter
{
  Q_OBJECT

#ifdef SIMPL_ENABLE_PYTHON
  PYB11_CREATE_BINDINGS(ComputeUmeyamaTransform SUPERCLASS AbstractFilter)
  PYB11_SHARED_POINTERS(ComputeUmeyamaTransform)
  PYB11_FILTER_NEW_MACRO(ComputeUmeyamaTransform)
  PYB11_FILTER_PARAMETER(DataArrayPath, SourcePointSet)
  PYB11_FILTER_PARAMETER(DataArrayPath, DestPointSet)
  PYB11_FILTER_PARAMETER(bool, UseScaling)
  PYB11_FILTER_PARAMETER(QString, TransformationAttributeMatrixName)
  PYB11_FILTER_PARAMETER(QString, TransformationMatrixName)
  PYB11_PROPERTY(DataArrayPath SourcePointSet READ getSourcePointSet WRITE setSourcePointSet)
  PYB11_PROPERTY(DataArrayPath DestPointSet READ getDestPointSet WRITE setDestPointSet)
  PYB11_PROPERTY(bool UseScaling READ getUseScaling WRITE setUseScaling)
  PYB11_PROPERTY(QString TransformationAttributeMatrixName READ getTransformationAttributeMatrixName WRITE setTransformationAttributeMatrixName)
  PYB11_PROPERTY(QString TransformationMatrixName READ getTransformationMatrixName WRITE setTransformationMatrixName)
#endif

public:
  using Self = ComputeUmeyamaTransform;
  using Pointer = std::shared_ptr<Self>;
  using ConstPointer = std::shared_ptr<const Self>;
  using WeakPointer = std::weak_ptr<Self>;
  using ConstWeakPointer = std::weak_ptr<const Self>;
  static Pointer NullPointer();

  static std::shared_ptr<ComputeUmeyamaTransform> New();

  /**
   * @brief Returns the name of the class for ComputeUmeyamaTransform
   */
  QString getNameOfClass() const override;
  /**
   * @brief Returns the name of the class for ComputeUmeyamaTransform
   */
  static QString ClassName();

  ~ComputeUmeyamaTransform() override;

  /**
   * @brief Setter property for SourcePointSet
   */
  void setSourcePointSet(const DataArrayPath& value);
  /**
   * @brief Getter property for SourcePointSet
   * @return Value of SourcePointSet
   */
  DataArrayPath getSourcePointSet() const;

  Q_PROPERTY(DataArrayPath SourcePointSet READ getSourcePointSet WRITE setSourcePointSet)

  /**
   * @brief Setter property for DestPointSet
   */
  void setDestPointSet(const DataArrayPath& value);
  /**
   * @brief Getter property for DestPointSet
   * @return Value of DestPointSet
   */
  DataArrayPath getDestPointSet() const;

  Q_PROPERTY(DataArrayPath DestPointSet READ getDestPointSet WRITE setDestPointSet)

  /**
   * @brief Setter property for UseScaling
   */
  void setUseScaling(bool value);
  /**
   * @brief Getter property for UseScaling
   * @return Value of UseScaling
   */
  bool getUseScaling() const;

  Q_PROPERTY(bool UseScaling READ getUseScaling WRITE setUseScaling)

  /**
   * @brief Setter property for TransformationAttributeMatrixName
   */
  void setTransformationAttributeMatrixName(const QString& value);
  /**
   * @brief Getter property for TransformationAttributeMatrixName
   * @return Value of TransformationAttributeMatrixName
   */
  QString getTransformationAttributeMatrixName() const;

  Q_PROPERTY(QString TransformationAttributeMatrixName READ getTransformationAttributeMatrixName WRITE setTransformationAttributeMatrixName)

  /**
   * @brief Setter property for TransformationMatrixName
   */
  void setTransformationMatrixName(const QString& value);
  /**
   * @brief Getter property for TransformationMatrixName
   * @return Value of TransformationMatrixName
   */
  QString getTransformationMatrixName() const;

  Q_PROPERTY(QString TransformationMatrixName READ getTransformationMatrixName WRITE setTransformationMatrixName)

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
  ComputeUmeyamaTransform();

  /**
   * @brief dataCheck Checks for the appropriate parameter values and availability of arrays
   */
  void dataCheck() override;

  /**
   * @brief Initializes all the private instance variables.
   */
  void initialize();

private:
  std::weak_ptr<DataArray<float>> m_TransformationMatrixPtr;
  float* m_TransformationMatrix = nullptr;

  DataArrayPath m_SourcePointSet = {};
  DataArrayPath m_DestPointSet = {};
  bool m_UseScaling = {};
  QString m_TransformationAttributeMatrixName = {};
  QString m_TransformationMatrixName = {};

public:
  ComputeUmeyamaTransform(const ComputeUmeyamaTransform&) = delete; // Copy Constructor Not Implemented
  ComputeUmeyamaTransform(ComputeUmeyamaTransform&&) = delete;      // Move Constructor Not Implemented
  ComputeUmeyamaTransform& operator=(const ComputeUmeyamaTransform&) = delete; // Copy Assignment Not Implemented
  ComputeUmeyamaTransform& operator=(ComputeUmeyamaTransform&&) = delete;      // Move Assignment Not Implemented
};

