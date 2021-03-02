/* ============================================================================
 * Copyright 2021 The University of Utah
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
 * Neither the name of BlueQuartz Software, the US Air Force, the University of Utah nor the names of its contributors may be used to endorse or promote products derived from this software without
 * specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY
 * OF SUCH DAMAGE.
 *
 * The code contained herein was partially funded by the following contracts:
 *
 *
 * This code contained herein is based upon work supported by the following grants:
 *    DOE Office of Nuclear Energy's Nuclear Energy University Program Grant No.: DE-NE0008799
 *    DOD Office of Economic Adjustment Grant No.: ST1605-19-03
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
#pragma once

#include "SIMPLib/SIMPLib.h"
#include "SIMPLib/DataArrays/DataArray.hpp"
#include "SIMPLib/Filtering/AbstractFilter.h"

#include "DREAM3DReview/DREAM3DReviewPlugin.h"

/**
 * @brief The ComputeFeatureEigenstrains class. See [Filter documentation](@ref computefeatureeigenstrains) for details.
 */
class DREAM3DReview_EXPORT ComputeFeatureEigenstrains : public AbstractFilter
{
  Q_OBJECT

  // clang-format off
  PYB11_BEGIN_BINDINGS(ComputeFeatureEigenstrains SUPERCLASS AbstractFilter)
  PYB11_FILTER()
  PYB11_SHARED_POINTERS(ComputeFeatureEigenstrains)
  PYB11_FILTER_NEW_MACRO(ComputeFeatureEigenstrains)

  PYB11_PROPERTY(float PoissonRatio READ getPoissonRatio WRITE setPoissonRatio)
  PYB11_PROPERTY(bool UseEllipsoidalGrains READ getUseEllipsoidalGrains WRITE setUseEllipsoidalGrains)
  PYB11_PROPERTY(bool UseCorrectionalMatrix READ getUseCorrectionalMatrix WRITE setUseCorrectionalMatrix)
  PYB11_PROPERTY(float Beta11 READ getBeta11 WRITE setBeta11)
  PYB11_PROPERTY(float Beta22 READ getBeta22 WRITE setBeta22)
  PYB11_PROPERTY(float Beta33 READ getBeta33 WRITE setBeta33)
  PYB11_PROPERTY(float Beta23 READ getBeta23 WRITE setBeta23)
  PYB11_PROPERTY(float Beta13 READ getBeta13 WRITE setBeta13)
  PYB11_PROPERTY(float Beta12 READ getBeta12 WRITE setBeta12)
  PYB11_PROPERTY(DataArrayPath AxisLengthsArrayPath READ getAxisLengthsArrayPath WRITE setAxisLengthsArrayPath)
  PYB11_PROPERTY(DataArrayPath AxisEulerAnglesArrayPath READ getAxisEulerAnglesArrayPath WRITE setAxisEulerAnglesArrayPath)
  PYB11_PROPERTY(DataArrayPath ElasticStrainsArrayPath READ getElasticStrainsArrayPath WRITE setElasticStrainsArrayPath)
  PYB11_PROPERTY(QString EigenstrainsArrayName READ getEigenstrainsArrayName WRITE setEigenstrainsArrayName)
  
  PYB11_END_BINDINGS()
  // clang-format on

public:
  using Self = ComputeFeatureEigenstrains;
  using Pointer = std::shared_ptr<Self>;
  using ConstPointer = std::shared_ptr<const Self>;
  using WeakPointer = std::weak_ptr<Self>;
  using ConstWeakPointer = std::weak_ptr<const Self>;
  static Pointer NullPointer();

  static Pointer New();

  /**
   * @brief Returns the name of the class for ComputeFeatureEigenstrains
   */
  QString getNameOfClass() const override;

  /**
   * @brief Returns the name of the class for ComputeFeatureEigenstrains
   */
  static QString ClassName();

  ~ComputeFeatureEigenstrains() override;

  /**
   * @brief Setter property for PoissonRatio
   */
  void setPoissonRatio(float value);
  /**
   * @brief Getter property for PoissonRatio
   * @return Value of PoissonRatio
   */
  float getPoissonRatio() const;
  Q_PROPERTY(float PoissonRatio READ getPoissonRatio WRITE setPoissonRatio)

  /**
   * @brief Setter property for UseEllipsoidalGrains
   */
  void setUseEllipsoidalGrains(bool value);
  /**
   * @brief Getter property for UseEllipsoidalGrains
   * @return Value of UseEllipsoidalGrains
   */
  bool getUseEllipsoidalGrains() const;
  Q_PROPERTY(bool UseEllipsoidalGrains READ getUseEllipsoidalGrains WRITE setUseEllipsoidalGrains)

  /**
   * @brief Setter property for UseCorrectionalMatrix
   */
  void setUseCorrectionalMatrix(bool value);
  /**
   * @brief Getter property for UseCorrectionalMatrix
   * @return Value of UseCorrectionalMatrix
   */
  bool getUseCorrectionalMatrix() const;
  Q_PROPERTY(bool UseCorrectionalMatrix READ getUseCorrectionalMatrix WRITE setUseCorrectionalMatrix)

  /**
   * @brief Setter property for Beta11
   */
  void setBeta11(float value);
  /**
   * @brief Getter property for Beta11
   * @return Value of Beta11
   */
  float getBeta11() const;
  Q_PROPERTY(float Beta11 READ getBeta11 WRITE setBeta11)

  /**
   * @brief Setter property for Beta22
   */
  void setBeta22(float value);
  /**
   * @brief Getter property for Beta22
   * @return Value of Beta22
   */
  float getBeta22() const;
  Q_PROPERTY(float Beta22 READ getBeta22 WRITE setBeta22)

  /**
   * @brief Setter property for Beta33
   */
  void setBeta33(float value);
  /**
   * @brief Getter property for Beta33
   * @return Value of Beta33
   */
  float getBeta33() const;
  Q_PROPERTY(float Beta33 READ getBeta33 WRITE setBeta33)

  /**
   * @brief Setter property for Beta23
   */
  void setBeta23(float value);
  /**
   * @brief Getter property for Beta23
   * @return Value of Beta23
   */
  float getBeta23() const;
  Q_PROPERTY(float Beta23 READ getBeta23 WRITE setBeta23)

  /**
   * @brief Setter property for Beta13
   */
  void setBeta13(float value);
  /**
   * @brief Getter property for Beta13
   * @return Value of Beta13
   */
  float getBeta13() const;
  Q_PROPERTY(float Beta13 READ getBeta13 WRITE setBeta13)

  /**
   * @brief Setter property for Beta12
   */
  void setBeta12(float value);
  /**
   * @brief Getter property for Beta12
   * @return Value of Beta12
   */
  float getBeta12() const;
  Q_PROPERTY(float Beta12 READ getBeta12 WRITE setBeta12)

  /**
   * @brief Setter property for AxisLengthsArrayPath
   */
  void setAxisLengthsArrayPath(const DataArrayPath& value);
  /**
   * @brief Getter property for AxisLengthsArrayPath
   * @return Value of AxisLengthsArrayPath
   */
  DataArrayPath getAxisLengthsArrayPath() const;
  Q_PROPERTY(DataArrayPath AxisLengthsArrayPath READ getAxisLengthsArrayPath WRITE setAxisLengthsArrayPath)

  /**
   * @brief Setter property for AxisEulerAnglesArrayPath
   */
  void setAxisEulerAnglesArrayPath(const DataArrayPath& value);
  /**
   * @brief Getter property for AxisEulerAnglesArrayPath
   * @return Value of AxisEulerAnglesArrayPath
   */
  DataArrayPath getAxisEulerAnglesArrayPath() const;
  Q_PROPERTY(DataArrayPath AxisEulerAnglesArrayPath READ getAxisEulerAnglesArrayPath WRITE setAxisEulerAnglesArrayPath)

  /**
   * @brief Setter property for ElasticStrainsArrayPath
   */
  void setElasticStrainsArrayPath(const DataArrayPath& value);
  /**
   * @brief Getter property for ElasticStrainsArrayPath
   * @return Value of ElasticStrainsArrayPath
   */
  DataArrayPath getElasticStrainsArrayPath() const;
  Q_PROPERTY(DataArrayPath ElasticStrainsArrayPath READ getElasticStrainsArrayPath WRITE setElasticStrainsArrayPath)

  /**
   * @brief Setter property for EigenstrainsArrayName
   */
  void setEigenstrainsArrayName(const QString& value);
  /**
   * @brief Getter property for EigenstrainsArrayName
   * @return Value of EigenstrainsArrayName
   */
  QString getEigenstrainsArrayName() const;
  Q_PROPERTY(QString EigenstrainsArrayName READ getEigenstrainsArrayName WRITE setEigenstrainsArrayName)

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
   * @brief execute Reimplemented from @see AbstractFilter class
   */
  void execute() override;

  /**
   * @brief preflight Reimplemented from @see AbstractFilter class
   */
  void preflight() override;

protected:
  ComputeFeatureEigenstrains();

  /**
   * @brief dataCheck Checks for the appropriate parameter values and availability of arrays
   */
  void dataCheck() override;

  /**
   * @brief Initializes all the private instance variables.
   */
  void initialize();

  /**
   * @brief find_eigenstrains Calculates the eigenstrains for each feature
   */
  void find_eigenstrains();

private:
  std::weak_ptr<FloatArrayType> m_AxisLengthsPtr;
  std::weak_ptr<FloatArrayType> m_AxisEulerAnglesPtr;
  std::weak_ptr<FloatArrayType> m_ElasticStrainsPtr;
  std::weak_ptr<FloatArrayType> m_EigenstrainsPtr;

  float m_PoissonRatio = {0.33f};
  bool m_UseEllipsoidalGrains = {true};
  bool m_UseCorrectionalMatrix = {false};
  float m_Beta11 = {1.0f};
  float m_Beta22 = {1.0f};
  float m_Beta33 = {1.0f};
  float m_Beta23 = {1.0f};
  float m_Beta13 = {1.0f};
  float m_Beta12 = {1.0f};
  DataArrayPath m_FeatureIdsArrayPath = {"", "", ""};
  DataArrayPath m_AxisLengthsArrayPath = {SIMPL::Defaults::ImageDataContainerName, SIMPL::Defaults::CellFeatureAttributeMatrixName, SIMPL::FeatureData::AxisLengths};
  DataArrayPath m_AxisEulerAnglesArrayPath = {SIMPL::Defaults::ImageDataContainerName, SIMPL::Defaults::CellFeatureAttributeMatrixName, SIMPL::FeatureData::AxisEulerAngles};
  DataArrayPath m_ElasticStrainsArrayPath = {SIMPL::Defaults::ImageDataContainerName, SIMPL::Defaults::CellFeatureAttributeMatrixName, "ElasticStrains"};
  QString m_EigenstrainsArrayName = {"Eigenstrains"};

public:
  ComputeFeatureEigenstrains(const ComputeFeatureEigenstrains&) = delete;            // Copy Constructor Not Implemented
  ComputeFeatureEigenstrains& operator=(const ComputeFeatureEigenstrains&) = delete; // Copy Assignment Not Implemented
  ComputeFeatureEigenstrains(ComputeFeatureEigenstrains&&) = delete;                 // Move Constructor Not Implemented
  ComputeFeatureEigenstrains& operator=(ComputeFeatureEigenstrains&&) = delete;      // Move Assignment Not Implemented
};
