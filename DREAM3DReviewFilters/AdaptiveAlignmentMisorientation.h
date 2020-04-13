/* ============================================================================
 * Copyright (c) 2016 Czech Academy of Sciences, Institute of Physics,
 * Group of Bulk Nanomaterials and Interfaces, http://ams.fzu.cz
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
 * Neither the name of the Czech Academy of Sciences, nor the names of its
 * contributors may be used to endorse or promote products derived from this
 * software without specific prior written permission.
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
 * The code contained herein was partially funded by the followig grants:
 *    Czech Science Foundation (GA CR), project no. GBP108/12/G043
 *    Czech Ministry of Education, Youth and Sports (MSMT), project no. LM2015087
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#pragma once

#include <memory>

#include "SIMPLib/SIMPLib.h"
#include "SIMPLib/Filtering/AbstractFilter.h"
#include "SIMPLib/DataArrays/DataArray.hpp"

#include "AdaptiveAlignment.h"

#include "DREAM3DReview/DREAM3DReviewDLLExport.h"

class LaueOps;
using LaueOpsShPtrType = std::shared_ptr<LaueOps>;
using LaueOpsContainer = std::vector<LaueOpsShPtrType>;
/**
 * @brief The AdaptiveAlignmentMisorientation class. See [Filter documentation](@ref AdaptiveAlignmentmisorientation) for details.
 */
class DREAM3DReview_EXPORT AdaptiveAlignmentMisorientation : public AdaptiveAlignment
{
  Q_OBJECT

  // Start Python bindings declarations
  PYB11_BEGIN_BINDINGS(AdaptiveAlignmentMisorientation SUPERCLASS AbstractFilter)
  PYB11_FILTER()
  PYB11_SHARED_POINTERS(AdaptiveAlignmentMisorientation)
  PYB11_FILTER_NEW_MACRO(AdaptiveAlignmentMisorientation)
  PYB11_PROPERTY(float MisorientationTolerance READ getMisorientationTolerance WRITE setMisorientationTolerance)
  PYB11_PROPERTY(bool UseGoodVoxels READ getUseGoodVoxels WRITE setUseGoodVoxels)
  PYB11_PROPERTY(DataArrayPath QuatsArrayPath READ getQuatsArrayPath WRITE setQuatsArrayPath)
  PYB11_PROPERTY(DataArrayPath CellPhasesArrayPath READ getCellPhasesArrayPath WRITE setCellPhasesArrayPath)
  PYB11_PROPERTY(DataArrayPath GoodVoxelsArrayPath READ getGoodVoxelsArrayPath WRITE setGoodVoxelsArrayPath)
  PYB11_PROPERTY(DataArrayPath CrystalStructuresArrayPath READ getCrystalStructuresArrayPath WRITE setCrystalStructuresArrayPath)
  PYB11_END_BINDINGS()
  // End Python bindings declarations

public:
  using Self = AdaptiveAlignmentMisorientation;
  using Pointer = std::shared_ptr<Self>;
  using ConstPointer = std::shared_ptr<const Self>;
  using WeakPointer = std::weak_ptr<Self>;
  using ConstWeakPointer = std::weak_ptr<const Self>;
  static Pointer NullPointer();

  static std::shared_ptr<AdaptiveAlignmentMisorientation> New();

  /**
   * @brief Returns the name of the class for AdaptiveAlignmentMisorientation
   */
  QString getNameOfClass() const override;
  /**
   * @brief Returns the name of the class for AdaptiveAlignmentMisorientation
   */
  static QString ClassName();

  virtual ~AdaptiveAlignmentMisorientation();

  /**
   * @brief Setter property for MisorientationTolerance
   */
  void setMisorientationTolerance(float value);
  /**
   * @brief Getter property for MisorientationTolerance
   * @return Value of MisorientationTolerance
   */
  float getMisorientationTolerance() const;
  Q_PROPERTY(float MisorientationTolerance READ getMisorientationTolerance WRITE setMisorientationTolerance)

  /**
   * @brief Setter property for UseGoodVoxels
   */
  void setUseGoodVoxels(bool value);
  /**
   * @brief Getter property for UseGoodVoxels
   * @return Value of UseGoodVoxels
   */
  bool getUseGoodVoxels() const;
  Q_PROPERTY(bool UseGoodVoxels READ getUseGoodVoxels WRITE setUseGoodVoxels)

  /**
   * @brief Setter property for QuatsArrayPath
   */
  void setQuatsArrayPath(const DataArrayPath& value);
  /**
   * @brief Getter property for QuatsArrayPath
   * @return Value of QuatsArrayPath
   */
  DataArrayPath getQuatsArrayPath() const;
  Q_PROPERTY(DataArrayPath QuatsArrayPath READ getQuatsArrayPath WRITE setQuatsArrayPath)

  /**
   * @brief Setter property for CellPhasesArrayPath
   */
  void setCellPhasesArrayPath(const DataArrayPath& value);
  /**
   * @brief Getter property for CellPhasesArrayPath
   * @return Value of CellPhasesArrayPath
   */
  DataArrayPath getCellPhasesArrayPath() const;
  Q_PROPERTY(DataArrayPath CellPhasesArrayPath READ getCellPhasesArrayPath WRITE setCellPhasesArrayPath)

  /**
   * @brief Setter property for GoodVoxelsArrayPath
   */
  void setGoodVoxelsArrayPath(const DataArrayPath& value);
  /**
   * @brief Getter property for GoodVoxelsArrayPath
   * @return Value of GoodVoxelsArrayPath
   */
  DataArrayPath getGoodVoxelsArrayPath() const;
  Q_PROPERTY(DataArrayPath GoodVoxelsArrayPath READ getGoodVoxelsArrayPath WRITE setGoodVoxelsArrayPath)

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
  AdaptiveAlignmentMisorientation();

  /**
   * @brief dataCheck Checks for the appropriate parameter values and availability of arrays
   */
  void dataCheck() override;

  /**
   * @brief Initializes all the private instance variables.
   */
  void initialize();

  /**
   * @brief find_shifts Reimplemented from @see AdaptiveAlignment class
   */
  virtual void find_shifts(std::vector<int64_t>& xshifts, std::vector<int64_t>& yshifts, std::vector<float>& xneedshifts, std::vector<float>& yneedshifts);

  /**
   * @brief compute_error1 Determines error between the current and desired shifts (based on discrepancy of slopes)
   */
  virtual float compute_error1(uint64_t iter, uint64_t index, float xneedtrend, float yneedtrend, std::vector<std::vector<int64_t>>& newxshift, std::vector<std::vector<int64_t>>& newyshift,
                               std::vector<uint64_t>& curindex);

  /**
   * @brief compute_error1 Determines error between the current and desired shifts (based on discrepancy from the estimates obtained from SEM images)
   */
  virtual float compute_error2(uint64_t iter, uint64_t index, std::vector<float>& xshiftsest, std::vector<float>& yshiftsest, std::vector<std::vector<int64_t>>& newxshift,
                               std::vector<std::vector<int64_t>>& newyshift, std::vector<uint64_t>& curindex);

private:
  std::weak_ptr<DataArray<float>> m_QuatsPtr;
  float* m_Quats = nullptr;
  std::weak_ptr<DataArray<int32_t>> m_CellPhasesPtr;
  int32_t* m_CellPhases = nullptr;
  std::weak_ptr<DataArray<bool>> m_GoodVoxelsPtr;
  bool* m_GoodVoxels = nullptr;
  std::weak_ptr<DataArray<uint32_t>> m_CrystalStructuresPtr;
  uint32_t* m_CrystalStructures = nullptr;

  float m_MisorientationTolerance = {};
  bool m_UseGoodVoxels = {};
  DataArrayPath m_QuatsArrayPath = {};
  DataArrayPath m_CellPhasesArrayPath = {};
  DataArrayPath m_GoodVoxelsArrayPath = {};
  DataArrayPath m_CrystalStructuresArrayPath = {};

  DataArray<uint8_t>::WeakPointer m_ImageDataPtr;
  DataArray<AnisotropyConstants::DefaultPixelType>::WeakPointer SelectedCellArrayPtr;
  DataArray<AnisotropyConstants::DefaultPixelType>::WeakPointer NewCellArrayPtr;

  LaueOpsContainer m_OrientationOps;
  uint64_t m_RandomSeed;

public:
  AdaptiveAlignmentMisorientation(const AdaptiveAlignmentMisorientation&) = delete;            // Copy Constructor Not Implemented
  AdaptiveAlignmentMisorientation(AdaptiveAlignmentMisorientation&&) = delete;                 // Move Constructor Not Implemented
  AdaptiveAlignmentMisorientation& operator=(const AdaptiveAlignmentMisorientation&) = delete; // Copy Assignment Not Implemented
  AdaptiveAlignmentMisorientation& operator=(AdaptiveAlignmentMisorientation&&) = delete;      // Move Assignment Not Implemented
};
