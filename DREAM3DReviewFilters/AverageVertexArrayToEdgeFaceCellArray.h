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

#include "SIMPLib/SIMPLib.h"
#include "SIMPLib/DataArrays/DataArray.hpp"
#include "SIMPLib/Filtering/AbstractFilter.h"

class IDataArray;
using IDataArrayWkPtrType = std::weak_ptr<IDataArray>;

#include "DREAM3DReview/DREAM3DReviewDLLExport.h"

/**
 * @brief The AverageVertexArrayToEdgeFaceCellArray class. See [Filter documentation](@ref averagevertexarraytoedgefacecellArray) for details.
 */
class DREAM3DReview_EXPORT AverageVertexArrayToEdgeFaceCellArray : public AbstractFilter
{
  Q_OBJECT

  // Start Python bindings declarations
  PYB11_BEGIN_BINDINGS(AverageVertexArrayToEdgeFaceCellArray SUPERCLASS AbstractFilter)
  PYB11_FILTER()
  PYB11_SHARED_POINTERS(AverageVertexArrayToEdgeFaceCellArray)
  PYB11_FILTER_NEW_MACRO(AverageVertexArrayToEdgeFaceCellArray)
  PYB11_PROPERTY(DataArrayPath SelectedArrayPath READ getSelectedArrayPath WRITE setSelectedArrayPath)
  PYB11_PROPERTY(DataArrayPath AverageCellArrayPath READ getAverageCellArrayPath WRITE setAverageCellArrayPath)
  PYB11_PROPERTY(bool WeightedAverage READ getWeightedAverage WRITE setWeightedAverage)
  PYB11_END_BINDINGS()
  // End Python bindings declarations

public:
  using Self = AverageVertexArrayToEdgeFaceCellArray;
  using Pointer = std::shared_ptr<Self>;
  using ConstPointer = std::shared_ptr<const Self>;
  using WeakPointer = std::weak_ptr<Self>;
  using ConstWeakPointer = std::weak_ptr<const Self>;
  static Pointer NullPointer();

  static std::shared_ptr<AverageVertexArrayToEdgeFaceCellArray> New();

  /**
   * @brief Returns the name of the class for AverageVertexArrayToEdgeFaceCellArray
   */
  QString getNameOfClass() const override;
  /**
   * @brief Returns the name of the class for AverageVertexArrayToEdgeFaceCellArray
   */
  static QString ClassName();

  ~AverageVertexArrayToEdgeFaceCellArray() override;

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
   * @brief Setter property for AverageCellArrayPath
   */
  void setAverageCellArrayPath(const DataArrayPath& value);
  /**
   * @brief Getter property for AverageCellArrayPath
   * @return Value of AverageCellArrayPath
   */
  DataArrayPath getAverageCellArrayPath() const;
  Q_PROPERTY(DataArrayPath AverageCellArrayPath READ getAverageCellArrayPath WRITE setAverageCellArrayPath)

  /**
   * @brief Setter property for WeightedAverage
   */
  void setWeightedAverage(bool value);
  /**
   * @brief Getter property for WeightedAverage
   * @return Value of WeightedAverage
   */
  bool getWeightedAverage() const;
  Q_PROPERTY(bool WeightedAverage READ getWeightedAverage WRITE setWeightedAverage)

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
  AverageVertexArrayToEdgeFaceCellArray();
  /**
   * @brief dataCheck Checks for the appropriate parameter values and availability of arrays
   */
  void dataCheck() override;

  /**
   * @brief Initializes all the private instance variables.
   */
  void initialize();

private:
  IDataArrayWkPtrType m_InVertexArrayPtr;
  void* m_InVertexArray = nullptr;

  std::weak_ptr<DataArray<float>> m_AverageCellArrayPtr;
  float* m_AverageCellArray = nullptr;

  DataArrayPath m_SelectedArrayPath = {"", "", ""};
  DataArrayPath m_AverageCellArrayPath = {"", "", ""};
  bool m_WeightedAverage = {false};

public:
  AverageVertexArrayToEdgeFaceCellArray(const AverageVertexArrayToEdgeFaceCellArray&) = delete;            // Copy Constructor Not Implemented
  AverageVertexArrayToEdgeFaceCellArray(AverageVertexArrayToEdgeFaceCellArray&&) = delete;                 // Move Constructor Not Implemented
  AverageVertexArrayToEdgeFaceCellArray& operator=(const AverageVertexArrayToEdgeFaceCellArray&) = delete; // Copy Assignment Not Implemented
  AverageVertexArrayToEdgeFaceCellArray& operator=(AverageVertexArrayToEdgeFaceCellArray&&) = delete;      // Move Assignment Not Implemented
};
