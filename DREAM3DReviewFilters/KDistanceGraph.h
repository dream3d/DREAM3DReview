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
#include "SIMPLib/DataArrays/DataArray.hpp"
#include "SIMPLib/DataArrays/IDataArray.h"
#include "SIMPLib/Filtering/AbstractFilter.h"

#include "DREAM3DReview/DREAM3DReviewDLLExport.h"

/**
 * @brief The KDistanceGraph class. See [Filter documentation](@ref kdistancegraph) for details.
 */
class DREAM3DReview_EXPORT KDistanceGraph : public AbstractFilter
{
  Q_OBJECT

  // Start Python bindings declarations
  PYB11_BEGIN_BINDINGS(KDistanceGraph SUPERCLASS AbstractFilter)
  PYB11_FILTER()
  PYB11_SHARED_POINTERS(KDistanceGraph)
  PYB11_FILTER_NEW_MACRO(KDistanceGraph)
  PYB11_PROPERTY(DataArrayPath SelectedArrayPath READ getSelectedArrayPath WRITE setSelectedArrayPath)
  PYB11_PROPERTY(bool UseMask READ getUseMask WRITE setUseMask)
  PYB11_PROPERTY(DataArrayPath MaskArrayPath READ getMaskArrayPath WRITE setMaskArrayPath)
  PYB11_PROPERTY(DataArrayPath KDistanceArrayPath READ getKDistanceArrayPath WRITE setKDistanceArrayPath)
  PYB11_PROPERTY(int MinDist READ getMinDist WRITE setMinDist)
  PYB11_PROPERTY(int DistanceMetric READ getDistanceMetric WRITE setDistanceMetric)
  PYB11_END_BINDINGS()
  // End Python bindings declarations

public:
  using Self = KDistanceGraph;
  using Pointer = std::shared_ptr<Self>;
  using ConstPointer = std::shared_ptr<const Self>;
  using WeakPointer = std::weak_ptr<Self>;
  using ConstWeakPointer = std::weak_ptr<const Self>;
  static Pointer NullPointer();

  static std::shared_ptr<KDistanceGraph> New();

  /**
   * @brief Returns the name of the class for KDistanceGraph
   */
  QString getNameOfClass() const override;
  /**
   * @brief Returns the name of the class for KDistanceGraph
   */
  static QString ClassName();

  ~KDistanceGraph() override;

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
   * @brief Setter property for KDistanceArrayPath
   */
  void setKDistanceArrayPath(const DataArrayPath& value);
  /**
   * @brief Getter property for KDistanceArrayPath
   * @return Value of KDistanceArrayPath
   */
  DataArrayPath getKDistanceArrayPath() const;
  Q_PROPERTY(DataArrayPath KDistanceArrayPath READ getKDistanceArrayPath WRITE setKDistanceArrayPath)

  /**
   * @brief Setter property for MinDist
   */
  void setMinDist(int value);
  /**
   * @brief Getter property for MinDist
   * @return Value of MinDist
   */
  int getMinDist() const;
  Q_PROPERTY(int MinDist READ getMinDist WRITE setMinDist)

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
  KDistanceGraph();

  /**
   * @brief dataCheck Checks for the appropriate parameter values and availability of arrays
   */
  void dataCheck() override;

  /**
   * @brief Initializes all the private instance variables.
   */
  void initialize();

private:
  std::weak_ptr<DataArray<bool>> m_MaskPtr;
  bool* m_Mask = nullptr;
  std::weak_ptr<DataArray<double>> m_KDistanceArrayPtr;
  double* m_KDistanceArray = nullptr;

  DataArrayPath m_SelectedArrayPath = {"", "", ""};
  bool m_UseMask = {false};
  DataArrayPath m_MaskArrayPath = {"", "", ""};
  DataArrayPath m_KDistanceArrayPath = {"", "", "KDistance"};
  int m_MinDist = {1};
  int m_DistanceMetric = {0};

  IDataArray::WeakPointer m_InDataPtr;

public:
  KDistanceGraph(const KDistanceGraph&) = delete;            // Copy Constructor Not Implemented
  KDistanceGraph(KDistanceGraph&&) = delete;                 // Move Constructor Not Implemented
  KDistanceGraph& operator=(const KDistanceGraph&) = delete; // Copy Assignment Not Implemented
  KDistanceGraph& operator=(KDistanceGraph&&) = delete;      // Move Assignment Not Implemented
};
