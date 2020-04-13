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
#include "SIMPLib/DataArrays/DataArray.hpp"

#include "DREAM3DReview/DREAM3DReviewDLLExport.h"

/**
 * @brief The NormalizeArrays class. See [Filter documentation](@ref normalizearrays) for details.
 */
class DREAM3DReview_EXPORT NormalizeArrays : public AbstractFilter
{
  Q_OBJECT

  // Start Python bindings declarations
  PYB11_BEGIN_BINDINGS(NormalizeArrays SUPERCLASS AbstractFilter)
  PYB11_FILTER()
  PYB11_SHARED_POINTERS(NormalizeArrays)
  PYB11_FILTER_NEW_MACRO(NormalizeArrays)
  PYB11_PROPERTY(QVector<DataArrayPath> SelectedDataArrayPaths READ getSelectedDataArrayPaths WRITE setSelectedDataArrayPaths)
  PYB11_PROPERTY(int NormalizeType READ getNormalizeType WRITE setNormalizeType)
  PYB11_PROPERTY(double RangeMin READ getRangeMin WRITE setRangeMin)
  PYB11_PROPERTY(double RangeMax READ getRangeMax WRITE setRangeMax)
  PYB11_PROPERTY(QString Postfix READ getPostfix WRITE setPostfix)
  PYB11_PROPERTY(bool UseMask READ getUseMask WRITE setUseMask)
  PYB11_PROPERTY(DataArrayPath MaskArrayPath READ getMaskArrayPath WRITE setMaskArrayPath)
  PYB11_PROPERTY(double DefaultValue READ getDefaultValue WRITE setDefaultValue)
  PYB11_END_BINDINGS()
  // End Python bindings declarations

public:
  using Self = NormalizeArrays;
  using Pointer = std::shared_ptr<Self>;
  using ConstPointer = std::shared_ptr<const Self>;
  using WeakPointer = std::weak_ptr<Self>;
  using ConstWeakPointer = std::weak_ptr<const Self>;
  static Pointer NullPointer();

  static std::shared_ptr<NormalizeArrays> New();

  /**
   * @brief Returns the name of the class for NormalizeArrays
   */
  QString getNameOfClass() const override;
  /**
   * @brief Returns the name of the class for NormalizeArrays
   */
  static QString ClassName();

  ~NormalizeArrays() override;

  /**
   * @brief Setter property for SelectedDataArrayPaths
   */
  void setSelectedDataArrayPaths(const QVector<DataArrayPath>& value);
  /**
   * @brief Getter property for SelectedDataArrayPaths
   * @return Value of SelectedDataArrayPaths
   */
  QVector<DataArrayPath> getSelectedDataArrayPaths() const;
  Q_PROPERTY(QVector<DataArrayPath> SelectedDataArrayPaths READ getSelectedDataArrayPaths WRITE setSelectedDataArrayPaths)

  /**
   * @brief Setter property for NormalizeType
   */
  void setNormalizeType(int value);
  /**
   * @brief Getter property for NormalizeType
   * @return Value of NormalizeType
   */
  int getNormalizeType() const;
  Q_PROPERTY(int NormalizeType READ getNormalizeType WRITE setNormalizeType)

  /**
   * @brief Setter property for RangeMin
   */
  void setRangeMin(double value);
  /**
   * @brief Getter property for RangeMin
   * @return Value of RangeMin
   */
  double getRangeMin() const;
  Q_PROPERTY(double RangeMin READ getRangeMin WRITE setRangeMin)

  /**
   * @brief Setter property for RangeMax
   */
  void setRangeMax(double value);
  /**
   * @brief Getter property for RangeMax
   * @return Value of RangeMax
   */
  double getRangeMax() const;
  Q_PROPERTY(double RangeMax READ getRangeMax WRITE setRangeMax)

  /**
   * @brief Setter property for Postfix
   */
  void setPostfix(const QString& value);
  /**
   * @brief Getter property for Postfix
   * @return Value of Postfix
   */
  QString getPostfix() const;
  Q_PROPERTY(QString Postfix READ getPostfix WRITE setPostfix)

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
   * @brief Setter property for DefaultValue
   */
  void setDefaultValue(double value);
  /**
   * @brief Getter property for DefaultValue
   * @return Value of DefaultValue
   */
  double getDefaultValue() const;
  Q_PROPERTY(double DefaultValue READ getDefaultValue WRITE setDefaultValue)

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
  NormalizeArrays();

  void updateDataArrays(double* copyPtr, DoubleArrayType::Pointer newArray, int32_t arrayIndex);

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

  QVector<DataArrayPath> m_SelectedDataArrayPaths = {};
  int m_NormalizeType = {};
  double m_RangeMin = {};
  double m_RangeMax = {};
  QString m_Postfix = {};
  bool m_UseMask = {};
  DataArrayPath m_MaskArrayPath = {};
  double m_DefaultValue = {};

  std::vector<IDataArray::WeakPointer> m_SelectedWeakPtrVector;
  std::vector<DoubleArrayType::Pointer> m_NormalizedArraysPtrVector;

public:
  NormalizeArrays(const NormalizeArrays&) = delete; // Copy Constructor Not Implemented
  NormalizeArrays(NormalizeArrays&&) = delete;      // Move Constructor Not Implemented
  NormalizeArrays& operator=(const NormalizeArrays&) = delete; // Copy Assignment Not Implemented
  NormalizeArrays& operator=(NormalizeArrays&&) = delete;      // Move Assignment Not Implemented
};

