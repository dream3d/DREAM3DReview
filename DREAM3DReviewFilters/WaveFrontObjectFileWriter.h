/* ============================================================================
 * Copyright (c) 2020 BlueQuartz Software, LLC
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
 * Neither the names of any of the BlueQuartz Software contributors
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
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
#pragma once

#include "SIMPLib/SIMPLib.h"
#include "SIMPLib/DataArrays/DataArray.hpp"
#include "SIMPLib/Filtering/AbstractFilter.h"

#include "DREAM3DReview/DREAM3DReviewPlugin.h"

/**
 * @brief The WaveFrontObjectFileWriter class. See [Filter documentation](@ref WaveFrontObjectFileWriter) for details.
 */
class DREAM3DReview_EXPORT WaveFrontObjectFileWriter : public AbstractFilter
{
  Q_OBJECT

  // Start Python bindings declarations
  PYB11_BEGIN_BINDINGS(WaveFrontObjectFileWriter SUPERCLASS AbstractFilter)
  PYB11_FILTER()
  PYB11_SHARED_POINTERS(WaveFrontObjectFileWriter)
  PYB11_FILTER_NEW_MACRO(WaveFrontObjectFileWriter)
  PYB11_PROPERTY(QString OutputWaveFrontFile READ getOutputWaveFrontFile WRITE setOutputWaveFrontFile)
  PYB11_PROPERTY(DataArrayPath TriangleGeometry READ getTriangleGeometry WRITE setTriangleGeometry)
  PYB11_END_BINDINGS()
  // End Python bindings declarations

public:
  using Self = WaveFrontObjectFileWriter;
  using Pointer = std::shared_ptr<Self>;
  using ConstPointer = std::shared_ptr<const Self>;
  using WeakPointer = std::weak_ptr<Self>;
  using ConstWeakPointer = std::weak_ptr<const Self>;
  static Pointer NullPointer();

  static Pointer New();

  /**
   * @brief Returns the name of the class for WaveFrontObjectFileWriter
   */
  QString getNameOfClass() const override;

  /**
   * @brief Returns the name of the class for WaveFrontObjectFileWriter
   */
  static QString ClassName();

  ~WaveFrontObjectFileWriter() override;

  /**
   * @brief Setter property for OutputWaveFrontFile
   */
  void setOutputWaveFrontFile(const QString& value);
  /**
   * @brief Getter property for OutputWaveFrontFile
   * @return Value of OutputWaveFrontFile
   */
  QString getOutputWaveFrontFile() const;
  Q_PROPERTY(QString OutputWaveFrontFile READ getOutputWaveFrontFile WRITE setOutputWaveFrontFile)

  /**
   * @brief Setter property for SurfaceMeshNodeTypeArrayPath
   */
  void setTriangleGeometry(const DataArrayPath& value);
  /**
   * @brief Getter property for SurfaceMeshNodeTypeArrayPath
   * @return Value of SurfaceMeshNodeTypeArrayPath
   */
  DataArrayPath getTriangleGeometry() const;
  Q_PROPERTY(DataArrayPath TriangleGeometry READ getTriangleGeometry WRITE setTriangleGeometry)

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

protected:
  WaveFrontObjectFileWriter();

  /**
   * @brief dataCheck Checks for the appropriate parameter values and availability of arrays
   */
  void dataCheck() override;

  /**
   * @brief Initializes all the private instance variables.
   */
  void initialize();

private:
  QString m_OutputWaveFrontFile = {};
  DataArrayPath m_TriangleGeometry = {};

public:
  WaveFrontObjectFileWriter(const WaveFrontObjectFileWriter&) = delete;            // Copy Constructor Not Implemented
  WaveFrontObjectFileWriter& operator=(const WaveFrontObjectFileWriter&) = delete; // Copy Assignment Not Implemented
  WaveFrontObjectFileWriter(WaveFrontObjectFileWriter&&) = delete;                 // Move Constructor Not Implemented
  WaveFrontObjectFileWriter& operator=(WaveFrontObjectFileWriter&&) = delete;      // Move Assignment Not Implemented
};
