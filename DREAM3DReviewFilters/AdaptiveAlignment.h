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

#include "DREAM3DReview/DREAM3DReviewConstants.h"

#include "DREAM3DReview/DREAM3DReviewDLLExport.h"

/**
 * @brief The AdaptiveAlignment class. This class serves as a superclass for other classes
 * in the Reconstruction plugin.
 */
class DREAM3DReview_EXPORT AdaptiveAlignment : public AbstractFilter
{
  Q_OBJECT

#ifdef SIMPL_ENABLE_PYTHON
  PYB11_CREATE_BINDINGS(AdaptiveAlignment SUPERCLASS AbstractFilter)
  PYB11_SHARED_POINTERS(AdaptiveAlignment)
  PYB11_FILTER_NEW_MACRO(AdaptiveAlignment)
  PYB11_FILTER_PARAMETER(int, GlobalCorrection)
  PYB11_FILTER_PARAMETER(QString, InputPath)
  PYB11_FILTER_PARAMETER(float, ShiftX)
  PYB11_FILTER_PARAMETER(float, ShiftY)
  PYB11_FILTER_PARAMETER(DataArrayPath, ImageDataArrayPath)
  PYB11_FILTER_PARAMETER(QString, NewCellArrayName)
  PYB11_FILTER_PARAMETER(float, MinRadius)
  PYB11_FILTER_PARAMETER(float, MaxRadius)
  PYB11_FILTER_PARAMETER(int, NumberCircles)
  PYB11_FILTER_PARAMETER(QVector<DataArrayPath>, IgnoredDataArrayPaths)
  PYB11_PROPERTY(int GlobalCorrection READ getGlobalCorrection WRITE setGlobalCorrection)
  PYB11_PROPERTY(QString InputPath READ getInputPath WRITE setInputPath)
  PYB11_PROPERTY(float ShiftX READ getShiftX WRITE setShiftX)
  PYB11_PROPERTY(float ShiftY READ getShiftY WRITE setShiftY)
  PYB11_PROPERTY(DataArrayPath ImageDataArrayPath READ getImageDataArrayPath WRITE setImageDataArrayPath)
  PYB11_PROPERTY(QString NewCellArrayName READ getNewCellArrayName WRITE setNewCellArrayName)
  PYB11_PROPERTY(float MinRadius READ getMinRadius WRITE setMinRadius)
  PYB11_PROPERTY(float MaxRadius READ getMaxRadius WRITE setMaxRadius)
  PYB11_PROPERTY(int NumberCircles READ getNumberCircles WRITE setNumberCircles)
#endif

public:
  using Self = AdaptiveAlignment;
  using Pointer = std::shared_ptr<Self>;
  using ConstPointer = std::shared_ptr<const Self>;
  using WeakPointer = std::weak_ptr<Self>;
  using ConstWeakPointer = std::weak_ptr<const Self>;
  static Pointer NullPointer();

  static std::shared_ptr<AdaptiveAlignment> New();

  /**
   * @brief Returns the name of the class for AdaptiveAlignment
   */
  QString getNameOfClass() const override;
  /**
   * @brief Returns the name of the class for AdaptiveAlignment
   */
  static QString ClassName();

  ~AdaptiveAlignment() override;

  /**
   * @brief Setter property for DataContainerName
   */
  void setDataContainerName(const QString& value);
  /**
   * @brief Getter property for DataContainerName
   * @return Value of DataContainerName
   */
  QString getDataContainerName() const;

  /**
   * @brief Setter property for CellAttributeMatrixName
   */
  void setCellAttributeMatrixName(const QString& value);
  /**
   * @brief Getter property for CellAttributeMatrixName
   * @return Value of CellAttributeMatrixName
   */
  QString getCellAttributeMatrixName() const;

  /**
   * @brief Setter property for WriteAlignmentShifts
   */
  void setWriteAlignmentShifts(bool value);
  /**
   * @brief Getter property for WriteAlignmentShifts
   * @return Value of WriteAlignmentShifts
   */
  bool getWriteAlignmentShifts() const;

  Q_PROPERTY(bool WriteAlignmentShifts READ getWriteAlignmentShifts WRITE setWriteAlignmentShifts)

  /**
   * @brief Setter property for AlignmentShiftFileName
   */
  void setAlignmentShiftFileName(const QString& value);
  /**
   * @brief Getter property for AlignmentShiftFileName
   * @return Value of AlignmentShiftFileName
   */
  QString getAlignmentShiftFileName() const;

  Q_PROPERTY(QString AlignmentShiftFileName READ getAlignmentShiftFileName WRITE setAlignmentShiftFileName)

  /////////////// new:

  /**
   * @brief Setter property for GlobalCorrection
   */
  void setGlobalCorrection(int value);
  /**
   * @brief Getter property for GlobalCorrection
   * @return Value of GlobalCorrection
   */
  int getGlobalCorrection() const;

  Q_PROPERTY(int GlobalCorrection READ getGlobalCorrection WRITE setGlobalCorrection)

  /**
   * @brief Setter property for InputPath
   */
  void setInputPath(const QString& value);
  /**
   * @brief Getter property for InputPath
   * @return Value of InputPath
   */
  QString getInputPath() const;

  Q_PROPERTY(QString InputPath READ getInputPath WRITE setInputPath)

  /**
   * @brief Setter property for ShiftX
   */
  void setShiftX(float value);
  /**
   * @brief Getter property for ShiftX
   * @return Value of ShiftX
   */
  float getShiftX() const;

  Q_PROPERTY(float ShiftX READ getShiftX WRITE setShiftX)

  /**
   * @brief Setter property for ShiftY
   */
  void setShiftY(float value);
  /**
   * @brief Getter property for ShiftY
   * @return Value of ShiftY
   */
  float getShiftY() const;

  Q_PROPERTY(float ShiftY READ getShiftY WRITE setShiftY)

  /**
   * @brief Setter property for ImageDataArrayPath
   */
  void setImageDataArrayPath(const DataArrayPath& value);
  /**
   * @brief Getter property for ImageDataArrayPath
   * @return Value of ImageDataArrayPath
   */
  DataArrayPath getImageDataArrayPath() const;

  Q_PROPERTY(DataArrayPath ImageDataArrayPath READ getImageDataArrayPath WRITE setImageDataArrayPath)

  /**
   * @brief Setter property for NewCellArrayName
   */
  void setNewCellArrayName(const QString& value);
  /**
   * @brief Getter property for NewCellArrayName
   * @return Value of NewCellArrayName
   */
  QString getNewCellArrayName() const;

  Q_PROPERTY(QString NewCellArrayName READ getNewCellArrayName WRITE setNewCellArrayName)

  /**
   * @brief Setter property for MinRadius
   */
  void setMinRadius(float value);
  /**
   * @brief Getter property for MinRadius
   * @return Value of MinRadius
   */
  float getMinRadius() const;

  Q_PROPERTY(float MinRadius READ getMinRadius WRITE setMinRadius)

  /**
   * @brief Setter property for MaxRadius
   */
  void setMaxRadius(float value);
  /**
   * @brief Getter property for MaxRadius
   * @return Value of MaxRadius
   */
  float getMaxRadius() const;

  Q_PROPERTY(float MaxRadius READ getMaxRadius WRITE setMaxRadius)

  /**
   * @brief Setter property for NumberCircles
   */
  void setNumberCircles(int value);
  /**
   * @brief Getter property for NumberCircles
   * @return Value of NumberCircles
   */
  int getNumberCircles() const;

  Q_PROPERTY(int NumberCircles READ getNumberCircles WRITE setNumberCircles)

  /**
   * @brief Setter property for IgnoredDataArrayPaths
   */
  void setIgnoredDataArrayPaths(const QVector<DataArrayPath>& value);
  /**
   * @brief Getter property for IgnoredDataArrayPaths
   * @return Value of IgnoredDataArrayPaths
   */
  QVector<DataArrayPath> getIgnoredDataArrayPaths() const;

  Q_PROPERTY(QVector<DataArrayPath> IgnoredDataArrayPaths READ getIgnoredDataArrayPaths WRITE setIgnoredDataArrayPaths)

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
  AdaptiveAlignment();

  /**
   * @brief dataCheck Checks for the appropriate parameter values and availability of arrays
   */
  void dataCheck() override;

  /**
   * @brief Initializes all the private instance variables.
   */
  void initialize();

  /**
   * @brief create_array_for_flattened_image creates a temporary array
   */
  virtual void create_array_for_flattened_image();

  /**
   * @brief delete_array_for_flattened_image deletes a temporary array
   */
  virtual void delete_array_for_flattened_image();

  /**
   * @brief flatten_image for converting an RGB image to grayscale
   */
  virtual void flatten_image();

  /**
   * @brief find_circles for identification of two circles for alignment in SEM images by Hough transform
   */
  virtual bool find_calibrating_circles();

  /**
   * @brief find_rectangles for identification of the green rectangle
   */
  virtual bool find_rectangles();

  /**
   * @brief find_rectangles for identification of the edge of the specimen
   */
  virtual bool find_interface_edges();

  /**
   * @brief estimate_shifts for estimation of the shifts from SEM images
   */
  virtual void estimate_shifts_from_images(std::vector<float>& xneedshifts, std::vector<float>& yneedshifts);

  /**
   * @brief find_shifts Determines the x and y shifts to register a stacked 3D volume
   * @param xshifts Vector of integer shifts in x direction
   * @param yshifts Vector of integer shifts in y direction
   */
  virtual void find_shifts(std::vector<int64_t>& xshifts, std::vector<int64_t>& yshifts, std::vector<float>& xneedshifts, std::vector<float>& yneedshifts);

private:
  std::weak_ptr<DataArray<AnisotropyConstants::DefaultPixelType>> m_ImageDataPtr;
  AnisotropyConstants::DefaultPixelType* m_ImageData = nullptr;
  std::weak_ptr<DataArray<AnisotropyConstants::DefaultPixelType>> m_FlatImageDataPtr;
  AnisotropyConstants::DefaultPixelType* m_FlatImageData = nullptr;

  QString m_DataContainerName = {};
  QString m_CellAttributeMatrixName = {};
  bool m_WriteAlignmentShifts = {};
  QString m_AlignmentShiftFileName = {};
  int m_GlobalCorrection = {};
  QString m_InputPath = {};
  float m_ShiftX = {};
  float m_ShiftY = {};
  DataArrayPath m_ImageDataArrayPath = {};
  QString m_NewCellArrayName = {};
  float m_MinRadius = {};
  float m_MaxRadius = {};
  int m_NumberCircles = {};
  QVector<DataArrayPath> m_IgnoredDataArrayPaths = {};

  DataArray<AnisotropyConstants::DefaultPixelType>::WeakPointer m_NewCellArrayPtr;

  std::vector<std::vector<uint64_t>> m_RectangleCorners;
  std::vector<std::vector<float>> m_CalibratingCircles;
  std::vector<uint64_t> m_InterfaceEdges;

public:
  AdaptiveAlignment(const AdaptiveAlignment&) = delete;            // Copy Constructor Not Implemented
  AdaptiveAlignment(AdaptiveAlignment&&) = delete;                 // Move Constructor Not Implemented
  AdaptiveAlignment& operator=(const AdaptiveAlignment&) = delete; // Copy Assignment Not Implemented
  AdaptiveAlignment& operator=(AdaptiveAlignment&&) = delete;      // Move Assignment Not Implemented
};
