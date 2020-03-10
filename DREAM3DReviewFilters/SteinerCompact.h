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

#include "DREAM3DReview/DREAM3DReviewDLLExport.h"

/**
 * @brief The SteinerCompact class. See [Filter documentation](@ref steinercompact) for details.
 */
class DREAM3DReview_EXPORT SteinerCompact : public AbstractFilter
{
  Q_OBJECT

#ifdef SIMPL_ENABLE_PYTHON
  PYB11_CREATE_BINDINGS(SteinerCompact SUPERCLASS AbstractFilter)
  PYB11_SHARED_POINTERS(SteinerCompact)
  PYB11_FILTER_NEW_MACRO(SteinerCompact)
  PYB11_FILTER_PARAMETER(bool, VtkOutput)
  PYB11_FILTER_PARAMETER(QString, VtkFileName)
  PYB11_FILTER_PARAMETER(bool, TxtOutput)
  PYB11_FILTER_PARAMETER(QString, TxtFileName)
  PYB11_FILTER_PARAMETER(DataArrayPath, FeatureIdsArrayPath)
  PYB11_FILTER_PARAMETER(DataArrayPath, CellPhasesArrayPath)
  PYB11_FILTER_PARAMETER(int, Plane)
  PYB11_FILTER_PARAMETER(int, Sites)
  PYB11_PROPERTY(bool VtkOutput READ getVtkOutput WRITE setVtkOutput)
  PYB11_PROPERTY(QString VtkFileName READ getVtkFileName WRITE setVtkFileName)
  PYB11_PROPERTY(bool TxtOutput READ getTxtOutput WRITE setTxtOutput)
  PYB11_PROPERTY(QString TxtFileName READ getTxtFileName WRITE setTxtFileName)
  PYB11_PROPERTY(DataArrayPath FeatureIdsArrayPath READ getFeatureIdsArrayPath WRITE setFeatureIdsArrayPath)
  PYB11_PROPERTY(DataArrayPath CellPhasesArrayPath READ getCellPhasesArrayPath WRITE setCellPhasesArrayPath)
  PYB11_PROPERTY(int Plane READ getPlane WRITE setPlane)
  PYB11_PROPERTY(int Sites READ getSites WRITE setSites)
#endif

public:
  using Self = SteinerCompact;
  using Pointer = std::shared_ptr<Self>;
  using ConstPointer = std::shared_ptr<const Self>;
  using WeakPointer = std::weak_ptr<Self>;
  using ConstWeakPointer = std::weak_ptr<const Self>;
  static Pointer NullPointer();

  static std::shared_ptr<SteinerCompact> New();

  /**
   * @brief Returns the name of the class for SteinerCompact
   */
  QString getNameOfClass() const override;
  /**
   * @brief Returns the name of the class for SteinerCompact
   */
  static QString ClassName();

  ~SteinerCompact() override;

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
   * @brief Setter property for VtkOutput
   */
  void setVtkOutput(bool value);
  /**
   * @brief Getter property for VtkOutput
   * @return Value of VtkOutput
   */
  bool getVtkOutput() const;

  Q_PROPERTY(bool VtkOutput READ getVtkOutput WRITE setVtkOutput)

  /**
   * @brief Setter property for VtkFileName
   */
  void setVtkFileName(const QString& value);
  /**
   * @brief Getter property for VtkFileName
   * @return Value of VtkFileName
   */
  QString getVtkFileName() const;

  Q_PROPERTY(QString VtkFileName READ getVtkFileName WRITE setVtkFileName)

  /**
   * @brief Setter property for TxtOutput
   */
  void setTxtOutput(bool value);
  /**
   * @brief Getter property for TxtOutput
   * @return Value of TxtOutput
   */
  bool getTxtOutput() const;

  Q_PROPERTY(bool TxtOutput READ getTxtOutput WRITE setTxtOutput)

  /**
   * @brief Setter property for TxtFileName
   */
  void setTxtFileName(const QString& value);
  /**
   * @brief Getter property for TxtFileName
   * @return Value of TxtFileName
   */
  QString getTxtFileName() const;

  Q_PROPERTY(QString TxtFileName READ getTxtFileName WRITE setTxtFileName)

  /**
   * @brief Setter property for FeatureIdsArrayPath
   */
  void setFeatureIdsArrayPath(const DataArrayPath& value);
  /**
   * @brief Getter property for FeatureIdsArrayPath
   * @return Value of FeatureIdsArrayPath
   */
  DataArrayPath getFeatureIdsArrayPath() const;

  Q_PROPERTY(DataArrayPath FeatureIdsArrayPath READ getFeatureIdsArrayPath WRITE setFeatureIdsArrayPath)

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
   * @brief Setter property for Plane
   */
  void setPlane(int value);
  /**
   * @brief Getter property for Plane
   * @return Value of Plane
   */
  int getPlane() const;

  Q_PROPERTY(int Plane READ getPlane WRITE setPlane)

  /**
   * @brief Setter property for Sites
   */
  void setSites(int value);
  /**
   * @brief Getter property for Sites
   * @return Value of Sites
   */
  int getSites() const;

  Q_PROPERTY(int Sites READ getSites WRITE setSites)

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
  SteinerCompact();
  /**
   * @brief dataCheck Checks for the appropriate parameter values and availability of arrays
   */
  void dataCheck() override;

  /**
   * @brief Initializes all the private instance variables.
   */
  void initialize();

  /**
   * @brief  rose_of_intersections Counts the numbers of intersections of randomly placed lines at fixed directions with grain boundaries
   */
  virtual void rose_of_intersections(std::vector<std::vector<float>>& ROI);

  /**
   * @brief find_intersections Counts the number of intersections of given line with grain boundaries
   */
  virtual void find_intersections(float line[3], int64_t z, int64_t dims[3], float res[3], std::vector<float>& numofintersections);

  /**
   * @brief line_rectangle_intersections Finds the intersections of the line (line[0] * x + line[1] * y + line[2] = 0) with the rectangle (x1, y1, x2, y2)
   */
  virtual uint64_t line_rectangle_intersections(float xintersections[2], float yintersections[2], float rectangle[4], float line[3]);

  /**
   * @brief find_one_site_vertices Identifies vertices of a site of the Steiner compact
   */
  virtual void find_one_site_vertices(std::vector<float> ROI, int64_t index, float vertices[4], float& length);

  /**
   * @brief find_all_vertices Identifies vertices of all sites of the Steiner compact
   */
  virtual void find_all_vertices(std::vector<std::vector<float>>& vertices_x, std::vector<std::vector<float>>& vertices_y, std::vector<std::vector<float>>& radii,
                                 std::vector<std::vector<float>>& ROI);

  /**
   * @brief output_vtk Writes the Steiner compact to a vtk file
   */
  virtual void output_vtk(std::vector<std::vector<float>>& vertices_x, std::vector<std::vector<float>>& vertices_y, std::vector<std::vector<float>>& radii, std::vector<std::vector<float>>& ROI);

  /**
   * @brief output_txt Writes the Steiner compact to a text file
   */
  virtual void output_txt(std::vector<std::vector<float>>& vertices_x, std::vector<std::vector<float>>& vertices_y, std::vector<std::vector<float>>& ROI);

private:
  std::weak_ptr<DataArray<int32_t>> m_FeatureIdsPtr;
  int32_t* m_FeatureIds = nullptr;
  std::weak_ptr<DataArray<int32_t>> m_CellPhasesPtr;
  int32_t* m_CellPhases = nullptr;

  QString m_DataContainerName = {};
  QString m_CellAttributeMatrixName = {};
  bool m_VtkOutput = {};
  QString m_VtkFileName = {};
  bool m_TxtOutput = {};
  QString m_TxtFileName = {};
  DataArrayPath m_FeatureIdsArrayPath = {};
  DataArrayPath m_CellPhasesArrayPath = {};
  int m_Plane = {};
  int m_Sites = {};

public:
  SteinerCompact(const SteinerCompact&) = delete;            // Copy Constructor Not Implemented
  SteinerCompact(SteinerCompact&&) = delete;                 // Move Constructor Not Implemented
  SteinerCompact& operator=(const SteinerCompact&) = delete; // Copy Assignment Not Implemented
  SteinerCompact& operator=(SteinerCompact&&) = delete;      // Move Assignment Not Implemented
};
