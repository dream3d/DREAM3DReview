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
 * @brief The PrincipalComponentAnalysis class. See [Filter documentation](@ref principalcomponentanalysis) for details.
 */
class DREAM3DReview_EXPORT PrincipalComponentAnalysis : public AbstractFilter
{
  Q_OBJECT

#ifdef SIMPL_ENABLE_PYTHON
  PYB11_CREATE_BINDINGS(PrincipalComponentAnalysis SUPERCLASS AbstractFilter)
  PYB11_SHARED_POINTERS(PrincipalComponentAnalysis)
  PYB11_FILTER_NEW_MACRO(PrincipalComponentAnalysis)
  PYB11_FILTER_PARAMETER(QVector<DataArrayPath>, SelectedDataArrayPaths)
  PYB11_FILTER_PARAMETER(QString, PCAttributeMatrixName)
  PYB11_FILTER_PARAMETER(QString, PCEigenvaluesName)
  PYB11_FILTER_PARAMETER(QString, PCEigenvectorsName)
  PYB11_FILTER_PARAMETER(int, MatrixApproach)
  PYB11_FILTER_PARAMETER(bool, ProjectDataSpace)
  PYB11_FILTER_PARAMETER(int, NumberOfDimensionsForProjection)
  PYB11_FILTER_PARAMETER(DataArrayPath, ProjectedDataSpaceArrayPath)
  PYB11_PROPERTY(QVector<DataArrayPath> SelectedDataArrayPaths READ getSelectedDataArrayPaths WRITE setSelectedDataArrayPaths)
  PYB11_PROPERTY(QString PCAttributeMatrixName READ getPCAttributeMatrixName WRITE setPCAttributeMatrixName)
  PYB11_PROPERTY(QString PCEigenvaluesName READ getPCEigenvaluesName WRITE setPCEigenvaluesName)
  PYB11_PROPERTY(QString PCEigenvectorsName READ getPCEigenvectorsName WRITE setPCEigenvectorsName)
  PYB11_PROPERTY(int MatrixApproach READ getMatrixApproach WRITE setMatrixApproach)
  PYB11_PROPERTY(bool ProjectDataSpace READ getProjectDataSpace WRITE setProjectDataSpace)
  PYB11_PROPERTY(int NumberOfDimensionsForProjection READ getNumberOfDimensionsForProjection WRITE setNumberOfDimensionsForProjection)
  PYB11_PROPERTY(DataArrayPath ProjectedDataSpaceArrayPath READ getProjectedDataSpaceArrayPath WRITE setProjectedDataSpaceArrayPath)
#endif

public:
  using Self = PrincipalComponentAnalysis;
  using Pointer = std::shared_ptr<Self>;
  using ConstPointer = std::shared_ptr<const Self>;
  using WeakPointer = std::weak_ptr<Self>;
  using ConstWeakPointer = std::weak_ptr<const Self>;
  static Pointer NullPointer();

  static std::shared_ptr<PrincipalComponentAnalysis> New();

  /**
   * @brief Returns the name of the class for PrincipalComponentAnalysis
   */
  QString getNameOfClass() const override;
  /**
   * @brief Returns the name of the class for PrincipalComponentAnalysis
   */
  static QString ClassName();

  ~PrincipalComponentAnalysis() override;

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
   * @brief Setter property for PCAttributeMatrixName
   */
  void setPCAttributeMatrixName(const QString& value);
  /**
   * @brief Getter property for PCAttributeMatrixName
   * @return Value of PCAttributeMatrixName
   */
  QString getPCAttributeMatrixName() const;

  Q_PROPERTY(QString PCAttributeMatrixName READ getPCAttributeMatrixName WRITE setPCAttributeMatrixName)

  /**
   * @brief Setter property for PCEigenvaluesName
   */
  void setPCEigenvaluesName(const QString& value);
  /**
   * @brief Getter property for PCEigenvaluesName
   * @return Value of PCEigenvaluesName
   */
  QString getPCEigenvaluesName() const;

  Q_PROPERTY(QString PCEigenvaluesName READ getPCEigenvaluesName WRITE setPCEigenvaluesName)

  /**
   * @brief Setter property for PCEigenvectorsName
   */
  void setPCEigenvectorsName(const QString& value);
  /**
   * @brief Getter property for PCEigenvectorsName
   * @return Value of PCEigenvectorsName
   */
  QString getPCEigenvectorsName() const;

  Q_PROPERTY(QString PCEigenvectorsName READ getPCEigenvectorsName WRITE setPCEigenvectorsName)

  /**
   * @brief Setter property for MatrixApproach
   */
  void setMatrixApproach(int value);
  /**
   * @brief Getter property for MatrixApproach
   * @return Value of MatrixApproach
   */
  int getMatrixApproach() const;

  Q_PROPERTY(int MatrixApproach READ getMatrixApproach WRITE setMatrixApproach)

  /**
   * @brief Setter property for ProjectDataSpace
   */
  void setProjectDataSpace(bool value);
  /**
   * @brief Getter property for ProjectDataSpace
   * @return Value of ProjectDataSpace
   */
  bool getProjectDataSpace() const;

  Q_PROPERTY(bool ProjectDataSpace READ getProjectDataSpace WRITE setProjectDataSpace)

  /**
   * @brief Setter property for NumberOfDimensionsForProjection
   */
  void setNumberOfDimensionsForProjection(int value);
  /**
   * @brief Getter property for NumberOfDimensionsForProjection
   * @return Value of NumberOfDimensionsForProjection
   */
  int getNumberOfDimensionsForProjection() const;

  Q_PROPERTY(int NumberOfDimensionsForProjection READ getNumberOfDimensionsForProjection WRITE setNumberOfDimensionsForProjection)

  /**
   * @brief Setter property for ProjectedDataSpaceArrayPath
   */
  void setProjectedDataSpaceArrayPath(const DataArrayPath& value);
  /**
   * @brief Getter property for ProjectedDataSpaceArrayPath
   * @return Value of ProjectedDataSpaceArrayPath
   */
  DataArrayPath getProjectedDataSpaceArrayPath() const;

  Q_PROPERTY(DataArrayPath ProjectedDataSpaceArrayPath READ getProjectedDataSpaceArrayPath WRITE setProjectedDataSpaceArrayPath)

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

  /**
  * @brief preflight Reimplemented from @see AbstractFilter class
  */
  void preflight() override;

protected:
  PrincipalComponentAnalysis();

  /**
   * @brief dataCheck Checks for the appropriate parameter values and availability of arrays
   */
  void dataCheck();

  /**
   * @brief Initializes all the private instance variables.
   */
  void initialize();

private:
  std::weak_ptr<DataArray<double>> m_PCEigenvaluesPtr;
  double* m_PCEigenvalues = nullptr;
  std::weak_ptr<DataArray<double>> m_PCEigenvectorsPtr;
  double* m_PCEigenvectors = nullptr;
  std::weak_ptr<DataArray<double>> m_ProjectedDataSpacePtr;
  double* m_ProjectedDataSpace = nullptr;

  QVector<DataArrayPath> m_SelectedDataArrayPaths = {};
  QString m_PCAttributeMatrixName = {};
  QString m_PCEigenvaluesName = {};
  QString m_PCEigenvectorsName = {};
  int m_MatrixApproach = {};
  bool m_ProjectDataSpace = {};
  int m_NumberOfDimensionsForProjection = {};
  DataArrayPath m_ProjectedDataSpaceArrayPath = {};

  std::vector<IDataArray::WeakPointer> m_SelectedWeakPtrVector;

public:
  PrincipalComponentAnalysis(const PrincipalComponentAnalysis&) = delete; // Copy Constructor Not Implemented
  PrincipalComponentAnalysis(PrincipalComponentAnalysis&&) = delete;      // Move Constructor Not Implemented
  PrincipalComponentAnalysis& operator=(const PrincipalComponentAnalysis&) = delete; // Copy Assignment Not Implemented
  PrincipalComponentAnalysis& operator=(PrincipalComponentAnalysis&&) = delete;      // Move Assignment Not Implemented
};

