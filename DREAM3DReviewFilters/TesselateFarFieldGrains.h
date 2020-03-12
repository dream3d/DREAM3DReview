/* ============================================================================
 * Copyright (c) 2011 Michael A. Jackson (BlueQuartz Software)
 * Copyright (c) 2011 Dr. Michael A. Groeber (US Air Force Research Laboratories)
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
 * Neither the name of Michael A. Groeber, Michael A. Jackson, the US Air Force,
 * BlueQuartz Software nor the names of its contributors may be used to endorse
 * or promote products derived from this software without specific prior written
 * permission.
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
 *  This code was written under United States Air Force Contract number
 *                           FA8650-07-D-5800
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#pragma once

#include <memory>

#include <QtCore/QString>
#include <map>
#include <vector>

#include "DREAM3DReview/DREAM3DReviewConstants.h"
#include "OrientationLib/LaueOps/OrthoRhombicOps.h"
#include "SIMPLib/SIMPLib.h"
#include "SIMPLib/DataArrays/DataArray.hpp"
#include "SIMPLib/DataArrays/IDataArray.h"
#include "SIMPLib/DataArrays/StatsDataArray.h"
#include "SIMPLib/DataContainers/DataContainer.h"
#include "SIMPLib/FilterParameters/FileListInfoFilterParameter.h"
#include "SIMPLib/FilterParameters/StackFileListInfo.h"
#include "SIMPLib/Filtering/AbstractFilter.h"
#include "SIMPLib/Geometry/ShapeOps/ShapeOps.h"
#include "SIMPLib/StatsData/StatsData.h"

#include "DREAM3DReview/DREAM3DReviewDLLExport.h"

/**
 * @class TesselateFarFieldGrains TesselateFarFieldGrains.h DREAM3DLib/SyntheticBuilderFilters/TesselateFarFieldGrains.h
 * @brief
 * @author
 * @date Nov 19, 2011
 * @version 1.0
 */
class DREAM3DReview_EXPORT TesselateFarFieldGrains : public AbstractFilter
{
  Q_OBJECT

#ifdef SIMPL_ENABLE_PYTHON
  PYB11_CREATE_BINDINGS(TesselateFarFieldGrains SUPERCLASS AbstractFilter)
  PYB11_SHARED_POINTERS(TesselateFarFieldGrains)
  PYB11_FILTER_NEW_MACRO(TesselateFarFieldGrains)
  PYB11_FILTER_PARAMETER(DataArrayPath, OutputCellAttributeMatrixName)
  PYB11_FILTER_PARAMETER(QString, FeatureIdsArrayName)
  PYB11_FILTER_PARAMETER(QString, CellPhasesArrayName)
  PYB11_FILTER_PARAMETER(QString, SlabIdArrayName)
  PYB11_FILTER_PARAMETER(QString, FeaturePhasesArrayName)
  PYB11_FILTER_PARAMETER(QString, FeatureEulerAnglesArrayName)
  PYB11_FILTER_PARAMETER(QString, ElasticStrainsArrayName)
  PYB11_FILTER_PARAMETER(QString, CentroidsArrayName)
  PYB11_FILTER_PARAMETER(QString, VolumesArrayName)
  PYB11_FILTER_PARAMETER(QString, AxisLengthsArrayName)
  PYB11_FILTER_PARAMETER(QString, AxisEulerAnglesArrayName)
  PYB11_FILTER_PARAMETER(QString, Omega3sArrayName)
  PYB11_FILTER_PARAMETER(QString, EquivalentDiametersArrayName)
  PYB11_FILTER_PARAMETER(QString, CrystalStructuresArrayName)
  PYB11_FILTER_PARAMETER(DataArrayPath, MaskArrayPath)
  PYB11_FILTER_PARAMETER(StackFileListInfo, FeatureInputFileListInfo)
  PYB11_PROPERTY(DataArrayPath OutputCellAttributeMatrixName READ getOutputCellAttributeMatrixName WRITE setOutputCellAttributeMatrixName)
  PYB11_PROPERTY(QString FeatureIdsArrayName READ getFeatureIdsArrayName WRITE setFeatureIdsArrayName)
  PYB11_PROPERTY(QString CellPhasesArrayName READ getCellPhasesArrayName WRITE setCellPhasesArrayName)
  PYB11_PROPERTY(QString SlabIdArrayName READ getSlabIdArrayName WRITE setSlabIdArrayName)
  PYB11_PROPERTY(QString FeaturePhasesArrayName READ getFeaturePhasesArrayName WRITE setFeaturePhasesArrayName)
  PYB11_PROPERTY(QString FeatureEulerAnglesArrayName READ getFeatureEulerAnglesArrayName WRITE setFeatureEulerAnglesArrayName)
  PYB11_PROPERTY(QString ElasticStrainsArrayName READ getElasticStrainsArrayName WRITE setElasticStrainsArrayName)
  PYB11_PROPERTY(QString CentroidsArrayName READ getCentroidsArrayName WRITE setCentroidsArrayName)
  PYB11_PROPERTY(QString VolumesArrayName READ getVolumesArrayName WRITE setVolumesArrayName)
  PYB11_PROPERTY(QString AxisLengthsArrayName READ getAxisLengthsArrayName WRITE setAxisLengthsArrayName)
  PYB11_PROPERTY(QString AxisEulerAnglesArrayName READ getAxisEulerAnglesArrayName WRITE setAxisEulerAnglesArrayName)
  PYB11_PROPERTY(QString Omega3sArrayName READ getOmega3sArrayName WRITE setOmega3sArrayName)
  PYB11_PROPERTY(QString EquivalentDiametersArrayName READ getEquivalentDiametersArrayName WRITE setEquivalentDiametersArrayName)
  PYB11_PROPERTY(QString CrystalStructuresArrayName READ getCrystalStructuresArrayName WRITE setCrystalStructuresArrayName)
  PYB11_PROPERTY(DataArrayPath MaskArrayPath READ getMaskArrayPath WRITE setMaskArrayPath)
  PYB11_PROPERTY(StackFileListInfo FeatureInputFileListInfo READ getFeatureInputFileListInfo WRITE setFeatureInputFileListInfo)
#endif

public:
  using Self = TesselateFarFieldGrains;
  using Pointer = std::shared_ptr<Self>;
  using ConstPointer = std::shared_ptr<const Self>;
  using WeakPointer = std::weak_ptr<Self>;
  using ConstWeakPointer = std::weak_ptr<const Self>;
  static Pointer NullPointer();

  static std::shared_ptr<TesselateFarFieldGrains> New();

  /**
   * @brief Returns the name of the class for TesselateFarFieldGrains
   */
  QString getNameOfClass() const override;
  /**
   * @brief Returns the name of the class for TesselateFarFieldGrains
   */
  static QString ClassName();

  ~TesselateFarFieldGrains() override;
  /**
   * @brief Setter property for OutputCellAttributeMatrixName
   */
  void setOutputCellAttributeMatrixName(const DataArrayPath& value);
  /**
   * @brief Getter property for OutputCellAttributeMatrixName
   * @return Value of OutputCellAttributeMatrixName
   */
  DataArrayPath getOutputCellAttributeMatrixName() const;

  Q_PROPERTY(DataArrayPath OutputCellAttributeMatrixName READ getOutputCellAttributeMatrixName WRITE setOutputCellAttributeMatrixName)

  /**
   * @brief Setter property for OutputCellFeatureAttributeMatrixName
   */
  void setOutputCellFeatureAttributeMatrixName(const QString& value);
  /**
   * @brief Getter property for OutputCellFeatureAttributeMatrixName
   * @return Value of OutputCellFeatureAttributeMatrixName
   */
  QString getOutputCellFeatureAttributeMatrixName() const;

  Q_PROPERTY(QString OutputCellFeatureAttributeMatrixName READ getOutputCellFeatureAttributeMatrixName WRITE setOutputCellFeatureAttributeMatrixName)

  /**
   * @brief Setter property for OutputCellEnsembleAttributeMatrixName
   */
  void setOutputCellEnsembleAttributeMatrixName(const QString& value);
  /**
   * @brief Getter property for OutputCellEnsembleAttributeMatrixName
   * @return Value of OutputCellEnsembleAttributeMatrixName
   */
  QString getOutputCellEnsembleAttributeMatrixName() const;

  Q_PROPERTY(QString OutputCellEnsembleAttributeMatrixName READ getOutputCellEnsembleAttributeMatrixName WRITE setOutputCellEnsembleAttributeMatrixName)

  /**
   * @brief Setter property for FeatureIdsArrayName
   */
  void setFeatureIdsArrayName(const QString& value);
  /**
   * @brief Getter property for FeatureIdsArrayName
   * @return Value of FeatureIdsArrayName
   */
  QString getFeatureIdsArrayName() const;

  Q_PROPERTY(QString FeatureIdsArrayName READ getFeatureIdsArrayName WRITE setFeatureIdsArrayName)

  /**
   * @brief Setter property for CellPhasesArrayName
   */
  void setCellPhasesArrayName(const QString& value);
  /**
   * @brief Getter property for CellPhasesArrayName
   * @return Value of CellPhasesArrayName
   */
  QString getCellPhasesArrayName() const;

  Q_PROPERTY(QString CellPhasesArrayName READ getCellPhasesArrayName WRITE setCellPhasesArrayName)

  /**
   * @brief Setter property for SlabIdArrayName
   */
  void setSlabIdArrayName(const QString& value);
  /**
   * @brief Getter property for SlabIdArrayName
   * @return Value of SlabIdArrayName
   */
  QString getSlabIdArrayName() const;

  /**
   * @brief Setter property for FeaturePhasesArrayName
   */
  void setFeaturePhasesArrayName(const QString& value);
  /**
   * @brief Getter property for FeaturePhasesArrayName
   * @return Value of FeaturePhasesArrayName
   */
  QString getFeaturePhasesArrayName() const;

  Q_PROPERTY(QString FeaturePhasesArrayName READ getFeaturePhasesArrayName WRITE setFeaturePhasesArrayName)

  /**
   * @brief Setter property for FeatureEulerAnglesArrayName
   */
  void setFeatureEulerAnglesArrayName(const QString& value);
  /**
   * @brief Getter property for FeatureEulerAnglesArrayName
   * @return Value of FeatureEulerAnglesArrayName
   */
  QString getFeatureEulerAnglesArrayName() const;

  Q_PROPERTY(QString FeatureEulerAnglesArrayName READ getFeatureEulerAnglesArrayName WRITE setFeatureEulerAnglesArrayName)

  /**
   * @brief Setter property for ElasticStrainsArrayName
   */
  void setElasticStrainsArrayName(const QString& value);
  /**
   * @brief Getter property for ElasticStrainsArrayName
   * @return Value of ElasticStrainsArrayName
   */
  QString getElasticStrainsArrayName() const;

  Q_PROPERTY(QString ElasticStrainsArrayName READ getElasticStrainsArrayName WRITE setElasticStrainsArrayName)

  /**
   * @brief Setter property for CentroidsArrayName
   */
  void setCentroidsArrayName(const QString& value);
  /**
   * @brief Getter property for CentroidsArrayName
   * @return Value of CentroidsArrayName
   */
  QString getCentroidsArrayName() const;

  Q_PROPERTY(QString CentroidsArrayName READ getCentroidsArrayName WRITE setCentroidsArrayName)

  /**
   * @brief Setter property for VolumesArrayName
   */
  void setVolumesArrayName(const QString& value);
  /**
   * @brief Getter property for VolumesArrayName
   * @return Value of VolumesArrayName
   */
  QString getVolumesArrayName() const;

  Q_PROPERTY(QString VolumesArrayName READ getVolumesArrayName WRITE setVolumesArrayName)

  /**
   * @brief Setter property for AxisLengthsArrayName
   */
  void setAxisLengthsArrayName(const QString& value);
  /**
   * @brief Getter property for AxisLengthsArrayName
   * @return Value of AxisLengthsArrayName
   */
  QString getAxisLengthsArrayName() const;

  Q_PROPERTY(QString AxisLengthsArrayName READ getAxisLengthsArrayName WRITE setAxisLengthsArrayName)

  /**
   * @brief Setter property for AxisEulerAnglesArrayName
   */
  void setAxisEulerAnglesArrayName(const QString& value);
  /**
   * @brief Getter property for AxisEulerAnglesArrayName
   * @return Value of AxisEulerAnglesArrayName
   */
  QString getAxisEulerAnglesArrayName() const;

  Q_PROPERTY(QString AxisEulerAnglesArrayName READ getAxisEulerAnglesArrayName WRITE setAxisEulerAnglesArrayName)

  /**
   * @brief Setter property for Omega3sArrayName
   */
  void setOmega3sArrayName(const QString& value);
  /**
   * @brief Getter property for Omega3sArrayName
   * @return Value of Omega3sArrayName
   */
  QString getOmega3sArrayName() const;

  Q_PROPERTY(QString Omega3sArrayName READ getOmega3sArrayName WRITE setOmega3sArrayName)

  /**
   * @brief Setter property for EquivalentDiametersArrayName
   */
  void setEquivalentDiametersArrayName(const QString& value);
  /**
   * @brief Getter property for EquivalentDiametersArrayName
   * @return Value of EquivalentDiametersArrayName
   */
  QString getEquivalentDiametersArrayName() const;

  Q_PROPERTY(QString EquivalentDiametersArrayName READ getEquivalentDiametersArrayName WRITE setEquivalentDiametersArrayName)

  /**
   * @brief Setter property for CrystalStructuresArrayName
   */
  void setCrystalStructuresArrayName(const QString& value);
  /**
   * @brief Getter property for CrystalStructuresArrayName
   * @return Value of CrystalStructuresArrayName
   */
  QString getCrystalStructuresArrayName() const;

  Q_PROPERTY(QString CrystalStructuresArrayName READ getCrystalStructuresArrayName WRITE setCrystalStructuresArrayName)

  typedef DataArray<float>::Pointer SharedFloatArray;
  typedef DataArray<int>::Pointer SharedIntArray;

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
   * @brief Setter property for FeatureInputFileListInfo
   */
  void setFeatureInputFileListInfo(const StackFileListInfo& value);
  /**
   * @brief Getter property for FeatureInputFileListInfo
   * @return Value of FeatureInputFileListInfo
   */
  StackFileListInfo getFeatureInputFileListInfo() const;

  Q_PROPERTY(StackFileListInfo FeatureInputFileListInfo READ getFeatureInputFileListInfo WRITE setFeatureInputFileListInfo)

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
   * @brief This method will read the options from a file
   * @param reader The reader that is used to read the options from a file
   */
  void readFilterParameters(AbstractFilterParametersReader* reader, int index) override;

  /**
   * @brief Reimplemented from @see AbstractFilter class
   */
  void execute() override;

protected:
  TesselateFarFieldGrains();

  /**
   * @brief dataCheck Checks for the appropriate parameter values and availability of arrays
   */
  void dataCheck() override;

  /**
   * @brief Initializes all the private instance variables.
   */
  void initialize();

  void load_features();
  void merge_twins();
  void assign_voxels();
  void assign_gaps_only();
  void assign_orientations();

private:
  std::weak_ptr<DataArray<int32_t>> m_FeatureIdsPtr;
  int32_t* m_FeatureIds = nullptr;
  std::weak_ptr<DataArray<int32_t>> m_CellPhasesPtr;
  int32_t* m_CellPhases = nullptr;
  std::weak_ptr<DataArray<bool>> m_MaskPtr;
  bool* m_Mask = nullptr;
  std::weak_ptr<DataArray<int32_t>> m_FeaturePhasesPtr;
  int32_t* m_FeaturePhases = nullptr;
  std::weak_ptr<DataArray<int32_t>> m_SlabIdPtr;
  int32_t* m_SlabId = nullptr;
  std::weak_ptr<DataArray<float>> m_FeatureEulerAnglesPtr;
  float* m_FeatureEulerAngles = nullptr;
  std::weak_ptr<DataArray<float>> m_ElasticStrainsPtr;
  float* m_ElasticStrains = nullptr;
  std::weak_ptr<DataArray<float>> m_CentroidsPtr;
  float* m_Centroids = nullptr;
  std::weak_ptr<DataArray<float>> m_VolumesPtr;
  float* m_Volumes = nullptr;
  std::weak_ptr<DataArray<float>> m_AxisLengthsPtr;
  float* m_AxisLengths = nullptr;
  std::weak_ptr<DataArray<float>> m_AxisEulerAnglesPtr;
  float* m_AxisEulerAngles = nullptr;
  std::weak_ptr<DataArray<float>> m_Omega3sPtr;
  float* m_Omega3s = nullptr;
  std::weak_ptr<DataArray<float>> m_EquivalentDiametersPtr;
  float* m_EquivalentDiameters = nullptr;
  std::weak_ptr<DataArray<unsigned int>> m_CrystalStructuresPtr;
  unsigned int* m_CrystalStructures = nullptr;

  DataArrayPath m_OutputCellAttributeMatrixName = {};
  QString m_OutputCellFeatureAttributeMatrixName = {};
  QString m_OutputCellEnsembleAttributeMatrixName = {};
  QString m_FeatureIdsArrayName = {};
  QString m_CellPhasesArrayName = {};
  QString m_SlabIdArrayName = {};
  QString m_FeaturePhasesArrayName = {};
  QString m_FeatureEulerAnglesArrayName = {};
  QString m_ElasticStrainsArrayName = {};
  QString m_CentroidsArrayName = {};
  QString m_VolumesArrayName = {};
  QString m_AxisLengthsArrayName = {};
  QString m_AxisEulerAnglesArrayName = {};
  QString m_Omega3sArrayName = {};
  QString m_EquivalentDiametersArrayName = {};
  QString m_CrystalStructuresArrayName = {};
  DataArrayPath m_MaskArrayPath = {};
  StackFileListInfo m_FeatureInputFileListInfo = {};

  int32_t* m_Neighbors;

  // Cell Data - make sure these are all initialized to nullptr in the constructor

  int8_t* m_BoundaryCells;

  // Feature Data - make sure these are all initialized to nullptr in the constructor

  // Ensemble Data - make sure these are all initialized to nullptr in the constructor

  // All other private variables
  ShapeOps::Pointer m_EllipsoidOps;

  OrthoRhombicOps::Pointer m_OrthoOps;

  unsigned long long int m_RandomSeed;

  float m_SizeX;
  float m_SizeY;
  float m_SizeZ;
  float m_TotalVol;

  std::vector<int> m_NewNames;
  std::vector<int> m_PackQualities;
  std::vector<int> m_GSizes;

  void updateFeatureInstancePointers();

  void updateEnsembleInstancePointers();

public:
  TesselateFarFieldGrains(const TesselateFarFieldGrains&) = delete;            // Copy Constructor Not Implemented
  TesselateFarFieldGrains(TesselateFarFieldGrains&&) = delete;                 // Move Constructor Not Implemented
  TesselateFarFieldGrains& operator=(const TesselateFarFieldGrains&) = delete; // Copy Assignment Not Implemented
  TesselateFarFieldGrains& operator=(TesselateFarFieldGrains&&) = delete;      // Move Assignment Not Implemented
};
