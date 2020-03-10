/* ============================================================================
 * Copyright (c) 2011 Michael A. Jackson (BlueQuartz Software)
 * Copyright (c) 2011 Dr. Michael A. Groeber (US Air Force Research Laboratories)
 * Copyright (c) 2014 Dr. Joseph C. Tucker (UES, Inc.)
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
 * Neither the name of Joseph C. Tucker, Michael A. Groeber, Michael A. Jackson,
 * UES, Inc., the US Air Force, BlueQuartz Software nor the names of its contributors
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
 *  This code was written under United States Air Force Contract number
 *                   FA8650-07-D-5800 and FA8650-10-D-5226
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
#pragma once

#include <memory>

#include <QtCore/QString>

#include "SIMPLib/DataArrays/IDataArray.h"
#include "SIMPLib/DataArrays/NeighborList.hpp"
#include "SIMPLib/DataContainers/DataContainer.h"
#include "SIMPLib/FilterParameters/FloatVec3FilterParameter.h"
#include "SIMPLib/Filtering/AbstractFilter.h"
#include "SIMPLib/SIMPLib.h"
#include "SIMPLib/DataArrays/DataArray.hpp"

#include "DREAM3DReview/DREAM3DReviewDLLExport.h"

/**
 * @class TiDwellFatigueCrystallographicAnalysis TiDwellFatigueCrystallographicAnalysis.h /TransformationPhase/TiDwellFatigueCrystallographicAnalysis.h
 * @brief
 * @author Joseph C. Tucker
 * @date Jan 23, 2015
 * @version 5.1
 */
class DREAM3DReview_EXPORT TiDwellFatigueCrystallographicAnalysis : public AbstractFilter
{
  Q_OBJECT
public:
  using Self = TiDwellFatigueCrystallographicAnalysis;
  using Pointer = std::shared_ptr<Self>;
  using ConstPointer = std::shared_ptr<const Self>;
  using WeakPointer = std::weak_ptr<Self>;
  using ConstWeakPointer = std::weak_ptr<const Self>;
  static Pointer NullPointer();

  static std::shared_ptr<TiDwellFatigueCrystallographicAnalysis> New();

  /**
   * @brief Returns the name of the class for TiDwellFatigueCrystallographicAnalysis
   */
  QString getNameOfClass() const override;
  /**
   * @brief Returns the name of the class for TiDwellFatigueCrystallographicAnalysis
   */
  static QString ClassName();

  ~TiDwellFatigueCrystallographicAnalysis() override;
  /**
   * @brief Setter property for DataContainerName
   */
  void setDataContainerName(const DataArrayPath& value);
  /**
   * @brief Getter property for DataContainerName
   * @return Value of DataContainerName
   */
  DataArrayPath getDataContainerName() const;

  Q_PROPERTY(DataArrayPath DataContainerName READ getDataContainerName WRITE setDataContainerName)

  /* Place your input parameters here. You can use some of the DREAM3D Macros if you want to */
  /**
   * @brief Setter property for AlphaGlobPhasePresent
   */
  void setAlphaGlobPhasePresent(bool value);
  /**
   * @brief Getter property for AlphaGlobPhasePresent
   * @return Value of AlphaGlobPhasePresent
   */
  bool getAlphaGlobPhasePresent() const;

  Q_PROPERTY(bool AlphaGlobPhasePresent READ getAlphaGlobPhasePresent WRITE setAlphaGlobPhasePresent)
  /**
   * @brief Setter property for AlphaGlobPhase
   */
  void setAlphaGlobPhase(int value);
  /**
   * @brief Getter property for AlphaGlobPhase
   * @return Value of AlphaGlobPhase
   */
  int getAlphaGlobPhase() const;

  Q_PROPERTY(int AlphaGlobPhase READ getAlphaGlobPhase WRITE setAlphaGlobPhase)
  /**
   * @brief Setter property for MTRPhase
   */
  void setMTRPhase(int value);
  /**
   * @brief Getter property for MTRPhase
   * @return Value of MTRPhase
   */
  int getMTRPhase() const;

  Q_PROPERTY(int MTRPhase READ getMTRPhase WRITE setMTRPhase)
  /**
   * @brief Setter property for LatticeParameterA
   */
  void setLatticeParameterA(float value);
  /**
   * @brief Getter property for LatticeParameterA
   * @return Value of LatticeParameterA
   */
  float getLatticeParameterA() const;

  Q_PROPERTY(float LatticeParameterA READ getLatticeParameterA WRITE setLatticeParameterA)
  /**
   * @brief Setter property for LatticeParameterC
   */
  void setLatticeParameterC(float value);
  /**
   * @brief Getter property for LatticeParameterC
   * @return Value of LatticeParameterC
   */
  float getLatticeParameterC() const;

  Q_PROPERTY(float LatticeParameterC READ getLatticeParameterC WRITE setLatticeParameterC)
  /**
   * @brief Setter property for StressAxis
   */
  void setStressAxis(const FloatVec3Type& value);
  /**
   * @brief Getter property for StressAxis
   * @return Value of StressAxis
   */
  FloatVec3Type getStressAxis() const;

  Q_PROPERTY(FloatVec3Type StressAxis READ getStressAxis WRITE setStressAxis)
  /**
   * @brief Setter property for SubsurfaceDistance
   */
  void setSubsurfaceDistance(int value);
  /**
   * @brief Getter property for SubsurfaceDistance
   * @return Value of SubsurfaceDistance
   */
  int getSubsurfaceDistance() const;

  Q_PROPERTY(int SubsurfaceDistance READ getSubsurfaceDistance WRITE setSubsurfaceDistance)
  Q_PROPERTY(float ConsiderationFraction READ getConsiderationFraction WRITE setConsiderationFraction)
  /**
   * @brief Setter property for ConsiderationFraction
   */
  void setConsiderationFraction(float value);
  /**
   * @brief Getter property for ConsiderationFraction
   * @return Value of ConsiderationFraction
   */
  float getConsiderationFraction() const;

  /**
   * @brief Setter property for DoNotAssumeInitiatorPresence
   */
  void setDoNotAssumeInitiatorPresence(bool value);
  /**
   * @brief Getter property for DoNotAssumeInitiatorPresence
   * @return Value of DoNotAssumeInitiatorPresence
   */
  bool getDoNotAssumeInitiatorPresence() const;

  Q_PROPERTY(bool DoNotAssumeInitiatorPresence READ getDoNotAssumeInitiatorPresence WRITE setDoNotAssumeInitiatorPresence)
  /**
   * @brief Setter property for InitiatorLowerThreshold
   */
  void setInitiatorLowerThreshold(float value);
  /**
   * @brief Getter property for InitiatorLowerThreshold
   * @return Value of InitiatorLowerThreshold
   */
  float getInitiatorLowerThreshold() const;

  Q_PROPERTY(float InitiatorLowerThreshold READ getInitiatorLowerThreshold WRITE setInitiatorLowerThreshold)
  /**
   * @brief Setter property for InitiatorUpperThreshold
   */
  void setInitiatorUpperThreshold(float value);
  /**
   * @brief Getter property for InitiatorUpperThreshold
   * @return Value of InitiatorUpperThreshold
   */
  float getInitiatorUpperThreshold() const;

  Q_PROPERTY(float InitiatorUpperThreshold READ getInitiatorUpperThreshold WRITE setInitiatorUpperThreshold)
  /**
   * @brief Setter property for HardFeatureLowerThreshold
   */
  void setHardFeatureLowerThreshold(float value);
  /**
   * @brief Getter property for HardFeatureLowerThreshold
   * @return Value of HardFeatureLowerThreshold
   */
  float getHardFeatureLowerThreshold() const;

  Q_PROPERTY(float HardFeatureLowerThreshold READ getHardFeatureLowerThreshold WRITE setHardFeatureLowerThreshold)
  /**
   * @brief Setter property for HardFeatureUpperThreshold
   */
  void setHardFeatureUpperThreshold(float value);
  /**
   * @brief Getter property for HardFeatureUpperThreshold
   * @return Value of HardFeatureUpperThreshold
   */
  float getHardFeatureUpperThreshold() const;

  Q_PROPERTY(float HardFeatureUpperThreshold READ getHardFeatureUpperThreshold WRITE setHardFeatureUpperThreshold)
  /**
   * @brief Setter property for SoftFeatureLowerThreshold
   */
  void setSoftFeatureLowerThreshold(float value);
  /**
   * @brief Getter property for SoftFeatureLowerThreshold
   * @return Value of SoftFeatureLowerThreshold
   */
  float getSoftFeatureLowerThreshold() const;

  Q_PROPERTY(float SoftFeatureLowerThreshold READ getSoftFeatureLowerThreshold WRITE setSoftFeatureLowerThreshold)
  /**
   * @brief Setter property for SoftFeatureUpperThreshold
   */
  void setSoftFeatureUpperThreshold(float value);
  /**
   * @brief Getter property for SoftFeatureUpperThreshold
   * @return Value of SoftFeatureUpperThreshold
   */
  float getSoftFeatureUpperThreshold() const;

  Q_PROPERTY(float SoftFeatureUpperThreshold READ getSoftFeatureUpperThreshold WRITE setSoftFeatureUpperThreshold)

  /**
   * @brief Setter property for NewCellFeatureAttributeMatrixName
   */
  void setNewCellFeatureAttributeMatrixName(const QString& value);
  /**
   * @brief Getter property for NewCellFeatureAttributeMatrixName
   * @return Value of NewCellFeatureAttributeMatrixName
   */
  QString getNewCellFeatureAttributeMatrixName() const;

  Q_PROPERTY(QString NewCellFeatureAttributeMatrixName READ getNewCellFeatureAttributeMatrixName WRITE setNewCellFeatureAttributeMatrixName)

  /**
   * @brief Setter property for SelectedFeaturesArrayName
   */
  void setSelectedFeaturesArrayName(const QString& value);
  /**
   * @brief Getter property for SelectedFeaturesArrayName
   * @return Value of SelectedFeaturesArrayName
   */
  QString getSelectedFeaturesArrayName() const;

  Q_PROPERTY(QString SelectedFeaturesArrayName READ getSelectedFeaturesArrayName WRITE setSelectedFeaturesArrayName)

  /**
   * @brief Setter property for InitiatorsArrayName
   */
  void setInitiatorsArrayName(const QString& value);
  /**
   * @brief Getter property for InitiatorsArrayName
   * @return Value of InitiatorsArrayName
   */
  QString getInitiatorsArrayName() const;

  Q_PROPERTY(QString InitiatorsArrayName READ getInitiatorsArrayName WRITE setInitiatorsArrayName)

  /**
   * @brief Setter property for HardFeaturesArrayName
   */
  void setHardFeaturesArrayName(const QString& value);
  /**
   * @brief Getter property for HardFeaturesArrayName
   * @return Value of HardFeaturesArrayName
   */
  QString getHardFeaturesArrayName() const;

  Q_PROPERTY(QString HardFeaturesArrayName READ getHardFeaturesArrayName WRITE setHardFeaturesArrayName)

  /**
   * @brief Setter property for SoftFeaturesArrayName
   */
  void setSoftFeaturesArrayName(const QString& value);
  /**
   * @brief Getter property for SoftFeaturesArrayName
   * @return Value of SoftFeaturesArrayName
   */
  QString getSoftFeaturesArrayName() const;

  Q_PROPERTY(QString SoftFeaturesArrayName READ getSoftFeaturesArrayName WRITE setSoftFeaturesArrayName)

  /**
   * @brief Setter property for HardSoftGroupsArrayName
   */
  void setHardSoftGroupsArrayName(const QString& value);
  /**
   * @brief Getter property for HardSoftGroupsArrayName
   * @return Value of HardSoftGroupsArrayName
   */
  QString getHardSoftGroupsArrayName() const;

  Q_PROPERTY(QString HardSoftGroupsArrayName READ getHardSoftGroupsArrayName WRITE setHardSoftGroupsArrayName)

  /**
   * @brief Setter property for CellFeatureAttributeMatrixName
   */
  void setCellFeatureAttributeMatrixName(const QString& value);
  /**
   * @brief Getter property for CellFeatureAttributeMatrixName
   * @return Value of CellFeatureAttributeMatrixName
   */
  QString getCellFeatureAttributeMatrixName() const;

  Q_PROPERTY(QString CellFeatureAttributeMatrixName READ getCellFeatureAttributeMatrixName WRITE setCellFeatureAttributeMatrixName)

  /**
   * @brief Setter property for CellFeatureAttributeMatrixPath
   */
  void setCellFeatureAttributeMatrixPath(const DataArrayPath& value);
  /**
   * @brief Getter property for CellFeatureAttributeMatrixPath
   * @return Value of CellFeatureAttributeMatrixPath
   */
  DataArrayPath getCellFeatureAttributeMatrixPath() const;

  Q_PROPERTY(DataArrayPath CellFeatureAttributeMatrixPath READ getCellFeatureAttributeMatrixPath WRITE setCellFeatureAttributeMatrixPath)

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
   * @brief Setter property for CellParentIdsArrayName
   */
  void setCellParentIdsArrayName(const QString& value);
  /**
   * @brief Getter property for CellParentIdsArrayName
   * @return Value of CellParentIdsArrayName
   */
  QString getCellParentIdsArrayName() const;

  Q_PROPERTY(QString CellParentIdsArrayName READ getCellParentIdsArrayName WRITE setCellParentIdsArrayName)

  /**
   * @brief Setter property for FeatureParentIdsArrayName
   */
  void setFeatureParentIdsArrayName(const QString& value);
  /**
   * @brief Getter property for FeatureParentIdsArrayName
   * @return Value of FeatureParentIdsArrayName
   */
  QString getFeatureParentIdsArrayName() const;

  Q_PROPERTY(QString FeatureParentIdsArrayName READ getFeatureParentIdsArrayName WRITE setFeatureParentIdsArrayName)

  /**
   * @brief Setter property for ActiveArrayName
   */
  void setActiveArrayName(const QString& value);
  /**
   * @brief Getter property for ActiveArrayName
   * @return Value of ActiveArrayName
   */
  QString getActiveArrayName() const;

  Q_PROPERTY(QString ActiveArrayName READ getActiveArrayName WRITE setActiveArrayName)

  /**
   * @brief Setter property for FeatureEulerAnglesArrayPath
   */
  void setFeatureEulerAnglesArrayPath(const DataArrayPath& value);
  /**
   * @brief Getter property for FeatureEulerAnglesArrayPath
   * @return Value of FeatureEulerAnglesArrayPath
   */
  DataArrayPath getFeatureEulerAnglesArrayPath() const;

  Q_PROPERTY(DataArrayPath FeatureEulerAnglesArrayPath READ getFeatureEulerAnglesArrayPath WRITE setFeatureEulerAnglesArrayPath)

  /**
   * @brief Setter property for FeaturePhasesArrayPath
   */
  void setFeaturePhasesArrayPath(const DataArrayPath& value);
  /**
   * @brief Getter property for FeaturePhasesArrayPath
   * @return Value of FeaturePhasesArrayPath
   */
  DataArrayPath getFeaturePhasesArrayPath() const;

  Q_PROPERTY(DataArrayPath FeaturePhasesArrayPath READ getFeaturePhasesArrayPath WRITE setFeaturePhasesArrayPath)

  /**
   * @brief Setter property for NeighborListArrayPath
   */
  void setNeighborListArrayPath(const DataArrayPath& value);
  /**
   * @brief Getter property for NeighborListArrayPath
   * @return Value of NeighborListArrayPath
   */
  DataArrayPath getNeighborListArrayPath() const;

  Q_PROPERTY(DataArrayPath NeighborListArrayPath READ getNeighborListArrayPath WRITE setNeighborListArrayPath)

  /**
   * @brief Setter property for CentroidsArrayPath
   */
  void setCentroidsArrayPath(const DataArrayPath& value);
  /**
   * @brief Getter property for CentroidsArrayPath
   * @return Value of CentroidsArrayPath
   */
  DataArrayPath getCentroidsArrayPath() const;

  Q_PROPERTY(DataArrayPath CentroidsArrayPath READ getCentroidsArrayPath WRITE setCentroidsArrayPath)

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
   * @brief readFilterParametersFromJson Reads the filter parameters from a file
   * @param reader Reader that is used to read the parameters from a file
   */
  void readFilterParameters(QJsonObject& obj) override;

  /**
   * @brief writeFilterParametersToJson Writes the filter parameters to a file
   * @param root The root json object
   */
  void writeFilterParameters(QJsonObject& obj) const override;

  /**
   * @brief execute Reimplemented from @see AbstractFilter class
   */
  void execute() override;

protected:
  TiDwellFatigueCrystallographicAnalysis();

  bool determine_subsurfacefeatures(int index);
  bool determine_hardfeatures(int index);
  void determine_initiators(int index);
  void determine_softfeatures(int index);
  void group_flaggedfeatures(int index);
  void assign_HardSoftGroups(int index);
  float find_angle(float g[3][3], float planeNormalU, float planeNormalV, float planeNormalW);

  /**
   * @brief dataCheck Checks for the appropriate parameter values and availability of arrays
   */
  void dataCheck() override;

  /**
   * @brief Initializes all the private instance variables.
   */
  void initialize();

  void updateFeatureInstancePointers();

private:
  std::weak_ptr<DataArray<bool>> m_SelectedFeaturesPtr;
  bool* m_SelectedFeatures = nullptr;
  std::weak_ptr<DataArray<bool>> m_InitiatorsPtr;
  bool* m_Initiators = nullptr;
  std::weak_ptr<DataArray<bool>> m_HardFeaturesPtr;
  bool* m_HardFeatures = nullptr;
  std::weak_ptr<DataArray<bool>> m_SoftFeaturesPtr;
  bool* m_SoftFeatures = nullptr;
  std::weak_ptr<DataArray<bool>> m_HardSoftGroupsPtr;
  bool* m_HardSoftGroups = nullptr;
  std::weak_ptr<DataArray<int32_t>> m_FeatureIdsPtr;
  int32_t* m_FeatureIds = nullptr;
  std::weak_ptr<DataArray<int32_t>> m_CellParentIdsPtr;
  int32_t* m_CellParentIds = nullptr;
  std::weak_ptr<DataArray<int32_t>> m_FeatureParentIdsPtr;
  int32_t* m_FeatureParentIds = nullptr;
  std::weak_ptr<DataArray<bool>> m_ActivePtr;
  bool* m_Active = nullptr;
  std::weak_ptr<DataArray<float>> m_FeatureEulerAnglesPtr;
  float* m_FeatureEulerAngles = nullptr;
  std::weak_ptr<DataArray<int32_t>> m_FeaturePhasesPtr;
  int32_t* m_FeaturePhases = nullptr;
  std::weak_ptr<DataArray<float>> m_CentroidsPtr;
  float* m_Centroids = nullptr;
  std::weak_ptr<DataArray<uint32_t>> m_CrystalStructuresPtr;
  uint32_t* m_CrystalStructures = nullptr;

  DataArrayPath m_DataContainerName = {};
  bool m_AlphaGlobPhasePresent = {};
  int m_AlphaGlobPhase = {};
  int m_MTRPhase = {};
  float m_LatticeParameterA = {};
  float m_LatticeParameterC = {};
  FloatVec3Type m_StressAxis = {};
  int m_SubsurfaceDistance = {};
  float m_ConsiderationFraction = {};
  bool m_DoNotAssumeInitiatorPresence = {};
  float m_InitiatorLowerThreshold = {};
  float m_InitiatorUpperThreshold = {};
  float m_HardFeatureLowerThreshold = {};
  float m_HardFeatureUpperThreshold = {};
  float m_SoftFeatureLowerThreshold = {};
  float m_SoftFeatureUpperThreshold = {};
  QString m_NewCellFeatureAttributeMatrixName = {};
  QString m_SelectedFeaturesArrayName = {};
  QString m_InitiatorsArrayName = {};
  QString m_HardFeaturesArrayName = {};
  QString m_SoftFeaturesArrayName = {};
  QString m_HardSoftGroupsArrayName = {};
  QString m_CellFeatureAttributeMatrixName = {};
  DataArrayPath m_CellFeatureAttributeMatrixPath = {};
  DataArrayPath m_FeatureIdsArrayPath = {};
  QString m_CellParentIdsArrayName = {};
  QString m_FeatureParentIdsArrayName = {};
  QString m_ActiveArrayName = {};
  DataArrayPath m_FeatureEulerAnglesArrayPath = {};
  DataArrayPath m_FeaturePhasesArrayPath = {};
  DataArrayPath m_NeighborListArrayPath = {};
  DataArrayPath m_CentroidsArrayPath = {};
  DataArrayPath m_CrystalStructuresArrayPath = {};

  // Feature Data - make sure these are all initialized to nullptr in the constructor

  NeighborList<int>::WeakPointer m_NeighborList;
  NeighborList<int>::WeakPointer m_ParentNeighborList;

  // Ensemble Data - make sure these are all initialized to nullptr in the constructor

public:
  TiDwellFatigueCrystallographicAnalysis(const TiDwellFatigueCrystallographicAnalysis&) = delete;            // Copy Constructor Not Implemented
  TiDwellFatigueCrystallographicAnalysis(TiDwellFatigueCrystallographicAnalysis&&) = delete;                 // Move Constructor Not Implemented
  TiDwellFatigueCrystallographicAnalysis& operator=(const TiDwellFatigueCrystallographicAnalysis&) = delete; // Copy Assignment Not Implemented
  TiDwellFatigueCrystallographicAnalysis& operator=(TiDwellFatigueCrystallographicAnalysis&&) = delete;      // Move Assignment Not Implemented
};
