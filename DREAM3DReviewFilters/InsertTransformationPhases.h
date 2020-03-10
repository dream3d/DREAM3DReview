/* ============================================================================
 * Copyright (c) 2011 Michael A. Jackson (BlueQuartz Software)
 * Copyright (c) 2011 Dr. Michael A. Groeber (US Air Force Research Laboratories)
 * Copyright (c) 2015 Dr. Joseph C. Tucker (UES, Inc.)
 * Copyright (c) 2014 Sudipto Mandal & Prof. Anthony D. Rollett (Carnegie Mellon
 * University)
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
 * Sudipto Mandal, Anthony D. Rollett, UES, Inc., the US Air Force, BlueQuartz
 * Software, Carnegie Mellon University nor the names of its contributors may be
 * used to endorse or promote products derived from this software without specific
 * prior written permission.
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
 *     FA8650-07-D-5800 and FA8650-10-D-5226 and The Boeing Company
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#pragma once

#include <memory>

#include <vector>

#include <QtCore/QString>

#include "SIMPLib/SIMPLib.h"
#include "SIMPLib/Common/PhaseType.h"
#include "SIMPLib/Common/ShapeType.h"
#include "SIMPLib/DataArrays/IDataArray.h"
#include "SIMPLib/DataContainers/DataContainer.h"
#include "SIMPLib/FilterParameters/FloatVec3FilterParameter.h"
#include "SIMPLib/Filtering/AbstractFilter.h"
#include "SIMPLib/DataArrays/DataArray.hpp"

#include "OrientationLib/LaueOps/LaueOps.h"

#include "DREAM3DReview/DREAM3DReviewDLLExport.h"

/**
 * @class InsertTransformationPhases InsertTransformationPhases.h DREAM3DLib/SyntheticBuildingFilters/InsertTransformationPhases.h
 * @brief
 * @author Joseph C. Tucker (UES, Inc.)
 * @author Sudipto Mandal (CMU)
 * @date Feb 6, 2015
 * @version 5.1
 */
class DREAM3DReview_EXPORT InsertTransformationPhases : public AbstractFilter
{
  Q_OBJECT

#ifdef SIMPL_ENABLE_PYTHON
  PYB11_CREATE_BINDINGS(InsertTransformationPhases SUPERCLASS AbstractFilter)
  PYB11_SHARED_POINTERS(InsertTransformationPhases)
  PYB11_FILTER_NEW_MACRO(InsertTransformationPhases)
  PYB11_FILTER_PARAMETER(int, ParentPhase)
  PYB11_FILTER_PARAMETER(int, TransCrystalStruct)
  PYB11_FILTER_PARAMETER(float, TransformationPhaseMisorientation)
  PYB11_FILTER_PARAMETER(FloatVec3Type, TransformationPhaseHabitPlane)
  PYB11_FILTER_PARAMETER(bool, DefineHabitPlane)
  PYB11_FILTER_PARAMETER(bool, UseAllVariants)
  PYB11_FILTER_PARAMETER(float, CoherentFrac)
  PYB11_FILTER_PARAMETER(float, TransformationPhaseThickness)
  PYB11_FILTER_PARAMETER(int, NumTransformationPhasesPerFeature)
  PYB11_FILTER_PARAMETER(float, PeninsulaFrac)
  PYB11_FILTER_PARAMETER(DataArrayPath, StatsGenCellEnsembleAttributeMatrixPath)
  PYB11_FILTER_PARAMETER(DataArrayPath, CellFeatureAttributeMatrixName)
  PYB11_FILTER_PARAMETER(DataArrayPath, FeatureIdsArrayPath)
  PYB11_FILTER_PARAMETER(DataArrayPath, CellEulerAnglesArrayPath)
  PYB11_FILTER_PARAMETER(DataArrayPath, CellPhasesArrayPath)
  PYB11_FILTER_PARAMETER(DataArrayPath, AvgQuatsArrayPath)
  PYB11_FILTER_PARAMETER(DataArrayPath, CentroidsArrayPath)
  PYB11_FILTER_PARAMETER(DataArrayPath, EquivalentDiametersArrayPath)
  PYB11_FILTER_PARAMETER(DataArrayPath, FeatureEulerAnglesArrayPath)
  PYB11_FILTER_PARAMETER(DataArrayPath, FeaturePhasesArrayPath)
  PYB11_FILTER_PARAMETER(QString, FeatureParentIdsArrayName)
  PYB11_FILTER_PARAMETER(DataArrayPath, NumFeaturesPerParentArrayPath)
  PYB11_FILTER_PARAMETER(DataArrayPath, CrystalStructuresArrayPath)
  PYB11_FILTER_PARAMETER(DataArrayPath, PhaseTypesArrayPath)
  PYB11_FILTER_PARAMETER(DataArrayPath, ShapeTypesArrayPath)
  PYB11_FILTER_PARAMETER(DataArrayPath, NumFeaturesArrayPath)
  PYB11_PROPERTY(int ParentPhase READ getParentPhase WRITE setParentPhase)
  PYB11_PROPERTY(int TransCrystalStruct READ getTransCrystalStruct WRITE setTransCrystalStruct)
  PYB11_PROPERTY(float TransformationPhaseMisorientation READ getTransformationPhaseMisorientation WRITE setTransformationPhaseMisorientation)
  PYB11_PROPERTY(FloatVec3Type TransformationPhaseHabitPlane READ getTransformationPhaseHabitPlane WRITE setTransformationPhaseHabitPlane)
  PYB11_PROPERTY(bool DefineHabitPlane READ getDefineHabitPlane WRITE setDefineHabitPlane)
  PYB11_PROPERTY(bool UseAllVariants READ getUseAllVariants WRITE setUseAllVariants)
  PYB11_PROPERTY(float CoherentFrac READ getCoherentFrac WRITE setCoherentFrac)
  PYB11_PROPERTY(float TransformationPhaseThickness READ getTransformationPhaseThickness WRITE setTransformationPhaseThickness)
  PYB11_PROPERTY(int NumTransformationPhasesPerFeature READ getNumTransformationPhasesPerFeature WRITE setNumTransformationPhasesPerFeature)
  PYB11_PROPERTY(float PeninsulaFrac READ getPeninsulaFrac WRITE setPeninsulaFrac)
  PYB11_PROPERTY(DataArrayPath StatsGenCellEnsembleAttributeMatrixPath READ getStatsGenCellEnsembleAttributeMatrixPath WRITE setStatsGenCellEnsembleAttributeMatrixPath)
  PYB11_PROPERTY(DataArrayPath CellFeatureAttributeMatrixName READ getCellFeatureAttributeMatrixName WRITE setCellFeatureAttributeMatrixName)
  PYB11_PROPERTY(DataArrayPath FeatureIdsArrayPath READ getFeatureIdsArrayPath WRITE setFeatureIdsArrayPath)
  PYB11_PROPERTY(DataArrayPath CellEulerAnglesArrayPath READ getCellEulerAnglesArrayPath WRITE setCellEulerAnglesArrayPath)
  PYB11_PROPERTY(DataArrayPath CellPhasesArrayPath READ getCellPhasesArrayPath WRITE setCellPhasesArrayPath)
  PYB11_PROPERTY(DataArrayPath AvgQuatsArrayPath READ getAvgQuatsArrayPath WRITE setAvgQuatsArrayPath)
  PYB11_PROPERTY(DataArrayPath CentroidsArrayPath READ getCentroidsArrayPath WRITE setCentroidsArrayPath)
  PYB11_PROPERTY(DataArrayPath EquivalentDiametersArrayPath READ getEquivalentDiametersArrayPath WRITE setEquivalentDiametersArrayPath)
  PYB11_PROPERTY(DataArrayPath FeatureEulerAnglesArrayPath READ getFeatureEulerAnglesArrayPath WRITE setFeatureEulerAnglesArrayPath)
  PYB11_PROPERTY(DataArrayPath FeaturePhasesArrayPath READ getFeaturePhasesArrayPath WRITE setFeaturePhasesArrayPath)
  PYB11_PROPERTY(QString FeatureParentIdsArrayName READ getFeatureParentIdsArrayName WRITE setFeatureParentIdsArrayName)
  PYB11_PROPERTY(DataArrayPath NumFeaturesPerParentArrayPath READ getNumFeaturesPerParentArrayPath WRITE setNumFeaturesPerParentArrayPath)
  PYB11_PROPERTY(DataArrayPath CrystalStructuresArrayPath READ getCrystalStructuresArrayPath WRITE setCrystalStructuresArrayPath)
  PYB11_PROPERTY(DataArrayPath PhaseTypesArrayPath READ getPhaseTypesArrayPath WRITE setPhaseTypesArrayPath)
  PYB11_PROPERTY(DataArrayPath ShapeTypesArrayPath READ getShapeTypesArrayPath WRITE setShapeTypesArrayPath)
  PYB11_PROPERTY(DataArrayPath NumFeaturesArrayPath READ getNumFeaturesArrayPath WRITE setNumFeaturesArrayPath)
#endif

public:
  using Self = InsertTransformationPhases;
  using Pointer = std::shared_ptr<Self>;
  using ConstPointer = std::shared_ptr<const Self>;
  using WeakPointer = std::weak_ptr<Self>;
  using ConstWeakPointer = std::weak_ptr<const Self>;
  static Pointer NullPointer();

  static std::shared_ptr<InsertTransformationPhases> New();

  /**
   * @brief Returns the name of the class for InsertTransformationPhases
   */
  QString getNameOfClass() const override;
  /**
   * @brief Returns the name of the class for InsertTransformationPhases
   */
  static QString ClassName();

  ~InsertTransformationPhases() override;

  /**
   * @brief Setter property for ParentPhase
   */
  void setParentPhase(int value);
  /**
   * @brief Getter property for ParentPhase
   * @return Value of ParentPhase
   */
  int getParentPhase() const;

  Q_PROPERTY(int ParentPhase READ getParentPhase WRITE setParentPhase)

  /**
   * @brief Setter property for TransCrystalStruct
   */
  void setTransCrystalStruct(int value);
  /**
   * @brief Getter property for TransCrystalStruct
   * @return Value of TransCrystalStruct
   */
  int getTransCrystalStruct() const;

  Q_PROPERTY(int TransCrystalStruct READ getTransCrystalStruct WRITE setTransCrystalStruct)

  /**
   * @brief Setter property for TransformationPhaseMisorientation
   */
  void setTransformationPhaseMisorientation(float value);
  /**
   * @brief Getter property for TransformationPhaseMisorientation
   * @return Value of TransformationPhaseMisorientation
   */
  float getTransformationPhaseMisorientation() const;

  Q_PROPERTY(float TransformationPhaseMisorientation READ getTransformationPhaseMisorientation WRITE setTransformationPhaseMisorientation)

  /**
   * @brief Setter property for TransformationPhaseHabitPlane
   */
  void setTransformationPhaseHabitPlane(const FloatVec3Type& value);
  /**
   * @brief Getter property for TransformationPhaseHabitPlane
   * @return Value of TransformationPhaseHabitPlane
   */
  FloatVec3Type getTransformationPhaseHabitPlane() const;

  Q_PROPERTY(FloatVec3Type TransformationPhaseHabitPlane READ getTransformationPhaseHabitPlane WRITE setTransformationPhaseHabitPlane)

  /**
   * @brief Setter property for DefineHabitPlane
   */
  void setDefineHabitPlane(bool value);
  /**
   * @brief Getter property for DefineHabitPlane
   * @return Value of DefineHabitPlane
   */
  bool getDefineHabitPlane() const;

  Q_PROPERTY(bool DefineHabitPlane READ getDefineHabitPlane WRITE setDefineHabitPlane)

  /**
   * @brief Setter property for UseAllVariants
   */
  void setUseAllVariants(bool value);
  /**
   * @brief Getter property for UseAllVariants
   * @return Value of UseAllVariants
   */
  bool getUseAllVariants() const;

  Q_PROPERTY(bool UseAllVariants READ getUseAllVariants WRITE setUseAllVariants)

  /**
   * @brief Setter property for CoherentFrac
   */
  void setCoherentFrac(float value);
  /**
   * @brief Getter property for CoherentFrac
   * @return Value of CoherentFrac
   */
  float getCoherentFrac() const;

  Q_PROPERTY(float CoherentFrac READ getCoherentFrac WRITE setCoherentFrac)

  /**
   * @brief Setter property for TransformationPhaseThickness
   */
  void setTransformationPhaseThickness(float value);
  /**
   * @brief Getter property for TransformationPhaseThickness
   * @return Value of TransformationPhaseThickness
   */
  float getTransformationPhaseThickness() const;

  Q_PROPERTY(float TransformationPhaseThickness READ getTransformationPhaseThickness WRITE setTransformationPhaseThickness)

  /**
   * @brief Setter property for NumTransformationPhasesPerFeature
   */
  void setNumTransformationPhasesPerFeature(int value);
  /**
   * @brief Getter property for NumTransformationPhasesPerFeature
   * @return Value of NumTransformationPhasesPerFeature
   */
  int getNumTransformationPhasesPerFeature() const;

  Q_PROPERTY(int NumTransformationPhasesPerFeature READ getNumTransformationPhasesPerFeature WRITE setNumTransformationPhasesPerFeature)

  /**
   * @brief Setter property for PeninsulaFrac
   */
  void setPeninsulaFrac(float value);
  /**
   * @brief Getter property for PeninsulaFrac
   * @return Value of PeninsulaFrac
   */
  float getPeninsulaFrac() const;

  Q_PROPERTY(float PeninsulaFrac READ getPeninsulaFrac WRITE setPeninsulaFrac)

  /**
   * @brief Setter property for StatsGenCellEnsembleAttributeMatrixPath
   */
  void setStatsGenCellEnsembleAttributeMatrixPath(const DataArrayPath& value);
  /**
   * @brief Getter property for StatsGenCellEnsembleAttributeMatrixPath
   * @return Value of StatsGenCellEnsembleAttributeMatrixPath
   */
  DataArrayPath getStatsGenCellEnsembleAttributeMatrixPath() const;

  Q_PROPERTY(DataArrayPath StatsGenCellEnsembleAttributeMatrixPath READ getStatsGenCellEnsembleAttributeMatrixPath WRITE setStatsGenCellEnsembleAttributeMatrixPath)

  /**
   * @brief Setter property for CellFeatureAttributeMatrixName
   */
  void setCellFeatureAttributeMatrixName(const DataArrayPath& value);
  /**
   * @brief Getter property for CellFeatureAttributeMatrixName
   * @return Value of CellFeatureAttributeMatrixName
   */
  DataArrayPath getCellFeatureAttributeMatrixName() const;

  Q_PROPERTY(DataArrayPath CellFeatureAttributeMatrixName READ getCellFeatureAttributeMatrixName WRITE setCellFeatureAttributeMatrixName)

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
   * @brief Setter property for CellEulerAnglesArrayPath
   */
  void setCellEulerAnglesArrayPath(const DataArrayPath& value);
  /**
   * @brief Getter property for CellEulerAnglesArrayPath
   * @return Value of CellEulerAnglesArrayPath
   */
  DataArrayPath getCellEulerAnglesArrayPath() const;

  Q_PROPERTY(DataArrayPath CellEulerAnglesArrayPath READ getCellEulerAnglesArrayPath WRITE setCellEulerAnglesArrayPath)

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
   * @brief Setter property for AvgQuatsArrayPath
   */
  void setAvgQuatsArrayPath(const DataArrayPath& value);
  /**
   * @brief Getter property for AvgQuatsArrayPath
   * @return Value of AvgQuatsArrayPath
   */
  DataArrayPath getAvgQuatsArrayPath() const;

  Q_PROPERTY(DataArrayPath AvgQuatsArrayPath READ getAvgQuatsArrayPath WRITE setAvgQuatsArrayPath)

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
   * @brief Setter property for EquivalentDiametersArrayPath
   */
  void setEquivalentDiametersArrayPath(const DataArrayPath& value);
  /**
   * @brief Getter property for EquivalentDiametersArrayPath
   * @return Value of EquivalentDiametersArrayPath
   */
  DataArrayPath getEquivalentDiametersArrayPath() const;

  Q_PROPERTY(DataArrayPath EquivalentDiametersArrayPath READ getEquivalentDiametersArrayPath WRITE setEquivalentDiametersArrayPath)

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
   * @brief Setter property for NumFeaturesPerParentArrayPath
   */
  void setNumFeaturesPerParentArrayPath(const DataArrayPath& value);
  /**
   * @brief Getter property for NumFeaturesPerParentArrayPath
   * @return Value of NumFeaturesPerParentArrayPath
   */
  DataArrayPath getNumFeaturesPerParentArrayPath() const;

  Q_PROPERTY(DataArrayPath NumFeaturesPerParentArrayPath READ getNumFeaturesPerParentArrayPath WRITE setNumFeaturesPerParentArrayPath)

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
   * @brief Setter property for PhaseTypesArrayPath
   */
  void setPhaseTypesArrayPath(const DataArrayPath& value);
  /**
   * @brief Getter property for PhaseTypesArrayPath
   * @return Value of PhaseTypesArrayPath
   */
  DataArrayPath getPhaseTypesArrayPath() const;

  Q_PROPERTY(DataArrayPath PhaseTypesArrayPath READ getPhaseTypesArrayPath WRITE setPhaseTypesArrayPath)

  /**
   * @brief Setter property for ShapeTypesArrayPath
   */
  void setShapeTypesArrayPath(const DataArrayPath& value);
  /**
   * @brief Getter property for ShapeTypesArrayPath
   * @return Value of ShapeTypesArrayPath
   */
  DataArrayPath getShapeTypesArrayPath() const;

  Q_PROPERTY(DataArrayPath ShapeTypesArrayPath READ getShapeTypesArrayPath WRITE setShapeTypesArrayPath)

  /**
   * @brief Setter property for NumFeaturesArrayPath
   */
  void setNumFeaturesArrayPath(const DataArrayPath& value);
  /**
   * @brief Getter property for NumFeaturesArrayPath
   * @return Value of NumFeaturesArrayPath
   */
  DataArrayPath getNumFeaturesArrayPath() const;

  Q_PROPERTY(DataArrayPath NumFeaturesArrayPath READ getNumFeaturesArrayPath WRITE setNumFeaturesArrayPath)

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
  InsertTransformationPhases();

  /**
   * @brief dataCheck Checks for the appropriate parameter values and availability of arrays
   */
  void dataCheck() override;

  /**
   * @brief Initializes all the private instance variables.
   */
  void initialize();

  /**
   * @brief insertTransformationPhases
   */
  void insertTransformationPhases();

  /**
   * @brief placeTransformationPhase
   * @param curFeature
   * @param sample111
   * @param totalFeatures
   * @param plateThickness
   * @param d
   * @param numFeatures
   * @param euler
   * @return
   */
  bool placeTransformationPhase(int32_t curFeature, float sample111[], int32_t totalFeatures, float plateThickness, float d, float* euler);

  /**
   * @brief peninsulaTransformationPhase
   * @param curFeature
   * @param totalFeatures
   */
  void peninsulaTransformationPhase(int32_t curFeature, int32_t totalFeatures);

  /**
   * @brief transferAttributes
   * @param totalFeatures
   * @param q
   * @param e
   * @param curFeature
   * @return
   */
  size_t transferAttributes(size_t totalFeatures, QuatF q, float e[], int32_t curFeature);

private:
  std::weak_ptr<DataArray<bool>> m_ActivePtr;

  std::weak_ptr<DataArray<int32_t>> m_FeatureIdsPtr;
  int32_t* m_FeatureIds = nullptr;
  std::weak_ptr<DataArray<float>> m_CellEulerAnglesPtr;
  float* m_CellEulerAngles = nullptr;
  std::weak_ptr<DataArray<int32_t>> m_CellPhasesPtr;
  int32_t* m_CellPhases = nullptr;
  std::weak_ptr<DataArray<float>> m_AvgQuatsPtr;
  float* m_AvgQuats = nullptr;
  std::weak_ptr<DataArray<float>> m_CentroidsPtr;
  float* m_Centroids = nullptr;
  std::weak_ptr<DataArray<float>> m_EquivalentDiametersPtr;
  float* m_EquivalentDiameters = nullptr;
  std::weak_ptr<DataArray<float>> m_FeatureEulerAnglesPtr;
  float* m_FeatureEulerAngles = nullptr;
  std::weak_ptr<DataArray<int32_t>> m_FeaturePhasesPtr;
  int32_t* m_FeaturePhases = nullptr;
  std::weak_ptr<DataArray<int32_t>> m_FeatureParentIdsPtr;
  int32_t* m_FeatureParentIds = nullptr;
  std::weak_ptr<DataArray<int32_t>> m_NumFeaturesPerParentPtr;
  int32_t* m_NumFeaturesPerParent = nullptr;
  std::weak_ptr<DataArray<uint32_t>> m_CrystalStructuresPtr;
  uint32_t* m_CrystalStructures = nullptr;
  std::weak_ptr<DataArray<PhaseType::EnumType>> m_PhaseTypesPtr;
  PhaseType::EnumType* m_PhaseTypes = nullptr;
  std::weak_ptr<DataArray<ShapeType::EnumType>> m_ShapeTypesPtr;
  ShapeType::EnumType* m_ShapeTypes = nullptr;
  std::weak_ptr<DataArray<int32_t>> m_NumFeaturesPtr;
  int32_t* m_NumFeatures = nullptr;

  int m_ParentPhase = {};
  int m_TransCrystalStruct = {};
  float m_TransformationPhaseMisorientation = {};
  FloatVec3Type m_TransformationPhaseHabitPlane = {};
  bool m_DefineHabitPlane = {};
  bool m_UseAllVariants = {};
  float m_CoherentFrac = {};
  float m_TransformationPhaseThickness = {};
  int m_NumTransformationPhasesPerFeature = {};
  float m_PeninsulaFrac = {};
  DataArrayPath m_StatsGenCellEnsembleAttributeMatrixPath = {};
  DataArrayPath m_CellFeatureAttributeMatrixName = {};
  DataArrayPath m_FeatureIdsArrayPath = {};
  DataArrayPath m_CellEulerAnglesArrayPath = {};
  DataArrayPath m_CellPhasesArrayPath = {};
  DataArrayPath m_AvgQuatsArrayPath = {};
  DataArrayPath m_CentroidsArrayPath = {};
  DataArrayPath m_EquivalentDiametersArrayPath = {};
  DataArrayPath m_FeatureEulerAnglesArrayPath = {};
  DataArrayPath m_FeaturePhasesArrayPath = {};
  QString m_FeatureParentIdsArrayName = {};
  DataArrayPath m_NumFeaturesPerParentArrayPath = {};
  DataArrayPath m_CrystalStructuresArrayPath = {};
  DataArrayPath m_PhaseTypesArrayPath = {};
  DataArrayPath m_ShapeTypesArrayPath = {};
  DataArrayPath m_NumFeaturesArrayPath = {};

  QVector<LaueOps::Pointer> m_OrientationOps;

  // Cell Data - make sure these are all initialized to nullptr in the constructor

  // Feature Data - make sure these are all initialized to nullptr in the constructor

  // Ensemble Data - make sure these are all initialized to nullptr in the constructor

  FloatVec3Type m_NormalizedTransformationPhaseHabitPlane;

  /**
   * @brief updateFeatureInstancePointers
   */
  void updateFeatureInstancePointers();

  /**
   * @brief updateStatsGenEnsembleInstancePointers
   */
  void updateStatsGenEnsembleInstancePointers();

  /**
   * @brief updateVolEnsembleInstancePointers
   */
  void updateVolEnsembleInstancePointers();

public:
  InsertTransformationPhases(const InsertTransformationPhases&) = delete;            // Copy Constructor Not Implemented
  InsertTransformationPhases(InsertTransformationPhases&&) = delete;                 // Move Constructor Not Implemented
  InsertTransformationPhases& operator=(const InsertTransformationPhases&) = delete; // Copy Assignment Not Implemented
  InsertTransformationPhases& operator=(InsertTransformationPhases&&) = delete;      // Move Assignment Not Implemented
};
