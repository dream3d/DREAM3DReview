/* ============================================================================
 * Copyright (c) 2011 Michael A. Jackson (BlueQuartz Software)
 * Copyright (c) 2011 Dr. Michael A. Groeber (US Air Force Research Laboratories)
 * Copyright (c) 2015 Dr. Joseph C. Tucker (UES, Inc.)
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
#include "TiDwellFatigueCrystallographicAnalysis.h"

#include <QtCore/QTextStream>

#include "SIMPLib/Common/Constants.h"
#include "SIMPLib/DataContainers/DataContainer.h"
#include "SIMPLib/DataContainers/DataContainerArray.h"
#include "SIMPLib/FilterParameters/AbstractFilterParametersReader.h"
#include "SIMPLib/FilterParameters/AttributeMatrixSelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/DataArraySelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/DataContainerSelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/FloatFilterParameter.h"
#include "SIMPLib/FilterParameters/FloatVec3FilterParameter.h"
#include "SIMPLib/FilterParameters/IntFilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedBooleanFilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedPathCreationFilterParameter.h"
#include "SIMPLib/FilterParameters/SeparatorFilterParameter.h"
#include "SIMPLib/Filtering/FilterFactory.hpp"
#include "SIMPLib/Filtering/FilterManager.h"
#include "SIMPLib/Geometry/ImageGeom.h"
#include "SIMPLib/Math/GeometryMath.h"
#include "SIMPLib/Math/MatrixMath.h"
#include "SIMPLib/Math/SIMPLibRandom.h"
#include "SIMPLib/Plugin/ISIMPLibPlugin.h"

#include "EbsdLib/Core/Orientation.hpp"
#include "EbsdLib/Core/OrientationTransformation.hpp"

#include "DREAM3DReview/DREAM3DReviewVersion.h"

/* Create Enumerations to allow the created Attribute Arrays to take part in renaming */
enum createdPathID : RenameDataPath::DataID_t
{
  AttributeMatrixID21 = 21,

  DataArrayID30 = 30,
  DataArrayID31 = 31,
  DataArrayID32 = 32,
  DataArrayID33 = 33,
  DataArrayID34 = 34,
  DataArrayID35 = 35,
  DataArrayID36 = 36,
  DataArrayID37 = 37,
};

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
TiDwellFatigueCrystallographicAnalysis::TiDwellFatigueCrystallographicAnalysis()
{
  m_StressAxis[0] = 0.0f;
  m_StressAxis[1] = 0.0f;
  m_StressAxis[2] = 1.0f;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
TiDwellFatigueCrystallographicAnalysis::~TiDwellFatigueCrystallographicAnalysis() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TiDwellFatigueCrystallographicAnalysis::setupFilterParameters()
{
  FilterParameterVectorType parameters;
  std::vector<QString> linkedProps1 = {"AlphaGlobPhase"};
  parameters.push_back(SIMPL_NEW_LINKED_BOOL_FP("Alpha Glob Phase Present", AlphaGlobPhasePresent, FilterParameter::Category::Parameter, TiDwellFatigueCrystallographicAnalysis, linkedProps1));
  parameters.push_back(SIMPL_NEW_INTEGER_FP("Alpha Glob Phase Number", AlphaGlobPhase, FilterParameter::Category::Parameter, TiDwellFatigueCrystallographicAnalysis));
  parameters.push_back(SIMPL_NEW_INTEGER_FP("Microtextured Region Phase Number", MTRPhase, FilterParameter::Category::Parameter, TiDwellFatigueCrystallographicAnalysis));
  parameters.push_back(SIMPL_NEW_FLOAT_FP("Lattice Parameter A", LatticeParameterA, FilterParameter::Category::Parameter, TiDwellFatigueCrystallographicAnalysis));
  parameters.push_back(SIMPL_NEW_FLOAT_FP("Lattice Parameter C", LatticeParameterC, FilterParameter::Category::Parameter, TiDwellFatigueCrystallographicAnalysis));
  parameters.push_back(SIMPL_NEW_FLOAT_VEC3_FP("Stress Axis", StressAxis, FilterParameter::Category::Parameter, TiDwellFatigueCrystallographicAnalysis));

  parameters.push_back(SIMPL_NEW_INTEGER_FP("Subsurface Feature Distance to Consider (Microns)", SubsurfaceDistance, FilterParameter::Category::Parameter, TiDwellFatigueCrystallographicAnalysis));
  parameters.push_back(SIMPL_NEW_FLOAT_FP("Fraction of Features to Consider", ConsiderationFraction, FilterParameter::Category::Parameter, TiDwellFatigueCrystallographicAnalysis));
  std::vector<QString> linkedProps2;
  linkedProps2.push_back("InitiatorLowerThreshold");
  linkedProps2.push_back("InitiatorUpperThreshold");
  parameters.push_back(
      SIMPL_NEW_LINKED_BOOL_FP("Do Not Assume Initiator Presence", DoNotAssumeInitiatorPresence, FilterParameter::Category::Parameter, TiDwellFatigueCrystallographicAnalysis, linkedProps2));
  parameters.push_back(SIMPL_NEW_FLOAT_FP("Initiator Lower Threshold (Degrees)", InitiatorLowerThreshold, FilterParameter::Category::Parameter, TiDwellFatigueCrystallographicAnalysis));
  parameters.push_back(SIMPL_NEW_FLOAT_FP("Initiator Upper Threshold (Degrees)", InitiatorUpperThreshold, FilterParameter::Category::Parameter, TiDwellFatigueCrystallographicAnalysis));
  parameters.push_back(SIMPL_NEW_FLOAT_FP("Hard Feature Lower Threshold (Degrees)", HardFeatureLowerThreshold, FilterParameter::Category::Parameter, TiDwellFatigueCrystallographicAnalysis));
  parameters.push_back(SIMPL_NEW_FLOAT_FP("Hard Feature Upper Threshold (Degrees)", HardFeatureUpperThreshold, FilterParameter::Category::Parameter, TiDwellFatigueCrystallographicAnalysis));
  parameters.push_back(SIMPL_NEW_FLOAT_FP("Soft Feature Lower Threshold (Degrees)", SoftFeatureLowerThreshold, FilterParameter::Category::Parameter, TiDwellFatigueCrystallographicAnalysis));
  parameters.push_back(SIMPL_NEW_FLOAT_FP("Soft Feature Upper Threshold (Degrees)", SoftFeatureUpperThreshold, FilterParameter::Category::Parameter, TiDwellFatigueCrystallographicAnalysis));

  {
    DataContainerSelectionFilterParameter::RequirementType req;
    parameters.push_back(SIMPL_NEW_DC_SELECTION_FP("Data Container", DataContainerName, FilterParameter::Category::RequiredArray, TiDwellFatigueCrystallographicAnalysis, req));
  }
  parameters.push_back(SeparatorFilterParameter::Create("Cell Data", FilterParameter::Category::RequiredArray));

  {
    DataArraySelectionFilterParameter::RequirementType req;
    parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("FeatureIds", FeatureIdsArrayPath, FilterParameter::Category::RequiredArray, TiDwellFatigueCrystallographicAnalysis, req));
  }

  parameters.push_back(SeparatorFilterParameter::Create("Cell Feature Data", FilterParameter::Category::RequiredArray));
  {
    AttributeMatrixSelectionFilterParameter::RequirementType req = AttributeMatrixSelectionFilterParameter::CreateRequirement(AttributeMatrix::Type::CellFeature, IGeometry::Type::Unknown);
    parameters.push_back(
        SIMPL_NEW_AM_SELECTION_FP("Cell Feature Attribute Matrix", CellFeatureAttributeMatrixPath, FilterParameter::Category::RequiredArray, TiDwellFatigueCrystallographicAnalysis, req));
  }
  {
    DataArraySelectionFilterParameter::RequirementType req = DataArraySelectionFilterParameter::CreateCategoryRequirement(SIMPL::Defaults::AnyPrimitive, 3, AttributeMatrix::Category::Feature);
    parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Average Euler Angles", FeatureEulerAnglesArrayPath, FilterParameter::Category::RequiredArray, TiDwellFatigueCrystallographicAnalysis, req));
  }
  {
    DataArraySelectionFilterParameter::RequirementType req = DataArraySelectionFilterParameter::CreateCategoryRequirement(SIMPL::Defaults::AnyPrimitive, 1, AttributeMatrix::Category::Feature);
    parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Phases", FeaturePhasesArrayPath, FilterParameter::Category::RequiredArray, TiDwellFatigueCrystallographicAnalysis, req));
  }
  {
    DataArraySelectionFilterParameter::RequirementType req = DataArraySelectionFilterParameter::CreateCategoryRequirement(SIMPL::Defaults::AnyPrimitive, 1, AttributeMatrix::Category::Feature);
    parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Neighbor List", NeighborListArrayPath, FilterParameter::Category::RequiredArray, TiDwellFatigueCrystallographicAnalysis, req));
  }
  {
    DataArraySelectionFilterParameter::RequirementType req = DataArraySelectionFilterParameter::CreateCategoryRequirement(SIMPL::Defaults::AnyPrimitive, 3, AttributeMatrix::Category::Feature);
    parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Centroids", CentroidsArrayPath, FilterParameter::Category::RequiredArray, TiDwellFatigueCrystallographicAnalysis, req));
  }

  parameters.push_back(SeparatorFilterParameter::Create("Cell Ensemble Data", FilterParameter::Category::RequiredArray));
  {
    DataArraySelectionFilterParameter::RequirementType req = DataArraySelectionFilterParameter::CreateCategoryRequirement(SIMPL::Defaults::AnyPrimitive, 1, AttributeMatrix::Category::Ensemble);
    parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Crystal Structures", CrystalStructuresArrayPath, FilterParameter::Category::RequiredArray, TiDwellFatigueCrystallographicAnalysis, req));
  }

  parameters.push_back(SeparatorFilterParameter::Create("Cell Data", FilterParameter::Category::CreatedArray));
  parameters.push_back(
      SIMPL_NEW_DA_WITH_LINKED_AM_FP("Parent Ids", CellParentIdsArrayName, FeatureIdsArrayPath, FeatureIdsArrayPath, FilterParameter::Category::CreatedArray, TiDwellFatigueCrystallographicAnalysis));

  parameters.push_back(SeparatorFilterParameter::Create("Cell Feature Data", FilterParameter::Category::CreatedArray));
  parameters.push_back(SIMPL_NEW_AM_WITH_LINKED_DC_FP("Cell Feature Attribute Matrix", NewCellFeatureAttributeMatrixName, FeatureIdsArrayPath, FilterParameter::Category::CreatedArray,
                                                      TiDwellFatigueCrystallographicAnalysis));
  parameters.push_back(SIMPL_NEW_DA_WITH_LINKED_AM_FP("Selected Features", SelectedFeaturesArrayName, CellFeatureAttributeMatrixPath, CellFeatureAttributeMatrixPath,
                                                      FilterParameter::Category::CreatedArray, TiDwellFatigueCrystallographicAnalysis));
  parameters.push_back(SIMPL_NEW_DA_WITH_LINKED_AM_FP("Initiators", InitiatorsArrayName, CellFeatureAttributeMatrixPath, CellFeatureAttributeMatrixPath, FilterParameter::Category::CreatedArray,
                                                      TiDwellFatigueCrystallographicAnalysis));
  parameters.push_back(SIMPL_NEW_DA_WITH_LINKED_AM_FP("Hard Features", HardFeaturesArrayName, CellFeatureAttributeMatrixPath, CellFeatureAttributeMatrixPath, FilterParameter::Category::CreatedArray,
                                                      TiDwellFatigueCrystallographicAnalysis));
  parameters.push_back(SIMPL_NEW_DA_WITH_LINKED_AM_FP("Soft Features", SoftFeaturesArrayName, CellFeatureAttributeMatrixPath, CellFeatureAttributeMatrixPath, FilterParameter::Category::CreatedArray,
                                                      TiDwellFatigueCrystallographicAnalysis));
  parameters.push_back(SIMPL_NEW_DA_WITH_LINKED_AM_FP("Hard-Soft Groups", HardSoftGroupsArrayName, CellFeatureAttributeMatrixPath, CellFeatureAttributeMatrixPath,
                                                      FilterParameter::Category::CreatedArray, TiDwellFatigueCrystallographicAnalysis));
  parameters.push_back(SIMPL_NEW_DA_WITH_LINKED_AM_FP("Parent Ids", FeatureParentIdsArrayName, CellFeatureAttributeMatrixPath, CellFeatureAttributeMatrixPath, FilterParameter::Category::CreatedArray,
                                                      TiDwellFatigueCrystallographicAnalysis));

  // parameters.push_back(SIMPL_NEW_STRING_FP("Active", ActiveArrayName, FilterParameter::Category::CreatedArray, TiDwellFatigueCrystallographicAnalysis));
  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TiDwellFatigueCrystallographicAnalysis::readFilterParameters(AbstractFilterParametersReader* reader, int index)
{
  reader->openFilterGroup(this, index);
  setDataContainerName(reader->readDataArrayPath("DataContainerName", getDataContainerName()));
  setAlphaGlobPhasePresent(reader->readValue("AlphaGlobPhasePresent", getAlphaGlobPhasePresent()));
  setAlphaGlobPhase(reader->readValue("AlphaGlobPhase", getAlphaGlobPhase()));
  setMTRPhase(reader->readValue("MTRPhase", getMTRPhase()));
  setLatticeParameterA(reader->readValue("LatticeParameterA", getLatticeParameterA()));
  setLatticeParameterC(reader->readValue("LatticeParameterC", getLatticeParameterC()));
  setStressAxis(reader->readFloatVec3("Stress Axis", getStressAxis()));
  setSubsurfaceDistance(reader->readValue("SubsurfaceDistance", getSubsurfaceDistance()));
  setConsiderationFraction(reader->readValue("ConsiderationFraction", getConsiderationFraction()));
  setDoNotAssumeInitiatorPresence(reader->readValue("DoNotAssumeInitiatorPresence", getDoNotAssumeInitiatorPresence()));
  setInitiatorLowerThreshold(reader->readValue("InitiatorLowerThreshold", getInitiatorLowerThreshold()));
  setInitiatorUpperThreshold(reader->readValue("InitiatorUpperThreshold", getInitiatorUpperThreshold()));
  setHardFeatureLowerThreshold(reader->readValue("HardFeatureLowerThreshold", getHardFeatureLowerThreshold()));
  setHardFeatureUpperThreshold(reader->readValue("HardFeatureUpperThreshold", getHardFeatureUpperThreshold()));
  setSoftFeatureLowerThreshold(reader->readValue("SoftFeatureLowerThreshold", getSoftFeatureLowerThreshold()));
  setSoftFeatureUpperThreshold(reader->readValue("SoftFeatureUpperThreshold", getSoftFeatureUpperThreshold()));
  setNewCellFeatureAttributeMatrixName(reader->readString("NewCellFeatureAttributeMatrixName", getNewCellFeatureAttributeMatrixName()));
  setSelectedFeaturesArrayName(reader->readString("SelectedFeaturesArrayName", getSelectedFeaturesArrayName()));
  setInitiatorsArrayName(reader->readString("InitiatorsArrayName", getInitiatorsArrayName()));
  setHardFeaturesArrayName(reader->readString("HardFeaturesArrayName", getHardFeaturesArrayName()));
  setSoftFeaturesArrayName(reader->readString("SoftFeaturesArrayName", getSoftFeaturesArrayName()));
  setHardSoftGroupsArrayName(reader->readString("HardSoftGroupsArrayName", getHardSoftGroupsArrayName()));
  setCellFeatureAttributeMatrixName(reader->readString("CellFeatureAttributeMatrixName", getCellFeatureAttributeMatrixName()));
  setCellFeatureAttributeMatrixPath(reader->readDataArrayPath("CellFeatureAttributeMatrixPath", getCellFeatureAttributeMatrixPath()));
  setActiveArrayName(reader->readString("ActiveArrayName", getActiveArrayName()));
  setFeatureIdsArrayPath(reader->readDataArrayPath("FeatureIdsArrayPath", getFeatureIdsArrayPath()));
  setFeatureParentIdsArrayName(reader->readString("FeatureParentIdsArrayName", getFeatureParentIdsArrayName()));
  setCellParentIdsArrayName(reader->readString("CellParentIdsArrayName", getCellParentIdsArrayName()));
  setFeatureEulerAnglesArrayPath(reader->readDataArrayPath("FeatureEulerAnglesArrayPath", getFeatureEulerAnglesArrayPath()));
  setFeaturePhasesArrayPath(reader->readDataArrayPath("FeaturePhasesArrayPath", getFeaturePhasesArrayPath()));
  setNeighborListArrayPath(reader->readDataArrayPath("NeighborListArrayPath", getNeighborListArrayPath()));
  setCentroidsArrayPath(reader->readDataArrayPath("CentroidsArrayPath", getCentroidsArrayPath()));
  setCrystalStructuresArrayPath(reader->readDataArrayPath("CrystalStructuresArrayPath", getCrystalStructuresArrayPath()));
  reader->closeFilterGroup();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TiDwellFatigueCrystallographicAnalysis::readFilterParameters(QJsonObject& obj)
{
  AbstractFilter::readFilterParameters(obj);
  setCellFeatureAttributeMatrixName(obj["CellFeatureAttributeMatrixName"].toString());
  setActiveArrayName(obj["ActiveArrayName"].toString());
}

// FP: Check why these values are not connected to a filter parameter!

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TiDwellFatigueCrystallographicAnalysis::writeFilterParameters(QJsonObject& obj) const
{
  AbstractFilter::writeFilterParameters(obj);
  obj["CellFeatureAttributeMatrixName"] = getCellFeatureAttributeMatrixName();
  obj["ActiveArrayName"] = getActiveArrayName();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TiDwellFatigueCrystallographicAnalysis::updateFeatureInstancePointers()
{
  clearErrorCode();
  clearWarningCode();

  if(nullptr != m_ActivePtr.lock())
  {
    m_Active = m_ActivePtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TiDwellFatigueCrystallographicAnalysis::initialize()
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TiDwellFatigueCrystallographicAnalysis::dataCheck()
{
  clearErrorCode();
  clearWarningCode();
  DataArrayPath tempPath;

  DataContainer::Pointer m = getDataContainerArray()->getPrereqDataContainer(this, m_FeatureIdsArrayPath.getDataContainerName(), false);
  if(getErrorCode() < 0 || nullptr == m.get())
  {
    return;
  }
  std::vector<size_t> tDims(1, 0);
  AttributeMatrix::Pointer newCellFeatureAttrMat = m->createNonPrereqAttributeMatrix(this, getNewCellFeatureAttributeMatrixName(), tDims, AttributeMatrix::Type::CellFeature, AttributeMatrixID21);
  if(getErrorCode() < 0)
  {
    return;
  }

  ImageGeom::Pointer image = m->getPrereqGeometry<ImageGeom>(this);
  if(getErrorCode() < 0 || nullptr == image.get())
  {
    return;
  }

  // Feature Data
  std::vector<size_t> dims(1, 1);
  m_FeatureIdsPtr = getDataContainerArray()->getPrereqArrayFromPath<DataArray<int32_t>>(this, getFeatureIdsArrayPath(), dims);
  if(nullptr != m_FeatureIdsPtr.lock())
  {
    m_FeatureIds = m_FeatureIdsPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */

  tempPath.update(getCellFeatureAttributeMatrixPath().getDataContainerName(), getCellFeatureAttributeMatrixPath().getAttributeMatrixName(), getSelectedFeaturesArrayName());
  m_SelectedFeaturesPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<bool>>(this, tempPath, false, dims, "", DataArrayID31);
  if(nullptr != m_SelectedFeaturesPtr.lock())
  {
    m_SelectedFeatures = m_SelectedFeaturesPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */

  tempPath.update(getCellFeatureAttributeMatrixPath().getDataContainerName(), getCellFeatureAttributeMatrixPath().getAttributeMatrixName(), getInitiatorsArrayName());
  m_InitiatorsPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<bool>>(this, tempPath, false, dims, "", DataArrayID32);
  if(nullptr != m_InitiatorsPtr.lock())
  {
    m_Initiators = m_InitiatorsPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */

  tempPath.update(getCellFeatureAttributeMatrixPath().getDataContainerName(), getCellFeatureAttributeMatrixPath().getAttributeMatrixName(), getHardFeaturesArrayName());
  m_HardFeaturesPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<bool>>(this, tempPath, false, dims, "", DataArrayID33);
  if(nullptr != m_HardFeaturesPtr.lock())
  {
    m_HardFeatures = m_HardFeaturesPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */

  tempPath.update(getCellFeatureAttributeMatrixPath().getDataContainerName(), getCellFeatureAttributeMatrixPath().getAttributeMatrixName(), getSoftFeaturesArrayName());
  m_SoftFeaturesPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<bool>>(this, tempPath, false, dims, "", DataArrayID34);
  if(nullptr != m_SoftFeaturesPtr.lock())
  {
    m_SoftFeatures = m_SoftFeaturesPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */

  tempPath.update(getCellFeatureAttributeMatrixPath().getDataContainerName(), getCellFeatureAttributeMatrixPath().getAttributeMatrixName(), getHardSoftGroupsArrayName());
  m_HardSoftGroupsPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<bool>>(this, tempPath, false, dims, "", DataArrayID35);
  if(nullptr != m_HardSoftGroupsPtr.lock())
  {
    m_HardSoftGroups = m_HardSoftGroupsPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */

  tempPath.update(getFeatureIdsArrayPath().getDataContainerName(), getFeatureIdsArrayPath().getAttributeMatrixName(), getCellParentIdsArrayName());
  m_CellParentIdsPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<int32_t>>(this, tempPath, -1, dims, "", DataArrayID36);
  if(nullptr != m_CellParentIdsPtr.lock())
  {
    m_CellParentIds = m_CellParentIdsPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */

  tempPath.update(getCellFeatureAttributeMatrixPath().getDataContainerName(), getCellFeatureAttributeMatrixPath().getAttributeMatrixName(), getFeatureParentIdsArrayName());
  m_FeatureParentIdsPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<int32_t>>(this, tempPath, -1, dims, "", DataArrayID37);
  if(nullptr != m_FeatureParentIdsPtr.lock())
  {
    m_FeatureParentIds = m_FeatureParentIdsPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */

  dims[0] = 3;
  m_FeatureEulerAnglesPtr = getDataContainerArray()->getPrereqArrayFromPath<DataArray<float>>(this, getFeatureEulerAnglesArrayPath(), dims);
  if(nullptr != m_FeatureEulerAnglesPtr.lock())
  {
    m_FeatureEulerAngles = m_FeatureEulerAnglesPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */

  dims[0] = 1;
  m_FeaturePhasesPtr = getDataContainerArray()->getPrereqArrayFromPath<DataArray<int32_t>>(this, getFeaturePhasesArrayPath(), dims);
  if(nullptr != m_FeaturePhasesPtr.lock())
  {
    m_FeaturePhases = m_FeaturePhasesPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */

  m_NeighborList = getDataContainerArray()->getPrereqArrayFromPath<NeighborList<int>>(this, getNeighborListArrayPath(), dims);

  dims[0] = 3;
  m_CentroidsPtr = getDataContainerArray()->getPrereqArrayFromPath<DataArray<float>>(this, getCentroidsArrayPath(), dims);
  if(nullptr != m_CentroidsPtr.lock())
  {
    m_Centroids = m_CentroidsPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */

  // New Feature Data

  // Ensemble Data
  dims[0] = 1;

  m_CrystalStructuresPtr = getDataContainerArray()->getPrereqArrayFromPath<DataArray<unsigned int>>(this, getCrystalStructuresArrayPath(), dims);
  if(nullptr != m_CrystalStructuresPtr.lock())
  {
    m_CrystalStructures = m_CrystalStructuresPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TiDwellFatigueCrystallographicAnalysis::execute()
{
  clearErrorCode();
  clearWarningCode();

  dataCheck();
  if(getErrorCode() < 0)
  {
    return;
  }

  SIMPL_RANDOMNG_NEW()

  size_t totalFeatures = m_FeaturePhasesPtr.lock()->getNumberOfTuples();
  size_t totalPoints = static_cast<size_t>(m_FeatureIdsPtr.lock()->getNumberOfTuples());

  bool subsurfaceFlag = false;
  bool hardfeatureFlag = false;
  //  bool initiatorFlag = false;
  //  bool softfeatureFlag = false;

  float random = 0.0f;

  // Normalize input stress axis
  MatrixMath::Normalize3x1(m_StressAxis[0], m_StressAxis[1], m_StressAxis[2]);

  for(size_t i = 1; i < totalFeatures; ++i)
  {
    //  int stop = 0;
    random = static_cast<float>(rg.genrand_res53());
    if(random < m_ConsiderationFraction)
    {
      m_SelectedFeatures[i] = true;
      subsurfaceFlag = determine_subsurfacefeatures(i);
    }
    if(subsurfaceFlag)
    {
      // Determine if it's a hard feature
      if(m_FeaturePhases[i] == m_MTRPhase)
      {
        hardfeatureFlag = determine_hardfeatures(i);
      }
      // Determine if it's a soft feature only if it's not a hard feature
      if(!hardfeatureFlag && m_FeaturePhases[i] == m_MTRPhase)
      {
        determine_softfeatures(i);
      }
      // Determine if it's an initiator only if we're assuming initiators are not necessarily present
      if(m_DoNotAssumeInitiatorPresence && m_FeaturePhases[i] == m_AlphaGlobPhase)
      {
        determine_initiators(i);
      }
      subsurfaceFlag = false;
    }
  }

  for(size_t i = 1; i < totalFeatures; ++i)
  {
    // Group neighboring hard features and soft features
    if(m_HardFeatures[i] || m_SoftFeatures[i])
    {
      group_flaggedfeatures(i);
    }
  }

  size_t fid = 0;

  // map grouped features to the cells
  for(size_t i = 0; i < totalPoints; ++i)
  {
    fid = m_FeatureParentIds[m_FeatureIds[i]];
    if(fid != -1)
    {
      m_CellParentIds[i] = fid;
    }
    else
    {
      m_CellParentIds[i] = m_FeatureIds[i];
    }
  }

  QString filtName = "FindNeighbors";
  FilterManager* fm = FilterManager::Instance();
  IFilterFactory::Pointer findNeighborFactory = fm->getFactoryFromClassName(filtName);

  DataArrayPath tempPath;

  if(nullptr != findNeighborFactory.get())
  {
    // If we get this far, the Factory is good so creating the filter should not fail unless something has
    // horribly gone wrong in which case the system is going to come down quickly after this.
    AbstractFilter::Pointer find_Neighbor = findNeighborFactory->create();

    // Connect up the Error/Warning/Progress object so the filter can report those things
    connect(find_Neighbor.get(), SIGNAL(messageGenerated(const AbstractMessage::Pointer&)), this, SIGNAL(messageGenerated(const AbstractMessage::Pointer&)));
    find_Neighbor->setDataContainerArray(getDataContainerArray()); // AbstractFilter implements this so no problem
    // Now set the filter parameters for the filter using QProperty System since we can not directly
    // instantiate the filter since it resides in a plugin. These calls are SLOW. DO NOT EVER do this in a
    // tight loop. Your filter will slow down by 10X.
    tempPath.update(getDataContainerName().getDataContainerName(), getCellFeatureAttributeMatrixName(), "");
    QVariant v;
    v.setValue(tempPath);
    bool propWasSet = find_Neighbor->setProperty("CellFeatureAttributeMatrixPath", v);
    if(!propWasSet)
    {
      QString ss = QObject::tr("Ti Dwell Fatigue Error Setting Property '%1' into filter '%2' which is a subfilter called by Ti Dwell Fatigue Error. The property was not set which could mean the "
                               "property was not exposed with a Q_PROPERTY macro. Please notify the developers.")
                       .arg("CellFeatureAttributeMatrixPath")
                       .arg(filtName);
      setErrorCondition(-109872, ss);
    }

    tempPath.update(getFeatureIdsArrayPath().getDataContainerName(), getFeatureIdsArrayPath().getAttributeMatrixName(), getCellParentIdsArrayName());
    v.setValue(tempPath);
    propWasSet = find_Neighbor->setProperty("FeatureIdsArrayPath", v);
    if(!propWasSet)
    {
      QString ss = QObject::tr("Ti Dwell Fatigue Error Setting Property '%1' into filter '%2' which is a subfilter called by Ti Dwell Fatigue Error. The property was not set which could mean the "
                               "property was not exposed with a Q_PROPERTY macro. Please notify the developers.")
                       .arg("FeatureIdsArrayPath")
                       .arg(filtName);
      setErrorCondition(-109873, ss);
    }

    QString tempString("BoundaryCells");
    v.setValue(tempString);
    propWasSet = find_Neighbor->setProperty("BoundaryCellsArrayName", v);
    if(!propWasSet)
    {
      QString ss = QObject::tr("Ti Dwell Fatigue Error Setting Property '%1' into filter '%2' which is a subfilter called by Ti Dwell Fatigue Error. The property was not set which could mean the "
                               "property was not exposed with a Q_PROPERTY macro. Please notify the developers.")
                       .arg("BoundaryCellsArrayName")
                       .arg(filtName);
      setErrorCondition(-109874, ss);
    }

    tempString = "SurfaceFeatures";
    v.setValue(tempString);
    propWasSet = find_Neighbor->setProperty("SurfaceFeaturesArrayName", v);
    if(!propWasSet)
    {
      QString ss = QObject::tr("Ti Dwell Fatigue Error Setting Property '%1' into filter '%2' which is a subfilter called by Ti Dwell Fatigue Error. The property was not set which could mean the "
                               "property was not exposed with a Q_PROPERTY macro. Please notify the developers.")
                       .arg("SurfaceFeaturesArrayName")
                       .arg(filtName);
      setErrorCondition(-109875, ss);
    }

    v.setValue(false);
    propWasSet = find_Neighbor->setProperty("StoreBoundaryCells", v);
    if(!propWasSet)
    {
      QString ss = QObject::tr("Ti Dwell Fatigue Error Setting Property '%1' into filter '%2' which is a subfilter called by Ti Dwell Fatigue Error. The property was not set which could mean the "
                               "property was not exposed with a Q_PROPERTY macro. Please notify the developers.")
                       .arg("StoreBoundaryCellsArrayName")
                       .arg(filtName);
      setErrorCondition(-109876, ss);
    }

    v.setValue(false);
    propWasSet = find_Neighbor->setProperty("StoreSurfaceFeatures", v);
    if(!propWasSet)
    {
      QString ss = QObject::tr("Ti Dwell Fatigue Error Setting Property '%1' into filter '%2' which is a subfilter called by Ti Dwell Fatigue Error. The property was not set which could mean the "
                               "property was not exposed with a Q_PROPERTY macro. Please notify the developers.")
                       .arg("StoreSurfaceFeaturesArrayName")
                       .arg(filtName);
      setErrorCondition(-109877, ss);
    }

    tempString = "ParentNumNeighbors";
    v.setValue(tempString);
    propWasSet = find_Neighbor->setProperty("NumNeighborsArrayName", v);
    if(!propWasSet)
    {
      QString ss = QObject::tr("Ti Dwell Fatigue Error Setting Property '%1' into filter '%2' which is a subfilter called by Ti Dwell Fatigue Error. The property was not set which could mean the "
                               "property was not exposed with a Q_PROPERTY macro. Please notify the developers.")
                       .arg("ParentNumNeighbors")
                       .arg(filtName);
      setErrorCondition(-109878, ss);
    }

    tempString = "ParentNeighborList";
    v.setValue(tempString);
    propWasSet = find_Neighbor->setProperty("NeighborListArrayName", v);
    if(!propWasSet)
    {
      QString ss = QObject::tr("Ti Dwell Fatigue Error Setting Property '%1' into filter '%2' which is a subfilter called by Ti Dwell Fatigue Error. The property was not set which could mean the "
                               "property was not exposed with a Q_PROPERTY macro. Please notify the developers.")
                       .arg("ParentNeighborList")
                       .arg(filtName);
      setErrorCondition(-109879, ss);
    }

    tempString = "ParentSharedSurfaceAreaList";
    v.setValue(tempString);
    propWasSet = find_Neighbor->setProperty("SharedSurfaceAreaListArrayName", v);
    if(!propWasSet)
    {
      QString ss = QObject::tr("Ti Dwell Fatigue Error Setting Property '%1' into filter '%2' which is a subfilter called by Ti Dwell Fatigue Error. The property was not set which could mean the "
                               "property was not exposed with a Q_PROPERTY macro. Please notify the developers.")
                       .arg("ParentSharedSurfaceAreaList")
                       .arg(filtName);
      setErrorCondition(-109880, ss);
    }
    find_Neighbor->execute();
  }

  std::vector<size_t> cDims(1, 1);
  tempPath.update(getFeatureEulerAnglesArrayPath().getDataContainerName(), getFeatureEulerAnglesArrayPath().getAttributeMatrixName(), "ParentNeighborList");
  m_ParentNeighborList = getDataContainerArray()->getPrereqArrayFromPath<NeighborList<int>>(this, tempPath, cDims);

  for(size_t i = 1; i < totalFeatures; ++i)
  {
    // only proceed if it's either an MTR or alpha glob
    if(m_FeaturePhases[i] == m_MTRPhase || (m_FeaturePhases[i] == m_AlphaGlobPhase && m_DoNotAssumeInitiatorPresence))
    {
      // Determine if it's a hard-soft pair only if there the current ID is either a hardfeature, initiator or soft feature
      if(m_Initiators[i] || m_HardFeatures[i] || m_SoftFeatures[i])
      {
        assign_HardSoftGroups(i);
      }
    }
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
bool TiDwellFatigueCrystallographicAnalysis::determine_subsurfacefeatures(int index)
{
  // using feature euler angles simply because it's available
  DataContainer::Pointer m = getDataContainerArray()->getDataContainer(m_FeatureEulerAnglesArrayPath.getDataContainerName());
  FloatVec3Type origin = m->getGeometryAs<ImageGeom>()->getOrigin();

  bool subsurfaceFlag = false;

  SizeVec3Type dims = m->getGeometryAs<ImageGeom>()->getDimensions();
  FloatVec3Type spacing = m->getGeometryAs<ImageGeom>()->getSpacing();
  FloatVec3Type xyzScaledDimension(origin[0] + dims[0] * spacing[0], origin[1] + dims[1] * spacing[1], origin[2] + dims[2] * spacing[2]);

  // check if current feature centroid is within the subsurface defined centroid
  if(m_Centroids[3 * index + 0] >= m_SubsurfaceDistance && m_Centroids[3 * index + 0] <= (xyzScaledDimension[0] - m_SubsurfaceDistance) && m_Centroids[3 * index + 1] >= m_SubsurfaceDistance &&
     m_Centroids[3 * index + 1] <= (xyzScaledDimension[1] - m_SubsurfaceDistance) && m_Centroids[3 * index + 2] >= m_SubsurfaceDistance &&
     m_Centroids[3 * index + 2] <= (xyzScaledDimension[2] - m_SubsurfaceDistance))
  {
    subsurfaceFlag = true;
  }
  return subsurfaceFlag;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
bool TiDwellFatigueCrystallographicAnalysis::determine_hardfeatures(int index)
{
  float g[3][3] = {{0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f}};
  float w = 0.0f;
  const int hardfeaturePlane[12][4] = {{-1, 0, 1, 7},  {-1, 1, 0, -7}, {0, -1, 1, 7},  {1, 0, -1, 7},  {1, -1, 0, -7}, {0, 1, -1, 7},
                                       {0, -1, 1, -7}, {0, 1, -1, -7}, {1, 0, -1, -7}, {-1, 0, 1, -7}, {-1, 1, 0, 7},  {1, -1, 0, 7}};
  float hardfeaturePlaneNormal[12][3];
  for(size_t i = 0; i < 12; i++)
  {
    for(size_t j = 0; j < 3; j++)
    {
      hardfeaturePlaneNormal[i][j] = 0.0f;
    }
  }
  const float m_OneOverA = 1 / m_LatticeParameterA;
  const float m_OneOverAxSqrtThree = 1 / (m_LatticeParameterA * sqrtf(3.0f));
  const float m_OneOverC = 1 / m_LatticeParameterC;
  bool hardfeatureFlag = false;

  // convert Miller-Bravais to unit normal
  for(int j = 0; j < 12; ++j)
  {
    hardfeaturePlaneNormal[j][0] = hardfeaturePlane[j][0] * m_OneOverA;
    hardfeaturePlaneNormal[j][1] = (2.0f * hardfeaturePlane[j][1] + hardfeaturePlane[j][0]) * m_OneOverAxSqrtThree;
    hardfeaturePlaneNormal[j][2] = hardfeaturePlane[j][3] * m_OneOverC;
  }

  OrientationTransformation::eu2om<OrientationF, OrientationF>(OrientationF(m_FeatureEulerAngles + 3 * index, 3)).toGMatrix(g);
  if(m_FeaturePhases[index] == m_MTRPhase)
  {
    for(int j = 0; j < 12; ++j)
    {
      w = find_angle(g, hardfeaturePlaneNormal[j][0], hardfeaturePlaneNormal[j][1], hardfeaturePlaneNormal[j][2]);
      if(w >= m_HardFeatureLowerThreshold && w <= m_HardFeatureUpperThreshold)
      {
        m_HardFeatures[index] = true;
        hardfeatureFlag = true;
        break;
      }
    }
  }
  return hardfeatureFlag;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TiDwellFatigueCrystallographicAnalysis::determine_initiators(int index)
{
  const float caxis[3] = {0.0f, 0.0f, 1.0f};
  bool initiatorFlag = false;
  float w = 0.0f;
  float g[3][3] = {{0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f}};

  if(m_FeaturePhases[index] == m_AlphaGlobPhase)
  {
    OrientationTransformation::eu2om<OrientationF, OrientationF>(OrientationF(m_FeatureEulerAngles + 3 * index, 3)).toGMatrix(g);

    w = find_angle(g, caxis[0], caxis[1], caxis[2]);
    if(w >= m_InitiatorLowerThreshold && w <= m_InitiatorUpperThreshold)
    {
      m_Initiators[index] = true;
      initiatorFlag = true;
    }
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TiDwellFatigueCrystallographicAnalysis::determine_softfeatures(int index)
{
  const float caxis[3] = {0.0f, 0.0f, 1.0f};
  float w = 0.0f;
  float g[3][3] = {{0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f}};
  // bool softfeatureFlag = false;

  if(m_FeaturePhases[index] == m_MTRPhase)
  {
    OrientationTransformation::eu2om<OrientationF, OrientationF>(OrientationF(m_FeatureEulerAngles + 3 * index, 3)).toGMatrix(g);

    w = find_angle(g, caxis[0], caxis[1], caxis[2]);
    if(w >= m_SoftFeatureLowerThreshold && w <= m_SoftFeatureUpperThreshold)
    {
      m_SoftFeatures[index] = true;
    }
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TiDwellFatigueCrystallographicAnalysis::group_flaggedfeatures(int index)
{
  // But since a pointer is difficult to use operators with we will now create a
  // reference variable to the pointer with the correct variable name that allows
  // us to use the same syntax as the "vector of vectors"
  NeighborList<int>& neighborlist = *(m_NeighborList.lock());

  for(int j = 0; j < neighborlist[index].size(); ++j)
  {
    if(m_FeatureParentIds[index] != neighborlist[index][j] && ((m_HardFeatures[index] && m_HardFeatures[neighborlist[index][j]]) || (m_SoftFeatures[index] && m_SoftFeatures[neighborlist[index][j]])))
    {
      if(m_FeatureParentIds[neighborlist[index][j]] == -1)
      {
        if(m_FeatureParentIds[index] == -1)
        {
          m_FeatureParentIds[neighborlist[index][j]] = index;
        }
        if(m_FeatureParentIds[index] != -1)
        {
          m_FeatureParentIds[neighborlist[index][j]] = m_FeatureParentIds[index];
        }
        group_flaggedfeatures(neighborlist[index][j]);
      }
    }
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TiDwellFatigueCrystallographicAnalysis::assign_HardSoftGroups(int index)
{
  // But since a pointer is difficult to use operators with we will now create a
  // reference variable to the pointer with the correct variable name that allows
  // us to use the same syntax as the "vector of vectors"
  NeighborList<int>& parentneighborlist = *(m_ParentNeighborList.lock());

  bool hardfeatureFlag = false;
  bool initiatorFlag = false;
  bool softfeatureFlag = false;
  int hardfeatureIndex = 0;
  int softfeatureIndex = 0;

  for(int j = 0; j < parentneighborlist[index].size(); ++j)
  {
    if((m_DoNotAssumeInitiatorPresence && m_FeaturePhases[index] == m_AlphaGlobPhase) || m_FeaturePhases[parentneighborlist[index][j]] == m_MTRPhase)
    {
      if(m_Initiators[index] && !initiatorFlag && m_DoNotAssumeInitiatorPresence && m_FeaturePhases[index] == m_AlphaGlobPhase)
      {
        initiatorFlag = true;
      }
      if(m_Initiators[parentneighborlist[index][j]] && !initiatorFlag && m_DoNotAssumeInitiatorPresence && m_FeaturePhases[parentneighborlist[index][j]] == m_AlphaGlobPhase)
      {
        initiatorFlag = true;
      }
      if(m_HardFeatures[index] && !hardfeatureFlag && m_FeaturePhases[index] == m_MTRPhase)
      {
        hardfeatureFlag = true;
        hardfeatureIndex = index;
      }
      if(m_HardFeatures[parentneighborlist[index][j]] && !hardfeatureFlag && m_FeaturePhases[parentneighborlist[index][j]] == m_MTRPhase)
      {
        hardfeatureFlag = true;
        hardfeatureIndex = parentneighborlist[index][j];
      }
      if(m_SoftFeatures[index] && !softfeatureFlag && m_FeaturePhases[index] == m_MTRPhase)
      {
        softfeatureFlag = true;
        softfeatureIndex = index;
      }
      if(m_SoftFeatures[parentneighborlist[index][j]] && !softfeatureFlag && m_FeaturePhases[parentneighborlist[index][j]] == m_MTRPhase)
      {
        softfeatureFlag = true;
        softfeatureIndex = parentneighborlist[index][j];
      }
      // only flag as a hard-soft pair if there's a neighboring group of initiator - hardfeature - soft feature and either the current index of the hardfeature or
      // soft feature have not already been flagged as a bad acting pair
      if((!m_DoNotAssumeInitiatorPresence || initiatorFlag) && hardfeatureFlag && softfeatureFlag && (!m_HardSoftGroups[hardfeatureIndex] || !m_HardSoftGroups[softfeatureIndex]))
      {
        m_HardSoftGroups[hardfeatureIndex] = true;
        m_HardSoftGroups[softfeatureIndex] = true;
        break;
      }
    }
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
float TiDwellFatigueCrystallographicAnalysis::find_angle(float g[3][3], float planeNormalU, float planeNormalV, float planeNormalW)
{
  float gt[3][3] = {{0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f}};
  float v[3] = {0};
  float sampleLoading[3] = {m_StressAxis[0], m_StressAxis[1], m_StressAxis[2]};
  float w = 0.0f;

  MatrixMath::Normalize3x1(planeNormalU, planeNormalV, planeNormalW);
  float planeNormal[3] = {planeNormalU, planeNormalV, planeNormalW};
  MatrixMath::Transpose3x3(g, gt);
  MatrixMath::Multiply3x3with3x1(gt, planeNormal, v);
  // Normalize so that the magnitude is 1
  MatrixMath::Normalize3x1(v);
  if(v[2] < 0)
  {
    MatrixMath::Multiply3x1withConstant(v, -1.0f);
  }
  w = GeometryMath::CosThetaBetweenVectors(v, sampleLoading);
  w = acos(w);
  // Convert from radian to degrees
  return w *= SIMPLib::Constants::k_180OverPiD;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AbstractFilter::Pointer TiDwellFatigueCrystallographicAnalysis::newFilterInstance(bool copyFilterParameters) const
{
  TiDwellFatigueCrystallographicAnalysis::Pointer filter = TiDwellFatigueCrystallographicAnalysis::New();
  if(copyFilterParameters)
  {
    copyFilterParameterInstanceVariables(filter.get());
  }
  return filter;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString TiDwellFatigueCrystallographicAnalysis::getCompiledLibraryName() const
{
  return DREAM3DReviewConstants::DREAM3DReviewBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString TiDwellFatigueCrystallographicAnalysis::getBrandingString() const
{
  return TransformationPhaseConstants::TransformationPhasePluginDisplayName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString TiDwellFatigueCrystallographicAnalysis::getFilterVersion() const
{
  QString version;
  QTextStream vStream(&version);
  vStream << DREAM3DReview::Version::Major() << "." << DREAM3DReview::Version::Minor() << "." << DREAM3DReview::Version::Patch();
  return version;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString TiDwellFatigueCrystallographicAnalysis::getGroupName() const
{
  return SIMPL::FilterGroups::Unsupported;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QUuid TiDwellFatigueCrystallographicAnalysis::getUuid() const
{
  return QUuid("{81e94b15-efb6-5bae-9ab1-c74099136174}");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString TiDwellFatigueCrystallographicAnalysis::getSubGroupName() const
{
  return SIMPL::FilterSubGroups::CrystallographyFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString TiDwellFatigueCrystallographicAnalysis::getHumanLabel() const
{
  return "Ti Dwell Fatigue Crystallographic Analysis";
}

// -----------------------------------------------------------------------------
TiDwellFatigueCrystallographicAnalysis::Pointer TiDwellFatigueCrystallographicAnalysis::NullPointer()
{
  return Pointer(static_cast<Self*>(nullptr));
}

// -----------------------------------------------------------------------------
std::shared_ptr<TiDwellFatigueCrystallographicAnalysis> TiDwellFatigueCrystallographicAnalysis::New()
{
  struct make_shared_enabler : public TiDwellFatigueCrystallographicAnalysis
  {
  };
  std::shared_ptr<make_shared_enabler> val = std::make_shared<make_shared_enabler>();
  val->setupFilterParameters();
  return val;
}

// -----------------------------------------------------------------------------
QString TiDwellFatigueCrystallographicAnalysis::getNameOfClass() const
{
  return QString("TiDwellFatigueCrystallographicAnalysis");
}

// -----------------------------------------------------------------------------
QString TiDwellFatigueCrystallographicAnalysis::ClassName()
{
  return QString("TiDwellFatigueCrystallographicAnalysis");
}

// -----------------------------------------------------------------------------
void TiDwellFatigueCrystallographicAnalysis::setDataContainerName(const DataArrayPath& value)
{
  m_DataContainerName = value;
}

// -----------------------------------------------------------------------------
DataArrayPath TiDwellFatigueCrystallographicAnalysis::getDataContainerName() const
{
  return m_DataContainerName;
}

// -----------------------------------------------------------------------------
void TiDwellFatigueCrystallographicAnalysis::setAlphaGlobPhasePresent(bool value)
{
  m_AlphaGlobPhasePresent = value;
}

// -----------------------------------------------------------------------------
bool TiDwellFatigueCrystallographicAnalysis::getAlphaGlobPhasePresent() const
{
  return m_AlphaGlobPhasePresent;
}

// -----------------------------------------------------------------------------
void TiDwellFatigueCrystallographicAnalysis::setAlphaGlobPhase(int value)
{
  m_AlphaGlobPhase = value;
}

// -----------------------------------------------------------------------------
int TiDwellFatigueCrystallographicAnalysis::getAlphaGlobPhase() const
{
  return m_AlphaGlobPhase;
}

// -----------------------------------------------------------------------------
void TiDwellFatigueCrystallographicAnalysis::setMTRPhase(int value)
{
  m_MTRPhase = value;
}

// -----------------------------------------------------------------------------
int TiDwellFatigueCrystallographicAnalysis::getMTRPhase() const
{
  return m_MTRPhase;
}

// -----------------------------------------------------------------------------
void TiDwellFatigueCrystallographicAnalysis::setLatticeParameterA(float value)
{
  m_LatticeParameterA = value;
}

// -----------------------------------------------------------------------------
float TiDwellFatigueCrystallographicAnalysis::getLatticeParameterA() const
{
  return m_LatticeParameterA;
}

// -----------------------------------------------------------------------------
void TiDwellFatigueCrystallographicAnalysis::setLatticeParameterC(float value)
{
  m_LatticeParameterC = value;
}

// -----------------------------------------------------------------------------
float TiDwellFatigueCrystallographicAnalysis::getLatticeParameterC() const
{
  return m_LatticeParameterC;
}

// -----------------------------------------------------------------------------
void TiDwellFatigueCrystallographicAnalysis::setStressAxis(const FloatVec3Type& value)
{
  m_StressAxis = value;
}

// -----------------------------------------------------------------------------
FloatVec3Type TiDwellFatigueCrystallographicAnalysis::getStressAxis() const
{
  return m_StressAxis;
}

// -----------------------------------------------------------------------------
void TiDwellFatigueCrystallographicAnalysis::setSubsurfaceDistance(int value)
{
  m_SubsurfaceDistance = value;
}

// -----------------------------------------------------------------------------
int TiDwellFatigueCrystallographicAnalysis::getSubsurfaceDistance() const
{
  return m_SubsurfaceDistance;
}

// -----------------------------------------------------------------------------
void TiDwellFatigueCrystallographicAnalysis::setConsiderationFraction(float value)
{
  m_ConsiderationFraction = value;
}

// -----------------------------------------------------------------------------
float TiDwellFatigueCrystallographicAnalysis::getConsiderationFraction() const
{
  return m_ConsiderationFraction;
}

// -----------------------------------------------------------------------------
void TiDwellFatigueCrystallographicAnalysis::setDoNotAssumeInitiatorPresence(bool value)
{
  m_DoNotAssumeInitiatorPresence = value;
}

// -----------------------------------------------------------------------------
bool TiDwellFatigueCrystallographicAnalysis::getDoNotAssumeInitiatorPresence() const
{
  return m_DoNotAssumeInitiatorPresence;
}

// -----------------------------------------------------------------------------
void TiDwellFatigueCrystallographicAnalysis::setInitiatorLowerThreshold(float value)
{
  m_InitiatorLowerThreshold = value;
}

// -----------------------------------------------------------------------------
float TiDwellFatigueCrystallographicAnalysis::getInitiatorLowerThreshold() const
{
  return m_InitiatorLowerThreshold;
}

// -----------------------------------------------------------------------------
void TiDwellFatigueCrystallographicAnalysis::setInitiatorUpperThreshold(float value)
{
  m_InitiatorUpperThreshold = value;
}

// -----------------------------------------------------------------------------
float TiDwellFatigueCrystallographicAnalysis::getInitiatorUpperThreshold() const
{
  return m_InitiatorUpperThreshold;
}

// -----------------------------------------------------------------------------
void TiDwellFatigueCrystallographicAnalysis::setHardFeatureLowerThreshold(float value)
{
  m_HardFeatureLowerThreshold = value;
}

// -----------------------------------------------------------------------------
float TiDwellFatigueCrystallographicAnalysis::getHardFeatureLowerThreshold() const
{
  return m_HardFeatureLowerThreshold;
}

// -----------------------------------------------------------------------------
void TiDwellFatigueCrystallographicAnalysis::setHardFeatureUpperThreshold(float value)
{
  m_HardFeatureUpperThreshold = value;
}

// -----------------------------------------------------------------------------
float TiDwellFatigueCrystallographicAnalysis::getHardFeatureUpperThreshold() const
{
  return m_HardFeatureUpperThreshold;
}

// -----------------------------------------------------------------------------
void TiDwellFatigueCrystallographicAnalysis::setSoftFeatureLowerThreshold(float value)
{
  m_SoftFeatureLowerThreshold = value;
}

// -----------------------------------------------------------------------------
float TiDwellFatigueCrystallographicAnalysis::getSoftFeatureLowerThreshold() const
{
  return m_SoftFeatureLowerThreshold;
}

// -----------------------------------------------------------------------------
void TiDwellFatigueCrystallographicAnalysis::setSoftFeatureUpperThreshold(float value)
{
  m_SoftFeatureUpperThreshold = value;
}

// -----------------------------------------------------------------------------
float TiDwellFatigueCrystallographicAnalysis::getSoftFeatureUpperThreshold() const
{
  return m_SoftFeatureUpperThreshold;
}

// -----------------------------------------------------------------------------
void TiDwellFatigueCrystallographicAnalysis::setNewCellFeatureAttributeMatrixName(const QString& value)
{
  m_NewCellFeatureAttributeMatrixName = value;
}

// -----------------------------------------------------------------------------
QString TiDwellFatigueCrystallographicAnalysis::getNewCellFeatureAttributeMatrixName() const
{
  return m_NewCellFeatureAttributeMatrixName;
}

// -----------------------------------------------------------------------------
void TiDwellFatigueCrystallographicAnalysis::setSelectedFeaturesArrayName(const QString& value)
{
  m_SelectedFeaturesArrayName = value;
}

// -----------------------------------------------------------------------------
QString TiDwellFatigueCrystallographicAnalysis::getSelectedFeaturesArrayName() const
{
  return m_SelectedFeaturesArrayName;
}

// -----------------------------------------------------------------------------
void TiDwellFatigueCrystallographicAnalysis::setInitiatorsArrayName(const QString& value)
{
  m_InitiatorsArrayName = value;
}

// -----------------------------------------------------------------------------
QString TiDwellFatigueCrystallographicAnalysis::getInitiatorsArrayName() const
{
  return m_InitiatorsArrayName;
}

// -----------------------------------------------------------------------------
void TiDwellFatigueCrystallographicAnalysis::setHardFeaturesArrayName(const QString& value)
{
  m_HardFeaturesArrayName = value;
}

// -----------------------------------------------------------------------------
QString TiDwellFatigueCrystallographicAnalysis::getHardFeaturesArrayName() const
{
  return m_HardFeaturesArrayName;
}

// -----------------------------------------------------------------------------
void TiDwellFatigueCrystallographicAnalysis::setSoftFeaturesArrayName(const QString& value)
{
  m_SoftFeaturesArrayName = value;
}

// -----------------------------------------------------------------------------
QString TiDwellFatigueCrystallographicAnalysis::getSoftFeaturesArrayName() const
{
  return m_SoftFeaturesArrayName;
}

// -----------------------------------------------------------------------------
void TiDwellFatigueCrystallographicAnalysis::setHardSoftGroupsArrayName(const QString& value)
{
  m_HardSoftGroupsArrayName = value;
}

// -----------------------------------------------------------------------------
QString TiDwellFatigueCrystallographicAnalysis::getHardSoftGroupsArrayName() const
{
  return m_HardSoftGroupsArrayName;
}

// -----------------------------------------------------------------------------
void TiDwellFatigueCrystallographicAnalysis::setCellFeatureAttributeMatrixName(const QString& value)
{
  m_CellFeatureAttributeMatrixName = value;
}

// -----------------------------------------------------------------------------
QString TiDwellFatigueCrystallographicAnalysis::getCellFeatureAttributeMatrixName() const
{
  return m_CellFeatureAttributeMatrixName;
}

// -----------------------------------------------------------------------------
void TiDwellFatigueCrystallographicAnalysis::setCellFeatureAttributeMatrixPath(const DataArrayPath& value)
{
  m_CellFeatureAttributeMatrixPath = value;
}

// -----------------------------------------------------------------------------
DataArrayPath TiDwellFatigueCrystallographicAnalysis::getCellFeatureAttributeMatrixPath() const
{
  return m_CellFeatureAttributeMatrixPath;
}

// -----------------------------------------------------------------------------
void TiDwellFatigueCrystallographicAnalysis::setFeatureIdsArrayPath(const DataArrayPath& value)
{
  m_FeatureIdsArrayPath = value;
}

// -----------------------------------------------------------------------------
DataArrayPath TiDwellFatigueCrystallographicAnalysis::getFeatureIdsArrayPath() const
{
  return m_FeatureIdsArrayPath;
}

// -----------------------------------------------------------------------------
void TiDwellFatigueCrystallographicAnalysis::setCellParentIdsArrayName(const QString& value)
{
  m_CellParentIdsArrayName = value;
}

// -----------------------------------------------------------------------------
QString TiDwellFatigueCrystallographicAnalysis::getCellParentIdsArrayName() const
{
  return m_CellParentIdsArrayName;
}

// -----------------------------------------------------------------------------
void TiDwellFatigueCrystallographicAnalysis::setFeatureParentIdsArrayName(const QString& value)
{
  m_FeatureParentIdsArrayName = value;
}

// -----------------------------------------------------------------------------
QString TiDwellFatigueCrystallographicAnalysis::getFeatureParentIdsArrayName() const
{
  return m_FeatureParentIdsArrayName;
}

// -----------------------------------------------------------------------------
void TiDwellFatigueCrystallographicAnalysis::setActiveArrayName(const QString& value)
{
  m_ActiveArrayName = value;
}

// -----------------------------------------------------------------------------
QString TiDwellFatigueCrystallographicAnalysis::getActiveArrayName() const
{
  return m_ActiveArrayName;
}

// -----------------------------------------------------------------------------
void TiDwellFatigueCrystallographicAnalysis::setFeatureEulerAnglesArrayPath(const DataArrayPath& value)
{
  m_FeatureEulerAnglesArrayPath = value;
}

// -----------------------------------------------------------------------------
DataArrayPath TiDwellFatigueCrystallographicAnalysis::getFeatureEulerAnglesArrayPath() const
{
  return m_FeatureEulerAnglesArrayPath;
}

// -----------------------------------------------------------------------------
void TiDwellFatigueCrystallographicAnalysis::setFeaturePhasesArrayPath(const DataArrayPath& value)
{
  m_FeaturePhasesArrayPath = value;
}

// -----------------------------------------------------------------------------
DataArrayPath TiDwellFatigueCrystallographicAnalysis::getFeaturePhasesArrayPath() const
{
  return m_FeaturePhasesArrayPath;
}

// -----------------------------------------------------------------------------
void TiDwellFatigueCrystallographicAnalysis::setNeighborListArrayPath(const DataArrayPath& value)
{
  m_NeighborListArrayPath = value;
}

// -----------------------------------------------------------------------------
DataArrayPath TiDwellFatigueCrystallographicAnalysis::getNeighborListArrayPath() const
{
  return m_NeighborListArrayPath;
}

// -----------------------------------------------------------------------------
void TiDwellFatigueCrystallographicAnalysis::setCentroidsArrayPath(const DataArrayPath& value)
{
  m_CentroidsArrayPath = value;
}

// -----------------------------------------------------------------------------
DataArrayPath TiDwellFatigueCrystallographicAnalysis::getCentroidsArrayPath() const
{
  return m_CentroidsArrayPath;
}

// -----------------------------------------------------------------------------
void TiDwellFatigueCrystallographicAnalysis::setCrystalStructuresArrayPath(const DataArrayPath& value)
{
  m_CrystalStructuresArrayPath = value;
}

// -----------------------------------------------------------------------------
DataArrayPath TiDwellFatigueCrystallographicAnalysis::getCrystalStructuresArrayPath() const
{
  return m_CrystalStructuresArrayPath;
}
