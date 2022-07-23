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
#include "InsertTransformationPhases.h"

#include <QtCore/QTextStream>

#include "SIMPLib/Common/Constants.h"
#include "SIMPLib/DataContainers/DataContainer.h"
#include "SIMPLib/DataContainers/DataContainerArray.h"
#include "SIMPLib/FilterParameters/AbstractFilterParametersReader.h"
#include "SIMPLib/FilterParameters/AttributeMatrixSelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/ChoiceFilterParameter.h"
#include "SIMPLib/FilterParameters/DataArrayCreationFilterParameter.h"
#include "SIMPLib/FilterParameters/DataArraySelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/FloatFilterParameter.h"
#include "SIMPLib/FilterParameters/FloatVec3FilterParameter.h"
#include "SIMPLib/FilterParameters/IntFilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedBooleanFilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedPathCreationFilterParameter.h"
#include "SIMPLib/FilterParameters/SeparatorFilterParameter.h"
#include "SIMPLib/Geometry/ImageGeom.h"
#include "SIMPLib/Math/MatrixMath.h"
#include "SIMPLib/Math/SIMPLibRandom.h"

#include "EbsdLib/Core/Orientation.hpp"
#include "EbsdLib/LaueOps/LaueOps.h"

#include "DREAM3DReview/DREAM3DReviewConstants.h"
#include "DREAM3DReview/DREAM3DReviewVersion.h"

/* Create Enumerations to allow the created Attribute Arrays to take part in renaming */
enum createdPathID : RenameDataPath::DataID_t
{
  DataArrayID30 = 30,
  DataArrayID31 = 31,
  DataArrayID32 = 32,
};

const static float m_pi = static_cast<float>(M_PI);

QTextStream& operator<<(QTextStream& os, const ShapeType::Type& value)
{
  os << static_cast<ShapeType::EnumType>(value);
  return os;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
InsertTransformationPhases::InsertTransformationPhases()
{
  m_OrientationOps = LaueOps::GetAllOrientationOps();
  m_TransformationPhaseHabitPlane[0] = 1.0f;
  m_TransformationPhaseHabitPlane[1] = 1.0f;
  m_TransformationPhaseHabitPlane[2] = 1.0f;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
InsertTransformationPhases::~InsertTransformationPhases() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void InsertTransformationPhases::setupFilterParameters()
{
  FilterParameterVectorType parameters;
  parameters.push_back(SIMPL_NEW_INTEGER_FP("Parent Phase", ParentPhase, FilterParameter::Category::Parameter, InsertTransformationPhases));
  {
    ChoiceFilterParameter::Pointer option = ChoiceFilterParameter::New();
    option->setHumanLabel("Transformation Phase Crystal Structure");
    option->setPropertyName("TransCrystalStruct");
    option->setSetterCallback(SIMPL_BIND_SETTER(InsertTransformationPhases, this, TransCrystalStruct));
    option->setGetterCallback(SIMPL_BIND_GETTER(InsertTransformationPhases, this, TransCrystalStruct));
    // The choices here are IN ORDER of the enumerations from the EBSDLib. DO NOT CHANGE THE ORDER.
    std::vector<std::string> names = LaueOps::GetLaueNames();
    std::vector<QString> choices;
    for(const auto& name : names)
    {
      choices.push_back(QString::fromStdString(name));
    }

    choices.pop_back(); // Remove the last name because we don't need it.
    option->setCategory(FilterParameter::Category::Parameter);
    option->setChoices(choices);
    parameters.push_back(option);
  }
  parameters.push_back(SIMPL_NEW_FLOAT_FP("Transformation Phase Misorientation", TransformationPhaseMisorientation, FilterParameter::Category::Parameter, InsertTransformationPhases));
  std::vector<QString> linkedProps1;
  linkedProps1.push_back("TransformationPhaseHabitPlane");
  linkedProps1.push_back("UseAllVariants");
  parameters.push_back(SIMPL_NEW_LINKED_BOOL_FP("Define Habit Plane", DefineHabitPlane, FilterParameter::Category::Parameter, InsertTransformationPhases, linkedProps1));
  parameters.push_back(SIMPL_NEW_FLOAT_VEC3_FP("Transformation Phase Habit Plane", TransformationPhaseHabitPlane, FilterParameter::Category::Parameter, InsertTransformationPhases));

  std::vector<QString> linkedProps2 = {"CoherentFrac"};
  parameters.push_back(SIMPL_NEW_LINKED_BOOL_FP("Use All Variants", UseAllVariants, FilterParameter::Category::Parameter, InsertTransformationPhases, linkedProps2));
  parameters.push_back(SIMPL_NEW_FLOAT_FP("Coherent Fraction", CoherentFrac, FilterParameter::Category::Parameter, InsertTransformationPhases));
  parameters.push_back(SIMPL_NEW_FLOAT_FP("Transformation Phase Thickness", TransformationPhaseThickness, FilterParameter::Category::Parameter, InsertTransformationPhases));
  parameters.push_back(
      SIMPL_NEW_INTEGER_FP("Average Number Of Transformation Phases Per Feature", NumTransformationPhasesPerFeature, FilterParameter::Category::Parameter, InsertTransformationPhases));
  parameters.push_back(SIMPL_NEW_FLOAT_FP("Peninsula Transformation Phase Fraction", PeninsulaFrac, FilterParameter::Category::Parameter, InsertTransformationPhases));

  parameters.push_back(SeparatorFilterParameter::Create("Cell Data", FilterParameter::Category::RequiredArray));
  {
    DataArraySelectionFilterParameter::RequirementType req = DataArraySelectionFilterParameter::CreateRequirement(SIMPL::TypeNames::Int32, 1, AttributeMatrix::Type::Cell, IGeometry::Type::Any);
    parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Feature Ids", FeatureIdsArrayPath, FilterParameter::Category::RequiredArray, InsertTransformationPhases, req));
  }
  {
    DataArraySelectionFilterParameter::RequirementType req = DataArraySelectionFilterParameter::CreateRequirement(SIMPL::TypeNames::Float, 3, AttributeMatrix::Type::Cell, IGeometry::Type::Any);
    parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Euler Angles", CellEulerAnglesArrayPath, FilterParameter::Category::RequiredArray, InsertTransformationPhases, req));
  }
  {
    DataArraySelectionFilterParameter::RequirementType req = DataArraySelectionFilterParameter::CreateRequirement(SIMPL::TypeNames::Int32, 1, AttributeMatrix::Type::Cell, IGeometry::Type::Any);
    parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Phases", CellPhasesArrayPath, FilterParameter::Category::RequiredArray, InsertTransformationPhases, req));
  }

  parameters.push_back(SeparatorFilterParameter::Create("Cell Feature Data", FilterParameter::Category::RequiredArray));
  {
    AttributeMatrixSelectionFilterParameter::RequirementType req = AttributeMatrixSelectionFilterParameter::CreateRequirement(AttributeMatrix::Type::CellFeature, IGeometry::Type::Any);
    parameters.push_back(SIMPL_NEW_AM_SELECTION_FP("Cell Feature Attribute Matrix", CellFeatureAttributeMatrixName, FilterParameter::Category::RequiredArray, InsertTransformationPhases, req));
  }
  {
    DataArraySelectionFilterParameter::RequirementType req = DataArraySelectionFilterParameter::CreateCategoryRequirement(SIMPL::Defaults::AnyPrimitive, 3, AttributeMatrix::Category::Feature);
    parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Average Euler Angles", FeatureEulerAnglesArrayPath, FilterParameter::Category::RequiredArray, InsertTransformationPhases, req));
  }
  {
    DataArraySelectionFilterParameter::RequirementType req = DataArraySelectionFilterParameter::CreateCategoryRequirement(SIMPL::Defaults::AnyPrimitive, 4, AttributeMatrix::Category::Feature);
    parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Average Quaternions", AvgQuatsArrayPath, FilterParameter::Category::RequiredArray, InsertTransformationPhases, req));
  }
  {
    DataArraySelectionFilterParameter::RequirementType req = DataArraySelectionFilterParameter::CreateCategoryRequirement(SIMPL::Defaults::AnyPrimitive, 3, AttributeMatrix::Category::Feature);
    parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Centroids", CentroidsArrayPath, FilterParameter::Category::RequiredArray, InsertTransformationPhases, req));
  }
  {
    DataArraySelectionFilterParameter::RequirementType req = DataArraySelectionFilterParameter::CreateCategoryRequirement(SIMPL::Defaults::AnyPrimitive, 1, AttributeMatrix::Category::Feature);
    parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Equivalent Diameters", EquivalentDiametersArrayPath, FilterParameter::Category::RequiredArray, InsertTransformationPhases, req));
  }
  {
    DataArraySelectionFilterParameter::RequirementType req = DataArraySelectionFilterParameter::CreateCategoryRequirement(SIMPL::Defaults::AnyPrimitive, 1, AttributeMatrix::Category::Feature);
    parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Phases", FeaturePhasesArrayPath, FilterParameter::Category::RequiredArray, InsertTransformationPhases, req));
  }

  parameters.push_back(SeparatorFilterParameter::Create("Cell Ensemble Data", FilterParameter::Category::RequiredArray));
  {
    AttributeMatrixSelectionFilterParameter::RequirementType req = AttributeMatrixSelectionFilterParameter::CreateRequirement(AttributeMatrix::Type::CellEnsemble, IGeometry::Type::Any);
    parameters.push_back(
        SIMPL_NEW_AM_SELECTION_FP("Cell Ensemble Attribute Matrix", StatsGenCellEnsembleAttributeMatrixPath, FilterParameter::Category::RequiredArray, InsertTransformationPhases, req));
  }
  {
    DataArraySelectionFilterParameter::RequirementType req = DataArraySelectionFilterParameter::CreateCategoryRequirement(SIMPL::Defaults::AnyPrimitive, 1, AttributeMatrix::Category::Ensemble);
    parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Crystal Structures", CrystalStructuresArrayPath, FilterParameter::Category::RequiredArray, InsertTransformationPhases, req));
  }
  {
    DataArraySelectionFilterParameter::RequirementType req = DataArraySelectionFilterParameter::CreateCategoryRequirement(SIMPL::Defaults::AnyPrimitive, 1, AttributeMatrix::Category::Ensemble);
    parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Phase Types", PhaseTypesArrayPath, FilterParameter::Category::RequiredArray, InsertTransformationPhases, req));
  }
  {
    DataArraySelectionFilterParameter::RequirementType req = DataArraySelectionFilterParameter::CreateCategoryRequirement(SIMPL::Defaults::AnyPrimitive, 1, AttributeMatrix::Category::Ensemble);
    parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Shape Types", ShapeTypesArrayPath, FilterParameter::Category::RequiredArray, InsertTransformationPhases, req));
  }
  {
    DataArraySelectionFilterParameter::RequirementType req = DataArraySelectionFilterParameter::CreateCategoryRequirement(SIMPL::Defaults::AnyPrimitive, 1, AttributeMatrix::Category::Ensemble);
    parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Number of Features", NumFeaturesArrayPath, FilterParameter::Category::RequiredArray, InsertTransformationPhases, req));
  }

  parameters.push_back(SeparatorFilterParameter::Create("Cell Feature Data", FilterParameter::Category::CreatedArray));
  parameters.push_back(SIMPL_NEW_DA_WITH_LINKED_AM_FP("Parent Ids", FeatureParentIdsArrayName, CellFeatureAttributeMatrixName, CellFeatureAttributeMatrixName, FilterParameter::Category::CreatedArray,
                                                      InsertTransformationPhases));
  {
    DataArrayCreationFilterParameter::RequirementType req = DataArrayCreationFilterParameter::CreateRequirement(AttributeMatrix::Category::Feature);
    parameters.push_back(SIMPL_NEW_DA_CREATION_FP("Number of Features Per Parent", NumFeaturesPerParentArrayPath, FilterParameter::Category::CreatedArray, InsertTransformationPhases, req));
  }

  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void InsertTransformationPhases::readFilterParameters(AbstractFilterParametersReader* reader, int index)
{
  reader->openFilterGroup(this, index);
  setParentPhase(reader->readValue("ParentPhase", getParentPhase()));
  setTransCrystalStruct(reader->readValue("TransCrystalStruct", getTransCrystalStruct()));
  setTransformationPhaseMisorientation(reader->readValue("TransformationPhaseMisorientation", getTransformationPhaseMisorientation()));
  setDefineHabitPlane(reader->readValue("DefineHabitPlane", getDefineHabitPlane()));
  setTransformationPhaseHabitPlane(reader->readFloatVec3("TransformationPhaseHabitPlane", getTransformationPhaseHabitPlane()));
  setUseAllVariants(reader->readValue("UseAllVariants", getUseAllVariants()));
  setCoherentFrac(reader->readValue("CoherentFrac", getCoherentFrac()));
  setTransformationPhaseThickness(reader->readValue("TransformationPhaseThickness", getTransformationPhaseThickness()));
  setNumTransformationPhasesPerFeature(reader->readValue("NumTransformationPhasesPerFeature", getNumTransformationPhasesPerFeature()));
  setPeninsulaFrac(reader->readValue("PeninsulaFrac", getPeninsulaFrac()));

  setStatsGenCellEnsembleAttributeMatrixPath(reader->readDataArrayPath("StatsGenCellEnsembleAttributeMatrixPath", getStatsGenCellEnsembleAttributeMatrixPath()));
  setCellFeatureAttributeMatrixName(reader->readDataArrayPath("CellFeatureAttributeMatrixName", getCellFeatureAttributeMatrixName()));
  setFeatureIdsArrayPath(reader->readDataArrayPath("FeatureIdsArrayPath", getFeatureIdsArrayPath()));
  setCellEulerAnglesArrayPath(reader->readDataArrayPath("CellEulerAnglesArrayPath", getCellEulerAnglesArrayPath()));
  setCentroidsArrayPath(reader->readDataArrayPath("CentroidsArrayPath", getCentroidsArrayPath()));
  setEquivalentDiametersArrayPath(reader->readDataArrayPath("EquivalentDiametersArrayPath", getEquivalentDiametersArrayPath()));
  setFeatureEulerAnglesArrayPath(reader->readDataArrayPath("FeatureEulerAnglesArrayPath", getFeatureEulerAnglesArrayPath()));
  setFeaturePhasesArrayPath(reader->readDataArrayPath("FeaturePhasesArrayPath", getFeaturePhasesArrayPath()));
  setFeatureParentIdsArrayName(reader->readString("FeatureParentIdsArrayName", getFeatureParentIdsArrayName()));
  setNumFeaturesPerParentArrayPath(reader->readDataArrayPath("NumFeaturesPerParentArrayPath", getNumFeaturesPerParentArrayPath()));
  setCrystalStructuresArrayPath(reader->readDataArrayPath("CrystalStructuresArrayPath", getCrystalStructuresArrayPath()));
  setPhaseTypesArrayPath(reader->readDataArrayPath("PhaseTypesArrayPath", getPhaseTypesArrayPath()));
  setShapeTypesArrayPath(reader->readDataArrayPath("ShapeTypesArrayPath", getShapeTypesArrayPath()));
  setNumFeaturesArrayPath(reader->readDataArrayPath("NumFeaturesArrayPath", getNumFeaturesArrayPath()));
  reader->closeFilterGroup();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void InsertTransformationPhases::updateFeatureInstancePointers()
{
  clearErrorCode();
  clearWarningCode();

  if(nullptr != m_AvgQuatsPtr.lock())
  {
    m_AvgQuats = m_AvgQuatsPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */
  if(nullptr != m_CentroidsPtr.lock())
  {
    m_Centroids = m_CentroidsPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */
  if(nullptr != m_FeatureEulerAnglesPtr.lock())
  {
    m_FeatureEulerAngles = m_FeatureEulerAnglesPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */
  if(nullptr != m_EquivalentDiametersPtr.lock())
  {
    m_EquivalentDiameters = m_EquivalentDiametersPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */
  if(nullptr != m_FeaturePhasesPtr.lock())
  {
    m_FeaturePhases = m_FeaturePhasesPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */
  if(nullptr != m_FeatureParentIdsPtr.lock())
  {
    m_FeatureParentIds = m_FeatureParentIdsPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */
  if(nullptr != m_NumFeaturesPerParentPtr.lock())
  {
    m_NumFeaturesPerParent = m_NumFeaturesPerParentPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void InsertTransformationPhases::updateStatsGenEnsembleInstancePointers()
{
  clearErrorCode();
  clearWarningCode();

  if(nullptr != m_CrystalStructuresPtr.lock())
  {
    m_CrystalStructures = m_CrystalStructuresPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */
  if(nullptr != m_PhaseTypesPtr.lock())
  {
    m_PhaseTypes = m_PhaseTypesPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */
  if(nullptr != m_ShapeTypesPtr.lock())
  {
    m_ShapeTypes = m_ShapeTypesPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void InsertTransformationPhases::updateVolEnsembleInstancePointers()
{
  clearErrorCode();
  clearWarningCode();

  if(nullptr != m_NumFeaturesPtr.lock())
  {
    m_NumFeatures = m_NumFeaturesPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void InsertTransformationPhases::initialize()
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void InsertTransformationPhases::dataCheck()
{

  clearErrorCode();
  clearWarningCode();
  DataArrayPath tempPath;

  DataContainer::Pointer m = getDataContainerArray()->getPrereqDataContainer(this, m_FeatureIdsArrayPath.getDataContainerName(), false);
  if(getErrorCode() < 0 || m == nullptr)
  {
    return;
  }

  AttributeMatrix::Pointer statsGenAttrMat = m->getPrereqAttributeMatrix(this, getStatsGenCellEnsembleAttributeMatrixPath().getAttributeMatrixName(), -301);
  if(getErrorCode() < 0 || statsGenAttrMat == nullptr)
  {
    return;
  }
  AttributeMatrix::Pointer volAttrMat = m->getPrereqAttributeMatrix(this, getNumFeaturesPerParentArrayPath().getAttributeMatrixName(), -301);
  if(getErrorCode() < 0 || volAttrMat == nullptr)
  {
    return;
  }

  ImageGeom::Pointer image = m->getPrereqGeometry<ImageGeom>(this);
  if(getErrorCode() < 0 || nullptr == image.get())
  {
    return;
  }

  std::vector<size_t> dims(1, 1);
  // Cell Data
  m_FeatureIdsPtr = getDataContainerArray()->getPrereqArrayFromPath<DataArray<int32_t>>(this, getFeatureIdsArrayPath(), dims);
  if(nullptr != m_FeatureIdsPtr.lock())
  {
    m_FeatureIds = m_FeatureIdsPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */

  dims[0] = 3;
  m_CellEulerAnglesPtr = getDataContainerArray()->getPrereqArrayFromPath<DataArray<float>>(this, getCellEulerAnglesArrayPath(), dims);
  if(nullptr != m_CellEulerAnglesPtr.lock())
  {
    m_CellEulerAngles = m_CellEulerAnglesPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */

  dims[0] = 1;
  m_CellPhasesPtr = getDataContainerArray()->getPrereqArrayFromPath<Int32ArrayType>(this, getCellPhasesArrayPath(), dims);
  if(nullptr != m_CellEulerAnglesPtr.lock())
  {
    m_CellPhases = m_CellPhasesPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */

  // Feature Data
  dims[0] = 4;
  m_AvgQuatsPtr = getDataContainerArray()->getPrereqArrayFromPath<DataArray<float>>(this, getAvgQuatsArrayPath(), dims);
  if(nullptr != m_AvgQuatsPtr.lock())
  {
    m_AvgQuats = m_AvgQuatsPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */

  dims[0] = 1;
  m_EquivalentDiametersPtr = getDataContainerArray()->getPrereqArrayFromPath<DataArray<float>>(this, getEquivalentDiametersArrayPath(), dims);
  if(nullptr != m_EquivalentDiametersPtr.lock())
  {
    m_EquivalentDiameters = m_EquivalentDiametersPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */

  dims[0] = 3;
  m_CentroidsPtr = getDataContainerArray()->getPrereqArrayFromPath<DataArray<float>>(this, getCentroidsArrayPath(), dims);
  if(nullptr != m_CentroidsPtr.lock())
  {
    m_Centroids = m_CentroidsPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */

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

  // New Feature Data

  tempPath.update(getCellFeatureAttributeMatrixName().getDataContainerName(), getCellFeatureAttributeMatrixName().getAttributeMatrixName(), getFeatureParentIdsArrayName());
  m_FeatureParentIdsPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<int32_t>>(this, tempPath, -1, dims, "", DataArrayID31);
  if(nullptr != m_FeatureParentIdsPtr.lock())
  {
    m_FeatureParentIds = m_FeatureParentIdsPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */

  tempPath.update(getCellFeatureAttributeMatrixName().getDataContainerName(), getCellFeatureAttributeMatrixName().getAttributeMatrixName(), getNumFeaturesPerParentArrayPath().getDataArrayName());
  m_NumFeaturesPerParentPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<int32_t>>(this, tempPath, 0, dims, "", DataArrayID32);
  if(nullptr != m_NumFeaturesPerParentPtr.lock())
  {
    m_NumFeaturesPerParent = m_NumFeaturesPerParentPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */

  // Ensemble Data
  m_CrystalStructuresPtr = getDataContainerArray()->getPrereqArrayFromPath<DataArray<unsigned int>>(this, getCrystalStructuresArrayPath(), dims);
  if(nullptr != m_CrystalStructuresPtr.lock())
  {
    m_CrystalStructures = m_CrystalStructuresPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */

  m_PhaseTypesPtr = getDataContainerArray()->getPrereqArrayFromPath<DataArray<unsigned int>>(this, getPhaseTypesArrayPath(), dims);
  if(nullptr != m_PhaseTypesPtr.lock())
  {
    m_PhaseTypes = m_PhaseTypesPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */

  m_ShapeTypesPtr = getDataContainerArray()->getPrereqArrayFromPath<DataArray<ShapeType::EnumType>>(this, getShapeTypesArrayPath(), dims);
  if(nullptr != m_ShapeTypesPtr.lock())
  {
    m_ShapeTypes = m_ShapeTypesPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */

  m_NumFeaturesPtr = getDataContainerArray()->getPrereqArrayFromPath<DataArray<int32_t>>(this, getNumFeaturesArrayPath(), dims);
  if(nullptr != m_NumFeaturesPtr.lock())
  {
    m_NumFeatures = m_NumFeaturesPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void InsertTransformationPhases::execute()
{
  dataCheck();
  if(getErrorCode() < 0)
  {
    return;
  }

  QString msg;
  QTextStream ss(&msg);
  ss << getHumanLabel() << " is an experimental filter that may or may *not* give the correct results. Use at your own risk";
  setWarningCondition(1, msg);

  m_NormalizedTransformationPhaseHabitPlane = m_TransformationPhaseHabitPlane;

  float denom = std::sqrt(m_NormalizedTransformationPhaseHabitPlane[0] * m_NormalizedTransformationPhaseHabitPlane[0] +
                          m_NormalizedTransformationPhaseHabitPlane[1] * m_NormalizedTransformationPhaseHabitPlane[1] +
                          m_NormalizedTransformationPhaseHabitPlane[2] * m_NormalizedTransformationPhaseHabitPlane[2]);
  m_NormalizedTransformationPhaseHabitPlane[0] = m_NormalizedTransformationPhaseHabitPlane[0] / denom;
  m_NormalizedTransformationPhaseHabitPlane[1] = m_NormalizedTransformationPhaseHabitPlane[1] / denom;
  m_NormalizedTransformationPhaseHabitPlane[2] = m_NormalizedTransformationPhaseHabitPlane[2] / denom;

  // DataContainer::Pointer m = getDataContainerArray()->getDataContainer(getFeatureIdsArrayPath().getDataContainerName());
  DataContainerArray::Pointer dca = getDataContainerArray();

  // defining separate ensemble attribute matrix for statsgen & vol
  // AttributeMatrix::Pointer statsGenAttrMat = dca->getAttributeMatrix(getStatsGenCellEnsembleAttributeMatrixPath());
  AttributeMatrix::Pointer volAttrMat = dca->getAttributeMatrix(getNumFeaturesPerParentArrayPath());

  size_t totalFeatures = m_FeaturePhasesPtr.lock()->getNumberOfTuples();
  size_t numensembles = m_PhaseTypesPtr.lock()->getNumberOfTuples();

  // resizing statsgen attribute matrix to add new phase info
  //    std::vector<size_t> tDims(1, numensembles + 1);
  //  statsGenAttrMat->resizeAttributeArrays(tDims);
  updateStatsGenEnsembleInstancePointers();

  // hard-coded to the below stats for now
  m_CrystalStructures[numensembles] = static_cast<uint32_t>(m_TransCrystalStruct);
  m_PhaseTypes[numensembles] = static_cast<PhaseType::EnumType>(PhaseType::Type::Transformation);
  m_ShapeTypes[numensembles] = static_cast<ShapeType::EnumType>(ShapeType::Type::Ellipsoid);

  // start insert transformation phases routine
  insertTransformationPhases();

  notifyStatusMessage("Placement Complete");

  // updating number of ensembles and features
  totalFeatures = m_FeaturePhasesPtr.lock()->getNumberOfTuples();
  numensembles = m_PhaseTypesPtr.lock()->getNumberOfTuples();

  // resizing vol attribute matrix
  volAttrMat->resizeAttributeArrays({numensembles + 1});
  updateVolEnsembleInstancePointers();

  // finding ensemble level number of features per phase
  for(size_t i = 1; i < numensembles; ++i)
  {
    m_NumFeatures[i] = 0;
  }
  for(size_t i = 1; i < totalFeatures; ++i)
  {
    ++m_NumFeatures[m_FeaturePhases[i]];
  }

  // filter_calls();

  notifyStatusMessage("Execute Complete");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void InsertTransformationPhases::insertTransformationPhases()
{
  DataContainer::Pointer m = getDataContainerArray()->getDataContainer(getFeatureIdsArrayPath().getDataContainerName());
  SIMPL_RANDOMNG_NEW()

  size_t totalFeatures = m_FeaturePhasesPtr.lock()->getNumberOfTuples();
  // int64_t totalPoints = static_cast<size_t>(m_FeatureIdsPtr.lock()->getNumberOfTuples());
  std::vector<size_t> tDims(1, 1);

  // find the minimum resolution
  FloatVec3Type spacing = m->getGeometryAs<ImageGeom>()->getSpacing();
  float minRes = spacing[0];
  if(minRes > spacing[1])
  {
    minRes = spacing[1];
  }
  if(minRes > spacing[2])
  {
    minRes = spacing[2];
  }

  float sampleHabitPlane[3] = {0.0f};
  // NOTE: incoherent habit planes e[-20,20] & are ints
  float crystalHabitPlane[3] = {static_cast<float>(rg.genrand_res53() * 20), static_cast<float>(rg.genrand_res53() * 20), static_cast<float>(rg.genrand_res53() * 20)};

  float g[3][3], gT[3][3], rotMat[3][3], newMat[3][3], newMatCopy[3][3], symMat[3][3];
  float plateThickness = 0.0f;
  size_t numFeatures = totalFeatures;
  float d, random = 0.0f, random2 = 0.0f;
  float sig3 = m_TransformationPhaseMisorientation * (m_pi / 180.0f);
  float e[3];
  float traceMin = -1.0f;
  float trace = 0.0f;
  int minPos = 0;
  std::vector<float> shifts;
  size_t numTransformationPhases = 0;
  bool createdTransformationPhase = false;

  float voxelDiagonal = sqrtf(spacing[0] * spacing[0] + spacing[1] * spacing[1] + spacing[2] * spacing[2]);

  for(size_t curFeature = 1; curFeature < numFeatures; ++curFeature)
  {
    if(m_FeaturePhases[curFeature] == m_ParentPhase)
    {
      // set the grain Id to the parent Id for if/when the features get uniquely renumbered
      m_FeatureParentIds[curFeature] = static_cast<int32_t>(curFeature);

      if(m_DefineHabitPlane)
      {
        crystalHabitPlane[0] = m_NormalizedTransformationPhaseHabitPlane[0];
        crystalHabitPlane[1] = m_NormalizedTransformationPhaseHabitPlane[1];
        crystalHabitPlane[2] = m_NormalizedTransformationPhaseHabitPlane[2];
      }
      // pick a habit plane variant if user desires
      if(m_UseAllVariants)
      {
        random = static_cast<float>(rg.genrand_res53());
        for(int i = 0; i < 3; ++i)
        {
          // decide whether to make it coherent or incoherent &
          // decide whether its negative or positive
          if(random < m_CoherentFrac && random2 < 0.5f)
          {
            crystalHabitPlane[i] *= -1.0f;
          }
        }
        // jumble the i,j,k order
        float tempCrystalHabitPlane[3] = {crystalHabitPlane[0], crystalHabitPlane[1], crystalHabitPlane[2]};
        random = static_cast<float>(rg.genrand_res53());
        if(random < 0.16667f)
        {
          crystalHabitPlane[0] = tempCrystalHabitPlane[0];
          crystalHabitPlane[1] = tempCrystalHabitPlane[1];
          crystalHabitPlane[2] = tempCrystalHabitPlane[2];
        }
        else if(random >= 0.16667f && random < 0.33333f)
        {
          crystalHabitPlane[0] = tempCrystalHabitPlane[0];
          crystalHabitPlane[1] = tempCrystalHabitPlane[2];
          crystalHabitPlane[2] = tempCrystalHabitPlane[1];
        }
        else if(random >= 0.33333f && random < 0.5f)
        {
          crystalHabitPlane[0] = tempCrystalHabitPlane[1];
          crystalHabitPlane[1] = tempCrystalHabitPlane[0];
          crystalHabitPlane[2] = tempCrystalHabitPlane[2];
        }
        else if(random >= 0.5f && random < 0.66667f)
        {
          crystalHabitPlane[0] = tempCrystalHabitPlane[1];
          crystalHabitPlane[1] = tempCrystalHabitPlane[2];
          crystalHabitPlane[2] = tempCrystalHabitPlane[0];
        }
        else if(random >= 0.66667f && random < 0.83333f)
        {
          crystalHabitPlane[0] = tempCrystalHabitPlane[2];
          crystalHabitPlane[1] = tempCrystalHabitPlane[0];
          crystalHabitPlane[2] = tempCrystalHabitPlane[1];
        }
        else
        {
          crystalHabitPlane[0] = tempCrystalHabitPlane[2];
          crystalHabitPlane[1] = tempCrystalHabitPlane[1];
          crystalHabitPlane[2] = tempCrystalHabitPlane[0];
        }
      }
      // find where the habit plane points
      float* currentAvgQuatPtr = m_AvgQuats + curFeature * 4;
      {
        QuatF q1(currentAvgQuatPtr[0], currentAvgQuatPtr[1], currentAvgQuatPtr[2], currentAvgQuatPtr[3]);
        OrientationTransformation::qu2om<QuatF, OrientationF>(q1).toGMatrix(g);
      }
      MatrixMath::Transpose3x3(g, gT);
      MatrixMath::Multiply3x3with3x1(gT, crystalHabitPlane, sampleHabitPlane);

      //    rotMat[0][0]=0.5774;
      //    rotMat[0][1]=0.5774;
      //    rotMat[0][2]=0.5774;
      //    rotMat[1][0]=-0.4082;
      //    rotMat[1][1]=-0.4082;
      //    rotMat[1][2]=0.8165;
      //    rotMat[2][0]=0.7071;
      //    rotMat[2][1]=-0.7071;
      //    rotMat[2][2]=0;

      // generate transformation phase orientation with a user-defined rotation about the habit plane
      OrientationTransformation::ax2om<OrientationF, OrientationF>(OrientationF(crystalHabitPlane[0], crystalHabitPlane[1], crystalHabitPlane[2], sig3)).toGMatrix(rotMat);
      MatrixMath::Multiply3x3with3x3(rotMat, g, newMat);

      // find the minimum angle
      MatrixMath::Copy3x3(newMat, newMatCopy);
      // Get our LaueOps pointer for the selected crystal structure
      LaueOps::Pointer orientOps = m_OrientationOps[m_TransCrystalStruct];

      // get number of symmetry operators
      traceMin = -1.0f;
      minPos = 0;
      int n_sym = orientOps->getNumSymOps();
      for(int i = 0; i < n_sym; ++i)
      {
        orientOps->getMatSymOp(i, symMat);
        MatrixMath::Multiply3x3with3x3(symMat, newMatCopy, newMat);
        trace = newMat[0][0] + newMat[1][1] + newMat[2][2];
        if(trace > traceMin)
        {
          traceMin = trace;
          minPos = i;
        }
      }
      // assign our symmetry matrix to that which produced the minimum angle
      orientOps->getMatSymOp(minPos, symMat);
      MatrixMath::Multiply3x3with3x3(symMat, newMatCopy, newMat);
      OrientationF eOut = OrientationTransformation::om2eu<OrientationF, OrientationF>(OrientationF(symMat));

      QuatF q2 = OrientationTransformation::eu2qu<OrientationF, QuatF>(eOut);

      // define plate = user input fraction of eq dia centered at centroid
      // NOTE: we multiply by 0.5 because the transformation phase thickness will be established by
      // a search from both sides
      plateThickness = m_EquivalentDiameters[curFeature] * m_TransformationPhaseThickness * 0.5f;
      // if the plate thickness is less than the minimum dimension resolution,
      // then don't try to insert a transformation phase because the grain is too small and
      // the transformation phase thickness will be too small
      if(plateThickness > minRes)
      {
        // select number of transformation phases to insert per grain from 2x user input to zero
        numTransformationPhases = static_cast<size_t>(rg.genrand_res53() * double(2 * m_NumTransformationPhasesPerFeature + 1));
        shifts.resize(numTransformationPhases);
        int attempt = 0;

        // each loop represents an attempt to insert a single transformation phase (all in the same grain)
        for(size_t i = 0; i < static_cast<size_t>(numTransformationPhases); ++i)
        {
          // define the shift placement from centroid
          random = static_cast<float>(rg.genrand_res53());
          // shift defined as random (float) e[0,1] x grain's eq dia x 0.5
          // NOTE: we multiply by 0.5 because it starts at the center and we don't
          // want to run off the end of the grain
          shifts[i] = random * m_EquivalentDiameters[curFeature] * 0.5f;
          random = static_cast<float>(rg.genrand_res53());
          if(random < 0.5f)
          {
            shifts[i] = -shifts[i];
          }

          // check if new transformation phase will overlap an old transformation phase
          size_t ii = i;
          for(size_t j = 0; j <= ii; ++j)
          {
            // if adding more than one transformation phase in a grain, check the current shift from
            // center against the previous ones to make sure the new transformation phase does not
            // overlap
            // NOTE: that the transformation phases can touch since we're using equiv dia,
            // they can be high aspect ratio grains inserted in the low aspect direction
            // which skews the calculation (FIX)
            // tolerance is set to 3x the plate thickness + 2x the voxel diagonal -- this empirically
            // keeps the transformation phases separated enough
            if(std::fabs(shifts[i] - shifts[j]) <= (plateThickness * 3.0f + voxelDiagonal * 2.0f) && ii != j)
            {
              // if the attempted placement of a second, third, fourth, ... transformation phase is too close to a
              // previously inserted transformation phase in the current grain, then increment the attempts counter
              // and try a different shift
              ++attempt;
              --i;
              break;
            }

            // set the origin of the plane
            d = -sampleHabitPlane[0] * (m_Centroids[3 * curFeature + 0] + shifts[i]) - sampleHabitPlane[1] * (m_Centroids[3 * curFeature + 1] + shifts[i]) -
                sampleHabitPlane[2] * (m_Centroids[3 * curFeature + 2] + shifts[i]);

            createdTransformationPhase = placeTransformationPhase(curFeature, sampleHabitPlane, totalFeatures, plateThickness, d, e);
            // int stop = 0;
            // change an isthmus transformation phase to a peninsula transformation phase
            random = static_cast<float>(rg.genrand_res53());
            if(createdTransformationPhase)
            {
              if(random < m_PeninsulaFrac)
              {
                peninsulaTransformationPhase(curFeature, totalFeatures);
              }

              // filling in transformation phase stats that already exist for parents
              tDims[0] = totalFeatures + 1;
              m->getAttributeMatrix(getFeaturePhasesArrayPath().getAttributeMatrixName())->resizeAttributeArrays(tDims);
              updateFeatureInstancePointers();
              totalFeatures = transferAttributes(totalFeatures, q2, e, curFeature);
              ++m_NumFeaturesPerParent[curFeature];
            }
          }
          // try 10 times to insert 2+ transformation phases, if it can't, the grain is probably too small
          if(attempt == 10)
          {
            break;
          }
        }
      }
    }
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
bool InsertTransformationPhases::placeTransformationPhase(int32_t curFeature, float sampleHabitPlane[3], int32_t totalFeatures, float plateThickness, float d, float* euler)
{
  DataContainer::Pointer m = getDataContainerArray()->getDataContainer(getFeatureIdsArrayPath().getDataContainerName());
  SIMPL_RANDOMNG_NEW()

  // size_t totalPoints = static_cast<size_t>(m_FeatureIdsPtr.lock()->getNumberOfTuples());
  int xPoints = static_cast<int>(m->getGeometryAs<ImageGeom>()->getXPoints());
  int yPoints = static_cast<int>(m->getGeometryAs<ImageGeom>()->getYPoints());
  int zPoints = static_cast<int>(m->getGeometryAs<ImageGeom>()->getZPoints());
  FloatVec3Type spacing = m->getGeometryAs<ImageGeom>()->getSpacing();
  bool flag = false;
  float x, y, z, D;
  bool firstVoxel = true;

  // loop through all cells to find matching grain IDs
  int zStride, yStride;
  for(int i = 0; i < zPoints; i++)
  {
    zStride = i * xPoints * yPoints;
    for(int j = 0; j < yPoints; j++)
    {
      yStride = j * xPoints;
      for(int k = 0; k < xPoints; k++)
      {
        int32_t gnum = m_FeatureIds[zStride + yStride + k];
        x = float(k) * spacing[0];
        y = float(j) * spacing[1];
        z = float(i) * spacing[2];

        // if the grain IDs match...
        if(gnum == curFeature)
        {
          // calculate the distance between the cell and the plane
          float denom = 1 / sqrtf(sampleHabitPlane[0] * sampleHabitPlane[0] + sampleHabitPlane[1] * sampleHabitPlane[1] + sampleHabitPlane[2] * sampleHabitPlane[2]);

          D = sampleHabitPlane[0] * x * denom + sampleHabitPlane[1] * y * denom + sampleHabitPlane[2] * z * denom + d * denom;

          // if the cell-plane distance is less than the plate thickness then place a transformation phase voxel
          // and that a neighboring cell is not a transformation phase cell
          if(std::fabs(D) < plateThickness)
          ////      && (k == 0 || m_FeatureIds[zStride+yStride+k-1] < numFeatures || m_FeatureIds[zStride+yStride+k-1] == totalFeatures)
          ////      && (k == (m->getGeometryAs<ImageGeom>()->getXPoints() - 1) || m_FeatureIds[zStride+yStride+k+1] < numFeatures || m_FeatureIds[zStride+yStride+k+1] == totalFeatures)
          ////      && (j == 0 || m_FeatureIds[zStride+yStride+k-xPoints] < numFeatures || m_FeatureIds[zStride+yStride+k-xPoints] == totalFeatures)
          ////      && (j == (m->getGeometryAs<ImageGeom>()->getYPoints() - 1) || m_FeatureIds[zStride+yStride+k+xPoints] < numFeatures || m_FeatureIds[zStride+yStride+k+xPoints] == totalFeatures)
          ////      && (i == 0 || m_FeatureIds[zStride+yStride+k-xPoints*yPoints] < numFeatures || m_FeatureIds[zStride+yStride+k-xPoints*yPoints] == totalFeatures)
          ////      && (i == (m->getGeometryAs<ImageGeom>()->getZPoints() - 1) || m_FeatureIds[zStride+yStride+k+xPoints*yPoints] < numFeatures || m_FeatureIds[zStride+yStride+k+xPoints*yPoints] ==
          /// totalFeatures))
          //      && (k == 0 || j == 0 || m_FeatureIds[zStride+yStride+k-xPoints-1] < numFeatures || m_FeatureIds[zStride+yStride+k-xPoints-1] == totalFeatures)
          //      && (k == (m->getGeometryAs<ImageGeom>()->getXPoints() - 1) || j == 0 || m_FeatureIds[zStride+yStride+k-xPoints+1] < numFeatures || m_FeatureIds[zStride+yStride+k-xPoints+1] ==
          //      totalFeatures)
          //      && (k == 0 || j == (m->getGeometryAs<ImageGeom>()->getYPoints() - 1) || m_FeatureIds[zStride+yStride+k+xPoints-1] < numFeatures || m_FeatureIds[zStride+yStride+k+xPoints-1] ==
          //      totalFeatures)
          //      && (k == (m->getGeometryAs<ImageGeom>()->getXPoints() - 1) || j == (m->getGeometryAs<ImageGeom>()->getYPoints() - 1) || m_FeatureIds[zStride+yStride+k+xPoints+1] < numFeatures ||
          //      m_FeatureIds[zStride+yStride+k+xPoints+1] == totalFeatures)
          //      && (k == 0 || j == 0 || i == (m->getGeometryAs<ImageGeom>()->getZPoints() - 1) || m_FeatureIds[zStride+yStride+k+xPoints*yPoints-xPoints-1] < numFeatures ||
          //      m_FeatureIds[zStride+yStride+k+xPoints*yPoints-xPoints-1] == totalFeatures)
          //      && (k == (m->getGeometryAs<ImageGeom>()->getXPoints() - 1) || j == 0 || i == (m->getGeometryAs<ImageGeom>()->getZPoints() - 1) ||
          //      m_FeatureIds[zStride+yStride+k+xPoints*yPoints-xPoints+1] < numFeatures || m_FeatureIds[zStride+yStride+k+xPoints*yPoints-xPoints+1] == totalFeatures)
          //      && (k == 0 || j == (m->getGeometryAs<ImageGeom>()->getYPoints() - 1) || i == (m->getGeometryAs<ImageGeom>()->getZPoints() - 1) ||
          //      m_FeatureIds[zStride+yStride+k+xPoints*yPoints+xPoints-1] < numFeatures || m_FeatureIds[zStride+yStride+k+xPoints*yPoints+xPoints-1] == totalFeatures)
          //      && (k == (m->getGeometryAs<ImageGeom>()->getXPoints() - 1) || j == (m->getGeometryAs<ImageGeom>()->getYPoints() - 1) || i == (m->getGeometryAs<ImageGeom>()->getZPoints() - 1) ||
          //      m_FeatureIds[zStride+yStride+k+xPoints*yPoints+xPoints+1] < numFeatures || m_FeatureIds[zStride+yStride+k+xPoints*yPoints+xPoints+1] == totalFeatures)
          //      && (k == 0 || j == 0 || i == 0 || m_FeatureIds[zStride+yStride+k-xPoints*yPoints-xPoints-1] < numFeatures || m_FeatureIds[zStride+yStride+k-xPoints*yPoints-xPoints-1] ==
          //      totalFeatures)
          //      && (k == (m->getGeometryAs<ImageGeom>()->getXPoints() - 1) || j == 0 || i == 0 || m_FeatureIds[zStride+yStride+k-xPoints*yPoints-xPoints+1] < numFeatures ||
          //      m_FeatureIds[zStride+yStride+k-xPoints*yPoints-xPoints+1] == totalFeatures)
          //      && (k == (m->getGeometryAs<ImageGeom>()->getXPoints() - 1) || j == (m->getGeometryAs<ImageGeom>()->getYPoints() - 1) || i == 0 ||
          //      m_FeatureIds[zStride+yStride+k-xPoints*yPoints+xPoints-1] < numFeatures || m_FeatureIds[zStride+yStride+k-xPoints*yPoints+xPoints-1] == totalFeatures)
          //      && (k == (m->getGeometryAs<ImageGeom>()->getXPoints() - 1) || j == (m->getGeometryAs<ImageGeom>()->getYPoints() - 1) || i == 0 ||
          //      m_FeatureIds[zStride+yStride+k-xPoints*yPoints+xPoints+1] < numFeatures || m_FeatureIds[zStride+yStride+k-xPoints*yPoints+xPoints+1] == totalFeatures))
          {
            // check if an "island" transformation phase voxel will be placed (excluding the first voxel placement)
            /*
            if ((k == 0 || m_FeatureIds[zStride+yStride+k-1] != totalFeatures)
            && (k == (xPoints - 1) || m_FeatureIds[zStride+yStride+k+1] != totalFeatures)
            && (j == 0 || m_FeatureIds[zStride+yStride+k-xPoints] != totalFeatures)
            && (j == (yPoints - 1) || m_FeatureIds[zStride+yStride+k+xPoints] != totalFeatures)
            && (i == 0 || m_FeatureIds[zStride+yStride+k-xPoints*yPoints] != totalFeatures)
            && (i == (zPoints - 1) || m_FeatureIds[zStride+yStride+k+xPoints*yPoints] != totalFeatures)
            && firstVoxel == false)
            {
            // if an "island" transformation phase voxel is inserted flip all the transformation phase voxels back to the grain ID
            // and this will go in the books as a failed transformation phase insertion attempt
            for (int64_t l = 0; l < totalPoints; ++l)
            {
            if (m_FeatureIds[l] == totalFeatures)
            {
            //        int goal = m_FeatureIds[l];
            m_FeatureIds[l] = curFeature;
            //        int one = m_FeatureIds[zStride+yStride+k+xPoints*yPoints];
            //        int two =  m_FeatureIds[zStride+yStride+k-xPoints*yPoints];
            //        int three = m_FeatureIds[zStride+yStride+k+xPoints];
            //        int four = m_FeatureIds[zStride+yStride+k-xPoints];
            //        int five = m_FeatureIds[zStride+yStride+k+1];
            //        int six = m_FeatureIds[zStride+yStride+k-1];
            int stop = 0;
            }
            }
            flag = false;
            }
            else
            {
            */
            size_t cellIndex = static_cast<size_t>(zStride + yStride + k);
            m_FeatureIds[cellIndex] = static_cast<int32_t>(totalFeatures);
            m_CellEulerAngles[cellIndex * 3] = euler[0];
            m_CellEulerAngles[cellIndex * 3 + 1] = euler[1];
            m_CellEulerAngles[cellIndex * 3 + 2] = euler[2];
            m_CellPhases[cellIndex] = static_cast<int32_t>(m_PhaseTypesPtr.lock()->getNumberOfTuples() - 1);
            flag = true;
            firstVoxel = false;
            //    }
          }
        }
      }
    }
  }
  return flag;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void InsertTransformationPhases::peninsulaTransformationPhase(int32_t curFeature, int32_t totalFeatures)
{
  DataContainer::Pointer m = getDataContainerArray()->getDataContainer(getFeatureIdsArrayPath().getDataContainerName());
  SIMPL_RANDOMNG_NEW()
  int32_t xPoints = static_cast<int>(m->getGeometryAs<ImageGeom>()->getXPoints());
  int32_t yPoints = static_cast<int>(m->getGeometryAs<ImageGeom>()->getYPoints());
  int32_t zPoints = static_cast<int>(m->getGeometryAs<ImageGeom>()->getZPoints());

  int32_t x1 = -1, x2 = -1, y1 = -1, y2 = -1, z1 = -1, z2 = -1;
  float transformationPhaseLength = 0.0f, random = 0.0f, fractionKept = 0.0f, currentDistance = 0.0f;
  int32_t zStride, yStride;

  // loop through all cells to find matching grain IDs
  for(int32_t i = 0; i < zPoints; i++)
  {
    zStride = i * xPoints * yPoints;
    for(int32_t j = 0; j < yPoints; j++)
    {
      yStride = j * xPoints;
      for(int32_t k = 0; k < xPoints; k++)
      {
        int32_t gnum = m_FeatureIds[zStride + yStride + k];
        // if the grain IDs match...
        if(gnum == totalFeatures)
        {
          if(x1 == -1)
          {
            // establish one extremum of the transformation phase
            x1 = i;
            y1 = j;
            z1 = k;
          }
          // establish the other extremum of the transformation phase
          x2 = i;
          y2 = j;
          z2 = k;
        }
      }
    }
  }

  // calculate the distance between the extrema
  transformationPhaseLength = sqrtf((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1) + (z2 - z1) * (z2 - z1));
  // if the transformation phase is 1 voxel, set the length to 2 so it makes it through the below loop
  // but won't change anything
  if(transformationPhaseLength == 0.0f)
  {
    transformationPhaseLength = 2.0f;
  }

  // choose which extremum to start from
  random = static_cast<float>(rg.genrand_res53());
  if(random < 0.5f)
  {
    x2 = x1;
    y2 = y1;
    z2 = z1;
  }

  // choose how much of the transformation phase to keep (has to be at least one voxel)
  while(static_cast<int32_t>(fractionKept * transformationPhaseLength) < 1)
  {
    fractionKept = static_cast<float>(rg.genrand_res53());
  }

  // loop through again to decide which transformation phase Ids get flipped back to grain Ids
  for(int32_t i = 0; i < zPoints; i++)
  {
    zStride = i * xPoints * yPoints;
    for(int32_t j = 0; j < yPoints; j++)
    {
      yStride = j * xPoints;
      for(int32_t k = 0; k < xPoints; k++)
      {
        int gnum = m_FeatureIds[zStride + yStride + k];
        // if the grain IDs match...
        if(gnum == totalFeatures)
        {
          currentDistance = sqrtf((i - x1) * (i - x1) + (j - y1) * (j - y1) + (k - z1) * (k - z1));
          // if the distance is longer than the transformation phase length, flip back to the parent ID
          if(currentDistance > transformationPhaseLength * fractionKept)
          {
            m_FeatureIds[zStride + yStride + k] = curFeature;
          }
        }
      }
    }
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
size_t InsertTransformationPhases::transferAttributes(size_t totalFeatures, QuatF q, float e[], int32_t curFeature)
{
  // DataContainer::Pointer m = getDataContainerArray()->getDataContainer(getFeatureIdsArrayPath().getDataContainerName());

  m_AvgQuats[4 * totalFeatures + 0] = q[0];
  m_AvgQuats[4 * totalFeatures + 1] = q[1];
  m_AvgQuats[4 * totalFeatures + 2] = q[2];
  m_AvgQuats[4 * totalFeatures + 3] = q[3];
  m_Centroids[3 * totalFeatures + 0] = 1.0f;
  m_Centroids[3 * totalFeatures + 1] = 2.0f;
  m_Centroids[3 * totalFeatures + 2] = 3.0f;
  m_FeatureEulerAngles[3 * totalFeatures + 0] = e[0];
  m_FeatureEulerAngles[3 * totalFeatures + 1] = e[1];
  m_FeatureEulerAngles[3 * totalFeatures + 2] = e[2];
  m_EquivalentDiameters[totalFeatures] = 4.0f;
  m_FeaturePhases[totalFeatures] = static_cast<int32_t>(m_PhaseTypesPtr.lock()->getNumberOfTuples() - 1);
  m_FeatureParentIds[totalFeatures] = curFeature;
  m_NumFeaturesPerParent[totalFeatures] = 0;
  return ++totalFeatures;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AbstractFilter::Pointer InsertTransformationPhases::newFilterInstance(bool copyFilterParameters) const
{
  /*
   * ApplyToVoxelVolume
   * ApplyToSurfaceMesh
   * Origin
   */
  InsertTransformationPhases::Pointer filter = InsertTransformationPhases::New();
  if(copyFilterParameters)
  {
    copyFilterParameterInstanceVariables(filter.get());
  }
  return filter;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString InsertTransformationPhases::getCompiledLibraryName() const
{
  return DREAM3DReviewConstants::DREAM3DReviewBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString InsertTransformationPhases::getBrandingString() const
{
  return TransformationPhaseConstants::TransformationPhasePluginDisplayName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString InsertTransformationPhases::getFilterVersion() const
{
  QString version;
  QTextStream vStream(&version);
  vStream << DREAM3DReview::Version::Major() << "." << DREAM3DReview::Version::Minor() << "." << DREAM3DReview::Version::Patch();
  return version;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString InsertTransformationPhases::getGroupName() const
{
  return SIMPL::FilterGroups::Unsupported;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QUuid InsertTransformationPhases::getUuid() const
{
  return QUuid("{fd6da27d-0f2c-5c35-802f-7d6ce1ad0ca1}");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString InsertTransformationPhases::getSubGroupName() const
{
  return SIMPL::FilterSubGroups::PackingFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString InsertTransformationPhases::getHumanLabel() const
{
  return "Insert Transformation Phases";
}

// -----------------------------------------------------------------------------
InsertTransformationPhases::Pointer InsertTransformationPhases::NullPointer()
{
  return Pointer(static_cast<Self*>(nullptr));
}

// -----------------------------------------------------------------------------
std::shared_ptr<InsertTransformationPhases> InsertTransformationPhases::New()
{
  struct make_shared_enabler : public InsertTransformationPhases
  {
  };
  std::shared_ptr<make_shared_enabler> val = std::make_shared<make_shared_enabler>();
  val->setupFilterParameters();
  return val;
}

// -----------------------------------------------------------------------------
QString InsertTransformationPhases::getNameOfClass() const
{
  return QString("InsertTransformationPhases");
}

// -----------------------------------------------------------------------------
QString InsertTransformationPhases::ClassName()
{
  return QString("InsertTransformationPhases");
}

// -----------------------------------------------------------------------------
void InsertTransformationPhases::setParentPhase(int value)
{
  m_ParentPhase = value;
}

// -----------------------------------------------------------------------------
int InsertTransformationPhases::getParentPhase() const
{
  return m_ParentPhase;
}

// -----------------------------------------------------------------------------
void InsertTransformationPhases::setTransCrystalStruct(int value)
{
  m_TransCrystalStruct = value;
}

// -----------------------------------------------------------------------------
int InsertTransformationPhases::getTransCrystalStruct() const
{
  return m_TransCrystalStruct;
}

// -----------------------------------------------------------------------------
void InsertTransformationPhases::setTransformationPhaseMisorientation(float value)
{
  m_TransformationPhaseMisorientation = value;
}

// -----------------------------------------------------------------------------
float InsertTransformationPhases::getTransformationPhaseMisorientation() const
{
  return m_TransformationPhaseMisorientation;
}

// -----------------------------------------------------------------------------
void InsertTransformationPhases::setTransformationPhaseHabitPlane(const FloatVec3Type& value)
{
  m_TransformationPhaseHabitPlane = value;
}

// -----------------------------------------------------------------------------
FloatVec3Type InsertTransformationPhases::getTransformationPhaseHabitPlane() const
{
  return m_TransformationPhaseHabitPlane;
}

// -----------------------------------------------------------------------------
void InsertTransformationPhases::setDefineHabitPlane(bool value)
{
  m_DefineHabitPlane = value;
}

// -----------------------------------------------------------------------------
bool InsertTransformationPhases::getDefineHabitPlane() const
{
  return m_DefineHabitPlane;
}

// -----------------------------------------------------------------------------
void InsertTransformationPhases::setUseAllVariants(bool value)
{
  m_UseAllVariants = value;
}

// -----------------------------------------------------------------------------
bool InsertTransformationPhases::getUseAllVariants() const
{
  return m_UseAllVariants;
}

// -----------------------------------------------------------------------------
void InsertTransformationPhases::setCoherentFrac(float value)
{
  m_CoherentFrac = value;
}

// -----------------------------------------------------------------------------
float InsertTransformationPhases::getCoherentFrac() const
{
  return m_CoherentFrac;
}

// -----------------------------------------------------------------------------
void InsertTransformationPhases::setTransformationPhaseThickness(float value)
{
  m_TransformationPhaseThickness = value;
}

// -----------------------------------------------------------------------------
float InsertTransformationPhases::getTransformationPhaseThickness() const
{
  return m_TransformationPhaseThickness;
}

// -----------------------------------------------------------------------------
void InsertTransformationPhases::setNumTransformationPhasesPerFeature(int value)
{
  m_NumTransformationPhasesPerFeature = value;
}

// -----------------------------------------------------------------------------
int InsertTransformationPhases::getNumTransformationPhasesPerFeature() const
{
  return m_NumTransformationPhasesPerFeature;
}

// -----------------------------------------------------------------------------
void InsertTransformationPhases::setPeninsulaFrac(float value)
{
  m_PeninsulaFrac = value;
}

// -----------------------------------------------------------------------------
float InsertTransformationPhases::getPeninsulaFrac() const
{
  return m_PeninsulaFrac;
}

// -----------------------------------------------------------------------------
void InsertTransformationPhases::setStatsGenCellEnsembleAttributeMatrixPath(const DataArrayPath& value)
{
  m_StatsGenCellEnsembleAttributeMatrixPath = value;
}

// -----------------------------------------------------------------------------
DataArrayPath InsertTransformationPhases::getStatsGenCellEnsembleAttributeMatrixPath() const
{
  return m_StatsGenCellEnsembleAttributeMatrixPath;
}

// -----------------------------------------------------------------------------
void InsertTransformationPhases::setCellFeatureAttributeMatrixName(const DataArrayPath& value)
{
  m_CellFeatureAttributeMatrixName = value;
}

// -----------------------------------------------------------------------------
DataArrayPath InsertTransformationPhases::getCellFeatureAttributeMatrixName() const
{
  return m_CellFeatureAttributeMatrixName;
}

// -----------------------------------------------------------------------------
void InsertTransformationPhases::setFeatureIdsArrayPath(const DataArrayPath& value)
{
  m_FeatureIdsArrayPath = value;
}

// -----------------------------------------------------------------------------
DataArrayPath InsertTransformationPhases::getFeatureIdsArrayPath() const
{
  return m_FeatureIdsArrayPath;
}

// -----------------------------------------------------------------------------
void InsertTransformationPhases::setCellEulerAnglesArrayPath(const DataArrayPath& value)
{
  m_CellEulerAnglesArrayPath = value;
}

// -----------------------------------------------------------------------------
DataArrayPath InsertTransformationPhases::getCellEulerAnglesArrayPath() const
{
  return m_CellEulerAnglesArrayPath;
}

// -----------------------------------------------------------------------------
void InsertTransformationPhases::setCellPhasesArrayPath(const DataArrayPath& value)
{
  m_CellPhasesArrayPath = value;
}

// -----------------------------------------------------------------------------
DataArrayPath InsertTransformationPhases::getCellPhasesArrayPath() const
{
  return m_CellPhasesArrayPath;
}

// -----------------------------------------------------------------------------
void InsertTransformationPhases::setAvgQuatsArrayPath(const DataArrayPath& value)
{
  m_AvgQuatsArrayPath = value;
}

// -----------------------------------------------------------------------------
DataArrayPath InsertTransformationPhases::getAvgQuatsArrayPath() const
{
  return m_AvgQuatsArrayPath;
}

// -----------------------------------------------------------------------------
void InsertTransformationPhases::setCentroidsArrayPath(const DataArrayPath& value)
{
  m_CentroidsArrayPath = value;
}

// -----------------------------------------------------------------------------
DataArrayPath InsertTransformationPhases::getCentroidsArrayPath() const
{
  return m_CentroidsArrayPath;
}

// -----------------------------------------------------------------------------
void InsertTransformationPhases::setEquivalentDiametersArrayPath(const DataArrayPath& value)
{
  m_EquivalentDiametersArrayPath = value;
}

// -----------------------------------------------------------------------------
DataArrayPath InsertTransformationPhases::getEquivalentDiametersArrayPath() const
{
  return m_EquivalentDiametersArrayPath;
}

// -----------------------------------------------------------------------------
void InsertTransformationPhases::setFeatureEulerAnglesArrayPath(const DataArrayPath& value)
{
  m_FeatureEulerAnglesArrayPath = value;
}

// -----------------------------------------------------------------------------
DataArrayPath InsertTransformationPhases::getFeatureEulerAnglesArrayPath() const
{
  return m_FeatureEulerAnglesArrayPath;
}

// -----------------------------------------------------------------------------
void InsertTransformationPhases::setFeaturePhasesArrayPath(const DataArrayPath& value)
{
  m_FeaturePhasesArrayPath = value;
}

// -----------------------------------------------------------------------------
DataArrayPath InsertTransformationPhases::getFeaturePhasesArrayPath() const
{
  return m_FeaturePhasesArrayPath;
}

// -----------------------------------------------------------------------------
void InsertTransformationPhases::setFeatureParentIdsArrayName(const QString& value)
{
  m_FeatureParentIdsArrayName = value;
}

// -----------------------------------------------------------------------------
QString InsertTransformationPhases::getFeatureParentIdsArrayName() const
{
  return m_FeatureParentIdsArrayName;
}

// -----------------------------------------------------------------------------
void InsertTransformationPhases::setNumFeaturesPerParentArrayPath(const DataArrayPath& value)
{
  m_NumFeaturesPerParentArrayPath = value;
}

// -----------------------------------------------------------------------------
DataArrayPath InsertTransformationPhases::getNumFeaturesPerParentArrayPath() const
{
  return m_NumFeaturesPerParentArrayPath;
}

// -----------------------------------------------------------------------------
void InsertTransformationPhases::setCrystalStructuresArrayPath(const DataArrayPath& value)
{
  m_CrystalStructuresArrayPath = value;
}

// -----------------------------------------------------------------------------
DataArrayPath InsertTransformationPhases::getCrystalStructuresArrayPath() const
{
  return m_CrystalStructuresArrayPath;
}

// -----------------------------------------------------------------------------
void InsertTransformationPhases::setPhaseTypesArrayPath(const DataArrayPath& value)
{
  m_PhaseTypesArrayPath = value;
}

// -----------------------------------------------------------------------------
DataArrayPath InsertTransformationPhases::getPhaseTypesArrayPath() const
{
  return m_PhaseTypesArrayPath;
}

// -----------------------------------------------------------------------------
void InsertTransformationPhases::setShapeTypesArrayPath(const DataArrayPath& value)
{
  m_ShapeTypesArrayPath = value;
}

// -----------------------------------------------------------------------------
DataArrayPath InsertTransformationPhases::getShapeTypesArrayPath() const
{
  return m_ShapeTypesArrayPath;
}

// -----------------------------------------------------------------------------
void InsertTransformationPhases::setNumFeaturesArrayPath(const DataArrayPath& value)
{
  m_NumFeaturesArrayPath = value;
}

// -----------------------------------------------------------------------------
DataArrayPath InsertTransformationPhases::getNumFeaturesArrayPath() const
{
  return m_NumFeaturesArrayPath;
}
