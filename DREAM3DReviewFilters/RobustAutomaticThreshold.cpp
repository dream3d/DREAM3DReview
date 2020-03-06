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
#include "RobustAutomaticThreshold.h"

#include <memory>

#include <QtCore/QTextStream>

#include "SIMPLib/Common/Constants.h"
#include "SIMPLib/Common/TemplateHelpers.h"
#include "SIMPLib/FilterParameters/AbstractFilterParametersReader.h"
#include "SIMPLib/FilterParameters/DataArrayCreationFilterParameter.h"
#include "SIMPLib/FilterParameters/DataArraySelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/IntFilterParameter.h"
#include "SIMPLib/DataContainers/DataContainerArray.h"
#include "SIMPLib/ITK/Dream3DTemplateAliasMacro.h"

#include "DREAM3DReview/DREAM3DReviewConstants.h"
#include "DREAM3DReview/DREAM3DReviewVersion.h"

/* Create Enumerations to allow the created Attribute Arrays to take part in renaming */
enum createdPathID : RenameDataPath::DataID_t
{
  DataArrayID30 = 30,
  DataArrayID31 = 31,
};

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
RobustAutomaticThreshold::RobustAutomaticThreshold()
: m_InputArrayPath("", "", "")
, m_FeatureIdsArrayPath("", "", "Mask")
, m_GradientMagnitudeArrayPath("", "", "")
{
  initialize();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
RobustAutomaticThreshold::~RobustAutomaticThreshold() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void RobustAutomaticThreshold::setupFilterParameters()
{
  FilterParameterVectorType parameters;
  DataArraySelectionFilterParameter::RequirementType dasReq = DataArraySelectionFilterParameter::CreateRequirement(SIMPL::Defaults::AnyPrimitive, 1, AttributeMatrix::Type::Any, IGeometry::Type::Any);
  parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Attribute Array to Threshold", InputArrayPath, FilterParameter::RequiredArray, RobustAutomaticThreshold, dasReq));
  dasReq = DataArraySelectionFilterParameter::CreateRequirement(SIMPL::TypeNames::Float, 1, AttributeMatrix::Type::Any, IGeometry::Type::Any);
  parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Gradient Magnitude", GradientMagnitudeArrayPath, FilterParameter::RequiredArray, RobustAutomaticThreshold, dasReq));
  DataArrayCreationFilterParameter::RequirementType dacReq = DataArrayCreationFilterParameter::CreateRequirement(AttributeMatrix::Type::Any, IGeometry::Type::Any);
  parameters.push_back(SIMPL_NEW_DA_CREATION_FP("Mask", FeatureIdsArrayPath, FilterParameter::RequiredArray, RobustAutomaticThreshold, dacReq));
  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void RobustAutomaticThreshold::readFilterParameters(AbstractFilterParametersReader* reader, int index)
{
  reader->openFilterGroup(this, index);
  setInputArrayPath(reader->readDataArrayPath("InputArrayPath", getInputArrayPath()));
  setFeatureIdsArrayPath(reader->readDataArrayPath("FeatureIdsArrayPath", getFeatureIdsArrayPath()));
  setGradientMagnitudeArrayPath(reader->readDataArrayPath("GradientMagnitudeArrayPath", getGradientMagnitudeArrayPath()));
  reader->closeFilterGroup();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void RobustAutomaticThreshold::initialize()
{
  clearErrorCode();
  clearWarningCode();
  setCancel(false);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void RobustAutomaticThreshold::dataCheck()
{
  clearErrorCode();
  clearWarningCode();

  QVector<DataArrayPath> dataArrayPaths;

  std::vector<size_t> cDims(1, 1);

  m_InputArrayPtr = getDataContainerArray()->getPrereqIDataArrayFromPath(this, getInputArrayPath());
  if(getErrorCode() < 0)
  {
    return;
  }

  if(m_InputArrayPtr.lock()->getTypeAsString() == SIMPL::TypeNames::Bool)
  {
    QString ss = QObject::tr("Input Attribute Array to threshold cannot be of type bool");
    setErrorCondition(-11001, ss);
  }
  if(getErrorCode() >= 0)
  {
    dataArrayPaths.push_back(getInputArrayPath());
  }

  m_GradientMagnitudePtr = getDataContainerArray()->getPrereqArrayFromPath<DataArray<float>>(this, getGradientMagnitudeArrayPath(), cDims);
  if(m_GradientMagnitudePtr.lock())
  {
    m_GradientMagnitude = m_GradientMagnitudePtr.lock()->getPointer(0);
  }
  if(getErrorCode() >= 0)
  {
    dataArrayPaths.push_back(getGradientMagnitudeArrayPath());
  }

  m_FeatureIdsPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<bool>>(this, getFeatureIdsArrayPath(), false, cDims, "", DataArrayID31);
  if(m_FeatureIdsPtr.lock())
  {
    m_FeatureIds = m_FeatureIdsPtr.lock()->getPointer(0);
  }
  if(getErrorCode() >= 0)
  {
    dataArrayPaths.push_back(getFeatureIdsArrayPath());
  }

  getDataContainerArray()->validateNumberOfTuples(this, dataArrayPaths);
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
template <typename T> void findThreshold(IDataArray::Pointer inputPtr, const FloatArrayType::Pointer& gradMagPtr, const BoolArrayType::Pointer& maskPtr)
{
  typename DataArray<T>::Pointer input = std::dynamic_pointer_cast<DataArray<T>>(inputPtr);
  T* iPtr = input->getPointer(0);
  float* gradMag = gradMagPtr->getPointer(0);
  bool* mask = maskPtr->getPointer(0);

  size_t numTuples = input->getNumberOfTuples();
  float numerator = 0;
  float denominator = 0;

  for(size_t i = 0; i < numTuples; i++)
  {
    numerator += (iPtr[i] * gradMag[i]);
    denominator += gradMag[i];
  }

  float threshold = numerator / denominator;

  for(size_t i = 0; i < numTuples; i++)
  {
    if(iPtr[i] < threshold)
    {
      mask[i] = false;
    }
    else
    {
      mask[i] = true;
    }
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void RobustAutomaticThreshold::execute()
{
  initialize();
  dataCheck();
  if(getErrorCode() < 0)
  {
    return;
  }

  //float threshold = 0.0f;

  EXECUTE_FUNCTION_TEMPLATE_NO_BOOL(this, findThreshold, m_InputArrayPtr.lock(), m_InputArrayPtr.lock(), m_GradientMagnitudePtr.lock(), m_FeatureIdsPtr.lock());

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AbstractFilter::Pointer RobustAutomaticThreshold::newFilterInstance(bool copyFilterParameters) const
{
  RobustAutomaticThreshold::Pointer filter = RobustAutomaticThreshold::New();
  if(copyFilterParameters)
  {
    copyFilterParameterInstanceVariables(filter.get());
  }
  return filter;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString RobustAutomaticThreshold::getCompiledLibraryName() const
{
  return DREAM3DReviewConstants::DREAM3DReviewBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString RobustAutomaticThreshold::getBrandingString() const
{
  return "DREAM3DReview";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString RobustAutomaticThreshold::getFilterVersion() const
{
  QString version;
  QTextStream vStream(&version);
  vStream << DREAM3DReview::Version::Major() << "." << DREAM3DReview::Version::Minor() << "." << DREAM3DReview::Version::Patch();
  return version;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString RobustAutomaticThreshold::getGroupName() const
{
  return DREAM3DReviewConstants::FilterGroups::DREAM3DReviewFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QUuid RobustAutomaticThreshold::getUuid() const
{
  return QUuid("{3062fc2c-76b2-5c50-92b7-edbbb424c41d}");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString RobustAutomaticThreshold::getSubGroupName() const
{
  return DREAM3DReviewConstants::FilterSubGroups::ThresholdFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString RobustAutomaticThreshold::getHumanLabel() const
{
  return "Robust Automatic Threshold";
}

// -----------------------------------------------------------------------------
RobustAutomaticThreshold::Pointer RobustAutomaticThreshold::NullPointer()
{
  return Pointer(static_cast<Self*>(nullptr));
}

// -----------------------------------------------------------------------------
std::shared_ptr<RobustAutomaticThreshold> RobustAutomaticThreshold::New()
{
  struct make_shared_enabler : public RobustAutomaticThreshold
  {
  };
  std::shared_ptr<make_shared_enabler> val = std::make_shared<make_shared_enabler>();
  val->setupFilterParameters();
  return val;
}

// -----------------------------------------------------------------------------
QString RobustAutomaticThreshold::getNameOfClass() const
{
  return QString("RobustAutomaticThreshold");
}

// -----------------------------------------------------------------------------
QString RobustAutomaticThreshold::ClassName()
{
  return QString("RobustAutomaticThreshold");
}

// -----------------------------------------------------------------------------
void RobustAutomaticThreshold::setInputArrayPath(const DataArrayPath& value)
{
  m_InputArrayPath = value;
}

// -----------------------------------------------------------------------------
DataArrayPath RobustAutomaticThreshold::getInputArrayPath() const
{
  return m_InputArrayPath;
}

// -----------------------------------------------------------------------------
void RobustAutomaticThreshold::setFeatureIdsArrayPath(const DataArrayPath& value)
{
  m_FeatureIdsArrayPath = value;
}

// -----------------------------------------------------------------------------
DataArrayPath RobustAutomaticThreshold::getFeatureIdsArrayPath() const
{
  return m_FeatureIdsArrayPath;
}

// -----------------------------------------------------------------------------
void RobustAutomaticThreshold::setGradientMagnitudeArrayPath(const DataArrayPath& value)
{
  m_GradientMagnitudeArrayPath = value;
}

// -----------------------------------------------------------------------------
DataArrayPath RobustAutomaticThreshold::getGradientMagnitudeArrayPath() const
{
  return m_GradientMagnitudeArrayPath;
}
