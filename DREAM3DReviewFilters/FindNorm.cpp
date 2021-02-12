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
#include "FindNorm.h"

#include <cmath>

#include <QtCore/QTextStream>

#include "SIMPLib/Common/Constants.h"
#include "SIMPLib/Common/TemplateHelpers.h"
#include "SIMPLib/DataContainers/DataContainerArray.h"
#include "SIMPLib/FilterParameters/AbstractFilterParametersReader.h"
#include "SIMPLib/FilterParameters/DataArrayCreationFilterParameter.h"
#include "SIMPLib/FilterParameters/DataArraySelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/FloatFilterParameter.h"

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
FindNorm::FindNorm() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
FindNorm::~FindNorm() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FindNorm::setupFilterParameters()
{
  FilterParameterVectorType parameters;
  parameters.push_back(SIMPL_NEW_FLOAT_FP("p-Space Value", PSpace, FilterParameter::Category::Parameter, FindNorm));
  DataArraySelectionFilterParameter::RequirementType dasReq =
      DataArraySelectionFilterParameter::CreateRequirement(SIMPL::Defaults::AnyPrimitive, SIMPL::Defaults::AnyComponentSize, AttributeMatrix::Type::Any, IGeometry::Type::Any);
  parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Input Attribute Array", SelectedArrayPath, FilterParameter::Category::RequiredArray, FindNorm, dasReq));
  DataArrayCreationFilterParameter::RequirementType dacReq = DataArrayCreationFilterParameter::CreateRequirement(AttributeMatrix::Type::Any, IGeometry::Type::Any);
  parameters.push_back(SIMPL_NEW_DA_CREATION_FP("Norm", NormArrayPath, FilterParameter::Category::CreatedArray, FindNorm, dacReq));
  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FindNorm::readFilterParameters(AbstractFilterParametersReader* reader, int index)
{
  reader->openFilterGroup(this, index);
  setSelectedArrayPath(reader->readDataArrayPath("SelectedArrayPath", getSelectedArrayPath()));
  setNormArrayPath(reader->readDataArrayPath("NormArrayPath", getNormArrayPath()));
  setPSpace(reader->readValue("PSpace", getPSpace()));
  reader->closeFilterGroup();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FindNorm::initialize()
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FindNorm::dataCheck()
{
  clearErrorCode();
  clearWarningCode();

  if(getPSpace() < 0)
  {
    QString ss = QObject::tr("p-space value must be greater than or equal to 0");
    setErrorCondition(-11002, ss);
  }

  QVector<DataArrayPath> dataArrayPaths;

  m_InArrayPtr = getDataContainerArray()->getPrereqIDataArrayFromPath(this, getSelectedArrayPath());
  if(getErrorCode() >= 0)
  {
    dataArrayPaths.push_back(getSelectedArrayPath());
  }

  std::vector<size_t> cDims(1, 1);

  m_NormPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<float>>(this, getNormArrayPath(), 0, cDims, "", DataArrayID31);
  if(m_NormPtr.lock())
  {
    m_Norm = m_NormPtr.lock()->getPointer(0);
  }
  if(getErrorCode() >= 0)
  {
    dataArrayPaths.push_back(getNormArrayPath());
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
template <typename T>
void findPthNorm(IDataArray::Pointer inDataPtr, const FloatArrayType::Pointer& normPtr, float p)
{
  typename DataArray<T>::Pointer inputDataPtr = std::dynamic_pointer_cast<DataArray<T>>(inDataPtr);
  T* inData = inputDataPtr->getPointer(0);
  float* norm = normPtr->getPointer(0);

  int32_t nDims = inputDataPtr->getNumberOfComponents();
  size_t nTuples = inputDataPtr->getNumberOfTuples();

  for(size_t i = 0; i < nTuples; i++)
  {
    float normTmp = 0.0f;
    //    float valTmp = 0.0f;
    for(int32_t j = 0; j < nDims; j++)
    {
      normTmp += std::pow(static_cast<float>(inData[nDims * i + j]), p);
    }
    norm[i] = std::pow(normTmp, 1 / p);
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FindNorm::execute()
{
  dataCheck();
  if(getErrorCode() < 0)
  {
    return;
  }

  EXECUTE_FUNCTION_TEMPLATE(this, findPthNorm, m_InArrayPtr.lock(), m_InArrayPtr.lock(), m_NormPtr.lock(), m_PSpace);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AbstractFilter::Pointer FindNorm::newFilterInstance(bool copyFilterParameters) const
{
  FindNorm::Pointer filter = FindNorm::New();
  if(copyFilterParameters)
  {
    copyFilterParameterInstanceVariables(filter.get());
  }
  return filter;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString FindNorm::getCompiledLibraryName() const
{
  return DREAM3DReviewConstants::DREAM3DReviewBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString FindNorm::getBrandingString() const
{
  return "DREAM3DReview";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString FindNorm::getFilterVersion() const
{
  QString version;
  QTextStream vStream(&version);
  vStream << DREAM3DReview::Version::Major() << "." << DREAM3DReview::Version::Minor() << "." << DREAM3DReview::Version::Patch();
  return version;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString FindNorm::getGroupName() const
{
  return DREAM3DReviewConstants::FilterGroups::DREAM3DReviewFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QUuid FindNorm::getUuid() const
{
  return QUuid("{5d0cd577-3e3e-57b8-a36d-b215b834251f}");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString FindNorm::getSubGroupName() const
{
  return DREAM3DReviewConstants::FilterSubGroups::StatisticsFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString FindNorm::getHumanLabel() const
{
  return "Find Norm";
}

// -----------------------------------------------------------------------------
FindNorm::Pointer FindNorm::NullPointer()
{
  return Pointer(static_cast<Self*>(nullptr));
}

// -----------------------------------------------------------------------------
std::shared_ptr<FindNorm> FindNorm::New()
{
  struct make_shared_enabler : public FindNorm
  {
  };
  std::shared_ptr<make_shared_enabler> val = std::make_shared<make_shared_enabler>();
  val->setupFilterParameters();
  return val;
}

// -----------------------------------------------------------------------------
QString FindNorm::getNameOfClass() const
{
  return QString("FindNorm");
}

// -----------------------------------------------------------------------------
QString FindNorm::ClassName()
{
  return QString("FindNorm");
}

// -----------------------------------------------------------------------------
void FindNorm::setSelectedArrayPath(const DataArrayPath& value)
{
  m_SelectedArrayPath = value;
}

// -----------------------------------------------------------------------------
DataArrayPath FindNorm::getSelectedArrayPath() const
{
  return m_SelectedArrayPath;
}

// -----------------------------------------------------------------------------
void FindNorm::setNormArrayPath(const DataArrayPath& value)
{
  m_NormArrayPath = value;
}

// -----------------------------------------------------------------------------
DataArrayPath FindNorm::getNormArrayPath() const
{
  return m_NormArrayPath;
}

// -----------------------------------------------------------------------------
void FindNorm::setPSpace(float value)
{
  m_PSpace = value;
}

// -----------------------------------------------------------------------------
float FindNorm::getPSpace() const
{
  return m_PSpace;
}
