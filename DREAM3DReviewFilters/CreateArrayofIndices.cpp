/*
 * Your License or Copyright can go here
 */

#include "CreateArrayofIndices.h"

#include <numeric>

#include <QtCore/QTextStream>

#include "SIMPLib/Common/Constants.h"
#include "SIMPLib/FilterParameters/DataArrayCreationFilterParameter.h"

#include "DREAM3DReview/DREAM3DReviewConstants.h"
#include "DREAM3DReview/DREAM3DReviewVersion.h"
#include "SIMPLib/DataContainers/DataContainer.h"
#include "SIMPLib/DataContainers/DataContainerArray.h"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
CreateArrayofIndices::CreateArrayofIndices()
: m_IndexArrayPath("IndexArray")
{
  initialize();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
CreateArrayofIndices::~CreateArrayofIndices() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void CreateArrayofIndices::initialize()
{
  clearErrorCode();
  clearWarningCode();
  setCancel(false);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void CreateArrayofIndices::setupFilterParameters()
{
  FilterParameterVectorType parameters;
  DataArrayCreationFilterParameter::RequirementType dacReq;
  parameters.push_back(SIMPL_NEW_DA_CREATION_FP("Index Array Path", IndexArrayPath, FilterParameter::CreatedArray, CreateArrayofIndices, dacReq));
  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void CreateArrayofIndices::dataCheck()
{
  clearErrorCode();
  clearWarningCode();

  std::vector<size_t> cDims(1, 1);
  m_IndicesPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<size_t>>(this, m_IndexArrayPath, 0, cDims);
  if(nullptr != m_IndicesPtr.lock())
  {
    m_Indices = m_IndicesPtr.lock()->getPointer(0);
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void CreateArrayofIndices::execute()
{
  initialize();
  dataCheck();
  if(getErrorCode() < 0)
  {
    return;
  }

  if(getCancel())
  {
    return;
  }

  SizeTArrayType& array = *(m_IndicesPtr.lock());
  std::iota(array.begin(), array.end(), 0);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AbstractFilter::Pointer CreateArrayofIndices::newFilterInstance(bool copyFilterParameters) const
{
  CreateArrayofIndices::Pointer filter = CreateArrayofIndices::New();
  if(copyFilterParameters)
  {
    copyFilterParameterInstanceVariables(filter.get());
  }
  return filter;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString CreateArrayofIndices::getCompiledLibraryName() const
{
  return DREAM3DReviewConstants::DREAM3DReviewBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString CreateArrayofIndices::getBrandingString() const
{
  return "DREAM3DReview";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString CreateArrayofIndices::getFilterVersion() const
{
  QString version;
  QTextStream vStream(&version);
  vStream << DREAM3DReview::Version::Major() << "." << DREAM3DReview::Version::Minor() << "." << DREAM3DReview::Version::Patch();
  return version;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString CreateArrayofIndices::getGroupName() const
{
  return SIMPL::FilterGroups::Unsupported;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString CreateArrayofIndices::getSubGroupName() const
{
  return "DREAM3DReview";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString CreateArrayofIndices::getHumanLabel() const
{
  return "Create Array of Indices";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QUuid CreateArrayofIndices::getUuid() const
{
  return QUuid("{3dd5d1cb-d39d-5bb6-aab5-4f64c0892e24}");
}

// -----------------------------------------------------------------------------
CreateArrayofIndices::Pointer CreateArrayofIndices::NullPointer()
{
  return Pointer(static_cast<Self*>(nullptr));
}

// -----------------------------------------------------------------------------
std::shared_ptr<CreateArrayofIndices> CreateArrayofIndices::New()
{
  struct make_shared_enabler : public CreateArrayofIndices
  {
  };
  std::shared_ptr<make_shared_enabler> val = std::make_shared<make_shared_enabler>();
  val->setupFilterParameters();
  return val;
}

// -----------------------------------------------------------------------------
QString CreateArrayofIndices::getNameOfClass() const
{
  return QString("CreateArrayofIndices");
}

// -----------------------------------------------------------------------------
QString CreateArrayofIndices::ClassName()
{
  return QString("CreateArrayofIndices");
}

// -----------------------------------------------------------------------------
void CreateArrayofIndices::setIndexArrayPath(const DataArrayPath& value)
{
  m_IndexArrayPath = value;
}
// -----------------------------------------------------------------------------
DataArrayPath CreateArrayofIndices::getIndexArrayPath() const
{
  return m_IndexArrayPath;
}
