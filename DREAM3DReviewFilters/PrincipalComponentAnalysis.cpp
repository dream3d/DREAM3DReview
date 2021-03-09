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
#include "PrincipalComponentAnalysis.h"

#include <Eigen/Core>
#include <Eigen/Eigen>

#include <QtCore/QTextStream>

#include "SIMPLib/Common/Constants.h"
#include "SIMPLib/Common/TemplateHelpers.h"
#include "SIMPLib/DataArrays/IDataArray.h"
#include "SIMPLib/DataContainers/DataContainer.h"
#include "SIMPLib/DataContainers/DataContainerArray.h"
#include "SIMPLib/FilterParameters/AbstractFilterParametersReader.h"
#include "SIMPLib/FilterParameters/BooleanFilterParameter.h"
#include "SIMPLib/FilterParameters/ChoiceFilterParameter.h"
#include "SIMPLib/FilterParameters/DataArrayCreationFilterParameter.h"
#include "SIMPLib/FilterParameters/IntFilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedBooleanFilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedPathCreationFilterParameter.h"
#include "SIMPLib/FilterParameters/MultiDataArraySelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/SeparatorFilterParameter.h"
#include "SIMPLib/FilterParameters/StringFilterParameter.h"

#include "DREAM3DReview/DREAM3DReviewConstants.h"
#include "DREAM3DReview/DREAM3DReviewVersion.h"

/* Create Enumerations to allow the created Attribute Arrays to take part in renaming */
enum createdPathID : RenameDataPath::DataID_t
{
  AttributeMatrixID21 = 21,

  DataArrayID30 = 30,
  DataArrayID31 = 31,
  DataArrayID32 = 32,
  DataArrayID33 = 33,
};

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
PrincipalComponentAnalysis::PrincipalComponentAnalysis() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
PrincipalComponentAnalysis::~PrincipalComponentAnalysis() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void PrincipalComponentAnalysis::setupFilterParameters()
{
  FilterParameterVectorType parameters;
  MultiDataArraySelectionFilterParameter::RequirementType mdaReq =
      MultiDataArraySelectionFilterParameter::CreateRequirement(SIMPL::Defaults::AnyPrimitive, 1, AttributeMatrix::Type::Any, IGeometry::Type::Any);
  parameters.push_back(
      SIMPL_NEW_MDA_SELECTION_FP("Attribute Arrays for Computing Principal Components", SelectedDataArrayPaths, FilterParameter::Category::RequiredArray, PrincipalComponentAnalysis, mdaReq));
  {
    ChoiceFilterParameter::Pointer choices = ChoiceFilterParameter::New();
    choices->setHumanLabel("Matrix Approach");
    choices->setPropertyName("MatrixApproach");
    choices->setSetterCallback(SIMPL_BIND_SETTER(PrincipalComponentAnalysis, this, MatrixApproach));
    choices->setGetterCallback(SIMPL_BIND_GETTER(PrincipalComponentAnalysis, this, MatrixApproach));
    std::vector<QString> approaches = {"Correlation", "Covariance"};
    choices->setChoices(approaches);
    choices->setCategory(FilterParameter::Category::Parameter);
    parameters.push_back(choices);
  }
  std::vector<QString> linkedProps = {"NumberOfDimensionsForProjection", "ProjectedDataSpaceArrayPath"};
  parameters.push_back(SIMPL_NEW_LINKED_BOOL_FP("Project Data Space", ProjectDataSpace, FilterParameter::Category::Parameter, PrincipalComponentAnalysis, linkedProps));
  parameters.push_back(SIMPL_NEW_INTEGER_FP("Number of Dimensions for Projection", NumberOfDimensionsForProjection, FilterParameter::Category::Parameter, PrincipalComponentAnalysis));
  DataArrayCreationFilterParameter::RequirementType dacReq = DataArrayCreationFilterParameter::CreateRequirement(AttributeMatrix::Type::Any, IGeometry::Type::Any);
  parameters.push_back(SIMPL_NEW_DA_CREATION_FP("Projected Data Space", ProjectedDataSpaceArrayPath, FilterParameter::Category::CreatedArray, PrincipalComponentAnalysis, dacReq));
  parameters.push_back(SeparatorFilterParameter::Create("Principal Component Data", FilterParameter::Category::CreatedArray));
  parameters.push_back(SIMPL_NEW_STRING_FP("Principal Component Attribute Matrix", PCAttributeMatrixName, FilterParameter::Category::CreatedArray, PrincipalComponentAnalysis));
  parameters.push_back(SIMPL_NEW_STRING_FP("Principal Component Eigenvalues", PCEigenvaluesName, FilterParameter::Category::CreatedArray, PrincipalComponentAnalysis));
  parameters.push_back(SIMPL_NEW_STRING_FP("Principal Component Eigenvectors", PCEigenvectorsName, FilterParameter::Category::CreatedArray, PrincipalComponentAnalysis));
  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void PrincipalComponentAnalysis::initialize()
{
  m_SelectedWeakPtrVector.clear();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void PrincipalComponentAnalysis::dataCheck()
{
  clearErrorCode();
  clearWarningCode();
  initialize();

  if(getSelectedDataArrayPaths().size() < 2)
  {
    QString ss = QObject::tr("At least two Attribute Array must be selected");
    setErrorCondition(-11001, ss);
    return;
  }

  std::vector<DataArrayPath> paths = getSelectedDataArrayPaths();

  if(!DataArrayPath::ValidateVector(paths))
  {
    QString ss = QObject::tr("There are Attribute Arrays selected that are not contained in the same Attribute Matrix; all selected Attribute Arrays must belong to the same Attribute Matrix");
    setErrorCondition(-11002, ss);
  }

  for(auto&& path : paths)
  {
    IDataArray::WeakPointer ptr = getDataContainerArray()->getPrereqIDataArrayFromPath(this, path);
    if(ptr.lock())
    {
      m_SelectedWeakPtrVector.push_back(ptr);
      int32_t numComps = ptr.lock()->getNumberOfComponents();
      if(numComps != 1)
      {
        QString ss = QObject::tr("Attribute Arrays must be scalar arrays, but %1 has %2 total components").arg(ptr.lock()->getName()).arg(numComps);
        setErrorCondition(-11003, ss);
      }
    }
  }

  if(getErrorCode() < 0)
  {
    return;
  }

  std::vector<size_t> tDims(1, paths.size());

  DataContainer::Pointer m = getDataContainerArray()->getDataContainer(getSelectedDataArrayPaths().at(0).getDataContainerName());
  m->createNonPrereqAttributeMatrix(this, getPCAttributeMatrixName(), tDims, AttributeMatrix::Type::Generic, AttributeMatrixID21);

  DataArrayPath tempPath;
  std::vector<size_t> cDims(1, 1);

  tempPath.update(getSelectedDataArrayPaths().at(0).getDataContainerName(), getPCAttributeMatrixName(), getPCEigenvaluesName());
  m_PCEigenvaluesPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<double>>(this, tempPath, 0, cDims, "", DataArrayID31);
  if(m_PCEigenvaluesPtr.lock())
  {
    m_PCEigenvalues = m_PCEigenvaluesPtr.lock()->getPointer(0);
  }

  cDims[0] = paths.size();

  tempPath.update(getSelectedDataArrayPaths().at(0).getDataContainerName(), getPCAttributeMatrixName(), getPCEigenvectorsName());
  m_PCEigenvectorsPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<double>>(this, tempPath, 0, cDims, "", DataArrayID32);
  if(m_PCEigenvectorsPtr.lock())
  {
    m_PCEigenvectors = m_PCEigenvectorsPtr.lock()->getPointer(0);
  }

  if(getProjectDataSpace())
  {
    if(getNumberOfDimensionsForProjection() <= 0)
    {
      QString ss = QObject::tr("Number of dimensions for the projected space (%1) must be greater than 0").arg(getNumberOfDimensionsForProjection());
      setErrorCondition(-11004, ss);
    }

    if(getNumberOfDimensionsForProjection() > paths.size())
    {
      QString ss = QObject::tr("Number of dimensions for the projected space (%1) must be less than or equal to the number of selected Attribute Arrays (%2)")
                       .arg(getNumberOfDimensionsForProjection())
                       .arg(paths.size());
      setErrorCondition(-11005, ss);
    }

    cDims[0] = getNumberOfDimensionsForProjection();

    m_ProjectedDataSpacePtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<double>>(this, getProjectedDataSpaceArrayPath(), 0, cDims, "", DataArrayID33);
    if(m_ProjectedDataSpacePtr.lock())
    {
      m_ProjectedDataSpace = m_ProjectedDataSpacePtr.lock()->getPointer(0);
    }
    if(getErrorCode() >= 0)
    {
      paths.push_back(getProjectedDataSpaceArrayPath());
    }
  }

  getDataContainerArray()->validateNumberOfTuples(this, paths);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
template <typename T>
void copyDataArrays(IDataArray::Pointer dataPtr, std::vector<double>& copy, int32_t arrayIndex)
{
  typename DataArray<T>::Pointer inDataPtr = std::dynamic_pointer_cast<DataArray<T>>(dataPtr);
  T* dPtr = inDataPtr->getPointer(0);
  size_t numTuples = inDataPtr->getNumberOfTuples();

  for(size_t i = 0; i < numTuples; i++)
  {
    copy[numTuples * arrayIndex + i] = static_cast<double>(dPtr[i]);
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void PrincipalComponentAnalysis::execute()
{
  dataCheck();
  if(getErrorCode() < 0)
  {
    return;
  }

  if(m_SelectedDataArrayPaths.size() != m_SelectedWeakPtrVector.size())
  {
    QString ss = QObject::tr("The number of selected Attribute Arrays does not equal the number of internal weak pointers");
    setErrorCondition(-11008, ss);
    return;
  }

  // Set up our Eigen variables ; we use a dynamic double matrix since we don't know
  // the number of arrays or their tuples at compile time ; could have just used the
  // MatrixXd typedef instead of our own, but this is explicit about the ordering
  // being column major
  auto numArrays = m_SelectedWeakPtrVector.size();
  size_t numTuples = m_SelectedWeakPtrVector[0].lock()->getNumberOfTuples();

  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> CovarianceMatrix;

  // Copy the incoming data into a contiguous chunk of memory, casting everything to doubles
  std::vector<double> inDataVector(numTuples * numArrays, 0);

  for(auto i = 0; i < numArrays; i++)
  {
    EXECUTE_FUNCTION_TEMPLATE(this, copyDataArrays, m_SelectedWeakPtrVector[i].lock(), m_SelectedWeakPtrVector[i].lock(), inDataVector, i)
  }

  // Interface our CovarianceMatrix typedef with our contiguous array in memory
  // Using a colum major matrix, so each array is a column, so the arrays should
  // be layed out end-to-end in memory ; initialization for sizes in Eigen is:
  // numRows x numColumns
  Eigen::Map<CovarianceMatrix> dataMat(inDataVector.data(), numTuples, numArrays);

  // If the correlation approach is being used, the data must be standardize to have
  // mean 0 and unit variance
  if(m_MatrixApproach == 0)
  {
    // Standardize data to have mean 0 and unit variance
    dataMat = (dataMat.rowwise() - dataMat.colwise().mean());
    Eigen::RowVectorXd stdDev = dataMat.colwise().squaredNorm();
    for(auto i = 0; i < stdDev.cols(); i++)
    {
      stdDev(i) /= dataMat.rows();
      stdDev(i) = sqrt(stdDev(i));
    }
    dataMat = dataMat.array().rowwise() / stdDev.array();
  }

  // Calculate the covariance matrix using the matrix formulation, checking if
  // tuples are 1 to avoid division by zero
  // If correlation was chosen, then this is technically the correlation matrix,
  // but since the means are zero in that case the centering operation is trivial
  Eigen::MatrixXd centeredMat = dataMat.rowwise() - dataMat.colwise().mean();
  Eigen::MatrixXd covMat;

  if(numTuples == 1)
  {
    covMat = (centeredMat.adjoint() * centeredMat) / (1.0);
  }
  else
  {
    covMat = (centeredMat.adjoint() * centeredMat) / (double(numTuples - 1.0));
  }

  // Perform the eigen decomposition to get the eigenvectors and eigenvalues
  // of the covariance matrix, then copy those values into the primary
  // DataArray pointers
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> pca(covMat);

  for(auto i = 0; i < pca.eigenvalues().size(); i++)
  {
    m_PCEigenvalues[i] = pca.eigenvalues()(i);
  }

  for(auto i = 0; i < pca.eigenvectors().cols(); i++)
  {
    for(auto j = 0; j < pca.eigenvectors().rows(); j++)
    {
      m_PCEigenvectors[pca.eigenvectors().rows() * i + j] = pca.eigenvectors()(j, i);
    }
  }

  if(m_ProjectDataSpace)
  {
    // Extract the projective transform
    // Eigen orders the eigenvalues/eigenvectors in ascending order, so just grab
    // the rightmost columns equal to the number of projective dimensions
    Eigen::MatrixXd transform = pca.eigenvectors().rightCols(m_NumberOfDimensionsForProjection);

    // Multiply each row in the centered data matrix by transform and copy the projected
    // points into the DataArray pointer
    for(auto i = 0; i < centeredMat.rows(); i++)
    {
      Eigen::VectorXd tmp = centeredMat.row(i) * transform;
      for(int32_t j = 0; j < m_NumberOfDimensionsForProjection; j++)
      {
        m_ProjectedDataSpace[m_NumberOfDimensionsForProjection * i + j] = tmp(j);
      }
    }
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AbstractFilter::Pointer PrincipalComponentAnalysis::newFilterInstance(bool copyFilterParameters) const
{
  PrincipalComponentAnalysis::Pointer filter = PrincipalComponentAnalysis::New();
  if(copyFilterParameters)
  {
    copyFilterParameterInstanceVariables(filter.get());
  }
  return filter;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString PrincipalComponentAnalysis::getCompiledLibraryName() const
{
  return DREAM3DReviewConstants::DREAM3DReviewBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString PrincipalComponentAnalysis::getBrandingString() const
{
  return "DREAM3DReview";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString PrincipalComponentAnalysis::getFilterVersion() const
{
  QString version;
  QTextStream vStream(&version);
  vStream << DREAM3DReview::Version::Major() << "." << DREAM3DReview::Version::Minor() << "." << DREAM3DReview::Version::Patch();
  return version;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString PrincipalComponentAnalysis::getGroupName() const
{
  return DREAM3DReviewConstants::FilterGroups::DREAM3DReviewFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QUuid PrincipalComponentAnalysis::getUuid() const
{
  return QUuid("{ec163736-39c8-5c69-9a56-61940a337c07}");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString PrincipalComponentAnalysis::getSubGroupName() const
{
  return DREAM3DReviewConstants::FilterSubGroups::DimensionalityReductionFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString PrincipalComponentAnalysis::getHumanLabel() const
{
  return "Principal Component Analysis";
}

// -----------------------------------------------------------------------------
PrincipalComponentAnalysis::Pointer PrincipalComponentAnalysis::NullPointer()
{
  return Pointer(static_cast<Self*>(nullptr));
}

// -----------------------------------------------------------------------------
std::shared_ptr<PrincipalComponentAnalysis> PrincipalComponentAnalysis::New()
{
  struct make_shared_enabler : public PrincipalComponentAnalysis
  {
  };
  std::shared_ptr<make_shared_enabler> val = std::make_shared<make_shared_enabler>();
  val->setupFilterParameters();
  return val;
}

// -----------------------------------------------------------------------------
QString PrincipalComponentAnalysis::getNameOfClass() const
{
  return QString("PrincipalComponentAnalysis");
}

// -----------------------------------------------------------------------------
QString PrincipalComponentAnalysis::ClassName()
{
  return QString("PrincipalComponentAnalysis");
}

// -----------------------------------------------------------------------------
void PrincipalComponentAnalysis::setSelectedDataArrayPaths(const std::vector<DataArrayPath>& value)
{
  m_SelectedDataArrayPaths = value;
}

// -----------------------------------------------------------------------------
std::vector<DataArrayPath> PrincipalComponentAnalysis::getSelectedDataArrayPaths() const
{
  return m_SelectedDataArrayPaths;
}

// -----------------------------------------------------------------------------
void PrincipalComponentAnalysis::setPCAttributeMatrixName(const QString& value)
{
  m_PCAttributeMatrixName = value;
}

// -----------------------------------------------------------------------------
QString PrincipalComponentAnalysis::getPCAttributeMatrixName() const
{
  return m_PCAttributeMatrixName;
}

// -----------------------------------------------------------------------------
void PrincipalComponentAnalysis::setPCEigenvaluesName(const QString& value)
{
  m_PCEigenvaluesName = value;
}

// -----------------------------------------------------------------------------
QString PrincipalComponentAnalysis::getPCEigenvaluesName() const
{
  return m_PCEigenvaluesName;
}

// -----------------------------------------------------------------------------
void PrincipalComponentAnalysis::setPCEigenvectorsName(const QString& value)
{
  m_PCEigenvectorsName = value;
}

// -----------------------------------------------------------------------------
QString PrincipalComponentAnalysis::getPCEigenvectorsName() const
{
  return m_PCEigenvectorsName;
}

// -----------------------------------------------------------------------------
void PrincipalComponentAnalysis::setMatrixApproach(int value)
{
  m_MatrixApproach = value;
}

// -----------------------------------------------------------------------------
int PrincipalComponentAnalysis::getMatrixApproach() const
{
  return m_MatrixApproach;
}

// -----------------------------------------------------------------------------
void PrincipalComponentAnalysis::setProjectDataSpace(bool value)
{
  m_ProjectDataSpace = value;
}

// -----------------------------------------------------------------------------
bool PrincipalComponentAnalysis::getProjectDataSpace() const
{
  return m_ProjectDataSpace;
}

// -----------------------------------------------------------------------------
void PrincipalComponentAnalysis::setNumberOfDimensionsForProjection(int value)
{
  m_NumberOfDimensionsForProjection = value;
}

// -----------------------------------------------------------------------------
int PrincipalComponentAnalysis::getNumberOfDimensionsForProjection() const
{
  return m_NumberOfDimensionsForProjection;
}

// -----------------------------------------------------------------------------
void PrincipalComponentAnalysis::setProjectedDataSpaceArrayPath(const DataArrayPath& value)
{
  m_ProjectedDataSpaceArrayPath = value;
}

// -----------------------------------------------------------------------------
DataArrayPath PrincipalComponentAnalysis::getProjectedDataSpaceArrayPath() const
{
  return m_ProjectedDataSpaceArrayPath;
}
