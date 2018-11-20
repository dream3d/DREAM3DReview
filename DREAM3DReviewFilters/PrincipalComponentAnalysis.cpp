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

#include "SIMPLib/Common/Constants.h"
#include "SIMPLib/Common/TemplateHelpers.h"
#include "SIMPLib/FilterParameters/AbstractFilterParametersReader.h"
#include "SIMPLib/FilterParameters/BooleanFilterParameter.h"
#include "SIMPLib/FilterParameters/ChoiceFilterParameter.h"
#include "SIMPLib/FilterParameters/DataArrayCreationFilterParameter.h"
#include "SIMPLib/FilterParameters/IntFilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedBooleanFilterParameter.h"
#include "SIMPLib/FilterParameters/MultiDataArraySelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/SeparatorFilterParameter.h"
#include "SIMPLib/FilterParameters/StringFilterParameter.h"

#include "DREAM3DReview/DREAM3DReviewConstants.h"
#include "DREAM3DReview/DREAM3DReviewVersion.h"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
PrincipalComponentAnalysis::PrincipalComponentAnalysis()
: m_SelectedDataArrayPaths(QVector<DataArrayPath>())
, m_PCAttributeMatrixName("PrincipalComponentAnalysis")
, m_PCEigenvaluesName("PrincipalComponentEigenvalues")
, m_PCEigenvectorsName("PrincipalComponentEigenvectors")
, m_MatrixApproach(0)
, m_ProjectDataSpace(false)
, m_NumberOfDimensionsForProjection(0)
, m_ProjectedDataSpaceArrayPath("", "", "ProjectedDataSpace")
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
PrincipalComponentAnalysis::~PrincipalComponentAnalysis()
= default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void PrincipalComponentAnalysis::setupFilterParameters()
{
  FilterParameterVector parameters;
  MultiDataArraySelectionFilterParameter::RequirementType mdaReq =
      MultiDataArraySelectionFilterParameter::CreateRequirement(SIMPL::Defaults::AnyPrimitive, 1, AttributeMatrix::Type::Any, IGeometry::Type::Any);
  parameters.push_back(SIMPL_NEW_MDA_SELECTION_FP("Attribute Arrays for Computing Principal Components", SelectedDataArrayPaths, FilterParameter::RequiredArray, PrincipalComponentAnalysis, mdaReq));
  {
    ChoiceFilterParameter::Pointer choices = ChoiceFilterParameter::New();
    choices->setHumanLabel("Matrix Approach");
    choices->setPropertyName("MatrixApproach");
    choices->setSetterCallback(SIMPL_BIND_SETTER(PrincipalComponentAnalysis, this, MatrixApproach));
    choices->setGetterCallback(SIMPL_BIND_GETTER(PrincipalComponentAnalysis, this, MatrixApproach));
    QVector<QString> approaches = {"Correlation", "Covariance"};
    choices->setChoices(approaches);
    choices->setCategory(FilterParameter::Parameter);
    parameters.push_back(choices);
  }
  QStringList linkedProps = {"NumberOfDimensionsForProjection", "ProjectedDataSpaceArrayPath"};
  parameters.push_back(SIMPL_NEW_LINKED_BOOL_FP("Project Data Space", ProjectDataSpace, FilterParameter::Parameter, PrincipalComponentAnalysis, linkedProps));
  parameters.push_back(SIMPL_NEW_INTEGER_FP("Number of Dimensions for Projection", NumberOfDimensionsForProjection, FilterParameter::Parameter, PrincipalComponentAnalysis));
  DataArrayCreationFilterParameter::RequirementType dacReq = DataArrayCreationFilterParameter::CreateRequirement(AttributeMatrix::Type::Any, IGeometry::Type::Any);
  parameters.push_back(SIMPL_NEW_DA_CREATION_FP("Projected Data Space", ProjectedDataSpaceArrayPath, FilterParameter::CreatedArray, PrincipalComponentAnalysis, dacReq));
  parameters.push_back(SeparatorFilterParameter::New("Principal Component Data", FilterParameter::CreatedArray));
  parameters.push_back(SIMPL_NEW_STRING_FP("Principal Component Attribute Matrix", PCAttributeMatrixName, FilterParameter::CreatedArray, PrincipalComponentAnalysis));
  parameters.push_back(SIMPL_NEW_STRING_FP("Principal Component Eigenvalues", PCEigenvaluesName, FilterParameter::CreatedArray, PrincipalComponentAnalysis));
  parameters.push_back(SIMPL_NEW_STRING_FP("Principal Component Eigenvectors", PCEigenvectorsName, FilterParameter::CreatedArray, PrincipalComponentAnalysis));
  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void PrincipalComponentAnalysis::readFilterParameters(AbstractFilterParametersReader* reader, int index)
{
  reader->openFilterGroup(this, index);
  setSelectedDataArrayPaths(reader->readDataArrayPathVector("SelectedDataArrayPaths", getSelectedDataArrayPaths()));
  setPCAttributeMatrixName(reader->readString("PCAttributeMatrixName", getPCAttributeMatrixName()));
  setPCEigenvaluesName(reader->readString("PCEigenvaluesName", getPCEigenvaluesName()));
  setPCEigenvectorsName(reader->readString("PCEigenvectorsName", getPCEigenvectorsName()));
  setMatrixApproach(reader->readValue("MatrixApproach", getMatrixApproach()));
  setProjectDataSpace(reader->readValue("ProjectDataSpace", getProjectDataSpace()));
  setNumberOfDimensionsForProjection(reader->readValue("NumberOfDimensionsForProjection", getNumberOfDimensionsForProjection()));
  setProjectedDataSpaceArrayPath(reader->readDataArrayPath("ProjectedDataSpaceArrayPath", getProjectedDataSpaceArrayPath()));
  reader->closeFilterGroup();
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
  setErrorCondition(0);
  setWarningCondition(0);
  initialize();

  if(getSelectedDataArrayPaths().size() < 2)
  {
    setErrorCondition(-11001);
    QString ss = QObject::tr("At least two Attribute Array must be selected");
    notifyErrorMessage(getHumanLabel(), ss, getErrorCondition());
    return;
  }

  QVector<DataArrayPath> paths = getSelectedDataArrayPaths();

  if(!DataArrayPath::ValidateVector(paths))
  {
    setErrorCondition(-11002);
    QString ss = QObject::tr("There are Attribute Arrays selected that are not contained in the same Attribute Matrix; all selected Attribute Arrays must belong to the same Attribute Matrix");
    notifyErrorMessage(getHumanLabel(), ss, getErrorCondition());
  }

  for(auto&& path : paths)
  {
    IDataArray::WeakPointer ptr = getDataContainerArray()->getPrereqIDataArrayFromPath<IDataArray, AbstractFilter>(this, path);
    if(ptr.lock())
    {
      m_SelectedWeakPtrVector.push_back(ptr);
      int32_t numComps = ptr.lock()->getNumberOfComponents();
      if(numComps != 1)
      {
        setErrorCondition(-11003);
        QString ss = QObject::tr("Attribute Arrays must be scalar arrays, but %1 has %2 total components").arg(ptr.lock()->getName()).arg(numComps);
        notifyErrorMessage(getHumanLabel(), ss, getErrorCondition());
      }
    }
  }

  if(getErrorCondition() < 0)
  {
    return;
  }

  QVector<size_t> tDims(1, paths.size());

  DataContainer::Pointer m = getDataContainerArray()->getDataContainer(getSelectedDataArrayPaths()[0].getDataContainerName());
  m->createNonPrereqAttributeMatrix(this, getPCAttributeMatrixName(), tDims, AttributeMatrix::Type::Generic);

  DataArrayPath tempPath;
  QVector<size_t> cDims(1, 1);

  tempPath.update(getSelectedDataArrayPaths()[0].getDataContainerName(), getPCAttributeMatrixName(), getPCEigenvaluesName());
  m_PCEigenvaluesPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<double>, AbstractFilter, double>(this, tempPath, 0, cDims);
  if(m_PCEigenvaluesPtr.lock())
  {
    m_PCEigenvalues = m_PCEigenvaluesPtr.lock()->getPointer(0);
  }

  cDims[0] = paths.size();

  tempPath.update(getSelectedDataArrayPaths()[0].getDataContainerName(), getPCAttributeMatrixName(), getPCEigenvectorsName());
  m_PCEigenvectorsPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<double>, AbstractFilter, double>(this, tempPath, 0, cDims);
  if(m_PCEigenvectorsPtr.lock())
  {
    m_PCEigenvectors = m_PCEigenvectorsPtr.lock()->getPointer(0);
  }

  if(getProjectDataSpace())
  {
    if(getNumberOfDimensionsForProjection() <= 0)
    {
      setErrorCondition(-11004);
      QString ss = QObject::tr("Number of dimensions for the projected space (%1) must be greater than 0").arg(getNumberOfDimensionsForProjection());
      notifyErrorMessage(getHumanLabel(), ss, getErrorCondition());
    }

    if(getNumberOfDimensionsForProjection() > paths.size())
    {
      setErrorCondition(-11005);
      QString ss = QObject::tr("Number of dimensions for the projected space (%1) must be less than or equal to the number of selected Attribute Arrays (%2)")
                       .arg(getNumberOfDimensionsForProjection())
                       .arg(paths.size());
      notifyErrorMessage(getHumanLabel(), ss, getErrorCondition());
    }

    cDims[0] = getNumberOfDimensionsForProjection();

    m_ProjectedDataSpacePtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<double>, AbstractFilter, double>(this, getProjectedDataSpaceArrayPath(), 0, cDims);
    if(m_ProjectedDataSpacePtr.lock())
    {
      m_ProjectedDataSpace = m_ProjectedDataSpacePtr.lock()->getPointer(0);
    }
    if(getErrorCondition() >= 0)
    {
      paths.push_back(getProjectedDataSpaceArrayPath());
    }
  }

  getDataContainerArray()->validateNumberOfTuples<AbstractFilter>(this, paths);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void PrincipalComponentAnalysis::preflight()
{
  setInPreflight(true);
  emit preflightAboutToExecute();
  emit updateFilterParameters(this);
  dataCheck();
  emit preflightExecuted();
  setInPreflight(false);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
template <typename T> void copyDataArrays(IDataArray::Pointer dataPtr, std::vector<double>& copy, int32_t arrayIndex)
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
  setErrorCondition(0);
  setWarningCondition(0);
  dataCheck();
  if(getErrorCondition() < 0)
  {
    return;
  }

  if(m_SelectedDataArrayPaths.size() != m_SelectedWeakPtrVector.size())
  {
    QString ss = QObject::tr("The number of selected Attribute Arrays does not equal the number of internal weak pointers");
    setErrorCondition(-11008);
    notifyErrorMessage(getHumanLabel(), ss, getErrorCondition());
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
    EXECUTE_FUNCTION_TEMPLATE(this, copyDataArrays, m_SelectedWeakPtrVector[i].lock(), m_SelectedWeakPtrVector[i].lock(), inDataVector, i);
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

  notifyStatusMessage(getHumanLabel(), "Complete");
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
const QString PrincipalComponentAnalysis::getCompiledLibraryName() const
{
  return DREAM3DReviewConstants::DREAM3DReviewBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString PrincipalComponentAnalysis::getBrandingString() const
{
  return "DREAM3DReview";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString PrincipalComponentAnalysis::getFilterVersion() const
{
  QString version;
  QTextStream vStream(&version);
  vStream << DREAM3DReview::Version::Major() << "." << DREAM3DReview::Version::Minor() << "." << DREAM3DReview::Version::Patch();
  return version;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString PrincipalComponentAnalysis::getGroupName() const
{
  return DREAM3DReviewConstants::FilterGroups::DREAM3DReviewFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QUuid PrincipalComponentAnalysis::getUuid()
{
  return QUuid("{ec163736-39c8-5c69-9a56-61940a337c07}");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString PrincipalComponentAnalysis::getSubGroupName() const
{
  return DREAM3DReviewConstants::FilterSubGroups::DimensionalityReductionFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString PrincipalComponentAnalysis::getHumanLabel() const
{
  return "Principal Component Analysis";
}
