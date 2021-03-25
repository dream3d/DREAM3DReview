/* ============================================================================
 * Copyright 2021 The University of Utah
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
 * Neither the name of BlueQuartz Software, the US Air Force, the University of Utah nor the names of its contributors may be used to endorse or promote products derived from this software without
 * specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY
 * OF SUCH DAMAGE.
 *
 * The code contained herein was partially funded by the following contracts:
 *
 *
 * This code contained herein is based upon work supported by the following grants:
 *    DOE Office of Nuclear Energy's Nuclear Energy University Program Grant No.: DE-NE0008799
 *    DOD Office of Economic Adjustment Grant No.: ST1605-19-03
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include "ComputeFeatureEigenstrains.h"

#include <cmath>

#include <QtCore/QTextStream>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigen>

#include "EbsdLib/Core/Orientation.hpp"
#include "EbsdLib/Core/OrientationTransformation.hpp"

#include "SIMPLib/Common/Constants.h"
#include "SIMPLib/DataContainers/DataContainer.h"
#include "SIMPLib/DataContainers/DataContainerArray.h"
#include "SIMPLib/FilterParameters/AbstractFilterParametersReader.h"
#include "SIMPLib/FilterParameters/DataArraySelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/FloatFilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedBooleanFilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedPathCreationFilterParameter.h"
#include "SIMPLib/FilterParameters/SeparatorFilterParameter.h"
#include "SIMPLib/Geometry/ImageGeom.h"
#include "SIMPLib/Math/SIMPLibMath.h"

#include "DREAM3DReview/DREAM3DReviewConstants.h"
#include "DREAM3DReview/DREAM3DReviewVersion.h"

#include "DREAM3DReview/DREAM3DReviewFilters/util/EigenstrainsHelper.hpp"

namespace SIMPLMath = SIMPLib::Constants;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
ComputeFeatureEigenstrains::ComputeFeatureEigenstrains() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
ComputeFeatureEigenstrains::~ComputeFeatureEigenstrains() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ComputeFeatureEigenstrains::initialize()
{
  clearErrorCode();
  clearWarningCode();
  setCancel(false);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ComputeFeatureEigenstrains::setupFilterParameters()
{
  FilterParameterVectorType parameters;

  // Poisson's ratio user input
  parameters.push_back(SIMPL_NEW_FLOAT_FP("Poisson's Ratio", PoissonRatio, FilterParameter::Category::Parameter, ComputeFeatureEigenstrains));

  std::vector<QString> linkedProps;
  linkedProps.push_back("AxisLengthsArrayPath");
  linkedProps.push_back("AxisEulerAnglesArrayPath");

  parameters.push_back(
      SIMPL_NEW_LINKED_BOOL_FP("Use Ellipsoidal Grains (versus spherical assumption)", UseEllipsoidalGrains, FilterParameter::Category::Parameter, ComputeFeatureEigenstrains, linkedProps));

  // Correctional matrix beta user inputs
  linkedProps.clear();
  linkedProps.push_back("Beta11");
  linkedProps.push_back("Beta22");
  linkedProps.push_back("Beta33");
  linkedProps.push_back("Beta23");
  linkedProps.push_back("Beta13");
  linkedProps.push_back("Beta12");
  parameters.push_back(SIMPL_NEW_LINKED_BOOL_FP("Use Correctional Matrix", UseCorrectionalMatrix, FilterParameter::Category::Parameter, ComputeFeatureEigenstrains, linkedProps));

  parameters.push_back(SIMPL_NEW_FLOAT_FP("Beta11", Beta11, FilterParameter::Category::Parameter, ComputeFeatureEigenstrains));
  parameters.push_back(SIMPL_NEW_FLOAT_FP("Beta22", Beta22, FilterParameter::Category::Parameter, ComputeFeatureEigenstrains));
  parameters.push_back(SIMPL_NEW_FLOAT_FP("Beta33", Beta33, FilterParameter::Category::Parameter, ComputeFeatureEigenstrains));
  parameters.push_back(SIMPL_NEW_FLOAT_FP("Beta23", Beta23, FilterParameter::Category::Parameter, ComputeFeatureEigenstrains));
  parameters.push_back(SIMPL_NEW_FLOAT_FP("Beta13", Beta13, FilterParameter::Category::Parameter, ComputeFeatureEigenstrains));
  parameters.push_back(SIMPL_NEW_FLOAT_FP("Beta12", Beta12, FilterParameter::Category::Parameter, ComputeFeatureEigenstrains));

  // Axis lengths and euler angles 3xN feature arrays
  parameters.push_back(SeparatorFilterParameter::Create("Cell Feature Data", FilterParameter::Category::RequiredArray));

  {
    DataArraySelectionFilterParameter::RequirementType req =
        DataArraySelectionFilterParameter::CreateRequirement(SIMPL::TypeNames::Float, 3, AttributeMatrix::Type::CellFeature, IGeometry::Type::Image);
    parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Axis Lengths", AxisLengthsArrayPath, FilterParameter::Category::RequiredArray, ComputeFeatureEigenstrains, req));
  }

  {
    DataArraySelectionFilterParameter::RequirementType req =
        DataArraySelectionFilterParameter::CreateRequirement(SIMPL::TypeNames::Float, 3, AttributeMatrix::Type::CellFeature, IGeometry::Type::Image);
    parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Axis Euler Angles", AxisEulerAnglesArrayPath, FilterParameter::Category::RequiredArray, ComputeFeatureEigenstrains, req));
  }

  // Elastic strain 6xN feature array
  {
    DataArraySelectionFilterParameter::RequirementType req =
        DataArraySelectionFilterParameter::CreateRequirement(SIMPL::TypeNames::Float, 6, AttributeMatrix::Type::CellFeature, IGeometry::Type::Image);
    parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Elastic Strains (Voigt Notation)", ElasticStrainsArrayPath, FilterParameter::Category::RequiredArray, ComputeFeatureEigenstrains, req));
  }

  // Output 6xN eigenstrain feature array
  parameters.push_back(SeparatorFilterParameter::Create("Cell Feature Data", FilterParameter::Category::CreatedArray));
  parameters.push_back(
      SIMPL_NEW_DA_WITH_LINKED_AM_FP("Eigenstrains", EigenstrainsArrayName, ElasticStrainsArrayPath, ElasticStrainsArrayPath, FilterParameter::Category::CreatedArray, ComputeFeatureEigenstrains));

  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ComputeFeatureEigenstrains::dataCheck()
{
  clearErrorCode();
  clearWarningCode();
  initialize();
  DataArrayPath tempPath;

  // Negative Poisson's ratio works in the calculations but warn
  if(getPoissonRatio() < 0.0f)
  {
    QString ss = QObject::tr("Poisson's ratio is negative");
    setWarningCondition(-94001, ss);
  }

  // Incompressible v=0.5 results in singular matrix
  if(getPoissonRatio() > 0.49999999f)
  {
    QString ss = QObject::tr("Poisson's ratio cannot be 0.5 or greater");
    setErrorCondition(-94002, ss);
  }

  // Check correction values and warn if they are large (potential typos)
  if(m_UseCorrectionalMatrix)
  {
    if(getBeta11() < 0.5f || getBeta11() > 2.0f)
    {
      QString ss = QObject::tr("Beta11 correction is pretty large; this may be a typo");
      setWarningCondition(-94003, ss);
    }

    if(getBeta22() < 0.5f || getBeta22() > 2.0f)
    {
      QString ss = QObject::tr("Beta22 correction is pretty large; this may be a typo");
      setWarningCondition(-94004, ss);
    }

    if(getBeta33() < 0.5f || getBeta33() > 2.0f)
    {
      QString ss = QObject::tr("Beta33 correction is pretty large; this may be a typo");
      setWarningCondition(-94005, ss);
    }

    if(getBeta23() < 0.5f || getBeta23() > 2.0f)
    {
      QString ss = QObject::tr("Beta23 correction is pretty large; this may be a typo");
      setWarningCondition(-94006, ss);
    }

    if(getBeta13() < 0.5f || getBeta13() > 2.0f)
    {
      QString ss = QObject::tr("Beta13 correction is pretty large; this may be a typo");
      setWarningCondition(-94007, ss);
    }

    if(getBeta12() < 0.5f || getBeta12() > 2.0f)
    {
      QString ss = QObject::tr("Beta12 correction is pretty large; this may be a typo");
      setWarningCondition(-94008, ss);
    }
  }

  // Check Required Objects
  std::vector<size_t> cDims(1, 1);
  if(m_UseEllipsoidalGrains)
  {
    cDims[0] = 3;
    m_AxisLengthsPtr = getDataContainerArray()->getPrereqArrayFromPath<FloatArrayType>(this, getAxisLengthsArrayPath(), cDims);
    if(getErrorCode() < 0)
    {
      return;
    }

    cDims[0] = 3;
    m_AxisEulerAnglesPtr = getDataContainerArray()->getPrereqArrayFromPath<FloatArrayType>(this, getAxisEulerAnglesArrayPath(), cDims);
    if(getErrorCode() < 0)
    {
      return;
    }
  }

  cDims[0] = 6;
  m_ElasticStrainsPtr = getDataContainerArray()->getPrereqArrayFromPath<FloatArrayType>(this, getElasticStrainsArrayPath(), cDims);
  if(getErrorCode() < 0)
  {
    return;
  }

  // Check Eigenstrain output
  cDims[0] = 6;
  tempPath.update(getElasticStrainsArrayPath().getDataContainerName(), getElasticStrainsArrayPath().getAttributeMatrixName(), getEigenstrainsArrayName());
  m_EigenstrainsPtr = getDataContainerArray()->createNonPrereqArrayFromPath<FloatArrayType>(this, tempPath, 0, cDims, "", 1);
  if(getErrorCode() < 0)
  {
    return;
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ComputeFeatureEigenstrains::execute()
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

  find_eigenstrains();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ComputeFeatureEigenstrains::find_eigenstrains()
{
  size_t numfeatures = m_ElasticStrainsPtr.lock()->getNumberOfTuples();

  double phi1 = 0.0;
  double theta = 0.0;
  double phi2 = 0.0;

  double semiAxisA = 0.0;
  double semiAxisB = 0.0;
  double semiAxisC = 0.0;

  double E11 = 0.0;
  double E22 = 0.0;
  double E33 = 0.0;
  double E23 = 0.0;
  double E13 = 0.0;
  double E12 = 0.0;

  Eigen::Matrix<double, 3, 3> beta;
  beta.setOnes(3, 3); // Default no correction
  if(m_UseCorrectionalMatrix)
  {
    // clang-format off
    beta << m_Beta11, m_Beta12, m_Beta13,
            m_Beta12, m_Beta22, m_Beta23,
            m_Beta13, m_Beta23, m_Beta33;
    // clang-format on
  }

  OrientationD orientationMatrix;
  Eigen::Matrix<double, 3, 3> OM;
  Eigen::Matrix<double, 3, 3> OMT;
  OM.setIdentity(3, 3); // Defaults to no orientation identity
  OMT.setIdentity(3, 3);

  Eigen::Matrix<double, 3, 3> elasticStrainTensor;
  Eigen::Matrix<double, 3, 3> elasticStrainTensorRot;

  EigenstrainsHelper::Tensor4DType eshelbyTensor;
  EigenstrainsHelper::Tensor4DType eshelbyInverse;

  Eigen::Matrix<double, 9, 9> eshelbyTensor99;
  Eigen::Matrix<double, 9, 9> eshelbyInverse99;
  Eigen::Matrix<double, 9, 9> I9;
  I9.setIdentity(9, 9);

  Eigen::Matrix<double, 3, 3> eigenstrainTensorRot;
  Eigen::Matrix<double, 3, 3> eigenstrainTensor;
  Eigen::Matrix<double, 3, 3> eigenstrainTensorCorrected;

  FloatArrayType& axisEulerAngles = *(m_AxisEulerAnglesPtr.lock().get());
  FloatArrayType& axisLengths = *(m_AxisLengthsPtr.lock().get());
  FloatArrayType& elasticStrains = *(m_ElasticStrainsPtr.lock().get());
  FloatArrayType& eigenstrains = *(m_EigenstrainsPtr.lock().get());

  // small eps term added to euler angle bounds check as FindFeatureShapes has some small noise outside the bounds
  double eps = 0.001;
  double eulerPhiMax = 2 * SIMPLMath::k_PiD + eps;
  double eulerThetaMax = SIMPLMath::k_PiD + eps;
  double eulerMin = 0 - eps;

  for(size_t feature = 0; feature < numfeatures; feature++)
  {
    if(m_UseEllipsoidalGrains)
    {
      phi1 = axisEulerAngles[feature * 3 + 0];
      theta = axisEulerAngles[feature * 3 + 1];
      phi2 = axisEulerAngles[feature * 3 + 2];

      if(std::isnan(phi1) || std::isnan(theta) || std::isnan(phi2))
      {
        QString ss = QObject::tr("NaN Axis Euler angle found in feature ID #%1, skipping").arg(feature);
        notifyStatusMessage(ss);
        eigenstrains[feature * 6 + 0] = 0;
        eigenstrains[feature * 6 + 1] = 0;
        eigenstrains[feature * 6 + 2] = 0;
        eigenstrains[feature * 6 + 3] = 0;
        eigenstrains[feature * 6 + 4] = 0;
        eigenstrains[feature * 6 + 5] = 0;
        continue;
      }

      if(phi1 > eulerPhiMax || phi1 < eulerMin)
      {
        QString ss = QObject::tr("Feature %1 euler angle phi1=%2 out of bounds 2pi >= phi1 >= 0. Euler angles may be in degrees").arg(feature).arg(phi1);
        setErrorCondition(-94000, ss);
        return;
      }

      if(theta > eulerThetaMax || theta < eulerMin)
      {
        QString ss = QObject::tr("Feature %1 euler angle theta=%2 out of bounds pi >= theta >= 0. Euler angles may be in degrees").arg(feature).arg(theta);
        setErrorCondition(-94000, ss);
        return;
      }

      if(phi2 > eulerPhiMax || phi2 < eulerMin)
      {
        QString ss = QObject::tr("Feature %1 euler angle phi2=%2 out of bounds 2pi >= phi2 >= 0. Euler angles may be in degrees").arg(feature).arg(phi2);
        setErrorCondition(-94000, ss);
        return;
      }

      semiAxisA = axisLengths[feature * 3 + 0];
      semiAxisB = axisLengths[feature * 3 + 1];
      semiAxisC = axisLengths[feature * 3 + 2];

      if(std::isnan(semiAxisA) || std::isnan(semiAxisB) || std::isnan(semiAxisC))
      {
        QString ss = QObject::tr("NaN Axis length found in feature ID #%1, skipping").arg(feature);
        notifyStatusMessage(ss);
        eigenstrains[feature * 6 + 0] = 0;
        eigenstrains[feature * 6 + 1] = 0;
        eigenstrains[feature * 6 + 2] = 0;
        eigenstrains[feature * 6 + 3] = 0;
        eigenstrains[feature * 6 + 4] = 0;
        eigenstrains[feature * 6 + 5] = 0;
        continue;
      }

      if(semiAxisB > semiAxisA)
      {
        QString ss = QObject::tr("Feature %1 semi-axis b=%2 is greater than semi-axis a=%3. Criteria a>=b>=c must be satisfied").arg(feature).arg(semiAxisB).arg(semiAxisA);
        setErrorCondition(-94000, ss);
        return;
      }

      if(semiAxisC > semiAxisB)
      {
        QString ss = QObject::tr("Feature %1 semi-axis c=%2 is greater than semi-axis b=%3. Criteria a>=b>=c must be satisfied").arg(feature).arg(semiAxisC).arg(semiAxisB);
        setErrorCondition(-94000, ss);
        return;
      }
    }

    E11 = elasticStrains[feature * 6 + 0];
    E22 = elasticStrains[feature * 6 + 1];
    E33 = elasticStrains[feature * 6 + 2];
    E23 = elasticStrains[feature * 6 + 3];
    E13 = elasticStrains[feature * 6 + 4];
    E12 = elasticStrains[feature * 6 + 5];

    // clang-format off
    elasticStrainTensor << E11, E12, E13,
                           E12, E22, E23,
                           E13, E23, E33;
    // clang-format on

    // Check if the elastic strains are zero (or negligible) then so are the eigenstrains
    if(elasticStrainTensor.isMuchSmallerThan(1e-10))
    {
      eigenstrains[feature * 6 + 0] = 0;
      eigenstrains[feature * 6 + 1] = 0;
      eigenstrains[feature * 6 + 2] = 0;
      eigenstrains[feature * 6 + 3] = 0;
      eigenstrains[feature * 6 + 4] = 0;
      eigenstrains[feature * 6 + 5] = 0;
      continue;
    }

    if(m_UseEllipsoidalGrains)
    {
      orientationMatrix = OrientationTransformation::eu2om<OrientationD, OrientationD>({phi1, theta, phi2});

      // clang-format off
      OM << orientationMatrix[0], orientationMatrix[1], orientationMatrix[2],
            orientationMatrix[3], orientationMatrix[4], orientationMatrix[5],
            orientationMatrix[6], orientationMatrix[7], orientationMatrix[8];
      OMT = OM.transpose();
      // clang-format on
    }

    // Change basis into ellipsoid reference frame | e' = Q e Q^T
    elasticStrainTensorRot = OM * elasticStrainTensor * OMT;

    // Calculate Eshelby tensor | S = f(a, b, c, nu)
    eshelbyTensor = EigenstrainsHelper::find_eshelby(semiAxisA, semiAxisB, semiAxisC, m_PoissonRatio, m_UseEllipsoidalGrains);

    // Map Eshelby tensor into 9x9 matrix
    size_t col = 0;
    size_t row = 0;
    for(size_t i = 0; i < 3; i++)
    {
      for(size_t j = 0; j < 3; j++)
      {
        for(size_t k = 0; k < 3; k++)
        {
          for(size_t l = 0; l < 3; l++)
          {
            eshelbyTensor99(col, row) = eshelbyTensor(i, j, k, l);
            col++;
          }
        }
        row++;
        col = 0;
      }
    }

    // Calculate inverse | (S-I)^-1
    eshelbyInverse99 = (eshelbyTensor99 - I9).inverse();

    // Remap inverse back into a 3x3x3x3 tensor
    col = 0;
    row = 0;
    for(size_t i = 0; i < 3; i++)
    {
      for(size_t j = 0; j < 3; j++)
      {
        for(size_t k = 0; k < 3; k++)
        {
          for(size_t l = 0; l < 3; l++)
          {
            eshelbyInverse(i, j, k, l) = eshelbyInverse99(col, row);
            col++;
          }
        }
        row++;
        col = 0;
      }
    }

    // Calculate eigenstrain tensor | e*' = (S-I)^-1 e'
    eigenstrainTensorRot.setZero(3, 3);
    for(size_t i = 0; i < 3; i++)
    {
      for(size_t j = 0; j < 3; j++)
      {
        for(size_t k = 0; k < 3; k++)
        {
          for(size_t l = 0; l < 3; l++)
          {
            eigenstrainTensorRot(i, j) += eshelbyInverse(i, j, k, l) * elasticStrainTensorRot(k, l);
          }
        }
      }
    }

    // Change basis back to global reference frame | e* = Q^T e*' Q
    eigenstrainTensor = OMT * eigenstrainTensorRot * OM;

    // Add correction | e*c = B e*
    eigenstrainTensorCorrected = eigenstrainTensor.cwiseProduct(beta);

    eigenstrains[feature * 6 + 0] = eigenstrainTensorCorrected(0, 0);
    eigenstrains[feature * 6 + 1] = eigenstrainTensorCorrected(1, 1);
    eigenstrains[feature * 6 + 2] = eigenstrainTensorCorrected(2, 2);
    eigenstrains[feature * 6 + 3] = eigenstrainTensorCorrected(1, 2);
    eigenstrains[feature * 6 + 4] = eigenstrainTensorCorrected(0, 2);
    eigenstrains[feature * 6 + 5] = eigenstrainTensorCorrected(0, 1);
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AbstractFilter::Pointer ComputeFeatureEigenstrains::newFilterInstance(bool copyFilterParameters) const
{
  ComputeFeatureEigenstrains::Pointer filter = ComputeFeatureEigenstrains::New();
  if(copyFilterParameters)
  {
    copyFilterParameterInstanceVariables(filter.get());
  }
  return filter;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ComputeFeatureEigenstrains::getCompiledLibraryName() const
{
  return DREAM3DReviewConstants::DREAM3DReviewBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ComputeFeatureEigenstrains::getBrandingString() const
{
  return "DREAM3DReview";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ComputeFeatureEigenstrains::getFilterVersion() const
{
  QString version;
  QTextStream vStream(&version);
  vStream << DREAM3DReview::Version::Major() << "." << DREAM3DReview::Version::Minor() << "." << DREAM3DReview::Version::Patch();
  return version;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ComputeFeatureEigenstrains::getGroupName() const
{
  return DREAM3DReviewConstants::FilterGroups::DREAM3DReviewFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ComputeFeatureEigenstrains::getSubGroupName() const
{
  return DREAM3DReviewConstants::FilterSubGroups::RegistrationFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ComputeFeatureEigenstrains::getHumanLabel() const
{
  return "Compute Eigenstrains by Feature (Grain/Inclusion)";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QUuid ComputeFeatureEigenstrains::getUuid() const
{
  return QUuid("{879e1eb8-40dc-5a5b-abe5-7e0baa77ed73}");
}

// -----------------------------------------------------------------------------
ComputeFeatureEigenstrains::Pointer ComputeFeatureEigenstrains::NullPointer()
{
  return Pointer(static_cast<Self*>(nullptr));
}

// -----------------------------------------------------------------------------
std::shared_ptr<ComputeFeatureEigenstrains> ComputeFeatureEigenstrains::New()
{
  struct make_shared_enabler : public ComputeFeatureEigenstrains
  {
  };
  std::shared_ptr<make_shared_enabler> val = std::make_shared<make_shared_enabler>();
  val->setupFilterParameters();
  return val;
}

// -----------------------------------------------------------------------------
QString ComputeFeatureEigenstrains::getNameOfClass() const
{
  return QString("ComputeFeatureEigenstrains");
}

// -----------------------------------------------------------------------------
QString ComputeFeatureEigenstrains::ClassName()
{
  return QString("ComputeFeatureEigenstrains");
}

// -----------------------------------------------------------------------------
void ComputeFeatureEigenstrains::setPoissonRatio(float value)
{
  m_PoissonRatio = value;
}

// -----------------------------------------------------------------------------
float ComputeFeatureEigenstrains::getPoissonRatio() const
{
  return m_PoissonRatio;
}

// -----------------------------------------------------------------------------
void ComputeFeatureEigenstrains::setUseEllipsoidalGrains(bool value)
{
  m_UseEllipsoidalGrains = value;
}

// -----------------------------------------------------------------------------
bool ComputeFeatureEigenstrains::getUseEllipsoidalGrains() const
{
  return m_UseEllipsoidalGrains;
}

// -----------------------------------------------------------------------------
void ComputeFeatureEigenstrains::setUseCorrectionalMatrix(bool value)
{
  m_UseCorrectionalMatrix = value;
}

// -----------------------------------------------------------------------------
bool ComputeFeatureEigenstrains::getUseCorrectionalMatrix() const
{
  return m_UseCorrectionalMatrix;
}

// -----------------------------------------------------------------------------
void ComputeFeatureEigenstrains::setBeta11(float value)
{
  m_Beta11 = value;
}

// -----------------------------------------------------------------------------
float ComputeFeatureEigenstrains::getBeta11() const
{
  return m_Beta11;
}

// -----------------------------------------------------------------------------
void ComputeFeatureEigenstrains::setBeta22(float value)
{
  m_Beta22 = value;
}

// -----------------------------------------------------------------------------
float ComputeFeatureEigenstrains::getBeta22() const
{
  return m_Beta22;
}

// -----------------------------------------------------------------------------
void ComputeFeatureEigenstrains::setBeta33(float value)
{
  m_Beta33 = value;
}

// -----------------------------------------------------------------------------
float ComputeFeatureEigenstrains::getBeta33() const
{
  return m_Beta33;
}

// -----------------------------------------------------------------------------
void ComputeFeatureEigenstrains::setBeta23(float value)
{
  m_Beta23 = value;
}

// -----------------------------------------------------------------------------
float ComputeFeatureEigenstrains::getBeta23() const
{
  return m_Beta23;
}

// -----------------------------------------------------------------------------
void ComputeFeatureEigenstrains::setBeta13(float value)
{
  m_Beta13 = value;
}

// -----------------------------------------------------------------------------
float ComputeFeatureEigenstrains::getBeta13() const
{
  return m_Beta13;
}

// -----------------------------------------------------------------------------
void ComputeFeatureEigenstrains::setBeta12(float value)
{
  m_Beta12 = value;
}

// -----------------------------------------------------------------------------
float ComputeFeatureEigenstrains::getBeta12() const
{
  return m_Beta12;
}

// -----------------------------------------------------------------------------
void ComputeFeatureEigenstrains::setAxisLengthsArrayPath(const DataArrayPath& value)
{
  m_AxisLengthsArrayPath = value;
}

// -----------------------------------------------------------------------------
DataArrayPath ComputeFeatureEigenstrains::getAxisLengthsArrayPath() const
{
  return m_AxisLengthsArrayPath;
}

// -----------------------------------------------------------------------------
void ComputeFeatureEigenstrains::setAxisEulerAnglesArrayPath(const DataArrayPath& value)
{
  m_AxisEulerAnglesArrayPath = value;
}

// -----------------------------------------------------------------------------
DataArrayPath ComputeFeatureEigenstrains::getAxisEulerAnglesArrayPath() const
{
  return m_AxisEulerAnglesArrayPath;
}

// -----------------------------------------------------------------------------
void ComputeFeatureEigenstrains::setElasticStrainsArrayPath(const DataArrayPath& value)
{
  m_ElasticStrainsArrayPath = value;
}

// -----------------------------------------------------------------------------
DataArrayPath ComputeFeatureEigenstrains::getElasticStrainsArrayPath() const
{
  return m_ElasticStrainsArrayPath;
}

// -----------------------------------------------------------------------------
void ComputeFeatureEigenstrains::setEigenstrainsArrayName(const QString& value)
{
  m_EigenstrainsArrayName = value;
}

// -----------------------------------------------------------------------------
QString ComputeFeatureEigenstrains::getEigenstrainsArrayName() const
{
  return m_EigenstrainsArrayName;
}
