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

#include "ApplyTransformationToGeometry.h"

#include <Eigen/Dense>

#include "SIMPLib/Common/Constants.h"
#include "SIMPLib/FilterParameters/AbstractFilterParametersReader.h"
#include "SIMPLib/FilterParameters/DataArraySelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/DataContainerSelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/DynamicTableFilterParameter.h"
#include "SIMPLib/FilterParameters/FloatFilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedChoicesFilterParameter.h"
#include "SIMPLib/Geometry/EdgeGeom.h"
#include "SIMPLib/Geometry/IGeometry2D.h"
#include "SIMPLib/Geometry/IGeometry3D.h"
#include "SIMPLib/Geometry/VertexGeom.h"

#include "OrientationLib/OrientationMath/OrientationTransforms.hpp"

#include "DREAM3DReview/DREAM3DReviewConstants.h"
#include "DREAM3DReview/DREAM3DReviewVersion.h"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
ApplyTransformationToGeometry::ApplyTransformationToGeometry()
: m_ComputedTransformationMatrix("", "", "TransformationMatrix")
, m_GeometryToTransform("")
, m_TransformationMatrixType(1)
{
  m_RotationAngle = 0.0f;
  m_RotationAxis.x = 0.0f;
  m_RotationAxis.y = 0.0f;
  m_RotationAxis.z = 1.0f;

  m_Translation.x = 0.0f;
  m_Translation.y = 0.0f;
  m_Translation.z = 0.0f;

  m_Scale.x = 0.0f;
  m_Scale.y = 0.0f;
  m_Scale.z = 0.0f;

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
ApplyTransformationToGeometry::~ApplyTransformationToGeometry()
= default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ApplyTransformationToGeometry::setupFilterParameters()
{
  FilterParameterVector parameters;
  {
    LinkedChoicesFilterParameter::Pointer parameter = LinkedChoicesFilterParameter::New();
    parameter->setHumanLabel("Transformation Type");
    parameter->setPropertyName("TransformationMatrixType");
    parameter->setSetterCallback(SIMPL_BIND_SETTER(ApplyTransformationToGeometry, this, TransformationMatrixType));
    parameter->setGetterCallback(SIMPL_BIND_GETTER(ApplyTransformationToGeometry, this, TransformationMatrixType));
    QVector<QString> choices;
    choices.push_back("No Transformation");
    choices.push_back("Pre-Computed Transformation Matrix");
    choices.push_back("Manual Transformation Matrix");
    choices.push_back("Rotation");
    choices.push_back("Translation");
    choices.push_back("Scale");
    parameter->setChoices(choices);
    QStringList linkedProps = {"ComputedTransformationMatrix", "ManualTransformationMatrix", "RotationAngle", "RotationAxis", "Translation", "Scale"};
    parameter->setLinkedProperties(linkedProps);
    parameter->setEditable(false);
    parameter->setCategory(FilterParameter::Parameter);
    parameters.push_back(parameter);
  }
  {
    QStringList rHeaders, cHeaders;
    std::vector<std::vector<double>> defaultTable;
    for(size_t i = 0; i < 4; i++)
    {
      std::vector<double> row(4, 0);
      row[i] = 1.0;
      defaultTable.push_back(row);
    }
    m_ManualTransformationMatrix.setTableData(defaultTable);
    parameters.push_back(SIMPL_NEW_DYN_TABLE_FP("Transformation Matrix", ManualTransformationMatrix, FilterParameter::Parameter, ApplyTransformationToGeometry, 2));
  }
  parameters.push_back(SIMPL_NEW_FLOAT_FP("Rotation Angle (Degrees)", RotationAngle, FilterParameter::Parameter, ApplyTransformationToGeometry, 3));
  parameters.push_back(SIMPL_NEW_FLOAT_VEC3_FP("Rotation Axis (ijk)", RotationAxis, FilterParameter::Parameter, ApplyTransformationToGeometry, 3));
  parameters.push_back(SIMPL_NEW_FLOAT_VEC3_FP("Translation", Translation, FilterParameter::Parameter, ApplyTransformationToGeometry, 4));
  parameters.push_back(SIMPL_NEW_FLOAT_VEC3_FP("Scale", Scale, FilterParameter::Parameter, ApplyTransformationToGeometry, 5));
  DataContainerSelectionFilterParameter::RequirementType dcReq;
  IGeometry::Types geomTypes = {IGeometry::Type::Vertex, IGeometry::Type::Edge, IGeometry::Type::Triangle, IGeometry::Type::Quad, IGeometry::Type::Tetrahedral};
  dcReq.dcGeometryTypes = geomTypes;
  parameters.push_back(SIMPL_NEW_DC_SELECTION_FP("Geometry to Transform", GeometryToTransform, FilterParameter::RequiredArray, ApplyTransformationToGeometry, dcReq));
  {
    DataArraySelectionFilterParameter::RequirementType dasReq =
        DataArraySelectionFilterParameter::CreateRequirement(SIMPL::TypeNames::Float, SIMPL::Defaults::AnyComponentSize, AttributeMatrix::Type::Generic, IGeometry::Type::Any);
    parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Transformation Matrix", ComputedTransformationMatrix, FilterParameter::RequiredArray, ApplyTransformationToGeometry, dasReq, 1));
  }
  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ApplyTransformationToGeometry::readFilterParameters(AbstractFilterParametersReader* reader, int index)
{
  reader->openFilterGroup(this, index);
  setManualTransformationMatrix(reader->readDynamicTableData("ManualTransformationMatrix", getManualTransformationMatrix()));
  setComputedTransformationMatrix(reader->readDataArrayPath("ComputedTransformationMatrix", getComputedTransformationMatrix()));
  setGeometryToTransform(reader->readString("GeometryToTransform", getGeometryToTransform()));
  setTransformationMatrixType(reader->readValue("TransformationMatrixType", getTransformationMatrixType()));
  setRotationAxis(reader->readFloatVec3("RotationAxis", getRotationAxis()));
  setRotationAngle(reader->readValue("RotationAngle", getRotationAngle()));
  setTranslation(reader->readFloatVec3("Translation", getTranslation()));
  setScale(reader->readFloatVec3("Scale", getScale()));
  reader->closeFilterGroup();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ApplyTransformationToGeometry::dataCheck()
{
  setErrorCondition(0);
  setWarningCondition(0);

  IGeometry::Pointer igeom = getDataContainerArray()->getPrereqGeometryFromDataContainer<IGeometry, AbstractFilter>(this, getGeometryToTransform());

  if(getErrorCondition() < 0)
  {
    return;
  }

  if(!std::dynamic_pointer_cast<IGeometry2D>(igeom) && !std::dynamic_pointer_cast<IGeometry3D>(igeom) && !std::dynamic_pointer_cast<VertexGeom>(igeom) && !std::dynamic_pointer_cast<EdgeGeom>(igeom))
  {
    setErrorCondition(-702);
    QString ss =
        QObject::tr("Geometry to transform must be an unstructured geometry (Vertex, Edge, Triangle, Quadrilateral, or Tetrahedral), but the type is %1").arg(igeom->getGeometryTypeAsString());
    notifyErrorMessage(getHumanLabel(), ss, getErrorCondition());
  }

  QVector<size_t> cDims = {4, 4};

  switch(getTransformationMatrixType())
  {
  case 0: // No-Op
  {
    QString ss = QObject::tr("No transformation has been selected, so this filter will perform no operations");
    setWarningCondition(-701);
    notifyWarningMessage(getHumanLabel(), ss, getErrorCondition());
  }
  case 1: // Transformation matrix from array
  {
    m_TransformationMatrixPtr = getDataContainerArray()->getPrereqArrayFromPath<DataArray<float>, AbstractFilter>(this, getComputedTransformationMatrix(), cDims);
    if(m_TransformationMatrixPtr.lock())
    {
      m_TransformationMatrix = m_TransformationMatrixPtr.lock()->getPointer(0);
    }
    break;
  }
  case 2: // Manual transformation matrix
  {
    if(getManualTransformationMatrix().getNumRows() != 4)
    {
      setErrorCondition(-702);
      QString ss = QObject::tr("Manually entered transformation matrix must have exactly 4 rows");
      notifyErrorMessage(getHumanLabel(), ss, getErrorCondition());
      return;
    }
    if(getManualTransformationMatrix().getNumCols() != 4)
    {
      setErrorCondition(-703);
      QString ss = QObject::tr("Manually entered transformation matrix must have exactly 4 columns");
      notifyErrorMessage(getHumanLabel(), ss, getErrorCondition());
      return;
    }
    std::vector<std::vector<double>> tableData = getManualTransformationMatrix().getTableData();
    m_TransformationReference = FloatArrayType::CreateArray(1, cDims, "_INTERNAL_USE_ONLY_ManualTransformationMatrix");
    m_TransformationReference->initializeWithZeros();
    m_TransformationMatrixPtr = m_TransformationReference;
    if(m_TransformationMatrixPtr.lock())
    {
      m_TransformationMatrix = m_TransformationMatrixPtr.lock()->getPointer(0);
      for(size_t i = 0; i < tableData.size(); i++)
      {
        std::vector<double> row = tableData[i];
        for(size_t j = 0; j < row.size(); j++)
        {
          m_TransformationMatrix[4 * i + j] = static_cast<float>(row[j]);
        }
      }
    }
    break;
  }
  case 3: // Rotation via axis-angle
  {
    float rotAngle = m_RotationAngle * SIMPLib::Constants::k_Pi / 180.0;
    FOrientArrayType om(9);
    FOrientTransformsType::ax2om(FOrientArrayType(m_RotationAxis.x, m_RotationAxis.y, m_RotationAxis.z, rotAngle), om);

    m_TransformationReference = FloatArrayType::CreateArray(1, cDims, "_INTERNAL_USE_ONLY_ManualTransformationMatrix");
    m_TransformationReference->initializeWithZeros();
    m_TransformationMatrixPtr = m_TransformationReference;
    if(m_TransformationMatrixPtr.lock())
    {
      m_TransformationMatrix = m_TransformationMatrixPtr.lock()->getPointer(0);
      for(size_t i = 0; i < 3; i++)
      {
        m_TransformationMatrix[4 * i + 0] = om[3 * i + 0];
        m_TransformationMatrix[4 * i + 1] = om[3 * i + 1];
        m_TransformationMatrix[4 * i + 2] = om[3 * i + 2];
        m_TransformationMatrix[4 * i + 3] = 0.0f;
      }
      m_TransformationMatrix[4 * 3 + 3] = 1.0f;
    }
    break;
  }
  case 4: // Translation
  {
    m_TransformationReference = FloatArrayType::CreateArray(1, cDims, "_INTERNAL_USE_ONLY_ManualTransformationMatrix");
    m_TransformationReference->initializeWithZeros();
    m_TransformationMatrixPtr = m_TransformationReference;
    if(m_TransformationMatrixPtr.lock())
    {
      m_TransformationMatrix = m_TransformationMatrixPtr.lock()->getPointer(0);
      m_TransformationMatrix[4 * 0 + 0] = 1.0f;
      m_TransformationMatrix[4 * 1 + 1] = 1.0f;
      m_TransformationMatrix[4 * 2 + 2] = 1.0f;
      m_TransformationMatrix[4 * 0 + 3] = m_Translation.x;
      m_TransformationMatrix[4 * 1 + 3] = m_Translation.y;
      m_TransformationMatrix[4 * 2 + 3] = m_Translation.z;
      m_TransformationMatrix[4 * 3 + 3] = 1.0f;
    }
    break;
  }
  case 5: // Scale
  {
    m_TransformationReference = FloatArrayType::CreateArray(1, cDims, "_INTERNAL_USE_ONLY_ManualTransformationMatrix");
    m_TransformationReference->initializeWithZeros();
    m_TransformationMatrixPtr = m_TransformationReference;
    if(m_TransformationMatrixPtr.lock())
    {
      m_TransformationMatrix = m_TransformationMatrixPtr.lock()->getPointer(0);
      m_TransformationMatrix[4 * 0 + 0] = m_Scale.x;
      m_TransformationMatrix[4 * 1 + 1] = m_Scale.y;
      m_TransformationMatrix[4 * 2 + 2] = m_Scale.z;
      m_TransformationMatrix[4 * 3 + 3] = 1.0f;
    }
    break;
  }
  default:
  {
    QString ss = QObject::tr("Invalid selection for transformation type");
    setErrorCondition(-701);
    notifyErrorMessage(getHumanLabel(), ss, getErrorCondition());
    break;
  }
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ApplyTransformationToGeometry::preflight()
{
  // These are the REQUIRED lines of CODE to make sure the filter behaves correctly
  setInPreflight(true);              // Set the fact that we are preflighting.
  Q_EMIT preflightAboutToExecute();    // Emit this signal so that other widgets can do one file update
  Q_EMIT updateFilterParameters(this); // Emit this signal to have the widgets push their values down to the filter
  dataCheck();                       // Run our DataCheck to make sure everthing is setup correctly
  Q_EMIT preflightExecuted();          // We are done preflighting this filter
  setInPreflight(false);             // Inform the system this filter is NOT in preflight mode anymore.
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ApplyTransformationToGeometry::applyTransformation()
{
  IGeometry::Pointer igeom = getDataContainerArray()->getDataContainer(m_GeometryToTransform)->getGeometry();

  int64_t numVertices = 0;
  float* vertices = nullptr;

  if(IGeometry2D::Pointer igeom2D = std::dynamic_pointer_cast<IGeometry2D>(igeom))
  {
    numVertices = igeom2D->getNumberOfVertices();
    vertices = igeom2D->getVertexPointer(0);
  }
  else if(IGeometry3D::Pointer igeom3D = std::dynamic_pointer_cast<IGeometry3D>(igeom))
  {
    numVertices = igeom3D->getNumberOfVertices();
    vertices = igeom3D->getVertexPointer(0);
  }
  else if(VertexGeom::Pointer vertex = std::dynamic_pointer_cast<VertexGeom>(igeom))
  {
    numVertices = vertex->getNumberOfVertices();
    vertices = vertex->getVertexPointer(0);
  }
  else if(EdgeGeom::Pointer edge = std::dynamic_pointer_cast<EdgeGeom>(igeom))
  {
    numVertices = edge->getNumberOfVertices();
    vertices = vertex->getVertexPointer(0);
  }
  else
  {
    return;
  }

  typedef Eigen::Matrix<float, 4, 4, Eigen::RowMajor> ProjectiveMatrix;
  Eigen::Map<ProjectiveMatrix> transformation(m_TransformationMatrix);

  int64_t progIncrement = numVertices / 100;
  int64_t prog = 1;
  int64_t progressInt = 0;
  int64_t counter = 0;

  for(int64_t i = 0; i < numVertices; i++)
  {
    if(getCancel())
    {
      return;
    }
    Eigen::Vector4f position(vertices[3 * i + 0], vertices[3 * i + 1], vertices[3 * i + 2], 1);
    Eigen::Vector4f transformedPosition = transformation * position;
    std::memcpy(vertices + (3 * i), transformedPosition.data(), sizeof(float) * 3);

    if(counter > prog)
    {
      progressInt = static_cast<int64_t>((static_cast<float>(counter) / numVertices) * 100.0f);
      QString ss = QObject::tr("Transforming Geometry || %1% Completed").arg(progressInt);
      notifyStatusMessage(getMessagePrefix(), getHumanLabel(), ss);
      prog = prog + progIncrement;
    }
    counter++;
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ApplyTransformationToGeometry::execute()
{
  setErrorCondition(0);
  setWarningCondition(0);
  dataCheck();
  if(getErrorCondition() < 0)
  {
    return;
  }
  if(getWarningCondition() < 0)
  {
    return;
  }

  if(m_TransformationMatrixType == 0)
  {
    return;
  }
  
  
    applyTransformation();
  

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AbstractFilter::Pointer ApplyTransformationToGeometry::newFilterInstance(bool copyFilterParameters) const
{
  ApplyTransformationToGeometry::Pointer filter = ApplyTransformationToGeometry::New();
  if(copyFilterParameters)
  {
    copyFilterParameterInstanceVariables(filter.get());
  }
  return filter;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString ApplyTransformationToGeometry::getCompiledLibraryName() const
{
  return DREAM3DReviewConstants::DREAM3DReviewBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString ApplyTransformationToGeometry::getBrandingString() const
{
  return "DREAM3DReview";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString ApplyTransformationToGeometry::getFilterVersion() const
{
  QString version;
  QTextStream vStream(&version);
  vStream << DREAM3DReview::Version::Major() << "." << DREAM3DReview::Version::Minor() << "." << DREAM3DReview::Version::Patch();
  return version;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString ApplyTransformationToGeometry::getGroupName() const
{
  return DREAM3DReviewConstants::FilterGroups::DREAM3DReviewFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QUuid ApplyTransformationToGeometry::getUuid()
{
  return QUuid("{c681caf4-22f2-5885-bbc9-a0476abc72eb}");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString ApplyTransformationToGeometry::getSubGroupName() const
{
  return DREAM3DReviewConstants::FilterSubGroups::RotationTransformationFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString ApplyTransformationToGeometry::getHumanLabel() const
{
  return "Apply Transformation to Geometry";
}
