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

#include "ComputeUmeyamaTransform.h"

#include "SIMPLib/Common/Constants.h"
#include "SIMPLib/FilterParameters/AbstractFilterParametersReader.h"
#include "SIMPLib/FilterParameters/BooleanFilterParameter.h"
#include "SIMPLib/FilterParameters/DataContainerSelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/SeparatorFilterParameter.h"
#include "SIMPLib/FilterParameters/StringFilterParameter.h"
#include "SIMPLib/Geometry/EdgeGeom.h"
#include "SIMPLib/Geometry/IGeometry2D.h"
#include "SIMPLib/Geometry/IGeometry3D.h"
#include "SIMPLib/Geometry/VertexGeom.h"

#include <Eigen/Dense>
#include <Eigen/Geometry>

#include "DREAM3DReview/DREAM3DReviewConstants.h"
#include "DREAM3DReview/DREAM3DReviewVersion.h"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
ComputeUmeyamaTransform::ComputeUmeyamaTransform()
: m_SourcePointSet("")
, m_DestPointSet("")
, m_UseScaling(false)
, m_TransformationAttributeMatrixName("TransformationData")
, m_TransformationMatrixName("TransformationMatrix")
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
ComputeUmeyamaTransform::~ComputeUmeyamaTransform()
= default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ComputeUmeyamaTransform::setupFilterParameters()
{
  FilterParameterVector parameters;
  parameters.push_back(SIMPL_NEW_BOOL_FP("Use Scaling", UseScaling, FilterParameter::Parameter, ComputeUmeyamaTransform));
  DataContainerSelectionFilterParameter::RequirementType dcReq;
  dcReq.dcGeometryTypes = {IGeometry::Type::Vertex, IGeometry::Type::Edge, IGeometry::Type::Triangle, IGeometry::Type::Quad, IGeometry::Type::Tetrahedral};
  parameters.push_back(SIMPL_NEW_DC_SELECTION_FP("Moving Geometry", SourcePointSet, FilterParameter::RequiredArray, ComputeUmeyamaTransform, dcReq));
  parameters.push_back(SIMPL_NEW_DC_SELECTION_FP("Fixed Geometry", DestPointSet, FilterParameter::RequiredArray, ComputeUmeyamaTransform, dcReq));
  parameters.push_back(SeparatorFilterParameter::New("Transformation Data", FilterParameter::CreatedArray));
  parameters.push_back(SIMPL_NEW_STRING_FP("Transformation Attribute Matrix", TransformationAttributeMatrixName, FilterParameter::CreatedArray, ComputeUmeyamaTransform));
  parameters.push_back(SIMPL_NEW_STRING_FP("Transformation Matrix", TransformationMatrixName, FilterParameter::CreatedArray, ComputeUmeyamaTransform));
  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ComputeUmeyamaTransform::readFilterParameters(AbstractFilterParametersReader* reader, int index)
{
  reader->openFilterGroup(this, index);
  setSourcePointSet(reader->readString("SourcePointSet", getSourcePointSet()));
  setDestPointSet(reader->readString("DestPointSet", getDestPointSet()));
  setUseScaling(reader->readValue("UseScaling", getUseScaling()));
  setTransformationAttributeMatrixName(reader->readString("TransformationAttributeMatrixName", getTransformationAttributeMatrixName()));
  setTransformationMatrixName(reader->readString("TransformationMatrixName", getTransformationMatrixName()));
  reader->closeFilterGroup();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ComputeUmeyamaTransform::initialize()
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ComputeUmeyamaTransform::dataCheck()
{
  setErrorCondition(0);
  setWarningCondition(0);

  IGeometry::Pointer movingGeom = getDataContainerArray()->getPrereqGeometryFromDataContainer<IGeometry, AbstractFilter>(this, getSourcePointSet());
  IGeometry::Pointer fixedGeom = getDataContainerArray()->getPrereqGeometryFromDataContainer<IGeometry, AbstractFilter>(this, getDestPointSet());

  if(getErrorCondition() < 0)
  {
    return;
  }

  if(!std::dynamic_pointer_cast<IGeometry2D>(movingGeom) && !std::dynamic_pointer_cast<IGeometry3D>(movingGeom) && !std::dynamic_pointer_cast<VertexGeom>(movingGeom) &&
     !std::dynamic_pointer_cast<EdgeGeom>(movingGeom))
  {
    setErrorCondition(-702);
    QString ss = QObject::tr("Moving Geometry must be an unstructured geometry (Vertex, Edge, Triangle, Quadrilateral, or Tetrahedral), but the type is %1").arg(movingGeom->getGeometryTypeAsString());
    notifyErrorMessage(getHumanLabel(), ss, getErrorCondition());
  }

  if(!std::dynamic_pointer_cast<IGeometry2D>(fixedGeom) && !std::dynamic_pointer_cast<IGeometry3D>(fixedGeom) && !std::dynamic_pointer_cast<VertexGeom>(fixedGeom) &&
     !std::dynamic_pointer_cast<EdgeGeom>(fixedGeom))
  {
    setErrorCondition(-702);
    QString ss = QObject::tr("Fixed Geometry must be an unstructured geometry (Vertex, Edge, Triangle, Quadrilateral, or Tetrahedral), but the type is %1").arg(fixedGeom->getGeometryTypeAsString());
    notifyErrorMessage(getHumanLabel(), ss, getErrorCondition());
  }

  size_t numMovingVertices = 0;
  size_t numFixedVertices = 0;

  if(IGeometry2D::Pointer igeom2D = std::dynamic_pointer_cast<IGeometry2D>(movingGeom))
  {
    numMovingVertices = igeom2D->getNumberOfVertices();
  }
  else if(IGeometry3D::Pointer igeom3D = std::dynamic_pointer_cast<IGeometry3D>(movingGeom))
  {
    numMovingVertices = igeom3D->getNumberOfVertices();
  }
  else if(VertexGeom::Pointer vertex = std::dynamic_pointer_cast<VertexGeom>(movingGeom))
  {
    numMovingVertices = vertex->getNumberOfVertices();
  }
  else if(EdgeGeom::Pointer edge = std::dynamic_pointer_cast<EdgeGeom>(movingGeom))
  {
    numMovingVertices = edge->getNumberOfVertices();
  }

  if(IGeometry2D::Pointer igeom2D = std::dynamic_pointer_cast<IGeometry2D>(fixedGeom))
  {
    numFixedVertices = igeom2D->getNumberOfVertices();
  }
  else if(IGeometry3D::Pointer igeom3D = std::dynamic_pointer_cast<IGeometry3D>(fixedGeom))
  {
    numFixedVertices = igeom3D->getNumberOfVertices();
  }
  else if(VertexGeom::Pointer vertex = std::dynamic_pointer_cast<VertexGeom>(fixedGeom))
  {
    numFixedVertices = vertex->getNumberOfVertices();
  }
  else if(EdgeGeom::Pointer edge = std::dynamic_pointer_cast<EdgeGeom>(fixedGeom))
  {
    numFixedVertices = edge->getNumberOfVertices();
  }

  if(numMovingVertices != numFixedVertices)
  {
    setErrorCondition(-11000);
    QString ss = QObject::tr("The moving and fixed Geometries must have the same number of Vertices; the number of moving Vertices is %1 and the number of fixed Vertices is %2")
                     .arg(numMovingVertices)
                     .arg(numFixedVertices);
    notifyErrorMessage(getHumanLabel(), ss, getErrorCondition());
  }

  DataContainer::Pointer m = getDataContainerArray()->getPrereqDataContainer(this, getSourcePointSet());

  if(getErrorCondition() < 0)
  {
    return;
  }

  QVector<size_t> tDims(1, 1);
  m->createNonPrereqAttributeMatrix(this, getTransformationAttributeMatrixName(), tDims, AttributeMatrix::Type::Generic);

  QVector<size_t> cDims(2, 4);
  DataArrayPath path(getSourcePointSet(), getTransformationAttributeMatrixName(), getTransformationMatrixName());

  m_TransformationMatrixPtr =
      getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<float>, AbstractFilter, float>(this, path, 0, cDims); /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
  if(nullptr != m_TransformationMatrixPtr.lock()) /* Validate the Weak Pointer wraps a non-nullptr pointer to a DataArray<T> object */
  {
    m_TransformationMatrix = m_TransformationMatrixPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ComputeUmeyamaTransform::preflight()
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
void ComputeUmeyamaTransform::execute()
{
  setErrorCondition(0);
  setWarningCondition(0);
  dataCheck();
  if(getErrorCondition() < 0)
  {
    return;
  }

  IGeometry::Pointer movingGeom = getDataContainerArray()->getDataContainer(m_SourcePointSet)->getGeometry();
  IGeometry::Pointer fixedGeom = getDataContainerArray()->getDataContainer(m_DestPointSet)->getGeometry();

  float* movingPointsPtr = nullptr;
  float* fixedPointsPtr = nullptr;
  size_t numMovingVertices = 0;
  size_t numFixedVertices = 0;

  if(IGeometry2D::Pointer igeom2D = std::dynamic_pointer_cast<IGeometry2D>(movingGeom))
  {
    movingPointsPtr = igeom2D->getVertexPointer(0);
    numMovingVertices = igeom2D->getNumberOfVertices();
  }
  else if(IGeometry3D::Pointer igeom3D = std::dynamic_pointer_cast<IGeometry3D>(movingGeom))
  {
    movingPointsPtr = igeom3D->getVertexPointer(0);
    numMovingVertices = igeom3D->getNumberOfVertices();
  }
  else if(VertexGeom::Pointer vertex = std::dynamic_pointer_cast<VertexGeom>(movingGeom))
  {
    movingPointsPtr = vertex->getVertexPointer(0);
    numMovingVertices = vertex->getNumberOfVertices();
  }
  else if(EdgeGeom::Pointer edge = std::dynamic_pointer_cast<EdgeGeom>(movingGeom))
  {
    movingPointsPtr = edge->getVertexPointer(0);
    numMovingVertices = edge->getNumberOfVertices();
  }

  if(IGeometry2D::Pointer igeom2D = std::dynamic_pointer_cast<IGeometry2D>(fixedGeom))
  {
    fixedPointsPtr = igeom2D->getVertexPointer(0);
    numFixedVertices = igeom2D->getNumberOfVertices();
  }
  else if(IGeometry3D::Pointer igeom3D = std::dynamic_pointer_cast<IGeometry3D>(fixedGeom))
  {
    fixedPointsPtr = igeom3D->getVertexPointer(0);
    numFixedVertices = igeom3D->getNumberOfVertices();
  }
  else if(VertexGeom::Pointer vertex = std::dynamic_pointer_cast<VertexGeom>(fixedGeom))
  {
    fixedPointsPtr = vertex->getVertexPointer(0);
    numFixedVertices = vertex->getNumberOfVertices();
  }
  else if(EdgeGeom::Pointer edge = std::dynamic_pointer_cast<EdgeGeom>(fixedGeom))
  {
    fixedPointsPtr = edge->getVertexPointer(0);
    numFixedVertices = edge->getNumberOfVertices();
  }

  // Eigen does best with things in column major, but we'll ultimately want
  // the matrix as row major; perform the actual computation in column
  // major and then transpose when we're ready to move it into a DataArray
  typedef Eigen::Matrix<float, 3, Eigen::Dynamic, Eigen::ColMajor> PointCloud;
  typedef Eigen::Matrix<float, 4, 4, Eigen::ColMajor> UmeyamaTransform;

  // Map our DataArray pointers directly to Eigen objects in memory to avoid making copies
  Eigen::Map<PointCloud> moving(movingPointsPtr, 3, numMovingVertices);
  Eigen::Map<PointCloud> fixed(fixedPointsPtr, 3, numFixedVertices);

  UmeyamaTransform transformMatrix = Eigen::umeyama(moving, fixed, m_UseScaling);

  // The matrix is in column major, transpose it in place to make it row major
  transformMatrix.transposeInPlace();
  float* umeyamaArray = &transformMatrix(0);

  // Copy the transformation matrix over to our DataArray
  // We know the size will be exactly sixteen...
  for(size_t i = 0; i < 16; i++)
  {
    m_TransformationMatrix[i] = umeyamaArray[i];
  }

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AbstractFilter::Pointer ComputeUmeyamaTransform::newFilterInstance(bool copyFilterParameters) const
{
  ComputeUmeyamaTransform::Pointer filter = ComputeUmeyamaTransform::New();
  if(copyFilterParameters)
  {
    copyFilterParameterInstanceVariables(filter.get());
  }
  return filter;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString ComputeUmeyamaTransform::getCompiledLibraryName() const
{
  return DREAM3DReviewConstants::DREAM3DReviewBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString ComputeUmeyamaTransform::getBrandingString() const
{
  return "DREAM3DReview";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString ComputeUmeyamaTransform::getFilterVersion() const
{
  QString version;
  QTextStream vStream(&version);
  vStream << DREAM3DReview::Version::Major() << "." << DREAM3DReview::Version::Minor() << "." << DREAM3DReview::Version::Patch();
  return version;
}
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString ComputeUmeyamaTransform::getGroupName() const
{
  return DREAM3DReviewConstants::FilterGroups::DREAM3DReviewFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QUuid ComputeUmeyamaTransform::getUuid()
{
  return QUuid("{3192d494-d1ec-5ee7-a345-e9963f02aaab}");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString ComputeUmeyamaTransform::getSubGroupName() const
{
  return DREAM3DReviewConstants::FilterSubGroups::RegistrationFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString ComputeUmeyamaTransform::getHumanLabel() const
{
  return "Compute Umeyama Transform";
}
