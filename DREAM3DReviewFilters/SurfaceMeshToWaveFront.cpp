/* ============================================================================
 * Copyright (c) 2020 BlueQuartz Software, LLC
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
 * Neither the names of any of the BlueQuartz Software contributors
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
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include "SurfaceMeshToWaveFront.h"

#include <fstream>

#include "SIMPLib/Common/Constants.h"
#include "SIMPLib/DataContainers/DataContainer.h"
#include "SIMPLib/DataContainers/DataContainerArray.h"
#include "SIMPLib/FilterParameters/DataContainerSelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/DataArraySelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/MultiDataArraySelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/SeparatorFilterParameter.h"
#include "SIMPLib/FilterParameters/OutputFileFilterParameter.h"
#include "SIMPLib/Geometry/TriangleGeom.h"
#include "SIMPLib/Utilities/FileSystemPathHelper.h"

#include "DREAM3DReview/DREAM3DReviewConstants.h"
#include "DREAM3DReview/DREAM3DReviewVersion.h"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
SurfaceMeshToWaveFront::SurfaceMeshToWaveFront()
: AbstractFilter()
{
  initialize();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
SurfaceMeshToWaveFront::~SurfaceMeshToWaveFront() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void SurfaceMeshToWaveFront::initialize()
{
  clearErrorCode();
  clearWarningCode();
  setCancel(false);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void SurfaceMeshToWaveFront::setupFilterParameters()
{
  FilterParameterVectorType parameters;

  parameters.push_back(SIMPL_NEW_OUTPUT_FILE_FP("Output Wavefront File", OutputWaveFrontFile, FilterParameter::Parameter, SurfaceMeshToWaveFront, "*.obj", "Wavefront Object File"));

  DataArraySelectionFilterParameter::RequirementType req;
  req.dcGeometryTypes = IGeometry::Types(1, IGeometry::Type::Triangle);
  req.amTypes = AttributeMatrix::Types(1, AttributeMatrix::Type::Vertex);
  req.daTypes = QVector<QString>(1, SIMPL::TypeNames::Int8);
  std::vector<std::vector<size_t>> comp;
  comp.push_back(std::vector<size_t>(1, 1));
  req.componentDimensions = comp;
  parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Surface Mesh Node Type", SurfaceMeshNodeTypeArrayPath, FilterParameter::RequiredArray, SurfaceMeshToWaveFront, req));

  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void SurfaceMeshToWaveFront::dataCheck()
{
  clearErrorCode();
  clearWarningCode();

  FileSystemPathHelper::CheckOutputFile(this, "Output Wavefront File", getOutputWaveFrontFile(), true);

  std::vector<size_t> dims = {1};
  m_SurfaceMeshNodeTypePtr = getDataContainerArray()->getPrereqArrayFromPath<Int8ArrayType>(this, getSurfaceMeshNodeTypeArrayPath(), dims);
  if(nullptr != m_SurfaceMeshNodeTypePtr.lock())
  {
    m_SurfaceMeshNodeType = m_SurfaceMeshNodeTypePtr.lock()->getPointer(0);
  }

  DataContainer::Pointer sm = getDataContainerArray()->getPrereqDataContainer(this, getSurfaceMeshNodeTypeArrayPath().getDataContainerName(), false);
  if(getErrorCode() < 0)
  {
    return;
  }

  TriangleGeom::Pointer triangleGeom = sm->getPrereqGeometry<TriangleGeom>(this);
  if(getErrorCode() < 0)
  {
    return;
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void SurfaceMeshToWaveFront::execute()
{
  initialize();
  dataCheck();
  if(getErrorCode() < 0)
  {
    return;
  }

  DataContainer::Pointer dc = getDataContainerArray()->getDataContainer(getSurfaceMeshNodeTypeArrayPath().getDataContainerName());
  TriangleGeom::Pointer triangleGeom = dc->getGeometryAs<TriangleGeom>();
  MeshIndexType numberOfTriangles = triangleGeom->getNumberOfTris();
  SharedTriList::Pointer triangles = triangleGeom->getTriangles();
  SharedVertexList::Pointer vertices = triangleGeom->getVertices();
  size_t numberOfVertices = triangleGeom->getNumberOfVertices();

  std::ofstream ss;
  ss.open(getOutputWaveFrontFile().toStdString().c_str());
  if(!ss.is_open())
  {
    QString ss = QObject::tr("Could not open output file %1 for writing.").arg(getOutputWaveFrontFile());
    setErrorCondition(-2000, ss);
    return;
  }

  ss << "# Vertices\n";

  // Dump the vertices
  for(size_t i = 0; i < numberOfVertices; i++)
  {
    float c0 = vertices->getComponent(i, 0);
    float c1 = vertices->getComponent(i, 1);
    float c2 = vertices->getComponent(i, 2);

    ss << "v " << c0 << " " << c1 << " " << c2 << "\n";
  }

  ss << "\n# Faces\n";

  // Dump the triangle faces
  for(size_t i = 0; i < numberOfTriangles; i++)
  {
    size_t c0 = triangles->getComponent(i, 0);
    size_t c1 = triangles->getComponent(i, 1);
    size_t c2 = triangles->getComponent(i, 2);

    // These vertex values that make up the face must be 1-based
    ss << "f " << c0 + 1 << " " << c1 + 1 << " " << c2 + 1 << "\n";
  }

  ss.close();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AbstractFilter::Pointer SurfaceMeshToWaveFront::newFilterInstance(bool copyFilterParameters) const
{
  SurfaceMeshToWaveFront::Pointer filter = SurfaceMeshToWaveFront::New();
  if(copyFilterParameters)
  {
    copyFilterParameterInstanceVariables(filter.get());
  }
  return filter;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString SurfaceMeshToWaveFront::getCompiledLibraryName() const
{
  return DREAM3DReviewConstants::DREAM3DReviewBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString SurfaceMeshToWaveFront::getBrandingString() const
{
  return "DREAM3DReview";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString SurfaceMeshToWaveFront::getFilterVersion() const
{
  QString version;
  QTextStream vStream(&version);
  vStream << DREAM3DReview::Version::Major() << "." << DREAM3DReview::Version::Minor() << "." << DREAM3DReview::Version::Patch();
  return version;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString SurfaceMeshToWaveFront::getGroupName() const
{
  return SIMPL::FilterGroups::Unsupported;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString SurfaceMeshToWaveFront::getSubGroupName() const
{
  return "DREAM3DReview";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString SurfaceMeshToWaveFront::getHumanLabel() const
{
  return "Surface Mesh To Wavefront";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QUuid SurfaceMeshToWaveFront::getUuid() const
{
  return QUuid("{5d4c38fe-af0a-5cb3-b7fa-fb9e57b2097f}");
}

// -----------------------------------------------------------------------------
SurfaceMeshToWaveFront::Pointer SurfaceMeshToWaveFront::NullPointer()
{
  return Pointer(static_cast<Self*>(nullptr));
}

// -----------------------------------------------------------------------------
std::shared_ptr<SurfaceMeshToWaveFront> SurfaceMeshToWaveFront::New()
{
  struct make_shared_enabler : public SurfaceMeshToWaveFront
  {
  };
  std::shared_ptr<make_shared_enabler> val = std::make_shared<make_shared_enabler>();
  val->setupFilterParameters();
  return val;
}

// -----------------------------------------------------------------------------
QString SurfaceMeshToWaveFront::getNameOfClass() const
{
  return QString("SurfaceMeshToWaveFront");
}

// -----------------------------------------------------------------------------
QString SurfaceMeshToWaveFront::ClassName()
{
  return QString("SurfaceMeshToWaveFront");
}

// -----------------------------------------------------------------------------
void SurfaceMeshToWaveFront::setOutputWaveFrontFile(const QString& value)
{
  m_OutputWaveFrontFile = value;
}

// -----------------------------------------------------------------------------
QString SurfaceMeshToWaveFront::getOutputWaveFrontFile() const
{
  return m_OutputWaveFrontFile;
}

// -----------------------------------------------------------------------------
void SurfaceMeshToWaveFront::setSurfaceMeshNodeTypeArrayPath(const DataArrayPath& value)
{
  m_SurfaceMeshNodeTypeArrayPath = value;
}

// -----------------------------------------------------------------------------
DataArrayPath SurfaceMeshToWaveFront::getSurfaceMeshNodeTypeArrayPath() const
{
  return m_SurfaceMeshNodeTypeArrayPath;
}
