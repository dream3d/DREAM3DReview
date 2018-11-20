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

#include "PointSampleTriangleGeometry.h"

#include <cassert>
#include <cstring>
#include <ctime>

#include "SIMPLib/Common/Constants.h"
#include "SIMPLib/Common/TemplateHelpers.h"
#include "SIMPLib/FilterParameters/AbstractFilterParametersReader.h"
#include "SIMPLib/FilterParameters/DataArraySelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/DataContainerSelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/IntFilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedBooleanFilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedChoicesFilterParameter.h"
#include "SIMPLib/FilterParameters/MultiDataArraySelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/SeparatorFilterParameter.h"
#include "SIMPLib/FilterParameters/StringFilterParameter.h"
#include "SIMPLib/Geometry/EdgeGeom.h"
#include "SIMPLib/Geometry/IGeometry2D.h"
#include "SIMPLib/Geometry/IGeometry3D.h"
#include "SIMPLib/Geometry/IGeometryGrid.h"
#include "SIMPLib/Geometry/TriangleGeom.h"

#include "DREAM3DReview/DREAM3DReviewConstants.h"
#include "DREAM3DReview/DREAM3DReviewVersion.h"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
PointSampleTriangleGeometry::PointSampleTriangleGeometry()
: m_SamplesNumberType(0)
, m_TriangleGeometry(SIMPL::Defaults::TriangleDataContainerName)
, m_VertexGeometry(SIMPL::Defaults::VertexDataContainerName)
, m_VertexAttributeMatrixName(SIMPL::Defaults::VertexAttributeMatrixName)
, m_NumberOfSamples(100000)
, m_ParentGeometry("")
, m_TriangleAreasArrayPath(SIMPL::Defaults::TriangleDataContainerName, SIMPL::Defaults::FaceAttributeMatrixName, SIMPL::FaceData::SurfaceMeshFaceAreas)
, m_UseMask(false)
, m_MaskArrayPath("", "", "")
, m_SelectedDataArrayPaths(QVector<DataArrayPath>())
, m_NumSamples(0)
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
PointSampleTriangleGeometry::~PointSampleTriangleGeometry()
= default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void PointSampleTriangleGeometry::setupFilterParameters()
{
  FilterParameterVector parameters;
  {
    LinkedChoicesFilterParameter::Pointer parameter = LinkedChoicesFilterParameter::New();
    parameter->setHumanLabel("Source for Number of Samples");
    parameter->setPropertyName("SamplesNumberType");
    parameter->setSetterCallback(SIMPL_BIND_SETTER(PointSampleTriangleGeometry, this, SamplesNumberType));
    parameter->setGetterCallback(SIMPL_BIND_GETTER(PointSampleTriangleGeometry, this, SamplesNumberType));
    QVector<QString> choices;
    choices.push_back("Manual");
    choices.push_back("Other Geometry");
    parameter->setChoices(choices);
    QStringList linkedProps = {"NumberOfSamples", "ParentGeometry"};
    parameter->setLinkedProperties(linkedProps);
    parameter->setEditable(false);
    parameter->setCategory(FilterParameter::Parameter);
    parameters.push_back(parameter);
  }
  parameters.push_back(SIMPL_NEW_INTEGER_FP("Number of Sample Points", NumberOfSamples, FilterParameter::Parameter, PointSampleTriangleGeometry, 0));
  QStringList linkedProps;
  linkedProps << "MaskArrayPath";
  parameters.push_back(SIMPL_NEW_LINKED_BOOL_FP("Use Mask", UseMask, FilterParameter::Parameter, PointSampleTriangleGeometry, linkedProps));
  DataContainerSelectionFilterParameter::RequirementType dcsReq;
  IGeometry::Types geomTypes = {IGeometry::Type::Triangle};
  dcsReq.dcGeometryTypes = geomTypes;
  parameters.push_back(SIMPL_NEW_DC_SELECTION_FP("Triangle Geometry to Sample", TriangleGeometry, FilterParameter::RequiredArray, PointSampleTriangleGeometry, dcsReq));
  geomTypes = {IGeometry::Type::Image, IGeometry::Type::RectGrid, IGeometry::Type::Vertex, IGeometry::Type::Edge, IGeometry::Type::Triangle, IGeometry::Type::Quad, IGeometry::Type::Tetrahedral};
  dcsReq.dcGeometryTypes = geomTypes;
  parameters.push_back(SIMPL_NEW_DC_SELECTION_FP("Source Geometry for Number of Sample Points", ParentGeometry, FilterParameter::RequiredArray, PointSampleTriangleGeometry, dcsReq, 1));
  parameters.push_back(SeparatorFilterParameter::New("Face Data", FilterParameter::RequiredArray));
  DataArraySelectionFilterParameter::RequirementType dasReq = DataArraySelectionFilterParameter::CreateRequirement(SIMPL::TypeNames::Double, 1, AttributeMatrix::Type::Face, IGeometry::Type::Triangle);
  parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Face Areas", TriangleAreasArrayPath, FilterParameter::RequiredArray, PointSampleTriangleGeometry, dasReq));
  dasReq = DataArraySelectionFilterParameter::CreateRequirement(SIMPL::TypeNames::Bool, 1, AttributeMatrix::Type::Face, IGeometry::Type::Triangle);
  parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Mask", MaskArrayPath, FilterParameter::RequiredArray, PointSampleTriangleGeometry, dasReq));
  MultiDataArraySelectionFilterParameter::RequirementType mdaReq =
      MultiDataArraySelectionFilterParameter::CreateRequirement(SIMPL::Defaults::AnyPrimitive, SIMPL::Defaults::AnyComponentSize, AttributeMatrix::Type::Face, IGeometry::Type::Any);
  geomTypes = {IGeometry::Type::Triangle};
  mdaReq.dcGeometryTypes = geomTypes;
  parameters.push_back(SIMPL_NEW_MDA_SELECTION_FP("Attribute Arrays to Transfer", SelectedDataArrayPaths, FilterParameter::RequiredArray, PointSampleTriangleGeometry, mdaReq));
  parameters.push_back(SIMPL_NEW_STRING_FP("Vertex Geometry", VertexGeometry, FilterParameter::CreatedArray, PointSampleTriangleGeometry));
  parameters.push_back(SIMPL_NEW_STRING_FP("Vertex Attribute Matrix", VertexAttributeMatrixName, FilterParameter::CreatedArray, PointSampleTriangleGeometry));
  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void PointSampleTriangleGeometry::readFilterParameters(AbstractFilterParametersReader* reader, int index)
{
  reader->openFilterGroup(this, index);
  setSamplesNumberType(reader->readValue("SamplesNumberType", getSamplesNumberType()));
  setTriangleGeometry(reader->readString("TriangleGeometry", getTriangleGeometry()));
  setVertexGeometry(reader->readString("VertexGeometry", getVertexGeometry()));
  setVertexAttributeMatrixName(reader->readString("VertexAttributeMatrixName", getVertexAttributeMatrixName()));
  setNumberOfSamples(reader->readValue("NumberOfSamples", getNumberOfSamples()));
  setParentGeometry(reader->readString("ParentGeometry", getParentGeometry()));
  setTriangleAreasArrayPath(reader->readDataArrayPath("TriangleAreasArrayPath", getTriangleAreasArrayPath()));
  setUseMask(reader->readValue("UseMask", getUseMask()));
  setMaskArrayPath(reader->readDataArrayPath("MaskArrayPath", getMaskArrayPath()));
  reader->closeFilterGroup();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void PointSampleTriangleGeometry::initialize()
{
  m_SelectedWeakPtrVector.clear();
  m_CreatedWeakPtrVector.clear();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void PointSampleTriangleGeometry::dataCheck()
{
  setErrorCondition(0);
  setWarningCondition(0);
  initialize();

  switch(getSamplesNumberType())
  {
  case 0: // Manual
  {
    m_NumSamples = getNumberOfSamples();
    break;
  }
  case 1: // From other Geometry
  {
    IGeometry::Pointer igeom = getDataContainerArray()->getPrereqGeometryFromDataContainer<IGeometry, AbstractFilter>(this, getParentGeometry());
    if(getErrorCondition() < 0)
    {
      return;
    }

    if(VertexGeom::Pointer vertex = std::dynamic_pointer_cast<VertexGeom>(igeom))
    {
      m_NumSamples = vertex->getNumberOfVertices();
    }
    else if(EdgeGeom::Pointer edge = std::dynamic_pointer_cast<EdgeGeom>(igeom))
    {
      m_NumSamples = edge->getNumberOfVertices();
    }
    else if(IGeometryGrid::Pointer grid = std::dynamic_pointer_cast<IGeometryGrid>(igeom))
    {
      m_NumSamples = grid->getNumberOfElements();
    }
    else if(IGeometry2D::Pointer igeom2D = std::dynamic_pointer_cast<IGeometry2D>(igeom))
    {
      m_NumSamples = igeom2D->getNumberOfVertices();
    }
    else if(IGeometry3D::Pointer igeom3D = std::dynamic_pointer_cast<IGeometry3D>(igeom))
    {
      m_NumSamples = igeom3D->getNumberOfVertices();
    }
    else
    {
      QString ss =
          QObject::tr("Source Geometry must be of type Image, RectilinearGrid, Vertex, Edge, Triangle, Quadrilateral, or Tetrahedral, but the type is %1").arg(igeom->getGeometryTypeAsString());
      setErrorCondition(-701);
      notifyErrorMessage(getHumanLabel(), ss, getErrorCondition());
    }
    break;
  }
  default:
  {
    QString ss = QObject::tr("Invalid selection for determining the number of samples");
    setErrorCondition(-701);
    notifyErrorMessage(getHumanLabel(), ss, getErrorCondition());
    break;
  }
  }

  if(getErrorCondition() < 0)
  {
    return;
  }

  if(getSamplesNumberType() == 0 && getNumberOfSamples() <= 0)
  {
    QString ss = QObject::tr("Number of sample points must be greater than 0");
    setErrorCondition(-700);
    notifyErrorMessage(getHumanLabel(), ss, getErrorCondition());
  }

  QVector<IDataArray::Pointer> dataArrays;

  TriangleGeom::Pointer triangle = getDataContainerArray()->getPrereqGeometryFromDataContainer<TriangleGeom, AbstractFilter>(this, getTriangleGeometry());
  DataContainer::Pointer m = getDataContainerArray()->createNonPrereqDataContainer<AbstractFilter>(this, getVertexGeometry());

  if(getErrorCondition() < 0)
  {
    return;
  }

  dataArrays.push_back(triangle->getTriangles());

  VertexGeom::Pointer vertex = VertexGeom::CreateGeometry(m_NumSamples, SIMPL::Geometry::VertexGeometry, !getInPreflight());
  m->setGeometry(vertex);

  QVector<size_t> tDims(1, m_NumSamples);
  QVector<size_t> cDims(1, 1);

  m->createNonPrereqAttributeMatrix(this, getVertexAttributeMatrixName(), tDims, AttributeMatrix::Type::Vertex);

  m_TriangleAreasPtr = getDataContainerArray()->getPrereqArrayFromPath<DataArray<double>, AbstractFilter>(this, getTriangleAreasArrayPath(), cDims);
  if(m_TriangleAreasPtr.lock())
  {
    m_TriangleAreas = m_TriangleAreasPtr.lock()->getPointer(0);
  }
  if(getErrorCondition() >= 0)
  {
    dataArrays.push_back(m_TriangleAreasPtr.lock());
  }

  if(getUseMask())
  {
    m_MaskPtr = getDataContainerArray()->getPrereqArrayFromPath<DataArray<bool>, AbstractFilter>(this, getMaskArrayPath(), cDims);
    if(m_MaskPtr.lock())
    {
      m_Mask = m_MaskPtr.lock()->getPointer(0);
    }
    if(getErrorCondition() >= 0)
    {
      dataArrays.push_back(m_MaskPtr.lock());
    }
  }

  if(getErrorCondition() < 0)
  {
    return;
  }

  QVector<DataArrayPath> paths = getSelectedDataArrayPaths();

  if(!DataArrayPath::ValidateVector(paths))
  {
    setErrorCondition(-11004);
    QString ss = QObject::tr("There are Attribute Arrays selected that are not contained in the same Attribute Matrix; all selected Attribute Arrays must belong to the same Attribute Matrix");
    notifyErrorMessage(getHumanLabel(), ss, getErrorCondition());
  }

  for(auto&& path : paths)
  {
    IDataArray::WeakPointer ptr = getDataContainerArray()->getPrereqIDataArrayFromPath<IDataArray, AbstractFilter>(this, path);
    if(getErrorCondition() >= 0)
    {
      dataArrays.push_back(ptr.lock());
      m_SelectedWeakPtrVector.push_back(ptr);
      QVector<size_t> cDims = ptr.lock()->getComponentDimensions();
      DataArrayPath newPath(getVertexGeometry(), getVertexAttributeMatrixName(), ptr.lock()->getName());
      IDataArray::WeakPointer newPtr = TemplateHelpers::CreateNonPrereqArrayFromArrayType()(this, newPath, cDims, ptr.lock());
      m_CreatedWeakPtrVector.push_back(newPtr);
    }
  }

  getDataContainerArray()->validateNumberOfTuples<AbstractFilter>(this, dataArrays);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void PointSampleTriangleGeometry::preflight()
{
  // These are the REQUIRED lines of CODE to make sure the filter behaves correctly
  setInPreflight(true);              // Set the fact that we are preflighting.
  emit preflightAboutToExecute();    // Emit this signal so that other widgets can do one file update
  emit updateFilterParameters(this); // Emit this signal to have the widgets push their values down to the filter
  dataCheck();                       // Run our DataCheck to make sure everthing is setup correctly
  emit preflightExecuted();          // We are done preflighting this filter
  setInPreflight(false);             // Inform the system this filter is NOT in preflight mode anymore.
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
template <typename T> void copyDataToPoints(IDataArray::Pointer source, IDataArray::Pointer dest, int64_t sourceIdx, int64_t destIdx)
{
  typename DataArray<T>::Pointer sourcePtr = std::dynamic_pointer_cast<DataArray<T>>(source);
  typename DataArray<T>::Pointer destPtr = std::dynamic_pointer_cast<DataArray<T>>(dest);
  T* sPtr = sourcePtr->getPointer(0);
  T* dPtr = destPtr->getPointer(0);

  assert(sourcePtr->getNumberOfComponents() == destPtr->getNumberOfComponents());

  std::memcpy(dPtr + (destIdx * destPtr->getNumberOfComponents()), sPtr + (sourceIdx * sourcePtr->getNumberOfComponents()), sizeof(T) * sourcePtr->getNumberOfComponents());
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void PointSampleTriangleGeometry::sampleTriangle(float a[3], float b[3], float c[3], 
                                                  int64_t curVertex, VertexGeom::Pointer vertex, 
                                                  int64_t tri, std::mt19937_64& gen, 
                                                  std::uniform_real_distribution<>& dist)
{
  float* vertices = vertex->getVertexPointer(curVertex);

  float r1 = dist(gen);
  float r2 = dist(gen);

  float prefactorA = 1.0f - sqrtf(r1);
  float prefactorB = sqrtf(r1) * (1 - r2);
  float prefactorC = sqrtf(r1) * r2;

  vertices[0] = (prefactorA * a[0]) + (prefactorB * b[0]) + (prefactorC * c[0]);
  vertices[1] = (prefactorA * a[1]) + (prefactorB * b[1]) + (prefactorC * c[1]);
  vertices[2] = (prefactorA * a[2]) + (prefactorB * b[2]) + (prefactorC * c[2]);

  for(std::vector<IDataArray::WeakPointer>::size_type i = 0; i < m_SelectedWeakPtrVector.size(); i++)
  {
    EXECUTE_FUNCTION_TEMPLATE(this, copyDataToPoints, m_SelectedWeakPtrVector[i].lock(), m_SelectedWeakPtrVector[i].lock(), m_CreatedWeakPtrVector[i].lock(), tri, curVertex);
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void PointSampleTriangleGeometry::execute()
{
  setErrorCondition(0);
  setWarningCondition(0);
  dataCheck();
  if(getErrorCondition() < 0)
  {
    return;
  }

  TriangleGeom::Pointer triangle = getDataContainerArray()->getDataContainer(m_TriangleGeometry)->getGeometryAs<TriangleGeom>();
  int64_t numTris = triangle->getNumberOfTris();

  AttributeMatrix::Pointer attrMat = getDataContainerArray()->getDataContainer(m_VertexGeometry)->getAttributeMatrix(m_VertexAttributeMatrixName);
  QVector<size_t> tDims(1, m_NumSamples);
  attrMat->resizeAttributeArrays(tDims);

  VertexGeom::Pointer vertex = getDataContainerArray()->getDataContainer(m_VertexGeometry)->getGeometryAs<VertexGeom>();
  vertex->resizeVertexList(m_NumSamples);

  std::mt19937_64::result_type seed = static_cast<std::mt19937_64::result_type>(std::chrono::steady_clock::now().time_since_epoch().count());
  std::mt19937_64 generator(seed);
  std::uniform_real_distribution<> distribution(0.0f, 1.0f);

  std::vector<int64_t> triangleIds(numTris);
  std::iota(std::begin(triangleIds), std::end(triangleIds), 0);

  std::vector<double> triangleWeights(numTris);
  std::memcpy(triangleWeights.data(), m_TriangleAreas, sizeof(double) * numTris);

  // VS2013 does not implement the iterator contructor for std::discrete_distribution<>, which
  // really is a massively idiotic oversight; hack the equivalent using the unary_op constructor
  std::discrete_distribution<size_t> triangle_distribution(triangleWeights.size(), -0.5, -0.5 + triangleWeights.size(),
                                                           [&triangleWeights](double index) { return triangleWeights[static_cast<size_t>(index)]; });

  int64_t progIncrement = m_NumSamples / 100;
  int64_t prog = 1;
  int64_t progressInt = 0;
  int64_t counter = 0;

  float a[3] = {0.0f, 0.0f, 0.0f};
  float b[3] = {0.0f, 0.0f, 0.0f};
  float c[3] = {0.0f, 0.0f, 0.0f};

  assert(m_SelectedWeakPtrVector.size() == m_CreatedWeakPtrVector.size());

  for(int64_t i = 0; i < m_NumSamples; i++)
  {
    if(getCancel())
    {
      return;
    }

    int64_t randomTri = 0;

    if(m_UseMask)
    {
      while(true)
      {
        randomTri = triangle_distribution(generator);
        assert(randomTri < numTris && randomTri >= 0);
        if(m_Mask[randomTri])
        {
          break;
        }
      }
    }
    else
    {
      randomTri = triangle_distribution(generator);
      assert(randomTri < numTris && randomTri >= 0);
    }

    triangle->getVertCoordsAtTri(randomTri, a, b, c);
    sampleTriangle(a, b, c, i, vertex, randomTri, generator, distribution);

    if(counter > prog)
    {
      progressInt = static_cast<int64_t>((static_cast<float>(counter) / m_NumSamples) * 100.0f);
      QString ss = QObject::tr("Sampling Triangles || %1% Completed").arg(progressInt);
      notifyStatusMessage(getMessagePrefix(), getHumanLabel(), ss);
      prog = prog + progIncrement;
    }
    counter++;
  }

  notifyStatusMessage(getHumanLabel(), "Complete");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AbstractFilter::Pointer PointSampleTriangleGeometry::newFilterInstance(bool copyFilterParameters) const
{
  PointSampleTriangleGeometry::Pointer filter = PointSampleTriangleGeometry::New();
  if(copyFilterParameters)
  {
    copyFilterParameterInstanceVariables(filter.get());
  }
  return filter;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString PointSampleTriangleGeometry::getCompiledLibraryName() const
{
  return DREAM3DReviewConstants::DREAM3DReviewBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString PointSampleTriangleGeometry::getBrandingString() const
{
  return "DREAM3DReview";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString PointSampleTriangleGeometry::getFilterVersion() const
{
  QString version;
  QTextStream vStream(&version);
  vStream << DREAM3DReview::Version::Major() << "." << DREAM3DReview::Version::Minor() << "." << DREAM3DReview::Version::Patch();
  return version;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString PointSampleTriangleGeometry::getGroupName() const
{
  return DREAM3DReviewConstants::FilterGroups::DREAM3DReviewFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QUuid PointSampleTriangleGeometry::getUuid()
{
  return QUuid("{119861c5-e303-537e-b210-2e62936222e9}");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString PointSampleTriangleGeometry::getSubGroupName() const
{
  return DREAM3DReviewConstants::FilterSubGroups::GeometryFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString PointSampleTriangleGeometry::getHumanLabel() const
{
  return "Point Sample Triangle Geometry";
}
