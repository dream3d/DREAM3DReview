/* ============================================================================
* Software developed by US federal government employees (including military personnel) 
* as part of their official duties is not subject to copyright protection and is 
* considered “public domain” (see 17 USC Section 105). Public domain software can be used 
* by anyone for any purpose, and cannot be released under a copyright license 
* (including typical open source software licenses).
* 
* This source code file was originally written by United States DoD employees. The
* original source code files are released into the Public Domain.
* 
* Subsequent changes to the codes by others may elect to add a copyright and license
* for those changes.
* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
#include "ExtractTripleLinesFromTriangleGeometry.h"

#include <array>
#include <random>
#include <unordered_map>
#include <unordered_set>

#include "SIMPLib/Common/Constants.h"
#include "SIMPLib/FilterParameters/AbstractFilterParametersReader.h"
#include "SIMPLib/FilterParameters/BooleanFilterParameter.h"
#include "SIMPLib/FilterParameters/DataArraySelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/SeparatorFilterParameter.h"
#include "SIMPLib/FilterParameters/StringFilterParameter.h"
#include "SIMPLib/Geometry/EdgeGeom.h"
#include "SIMPLib/Geometry/GeometryHelpers.h"
#include "SIMPLib/Geometry/TriangleGeom.h"
#include "SIMPLib/Geometry/VertexGeom.h"

#include "DREAM3DReview/DREAM3DReviewConstants.h"
#include "DREAM3DReview/DREAM3DReviewVersion.h"

namespace
{
template <class T> inline void hashCombine(size_t& seed, const T& obj)
{
  std::hash<T> hasher;
  seed ^= hasher(obj) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

using Vertex = std::array<float, 3>;
using Edge = std::array<int64_t, 2>;

struct VertexHasher
{
  size_t operator()(const Vertex& vert) const
  {
    size_t hash = std::hash<float>()(vert[0]);
    hashCombine(hash, vert[1]);
    hashCombine(hash, vert[2]);
    return hash;
  }
};

struct EdgeHasher
{
  size_t operator()(const Edge& edge) const
  {
    size_t hash = std::hash<int64_t>()(edge[0]);
    hashCombine(hash, edge[1]);
    return hash;
  }
};

using VertexMap = std::unordered_map<Vertex, int64_t, VertexHasher>;
using EdgeMap = std::unordered_map<Edge, int64_t, EdgeHasher>;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
ExtractTripleLinesFromTriangleGeometry::ExtractTripleLinesFromTriangleGeometry()
: m_NodeTypesArrayPath(SIMPL::Defaults::TriangleDataContainerName, SIMPL::Defaults::VertexAttributeMatrixName, SIMPL::VertexData::SurfaceMeshNodeType)
, m_EdgeGeometry(SIMPL::Defaults::EdgeDataContainerName)
, m_VertexAttributeMatrixName(SIMPL::Defaults::VertexAttributeMatrixName)
, m_EdgeAttributeMatrixName(SIMPL::Defaults::EdgeAttributeMatrixName)
, m_NodeTypesArrayName(SIMPL::VertexData::SurfaceMeshNodeType)
, m_SmoothTripleLines(false)
, m_NodeTypes(nullptr)
, m_TripleLineNodeTypes(nullptr)
{
  initialize();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
ExtractTripleLinesFromTriangleGeometry::~ExtractTripleLinesFromTriangleGeometry() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ExtractTripleLinesFromTriangleGeometry::initialize()
{
  clearErrorCode();
  setCancel(false);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ExtractTripleLinesFromTriangleGeometry::setupFilterParameters()
{
  FilterParameterVectorType parameters;
  parameters.push_back(SIMPL_NEW_BOOL_FP("Compactify Triple Lines", SmoothTripleLines, FilterParameter::Parameter, ExtractTripleLinesFromTriangleGeometry));
  DataArraySelectionFilterParameter::RequirementType dasReq = DataArraySelectionFilterParameter::CreateRequirement(SIMPL::TypeNames::Int8, 1, AttributeMatrix::Type::Vertex, IGeometry::Type::Triangle);
  parameters.push_back(SeparatorFilterParameter::New("Vertex Data", FilterParameter::RequiredArray));
  parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Node Types", NodeTypesArrayPath, FilterParameter::RequiredArray, ExtractTripleLinesFromTriangleGeometry, dasReq));
  parameters.push_back(SIMPL_NEW_STRING_FP("Edge Geometry", EdgeGeometry, FilterParameter::CreatedArray, ExtractTripleLinesFromTriangleGeometry));
  parameters.push_back(SeparatorFilterParameter::New("Vertex Data", FilterParameter::CreatedArray));
  parameters.push_back(SIMPL_NEW_STRING_FP("Vertex Attribute Matrix", VertexAttributeMatrixName, FilterParameter::CreatedArray, ExtractTripleLinesFromTriangleGeometry));
  parameters.push_back(SIMPL_NEW_STRING_FP("Node Types", NodeTypesArrayName, FilterParameter::CreatedArray, ExtractTripleLinesFromTriangleGeometry));
  parameters.push_back(SeparatorFilterParameter::New("Edge Data", FilterParameter::CreatedArray));
  parameters.push_back(SIMPL_NEW_STRING_FP("Edge Attribute Matrix", EdgeAttributeMatrixName, FilterParameter::CreatedArray, ExtractTripleLinesFromTriangleGeometry));
  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ExtractTripleLinesFromTriangleGeometry::dataCheck()
{
  clearErrorCode();

  QVector<IDataArray::Pointer> dataArrays;

  TriangleGeom::Pointer triangle = getDataContainerArray()->getPrereqGeometryFromDataContainer<TriangleGeom, AbstractFilter>(this, getNodeTypesArrayPath().getDataContainerName());

  if(getErrorCode() < 0)
  {
    return;
  }

  dataArrays.push_back(triangle->getVertices());

  QVector<size_t> cDims(1, 1);
  m_NodeTypesPtr = getDataContainerArray()->getPrereqArrayFromPath<DataArray<int8_t>, AbstractFilter>(this, getNodeTypesArrayPath(),
                                                                                                      cDims); /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
  if(nullptr != m_NodeTypesPtr.lock().get())                                                                  /* Validate the Weak Pointer wraps a non-nullptr pointer to a DataArray<T> object */
  {
    m_NodeTypes = m_NodeTypesPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */
  if(getErrorCode() >= 0)
  {
    dataArrays.push_back(m_NodeTypesPtr.lock());
  }

  getDataContainerArray()->validateNumberOfTuples<AbstractFilter>(this, dataArrays);

  DataContainer::Pointer edc = getDataContainerArray()->createNonPrereqDataContainer<AbstractFilter>(this, getEdgeGeometry());

  if(getErrorCode() < 0)
  {
    return;
  }

  SharedVertexList::Pointer vertex = EdgeGeom::CreateSharedVertexList(0, !getInPreflight());
  EdgeGeom::Pointer edge = EdgeGeom::CreateGeometry(0, vertex, SIMPL::Geometry::EdgeGeometry, !getInPreflight());
  edc->setGeometry(edge);

  QVector<size_t> tDims(1, 0);
  edc->createNonPrereqAttributeMatrix(this, getVertexAttributeMatrixName(), tDims, AttributeMatrix::Type::Vertex);
  edc->createNonPrereqAttributeMatrix(this, getEdgeAttributeMatrixName(), tDims, AttributeMatrix::Type::Edge);

  DataArrayPath path(getEdgeGeometry(), getVertexAttributeMatrixName(), getNodeTypesArrayName());
  m_TripleLineNodeTypesPtr = getDataContainerArray()->createNonPrereqArrayFromPath<Int8ArrayType, AbstractFilter, int8_t>(this, path, 0, cDims);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ExtractTripleLinesFromTriangleGeometry::preflight()
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
bool checkBoxEdge(float minExtents[3], float maxExtents[3], float* coords)
{
  if(coords[0] == minExtents[0] && coords[1] == minExtents[1])
  {
    return true;
  }
  if(coords[0] == minExtents[0] && coords[1] == maxExtents[1])
  {
    return true;
  }
  if(coords[0] == minExtents[0] && coords[2] == minExtents[2])
  {
    return true;
  }
  if(coords[0] == minExtents[0] && coords[2] == maxExtents[2])
  {
    return true;
  }

  if(coords[0] == maxExtents[0] && coords[1] == minExtents[1])
  {
    return true;
  }
  if(coords[0] == maxExtents[0] && coords[1] == maxExtents[1])
  {
    return true;
  }
  if(coords[0] == maxExtents[0] && coords[2] == minExtents[2])
  {
    return true;
  }
  if(coords[0] == maxExtents[0] && coords[2] == maxExtents[2])
  {
    return true;
  }

  if(coords[1] == minExtents[1] && coords[2] == minExtents[2])
  {
    return true;
  }
  if(coords[1] == minExtents[1] && coords[2] == maxExtents[2])
  {
    return true;
  }
  if(coords[1] == maxExtents[1] && coords[2] == minExtents[2])
  {
    return true;
  }
  if(coords[1] == maxExtents[1] && coords[2] == maxExtents[2])
  {
    return true;
  }

  return false;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ExtractTripleLinesFromTriangleGeometry::extractTripleLines()
{
  TriangleGeom::Pointer triangle = getDataContainerArray()->getDataContainer(m_NodeTypesArrayPath.getDataContainerName())->getGeometryAs<TriangleGeom>();
  float* triVerts = triangle->getVertexPointer(0);
  int64_t numVerts = triangle->getNumberOfVertices();

  if(!triangle->getEdges())
  {
    triangle->findEdges();
  }

  SharedEdgeList::Pointer edges = triangle->getEdges();
  int64_t numEdges = triangle->getNumberOfEdges();
  int64_t* edgePtr = edges->getPointer(0);
  ElementDynamicList::Pointer edgesContainingVert = ElementDynamicList::New();
  GeometryHelpers::Connectivity::FindElementsContainingVert<uint16_t, int64_t>(edges, edgesContainingVert, numVerts);

  VertexMap vertexMap;
  EdgeMap edgeMap;
  std::vector<Vertex> tmpVerts;
  std::vector<Edge> tmpEdges;
  std::vector<int8_t> tmpNodeTypes;
  tmpVerts.reserve(numVerts);
  tmpEdges.reserve(numEdges);
  tmpNodeTypes.reserve(numVerts);
  int64_t vertCounter = 0;
  int64_t edgeCounter = 0;
  float minExtents[3] = {std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max()};
  float maxExtents[3] = {std::numeric_limits<float>::lowest(), std::numeric_limits<float>::lowest(), std::numeric_limits<float>::lowest()};

  for(int64_t i = 0; i < numVerts; i++)
  {
    if(triVerts[3 * i + 0] < minExtents[0])
    {
      minExtents[0] = triVerts[3 * i + 0];
    }
    if(triVerts[3 * i + 0] > maxExtents[0])
    {
      maxExtents[0] = triVerts[3 * i + 0];
    }
    if(triVerts[3 * i + 1] < minExtents[1])
    {
      minExtents[1] = triVerts[3 * i + 1];
    }
    if(triVerts[3 * i + 1] > maxExtents[1])
    {
      maxExtents[1] = triVerts[3 * i + 1];
    }
    if(triVerts[3 * i + 2] < minExtents[2])
    {
      minExtents[2] = triVerts[3 * i + 2];
    }
    if(triVerts[3 * i + 2] > maxExtents[2])
    {
      maxExtents[2] = triVerts[3 * i + 2];
    }
  }

  // TODO: generalize assumption that the surface of the mesh is an axis-aligned box
  std::vector<Vertex> corners;
  Vertex corner0 = {minExtents[0], minExtents[1], minExtents[2]};
  corners.push_back(corner0);
  Vertex corner1 = {maxExtents[0], minExtents[1], minExtents[2]};
  corners.push_back(corner1);
  Vertex corner2 = {maxExtents[0], maxExtents[1], minExtents[2]};
  corners.push_back(corner2);
  Vertex corner3 = {minExtents[0], maxExtents[1], minExtents[2]};
  corners.push_back(corner3);
  Vertex corner4 = {minExtents[0], minExtents[1], maxExtents[2]};
  corners.push_back(corner4);
  Vertex corner5 = {maxExtents[0], minExtents[1], maxExtents[2]};
  corners.push_back(corner5);
  Vertex corner6 = {maxExtents[0], maxExtents[1], maxExtents[2]};
  corners.push_back(corner6);
  Vertex corner7 = {minExtents[0], maxExtents[1], maxExtents[2]};
  corners.push_back(corner7);

  for(auto vert : corners)
  {
    tmpVerts.push_back(vert);
    tmpNodeTypes.push_back(12);
    vertexMap[vert] = vertCounter;
    vertCounter++;
  }

  int64_t progIncrement = numVerts / 100;
  int64_t prog = 1;
  int64_t progressInt = 0;

  for(int64_t i = 0; i < numVerts; i++)
  {
    if(m_NodeTypes[i] == 12)
    {
      if(checkBoxEdge(minExtents, maxExtents, triVerts + (3 * i)))
      {
        uint16_t numEdgesToVert = edgesContainingVert->getNumberOfElements(i);
        int64_t* edgesAtVert = edgesContainingVert->getElementListPointer(i);

        Vertex tmpCoords0 = {triVerts[3 * i + 0], triVerts[3 * i + 1], triVerts[3 * i + 2]};
        auto iter = vertexMap.find(tmpCoords0);
        if(iter == vertexMap.end())
        {
          tmpVerts.push_back(tmpCoords0);
          tmpNodeTypes.push_back(m_NodeTypes[i]);
          vertexMap[tmpCoords0] = vertCounter;
          vertCounter++;
        }

        for(uint16_t j = 0; j < numEdgesToVert; j++)
        {
          std::vector<std::pair<Vertex, uint8_t>> vertsToCheck;
          if(m_NodeTypes[edgePtr[2 * edgesAtVert[j] + 0]] == 3 || m_NodeTypes[edgePtr[2 * edgesAtVert[j] + 0]] == 4 || m_NodeTypes[edgePtr[2 * edgesAtVert[j] + 0]] == 13 ||
             m_NodeTypes[edgePtr[2 * edgesAtVert[j] + 0]] == 14 ||
             (m_NodeTypes[edgePtr[2 * edgesAtVert[j] + 0]] == 12 && checkBoxEdge(minExtents, maxExtents, triVerts + (3 * edgePtr[2 * edgesAtVert[j] + 0]))))
          {
            Vertex vert = {triVerts[3 * edgePtr[2 * edgesAtVert[j] + 0] + 0], triVerts[3 * edgePtr[2 * edgesAtVert[j] + 0] + 1], triVerts[3 * edgePtr[2 * edgesAtVert[j] + 0] + 2]};
            vertsToCheck.push_back({vert, 0});
          }
          if(m_NodeTypes[edgePtr[2 * edgesAtVert[j] + 1]] == 3 || m_NodeTypes[edgePtr[2 * edgesAtVert[j] + 1]] == 4 || m_NodeTypes[edgePtr[2 * edgesAtVert[j] + 1]] == 13 ||
             m_NodeTypes[edgePtr[2 * edgesAtVert[j] + 1]] == 14 ||
             (m_NodeTypes[edgePtr[2 * edgesAtVert[j] + 1]] == 12 && checkBoxEdge(minExtents, maxExtents, triVerts + (3 * edgePtr[2 * edgesAtVert[j] + 1]))))
          {
            Vertex vert = {triVerts[3 * edgePtr[2 * edgesAtVert[j] + 1] + 0], triVerts[3 * edgePtr[2 * edgesAtVert[j] + 1] + 1], triVerts[3 * edgePtr[2 * edgesAtVert[j] + 1] + 2]};
            vertsToCheck.push_back({vert, 1});
          }
          for(auto vert : vertsToCheck)
          {
            auto viter = vertexMap.find(vert.first);
            if(viter == vertexMap.end())
            {
              tmpVerts.push_back(vert.first);
              tmpNodeTypes.push_back(m_NodeTypes[edgePtr[2 * edgesAtVert[j] + vert.second]]);
              vertexMap[vert.first] = vertCounter;
              vertCounter++;
            }

            int64_t startVert = vertexMap[tmpCoords0];
            int64_t endVert = vertexMap[vert.first];
            if(startVert != endVert)
            {
              if(startVert > endVert)
              {
                std::swap(startVert, endVert);
              }
              Edge tmpEdge = {startVert, endVert};
              auto eiter = edgeMap.find(tmpEdge);
              if(eiter == edgeMap.end())
              {
                tmpEdges.push_back(tmpEdge);
                edgeMap[tmpEdge] = edgeCounter;
                edgeCounter++;
              }
            }
          }
        }
      }
    }

    if(m_NodeTypes[i] == 3 || m_NodeTypes[i] == 4 || m_NodeTypes[i] == 13 || m_NodeTypes[i] == 14)
    {
      uint16_t numEdgesToVert = edgesContainingVert->getNumberOfElements(i);
      int64_t* edgesAtVert = edgesContainingVert->getElementListPointer(i);

      Vertex tmpCoords0 = {triVerts[3 * i + 0], triVerts[3 * i + 1], triVerts[3 * i + 2]};
      auto iter = vertexMap.find(tmpCoords0);
      if(iter == vertexMap.end())
      {
        tmpVerts.push_back(tmpCoords0);
        tmpNodeTypes.push_back(m_NodeTypes[i]);
        vertexMap[tmpCoords0] = vertCounter;
        vertCounter++;
      }

      for(uint16_t j = 0; j < numEdgesToVert; j++)
      {
        std::vector<std::pair<Vertex, uint8_t>> vertsToCheck;
        if(m_NodeTypes[edgePtr[2 * edgesAtVert[j] + 0]] == 3 || m_NodeTypes[edgePtr[2 * edgesAtVert[j] + 0]] == 4 || m_NodeTypes[edgePtr[2 * edgesAtVert[j] + 0]] == 13 ||
           m_NodeTypes[edgePtr[2 * edgesAtVert[j] + 0]] == 14)
        {
          Vertex vert = {triVerts[3 * edgePtr[2 * edgesAtVert[j] + 0] + 0], triVerts[3 * edgePtr[2 * edgesAtVert[j] + 0] + 1], triVerts[3 * edgePtr[2 * edgesAtVert[j] + 0] + 2]};
          vertsToCheck.push_back({vert, 0});
        }
        if(m_NodeTypes[edgePtr[2 * edgesAtVert[j] + 1]] == 3 || m_NodeTypes[edgePtr[2 * edgesAtVert[j] + 1]] == 4 || m_NodeTypes[edgePtr[2 * edgesAtVert[j] + 1]] == 13 ||
           m_NodeTypes[edgePtr[2 * edgesAtVert[j] + 1]] == 14)
        {
          Vertex vert = {triVerts[3 * edgePtr[2 * edgesAtVert[j] + 1] + 0], triVerts[3 * edgePtr[2 * edgesAtVert[j] + 1] + 1], triVerts[3 * edgePtr[2 * edgesAtVert[j] + 1] + 2]};
          vertsToCheck.push_back({vert, 1});
        }
        for(auto vert : vertsToCheck)
        {
          auto viter = vertexMap.find(vert.first);
          if(viter == vertexMap.end())
          {
            tmpVerts.push_back(vert.first);
            tmpNodeTypes.push_back(m_NodeTypes[edgePtr[2 * edgesAtVert[j] + vert.second]]);
            vertexMap[vert.first] = vertCounter;
            vertCounter++;
          }

          int64_t startVert = vertexMap[tmpCoords0];
          int64_t endVert = vertexMap[vert.first];
          if(startVert != endVert)
          {
            if(startVert > endVert)
            {
              std::swap(startVert, endVert);
            }
            Edge tmpEdge = {startVert, endVert};
            auto eiter = edgeMap.find(tmpEdge);
            if(eiter == edgeMap.end())
            {
              tmpEdges.push_back(tmpEdge);
              edgeMap[tmpEdge] = edgeCounter;
              edgeCounter++;
            }
          }
        }
      }
    }

    if(i > prog)
    {
      progressInt = static_cast<int64_t>((static_cast<float>(i) / numVerts) * 100.0f);
      QString ss = QObject::tr("Extracting Triple Lines || %1% Completed").arg(progressInt);
      notifyStatusMessage(ss);
      prog = prog + progIncrement;
    }
  }

  EdgeGeom::Pointer tripleLineEdge = getDataContainerArray()->getDataContainer(m_EdgeGeometry)->getGeometryAs<EdgeGeom>();
  tripleLineEdge->resizeEdgeList(tmpEdges.size());
  tripleLineEdge->resizeVertexList(tmpVerts.size());
  float* tripleLineVerts = tripleLineEdge->getVertexPointer(0);
  int64_t* tripleLineEdges = tripleLineEdge->getEdgePointer(0);
  int64_t numTripleLineVerts = tripleLineEdge->getNumberOfVertices();
  int64_t numTripleLineEdges = tripleLineEdge->getNumberOfEdges();

  QString ss = QObject::tr("Building Edge Geoemetry...");
  notifyStatusMessage(ss);

  for(int64_t i = 0; i < numTripleLineVerts; i++)
  {
    tripleLineVerts[3 * i + 0] = tmpVerts[i][0];
    tripleLineVerts[3 * i + 1] = tmpVerts[i][1];
    tripleLineVerts[3 * i + 2] = tmpVerts[i][2];
  }

  for(int64_t i = 0; i < numTripleLineEdges; i++)
  {
    tripleLineEdges[2 * i + 0] = tmpEdges[i][0];
    tripleLineEdges[2 * i + 1] = tmpEdges[i][1];
  }

  AttributeMatrix::Pointer attrMat = getDataContainerArray()->getDataContainer(m_EdgeGeometry)->getAttributeMatrix(m_VertexAttributeMatrixName);
  QVector<size_t> tDims(1, numTripleLineVerts);
  attrMat->resizeAttributeArrays(tDims);
  attrMat = getDataContainerArray()->getDataContainer(m_EdgeGeometry)->getAttributeMatrix(m_EdgeAttributeMatrixName);
  tDims[0] = numTripleLineEdges;
  attrMat->resizeAttributeArrays(tDims);

  m_TripleLineNodeTypes = m_TripleLineNodeTypesPtr.lock()->getPointer(0);
  std::memcpy(m_TripleLineNodeTypes, tmpNodeTypes.data(), sizeof(int8_t) * numTripleLineVerts);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void burnVertices(std::vector<int64_t> vertsToBurn, int64_t seedVert, std::unordered_set<int64_t>& burnedVerts, int64_t& vertCounter, int64_t& edgeCounter, float* vertPtr, int64_t* edgePtr,
                  int8_t* nodeTypes, ElementDynamicList::Pointer edgesContainingVert, VertexMap& vertexMap, EdgeMap& edgeMap, std::vector<Vertex>& tmpVerts, std::vector<Edge>& tmpEdges,
                  std::vector<int8_t>& tmpNodeTypes, AbstractFilter* filter)
{
  if(vertsToBurn.empty())
  {
    return;
  }

  std::vector<int64_t> chain;
  int64_t startVert = 0;
  int64_t endVert = 0;
  for(auto vert : vertsToBurn)
  {
    burnedVerts.insert(vert);

    Vertex headVert = {vertPtr[3 * vert + 0], vertPtr[3 * vert + 1], vertPtr[3 * vert + 2]};
    auto hiter = vertexMap.find(headVert);
    if(hiter == vertexMap.end())
    {
      tmpVerts.push_back(headVert);
      tmpNodeTypes.push_back(nodeTypes[vert]);
      vertexMap[headVert] = vertCounter;
      vertCounter++;
    }

    uint16_t numEdgesToVert = edgesContainingVert->getNumberOfElements(vert);
    int64_t* edgesAtVert = edgesContainingVert->getElementListPointer(vert);
    for(uint16_t i = 0; i < numEdgesToVert; i++)
    {
      startVert = vert;
      endVert = (edgePtr[2 * edgesAtVert[i] + 0] != startVert) ? edgePtr[2 * edgesAtVert[i] + 0] : edgePtr[2 * edgesAtVert[i] + 1];
      if((nodeTypes[endVert] == 4 || nodeTypes[endVert] == 14) && burnedVerts.count(endVert) == 0)
      {
        Vertex tailVert = {vertPtr[3 * endVert + 0], vertPtr[3 * endVert + 1], vertPtr[3 * endVert + 2]};
        auto titer = vertexMap.find(tailVert);
        if(titer == vertexMap.end())
        {
          tmpVerts.push_back(tailVert);
          tmpNodeTypes.push_back(nodeTypes[endVert]);
          vertexMap[tailVert] = vertCounter;
          vertCounter++;
        }

        int64_t head = vertexMap[headVert];
        int64_t tail = vertexMap[tailVert];
        if(head > tail)
        {
          std::swap(head, tail);
        }

        Edge tmpEdge = {head, tail};
        auto eiter = edgeMap.find(tmpEdge);
        if(eiter == edgeMap.end())
        {
          tmpEdges.push_back(tmpEdge);
          edgeMap[tmpEdge] = edgeCounter;
          edgeCounter++;
        }
        return;
      }
    }

    for(uint16_t i = 0; i < numEdgesToVert; i++)
    {
      startVert = vert;
      endVert = (edgePtr[2 * edgesAtVert[i] + 0] != startVert) ? edgePtr[2 * edgesAtVert[i] + 0] : edgePtr[2 * edgesAtVert[i] + 1];
      if(startVert != endVert)
      {
        if(burnedVerts.count(endVert) > 0)
        {
          continue;
        }
        else
        {
          Vertex tailVert = {vertPtr[3 * endVert + 0], vertPtr[3 * endVert + 1], vertPtr[3 * endVert + 2]};
          auto titer = vertexMap.find(tailVert);
          if(titer == vertexMap.end())
          {
            tmpVerts.push_back(tailVert);
            tmpNodeTypes.push_back(nodeTypes[endVert]);
            vertexMap[tailVert] = vertCounter;
            vertCounter++;
          }

          int64_t head = vertexMap[headVert];
          int64_t tail = vertexMap[tailVert];
          if(head > tail)
          {
            std::swap(head, tail);
          }

          Edge tmpEdge = {head, tail};
          auto eiter = edgeMap.find(tmpEdge);
          if(eiter == edgeMap.end())
          {
            tmpEdges.push_back(tmpEdge);
            edgeMap[tmpEdge] = edgeCounter;
            edgeCounter++;
          }

          if(!(nodeTypes[endVert] == 4 || nodeTypes[endVert] == 14))
          {
            chain.push_back(endVert);
          }
        }
      }
    }
  }

  QString ss = QObject::tr("Burned %1 Vertices").arg(burnedVerts.size());
  filter->notifyStatusMessage(ss);

  burnVertices(chain, seedVert, burnedVerts, vertCounter, edgeCounter, vertPtr, edgePtr, nodeTypes, edgesContainingVert, vertexMap, edgeMap, tmpVerts, tmpEdges, tmpNodeTypes, filter);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ExtractTripleLinesFromTriangleGeometry::smoothTripleLines()
{
  EdgeGeom::Pointer tripleLineEdge = getDataContainerArray()->getDataContainer(m_EdgeGeometry)->getGeometryAs<EdgeGeom>();
  tripleLineEdge->findElementsContainingVert();
  float* vertPtr = tripleLineEdge->getVertexPointer(0);
  int64_t* edgePtr = tripleLineEdge->getEdgePointer(0);
  ElementDynamicList::Pointer edgesContainingVert = tripleLineEdge->getElementsContainingVert();
  int64_t numVerts = tripleLineEdge->getNumberOfVertices();

  VertexMap vertexMap;
  EdgeMap edgeMap;
  std::vector<Vertex> tmpVerts;
  std::vector<Edge> tmpEdges;
  std::vector<int8_t> tmpNodeTypes;
  std::vector<int32_t> tmpTripleLineIds;
  int64_t vertCounter = 0;
  int64_t edgeCounter = 0;
  // int32_t tripleLineCounter = 1;
  int64_t seedVert = -1;

  // for(int64_t i = 0; i < numVerts; i++)
  //{
  //  if(m_TripleLineNodeTypes[i] == 4 || m_TripleLineNodeTypes[i] == 14)
  //  {
  //    bool valid = true;
  //    uint16_t numEdgesToVert = edgesContainingVert->getNumberOfElements(i);
  //    int64_t* edgesAtVert = edgesContainingVert->getElementListPointer(i);
  //    for(uint16_t j = 0; j < numEdgesToVert; j++)
  //    {
  //      if(((m_TripleLineNodeTypes[edgePtr[2 * edgesAtVert[j] + 0]] == 4 || m_TripleLineNodeTypes[edgePtr[2 * edgesAtVert[j] + 0]] == 14) && edgePtr[2 * edgesAtVert[j] + 0] != i) ||
  //         ((m_TripleLineNodeTypes[edgePtr[2 * edgesAtVert[j] + 1]] == 4 || m_TripleLineNodeTypes[edgePtr[2 * edgesAtVert[j] + 1]] == 14) && edgePtr[2 * edgesAtVert[j] + 1] != i))
  //      {
  //        valid = false;
  //        break;
  //      }
  //    }
  //    if(valid)
  //    {
  //      seedVert = i;
  //      break;
  //    }
  //  }
  //}

  // if(seedVert < 0)
  //{
  //  setErrorCondition(-1, "Unable to find isolated seed quad point");
  //  return;
  //}

  float minExtents[3] = {std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max()};
  float maxExtents[3] = {std::numeric_limits<float>::lowest(), std::numeric_limits<float>::lowest(), std::numeric_limits<float>::lowest()};
  for(int64_t i = 0; i < numVerts; i++)
  {
    if(vertPtr[3 * i + 0] < minExtents[0])
    {
      minExtents[0] = vertPtr[3 * i + 0];
    }
    if(vertPtr[3 * i + 0] > maxExtents[0])
    {
      vertPtr[0] = vertPtr[3 * i + 0];
    }
    if(vertPtr[3 * i + 1] < minExtents[1])
    {
      minExtents[1] = vertPtr[3 * i + 1];
    }
    if(vertPtr[3 * i + 1] > maxExtents[1])
    {
      maxExtents[1] = vertPtr[3 * i + 1];
    }
    if(vertPtr[3 * i + 2] < minExtents[2])
    {
      minExtents[2] = vertPtr[3 * i + 2];
    }
    if(vertPtr[3 * i + 2] > maxExtents[2])
    {
      maxExtents[2] = vertPtr[3 * i + 2];
    }
  }

  Vertex vert = {minExtents[0], minExtents[1], minExtents[2]};
  tmpVerts.push_back(vert);
  tmpNodeTypes.push_back(m_TripleLineNodeTypes[seedVert]);
  vertexMap[vert] = vertCounter;
  vertCounter++;

  uint16_t numEdgesToVert = edgesContainingVert->getNumberOfElements(seedVert);
  int64_t* edgesAtVert = edgesContainingVert->getElementListPointer(seedVert);
  std::unordered_set<int64_t> burnedVerts;

  for(uint16_t i = 0; i < numEdgesToVert; i++)
  {
    burnedVerts.insert(edgePtr[2 * edgesAtVert[i] + 0]);
    burnedVerts.insert(edgePtr[2 * edgesAtVert[i] + 1]);
  }

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dis(0, numEdgesToVert - 1);
  uint16_t idx = dis(gen);
  int64_t begin = (edgePtr[2 * edgesAtVert[idx] + 0] != seedVert) ? edgePtr[2 * edgesAtVert[idx] + 0] : edgePtr[2 * edgesAtVert[idx] + 1];

  Edge tmpEdge = {0, 1};
  tmpEdges.push_back(tmpEdge);
  edgeMap[tmpEdge] = edgeCounter;
  edgeCounter++;

  std::vector<int64_t> vertsToBurn = {begin};

  burnVertices(vertsToBurn, seedVert, burnedVerts, vertCounter, edgeCounter, vertPtr, edgePtr, m_TripleLineNodeTypes, edgesContainingVert, vertexMap, edgeMap, tmpVerts, tmpEdges, tmpNodeTypes, this);

  std::vector<std::vector<float>> controlPoints(tmpVerts.size(), std::vector<float>(3, 0.0f));
  for(auto i = 0; i < controlPoints.size(); i++)
  {
    controlPoints[i][0] = tmpVerts[i][0];
    controlPoints[i][1] = tmpVerts[i][1];
    controlPoints[i][2] = tmpVerts[i][2];
  }

  std::vector<float> smoothPoints = computeBsplineInterpolation(controlPoints, 0, 1);
  VertexGeom::Pointer vertex = VertexGeom::CreateGeometry(smoothPoints.size() / 3, SIMPL::Geometry::VertexGeometry);
  int64_t smoothNumVerts = vertex->getNumberOfVertices();
  float* smoothPtr = vertex->getVertexPointer(0);
  std::memcpy(smoothPtr, smoothPoints.data(), sizeof(float) * smoothNumVerts * 3);
  DataContainer::Pointer dc = DataContainer::New("TEST_SMOOTHING");
  dc->setGeometry(vertex);
  getDataContainerArray()->addOrReplaceDataContainer(dc);

  tripleLineEdge->resizeVertexList(tmpVerts.size());
  tripleLineEdge->resizeEdgeList(tmpEdges.size());
  vertPtr = tripleLineEdge->getVertexPointer(0);
  edgePtr = tripleLineEdge->getEdgePointer(0);
  for(auto i = 0; i < tmpVerts.size(); i++)
  {
    vertPtr[3 * i + 0] = tmpVerts[i][0];
    vertPtr[3 * i + 1] = tmpVerts[i][1];
    vertPtr[3 * i + 2] = tmpVerts[i][2];
  }
  for(auto i = 0; i < tmpEdges.size(); i++)
  {
    edgePtr[2 * i + 0] = tmpEdges[i][0];
    edgePtr[2 * i + 1] = tmpEdges[i][1];
  }

  AttributeMatrix::Pointer attrMat = getDataContainerArray()->getDataContainer(m_EdgeGeometry)->getAttributeMatrix(m_VertexAttributeMatrixName);
  QVector<size_t> tDims(1, tripleLineEdge->getNumberOfVertices());
  attrMat->resizeAttributeArrays(tDims);
  attrMat = getDataContainerArray()->getDataContainer(m_EdgeGeometry)->getAttributeMatrix(m_EdgeAttributeMatrixName);
  tDims[0] = tripleLineEdge->getNumberOfEdges();
  attrMat->resizeAttributeArrays(tDims);

  m_TripleLineNodeTypes = m_TripleLineNodeTypesPtr.lock()->getPointer(0);
  std::memcpy(m_TripleLineNodeTypes, tmpNodeTypes.data(), sizeof(int8_t) * tripleLineEdge->getNumberOfVertices());
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
std::vector<float> ExtractTripleLinesFromTriangleGeometry::computeBsplineInterpolation(std::vector<std::vector<float>> points, float minVal, float maxVal)
{
  std::vector<float> result;
  size_t counter = 1;
  for(float value = minVal; value <= maxVal; value += 0.01f)
  {
    size_t degree = 4;
    auto numKnots = points.size() + degree + 1;
    std::vector<float> knots(numKnots);
    std::vector<float> weights(points.size(), 1.0f);
    std::iota(std::begin(knots), std::end(knots), 0);
    float domainLow = static_cast<float>(degree);
    float domainHigh = static_cast<float>(knots.size() - 1 - degree);
    float low = knots[domainLow];
    float high = knots[domainHigh];
    float t = value * (high - low) + low;
    size_t s;
    for(s = domainLow; s < domainHigh; s++)
    {
      if(t >= knots[s] && t <= knots[s + 1])
      {
        break;
      }
    }

    std::vector<std::vector<float>> v(points);
    for(auto&& point : v)
    {
      point.push_back(1.0f);
    }

    float alpha;
    for(size_t l = 1; l <= degree; l++)
    {
      for(size_t i = s; i > (s - degree - 1 + 1); i--)
      {
        alpha = (t - knots[i]) / (knots[i + degree + 1 - 1] - knots[i]);
        for(size_t j = 0; j < (3 + 1); j++)
        {
          v[i][j] = (1 - alpha) * v[i - 1][j] + alpha * v[i][j];
        }
      }
    }

    for(size_t i = 0; i < 3; i++)
    {
      result.push_back(v[s][i] / v[s][3]);
    }

    QString ss = QObject::tr("Smoothing Triple Line || Iteration %1").arg(counter);
    notifyStatusMessage(ss);
    counter++;
  }

  return result;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ExtractTripleLinesFromTriangleGeometry::execute()
{
  initialize();
  dataCheck();
  if(getErrorCode() < 0)
  {
    return;
  }

  extractTripleLines();

  if(getErrorCode() < 0)
  {
    return;
  }

  if(m_SmoothTripleLines)
  {
    smoothTripleLines();
  }

  notifyStatusMessage("Complete");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AbstractFilter::Pointer ExtractTripleLinesFromTriangleGeometry::newFilterInstance(bool copyFilterParameters) const
{
  ExtractTripleLinesFromTriangleGeometry::Pointer filter = ExtractTripleLinesFromTriangleGeometry::New();
  if(copyFilterParameters)
  {
    copyFilterParameterInstanceVariables(filter.get());
  }
  return filter;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString ExtractTripleLinesFromTriangleGeometry::getCompiledLibraryName() const
{
  return DREAM3DReviewConstants::DREAM3DReviewBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString ExtractTripleLinesFromTriangleGeometry::getBrandingString() const
{
  return "MeshUtilities";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString ExtractTripleLinesFromTriangleGeometry::getFilterVersion() const
{
  QString version;
  QTextStream vStream(&version);
  vStream << DREAM3DReview::Version::Major() << "." << DREAM3DReview::Version::Minor() << "." << DREAM3DReview::Version::Patch();
  return version;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString ExtractTripleLinesFromTriangleGeometry::getGroupName() const
{
  return SIMPL::FilterGroups::SamplingFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString ExtractTripleLinesFromTriangleGeometry::getSubGroupName() const
{
  return SIMPL::FilterSubGroups::GeometryFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString ExtractTripleLinesFromTriangleGeometry::getHumanLabel() const
{
  return "Extract Triple Lines from Triangle Geometry";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QUuid ExtractTripleLinesFromTriangleGeometry::getUuid()
{
  return QUuid("{1422a4d7-ecd3-5e84-8883-8d2c8c551675}");
}
