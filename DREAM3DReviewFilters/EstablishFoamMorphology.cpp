/*
 * This material is based upon work supported by the National Science Foundation
 * under Grant No. 1629660.  Any opinions, findings, and conclusions or
 * recommendations expressed in this material are those of the author(s) and do not
 * necessarily reflect the views of the National Science Foundation.
 *
 *
 * The code is a derivative work based on the filter PackPrimaryPhases,
 * InsertPrecipitatePhases and EuclideanDistanceMap which were developed under DoD
 * contract:   United States Air Force Prime Contract FA8650-07-D-5800. See those files for their respective copyright
 * and distribution license.
 */

#include "EstablishFoamMorphology.h"

#include <fstream>

#include <QtCore/QDir>
#include <QtCore/QFile>
#include <QtCore/QFileInfo>

#include "SIMPLib/SIMPLib.h"
#include "SIMPLib/Common/Constants.h"
#include "SIMPLib/Common/ShapeType.h"
#include "SIMPLib/DataArrays/IDataArray.h"
#include "SIMPLib/DataArrays/NeighborList.hpp"
#include "SIMPLib/DataContainers/DataContainer.h"
#include "SIMPLib/DataContainers/DataContainerArray.h"
#include "SIMPLib/FilterParameters/BooleanFilterParameter.h"
#include "SIMPLib/FilterParameters/DataArraySelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/DataContainerSelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/DoubleFilterParameter.h"
#include "SIMPLib/FilterParameters/FloatFilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedBooleanFilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedChoicesFilterParameter.h"
#include "SIMPLib/FilterParameters/OutputFileFilterParameter.h"
#include "SIMPLib/FilterParameters/SeparatorFilterParameter.h"
#include "SIMPLib/FilterParameters/StringFilterParameter.h"
#include "SIMPLib/Geometry/ImageGeom.h"
#include "SIMPLib/Math/MatrixMath.h"
#include "SIMPLib/Math/SIMPLibMath.h"
#include "SIMPLib/Math/SIMPLibRandom.h"
#include "SIMPLib/StatsData/PrecipitateStatsData.h"
#include "SIMPLib/Utilities/TimeUtilities.h"

#include "EbsdLib/Core/Orientation.hpp"
#include "EbsdLib/Core/OrientationTransformation.hpp"
#include "EbsdLib/Core/Quaternion.hpp"
#include "EbsdLib/LaueOps/OrthoRhombicOps.h"

#include "DREAM3DReview/DREAM3DReviewConstants.h"
#include "DREAM3DReview/DREAM3DReviewVersion.h"

// SIMPLib.h MUST be included before this or the guard will block the include but not its uses below.
// This is consistent with previous behavior, only earlier parallelization split the includes between
// the corresponding .h and .cpp files.
#ifdef SIMPL_USE_PARALLEL_ALGORITHMS
// clang-format off
#include <tbb/blocked_range.h>
#include <tbb/blocked_range3d.h>
#include <tbb/parallel_for.h>
#include <tbb/partitioner.h>
#include <tbb/task_group.h>
// clang-format on
#endif

/**
 * @brief The FindEuclideanMap2 class implements a threaded algorithm that computes the Euclidean distance map
 * for each point in the supplied volume
 */
class FindEuclideanMap2
{
  DataContainer::Pointer m_DataContainer;
  int32_t* m_FeatureIds;
  int32_t* m_NearestNeighbors;
  float* m_GBEuclideanDistances;
  float* m_TJEuclideanDistances;
  float* m_QPEuclideanDistances;
  int32_t mapType;

public:
  FindEuclideanMap2(const DataContainer::Pointer& datacontainer, int32_t* fIds, int32_t* nearNeighs, float* gbDists, float* tjDists, float* qpDists, int32_t type)
  : m_DataContainer(datacontainer)
  , m_FeatureIds(fIds)
  , m_NearestNeighbors(nearNeighs)
  , m_GBEuclideanDistances(gbDists)
  , m_TJEuclideanDistances(tjDists)
  , m_QPEuclideanDistances(qpDists)
  , mapType(type)
  {
  }

  ~FindEuclideanMap2() = default;

  void operator()() const
  {
    size_t totalPoints = m_DataContainer->getGeometryAs<ImageGeom>()->getNumberOfElements();
    double euclideanDistance = 0.0f;
    size_t count = 1;
    size_t changed = 1;
    size_t neighpoint = 0;
    int64_t nearestneighbor;
    std::array<int64_t, 6> neighbors = {0, 0, 0, 0, 0, 0};
    ImageGeom::Pointer imageGeom = m_DataContainer->getGeometryAs<ImageGeom>();
    int64_t xpoints = static_cast<int64_t>(imageGeom->getXPoints());
    int64_t ypoints = static_cast<int64_t>(imageGeom->getYPoints());
    int64_t zpoints = static_cast<int64_t>(imageGeom->getZPoints());

    FloatVec3Type spacing = imageGeom->getSpacing();

    double resx = static_cast<double>(spacing[0]);
    double resy = static_cast<double>(spacing[1]);
    double resz = static_cast<double>(spacing[2]);

    neighbors[0] = -xpoints * ypoints;
    neighbors[1] = -xpoints;
    neighbors[2] = -1;
    neighbors[3] = 1;
    neighbors[4] = xpoints;
    neighbors[5] = xpoints * ypoints;

    // Use a std::vector to get an auto cleaned up array thus not needing the 'delete' keyword later on.
    std::vector<int64_t> voxNN(totalPoints, 0);
    int64_t* voxel_NearestNeighbor = &(voxNN.front());
    std::vector<double> voxEDist(totalPoints, 0.0);
    double* voxel_EuclideanDistance = &(voxEDist.front());

    euclideanDistance = 0;
    for(size_t a = 0; a < totalPoints; ++a)
    {
      if(m_NearestNeighbors[a * 3 + mapType] >= 0)
      {
        voxel_NearestNeighbor[a] = a;
      } // if voxel is boundary voxel, then want to use itself as nearest boundary voxel
      else
      {
        voxel_NearestNeighbor[a] = -1;
      }
      if(mapType == 0)
      {
        voxel_EuclideanDistance[a] = static_cast<double>(m_GBEuclideanDistances[a]);
      }
      else if(mapType == 1)
      {
        voxel_EuclideanDistance[a] = static_cast<double>(m_TJEuclideanDistances[a]);
      }
      else if(mapType == 2)
      {
        voxel_EuclideanDistance[a] = static_cast<double>(m_QPEuclideanDistances[a]);
      }
    }
    count = 1;
    changed = 1;
    int64_t i = 0;
    int64_t zBlock = xpoints * ypoints;
    int64_t zStride = 0;
    int64_t yStride = 0;
    std::array<char, 6> mask = {0, 0, 0, 0, 0, 0};
    while(count > 0 && changed > 0)
    {
      count = 0;
      changed = 0;
      euclideanDistance++;

      for(int64_t z = 0; z < zpoints; ++z)
      {
        zStride = z * zBlock;
        mask[0] = mask[5] = 1;
        if(z == 0)
        {
          mask[0] = 0;
        }
        if(z == zpoints - 1)
        {
          mask[5] = 0;
        }

        for(int64_t y = 0; y < ypoints; ++y)
        {
          yStride = y * xpoints;
          mask[1] = mask[4] = 1;
          if(y == 0)
          {
            mask[1] = 0;
          }
          if(y == ypoints - 1)
          {
            mask[4] = 0;
          }

          for(int64_t x = 0; x < xpoints; ++x)
          {
            mask[2] = mask[3] = 1;
            if(x == 0)
            {
              mask[2] = 0;
            }
            if(x == xpoints - 1)
            {
              mask[3] = 0;
            }

            i = zStride + yStride + x;
            if(voxel_NearestNeighbor[i] == -1 && m_FeatureIds[i] > 0)
            {
              count++;
              for(int32_t j = 0; j < 6; j++)
              {
                neighpoint = i + neighbors[j];
                if(mask[j] == 1)
                {
                  if(voxel_EuclideanDistance[neighpoint] != -1.0)
                  {
                    voxel_NearestNeighbor[i] = voxel_NearestNeighbor[neighpoint];
                  }
                }
              }
            }
          }
        }
      }
      for(size_t j = 0; j < totalPoints; ++j)
      {
        if(voxel_NearestNeighbor[j] != -1 && voxel_EuclideanDistance[j] == -1.0 && m_FeatureIds[j] > 0)
        {
          changed++;
          voxel_EuclideanDistance[j] = euclideanDistance;
        }
      }
    }
    double x1 = 0.0;
    double x2 = 0.0;
    double y1 = 0.0;
    double y2 = 0.0;
    double z1 = 0.0;
    double z2 = 0.0;
    double dist = 0.0;
    double oneOverzBlock = 1.0 / double(zBlock);
    double oneOverxpoints = 1.0 / double(xpoints);
    for(int64_t m = 0; m < zpoints; m++)
    {
      zStride = m * zBlock;
      for(int64_t n = 0; n < ypoints; n++)
      {
        yStride = n * xpoints;
        for(int64_t p = 0; p < xpoints; p++)
        {
          x1 = double(p) * resx;
          y1 = double(n) * resy;
          z1 = double(m) * resz;
          nearestneighbor = voxel_NearestNeighbor[zStride + yStride + p];
          if(nearestneighbor >= 0)
          {
            x2 = resx * double(nearestneighbor % xpoints);                           // find_xcoord(nearestneighbor);
            y2 = resy * double(int64_t(nearestneighbor * oneOverxpoints) % ypoints); // find_ycoord(nearestneighbor);
            z2 = resz * floor(nearestneighbor * oneOverzBlock);                      // find_zcoord(nearestneighbor);
            dist = ((x1 - x2) * (x1 - x2)) + ((y1 - y2) * (y1 - y2)) + ((z1 - z2) * (z1 - z2));
            dist = sqrt(dist);
            voxel_EuclideanDistance[zStride + yStride + p] = dist;
          }
        }
      }
    }
    for(size_t a = 0; a < totalPoints; ++a)
    {
      m_NearestNeighbors[a * 3 + mapType] = voxel_NearestNeighbor[a];
      if(mapType == 0)
      {
        m_GBEuclideanDistances[a] = static_cast<float>(voxel_EuclideanDistance[a]);
      }
      else if(mapType == 1)
      {
        m_TJEuclideanDistances[a] = static_cast<float>(voxel_EuclideanDistance[a]);
      }
      else if(mapType == 2)
      {
        m_QPEuclideanDistances[a] = static_cast<float>(voxel_EuclideanDistance[a]);
      }
    }
  }
};

/**
 * @brief The FoamAssignVoxelsGapsImpl class implements a threaded algorithm that assigns all the voxels
 * in the volume to a unique Feature.
 */
class FoamAssignVoxelsGapsImpl
{
  std::array<size_t, 3> m_UDims;
  std::array<float, 3> m_Invradcur;
  std::array<float, 3> m_Spacing;

  float ga[3][3];
  int32_t curFeature = 0;

  float xc = 0.0f;
  float yc = 0.0f;
  float zc = 0.0f;
  // int32_t* m_FeatureIds = nullptr;
  ShapeOps* m_ShapeOps = nullptr;

  Int32ArrayType::Pointer m_NewOwnersPtr;
  FloatArrayType::Pointer mEllipFuncsPtr;

public:
  FoamAssignVoxelsGapsImpl(SizeVec3Type& dimensions, FloatVec3Type& resolution, int32_t* featureIds, const float* radCur, const float* xx, ShapeOps* shapeOps, float gA[3][3], const float* size,
                           int32_t cur_feature, const Int32ArrayType::Pointer& newowners, const FloatArrayType::Pointer& ellipfuncs)
  : curFeature(cur_feature)
  , m_ShapeOps(shapeOps)
  {
    m_UDims[0] = dimensions[0];
    m_UDims[1] = dimensions[1];
    m_UDims[2] = dimensions[2];

    m_Invradcur[0] = 1.0 / radCur[0];
    m_Invradcur[1] = 1.0 / radCur[1];
    m_Invradcur[2] = 1.0 / radCur[2];

    m_Spacing[0] = resolution[0];
    m_Spacing[1] = resolution[1];
    m_Spacing[2] = resolution[2];

    xc = xx[0];
    yc = xx[1];
    zc = xx[2];

    ga[0][0] = gA[0][0];
    ga[0][1] = gA[0][1];
    ga[0][2] = gA[0][2];
    ga[1][0] = gA[1][0];
    ga[1][1] = gA[1][1];
    ga[1][2] = gA[1][2];
    ga[2][0] = gA[2][0];
    ga[2][1] = gA[2][1];
    ga[2][2] = gA[2][2];

    m_NewOwnersPtr = newowners;
    mEllipFuncsPtr = ellipfuncs;
  }

  ~FoamAssignVoxelsGapsImpl() = default;

  // -----------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------
  void convert(int64_t zStart, int64_t zEnd, int64_t yStart, int64_t yEnd, int64_t xStart, int64_t xEnd) const
  {
    int64_t column = 0;
    int64_t row = 0;
    int64_t plane = 0;
    int64_t index = 0;
    std::array<float, 3> coords = {0.0f, 0.0f, 0.0f};
    float inside = 0.0f;
    std::array<float, 3> coordsRotated = {0.0f, 0.0f, 0.0f};
    int32_t* newowners = m_NewOwnersPtr->getPointer(0);
    float* ellipfuncs = mEllipFuncsPtr->getPointer(0);

    const std::array<int64_t, 3> dims = {
        static_cast<int64_t>(m_UDims[0]),
        static_cast<int64_t>(m_UDims[1]),
        static_cast<int64_t>(m_UDims[2]),
    };

    int64_t dim0_dim_1 = m_UDims[0] * m_UDims[1];
    for(int64_t iter1 = xStart; iter1 < xEnd; iter1++)
    {
      column = iter1;
      if(iter1 < 0)
      {
        column = iter1 + dims[0];
      }
      else if(iter1 > dims[0] - 1)
      {
        column = iter1 - dims[0];
      }

      for(int64_t iter2 = yStart; iter2 < yEnd; iter2++)
      {
        row = iter2;
        if(iter2 < 0)
        {
          row = iter2 + dims[1];
        }
        else if(iter2 > dims[1] - 1)
        {
          row = iter2 - dims[1];
        }
        size_t row_dim = row * dims[0];

        for(int64_t iter3 = zStart; iter3 < zEnd; iter3++)
        {
          plane = iter3;
          if(iter3 < 0)
          {
            plane = iter3 + dims[2];
          }
          else if(iter3 > dims[2] - 1)
          {
            plane = iter3 - dims[2];
          }

          index = static_cast<int64_t>((plane * dim0_dim_1) + (row_dim) + column);

          inside = -1.0f;
          coords[0] = float(iter1) * m_Spacing[0];
          coords[1] = float(iter2) * m_Spacing[1];
          coords[2] = float(iter3) * m_Spacing[2];

          coords[0] = coords[0] - xc;
          coords[1] = coords[1] - yc;
          coords[2] = coords[2] - zc;
          MatrixMath::Multiply3x3with3x1(ga, coords.data(), coordsRotated.data());
          float axis1comp = coordsRotated[0] * m_Invradcur[0];
          float axis2comp = coordsRotated[1] * m_Invradcur[1];
          float axis3comp = coordsRotated[2] * m_Invradcur[2];
          inside = m_ShapeOps->inside(axis1comp, axis2comp, axis3comp);

          if(inside >= 0 && newowners[index] > 0 && inside > ellipfuncs[index])
          {
            newowners[index] = curFeature;
            ellipfuncs[index] = inside;
          }
          else if(inside >= 0 && newowners[index] == -1)
          {
            newowners[index] = curFeature;
            ellipfuncs[index] = inside;
          }
        }
      }
    }
  }

#ifdef SIMPL_USE_PARALLEL_ALGORITHMS
  void operator()(const tbb::blocked_range3d<int64_t, int64_t, int64_t>& r) const
  {
    convert(r.pages().begin(), r.pages().end(), r.rows().begin(), r.rows().end(), r.cols().begin(), r.cols().end());
  }
#endif
};

// Include the MOC generated file for this class
#include "moc_EstablishFoamMorphology.cpp"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
EstablishFoamMorphology::EstablishFoamMorphology()
{
  initialize();
  setupFilterParameters();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
EstablishFoamMorphology::~EstablishFoamMorphology() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void EstablishFoamMorphology::initialize()
{
  m_Neighbors = nullptr;
  m_BoundaryCells = nullptr;

  m_StatsDataArray = StatsDataArray::NullPointer();

  m_ShapeOps = ShapeOps::getShapeOpsQVector();
  m_UnknownShapeOps = ShapeOps::NullPointer();
  m_CubicOctohedronOps = ShapeOps::NullPointer();
  m_CylinderOps = ShapeOps::NullPointer();
  m_EllipsoidOps = ShapeOps::NullPointer();
  m_SuperEllipsoidOps = ShapeOps::NullPointer();
  m_OrthoOps = OrthoRhombicOps::New();

  m_ColumnList.clear();
  m_RowList.clear();
  m_PlaneList.clear();
  m_EllipFuncList.clear();

  m_PointsToAdd.clear();
  m_PointsToRemove.clear();
  m_Seed = QDateTime::currentMSecsSinceEpoch();
  m_FirstFoamFeature = 1;
  m_SizeX = m_SizeY = m_SizeZ = m_TotalVol = 0.0f;
  m_TotalVol = 1.0f;
  m_HalfPackingRes[0] = m_HalfPackingRes[1] = m_HalfPackingRes[2] = 1.0f;
  m_OneOverPackingRes[0] = m_OneOverPackingRes[1] = m_OneOverPackingRes[2] = 1.0f;
  m_OneOverHalfPackingRes[0] = m_OneOverHalfPackingRes[1] = m_OneOverHalfPackingRes[2] = 1.0f;
  m_PackingRes[0] = m_PackingRes[1] = m_PackingRes[2] = 1.0f;

  m_PackingPoints[0] = m_PackingPoints[1] = m_PackingPoints[2] = 1;
  m_TotalPackingPoints = 1;

  m_FeatureSizeDist.clear();
  m_SimFeatureSizeDist.clear();
  m_NeighborDist.clear();
  m_SimNeighborDist.clear();

  m_FeatureSizeDistStep.clear();
  m_NeighborDistStep.clear();

  m_PackQualities.clear();
  m_GSizes.clear();

  m_PrecipitatePhases.clear();
  m_PrecipitatePhaseFractions.clear();

  m_AvailablePointsCount = 1;
  m_FillingError = m_OldFillingError = 0.0f;
  m_CurrentNeighborhoodError = m_OldNeighborhoodError = 0.0f;
  m_CurrentSizeDistError = m_OldSizeDistError = 0.0f;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void EstablishFoamMorphology::setupFilterParameters()
{
  FilterParameterVectorType parameters;
  parameters.push_back(SIMPL_NEW_BOOL_FP("Periodic Boundaries", PeriodicBoundaries, FilterParameter::Category::Parameter, EstablishFoamMorphology));

  parameters.push_back(SeparatorFilterParameter::Create("Data Container", FilterParameter::Category::RequiredArray));
  {
    DataContainerSelectionFilterParameter::RequirementType req;
    parameters.push_back(SIMPL_NEW_DC_SELECTION_FP("Data Container", DataContainerName, FilterParameter::Category::RequiredArray, EstablishFoamMorphology, req));
  }

  parameters.push_back(SeparatorFilterParameter::Create("Cell Ensemble Data", FilterParameter::Category::RequiredArray));
  {
    DataArraySelectionFilterParameter::RequirementType req =
        DataArraySelectionFilterParameter::CreateRequirement(SIMPL::TypeNames::StatsDataArray, 1, AttributeMatrix::Type::CellEnsemble, IGeometry::Type::Any);
    IGeometry::Types geomTypes;
    geomTypes.push_back(IGeometry::Type::Image);
    geomTypes.push_back(IGeometry::Type::Unknown);
    req.dcGeometryTypes = geomTypes;
    parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Statistics", InputStatsArrayPath, FilterParameter::Category::RequiredArray, EstablishFoamMorphology, req));
  }
  {
    DataArraySelectionFilterParameter::RequirementType req =
        DataArraySelectionFilterParameter::CreateRequirement(SIMPL::TypeNames::UInt32, 1, AttributeMatrix::Type::CellEnsemble, IGeometry::Type::Any);
    IGeometry::Types geomTypes;
    req.dcGeometryTypes = geomTypes;
    parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Phase Types", InputPhaseTypesArrayPath, FilterParameter::Category::RequiredArray, EstablishFoamMorphology, req));
  }
  {
    DataArraySelectionFilterParameter::RequirementType req =
        DataArraySelectionFilterParameter::CreateRequirement(SIMPL::TypeNames::StringArray, 1, AttributeMatrix::Type::CellEnsemble, IGeometry::Type::Any);
    IGeometry::Types geomTypes;
    req.dcGeometryTypes = geomTypes;
    parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Phase Names", InputPhaseNamesArrayPath, FilterParameter::Category::RequiredArray, EstablishFoamMorphology, req));
  }
  {
    DataArraySelectionFilterParameter::RequirementType req =
        DataArraySelectionFilterParameter::CreateRequirement(SIMPL::TypeNames::UInt32, 1, AttributeMatrix::Type::CellEnsemble, IGeometry::Type::Any);
    IGeometry::Types geomTypes;
    geomTypes.push_back(IGeometry::Type::Image);
    geomTypes.push_back(IGeometry::Type::Unknown);
    req.dcGeometryTypes = geomTypes;
    parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Shape Types", InputShapeTypesArrayPath, FilterParameter::Category::RequiredArray, EstablishFoamMorphology, req));
  }

  parameters.push_back(SeparatorFilterParameter::Create("Cell Data", FilterParameter::Category::RequiredArray));
  {
    DataArraySelectionFilterParameter::RequirementType req = DataArraySelectionFilterParameter::CreateRequirement(SIMPL::TypeNames::Int32, 1, AttributeMatrix::Type::Cell, IGeometry::Type::Image);
    parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Feature Ids", InputCellFeatureIdsArrayPath, FilterParameter::Category::RequiredArray, EstablishFoamMorphology, req, {1}));
  }

  parameters.push_back(SeparatorFilterParameter::Create("Cell Ensemble Data", FilterParameter::Category::CreatedArray));
  parameters.push_back(SIMPL_NEW_STRING_FP("Cell Ensemble Attribute Matrix", OutputCellEnsembleAttributeMatrixName, FilterParameter::Category::CreatedArray, EstablishFoamMorphology));
  parameters.push_back(SIMPL_NEW_STRING_FP("Number of Features", NumFeaturesArrayName, FilterParameter::Category::CreatedArray, EstablishFoamMorphology));

  parameters.push_back(SeparatorFilterParameter::Create("Cell Feature Data", FilterParameter::Category::CreatedArray));
  parameters.push_back(SIMPL_NEW_STRING_FP("Cell Feature Attribute Matrix", OutputCellFeatureAttributeMatrixName, FilterParameter::Category::CreatedArray, EstablishFoamMorphology));

  parameters.push_back(SeparatorFilterParameter::Create("Cell Data", FilterParameter::Category::CreatedArray));
  parameters.push_back(SIMPL_NEW_STRING_FP("Feature Ids", FeatureIdsArrayName, FilterParameter::Category::CreatedArray, EstablishFoamMorphology, {0}));
  parameters.push_back(SIMPL_NEW_STRING_FP("Mask", MaskArrayName, FilterParameter::Category::CreatedArray, EstablishFoamMorphology));
  parameters.push_back(SIMPL_NEW_STRING_FP("Phases", CellPhasesArrayName, FilterParameter::Category::CreatedArray, EstablishFoamMorphology));
  parameters.push_back(SIMPL_NEW_STRING_FP("Quadruple Point Euclidean Distances", QPEuclideanDistancesArrayName, FilterParameter::Category::CreatedArray, EstablishFoamMorphology));
  parameters.push_back(SIMPL_NEW_STRING_FP("Triple Junction Euclidean Distances", TJEuclideanDistancesArrayName, FilterParameter::Category::CreatedArray, EstablishFoamMorphology));
  parameters.push_back(SIMPL_NEW_STRING_FP("Grain Boundary Euclidean Distances", GBEuclideanDistancesArrayName, FilterParameter::Category::CreatedArray, EstablishFoamMorphology));

  std::vector<QString> linkedProps = {"CsvOutputFile"};
  linkedProps.push_back("CsvOutputFile");
  parameters.push_back(SIMPL_NEW_LINKED_BOOL_FP("Write Goal Attributes", WriteGoalAttributes, FilterParameter::Category::Parameter, EstablishFoamMorphology, linkedProps));
  parameters.push_back(SIMPL_NEW_OUTPUT_FILE_FP("Goal Attribute CSV File", CsvOutputFile, FilterParameter::Category::Parameter, EstablishFoamMorphology, "*.csv", "Comma Separated Data"));
  linkedProps.clear();
  linkedProps.push_back("FeatureIdsArrayName");
  linkedProps.push_back("InputCellFeatureIdsArrayPath");
  {
    LinkedChoicesFilterParameter::Pointer parameter = LinkedChoicesFilterParameter::New();
    parameter->setHumanLabel("Already Have Features");
    parameter->setPropertyName("HaveFeatures");
    parameter->setSetterCallback(SIMPL_BIND_SETTER(EstablishFoamMorphology, this, HaveFeatures));
    parameter->setGetterCallback(SIMPL_BIND_GETTER(EstablishFoamMorphology, this, HaveFeatures));

    parameter->setDefaultValue(getHaveFeatures()); // Just set the first index
    std::vector<QString> choices;
    choices.push_back("No");
    choices.push_back("Yes");
    parameter->setChoices(choices);
    parameter->setLinkedProperties(linkedProps);
    parameter->setEditable(false);
    parameter->setCategory(FilterParameter::Category::Parameter);
    parameters.push_back(parameter);
  }
  parameters.push_back(SIMPL_NEW_DOUBLE_FP("Minimum Strut Thickness", MinStrutThickness, FilterParameter::Category::Parameter, EstablishFoamMorphology));
  parameters.push_back(SIMPL_NEW_FLOAT_FP("Strut Thickness Variability Factor", StrutThicknessVariability, FilterParameter::Category::Parameter, EstablishFoamMorphology));
  parameters.push_back(SIMPL_NEW_FLOAT_FP("Strut Cross Section Shape Factor", StrutShapeVariability, FilterParameter::Category::Parameter, EstablishFoamMorphology));

  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
// void EstablishFoamMorphology::readFilterParameters(AbstractFilterParametersReader* reader, int index)
//{
//	reader->openFilterGroup(this, index);
//	setOutputCellAttributeMatrixPath(reader->readDataArrayPath("OutputCellAttributeMatrixPath", getOutputCellAttributeMatrixPath()));
//	setInputCellFeatureIdsArrayPath(reader->readDataArrayPath("InputCellFeatureIdsArrayPath", getInputCellFeatureIdsArrayPath()));
//	setOutputCellEnsembleAttributeMatrixName(reader->readString("OutputCellEnsembleAttributeMatrixName", getOutputCellEnsembleAttributeMatrixName()));
//	setOutputCellFeatureAttributeMatrixName(reader->readString("OutputCellFeatureAttributeMatrixName", getOutputCellFeatureAttributeMatrixName()));
//	setFeatureIdsArrayName(reader->readString("FeatureIdsArrayName", getFeatureIdsArrayName()));
//	setMaskArrayName(reader->readString("MaskArrayName", getMaskArrayName()));
//	setCellPhasesArrayName(reader->readString("CellPhasesArrayName", getCellPhasesArrayName()));
//	setNumFeaturesArrayName(reader->readString("NumFeaturesArrayName", getNumFeaturesArrayName()));
//	setQPEuclideanDistancesArrayName(reader->readString("QPEuclideanDistancesArrayName", getQPEuclideanDistancesArrayName()));
//	setTJEuclideanDistancesArrayName(reader->readString("TJEuclideanDistancesArrayName", getTJEuclideanDistancesArrayName()));
//	setGBEuclideanDistancesArrayName(reader->readString("GBEuclideanDistancesArrayName", getGBEuclideanDistancesArrayName()));
//	setPeriodicBoundaries(reader->readValue("PeriodicBoundaries", false));
//	setWriteGoalAttributes(reader->readValue("WriteGoalAttributes", false));
//	setCsvOutputFile(reader->readString("CsvOutputFile", getCsvOutputFile()));
//	setHaveFeatures(reader->readValue("HaveFeatures", getHaveFeatures()));
//	setInputStatsArrayPath(reader->readDataArrayPath("InputStatsArrayPath", getInputStatsArrayPath()));
//	setInputPhaseTypesArrayPath(reader->readDataArrayPath("InputPhaseTypesArrayPath", getInputPhaseTypesArrayPath()));
//	setInputShapeTypesArrayPath(reader->readDataArrayPath("InputShapeTypesArrayPath", getInputShapeTypesArrayPath()));
//	reader->closeFilterGroup();
//}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void EstablishFoamMorphology::updateFeatureInstancePointers()
{
  clearErrorCode();
  clearWarningCode();
  if(nullptr != m_FeaturePhasesPtr.lock())
  {
    m_FeaturePhases = m_FeaturePhasesPtr.lock()->getPointer(0);
  }
  if(nullptr != m_NeighborhoodsPtr.lock())
  {
    m_Neighborhoods = m_NeighborhoodsPtr.lock()->getPointer(0);
  }
  if(nullptr != m_EquivalentDiametersPtr.lock())
  {
    m_EquivalentDiameters = m_EquivalentDiametersPtr.lock()->getPointer(0);
  }
  if(nullptr != m_VolumesPtr.lock())
  {
    m_Volumes = m_VolumesPtr.lock()->getPointer(0);
  }
  if(nullptr != m_Omega3sPtr.lock())
  {
    m_Omega3s = m_Omega3sPtr.lock()->getPointer(0);
  }
  if(nullptr != m_CentroidsPtr.lock())
  {
    m_Centroids = m_CentroidsPtr.lock()->getPointer(0);
  }
  if(nullptr != m_AxisEulerAnglesPtr.lock())
  {
    m_AxisEulerAngles = m_AxisEulerAnglesPtr.lock()->getPointer(0);
  }
  if(nullptr != m_AxisLengthsPtr.lock())
  {
    m_AxisLengths = m_AxisLengthsPtr.lock()->getPointer(0);
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void EstablishFoamMorphology::dataCheck()
{
  clearErrorCode();
  clearWarningCode();
  DataArrayPath tempPath;

  // Make sure we have our input DataContainer with the proper Ensemble data
  getDataContainerArray()->getPrereqGeometryFromDataContainer<ImageGeom>(this, getOutputCellAttributeMatrixPath().getDataContainerName());

  QVector<DataArrayPath> cellDataArrayPaths;
  QVector<DataArrayPath> ensembleDataArrayPaths;

  // Input Ensemble Data that we require
  std::vector<size_t> cDims(1, 1);
  m_PhaseTypesPtr = getDataContainerArray()->getPrereqArrayFromPath<DataArray<uint32_t>>(this, getInputPhaseTypesArrayPath(), cDims);
  if(nullptr != m_PhaseTypesPtr.lock())
  {
    m_PhaseTypes = m_PhaseTypesPtr.lock()->getPointer(0);
  }
  if(getErrorCode() >= 0)
  {
    ensembleDataArrayPaths.push_back(getInputPhaseTypesArrayPath());
  }

  m_PhaseNamesPtr = getDataContainerArray()->getPrereqArrayFromPath<StringDataArray>(this, getInputPhaseNamesArrayPath(), cDims);
  if(getErrorCode() >= 0)
  {
    ensembleDataArrayPaths.push_back(getInputPhaseNamesArrayPath());
  }

  m_ShapeTypesPtr = getDataContainerArray()->getPrereqArrayFromPath<DataArray<uint32_t>>(this, getInputShapeTypesArrayPath(), cDims);
  if(nullptr != m_ShapeTypesPtr.lock())
  {
    m_ShapeTypes = m_ShapeTypesPtr.lock()->getPointer(0);
  }
  if(getErrorCode() >= 0)
  {
    ensembleDataArrayPaths.push_back(getInputShapeTypesArrayPath());
  }

  m_StatsDataArray = getDataContainerArray()->getPrereqArrayFromPath<StatsDataArray>(this, getInputStatsArrayPath(), cDims);
  if(m_StatsDataArray.lock() == nullptr)
  {
    QString ss = QObject::tr("Statistics array is not initialized correctly. The path is %1").arg(getInputStatsArrayPath().serialize());
    setErrorCondition(-308, ss);
  }
  if(getErrorCode() >= 0)
  {
    ensembleDataArrayPaths.push_back(getInputStatsArrayPath());
  }

  cDims[0] = 1;
  // Cell Data
  if(m_HaveFeatures == 1)
  {
    m_FeatureIdsPtr = getDataContainerArray()->getPrereqArrayFromPath<DataArray<int32_t>>(this, getInputCellFeatureIdsArrayPath(), cDims);
    if(nullptr != m_FeatureIdsPtr.lock())
    {
      m_FeatureIds = m_FeatureIdsPtr.lock()->getPointer(0);
    }
  }

  if(m_HaveFeatures == 0)
  {
    tempPath.update(getOutputCellAttributeMatrixPath().getDataContainerName(), getOutputCellAttributeMatrixPath().getAttributeMatrixName(), getFeatureIdsArrayName());
    m_FeatureIdsPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<int32_t>>(
        this, tempPath, -1, cDims); /* Assigns the shared_ptr<>(this, tempPath, -1, udims); Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
    if(nullptr != m_FeatureIdsPtr.lock())
    {
      m_FeatureIds = m_FeatureIdsPtr.lock()->getPointer(0);
    }
  }

  tempPath.update(getOutputCellAttributeMatrixPath().getDataContainerName(), getOutputCellAttributeMatrixPath().getAttributeMatrixName(), getMaskArrayName());
  m_MaskPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<bool>>(
      this, tempPath, false, cDims); /* Assigns the shared_ptr<>(this, tempPath, -1, udims); Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
  if(nullptr != m_MaskPtr.lock())
  {
    m_Mask = m_MaskPtr.lock()->getPointer(0);
  }

  tempPath.update(getOutputCellAttributeMatrixPath().getDataContainerName(), getOutputCellAttributeMatrixPath().getAttributeMatrixName(), getCellPhasesArrayName());
  m_CellPhasesPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<int32_t>>(
      this, tempPath, 0, cDims); /* Assigns the shared_ptr<>(this, tempPath, 0, udims); Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
  if(nullptr != m_CellPhasesPtr.lock())
  {
    m_CellPhases = m_CellPhasesPtr.lock()->getPointer(0);
  }

  if(getErrorCode() < 0)
  {
    return;
  }

  DataContainer::Pointer m = getDataContainerArray()->getDataContainer(getOutputCellAttributeMatrixPath().getDataContainerName());

  std::vector<size_t> tDims(1, 0);
  m->createNonPrereqAttributeMatrix(this, getOutputCellFeatureAttributeMatrixName(), tDims, AttributeMatrix::Type::CellFeature);

  AttributeMatrix::Pointer outEnsembleAttrMat = AttributeMatrix::NullPointer();
  if(m->doesAttributeMatrixExist(getOutputCellEnsembleAttributeMatrixName()))
  {
    outEnsembleAttrMat = m->getPrereqAttributeMatrix(this, getOutputCellEnsembleAttributeMatrixName(), -350);
  }
  else
  {
    tDims[0] = m_PhaseTypesPtr.lock()->getNumberOfTuples();
    outEnsembleAttrMat = m->createNonPrereqAttributeMatrix(this, getOutputCellEnsembleAttributeMatrixName(), tDims, AttributeMatrix::Type::CellEnsemble);
  }

  tempPath = getOutputCellAttributeMatrixPath();
  tempPath.setAttributeMatrixName(getOutputCellEnsembleAttributeMatrixName());
  tempPath.setDataArrayName(SIMPL::EnsembleData::PhaseName);
  m_PhaseNamesPtr = getDataContainerArray()->createNonPrereqArrayFromPath<StringDataArray>(this, tempPath, nullptr, cDims);

  tempPath.update(getOutputCellAttributeMatrixPath().getDataContainerName(), getOutputCellAttributeMatrixPath().getAttributeMatrixName(), m_GBEuclideanDistancesArrayName);
  m_GBEuclideanDistancesPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<float>>(this, tempPath, -1, cDims);
  if(nullptr != m_GBEuclideanDistancesPtr.lock())
  {
    m_GBEuclideanDistances = m_GBEuclideanDistancesPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */

  tempPath.update(getOutputCellAttributeMatrixPath().getDataContainerName(), getOutputCellAttributeMatrixPath().getAttributeMatrixName(), getTJEuclideanDistancesArrayName());
  m_TJEuclideanDistancesPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<float>>(this, tempPath, -1, cDims);
  if(nullptr != m_TJEuclideanDistancesPtr.lock())
  {
    m_TJEuclideanDistances = m_TJEuclideanDistancesPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */

  tempPath.update(getOutputCellAttributeMatrixPath().getDataContainerName(), getOutputCellAttributeMatrixPath().getAttributeMatrixName(), getQPEuclideanDistancesArrayName());
  m_QPEuclideanDistancesPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<float>>(this, tempPath, -1, cDims);
  if(nullptr != m_QPEuclideanDistancesPtr.lock())
  {
    m_QPEuclideanDistances = m_QPEuclideanDistancesPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */

  cDims[0] = 3;
  tempPath.update(getOutputCellAttributeMatrixPath().getDataContainerName(), getOutputCellAttributeMatrixPath().getAttributeMatrixName(), m_NearestNeighborsArrayName);
  m_NearestNeighborsPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<int32_t>>(this, tempPath, 0, cDims);
  if(nullptr != m_NearestNeighborsPtr.lock())
  {
    m_NearestNeighbors = m_NearestNeighborsPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */

  // Feature Data
  tempPath.update(getOutputCellAttributeMatrixPath().getDataContainerName(), getOutputCellFeatureAttributeMatrixName(), m_FeaturePhasesArrayName);
  m_FeaturePhasesPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<int32_t>>(this, tempPath, 0, cDims);
  if(nullptr != m_FeaturePhasesPtr.lock())
  {
    m_FeaturePhases = m_FeaturePhasesPtr.lock()->getPointer(0);
  }

  tempPath.update(getOutputCellAttributeMatrixPath().getDataContainerName(), getOutputCellFeatureAttributeMatrixName(), m_NeighborhoodsArrayName);
  m_NeighborhoodsPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<int32_t>>(this, tempPath, 0, cDims);
  if(nullptr != m_NeighborhoodsPtr.lock())
  {
    m_Neighborhoods = m_NeighborhoodsPtr.lock()->getPointer(0);
  }

  tempPath.update(getOutputCellAttributeMatrixPath().getDataContainerName(), getOutputCellFeatureAttributeMatrixName(), m_EquivalentDiametersArrayName);
  m_EquivalentDiametersPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<float>>(this, tempPath, 0, cDims);
  if(nullptr != m_EquivalentDiametersPtr.lock())
  {
    m_EquivalentDiameters = m_EquivalentDiametersPtr.lock()->getPointer(0);
  }

  tempPath.update(getOutputCellAttributeMatrixPath().getDataContainerName(), getOutputCellFeatureAttributeMatrixName(), m_VolumesArrayName);
  m_VolumesPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<float>>(this, tempPath, 0, cDims);
  if(nullptr != m_VolumesPtr.lock())
  {
    m_Volumes = m_VolumesPtr.lock()->getPointer(0);
  }

  tempPath.update(getOutputCellAttributeMatrixPath().getDataContainerName(), getOutputCellFeatureAttributeMatrixName(), m_Omega3sArrayName);
  m_Omega3sPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<float>>(this, tempPath, 0, cDims);
  if(nullptr != m_Omega3sPtr.lock())
  {
    m_Omega3s = m_Omega3sPtr.lock()->getPointer(0);
  }

  cDims[0] = 3;
  tempPath.update(getOutputCellAttributeMatrixPath().getDataContainerName(), getOutputCellFeatureAttributeMatrixName(), m_CentroidsArrayName);
  m_CentroidsPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<float>>(this, tempPath, 0, cDims);
  if(nullptr != m_CentroidsPtr.lock())
  {
    m_Centroids = m_CentroidsPtr.lock()->getPointer(0);
  }

  tempPath.update(getOutputCellAttributeMatrixPath().getDataContainerName(), getOutputCellFeatureAttributeMatrixName(), m_AxisEulerAnglesArrayName);
  m_AxisEulerAnglesPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<float>>(this, tempPath, 0, cDims);
  if(nullptr != m_AxisEulerAnglesPtr.lock())
  {
    m_AxisEulerAngles = m_AxisEulerAnglesPtr.lock()->getPointer(0);
  }

  tempPath.update(getOutputCellAttributeMatrixPath().getDataContainerName(), getOutputCellFeatureAttributeMatrixName(), m_AxisLengthsArrayName);
  m_AxisLengthsPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<float>>(this, tempPath, 0, cDims);
  if(nullptr != m_AxisLengthsPtr.lock())
  {
    m_AxisLengths = m_AxisLengthsPtr.lock()->getPointer(0);
  }

  // Ensemble Data
  cDims[0] = 1;
  tempPath.update(getOutputCellAttributeMatrixPath().getDataContainerName(), getOutputCellEnsembleAttributeMatrixName(), getNumFeaturesArrayName());
  m_NumFeaturesPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<int32_t>>(this, tempPath, 0, cDims);
  if(nullptr != m_NumFeaturesPtr.lock())
  {
    m_NumFeatures = m_NumFeaturesPtr.lock()->getPointer(0);
  }

  if(m_WriteGoalAttributes && getCsvOutputFile().isEmpty())
  {
    QString ss = QObject::tr("The Csv output file must be set").arg(ClassName());
    setErrorCondition(-1, ss);
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void EstablishFoamMorphology::preflight()
{
  setInPreflight(true);
  Q_EMIT preflightAboutToExecute();
  Q_EMIT updateFilterParameters(this);
  dataCheck();
  Q_EMIT preflightExecuted();

  DataContainer::Pointer dc = getDataContainerArray()->getDataContainer(getOutputCellAttributeMatrixPath());
  if(dc == nullptr)
  {
    setInPreflight(false);
    return;
  }
  AttributeMatrix::Pointer featureAttrMat = dc->getAttributeMatrix(getOutputCellFeatureAttributeMatrixName());
  if(featureAttrMat == nullptr)
  {
    setInPreflight(false);
    return;
  }
  AttributeMatrix::Pointer cellAttrMat = dc->getAttributeMatrix(getOutputCellAttributeMatrixPath().getAttributeMatrixName());
  if(cellAttrMat == nullptr)
  {
    setInPreflight(false);
    return;
  }
  cellAttrMat->removeAttributeArray(m_NearestNeighborsArrayName);
  featureAttrMat->removeAttributeArray(m_AxisEulerAnglesArrayName);
  featureAttrMat->removeAttributeArray(m_AxisLengthsArrayName);
  featureAttrMat->removeAttributeArray(m_CentroidsArrayName);
  featureAttrMat->removeAttributeArray(m_EquivalentDiametersArrayName);
  featureAttrMat->removeAttributeArray(m_NeighborhoodsArrayName);
  featureAttrMat->removeAttributeArray(m_Omega3sArrayName);
  featureAttrMat->removeAttributeArray(m_FeaturePhasesArrayName);
  featureAttrMat->removeAttributeArray(m_VolumesArrayName);

  setInPreflight(false);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void EstablishFoamMorphology::execute()
{
  initialize();

  clearErrorCode();
  dataCheck();
  if(getErrorCode() < 0)
  {
    return;
  }

  if(m_HaveFeatures == 0)
  {
    notifyStatusMessage("Packing Features || Initializing Volume");
    // this initializes the arrays to hold the details of the locations of all of the features during packing
    Int32ArrayType::Pointer featureOwnersPtr = initialize_packinggrid();
    if(getErrorCode() < 0)
    {
      return;
    }
    notifyStatusMessage("Packing Features || Placing Features");
    place_features(featureOwnersPtr);
    if(getErrorCode() < 0)
    {
      return;
    }
    if(getCancel())
    {
      return;
    }

    notifyStatusMessage("Packing Features || Assigning Voxels");
    assign_voxels();
    if(getErrorCode() < 0)
    {
      return;
    }
    if(getCancel())
    {
      return;
    }

    notifyStatusMessage("Packing Features || Assigning Gaps");
    assign_gaps_only();
    if(getErrorCode() < 0)
    {
      return;
    }
    if(getCancel())
    {
      return;
    }
  }

  notifyStatusMessage("Foam Operations || Determining Strut Formation Options");
  find_euclidean_dist_map();
  if(getErrorCode() < 0)
  {
    return;
  }
  if(getCancel())
  {
    return;
  }

  notifyStatusMessage("Foam Operations || Forming Struts");
  form_struts();
  if(getErrorCode() < 0)
  {
    return;
  }
  if(getCancel())
  {
    return;
  }

  size_t totalPoints = m_FeatureIdsPtr.lock()->getNumberOfTuples();
  int32_t featurename = 0;
  int32_t feature = 0;
  for(size_t i = 0; i < totalPoints; i++)
  {
    featurename = m_FeatureIds[i];
    if(featurename == -1)
    {
      featurename = 0;
      m_FeatureIds[i] = feature;
    }
  }

  if(getCancel())
  {
    return;
  }

  if(m_WriteGoalAttributes)
  {
    write_goal_attributes();
  }
  if(getErrorCode() < 0)
  {
    return;
  }

  DataContainer::Pointer m = getDataContainerArray()->getDataContainer(getOutputCellAttributeMatrixPath().getDataContainerName());
  AttributeMatrix::Pointer cellAttrMat = m->getAttributeMatrix(getOutputCellAttributeMatrixPath().getAttributeMatrixName());
  AttributeMatrix::Pointer featureAttrMat = m->getAttributeMatrix(getOutputCellFeatureAttributeMatrixName());
  cellAttrMat->removeAttributeArray(m_NearestNeighborsArrayName);
  featureAttrMat->removeAttributeArray(m_AxisEulerAnglesArrayName);
  featureAttrMat->removeAttributeArray(m_AxisLengthsArrayName);
  featureAttrMat->removeAttributeArray(m_CentroidsArrayName);
  featureAttrMat->removeAttributeArray(m_EquivalentDiametersArrayName);
  featureAttrMat->removeAttributeArray(m_NeighborhoodsArrayName);
  featureAttrMat->removeAttributeArray(m_Omega3sArrayName);
  featureAttrMat->removeAttributeArray(m_FeaturePhasesArrayName);
  featureAttrMat->removeAttributeArray(m_VolumesArrayName);

  AttributeMatrix::Pointer ensembleAttrMat = getDataContainerArray()->getAttributeMatrix(getInputPhaseNamesArrayPath());
  IDataArray::Pointer inputPhaseNames = ensembleAttrMat->getAttributeArray(getInputPhaseNamesArrayPath().getDataArrayName());
  if(inputPhaseNames.get() != nullptr)
  {
    AttributeMatrix::Pointer cellEnsembleAttrMat = m->getAttributeMatrix(m_OutputCellEnsembleAttributeMatrixName);
    IDataArray::Pointer outputPhaseNames = inputPhaseNames->deepCopy();
    cellEnsembleAttrMat->addOrReplaceAttributeArray(outputPhaseNames);
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void EstablishFoamMorphology::place_features(const Int32ArrayType::Pointer& featureOwnersPtr)
{
  bool writeErrorFile = false;
  std::ofstream outFile;
  if(!m_ErrorOutputFile.isEmpty())
  {
    outFile.open(m_ErrorOutputFile.toLatin1().data(), std::ios_base::binary);
    writeErrorFile = outFile.is_open();
  }

  clearErrorCode();
  m_Seed = QDateTime::currentMSecsSinceEpoch();
  SIMPL_RANDOMNG_NEW_SEEDED(m_Seed);

  DataContainer::Pointer m = getDataContainerArray()->getDataContainer(getOutputCellAttributeMatrixPath().getDataContainerName());

  size_t totalFeatures = m_FeaturePhasesPtr.lock()->getNumberOfTuples();

  // Get the number of input ensembles from one of the input arrays that are located in the Input Ensemble AttributeMatrix
  size_t totalEnsembles = m_PhaseTypesPtr.lock()->getNumberOfTuples();

  StatsDataArray& statsDataArray = *(m_StatsDataArray.lock());

  SizeVec3Type udims = m->getGeometryAs<ImageGeom>()->getDimensions();

  FloatVec3Type res = m->getGeometryAs<ImageGeom>()->getSpacing();

  m_SizeX = static_cast<float>(udims[0] * res[0]);
  m_SizeY = static_cast<float>(udims[1] * res[1]);
  m_SizeZ = static_cast<float>(udims[2] * res[2]);
  m_TotalVol = m_SizeX * m_SizeY * m_SizeZ;

  // Making a double to prevent float overflow on incrementing
  double totalprecipitatevolTEMP = 0;
  size_t totalVox = udims[0] * udims[1] * udims[2];
  for(size_t i = 0; i < totalVox; i++)
  {
    if(m_FeatureIds[i] <= 0)
    {
      totalprecipitatevolTEMP++;
    }
  }
  float totalprecipitatevol = static_cast<float>(totalprecipitatevolTEMP);
  totalprecipitatevol = totalprecipitatevol * (res[0] * res[1] * res[2]);

  float change = 0.0f;
  int32_t phase = 0;
  int32_t randomfeature = 0;
  float xc = 0.0f;
  float yc = 0.0f;
  float zc = 0.0f;
  float oldxc = 0.0f;
  float oldyc = 0.0f;
  float oldzc = 0.0f;
  m_OldFillingError = 0.0f;
  //	m_CurrentNeighborhoodError = 0.0f, m_OldNeighborhoodError = 0.0f;
  //	m_CurrentSizeDistError = 0.0f, m_OldSizeDistError = 0.0f;
  int32_t acceptedmoves = 0;
  float totalprecipitatefractions = 0.0f;

  // find which phases are precipitate phases
  for(size_t i = 1; i < totalEnsembles; ++i)
  {
    if(m_PhaseTypes[i] == static_cast<PhaseType::EnumType>(PhaseType::Type::Precipitate))
    {
      PrecipitateStatsData::Pointer pp = std::dynamic_pointer_cast<PrecipitateStatsData>(statsDataArray[i]);
      if(nullptr == pp)
      {
        QString ss = QObject::tr("Tried to cast a statsDataArray[%1].get() to a PrecipitateStatsData* "
                                 "pointer but this resulted in a nullptr pointer. The value at m_PhaseTypes[%2] = %3 does not match up "
                                 "with the type of pointer stored in the StatsDataArray (PrecipitateStatsData)\n")
                         .arg(i)
                         .arg(i)
                         .arg(m_PhaseTypes[i]);
        setErrorCondition(-666, ss);
        return;
      }
      m_PrecipitatePhases.push_back(static_cast<int32_t>(i));
      m_PrecipitatePhaseFractions.push_back(pp->getPhaseFraction());
      totalprecipitatefractions = totalprecipitatefractions + pp->getPhaseFraction();
    }
  }
  // scale the precipitate phase fractions to total to 1
  for(auto& precipPhaseFraction : m_PrecipitatePhaseFractions)
  {
    precipPhaseFraction = precipPhaseFraction / totalprecipitatefractions;
  }

  std::vector<size_t> cDim(1, 1);
  Int32ArrayType::Pointer exclusionOwnersPtr = Int32ArrayType::CreateArray(m_TotalPackingPoints, cDim, "_INTERNAL_USE_ONLY_PackPrecipitateFeatures::exclusions_owners", true);
  exclusionOwnersPtr->initializeWithValue(0);

  // This is the set that we are going to keep updated with the points that are not in an exclusion zone
  std::map<size_t, size_t> availablePoints;
  std::map<size_t, size_t> availablePointsInv;

  // Get a pointer to the Feature Owners that was just initialized in the initialize_packinggrid() method
  int32_t* featureOwners = featureOwnersPtr->getPointer(0);
  int32_t* exclusionOwners = exclusionOwnersPtr->getPointer(0);
  int64_t featureOwnersIdx = 0;

  // determine initial set of available points
  m_AvailablePointsCount = 0;
  for(int64_t i = 0; i < m_TotalPackingPoints; i++)
  {
    if(exclusionOwners[i] == 0)
    {
      availablePoints[i] = m_AvailablePointsCount;
      availablePointsInv[m_AvailablePointsCount] = i;
      m_AvailablePointsCount++;
    }
  }
  // and clear the pointsToRemove and pointsToAdd vectors from the initial packing
  m_PointsToRemove.clear();
  m_PointsToAdd.clear();

  // initialize the sim and goal size distributions for the precipitate phases
  m_FeatureSizeDist.resize(m_PrecipitatePhases.size());
  m_SimFeatureSizeDist.resize(m_PrecipitatePhases.size());
  m_FeatureSizeDistStep.resize(m_PrecipitatePhases.size());
  size_t numPrecipitatePhases = m_PrecipitatePhases.size();
  for(size_t i = 0; i < numPrecipitatePhases; i++)
  {
    phase = m_PrecipitatePhases[i];
    PrecipitateStatsData::Pointer pp = std::dynamic_pointer_cast<PrecipitateStatsData>(statsDataArray[phase]);
    if(pp == nullptr)
    {
    }
    m_FeatureSizeDist[i].resize(40);
    m_SimFeatureSizeDist[i].resize(40);
    m_FeatureSizeDistStep[i] = static_cast<float>(((2 * pp->getMaxFeatureDiameter()) - (pp->getMinFeatureDiameter() / 2.0f)) / m_FeatureSizeDist[i].size());
    float input = 0.0f;
    float previoustotal = 0.0f;
    VectorOfFloatArray GSdist = pp->getFeatureSizeDistribution();
    float avg = GSdist[0]->getValue(0);
    float stdev = GSdist[1]->getValue(0);
    float denominatorConst = 1.0f / sqrtf(2.0f * stdev * stdev); // Calculate it here rather than calculating the same thing multiple times below
    size_t numFeatureSizeDist = m_FeatureSizeDist[i].size();
    for(size_t j = 0; j < numFeatureSizeDist; j++)
    {
      input = (float(j + 1) * m_FeatureSizeDistStep[i]) + (pp->getMinFeatureDiameter() / 2.0f);
      float logInput = logf(input);
      if(logInput <= avg)
      {
        m_FeatureSizeDist[i][j] = 0.5f - 0.5f * (SIMPLibMath::erf((avg - logInput) * denominatorConst)) - previoustotal;
      }
      if(logInput > avg)
      {
        m_FeatureSizeDist[i][j] = 0.5f + 0.5f * (SIMPLibMath::erf((logInput - avg) * denominatorConst)) - previoustotal;
      }
      previoustotal = previoustotal + m_FeatureSizeDist[i][j];
    }
  }

  if(getCancel())
  {
    return;
  }

  // generate the features and monitor the size distribution error while doing so. After features are generated, no new features can enter or leave the structure.
  Feature_t feature;

  // Estimate the total Number of features here
  int32_t estNumFeatures = estimate_numfeatures(udims[0], udims[1], udims[2], res[0], res[1], res[2]);
  std::vector<size_t> tDims(1, estNumFeatures);
  m->getAttributeMatrix(getOutputCellFeatureAttributeMatrixName())->resizeAttributeArrays(tDims);
  // need to update pointers after resize, buut do not need to run full data check because pointers are still valid
  updateFeatureInstancePointers();

  int32_t gid = 1;
  m_FirstFoamFeature = gid;
  std::vector<float> curphasevol;
  curphasevol.resize(m_PrecipitatePhases.size());
  float factor = 1.0f;
  size_t iter = 0;
  for(size_t j = 0; j < numPrecipitatePhases; ++j)
  {
    curphasevol[j] = 0;
    float curphasetotalvol = totalprecipitatevol * m_PrecipitatePhaseFractions[j];
    while(curphasevol[j] < (factor * curphasetotalvol))
    {
      iter++;
      m_Seed++;
      phase = m_PrecipitatePhases[j];
      generate_feature(phase, &feature, m_ShapeTypes[phase]);
      m_CurrentSizeDistError = check_sizedisterror(&feature);
      change = (m_CurrentSizeDistError) - (m_OldSizeDistError);
      if(change > 0.0f || m_CurrentSizeDistError > (1.0f - (float(iter) * 0.001f)) || curphasevol[j] < (0.75f * factor * curphasetotalvol))
      {
        if(gid % 100 == 0)
        {
          QString ss = QObject::tr("Packing Features (1/2) || Generating Feature #%1").arg(gid);
          notifyStatusMessage(ss);
        }
        if(gid + 1 >= static_cast<int32_t>(m->getAttributeMatrix(getOutputCellFeatureAttributeMatrixName())->getNumberOfTuples()))
        {
          tDims[0] = static_cast<size_t>(gid + 1);
          m->getAttributeMatrix(getOutputCellFeatureAttributeMatrixName())->resizeAttributeArrays(tDims);
          // need to update pointers after resize, buut do not need to run full data check because pointers are still valid
          updateFeatureInstancePointers();
        }

        transfer_attributes(gid, &feature);
        m_OldSizeDistError = m_CurrentSizeDistError;
        curphasevol[j] = curphasevol[j] + m_Volumes[gid];
        iter = 0;
        gid++;
      }
      if(getCancel())
      {
        return;
      }
    }
  }

  if(m_PeriodicBoundaries)
  {
    iter = 0;
    int32_t xfeatures = 0;
    int32_t yfeatures = 0;
    int32_t zfeatures = 0;
    xfeatures = static_cast<int32_t>(powf((m->getAttributeMatrix(getOutputCellFeatureAttributeMatrixName())->getNumberOfTuples() * (m_SizeX / m_SizeY) * (m_SizeX / m_SizeZ)), (1.0f / 3.0f)) + 1);
    yfeatures = static_cast<int32_t>(xfeatures * (m_SizeY / m_SizeX) + 1);
    zfeatures = static_cast<int32_t>(xfeatures * (m_SizeZ / m_SizeX) + 1);
    factor = 0.25f * (1.0f - (float((xfeatures - 2) * (yfeatures - 2) * (zfeatures - 2)) / float(xfeatures * yfeatures * zfeatures)));
    for(size_t j = 0; j < numPrecipitatePhases; ++j)
    {
      float curphasetotalvol = totalprecipitatevol * m_PrecipitatePhaseFractions[j];
      while(curphasevol[j] < ((1 + factor) * curphasetotalvol))
      {
        iter++;
        m_Seed++;
        phase = m_PrecipitatePhases[j];
        generate_feature(phase, &feature, m_ShapeTypes[phase]);
        m_CurrentSizeDistError = check_sizedisterror(&feature);
        change = (m_CurrentSizeDistError) - (m_OldSizeDistError);
        if(change > 0 || m_CurrentSizeDistError > (1.0f - (iter * 0.001f)) || curphasevol[j] < (0.75f * factor * curphasetotalvol))
        {
          QString ss = QObject::tr("Packing Features (2/2) || Generating Feature #%1").arg(gid);
          notifyStatusMessage(ss);
          if(gid + 1 >= static_cast<int32_t>(m->getAttributeMatrix(getOutputCellFeatureAttributeMatrixName())->getNumberOfTuples()))
          {
            tDims[0] = static_cast<size_t>(gid + 1);
            m->getAttributeMatrix(getOutputCellFeatureAttributeMatrixName())->resizeAttributeArrays(tDims);
            // need to update pointers after resize, but do not need to run full data check because pointers are still valid
            updateFeatureInstancePointers();
          }
          transfer_attributes(gid, &feature);
          m_OldSizeDistError = m_CurrentSizeDistError;
          curphasevol[j] = curphasevol[j] + m_Volumes[gid];
          iter = 0;
          gid++;
        }
        if(getCancel())
        {
          return;
        }
      }
    }
  }

  QString ss = QObject::tr("Packing Features || Starting Feature Placement...");
  notifyStatusMessage(ss);

  tDims[0] = static_cast<size_t>(gid);
  m->getAttributeMatrix(getOutputCellFeatureAttributeMatrixName())->resizeAttributeArrays(tDims);
  totalFeatures = m->getAttributeMatrix(getOutputCellFeatureAttributeMatrixName())->getNumberOfTuples();
  // need to update pointers after resize, buut do not need to run full data check because pointers are still valid
  updateFeatureInstancePointers();

  if(getCancel())
  {
    return;
  }

  // initialize the sim and goal neighbor distribution for the precipitate phases
  m_NeighborDist.resize(m_PrecipitatePhases.size());
  m_SimNeighborDist.resize(m_PrecipitatePhases.size());
  m_NeighborDistStep.resize(m_PrecipitatePhases.size());
  for(size_t i = 0; i < numPrecipitatePhases; i++)
  {
    phase = m_PrecipitatePhases[i];
    PrecipitateStatsData::Pointer pp = std::dynamic_pointer_cast<PrecipitateStatsData>(statsDataArray[phase]);
    m_NeighborDist[i].resize(pp->getBinNumbers()->getSize());
    m_SimNeighborDist[i].resize(pp->getBinNumbers()->getSize());
  }

  if(getCancel())
  {
    return;
  }

  m_ColumnList.resize(totalFeatures);
  m_RowList.resize(totalFeatures);
  m_PlaneList.resize(totalFeatures);
  m_EllipFuncList.resize(totalFeatures);
  m_PackQualities.resize(totalFeatures);
  m_FillingError = 1.0f;

  int64_t count = 0;
  int64_t column = 0;
  int64_t row = 0;
  int64_t plane = 0;
  int32_t progFeature = 0;
  int32_t progFeatureInc = static_cast<int32_t>(totalFeatures * 0.01f);
  for(size_t i = m_FirstFoamFeature; i < totalFeatures; i++)
  {
    if(getCancel())
    {
      return;
    }

    if((int32_t)i > progFeature + progFeatureInc)
    {
      ss = QObject::tr("Placing Feature #%1/%2").arg(i).arg(totalFeatures);
      notifyStatusMessage(ss);
      progFeature = i;
    }

    if(i == (totalFeatures - 1))
    {
      ss = QObject::tr("Placing Feature #%1/%2").arg(i + 1).arg(totalFeatures);
      notifyStatusMessage(ss);
    }

    // we always put the feature in the center of the box to make sure the feature has the optimal chance to not touch the edge of the box
    // this is because in the insert function below we are just determining which voxels relative to the centroid are in the feature
    xc = m_SizeX * 0.5f;
    yc = m_SizeY * 0.5f;
    zc = m_SizeZ * 0.5f;
    m_Centroids[3 * i] = xc;
    m_Centroids[3 * i + 1] = yc;
    m_Centroids[3 * i + 2] = zc;
    insert_feature(i);
    if(getErrorCode() < 0)
    {
      return;
    }
    count = 0;
    // now we randomly pick a place to try to place the feature
    xc = static_cast<float>(rg.genrand_res53() * m_SizeX);
    yc = static_cast<float>(rg.genrand_res53() * m_SizeY);
    zc = static_cast<float>(rg.genrand_res53() * m_SizeZ);
    column = static_cast<int64_t>((xc - (m_HalfPackingRes[0])) * m_OneOverPackingRes[0]);
    row = static_cast<int64_t>((yc - (m_HalfPackingRes[1])) * m_OneOverPackingRes[1]);
    plane = static_cast<int64_t>((zc - (m_HalfPackingRes[2])) * m_OneOverPackingRes[2]);
    featureOwnersIdx = (m_PackingPoints[0] * m_PackingPoints[1] * plane) + (m_PackingPoints[0] * row) + column;
    // now we walk til we find a point that is not in an exclusion zone
    while(exclusionOwners[featureOwnersIdx] > 0 && count < m_TotalPackingPoints)
    {
      featureOwnersIdx++;
      if(featureOwnersIdx >= m_TotalPackingPoints)
      {
        featureOwnersIdx = 0;
      }
      count++;
    }
    column = static_cast<int64_t>(featureOwnersIdx % m_PackingPoints[0]);
    row = static_cast<int64_t>(featureOwnersIdx / m_PackingPoints[0]) % m_PackingPoints[1];
    plane = static_cast<int64_t>(featureOwnersIdx / (m_PackingPoints[0] * m_PackingPoints[1]));
    xc = static_cast<float>((column * m_PackingRes[0]) + (m_PackingRes[0] * 0.5));
    yc = static_cast<float>((row * m_PackingRes[1]) + (m_PackingRes[1] * 0.5));
    zc = static_cast<float>((plane * m_PackingRes[2]) + (m_PackingRes[2] * 0.5));
    move_feature(i, xc, yc, zc);
    m_FillingError = check_fillingerror(i, -1000, featureOwnersPtr, exclusionOwnersPtr);
  }

  progFeature = 0;
  progFeatureInc = static_cast<int32_t>(totalFeatures * 0.01f);
  uint64_t millis = QDateTime::currentMSecsSinceEpoch();
  uint64_t currentMillis = millis;
  uint64_t startMillis = millis;
  uint64_t estimatedTime = 0;
  float timeDiff = 0.0f;

  // begin swaping/moving/adding/removing features to try to improve packing
  int32_t totalAdjustments = static_cast<int32_t>(100 * (totalFeatures - 1));

  // determine initial set of available points
  m_AvailablePointsCount = 0;
  for(int64_t i = 0; i < m_TotalPackingPoints; i++)
  {
    if(exclusionOwners[i] == 0)
    {
      availablePoints[i] = m_AvailablePointsCount;
      availablePointsInv[m_AvailablePointsCount] = i;
      m_AvailablePointsCount++;
    }
  }

  // and clear the pointsToRemove and pointsToAdd vectors from the initial packing
  m_PointsToRemove.clear();
  m_PointsToAdd.clear();

  millis = QDateTime::currentMSecsSinceEpoch();
  startMillis = millis;
  bool good = false;
  size_t key = 0;
  float xshift = 0.0f;
  float yshift = 0.0f;
  float zshift = 0.0f;
  int32_t lastIteration = 0;
  for(int32_t iteration = 0; iteration < totalAdjustments; ++iteration)
  {
    currentMillis = QDateTime::currentMSecsSinceEpoch();
    if(currentMillis - millis > 1000)
    {
      ss = QObject::tr("Swapping/Moving/Adding/Removing Features Iteration %1/%2").arg(iteration).arg(totalAdjustments);
      timeDiff = ((float)iteration / (float)(currentMillis - startMillis));
      estimatedTime = (float)(totalAdjustments - iteration) / timeDiff;

      ss += QObject::tr(" || Est. Time Remain: %1 || Iterations/Sec: %2").arg(DREAM3D::convertMillisToHrsMinSecs(estimatedTime)).arg(timeDiff * 1000);
      notifyStatusMessage(ss);

      millis = QDateTime::currentMSecsSinceEpoch();
      lastIteration = iteration;
    }

    if(getCancel())
    {
      return;
    }

    int32_t option = iteration % 2;

    if(writeErrorFile && iteration % 25 == 0)
    {
      outFile << iteration << " " << m_FillingError << "  " << availablePoints.size() << "  " << m_AvailablePointsCount << " " << totalFeatures << " " << acceptedmoves << "\n";
    }

    // JUMP - this option moves one feature to a random spot in the volume
    if(option == 0)
    {
      randomfeature = m_FirstFoamFeature + int32_t(rg.genrand_res53() * (totalFeatures - m_FirstFoamFeature));
      good = false;
      count = 0;
      while(!good && count < static_cast<int32_t>((totalFeatures - m_FirstFoamFeature)))
      {
        xc = m_Centroids[3 * randomfeature];
        yc = m_Centroids[3 * randomfeature + 1];
        zc = m_Centroids[3 * randomfeature + 2];
        column = static_cast<int64_t>((xc - (m_HalfPackingRes[0])) * m_OneOverPackingRes[0]);
        row = static_cast<int64_t>((yc - (m_HalfPackingRes[1])) * m_OneOverPackingRes[1]);
        plane = static_cast<int64_t>((zc - (m_HalfPackingRes[2])) * m_OneOverPackingRes[2]);
        featureOwnersIdx = (m_PackingPoints[0] * m_PackingPoints[1] * plane) + (m_PackingPoints[0] * row) + column;
        if(featureOwners[featureOwnersIdx] > 1)
        {
          good = true;
        }
        else
        {
          randomfeature++;
        }
        if(static_cast<size_t>(randomfeature) >= totalFeatures)
        {
          randomfeature = m_FirstFoamFeature;
        }
        count++;
      }
      m_Seed++;

      if(!availablePoints.empty())
      {
        key = static_cast<size_t>(rg.genrand_res53() * (m_AvailablePointsCount - 1));
        featureOwnersIdx = availablePointsInv[key];
      }
      else
      {
        featureOwnersIdx = static_cast<size_t>(rg.genrand_res53() * m_TotalPackingPoints);
      }

      // find the column row and plane of that point
      column = static_cast<int64_t>(featureOwnersIdx % m_PackingPoints[0]);
      row = static_cast<int64_t>(featureOwnersIdx / m_PackingPoints[0]) % m_PackingPoints[1];
      plane = static_cast<int64_t>(featureOwnersIdx / (m_PackingPoints[0] * m_PackingPoints[1]));
      xc = static_cast<float>((column * m_PackingRes[0]) + (m_PackingRes[0] * 0.5));
      yc = static_cast<float>((row * m_PackingRes[1]) + (m_PackingRes[1] * 0.5));
      zc = static_cast<float>((plane * m_PackingRes[2]) + (m_PackingRes[2] * 0.5));
      oldxc = m_Centroids[3 * randomfeature];
      oldyc = m_Centroids[3 * randomfeature + 1];
      oldzc = m_Centroids[3 * randomfeature + 2];
      m_OldFillingError = m_FillingError;
      m_FillingError = check_fillingerror(-1000, static_cast<int32_t>(randomfeature), featureOwnersPtr, exclusionOwnersPtr);
      move_feature(randomfeature, xc, yc, zc);
      m_FillingError = check_fillingerror(static_cast<int32_t>(randomfeature), -1000, featureOwnersPtr, exclusionOwnersPtr);

      if(m_FillingError > m_OldFillingError)
      {
        m_FillingError = check_fillingerror(-1000, static_cast<int32_t>(randomfeature), featureOwnersPtr, exclusionOwnersPtr);
        move_feature(randomfeature, oldxc, oldyc, oldzc);
        m_FillingError = check_fillingerror(static_cast<int32_t>(randomfeature), -1000, featureOwnersPtr, exclusionOwnersPtr);
        m_PointsToRemove.clear();
        m_PointsToAdd.clear();
      }
    }

    // NUDGE - this option moves one feature to a spot close to its current centroid
    if(option == 1)
    {
      randomfeature = m_FirstFoamFeature + int32_t(rg.genrand_res53() * (totalFeatures - m_FirstFoamFeature));
      good = false;
      count = 0;
      while(!good && count < static_cast<int32_t>((totalFeatures - m_FirstFoamFeature)))
      {
        xc = m_Centroids[3 * randomfeature];
        yc = m_Centroids[3 * randomfeature + 1];
        zc = m_Centroids[3 * randomfeature + 2];
        column = static_cast<int64_t>((xc - (m_HalfPackingRes[0])) * m_OneOverPackingRes[0]);
        row = static_cast<int64_t>((yc - (m_HalfPackingRes[1])) * m_OneOverPackingRes[1]);
        plane = static_cast<int64_t>((zc - (m_HalfPackingRes[2])) * m_OneOverPackingRes[2]);
        featureOwnersIdx = (m_PackingPoints[0] * m_PackingPoints[1] * plane) + (m_PackingPoints[0] * row) + column;
        if(featureOwners[featureOwnersIdx] > 1)
        {
          good = true;
        }
        else
        {
          randomfeature++;
        }
        if(static_cast<size_t>(randomfeature) >= totalFeatures)
        {
          randomfeature = m_FirstFoamFeature;
        }
        count++;
      }
      m_Seed++;
      oldxc = m_Centroids[3 * randomfeature];
      oldyc = m_Centroids[3 * randomfeature + 1];
      oldzc = m_Centroids[3 * randomfeature + 2];
      xshift = static_cast<float>(((2.0f * (rg.genrand_res53() - 0.5f)) * (2.0f * m_PackingRes[0])));
      yshift = static_cast<float>(((2.0f * (rg.genrand_res53() - 0.5f)) * (2.0f * m_PackingRes[1])));
      zshift = static_cast<float>(((2.0f * (rg.genrand_res53() - 0.5f)) * (2.0f * m_PackingRes[2])));
      if((oldxc + xshift) < m_SizeX && (oldxc + xshift) > 0)
      {
        xc = oldxc + xshift;
      }
      else
      {
        xc = oldxc;
      }
      if((oldyc + yshift) < m_SizeY && (oldyc + yshift) > 0)
      {
        yc = oldyc + yshift;
      }
      else
      {
        yc = oldyc;
      }
      if((oldzc + zshift) < m_SizeZ && (oldzc + zshift) > 0)
      {
        zc = oldzc + zshift;
      }
      else
      {
        zc = oldzc;
      }
      m_OldFillingError = m_FillingError;
      m_FillingError = check_fillingerror(-1000, static_cast<int32_t>(randomfeature), featureOwnersPtr, exclusionOwnersPtr);
      move_feature(randomfeature, xc, yc, zc);
      m_FillingError = check_fillingerror(static_cast<int32_t>(randomfeature), -1000, featureOwnersPtr, exclusionOwnersPtr);

      if(m_FillingError > m_OldFillingError)
      {
        m_FillingError = check_fillingerror(-1000, static_cast<int>(randomfeature), featureOwnersPtr, exclusionOwnersPtr);
        move_feature(randomfeature, oldxc, oldyc, oldzc);
        m_FillingError = check_fillingerror(static_cast<int>(randomfeature), -1000, featureOwnersPtr, exclusionOwnersPtr);
        m_PointsToRemove.clear();
        m_PointsToAdd.clear();
      }
    }
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
Int32ArrayType::Pointer EstablishFoamMorphology::initialize_packinggrid()
{
  DataContainer::Pointer m = getDataContainerArray()->getDataContainer(getOutputCellAttributeMatrixPath().getDataContainerName());

  m_PackingRes = m->getGeometryAs<ImageGeom>()->getSpacing();
  m_PackingRes[0] *= 2.0f;
  m_PackingRes[1] *= 2.0f;
  m_PackingRes[2] *= 2.0f;

  m_HalfPackingRes[0] = m_PackingRes[0] * 0.5f;
  m_HalfPackingRes[1] = m_PackingRes[1] * 0.5f;
  m_HalfPackingRes[2] = m_PackingRes[2] * 0.5f;

  m_OneOverHalfPackingRes[0] = 1.0f / m_HalfPackingRes[0];
  m_OneOverHalfPackingRes[1] = 1.0f / m_HalfPackingRes[1];
  m_OneOverHalfPackingRes[2] = 1.0f / m_HalfPackingRes[2];

  m_OneOverPackingRes[0] = 1.0f / m_PackingRes[0];
  m_OneOverPackingRes[1] = 1.0f / m_PackingRes[1];
  m_OneOverPackingRes[2] = 1.0f / m_PackingRes[2];

  m_PackingPoints[0] = static_cast<int64_t>(m->getGeometryAs<ImageGeom>()->getXPoints() / 2);
  m_PackingPoints[1] = static_cast<int64_t>(m->getGeometryAs<ImageGeom>()->getYPoints() / 2);
  m_PackingPoints[2] = static_cast<int64_t>(m->getGeometryAs<ImageGeom>()->getZPoints() / 2);

  for(auto&& packingValue : m_PackingPoints)
  {
    if(packingValue == 0)
    {
      packingValue = 1;
    }
  }

  m_TotalPackingPoints = m_PackingPoints[0] * m_PackingPoints[1] * m_PackingPoints[2];

  Int32ArrayType::Pointer featureOwnersPtr = Int32ArrayType::CreateArray(m_TotalPackingPoints, std::string("_INTERNAL_USE_ONLY_PackPrecipitateFeatures::feature_owners"), true);
  featureOwnersPtr->initializeWithZeros();

  return featureOwnersPtr;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void EstablishFoamMorphology::generate_feature(int32_t phase, Feature_t* feature, uint32_t shapeclass)
{
  SIMPL_RANDOMNG_NEW_SEEDED(m_Seed)
  std::array<double, 3> randx3 = {rg.genrand_res53(), rg.genrand_res53(), rg.genrand_res53()};

  StatsDataArray& statsDataArray = *(m_StatsDataArray.lock());

  float r1 = 1.0f;
  float a2 = 0.0f;
  float a3 = 0.0f;
  float b2 = 0.0f;
  float b3 = 0.0f;
  float diam = 0.0f;
  float vol = 0.0f;
  bool volgood = false;
  float fourThirdsPiOverEight = static_cast<float>(((4.0f / 3.0f) * (SIMPLib::Constants::k_PiD)) / 8.0f);
  PrecipitateStatsData::Pointer pp = std::dynamic_pointer_cast<PrecipitateStatsData>(statsDataArray[phase]);
  VectorOfFloatArray GSdist = pp->getFeatureSizeDistribution();
  float avg = GSdist[0]->getValue(0);
  float stdev = GSdist[1]->getValue(0);
  while(!volgood)
  {
    volgood = true;
    diam = static_cast<float>(rg.genrand_norm(avg, stdev));
    diam = expf(diam);
    if(diam >= pp->getMaxFeatureDiameter())
    {
      volgood = false;
    }
    if(diam < pp->getMinFeatureDiameter())
    {
      volgood = false;
    }
    vol = fourThirdsPiOverEight * (diam * diam * diam);
  }
  int32_t diameter = int32_t((diam - pp->getMinFeatureDiameter()) / pp->getBinStepSize());
  float r2 = 0.0f;
  float r3 = 1.0f;
  VectorOfFloatArray bovera = pp->getFeatureSize_BOverA();
  VectorOfFloatArray covera = pp->getFeatureSize_COverA();
  if(diameter >= static_cast<int32_t>(bovera[0]->getSize()))
  {
    diameter = static_cast<int32_t>(bovera[0]->getSize()) - 1;
  }
  while(r2 < r3)
  {
    r2 = 0.0f, r3 = 0.0f;
    a2 = bovera[0]->getValue(diameter);
    b2 = bovera[1]->getValue(diameter);
    a3 = covera[0]->getValue(diameter);
    b3 = covera[1]->getValue(diameter);
    int32_t tmpDiameter = diameter;
    int32_t increment = -1;
    while(a2 == 0 || b2 == 0 || a3 == 0 || b3 == 0)
    {
      tmpDiameter += increment;
      if(tmpDiameter < 0)
      {
        tmpDiameter = diameter + 1;
        increment = 1;
      }
      if(tmpDiameter >= static_cast<int32_t>(bovera[0]->getSize()))
      {
        a2 = 1.0f;
        b2 = 0.0f;
        a3 = 1.0f;
        b3 = 0.0f;
        break;
      }
      a2 = bovera[0]->getValue(tmpDiameter);
      b2 = bovera[1]->getValue(tmpDiameter);
      a3 = covera[0]->getValue(tmpDiameter);
      b3 = covera[1]->getValue(tmpDiameter);
    }
    r2 = static_cast<float>(rg.genrand_beta(a2, b2));
    r3 = static_cast<float>(rg.genrand_beta(a3, b3));
  }
  FloatArrayType::Pointer axisodf = pp->getAxisOrientation();
  int32_t numbins = axisodf->getNumberOfTuples();
  float random = static_cast<float>(rg.genrand_res53());
  float totaldensity = 0.0f;
  int32_t bin = 0;
  for(int32_t j = 0; j < numbins; j++)
  {
    float density = axisodf->getValue(j);
    float td1 = totaldensity;
    totaldensity = totaldensity + density;
    if(random < totaldensity && random >= td1)
    {
      bin = j;
      break;
    }
  }
  OrientationD eulers = m_OrthoOps->determineEulerAngles(randx3.data(), bin);
  VectorOfFloatArray omega3 = pp->getFeatureSize_Omegas();
  float mf = omega3[0]->getValue(diameter);
  float s = omega3[1]->getValue(diameter);
  float omega3f = static_cast<float>(rg.genrand_beta(mf, s));
  if(shapeclass == static_cast<ShapeType::EnumType>(ShapeType::Type::Ellipsoid))
  {
    omega3f = 1;
  }

  feature->m_Volumes = vol;
  feature->m_EquivalentDiameters = diam;
  feature->m_AxisLengths[0] = r1;
  feature->m_AxisLengths[1] = r2;
  feature->m_AxisLengths[2] = r3;
  feature->m_AxisEulerAngles[0] = eulers[0];
  feature->m_AxisEulerAngles[1] = eulers[1];
  feature->m_AxisEulerAngles[2] = eulers[2];
  feature->m_Omega3s = omega3f;
  feature->m_FeaturePhases = phase;
  feature->m_Neighborhoods = 0;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void EstablishFoamMorphology::transfer_attributes(int32_t gnum, Feature_t* feature)
{
  m_Volumes[gnum] = feature->m_Volumes;
  m_EquivalentDiameters[gnum] = feature->m_EquivalentDiameters;
  m_AxisLengths[3 * gnum + 0] = feature->m_AxisLengths[0];
  m_AxisLengths[3 * gnum + 1] = feature->m_AxisLengths[1];
  m_AxisLengths[3 * gnum + 2] = feature->m_AxisLengths[2];
  m_AxisEulerAngles[3 * gnum + 0] = feature->m_AxisEulerAngles[0];
  m_AxisEulerAngles[3 * gnum + 1] = feature->m_AxisEulerAngles[1];
  m_AxisEulerAngles[3 * gnum + 2] = feature->m_AxisEulerAngles[2];
  m_Omega3s[gnum] = feature->m_Omega3s;
  m_FeaturePhases[gnum] = feature->m_FeaturePhases;
  m_Neighborhoods[gnum] = feature->m_Neighborhoods;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void EstablishFoamMorphology::move_feature(size_t gnum, float xc, float yc, float zc)
{
  int64_t occolumn = 0;
  int64_t ocrow = 0;
  int64_t ocplane = 0;
  int64_t nccolumn = 0;
  int64_t ncrow = 0;
  int64_t ncplane = 0;
  int64_t shiftcolumn = 0;
  int64_t shiftrow = 0;
  int64_t shiftplane = 0;
  float oxc = m_Centroids[3 * gnum];
  float oyc = m_Centroids[3 * gnum + 1];
  float ozc = m_Centroids[3 * gnum + 2];
  occolumn = static_cast<int64_t>((oxc - (m_HalfPackingRes[0])) * m_OneOverPackingRes[0]);
  ocrow = static_cast<int64_t>((oyc - (m_HalfPackingRes[1])) * m_OneOverPackingRes[1]);
  ocplane = static_cast<int64_t>((ozc - (m_HalfPackingRes[2])) * m_OneOverPackingRes[2]);
  nccolumn = static_cast<int64_t>((xc - (m_HalfPackingRes[0])) * m_OneOverPackingRes[0]);
  ncrow = static_cast<int64_t>((yc - (m_HalfPackingRes[1])) * m_OneOverPackingRes[1]);
  ncplane = static_cast<int64_t>((zc - (m_HalfPackingRes[2])) * m_OneOverPackingRes[2]);
  shiftcolumn = nccolumn - occolumn;
  shiftrow = ncrow - ocrow;
  shiftplane = ncplane - ocplane;
  m_Centroids[3 * gnum] = xc;
  m_Centroids[3 * gnum + 1] = yc;
  m_Centroids[3 * gnum + 2] = zc;
  size_t size = m_ColumnList[gnum].size();

  for(size_t i = 0; i < size; i++)
  {
    int64_t& cl = m_ColumnList[gnum][i];
    cl += shiftcolumn;
    int64_t& rl = m_RowList[gnum][i];
    rl += shiftrow;
    int64_t& pl = m_PlaneList[gnum][i];
    pl += shiftplane;
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void EstablishFoamMorphology::determine_neighbors(size_t gnum, bool add)
{
  DataContainer::Pointer m = getDataContainerArray()->getDataContainer(getOutputCellAttributeMatrixPath().getDataContainerName());

  float x = 0.0f;
  float y = 0.0f;
  float z = 0.0f;
  float xn = 0.0f;
  float yn = 0.0f;
  float zn = 0.0f;
  float dia = 0.0f;
  float dia2 = 0.0f;
  float dx = 0.0f;
  float dy = 0.0f;
  float dz = 0.0f;
  x = m_Centroids[3 * gnum];
  y = m_Centroids[3 * gnum + 1];
  z = m_Centroids[3 * gnum + 2];
  dia = m_EquivalentDiameters[gnum];
  size_t totalFeatures = m->getAttributeMatrix(getOutputCellFeatureAttributeMatrixName())->getNumberOfTuples();
  int32_t increment = 0;
  if(add)
  {
    increment = 1;
  }
  else
  {
    increment = -1;
  }
  for(size_t n = m_FirstFoamFeature; n < totalFeatures; n++)
  {
    xn = m_Centroids[3 * n];
    yn = m_Centroids[3 * n + 1];
    zn = m_Centroids[3 * n + 2];
    dia2 = m_EquivalentDiameters[n];
    dx = fabs(x - xn);
    dy = fabs(y - yn);
    dz = fabs(z - zn);
    if(dx < dia && dy < dia && dz < dia)
    {
      m_Neighborhoods[gnum] = m_Neighborhoods[gnum] + increment;
    }
    if(dx < dia2 && dy < dia2 && dz < dia2)
    {
      m_Neighborhoods[n] = m_Neighborhoods[n] + increment;
    }
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
float EstablishFoamMorphology::check_neighborhooderror(int32_t gadd, int32_t gremove)
{
  // Optimized Code
  DataContainer::Pointer m = getDataContainerArray()->getDataContainer(getOutputCellAttributeMatrixPath().getDataContainerName());

  StatsDataArray& statsDataArray = *(m_StatsDataArray.lock());

  float neighborerror = 0.0f;
  float bhattdist = 0.0f;
  float dia = 0.0f;
  int32_t nnum = 0;
  size_t diabin = 0;
  size_t nnumbin = 0;
  int32_t index = 0;
  int32_t phase = 0;

  using VectOfVectFloat_t = std::vector<std::vector<float>>;
  size_t numPhases = m_SimNeighborDist.size();
  for(size_t iter = 0; iter < numPhases; ++iter)
  {
    phase = m_PrecipitatePhases[iter];
    PrecipitateStatsData::Pointer pp = std::dynamic_pointer_cast<PrecipitateStatsData>(statsDataArray[phase]);
    VectOfVectFloat_t& curSimNeighborDist = m_SimNeighborDist[iter];
    size_t curSImNeighborDist_Size = curSimNeighborDist.size();
    float oneOverNeighborDistStep = 1.0f / m_NeighborDistStep[iter];

    std::vector<int32_t> count(curSImNeighborDist_Size, 0);
    for(size_t i = 0; i < curSImNeighborDist_Size; i++)
    {
      curSimNeighborDist[i].resize(40);
      for(size_t j = 0; j < 40; j++)
      {
        curSimNeighborDist[i][j] = 0;
      }
    }
    if(gadd > 0 && m_FeaturePhases[gadd] == phase)
    {
      determine_neighbors(gadd, true);
    }
    if(gremove > 0 && m_FeaturePhases[gremove] == phase)
    {
      determine_neighbors(gremove, false);
    }

    float maxFeatureDia = pp->getMaxFeatureDiameter();
    float minFeatureDia = pp->getMinFeatureDiameter();
    float oneOverBinStepSize = 1.0f / pp->getBinStepSize();

    size_t totalFeatures = m->getAttributeMatrix(getOutputCellFeatureAttributeMatrixName())->getNumberOfTuples();
    for(size_t i = m_FirstFoamFeature; i < totalFeatures; i++)
    {
      nnum = 0;
      index = static_cast<int32_t>(i);
      if(index != gremove && m_FeaturePhases[index] == phase)
      {
        dia = m_EquivalentDiameters[index];
        if(dia > maxFeatureDia)
        {
          dia = maxFeatureDia;
        }
        if(dia < minFeatureDia)
        {
          dia = minFeatureDia;
        }
        diabin = static_cast<size_t>(((dia - minFeatureDia) * oneOverBinStepSize));
        if(diabin >= curSImNeighborDist_Size)
        {
          diabin = curSImNeighborDist_Size - 1;
        }
        nnum = m_Neighborhoods[index];
        nnumbin = static_cast<size_t>(nnum * oneOverNeighborDistStep);
        if(nnumbin >= 40)
        {
          nnumbin = 39;
        }
        curSimNeighborDist[diabin][nnumbin]++;
        count[diabin]++;
      }
    }
    if(gadd > 0 && m_FeaturePhases[gadd] == phase)
    {
      dia = m_EquivalentDiameters[gadd];
      if(dia > maxFeatureDia)
      {
        dia = maxFeatureDia;
      }
      if(dia < minFeatureDia)
      {
        dia = minFeatureDia;
      }
      diabin = static_cast<size_t>(((dia - minFeatureDia) * oneOverBinStepSize));
      if(diabin >= curSImNeighborDist_Size)
      {
        diabin = curSImNeighborDist_Size - 1;
      }
      nnum = m_Neighborhoods[gadd];
      nnumbin = static_cast<size_t>(nnum * oneOverNeighborDistStep);
      if(nnumbin >= 40)
      {
        nnumbin = 39;
      }
      curSimNeighborDist[diabin][nnumbin]++;
      count[diabin]++;
    }
    float runningtotal = 0.0f;

    for(size_t i = 0; i < curSImNeighborDist_Size; i++)
    {
      if(count[i] == 0)
      {
        for(size_t j = 0; j < 40; j++)
        {
          curSimNeighborDist[i][j] = 0.0f;
        }
      }
      else
      {
        float oneOverCount = 1.0f / (float)(count[i]);
        for(size_t j = 0; j < 40; j++)
        {
          curSimNeighborDist[i][j] = curSimNeighborDist[i][j] * oneOverCount;
          runningtotal = runningtotal + curSimNeighborDist[i][j];
        }
      }
    }

    runningtotal = 1.0f / runningtotal; // Flip this so that we have a multiply instead of a divide.
    for(size_t i = 0; i < curSImNeighborDist_Size; i++)
    {
      for(size_t j = 0; j < 40; j++)
      {
        curSimNeighborDist[i][j] = curSimNeighborDist[i][j] * runningtotal;
      }
    }

    if(gadd > 0 && m_FeaturePhases[gadd] == phase)
    {
      determine_neighbors(gadd, false);
    }

    if(gremove > 0 && m_FeaturePhases[gremove] == phase)
    {
      determine_neighbors(gremove, true);
    }
  }
  compare_3Ddistributions(m_SimNeighborDist, m_NeighborDist, bhattdist);
  neighborerror = bhattdist;
  return neighborerror;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void EstablishFoamMorphology::compare_1Ddistributions(std::vector<float> array1, std::vector<float> array2, float& bhattdist)
{
  bhattdist = 0.0f;
  size_t array1Size = array1.size();
  for(size_t i = 0; i < array1Size; i++)
  {
    bhattdist = bhattdist + sqrtf((array1[i] * array2[i]));
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void EstablishFoamMorphology::compare_2Ddistributions(std::vector<std::vector<float>> array1, std::vector<std::vector<float>> array2, float& bhattdist)
{
  bhattdist = 0.0f;
  size_t array1Size = array1.size();
  for(size_t i = 0; i < array1Size; i++)
  {
    size_t array2Size = array1[i].size();
    for(size_t j = 0; j < array2Size; j++)
    {
      bhattdist = bhattdist + sqrtf((array1[i][j] * array2[i][j]));
    }
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void EstablishFoamMorphology::compare_3Ddistributions(std::vector<std::vector<std::vector<float>>> array1, std::vector<std::vector<std::vector<float>>> array2, float& bhattdist)
{
  bhattdist = 0.0f;
  size_t array1Size = array1.size();
  for(size_t i = 0; i < array1Size; i++)
  {
    size_t array2Size = array1[i].size();
    for(size_t j = 0; j < array2Size; j++)
    {
      size_t array3Size = array1[i][j].size();
      for(size_t k = 0; k < array3Size; k++)
      {
        bhattdist = bhattdist + sqrtf((array1[i][j][k] * array2[i][j][k]));
      }
    }
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
float EstablishFoamMorphology::check_sizedisterror(Feature_t* feature)
{
  DataContainer::Pointer m = getDataContainerArray()->getDataContainer(getOutputCellAttributeMatrixPath().getDataContainerName());

  StatsDataArray& statsDataArray = *(m_StatsDataArray.lock());

  float dia = 0.0f;
  float sizedisterror = 0.0f;
  float bhattdist = 0.0f;
  int32_t index = 0;
  int32_t count = 0;
  int32_t phase = 0;
  size_t featureSizeDist_Size = m_FeatureSizeDist.size();
  for(size_t iter = 0; iter < featureSizeDist_Size; ++iter)
  {
    phase = m_PrecipitatePhases[iter];
    PrecipitateStatsData::Pointer pp = std::dynamic_pointer_cast<PrecipitateStatsData>(statsDataArray[phase]);
    count = 0;
    std::vector<float>& curFeatureSizeDist = m_FeatureSizeDist[iter];
    std::vector<float>::size_type curFeatureSizeDistSize = curFeatureSizeDist.size();
    std::vector<float>& curSimFeatureSizeDist = m_SimFeatureSizeDist[iter];
    // Initialize all Values to Zero
    for(size_t i = 0; i < curFeatureSizeDistSize; i++)
    {
      curSimFeatureSizeDist[i] = 0.0f;
    }
    int64_t totalFeatures = m->getAttributeMatrix(getOutputCellFeatureAttributeMatrixName())->getNumberOfTuples();
    float oneOverCurFeatureSizeDistStep = 1.0f / m_FeatureSizeDistStep[iter];
    float halfMinFeatureDiameter = pp->getMinFeatureDiameter() * 0.5f;
    for(int32_t b = m_FirstFoamFeature; b < totalFeatures; b++)
    {
      index = b;
      if(m_FeaturePhases[index] == phase)
      {
        dia = m_EquivalentDiameters[index];
        dia = (dia - halfMinFeatureDiameter) * oneOverCurFeatureSizeDistStep;
        if(dia < 0)
        {
          dia = 0.0f;
        }
        if(dia > curFeatureSizeDistSize - 1.0f)
        {
          dia = curFeatureSizeDistSize - 1.0f;
        }
        curSimFeatureSizeDist[int32_t(dia)]++;
        count++;
      }
    }

    if(feature->m_FeaturePhases == phase)
    {
      dia = feature->m_EquivalentDiameters;
      dia = (dia - halfMinFeatureDiameter) * oneOverCurFeatureSizeDistStep;
      if(dia < 0)
      {
        dia = 0.0f;
      }
      if(dia > curFeatureSizeDistSize - 1.0f)
      {
        dia = curFeatureSizeDistSize - 1.0f;
      }
      curSimFeatureSizeDist[int32_t(dia)]++;
      count++;
    }
    float oneOverCount = 1.0f / count;

    if(count == 0)
    {
      for(size_t i = 0; i < curFeatureSizeDistSize; i++)
      {
        curSimFeatureSizeDist[i] = 0.0;
      }
    }
    else
    {
      for(size_t i = 0; i < curFeatureSizeDistSize; i++)
      {
        curSimFeatureSizeDist[i] = curSimFeatureSizeDist[i] * oneOverCount;
      }
    }
  }
  compare_2Ddistributions(m_SimFeatureSizeDist, m_FeatureSizeDist, bhattdist);
  sizedisterror = bhattdist;
  return sizedisterror;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
float EstablishFoamMorphology::check_fillingerror(int32_t gadd, int32_t gremove, const Int32ArrayType::Pointer& featureOwnersPtr, const Int32ArrayType::Pointer& exclusionOwnersPtr)
{
  size_t featureOwnersIdx = 0;
  int32_t* featureOwners = featureOwnersPtr->getPointer(0);
  int32_t* exclusionOwners = exclusionOwnersPtr->getPointer(0);

  m_FillingError = m_FillingError * float(m_TotalPackingPoints);
  int64_t col = 0;
  int64_t row = 0;
  int64_t plane = 0;
  int32_t k1 = 0;
  int32_t k2 = 0;
  int32_t k3 = 0;
  if(gadd > 0)
  {
    k1 = 2;
    k2 = -1;
    k3 = 1;
    size_t numVoxelsForCurrentGrain = m_ColumnList[gadd].size();
    std::vector<int64_t>& cl = m_ColumnList[gadd];
    std::vector<int64_t>& rl = m_RowList[gadd];
    std::vector<int64_t>& pl = m_PlaneList[gadd];
    std::vector<float>& efl = m_EllipFuncList[gadd];
    float packquality = 0;
    for(size_t i = 0; i < numVoxelsForCurrentGrain; i++)
    {
      col = cl[i];
      row = rl[i];
      plane = pl[i];
      if(m_PeriodicBoundaries)
      {
        if(col < 0)
        {
          col = col + m_PackingPoints[0];
        }
        if(col > m_PackingPoints[0] - 1)
        {
          col = col - m_PackingPoints[0];
        }
        if(row < 0)
        {
          row = row + m_PackingPoints[1];
        }
        if(row > m_PackingPoints[1] - 1)
        {
          row = row - m_PackingPoints[1];
        }
        if(plane < 0)
        {
          plane = plane + m_PackingPoints[2];
        }
        if(plane > m_PackingPoints[2] - 1)
        {
          plane = plane - m_PackingPoints[2];
        }
        featureOwnersIdx = (m_PackingPoints[0] * m_PackingPoints[1] * plane) + (m_PackingPoints[0] * row) + col;
        int32_t currentFeatureOwner = featureOwners[featureOwnersIdx];
        if(efl[i] > 0.1f)
        {
          if(exclusionOwners[featureOwnersIdx] == 0)
          {
            m_PointsToRemove.push_back(featureOwnersIdx);
          }
          exclusionOwners[featureOwnersIdx]++;
        }
        m_FillingError = static_cast<float>(m_FillingError + ((k1 * currentFeatureOwner + k2)));
        //        fillingerror = fillingerror + (multiplier * (k1 * currentFeatureOwner  + k2));
        featureOwners[featureOwnersIdx] = currentFeatureOwner + k3;
        packquality = static_cast<float>(packquality + ((currentFeatureOwner) * (currentFeatureOwner)));
      }
      else
      {
        if(col >= 0 && col < m_PackingPoints[0] && row >= 0 && row < m_PackingPoints[1] && plane >= 0 && plane < m_PackingPoints[2])
        {
          featureOwnersIdx = (m_PackingPoints[0] * m_PackingPoints[1] * plane) + (m_PackingPoints[0] * row) + col;
          int32_t currentFeatureOwner = featureOwners[featureOwnersIdx];
          if(exclusionOwners[featureOwnersIdx] > 0)
          {
          }
          if(efl[i] > 0.1f)
          {
            if(exclusionOwners[featureOwnersIdx] == 0)
            {
              m_PointsToRemove.push_back(featureOwnersIdx);
            }
            exclusionOwners[featureOwnersIdx]++;
          }
          m_FillingError = static_cast<float>(m_FillingError + ((k1 * currentFeatureOwner + k2)));
          //        fillingerror = fillingerror + (multiplier * (k1 * currentFeatureOwner  + k2));
          featureOwners[featureOwnersIdx] = currentFeatureOwner + k3;
          packquality = static_cast<float>(packquality + ((currentFeatureOwner) * (currentFeatureOwner)));
        }
      }
    }
    m_PackQualities[gadd] = static_cast<int64_t>(packquality / float(numVoxelsForCurrentGrain));
  }
  if(gremove > 0)
  {
    k1 = -2;
    k2 = 3;
    k3 = -1;
    size_t size = m_ColumnList[gremove].size();
    std::vector<int64_t>& cl = m_ColumnList[gremove];
    std::vector<int64_t>& rl = m_RowList[gremove];
    std::vector<int64_t>& pl = m_PlaneList[gremove];
    std::vector<float>& efl = m_EllipFuncList[gremove];
    for(size_t i = 0; i < size; i++)
    {
      col = cl[i];
      row = rl[i];
      plane = pl[i];
      if(m_PeriodicBoundaries)
      {
        if(col < 0)
        {
          col = col + m_PackingPoints[0];
        }
        if(col > m_PackingPoints[0] - 1)
        {
          col = col - m_PackingPoints[0];
        }
        if(row < 0)
        {
          row = row + m_PackingPoints[1];
        }
        if(row > m_PackingPoints[1] - 1)
        {
          row = row - m_PackingPoints[1];
        }
        if(plane < 0)
        {
          plane = plane + m_PackingPoints[2];
        }
        if(plane > m_PackingPoints[2] - 1)
        {
          plane = plane - m_PackingPoints[2];
        }
        featureOwnersIdx = (m_PackingPoints[0] * m_PackingPoints[1] * plane) + (m_PackingPoints[0] * row) + col;
        int32_t currentFeatureOwner = featureOwners[featureOwnersIdx];
        if(efl[i] > 0.1f)
        {
          exclusionOwners[featureOwnersIdx]--;
          if(exclusionOwners[featureOwnersIdx] == 0)
          {
            m_PointsToAdd.push_back(featureOwnersIdx);
          }
        }
        m_FillingError = static_cast<float>(m_FillingError + ((k1 * currentFeatureOwner + k2)));
        //        fillingerror = fillingerror + (multiplier * (k1 * currentFeatureOwner  + k2));
        featureOwners[featureOwnersIdx] = currentFeatureOwner + k3;
      }
      else
      {
        if(col >= 0 && col < m_PackingPoints[0] && row >= 0 && row < m_PackingPoints[1] && plane >= 0 && plane < m_PackingPoints[2])
        {
          featureOwnersIdx = (m_PackingPoints[0] * m_PackingPoints[1] * plane) + (m_PackingPoints[0] * row) + col;
          int32_t currentFeatureOwner = featureOwners[featureOwnersIdx];
          if(efl[i] > 0.1f)
          {
            exclusionOwners[featureOwnersIdx]--;
            if(exclusionOwners[featureOwnersIdx] == 0)
            {
              m_PointsToAdd.push_back(featureOwnersIdx);
            }
          }
          m_FillingError = static_cast<float>(m_FillingError + ((k1 * currentFeatureOwner + k2)));
          //          fillingerror = fillingerror + (multiplier * (k1 * currentFeatureOwner  + k2));
          featureOwners[featureOwnersIdx] = currentFeatureOwner + k3;
        }
      }
    }
  }
  m_FillingError = m_FillingError / float(m_TotalPackingPoints);
  return m_FillingError;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void EstablishFoamMorphology::update_availablepoints(std::map<size_t, size_t>& availablePoints, std::map<size_t, size_t>& availablePointsInv)
{
  size_t removeSize = m_PointsToRemove.size();
  size_t addSize = m_PointsToAdd.size();
  size_t featureOwnersIdx = 0;
  size_t key = 0;
  size_t val = 0;
  for(size_t i = 0; i < removeSize; i++)
  {
    featureOwnersIdx = m_PointsToRemove[i];
    key = availablePoints[featureOwnersIdx];
    //  availablePoints.erase(featureOwnersIdx);
    val = availablePointsInv[m_AvailablePointsCount - 1];
    //  availablePointsInv.erase(availablePointsCount-1);
    if(key < m_AvailablePointsCount - 1)
    {
      availablePointsInv[key] = val;
      availablePoints[val] = key;
    }
    m_AvailablePointsCount--;
  }
  for(size_t i = 0; i < addSize; i++)
  {
    featureOwnersIdx = m_PointsToAdd[i];
    availablePoints[featureOwnersIdx] = m_AvailablePointsCount;
    availablePointsInv[m_AvailablePointsCount] = featureOwnersIdx;
    m_AvailablePointsCount++;
  }
  m_PointsToRemove.clear();
  m_PointsToAdd.clear();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void EstablishFoamMorphology::insert_feature(size_t gnum)
{
  SIMPL_RANDOMNG_NEW();

  float inside = -1.0f;
  int64_t column = 0;
  int64_t row = 0;
  int64_t plane = 0;
  int64_t centercolumn = 0;
  int64_t centerrow = 0;
  int64_t centerplane = 0;
  int64_t xmin = 0;
  int64_t xmax = 0;
  int64_t ymin = 0;
  int64_t ymax = 0;
  int64_t zmin = 0;
  int64_t zmax = 0;
  float xc = 0.0f;
  float yc = 0.0f;
  float zc = 0.0f;
  std::array<float, 3> coordsRotated = {0.0f, 0.0f, 0.0f};
  std::array<float, 3> coords = {0.0f, 0.0f, 0.0f};
  float volcur = m_Volumes[gnum];
  float bovera = m_AxisLengths[3 * gnum + 1];
  float covera = m_AxisLengths[3 * gnum + 2];
  float omega3 = m_Omega3s[gnum];
  float radcur1 = 1.0f;
  uint32_t shapeclass = m_ShapeTypes[m_FeaturePhases[gnum]];

  // Bail if the shapeclass is not one of our enumerated types
  if(shapeclass >= static_cast<ShapeType::EnumType>(ShapeType::Type::ShapeTypeEnd))
  {
    QString ss = QObject::tr("Undefined shape class in shape types array with path %1").arg(m_InputShapeTypesArrayPath.serialize());
    setErrorCondition(-666, ss);
    return;
  }

  // init any values for each of the Shape Ops
  // init any values for each of the Shape Ops
  for(auto& shapeOp : m_ShapeOps)
  {
    shapeOp->init();
  }
  // Create our Argument Map
  QMap<ShapeOps::ArgName, float> shapeArgMap;
  shapeArgMap[ShapeOps::Omega3] = omega3;
  shapeArgMap[ShapeOps::VolCur] = volcur;
  shapeArgMap[ShapeOps::B_OverA] = bovera;
  shapeArgMap[ShapeOps::C_OverA] = covera;

  radcur1 = m_ShapeOps[shapeclass]->radcur1(shapeArgMap);

  float radcur2 = (radcur1 * bovera);
  float radcur3 = (radcur1 * covera);
  float phi1 = m_AxisEulerAngles[3 * gnum];
  float PHI = m_AxisEulerAngles[3 * gnum + 1];
  float phi2 = m_AxisEulerAngles[3 * gnum + 2];
  float ga[3][3] = {{0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f}};

  OrientationTransformation::eu2om<OrientationF, OrientationF>(OrientationF(phi1, PHI, phi2)).toGMatrix(ga);

  xc = m_Centroids[3 * gnum];
  yc = m_Centroids[3 * gnum + 1];
  zc = m_Centroids[3 * gnum + 2];
  centercolumn = static_cast<int64_t>((xc - (m_PackingRes[0] * 0.5f)) * m_OneOverPackingRes[0]);
  centerrow = static_cast<int64_t>((yc - (m_PackingRes[1] * 0.5f)) * m_OneOverPackingRes[1]);
  centerplane = static_cast<int64_t>((zc - (m_PackingRes[2] * 0.5f)) * m_OneOverPackingRes[2]);
  xmin = static_cast<int64_t>(centercolumn - ((radcur1 * m_OneOverPackingRes[0]) + 1));
  xmax = static_cast<int64_t>(centercolumn + ((radcur1 * m_OneOverPackingRes[0]) + 1));
  ymin = static_cast<int64_t>(centerrow - ((radcur1 * m_OneOverPackingRes[1]) + 1));
  ymax = static_cast<int64_t>(centerrow + ((radcur1 * m_OneOverPackingRes[1]) + 1));
  zmin = static_cast<int64_t>(centerplane - ((radcur1 * m_OneOverPackingRes[2]) + 1));
  zmax = static_cast<int64_t>(centerplane + ((radcur1 * m_OneOverPackingRes[2]) + 1));
  if(xmin < -m_PackingPoints[0])
  {
    xmin = -m_PackingPoints[0];
  }
  if(xmax > 2 * m_PackingPoints[0] - 1)
  {
    xmax = (2 * m_PackingPoints[0] - 1);
  }
  if(ymin < -m_PackingPoints[1])
  {
    ymin = -m_PackingPoints[1];
  }
  if(ymax > 2 * m_PackingPoints[1] - 1)
  {
    ymax = (2 * m_PackingPoints[1] - 1);
  }
  if(zmin < -m_PackingPoints[2])
  {
    zmin = -m_PackingPoints[2];
  }
  if(zmax > 2 * m_PackingPoints[2] - 1)
  {
    zmax = (2 * m_PackingPoints[2] - 1);
  }

  float OneOverRadcur1 = 1.0f / radcur1;
  float OneOverRadcur2 = 1.0f / radcur2;
  float OneOverRadcur3 = 1.0f / radcur3;
  for(int64_t iter1 = xmin; iter1 < xmax + 1; iter1++)
  {
    for(int64_t iter2 = ymin; iter2 < ymax + 1; iter2++)
    {
      for(int64_t iter3 = zmin; iter3 < zmax + 1; iter3++)
      {
        column = iter1;
        row = iter2;
        plane = iter3;
        coords[0] = float(column) * m_PackingRes[0];
        coords[1] = float(row) * m_PackingRes[1];
        coords[2] = float(plane) * m_PackingRes[2];
        inside = -1.0f;
        coords[0] = coords[0] - xc;
        coords[1] = coords[1] - yc;
        coords[2] = coords[2] - zc;
        MatrixMath::Multiply3x3with3x1(ga, coords.data(), coordsRotated.data());
        float axis1comp = coordsRotated[0] * OneOverRadcur1;
        float axis2comp = coordsRotated[1] * OneOverRadcur2;
        float axis3comp = coordsRotated[2] * OneOverRadcur3;
        inside = m_ShapeOps[shapeclass]->inside(axis1comp, axis2comp, axis3comp);
        if(inside >= 0)
        {
          m_ColumnList[gnum].push_back(column);
          m_RowList[gnum].push_back(row);
          m_PlaneList[gnum].push_back(plane);
          m_EllipFuncList[gnum].push_back(inside);
        }
      }
    }
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void EstablishFoamMorphology::assign_voxels()
{

  DataContainer::Pointer m = getDataContainerArray()->getDataContainer(getOutputCellAttributeMatrixPath().getDataContainerName());

  size_t totalPoints = m->getAttributeMatrix(m_OutputCellAttributeMatrixPath.getAttributeMatrixName())->getNumberOfTuples();

  SizeVec3Type udims = m->getGeometryAs<ImageGeom>()->getDimensions();

  const std::array<int64_t, 3> dims = {
      static_cast<int64_t>(udims[0]),
      static_cast<int64_t>(udims[1]),
      static_cast<int64_t>(udims[2]),
  };

  int64_t column = 0;
  int64_t row = 0;
  int64_t plane = 0;
  float xc = 0.0f;
  float yc = 0.0f;
  float zc = 0.0f;
  std::array<float, 3> size = {m_SizeX, m_SizeY, m_SizeZ};

  int64_t xmin = 0;
  int64_t xmax = 0;
  int64_t ymin = 0;
  int64_t ymax = 0;
  int64_t zmin = 0;
  int64_t zmax = 0;

  FloatVec3Type res = m->getGeometryAs<ImageGeom>()->getSpacing();

  Int32ArrayType::Pointer newownersPtr = Int32ArrayType::CreateArray(totalPoints, std::string("_INTERNAL_USE_ONLY_newowners"), true);
  newownersPtr->initializeWithValue(-1);
  int32_t* newowners = newownersPtr->getPointer(0);

  FloatArrayType::Pointer ellipfuncsPtr = FloatArrayType::CreateArray(totalPoints, std::string("_INTERNAL_USE_ONLY_ellipfuncs"), true);
  ellipfuncsPtr->initializeWithValue(-1);
  float* ellipfuncs = ellipfuncsPtr->getPointer(0);

  float featuresPerTime = 0;
  uint64_t millis = QDateTime::currentMSecsSinceEpoch();
  uint64_t currentMillis = millis;

  int64_t totalFeatures = m->getAttributeMatrix(getOutputCellFeatureAttributeMatrixName())->getNumberOfTuples();
  for(int64_t i = m_FirstFoamFeature; i < totalFeatures; i++)
  {
    featuresPerTime++;
    currentMillis = QDateTime::currentMSecsSinceEpoch();
    if(currentMillis - millis > 1000)
    {
      float rate = featuresPerTime / ((float)(currentMillis - millis)) * 1000.0f;

      QString ss = QObject::tr("Assign Voxels & Gaps|| Features Checked: %1 || Features/Second: %2").arg(i).arg((int)rate);
      notifyStatusMessage(ss);
      featuresPerTime = 0;
      millis = QDateTime::currentMSecsSinceEpoch();
    }
    float volcur = m_Volumes[i];
    float bovera = m_AxisLengths[3 * i + 1];
    float covera = m_AxisLengths[3 * i + 2];
    float omega3 = m_Omega3s[i];
    xc = m_Centroids[3 * i];
    yc = m_Centroids[3 * i + 1];
    zc = m_Centroids[3 * i + 2];
    float radcur1 = 0.0f;
    // Unbounded Check for the size of shapeTypes. We assume a 1:1 with phase
    uint32_t shapeclass = m_ShapeTypes[m_FeaturePhases[i]];

    if(shapeclass != 0 && shapeclass != 1 && shapeclass != 2 && shapeclass != 3)
    {
      QString ss = QObject::tr("Undefined shape class in shape types array with path %1").arg(m_InputShapeTypesArrayPath.serialize());
      setErrorCondition(-666, ss);
      return;
    }

    // init any values for each of the Shape Ops
    for(auto& shapeOp : m_ShapeOps)
    {
      shapeOp->init();
    }
    // Create our Argument Map
    QMap<ShapeOps::ArgName, float> shapeArgMap;
    shapeArgMap[ShapeOps::Omega3] = omega3;
    shapeArgMap[ShapeOps::VolCur] = volcur;
    shapeArgMap[ShapeOps::B_OverA] = bovera;
    shapeArgMap[ShapeOps::C_OverA] = covera;

    radcur1 = m_ShapeOps[shapeclass]->radcur1(shapeArgMap);

    float radcur2 = (radcur1 * bovera);
    float radcur3 = (radcur1 * covera);
    float phi1 = m_AxisEulerAngles[3 * i];
    float PHI = m_AxisEulerAngles[3 * i + 1];
    float phi2 = m_AxisEulerAngles[3 * i + 2];
    float ga[3][3] = {{0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f}};
    OrientationTransformation::eu2om<OrientationF, OrientationF>(OrientationF(phi1, PHI, phi2)).toGMatrix(ga);
    column = static_cast<int64_t>(xc / res[0]);
    row = static_cast<int64_t>(yc / res[1]);
    plane = static_cast<int64_t>(zc / res[2]);
    xmin = static_cast<int64_t>(column - ((radcur1 / res[0]) + 1));
    xmax = static_cast<int64_t>(column + ((radcur1 / res[0]) + 1));
    ymin = static_cast<int64_t>(row - ((radcur1 / res[1]) + 1));
    ymax = static_cast<int64_t>(row + ((radcur1 / res[1]) + 1));
    zmin = static_cast<int64_t>(plane - ((radcur1 / res[2]) + 1));
    zmax = static_cast<int64_t>(plane + ((radcur1 / res[2]) + 1));

    if(m_PeriodicBoundaries)
    {
      if(xmin < -dims[0])
      {
        xmin = -dims[0];
      }
      if(xmax > 2 * dims[0] - 1)
      {
        xmax = (2 * dims[0] - 1);
      }
      if(ymin < -dims[1])
      {
        ymin = -dims[1];
      }
      if(ymax > 2 * dims[1] - 1)
      {
        ymax = (2 * dims[1] - 1);
      }
      if(zmin < -dims[2])
      {
        zmin = -dims[2];
      }
      if(zmax > 2 * dims[2] - 1)
      {
        zmax = (2 * dims[2] - 1);
      }
    }
    else
    {
      if(xmin < 0)
      {
        xmin = 0;
      }
      if(xmax > dims[0] - 1)
      {
        xmax = dims[0] - 1;
      }
      if(ymin < 0)
      {
        ymin = 0;
      }
      if(ymax > dims[1] - 1)
      {
        ymax = dims[1] - 1;
      }
      if(zmin < 0)
      {
        zmin = 0;
      }
      if(zmax > dims[2] - 1)
      {
        zmax = dims[2] - 1;
      }
    }

    std::array<float, 3> radCur = {radcur1, radcur2, radcur3};
    std::array<float, 3> xx = {xc, yc, zc};
    ShapeOps* shapeOps = m_ShapeOps[shapeclass].get();

#ifdef SIMPL_USE_PARALLEL_ALGORITHMS
    bool doParallel = true;

    if(doParallel)
    {
      tbb::parallel_for(tbb::blocked_range3d<int64_t, int64_t, int64_t>(zmin, zmax + 1, ymin, ymax + 1, xmin, xmax + 1),
                        FoamAssignVoxelsGapsImpl(udims, res, m_FeatureIds, radCur.data(), xx.data(), shapeOps, ga, size.data(), i, newownersPtr, ellipfuncsPtr), tbb::auto_partitioner());
    }
    else
#endif
    {
      FoamAssignVoxelsGapsImpl serial(udims, res, m_FeatureIds, radCur.data(), xx.data(), shapeOps, ga, size.data(), i, newownersPtr, ellipfuncsPtr);
      serial.convert(zmin, zmax + 1, ymin, ymax + 1, xmin, xmax + 1);
    }
  }

  QVector<bool> activeObjects(totalFeatures, false);
  int32_t gnum = 0;
  for(size_t i = 0; i < totalPoints; i++)
  {
    //    if(ellipfuncs[i] >= 0) { m_FeatureIds[i] = newowners[i]; }
    if(ellipfuncs[i] >= 0)
    {
      m_FeatureIds[i] = newowners[i];
      //	m_Mask[i] = false;
    }
    gnum = m_FeatureIds[i];
    if(gnum >= 0)
    {
      activeObjects[gnum] = true;
    }
    newowners[i] = -1;
    ellipfuncs[i] = -1.0f;
  }

  AttributeMatrix::Pointer cellFeatureAttrMat = m->getAttributeMatrix(getOutputCellFeatureAttributeMatrixName());
  cellFeatureAttrMat->removeInactiveObjects(activeObjects, m_FeatureIdsPtr.lock().get());
  // need to update pointers after removing inactive objects
  updateFeatureInstancePointers();
  totalFeatures = cellFeatureAttrMat->getNumberOfTuples();
  // counting the number of features for each phase
  for(int64_t i = 1; i < totalFeatures; i++)
  {
    int32_t phase = m_FeaturePhases[i];
    if(phase >= 0)
    {
      m_NumFeatures[phase]++;
    }
  }

  // need to update pointers after resize, but do not need to run full data check because pointers are still valid
  updateFeatureInstancePointers();
  if(getCancel())
  {
    return;
  }

  for(size_t i = 0; i < totalPoints; i++)
  {
    gnum = m_FeatureIds[i];
    if(gnum >= 0)
    {
      m_CellPhases[i] = m_FeaturePhases[gnum];
    }
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void EstablishFoamMorphology::assign_gaps_only()
{
  DataContainer::Pointer m = getDataContainerArray()->getDataContainer(getOutputCellAttributeMatrixPath().getDataContainerName());

  int32_t featurename = 0;
  int32_t feature = 0;
  int32_t current = 0;
  int32_t most = 0;
  int64_t gapVoxelCount = 1;
  int64_t previousGapVoxelCount = 0;
  int32_t iterationCounter = 0;
  int64_t neighpoint;
  bool good = false;
  int32_t neighbor = 0;

  int64_t xPoints = static_cast<int64_t>(m->getGeometryAs<ImageGeom>()->getXPoints());
  int64_t yPoints = static_cast<int64_t>(m->getGeometryAs<ImageGeom>()->getYPoints());
  int64_t zPoints = static_cast<int64_t>(m->getGeometryAs<ImageGeom>()->getZPoints());
  size_t totalPoints = m->getAttributeMatrix(m_OutputCellAttributeMatrixPath.getAttributeMatrixName())->getNumberOfTuples();
  size_t totalFeatures = m->getAttributeMatrix(getOutputCellFeatureAttributeMatrixName())->getNumberOfTuples();

  std::array<int64_t, 6> neighpoints = {0, 0, 0, 0, 0, 0};
  neighpoints[0] = -xPoints * yPoints;
  neighpoints[1] = -xPoints;
  neighpoints[2] = -1;
  neighpoints[3] = 1;
  neighpoints[4] = xPoints;
  neighpoints[5] = xPoints * yPoints;

  Int64ArrayType::Pointer neighborsPtr = Int64ArrayType::CreateArray(m->getGeometryAs<ImageGeom>()->getNumberOfElements(), std::string("_INTERNAL_USE_ONLY_Neighbors"), true);
  neighborsPtr->initializeWithValue(-1);
  m_Neighbors = neighborsPtr->getPointer(0);

  std::vector<int32_t> n(totalFeatures + 1, 0);

  while(gapVoxelCount != 0 && gapVoxelCount != previousGapVoxelCount)
  {
    iterationCounter++;
    previousGapVoxelCount = gapVoxelCount;
    gapVoxelCount = 0;
    int64_t zStride;
    int64_t yStride;
    for(int64_t i = 0; i < zPoints; i++)
    {
      zStride = i * xPoints * yPoints;
      for(int64_t j = 0; j < yPoints; j++)
      {
        yStride = j * xPoints;
        for(int64_t k = 0; k < xPoints; k++)
        {
          featurename = m_FeatureIds[zStride + yStride + k];
          if(featurename < 0)
          {
            gapVoxelCount++;
            current = 0;
            most = 0;
            for(int32_t l = 0; l < 6; l++)
            {
              good = true;
              neighpoint = zStride + yStride + k + neighpoints[l];
              if(l == 0 && i == 0)
              {
                good = false;
              }
              if(l == 5 && i == (zPoints - 1))
              {
                good = false;
              }
              if(l == 1 && j == 0)
              {
                good = false;
              }
              if(l == 4 && j == (yPoints - 1))
              {
                good = false;
              }
              if(l == 2 && k == 0)
              {
                good = false;
              }
              if(l == 3 && k == (xPoints - 1))
              {
                good = false;
              }
              if(good)
              {
                feature = m_FeatureIds[neighpoint];
                if(feature > 0)
                {
                  n[feature]++;
                  current = n[feature];
                  if(current > most)
                  {
                    most = current;
                    m_Neighbors[zStride + yStride + k] = neighpoint;
                  }
                }
              }
            }
            for(int32_t l = 0; l < 6; l++)
            {
              good = true;
              neighpoint = zStride + yStride + k + neighpoints[l];
              if(l == 0 && i == 0)
              {
                good = false;
              }
              if(l == 5 && i == (zPoints - 1))
              {
                good = false;
              }
              if(l == 1 && j == 0)
              {
                good = false;
              }
              if(l == 4 && j == (yPoints - 1))
              {
                good = false;
              }
              if(l == 2 && k == 0)
              {
                good = false;
              }
              if(l == 3 && k == (xPoints - 1))
              {
                good = false;
              }
              if(good)
              {
                feature = m_FeatureIds[neighpoint];
                if(feature > 0)
                {
                  n[feature] = 0;
                }
              }
            }
          }
        }
      }
    }
    for(size_t j = 0; j < totalPoints; j++)
    {
      featurename = m_FeatureIds[j];
      neighbor = m_Neighbors[j];
      if(featurename < 0 && neighbor != -1 && m_FeatureIds[neighbor] > 0)
      {
        m_FeatureIds[j] = m_FeatureIds[neighbor];
        m_CellPhases[j] = m_FeaturePhases[m_FeatureIds[neighbor]];
      }
    }
    if(iterationCounter >= 1)
    {
      QString ss = QObject::tr("Assign Gaps || Cycle#: %1 || Remaining Unassigned Voxel Count: %2").arg(iterationCounter).arg(gapVoxelCount);
      notifyStatusMessage(ss);
    }
    if(getCancel())
    {
      return;
    }
  }
  if(gapVoxelCount != 0)
  {
    for(size_t j = 0; j < totalPoints; j++)
    {
      if(m_FeatureIds[j] < 0)
      {
        m_FeatureIds[j] = 0;
      }
    }
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int32_t EstablishFoamMorphology::estimate_numfeatures(size_t xpoints, size_t ypoints, size_t zpoints, float xRes, float yRes, float zRes)
{
  float ptotalvol = 0.0f;
  int32_t phase = 0;

  ptotalvol = static_cast<float>((xpoints * xRes) * (ypoints * yRes) * (zpoints * zRes));
  if(ptotalvol == 0.0)
  {
    return 1;
  }

  // This is for convenience
  DataContainerArray::Pointer dca = getDataContainerArray();

  // DataContainer::Pointer m = getDataContainerArray()->getDataContainer(getOutputCellAttributeMatrixPath().getDataContainerName());
  std::vector<size_t> cDims(1, 1);
  m_PhaseTypesPtr = dca->getPrereqArrayFromPath<DataArray<uint32_t>>(this, getInputPhaseTypesArrayPath(), cDims);
  DataArray<uint32_t>* phaseType = m_PhaseTypesPtr.lock().get();

  StatsDataArray::Pointer statsPtr = dca->getPrereqArrayFromPath<StatsDataArray>(this, getInputStatsArrayPath(), cDims);
  m_StatsDataArray = std::dynamic_pointer_cast<StatsDataArray>(statsPtr);
  if(m_StatsDataArray.lock() == nullptr)
  {
    QString ss = QObject::tr("Stats Array Not Initialized correctly");
    setErrorCondition(-308, ss);
    return 1;
  }

#if 0

	IDataArray::Pointer iPtr = m->getAttributeMatrix(m_CellEnsembleAttributeMatrixName)->getAttributeArray(SIMPL::EnsembleData::PhaseTypes);
	// Get the PhaseTypes - Remember there is a Dummy PhaseType in the first slot of the array
	DataArray<uint32_t>* phaseType = DataArray<uint32_t>::SafePointerDownCast(iPtr.get());

	iPtr = m->getAttributeMatrix(m_CellEnsembleAttributeMatrixName)->getAttributeArray(SIMPL::EnsembleData::Statistics);
	StatsDataArray* statsDataArrayPtr = StatsDataArray::SafePointerDownCast(iPtr.get());
	if (nullptr == statsDataArrayPtr)
	{
		return 1;
	}
#endif
  // Create a Reference Variable so we can use the [] syntax
  StatsDataArray& statsDataArray = *(m_StatsDataArray.lock());

  SIMPL_RANDOMNG_NEW()

  std::vector<int32_t> precipitatePhasesLocal;
  std::vector<double> precipitatePhaseFractionsLocal;
  double totalprecipitatefractions = 0.0;
  // find which phases are precipitate phases
  size_t numPhases = phaseType->getNumberOfTuples();
  for(size_t i = 1; i < numPhases; ++i)
  {
    if(phaseType->getValue(i) == static_cast<PhaseType::EnumType>(PhaseType::Type::Precipitate))
    {
      PrecipitateStatsData::Pointer pp = std::dynamic_pointer_cast<PrecipitateStatsData>(statsDataArray[i]);
      precipitatePhasesLocal.push_back(i);
      precipitatePhaseFractionsLocal.push_back(pp->getPhaseFraction());
      totalprecipitatefractions = totalprecipitatefractions + pp->getPhaseFraction();
    }
  }

  // scale the precipitate phase fractions to total to 1
  for(auto& precipPhaseFraction : precipitatePhaseFractionsLocal)
  {
    precipPhaseFraction = precipPhaseFraction / totalprecipitatefractions;
  }

  // generate the features
  int32_t gid = 1;

  float currentvol = 0.0f;
  float vol = 0.0f;
  float diam = 0.0f;
  bool volgood = false;
  size_t numLocalPrecipitatePhases = precipitatePhasesLocal.size();

  for(size_t j = 0; j < numLocalPrecipitatePhases; ++j)
  {
    float curphasetotalvol = m_TotalVol * precipitatePhaseFractionsLocal[j];
    while(currentvol < curphasetotalvol)
    {
      volgood = false;
      phase = precipitatePhasesLocal[j];
      PrecipitateStatsData::Pointer pp = std::dynamic_pointer_cast<PrecipitateStatsData>(statsDataArray[phase]);
      while(!volgood)
      {
        volgood = true;
        if(pp->getFeatureSize_DistType() == SIMPL::DistributionType::LogNormal)
        {
          float avgdiam = pp->getFeatureSizeDistribution().at(0)->getValue(0);
          float sddiam = pp->getFeatureSizeDistribution().at(1)->getValue(0);
          diam = rg.genrand_norm(avgdiam, sddiam);
          diam = expf(diam);
          if(diam >= pp->getMaxFeatureDiameter())
          {
            volgood = false;
          }
          if(diam < pp->getMinFeatureDiameter())
          {
            volgood = false;
          }
          vol = (4.0f / 3.0f) * (M_PI) * ((diam * 0.5f) * (diam * 0.5f) * (diam * 0.5f));
        }
      }
      currentvol = currentvol + vol;
      gid++;
    }
  }
  return gid;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void EstablishFoamMorphology::write_goal_attributes()
{
  clearErrorCode();
  DataContainer::Pointer m = getDataContainerArray()->getDataContainer(getOutputCellAttributeMatrixPath().getDataContainerName());

  size_t totalFeatures = m->getAttributeMatrix(getOutputCellFeatureAttributeMatrixName())->getNumberOfTuples();

  // Make sure any directory path is also available as the user may have just typed
  // in a path without actually creating the full path
  QFileInfo fi(m_CsvOutputFile);
  QDir parentPath = fi.path();
  if(!parentPath.mkpath("."))
  {
    QString ss = QObject::tr("Error creating parent path '%1'").arg(parentPath.absolutePath());
    setErrorCondition(-1, ss);
    return;
  }

  QFile outFile(getCsvOutputFile());
  if(!outFile.open(QIODevice::WriteOnly))
  {
    QString msg = QObject::tr("CSV Output file could not be opened: %1").arg(getCsvOutputFile());
    setErrorCondition(-200, msg);
    return;
  }

  QTextStream dStream(&outFile);

  char space = ',';
  // Write the total number of features
  dStream << totalFeatures - m_FirstFoamFeature;
  // Get all the names of the arrays from the Data Container
  QList<QString> headers = m->getAttributeMatrix(getOutputCellFeatureAttributeMatrixName())->getAttributeArrayNames();

  QVector<IDataArray::Pointer> data;

  // For checking if an array is a neighborlist
  NeighborList<int32_t>::Pointer neighborlistPtr = NeighborList<int32_t>::New();

  // Print the FeatureIds Header before the rest of the headers
  dStream << SIMPL::FeatureData::FeatureID;
  // Loop throught the list and print the rest of the headers, ignoring those we don't want
  // for(QList<QString>::iterator iter = headers.begin(); iter != headers.end(); ++iter)
  for(const auto& header : headers)
  {
    // Only get the array if the name does NOT match those listed
    IDataArray::Pointer p = m->getAttributeMatrix(getOutputCellFeatureAttributeMatrixName())->getAttributeArray(header);
    if(p->getNameOfClass().compare(neighborlistPtr->getNameOfClass()) != 0)
    {
      if(p->getNumberOfComponents() == 1)
      {
        dStream << space << header;
      }
      else // There are more than a single component so we need to add multiple header values
      {
        for(int32_t k = 0; k < p->getNumberOfComponents(); ++k)
        {
          dStream << space << header << "_" << k;
        }
      }
      // Get the IDataArray from the DataContainer
      data.push_back(p);
    }
  }
  dStream << "\n";

  // Get the number of tuples in the arrays
  size_t numTuples = data[0]->getNumberOfTuples();

  float threshold = 0.0f;

  // Skip the first feature
  for(size_t i = m_FirstFoamFeature; i < numTuples; ++i)
  {
    if(((float)i / numTuples) * 100.0f > threshold)
    {

      QString ss = QObject::tr("Writing Feature Data || %1% Complete").arg(((float)i / numTuples) * 100);
      notifyStatusMessage(ss);
      threshold = threshold + 5.0f;
      if(threshold < ((float)i / numTuples) * 100.0f)
      {
        threshold = ((float)i / numTuples) * 100.0f;
      }
    }

    // Print the feature id
    dStream << i;
    // Print a row of data
    for(const auto& d : data)
    {
      dStream << space;
      d->printTuple(dStream, i, space);
    }
    dStream << "\n";
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void EstablishFoamMorphology::form_struts()
{
  DataContainer::Pointer m = getDataContainerArray()->getDataContainer(getOutputCellAttributeMatrixPath().getDataContainerName());

  SizeVec3Type udims = m->getGeometryAs<ImageGeom>()->getDimensions();

  size_t point = 0;
  size_t kstride = 0;
  size_t jstride = 0;

  float QPpoint = 0.0f;
  float TJpoint = 0.0f;
  float GBpoint = 0.0f;
  float QPTJpoint = 0.0f;
  float TJGBpoint = 0.0f;

  // index into the flat array of points
  for(size_t k = 0; k < udims[2]; k++)
  {
    kstride = udims[0] * udims[1] * k;
    for(size_t j = 0; j < udims[1]; j++)
    {
      jstride = udims[0] * j;
      for(size_t i = 0; i < udims[0]; i++)
      {
        point = kstride + jstride + i;
        QPpoint = m_QPEuclideanDistances[point];
        TJpoint = m_TJEuclideanDistances[point];
        GBpoint = m_GBEuclideanDistances[point];
        QPTJpoint = QPpoint * TJpoint;
        TJGBpoint = TJpoint * GBpoint;
        if(QPpoint == 0.0f || TJpoint <= m_MinStrutThickness || (QPTJpoint < m_StrutThicknessVariability && TJGBpoint < m_StrutShapeVariability))
        {
          m_FeatureIds[point] = 0;
          m_Mask[point] = true;
        }
      }
    }
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void EstablishFoamMorphology::find_euclidean_dist_map()
{
  if(getErrorCode() < 0)
  {
    return;
  }

  find_euclideandistmap();

  DataContainer::Pointer m = getDataContainerArray()->getDataContainer(getOutputCellAttributeMatrixPath().getDataContainerName());
  AttributeMatrix::Pointer attrMat = m->getAttributeMatrix(getOutputCellAttributeMatrixPath().getAttributeMatrixName());
  attrMat->removeAttributeArray(m_NearestNeighborsArrayName);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void EstablishFoamMorphology::find_euclideandistmap()
{
  DataContainer::Pointer m = getDataContainerArray()->getDataContainer(getOutputCellAttributeMatrixPath().getDataContainerName());

  size_t totalPoints = m_FeatureIdsPtr.lock()->getNumberOfTuples();

  for(size_t i = 0; i < totalPoints; i++)
  {
    m_GBEuclideanDistances[i] = -1;
    m_TJEuclideanDistances[i] = -1;
    m_QPEuclideanDistances[i] = -1;
  }

  int64_t column = 0;
  int64_t row = 0;
  int64_t plane = 0;
  bool good = false;
  bool add = true;
  int32_t feature = 0;
  std::vector<int32_t> coordination;

  SizeVec3Type udims = m->getGeometryAs<ImageGeom>()->getDimensions();

  std::array<int64_t, 3> dims = {
      static_cast<int64_t>(udims[0]),
      static_cast<int64_t>(udims[1]),
      static_cast<int64_t>(udims[2]),
  };

  int64_t neighbor = 0;
  std::array<int64_t, 6> neighbors = {0, 0, 0, 0, 0, 0};
  neighbors[0] = -dims[0] * dims[1];
  neighbors[1] = -dims[0];
  neighbors[2] = -1;
  neighbors[3] = 1;
  neighbors[4] = udims[0];
  neighbors[5] = udims[0] * udims[1];

  int64_t xPoints = static_cast<int64_t>(m->getGeometryAs<ImageGeom>()->getXPoints());
  int64_t yPoints = static_cast<int64_t>(m->getGeometryAs<ImageGeom>()->getYPoints());
  int64_t zPoints = static_cast<int64_t>(m->getGeometryAs<ImageGeom>()->getZPoints());

  for(size_t a = 0; a < totalPoints; ++a)
  {
    feature = m_FeatureIds[a];
    if(feature > 0)
    {
      column = static_cast<int64_t>(a % xPoints);
      row = static_cast<int64_t>((a / xPoints) % yPoints);
      plane = static_cast<int64_t>(a / (xPoints * yPoints));
      for(int32_t k = 0; k < 6; k++)
      {
        good = true;
        neighbor = static_cast<int64_t>(a + neighbors[k]);
        if(k == 0 && plane == 0)
        {
          good = false;
        }
        if(k == 5 && plane == (zPoints - 1))
        {
          good = false;
        }
        if(k == 1 && row == 0)
        {
          good = false;
        }
        if(k == 4 && row == (yPoints - 1))
        {
          good = false;
        }
        if(k == 2 && column == 0)
        {
          good = false;
        }
        if(k == 3 && column == (xPoints - 1))
        {
          good = false;
        }
        if(good && m_FeatureIds[neighbor] != feature && m_FeatureIds[neighbor] >= 0)
        {
          add = true;
          for(const auto& coordinationValue : coordination)
          {
            if(m_FeatureIds[neighbor] == coordinationValue)
            {
              add = false;
            }
          }
          if(add)
          {
            coordination.push_back(m_FeatureIds[neighbor]);
          }
        }
      }
      if(coordination.empty())
      {
        m_NearestNeighbors[a * 3 + 0] = -1, m_NearestNeighbors[a * 3 + 1] = -1, m_NearestNeighbors[a * 3 + 2] = -1;
      }
      if(!coordination.empty())
      {
        m_GBEuclideanDistances[a] = 0.0f, m_NearestNeighbors[a * 3 + 0] = coordination[0], m_NearestNeighbors[a * 3 + 1] = -1, m_NearestNeighbors[a * 3 + 2] = -1;
      }
      if(coordination.size() >= 2)
      {
        m_TJEuclideanDistances[a] = 0.0f, m_NearestNeighbors[a * 3 + 0] = coordination[0], m_NearestNeighbors[a * 3 + 1] = coordination[0], m_NearestNeighbors[a * 3 + 2] = -1;
      }
      if(coordination.size() > 2)
      {
        m_QPEuclideanDistances[a] = 0.0f, m_NearestNeighbors[a * 3 + 0] = coordination[0], m_NearestNeighbors[a * 3 + 1] = coordination[0], m_NearestNeighbors[a * 3 + 2] = coordination[0];
      }
      coordination.resize(0);
    }
  }

#ifdef SIMPL_USE_PARALLEL_ALGORITHMS
  bool doParallel = true;

  if(doParallel)
  {
    tbb::task_group* g = new tbb::task_group;
    g->run(FindEuclideanMap2(m, m_FeatureIds, m_NearestNeighbors, m_GBEuclideanDistances, m_TJEuclideanDistances, m_QPEuclideanDistances, 0));
    g->run(FindEuclideanMap2(m, m_FeatureIds, m_NearestNeighbors, m_GBEuclideanDistances, m_TJEuclideanDistances, m_QPEuclideanDistances, 1));
    g->run(FindEuclideanMap2(m, m_FeatureIds, m_NearestNeighbors, m_GBEuclideanDistances, m_TJEuclideanDistances, m_QPEuclideanDistances, 2));
    g->wait();
    delete g;
  }
  else
#endif
  {
    for(int32_t i = 0; i < 3; i++)
    {
      if((i == 0) || (i == 1) || (i == 2))
      {
        FindEuclideanMap2 f(m, m_FeatureIds, m_NearestNeighbors, m_GBEuclideanDistances, m_TJEuclideanDistances, m_QPEuclideanDistances, i);
        f();
      }
    }
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AbstractFilter::Pointer EstablishFoamMorphology::newFilterInstance(bool copyFilterParameters) const
{
  EstablishFoamMorphology::Pointer filter = EstablishFoamMorphology::New();
  if(copyFilterParameters)
  {
    copyFilterParameterInstanceVariables(filter.get());
  }
  return filter;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString EstablishFoamMorphology::getCompiledLibraryName() const
{
  return DREAM3DReviewConstants::DREAM3DReviewBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString EstablishFoamMorphology::getBrandingString() const
{
  return "Univ. of Utah DMREF";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString EstablishFoamMorphology::getFilterVersion() const
{
  QString version;
  QTextStream vStream(&version);
  vStream << DREAM3DReview::Version::Major() << "." << DREAM3DReview::Version::Minor() << "." << DREAM3DReview::Version::Patch();
  return version;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString EstablishFoamMorphology::getGroupName() const
{
  return SIMPL::FilterGroups::SyntheticBuildingFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QUuid EstablishFoamMorphology::getUuid() const
{
  return QUuid("{14ad5792-a186-5e96-8e9e-4d623919dabc}");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString EstablishFoamMorphology::getSubGroupName() const
{
  return SIMPL::FilterSubGroups::PackingFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString EstablishFoamMorphology::getHumanLabel() const
{
  return "Establish Foam Morphology";
}

// -----------------------------------------------------------------------------
EstablishFoamMorphology::Pointer EstablishFoamMorphology::NullPointer()
{
  return Pointer(static_cast<Self*>(nullptr));
}

// -----------------------------------------------------------------------------
EstablishFoamMorphology::Pointer EstablishFoamMorphology::New()
{
  Pointer sharedPtr(new(EstablishFoamMorphology));
  return sharedPtr;
}

// -----------------------------------------------------------------------------
QString EstablishFoamMorphology::getNameOfClass() const
{
  return QString("EstablishFoamMorphology");
}

// -----------------------------------------------------------------------------
QString EstablishFoamMorphology::ClassName()
{
  return QString("EstablishFoamMorphology");
}

// -----------------------------------------------------------------------------
void EstablishFoamMorphology::setDataContainerName(const DataArrayPath& value)
{
  m_DataContainerName = value;
}

// -----------------------------------------------------------------------------
DataArrayPath EstablishFoamMorphology::getDataContainerName() const
{
  return m_DataContainerName;
}

// -----------------------------------------------------------------------------
void EstablishFoamMorphology::setOutputCellAttributeMatrixPath(const DataArrayPath& value)
{
  m_OutputCellAttributeMatrixPath = value;
}

// -----------------------------------------------------------------------------
DataArrayPath EstablishFoamMorphology::getOutputCellAttributeMatrixPath() const
{
  return m_OutputCellAttributeMatrixPath;
}

// -----------------------------------------------------------------------------
void EstablishFoamMorphology::setInputCellFeatureIdsArrayPath(const DataArrayPath& value)
{
  m_InputCellFeatureIdsArrayPath = value;
}

// -----------------------------------------------------------------------------
DataArrayPath EstablishFoamMorphology::getInputCellFeatureIdsArrayPath() const
{
  return m_InputCellFeatureIdsArrayPath;
}

// -----------------------------------------------------------------------------
void EstablishFoamMorphology::setOutputCellEnsembleAttributeMatrixName(const QString& value)
{
  m_OutputCellEnsembleAttributeMatrixName = value;
}

// -----------------------------------------------------------------------------
QString EstablishFoamMorphology::getOutputCellEnsembleAttributeMatrixName() const
{
  return m_OutputCellEnsembleAttributeMatrixName;
}

// -----------------------------------------------------------------------------
void EstablishFoamMorphology::setOutputCellFeatureAttributeMatrixName(const QString& value)
{
  m_OutputCellFeatureAttributeMatrixName = value;
}

// -----------------------------------------------------------------------------
QString EstablishFoamMorphology::getOutputCellFeatureAttributeMatrixName() const
{
  return m_OutputCellFeatureAttributeMatrixName;
}

// -----------------------------------------------------------------------------
void EstablishFoamMorphology::setFeatureIdsArrayName(const QString& value)
{
  m_FeatureIdsArrayName = value;
}

// -----------------------------------------------------------------------------
QString EstablishFoamMorphology::getFeatureIdsArrayName() const
{
  return m_FeatureIdsArrayName;
}

// -----------------------------------------------------------------------------
void EstablishFoamMorphology::setMaskArrayName(const QString& value)
{
  m_MaskArrayName = value;
}

// -----------------------------------------------------------------------------
QString EstablishFoamMorphology::getMaskArrayName() const
{
  return m_MaskArrayName;
}

// -----------------------------------------------------------------------------
void EstablishFoamMorphology::setCellPhasesArrayName(const QString& value)
{
  m_CellPhasesArrayName = value;
}

// -----------------------------------------------------------------------------
QString EstablishFoamMorphology::getCellPhasesArrayName() const
{
  return m_CellPhasesArrayName;
}

// -----------------------------------------------------------------------------
void EstablishFoamMorphology::setNumFeaturesArrayName(const QString& value)
{
  m_NumFeaturesArrayName = value;
}

// -----------------------------------------------------------------------------
QString EstablishFoamMorphology::getNumFeaturesArrayName() const
{
  return m_NumFeaturesArrayName;
}

// -----------------------------------------------------------------------------
void EstablishFoamMorphology::setQPEuclideanDistancesArrayName(const QString& value)
{
  m_QPEuclideanDistancesArrayName = value;
}

// -----------------------------------------------------------------------------
QString EstablishFoamMorphology::getQPEuclideanDistancesArrayName() const
{
  return m_QPEuclideanDistancesArrayName;
}

// -----------------------------------------------------------------------------
void EstablishFoamMorphology::setTJEuclideanDistancesArrayName(const QString& value)
{
  m_TJEuclideanDistancesArrayName = value;
}

// -----------------------------------------------------------------------------
QString EstablishFoamMorphology::getTJEuclideanDistancesArrayName() const
{
  return m_TJEuclideanDistancesArrayName;
}

// -----------------------------------------------------------------------------
void EstablishFoamMorphology::setGBEuclideanDistancesArrayName(const QString& value)
{
  m_GBEuclideanDistancesArrayName = value;
}

// -----------------------------------------------------------------------------
QString EstablishFoamMorphology::getGBEuclideanDistancesArrayName() const
{
  return m_GBEuclideanDistancesArrayName;
}

// -----------------------------------------------------------------------------
void EstablishFoamMorphology::setInputStatsArrayPath(const DataArrayPath& value)
{
  m_InputStatsArrayPath = value;
}

// -----------------------------------------------------------------------------
DataArrayPath EstablishFoamMorphology::getInputStatsArrayPath() const
{
  return m_InputStatsArrayPath;
}

// -----------------------------------------------------------------------------
void EstablishFoamMorphology::setInputPhaseTypesArrayPath(const DataArrayPath& value)
{
  m_InputPhaseTypesArrayPath = value;
}

// -----------------------------------------------------------------------------
DataArrayPath EstablishFoamMorphology::getInputPhaseTypesArrayPath() const
{
  return m_InputPhaseTypesArrayPath;
}

// -----------------------------------------------------------------------------
void EstablishFoamMorphology::setInputPhaseNamesArrayPath(const DataArrayPath& value)
{
  m_InputPhaseNamesArrayPath = value;
}

// -----------------------------------------------------------------------------
DataArrayPath EstablishFoamMorphology::getInputPhaseNamesArrayPath() const
{
  return m_InputPhaseNamesArrayPath;
}

// -----------------------------------------------------------------------------
void EstablishFoamMorphology::setInputShapeTypesArrayPath(const DataArrayPath& value)
{
  m_InputShapeTypesArrayPath = value;
}

// -----------------------------------------------------------------------------
DataArrayPath EstablishFoamMorphology::getInputShapeTypesArrayPath() const
{
  return m_InputShapeTypesArrayPath;
}

// -----------------------------------------------------------------------------
void EstablishFoamMorphology::setCsvOutputFile(const QString& value)
{
  m_CsvOutputFile = value;
}

// -----------------------------------------------------------------------------
QString EstablishFoamMorphology::getCsvOutputFile() const
{
  return m_CsvOutputFile;
}

// -----------------------------------------------------------------------------
void EstablishFoamMorphology::setPeriodicBoundaries(const bool& value)
{
  m_PeriodicBoundaries = value;
}

// -----------------------------------------------------------------------------
bool EstablishFoamMorphology::getPeriodicBoundaries() const
{
  return m_PeriodicBoundaries;
}

// -----------------------------------------------------------------------------
void EstablishFoamMorphology::setWriteGoalAttributes(const bool& value)
{
  m_WriteGoalAttributes = value;
}

// -----------------------------------------------------------------------------
bool EstablishFoamMorphology::getWriteGoalAttributes() const
{
  return m_WriteGoalAttributes;
}

// -----------------------------------------------------------------------------
void EstablishFoamMorphology::setHaveFeatures(const int& value)
{
  m_HaveFeatures = value;
}

// -----------------------------------------------------------------------------
int EstablishFoamMorphology::getHaveFeatures() const
{
  return m_HaveFeatures;
}

// -----------------------------------------------------------------------------
void EstablishFoamMorphology::setMinStrutThickness(const double& value)
{
  m_MinStrutThickness = value;
}

// -----------------------------------------------------------------------------
double EstablishFoamMorphology::getMinStrutThickness() const
{
  return m_MinStrutThickness;
}

// -----------------------------------------------------------------------------
void EstablishFoamMorphology::setStrutThicknessVariability(const float& value)
{
  m_StrutThicknessVariability = value;
}

// -----------------------------------------------------------------------------
float EstablishFoamMorphology::getStrutThicknessVariability() const
{
  return m_StrutThicknessVariability;
}

// -----------------------------------------------------------------------------
void EstablishFoamMorphology::setStrutShapeVariability(const float& value)
{
  m_StrutShapeVariability = value;
}

// -----------------------------------------------------------------------------
float EstablishFoamMorphology::getStrutShapeVariability() const
{
  return m_StrutShapeVariability;
}

// -----------------------------------------------------------------------------
void EstablishFoamMorphology::setSmoothStruts(const bool& value)
{
  m_SmoothStruts = value;
}

// -----------------------------------------------------------------------------
bool EstablishFoamMorphology::getSmoothStruts() const
{
  return m_SmoothStruts;
}
