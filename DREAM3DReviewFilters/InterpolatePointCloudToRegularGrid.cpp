/* ============================================================================
 * Software developed by US federal government employees (including military personnel)
 * as part of their official duties is not subject to copyright protection and is
 * considered "public domain" (see 17 USC Section 105). Public domain software can be used
 * by anyone for any purpose, and cannot be released under a copyright license
 * (including typical open source software licenses).
 *
 * This source code file was originally written by United States DoD employees. The
 * original source code files are released into the Public Domain.
 *
 * Subsequent changes to the codes by others may elect to add a copyright and license
 * for those changes.
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include "InterpolatePointCloudToRegularGrid.h"

#include <memory>

#include <QtCore/QTextStream>

#include "SIMPLib/Common/Constants.h"
#include "SIMPLib/Common/TemplateHelpers.h"
#include "SIMPLib/FilterParameters/AbstractFilterParametersReader.h"
#include "SIMPLib/FilterParameters/DataArraySelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/DataContainerSelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedBooleanFilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedChoicesFilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedPathCreationFilterParameter.h"
#include "SIMPLib/FilterParameters/MultiDataArraySelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/SeparatorFilterParameter.h"
#include "SIMPLib/FilterParameters/StringFilterParameter.h"
#include "SIMPLib/Geometry/VertexGeom.h"
#include "SIMPLib/Utilities/TimeUtilities.h"
#include "SIMPLib/DataContainers/DataContainerArray.h"
#include "SIMPLib/DataArrays/IDataArray.h"
#include "SIMPLib/DataContainers/DataContainer.h"
#include "SIMPLib/ITK/Dream3DTemplateAliasMacro.h"

#include "DREAM3DReview/DREAM3DReviewConstants.h"
#include "DREAM3DReview/DREAM3DReviewVersion.h"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
InterpolatePointCloudToRegularGrid::InterpolatePointCloudToRegularGrid()
: m_DataContainerName("")
, m_ArraysToInterpolate(QVector<DataArrayPath>())
, m_ArraysToCopy(QVector<DataArrayPath>())
, m_VoxelIndicesArrayPath("", "", "VoxelIndices")
, m_InterpolatedDataContainerName("InterpolatedDataContainer")
, m_InterpolatedAttributeMatrixName("InterpolatedAttributeMatrix")
, m_MaskArrayPath("", "", "")
, m_KernelDistancesArrayName("KernelDistances")
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
InterpolatePointCloudToRegularGrid::~InterpolatePointCloudToRegularGrid() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void InterpolatePointCloudToRegularGrid::setupFilterParameters()
{
  FilterParameterVectorType parameters;
  QStringList linkedProps("MaskArrayPath");
  parameters.push_back(SIMPL_NEW_LINKED_BOOL_FP("Use Mask", UseMask, FilterParameter::Parameter, InterpolatePointCloudToRegularGrid, linkedProps));
  linkedProps.clear();
  linkedProps << "KernelDistancesArrayName";
  parameters.push_back(SIMPL_NEW_LINKED_BOOL_FP("Store Kernel Distances", StoreKernelDistances, FilterParameter::Parameter, InterpolatePointCloudToRegularGrid, linkedProps));
  {
    QVector<QString> choices;
    choices.push_back("Uniform");
    choices.push_back("Gaussian");
    QStringList linkedProps;
    linkedProps << "Sigmas";
    LinkedChoicesFilterParameter::Pointer parameter = LinkedChoicesFilterParameter::New();
    parameter->setHumanLabel("Interpolation Technique");
    parameter->setPropertyName("InterpolationTechnique");
    parameter->setSetterCallback(SIMPL_BIND_SETTER(InterpolatePointCloudToRegularGrid, this, InterpolationTechnique));
    parameter->setGetterCallback(SIMPL_BIND_GETTER(InterpolatePointCloudToRegularGrid, this, InterpolationTechnique));
    parameter->setChoices(choices);
    parameter->setLinkedProperties(linkedProps);
    parameter->setCategory(FilterParameter::Parameter);
    parameters.push_back(parameter);
  }
  parameters.push_back(SIMPL_NEW_FLOAT_VEC3_FP("Kernel Size", KernelSize, FilterParameter::Parameter, InterpolatePointCloudToRegularGrid));

  parameters.push_back(SIMPL_NEW_FLOAT_VEC3_FP("Gaussian Sigmas", Sigmas, FilterParameter::Parameter, InterpolatePointCloudToRegularGrid, 1));
  {
    DataContainerSelectionFilterParameter::RequirementType req;
    IGeometry::Types reqGeom = {IGeometry::Type::Vertex};
    req.dcGeometryTypes = reqGeom;
    parameters.push_back(SIMPL_NEW_DC_SELECTION_FP("Data Container to Interpolate", DataContainerName, FilterParameter::RequiredArray, InterpolatePointCloudToRegularGrid, req));
  }
  {
    DataContainerSelectionFilterParameter::RequirementType req;
    IGeometry::Types reqGeom = {IGeometry::Type::Image};
    req.dcGeometryTypes = reqGeom;
    parameters.push_back(SIMPL_NEW_DC_SELECTION_FP("Interpolated Data Container", InterpolatedDataContainerName, FilterParameter::RequiredArray, InterpolatePointCloudToRegularGrid, req));
  }
  parameters.push_back(SeparatorFilterParameter::New("Vertex Data", FilterParameter::RequiredArray));
  {
    DataArraySelectionFilterParameter::RequirementType req = DataArraySelectionFilterParameter::CreateRequirement(SIMPL::TypeNames::Int64, 1, AttributeMatrix::Type::Vertex, IGeometry::Type::Vertex);
    parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Voxel Indices", VoxelIndicesArrayPath, FilterParameter::RequiredArray, InterpolatePointCloudToRegularGrid, req));
  }
  {
    DataArraySelectionFilterParameter::RequirementType req = DataArraySelectionFilterParameter::CreateRequirement(SIMPL::TypeNames::Bool, 1, AttributeMatrix::Type::Vertex, IGeometry::Type::Vertex);
    parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Mask", MaskArrayPath, FilterParameter::RequiredArray, InterpolatePointCloudToRegularGrid, req));
  }
  {
    MultiDataArraySelectionFilterParameter::RequirementType req =
        MultiDataArraySelectionFilterParameter::CreateRequirement(SIMPL::Defaults::AnyPrimitive, SIMPL::Defaults::AnyComponentSize, AttributeMatrix::Type::Vertex, IGeometry::Type::Vertex);
    parameters.push_back(SIMPL_NEW_MDA_SELECTION_FP("Attribute Arrays to Interpolate", ArraysToInterpolate, FilterParameter::RequiredArray, InterpolatePointCloudToRegularGrid, req));
  }
  {
    MultiDataArraySelectionFilterParameter::RequirementType req =
        MultiDataArraySelectionFilterParameter::CreateRequirement(SIMPL::Defaults::AnyPrimitive, SIMPL::Defaults::AnyComponentSize, AttributeMatrix::Type::Vertex, IGeometry::Type::Vertex);
    parameters.push_back(SIMPL_NEW_MDA_SELECTION_FP("Attribute Arrays to Copy", ArraysToCopy, FilterParameter::RequiredArray, InterpolatePointCloudToRegularGrid, req));
  }
  parameters.push_back(SeparatorFilterParameter::New("Cell Data", FilterParameter::CreatedArray));
  parameters.push_back(SIMPL_NEW_AM_WITH_LINKED_DC_FP("Interpolated Attribute Matrix", InterpolatedAttributeMatrixName, InterpolatedDataContainerName, FilterParameter::CreatedArray, InterpolatePointCloudToRegularGrid));
  parameters.push_back(SIMPL_NEW_DA_WITH_LINKED_AM_FP("Kernel Distances", KernelDistancesArrayName, InterpolatedDataContainerName, InterpolatedAttributeMatrixName, FilterParameter::CreatedArray, InterpolatePointCloudToRegularGrid));
  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void InterpolatePointCloudToRegularGrid::readFilterParameters(AbstractFilterParametersReader* reader, int index)
{
  reader->openFilterGroup(this, index);
  setDataContainerName(reader->readDataArrayPath("DataContainerName", getDataContainerName()));
  setArraysToInterpolate(reader->readDataArrayPathVector("ArraysToInterpolate", getArraysToInterpolate()));
  setArraysToCopy(reader->readDataArrayPathVector("ArraysToCopy", getArraysToCopy()));
  setVoxelIndicesArrayPath(reader->readDataArrayPath("VoxelIndicesArrayPath", getVoxelIndicesArrayPath()));
  setInterpolatedDataContainerName(reader->readDataArrayPath("InterpolatedDataContainerName", getInterpolatedDataContainerName()));
  setInterpolatedAttributeMatrixName(reader->readString("InterpolatedAttributeMatrixName", getInterpolatedAttributeMatrixName()));
  setInterpolationTechnique(reader->readValue("InterpolationTechnique", getInterpolationTechnique()));
  setKernelSize(reader->readFloatVec3("KernelSize", getKernelSize()));
  setSigmas(reader->readFloatVec3("Sigmas", getSigmas()));
  setUseMask(reader->readValue("UseMask", getUseMask()));
  setMaskArrayPath(reader->readDataArrayPath("MaskArrayPath", getMaskArrayPath()));
  setStoreKernelDistances(reader->readValue("StoreKernelDistances", getStoreKernelDistances()));
  setKernelDistancesArrayName(reader->readString("KernelDistancesArrayName", getKernelDistancesArrayName()));
  reader->closeFilterGroup();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
template <typename T>
void createCompatibleNeighborList(AbstractFilter* filter, DataArrayPath path, std::vector<size_t> cDims, std::vector<IDataArray::WeakPointer>& dynamicArrays)
{
  IDataArray::WeakPointer ptr = filter->getDataContainerArray()->createNonPrereqArrayFromPath<NeighborList<T>>(filter, path, 0, cDims);
  dynamicArrays.push_back(ptr);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void InterpolatePointCloudToRegularGrid::initialize()
{
  m_KernelDistances = NeighborList<float>::NullPointer();
  m_AttrMatList.clear();
  m_SourceArraysToInterpolate.clear();
  m_SourceArraysToCopy.clear();
  m_DynamicArraysToInterpolate.clear();
  m_DynamicArraysToCopy.clear();
  m_Kernel.clear();
  m_KernelValDistances.clear();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void InterpolatePointCloudToRegularGrid::dataCheck()
{

  clearErrorCode();
  clearWarningCode();
  initialize();

  QVector<IDataArray::Pointer> dataArrays;

  VertexGeom::Pointer vertices = getDataContainerArray()->getPrereqGeometryFromDataContainer<VertexGeom>(this, getDataContainerName());
  ImageGeom::Pointer image = getDataContainerArray()->getPrereqGeometryFromDataContainer<ImageGeom>(this, getInterpolatedDataContainerName());

  if(getErrorCode() < 0)
  {
    return;
  }

  dataArrays.push_back(vertices->getVertices());

  if(getInterpolationTechnique() < 0 || getInterpolationTechnique() > 1)
  {
    QString ss = QObject::tr("Invalid selection for interpolation technique");
    setErrorCondition(-11000, ss);
  }

  if(getKernelSize()[0] <= 0 || getKernelSize()[1] <= 0 || getKernelSize()[2] <= 0)
  {
    QString ss = QObject::tr("All kernel dimensions must be positive.\n "
                             "Current kernel dimensions:\n x = %1\n y = %2\n z = %3\n")
                     .arg(getKernelSize()[0])
                     .arg(getKernelSize()[1])
                     .arg(getKernelSize()[2]);
    setErrorCondition(-11000, ss);
  }

  if(getSigmas()[0] <= 0 || getSigmas()[1] <= 0 || getSigmas()[2] <= 0)
  {
    QString ss = QObject::tr("All sigmas must be positive.\n "
                             "Current sigmas:\n x = %1\n y = %2\n z = %3\n")
                     .arg(getKernelSize()[0])
                     .arg(getKernelSize()[1])
                     .arg(getKernelSize()[2]);
    setErrorCondition(-11000, ss);
  }

  // Get the list of all attribute matrices in the input data container and add them to the regular grid data container
  DataContainer::Pointer m = getDataContainerArray()->getDataContainer(getDataContainerName());
  DataContainer::Pointer interpolatedDC = getDataContainerArray()->getDataContainer(getInterpolatedDataContainerName());
  m_AttrMatList = m->getAttributeMatrixNames();
  SizeVec3Type dims = image->getDimensions();
  std::vector<size_t> tDims = {dims[0], dims[1], dims[2]};
  QList<QString> tempDataArrayList;
  DataArrayPath tempPath;
  QString tempAttrMatName;

  // All vertex data in the original point cloud will be interpolated to the regular grid's cells
  // Feature/Ensemble attribute matrices will remain unchanged and are deep copied to the new data container below
  // Create the attribute matrix where all the interpolated data will be stored
  interpolatedDC->createNonPrereqAttributeMatrix(this, getInterpolatedAttributeMatrixName(), tDims, AttributeMatrix::Type::Cell);

  QVector<DataArrayPath> interpolatePaths = getArraysToInterpolate();
  QVector<DataArrayPath> copyPaths = getArraysToCopy();

  if(!DataArrayPath::ValidateVector(interpolatePaths) || !DataArrayPath::ValidateVector(copyPaths))
  {
    QString ss = QObject::tr("There are Attribute Arrays selected that are not contained in the same Attribute Matrix. All selected Attribute Arrays must belong to the same Attribute Matrix");
    setErrorCondition(-11002, ss);
  }

  DataArrayPath path(getDataContainerName().getDataContainerName(), "", "");

  if(!interpolatePaths.isEmpty())
  {
    if(!interpolatePaths[0].hasSameDataContainer(path))
    {
      QString ss = QObject::tr("Attribute Arrays selected for interpolating must belong to the same Data Container selected for interpolation");
      setErrorCondition(-11002, ss);
    }
  }

  if(!copyPaths.isEmpty())
  {
    if(!copyPaths[0].hasSameDataContainer(path))
    {
      QString ss = QObject::tr("Attribute Arrays selected for copying must belong to the same Data Container selected for interpolation");
      setErrorCondition(-11002, ss);
    }
  }

  if(!interpolatePaths.isEmpty() && !copyPaths.isEmpty())
  {
    if(!interpolatePaths[0].hasSameAttributeMatrixPath(copyPaths[0]))
    {
      QString ss = QObject::tr("Attribute Arrays selected for interpolation and copying must belong to the same Data Container and Attribute Matrix");
      setErrorCondition(-11002, ss);
    }
  }

  if(getErrorCode() < 0)
  {
    return;
  }

  QString attrMatName;

  if(attrMatName.isEmpty() && !interpolatePaths.isEmpty())
  {
    attrMatName = interpolatePaths[0].getAttributeMatrixName();
  }

  if(attrMatName.isEmpty() && !copyPaths.isEmpty())
  {
    attrMatName = copyPaths[0].getAttributeMatrixName();
  }

  // Loop through all the attribute matrices in the original data container
  // If we are in a vertex attribute matrix, create data arrays for all in the new interpolated data attribute matrix
  // Else, we are in a feature/ensemble attribute matrix, and just deep copy it into the new data container

  std::vector<size_t> cDims(1, 1);

  if(!getDataContainerName().isEmpty())
  {
    for(QList<QString>::Iterator it = m_AttrMatList.begin(); it != m_AttrMatList.end(); ++it)
    {
      AttributeMatrix::Pointer tmpAttrMat = m->getPrereqAttributeMatrix(this, *it, -301);
      if(getErrorCode() >= 0)
      {
        tempAttrMatName = tmpAttrMat->getName();
        if(tempAttrMatName != attrMatName)
        {
          AttributeMatrix::Pointer attrMat = tmpAttrMat->deepCopy();
          interpolatedDC->addOrReplaceAttributeMatrix(attrMat);
        }
        else
        {
          for(int32_t i = 0; i < interpolatePaths.size(); i++)
          {
            tempPath.update(getInterpolatedDataContainerName().getDataContainerName(), getInterpolatedAttributeMatrixName(), interpolatePaths[i].getDataArrayName() + "Interpolation");
            IDataArray::Pointer tmpDataArray = tmpAttrMat->getPrereqIDataArray(this, interpolatePaths[i].getDataArrayName(), -90002);
            if(getErrorCode() >= 0)
            {
              m_SourceArraysToInterpolate.push_back(tmpDataArray);
              dataArrays.push_back(tmpDataArray);
              std::vector<size_t> tmpDims = tmpDataArray->getComponentDimensions();
              if(tmpDims != cDims)
              {
                QString ss = QObject::tr("Attribute Arrays selected for interpolation must be scalar arrays");
                setErrorCondition(-11002, ss);
                return;
              }
              EXECUTE_FUNCTION_TEMPLATE_NO_BOOL(this, createCompatibleNeighborList, tmpDataArray, this, tempPath, cDims, m_DynamicArraysToInterpolate)
            }
          }

          for(int32_t i = 0; i < copyPaths.size(); i++)
          {
            tempPath.update(getInterpolatedDataContainerName().getDataContainerName(), getInterpolatedAttributeMatrixName(), copyPaths[i].getDataArrayName() + "Copy");
            IDataArray::Pointer tmpDataArray = tmpAttrMat->getPrereqIDataArray(this, copyPaths[i].getDataArrayName(), -90002);
            if(getErrorCode() >= 0)
            {
              m_SourceArraysToCopy.push_back(tmpDataArray);
              dataArrays.push_back(tmpDataArray);
              std::vector<size_t> tmpDims = tmpDataArray->getComponentDimensions();
              if(tmpDims != cDims)
              {
                QString ss = QObject::tr("Attribute Arrays selected for copying must be scalar arrays");
                setErrorCondition(-11002, ss);
                return;
              }
              EXECUTE_FUNCTION_TEMPLATE_NO_BOOL(this, createCompatibleNeighborList, tmpDataArray, this, tempPath, cDims, m_DynamicArraysToCopy)
            }
          }
        }
      }
    }
  }

  m_VoxelIndicesPtr = getDataContainerArray()->getPrereqArrayFromPath<DataArray<size_t>>(this, getVoxelIndicesArrayPath(),
                                                                                                         cDims);  /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
  if(nullptr != m_VoxelIndicesPtr.lock().get())                                                                   /* Validate the Weak Pointer wraps a non-nullptr pointer to a DataArray<T> object */
  {
    m_VoxelIndices = m_VoxelIndicesPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */
  if(getErrorCode() >= 0)
  {
    dataArrays.push_back(m_VoxelIndicesPtr.lock());
  }

  if(getUseMask())
  {
    m_MaskPtr =
        getDataContainerArray()->getPrereqArrayFromPath<DataArray<bool>>(this, getMaskArrayPath(), cDims); /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
    if(nullptr != m_MaskPtr.lock().get()) /* Validate the Weak Pointer wraps a non-nullptr pointer to a DataArray<T> object */
    {
      m_Mask = m_MaskPtr.lock()->getPointer(0);
    } /* Now assign the raw pointer to data from the DataArray<T> object */
    if(getErrorCode() >= 0)
    {
      dataArrays.push_back(m_MaskPtr.lock());
    }
  }

  path.update(getInterpolatedDataContainerName().getDataContainerName(), getInterpolatedAttributeMatrixName(), getKernelDistancesArrayName());

  if(getStoreKernelDistances())
  {
    m_KernelDistances = getDataContainerArray()->createNonPrereqArrayFromPath<NeighborList<float>>(this, path, 0, cDims);
  }

  getDataContainerArray()->validateNumberOfTuples(this, dataArrays);
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
template <typename T>
void mapPointCloudDataByKernel(IDataArray::Pointer source, IDataArray::Pointer dynamic, std::vector<float>& kernelVals, int64_t kernel[3], size_t dims[3], size_t curX, size_t curY, size_t curZ,
                               size_t vertIdx)
{
  typename DataArray<T>::Pointer inputDataPtr = std::dynamic_pointer_cast<DataArray<T>>(source);
  T* inputData = static_cast<T*>(inputDataPtr->getPointer(0));
  typename NeighborList<T>::Pointer interpolatedDataPtr = std::dynamic_pointer_cast<NeighborList<T>>(dynamic);

  size_t index = 0;
  int64_t startKernel[3] = {0, 0, 0};
  int64_t endKernel[3] = {0, 0, 0};
  size_t counter = 0;

  kernel[0] > int64_t(curX) ? startKernel[0] = 0 : startKernel[0] = curX - kernel[0];
  kernel[1] > int64_t(curY) ? startKernel[1] = 0 : startKernel[1] = curY - kernel[1];
  kernel[2] > int64_t(curZ) ? startKernel[2] = 0 : startKernel[2] = curZ - kernel[2];

  int64_t(curX) + kernel[0] >= int64_t(dims[0]) ? endKernel[0] = dims[0] - 1 : endKernel[0] = curX + kernel[0];
  int64_t(curY) + kernel[1] >= int64_t(dims[1]) ? endKernel[1] = dims[1] - 1 : endKernel[1] = curY + kernel[1];
  endKernel[2] = int64_t(curZ);

  for(int64_t z = startKernel[2]; z <= endKernel[2]; z++)
  {
    for(int64_t y = startKernel[1]; y <= endKernel[1]; y++)
    {
      for(int64_t x = startKernel[0]; x <= endKernel[0]; x++)
      {
        if(kernelVals[counter] == 0.0f)
        {
          continue;
        }
        index = (z * dims[1] * dims[0]) + (y * dims[0]) + x;
        interpolatedDataPtr->addEntry(index, kernelVals[counter] * inputData[vertIdx]);
        counter++;
      }
    }
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void InterpolatePointCloudToRegularGrid::determineKernel(int64_t kernelNumVoxels[3])
{
  size_t counter = 0;

  for(int64_t z = -kernelNumVoxels[2]; z <= kernelNumVoxels[2]; z++)
  {
    for(int64_t y = -kernelNumVoxels[1]; y <= kernelNumVoxels[1]; y++)
    {
      for(int64_t x = -kernelNumVoxels[0]; x <= kernelNumVoxels[0]; x++)
      {
        if(m_InterpolationTechnique == 0)
        {
          m_Kernel[counter] = 1.0f;
        }
        else if(m_InterpolationTechnique == 1)
        {
          m_Kernel[counter] = expf(-((x * x) / (2 * m_Sigmas[0] * m_Sigmas[0]) + (y * y) / (2 * m_Sigmas[1] * m_Sigmas[1]) + (z * z) / (2 * m_Sigmas[2] * m_Sigmas[2])));
        }
        counter++;
      }
    }
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void InterpolatePointCloudToRegularGrid::determineKernelDistances(int64_t kernelNumVoxels[3], FloatVec3Type res)
{
  size_t counter = 0;

  for(int64_t z = -kernelNumVoxels[2]; z <= kernelNumVoxels[2]; z++)
  {
    for(int64_t y = -kernelNumVoxels[1]; y <= kernelNumVoxels[1]; y++)
    {
      for(int64_t x = -kernelNumVoxels[0]; x <= kernelNumVoxels[0]; x++)
      {
        m_KernelValDistances[counter] = (x * x * res[0] * res[0]) + (y * y * res[1] * res[1]) + (z * z * res[2] * res[2]);
        m_KernelValDistances[counter] = sqrtf(m_KernelValDistances[counter]);
        counter++;
      }
    }
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void InterpolatePointCloudToRegularGrid::mapKernelDistances(int64_t kernel[3], size_t dims[3], size_t curX, size_t curY, size_t curZ)
{
  size_t index = 0;
  int64_t startKernel[3] = {0, 0, 0};
  int64_t endKernel[3] = {0, 0, 0};
  size_t counter = 0;

  kernel[0] > int64_t(curX) ? startKernel[0] = 0 : startKernel[0] = curX - kernel[0];
  kernel[1] > int64_t(curY) ? startKernel[1] = 0 : startKernel[1] = curY - kernel[1];
  kernel[2] > int64_t(curZ) ? startKernel[2] = 0 : startKernel[2] = curZ - kernel[2];

  int64_t(curX) + kernel[0] >= int64_t(dims[0]) ? endKernel[0] = dims[0] - 1 : endKernel[0] = curX + kernel[0];
  int64_t(curY) + kernel[1] >= int64_t(dims[1]) ? endKernel[1] = dims[1] - 1 : endKernel[1] = curY + kernel[1];
  endKernel[2] = int64_t(curZ);

  for(int64_t z = startKernel[2]; z <= endKernel[2]; z++)
  {
    for(int64_t y = startKernel[1]; y <= endKernel[1]; y++)
    {
      for(int64_t x = startKernel[0]; x <= endKernel[0]; x++)
      {
        if(m_Kernel[counter] == 0.0f)
        {
          continue;
        }
        index = (z * dims[1] * dims[0]) + (y * dims[0]) + x;
        m_KernelDistances.lock()->addEntry(index, m_KernelValDistances[counter]);
        counter++;
      }
    }
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void InterpolatePointCloudToRegularGrid::execute()
{
  clearErrorCode();
  clearWarningCode();
  dataCheck();
  if(getErrorCode() < 0)
  {
    return;
  }

  DataContainer::Pointer m = getDataContainerArray()->getDataContainer(getDataContainerName());
  DataContainer::Pointer interpolatedDC = getDataContainerArray()->getDataContainer(getInterpolatedDataContainerName());

  VertexGeom::Pointer vertices = m->getGeometryAs<VertexGeom>();
  ImageGeom::Pointer image = interpolatedDC->getGeometryAs<ImageGeom>();
  SizeVec3Type dims = image->getDimensions();
  FloatVec3Type res = image->getSpacing();
  int64_t kernelNumVoxels[3] = {0, 0, 0};

  MeshIndexType numVerts = vertices->getNumberOfVertices();
  size_t index = 0;
  size_t x = 0;
  size_t y = 0;
  size_t z = 0;

  size_t maxImageIndex = ((dims[2] - 1) * dims[0] * dims[1]) + ((dims[1] - 1) * dims[0]) + (dims[0] - 1);

  kernelNumVoxels[0] = int64_t(ceil((m_KernelSize[0] / res[0]) * 0.5f));
  kernelNumVoxels[1] = int64_t(ceil((m_KernelSize[1] / res[1]) * 0.5f));
  kernelNumVoxels[2] = int64_t(ceil((m_KernelSize[2] / res[2]) * 0.5f));

  if(m_KernelSize[0] < res[0])
  {
    kernelNumVoxels[0] = 0;
  }
  if(m_KernelSize[1] < res[1])
  {
    kernelNumVoxels[1] = 0;
  }
  if(m_KernelSize[2] < res[2])
  {
    kernelNumVoxels[2] = 0;
  }

  int64_t kernelSize[3] = {1, 1, 1};
  int64_t totalKernel = 0;

  for(size_t i = 0; i < 3; i++)
  {
    kernelSize[i] *= (kernelNumVoxels[i] * 2) + 1;
  }

  totalKernel = kernelSize[0] * kernelSize[1] * kernelSize[2];

  m_Kernel.resize(totalKernel);
  std::fill(m_Kernel.begin(), m_Kernel.end(), 0.0f);
  determineKernel(kernelNumVoxels);

  std::vector<float> uniformKernel(totalKernel, 1.0f);

  if(m_StoreKernelDistances)
  {
    m_KernelValDistances.resize(totalKernel);
    std::fill(m_KernelValDistances.begin(), m_KernelValDistances.end(), 0.0f);
    determineKernelDistances(kernelNumVoxels, res);
  }

  size_t progIncrement = numVerts / 100;
  size_t prog = 1;
  size_t progressInt = 0;

  for(size_t i = 0; i < numVerts; i++)
  {
    if(getCancel())
    {
      break;
    }
    if(m_UseMask)
    {
      if(!m_Mask[i])
      {
        continue;
      }
    }
    index = m_VoxelIndices[i];
    if(index > maxImageIndex)
    {
      QString ss = QObject::tr("Index present in the selected Voxel Indices array that falls outside the selected Image Geometry for interpolation.\n Index = %1\n Max Image Index = %2\n")
                       .arg(index)
                       .arg(maxImageIndex);
      setErrorCondition(-1, ss);
    }
    x = index % dims[0];
    y = (index / dims[0]) % dims[1];
    z = index / (dims[0] * dims[1]);

    if(!m_SourceArraysToInterpolate.empty())
    {
      for(std::vector<IDataArray::WeakPointer>::size_type j = 0; j < m_SourceArraysToInterpolate.size(); j++)
      {
        EXECUTE_FUNCTION_TEMPLATE_NO_BOOL(this, mapPointCloudDataByKernel, m_SourceArraysToInterpolate[j].lock(), m_SourceArraysToInterpolate[j].lock(), m_DynamicArraysToInterpolate[j].lock(),
                                          m_Kernel, kernelNumVoxels, dims.data(), x, y, z, i)
      }
    }

    if(m_SourceArraysToCopy.size() > 0)
    {
      for(std::vector<IDataArray::WeakPointer>::size_type j = 0; j < m_SourceArraysToCopy.size(); j++)
      {
        EXECUTE_FUNCTION_TEMPLATE_NO_BOOL(this, mapPointCloudDataByKernel, m_SourceArraysToCopy[j].lock(), m_SourceArraysToCopy[j].lock(), m_DynamicArraysToCopy[j].lock(), uniformKernel,
                                          kernelNumVoxels, dims.data(), x, y, z, i)
      }
    }

    if(m_StoreKernelDistances)
    {
      mapKernelDistances(kernelNumVoxels, dims.data(), x, y, z);
    }

    if(i > prog)
    {
      progressInt = static_cast<int64_t>((static_cast<float>(i) / numVerts) * 100.0f);
      QString ss = QObject::tr("Interpolating Point Cloud || %1% Completed").arg(progressInt);
      notifyStatusMessage(ss);
      prog = prog + progIncrement;
    }
  }

  notifyStatusMessage("Complete");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AbstractFilter::Pointer InterpolatePointCloudToRegularGrid::newFilterInstance(bool copyFilterParameters) const
{
  InterpolatePointCloudToRegularGrid::Pointer filter = InterpolatePointCloudToRegularGrid::New();
  if(copyFilterParameters)
  {
    copyFilterParameterInstanceVariables(filter.get());
  }
  return filter;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString InterpolatePointCloudToRegularGrid::getCompiledLibraryName() const
{
  return DREAM3DReviewConstants::DREAM3DReviewBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString InterpolatePointCloudToRegularGrid::getBrandingString() const
{
  return "DREAM3DReview";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString InterpolatePointCloudToRegularGrid::getFilterVersion() const
{
  QString version;
  QTextStream vStream(&version);
  vStream << DREAM3DReview::Version::Major() << "." << DREAM3DReview::Version::Minor() << "." << DREAM3DReview::Version::Patch();
  return version;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString InterpolatePointCloudToRegularGrid::getGroupName() const
{
  return SIMPL::FilterGroups::SamplingFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString InterpolatePointCloudToRegularGrid::getSubGroupName() const
{
  return DREAM3DReviewConstants::FilterSubGroups::InterpolationFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString InterpolatePointCloudToRegularGrid::getHumanLabel() const
{
  return "Interpolate Point Cloud to Regular Grid";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QUuid InterpolatePointCloudToRegularGrid::getUuid() const
{
  return QUuid("{4b551c15-418d-5081-be3f-d3aeb62408e5}");
}

// -----------------------------------------------------------------------------
InterpolatePointCloudToRegularGrid::Pointer InterpolatePointCloudToRegularGrid::NullPointer()
{
  return Pointer(static_cast<Self*>(nullptr));
}

// -----------------------------------------------------------------------------
std::shared_ptr<InterpolatePointCloudToRegularGrid> InterpolatePointCloudToRegularGrid::New()
{
  struct make_shared_enabler : public InterpolatePointCloudToRegularGrid
  {
  };
  std::shared_ptr<make_shared_enabler> val = std::make_shared<make_shared_enabler>();
  val->setupFilterParameters();
  return val;
}

// -----------------------------------------------------------------------------
QString InterpolatePointCloudToRegularGrid::getNameOfClass() const
{
  return QString("InterpolatePointCloudToRegularGrid");
}

// -----------------------------------------------------------------------------
QString InterpolatePointCloudToRegularGrid::ClassName()
{
  return QString("InterpolatePointCloudToRegularGrid");
}

// -----------------------------------------------------------------------------
void InterpolatePointCloudToRegularGrid::setDataContainerName(const DataArrayPath& value)
{
  m_DataContainerName = value;
}

// -----------------------------------------------------------------------------
DataArrayPath InterpolatePointCloudToRegularGrid::getDataContainerName() const
{
  return m_DataContainerName;
}

// -----------------------------------------------------------------------------
void InterpolatePointCloudToRegularGrid::setArraysToInterpolate(const QVector<DataArrayPath>& value)
{
  m_ArraysToInterpolate = value;
}

// -----------------------------------------------------------------------------
QVector<DataArrayPath> InterpolatePointCloudToRegularGrid::getArraysToInterpolate() const
{
  return m_ArraysToInterpolate;
}

// -----------------------------------------------------------------------------
void InterpolatePointCloudToRegularGrid::setArraysToCopy(const QVector<DataArrayPath>& value)
{
  m_ArraysToCopy = value;
}

// -----------------------------------------------------------------------------
QVector<DataArrayPath> InterpolatePointCloudToRegularGrid::getArraysToCopy() const
{
  return m_ArraysToCopy;
}

// -----------------------------------------------------------------------------
void InterpolatePointCloudToRegularGrid::setVoxelIndicesArrayPath(const DataArrayPath& value)
{
  m_VoxelIndicesArrayPath = value;
}

// -----------------------------------------------------------------------------
DataArrayPath InterpolatePointCloudToRegularGrid::getVoxelIndicesArrayPath() const
{
  return m_VoxelIndicesArrayPath;
}

// -----------------------------------------------------------------------------
void InterpolatePointCloudToRegularGrid::setInterpolatedDataContainerName(const DataArrayPath& value)
{
  m_InterpolatedDataContainerName = value;
}

// -----------------------------------------------------------------------------
DataArrayPath InterpolatePointCloudToRegularGrid::getInterpolatedDataContainerName() const
{
  return m_InterpolatedDataContainerName;
}

// -----------------------------------------------------------------------------
void InterpolatePointCloudToRegularGrid::setInterpolatedAttributeMatrixName(const QString& value)
{
  m_InterpolatedAttributeMatrixName = value;
}

// -----------------------------------------------------------------------------
QString InterpolatePointCloudToRegularGrid::getInterpolatedAttributeMatrixName() const
{
  return m_InterpolatedAttributeMatrixName;
}

// -----------------------------------------------------------------------------
void InterpolatePointCloudToRegularGrid::setInterpolationTechnique(int value)
{
  m_InterpolationTechnique = value;
}

// -----------------------------------------------------------------------------
int InterpolatePointCloudToRegularGrid::getInterpolationTechnique() const
{
  return m_InterpolationTechnique;
}

// -----------------------------------------------------------------------------
void InterpolatePointCloudToRegularGrid::setKernelSize(const FloatVec3Type& value)
{
  m_KernelSize = value;
}

// -----------------------------------------------------------------------------
FloatVec3Type InterpolatePointCloudToRegularGrid::getKernelSize() const
{
  return m_KernelSize;
}

// -----------------------------------------------------------------------------
void InterpolatePointCloudToRegularGrid::setSigmas(const FloatVec3Type& value)
{
  m_Sigmas = value;
}

// -----------------------------------------------------------------------------
FloatVec3Type InterpolatePointCloudToRegularGrid::getSigmas() const
{
  return m_Sigmas;
}

// -----------------------------------------------------------------------------
void InterpolatePointCloudToRegularGrid::setUseMask(bool value)
{
  m_UseMask = value;
}

// -----------------------------------------------------------------------------
bool InterpolatePointCloudToRegularGrid::getUseMask() const
{
  return m_UseMask;
}

// -----------------------------------------------------------------------------
void InterpolatePointCloudToRegularGrid::setMaskArrayPath(const DataArrayPath& value)
{
  m_MaskArrayPath = value;
}

// -----------------------------------------------------------------------------
DataArrayPath InterpolatePointCloudToRegularGrid::getMaskArrayPath() const
{
  return m_MaskArrayPath;
}

// -----------------------------------------------------------------------------
void InterpolatePointCloudToRegularGrid::setStoreKernelDistances(bool value)
{
  m_StoreKernelDistances = value;
}

// -----------------------------------------------------------------------------
bool InterpolatePointCloudToRegularGrid::getStoreKernelDistances() const
{
  return m_StoreKernelDistances;
}

// -----------------------------------------------------------------------------
void InterpolatePointCloudToRegularGrid::setKernelDistancesArrayName(const QString& value)
{
  m_KernelDistancesArrayName = value;
}

// -----------------------------------------------------------------------------
QString InterpolatePointCloudToRegularGrid::getKernelDistancesArrayName() const
{
  return m_KernelDistancesArrayName;
}
