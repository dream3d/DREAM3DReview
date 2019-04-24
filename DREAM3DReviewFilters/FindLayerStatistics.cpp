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

#include "FindLayerStatistics.h"

#ifdef SIMPLib_USE_PARALLEL_ALGORITHMS
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/partitioner.h>
#include <tbb/task_scheduler_init.h>
#endif

#include "SIMPLib/Common/Constants.h"
#include "SIMPLib/Common/TemplateHelpers.h"
#include "SIMPLib/FilterParameters/AbstractFilterParametersReader.h"
#include "SIMPLib/FilterParameters/ChoiceFilterParameter.h"
#include "SIMPLib/FilterParameters/DataArraySelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/SeparatorFilterParameter.h"
#include "SIMPLib/FilterParameters/StringFilterParameter.h"
#include "SIMPLib/Geometry/ImageGeom.h"
#include "SIMPLib/Math/SIMPLibMath.h"

#include "DREAM3DReview/DREAM3DReviewConstants.h"
#include "DREAM3DReview/DREAM3DReviewVersion.h"

/**
 * @brief The CalcProjectedStatsImpl class implements a templated threaded algorithm for
 * determining the projected image statistics of a given volume.
 */
template <typename T> class CalcLayerStatsImpl
{

public:
  CalcLayerStatsImpl(T* data, float* min, float* max, float* avg, float* std, float* var, int32_t* startPoints, size_t* stride, size_t* dims)
  : m_Data(data)
  , m_Min(min)
  , m_Max(max)
  , m_Avg(avg)
  , m_Std(std)
  , m_Var(var)
  , m_StartPoints(startPoints)
  , m_Stride(stride)
  , m_Dims(dims)
  {
  }
  virtual ~CalcLayerStatsImpl()
  {
  }

  void convert(size_t start, size_t end) const
  {
    size_t shift = 0;
    size_t count = 0;
    size_t point = 0, newPoint = 0;
    T val = static_cast<T>(0);
    for(size_t i = start; i < end; i++)
    {
      count = 0;
      point = m_StartPoints[i];
      for(size_t k = 0; k < m_Dims[1]; k++)
      {
        shift = k * m_Stride[1];
        for(size_t j = 0; j < m_Dims[0]; j++)
        {
          newPoint = static_cast<int32_t>(point + shift + (j * m_Stride[0]));
          val = m_Data[newPoint];
          if(val > 0.0)
          {
            if(val < m_Min[i])
            {
              m_Min[i] = val;
            }
            if(val > m_Max[i])
            {
              m_Max[i] = val;
            }
            m_Avg[i] += val;
            count++;
          }
        }
      }
      m_Avg[i] /= count;
      for(size_t k = 0; k < m_Dims[1]; k++)
      {
        shift = k * m_Stride[1];
        for(size_t j = 0; j < m_Dims[0]; j++)
        {
          newPoint = static_cast<int32_t>(point + shift + (j * m_Stride[0]));
          val = m_Data[newPoint];
          if(val > 0.0)
          {
            m_Std[i] += ((val - m_Avg[i]) * (val - m_Avg[i]));
          }
        }
      }
      m_Std[i] /= count;
      m_Var[i] = m_Std[i];
      m_Std[i] = sqrt(m_Std[i]);
    }
  }

#ifdef SIMPLib_USE_PARALLEL_ALGORITHMS
  void operator()(const tbb::blocked_range<size_t>& r) const
  {
    convert(r.begin(), r.end());
  }
#endif
private:
  T* m_Data;
  float* m_Min;
  float* m_Max;
  float* m_Avg;
  float* m_Std;
  float* m_Var;
  int32_t* m_StartPoints;
  size_t* m_Stride;
  size_t* m_Dims;
};

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
FindLayerStatistics::FindLayerStatistics()
: m_SelectedArrayPath("", "", "")
, m_Plane(0)
, m_LayerIDsArrayName("LayerIDs")
, m_LayerAttributeMatrixName("")
, m_LayerMinArrayName("LayerMin")
, m_LayerMaxArrayName("LayerMax")
, m_LayerAvgArrayName("LayerAvg")
, m_LayerStdArrayName("LayerStd")
, m_LayerVarArrayName("LayerVar")
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
FindLayerStatistics::~FindLayerStatistics() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FindLayerStatistics::setupFilterParameters()
{
  FilterParameterVectorType parameters;
  {
    ChoiceFilterParameter::Pointer parameter = ChoiceFilterParameter::New();
    parameter->setHumanLabel("Layer of Interest");
    parameter->setPropertyName("Plane");
    parameter->setSetterCallback(SIMPL_BIND_SETTER(FindLayerStatistics, this, Plane));
    parameter->setGetterCallback(SIMPL_BIND_GETTER(FindLayerStatistics, this, Plane));

    QVector<QString> choices;
    choices.push_back("XY");
    choices.push_back("XZ");
    choices.push_back("YZ");
    parameter->setChoices(choices);
    parameter->setCategory(FilterParameter::Parameter);
    parameters.push_back(parameter);
  }
  parameters.push_back(SeparatorFilterParameter::New("Cell Data", FilterParameter::RequiredArray));
  {
    DataArraySelectionFilterParameter::RequirementType req =
        DataArraySelectionFilterParameter::CreateRequirement(SIMPL::Defaults::AnyPrimitive, 1, AttributeMatrix::Type::Cell, IGeometry::Type::Image);
    QVector<QString> daTypes;
    daTypes.push_back(SIMPL::TypeNames::Int8);
    daTypes.push_back(SIMPL::TypeNames::Int16);
    daTypes.push_back(SIMPL::TypeNames::Int32);
    daTypes.push_back(SIMPL::TypeNames::Int64);
    daTypes.push_back(SIMPL::TypeNames::UInt8);
    daTypes.push_back(SIMPL::TypeNames::UInt16);
    daTypes.push_back(SIMPL::TypeNames::UInt32);
    daTypes.push_back(SIMPL::TypeNames::UInt64);
    daTypes.push_back(SIMPL::TypeNames::Float);
    daTypes.push_back(SIMPL::TypeNames::Double);
    req.daTypes = daTypes;
    parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Attribute Array to Quantify", SelectedArrayPath, FilterParameter::RequiredArray, FindLayerStatistics, req));
  }
  parameters.push_back(SIMPL_NEW_STRING_FP("Layer Attribute Matrix", LayerAttributeMatrixName, FilterParameter::CreatedArray, FindLayerStatistics));
  parameters.push_back(SeparatorFilterParameter::New("Cell Data", FilterParameter::CreatedArray));
  parameters.push_back(SIMPL_NEW_STRING_FP("Layer IDs", LayerIDsArrayName, FilterParameter::CreatedArray, FindLayerStatistics));
  parameters.push_back(SeparatorFilterParameter::New("Feature Data", FilterParameter::CreatedArray));
  parameters.push_back(SIMPL_NEW_STRING_FP("Layer Min", LayerMinArrayName, FilterParameter::CreatedArray, FindLayerStatistics));
  parameters.push_back(SIMPL_NEW_STRING_FP("Layer Max", LayerMaxArrayName, FilterParameter::CreatedArray, FindLayerStatistics));
  parameters.push_back(SIMPL_NEW_STRING_FP("Layer Avg", LayerAvgArrayName, FilterParameter::CreatedArray, FindLayerStatistics));
  parameters.push_back(SIMPL_NEW_STRING_FP("Layer Std", LayerStdArrayName, FilterParameter::CreatedArray, FindLayerStatistics));
  parameters.push_back(SIMPL_NEW_STRING_FP("Layer Var", LayerVarArrayName, FilterParameter::CreatedArray, FindLayerStatistics));
  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
void FindLayerStatistics::readFilterParameters(AbstractFilterParametersReader* reader, int index)
{
  reader->openFilterGroup(this, index);
  setLayerAttributeMatrixName(reader->readString("LayerAttributeMatrixName", getLayerAttributeMatrixName()));
  setLayerIDsArrayName(reader->readString("LayerIDsArrayName", getLayerIDsArrayName()));
  setLayerVarArrayName(reader->readString("LayerVarArrayName", getLayerVarArrayName()));
  setLayerStdArrayName(reader->readString("LayerStdArrayName", getLayerStdArrayName()));
  setLayerAvgArrayName(reader->readString("LayerAvgArrayName", getLayerAvgArrayName()));
  setLayerMaxArrayName(reader->readString("LayerMaxArrayName", getLayerMaxArrayName()));
  setLayerMinArrayName(reader->readString("LayerMinArrayName", getLayerMinArrayName()));
  setSelectedArrayPath(reader->readDataArrayPath("SelectedArrayPath", getSelectedArrayPath()));
  setPlane(reader->readValue("Plane", getPlane()));
  reader->closeFilterGroup();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FindLayerStatistics::initialize()
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FindLayerStatistics::dataCheck()
{
  clearErrorCode();
  clearWarningCode();

  DataArrayPath tempPath;

  m_InDataPtr = getDataContainerArray()->getPrereqIDataArrayFromPath<IDataArray, AbstractFilter>(this, getSelectedArrayPath());
  if(nullptr != m_InDataPtr.lock())
  {
    if(TemplateHelpers::CanDynamicCast<BoolArrayType>()(m_InDataPtr.lock()))
    {
      QString ss = QObject::tr("Selected array cannot be of type bool.  The path is %1").arg(getSelectedArrayPath().serialize());
      setErrorCondition(-11001, ss);
    }
  }

  ImageGeom::Pointer image = getDataContainerArray()->getPrereqGeometryFromDataContainer<ImageGeom, AbstractFilter>(this, getSelectedArrayPath().getDataContainerName());
  if(getErrorCode() < 0)
  {
    return;
  }

  if(image->getXPoints() <= 1 || image->getYPoints() <= 1 || image->getZPoints() <= 1)
  {
    setErrorCondition(-999, "The Image Geometry is not 3D and cannot be run through this Filter");
  }

  size_t numLayers = 0;
  if(m_Plane == 0)
  {
    numLayers = image->getZPoints();
  }
  else if(m_Plane == 1)
  {
    numLayers = image->getYPoints();
  }
  else if(m_Plane == 2)
  {
    numLayers = image->getXPoints();
  }

  QVector<size_t> tDims(1, numLayers);
  DataContainer::Pointer m = getDataContainerArray()->getDataContainer(getSelectedArrayPath().getDataContainerName());
  m->createNonPrereqAttributeMatrix(this, getLayerAttributeMatrixName(), tDims, AttributeMatrix::Type::CellFeature);

  QVector<size_t> cDims(1, 1);
  tempPath.update(getSelectedArrayPath().getDataContainerName(), getSelectedArrayPath().getAttributeMatrixName(), getLayerIDsArrayName());
  m_LayerIDsPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<int32_t>, AbstractFilter, int32_t>(this, tempPath, 0,
                                                                                                                     cDims); /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
  if(nullptr != m_LayerIDsPtr.lock()) /* Validate the Weak Pointer wraps a non-nullptr pointer to a DataArray<T> object */
  {
    m_LayerIDs = m_LayerIDsPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */

  tempPath.update(getSelectedArrayPath().getDataContainerName(), getLayerAttributeMatrixName(), getLayerMinArrayName());
  m_LayerMinPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<float>, AbstractFilter, float>(this, tempPath, 0,
                                                                                                                 cDims); /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
  if(nullptr != m_LayerMinPtr.lock()) /* Validate the Weak Pointer wraps a non-nullptr pointer to a DataArray<T> object */
  {
    m_LayerMin = m_LayerMinPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */

  tempPath.update(getSelectedArrayPath().getDataContainerName(), getLayerAttributeMatrixName(), getLayerMaxArrayName());
  m_LayerMaxPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<float>, AbstractFilter, float>(this, tempPath, 0,
                                                                                                                 cDims); /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
  if(nullptr != m_LayerMaxPtr.lock()) /* Validate the Weak Pointer wraps a non-nullptr pointer to a DataArray<T> object */
  {
    m_LayerMax = m_LayerMaxPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */

  tempPath.update(getSelectedArrayPath().getDataContainerName(), getLayerAttributeMatrixName(), getLayerAvgArrayName());
  m_LayerAvgPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<float>, AbstractFilter, float>(this, tempPath, 0,
                                                                                                                 cDims); /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
  if(nullptr != m_LayerAvgPtr.lock()) /* Validate the Weak Pointer wraps a non-nullptr pointer to a DataArray<T> object */
  {
    m_LayerAvg = m_LayerAvgPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */

  tempPath.update(getSelectedArrayPath().getDataContainerName(), getLayerAttributeMatrixName(), getLayerStdArrayName());
  m_LayerStdPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<float>, AbstractFilter, float>(this, tempPath, 0,
                                                                                                                 cDims); /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
  if(nullptr != m_LayerStdPtr.lock()) /* Validate the Weak Pointer wraps a non-nullptr pointer to a DataArray<T> object */
  {
    m_LayerStd = m_LayerStdPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */

  tempPath.update(getSelectedArrayPath().getDataContainerName(), getLayerAttributeMatrixName(), getLayerVarArrayName());
  m_LayerVarPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<float>, AbstractFilter, float>(this, tempPath, 0,
                                                                                                                 cDims); /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
  if(nullptr != m_LayerVarPtr.lock()) /* Validate the Weak Pointer wraps a non-nullptr pointer to a DataArray<T> object */
  {
    m_LayerVar = m_LayerVarPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FindLayerStatistics::preflight()
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
void FindLayerStatistics::execute()
{
  clearErrorCode();
  clearWarningCode();
  dataCheck();
  if(getErrorCode() < 0)
  {
    return;
  }

  DataContainer::Pointer m = getDataContainerArray()->getDataContainer(getSelectedArrayPath().getDataContainerName());

#ifdef SIMPLib_USE_PARALLEL_ALGORITHMS
  tbb::task_scheduler_init init;
  bool doParallel = true;
#endif

  SizeVec3Type dimsP = m->getGeometryAs<ImageGeom>()->getDimensions();

  Int32ArrayType::Pointer startingPoints = Int32ArrayType::CreateArray(0, "_INTERNAL_USE_ONLY_startingPoints");
  int32_t* startPoints = nullptr;
  size_t stride[2];
  size_t dims[2];
  // size_t yStride = 0;
  size_t depth = 0;
  if(m_Plane == 0)
  {
    startingPoints->resizeTuples(dimsP[2]);
    startPoints = startingPoints->getPointer(0);
    stride[0] = 1;
    stride[1] = dimsP[0];
    dims[0] = dimsP[0];
    dims[1] = dimsP[1];
    depth = dimsP[2];
    for(size_t i = 0; i < dimsP[2]; i++)
    {
      startPoints[i] = i * (dimsP[0] * dimsP[1]);
    }
  }
  if(m_Plane == 1)
  {
    startingPoints->resizeTuples(dimsP[1]);
    startPoints = startingPoints->getPointer(0);
    stride[0] = 1;
    stride[1] = (dimsP[0] * dimsP[1]);
    dims[0] = dimsP[0];
    dims[1] = dimsP[2];
    depth = dimsP[1];
    for(size_t i = 0; i < dimsP[1]; i++)
    {
      startPoints[i] = i * dimsP[0];
    }
  }
  if(m_Plane == 2)
  {
    startingPoints->resizeTuples(dimsP[0]);
    startPoints = startingPoints->getPointer(0);
    stride[0] = dimsP[0];
    stride[1] = (dimsP[0] * dimsP[1]);
    dims[0] = dimsP[1];
    dims[1] = dimsP[2];
    depth = dimsP[0];
    for(size_t i = 0; i < dimsP[0]; i++)
    {
      startPoints[i] = i;
    }
  }

  if(nullptr == startPoints)
  {
    QString ss = QObject::tr("Unable to establish starting location for supplied plane. The plane is %1").arg(m_Plane);
    setErrorCondition(-11001, ss);
    return;
  }

  if(TemplateHelpers::CanDynamicCast<Int8ArrayType>()(m_InDataPtr.lock()))
  {
    Int8ArrayType::Pointer cellArray = std::dynamic_pointer_cast<Int8ArrayType>(m_InDataPtr.lock());
    int8_t* cPtr = cellArray->getPointer(0);
#ifdef SIMPLib_USE_PARALLEL_ALGORITHMS
    if(doParallel == true)
    {
      tbb::parallel_for(tbb::blocked_range<size_t>(0, depth), CalcLayerStatsImpl<int8_t>(cPtr, m_LayerMin, m_LayerMax, m_LayerAvg, m_LayerStd, m_LayerVar, startPoints, stride, dims),
                        tbb::auto_partitioner());
    }
    else
#endif
    {
      CalcLayerStatsImpl<int8_t> serial(cPtr, m_LayerMin, m_LayerMax, m_LayerAvg, m_LayerStd, m_LayerVar, startPoints, stride, dims);
      serial.convert(0, depth);
    }
  }
  else if(TemplateHelpers::CanDynamicCast<UInt8ArrayType>()(m_InDataPtr.lock()))
  {
    UInt8ArrayType::Pointer cellArray = std::dynamic_pointer_cast<UInt8ArrayType>(m_InDataPtr.lock());
    uint8_t* cPtr = cellArray->getPointer(0);
#ifdef SIMPLib_USE_PARALLEL_ALGORITHMS
    if(doParallel == true)
    {
      tbb::parallel_for(tbb::blocked_range<size_t>(0, depth), CalcLayerStatsImpl<uint8_t>(cPtr, m_LayerMin, m_LayerMax, m_LayerAvg, m_LayerStd, m_LayerVar, startPoints, stride, dims),
                        tbb::auto_partitioner());
    }
    else
#endif
    {
      CalcLayerStatsImpl<uint8_t> serial(cPtr, m_LayerMin, m_LayerMax, m_LayerAvg, m_LayerStd, m_LayerVar, startPoints, stride, dims);
      serial.convert(0, depth);
    }
  }
  else if(TemplateHelpers::CanDynamicCast<Int16ArrayType>()(m_InDataPtr.lock()))
  {
    Int16ArrayType::Pointer cellArray = std::dynamic_pointer_cast<Int16ArrayType>(m_InDataPtr.lock());
    int16_t* cPtr = cellArray->getPointer(0);
#ifdef SIMPLib_USE_PARALLEL_ALGORITHMS
    if(doParallel == true)
    {
      tbb::parallel_for(tbb::blocked_range<size_t>(0, depth), CalcLayerStatsImpl<int16_t>(cPtr, m_LayerMin, m_LayerMax, m_LayerAvg, m_LayerStd, m_LayerVar, startPoints, stride, dims),
                        tbb::auto_partitioner());
    }
    else
#endif
    {
      CalcLayerStatsImpl<int16_t> serial(cPtr, m_LayerMin, m_LayerMax, m_LayerAvg, m_LayerStd, m_LayerVar, startPoints, stride, dims);
      serial.convert(0, depth);
    }
  }
  else if(TemplateHelpers::CanDynamicCast<UInt16ArrayType>()(m_InDataPtr.lock()))
  {
    UInt16ArrayType::Pointer cellArray = std::dynamic_pointer_cast<UInt16ArrayType>(m_InDataPtr.lock());
    uint16_t* cPtr = cellArray->getPointer(0);
#ifdef SIMPLib_USE_PARALLEL_ALGORITHMS
    if(doParallel == true)
    {
      tbb::parallel_for(tbb::blocked_range<size_t>(0, depth), CalcLayerStatsImpl<uint16_t>(cPtr, m_LayerMin, m_LayerMax, m_LayerAvg, m_LayerStd, m_LayerVar, startPoints, stride, dims),
                        tbb::auto_partitioner());
    }
    else
#endif
    {
      CalcLayerStatsImpl<uint16_t> serial(cPtr, m_LayerMin, m_LayerMax, m_LayerAvg, m_LayerStd, m_LayerVar, startPoints, stride, dims);
      serial.convert(0, depth);
    }
  }
  else if(TemplateHelpers::CanDynamicCast<Int32ArrayType>()(m_InDataPtr.lock()))
  {
    Int32ArrayType::Pointer cellArray = std::dynamic_pointer_cast<Int32ArrayType>(m_InDataPtr.lock());
    int32_t* cPtr = cellArray->getPointer(0);
#ifdef SIMPLib_USE_PARALLEL_ALGORITHMS
    if(doParallel == true)
    {
      tbb::parallel_for(tbb::blocked_range<size_t>(0, depth), CalcLayerStatsImpl<int32_t>(cPtr, m_LayerMin, m_LayerMax, m_LayerAvg, m_LayerStd, m_LayerVar, startPoints, stride, dims),
                        tbb::auto_partitioner());
    }
    else
#endif
    {
      CalcLayerStatsImpl<int32_t> serial(cPtr, m_LayerMin, m_LayerMax, m_LayerAvg, m_LayerStd, m_LayerVar, startPoints, stride, dims);
      serial.convert(0, depth);
    }
  }
  else if(TemplateHelpers::CanDynamicCast<UInt32ArrayType>()(m_InDataPtr.lock()))
  {
    UInt32ArrayType::Pointer cellArray = std::dynamic_pointer_cast<UInt32ArrayType>(m_InDataPtr.lock());
    uint32_t* cPtr = cellArray->getPointer(0);
#ifdef SIMPLib_USE_PARALLEL_ALGORITHMS
    if(doParallel == true)
    {
      tbb::parallel_for(tbb::blocked_range<size_t>(0, depth), CalcLayerStatsImpl<uint32_t>(cPtr, m_LayerMin, m_LayerMax, m_LayerAvg, m_LayerStd, m_LayerVar, startPoints, stride, dims),
                        tbb::auto_partitioner());
    }
    else
#endif
    {
      CalcLayerStatsImpl<uint32_t> serial(cPtr, m_LayerMin, m_LayerMax, m_LayerAvg, m_LayerStd, m_LayerVar, startPoints, stride, dims);
      serial.convert(0, depth);
    }
  }
  else if(TemplateHelpers::CanDynamicCast<Int64ArrayType>()(m_InDataPtr.lock()))
  {
    Int64ArrayType::Pointer cellArray = std::dynamic_pointer_cast<Int64ArrayType>(m_InDataPtr.lock());
    int64_t* cPtr = cellArray->getPointer(0);
#ifdef SIMPLib_USE_PARALLEL_ALGORITHMS
    if(doParallel == true)
    {
      tbb::parallel_for(tbb::blocked_range<size_t>(0, depth), CalcLayerStatsImpl<int64_t>(cPtr, m_LayerMin, m_LayerMax, m_LayerAvg, m_LayerStd, m_LayerVar, startPoints, stride, dims),
                        tbb::auto_partitioner());
    }
    else
#endif
    {
      CalcLayerStatsImpl<int64_t> serial(cPtr, m_LayerMin, m_LayerMax, m_LayerAvg, m_LayerStd, m_LayerVar, startPoints, stride, dims);
      serial.convert(0, depth);
    }
  }
  else if(TemplateHelpers::CanDynamicCast<UInt64ArrayType>()(m_InDataPtr.lock()))
  {
    UInt64ArrayType::Pointer cellArray = std::dynamic_pointer_cast<UInt64ArrayType>(m_InDataPtr.lock());
    uint64_t* cPtr = cellArray->getPointer(0);
#ifdef SIMPLib_USE_PARALLEL_ALGORITHMS
    if(doParallel == true)
    {
      tbb::parallel_for(tbb::blocked_range<size_t>(0, depth), CalcLayerStatsImpl<uint64_t>(cPtr, m_LayerMin, m_LayerMax, m_LayerAvg, m_LayerStd, m_LayerVar, startPoints, stride, dims),
                        tbb::auto_partitioner());
    }
    else
#endif
    {
      CalcLayerStatsImpl<uint64_t> serial(cPtr, m_LayerMin, m_LayerMax, m_LayerAvg, m_LayerStd, m_LayerVar, startPoints, stride, dims);
      serial.convert(0, depth);
    }
  }
  else if(TemplateHelpers::CanDynamicCast<FloatArrayType>()(m_InDataPtr.lock()))
  {
    FloatArrayType::Pointer cellArray = std::dynamic_pointer_cast<FloatArrayType>(m_InDataPtr.lock());
    float* cPtr = cellArray->getPointer(0);
#ifdef SIMPLib_USE_PARALLEL_ALGORITHMS
    if(doParallel == true)
    {
      tbb::parallel_for(tbb::blocked_range<size_t>(0, depth), CalcLayerStatsImpl<float>(cPtr, m_LayerMin, m_LayerMax, m_LayerAvg, m_LayerStd, m_LayerVar, startPoints, stride, dims),
                        tbb::auto_partitioner());
    }
    else
#endif
    {
      CalcLayerStatsImpl<float> serial(cPtr, m_LayerMin, m_LayerMax, m_LayerAvg, m_LayerStd, m_LayerVar, startPoints, stride, dims);
      serial.convert(0, depth);
    }
  }
  else if(TemplateHelpers::CanDynamicCast<DoubleArrayType>()(m_InDataPtr.lock()))
  {
    DoubleArrayType::Pointer cellArray = std::dynamic_pointer_cast<DoubleArrayType>(m_InDataPtr.lock());
    double* cPtr = cellArray->getPointer(0);
#ifdef SIMPLib_USE_PARALLEL_ALGORITHMS
    if(doParallel == true)
    {
      tbb::parallel_for(tbb::blocked_range<size_t>(0, depth), CalcLayerStatsImpl<double>(cPtr, m_LayerMin, m_LayerMax, m_LayerAvg, m_LayerStd, m_LayerVar, startPoints, stride, dims),
                        tbb::auto_partitioner());
    }
    else
#endif
    {
      CalcLayerStatsImpl<double> serial(cPtr, m_LayerMin, m_LayerMax, m_LayerAvg, m_LayerStd, m_LayerVar, startPoints, stride, dims);
      serial.convert(0, depth);
    }
  }
  else
  {
    QString ss = QObject::tr("Selected array is of unsupported type. The type is %1").arg(m_InDataPtr.lock()->getTypeAsString());
    setErrorCondition(-11001, ss);
    return;
  }

  notifyStatusMessage("Complete");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AbstractFilter::Pointer FindLayerStatistics::newFilterInstance(bool copyFilterParameters) const
{
  FindLayerStatistics::Pointer filter = FindLayerStatistics::New();
  if(copyFilterParameters)
  {
    copyFilterParameterInstanceVariables(filter.get());
  }
  return filter;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString FindLayerStatistics::getCompiledLibraryName() const
{
  return DREAM3DReviewConstants::DREAM3DReviewBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString FindLayerStatistics::getBrandingString() const
{
  return "DREAM3DReview";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString FindLayerStatistics::getFilterVersion() const
{
  QString version;
  QTextStream vStream(&version);
  vStream << DREAM3DReview::Version::Major() << "." << DREAM3DReview::Version::Minor() << "." << DREAM3DReview::Version::Patch();
  return version;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString FindLayerStatistics::getGroupName() const
{
  return SIMPL::FilterGroups::StatisticsFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString FindLayerStatistics::getSubGroupName() const
{
  return SIMPL::FilterSubGroups::ImageFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString FindLayerStatistics::getHumanLabel() const
{
  return "Find Layer Statistics";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QUuid FindLayerStatistics::getUuid()
{
  return QUuid("{c85277c7-9f02-5891-8141-04523658737d}");
}
