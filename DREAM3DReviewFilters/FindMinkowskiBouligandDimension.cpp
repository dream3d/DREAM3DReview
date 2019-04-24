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
#include "FindMinkowskiBouligandDimension.h"

#include <cmath>
#include <cstring>

#include "SIMPLib/Common/Constants.h"
#include "SIMPLib/FilterParameters/DataArraySelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/StringFilterParameter.h"
#include "SIMPLib/Geometry/ImageGeom.h"

#include "DREAM3DReview/DREAM3DReviewConstants.h"
#include "DREAM3DReview/DREAM3DReviewVersion.h"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
FindMinkowskiBouligandDimension::FindMinkowskiBouligandDimension()
: m_MaskArrayPath(SIMPL::Defaults::ImageDataContainerName, SIMPL::Defaults::CellAttributeMatrixName, SIMPL::CellData::Mask)
, m_AttributeMatrixName("FractalData")
, m_MinkowskiBouligandDimensionArrayName("MinkowskiBouligandDimension")
, m_Mask(nullptr)
, m_MinkowskiBouligandDimension(nullptr)
{
  initialize();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
FindMinkowskiBouligandDimension::~FindMinkowskiBouligandDimension() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FindMinkowskiBouligandDimension::initialize()
{
  clearErrorCode();
  clearWarningCode();
  setCancel(false);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FindMinkowskiBouligandDimension::setupFilterParameters()
{
  FilterParameterVectorType parameters;
  DataArraySelectionFilterParameter::RequirementType dasReq = DataArraySelectionFilterParameter::CreateRequirement(SIMPL::TypeNames::Bool, 1, AttributeMatrix::Type::Cell, IGeometry::Type::Image);
  parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Mask", MaskArrayPath, FilterParameter::RequiredArray, FindMinkowskiBouligandDimension, dasReq));
  parameters.push_back(SIMPL_NEW_STRING_FP("Fractal Attribute Matrix", AttributeMatrixName, FilterParameter::CreatedArray, FindMinkowskiBouligandDimension));
  parameters.push_back(SIMPL_NEW_STRING_FP("Minkowski-Bouligand Dimension", MinkowskiBouligandDimensionArrayName, FilterParameter::CreatedArray, FindMinkowskiBouligandDimension));
  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FindMinkowskiBouligandDimension::dataCheck()
{
  clearErrorCode();
  clearWarningCode();

  ImageGeom::Pointer image = getDataContainerArray()->getPrereqGeometryFromDataContainer<ImageGeom, AbstractFilter>(this, getMaskArrayPath().getDataContainerName());

  if(getErrorCode() < 0)
  {
    return;
  }

  FloatVec3Type res = image->getSpacing();
  if(std::adjacent_find(std::begin(res), std::end(res), std::not_equal_to<float>()) != std::end(res))
  {
    QString ss = QObject::tr("The input Image Geometry must have isotropic resolution");
    setErrorCondition(-1, ss);
  }

  QVector<size_t> cDims(1, 1);

  m_MaskPtr = getDataContainerArray()->getPrereqArrayFromPath<DataArray<bool>, AbstractFilter>(this, getMaskArrayPath(), cDims);
  if(nullptr != m_MaskPtr.lock().get())
  {
    m_Mask = m_MaskPtr.lock()->getPointer(0);
  }

  DataContainer::Pointer dc = getDataContainerArray()->getPrereqDataContainer<AbstractFilter>(this, getMaskArrayPath().getDataContainerName());

  if(getErrorCode() < 0)
  {
    return;
  }

  dc->createNonPrereqAttributeMatrix(this, getAttributeMatrixName(), QVector<size_t>(1, 1), AttributeMatrix::Type::CellFeature);

  DataArrayPath path(getMaskArrayPath().getDataContainerName(), getAttributeMatrixName(), getMinkowskiBouligandDimensionArrayName());

  m_MinkowskiBouligandDimensionPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<double>, AbstractFilter, double>(this, path, 0, cDims);
  if(nullptr != m_MinkowskiBouligandDimensionPtr.lock().get())
  {
    m_MinkowskiBouligandDimension = m_MinkowskiBouligandDimensionPtr.lock()->getPointer(0);
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FindMinkowskiBouligandDimension::preflight()
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
void FindMinkowskiBouligandDimension::findSourceIndices(std::vector<size_t>& indices, size_t sourceDim, size_t index, size_t box[3], size_t imageDim)
{
  size_t zloop = imageDim == 3 ? 2 : 1;
  for(size_t z = 0; z < zloop; z++)
  {
    size_t zStride = (z * sourceDim * sourceDim) + (box[2] * sourceDim * sourceDim * 2);
    for(size_t y = 0; y < 2; y++)
    {
      size_t yStride = (y * sourceDim) + (box[1] * sourceDim * 2);
      for(size_t x = 0; x < 2; x++)
      {
        size_t idx = (z * 2 * 2) + (y * 2) + x;
        size_t mov_idx = zStride + yStride + x + (box[0] * 2);
        indices[idx] = mov_idx;
      }
    }
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FindMinkowskiBouligandDimension::execute()
{
  initialize();
  dataCheck();
  if(getErrorCode() < 0)
  {
    return;
  }

  ImageGeom::Pointer image = getDataContainerArray()->getDataContainer(m_MaskArrayPath.getDataContainerName())->getGeometryAs<ImageGeom>();
  SizeVec3Type dims = image->getDimensions();
  // FloatVec3Type res = image->getSpacing();

  size_t maxDim = std::max({dims[0], dims[1], dims[2]});
  size_t imageDim = 3;

  if(std::find(std::begin(dims), std::end(dims), 1) != std::end(dims))
  {
    imageDim = 2;
    if(dims[0] == 1)
    {
      dims[0] = dims[1];
      dims[1] = dims[2];
    }
    else if(dims[1] == 1)
    {
      dims[1] = dims[2];
    }
    dims[2] = 1;
  }

  auto is_pow_2 = [](size_t x) -> bool { return !(x & (x - 1)); };

  auto next_pow_2 = [&](size_t x) -> size_t {
    if(is_pow_2(x))
    {
      return x;
    }
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    x |= x >> 32;
    return x + 1;
  };

  auto floor_log2 = [](size_t x) -> size_t {
    size_t exp = 0;
    while(x >>= 1)
    {
      ++exp;
    }
    return exp;
  };

  maxDim = next_pow_2(maxDim);
  size_t exponent = floor_log2(maxDim);
  size_t totalSize = std::pow(maxDim, imageDim);

  std::vector<bool> source(totalSize, 0);
  size_t offsets[3] = {(maxDim - dims[0]) / 2, (maxDim - dims[1]) / 2, (maxDim - dims[2]) / 2};

  if(imageDim == 2)
  {
    offsets[2] = 0;
  }

  for(size_t z = 0; z < dims[2]; z++)
  {
    size_t zStride = (z * dims[1] * dims[0]);
    for(size_t y = 0; y < dims[1]; y++)
    {
      size_t yStride = (y * dims[0]);
      for(size_t x = 0; x < dims[0]; x++)
      {
        size_t index = zStride + yStride + x;
        size_t gridIndex = ((z + offsets[2]) * maxDim * maxDim) + ((y + offsets[1]) * maxDim) + (x + offsets[0]);
        source[gridIndex] = m_Mask[index];
      }
    }
  }

  size_t edgeLength = maxDim;
  std::vector<size_t> boxDims(exponent + 1, edgeLength);
  std::generate(std::next(std::begin(boxDims)), std::end(boxDims), [&edgeLength]() {
    edgeLength /= 2;
    return edgeLength;
  });

  std::vector<size_t> covering(exponent + 1, 0);
  size_t nextDim = *std::next(std::begin(boxDims));
  std::vector<bool> sink;
  sink.reserve(std::pow(nextDim, imageDim));
  std::vector<size_t> indices(std::pow(2, imageDim), 0);
  auto dimsIter = std::next(std::begin(boxDims));
  auto coverIter = std::next(std::begin(covering));
  size_t box[3] = {0, 0, 0};

  for(auto&& box : source)
  {
    if(box)
    {
      (*std::begin(covering))++;
    }
  }

  for(; dimsIter != std::prev(std::end(boxDims)); ++dimsIter, ++coverIter)
  {
    sink.resize(std::pow(*dimsIter, imageDim));

    size_t zloop = imageDim == 3 ? *dimsIter : 1;
    for(size_t z = 0; z < zloop; z++)
    {
      for(size_t y = 0; y < *dimsIter; y++)
      {
        for(size_t x = 0; x < *dimsIter; x++)
        {
          size_t index = (z * (*dimsIter) * (*dimsIter)) + (y * (*dimsIter)) + x;
          box[0] = x;
          box[1] = y;
          box[2] = z;
          findSourceIndices(indices, *std::prev(dimsIter), index, box, imageDim);
          for(auto&& ind : indices)
          {
            if(source[ind])
            {
              sink[index] = true;
              break;
            }
          }
        }
      }
    }
    for(auto&& box : sink)
    {
      if(box)
      {
        (*coverIter)++;
      }
    }

    source = sink;
  }

  covering.back() = 1;

  std::vector<double> LnNumBoxes(covering.size(), 0);
  std::vector<double> LnOneOverE(boxDims.size(), 0);
#if 0
  std::transform(std::rbegin(covering), std::rend(covering), std::begin(LnNumBoxes), [](size_t& x) -> double {
    return std::log(x);
  });
#else
  typedef std::vector<size_t>::iterator iter_type;
  // ? 0 1 2 3 4 5 6 7 8 9 ?
  iter_type from(covering.begin());                 //   ^
                                                    //         ------>
  iter_type until(covering.end());                  //                       ^
                                                    //
  std::reverse_iterator<iter_type> rev_until(from); // ^
                                                    //         <------
  std::reverse_iterator<iter_type> rev_from(until); //                     ^
  std::vector<double>::size_type lnIdx = 0;
  while(rev_from != rev_until)
  {
    LnNumBoxes[lnIdx] = static_cast<double>(std::log(*rev_from));
    lnIdx++;
    rev_from++;
  }
#endif

  std::transform(std::begin(boxDims), std::end(boxDims), std::begin(LnOneOverE), [](size_t& x) -> double { return std::log(1.0 / x); });

  auto linear_regression_slope = [&]() -> double {
    double xmean = std::accumulate(std::begin(LnOneOverE), std::end(LnOneOverE), 0.0);
    double ymean = std::accumulate(std::begin(LnNumBoxes), std::end(LnNumBoxes), 0.0);
    xmean /= LnOneOverE.size();
    ymean /= LnNumBoxes.size();

    std::vector<double> xx;
    std::transform(std::begin(LnOneOverE), std::end(LnOneOverE), std::begin(LnOneOverE), std::back_inserter(xx), std::multiplies<double>());
    double sumxx = std::accumulate(std::begin(xx), std::end(xx), 0.0);
    double ssxx = sumxx - (xx.size() * xmean * xmean);

    std::vector<double> xy;
    std::transform(std::begin(LnOneOverE), std::end(LnOneOverE), std::begin(LnNumBoxes), std::back_inserter(xy), std::multiplies<double>());
    double sumxy = std::accumulate(std::begin(xy), std::end(xy), 0.0);
    double ssxy = sumxy - (xy.size() * xmean * ymean);

    return ssxy / ssxx;
  };

  m_MinkowskiBouligandDimension[0] = linear_regression_slope();

  notifyStatusMessage("Complete");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AbstractFilter::Pointer FindMinkowskiBouligandDimension::newFilterInstance(bool copyFilterParameters) const
{
  FindMinkowskiBouligandDimension::Pointer filter = FindMinkowskiBouligandDimension::New();
  if(copyFilterParameters)
  {
    copyFilterParameterInstanceVariables(filter.get());
  }
  return filter;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString FindMinkowskiBouligandDimension::getCompiledLibraryName() const
{
  return DREAM3DReviewConstants::DREAM3DReviewBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString FindMinkowskiBouligandDimension::getBrandingString() const
{
  return "DREAM3DReview";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString FindMinkowskiBouligandDimension::getFilterVersion() const
{
  QString version;
  QTextStream vStream(&version);
  vStream << DREAM3DReview::Version::Major() << "." << DREAM3DReview::Version::Minor() << "." << DREAM3DReview::Version::Patch();
  return version;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString FindMinkowskiBouligandDimension::getGroupName() const
{
  return SIMPL::FilterGroups::StatisticsFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString FindMinkowskiBouligandDimension::getSubGroupName() const
{
  return SIMPL::FilterSubGroups::GeometryFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString FindMinkowskiBouligandDimension::getHumanLabel() const
{
  return "Find Minkowski-Bouligand Dimension";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QUuid FindMinkowskiBouligandDimension::getUuid()
{
  return QUuid("{6cc3148c-0bad-53b4-9568-ee1971cadb00}");
}
