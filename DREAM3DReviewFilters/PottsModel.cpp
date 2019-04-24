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

#include "PottsModel.h"

#include <cassert>
#include <chrono>
#include <deque>
#include <random>
#include <utility>
#include <utility>

#include "SIMPLib/Common/Constants.h"

#include "SIMPLib/FilterParameters/BooleanFilterParameter.h"
#include "SIMPLib/FilterParameters/DataArraySelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/DataArraySelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/DoubleFilterParameter.h"
#include "SIMPLib/FilterParameters/IntFilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedBooleanFilterParameter.h"
#include "SIMPLib/FilterParameters/SeparatorFilterParameter.h"

#include "SIMPLib/Geometry/ImageGeom.h"

#include "DREAM3DReview/DREAM3DReviewConstants.h"
#include "DREAM3DReview/DREAM3DReviewVersion.h"

namespace
{
using Neighborhood = std::vector<std::array<int8_t, 3>>;

const double BOLTZMANN = 1.38064852e-23;

enum class Dimension : uint8_t
{
  Two,
  Three
};

class SpinLattice
{
public:
  SpinLattice(ImageGeom::Pointer image, double temperature, bool periodic, int32_t* fIds, bool* mask)
  : image_(std::move(image))
  , temperature_(temperature)
  , periodic_(periodic)
  , fIds_(fIds)
  , mask_(mask)
  {
    initialize();
  }

  virtual ~SpinLattice()= default;

  void attempt_flip(size_t index)
  {
    int32_t spin = fIds_[index];
    std::deque<size_t> neighbors = valid_neighbor_query(index);

    auto same = std::count_if(std::begin(neighbors), std::end(neighbors), [&](size_t& idx) { return spin == fIds_[idx]; });

    if(same == neighbors.size())
    {
      return;
    }

    size_t rand_neighbor = random_neighbor(neighbors);
    int32_t candidate = fIds_[rand_neighbor];

    if(candidate == 0)
    {
      return;
    }

    if(spin != candidate)
    {
      double dE = 0.0;
      for(auto&& neighbor : neighbors)
      {
        int32_t n = fIds_[neighbor];
        double a = (n != candidate) ? 1.0 : 0.0;
        double b = (n != spin) ? 1.0 : 0.0;
        dE += (a - b);
      }
      dE *= 0.5;
      if(dE <= 0.0)
      {
        fIds_[index] = candidate;
        total_flips_++;
      }
      else
      {
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        double r = dist(gen_);
        if(r < std::exp(-dE / kT_))
        {
          fIds_[index] = candidate;
          total_flips_++;
        }
      }
    }
  }

  size_t total_flips()
  {
    return total_flips_;
  }

  Dimension dimension()
  {
    return dim_type_;
  }

private:
  size_t random_neighbor(std::deque<size_t>& neighbors)
  {
    std::uniform_int_distribution<size_t> dist(0, neighbors.size() - 1);
    return neighbors[dist(gen_)];
  }

  std::deque<size_t> valid_neighbor_query(size_t index)
  {
    std::deque<size_t> neighbors;

    size_t x = index % dims_[0];
    size_t y = (index / dims_[0]) % dims_[1];
    size_t z = index / (dims_[0] * dims_[1]);

    for(auto&& neighbor : neighborhood_)
    {
      if(periodic_)
      {
        size_t modx = apply_modular_operation(x, neighbor[0], dims_[0]);
        size_t mody = apply_modular_operation(y, neighbor[1], dims_[1]);
        size_t modz = apply_modular_operation(z, neighbor[2], dims_[2]);
        size_t neigh = neighbor_index(modx, mody, modz);
        if(mask_ != nullptr)
        {
          if(mask_[neigh])
          {
            neighbors.push_back(neigh);
          }
        }
        else
        {
          neighbors.push_back(neigh);
        }
      }
      else
      {
        int64_t modx = static_cast<int64_t>(x) + neighbor[0];
        int64_t mody = static_cast<int64_t>(y) + neighbor[1];
        int64_t modz = static_cast<int64_t>(z) + neighbor[2];
        if((modx >= 0 && modx < static_cast<int64_t>(dims_[0])) && (mody >= 0 && mody < static_cast<int64_t>(dims_[1])) && (modz >= 0 && modz < static_cast<int64_t>(dims_[2])))
        {
          size_t neigh = neighbor_index(modx, mody, modz);
          if(mask_ != nullptr)
          {
            if(mask_[neigh])
            {
              neighbors.push_back(neigh);
            }
          }
          else
          {
            neighbors.push_back(neigh);
          }
        }
      }
    }
    return neighbors;
  }

  size_t modular_subtraction(size_t a, size_t b, size_t m)
  {
    if(a >= b)
    {
      return a - b;
    }
    
    
      return m - b + a;
    
  }

  size_t modular_addition(size_t a, size_t b, size_t m)
  {
    if(b == 0)
    {
      return a;
    }
    b = m - b;
    if(a >= b)
    {
      return a - b;
    }
    
    
      return m - b + a;
    
  }

  size_t apply_modular_operation(size_t a, int8_t b, size_t m)
  {
    if(b < 0)
    {
      b *= -1;
      //size_t c = static_cast<size_t>(b);
      return modular_subtraction(a, b, m);
    }
    
    
      return modular_addition(a, b, m);
    
  }

  template <typename T> size_t neighbor_index(T x, T y, T z)
  {
    size_t neigh = 0;
    switch(dim_type_)
    {
    case Dimension::Two:
    {
      neigh = (y * dims_[0]) + (x);
      break;
    }
    case Dimension::Three:
    {
      neigh = (z * dims_[1] * dims_[0]) + (y * dims_[0]) + (x);
      break;
    }
    default:
    {
      break;
    }
    }
    return neigh;
  }

  void initialize()
  {
    determine_dimensionality();
    generate_neighborhood();
    kT_ = BOLTZMANN * temperature_;
    total_flips_ = 0;
    seed_ = static_cast<std::mt19937_64::result_type>(std::chrono::steady_clock::now().time_since_epoch().count());
    gen_.seed(seed_);
  }

  void determine_dimensionality()
  {
    SizeVec3Type dims = image_->getDimensions();
    dims_[0] = dims[0];
    dims_[1] = dims[1];
    dims_[2] = dims[2];

    auto num_single_dim = std::count(std::begin(dims_), std::end(dims_), 1);
    assert(num_single_dim <= 1);

    dim_type_ = Dimension::Three;

    for(auto&& dim : dims_)
    {
      if(dim == 1)
      {
        dim_type_ = Dimension::Two;
      }
    }

    if(dim_type_ == Dimension::Two)
    {
      if(dims_[0] == 1)
      {
        dims_[0] = dims_[1];
        dims_[1] = dims_[2];
        dims_[2] = 1;
      }
      else if(dims_[1] == 1)
      {
        dims_[1] = dims_[2];
        dims_[2] = 1;
      }
    }
  }

  void generate_neighborhood()
  {
    switch(dim_type_)
    {
    case Dimension::Two:
    {
      neighborhood_ = {{{1, 0, 0}}, {{-1, 0, 0}}, {{0, 1, 0}}, {{0, -1, 0}}, {{1, 1, 0}}, {{-1, 1, 0}}, {{1, -1, 0}}, {{-1, -1, 0}}};
      break;
    }
    case Dimension::Three:
    {
      neighborhood_ = {{{1, 0, 0}},   {{-1, 0, 0}}, {{0, 1, 0}},  {{0, -1, 0}},  {{0, 0, 1}},   {{0, 0, -1}},  {{1, 1, 0}},   {{-1, 1, 0}},  {{1, -1, 0}},
                       {{-1, -1, 0}}, {{1, 0, 1}},  {{1, 0, -1}}, {{-1, 0, 1}},  {{-1, 0, -1}}, {{0, 1, 1}},   {{0, 1, -1}},  {{0, -1, 1}},  {{0, -1, -1}},
                       {{1, 1, 1}},   {{1, 1, -1}}, {{1, -1, 1}}, {{1, -1, -1}}, {{-1, 1, 1}},  {{-1, 1, -1}}, {{-1, -1, 1}}, {{-1, -1, -1}}};
      break;
    }
    default:
    {
      break;
    }
    }
  }

  ImageGeom::Pointer image_;
  double temperature_;
  double kT_{};
  bool periodic_;
  Neighborhood neighborhood_;
  int32_t* fIds_;
  bool* mask_;
  size_t dims_[3]{};
  Dimension dim_type_;
  size_t total_flips_{};
  std::mt19937_64::result_type seed_{};
  std::mt19937_64 gen_;
};
} // namespace

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
PottsModel::PottsModel()
: m_Iterations(100)
, m_Temperature(273.0)
, m_PeriodicBoundaries(false)
, m_UseMask(false)
, m_FeatureIdsArrayPath(SIMPL::Defaults::ImageDataContainerName, SIMPL::Defaults::CellAttributeMatrixName, SIMPL::CellData::FeatureIds)
, m_MaskArrayPath(SIMPL::Defaults::ImageDataContainerName, SIMPL::Defaults::CellAttributeMatrixName, SIMPL::CellData::Mask)
{
  initialize();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
PottsModel::~PottsModel()
= default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void PottsModel::initialize()
{
  clearErrorCode();
  clearWarningCode();
  setCancel(false);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void PottsModel::setupFilterParameters()
{
  FilterParameterVectorType parameters;
  parameters.push_back(SIMPL_NEW_INTEGER_FP("Iterations", Iterations, FilterParameter::Parameter, PottsModel));
  parameters.push_back(SIMPL_NEW_DOUBLE_FP("Temperature", Temperature, FilterParameter::Parameter, PottsModel));
  parameters.push_back(SIMPL_NEW_BOOL_FP("Periodic Boundaries", PeriodicBoundaries, FilterParameter::Parameter, PottsModel));
  QStringList linkedProps = {"MaskArrayPath"};
  parameters.push_back(SIMPL_NEW_LINKED_BOOL_FP("Use Mask", UseMask, FilterParameter::Parameter, PottsModel, linkedProps));
  parameters.push_back(SeparatorFilterParameter::New("Cell Data", FilterParameter::RequiredArray));
  DataArraySelectionFilterParameter::RequirementType dasReq = DataArraySelectionFilterParameter::CreateRequirement(SIMPL::TypeNames::Int32, 1, AttributeMatrix::Type::Cell, IGeometry::Type::Image);
  parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Feature Ids", FeatureIdsArrayPath, FilterParameter::RequiredArray, PottsModel, dasReq));
  dasReq = DataArraySelectionFilterParameter::CreateRequirement(SIMPL::TypeNames::Bool, 1, AttributeMatrix::Type::Cell, IGeometry::Type::Image);
  parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Mask", MaskArrayPath, FilterParameter::RequiredArray, PottsModel, dasReq));
  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void PottsModel::dataCheck()
{
  initialize();

  if(getIterations() < 1)
  {
    QString ss = QObject::tr("Number of iterations must be greater than 0");
    setErrorCondition(-5555, ss);
  }

  if(getTemperature() <= 0.0f)
  {
    QString ss = QObject::tr("Temperature must be greater than 0");
    setErrorCondition(-5555, ss);
  }

  getDataContainerArray()->getPrereqGeometryFromDataContainer<ImageGeom, AbstractFilter>(this, getFeatureIdsArrayPath().getDataContainerName());

  if(getErrorCode() < 0)
  {
    return;
  }

  QVector<size_t> cDims(1, 1);
  QVector<DataArrayPath> paths;

  m_FeatureIdsPtr = getDataContainerArray()->getPrereqArrayFromPath<DataArray<int32_t>, AbstractFilter>(this, getFeatureIdsArrayPath(), cDims);
  if(m_FeatureIdsPtr.lock())
  {
    m_FeatureIds = m_FeatureIdsPtr.lock()->getPointer(0);
  }
  if(getErrorCode() >= 0)
  {
    paths.push_back(getFeatureIdsArrayPath());
  }

  if(getUseMask())
  {
    getDataContainerArray()->getPrereqGeometryFromDataContainer<ImageGeom, AbstractFilter>(this, getMaskArrayPath().getDataContainerName());

    m_MaskPtr = getDataContainerArray()->getPrereqArrayFromPath<DataArray<bool>, AbstractFilter>(this, getMaskArrayPath(), cDims);
    if(m_MaskPtr.lock())
    {
      m_Mask = m_MaskPtr.lock()->getPointer(0);
    }
    if(getErrorCode() >= 0)
    {
      paths.push_back(getMaskArrayPath());
    }
  }

  getDataContainerArray()->validateNumberOfTuples<AbstractFilter>(this, paths);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void PottsModel::preflight()
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
void PottsModel::execute()
{
  initialize();
  dataCheck();
  if(getErrorCode() < 0)
  {
    return;
  }

  ImageGeom::Pointer image = getDataContainerArray()->getDataContainer(m_FeatureIdsArrayPath.getDataContainerName())->getGeometryAs<ImageGeom>();

  size_t numTuples = m_FeatureIdsPtr.lock()->getNumberOfTuples();
  size_t rangeMin = 0;
  size_t rangeMax = numTuples - 1;
  size_t numNeighbors = 5;
  std::mt19937_64::result_type seed = static_cast<std::mt19937_64::result_type>(std::chrono::steady_clock::now().time_since_epoch().count());
  std::mt19937_64 gen(seed);
  std::uniform_int_distribution<size_t> dist(rangeMin, rangeMax);
  std::uniform_int_distribution<size_t> ndist(rangeMin, numNeighbors);

  if(m_UseMask)
  {
    size_t trues = 0;
    for(size_t i = 0; i < numTuples; i++)
    {
      if(m_Mask[i])
      {
        trues++;
      }
    }
    numTuples = trues;
  }

  SpinLattice lattice(image, m_Temperature, m_PeriodicBoundaries, m_FeatureIds, m_Mask);

  int64_t iter = 0;
  int64_t progPercent = (lattice.dimension() == Dimension::Two) ? 10 : 100;
  int64_t progIncrement = numTuples / progPercent;

  while(iter < m_Iterations)
  {
    if(getCancel())
    {
      return;
    }

    int64_t prog = 1;
    int64_t progressInt = 0;
    int64_t counter = 0;

    for(size_t i = 0; i < numTuples; i++)
    {
      size_t current = 0;

      while(true)
      {
        current = dist(gen);
        if(m_FeatureIds[current] == 0)
        {
          continue;
        }
        if(m_UseMask)
        {
          if(!m_Mask[current])
          {
            continue;
          }
        }
        break;
      }

      lattice.attempt_flip(current);

      if(counter > prog)
      {
        progressInt = static_cast<int64_t>((static_cast<float>(counter) / numTuples) * static_cast<float>(progPercent));
        if(lattice.dimension() == Dimension::Two)
        {
          progressInt *= 10;
        }
        QString ss = QObject::tr("Iteration %1 of %2 || %3% Completed || %4 Total Flips").arg(iter + 1).arg(m_Iterations).arg(progressInt).arg(lattice.total_flips());
        notifyStatusMessage(ss);
        prog = prog + progIncrement;
      }
      counter++;
    }
    iter++;
  }

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AbstractFilter::Pointer PottsModel::newFilterInstance(bool copyFilterParameters) const
{
  PottsModel::Pointer filter = PottsModel::New();
  if(copyFilterParameters)
  {
    copyFilterParameterInstanceVariables(filter.get());
  }
  return filter;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString PottsModel::getCompiledLibraryName() const
{
  return DREAM3DReviewConstants::DREAM3DReviewBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString PottsModel::getBrandingString() const
{
  return "DREAM3DReview";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString PottsModel::getFilterVersion() const
{
  QString version;
  QTextStream vStream(&version);
  vStream << DREAM3DReview::Version::Major() << "." << DREAM3DReview::Version::Minor() << "." << DREAM3DReview::Version::Patch();
  return version;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString PottsModel::getGroupName() const
{
  return DREAM3DReviewConstants::FilterGroups::DREAM3DReviewFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QUuid PottsModel::getUuid()
{
  return QUuid("{e15ec84b-1e02-53a6-a830-59e0813775a1}");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString PottsModel::getSubGroupName() const
{
  return DREAM3DReviewConstants::FilterSubGroups::Coarsening;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString PottsModel::getHumanLabel() const
{
  return "Potts Model";
}
