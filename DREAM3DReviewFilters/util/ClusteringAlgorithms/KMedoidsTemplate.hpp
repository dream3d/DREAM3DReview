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
 * The code contained herein was partially funded by the following contracts:
 *    United States Air Force Prime Contract FA8650-07-D-5800
 *    United States Air Force Prime Contract FA8650-10-D-5210
 *    United States Prime Contract Navy N00173-07-C-2068
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#pragma once

#include <chrono>
#include <random>

#include "SIMPLib/SIMPLib.h"
#include "SIMPLib/DataArrays/DataArray.hpp"
#include "SIMPLib/Filtering/AbstractFilter.h"

#include "DREAM3DReview/DREAM3DReviewFilters/util/DistanceTemplate.hpp"

template <typename T>
class KMedoidsTemplate
{
public:
  using Self = KMedoidsTemplate;
  using Pointer = std::shared_ptr<Self>;
  using ConstPointer = std::shared_ptr<const Self>;
  using WeakPointer = std::weak_ptr<Self>;
  using ConstWeakPointer = std::weak_ptr<const Self>;
  static Pointer NullPointer()
  {
    return Pointer(static_cast<Self*>(nullptr));
  }

  /**
   * @brief Returns the name of the class for KMedoidsTemplate
   */
  /**
   * @brief Returns the name of the class for KMedoidsTemplate
   */
  QString getNameOfClass() const
  {
    return QString("KMedoidsTemplate");
  }

  /**
   * @brief Returns the name of the class for KMedoidsTemplate
   */
  QString ClassName()
  {
    return QString("KMedoidsTemplate");
  }

  KMedoidsTemplate() = default;
  virtual ~KMedoidsTemplate() = default;

  // -----------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------
  bool operator()(IDataArray::Pointer p)
  {
    return (std::dynamic_pointer_cast<DataArray<T>>(p).get() != nullptr);
  }

  // -----------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------
  void Execute(AbstractFilter* filter, IDataArray::Pointer inputIDataArray, IDataArray::Pointer outputIDataArray, BoolArrayType::Pointer maskDataArray, size_t numClusters,
               Int32ArrayType::Pointer fIds, int32_t distMetric, std::pair<bool, uint64_t> randomSeed)
  {
    typename DataArray<T>::Pointer inputDataPtr = std::dynamic_pointer_cast<DataArray<T>>(inputIDataArray);
    typename DataArray<T>::Pointer outputDataPtr = std::dynamic_pointer_cast<DataArray<T>>(outputIDataArray);
    T* inputData = inputDataPtr->getPointer(0);
    T* outputData = outputDataPtr->getPointer(0);

    size_t numTuples = inputDataPtr->getNumberOfTuples();
    int32_t numCompDims = inputDataPtr->getNumberOfComponents();

    size_t rangeMin = 0;
    size_t rangeMax = numTuples - 1;
    std::mt19937_64::result_type seed = static_cast<std::mt19937_64::result_type>(std::chrono::steady_clock::now().time_since_epoch().count());
    if(randomSeed.first)
    {
      seed = static_cast<std::mt19937_64::result_type>(randomSeed.second);
    }
    std::mt19937_64 gen(seed);
    std::uniform_int_distribution<size_t> dist(rangeMin, rangeMax);

    std::vector<size_t> clusterIdxs(numClusters);
    bool* mask = maskDataArray->getPointer(0);
    size_t clusterChoices = 0;

    while(clusterChoices < numClusters)
    {
      size_t index = dist(gen);
      if(mask[index])
      {
        clusterIdxs[clusterChoices] = index;
        clusterChoices++;
      }
    }

    for(size_t i = 0; i < numClusters; i++)
    {
      for(int32_t j = 0; j < numCompDims; j++)
      {
        outputData[numCompDims * (i + 1) + j] = inputData[numCompDims * clusterIdxs[i] + j];
      }
    }

    // size_t updateCheck = 0;
    int32_t* fPtr = fIds->getPointer(0);

    findClusters(filter, mask, inputData, outputData, fPtr, numTuples, numClusters, numCompDims, distMetric);

    std::vector<size_t> optClusterIdxs(clusterIdxs);
    std::vector<double> costs;

    costs = optimizeClusters(filter, mask, inputData, outputData, fPtr, numTuples, numClusters, numCompDims, clusterIdxs, distMetric);

    bool update = optClusterIdxs != clusterIdxs;
    size_t iteration = 1;

    while(update)
    {
      findClusters(filter, mask, inputData, outputData, fPtr, numTuples, numClusters, numCompDims, distMetric);

      optClusterIdxs = clusterIdxs;

      costs = optimizeClusters(filter, mask, inputData, outputData, fPtr, numTuples, numClusters, numCompDims, clusterIdxs, distMetric);

      update = optClusterIdxs != clusterIdxs;

      double sum = std::accumulate(std::begin(costs), std::end(costs), 0.0);
      QString ss = QObject::tr("Clustering Data || Iteration %1 || Total Cost: %2").arg(iteration).arg(sum);
      filter->notifyStatusMessage(ss);
      iteration++;
    }
  }

private:
  // -----------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------
  void findClusters(AbstractFilter* filter, bool* mask, T* input, T* medoids, int32_t* fIds, size_t tuples, int32_t clusters, int32_t dims, int32_t distMetric)
  {
    double dist = 0.0;

    for(size_t i = 0; i < tuples; i++)
    {
      if(filter->getCancel())
      {
        return;
      }
      if(mask[i])
      {
        double minDist = std::numeric_limits<double>::max();
        for(size_t j = 0; j < clusters; j++)
        {
          dist = DistanceTemplate::GetDistance<T, T, double>(input + (dims * i), medoids + (dims * (j + 1)), dims, distMetric);
          if(dist < minDist)
          {
            minDist = dist;
            fIds[i] = j + 1;
          }
        }
      }
    }
  }

  // -----------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------
  std::vector<double> optimizeClusters(AbstractFilter* filter, bool* mask, T* input, T* medoids, int32_t* fIds, size_t tuples, int32_t clusters, int32_t dims, std::vector<size_t>& clusterIdxs,
                                       int32_t distMetric)
  {
    double dist = 0.0;
    std::vector<double> minCosts(clusters, std::numeric_limits<double>::max());

    for(size_t i = 0; i < clusters; i++)
    {
      if(filter->getCancel())
      {
        return std::vector<double>();
      }
      for(size_t j = 0; j < tuples; j++)
      {
        if(filter->getCancel())
        {
          return std::vector<double>();
        }
        if(mask[j])
        {
          double cost = 0.0;
          if(fIds[j] == i + 1)
          {
            for(size_t k = 0; k < tuples; k++)
            {
              if(filter->getCancel())
              {
                return std::vector<double>();
              }
              if(fIds[k] == i + 1 && mask[k])
              {
                dist = DistanceTemplate::GetDistance<T, T, double>(input + (dims * k), input + (dims * j), dims, distMetric);
                cost += dist;
              }
            }

            if(cost < minCosts[i])
            {
              minCosts[i] = cost;
              clusterIdxs[i] = j;
            }
          }
        }
      }
    }

    for(size_t i = 0; i < clusters; i++)
    {
      for(int32_t j = 0; j < dims; j++)
      {
        medoids[dims * (i + 1) + j] = input[dims * clusterIdxs[i] + j];
      }
    }

    return minCosts;
  }

public:
  KMedoidsTemplate(const KMedoidsTemplate&) = delete;            // Copy Constructor Not Implemented
  KMedoidsTemplate(KMedoidsTemplate&&) = delete;                 // Move Constructor Not Implemented
  KMedoidsTemplate& operator=(const KMedoidsTemplate&) = delete; // Copy Assignment Not Implemented
  KMedoidsTemplate& operator=(KMedoidsTemplate&&) = delete;      // Move Assignment Not Implemented
};
