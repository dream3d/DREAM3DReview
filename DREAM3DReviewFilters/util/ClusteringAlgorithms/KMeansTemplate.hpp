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
#include "SIMPLib/Math/SIMPLibMath.h"

#include "DREAM3DReview/DREAM3DReviewFilters/util/DistanceTemplate.hpp"

template <typename T>
class KMeansTemplate
{
public:
  using Self = KMeansTemplate;
  using Pointer = std::shared_ptr<Self>;
  using ConstPointer = std::shared_ptr<const Self>;
  using WeakPointer = std::weak_ptr<Self>;
  using ConstWeakPointer = std::weak_ptr<const Self>;
  static Pointer NullPointer()
  {
    return Pointer(static_cast<Self*>(nullptr));
  }

  /**
   * @brief Returns the name of the class for KMeansTemplate
   */
  /**
   * @brief Returns the name of the class for KMeansTemplate
   */
  QString getNameOfClass() const
  {
    return QString("KMeansTemplate");
  }

  /**
   * @brief Returns the name of the class for KMeansTemplate
   */
  QString ClassName()
  {
    return QString("KMeansTemplate");
  }

  KMeansTemplate() = default;
  virtual ~KMeansTemplate() = default;

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
  void Execute(AbstractFilter* filter, IDataArray::Pointer inputIDataArray, DoubleArrayType::Pointer outputDataArray, BoolArrayType::Pointer maskDataArray, size_t numClusters,
               Int32ArrayType::Pointer fIds, int distMetric, std::pair<bool, uint64_t> randomSeed)
  {
    typename DataArray<T>::Pointer inputDataPtr = std::dynamic_pointer_cast<DataArray<T>>(inputIDataArray);
    T* inputData = inputDataPtr->getPointer(0);
    double* outputData = outputDataArray->getPointer(0);

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

    std::vector<size_t> initClusterIdxs(numClusters);
    bool* mask = maskDataArray->getPointer(0);
    size_t clusterChoices = 0;

    while(clusterChoices < numClusters)
    {
      size_t index = dist(gen);
      if(mask[index])
      {
        initClusterIdxs[clusterChoices] = index;
        clusterChoices++;
      }
    }

    for(size_t i = 0; i < numClusters; i++)
    {
      for(int32_t j = 0; j < numCompDims; j++)
      {
        outputData[numCompDims * (i + 1) + j] = inputData[numCompDims * initClusterIdxs[i] + j];
      }
    }

    size_t updateCheck = 0;
    std::vector<double> oldMeans(numClusters);
    std::vector<double> differences(numClusters);
    int32_t* fPtr = fIds->getPointer(0);
    size_t iteration = 1;

    while(true)
    {
      if(filter->getCancel())
      {
        return;
      }
      findClusters(filter, mask, inputData, outputData, fPtr, numTuples, numClusters, numCompDims, distMetric);

      for(size_t i = 0; i < numClusters; i++)
      {
        oldMeans[i] = outputData[i + 1];
      }

      findMeans(mask, inputData, outputData, fPtr, numTuples, numClusters, numCompDims);

      updateCheck = 0;
      for(size_t i = 0; i < numClusters; i++)
      {
        differences[i] = oldMeans[i] - outputData[i + 1];
        if(SIMPLibMath::closeEnough<double>(differences[i], 0.0))
        {
          updateCheck++;
        }
      }

      double sum = std::accumulate(std::begin(differences), std::end(differences), 0.0);
      QString ss = QObject::tr("Clustering Data || Iteration %1 || Total Mean Shift: %2").arg(iteration).arg(sum);
      filter->notifyStatusMessage(ss);
      iteration++;

      if(updateCheck == numClusters)
      {
        break;
      }
    }
  }

private:
  // -----------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------
  void findClusters(AbstractFilter* filter, bool* mask, T* input, double* averages, int32_t* fIds, size_t tuples, int32_t clusters, int32_t dims, int32_t distMetric)
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
          dist = DistanceTemplate::GetDistance<T, double, double>(input + (dims * i), averages + (dims * (j + 1)), dims, distMetric);
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
  void findMeans(bool* mask, T* input, double* averages, int32_t* fIds, size_t tuples, int32_t clusters, int32_t dims)
  {
    double value = 0.0;
    int32_t feature = 0;
    std::vector<size_t> counts(clusters + 1, 0);

    for(size_t i = 0; i <= clusters; i++)
    {
      for(size_t j = 0; j < dims; j++)
      {
        averages[dims * i + j] = 0.0;
      }
    }

    for(size_t i = 0; i < dims; i++)
    {
      for(size_t j = 0; j < tuples; j++)
      {
        value = static_cast<double>(input[dims * j + i]);
        feature = fIds[j];
        averages[dims * feature + i] += value;
        counts[feature] += 1;
      }
      for(size_t j = 0; j <= clusters; j++)
      {
        if(counts[j] == 0)
        {
          averages[dims * j + i] = 0.0;
        }
        else
        {
          averages[dims * j + i] /= static_cast<double>(counts[j]);
        }
      }
      std::fill(std::begin(counts), std::end(counts), 0);
    }
  }

  KMeansTemplate(const KMeansTemplate&); // Copy Constructor Not Implemented
  void operator=(const KMeansTemplate&); // Move assignment Not Implemented
};
