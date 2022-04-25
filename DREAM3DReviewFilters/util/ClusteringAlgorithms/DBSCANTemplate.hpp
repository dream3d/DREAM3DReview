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

#pragma once


#ifdef SIMPL_USE_PARALLEL_ALGORITHMS
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/partitioner.h>
#endif

#include "SIMPLib/Common/SIMPLibSetGetMacros.h"
#include "SIMPLib/Filtering/AbstractFilter.h"
#include "SIMPLib/Filtering/AbstractFilter.h"
#include "SIMPLib/SIMPLib.h"
#include "SIMPLib/SIMPLib.h"

#include "DREAM3DReview/DREAM3DReviewFilters/util/DistanceTemplate.hpp"

template<typename T>
class FindEpsilonNeighborhoodsImpl
{
  public:
  FindEpsilonNeighborhoodsImpl(AbstractFilter* filter, double epsilon, T* inputData, bool* mask,
                               size_t numCompDims, size_t numTuples, int32_t distMetric, std::vector<std::list<size_t>>& neighborhoods)
  : m_Filter(filter)
  , m_Epsilon(epsilon)
  , m_InputData(inputData)
  , m_Mask(mask)
  , m_NumCompDims(numCompDims)
  , m_NumTuples(numTuples)
  , m_DistMetric(distMetric)
  , m_Neighborhoods(neighborhoods)
  {}

  void compute(size_t start, size_t end) const
  {
    // int64_t counter = 0;
    // int64_t totalElements = end - start;
    // int64_t progIncrement = static_cast<int64_t>(totalElements / 100);

    for(size_t i = start; i < end; i++)
    {
      if(m_Mask[i])
      {
        std::list<size_t> neighbors = epsilon_neighbors(i);
        m_Neighborhoods[i] = neighbors;
      }
    }
  }

  std::list<size_t> epsilon_neighbors(size_t index) const
  {
    std::list<size_t> neighbors;

    for(size_t i = 0; i < m_NumTuples; i++)
    {
      if(m_Filter->getCancel()) { return std::list<size_t>(); }
      if(m_Mask[i])
      {
        double dist = DistanceTemplate::GetDistance<T, T, double>(m_InputData + (m_NumCompDims * index),
                                                                  m_InputData + (m_NumCompDims * i),
                                                                  m_NumCompDims,
                                                                  m_DistMetric);
        if(dist < m_Epsilon) { neighbors.push_back(i); }
      }
    }

    return neighbors;
  }

#ifdef SIMPL_USE_PARALLEL_ALGORITHMS
  void operator()(const tbb::blocked_range<size_t>& r) const
  {
    compute(r.begin(), r.end());
  }
#endif

  private:
  AbstractFilter* m_Filter;
  double m_Epsilon;
  T* m_InputData;
  bool* m_Mask;
  size_t m_NumCompDims;
  size_t m_NumTuples;
  int32_t m_DistMetric;
  std::vector<std::list<size_t>>& m_Neighborhoods;
};

template<typename T>
class DBSCANTemplate
{
public:
  SIMPL_SHARED_POINTERS(DBSCANTemplate)
  SIMPL_TYPE_MACRO(DBSCANTemplate)

  DBSCANTemplate() {}
  virtual ~DBSCANTemplate() {}

  // -----------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------
  bool operator()(IDataArray::Pointer p)
  {
    return (std::dynamic_pointer_cast<DataArray<T>>(p) != nullptr);
  }

  // -----------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------
  void Execute(AbstractFilter* filter, IDataArray::Pointer inputIDataArray, BoolArrayType::Pointer maskDataArray, 
               Int32ArrayType::Pointer fIds, float epsilon, int32_t minPnts, int32_t distMetric)
  {
    typename DataArray<T>::Pointer inputDataPtr = std::dynamic_pointer_cast<DataArray<T>>(inputIDataArray);
    T* inputData = inputDataPtr->getPointer(0);
    int32_t* fPtr = fIds->getPointer(0);

    size_t numTuples = inputDataPtr->getNumberOfTuples();
    size_t numCompDims = inputDataPtr->getNumberOfComponents();
    std::vector<bool> visited(numTuples, false);
    std::vector<bool> clustered(numTuples, false);

    double minDist = static_cast<double>(epsilon);
    int32_t cluster = 0;
    bool* mask = maskDataArray->getPointer(0);

    int64_t progIncrement = static_cast<int64_t>(numTuples / 100);
    int64_t prog = 1;
    int64_t progressInt = 0;
    int64_t counter = 0;

    std::vector<std::list<size_t>> epsilonNeighborhoods(numTuples);
#ifdef SIMPL_USE_PARALLEL_ALGORITHMS
        bool doParallel = true;
#endif

#ifdef SIMPL_USE_PARALLEL_ALGORITHMS
    if(doParallel == true)
    {
      tbb::parallel_for(tbb::blocked_range<size_t>(0, numTuples), FindEpsilonNeighborhoodsImpl<T>(filter, minDist, inputData, mask, 
        numCompDims, numTuples, distMetric, epsilonNeighborhoods), tbb::auto_partitioner());
    }
    else
#endif
    {
      FindEpsilonNeighborhoodsImpl<T> serial(filter, minDist, inputData, mask, 
                                             numCompDims, numTuples, distMetric, epsilonNeighborhoods);
      serial.compute(0, numTuples);
    }

    prog = 1;
    progressInt = 0;
    counter = 0;

    for(size_t i = 0; i < numTuples; i++)
    {
      if(filter->getCancel()) { return; }
      if(mask[i] && !visited[i])
      {
        visited[i] = true;

        if(counter > prog)
        {
          progressInt = static_cast<int64_t>((static_cast<float>(counter) / numTuples) * 100.0f);
          QString ss = QObject::tr("Scanning Data || Visited Point %1 of %2 || %3% Completed").arg(counter).arg(numTuples).arg(progressInt);
          filter->notifyStatusMessage(filter->getMessagePrefix(), filter->getHumanLabel(), ss);
          prog = prog + progIncrement;
        }
        counter++;

        std::list<size_t> neighbors = epsilonNeighborhoods[i];

        if(static_cast<int32_t>(neighbors.size()) < minPnts)
        {
          fPtr[i] = 0;
          clustered[i] = true;
        }
        else
        {
          cluster++;
          expand_cluster(filter, neighbors, fPtr, cluster, minPnts, visited, clustered,
                         minDist, i, inputData, mask,
                         numCompDims, numTuples, distMetric,
                         progIncrement, prog, progressInt, counter, epsilonNeighborhoods);
        }
      }
    }
  }

private:
  // -----------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------
  std::list<size_t> epsilon_neighbors(AbstractFilter* filter, double eps, size_t index, T* data, bool* mask, 
                                      size_t dims, size_t numTuples, int32_t metric)
  {
    std::list<size_t> neighbors;

    for(size_t i = 0; i < numTuples; i++)
    {
      if(filter->getCancel()) { return std::list<size_t>(); }
      if(mask[i])
      {
        double dist = DistanceTemplate::GetDistance<T, T, double>(data + (dims * index), 
                                                                  data + (dims * i),
                                                                  dims,
                                                                  metric);
        if(dist < eps) { neighbors.push_back(i); }
      }
    }

    return neighbors;
  }

  // -----------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------
  void expand_cluster(AbstractFilter* filter, std::list<size_t>& neighbors, int32_t* features, int32_t cluster,
                      int32_t minPnts, std::vector<bool>& visited, std::vector<bool>& clustered,
                      double eps, size_t index, T* data, bool* mask,
                      size_t dims, size_t numTuples, int32_t metric,
                      int64_t& progIncrement, int64_t& prog, int64_t& progressInt, int64_t& counter, std::vector<std::list<size_t>>& epsNeighbors)
  {
    features[index] = cluster;
    clustered[index] = true;

    for(auto&& idx : neighbors)
    {
      if(filter->getCancel()) { return; }
      if(mask[idx])
      {
        if(!visited[idx])
        {
          visited[idx] = true;

          if(counter > prog)
          {
            progressInt = static_cast<int64_t>((static_cast<float>(counter) / numTuples) * 100.0f);
            QString ss = QObject::tr("Scanning Data || Visited Point %1 of %2 || %3% Completed").arg(counter).arg(numTuples).arg(progressInt);
            filter->notifyStatusMessage(filter->getMessagePrefix(), filter->getHumanLabel(), ss);
            prog = prog + progIncrement;
          }
          counter++;

          std::list<size_t> neighbors_prime = epsNeighbors[idx];
          if(static_cast<int32_t>(neighbors_prime.size()) >= minPnts)
          {
            neighbors.splice(std::end(neighbors), neighbors_prime);
          }
        }
        if(!clustered[idx])
        {
          features[idx] = cluster;
          clustered[idx] = true;
        }
      }
    }
  }
  
  DBSCANTemplate(const DBSCANTemplate&); // Copy Constructor Not Implemented
  void operator=(const DBSCANTemplate&); // Move assignment Not Implemented
};

