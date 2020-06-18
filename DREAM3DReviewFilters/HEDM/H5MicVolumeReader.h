/* ============================================================================
 * Copyright (c) 2010, Michael A. Jackson (BlueQuartz Software)
 * Copyright (c) 2010, Dr. Michael A. Groeber (US Air Force Research Laboratories)
 * All rights reserved.
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
 * Neither the name of Michael A. Groeber, Michael A. Jackson, the US Air Force,
 * BlueQuartz Software nor the names of its contributors may be used to endorse
 * or promote products derived from this software without specific prior written
 * permission.
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
 *  This code was written under United States Air Force Contract number
 *                           FA8650-07-D-5800
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#pragma once

#include <vector>

//-- Ebsd Lib Includes
#include "EbsdLib/Core/EbsdLibConstants.h"
#include "EbsdLib/EbsdLib.h"
#include "EbsdLib/Core/EbsdSetGetMacros.h"
#include "EbsdLib/IO/H5EbsdVolumeReader.h"

#include "DREAM3DReview/DREAM3DReviewFilters/HEDM/MicPhase.h"

/**
 * @class H5MicVolumeReader H5MicVolumeReader.h Reconstruction/EbsdSupport/H5MicVolumeReader.h
 * @brief This class loads OIM data from an HDF5 based file.
 *
 * @date May 23, 2011
 * @version 1.0
 */
class H5MicVolumeReader : public H5EbsdVolumeReader
{
public:
  EBSD_SHARED_POINTERS(H5MicVolumeReader)
  EBSD_STATIC_NEW_SUPERCLASS(H5EbsdVolumeReader, H5MicVolumeReader)

  ~H5MicVolumeReader() override;

  EBSD_POINTER_PROPERTY(Euler1, Euler1, float)
  EBSD_POINTER_PROPERTY(Euler2, Euler2, float)
  EBSD_POINTER_PROPERTY(Euler3, Euler3, float)
  EBSD_POINTER_PROPERTY(X, X, float)
  EBSD_POINTER_PROPERTY(Y, Y, float)
  EBSD_POINTER_PROPERTY(Confidence, Conf, float)
  EBSD_POINTER_PROPERTY(Phase, Phase, int)

  /**
   * @brief This method does the actual loading of the OIM data from the data
   * source (files, streams, etc) into the data structures.
   * @param eulerMicles
   * @param phases
   * @param xpoints
   * @param ypoints
   * @param zpoints
   * @param filters
   * @return
   */
  int loadData(int64_t xpoints, int64_t ypoints, int64_t zpoints, uint32_t ZDir) override;

  /**
   * @brief
   * @return
   */
  std::vector<MicPhase::Pointer> getPhases();
  /**
   * @brief Returns the pointer to the data for a given feature
   * @param featureName The name of the feature to return the pointer to.
   */
  void* getPointerByName(const std::string& featureName) override;

  /**
   * @brief Returns an enumeration value that depicts the numerical
   * primitive type that the data is stored as (Int, Float, etc).
   * @param featureName The name of the feature.
   */
  EbsdLib::NumericTypes::Type getPointerType(const std::string& featureName) override;

  /** @brief Allocates the proper amount of memory (after reading the header portion of the file)
   * and then splats '0' across all the bytes of the memory allocation
   */
  void initPointers(size_t numElements) override;

  /** @brief 'free's the allocated memory and sets the pointer to nullptr
   */
  void deletePointers() override;

protected:
  H5MicVolumeReader();

private:
  std::vector<MicPhase::Pointer> m_Phases;

  /**
   * @brief Allocats a contiguous chunk of memory to store values from the .Mic file
   * @param numberOfElements The number of elements in the Array. This method can
   * also optionally produce SSE aligned memory for use with SSE intrinsics
   * @return Pointer to allocated memory
   */
  template <typename T>
  T* allocateArray(size_t numberOfElements)
  {
    T* buffer = nullptr;
    if(numberOfElements == 0)
    {
      return buffer;
    }
    buffer = new T[numberOfElements]();

    return buffer;
  }

  /**
   * @brief Deallocates memory that has been previously allocated. This will set the
   * value of the pointer passed in as the argument to nullptr.
   * @param ptr The pointer to be freed.
   */
  template <typename T>
  void deallocateArrayData(T*& ptr)
  {
    if(ptr != nullptr && getManageMemory())
    {
      delete[] ptr;
      ptr = nullptr;
    }
  }

public:
  H5MicVolumeReader(const H5MicVolumeReader&) = delete;            // Copy Constructor Not Implemented
  H5MicVolumeReader(H5MicVolumeReader&&) = delete;                 // Move Constructor Not Implemented
  H5MicVolumeReader& operator=(const H5MicVolumeReader&) = delete; // Copy Assignment Not Implemented
  H5MicVolumeReader& operator=(H5MicVolumeReader&&) = delete;      // Move Assignment Not Implemented
};
