/* ============================================================================
 * Copyright 2021 The University of Utah
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
 * Neither the name of BlueQuartz Software, the US Air Force, the University of Utah nor the names of its contributors may be used to endorse or promote products derived from this software without
 * specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY
 * OF SUCH DAMAGE.
 *
 * The code contained herein was partially funded by the following contracts:
 *
 *
 * This code contained herein is based upon work supported by the following grants:
 *    DOE Office of Nuclear Energy's Nuclear Energy University Program Grant No.: DE-NE0008799
 *    DOD Office of Economic Adjustment Grant No.: ST1605-19-03
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
#pragma once

#include <QtCore/QCoreApplication>
#include <QtCore/QFile>

#include "SIMPLib/SIMPLib.h"
#include "SIMPLib/DataArrays/DataArray.hpp"
#include "SIMPLib/Filtering/FilterFactory.hpp"
#include "SIMPLib/Filtering/FilterManager.h"
#include "SIMPLib/Filtering/FilterPipeline.h"
#include "SIMPLib/Filtering/QMetaObjectUtilities.h"
#include "SIMPLib/Plugin/ISIMPLibPlugin.h"
#include "SIMPLib/Plugin/SIMPLibPluginLoader.h"

#include "UnitTestSupport.hpp"

#include "DREAM3DReviewTestFileLocations.h"

#include "DREAM3DReview/DREAM3DReviewFilters/ComputeFeatureEigenstrains.h"

class ComputeFeatureEigenstrainsTest
{

public:
  ComputeFeatureEigenstrainsTest() = default;
  ~ComputeFeatureEigenstrainsTest() = default;
  ComputeFeatureEigenstrainsTest(const ComputeFeatureEigenstrainsTest&) = delete;            // Copy Constructor
  ComputeFeatureEigenstrainsTest(ComputeFeatureEigenstrainsTest&&) = delete;                 // Move Constructor
  ComputeFeatureEigenstrainsTest& operator=(const ComputeFeatureEigenstrainsTest&) = delete; // Copy Assignment
  ComputeFeatureEigenstrainsTest& operator=(ComputeFeatureEigenstrainsTest&&) = delete;      // Move Assignment

  // -----------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------
  int TestFilterAvailability()
  {
    // Now instantiate the ComputeFeatureEigenstrainsTest Filter from the FilterManager
    QString filtName = "ComputeFeatureEigenstrains";
    FilterManager* fm = FilterManager::Instance();
    IFilterFactory::Pointer filterFactory = fm->getFactoryFromClassName(filtName);
    if(nullptr == filterFactory.get())
    {
      std::stringstream ss;
      ss << "The ComputeFeatureEigenstrainsTest Requires the use of the " << filtName.toStdString() << " filter which is found in the DREAM3DReview Plugin";
      DREAM3D_TEST_THROW_EXCEPTION(ss.str())
    }
    return 0;
  }

  // -----------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------
  int TestFilterInputs()
  {
    bool propWasSet;
    int32_t err;

    ComputeFeatureEigenstrains::Pointer filter = ComputeFeatureEigenstrains::New();

    // Test setting impossible Poisson's ratio
    propWasSet = filter->setProperty("PoissonRatio", 0.5f);
    err = filter->getErrorCode();
    DREAM3D_REQUIRE(err < 0)

    // Test setting possible Poisson's ratio
    propWasSet = filter->setProperty("PoissonRatio", 0.12345f);
    DREAM3D_REQUIRE_EQUAL(propWasSet, true)
    err = filter->getErrorCode();
    DREAM3D_REQUIRE_EQUAL(err, 0)

    // Test setting UseElliposoidalGrains
    propWasSet = filter->setProperty("UseElliposoidalGrains", false);
    DREAM3D_REQUIRE_EQUAL(propWasSet, true)
    err = filter->getErrorCode();
    DREAM3D_REQUIRE_EQUAL(err, 0)

    // Test setting UseCorrectionalMatrix
    propWasSet = filter->setProperty("UseCorrectionalMatrix", true);
    DREAM3D_REQUIRE_EQUAL(propWasSet, true)
    err = filter->getErrorCode();
    DREAM3D_REQUIRE_EQUAL(err, 0)

    return EXIT_SUCCESS;
  }

  // -----------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------
  int TestComputeFeatureEigenstrainsTest()
  {

    return EXIT_SUCCESS;
  }

  // -----------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------
  void operator()()
  {
    int err = EXIT_SUCCESS;

    DREAM3D_REGISTER_TEST(TestFilterAvailability());
    DREAM3D_REGISTER_TEST(TestFilterInputs());
    DREAM3D_REGISTER_TEST(TestComputeFeatureEigenstrainsTest());
  }

private:
};
