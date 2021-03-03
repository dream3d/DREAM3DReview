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
  int TestFilterInputs()
  {
    int32_t err = 0;

    ComputeFeatureEigenstrains::Pointer filter = ComputeFeatureEigenstrains::New();
    // Test setting impossible Poisson's ratio
    filter->setPoissonRatio(0.5F);
    filter->preflight();
    err = filter->getErrorCode();
    DREAM3D_REQUIRED(err, <, 0)

    // Test setting possible Poisson's ratio
    filter->setPoissonRatio(0.12345f);
    filter->preflight();
    err = filter->getErrorCode();
    DREAM3D_REQUIRE_EQUAL(err, -80002)

    // filter->setUseEllipsoidalGrains(true);

    // Test setting UseCorrectionalMatrix
    filter->setUseCorrectionalMatrix(true);
    filter->preflight();
    err = filter->getErrorCode();
    DREAM3D_REQUIRE_EQUAL(err, -80002)

    return EXIT_SUCCESS;
  }

  // -----------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------
  int TestComputeFeatureEigenstrainsTest()
  {
    // Possible tests:
    // Check the following to make sure transition between solution types is valid
    // Eshelby tensor for a=1, b=1, c=1 should be very close to a=1.001, b=1. c=0.999
    // Eshelby tensor for a=3, b=3, c=1 should be very close to a=3, b=2.999, c=1
    // Eshelby tensor for a=3, b=1, c=1 should be very close to a=3, b=1, c=0.999
    // Calculate Eshelby tensor for variety of grain shapes
    // make sure a < b and b < c fails
    // make sure nu > 0.4999 fails
    // check that gauss integration code works for simple analytical functions
    // check that pipeline runs and gives correct result

    return EXIT_SUCCESS;
  }

  // -----------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------
  void operator()()
  {
    int err = EXIT_SUCCESS;

    DREAM3D_REGISTER_TEST(TestFilterInputs());
    DREAM3D_REGISTER_TEST(TestComputeFeatureEigenstrainsTest());
  }

private:
};
