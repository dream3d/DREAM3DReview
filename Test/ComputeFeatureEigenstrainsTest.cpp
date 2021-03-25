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
#include <QtCore/QFileInfo>
#include <QtCore/QString>

#include "SIMPLib/SIMPLib.h"
#include "SIMPLib/DataArrays/DataArray.hpp"
#include "SIMPLib/DataContainers/DataContainerArray.h"
#include "SIMPLib/FilterParameters/JsonFilterParametersReader.h"
#include "SIMPLib/Filtering/FilterFactory.hpp"
#include "SIMPLib/Filtering/FilterManager.h"
#include "SIMPLib/Filtering/FilterPipeline.h"
#include "SIMPLib/Plugin/ISIMPLibPlugin.h"
#include "SIMPLib/CoreFilters/DataContainerReader.h"

#include "UnitTestSupport.hpp"

#include "DREAM3DReviewTestFileLocations.h"

#include "DREAM3DReview/DREAM3DReviewFilters/ComputeFeatureEigenstrains.h"
#include "DREAM3DReview/DREAM3DReviewFilters/util/EigenstrainsHelper.hpp"
#include "SIMPLib/CoreFilters/DataContainerReader.h"

namespace SIMPLMath = SIMPLib::Constants;

namespace
{
void compare_tensors(EigenstrainsHelper::Tensor4DType L, EigenstrainsHelper::Tensor4DType R, double eps)
{
  for(size_t i = 0; i < 3; i++)
  {
    for(size_t j = 0; j < 3; j++)
    {
      for(size_t k = 0; k < 3; k++)
      {
        for(size_t l = 0; l < 3; l++)
        {
          if(!SIMPLibMath::closeEnough<>(L(i, j, k, l), R(i, j, k, l), eps))
          {
            QString buf;
            QTextStream ss(&buf);
            ss << "Your test required the following\n            '";
            ss << "SIMPLibMath::closeEnough<>(" << L(i, j, k, l) << ", " << R(i, j, k, l) << ", " << eps << "'\n             but this condition was not met with eps=" << eps << "\n";
            ss << "             S_" << i + 1 << j + 1 << k + 1 << l + 1 << ": " << L(i, j, k, l) << "==" << R(i, j, k, l);
            DREAM3D_TEST_THROW_EXCEPTION(buf.toStdString());
          }
        }
      }
    }
  }
}
} // namespace

class ComputeFeatureEigenstrainsTest
{
  const QString k_SyntheticVolumeDataContainer = {"SyntheticVolumeDataContainer"};
  const QString k_GrainData = {"Grain Data"};
  const QString k_Eigenstrains = {"Eigenstrains"};

public:
  ComputeFeatureEigenstrainsTest() = default;
  ~ComputeFeatureEigenstrainsTest() = default;
  ComputeFeatureEigenstrainsTest(const ComputeFeatureEigenstrainsTest&) = delete;            // Copy Constructor
  ComputeFeatureEigenstrainsTest(ComputeFeatureEigenstrainsTest&&) = delete;                 // Move Constructor
  ComputeFeatureEigenstrainsTest& operator=(const ComputeFeatureEigenstrainsTest&) = delete; // Copy Assignment
  ComputeFeatureEigenstrainsTest& operator=(ComputeFeatureEigenstrainsTest&&) = delete;      // Move Assignment


  // -----------------------------------------------------------------------------
  // Test compute feature eigenstrains example pipeline and compare the eigenstrains to a reference
  // -----------------------------------------------------------------------------
  int TestComputeFeatureEigenstrainsTest()
  {
    // ComputeFeatureEigenstrains example pipeline
    QString pipelineFile = UnitTest::PluginSourceDir + "/ExamplePipelines/" + UnitTest::PluginName + "/ComputeFeatureEigenstrainsExample.json";
    QFileInfo fi(pipelineFile);
    if(!fi.exists())
    {
      std::cout << "Compute feature eigenstrains example pipeline path: '" << pipelineFile.toStdString() << "' does not exist" << std::endl;
      return EXIT_FAILURE;
    }

    // Read example pipeline
    FilterPipeline::Pointer pipeline;
    JsonFilterParametersReader::Pointer jsonReader = JsonFilterParametersReader::New();
    pipeline = jsonReader->readPipelineFromFile(pipelineFile);
    if(nullptr == pipeline.get())
    {
      std::cout << "An error occurred trying to read the pipeline file. Exiting now." << std::endl;
      return EXIT_FAILURE;
    }

    // Create an Observer to report errors/progress from the executing pipeline
    std::cout << "Pipeline Count: " << pipeline->size() << std::endl;
    Observer obs;
    pipeline->addMessageReceiver(&obs);

    //Ensure we can find the input data because we don't know what directory the test will be run from.
    QString inputFilePath;
    QTextStream ss(&inputFilePath);
    ss << UnitTest::PluginSourceDir << "/Data/DREAM3DReview/SyntheticComputeEigenstrainData.dream3d";
    DataContainerReader::Pointer dataContainerReader = std::dynamic_pointer_cast<DataContainerReader>(pipeline->getFilterContainer()[0]);
    dataContainerReader->setInputFile(inputFilePath);

    // Preflight the pipeline
    int32_t err = pipeline->preflightPipeline();
    if(err < 0)
    {
      std::cout << "Errors preflighting the pipeline. Exiting Now." << std::endl;
    }
    DREAM3D_REQUIRE(err >= 0)

    // Execute the pipeline
    DataContainerArray::Pointer dca = pipeline->execute();
    err = pipeline->getErrorCode();
    if(err < 0)
    {
      std::cout << "Error Condition of Pipeline: " << err << std::endl;
    }
    DREAM3D_REQUIRE(err >= 0)

    // Get eigenstrains from output
    DataContainer::Pointer dc = dca->getDataContainer(k_SyntheticVolumeDataContainer);
    AttributeMatrix::Pointer featureAttrMat = dc->getAttributeMatrix(k_GrainData);
    FloatArrayType::Pointer EigenstrainsArray = std::dynamic_pointer_cast<FloatArrayType>(featureAttrMat->getAttributeArray(k_Eigenstrains));

    // Read reference eigenstrain HDF5 file
    QString referenceFile = "Data/DREAM3DReview/ReferenceOutputSyntheticComputeEigenstrainData.dream3d";
    DataContainerArray::Pointer dca2 = DataContainerArray::New();
    DataContainerReader::Pointer reader2 = DataContainerReader::New();
    reader2->setInputFile(referenceFile);
    reader2->setDataContainerArray(dca2);

    DataContainerArrayProxy dcaProxy = reader2->readDataContainerArrayStructure(referenceFile);
    reader2->setInputFileDataContainerArrayProxy(dcaProxy);
    reader2->execute();
    AttributeMatrix::Pointer featureAttrMatRef = dc->getAttributeMatrix(k_GrainData);
    FloatArrayType::Pointer EigenstrainsArrayRef = std::dynamic_pointer_cast<FloatArrayType>(featureAttrMat->getAttributeArray(k_Eigenstrains));

    // Compare calculated eigenstrains to reference eigenstrains
    size_t numTuples = EigenstrainsArray->getNumberOfTuples();
    size_t numTuplesRef = EigenstrainsArray->getNumberOfTuples();
    DREAM3D_REQUIRE_EQUAL(numTuples, numTuplesRef)

    FloatArrayType& eigenstrains = *(EigenstrainsArray.get());
    FloatArrayType& eigenstrainsRef = *(EigenstrainsArrayRef.get());
    double eps = 1e-5;
    double calc, ref, delta;
    for(size_t feature = 0; feature < numTuples; feature++)
    {
      for(size_t idx = 0; idx < 6; idx++)
      {
        calc = eigenstrains[feature * 6 + idx];
        ref = eigenstrainsRef[feature * 6 + idx];
        delta = std::abs(calc - ref);
        if(isnan(delta))
        {
          if(idx == 0)
          {
            std::cout << "Warning, featureID '" << feature << "' has NaN eigenstrains" << std::endl;
          }
        }
        else
        {
          DREAM3D_REQUIRED(delta, <, eps);
        }
      }
    }

    return EXIT_SUCCESS;
  }

  // -----------------------------------------------------------------------------
  // Test 32-point gauss integration on several functions
  // -----------------------------------------------------------------------------
  int GaussIntegrationTest()
  {
    double eps = 1e-6;
    double trueValue, numericalValue;

    // 1.) Integrate x^2
    auto function = [](double x) { return std::pow(x, 2); };

    // 0 to 1
    trueValue = 1.0 / 3;
    numericalValue = EigenstrainsHelper::gauss_integration(function, 0, 1);
    DREAM3D_REQUIRED(std::abs(numericalValue - trueValue), <, eps);

    // 0 to 4
    trueValue = 64.0 / 3;
    numericalValue = EigenstrainsHelper::gauss_integration(function, 0, 4);
    DREAM3D_REQUIRED(std::abs(numericalValue - trueValue), <, eps);

    // 3 to 7
    trueValue = 316.0 / 3;
    numericalValue = EigenstrainsHelper::gauss_integration(function, 3, 7);
    DREAM3D_REQUIRED(std::abs(numericalValue - trueValue), <, eps);

    // -4 to 5
    trueValue = 63.0;
    numericalValue = EigenstrainsHelper::gauss_integration(function, -4, 5);
    DREAM3D_REQUIRED(std::abs(numericalValue - trueValue), <, eps);

    // 2.) Integrate sin(x)
    auto function2 = [](double x) { return std::sin(x); };

    // 0 to pi
    trueValue = 2.0;
    numericalValue = EigenstrainsHelper::gauss_integration(function2, 0, SIMPLMath::k_PiD);
    DREAM3D_REQUIRED(std::abs(numericalValue - trueValue), <, eps);

    // -pi/4 to pi/4
    trueValue = 0.0;
    numericalValue = EigenstrainsHelper::gauss_integration(function2, -SIMPLMath::k_PiD / 4, SIMPLMath::k_PiD / 4);
    DREAM3D_REQUIRED(std::abs(numericalValue - trueValue), <, eps);

    // 3.) Integrate x / (1 + x)
    auto function3 = [](double x) { return x / (1 + x); };

    // 1 to 3
    trueValue = 10.0 - std::log(11);
    numericalValue = EigenstrainsHelper::gauss_integration(function3, 0, 10);
    DREAM3D_REQUIRED(std::abs(numericalValue - trueValue), <, eps);

    // 1 to 3
    trueValue = 2.0 - std::log(2);
    numericalValue = EigenstrainsHelper::gauss_integration(function3, 1, 3);
    DREAM3D_REQUIRED(std::abs(numericalValue - trueValue), <, eps);

    return EXIT_SUCCESS;
  }

  // -----------------------------------------------------------------------------
  // Test Eshelby tensor calculations to verify continuity and accuracy
  // -----------------------------------------------------------------------------
  int FindEshelbyTest()
  {
    EigenstrainsHelper::Tensor4DType eshelbyTensor1, eshelbyTensor2;
    double eps = 1e-3;
    double nu1 = 0.15;
    double nu2 = 0.33;
    double nu3 = 0.40;

    // Check that ellipsoid solution converges to spherical solution
    eshelbyTensor1 = EigenstrainsHelper::find_eshelby(1.0, 1.0, 1.0, nu1, true);
    eshelbyTensor2 = EigenstrainsHelper::find_eshelby(1.0 + eps, 1.0, 1.0 - eps, nu1, true);
    ::compare_tensors(eshelbyTensor1, eshelbyTensor2, eps);

    eshelbyTensor1 = EigenstrainsHelper::find_eshelby(1.0, 1.0, 1.0, nu2, true);
    eshelbyTensor2 = EigenstrainsHelper::find_eshelby(1.0 + eps, 1.0, 1.0 - eps, nu2, true);
    ::compare_tensors(eshelbyTensor1, eshelbyTensor2, eps);

    eshelbyTensor1 = EigenstrainsHelper::find_eshelby(1.0, 1.0, 1.0, nu3, true);
    eshelbyTensor2 = EigenstrainsHelper::find_eshelby(1.0 + eps, 1.0, 1.0 - eps, nu3, true);
    ::compare_tensors(eshelbyTensor1, eshelbyTensor2, eps);

    // Check that ellipsoid solution converges to oblate spheroid solution
    eshelbyTensor1 = EigenstrainsHelper::find_eshelby(2.0, 2.0, 1.0, nu1, true);
    eshelbyTensor2 = EigenstrainsHelper::find_eshelby(2.0 + eps, 2.0, 1.0, nu1, true);
    ::compare_tensors(eshelbyTensor1, eshelbyTensor2, eps);

    eshelbyTensor1 = EigenstrainsHelper::find_eshelby(2.0, 2.0, 1.0, nu2, true);
    eshelbyTensor2 = EigenstrainsHelper::find_eshelby(2.0 + eps, 2.0, 1.0, nu2, true);
    ::compare_tensors(eshelbyTensor1, eshelbyTensor2, eps);

    eshelbyTensor1 = EigenstrainsHelper::find_eshelby(2.0, 2.0, 1.0, nu3, true);
    eshelbyTensor2 = EigenstrainsHelper::find_eshelby(2.0 + eps, 2.0, 1.0, nu3, true);
    ::compare_tensors(eshelbyTensor1, eshelbyTensor2, eps);

    // Check that ellipsoid solution converges to prolate spheroid solution
    eshelbyTensor1 = EigenstrainsHelper::find_eshelby(2.0, 1.0, 1.0, nu1, true);
    eshelbyTensor2 = EigenstrainsHelper::find_eshelby(2.0, 1.0, 1.0 - eps, nu1, true);
    ::compare_tensors(eshelbyTensor1, eshelbyTensor2, eps);

    eshelbyTensor1 = EigenstrainsHelper::find_eshelby(2.0, 1.0, 1.0, nu2, true);
    eshelbyTensor2 = EigenstrainsHelper::find_eshelby(2.0, 1.0, 1.0 - eps, nu2, true);
    ::compare_tensors(eshelbyTensor1, eshelbyTensor2, eps);

    eshelbyTensor1 = EigenstrainsHelper::find_eshelby(2.0, 1.0, 1.0, nu3, true);
    eshelbyTensor2 = EigenstrainsHelper::find_eshelby(2.0, 1.0, 1.0 - eps, nu3, true);
    ::compare_tensors(eshelbyTensor1, eshelbyTensor2, eps);

    double S1111, S1122, S1313;

    // Test simple spherical case, probe for specific indices
    eps = 1e-5;
    eshelbyTensor1 = EigenstrainsHelper::find_eshelby(1.0, 1.0, 1.0, 0, true);
    S1111 = 7.0 / 15;
    DREAM3D_REQUIRED(std::abs(S1111 - eshelbyTensor1(0, 0, 0, 0)), <, eps);
    S1122 = -1.0 / 15;
    DREAM3D_REQUIRED(std::abs(S1122 - eshelbyTensor1(0, 0, 1, 1)), <, eps);
    S1313 = 4.0 / 15;
    DREAM3D_REQUIRED(std::abs(S1313 - eshelbyTensor1(0, 2, 0, 2)), <, eps);

    // Test elliptical values, probe for specific indices
    eshelbyTensor1 = EigenstrainsHelper::find_eshelby(2.0, 1.0, 0.5, 0.3, true);
    S1111 = 0.20843827;
    DREAM3D_REQUIRED(std::abs(S1111 - eshelbyTensor1(0, 0, 0, 0)), <, eps);
    S1122 = 0.00895464;
    DREAM3D_REQUIRED(std::abs(S1122 - eshelbyTensor1(0, 0, 1, 1)), <, eps);
    S1313 = 0.30071747;
    DREAM3D_REQUIRED(std::abs(S1313 - eshelbyTensor1(0, 2, 0, 2)), <, eps);

    return EXIT_SUCCESS;
  }

  // -----------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------
  void operator()()
  {
    std::cout << "###### ComputeFeatureEigenstrainsTest ######" << std::endl;
    int err = EXIT_SUCCESS;

    DREAM3D_REGISTER_TEST(TestComputeFeatureEigenstrainsTest());
    DREAM3D_REGISTER_TEST(GaussIntegrationTest());
    DREAM3D_REGISTER_TEST(FindEshelbyTest());
  }

private:
};
