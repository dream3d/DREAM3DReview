/* ============================================================================
 * Copyright (c) 2019-2019 BlueQuartz Software, LLC
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
 *    United States Air Force Prime Contract FA8650-10-D-5210
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include <QtCore/QFile>
#include <QtCore/QFileInfo>
#include <QtCore/QString>

#include "SIMPLib/Common/SIMPLibSetGetMacros.h"
#include "SIMPLib/DataArrays/DataArray.hpp"
#include "SIMPLib/FilterParameters/JsonFilterParametersReader.h"
#include "SIMPLib/Filtering/FilterFactory.hpp"
#include "SIMPLib/Filtering/FilterManager.h"
#include "SIMPLib/Filtering/FilterPipeline.h"
#include "SIMPLib/Filtering/QMetaObjectUtilities.h"
#include "SIMPLib/Plugin/ISIMPLibPlugin.h"
#include "SIMPLib/Plugin/SIMPLibPluginLoader.h"
#include "SIMPLib/SIMPLib.h"

#include "UnitTestSupport.hpp"

#include "DREAM3DReviewTestFileLocations.h"

// These values are dependant on the prebuilt example pipeline for this filter. If the names of any of the Data Structure
// elements in that file change, then these need to change.
namespace
{
const QString k_SyntheticVolumeDataContainer("SyntheticVolumeDataContainer");
const QString k_CellFeatureData("CellFeatureData");
const QString k_CellFeatureData2("CellFeatureData2");
const QString k_Centroids("Centroids");
const QString k_Phases2("Phases2");
} // namespace
class EstablishFoamMorphologyTest
{

public:
  EstablishFoamMorphologyTest() = default;
  virtual ~EstablishFoamMorphologyTest() = default;

  // -----------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------
  void RemoveTestFiles()
  {
#if REMOVE_TEST_FILES
      // QFile::remove(UnitTest::EstablishFoamMorphologyTest::TestFile1);
      // QFile::remove(UnitTest::EstablishFoamMorphologyTest::TestFile2);
#endif
  }

  // -----------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------
  int TestFilterAvailability()
  {
    // Now instantiate the EstablishFoamMorphologyTest Filter from the FilterManager
    QString filtName = "EstablishFoamMorphology";
    FilterManager* fm = FilterManager::Instance();
    IFilterFactory::Pointer filterFactory = fm->getFactoryFromClassName(filtName);
    if(nullptr == filterFactory.get())
    {
      std::stringstream ss;
      ss << "The EstablishFoamMorphologyTest Requires the use of the " << filtName.toStdString() << " filter which is found in the MASSIFUtilities Plugin";
      DREAM3D_TEST_THROW_EXCEPTION(ss.str())
    }
    return 0;
  }

  // -----------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------
  int TestEstablishFoamMorphologyTest()
  {
    QString pipelineFile = UnitTest::PluginSourceDir + "/" + "ExamplePipelines/" + UnitTest::PluginName + "/" + "Open-Cell-Foam-Example.json";

    // Sanity Check the filepath to make sure it exists, Report an error and bail if it does not
    QFileInfo fi(pipelineFile);
    if(!fi.exists())
    {
      std::cout << "The input file '" << pipelineFile.toStdString() << "' does not exist" << std::endl;
      return EXIT_FAILURE;
    }

    // Use the static method to read the Pipeline file and return a Filter Pipeline
    FilterPipeline::Pointer pipeline;
    JsonFilterParametersReader::Pointer jsonReader = JsonFilterParametersReader::New();
    pipeline = jsonReader->readPipelineFromFile(pipelineFile);

    if(nullptr == pipeline.get())
    {
      std::cout << "An error occurred trying to read the pipeline file. Exiting now." << std::endl;
      return EXIT_FAILURE;
    }

    std::cout << "Pipeline Count: " << pipeline->size() << std::endl;
    Observer obs; // Create an Observer to report errors/progress from the executing pipeline
    pipeline->addMessageReceiver(&obs);
    // Preflight the pipeline
    int32_t err = pipeline->preflightPipeline();
    if(err < 0)
    {
      std::cout << "Errors preflighting the pipeline. Exiting Now." << std::endl;
    }
    DREAM3D_REQUIRE(err >= 0)

    // Now actually execute the pipeline
    DataContainerArray::Pointer dca = pipeline->execute();
    err = pipeline->getErrorCode();
    if(err < 0)
    {
      std::cout << "Error Condition of Pipeline: " << err << std::endl;
    }
    DREAM3D_REQUIRE(err >= 0)

    DataContainer::Pointer dc = dca->getDataContainer(::k_SyntheticVolumeDataContainer);

    AttributeMatrix::Pointer featureAttrMat = dc->getAttributeMatrix(::k_CellFeatureData);

    FloatArrayType::Pointer centroids = std::dynamic_pointer_cast<FloatArrayType>(featureAttrMat->getAttributeArray(::k_Centroids));
    // In the example pipeline there should be somewhere around 6264 features created. Should be 61 but with the randomness
    // of the test we can use a tolerance instead of +-10.

    size_t numTuples = centroids->getNumberOfTuples();
    DREAM3D_REQUIRED(numTuples, >, 50)
    DREAM3D_REQUIRED(numTuples, <, 70)

    AttributeMatrix::Pointer featureAttrMat2 = dc->getAttributeMatrix(::k_CellFeatureData2);
    Int32ArrayType::Pointer phases2 = std::dynamic_pointer_cast<Int32ArrayType>(featureAttrMat2->getAttributeArray(::k_Phases2));
    numTuples = phases2->getNumberOfTuples();
    DREAM3D_REQUIRED(numTuples, >, 2600)
    DREAM3D_REQUIRED(numTuples, <, 2900)

    return EXIT_SUCCESS;
  }

  // -----------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------
  void operator()()
  {
    std::cout << "###### EstablishFoamMorphologyTest ######" << std::endl;
    int err = EXIT_SUCCESS;

    DREAM3D_REGISTER_TEST(TestFilterAvailability());

    DREAM3D_REGISTER_TEST(TestEstablishFoamMorphologyTest())

    DREAM3D_REGISTER_TEST(RemoveTestFiles())
  }

public:
  EstablishFoamMorphologyTest(const EstablishFoamMorphologyTest&) = delete;            // Copy Constructor Not Implemented
  EstablishFoamMorphologyTest(EstablishFoamMorphologyTest&&) = delete;                 // Move Constructor Not Implemented
  EstablishFoamMorphologyTest& operator=(const EstablishFoamMorphologyTest&) = delete; // Copy Assignment Not Implemented
  EstablishFoamMorphologyTest& operator=(EstablishFoamMorphologyTest&&) = delete;      // Move Assignment Not Implemented
};
