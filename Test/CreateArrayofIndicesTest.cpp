// -----------------------------------------------------------------------------
// Insert your license & copyright information here
// -----------------------------------------------------------------------------
#pragma once

#include <QtCore/QCoreApplication>
#include <QtCore/QFile>

#include "SIMPLib/SIMPLib.h"
#include "SIMPLib/DataArrays/DataArray.hpp"
#include "SIMPLib/Filtering/FilterPipeline.h"
#include "SIMPLib/Filtering/FilterManager.h"
#include "SIMPLib/Filtering/FilterFactory.hpp"
#include "SIMPLib/Plugin/ISIMPLibPlugin.h"
#include "SIMPLib/Plugin/SIMPLibPluginLoader.h"
#include "SIMPLib/Filtering/QMetaObjectUtilities.h"

#include "UnitTestSupport.hpp"

#include "DREAM3DReviewTestFileLocations.h"

class CreateArrayofIndicesTest
{

  public:
    CreateArrayofIndicesTest() = default;
    ~CreateArrayofIndicesTest() = default;
    CreateArrayofIndicesTest(const CreateArrayofIndicesTest&) = delete;            // Copy Constructor
    CreateArrayofIndicesTest(CreateArrayofIndicesTest&&) = delete;                 // Move Constructor
    CreateArrayofIndicesTest& operator=(const CreateArrayofIndicesTest&) = delete; // Copy Assignment
    CreateArrayofIndicesTest& operator=(CreateArrayofIndicesTest&&) = delete;      // Move Assignment

  // -----------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------
  int TestFilterAvailability()
  {
    // Now instantiate the CreateArrayofIndicesTest Filter from the FilterManager
    QString filtName = "CreateArrayofIndices";
    FilterManager* fm = FilterManager::Instance();
    IFilterFactory::Pointer filterFactory = fm->getFactoryFromClassName(filtName);
    if (nullptr == filterFactory.get())
    {
      std::stringstream ss;
      ss << "The CreateArrayofIndicesTest Requires the use of the " << filtName.toStdString()
         << " filter which is found in the DREAM3DReview Plugin";
      DREAM3D_TEST_THROW_EXCEPTION(ss.str())
    }
    return 0;
  }

  // -----------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------
  int TestCreateArrayofIndicesTest()
  {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   /* Please write CreateArrayofIndicesTest test code here.
    *
    * To create IO test files, please edit the file template at DREAM3DReview/Test/TestFileLocations.h.in.
    * Add a CreateArrayofIndicesTest namespace inside the UnitTest namespace, and add your test file paths to your new namespace.
    *
    * SIMPLib provides some macros that will throw exceptions when a test fails
    * and thus report that during testing. These macros are located in the
    * SIMPLib/Utilities/UnitTestSupport.hpp file. Some examples are:
    *
    * SIMPLib_REQUIRE_EQUAL(foo, 0)
    * This means that if the variable foo is NOT equal to Zero then test will fail
    * and the current test will exit immediately. If there are more tests registered
    * with the SIMPLib_REGISTER_TEST() macro, the next test will execute. There are
    * lots of examples in the SIMPLib/Test folder to look at.
    */
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    int foo = -1;
    DREAM3D_REQUIRE_EQUAL(foo, 0)

    return EXIT_SUCCESS;
  }

  // -----------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------
  void operator()()
  {
    int err = EXIT_SUCCESS;

    DREAM3D_REGISTER_TEST( TestFilterAvailability() );

    DREAM3D_REGISTER_TEST( TestCreateArrayofIndicesTest() )
  }

  private:


};

