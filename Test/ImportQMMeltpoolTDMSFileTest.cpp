// -----------------------------------------------------------------------------
// Insert your license & copyright information here
// -----------------------------------------------------------------------------
#pragma once

#include <QtCore/QCoreApplication>
#include <QtCore/QFile>

#include "SIMPLib/SIMPLib.h"
#include "SIMPLib/Common/SIMPLibSetGetMacros.h"
#include "SIMPLib/DataArrays/DataArray.hpp"
#include "SIMPLib/Filtering/FilterFactory.hpp"
#include "SIMPLib/Filtering/FilterManager.h"
#include "SIMPLib/Filtering/FilterPipeline.h"
#include "SIMPLib/Filtering/QMetaObjectUtilities.h"
#include "SIMPLib/Plugin/ISIMPLibPlugin.h"
#include "SIMPLib/Plugin/SIMPLibPluginLoader.h"

#include "UnitTestSupport.hpp"

#include "DREAM3DReview/DREAM3DReviewFilters/ImportQMMeltpoolTDMSFile.h"
#include "DREAM3DReviewTestFileLocations.h"

class ImportQMMeltpoolTDMSFileTest
{

public:
  ImportQMMeltpoolTDMSFileTest() = default;
  ~ImportQMMeltpoolTDMSFileTest() = default;
  ImportQMMeltpoolTDMSFileTest(const ImportQMMeltpoolTDMSFileTest&) = delete;            // Copy Constructor
  ImportQMMeltpoolTDMSFileTest(ImportQMMeltpoolTDMSFileTest&&) = delete;                 // Move Constructor
  ImportQMMeltpoolTDMSFileTest& operator=(const ImportQMMeltpoolTDMSFileTest&) = delete; // Copy Assignment
  ImportQMMeltpoolTDMSFileTest& operator=(ImportQMMeltpoolTDMSFileTest&&) = delete;      // Move Assignment

  const QString k_TestFile1 = UnitTest::TestTempDir + "/TestFile1.txt";
  const QString k_TestFile2 = UnitTest::TestTempDir + "/TestFile2.txt";

  // -----------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------
  void RemoveTestFiles()
  {
#if REMOVE_TEST_FILES
    QFile::remove(k_TestFile1);
    QFile::remove(k_TestFile2);
#endif
  }

  // -----------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------
  int TestFilterAvailability()
  {
    // Now instantiate the ImportQMMeltpoolTDMSFile Filter from the FilterManager
    ImportQMMeltpoolTDMSFile::Pointer filter = ImportQMMeltpoolTDMSFile::New();
    filter->preflight();
    int32_t err = filter->getErrorCode();
    DREAM3D_REQUIRE_EQUAL(err, -387)
    return 0;
  }

  // -----------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------
  int TestImportQMMeltpoolTDMSFileTest()
  {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    /* Please write ImportQMMeltpoolTDMSFileTest test code here.
     *
     * Your IO test files are:
     * UnitTest::ImportQMMeltpoolTDMSFileTest::TestFile1
     * UnitTest::ImportQMMeltpoolTDMSFileTest::TestFile2
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

    return EXIT_SUCCESS;
  }

  // -----------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------
  void operator()()
  {
    int err = EXIT_SUCCESS;

    DREAM3D_REGISTER_TEST(TestFilterAvailability());

    DREAM3D_REGISTER_TEST(TestImportQMMeltpoolTDMSFileTest())

    DREAM3D_REGISTER_TEST(RemoveTestFiles())
  }

private:
};
