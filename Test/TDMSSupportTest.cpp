

#include <QtCore/QString>

#include "DREAM3DReview/TDMSSupport/TDMSFileProxy.h"
#include "DREAM3DReview/TDMSSupport/TDMSObject.h"

#include "DREAM3DReviewTestFileLocations.h"
#include "UnitTestSupport.hpp"

class TDMSSupportTest
{
public:
  TDMSSupportTest() = default;
  ~TDMSSupportTest() = default;
  TDMSSupportTest(const TDMSSupportTest&) = delete;            // Copy Constructor Not Implemented
  TDMSSupportTest(TDMSSupportTest&&) = delete;                 // Move Constructor Not Implemented
  TDMSSupportTest& operator=(const TDMSSupportTest&) = delete; // Copy Assignment Not Implemented
  TDMSSupportTest& operator=(TDMSSupportTest&&) = delete;      // Move Assignment Not Implemented

  const QString k_TestFile = {""};

  // -----------------------------------------------------------------------------
  void RemoveTestFiles()
  {
#if REMOVE_TEST_FILES
    //    QFile::remove(k_TestFile);
#endif
  }

  // -----------------------------------------------------------------------------
  int ReadTest()
  {

    TDMSFileProxy::Pointer fileProxy = TDMSFileProxy::New("/Volumes/RAID-0/LockheedMartin/TDMS_200120_12-40_2020-01-20 ATRQ Build 2_Slice_00001_to_00040/Slice00001.tdms");
    fileProxy->readMetaData();
    using TDMSObjectPointerType = TDMSObject::Pointer;
    using TDMSObjectMapType = std::unordered_map<std::string, TDMSObjectPointerType>;

    TDMSObjectMapType objects = fileProxy->objects();
    for(const auto& object : objects)
    {
      std::cout << "Object Nmae: " << object.first << std::endl;
    }

    return EXIT_SUCCESS;
  }

  // -----------------------------------------------------------------------------
  void operator()()
  {
    int err = EXIT_SUCCESS;
    std::cout << "================ TDMSSupportTest =====================" << std::endl;
    DREAM3D_REGISTER_TEST(ReadTest());

    DREAM3D_REGISTER_TEST(RemoveTestFiles())
  }
};
