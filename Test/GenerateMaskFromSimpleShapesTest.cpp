// -----------------------------------------------------------------------------
// Insert your license & copyright information here
// -----------------------------------------------------------------------------
#pragma once

#include <QtCore/QCoreApplication>
#include <QtCore/QFile>

#include "SIMPLib/SIMPLib.h"
#include "SIMPLib/DataArrays/DataArray.hpp"
#include "SIMPLib/DataContainers/AttributeMatrix.h"
#include "SIMPLib/DataContainers/DataContainer.h"
#include "SIMPLib/DataContainers/DataContainerArray.h"
#include "SIMPLib/Filtering/FilterFactory.hpp"
#include "SIMPLib/Filtering/FilterManager.h"
#include "SIMPLib/Filtering/FilterPipeline.h"
#include "SIMPLib/Filtering/QMetaObjectUtilities.h"
#include "SIMPLib/Plugin/ISIMPLibPlugin.h"
#include "SIMPLib/Plugin/SIMPLibPluginLoader.h"

#include "DREAM3DReview/DREAM3DReviewFilters/GenerateMaskFromSimpleShapes.h"

#include "UnitTestSupport.hpp"

#include "DREAM3DReviewTestFileLocations.h"

class GenerateMaskFromSimpleShapesTest
{

public:
  GenerateMaskFromSimpleShapesTest() = default;
  ~GenerateMaskFromSimpleShapesTest() = default;
  GenerateMaskFromSimpleShapesTest(const GenerateMaskFromSimpleShapesTest&) = delete;            // Copy Constructor
  GenerateMaskFromSimpleShapesTest(GenerateMaskFromSimpleShapesTest&&) = delete;                 // Move Constructor
  GenerateMaskFromSimpleShapesTest& operator=(const GenerateMaskFromSimpleShapesTest&) = delete; // Copy Assignment
  GenerateMaskFromSimpleShapesTest& operator=(GenerateMaskFromSimpleShapesTest&&) = delete;      // Move Assignment

  // -----------------------------------------------------------------------------
  DataContainerArray::Pointer createDataStructure()
  {
    DataContainerArray::Pointer dca = DataContainerArray::New();
    DataContainer::Pointer dc = DataContainer::New("Test");
    AttributeMatrix::Pointer am = AttributeMatrix::New({100}, "AM", AttributeMatrix::Type::Cell);

    dca->addOrReplaceDataContainer(dc);
    dc->addOrReplaceAttributeMatrix(am);
    return dca;
  }

  // -----------------------------------------------------------------------------
  int TestGenerateMaskFromSimpleShapesTest()
  {

    DataContainerArray::Pointer dca = createDataStructure();
    GenerateMaskFromSimpleShapes::Pointer filter = GenerateMaskFromSimpleShapes::New();
    filter->setDataContainerArray(dca);
    filter->preflight();
    int32_t err = filter->getErrorCode();
    DREAM3D_REQUIRE_EQUAL(err, -999)

    return EXIT_SUCCESS;
  }

  // -----------------------------------------------------------------------------
  void operator()()
  {
    int err = EXIT_SUCCESS;

    DREAM3D_REGISTER_TEST(TestGenerateMaskFromSimpleShapesTest())
  }
};
