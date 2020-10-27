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

#include "UnitTestSupport.hpp"

#include "DREAM3DReviewTestFileLocations.h"

#include "DREAM3DReview/DREAM3DReviewFilters/CreateArrayofIndices.h"

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
  int TestCreateArrayofIndicesTest()
  {
    DataContainerArray::Pointer dca = createDataStructure();
    DataArrayPath indexArrayPath("Test", "AM", "Indices");

    CreateArrayofIndices::Pointer filter = CreateArrayofIndices::New();
    filter->setDataContainerArray(dca);
    filter->setIndexArrayPath(indexArrayPath);
    filter->preflight();
    int32_t err = filter->getErrorCode();
    DREAM3D_REQUIRE_EQUAL(err, 0)

    dca = createDataStructure();
    filter->setDataContainerArray(dca);
    filter->setIndexArrayPath(indexArrayPath);
    filter->execute();
    err = filter->getErrorCode();
    DREAM3D_REQUIRE_EQUAL(err, 0)

    SizeTArrayType& indices = *(dca->getAttributeMatrix(indexArrayPath)->getAttributeArrayAs<SizeTArrayType>("Indices"));
    for(size_t i = 0; i < indices.size(); i++)
    {
      DREAM3D_REQUIRE_EQUAL(i, indices[i]);
    }

    return EXIT_SUCCESS;
  }

  // -----------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------
  void operator()()
  {
    int err = EXIT_SUCCESS;

    DREAM3D_REGISTER_TEST(TestCreateArrayofIndicesTest())
  }

private:
};
