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
 * The code contained herein was partially funded by the following contracts:
 *    United States Air Force Prime Contract FA8650-15-D-5231
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

/*
TEST ONLY CHECKS FEATURES ADDED FOR IMAGE TRANSFORMATION BASED ON RotateSampleRefFrameTest
*/
#include <QtCore/QFile>

#include <Eigen/Dense>

#include "SIMPLib/SIMPLib.h"
#include "SIMPLib/CoreFilters/DataContainerReader.h"
#include "SIMPLib/DataArrays/DataArray.hpp"
#include "SIMPLib/DataContainers/DataContainerArray.h"
#include "SIMPLib/FilterParameters/DynamicTableData.h"
#include "SIMPLib/FilterParameters/FloatVec3FilterParameter.h"
#include "SIMPLib/Filtering/FilterFactory.hpp"
#include "SIMPLib/Filtering/FilterManager.h"
#include "SIMPLib/Filtering/FilterPipeline.h"
#include "SIMPLib/Filtering/QMetaObjectUtilities.h"
#include "SIMPLib/Geometry/ImageGeom.h"
#include "SIMPLib/Math/SIMPLibMath.h"
#include "SIMPLib/Plugin/ISIMPLibPlugin.h"
#include "SIMPLib/Plugin/SIMPLibPluginLoader.h"
#include "UnitTestSupport.hpp"

#include "DREAM3DReviewTestFileLocations.h"

class ApplyTransformationToGeometryTest
{
private:
  const QString k_FilterName = "ApplyTransformationToGeometry";
  const QString k_TransformationMatrixType = "TransformationMatrixType";
  const QString k_ManualTransformationMatrix = "ManualTransformationMatrix";
  const QString k_InterpolationType = "InterpolationType";
  const QString k_RotationAngle = "RotationAngle";
  const QString k_RotationAxis = "RotationAxis";
  const QString k_CellAttributeMatrixPath = "CellAttributeMatrixPath";

  template <class T>
  static bool dataArrayEqual(const DataArray<T>& a, const DataArray<T>& b)
  {
    if(a.getNumberOfComponents() != b.getNumberOfComponents())
    {
      return false;
    }

    if(a.getNumberOfTuples() != b.getNumberOfTuples())
    {
      return false;
    }

    return std::equal(a.begin(), a.end(), b.begin(), b.end());
  }

  static bool imageGeomEqual(const ImageGeom& a, const ImageGeom& b)
  {
    float min = .99;
    float max = 1.01;
    SizeVec3Type dimsA = a.getDimensions();
    FloatVec3Type originA = a.getOrigin();
    FloatVec3Type spacingA = a.getSpacing();

    SizeVec3Type dimsB = b.getDimensions();
    FloatVec3Type originB = b.getOrigin();
    FloatVec3Type spacingB = b.getSpacing();

    SizeVec3Type dimsBMax;
    SizeVec3Type dimsBMin;
    FloatVec3Type originBMax;
    FloatVec3Type originBMin;
    FloatVec3Type spacingBMax;
    FloatVec3Type spacingBMin;

    for(int i = 0; i < 3; i++)
    {
      spacingBMax[i] = spacingB[i] * max;
      spacingBMin[i] = spacingB[i] * min;
    }

    bool spacingValid = true;

    for(int i = 0; i < 3; i++)
    {
      if(!(spacingBMax[i] >= spacingA[i]))
      {
        spacingValid = false;
	  }
	  if(!(spacingBMin[i] <= spacingA[i]))
	  {
        spacingValid = false;
	  }
    }

	return (dimsA == dimsB) && (originA == originB) && spacingValid;
  }

  static bool matrixEqual(const AttributeMatrix& a, const AttributeMatrix& b)
  {
    return a.getTupleDimensions() == b.getTupleDimensions();
  }

  static void dataContainerEqual(const DataContainerArray& dca, const DataArrayPath& testPath, const DataArrayPath& expectedPath)
  {
    // Get test

    DataContainer::Pointer testContainer = dca.getDataContainer(testPath);
    DREAM3D_REQUIRE_VALID_POINTER(testContainer)

    ImageGeom::Pointer testGeom = testContainer->getGeometryAs<ImageGeom>();
    DREAM3D_REQUIRE_VALID_POINTER(testGeom)

    AttributeMatrix::Pointer testMatrix = testContainer->getAttributeMatrix(testPath);
    DREAM3D_REQUIRE_VALID_POINTER(testMatrix)

    UInt8ArrayType::Pointer testArray = testMatrix->getAttributeArrayAs<UInt8ArrayType>(testPath.getDataArrayName());
    DREAM3D_REQUIRE_VALID_POINTER(testArray)

    // Get expected

    DataContainer::Pointer expectedContainer = dca.getDataContainer(expectedPath);
    DREAM3D_REQUIRE_VALID_POINTER(expectedContainer)

    ImageGeom::Pointer expectedGeom = expectedContainer->getGeometryAs<ImageGeom>();
    DREAM3D_REQUIRE_VALID_POINTER(expectedGeom)

    AttributeMatrix::Pointer expectedMatrix = expectedContainer->getAttributeMatrix(expectedPath);
    DREAM3D_REQUIRE_VALID_POINTER(expectedMatrix)

    UInt8ArrayType::Pointer expectedArray = testMatrix->getAttributeArrayAs<UInt8ArrayType>(expectedPath.getDataArrayName());
    DREAM3D_REQUIRE_VALID_POINTER(expectedArray)

    // Compare

    bool geometryEqual = imageGeomEqual(*testGeom, *expectedGeom);
    DREAM3D_REQUIRE(geometryEqual)

    bool attributeMatrixEqual = matrixEqual(*testMatrix, *expectedMatrix);
    DREAM3D_REQUIRE(attributeMatrixEqual)

    bool arraysEqual = dataArrayEqual(*testArray, *expectedArray);
    DREAM3D_REQUIRE(arraysEqual)
  }

  static void resetGeometry(AttributeMatrix::Pointer matrix, ImageGeom::Pointer imageGeom, const std::vector<size_t>& tDims)
  {
    matrix->resizeAttributeArrays(tDims);
    imageGeom->setDimensions(tDims);
    imageGeom->setOrigin(0.0f, 0.0f, 0.0f);
    imageGeom->setSpacing(1.0f, 1.0f, 1.0f);
  }

  AbstractFilter::Pointer createFilter() const
  {
    FilterManager* fm = FilterManager::Instance();
    IFilterFactory::Pointer filterFactory = fm->getFactoryFromClassName(k_FilterName);
    DREAM3D_REQUIRE_VALID_POINTER(filterFactory)

    AbstractFilter::Pointer transformFilter = filterFactory->create();
    DREAM3D_REQUIRE_VALID_POINTER(transformFilter)

    return transformFilter;
  }

  template <class T>
  static void setProperty(AbstractFilter::Pointer filter, const QString& property, const T& value)
  {
    QVariant variant;
    variant.setValue(value);
    bool propWasSet = filter->setProperty(property.toStdString().c_str(), variant);
    DREAM3D_REQUIRE_EQUAL(propWasSet, true)
  }

public:
  ApplyTransformationToGeometryTest() = default;
  ~ApplyTransformationToGeometryTest() = default;
  ApplyTransformationToGeometryTest(const ApplyTransformationToGeometryTest&) = delete;            // Copy Constructor Not Implemented
  ApplyTransformationToGeometryTest(ApplyTransformationToGeometryTest&&) = delete;                 // Move Constructor Not Implemented
  ApplyTransformationToGeometryTest& operator=(const ApplyTransformationToGeometryTest&) = delete; // Copy Assignment Not Implemented
  ApplyTransformationToGeometryTest& operator=(ApplyTransformationToGeometryTest&&) = delete;      // Move Assignment Not Implemented

  // -----------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------
  void RemoveTestFiles()
  {
#if REMOVE_TEST_FILES
#endif
  }

  // -----------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------
  int TestFilterAvailability()
  {
    // Now instantiate the SampleSurfaceMeshSpecifiedPointsTest Filter from the FilterManager
    FilterManager* filterManager = FilterManager::Instance();
    IFilterFactory::Pointer filterFactory = filterManager->getFactoryFromClassName(k_FilterName);
    if(filterFactory == nullptr)
    {
      std::stringstream ss;
      ss << "The ApplyTransformationToGeometryTest Requires the use of the " << k_FilterName.toStdString() << " filter which is found in the Sampling Plugin";
      DREAM3D_TEST_THROW_EXCEPTION(ss.str())
    }
    return 0;
  }

  // -----------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------
  void TestFilterParameters()
  {
    const std::vector<size_t> tDims{10, 20, 30};
    const DataArrayPath path("DataContainer", "AttributeMatrix", "DataArray");

    DataContainerArray::Pointer dca = DataContainerArray::New();

    DataContainer::Pointer dc = DataContainer::New(path.getDataContainerName());

    ImageGeom::Pointer imageGeom = ImageGeom::New();

    imageGeom->setDimensions(tDims);

    dc->setGeometry(imageGeom);

    dca->addOrReplaceDataContainer(dc);

    AttributeMatrix::Pointer matrix = dc->createNonPrereqAttributeMatrix(nullptr, path.getAttributeMatrixName(), tDims, AttributeMatrix::Type::Cell);

    AbstractFilter::Pointer transformFilter = createFilter();
    transformFilter->setDataContainerArray(dca);

    setProperty(transformFilter, k_CellAttributeMatrixPath, matrix->getDataArrayPath());
    setProperty(transformFilter, k_TransformationMatrixType, 2);

    std::vector<std::vector<double>> dataMatrix = {{0, 0, 1, 0}, {0, 1, 0, 0}, {-1, 0, 0, 0}, {0, 0, 0, 1}};
    QString name = "_INTERNAL_USE_ONLY_TestFilterParameters_manMatrix";
    DynamicTableData manMatrix(dataMatrix);
    setProperty(transformFilter, k_ManualTransformationMatrix, manMatrix);
    setProperty(transformFilter, k_InterpolationType, 0);

    transformFilter->preflight();
    int error = transformFilter->getErrorCode();
    DREAM3D_REQUIRED(error, >=, 0)

    resetGeometry(matrix, imageGeom, tDims);

    // No op error
    setProperty(transformFilter, k_TransformationMatrixType, 0);

    transformFilter->preflight();
    error = transformFilter->getErrorCode();
    DREAM3D_REQUIRED(error, <, 0)

    resetGeometry(matrix, imageGeom, tDims);

    // Manual Transformation row number error
    setProperty(transformFilter, k_TransformationMatrixType, 2);
    std::vector<std::vector<double>> dataMatrixError = {{0, 0, 1, 0}, {0, 1, 0, 0}, {-1, 0, 0, 0}};
    DynamicTableData manMatrixError(dataMatrixError);
    setProperty(transformFilter, k_ManualTransformationMatrix, manMatrixError);
    transformFilter->preflight();
    error = transformFilter->getErrorCode();
    DREAM3D_REQUIRED(error, <, 0)

    resetGeometry(matrix, imageGeom, tDims);

    // Manual Transformation column number error
    std::vector<std::vector<double>> dataMatrixError2 = {{0, 0, 1}, {0, 1, 0}, {-1, 0, 0}, {0, 0, 0}};
    DynamicTableData manMatrixError2(dataMatrixError2);
    setProperty(transformFilter, k_ManualTransformationMatrix, manMatrixError2);
    transformFilter->preflight();
    error = transformFilter->getErrorCode();
    DREAM3D_REQUIRED(error, <, 0)

    resetGeometry(matrix, imageGeom, tDims);

    // Invalid Transformation Type
    setProperty(transformFilter, k_TransformationMatrixType, 9);
    setProperty(transformFilter, k_ManualTransformationMatrix, manMatrix);
    transformFilter->preflight();
    error = transformFilter->getErrorCode();
    DREAM3D_REQUIRED(error, <, 0)
    resetGeometry(matrix, imageGeom, tDims);
  }

  // -----------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------
  void TestApplyTransformationToGeometry()
  {
    const std::vector<size_t> tDims{10, 20, 30};
    const DataArrayPath basePath("", "CellData", "Data");

    // Read test file

    DataContainerArray::Pointer dca = DataContainerArray::New();

    DataContainerReader::Pointer dataContainerReader = DataContainerReader::New();
    dataContainerReader->setDataContainerArray(dca);
    dataContainerReader->setInputFile(UnitTest::ApplyTransformationToGeometryTest::TestFile1);
    DataContainerArrayProxy proxy = dataContainerReader->readDataContainerArrayStructure(UnitTest::ApplyTransformationToGeometryTest::TestFile1);
    dataContainerReader->setInputFileDataContainerArrayProxy(proxy);

    dataContainerReader->execute();
    int error = dataContainerReader->getErrorCode();
    DREAM3D_REQUIRED(error, >=, 0)

    QList<QString> dcNames = dca->getDataContainerNames();

    bool hasOriginal = dcNames.contains("Original");

    DREAM3D_REQUIRE(hasOriginal)

    AbstractFilter::Pointer transformFilter = createFilter();
    transformFilter->setDataContainerArray(dca);

    DataArrayPath originalPath = basePath;
    originalPath.setDataContainerName("Original");
    DataContainer::Pointer dcOriginal = dca->getDataContainer(originalPath);
    DREAM3D_REQUIRE_VALID_POINTER(dcOriginal)

    bool transformed = false;

    for(const auto& name : dcNames)
    {
      bool isTransformed = name.startsWith("Apply");
      transformed = isTransformed ? true : transformed;
      if(!isTransformed)
      {
        continue;
      }

      DataArrayPath expectedPath = basePath;
      expectedPath.setDataContainerName(name);

      QStringList parts = name.split("_");
      int size = parts.size();
      DREAM3D_REQUIRED(size, ==, 5)

      bool conversionSuccess = false;

      float x = parts[1].toFloat(&conversionSuccess);
      DREAM3D_REQUIRE(conversionSuccess)
      float y = parts[2].toFloat(&conversionSuccess);
      DREAM3D_REQUIRE(conversionSuccess)
      float z = parts[3].toFloat(&conversionSuccess);
      DREAM3D_REQUIRE(conversionSuccess)
      float angle = parts[4].toFloat(&conversionSuccess);
      DREAM3D_REQUIRE(conversionSuccess)

      {
        setProperty(transformFilter, k_TransformationMatrixType, 3);

        // Copy original data for rotation

        DataContainer::Pointer dcCopy = dcOriginal->deepCopy();
        dcCopy->setName(name + "_Copy");
        dca->addOrReplaceDataContainer(dcCopy);

        DataArrayPath testPath = basePath;
        testPath.setDataContainerName(dcCopy->getName());

        DataContainer::Pointer dcTest = dca->getDataContainer(testPath);
        DREAM3D_REQUIRE_VALID_POINTER(dcTest)

        setProperty(transformFilter, k_CellAttributeMatrixPath, testPath);

        // Set properties
        setProperty(transformFilter, k_RotationAngle, angle);
        setProperty(transformFilter, k_RotationAxis, FloatVec3Type(x, y, z));
        setProperty(transformFilter, k_InterpolationType, 0);

        // Execute
        transformFilter->execute();
        int error = transformFilter->getErrorCode();
        DREAM3D_REQUIRED(error, >=, 0)

        // Compare
        std::string arrname = testPath.getDataArrayName().toStdString();
        QJsonObject json;
        testPath.writeJson(json);
        dataContainerEqual(*dca, testPath, expectedPath);
      }

      {
        //// Set to rotation matrix representation
        setProperty(transformFilter, k_TransformationMatrixType, 2);

        // Copy original data for rotation

        DataContainer::Pointer dcCopy = dcOriginal->deepCopy();
        dcCopy->setName(name + "_Copy_Matrix");
        dca->addOrReplaceDataContainer(dcCopy);

        DataArrayPath testPath = basePath;
        testPath.setDataContainerName(dcCopy->getName());

        setProperty(transformFilter, k_CellAttributeMatrixPath, testPath);
        setProperty(transformFilter, k_InterpolationType, 0);

        // Set properties
        Eigen::Vector3f axis(x, y, z);
        float angleRadians = angle * SIMPLib::Constants::k_DegToRadF;
        Eigen::AngleAxisf axisAngle(angleRadians, axis);

        Eigen::Matrix3f rotationMatrix = axisAngle.toRotationMatrix();

        std::vector<std::vector<double>> dataMatrix;

        for(int i = 0; i <= rotationMatrix.rows(); i++)
        {
          std::vector<double> row;
          for(int j = 0; j <= rotationMatrix.cols(); j++)
          {
            if(j == rotationMatrix.cols())
            {
              double zero = 0.0;
              row.push_back(zero);
            }
            else
            {
              if(i < rotationMatrix.rows())
              {
                row.push_back(static_cast<double>(rotationMatrix(i, j)));
              }
            }
          }
          if(i == rotationMatrix.rows())
          {
            std::vector<double> vect = {0.0, 0.0, 0.0, 1.0};
            dataMatrix.push_back(vect);
          }
          else
          {
            dataMatrix.push_back(row);
          }
        }

        DynamicTableData tableData(dataMatrix);

        setProperty(transformFilter, k_ManualTransformationMatrix, tableData);

        // Execute
        transformFilter->execute();
        int error = transformFilter->getErrorCode();
        DREAM3D_REQUIRED(error, >=, 0)

        // Compare

        dataContainerEqual(*dca, testPath, expectedPath);
      }
    }

    DREAM3D_REQUIRE(transformed)
  }

  // -----------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------
  void operator()()
  {
    std::cout << "----Start ApplyTransformationToGeometryTest----\n";

    int err = EXIT_SUCCESS;
    DREAM3D_REGISTER_TEST(TestFilterAvailability())

    DREAM3D_REGISTER_TEST(TestFilterParameters())
    DREAM3D_REGISTER_TEST(TestApplyTransformationToGeometry())

    DREAM3D_REGISTER_TEST(RemoveTestFiles())

    std::cout << "----End ApplyTransformationToGeometryTest----\n";
  }
};
