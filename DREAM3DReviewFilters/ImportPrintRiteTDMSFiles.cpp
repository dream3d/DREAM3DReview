#include "ImportPrintRiteTDMSFiles.h"

#include <chrono>
#include <random>
#include <unordered_set>
#include <utility>

#include <QtCore/QDir>
#include <QtCore/QFileInfo>
#include <QtCore/QJsonArray>
#include <QtCore/QJsonDocument>

#include "H5Support/QH5Utilities.h"

using namespace H5Support;
#include "SIMPLib/Common/Constants.h"
#include "SIMPLib/Common/TemplateHelpers.h"
#include "SIMPLib/DataContainers/DataContainer.h"
#include "SIMPLib/DataContainers/DataContainerArray.h"
#include "SIMPLib/FilterParameters/BooleanFilterParameter.h"
#include "SIMPLib/FilterParameters/ChoiceFilterParameter.h"
#include "SIMPLib/FilterParameters/FloatFilterParameter.h"
#include "SIMPLib/FilterParameters/FloatVec3FilterParameter.h"
#include "SIMPLib/FilterParameters/InputFileFilterParameter.h"
#include "SIMPLib/FilterParameters/IntFilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedBooleanFilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedChoicesFilterParameter.h"
#include "SIMPLib/FilterParameters/OutputPathFilterParameter.h"
#include "SIMPLib/FilterParameters/SeparatorFilterParameter.h"
#include "SIMPLib/FilterParameters/StringFilterParameter.h"
#include "SIMPLib/Filtering/FilterFactory.hpp"
#include "SIMPLib/Filtering/FilterManager.h"
#include "SIMPLib/Geometry/EdgeGeom.h"
#include "SIMPLib/Geometry/TriangleGeom.h"
#include "SIMPLib/Utilities/FilePathGenerator.h"

#include "DREAM3DReview/DREAM3DReviewConstants.h"
#include "DREAM3DReview/DREAM3DReviewFilters/util/Delaunay2D.h"
#include "DREAM3DReview/DREAM3DReviewVersion.h"
#include "DREAM3DReview/TDMSSupport/TDMSExceptionHandler.h"
#include "DREAM3DReview/TDMSSupport/TDMSFileProxy.h"

namespace PRH = PrintRiteHelpers;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
ImportPrintRiteTDMSFiles::ImportPrintRiteTDMSFiles()
{
  m_InputFilesList.FileExtension = QString("tdms");
  m_InputFilesList.StartIndex = 0;
  m_InputFilesList.EndIndex = 0;
  m_InputFilesList.PaddingDigits = 0;
  m_InputFilesList.IncrementIndex = 1;

  initialize();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
ImportPrintRiteTDMSFiles::~ImportPrintRiteTDMSFiles() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ImportPrintRiteTDMSFiles::initialize()
{
  clearErrorCode();
  setCancel(false);
  m_PowerScalingCoefficients[0] = 1.0f;
  m_PowerScalingCoefficients[1] = 0.0f;
  m_TemperatureScalingCoefficients[0] = 1.0f;
  m_TemperatureScalingCoefficients[1] = 0.0f;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ImportPrintRiteTDMSFiles::setupFilterParameters()
{
  FilterParameterVectorType parameters;
  parameters.push_back(SeparatorFilterParameter::Create("Build Parameters", FilterParameter::Category::Parameter));
  parameters.push_back(SIMPL_NEW_FLOAT_FP("Layer Thickness", LayerThickness, FilterParameter::Category::Parameter, ImportPrintRiteTDMSFiles));
  {
    ChoiceFilterParameter::Pointer parameter = ChoiceFilterParameter::New();
    parameter->setHumanLabel("Array Used to Determine if Laser is On");
    parameter->setPropertyName("LaserOnArrayOption");
    parameter->setSetterCallback(SIMPL_BIND_SETTER(ImportPrintRiteTDMSFiles, this, LaserOnArrayOption));
    parameter->setGetterCallback(SIMPL_BIND_GETTER(ImportPrintRiteTDMSFiles, this, LaserOnArrayOption));
    std::vector<QString> choices;
    choices.push_back("Laser Drive");
    choices.push_back("On-Axis Photodiode");
    parameter->setChoices(choices);
    parameter->setEditable(false);
    parameter->setCategory(FilterParameter::Category::Parameter);
    parameters.push_back(parameter);
  }
  parameters.push_back(SIMPL_NEW_FLOAT_FP("Laser On Threshold", LaserOnThreshold, FilterParameter::Category::Parameter, ImportPrintRiteTDMSFiles));
  parameters.push_back(SIMPL_NEW_BOOL_FP("Downcast Raw Data", DowncastRawData, FilterParameter::Category::Parameter, ImportPrintRiteTDMSFiles));
  std::vector<QString> linkedProps = {"PowerScalingCoefficients"};
  parameters.push_back(SIMPL_NEW_LINKED_BOOL_FP("Scale Laser Power", ScaleLaserPower, FilterParameter::Category::Parameter, ImportPrintRiteTDMSFiles, linkedProps));
  parameters.push_back(SIMPL_NEW_FLOAT_VEC2_FP("Power Scaling Coefficients", PowerScalingCoefficients, FilterParameter::Category::Parameter, ImportPrintRiteTDMSFiles));

  linkedProps.clear();
  linkedProps.push_back("TemperatureScalingCoefficients");
  parameters.push_back(SIMPL_NEW_LINKED_BOOL_FP("Scale Pyrometer Tempeature", ScalePyrometerTemperature, FilterParameter::Category::Parameter, ImportPrintRiteTDMSFiles, linkedProps));
  parameters.push_back(SIMPL_NEW_FLOAT_VEC2_FP("Temperature Scaling Coefficients", TemperatureScalingCoefficients, FilterParameter::Category::Parameter, ImportPrintRiteTDMSFiles));
  parameters.push_back(SeparatorFilterParameter::Create("Spatial Parameters", FilterParameter::Category::Parameter));
  {
    LinkedChoicesFilterParameter::Pointer parameter = LinkedChoicesFilterParameter::New();
    parameter->setHumanLabel("Spatial Transform Options");
    parameter->setPropertyName("SpatialTransformOption");
    parameter->setSetterCallback(SIMPL_BIND_SETTER(ImportPrintRiteTDMSFiles, this, SpatialTransformOption));
    parameter->setGetterCallback(SIMPL_BIND_GETTER(ImportPrintRiteTDMSFiles, this, SpatialTransformOption));
    std::vector<QString> choices;
    choices.push_back("Do Not Compute Real Space Transformation");
    choices.push_back("Compute Real Space Transformation");
    choices.push_back("Import Spatial Transformation Coefficients From File");
    parameter->setChoices(choices);
    std::vector<QString> linkedProps2 = {"LayerForScaling", "SearchRadius", "SplitRegions1", "SplitRegions2", "STLFilePath1", "STLFilePath2", "InputSpatialTransformFilePath"};
    parameter->setLinkedProperties(linkedProps2);
    parameter->setEditable(false);
    parameter->setCategory(FilterParameter::Category::Parameter);
    parameters.push_back(parameter);
  }
  parameters.push_back(SIMPL_NEW_INTEGER_FP("Layer Index For Computing Real Space Transformation", LayerForScaling, FilterParameter::Category::Parameter, ImportPrintRiteTDMSFiles, {1}));
  parameters.push_back(SIMPL_NEW_FLOAT_FP("Search Radius for Finding Regions", SearchRadius, FilterParameter::Category::Parameter, ImportPrintRiteTDMSFiles, {1}));
  parameters.push_back(SIMPL_NEW_BOOL_FP("Split Contiguous Regions into Separate Files", SplitRegions1, FilterParameter::Category::Parameter, ImportPrintRiteTDMSFiles, {1}));
  parameters.push_back(SIMPL_NEW_BOOL_FP("Split Contiguous Regions into Separate Files", SplitRegions2, FilterParameter::Category::Parameter, ImportPrintRiteTDMSFiles, {2}));
  parameters.push_back(SeparatorFilterParameter::Create("Input File Parameters", FilterParameter::Category::Parameter));
  parameters.push_back(SIMPL_NEW_INPUT_FILE_FP("Build STL File", STLFilePath1, FilterParameter::Category::Parameter, ImportPrintRiteTDMSFiles, "*.stl", "Build STL File", {1}));
  parameters.push_back(SIMPL_NEW_INPUT_FILE_FP("Build STL File", STLFilePath2, FilterParameter::Category::Parameter, ImportPrintRiteTDMSFiles, "*.stl", "Build STL File", {2}));
  parameters.push_back(SIMPL_NEW_INPUT_FILE_FP("Spatial Transformation Coefficients File", InputSpatialTransformFilePath, FilterParameter::Category::Parameter, ImportPrintRiteTDMSFiles, "*.json",
                                               "Transform Coefficients", {2}));
  parameters.push_back(SIMPL_NEW_FILELISTINFO_FP("Input PrintRite Files", InputFilesList, FilterParameter::Category::Parameter, ImportPrintRiteTDMSFiles));
  parameters.push_back(SeparatorFilterParameter::Create("Output File Parameters", FilterParameter::Category::Parameter));
  parameters.push_back(SIMPL_NEW_OUTPUT_PATH_FP("Output File Directory", OutputDirectory, FilterParameter::Category::Parameter, ImportPrintRiteTDMSFiles));
  parameters.push_back(SIMPL_NEW_STRING_FP("Output File Prefix", OutputFilePrefix, FilterParameter::Category::Parameter, ImportPrintRiteTDMSFiles));
  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ImportPrintRiteTDMSFiles::dataCheck()
{
  clearErrorCode();
  clearWarningCode();

  if(getLayerThickness() <= 0.0f)
  {
    QString ss = QObject::tr("The layer thickness must be greater than zero");
    setErrorCondition(-389, ss);
  }

  if(getOutputDirectory().isEmpty())
  {
    QString ss = QObject::tr("The output directory must be set");
    setErrorCondition(-391, ss);
  }

  QDir dir(getOutputDirectory());

  if(getOutputDirectory().isEmpty())
  {
    QString ss = QObject::tr("The output directory must be set");
    setErrorCondition(-391, ss);
  }
  else if(!dir.exists())
  {
    QString ss = QObject::tr("The output directory does not exist; DREAM.3D will attempt to create this path during execution");
    setWarningCondition(-1, ss);
  }

  if(getOutputFilePrefix().isEmpty())
  {
    QString ss = QObject::tr("The output file prefix must be set");
    setErrorCondition(-392, ss);
  }

  if(getInputFilesList().InputPath.isEmpty())
  {
    QString ss = QObject::tr("The input directory must be set");
    setErrorCondition(-13, ss);
  }

  bool hasMissingFiles = false;
  bool orderAscending = false;

  if(getInputFilesList().Ordering == 0)
  {
    orderAscending = true;
  }
  else if(getInputFilesList().Ordering == 1)
  {
    orderAscending = false;
  }

  QVector<QString> fileList = FilePathGenerator::GenerateFileList(getInputFilesList().StartIndex, getInputFilesList().EndIndex, getInputFilesList().IncrementIndex, hasMissingFiles, orderAscending,
                                                                  getInputFilesList().InputPath, getInputFilesList().FilePrefix, getInputFilesList().FileSuffix, getInputFilesList().FileExtension,
                                                                  getInputFilesList().PaddingDigits);

  if(fileList.empty())
  {
    QString ss;
    QTextStream out(&ss);
    out << " No files have been selected for import. Have you set the input directory and other values so that input files will be generated?\n";
    out << "InputPath: " << getInputFilesList().InputPath << "\n";
    out << "FilePrefix: " << getInputFilesList().FilePrefix << "\n";
    out << "FileSuffix: " << getInputFilesList().FileSuffix << "\n";
    out << "FileExtension: " << getInputFilesList().FileExtension << "\n";
    out << "PaddingDigits: " << getInputFilesList().PaddingDigits << "\n";
    out << "StartIndex: " << getInputFilesList().StartIndex << "\n";
    out << "EndIndex: " << getInputFilesList().EndIndex << "\n";

    setErrorCondition(-11, ss);
  }

  QList<QString> missing_fnames;
  std::unordered_set<uint32_t> fileIndices;

  for(auto&& fname : fileList)
  {
    QFileInfo fi(fname);
    if(!fi.exists())
    {
      missing_fnames.push_back(fname);
    }
    else
    {
      QString baseName = fi.baseName();
      baseName.remove(getInputFilesList().FileSuffix);
      baseName.remove(getInputFilesList().FilePrefix);
      if(!baseName.isEmpty())
      {
        if(getInputFilesList().PaddingDigits > 1)
        {
          bool initialZero = (baseName.at(0) == "0");
          if(!initialZero)
          {
            fileIndices.insert(baseName.toUInt());
          }
          else
          {
            while(baseName.at(0) == "0")
            {
              baseName.remove(0, 1);
            }
            if(baseName.isEmpty())
            {
              fileIndices.insert(0);
            }
            else
            {
              fileIndices.insert(baseName.toUInt());
            }
          }
        }
        else
        {
          fileIndices.insert(baseName.toUInt());
        }
      }
    }
  }

  if(!missing_fnames.empty())
  {
    QString ss;
    QTextStream out(&ss);
    out << " File name options have resulted in missing files. Have you correctly set the start and end indices; file prefix, suffix, and extension; and padding digits for your files?\n";
    out << "InputPath: " << getInputFilesList().InputPath << "\n";
    out << "FilePrefix: " << getInputFilesList().FilePrefix << "\n";
    out << "FileSuffix: " << getInputFilesList().FileSuffix << "\n";
    out << "FileExtension: " << getInputFilesList().FileExtension << "\n";
    out << "PaddingDigits: " << getInputFilesList().PaddingDigits << "\n";
    out << "StartIndex: " << getInputFilesList().StartIndex << "\n";
    out << "EndIndex: " << getInputFilesList().EndIndex << "\n";
    out << "Missing Files:\n";

    for(auto&& fname : fileList)
    {
      out << fname << "\n";
    }

    setErrorCondition(-12, ss);
  }

  switch(getLaserOnArrayOption())
  {
  case 0:
  {
    m_LaserOnArrayName = "Laser Drive";
    break;
  }
  case 1:
  {
    m_LaserOnArrayName = "On-Axis PD";
    break;
  }
  default:
  {
    QString ss = QObject::tr("Invalid selection for array used to determine if laser is on");
    setErrorCondition(-701, ss);
    break;
  }
  }

  switch(getSpatialTransformOption())
  {
  case 0: // Not computing transform, this is a no-op
  {
    break;
  }
  case 1: // Computing transform, check other input parameters
  {
    QFileInfo stl(getSTLFilePath1());

    if(getSTLFilePath1().isEmpty())
    {
      QString ss = QObject::tr("The input build STL file must be set");
      setErrorCondition(-387, ss);
    }
    else if(!stl.exists() || !stl.isFile())
    {
      QString ss = QObject::tr("The input build STL file does not exist");
      setErrorCondition(-388, ss);
    }

    if(fileIndices.count(getLayerForScaling()) == 0)
    {
      QString ss = QObject::tr("Layer index selected for computing the spatial transform (%1) is not present in the supplied file list").arg(getLayerForScaling());
      setErrorCondition(-388, ss);
    }

    if(getSearchRadius() <= 0.0f)
    {
      QString ss = QObject::tr("The search radius must be greater than zero");
      setErrorCondition(-389, ss);
    }

    FilterManager* fm = FilterManager::Instance();
    IFilterFactory::Pointer factory = fm->getFactoryFromClassName("ReadStlFile");
    if(factory.get() == nullptr)
    {
      QString ss = QObject::tr("Unable to instantiate Filter with name 'ReadStlFile'\n"
                               "The 'ReadStlFile' Filter is needed to import the build STL file");
      setErrorCondition(-1, ss);
    }
    factory = fm->getFactoryFromClassName("DBSCAN");
    if(factory.get() == nullptr)
    {
      QString ss = QObject::tr("Unable to instantiate Filter with name 'DBSCAN'\n"
                               "The 'DBSCAN' Filter is needed to automatically segment regions in the PrintRite data");
      setErrorCondition(-1, ss);
    }
    break;
  }
  case 2: // Importing transform from file; need to check file exists + parse it for parameters
  {
    QFileInfo stl(getSTLFilePath2());

    if(getSTLFilePath2().isEmpty())
    {
      QString ss = QObject::tr("The input build STL file must be set");
      setErrorCondition(-387, ss);
    }
    else if(!stl.exists() || !stl.isFile())
    {
      QString ss = QObject::tr("The input build STL file does not exist");
      setErrorCondition(-388, ss);
    }

    QFileInfo fi(getInputSpatialTransformFilePath());

    if(getInputSpatialTransformFilePath().isEmpty())
    {
      QString ss = QObject::tr("The input spatial transform coefficients file must be set");
      setErrorCondition(-387, ss);
    }
    else if(!fi.exists() || !fi.isFile())
    {
      QString ss = QObject::tr("The input spatial transform coefficients file does not exist");
      setErrorCondition(-388, ss);
    }

    if(getErrorCode() < 0)
    {
      return;
    }

    QFile spatialTransforms(getInputSpatialTransformFilePath());

    if(!spatialTransforms.open(QIODevice::ReadOnly))
    {
      QString ss = QObject::tr("The input spatial transform coefficients file could not be opened");
      setErrorCondition(-388, ss);
    }

    QByteArray transformData = spatialTransforms.readAll();
    QJsonDocument jsonDoc(QJsonDocument::fromJson(transformData));

    if(jsonDoc.isNull())
    {
      QString ss = QObject::tr("The input spatial transform coefficients file is not a valid JSON file");
      setErrorCondition(-388, ss);
    }

    QJsonValue xCoefficients = jsonDoc.object()["X Coefficients"];
    QJsonValue yCoefficients = jsonDoc.object()["Y Coefficients"];
    QJsonArray xArray = xCoefficients.toArray();
    QJsonArray yArray = yCoefficients.toArray();
    Eigen::VectorXd xCoeffs;
    Eigen::VectorXd yCoeffs;
    xCoeffs.resize(10);
    yCoeffs.resize(10);
    for(auto i = 0; i < xArray.size(); i++)
    {
      xCoeffs(i) = xArray[i].toDouble();
      yCoeffs(i) = yArray[i].toDouble();
    }

    std::vector<Eigen::VectorXd> coefficients;
    coefficients.push_back(xCoeffs);
    coefficients.push_back(yCoeffs);
    m_Polynomial.setOrder(3);
    m_Polynomial.setCoefficients(coefficients);

    break;
  }
  default:
  {
    QString ss = QObject::tr("Invalid selection for spatial transform option");
    setErrorCondition(-701, ss);
    break;
  }
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ImportPrintRiteTDMSFiles::importLabelSliceSTL()
{
  QString ss = QObject::tr("Importing Build STL File...");
  notifyStatusMessage(ss);
  m_LocalStructure = DataContainerArray::New();
  FilterManager* fm = FilterManager::Instance();
  IFilterFactory::Pointer factory = fm->getFactoryFromClassName("ReadStlFile");
  if(factory.get() == nullptr)
  {
    ss = QObject::tr("Unable to instantiate Filter with name 'ReadStlFile'\n"
                     "The 'ReadStlFile' Filter is needed to import the build STL file");
    setErrorCondition(-1, ss);
    return;
  }
  AbstractFilter::Pointer importSTLFilter = factory->create();
  importSTLFilter->setDataContainerArray(m_LocalStructure);
  QVariant var;
  var.setValue(m_STLFilePath);
  importSTLFilter->setProperty("StlFilePath", var);
  importSTLFilter->execute();
  if(importSTLFilter->getErrorCode() < 0)
  {
    ss = QObject::tr("Error in filter 'ReadStlFile'");
    setErrorCondition(importSTLFilter->getErrorCode(), ss);
    return;
  }

  ss = QObject::tr("Identifying Individual Parts...");
  notifyStatusMessage(ss);
  factory = fm->getFactoryFromClassName("LabelTriangleGeometry");
  if(factory.get() == nullptr)
  {
    ss = QObject::tr("Unable to instantiate Filter with name 'LabelTriangleGeometry'\n"
                     "The 'LabelTriangleGeometry' Filter is needed to label the STL file regions");
    setErrorCondition(-1, ss);
    return;
  }
  AbstractFilter::Pointer labelFilter = factory->create();
  labelFilter->setDataContainerArray(m_LocalStructure);
  var.setValue(SIMPL::Defaults::TriangleDataContainerName);
  labelFilter->setProperty("CADDataContainerName", var);
  var.setValue(QString("_INTERNAL_USE_ONLY_RegionIdsAttributeMatrix"));
  labelFilter->setProperty("TriangleAttributeMatrixName", var);
  var.setValue(QString("_INTERNAL_USE_ONLY_RegionIds"));
  labelFilter->setProperty("RegionIdArrayName", var);
  labelFilter->execute();
  if(labelFilter->getErrorCode() < 0)
  {
    ss = QObject::tr("Error in filter 'LabelTriangleGeometry'");
    setErrorCondition(labelFilter->getErrorCode(), ss);
    return;
  }
  Int32ArrayType::Pointer regionIds = m_LocalStructure->getDataContainer(SIMPL::Defaults::TriangleDataContainerName)
                                          ->getAttributeMatrix("_INTERNAL_USE_ONLY_RegionIdsAttributeMatrix")
                                          ->getAttributeArrayAs<Int32ArrayType>("_INTERNAL_USE_ONLY_RegionIds");
  int32_t* regionIdsPtr = regionIds->getPointer(0);
  for(size_t i = 0; i < regionIds->getNumberOfTuples(); i++)
  {
    if(regionIdsPtr[i] > m_NumParts)
    {
      m_NumParts = regionIdsPtr[i];
    }
  }
  ss = QObject::tr("Found %1 Parts").arg(m_NumParts);
  notifyStatusMessage(ss);

  ss = QObject::tr("Finding Build Geometry Layers...");
  notifyStatusMessage(ss);
  factory = fm->getFactoryFromClassName("SliceTriangleGeometry");
  if(factory.get() == nullptr)
  {
    ss = QObject::tr("Unable to instantiate Filter with name 'SliceTriangleGeometry'\n"
                     "The 'SliceTriangleGeometry' Filter is needed to slice the STL file");
    setErrorCondition(-1, ss);
    return;
  }
  AbstractFilter::Pointer sliceFilter = factory->create();
  connect(sliceFilter.get(), SIGNAL(filterGeneratedMessage(const PipelineMessage&)), this, SLOT(broadcastPipelineMessage(const PipelineMessage&)));
  sliceFilter->setDataContainerArray(m_LocalStructure);
  FloatVec3Type sliceDir;
  sliceDir[0] = 0.0f;
  sliceDir[1] = 0.0f;
  sliceDir[2] = 1.0f;
  var.setValue(sliceDir);
  sliceFilter->setProperty("SliceDirection", var);
  var.setValue(m_LayerThickness);
  sliceFilter->setProperty("SliceResolution", var);
  var.setValue(true);
  sliceFilter->setProperty("HaveRegionIds", var);
  var.setValue(SIMPL::Defaults::TriangleDataContainerName);
  sliceFilter->setProperty("CADDataContainerName", var);
  DataArrayPath path(SIMPL::Defaults::TriangleDataContainerName, "_INTERNAL_USE_ONLY_RegionIdsAttributeMatrix", "_INTERNAL_USE_ONLY_RegionIds");
  var.setValue(path);
  sliceFilter->setProperty("RegionIdArrayPath", var);
  var.setValue(SIMPL::Defaults::EdgeDataContainerName);
  sliceFilter->setProperty("SliceDataContainerName", var);
  var.setValue(SIMPL::Defaults::EdgeAttributeMatrixName);
  sliceFilter->setProperty("EdgeAttributeMatrixName", var);
  var.setValue(QString("_INTERNAL_USE_ONLY_SliceData"));
  sliceFilter->setProperty("SliceAttributeMatrixName", var);
  var.setValue(QString("_INTERNAL_USE_ONLY_SliceIds"));
  sliceFilter->setProperty("SliceIdArrayName", var);
  sliceFilter->execute();
  if(sliceFilter->getErrorCode() < 0)
  {
    ss = QObject::tr("Error in filter 'SliceCADFile'");
    setErrorCondition(sliceFilter->getErrorCode(), ss);
    return;
  }

  Int32ArrayType::Pointer sliceIds = m_LocalStructure->getDataContainer(SIMPL::Defaults::EdgeDataContainerName)
                                         ->getAttributeMatrix(SIMPL::Defaults::EdgeAttributeMatrixName)
                                         ->getAttributeArrayAs<Int32ArrayType>("_INTERNAL_USE_ONLY_SliceIds");
  int32_t* sliceIdsPtr = sliceIds->getPointer(0);
  for(size_t i = 0; i < sliceIds->getNumberOfTuples(); i++)
  {
    if(sliceIdsPtr[i] > m_NumLayers)
    {
      m_NumLayers = sliceIdsPtr[i];
    }
  }
  ss = QObject::tr("Found %1 Total Layers in Build Geometry").arg(m_NumLayers);
  notifyStatusMessage(ss);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ImportPrintRiteTDMSFiles::createHDF5Files()
{
  QDir dir;
  if(!dir.mkpath(m_OutputDirectory))
  {
    QString ss = QObject::tr("Error creating output path '%1'").arg(m_OutputDirectory);
    setErrorCondition(-11110, ss);
    return;
  }

  if(m_SpatialTransformOption == 0 || (m_SpatialTransformOption == 1 && !m_SplitRegions1) || (m_SpatialTransformOption == 2 && !m_SplitRegions2))
  {
    QString fname = (m_OutputDirectory) + QDir::separator() + (m_OutputFilePrefix) + QString::number(1) + ".hdf5";
    hid_t fileId = QH5Utilities::createFile(fname);
    m_RegionFileMap[1] = fileId;
  }
  else
  {
    for(int32_t i = 0; i < m_NumParts + 1; i++)
    {
      QString fname = (m_OutputDirectory) + QDir::separator() + (m_OutputFilePrefix) + QString::number(i) + ".hdf5";
      hid_t fileId = QH5Utilities::createFile(fname);
      m_RegionFileMap[i] = fileId;
    }
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ImportPrintRiteTDMSFiles::closeHDF5Files()
{
  for(auto&& pair : m_RegionFileMap)
  {
    hid_t tmp = pair.second;
    if(tmp > 0)
    {
      H5Utilities::closeFile(tmp);
      pair.second = -1;
    }
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ImportPrintRiteTDMSFiles::writeMetaDataToHDF5Files()
{
  // TODO: also write h/l frequencies?
  std::vector<IDataArray::Pointer> buildMetaData;
  std::vector<IDataArray::Pointer> scalingMetaData;

  FloatArrayType::Pointer layerThickness = FloatArrayType::CreateArray(1, std::string("Layer Thickness"), true);
  layerThickness->setValue(0, m_LayerThickness);
  buildMetaData.push_back(layerThickness);

  FloatArrayType::Pointer laserDriveThreshold = FloatArrayType::CreateArray(1, std::string("Laser Drive Threshold"), true);
  laserDriveThreshold->setValue(0, m_LaserOnThreshold);
  buildMetaData.push_back(laserDriveThreshold);

  std::vector<size_t> cDims(1, 2);
  FloatArrayType::Pointer powerScales = FloatArrayType::CreateArray(1, cDims, "Laser Power Scaling Coefficients", true);
  powerScales->setComponent(0, 0, m_PowerScalingCoefficients[0]);
  powerScales->setComponent(0, 1, m_PowerScalingCoefficients[1]);
  scalingMetaData.push_back(powerScales);

  FloatArrayType::Pointer temperatureScales = FloatArrayType::CreateArray(1, cDims, "Pyrometer Temperature Scaling Coefficients", true);
  temperatureScales->setComponent(0, 0, m_TemperatureScalingCoefficients[0]);
  temperatureScales->setComponent(0, 1, m_TemperatureScalingCoefficients[1]);
  scalingMetaData.push_back(temperatureScales);

  DoubleArrayType::Pointer polyCoefficients = nullptr;
  if(m_SpatialTransformOption == 0)
  {
    cDims = {2};
    polyCoefficients = DoubleArrayType::CreateArray(10, cDims, "Spatial Scaling Coefficients", true);
    polyCoefficients->initializeWithZeros();
  }
  else
  {
    polyCoefficients = m_Polynomial.coefficientsAsDataArray();
  }
  scalingMetaData.push_back(polyCoefficients);

  for(auto&& pair : m_RegionFileMap)
  {
    hid_t buildMetaDataGroupId = QH5Utilities::createGroup(pair.second, "Build Meta Data");
    hid_t scalingMetaDataGroupId = QH5Utilities::createGroup(pair.second, "Scaling Meta Data");
    for(auto&& data : buildMetaData)
    {
      std::vector<size_t> tDims(1, data->getNumberOfTuples());
      data->writeH5Data(buildMetaDataGroupId, tDims);
    }
    for(auto&& data : scalingMetaData)
    {
      std::vector<size_t> tDims(1, data->getNumberOfTuples());
      data->writeH5Data(scalingMetaDataGroupId, tDims);
    }
    QH5Utilities::closeHDF5Object(buildMetaDataGroupId);
    QH5Utilities::closeHDF5Object(scalingMetaDataGroupId);
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ImportPrintRiteTDMSFiles::writeSpatialTransformToFile()
{
  QString fname = m_OutputDirectory + QDir::separator() + m_OutputFilePrefix + "SpatialTransformCoefficients.json";
  QFile file(fname);
  if(!file.open(QIODevice::WriteOnly))
  {
    QString ss = QObject::tr("Output file could not be opened: %1").arg(fname);
    setErrorCondition(-100, ss);
    return;
  }

  std::vector<Eigen::VectorXd> coefficients = m_Polynomial.coefficients();
  QJsonArray xCoefficients;
  QJsonArray yCoefficients;
  for(auto i = 0; i < coefficients[0].size(); i++)
  {
    xCoefficients.push_back(coefficients[0](i));
  }
  for(auto i = 0; i < coefficients[1].size(); i++)
  {
    yCoefficients.push_back(coefficients[1](i));
  }
  QJsonObject root;
  root["X Coefficients"] = xCoefficients;
  root["Y Coefficients"] = yCoefficients;

  QJsonDocument jsonDoc(root);
  file.write(jsonDoc.toJson());
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ImportPrintRiteTDMSFiles::computeSpatialTransformation(const QString& fname)
{
  TDMSFileProxy::Pointer proxy = nullptr;
  try
  {
    proxy = TDMSFileProxy::New(fname.toStdString());
    proxy->readMetaData();
    proxy->allocateObjects();
    proxy->readRawData();
  } catch(const FatalTDMSException& exc)
  {
    std::string tmp = exc.getMessage();
    QString msg = QString::fromStdString(exc.getMessage());
    setErrorCondition(-1, msg);
    return;
  }

  std::unordered_map<std::string, TDMSObject::Pointer> channels = proxy->channelObjects();

  if(channels.count("X Position") == 0)
  {
    QString ss = QObject::tr("Array named 'X Position' was not found in the following file:\n"
                             "%1\n"
                             "The x position array is required for determining the spatial layout of the build")
                     .arg(fname);
    setErrorCondition(-1, ss);
    return;
  }

  if(channels.count("Y Position") == 0)
  {
    QString ss = QObject::tr("Array named 'Y Position' was not found in the following file:\n"
                             "%1\n"
                             "The y position array is required for determining the spatial layout of the build")
                     .arg(fname);
    setErrorCondition(-1, ss);
    return;
  }

  if(channels.count(m_LaserOnArrayName) == 0)
  {
    QString ss = QObject::tr("Array named '%1' was not found in the following file:\n"
                             "%2\n"
                             "This array is required for determining the spatial layout of the build (i.e., where the laser is off)")
                     .arg(QString::fromStdString(m_LaserOnArrayName))
                     .arg(fname);
    setErrorCondition(-1, ss);
    return;
  }

  TDMSObject::Pointer xpos = channels["X Position"];
  TDMSObject::Pointer ypos = channels["Y Position"];
  TDMSObject::Pointer drive = channels[m_LaserOnArrayName];
  DoubleArrayType::Pointer xposArray = std::dynamic_pointer_cast<DoubleArrayType>(xpos->data());
  DoubleArrayType::Pointer yposArray = std::dynamic_pointer_cast<DoubleArrayType>(ypos->data());
  DoubleArrayType::Pointer driveArray = std::dynamic_pointer_cast<DoubleArrayType>(drive->data());
  double* xposPtr = xposArray->getPointer(0);
  double* yposPtr = yposArray->getPointer(0);
  double* drivePtr = driveArray->getPointer(0);

  if(xposArray->getNumberOfTuples() != yposArray->getNumberOfTuples())
  {
    QString ss = QObject::tr("The 'X Position' and 'Y Position' arrays do not have the same length\n"
                             "X Position Length: %1\n"
                             "Y Position Length: %2\n"
                             "File: %3")
                     .arg(xposArray->getNumberOfTuples())
                     .arg(yposArray->getNumberOfTuples())
                     .arg(fname);
    setErrorCondition(-1, ss);
    return;
  }

  if(xposArray->getNumberOfTuples() != driveArray->getNumberOfTuples())
  {
    QString ss = QObject::tr("The 'Laser Drive' and position arrays do not have the same length\n"
                             "Laser Drive Length: %1\n"
                             "Position Length: %2\n"
                             "File: %3")
                     .arg(driveArray->getNumberOfTuples())
                     .arg(xposArray->getNumberOfTuples())
                     .arg(fname);
    setErrorCondition(-1, ss);
    return;
  }

  size_t numBeamOnPoints = 0;
  for(size_t i = 0; i < xposArray->getNumberOfTuples(); i++)
  {
    if(drivePtr[i] > m_LaserOnThreshold)
    {
      numBeamOnPoints++;
    }
  }

  std::vector<size_t> cDims(1, 2);
  FloatArrayType::Pointer beamOnTDMS = FloatArrayType::CreateArray(numBeamOnPoints, cDims, "_INTERNAL_USE_ONLY_BeamOnTDMS", true);
  float* beamOnTDMSPtr = beamOnTDMS->getPointer(0);
  numBeamOnPoints = 0;
  for(size_t i = 0; i < xposArray->getNumberOfTuples(); i++)
  {
    if(drivePtr[i] > m_LaserOnThreshold)
    {
      beamOnTDMSPtr[2 * numBeamOnPoints + 0] = static_cast<float>(xposPtr[i]);
      beamOnTDMSPtr[2 * numBeamOnPoints + 1] = static_cast<float>(-yposPtr[i]);
      numBeamOnPoints++;
    }
  }

  float fractionForRegistration = 0.05f; // TODO: Make this a user option?
  size_t numPointsForRegistration = static_cast<size_t>(fractionForRegistration * numBeamOnPoints);
  BoolArrayType::Pointer mask = BoolArrayType::CreateArray(numBeamOnPoints, std::string("_INTERNAL_USE_ONLY_Mask"), true);
  mask->initializeWithValue(false);
  bool* maskPtr = mask->getPointer(0);
  size_t counter = 0;
  size_t rangeMin = 0;
  size_t rangeMax = numBeamOnPoints - 1;
  std::mt19937_64::result_type seed = static_cast<std::mt19937_64::result_type>(std::chrono::steady_clock::now().time_since_epoch().count());
  std::mt19937_64 gen(seed);
  std::uniform_int_distribution<size_t> dist(rangeMin, rangeMax);
  while(counter < numPointsForRegistration)
  {
    size_t index = dist(gen);
    if(!maskPtr[index])
    {
      maskPtr[index] = true;
      counter++;
    }
  }

  QString ss = QObject::tr("Finding regions in PrintRite data...");
  notifyStatusMessage(ss);

  FilterManager* fm = FilterManager::Instance();
  IFilterFactory::Pointer factory = fm->getFactoryFromClassName("DBSCAN");
  if(factory.get() == nullptr)
  {
    QString ss = QObject::tr("Unable to instantiate Filter with name 'DBSCAN'\n"
                             "The 'DBSCAN' Filter is needed to automatically segment regions in the PrintRite data");
    setErrorCondition(-1, ss);
    return;
  }
  DataContainer::Pointer dc = DataContainer::New("_INTERNAL_USE_ONLY_DBSCAN");
  std::vector<size_t> tDims(1, numBeamOnPoints);
  AttributeMatrix::Pointer am = AttributeMatrix::New(tDims, "_INTERNAL_USE_ONLY_DBSCANData", AttributeMatrix::Type::Vertex);
  am->addOrReplaceAttributeArray(mask);
  am->addOrReplaceAttributeArray(beamOnTDMS);
  dc->addOrReplaceAttributeMatrix(am);
  m_LocalStructure->addOrReplaceDataContainer(dc);
  AbstractFilter::Pointer dbscanFilter = factory->create();
  // TODO: fix progress output for DBSCAN
  // connect(dbscanFilter.get(), SIGNAL(filterGeneratedMessage(const PipelineMessage&)), this, SLOT(broadcastPipelineMessage(const PipelineMessage&)));
  dbscanFilter->setDataContainerArray(m_LocalStructure);
  QVariant var;
  var.setValue(m_SearchRadius);
  dbscanFilter->setProperty("Epsilon", var);
  var.setValue(20);
  dbscanFilter->setProperty("MinPnts", var);
  var.setValue(1);
  dbscanFilter->setProperty("DistanceMetric", var);
  var.setValue(true);
  dbscanFilter->setProperty("UseMask", var);
  DataArrayPath path("_INTERNAL_USE_ONLY_DBSCAN", "_INTERNAL_USE_ONLY_DBSCANData", "_INTERNAL_USE_ONLY_Mask");
  var.setValue(path);
  dbscanFilter->setProperty("MaskArrayPath", var);
  path.setDataArrayName("_INTERNAL_USE_ONLY_BeamOnTDMS");
  var.setValue(path);
  dbscanFilter->setProperty("SelectedArrayPath", var);
  dbscanFilter->execute();
  notifyStatusMessage(ss);
  if(dbscanFilter->getErrorCode() < 0)
  {
    QString ss = QObject::tr("Error in filter 'DBSCAN'");
    setErrorCondition(dbscanFilter->getErrorCode(), ss);
    return;
  }

  ss = QObject::tr("Determining rigid body transformation...");
  notifyStatusMessage(ss);

  Int32ArrayType::Pointer clusterIds =
      m_LocalStructure->getDataContainer("_INTERNAL_USE_ONLY_DBSCAN")->getAttributeMatrix("_INTERNAL_USE_ONLY_DBSCANData")->getAttributeArrayAs<Int32ArrayType>("ClusterIds");
  int32_t* clusterIdsPtr = clusterIds->getPointer(0);
  int32_t numClusters = std::numeric_limits<int32_t>::lowest();
  for(size_t i = 0; i < clusterIds->getNumberOfTuples(); i++)
  {
    if(clusterIdsPtr[i] > numClusters)
    {
      numClusters = clusterIdsPtr[i];
    }
  }
  ss = QObject::tr("Found %1 regions in layer %2").arg(numClusters).arg(m_LayerForScaling);
  notifyStatusMessage(ss);

  auto mostCentralCluster = findCentralCluster(beamOnTDMS, clusterIds);
  PRH::Polygons woundPolygons = extractLayerPolygons(m_LayerForScaling + m_Offset);
  auto mostCentralPolygon = findCentralPolygon(woundPolygons, m_LayerForScaling + m_Offset, fname);

  float extentsTDMS[2] = {mostCentralCluster.second.xmax - mostCentralCluster.second.xmin, mostCentralCluster.second.ymax - mostCentralCluster.second.ymin};
  float extentsSlice[2] = {mostCentralPolygon.second.xmax - mostCentralPolygon.second.xmin, mostCentralPolygon.second.ymax - mostCentralPolygon.second.ymin};
  float scale[2] = {extentsSlice[0] / extentsTDMS[0], extentsSlice[1] / extentsTDMS[1]};
  float centerTDMS[2] = {(((mostCentralCluster.second.xmax * scale[0]) - (mostCentralCluster.second.xmin * scale[0])) / 2) + (mostCentralCluster.second.xmin * scale[0]),
                         (((mostCentralCluster.second.ymax * scale[1]) - (mostCentralCluster.second.ymin * scale[1])) / 2) + (mostCentralCluster.second.ymin * scale[1])};
  float centerSlice[2] = {(extentsSlice[0] / 2) + mostCentralPolygon.second.xmin, (extentsSlice[1] / 2) + mostCentralPolygon.second.ymin};
  float translation[2] = {centerSlice[0] - centerTDMS[0], centerSlice[1] - centerTDMS[1]};

  FloatArrayType::Pointer rigidTransformBeamOnTDMS = std::dynamic_pointer_cast<FloatArrayType>(beamOnTDMS->deepCopy());
  float* rigidTransformBeamOnTDMSPtr = rigidTransformBeamOnTDMS->getPointer(0);
  for(size_t i = 0; i < numBeamOnPoints; i++)
  {
    rigidTransformBeamOnTDMSPtr[2 * i + 0] = beamOnTDMSPtr[2 * i + 0] * scale[0] + translation[0];
    rigidTransformBeamOnTDMSPtr[2 * i + 1] = beamOnTDMSPtr[2 * i + 1] * scale[1] + translation[1];
  }

  ss = QObject::tr("Found scale: [%1, %2]").arg(scale[0]).arg(scale[1]);
  notifyStatusMessage(ss);

  ss = QObject::tr("Found translation: [%1, %2]").arg(translation[0]).arg(translation[1]);
  notifyStatusMessage(ss);

  ss = QObject::tr("Associating points with polygons...");
  notifyStatusMessage(ss);

  auto associatedPointsPolys = associatePointsWithPolygons(woundPolygons, rigidTransformBeamOnTDMS);

  // TEST
  // DataContainer::Pointer tmpdc = DataContainer::New("TestRigidTransform");
  // VertexGeom::Pointer tmpverts = VertexGeom::CreateGeometry(rigidTransformBeamOnTDMS->getNumberOfTuples(), SIMPL::Geometry::VertexGeometry);
  // float* tmpvertsptr = tmpverts->getVertexPointer(0);
  // for(int64_t i = 0; i < tmpverts->getNumberOfVertices(); i++)
  //{
  //  tmpvertsptr[3 * i + 0] = rigidTransformBeamOnTDMSPtr[2 * i + 0];
  //  tmpvertsptr[3 * i + 1] = rigidTransformBeamOnTDMSPtr[2 * i + 1];
  //  tmpvertsptr[3 * i + 2] = 0.0f;
  //}
  // tmpdc->setGeometry(tmpverts);
  // getDataContainerArray()->addOrReplaceDataContainer(tmpdc);
  // TEST

  ss = QObject::tr("Computing polynomial warp...");
  notifyStatusMessage(ss);

  determinePointsForLeastSquares(woundPolygons, beamOnTDMS, associatedPointsPolys, mask);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ImportPrintRiteTDMSFiles::processLayers(const QVector<QString>& files)
{
  int32_t layerIndex = m_InputFilesList.StartIndex + m_Offset;
  size_t counter = 1;

  // TODO: create hf/lf/unknown groups + import lf data as well
  //      any arrays whose name is not recognized --> stick in unknown group

  for(auto&& fname : files)
  {
    QString ss = QObject::tr("Importing TDMS Layer %1 (%2 of %3)").arg(layerIndex).arg(counter).arg(m_NumLayersToImport);
    notifyStatusMessage(ss);

    TDMSFileProxy::Pointer proxy = nullptr;
    try
    {
      proxy = TDMSFileProxy::New(fname.toStdString());
      proxy->readMetaData();
      proxy->allocateObjects();
      proxy->readRawData();
    } catch(const FatalTDMSException& exc)
    {
      QString msg = QString::fromStdString(exc.getMessage());
      setErrorCondition(-1, msg);
      return;
    }

    std::unordered_map<std::string, TDMSObject::Pointer> channels = proxy->channelObjects();
    PRH::PrintRiteChannels printRiteChannels;
    std::vector<IDataArray::Pointer> hfArraysToWrite;
    std::vector<IDataArray::Pointer> lfArraysToWrite;
    std::vector<IDataArray::Pointer> unknownArraysToWrite;
    if(m_DowncastRawData)
    {
      std::vector<FloatArrayType::Pointer> downcastHfArrays = printRiteChannels.castChannelsTo<float>(PRH::PrintRiteChannels::ChannelType::HF, channels);
      std::vector<FloatArrayType::Pointer> downcastLfArrays = printRiteChannels.castChannelsTo<float>(PRH::PrintRiteChannels::ChannelType::LF, channels);
      std::vector<FloatArrayType::Pointer> downcastUnknownArrays = printRiteChannels.castChannelsTo<float>(PRH::PrintRiteChannels::ChannelType::Unknown, channels);
      hfArraysToWrite.insert(std::end(hfArraysToWrite), std::begin(downcastHfArrays), std::end(downcastHfArrays));
      lfArraysToWrite.insert(std::end(lfArraysToWrite), std::begin(downcastLfArrays), std::end(downcastLfArrays));
      unknownArraysToWrite.insert(std::end(unknownArraysToWrite), std::begin(downcastUnknownArrays), std::end(downcastUnknownArrays));
    }
    else
    {
      std::vector<IDataArray::Pointer> hfArrays = printRiteChannels.getChannelsOfType(PRH::PrintRiteChannels::ChannelType::HF, channels);
      std::vector<IDataArray::Pointer> lfArrays = printRiteChannels.getChannelsOfType(PRH::PrintRiteChannels::ChannelType::HF, channels);
      std::vector<IDataArray::Pointer> uknownArrays = printRiteChannels.getChannelsOfType(PRH::PrintRiteChannels::ChannelType::Unknown, channels);
      hfArraysToWrite.insert(std::end(hfArraysToWrite), std::begin(hfArrays), std::end(hfArrays));
      lfArraysToWrite.insert(std::end(lfArraysToWrite), std::begin(lfArrays), std::end(lfArrays));
      unknownArraysToWrite.insert(std::end(unknownArraysToWrite), std::begin(uknownArrays), std::end(uknownArrays));
    }
    std::vector<std::string> laserOnName = {m_LaserOnArrayName};
    std::vector<BoolArrayType::Pointer> thresholdHfArrays = printRiteChannels.thresholdHfChannels(channels, laserOnName, m_LaserOnThreshold);
    hfArraysToWrite.insert(std::end(hfArraysToWrite), std::begin(thresholdHfArrays), std::end(thresholdHfArrays));

    if(m_SpatialTransformOption == 0 || (m_SpatialTransformOption == 1 && !m_SplitRegions1) || (m_SpatialTransformOption == 2 && !m_SplitRegions2))
    {
      hid_t fileId = -1;
      try
      {
        fileId = m_RegionFileMap.at(1);
      } catch(const std::out_of_range& oor)
      {
        QString msg = oor.what();
        setErrorCondition(-1, msg);
        return;
      }

      hid_t layerDataGroupId = QH5Utilities::createGroup(fileId, "Layer Data");
      hid_t layerGroupId = QH5Utilities::createGroup(layerDataGroupId, QString::number(layerIndex));
      hid_t hfGroupId = QH5Utilities::createGroup(layerGroupId, "High Frequency Data");
      hid_t lfGroupId = QH5Utilities::createGroup(layerGroupId, "Low Frequency Data");
      hid_t unknownGroupId = QH5Utilities::createGroup(layerGroupId, "Unassociated Data");
      for(auto&& dataArray : hfArraysToWrite)
      {
        std::vector<size_t> tDims(1, dataArray->getNumberOfTuples());
        dataArray->writeH5Data(hfGroupId, tDims);
      }
      for(auto&& dataArray : lfArraysToWrite)
      {
        std::vector<size_t> tDims(1, dataArray->getNumberOfTuples());
        dataArray->writeH5Data(lfGroupId, tDims);
      }
      for(auto&& dataArray : unknownArraysToWrite)
      {
        std::vector<size_t> tDims(1, dataArray->getNumberOfTuples());
        dataArray->writeH5Data(unknownGroupId, tDims);
      }
      QH5Utilities::closeHDF5Object(layerDataGroupId);
      QH5Utilities::closeHDF5Object(layerGroupId);
      QH5Utilities::closeHDF5Object(hfGroupId);
      QH5Utilities::closeHDF5Object(lfGroupId);
      QH5Utilities::closeHDF5Object(unknownGroupId);
    }
    else
    {
      TDMSObject::Pointer xpos = channels["X Position"];
      TDMSObject::Pointer ypos = channels["Y Position"];
      DoubleArrayType::Pointer xposArray = std::dynamic_pointer_cast<DoubleArrayType>(xpos->data());
      DoubleArrayType::Pointer yposArray = std::dynamic_pointer_cast<DoubleArrayType>(ypos->data());
      double* xposPtr = xposArray->getPointer(0);
      double* yposPtr = yposArray->getPointer(0);

      std::vector<size_t> cDims(1, 2);
      FloatArrayType::Pointer tdmsPts = FloatArrayType::CreateArray(xposArray->getNumberOfTuples(), cDims, "_INTERNAL_USE_ONLY_TDMSPreScale", true);
      float* tdmsPtsPtr = tdmsPts->getPointer(0);
      for(size_t i = 0; i < xposArray->getNumberOfTuples(); i++)
      {
        float pt[2] = {static_cast<float>(xposPtr[i]), static_cast<float>(-yposPtr[i])};
        tdmsPtsPtr[2 * i + 0] = m_Polynomial.transformPoint(pt, 0);
        tdmsPtsPtr[2 * i + 1] = m_Polynomial.transformPoint(pt, 1);
      }

      PRH::Polygons woundPolygons = extractLayerPolygons(layerIndex);
      auto associatedPointsPolys = associatePointsWithPolygons(woundPolygons, tdmsPts);

      std::vector<size_t> numPointsForPoly(m_NumParts + 1, 0);
      Int32ArrayType::Pointer pointsToPolys = associatedPointsPolys.first;
      std::vector<bool> validPolys = associatedPointsPolys.second;
      validPolys.insert(std::begin(validPolys), true);
      int32_t* pointsToPolysPtr = pointsToPolys->getPointer(0);
      for(size_t i = 0; i < pointsToPolys->getNumberOfTuples(); i++)
      {
        numPointsForPoly[pointsToPolysPtr[i] + 1]++;
      }

      std::vector<std::vector<IDataArray::Pointer>> splitArraysToWrite(m_NumParts + 1);
      for(auto i = 0; i < splitArraysToWrite.size(); i++)
      {
        if(!validPolys[i])
        {
          continue;
        }
        for(auto&& it : hfArraysToWrite)
        {
          if(numPointsForPoly[i] > 0)
          {
            IDataArray::Pointer ptr = TemplateHelpers::CreateArrayFromArrayType()(this, numPointsForPoly[i], it->getComponentDimensions(), it->getName(), true, it);
            splitArraysToWrite[i].push_back(ptr);
          }
        }
      }

      std::vector<std::vector<size_t>> arrayCounters(m_NumParts + 1, std::vector<size_t>(hfArraysToWrite.size(), 0));
      for(size_t i = 0; i < pointsToPolys->getNumberOfTuples(); i++)
      {
        std::vector<IDataArray::Pointer> arraysToWriteForPoly = splitArraysToWrite[pointsToPolysPtr[i] + 1];
        for(auto j = 0; j < arraysToWriteForPoly.size(); j++)
        {
          IDataArray::Pointer source = hfArraysToWrite[j];
          IDataArray::Pointer destination = arraysToWriteForPoly[j];
          std::memcpy(destination->getVoidPointer(arrayCounters[pointsToPolysPtr[i] + 1][j]), source->getVoidPointer(i), source->getTypeSize());
          arrayCounters[pointsToPolysPtr[i] + 1][j]++;
        }
      }

      for(auto i = 0; i < splitArraysToWrite.size(); i++)
      {
        hid_t fileId = -1;
        try
        {
          fileId = m_RegionFileMap.at(i);
        } catch(const std::out_of_range& oor)
        {
          QString msg = oor.what();
          setErrorCondition(-1, msg);
          return;
        }

        hid_t layerDataGroupId = QH5Utilities::createGroup(fileId, "Layer Data");
        hid_t layerGroupId = QH5Utilities::createGroup(layerDataGroupId, QString::number(layerIndex));
        hid_t hfGroupId = QH5Utilities::createGroup(layerGroupId, "High Frequency Data");
        hid_t lfGroupId = QH5Utilities::createGroup(layerGroupId, "Low Frequency Data");
        hid_t unknownGroupId = QH5Utilities::createGroup(layerGroupId, "Unassociated Data");
        for(auto&& splitArray : splitArraysToWrite[i])
        {
          std::vector<size_t> tDims(1, splitArray->getNumberOfTuples());
          splitArray->writeH5Data(hfGroupId, tDims);
        }
        QH5Utilities::closeHDF5Object(layerDataGroupId);
        QH5Utilities::closeHDF5Object(layerGroupId);
        QH5Utilities::closeHDF5Object(hfGroupId);
        QH5Utilities::closeHDF5Object(lfGroupId);
        QH5Utilities::closeHDF5Object(unknownGroupId);
      }
    }

    layerIndex++;
    counter++;
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ImportPrintRiteTDMSFiles::determinePointsForLeastSquares(const PRH::Polygons& polygons, FloatArrayType::Pointer tdms, std::pair<Int32ArrayType::Pointer, std::vector<bool>> pointsToPolys,
                                                              BoolArrayType::Pointer mask)
{
  PRH::Polygons convexSTLPolys;
  PRH::Polygons convexTDMSPolys;
  size_t counter = 1;
  size_t numValidPolys = 0;
  for(auto&& poly : pointsToPolys.second)
  {
    if(poly)
    {
      numValidPolys++;
    }
  }

  for(auto p = 0; p < polygons.polygons.size(); p++)
  {
    if(pointsToPolys.second[p])
    {
      TriMesh::VertexCoordList vertexList(polygons.polygons[p].vertices.size(), std::vector<float>(3));
      for(auto v = 0; v < vertexList.size(); v++)
      {
        vertexList[v][0] = polygons.polygons[p].vertices[v].x;
        vertexList[v][1] = polygons.polygons[p].vertices[v].y;
        vertexList[v][2] = 0.0f;
      }
      Delaunay2D::Pointer delaunay = Delaunay2D::New(vertexList, 1.0, 1e-05, 0.0);
      QString ss = QObject::tr("Finding STL Region Convex Hulls || Region %1 of %2").arg(counter).arg(numValidPolys);
      notifyStatusMessage(ss);
      TriangleGeom::Pointer triangles = delaunay->triangulate();
      triangles->findUnsharedEdges();
      SharedEdgeList::Pointer boundaries = triangles->getUnsharedEdges();
      size_t* edge = boundaries->getPointer(0);
      float* verts = triangles->getVertexPointer(0);
      size_t numBoundaryEdges = boundaries->getNumberOfTuples();

      int64_t startVertex = edge[0];
      PRH::Polygon woundVertices;
      woundVertices.vertices.emplace_back(verts[3 * startVertex + 0], verts[3 * startVertex + 1]);
      std::unordered_set<size_t> visitedEdges;
      visitedEdges.insert(0);
      while(visitedEdges.size() < numBoundaryEdges)
      {
        for(size_t i = 1; i < numBoundaryEdges; i++)
        {
          if(visitedEdges.count(i) == 0)
          {
            if(edge[2 * i] == startVertex)
            {
              woundVertices.vertices.emplace_back(verts[3 * edge[2 * i + 1] + 0], verts[3 * edge[2 * i + 1] + 1]);
              startVertex = edge[2 * i + 1];
              visitedEdges.insert(i);
              break;
            }
            if(edge[2 * i + 1] == startVertex)
            {
              woundVertices.vertices.emplace_back(verts[3 * edge[2 * i] + 0], verts[3 * edge[2 * i] + 1]);
              startVertex = edge[2 * i];
              visitedEdges.insert(i);
              break;
            }
          }
        }
      }

      convexSTLPolys.polygons.push_back(woundVertices);
      counter++;
    }
  }

  float* tdmsPtr = tdms->getPointer(0);
  int32_t* pointsToPolysPtr = pointsToPolys.first->getPointer(0);
  bool* maskPtr = mask->getPointer(0);
  std::vector<TriMesh::VertexCoordList> vertexLists(polygons.polygons.size());
  counter = 1;
  for(size_t i = 0; i < tdms->getNumberOfTuples(); i++)
  {
    if(maskPtr[i] && pointsToPolysPtr[i] >= 0)
    {
      std::vector<float> vert = {tdmsPtr[2 * i + 0], tdmsPtr[2 * i + 1], 0.0f};
      vertexLists[pointsToPolysPtr[i]].push_back(vert);
    }
  }

  for(auto p = 0; p < vertexLists.size(); p++)
  {
    if(pointsToPolys.second[p])
    {
      Delaunay2D::Pointer delaunay = Delaunay2D::New(vertexLists[p], 1.0, 1e-05, 0.0, this);
      delaunay->setMessagePrefix(getHumanLabel());
      QString ss = QObject::tr("Finding TDMS Region Convex Hulls || Region %1 of %2").arg(counter).arg(numValidPolys);
      delaunay->setMessageTitle(ss);
      TriangleGeom::Pointer triangles = delaunay->triangulate();
      triangles->findUnsharedEdges();
      SharedEdgeList::Pointer boundaries = triangles->getUnsharedEdges();
      size_t* edge = boundaries->getPointer(0);
      float* verts = triangles->getVertexPointer(0);
      size_t numBoundaryEdges = boundaries->getNumberOfTuples();

      size_t startVertex = edge[0];
      PRH::Polygon woundVertices;
      woundVertices.vertices.emplace_back(verts[3 * startVertex + 0], verts[3 * startVertex + 1]);
      std::unordered_set<size_t> visitedEdges;
      visitedEdges.insert(0);
      while(visitedEdges.size() < numBoundaryEdges)
      {
        for(size_t i = 1; i < numBoundaryEdges; i++)
        {
          if(visitedEdges.count(i) == 0)
          {
            if(edge[2 * i] == startVertex)
            {
              woundVertices.vertices.emplace_back(verts[3 * edge[2 * i + 1] + 0], verts[3 * edge[2 * i + 1] + 1]);
              startVertex = edge[2 * i + 1];
              visitedEdges.insert(i);
              break;
            }
            if(edge[2 * i + 1] == startVertex)
            {
              woundVertices.vertices.emplace_back(verts[3 * edge[2 * i] + 0], verts[3 * edge[2 * i] + 1]);
              startVertex = edge[2 * i];
              visitedEdges.insert(i);
              break;
            }
          }
        }
      }

      convexTDMSPolys.polygons.push_back(woundVertices);
      counter++;
    }
  }

  std::vector<std::vector<float>> keyPointsSTL = findPolygonKeyPoints(convexSTLPolys);
  std::vector<std::vector<float>> keyPointsTDMS = findPolygonKeyPoints(convexTDMSPolys);

  m_Polynomial.setDependentVariables(keyPointsSTL);
  m_Polynomial.setIndependentVariables(keyPointsTDMS);
  m_Polynomial.setOrder(3);
  m_Polynomial.performRegression();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
std::vector<std::vector<float>> ImportPrintRiteTDMSFiles::findPolygonKeyPoints(const PRH::Polygons& polygons)
{
  std::vector<std::vector<float>> keyPoints(2, std::vector<float>(polygons.polygons.size() * 5));
  for(auto p = 0; p < polygons.polygons.size(); p++)
  {
    std::vector<std::vector<float>> keys = polygons.polygons[p].key_points();
    std::copy(std::begin(keys[0]), std::end(keys[0]), std::begin(keyPoints[0]) + 5 * p);
    std::copy(std::begin(keys[1]), std::end(keys[1]), std::begin(keyPoints[1]) + 5 * p);
  }
  return keyPoints;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
PRH::Polygons ImportPrintRiteTDMSFiles::extractLayerPolygons(int32_t layer)
{
  Int32ArrayType::Pointer regionIds = m_LocalStructure->getDataContainer(SIMPL::Defaults::EdgeDataContainerName)
                                          ->getAttributeMatrix(SIMPL::Defaults::EdgeAttributeMatrixName)
                                          ->getAttributeArrayAs<Int32ArrayType>("_INTERNAL_USE_ONLY_RegionIds");
  Int32ArrayType::Pointer sliceIds = m_LocalStructure->getDataContainer(SIMPL::Defaults::EdgeDataContainerName)
                                         ->getAttributeMatrix(SIMPL::Defaults::EdgeAttributeMatrixName)
                                         ->getAttributeArrayAs<Int32ArrayType>("_INTERNAL_USE_ONLY_SliceIds");
  EdgeGeom::Pointer edges = m_LocalStructure->getDataContainer(SIMPL::Defaults::EdgeDataContainerName)->getGeometryAs<EdgeGeom>();
  int32_t* sliceIdsPtr = sliceIds->getPointer(0);
  int32_t* regionIdsPtr = regionIds->getPointer(0);
  float* eVerts = edges->getVertexPointer(0);
  size_t* eEdges = edges->getEdgePointer(0);

  std::vector<std::vector<int64_t>> polyEdgeList(m_NumParts);
  PRH::Polygons woundPolygons;
  woundPolygons.polygons.resize(m_NumParts);

  for(size_t e = 0; e < edges->getNumberOfEdges(); e++)
  {
    if(sliceIdsPtr[e] == layer)
    {
      polyEdgeList[regionIdsPtr[e] - 1].push_back(e);
    }
  }

  for(auto p = 0; p < polyEdgeList.size(); p++)
  {
    if(!polyEdgeList[p].empty())
    {
      std::vector<bool> visited(polyEdgeList[p].size(), false);
      visited[0] = true;
      PRH::Vertex vert0(eVerts[3 * eEdges[2 * polyEdgeList[p][0] + 0] + 0], eVerts[3 * eEdges[2 * polyEdgeList[p][0] + 0] + 1]);
      PRH::Vertex vert1(eVerts[3 * eEdges[2 * polyEdgeList[p][0] + 1] + 0], eVerts[3 * eEdges[2 * polyEdgeList[p][0] + 1] + 1]);
      woundPolygons.polygons[p].vertices.push_back(vert0);
      woundPolygons.polygons[p].vertices.push_back(vert1);
      for(auto e = 0; e < polyEdgeList[p].size(); e++)
      {
        if(!visited[e])
        {
          PRH::Vertex vert2(eVerts[3 * eEdges[2 * polyEdgeList[p][e] + 0] + 0], eVerts[3 * eEdges[2 * polyEdgeList[p][e] + 0] + 1]);
          if(vert1 == vert2)
          {
            PRH::Vertex vert3(eVerts[3 * eEdges[2 * polyEdgeList[p][e] + 1] + 0], eVerts[3 * eEdges[2 * polyEdgeList[p][e] + 1] + 1]);
            woundPolygons.polygons[p].vertices.push_back(vert3);
            vert1 = vert3;
            visited[e] = true;
            e = 0;
          }
        }
      }
    }
  }

  for(auto&& poly : woundPolygons.polygons)
  {
    if(poly.valid())
    {
      poly.vertices.pop_back();
    }
  }

  return woundPolygons;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
std::pair<Int32ArrayType::Pointer, std::vector<bool>> ImportPrintRiteTDMSFiles::associatePointsWithPolygons(const PRH::Polygons& polygons, FloatArrayType::Pointer tdms)
{
  Int32ArrayType::Pointer pointsToPolygons = Int32ArrayType::CreateArray(tdms->getNumberOfTuples(), std::string("_INTERNAL_USE_ONLY_PointsToPolygons"), true);
  pointsToPolygons->initializeWithValue(-1);
  int32_t* pointsToPolygonsPtr = pointsToPolygons->getPointer(0);
  float* tdmsPtr = tdms->getPointer(0);
  std::vector<PRH::BoundingBox> boundingBoxes(polygons.polygons.size());

  for(auto p = 0; p < polygons.polygons.size(); p++)
  {
    if(polygons.exists(p))
    {
      PRH::BoundingBox bbox = polygons.polygons[p].bounding_box();
      bbox.pad(1);
      boundingBoxes[p] = bbox;
    }
  }

  std::vector<size_t> polyPointCounts(polygons.polygons.size(), 0);
  for(size_t i = 0; i < tdms->getNumberOfTuples(); i++)
  {
    for(auto b = 0; b < boundingBoxes.size(); b++)
    {
      if(polygons.exists(b))
      {
        if(boundingBoxes[b].inside(tdmsPtr[2 * i + 0], tdmsPtr[2 * i + 1]))
        {
          polyPointCounts[b]++;
          pointsToPolygonsPtr[i] = b;
          break;
        }
      }
    }
  }

  std::vector<bool> validPolygons(polygons.polygons.size(), true);
  for(auto i = 0; i < validPolygons.size(); i++)
  {
    if(polyPointCounts[i] < 1000) // TODO: what should this be really?
    {
      validPolygons[i] = false;
    }
  }

  return std::make_pair(pointsToPolygons, validPolygons);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
std::pair<std::vector<float>, PRH::BoundingBox> ImportPrintRiteTDMSFiles::findCentralCluster(FloatArrayType::Pointer tdms, Int32ArrayType::Pointer clusterIds)
{
  float* tdmsPtr = tdms->getPointer(0);
  float xMaxTDMS = std::numeric_limits<float>::lowest();
  float xMinTDMS = std::numeric_limits<float>::max();
  float yMaxTDMS = std::numeric_limits<float>::lowest();
  float yMinTDMS = std::numeric_limits<float>::max();
  for(size_t i = 0; i < tdms->getNumberOfTuples(); i++)
  {
    if(tdmsPtr[2 * i + 0] > xMaxTDMS)
    {
      xMaxTDMS = tdmsPtr[2 * i + 0];
    }
    if(tdmsPtr[2 * i + 0] < xMinTDMS)
    {
      xMinTDMS = tdmsPtr[2 * i + 0];
    }
    if(tdmsPtr[2 * i + 1] > yMaxTDMS)
    {
      yMaxTDMS = tdmsPtr[2 * i + 1];
    }
    if(tdmsPtr[2 * i + 1] < yMinTDMS)
    {
      yMinTDMS = tdmsPtr[2 * i + 1];
    }
  }

  // float extents[2] = {xMaxTDMS - xMinTDMS, yMaxTDMS - yMinTDMS};
  float center[2] = {((xMaxTDMS - xMinTDMS) / 2) + xMinTDMS, ((yMaxTDMS - yMinTDMS) / 2) + yMinTDMS};

  int32_t* clusterIdsPtr = clusterIds->getPointer(0);
  int32_t numClusters = 0;
  std::unordered_map<int32_t, uint32_t> clusterCounts;
  for(size_t i = 0; i < clusterIds->getNumberOfTuples(); i++)
  {
    if(clusterCounts.find(clusterIdsPtr[i]) == std::end(clusterCounts))
    {
      clusterCounts[clusterIdsPtr[i]] = 1;
    }
    else
    {
      clusterCounts[clusterIdsPtr[i]]++;
    }
  }
  for(auto&& pair : clusterCounts)
  {
    if(pair.first > numClusters)
    {
      numClusters = pair.first;
    }
  }

  std::vector<std::vector<float>> clusterCenters(numClusters + 1, std::vector<float>(2, 0.0f));
  for(size_t i = 0; i < tdms->getNumberOfTuples(); i++)
  {
    clusterCenters[clusterIdsPtr[i]][0] += tdmsPtr[2 * i + 0];
    clusterCenters[clusterIdsPtr[i]][1] += tdmsPtr[2 * i + 1];
  }
  for(auto&& pair : clusterCounts)
  {
    clusterCenters[pair.first][0] /= static_cast<float>(pair.second);
    clusterCenters[pair.first][1] /= static_cast<float>(pair.second);
  }

  int32_t mostCentralCluster = 0;
  float distance = std::numeric_limits<float>::max();
  for(auto c = 1; c < clusterCenters.size(); c++)
  {
    float d = std::sqrt((clusterCenters[c][0] - center[0]) * (clusterCenters[c][0] - center[0]) + (clusterCenters[c][1] - center[1]) * (clusterCenters[c][1] - center[1]));
    if(d < distance)
    {
      distance = d;
      mostCentralCluster = c;
    }
  }

  PRH::BoundingBox bbox;
  for(size_t i = 0; i < tdms->getNumberOfTuples(); i++)
  {
    if(clusterIdsPtr[i] == mostCentralCluster)
    {
      if(tdmsPtr[2 * i + 0] < bbox.xmin)
      {
        bbox.xmin = tdmsPtr[2 * i + 0];
      }
      if(tdmsPtr[2 * i + 0] > bbox.xmax)
      {
        bbox.xmax = tdmsPtr[2 * i + 0];
      }
      if(tdmsPtr[2 * i + 1] < bbox.ymin)
      {
        bbox.ymin = tdmsPtr[2 * i + 1];
      }
      if(tdmsPtr[2 * i + 1] > bbox.ymax)
      {
        bbox.ymax = tdmsPtr[2 * i + 1];
      }
    }
  }

  return std::make_pair(clusterCenters[mostCentralCluster], bbox);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
std::pair<std::vector<float>, PRH::BoundingBox> ImportPrintRiteTDMSFiles::findCentralPolygon(const PRH::Polygons& polygons, int32_t layer, const QString& fname)
{
  Int32ArrayType::Pointer sliceIds = m_LocalStructure->getDataContainer(SIMPL::Defaults::EdgeDataContainerName)
                                         ->getAttributeMatrix(SIMPL::Defaults::EdgeAttributeMatrixName)
                                         ->getAttributeArrayAs<Int32ArrayType>("_INTERNAL_USE_ONLY_SliceIds");
  EdgeGeom::Pointer edges = m_LocalStructure->getDataContainer(SIMPL::Defaults::EdgeDataContainerName)->getGeometryAs<EdgeGeom>();
  if(sliceIds->getNumberOfTuples() != edges->getNumberOfEdges())
  {
    QString ss = QObject::tr("The length of the internal layer Ids array does not match the number of edges in the sliced build geometry\n"
                             "Layer Ids Length: %1\n"
                             "Number of Edges: %2\n"
                             "Layer Index: %3\n"
                             "File: %3")
                     .arg(sliceIds->getNumberOfTuples())
                     .arg(edges->getNumberOfEdges())
                     .arg(layer)
                     .arg(fname);
    setErrorCondition(-1, ss);
    // return; TODO: fix return
  }

  int32_t* sliceIdsPtr = sliceIds->getPointer(0);
  float* sliceVerts = edges->getVertexPointer(0);
  size_t* sliceEdges = edges->getEdgePointer(0);
  float xMaxSlice = std::numeric_limits<float>::lowest();
  float xMinSlice = std::numeric_limits<float>::max();
  float yMaxSlice = std::numeric_limits<float>::lowest();
  float yMinSlice = std::numeric_limits<float>::max();
  for(size_t e = 0; e < edges->getNumberOfEdges(); e++)
  {
    if(sliceIdsPtr[e] == layer)
    {
      if(sliceVerts[3 * sliceEdges[2 * e + 0] + 0] > xMaxSlice)
      {
        xMaxSlice = sliceVerts[3 * sliceEdges[2 * e + 0] + 0];
      }
      if(sliceVerts[3 * sliceEdges[2 * e + 1] + 0] > xMaxSlice)
      {
        xMaxSlice = sliceVerts[3 * sliceEdges[2 * e + 1] + 0];
      }
      if(sliceVerts[3 * sliceEdges[2 * e + 0] + 0] < xMinSlice)
      {
        xMinSlice = sliceVerts[3 * sliceEdges[2 * e + 0] + 0];
      }
      if(sliceVerts[3 * sliceEdges[2 * e + 1] + 0] < xMinSlice)
      {
        xMinSlice = sliceVerts[3 * sliceEdges[2 * e + 1] + 0];
      }
      if(sliceVerts[3 * sliceEdges[2 * e + 0] + 1] > yMaxSlice)
      {
        yMaxSlice = sliceVerts[3 * sliceEdges[2 * e + 0] + 1];
      }
      if(sliceVerts[3 * sliceEdges[2 * e + 1] + 1] > yMaxSlice)
      {
        yMaxSlice = sliceVerts[3 * sliceEdges[2 * e + 1] + 1];
      }
      if(sliceVerts[3 * sliceEdges[2 * e + 0] + 1] < yMinSlice)
      {
        yMinSlice = sliceVerts[3 * sliceEdges[2 * e + 0] + 1];
      }
      if(sliceVerts[3 * sliceEdges[2 * e + 1] + 1] < yMinSlice)
      {
        yMinSlice = sliceVerts[3 * sliceEdges[2 * e + 1] + 1];
      }
    }
  }

  // float extents[2] = {xMaxSlice - xMinSlice, yMaxSlice - yMinSlice};
  float center[2] = {((xMaxSlice - xMinSlice) / 2) + xMinSlice, ((yMaxSlice - yMinSlice) / 2) + yMinSlice};

  std::vector<std::vector<float>> polyCenters(polygons.polygons.size());
  for(auto p = 0; p < polygons.polygons.size(); p++)
  {
    if(polygons.exists(p))
    {
      polyCenters[p] = polygons.polygons[p].centroid();
    }
  }

  int32_t mostCentralPolygon = 0;
  float distance = std::numeric_limits<float>::max();
  for(auto p = 0; p < polyCenters.size(); p++)
  {
    if(polygons.exists(p))
    {
      float d = std::sqrt((polyCenters[p][0] - center[0]) * (polyCenters[p][0] - center[0]) + (polyCenters[p][1] - center[1]) * (polyCenters[p][1] - center[1]));
      if(d < distance)
      {
        distance = d;
        mostCentralPolygon = p;
      }
    }
  }

  PRH::BoundingBox bbox = polygons.polygons[mostCentralPolygon].bounding_box();

  return std::make_pair(polyCenters[mostCentralPolygon], bbox);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ImportPrintRiteTDMSFiles::execute()
{
  initialize();
  dataCheck();
  if(getErrorCode() < 0)
  {
    return;
  }

  bool hasMissingFiles = false;
  bool orderAscending = false;

  if(getInputFilesList().Ordering == 0)
  {
    orderAscending = true;
  }
  else if(getInputFilesList().Ordering == 1)
  {
    orderAscending = false;
  }

  QVector<QString> fileList =
      FilePathGenerator::GenerateFileList(m_InputFilesList.StartIndex, m_InputFilesList.EndIndex, m_InputFilesList.IncrementIndex, hasMissingFiles, orderAscending, m_InputFilesList.InputPath,
                                          m_InputFilesList.FilePrefix, m_InputFilesList.FileSuffix, m_InputFilesList.FileExtension, m_InputFilesList.PaddingDigits);

  QString ss;

  if(fileList.empty())
  {
    ss.clear();
    QTextStream out(&ss);
    out << " No files have been selected for import. Have you set the input directory and other values so that input files will be generated?\n";
    out << "InputPath: " << m_InputFilesList.InputPath << "\n";
    out << "FilePrefix: " << m_InputFilesList.FilePrefix << "\n";
    out << "FileSuffix: " << m_InputFilesList.FileSuffix << "\n";
    out << "FileExtension: " << m_InputFilesList.FileExtension << "\n";
    out << "PaddingDigits: " << m_InputFilesList.PaddingDigits << "\n";
    out << "StartIndex: " << m_InputFilesList.StartIndex << "\n";
    out << "EndIndex: " << m_InputFilesList.EndIndex << "\n";

    setErrorCondition(-11, ss);
    return;
  }

  QVector<QString> zeroIndexFile = FilePathGenerator::GenerateFileList(0, 0, m_InputFilesList.IncrementIndex, hasMissingFiles, orderAscending, m_InputFilesList.InputPath, m_InputFilesList.FilePrefix,
                                                                       m_InputFilesList.FileSuffix, m_InputFilesList.FileExtension, m_InputFilesList.PaddingDigits);

  if(zeroIndexFile.size() != 1)
  {
    QString ss = QObject::tr("Unable to construct unique zero-indexed file name");
    setErrorCondition(-1, ss);
    return;
  }

  QFileInfo zeroIndexFname(zeroIndexFile[0]);

  if(m_InputFilesList.StartIndex == 0 || zeroIndexFname.exists())
  {
    m_Offset = 1;
  }

  m_NumLayersToImport = m_InputFilesList.EndIndex - m_InputFilesList.StartIndex + 1;

  ss = QObject::tr("Will Import PrintRite Data for Layers %1 to %2 (%3 Total Layer(s))").arg(m_InputFilesList.StartIndex).arg(m_InputFilesList.EndIndex).arg(m_NumLayersToImport);
  notifyStatusMessage(ss);

  // TODO: for all functions, return out on errors to bail; also add getCancel() checks

  if(m_SpatialTransformOption == 1 || m_SpatialTransformOption == 2)
  {
    m_STLFilePath = m_SpatialTransformOption == 1 ? m_STLFilePath1 : m_STLFilePath2;
    importLabelSliceSTL();

    if(getErrorCode() < 0)
    {
      return;
    }
  }

  if(m_SpatialTransformOption == 1)
  {
    QString fileForScaling;
    for(auto i = 0; i < fileList.size(); i++)
    {
      QFileInfo fi(fileList[i]);
      QString baseName = fi.baseName();
      baseName.remove(getInputFilesList().FileSuffix);
      baseName.remove(getInputFilesList().FilePrefix);
      if(!baseName.isEmpty())
      {
        if(getInputFilesList().PaddingDigits > 1)
        {
          bool initialZero = (baseName.at(0) == "0");
          if(!initialZero)
          {
            if(baseName.toUInt() == m_LayerForScaling)
            {
              fileForScaling = fileList[i];
              break;
            }
          }
          else
          {
            while(baseName.at(0) == "0")
            {
              baseName.remove(0, 1);
            }
            if(baseName.isEmpty() && m_LayerForScaling == 0)
            {
              fileForScaling = fileList[i];
              break;
            }
            if(baseName.toUInt() == m_LayerForScaling)
            {
              fileForScaling = fileList[i];
              break;
            }
          }
        }
        else
        {
          if(baseName.toUInt() == m_LayerForScaling)
          {
            fileForScaling = fileList[i];
            break;
          }
        }
      }
    }

    if(fileForScaling.isEmpty())
    {
      QString ss = QObject::tr("Unable to determine layer file to compute spatial transformation");
      setErrorCondition(-1, ss);
      return;
    }

    computeSpatialTransformation(fileForScaling);

    if(getErrorCode() < 0)
    {
      return;
    }
  }

  createHDF5Files();

  if(getErrorCode() < 0)
  {
    return;
  }

  processLayers(fileList);

  if(getErrorCode() < 0)
  {
    return;
  }

  writeMetaDataToHDF5Files();

  if(getErrorCode() < 0)
  {
    return;
  }

  closeHDF5Files();

  if(getErrorCode() < 0)
  {
    return;
  }

  if(m_SpatialTransformOption == 1)
  {
    writeSpatialTransformToFile();
  }

  if(getErrorCode() < 0)
  {
    return;
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AbstractFilter::Pointer ImportPrintRiteTDMSFiles::newFilterInstance(bool copyFilterParameters) const
{
  ImportPrintRiteTDMSFiles::Pointer filter = ImportPrintRiteTDMSFiles::New();
  if(copyFilterParameters)
  {
    copyFilterParameterInstanceVariables(filter.get());
  }
  return filter;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ImportPrintRiteTDMSFiles::getCompiledLibraryName() const
{
  return DREAM3DReviewConstants::DREAM3DReviewBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ImportPrintRiteTDMSFiles::getBrandingString() const
{
  return "DREAM3DReview";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ImportPrintRiteTDMSFiles::getFilterVersion() const
{
  QString version;
  QTextStream vStream(&version);
  vStream << DREAM3DReview::Version::Major() << "." << DREAM3DReview::Version::Minor() << "." << DREAM3DReview::Version::Patch();
  return version;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ImportPrintRiteTDMSFiles::getGroupName() const
{
  return SIMPL::FilterGroups::IOFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ImportPrintRiteTDMSFiles::getSubGroupName() const
{
  return SIMPL::FilterSubGroups::InputFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ImportPrintRiteTDMSFiles::getHumanLabel() const
{
  return "Import PrintRite TDMS File(s) to HDF5";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QUuid ImportPrintRiteTDMSFiles::getUuid() const
{
  return QUuid("{38fde424-0292-5c42-b3b4-18d80c95524d}");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
ImportPrintRiteTDMSFiles::Pointer ImportPrintRiteTDMSFiles::NullPointer()
{
  return Pointer(static_cast<Self*>(nullptr));
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
std::shared_ptr<ImportPrintRiteTDMSFiles> ImportPrintRiteTDMSFiles::New()
{
  struct make_shared_enabler : public ImportPrintRiteTDMSFiles
  {
  };
  std::shared_ptr<make_shared_enabler> val = std::make_shared<make_shared_enabler>();
  val->setupFilterParameters();
  return val;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ImportPrintRiteTDMSFiles::getNameOfClass() const
{
  return QString("ImportPrintRiteTDMSFiles");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ImportPrintRiteTDMSFiles::ClassName()
{
  return QString("ImportPrintRiteTDMSFiles");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ImportPrintRiteTDMSFiles::setInputFilesList(const StackFileListInfo& value)
{
  m_InputFilesList = value;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
StackFileListInfo ImportPrintRiteTDMSFiles::getInputFilesList() const
{
  return m_InputFilesList;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ImportPrintRiteTDMSFiles::setSTLFilePath1(const QString& value)
{
  m_STLFilePath1 = value;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ImportPrintRiteTDMSFiles::getSTLFilePath1() const
{
  return m_STLFilePath1;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ImportPrintRiteTDMSFiles::setSTLFilePath2(const QString& value)
{
  m_STLFilePath2 = value;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ImportPrintRiteTDMSFiles::getSTLFilePath2() const
{
  return m_STLFilePath2;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ImportPrintRiteTDMSFiles::setLayerThickness(const float& value)
{
  m_LayerThickness = value;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
float ImportPrintRiteTDMSFiles::getLayerThickness() const
{
  return m_LayerThickness;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ImportPrintRiteTDMSFiles::setLaserOnThreshold(const float& value)
{
  m_LaserOnThreshold = value;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
float ImportPrintRiteTDMSFiles::getLaserOnThreshold() const
{
  return m_LaserOnThreshold;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ImportPrintRiteTDMSFiles::setLaserOnArrayOption(const int& value)
{
  m_LaserOnArrayOption = value;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int ImportPrintRiteTDMSFiles::getLaserOnArrayOption() const
{
  return m_LaserOnArrayOption;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ImportPrintRiteTDMSFiles::setOutputDirectory(const QString& value)
{
  m_OutputDirectory = value;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ImportPrintRiteTDMSFiles::getOutputDirectory() const
{
  return m_OutputDirectory;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ImportPrintRiteTDMSFiles::setOutputFilePrefix(const QString& value)
{
  m_OutputFilePrefix = value;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ImportPrintRiteTDMSFiles::getOutputFilePrefix() const
{
  return m_OutputFilePrefix;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ImportPrintRiteTDMSFiles::setSpatialTransformOption(const int& value)
{
  m_SpatialTransformOption = value;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int ImportPrintRiteTDMSFiles::getSpatialTransformOption() const
{
  return m_SpatialTransformOption;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ImportPrintRiteTDMSFiles::setDowncastRawData(const bool& value)
{
  m_DowncastRawData = value;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
bool ImportPrintRiteTDMSFiles::getDowncastRawData() const
{
  return m_DowncastRawData;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ImportPrintRiteTDMSFiles::setScaleLaserPower(const bool& value)
{
  m_ScaleLaserPower = value;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
bool ImportPrintRiteTDMSFiles::getScaleLaserPower() const
{
  return m_ScaleLaserPower;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ImportPrintRiteTDMSFiles::setPowerScalingCoefficients(const FloatVec2Type& value)
{
  m_PowerScalingCoefficients = value;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
FloatVec2Type ImportPrintRiteTDMSFiles::getPowerScalingCoefficients() const
{
  return m_PowerScalingCoefficients;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ImportPrintRiteTDMSFiles::setScalePyrometerTemperature(const bool& value)
{
  m_ScalePyrometerTemperature = value;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
bool ImportPrintRiteTDMSFiles::getScalePyrometerTemperature() const
{
  return m_ScalePyrometerTemperature;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ImportPrintRiteTDMSFiles::setTemperatureScalingCoefficients(const FloatVec2Type& value)
{
  m_TemperatureScalingCoefficients = value;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
FloatVec2Type ImportPrintRiteTDMSFiles::getTemperatureScalingCoefficients() const
{
  return m_TemperatureScalingCoefficients;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ImportPrintRiteTDMSFiles::setSplitRegions1(const bool& value)
{
  m_SplitRegions1 = value;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
bool ImportPrintRiteTDMSFiles::getSplitRegions1() const
{
  return m_SplitRegions1;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ImportPrintRiteTDMSFiles::setSplitRegions2(const bool& value)
{
  m_SplitRegions2 = value;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
bool ImportPrintRiteTDMSFiles::getSplitRegions2() const
{
  return m_SplitRegions2;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ImportPrintRiteTDMSFiles::setLayerForScaling(const int& value)
{
  m_LayerForScaling = value;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int ImportPrintRiteTDMSFiles::getLayerForScaling() const
{
  return m_LayerForScaling;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ImportPrintRiteTDMSFiles::setInputSpatialTransformFilePath(const QString& value)
{
  m_InputSpatialTransformFilePath = value;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ImportPrintRiteTDMSFiles::getInputSpatialTransformFilePath() const
{
  return m_InputSpatialTransformFilePath;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ImportPrintRiteTDMSFiles::setSearchRadius(const float& value)
{
  m_SearchRadius = value;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
float ImportPrintRiteTDMSFiles::getSearchRadius() const
{
  return m_SearchRadius;
}
