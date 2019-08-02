/*
 * Your License or Copyright can go here
 */
#include "ITKResampleImage.h"

#include <map>

#include "itkAffineTransform.h"
#include "itkBSplineTransform.h"
#include "itkBSplineTransformInitializer.h"
#include "itkCastImageFilter.h"
#include "itkCenteredTransformInitializer.h"
#include "itkEuler2DTransform.h"
#include "itkImageFileWriter.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkResampleImageFilter.h"

#include <QtCore/QFile>
#include <QtCore/QFileInfo>

#include "H5Support/H5ScopedSentinel.h"
#include "H5Support/QH5Utilities.h"

#include "SIMPLib/Common/Constants.h"
#include "SIMPLib/Common/TemplateHelpers.h"
#include "SIMPLib/FilterParameters/ChoiceFilterParameter.h"
#include "SIMPLib/FilterParameters/DataArraySelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/DataContainerCreationFilterParameter.h"
#include "SIMPLib/FilterParameters/FileListInfoFilterParameter.h"
#include "SIMPLib/FilterParameters/InputFileFilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedChoicesFilterParameter.h"
//#include "SIMPLib/FilterParameters/LinkedPathCreationFilterParameter.h"
//#include "SIMPLib/FilterParameters/AttributeMatrixCreationFilterParameter.h"
//#include "SIMPLib/FilterParamters/DataArrayCreationFilterParameter.h"
#include "SIMPLib/FilterParameters/OutputPathFilterParameter.h"
#include "SIMPLib/FilterParameters/SeparatorFilterParameter.h"
#include "SIMPLib/FilterParameters/StringFilterParameter.h"
#include "SIMPLib/Filtering/FilterManager.h"
#include "SIMPLib/Geometry/TransformContainer.h"
#include "SIMPLib/ITK/itkInPlaceDream3DDataToImageFilter.h"
#include "SIMPLib/ITK/itkInPlaceImageToDream3DDataFilter.h"
#include "SIMPLib/Utilities/FilePathGenerator.h"
#include "SIMPLib/Utilities/FileSystemPathHelper.h"

#include "EbsdLib/EbsdConstants.h"
#include "EbsdLib/EbsdLib.h"
#include "EbsdLib/HKL/CtfReader.h"

#include "DREAM3DReview/DREAM3DReviewConstants.h"
#include "DREAM3DReview/DREAM3DReviewFilters/util/EBSDWriterFactory.h"
#include "DREAM3DReview/DREAM3DReviewFilters/util/ITKTransformHelpers.h"
#include "DREAM3DReview/DREAM3DReviewVersion.h"

// enum createdPathID : RenameDataPath::DataID_t
//{
//  AttributeMatrixID21 = 21,
//
//  DataArrayID31 = 31,
//
//  DataContainerID = 1
//};

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
ITKResampleImage::ITKResampleImage()
: m_OutputPath("")
, m_FileNamePrefix("")
, m_EBSDFileNamePrefix("")
, m_OutputPathEBSD("")
, m_TransformFileName("")
, m_InterpolationType(0)
, m_MovingImageArrayPath("", "", "")
, m_DataContainerName("ResampledDataDC")
, m_CellAttributeMatrixName("ResampledDataAM")
, m_ImageDataArrayName("ResampledData")
, m_OperationMode(0)
{
  m_ctfDataArrayNames.push_back("Phases");
  m_ctfDataArrayNames.push_back("X");
  m_ctfDataArrayNames.push_back("Y");
  m_ctfDataArrayNames.push_back("Bands");
  m_ctfDataArrayNames.push_back("Error");
  m_ctfDataArrayNames.push_back("EulerAngles");
  m_ctfDataArrayNames.push_back("MAD");
  m_ctfDataArrayNames.push_back("BC");
  m_ctfDataArrayNames.push_back("BS");

  m_ImageFileListInfo.FileExtension = QString("tif");
  m_ImageFileListInfo.StartIndex = 0;
  m_ImageFileListInfo.EndIndex = 0;
  m_ImageFileListInfo.PaddingDigits = 0;

  initialize();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
ITKResampleImage::~ITKResampleImage() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ITKResampleImage::initialize()
{
  setErrorCondition(0);
  setWarningCondition(0);
  setCancel(false);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ITKResampleImage::setupFilterParameters()
{
  FilterParameterVector parameters;

  /// MODE TYPE: Either a single pair or series
  {
    LinkedChoicesFilterParameter::Pointer parameter = LinkedChoicesFilterParameter::New();
    parameter->setHumanLabel("Operation Mode");
    parameter->setPropertyName("OperationMode");
    parameter->setSetterCallback(SIMPL_BIND_SETTER(ITKResampleImage, this, OperationMode));
    parameter->setGetterCallback(SIMPL_BIND_GETTER(ITKResampleImage, this, OperationMode));

    QVector<QString> choices;
    choices.push_back("Single Data Array");
    choices.push_back("Series of Images");
    choices.push_back("Series of CTFs");

    parameter->setChoices(choices);

    QStringList linkedProps;
    linkedProps << "ImageDataArrayName"
                << "CellAttributeMatrixName"
                << "DataContainerName"
                << "MovingImageArrayPath"
                << "ImageFileListInfo"
                << "OrientationFileListInfo"
                << "OutputPath"
                << "OutputPathEBSD"
                << "FileNamePrefix"
                << "EBSDFileNamePrefix";
    parameter->setLinkedProperties(linkedProps);
    parameter->setEditable(false);
    parameter->setCategory(FilterParameter::Parameter);
    parameters.push_back(parameter);
  }

  parameters.push_back(SIMPL_NEW_INPUT_FILE_FP("Transform File Name", TransformFileName, FilterParameter::Parameter, ITKResampleImage, "*.hdf5"));
  {
    ChoiceFilterParameter::Pointer parameter = ChoiceFilterParameter::New();
    parameter->setHumanLabel("Interpolation Type");
    parameter->setPropertyName("InterpolationType");

    QVector<QString> choices; // Please add choices to the choices QVector to finish this widget
    choices.push_back("Linear");
    choices.push_back("Nearest Neighbor");
    parameter->setChoices(choices);
    parameter->setChoices(choices);
    parameter->setCategory(FilterParameter::Parameter);
    parameter->setSetterCallback(SIMPL_BIND_SETTER(ITKResampleImage, this, InterpolationType));
    parameter->setGetterCallback(SIMPL_BIND_GETTER(ITKResampleImage, this, InterpolationType));
    parameters.push_back(parameter);
  }

  DataArraySelectionFilterParameter::RequirementType dasReq =
      DataArraySelectionFilterParameter::CreateRequirement(SIMPL::Defaults::AnyPrimitive, SIMPL::Defaults::AnyComponentSize, AttributeMatrix::Type::Cell, IGeometry::Type::Image);
  parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Moving Image", MovingImageArrayPath, FilterParameter::RequiredArray, ITKResampleImage, dasReq, 0));

  parameters.push_back(SIMPL_NEW_STRING_FP("Data Container", DataContainerName, FilterParameter::CreatedArray, ITKResampleImage, 0));
  parameters.push_back(SeparatorFilterParameter::New("Cell Data", FilterParameter::CreatedArray));
  parameters.push_back(SIMPL_NEW_STRING_FP("Cell Attribute Matrix", CellAttributeMatrixName, FilterParameter::CreatedArray, ITKResampleImage, 0));
  parameters.push_back(SIMPL_NEW_STRING_FP("Image Data", ImageDataArrayName, FilterParameter::CreatedArray, ITKResampleImage, 0));
  parameters.push_back(SIMPL_NEW_OUTPUT_PATH_FP("Output Path for Resampled Images", OutputPath, FilterParameter::Parameter, ITKResampleImage, "", "", 1));
  parameters.push_back(SIMPL_NEW_STRING_FP("Output Image File Name Prefix", FileNamePrefix, FilterParameter::Parameter, ITKResampleImage, 1));
  parameters.push_back(SIMPL_NEW_STRING_FP("Output EBSD File Name Prefix", EBSDFileNamePrefix, FilterParameter::Parameter, ITKResampleImage, 2));

  parameters.push_back(SIMPL_NEW_OUTPUT_PATH_FP("Output Path for Resampled Orientation Data", OutputPathEBSD, FilterParameter::Parameter, ITKResampleImage, "", "", 2));

  parameters.push_back(SIMPL_NEW_FILELISTINFO_FP("Image List", ImageFileListInfo, FilterParameter::Parameter, ITKResampleImage));
  FileListInfoFilterParameter::Pointer imageFileList = std::dynamic_pointer_cast<FileListInfoFilterParameter>(parameters.back());
  imageFileList->setGroupIndex(1);

  parameters.push_back(SIMPL_NEW_FILELISTINFO_FP("EBSD File List", OrientationFileListInfo, FilterParameter::Parameter, ITKResampleImage));
  FileListInfoFilterParameter::Pointer orientationFileList = std::dynamic_pointer_cast<FileListInfoFilterParameter>(parameters.back());
  orientationFileList->setGroupIndex(2);

  setFilterParameters(parameters);
}

QVector<QString> ITKResampleImage::getFileList(FileListInfo_t inputFileListInfo)
{
  bool hasMissingFiles = false;
  bool orderAscending = false;

  if(inputFileListInfo.Ordering == 0)
  {
    orderAscending = true;
  }
  else if(inputFileListInfo.Ordering == 1)
  {
    orderAscending = false;
  }

  // Now generate all the file names the user is asking for and populate the table
  return FilePathGenerator::GenerateFileList(inputFileListInfo.StartIndex, inputFileListInfo.EndIndex, inputFileListInfo.IncrementIndex, hasMissingFiles, orderAscending, inputFileListInfo.InputPath,
                                             inputFileListInfo.FilePrefix, inputFileListInfo.FileSuffix, inputFileListInfo.FileExtension, inputFileListInfo.PaddingDigits);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int32_t ITKResampleImage::checkInputFileList(FileListInfo_t inputFileListInfo)
{
  DataArrayPath tempPath;
  QString ss;

  if(inputFileListInfo.InputPath.isEmpty())
  {
    ss = QObject::tr("The moving image input directory must be set");
    setErrorCondition(-64500);
    notifyErrorMessage(getHumanLabel(), ss, getErrorCondition());
  }

  bool orderAscending = false;

  if(inputFileListInfo.Ordering == 0)
  {
    orderAscending = true;
  }
  else if(inputFileListInfo.Ordering == 1)
  {
    orderAscending = false;
  }

  // Now generate all the file names the user is asking for and populate the table
  const QVector<QString> fileList = this->getFileList(inputFileListInfo);
  if(fileList.empty())
  {
    ss.clear();
    QTextStream out(&ss);
    out << " No files have been selected for import. Have you set the input directory and other values so that input files will be generated?\n";
    out << "InputPath: " << inputFileListInfo.InputPath << "\n";
    out << "FilePrefix: " << inputFileListInfo.FilePrefix << "\n";
    out << "FileSuffix: " << inputFileListInfo.FileSuffix << "\n";
    out << "FileExtension: " << inputFileListInfo.FileExtension << "\n";
    out << "PaddingDigits: " << inputFileListInfo.PaddingDigits << "\n";
    out << "StartIndex: " << inputFileListInfo.StartIndex << "\n";
    out << "EndIndex: " << inputFileListInfo.EndIndex << "\n";
    setErrorCondition(-64501);
    notifyErrorMessage(getHumanLabel(), ss, getErrorCondition());

    return -1;
  }

  // Validate all the files in the list. Throw an error for each one if it does not exist
  for(const auto& filePath : fileList)
  {
    QFileInfo fi(filePath);
    if(!fi.exists())
    {
      QString errorMessage = QString("File does not exist: %1").arg(filePath);
      setErrorCondition(-64502);
      notifyErrorMessage(getHumanLabel(), errorMessage, getErrorCondition());
    }
  }
  if(getErrorCondition() < 0)
  {
    return -1;
  }

  // Create a subfilter to read each image, although for preflight we are going to read the first image in the
  // list and hope the rest are correct.
  if(m_OperationMode == 1) ///////FOR IMAGE DATA
  {
    FilterManager* fm = FilterManager::Instance();
    IFilterFactory::Pointer factory = fm->getFactoryFromClassName("ITKImageReader");
    if(factory.get() == nullptr)
    {
      QString ss = QObject::tr("Unable to instantiate Filter with name 'ITKImageReader'\n"
                               "The 'ITKImageReader' Filter is needed to import the image");
      setErrorCondition(-1);
      notifyErrorMessage(getHumanLabel(), ss, getErrorCondition());
    }
    AbstractFilter::Pointer itkImageReader = factory->create();
    DataContainerArray::Pointer dca = DataContainerArray::New();
    itkImageReader->setDataContainerArray(dca);
    QVariant var;
    var.setValue(fileList[0]);
    itkImageReader->setProperty("FileName", var);
    itkImageReader->preflight();
    if(itkImageReader->getErrorCondition() < 0)
    {
      setErrorCondition(itkImageReader->getErrorCondition());
      notifyErrorMessage(getHumanLabel(), "Error Reading Input Image.", getErrorCondition());

      return -1;
    }
  }

  if(m_OperationMode == 2) ////////// FOR ORIENATION DATA
  {
    // Based on the type of file (.ang or .ctf) get the list of arrays that would be created
    QFileInfo fi(fileList.front());
    QString ext = fi.suffix();
    if(ext.compare("ang") == 0)
    {
      // ebsdFeatures = new AngFields;
    }
    else if(ext.compare("ctf") == 0)
    {
      // ebsdFeatures = new CtfFields;
    }
    else
    {
      ss = QObject::tr("The file extension '%1' was not recognized. Currently .ang or .ctf are the only recognized file extensions").arg(ext);
      setErrorCondition(-997);
      notifyErrorMessage(getHumanLabel(), ss, getErrorCondition());

      return -1;
    }
  }

  return fileList.size();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ITKResampleImage::dataCheck()
{
  setErrorCondition(0);
  setWarningCondition(0);

  if(m_OperationMode == 0)
  {
    m_MovingImagePtr = getDataContainerArray()->getPrereqIDataArrayFromPath<IDataArray, AbstractFilter>(this, getMovingImageArrayPath());
    if(getErrorCondition() < 0)
    {
      return;
    }

    DataContainer::Pointer m = getDataContainerArray()->createNonPrereqDataContainer<AbstractFilter>(this, getDataContainerName());
    if(getErrorCondition() < 0)
    {
      return;
    }

    // Create the Image Geometry
    ImageGeom::Pointer image = ImageGeom::CreateGeometry(SIMPL::Geometry::ImageGeometry);
    m->setGeometry(image);

    QVector<size_t> tDims(3, 0);
    AttributeMatrix::Pointer cellAttrMat = m->createNonPrereqAttributeMatrix(this, getCellAttributeMatrixName(), tDims, AttributeMatrix::Type::Cell);
    if(getErrorCondition() < 0)
    {
      return;
    }

    EXECUTE_FUNCTION_TEMPLATE_NO_BOOL(this, createCompatibleArrays, m_MovingImagePtr.lock())
  }
  else if(m_OperationMode == 1)
  {
    hid_t fileId = QH5Utilities::openFile(m_TransformFileName, true);
    QStringList groupObjects;
    QH5Utilities::getGroupObjects(fileId, H5Utilities::H5Support_GROUP, groupObjects);
    int32_t numTransformObjects = groupObjects.size();

    int32_t numImages = checkInputFileList(m_ImageFileListInfo);

    if(numTransformObjects != numImages)
    {
      setErrorCondition(-64505);
      notifyErrorMessage(getHumanLabel(), "The number of images to transform and transform objects in the HDF5 file are required to be the same", getErrorCondition());
    }
  }
  else if(m_OperationMode == 2)
  {
    hid_t fileId = QH5Utilities::openFile(m_TransformFileName, true);
    QStringList groupObjects;
    QH5Utilities::getGroupObjects(fileId, H5Utilities::H5Support_GROUP, groupObjects);
    int32_t numTransformObjects = groupObjects.size();

    int32_t numOrientationFiles = checkInputFileList(m_OrientationFileListInfo);

    if(numTransformObjects != numOrientationFiles)
    {
      setErrorCondition(-64506);
      notifyErrorMessage(getHumanLabel(), "The number of images to transform and transform objects in the HDF5 file are required to be the same", getErrorCondition());
    }
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ITKResampleImage::preflight()
{
  // These are the REQUIRED lines of CODE to make sure the filter behaves correctly
  setInPreflight(true);              // Set the fact that we are preflighting.
  emit preflightAboutToExecute();    // Emit this signal so that other widgets can do one file update
  emit updateFilterParameters(this); // Emit this signal to have the widgets push their values down to the filter
  dataCheck();                       // Run our DataCheck to make sure everthing is setup correctly
  emit preflightExecuted();          // We are done preflighting this filter
  setInPreflight(false);             // Inform the system this filter is NOT in preflight mode anymore.
}
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
template <typename T>
void ITKResampleImage::createCompatibleArrays()
{

  DataArrayPath path(m_DataContainerName, m_CellAttributeMatrixName, m_ImageDataArrayName);

  m_FinalImagePtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<T>, AbstractFilter, T>(this, path, 0, m_MovingImagePtr.lock()->getComponentDimensions());
}

// -----------------------------------------------------------------------------
// BEGIN MACROS
// -----------------------------------------------------------------------------

#define RESAMPLE2D_VARIABLECOMPONENT(sliceNo)                                                                                                                                                          \
  using ToITKType = itk::InPlaceDream3DDataToImageFilter<PixelType, ImageDimension>;                                                                                                                   \
  typename ToITKType::Pointer movingtoITK = ToITKType::New();                                                                                                                                          \
  movingtoITK->SetDataArrayName(m_MovingImageArrayPath.getDataArrayName().toStdString());                                                                                                              \
  movingtoITK->SetAttributeMatrixArrayName(m_MovingImageArrayPath.getAttributeMatrixName().toStdString());                                                                                             \
  DataContainer::Pointer movingDC = getDataContainerArray()->getDataContainer(m_MovingImageArrayPath.getDataContainerName());                                                                          \
  movingtoITK->SetInput(movingDC);                                                                                                                                                                     \
  movingtoITK->InPlaceOn();                                                                                                                                                                            \
  movingtoITK->Update();                                                                                                                                                                               \
  using ImageType = itk::Dream3DImage<PixelType, ImageDimension>;                                                                                                                                      \
  typename ImageType::Pointer itkMovingImage = movingtoITK->GetOutput();                                                                                                                               \
  using ResampleFilterType = itk::ResampleImageFilter<ImageType, ImageType>;                                                                                                                           \
  typename ResampleFilterType::Pointer resample = ResampleFilterType::New();                                                                                                                           \
  ITKTransformHelpers transformhelper = getTransformAndFixedParams(sliceNo);                                                                                                                           \
  typename TransformContainer::Pointer d3dtransform = transformhelper.transform;                                                                                                                       \
  d3dtransform->getFixedParameters();                                                                                                                                                                  \
  QString stringTransformType = QString::fromStdString(d3dtransform->getTransformTypeAsString()).split("_")[0];                                                                                        \
  if(stringTransformType == "BSplineTransform")                                                                                                                                                        \
  {                                                                                                                                                                                                    \
    QString transformOrderString = QString::fromStdString(d3dtransform->getTransformTypeAsString());                                                                                                   \
    if(transformOrderString == "BSplineTransform_double_2_2")                                                                                                                                          \
    {                                                                                                                                                                                                  \
      const unsigned int SplineOrder = 3;                                                                                                                                                              \
      BSPLINESETUP()                                                                                                                                                                                   \
    }                                                                                                                                                                                                  \
    if(transformOrderString == "BSplineTransform_double_2_2_2")                                                                                                                                        \
    {                                                                                                                                                                                                  \
      const unsigned int SplineOrder = 2;                                                                                                                                                              \
      BSPLINESETUP()                                                                                                                                                                                   \
    }                                                                                                                                                                                                  \
    if(transformOrderString == "BSplineTransform_double_2_2_1")                                                                                                                                        \
    {                                                                                                                                                                                                  \
      const unsigned int SplineOrder = 1;                                                                                                                                                              \
      BSPLINESETUP()                                                                                                                                                                                   \
    }                                                                                                                                                                                                  \
  }                                                                                                                                                                                                    \
  if(stringTransformType == "Affine")                                                                                                                                                                  \
  {                                                                                                                                                                                                    \
    using TransformType = itk::AffineTransform<double, ImageDimension>;                                                                                                                                \
    typename TransformType::Pointer transform = TransformType::New();                                                                                                                                  \
    int32_t numFixedElements = d3dtransform->getFixedParameters().size();                                                                                                                              \
    typename TransformType::FixedParametersType fixedParameters(numFixedElements);                                                                                                                     \
    for(int32_t i = 0; i < numFixedElements; i++)                                                                                                                                                      \
    {                                                                                                                                                                                                  \
      fixedParameters.SetElement(i, d3dtransform->getFixedParameters()[i]);                                                                                                                            \
    }                                                                                                                                                                                                  \
    transform->SetFixedParameters(fixedParameters);                                                                                                                                                    \
    int32_t numLearnedElements = d3dtransform->getParameters().size();                                                                                                                                 \
    typename TransformType::ParametersType learnedParameters(numLearnedElements);                                                                                                                      \
    for(int32_t i = 0; i < numLearnedElements; i++)                                                                                                                                                    \
    {                                                                                                                                                                                                  \
      learnedParameters.SetElement(i, d3dtransform->getParameters()[i]);                                                                                                                               \
    }                                                                                                                                                                                                  \
    transform->SetParameters(learnedParameters);                                                                                                                                                       \
    resample->SetTransform(transform);                                                                                                                                                                 \
  }                                                                                                                                                                                                    \
  if(stringTransformType == "Euler2DTransform")                                                                                                                                                        \
  {                                                                                                                                                                                                    \
    using TransformType = itk::Euler2DTransform<double>;                                                                                                                                               \
    typename TransformType::Pointer transform = TransformType::New();                                                                                                                                  \
    int32_t numFixedElements = d3dtransform->getFixedParameters().size();                                                                                                                              \
    typename TransformType::FixedParametersType fixedParameters(numFixedElements);                                                                                                                     \
    for(int32_t i = 0; i < numFixedElements; i++)                                                                                                                                                      \
    {                                                                                                                                                                                                  \
      fixedParameters.SetElement(i, d3dtransform->getFixedParameters()[i]);                                                                                                                            \
    }                                                                                                                                                                                                  \
    transform->SetFixedParameters(fixedParameters);                                                                                                                                                    \
    int32_t numLearnedElements = d3dtransform->getParameters().size();                                                                                                                                 \
    typename TransformType::ParametersType learnedParameters(numLearnedElements);                                                                                                                      \
    for(int32_t i = 0; i < numLearnedElements; i++)                                                                                                                                                    \
    {                                                                                                                                                                                                  \
      learnedParameters.SetElement(i, d3dtransform->getParameters()[i]);                                                                                                                               \
    }                                                                                                                                                                                                  \
    transform->SetParameters(learnedParameters);                                                                                                                                                       \
    resample->SetTransform(transform);                                                                                                                                                                 \
  }                                                                                                                                                                                                    \
  resample->SetInput(itkMovingImage);                                                                                                                                                                  \
  resample->SetSize(transformhelper.FixedSize);                                                                                                                                                        \
  resample->SetOutputOrigin(transformhelper.FixedOrigin);                                                                                                                                              \
  resample->SetOutputSpacing(transformhelper.FixedSpacing);                                                                                                                                            \
  resample->SetOutputDirection(transformhelper.FixedDirection);                                                                                                                                        \
  if(m_InterpolationType == 1)                                                                                                                                                                         \
  {                                                                                                                                                                                                    \
    using InterpolatorType = itk::NearestNeighborInterpolateImageFunction<ImageType, double>;                                                                                                          \
    typename InterpolatorType::Pointer interpolator = InterpolatorType::New();                                                                                                                         \
    resample->SetInterpolator(interpolator);                                                                                                                                                           \
  }                                                                                                                                                                                                    \
  resample->Update();                                                                                                                                                                                  \
  typename ImageType::Pointer output = resample->GetOutput();                                                                                                                                          \
  typename ImageType::SizeType size = output->GetLargestPossibleRegion().GetSize();                                                                                                                    \
  using ToD3DType = itk::InPlaceImageToDream3DDataFilter<PixelType, ImageDimension>;                                                                                                                   \
  typename ToD3DType::Pointer movedD3D = ToD3DType::New();                                                                                                                                             \
  DataContainer::Pointer dc = getDataContainerArray()->getDataContainer(getDataContainerName());                                                                                                       \
  AttributeMatrix::Pointer am = dc->getAttributeMatrix(m_CellAttributeMatrixName);                                                                                                                     \
  QVector<size_t> tdims = {size[0], size[1], 1};                                                                                                                                                       \
  am->setTupleDimensions(tdims);                                                                                                                                                                       \
  movedD3D->SetDataContainer(dc);                                                                                                                                                                      \
  movedD3D->SetDataArrayName(m_ImageDataArrayName.toStdString());                                                                                                                                      \
  movedD3D->SetAttributeMatrixArrayName(m_CellAttributeMatrixName.toStdString());                                                                                                                      \
  movedD3D->InPlaceOn();                                                                                                                                                                               \
  movedD3D->SetInput(output);                                                                                                                                                                          \
  movedD3D->Update();

#define BSPLINESETUP()                                                                                                                                                                                 \
  using TransformType = itk::BSplineTransform<double, ImageDimension, SplineOrder>;                                                                                                                    \
  typename TransformType::Pointer transform = TransformType::New();                                                                                                                                    \
  typename TransformType::MeshSizeType meshSize;                                                                                                                                                       \
  meshSize.Fill(3);                                                                                                                                                                                    \
  transform->SetTransformDomainMeshSize(meshSize);                                                                                                                                                     \
  int32_t numFixedElements = d3dtransform->getFixedParameters().size();                                                                                                                                \
  typename TransformType::FixedParametersType fixedParameters(numFixedElements);                                                                                                                       \
  for(int32_t i = 0; i < numFixedElements; i++)                                                                                                                                                        \
  {                                                                                                                                                                                                    \
    fixedParameters.SetElement(i, d3dtransform->getFixedParameters()[i]);                                                                                                                              \
  }                                                                                                                                                                                                    \
  transform->SetFixedParameters(fixedParameters);                                                                                                                                                      \
  int32_t numLearnedElements = d3dtransform->getParameters().size();                                                                                                                                   \
  typename TransformType::ParametersType learnedParameters(numLearnedElements);                                                                                                                        \
  for(int32_t i = 0; i < numLearnedElements; i++)                                                                                                                                                      \
  {                                                                                                                                                                                                    \
    learnedParameters.SetElement(i, d3dtransform->getParameters()[i]);                                                                                                                                 \
  }                                                                                                                                                                                                    \
  transform->SetParameters(learnedParameters);                                                                                                                                                         \
  resample->SetTransform(transform);

// -----------------------------------------------------------------------------
// END MACROS
// -----------------------------------------------------------------------------

template <typename T>
void ITKResampleImage::Resample2D(const QString& sliceNo)
{
  const unsigned int ImageDimension = 2;

  if(m_MovingImagePtr.lock()->getComponentDimensions()[0] == 1)
  {
    using PixelType = itk::Vector<T, 1>;
    RESAMPLE2D_VARIABLECOMPONENT(sliceNo)
  }
  else if(m_MovingImagePtr.lock()->getComponentDimensions()[0] == 2)
  {
    using PixelType = itk::Vector<T, 2>;
    RESAMPLE2D_VARIABLECOMPONENT(sliceNo)
  }

  else if(m_MovingImagePtr.lock()->getComponentDimensions()[0] == 3)
  {
    using PixelType = itk::Vector<T, 3>;
    RESAMPLE2D_VARIABLECOMPONENT(sliceNo)
  }

  else if(m_MovingImagePtr.lock()->getComponentDimensions()[0] == 4)
  {
    using PixelType = itk::Vector<T, 4>;
    RESAMPLE2D_VARIABLECOMPONENT(sliceNo)
  }
  else if(m_MovingImagePtr.lock()->getComponentDimensions()[0] == 10)
  {
    using PixelType = itk::Vector<T, 10>;
    RESAMPLE2D_VARIABLECOMPONENT(sliceNo)
  }
  else if(m_MovingImagePtr.lock()->getComponentDimensions()[0] == 11)
  {
    using PixelType = itk::Vector<T, 11>;
    RESAMPLE2D_VARIABLECOMPONENT(sliceNo)
  }
}
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------

ITKTransformHelpers ITKResampleImage::getTransformAndFixedParams(QString sliceNo)
{
  hid_t fileId = QH5Utilities::openFile(m_TransformFileName, true);
  ITKTransformHelpers transformhelper(fileId, sliceNo.toStdString());
  transformhelper.readFixedInformation();
  transformhelper.readTransform();

  return transformhelper;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ITKResampleImage::SeriesResampling()
{
  QVector<QString> inputImageList = getFileList(m_ImageFileListInfo);

  FilterManager* fm = FilterManager::Instance();
  IFilterFactory::Pointer factory = fm->getFactoryFromClassName("ITKImageReader");
  IFilterFactory::Pointer factory2 = fm->getFactoryFromClassName("ITKImageWriter");
  if(factory.get() == nullptr)
  {
    QString ss = QObject::tr("Unable to instantiate Filter with name 'ITKImageReader'\n"
                             "The 'ITKImageReader' Filter is needed to import the image");
    setErrorCondition(-1);
    notifyErrorMessage(getHumanLabel(), ss, getErrorCondition());
  }
  if(factory2.get() == nullptr)
  {
    QString ss = QObject::tr("Unable to instantiate Filter with name 'ITKImageWriter'\n"
                             "The 'ITKImageWriter' Filter is needed to import the image");
    setErrorCondition(-1);
    notifyErrorMessage(getHumanLabel(), ss, getErrorCondition());
  }
  AbstractFilter::Pointer itkImageReader = factory->create();
  AbstractFilter::Pointer itkImageWriter = factory2->create();

  size_t startIndex = m_ImageFileListInfo.StartIndex;
  for(int32_t i = 0; i < inputImageList.size(); i++)
  {
    DataArrayPath movingpath("_INTERNAL_USE_ONLY_MovingImageDataContainerName", "_INTERNAL_USE_ONLY_attributeMatrixName", "_INTERNAL_USE_ONLY_imageDataArrayName");
    setMovingImageArrayPath(movingpath);

    itkImageReader->setDataContainerArray(getDataContainerArray());

    QVariant dcName;
    QString pathname = "_INTERNAL_USE_ONLY_MovingImageDataContainerName";
    dcName.setValue(pathname);
    itkImageReader->setProperty("DataContainerName", dcName);

    QVariant amName;
    amName.setValue(m_MovingImageArrayPath.getAttributeMatrixName());
    itkImageReader->setProperty("CellAttributeMatrixName", amName);

    QVariant imDAName;
    imDAName.setValue(m_MovingImageArrayPath.getDataArrayName());
    itkImageReader->setProperty("ImageDataArrayName", imDAName);

    QVariant var;
    var.setValue(inputImageList[i]);
    itkImageReader->setProperty("FileName", var);

    itkImageReader->execute();

    if(itkImageReader->getErrorCondition() < 0)
    {
      setErrorCondition(itkImageReader->getErrorCondition());
      notifyErrorMessage(getHumanLabel(), "Error Reading Input Image.", getErrorCondition());

      return;
    }

    DataContainer::Pointer m = getDataContainerArray()->createNonPrereqDataContainer<AbstractFilter>(this, getDataContainerName());
    if(getErrorCondition() < 0)
    {
      return;
    }

    // Create the Image Geometry
    ImageGeom::Pointer image = ImageGeom::CreateGeometry(SIMPL::Geometry::ImageGeometry);
    m->setGeometry(image);

    QVector<size_t> tDims(3, 0);
    AttributeMatrix::Pointer cellAttrMat = m->createNonPrereqAttributeMatrix(this, getCellAttributeMatrixName(), tDims, AttributeMatrix::Type::Cell);
    if(getErrorCondition() < 0)
    {
      return;
    }

    m_MovingImagePtr =
        getDataContainerArray()->getPrereqIDataArrayFromPath<IDataArray, AbstractFilter>(this, getMovingImageArrayPath()); /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */

    EXECUTE_FUNCTION_TEMPLATE_NO_BOOL(this, createCompatibleArrays, m_MovingImagePtr.lock())

    EXECUTE_FUNCTION_TEMPLATE_NO_BOOL(this, Resample2D, m_MovingImagePtr.lock(), QString::number(i))

    itkImageWriter->setDataContainerArray(getDataContainerArray());

    QVariant planeValue;
    planeValue.setValue(0);
    itkImageWriter->setProperty("Plane", planeValue);

    QString path = m_OutputPath + "/" + m_FileNamePrefix + QString::number(startIndex) + "." + m_ImageFileListInfo.FileExtension;

    itkImageWriter->setProperty("FileName", path);

    QVariant imArrayPath;
    DataArrayPath outpathname(m_DataContainerName, m_CellAttributeMatrixName, m_ImageDataArrayName);
    imArrayPath.setValue(outpathname);
    itkImageWriter->setProperty("ImageArrayPath", imArrayPath);
    itkImageWriter->execute();

    startIndex++;

    getDataContainerArray()->removeDataContainer("_INTERNAL_USE_ONLY_MovingImageDataContainerName");
    getDataContainerArray()->removeDataContainer(m_DataContainerName);
  }
}
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString ITKResampleImage::GetAndModifyHeader(QFile& reader, size_t dims[3], float spacing[3])
{
  QString ctfheader("");
  QString altLine("");
  QList<QByteArray> headerLines;
  QByteArray buf;
  bool ok = false;
  int32_t numPhases = -1;
  bool headerComplete = false;
  while(!reader.atEnd() && !headerComplete)
  {
    buf = reader.readLine();
    if(buf.startsWith("XStep"))
    {
      altLine = "XStep	" + QString::number(spacing[0], 'f', 13) + '\n';
      ctfheader.append(altLine);
    }
    else if(buf.startsWith("YStep"))
    {
      altLine = "YStep	" + QString::number(spacing[0], 'f', 13) + '\n';
      //+ QString::number(spacing[1]) +
      ctfheader.append(altLine);
    }
    else if(buf.startsWith("XCells"))
    {
      altLine = "XCells	" + QString::number(dims[0]) + '\n';
      ctfheader.append(altLine);
    }
    else if(buf.startsWith("YCells"))
    {
      altLine = "YCells	" + QString::number(dims[1]) + '\n';
      ctfheader.append(altLine);
    }
    else
    {
      ctfheader.append(QString(buf));
    }
    // Append the line to the complete header

    // remove the newline at the end of the line
    buf.chop(1);
    headerLines.push_back(buf);
    if(buf.startsWith("Phases"))
    {
      QList<QByteArray> tokens = buf.split('\t');
      numPhases = tokens.at(1).toInt(&ok, 10);
      break; //
    }
  }

  // Now read the phases line
  for(int32_t p = 0; p < numPhases; ++p)
  {
    buf = reader.readLine();
    ctfheader.append(QString(buf));

    // remove the newline at the end of the line
    buf.chop(1);
    headerLines.push_back(buf);
  }

  // one more line for column headings
  buf = reader.readLine();
  ctfheader.append(QString(buf));

  headerComplete = true;

  return ctfheader;
}
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ITKResampleImage::WriteResampledCTFfile(const QString& filename, std::vector<IDataArray::Pointer> ctfArrays, size_t dims[3], float spacing[3], size_t index)
{
  QFile in(filename);
  if(!in.open(QIODevice::ReadOnly | QIODevice::Text))
  {
    QString msg = QString("Ctf file could not be opened: ") + filename;
    setErrorCondition(-100);
    notifyErrorMessage(getHumanLabel(), msg, getErrorCondition());

    return;
  }

  QString header = GetAndModifyHeader(in, dims, spacing);
  QString outfile = m_OutputPathEBSD + "/" + m_EBSDFileNamePrefix + QString::number(index) + "." + m_OrientationFileListInfo.FileExtension;

  QFile file(outfile);
  if(!file.open(QIODevice::WriteOnly | QIODevice::Text))
  {
    return;
  }

  QTextStream out(&file);
  out << header;

  EBSDWriterFactory* factory = EBSDWriterFactory::Instance();

  int32_t arrayIndex = 0;
  for(int32_t i = 0; i < dims[0]; i++)
  {
    for(int32_t j = 0; j < dims[1]; j++)
    {
      QString rowData("");
      for(int32_t k = 0; k < ctfArrays.size(); k++)
      {
        std::function<QString(QString, IDataArray::Pointer, uint64_t)> writer = factory->getWriter(m_ctfDataArrayNames[k]);
        rowData = writer(rowData, ctfArrays[k], arrayIndex);
        if(k == (ctfArrays.size() - 1))
        {
          rowData.append('\n');
        }
        else
        {
          rowData.append("	");
        }
      }
      out << rowData;
      arrayIndex++;
    }
  }
}
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ITKResampleImage::EBSDSeriesResampling()
{
  QVector<QString> inputEBSDFileList = getFileList(m_OrientationFileListInfo);

  FilterManager* fm = FilterManager::Instance();

  IFilterFactory::Pointer factory = fm->getFactoryFromClassName("ReadCtfData");
  if(factory.get() == nullptr)
  {
    QString ss = QObject::tr("Unable to instantiate Filter with name 'ReadCtfData'\n"
                             "The 'ReadCtfData' Filter is needed to import EBSD data");
    setErrorCondition(-1);
    notifyErrorMessage(getHumanLabel(), ss, getErrorCondition());
  }
  AbstractFilter::Pointer ctfReader = factory->create();

  size_t index = m_OrientationFileListInfo.StartIndex;
  for(int32_t i = 0; i < inputEBSDFileList.size(); i++)
  {
    ctfReader->setDataContainerArray(getDataContainerArray());

    QVariant dcName;
    QString pathname = "_INTERNAL_USE_ONLY_MovingEBSDDataContainerName";
    dcName.setValue(pathname);
    ctfReader->setProperty("DataContainerName", dcName);

    QVariant amName;
    amName.setValue(QString("_INTERNAL_USE_ONLY_attributeMatrixName"));
    ctfReader->setProperty("CellAttributeMatrixName", amName);

    QVariant ensembleName;
    ensembleName.setValue(QString("_INTERNAL_USE_ONLY_ensembleAttributeMatrixName"));
    ctfReader->setProperty("CellEnsembleAttributeMatrixName", ensembleName);

    QVariant var;
    var.setValue(inputEBSDFileList[i]);
    ctfReader->setProperty("InputFile", var);

    QVariant degtorad;
    degtorad.setValue(false);
    ctfReader->setProperty("DegreesToRadians", degtorad);

    QVariant hexalign;
    hexalign.setValue(false);
    ctfReader->setProperty("EdaxHexagonalAlignment", hexalign);

    ctfReader->execute();

    DataArrayPath movingpath("_INTERNAL_USE_ONLY_MovingEBSDDataContainerName", "_INTERNAL_USE_ONLY_attributeMatrixName", "");
    setMovingImageArrayPath(movingpath);

    std::vector<IDataArray::Pointer> ctfArrays(m_ctfDataArrayNames.size(), nullptr);

    ImageGeom::Pointer finalImage;
    size_t dims[3];
    float spacing[3];

    for(int32_t j = 0; j < m_ctfDataArrayNames.size(); j++)
    {
      movingpath.setDataArrayName(m_ctfDataArrayNames[j]);
      setMovingImageArrayPath(movingpath);

      DataContainer::Pointer m = getDataContainerArray()->createNonPrereqDataContainer<AbstractFilter>(this, getDataContainerName());
      if(getErrorCondition() < 0)
      {
        return;
      }

      // Create the Image Geometry
      ImageGeom::Pointer image = ImageGeom::CreateGeometry(SIMPL::Geometry::ImageGeometry);
      m->setGeometry(image);

      QVector<size_t> tDims(3, 0);
      AttributeMatrix::Pointer cellAttrMat = m->createNonPrereqAttributeMatrix(this, getCellAttributeMatrixName(), tDims, AttributeMatrix::Type::Cell);
      if(getErrorCondition() < 0)
      {
        return;
      }

      m_MovingImagePtr =
          getDataContainerArray()->getPrereqIDataArrayFromPath<IDataArray, AbstractFilter>(this, getMovingImageArrayPath()); /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
      EXECUTE_FUNCTION_TEMPLATE_NO_BOOL(this, createCompatibleArrays, m_MovingImagePtr.lock())
      EXECUTE_FUNCTION_TEMPLATE_NO_BOOL(this, Resample2D, m_MovingImagePtr.lock(), QString::number(i))

      finalImage = getDataContainerArray()->getDataContainer(m_DataContainerName)->getGeometryAs<ImageGeom>();
      std::tie(dims[0], dims[1], dims[2]) = finalImage->getDimensions();
      finalImage->getResolution(spacing);

      if(m_ctfDataArrayNames[j] == "X")
      {
        IDataArray::Pointer X = getDataContainerArray()->getDataContainer(m_DataContainerName)->getAttributeMatrix(m_CellAttributeMatrixName)->getAttributeArray(m_ImageDataArrayName);
        DataArray<float>::Pointer XFloat = std::dynamic_pointer_cast<DataArray<float>>(X);

        int32_t xIndex = 0;
        for(int32_t yVals = 0; yVals < dims[1]; yVals++)
        {
          for(int32_t xVals = 0; xVals < dims[0]; xVals++)
          {
            XFloat->setValue(xIndex, xVals * spacing[0]);
            xIndex++;
          }
        }
      }

      if(m_ctfDataArrayNames[j] == "Y")
      {

        IDataArray::Pointer Y = getDataContainerArray()->getDataContainer(m_DataContainerName)->getAttributeMatrix(m_CellAttributeMatrixName)->getAttributeArray(m_ImageDataArrayName);
        DataArray<float>::Pointer YFloat = std::dynamic_pointer_cast<DataArray<float>>(Y);

        int32_t yIndex = 0;
        for(int32_t yVals = 0; yVals < dims[1]; yVals++)
        {
          for(int32_t xVals = 0; xVals < dims[0]; xVals++)
          {
            YFloat->setValue(yIndex, yVals * spacing[0]);
            yIndex++;
          }
        }
      }
      // This is a vector of all the arrays stored from reading in the CTF file
      ctfArrays[j] = getDataContainerArray()->getDataContainer(m_DataContainerName)->getAttributeMatrix(m_CellAttributeMatrixName)->getAttributeArray(m_ImageDataArrayName)->deepCopy();

      getDataContainerArray()->removeDataContainer(m_DataContainerName);
    }
    getDataContainerArray()->removeDataContainer("_INTERNAL_USE_ONLY_MovingEBSDDataContainerName");

    WriteResampledCTFfile(inputEBSDFileList[i], ctfArrays, dims, spacing, index);

    index++;
  }
}
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ITKResampleImage::execute()
{
  initialize();
  dataCheck();
  if(getErrorCondition() < 0)
  {
    return;
  }

  if(getCancel())
  {
    return;
  }

  if(getWarningCondition() < 0)
  {
    QString ss = QObject::tr("Some warning message");
    setWarningCondition(-88888888);
    notifyWarningMessage(getHumanLabel(), ss, getErrorCondition());
  }

  if(getErrorCondition() < 0)
  {
    QString ss = QObject::tr("Some error message");
    setErrorCondition(-99999999);
    notifyErrorMessage(getHumanLabel(), ss, getErrorCondition());

    return;
  }

  if(m_OperationMode == 0)
  {
    m_MovingImagePtr =
        getDataContainerArray()->getPrereqIDataArrayFromPath<IDataArray, AbstractFilter>(this, getMovingImageArrayPath()); /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
    EXECUTE_FUNCTION_TEMPLATE_NO_BOOL(this, Resample2D, m_MovingImagePtr.lock(), "0")
  }
  else if(m_OperationMode == 1)
  {
    SeriesResampling();
  }
  else if(m_OperationMode == 2)
  {
    EBSDSeriesResampling();
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AbstractFilter::Pointer ITKResampleImage::newFilterInstance(bool copyFilterParameters) const
{
  ITKResampleImage::Pointer filter = ITKResampleImage::New();
  if(copyFilterParameters)
  {
    copyFilterParameterInstanceVariables(filter.get());
  }
  return filter;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString ITKResampleImage::getCompiledLibraryName() const
{
  return DREAM3DReviewConstants::DREAM3DReviewBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString ITKResampleImage::getBrandingString() const
{
  return "DREAM3DReview";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString ITKResampleImage::getFilterVersion() const
{
  QString version;
  QTextStream vStream(&version);
  vStream << DREAM3DReview::Version::Major() << "." << DREAM3DReview::Version::Minor() << "." << DREAM3DReview::Version::Patch();
  return version;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString ITKResampleImage::getGroupName() const
{
  return SIMPL::FilterGroups::ReconstructionFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString ITKResampleImage::getSubGroupName() const
{
  return DREAM3DReviewConstants::FilterSubGroups::RegistrationFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString ITKResampleImage::getHumanLabel() const
{
  return "ITK::Resample Image";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QUuid ITKResampleImage::getUuid()
{
  return QUuid("{7d478bf6-1acc-5e49-86b6-45c94776bc48}");
}
