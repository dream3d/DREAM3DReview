/*
 * Your License or Copyright can go here
 */
#include "ITKPairwiseImageRegistration.h"

#include <QtCore/QFileInfo>

#include "itkAffineTransform.h"
#include "itkBSplineTransform.h"
#include "itkBSplineTransformInitializer.h"
#include "itkCastImageFilter.h"
#include "itkCenteredTransformInitializer.h"
#include "itkEuler2DTransform.h"
#include "itkFFTNormalizedCorrelationImageFilter.h"
#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "itkMaskedFFTNormalizedCorrelationImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkResampleImageFilter.h"
#include "itkTransformFileWriter.h"
// itkVersion 4 specific filters
#include "itkAmoebaOptimizerv4.h"
#include "itkCorrelationImageToImageMetricv4.h"
#include "itkGradientDescentLineSearchOptimizerv4.h"
#include "itkImageRegistrationMethodv4.h"
#include "itkLBFGSOptimizerv4.h"
#include "itkMattesMutualInformationImageToImageMetricv4.h"
#include "itkMeanSquaresImageToImageMetricv4.h"

#include "H5Support/H5ScopedSentinel.h"
#include "H5Support/QH5Utilities.h"

#include "SIMPLib/Common/Constants.h"
#include "SIMPLib/Common/TemplateHelpers.h"
#include "SIMPLib/FilterParameters/BooleanFilterParameter.h"
#include "SIMPLib/FilterParameters/ChoiceFilterParameter.h"
#include "SIMPLib/FilterParameters/DataArraySelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/FileListInfoFilterParameter.h"
#include "SIMPLib/FilterParameters/FloatFilterParameter.h"
#include "SIMPLib/FilterParameters/FloatVec3FilterParameter.h"
#include "SIMPLib/FilterParameters/IntFilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedBooleanFilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedChoicesFilterParameter.h"
#include "SIMPLib/FilterParameters/OutputFileFilterParameter.h"
#include "SIMPLib/FilterParameters/SeparatorFilterParameter.h"
#include "SIMPLib/Filtering/FilterFactory.hpp"
#include "SIMPLib/Filtering/FilterManager.h"
#include "SIMPLib/Geometry/ImageGeom.h"
#include "SIMPLib/Geometry/TransformContainer.h"
#include "SIMPLib/ITK/itkInPlaceDream3DDataToImageFilter.h"
#include "SIMPLib/ITK/itkInPlaceImageToDream3DDataFilter.h"
#include "SIMPLib/Utilities/FilePathGenerator.h"
#include "SIMPLib/Utilities/FileSystemPathHelper.h"

#include "DREAM3DReview/DREAM3DReviewConstants.h"
#include "DREAM3DReview/DREAM3DReviewVersion.h"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
ITKPairwiseImageRegistration::ITKPairwiseImageRegistration()
: m_ChangeOriginResolution(false)
, m_InitializeRigidWithFFT(true)
, m_InitializeAffineWithFFT(true)
, m_TransformFile("")
, m_FixedImageArrayPath("", "", "")
, m_MovingImageArrayPath("", "", "")
, m_OperationMode(0)
, m_TransformType(0)
, m_Metric(0)
, m_Optimizer(0)
, m_BSplineOrder(3)
, m_MeshSize(3)
, m_MaximumNumberOfFunctionEvaluations(1000)
, m_NumberOfIterations(100)
, m_GradientConvergenceTolerance(1e-5)
, m_LearningRate(0.5)
, m_ReqFracOverlapPxlsRigid(0.5)
, m_ReqFracOverlapPxlsAffine(0.5)
{

  m_MovingFileListInfo.FileExtension = QString("tif");
  m_MovingFileListInfo.StartIndex = 0;
  m_MovingFileListInfo.EndIndex = 0;
  m_MovingFileListInfo.PaddingDigits = 0;

  m_FixedFileListInfo.FileExtension = QString("tif");
  m_FixedFileListInfo.StartIndex = 0;
  m_FixedFileListInfo.EndIndex = 0;
  m_FixedFileListInfo.PaddingDigits = 0;

  m_MovingOrigin.x = 0;
  m_MovingOrigin.y = 0;
  m_MovingOrigin.z = 0;

  m_MovingSpacing.x = 1.0;
  m_MovingSpacing.y = 1.0;
  m_MovingSpacing.z = 1.0;

  m_FixedOrigin.x = 0;
  m_FixedOrigin.y = 0;
  m_FixedOrigin.z = 0;

  m_FixedSpacing.x = 1.0;
  m_FixedSpacing.y = 1.0;
  m_FixedSpacing.z = 1.0;

  initialize();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
ITKPairwiseImageRegistration::~ITKPairwiseImageRegistration() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ITKPairwiseImageRegistration::initialize()
{
  setErrorCondition(0);
  setWarningCondition(0);
  setCancel(false);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ITKPairwiseImageRegistration::setupFilterParameters()
{
  FilterParameterVector parameters;

  /// MODE TYPE: Either a single pair or series
  {
    LinkedChoicesFilterParameter::Pointer parameter = LinkedChoicesFilterParameter::New();
    parameter->setHumanLabel("Operation Mode");
    parameter->setPropertyName("OperationMode");
    parameter->setSetterCallback(SIMPL_BIND_SETTER(ITKPairwiseImageRegistration, this, OperationMode));
    parameter->setGetterCallback(SIMPL_BIND_GETTER(ITKPairwiseImageRegistration, this, OperationMode));

    QVector<QString> choices;
    choices.push_back("Single Pair of Images");
    choices.push_back("Series of Pairs");

    parameter->setChoices(choices);

    QStringList linkedProps;
    linkedProps << "FixedFileListInfo"
                << "MovingFileListInfo"
                << "FixedImageArrayPath"
                << "MovingImageArrayPath";
    parameter->setLinkedProperties(linkedProps);
    parameter->setEditable(false);
    parameter->setCategory(FilterParameter::Parameter);
    parameters.push_back(parameter);
  }

  /// TRANSFORMATION TYPE: B-spline, affine or rigid
  {
    LinkedChoicesFilterParameter::Pointer parameter = LinkedChoicesFilterParameter::New();
    parameter->setHumanLabel("Transformation Type");
    parameter->setPropertyName("TransformType");
    parameter->setSetterCallback(SIMPL_BIND_SETTER(ITKPairwiseImageRegistration, this, TransformType));
    parameter->setGetterCallback(SIMPL_BIND_GETTER(ITKPairwiseImageRegistration, this, TransformType));

    QVector<QString> choices;
    choices.push_back("BSpline");
    choices.push_back("Rigid");
    choices.push_back("Affine");
    parameter->setChoices(choices);

    QStringList linkedProps;
    linkedProps << "MeshSize"
                << "BSplineOrder"
                << "InitializeRigidWithFFT"
                << "InitializeAffineWithFFT"
                << "ReqFracOverlapPxlsRigid"
                << "ReqFracOverlapPxlsAffine";
    parameter->setLinkedProperties(linkedProps);
    parameter->setEditable(false);
    parameter->setCategory(FilterParameter::Parameter);
    parameters.push_back(parameter);
  }

  ChoiceFilterParameter::Pointer metricparameter = ChoiceFilterParameter::New();
  metricparameter->setHumanLabel("Evaluation Metric");
  metricparameter->setPropertyName("Metric");
  metricparameter->setSetterCallback(SIMPL_BIND_SETTER(ITKPairwiseImageRegistration, this, Metric));
  metricparameter->setGetterCallback(SIMPL_BIND_GETTER(ITKPairwiseImageRegistration, this, Metric));

  QVector<QString> metricchoices;
  metricchoices.push_back("Mattes Mutual Information");
  metricchoices.push_back("Mean Squares");
  metricchoices.push_back("Normalized Cross Correlation");
  metricparameter->setChoices(metricchoices);
  metricparameter->setEditable(false);
  metricparameter->setCategory(FilterParameter::Parameter);
  parameters.push_back(metricparameter);

  LinkedChoicesFilterParameter::Pointer optimizerparameter = LinkedChoicesFilterParameter::New();
  optimizerparameter->setHumanLabel("Optimizer");
  optimizerparameter->setPropertyName("Optimizer");
  optimizerparameter->setSetterCallback(SIMPL_BIND_SETTER(ITKPairwiseImageRegistration, this, Optimizer));
  optimizerparameter->setGetterCallback(SIMPL_BIND_GETTER(ITKPairwiseImageRegistration, this, Optimizer));

  QVector<QString> optimizerchoices;
  optimizerchoices.push_back("LBFGS");
  optimizerchoices.push_back("Amoeba");
  optimizerchoices.push_back("Gradient Descent");
  optimizerparameter->setChoices(optimizerchoices);

  QStringList optimizerlinkedProps;
  optimizerlinkedProps << "GradientConvergenceTolerance"
                       << "MaximumNumberOfFunctionEvaluations"
                       << "LearningRate";
  optimizerparameter->setLinkedProperties(optimizerlinkedProps);
  optimizerparameter->setEditable(false);
  optimizerparameter->setCategory(FilterParameter::Parameter);
  parameters.push_back(optimizerparameter);

  QStringList linkedProps;
  linkedProps << "MovingOrigin"
              << "FixedOrigin"
              << "FixedSpacing"
              << "MovingSpacing";
  parameters.push_back(SIMPL_NEW_LINKED_BOOL_FP("Origin/Resolution for Moving and Fixed are Different", ChangeOriginResolution, FilterParameter::Parameter, ITKPairwiseImageRegistration, linkedProps));
  parameters.push_back(SIMPL_NEW_FLOAT_VEC3_FP("Moving Origin", MovingOrigin, FilterParameter::Parameter, ITKPairwiseImageRegistration));
  parameters.push_back(SIMPL_NEW_FLOAT_VEC3_FP("Moving Spacing", MovingSpacing, FilterParameter::Parameter, ITKPairwiseImageRegistration));
  parameters.push_back(SIMPL_NEW_FLOAT_VEC3_FP("Fixed Origin", FixedOrigin, FilterParameter::Parameter, ITKPairwiseImageRegistration));
  parameters.push_back(SIMPL_NEW_FLOAT_VEC3_FP("Fixed Spacing", FixedSpacing, FilterParameter::Parameter, ITKPairwiseImageRegistration));

  DataArraySelectionFilterParameter::RequirementType dasReq =
      DataArraySelectionFilterParameter::CreateRequirement(SIMPL::Defaults::AnyPrimitive, 1, AttributeMatrix::Type::Cell, IGeometry::Type::Image);
  parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Fixed Image", FixedImageArrayPath, FilterParameter::RequiredArray, ITKPairwiseImageRegistration, dasReq, 0));
  parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Moving Image", MovingImageArrayPath, FilterParameter::RequiredArray, ITKPairwiseImageRegistration, dasReq, 0));
  parameters.push_back(SIMPL_NEW_INTEGER_FP("B-Spline Order", BSplineOrder, FilterParameter::Parameter, ITKPairwiseImageRegistration, 0));
  parameters.push_back(SIMPL_NEW_INTEGER_FP("Mesh Size", MeshSize, FilterParameter::Parameter, ITKPairwiseImageRegistration, 0));
  parameters.push_back(SIMPL_NEW_FLOAT_FP("Gradient Convergence Tolerance", GradientConvergenceTolerance, FilterParameter::Parameter, ITKPairwiseImageRegistration, 0));
  parameters.push_back(SIMPL_NEW_INTEGER_FP("Num Max Function Evaluations", MaximumNumberOfFunctionEvaluations, FilterParameter::Parameter, ITKPairwiseImageRegistration, 0));
  parameters.push_back(SIMPL_NEW_INTEGER_FP("Num Iterations", NumberOfIterations, FilterParameter::Parameter, ITKPairwiseImageRegistration));
  parameters.push_back(SIMPL_NEW_FLOAT_FP("Learning Rate", LearningRate, FilterParameter::Parameter, ITKPairwiseImageRegistration, 2));
  parameters.push_back(SIMPL_NEW_OUTPUT_FILE_FP("Output Transform File", TransformFile, FilterParameter::Parameter, ITKPairwiseImageRegistration));
  parameters.push_back(SIMPL_NEW_BOOL_FP("Initialize Translation With FFT", InitializeRigidWithFFT, FilterParameter::Parameter, ITKPairwiseImageRegistration, 1));
  parameters.push_back(SIMPL_NEW_BOOL_FP("Initialize Translation With FFT", InitializeAffineWithFFT, FilterParameter::Parameter, ITKPairwiseImageRegistration, 2));
  parameters.push_back(SIMPL_NEW_FLOAT_FP("Required Fraction of Overlapping Pixels for FFT", ReqFracOverlapPxlsRigid, FilterParameter::Parameter, ITKPairwiseImageRegistration, 1));
  parameters.push_back(SIMPL_NEW_FLOAT_FP("Required Fraction of Overlapping Pixels for FFT", ReqFracOverlapPxlsAffine, FilterParameter::Parameter, ITKPairwiseImageRegistration, 2));

  parameters.push_back(SeparatorFilterParameter::New("Fixed Image List", FilterParameter::Parameter));
  parameters.push_back(SIMPL_NEW_FILELISTINFO_FP("Fixed Image File List", FixedFileListInfo, FilterParameter::Parameter, ITKPairwiseImageRegistration));
  FileListInfoFilterParameter::Pointer fixedFileList = std::dynamic_pointer_cast<FileListInfoFilterParameter>(parameters.back());
  fixedFileList->setGroupIndex(1);

  parameters.push_back(SeparatorFilterParameter::New("Moving Image List", FilterParameter::Parameter));
  parameters.push_back(SIMPL_NEW_FILELISTINFO_FP("Moving Image File List", MovingFileListInfo, FilterParameter::Parameter, ITKPairwiseImageRegistration));
  FileListInfoFilterParameter::Pointer movingFileList = std::dynamic_pointer_cast<FileListInfoFilterParameter>(parameters.back());
  movingFileList->setGroupIndex(1);

  setFilterParameters(parameters);
}

QVector<QString> ITKPairwiseImageRegistration::getFileList(FileListInfo_t inputFileListInfo)
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

int ITKPairwiseImageRegistration::checkInputFileList(FileListInfo_t inputFileListInfo)
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

  return fileList.size();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ITKPairwiseImageRegistration::dataCheck()
{
  setErrorCondition(0);
  setWarningCondition(0);

  FileSystemPathHelper::CheckOutputFile(this, "Output Transform File", getTransformFile(), true);

  QVector<size_t> cDims(1, 1);
  if(m_OperationMode == 0)
  {
    ImageGeom::Pointer fixedImage = getDataContainerArray()->getPrereqGeometryFromDataContainer<ImageGeom, AbstractFilter>(this, getFixedImageArrayPath().getDataContainerName());
    ImageGeom::Pointer movingImage = getDataContainerArray()->getPrereqGeometryFromDataContainer<ImageGeom, AbstractFilter>(this, getMovingImageArrayPath().getDataContainerName());

    m_FixedImagePtr = getDataContainerArray()->getPrereqIDataArrayFromPath<IDataArray, AbstractFilter>(this, getFixedImageArrayPath());
    m_MovingImagePtr =
        getDataContainerArray()->getPrereqIDataArrayFromPath<IDataArray, AbstractFilter>(this, getMovingImageArrayPath()); /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
  }
  else if(m_OperationMode == 1)
  {
    int numFixedImages = checkInputFileList(m_FixedFileListInfo);
    int numMovingImages = checkInputFileList(m_MovingFileListInfo);

    if(numFixedImages != numMovingImages)
    {
      QString ss = "The number of fixed images and moving images are required to be the same";
      setErrorCondition(-64503);
      notifyErrorMessage(getHumanLabel(), ss, getErrorCondition());
    }
  }

  if(m_BSplineOrder > 3)
  {
    QString ss = QObject::tr("The largest B-spline order is 3rd Order");
    setErrorCondition(-5555);
    notifyErrorMessage(getHumanLabel(), ss, getErrorCondition());
  }

  if(m_MeshSize < m_BSplineOrder)
  {
    QString ss = QObject::tr("The mesh size cannot be smaller than the B-spline order");
    setErrorCondition(-5556);
    notifyErrorMessage(getHumanLabel(), ss, getErrorCondition());
  }

  if(m_TransformType == 1)
  {
    if((m_FixedSpacing.x != m_MovingSpacing.x || m_FixedSpacing.y != m_MovingSpacing.y) && m_InitializeRigidWithFFT)
    {
      QString ss = QObject::tr("The spacing must be the same for fixed and moving to initialize with FFT");
      setErrorCondition(-5557);
      notifyErrorMessage(getHumanLabel(), ss, getErrorCondition());
    }

    if(m_InitializeRigidWithFFT && (m_ReqFracOverlapPxlsRigid > 1 || m_ReqFracOverlapPxlsRigid < 0))
    {
      QString ss = QObject::tr("The pixel fraction overlap must be between 0 and 1");
      setErrorCondition(-5559);
      notifyErrorMessage(getHumanLabel(), ss, getErrorCondition());
    }
  }

  if(m_TransformType == 2)
  {
    if((m_FixedSpacing.x != m_MovingSpacing.x || m_FixedSpacing.y != m_MovingSpacing.y) && m_InitializeAffineWithFFT)
    {
      QString ss = QObject::tr("The spacing must be the same for fixed and moving to initialize with FFT");
      setErrorCondition(-5558);
      notifyErrorMessage(getHumanLabel(), ss, getErrorCondition());
    }

    if(m_InitializeAffineWithFFT && (m_ReqFracOverlapPxlsAffine > 1 || m_ReqFracOverlapPxlsAffine < 0))
    {
      QString ss = QObject::tr("The pixel fraction overlap must be between 0 and 1");
      setErrorCondition(-5559);
      notifyErrorMessage(getHumanLabel(), ss, getErrorCondition());
    }
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ITKPairwiseImageRegistration::preflight()
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
void ITKPairwiseImageRegistration::convertToDouble(std::weak_ptr<IDataArray> imagePtr, const DataArrayPath& dataArrayPath, const QString& dataArrayName, const DataContainerArray::Pointer& dca)
{
  typename DataArray<T>::Pointer FixedImagePtr = std::dynamic_pointer_cast<DataArray<T>>(imagePtr.lock());

  DataArray<double>::Pointer FixedImageDoublePtr = DataArray<double>::CreateArray(FixedImagePtr->getNumberOfTuples(), dataArrayName, true);

  T* FixedImageDataPtr = FixedImagePtr->getPointer(0);

  double* FixedImageDoubleDataPtr = FixedImageDoublePtr->getPointer(0);

  for(size_t i = 0; i < FixedImagePtr->getNumberOfTuples(); i++)
  {
    FixedImageDoubleDataPtr[i] = static_cast<double>(FixedImageDataPtr[i]);
  }

  dca->getDataContainer(dataArrayPath.getDataContainerName())->getAttributeMatrix(dataArrayPath.getAttributeMatrixName())->addAttributeArray(dataArrayName, FixedImageDoublePtr);
}

//////////////BEGIN MACROS/////////////////////////////////

#define CHOOSE_LBFGS_OPTIMIZER()                                                                                                                                                                       \
  typedef itk::LBFGSOptimizerv4 OptimizerType;                                                                                                                                                         \
  OptimizerType::Pointer optimizer = OptimizerType::New();                                                                                                                                             \
  registration->SetOptimizer(optimizer);                                                                                                                                                               \
  typedef itk::RegistrationParameterScalesFromPhysicalShift<MetricType> ScalesEstimatorType;                                                                                                           \
  ScalesEstimatorType::Pointer scalesEstimator = ScalesEstimatorType::New();                                                                                                                           \
  scalesEstimator->SetMetric(metric);                                                                                                                                                                  \
  scalesEstimator->SetTransformForward(true);                                                                                                                                                          \
  scalesEstimator->SetSmallParameterVariation(1.0);                                                                                                                                                    \
  optimizer->SetGradientConvergenceTolerance(m_GradientConvergenceTolerance);                                                                                                                          \
  optimizer->SetMaximumNumberOfFunctionEvaluations(m_MaximumNumberOfFunctionEvaluations);                                                                                                              \
  optimizer->SetScalesEstimator(scalesEstimator);                                                                                                                                                      \
  optimizer->SetNumberOfIterations(m_NumberOfIterations);

#define CHOOSE_AMOEBA_OPTIMIZER()                                                                                                                                                                      \
  typedef itk::AmoebaOptimizerv4 OptimizerType;                                                                                                                                                        \
  OptimizerType::Pointer optimizer = OptimizerType::New();                                                                                                                                             \
  registration->SetOptimizer(optimizer);                                                                                                                                                               \
  typedef itk::RegistrationParameterScalesFromPhysicalShift<MetricType> ScalesEstimatorType;                                                                                                           \
  ScalesEstimatorType::Pointer scalesEstimator = ScalesEstimatorType::New();                                                                                                                           \
  scalesEstimator->SetMetric(metric);                                                                                                                                                                  \
  scalesEstimator->SetTransformForward(true);                                                                                                                                                          \
  scalesEstimator->SetSmallParameterVariation(1.0);                                                                                                                                                    \
  optimizer->SetScalesEstimator(scalesEstimator);                                                                                                                                                      \
  optimizer->SetNumberOfIterations(m_NumberOfIterations);

#define CHOOSE_GRADIENT_LINE_SEARCH_OPTIMIZER()                                                                                                                                                        \
  typedef itk::GradientDescentLineSearchOptimizerv4 OptimizerType;                                                                                                                                     \
  OptimizerType::Pointer optimizer = OptimizerType::New();                                                                                                                                             \
  registration->SetOptimizer(optimizer);                                                                                                                                                               \
  typedef itk::RegistrationParameterScalesFromPhysicalShift<MetricType> ScalesEstimatorType;                                                                                                           \
  ScalesEstimatorType::Pointer scalesEstimator = ScalesEstimatorType::New();                                                                                                                           \
  scalesEstimator->SetMetric(metric);                                                                                                                                                                  \
  scalesEstimator->SetTransformForward(true);                                                                                                                                                          \
  scalesEstimator->SetSmallParameterVariation(1.0);                                                                                                                                                    \
  optimizer->SetScalesEstimator(scalesEstimator);                                                                                                                                                      \
  optimizer->SetNumberOfIterations(m_NumberOfIterations);

#define CHOOSE_METRIC(metricv4)                                                                                                                                                                        \
  typedef itk::metricv4<ImageType, ImageType> MetricType;                                                                                                                                              \
  MetricType::Pointer metric = MetricType::New();                                                                                                                                                      \
  registration->SetMetric(metric);                                                                                                                                                                     \
  if(m_Optimizer == 0)                                                                                                                                                                                 \
  {                                                                                                                                                                                                    \
    CHOOSE_LBFGS_OPTIMIZER();                                                                                                                                                                          \
  }                                                                                                                                                                                                    \
  if(m_Optimizer == 1)                                                                                                                                                                                 \
  {                                                                                                                                                                                                    \
    CHOOSE_AMOEBA_OPTIMIZER();                                                                                                                                                                         \
  }                                                                                                                                                                                                    \
  if(m_Optimizer == 2)                                                                                                                                                                                 \
  {                                                                                                                                                                                                    \
    CHOOSE_GRADIENT_LINE_SEARCH_OPTIMIZER();                                                                                                                                                           \
  }

#define BSPLINESETUP()                                                                                                                                                                                 \
  using TransformType = itk::BSplineTransform<CoordinateRepType, ImageDimension, SplineOrder>;                                                                                                         \
  typename TransformType::Pointer transform = TransformType::New();                                                                                                                                    \
  using TransformInitializerType = itk::BSplineTransformInitializer<TransformType, ImageType>;                                                                                                         \
  typename TransformInitializerType::Pointer transformInitializer = TransformInitializerType::New();                                                                                                   \
  typename TransformType::MeshSizeType meshSize;                                                                                                                                                       \
  meshSize.Fill(m_MeshSize);                                                                                                                                                                           \
  transformInitializer->SetTransform(transform);                                                                                                                                                       \
  transformInitializer->SetImage(itkFixedImage);                                                                                                                                                       \
  transformInitializer->SetTransformDomainMeshSize(meshSize);                                                                                                                                          \
  transformInitializer->InitializeTransform();                                                                                                                                                         \
  transform->SetIdentity();                                                                                                                                                                            \
  registration->SetInitialTransform(transform);

//////////////END MACROS/////////////////////////////////

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ITKPairwiseImageRegistration::registerImagePair2D(DataContainerArray::Pointer dca, size_t sliceNo, hid_t fileID)
{
  const unsigned int ImageDimension = 2;

  using ImageType = itk::Dream3DImage<double, ImageDimension>;

  using ToITKType = itk::InPlaceDream3DDataToImageFilter<double, ImageDimension>;
  ToITKType::Pointer fixedtoITK = ToITKType::New();
  ToITKType::Pointer movingtoITK = ToITKType::New();

  fixedtoITK->SetDataArrayName("_INTERNAL_USE_ONLY_FixedImageDouble");
  fixedtoITK->SetAttributeMatrixArrayName(m_FixedImageArrayPath.getAttributeMatrixName().toStdString());
  DataContainer::Pointer dcFixed = dca->getDataContainer(m_FixedImageArrayPath.getDataContainerName());
  fixedtoITK->SetInput(dcFixed);
  fixedtoITK->InPlaceOn();
  fixedtoITK->Update();

  ImageType::Pointer itkFixedImage = fixedtoITK->GetOutput();

  movingtoITK->SetDataArrayName("_INTERNAL_USE_ONLY_MovingImageDouble");
  movingtoITK->SetAttributeMatrixArrayName(m_MovingImageArrayPath.getAttributeMatrixName().toStdString());
  DataContainer::Pointer dcMoving = dca->getDataContainer(m_MovingImageArrayPath.getDataContainerName());
  movingtoITK->SetInput(dcMoving);
  movingtoITK->InPlaceOn();
  movingtoITK->Update();

  ImageType::Pointer itkMovingImage = movingtoITK->GetOutput();

  //////////////////////////////// SETUP ORIGIN _ RESOLUTION //////////////////////////////////////////////////

  // if the user selects the origin and resolution to be different for the moving and fixed image, then they must enter in the specific values.
  // if the user says the origin and resolution are the same for the moving and fixed, just set the origin to (0,0) and resolution (1.0, 1.0) for both

  if(!m_ChangeOriginResolution)
  {
    ImageType::PointType itkfixedorigin(0.0f);
    ImageType::PointType itkmovingorigin(0.0f);

    ImageType::SpacingType itkfixedspacing(1.0f);
    ImageType::SpacingType itkmovingspacing(1.0f);

    itkFixedImage->SetSpacing(itkfixedspacing);
    itkFixedImage->SetOrigin(itkfixedorigin);

    itkMovingImage->SetSpacing(itkmovingspacing);
    itkMovingImage->SetOrigin(itkmovingorigin);
  }
  else
  {
    ImageType::PointType itkfixedorigin;
    itkfixedorigin[0] = m_FixedOrigin.x; 
    itkfixedorigin[1] = m_FixedOrigin.y; 
    ImageType::PointType itkmovingorigin;
    itkmovingorigin[0] = m_MovingOrigin.x;
    itkmovingorigin[1] = m_MovingOrigin.y;

    ImageType::SpacingType itkfixedspacing;
    itkfixedspacing[0] = m_FixedSpacing.x;
    itkfixedspacing[1] = m_FixedSpacing.y;
    ImageType::SpacingType itkmovingspacing;
    itkmovingspacing[0] = m_MovingSpacing.x;
    itkmovingspacing[1] = m_MovingSpacing.y;

    itkFixedImage->SetSpacing(itkfixedspacing);
    itkFixedImage->SetOrigin(itkfixedorigin);

    itkMovingImage->SetSpacing(itkmovingspacing);
    itkMovingImage->SetOrigin(itkmovingorigin);
  }

  //////////////////////////////// REGISTRATION //////////////////////////////////////////////////

  using ImageRegistrationType = itk::ImageRegistrationMethodv4<ImageType, ImageType>;
  ImageRegistrationType::Pointer registration = ImageRegistrationType::New();
  registration->InPlaceOn();
  registration->SetFixedImage(itkFixedImage);
  registration->SetMovingImage(itkMovingImage);

  //////////////////////////////// TRANSFORM //////////////////////////////////////////////////

  using CoordinateRepType = double;
  using ImageType = itk::Dream3DImage<double, ImageDimension>;

  if(m_TransformType == 0)
  {
    if(m_BSplineOrder == 3)
    {
      const unsigned int SplineOrder = 3;
      BSPLINESETUP()
    }

    if(m_BSplineOrder == 2)
    {
      const unsigned int SplineOrder = 2;
      BSPLINESETUP()
    }

    if(m_BSplineOrder == 1)
    {
      const unsigned int SplineOrder = 1;
      BSPLINESETUP()
    }
  }

  if(m_TransformType == 1)
  {
    using TransformType = itk::Euler2DTransform<double>;
    TransformType::Pointer transform = TransformType::New();
    using TranformInitializerType = itk::CenteredTransformInitializer<TransformType, ImageType, ImageType>;
    TranformInitializerType::Pointer transformInitializer = TranformInitializerType::New();

    if(m_InitializeRigidWithFFT)
    {
      using FloatImageType = itk::Image<float, 2>;
      using FFTNXcorrType = itk::FFTNormalizedCorrelationImageFilter<ImageType, FloatImageType>;
      FFTNXcorrType::Pointer fftfilter = FFTNXcorrType::New();
      fftfilter->SetFixedImage(itkFixedImage);
      fftfilter->SetMovingImage(itkMovingImage);
      fftfilter->SetRequiredFractionOfOverlappingPixels(m_ReqFracOverlapPxlsRigid);
      // fftfilter->SetRequiredNumberOfOverlappingPixels(600000);
      fftfilter->Update();

      using MinimumMaximumImageCalculatorType = itk::MinimumMaximumImageCalculator<FloatImageType>;
      MinimumMaximumImageCalculatorType::Pointer minimumMaximumImageCalculatorFilter = MinimumMaximumImageCalculatorType::New();
      minimumMaximumImageCalculatorFilter->SetImage(fftfilter->GetOutput());
      minimumMaximumImageCalculatorFilter->Compute();
      itk::Index<2> maximumCorrelationPatchCenter = minimumMaximumImageCalculatorFilter->GetIndexOfMaximum();
      itk::Size<2> movingSize = itkMovingImage->GetLargestPossibleRegion().GetSize();

      itk::Index<2> maximumCorrelationPatchCenterFixed;
      maximumCorrelationPatchCenterFixed[0] = movingSize[0] - maximumCorrelationPatchCenter[0];
      maximumCorrelationPatchCenterFixed[1] = movingSize[1] - maximumCorrelationPatchCenter[1];
      TransformType::InputVectorType translation;

      translation[0] = maximumCorrelationPatchCenterFixed[0];
      translation[1] = maximumCorrelationPatchCenterFixed[1];

      transform->SetTranslation(translation);
    }

    transformInitializer->SetTransform(transform);
    transformInitializer->GeometryOn();
    transformInitializer->SetFixedImage(itkFixedImage);
    transformInitializer->SetMovingImage(itkMovingImage);
    // transformInitializer->GeometryOn();
    // transformInitializer->InitializeTransform();

    registration->SetInitialTransform(transform);
  }

  if(m_TransformType == 2)
  {
    using TransformType = itk::AffineTransform<double, ImageDimension>;
    TransformType::Pointer transform = TransformType::New();
    using TranformInitializerType = itk::CenteredTransformInitializer<TransformType, ImageType, ImageType>;
    TranformInitializerType::Pointer transformInitializer = TranformInitializerType::New();

    if(m_InitializeAffineWithFFT)
    {
      using FloatImageType = itk::Image<float, 2>;
      using FFTNXcorrType = itk::FFTNormalizedCorrelationImageFilter<ImageType, FloatImageType>;
      FFTNXcorrType::Pointer fftfilter = FFTNXcorrType::New();
      fftfilter->SetFixedImage(itkFixedImage);
      fftfilter->SetMovingImage(itkMovingImage);
      fftfilter->SetRequiredFractionOfOverlappingPixels(m_ReqFracOverlapPxlsAffine);
      fftfilter->Update();

      using MinimumMaximumImageCalculatorType = itk::MinimumMaximumImageCalculator<FloatImageType>;
      MinimumMaximumImageCalculatorType::Pointer minimumMaximumImageCalculatorFilter = MinimumMaximumImageCalculatorType::New();
      minimumMaximumImageCalculatorFilter->SetImage(fftfilter->GetOutput());
      minimumMaximumImageCalculatorFilter->Compute();
      itk::Index<2> maximumCorrelationPatchCenter = minimumMaximumImageCalculatorFilter->GetIndexOfMaximum();
      itk::Size<2> outputSize = fftfilter->GetOutput()->GetLargestPossibleRegion().GetSize();

      itk::Index<2> maximumCorrelationPatchCenterFixed;
      maximumCorrelationPatchCenterFixed[0] = outputSize[0] / 2 - maximumCorrelationPatchCenter[0];
      maximumCorrelationPatchCenterFixed[1] = outputSize[1] / 2 - maximumCorrelationPatchCenter[1];
      TransformType::InputVectorType translation;

      translation[0] = maximumCorrelationPatchCenterFixed[0];
      translation[1] = maximumCorrelationPatchCenterFixed[1];

      transform->SetTranslation(translation);
    }

    transformInitializer->SetTransform(transform);
    transformInitializer->GeometryOn();
    transformInitializer->SetFixedImage(itkFixedImage);
    transformInitializer->SetMovingImage(itkMovingImage);
    transformInitializer->GeometryOn();
    transformInitializer->InitializeTransform();

    registration->SetInitialTransform(transform);
  }

  //////////////////////////////// METRIC //////////////////////////////////////////////////

  // For each metric, an optimizer needs to be set. The optimizer needs to know the type of the metric, so there is no easy way to abstract this.
  // Instead each case must be explicitly written out. Here we use macros to help.

  using ImageType = itk::Dream3DImage<double, 2>;

  if(m_Metric == 0)
  {
    CHOOSE_METRIC(MattesMutualInformationImageToImageMetricv4);
  }

  if(m_Metric == 1)
  {
    CHOOSE_METRIC(MeanSquaresImageToImageMetricv4);
  }

  if(m_Metric == 2)
  {
    CHOOSE_METRIC(CorrelationImageToImageMetricv4);
  }

  //////////////////////////////// MUTLI-RESOLUTION //////////////////////////////////////////////////

  int numberOfLevels = 3;
  ImageRegistrationType::ShrinkFactorsArrayType shrinkFactorsPerLevel;
  shrinkFactorsPerLevel.SetSize(3);
  shrinkFactorsPerLevel[0] = 6;
  shrinkFactorsPerLevel[1] = 2;
  shrinkFactorsPerLevel[2] = 1;
  ImageRegistrationType::SmoothingSigmasArrayType smoothingSigmasPerLevel;
  smoothingSigmasPerLevel.SetSize(3);
  smoothingSigmasPerLevel[0] = 6;
  smoothingSigmasPerLevel[1] = 2;
  smoothingSigmasPerLevel[2] = 1;
  registration->SetNumberOfLevels(numberOfLevels);
  registration->SetSmoothingSigmasPerLevel(smoothingSigmasPerLevel);
  registration->SetShrinkFactorsPerLevel(shrinkFactorsPerLevel);

  //////////////////////////////// RUN REGISTRATION //////////////////////////////////////////////////

  registration->Update();

  /////////////////////////////// WRITE OUT TRANSFORM //////////////////////////////////////////////////
  int numFixedElements = registration->GetTransform()->GetFixedParameters().GetNumberOfElements();
  std::vector<double> fixedParameters(numFixedElements);
  for(int i = 0; i < numFixedElements; i++)
  {
    fixedParameters[i] = registration->GetTransform()->GetFixedParameters().GetElement(i);
  }

  int numElements = registration->GetTransform()->GetParameters().GetNumberOfElements();
  std::vector<double> learnedParameters(numElements);
  for(int i = 0; i < numElements; i++)
  {
    learnedParameters[i] = registration->GetTransform()->GetParameters().GetElement(i);
  }

  std::string transformTypeString = registration->GetTransform()->GetTransformTypeAsString();

  TransformContainer::Pointer transformContainer = TransformContainer::New();
  transformContainer->setFixedParameters(fixedParameters);
  transformContainer->setParameters(learnedParameters);
  transformContainer->setTransformTypeAsString(transformTypeString);

  ImageGeom::Pointer movingGeometry = dca->getDataContainer(m_MovingImageArrayPath.getDataContainerName())->getGeometryAs<ImageGeom>();

  movingGeometry->setTransformContainer(transformContainer);

  std::string groupName = std::to_string(sliceNo);
  transformContainer->writeTransformContainerToHDF5(fileID, groupName);
  writeFixedImageInfo2DtoHDF5(fileID, sliceNo, itkFixedImage);
}
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ITKPairwiseImageRegistration::SinglePairRegistration()
{
  EXECUTE_FUNCTION_TEMPLATE_NO_BOOL(this, convertToDouble, m_FixedImagePtr.lock(), m_FixedImagePtr, m_FixedImageArrayPath, "_INTERNAL_USE_ONLY_FixedImageDouble", getDataContainerArray());
  EXECUTE_FUNCTION_TEMPLATE_NO_BOOL(this, convertToDouble, m_MovingImagePtr.lock(), m_MovingImagePtr, m_MovingImageArrayPath, "_INTERNAL_USE_ONLY_MovingImageDouble", getDataContainerArray());
  hid_t fileId = QH5Utilities::createFile(m_TransformFile);
  registerImagePair2D(getDataContainerArray(), 0, fileId);
  H5Utilities::closeFile(fileId);
}
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int ITKPairwiseImageRegistration::writeFixedImageInfo2DtoHDF5(hid_t parentId, size_t sliceNo, itk::Dream3DImage<double, 2>::Pointer itkFixedImage)
{
  using ImageType = itk::Dream3DImage<double, 2>;

  ImageType::SizeType fixedsize = itkFixedImage->GetLargestPossibleRegion().GetSize();
  ImageType::PointType fixedorigin = itkFixedImage->GetOrigin();
  ImageType::SpacingType fixedspacing = itkFixedImage->GetSpacing();
  ImageType::DirectionType fixeddirection = itkFixedImage->GetDirection();

  std::vector<unsigned long long> sizeVector{fixedsize.GetElement(0), fixedsize.GetElement(1)};
  std::vector<double> originVector{fixedorigin.GetElement(0), fixedorigin.GetElement(1)};
  std::vector<double> spacingVector{fixedspacing.GetElement(0), fixedspacing.GetElement(1)};
  std::vector<double> directionVector{fixeddirection(0, 0), fixeddirection(0, 1), fixeddirection(1, 0), fixeddirection(1, 1)};

  std::string transformName = std::to_string(sliceNo);
  herr_t err = 0;

  hid_t transformId = H5Gopen(parentId, transformName.c_str(), H5P_DEFAULT);
  if(transformId < 0)
  {
    return -1;
  }
  H5ScopedGroupSentinel gSentinel(&transformId, false);
  std::vector<hsize_t> dims(1);

  dims[0] = 2;
  QString sizeString = "fixedSize";
  err = H5Lite::writeVectorDataset(transformId, sizeString.toLatin1().data(), dims, sizeVector);
  if(err < 0)
  {
    return err;
  }

  dims[0] = 2;
  QString originString = "fixedOrigin";
  err = H5Lite::writeVectorDataset(transformId, originString.toLatin1().data(), dims, originVector);
  if(err < 0)
  {
    return err;
  }

  dims[0] = 2;
  QString spacingString = "fixedSpacing";
  err = H5Lite::writeVectorDataset(transformId, spacingString.toLatin1().data(), dims, spacingVector);
  if(err < 0)
  {
    return err;
  }

  dims[0] = 4;
  QString directionString = "fixedDirection";
  err = H5Lite::writeVectorDataset(transformId, directionString.toLatin1().data(), dims, directionVector);
  if(err < 0)
  {
    return err;
  }

  return 0;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------

void ITKPairwiseImageRegistration::SeriesPairRegistration()
{

  QVector<QString> movingFileList = getFileList(m_MovingFileListInfo);
  QVector<QString> fixedFileList = getFileList(m_FixedFileListInfo);

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

  hid_t fileId = QH5Utilities::createFile(m_TransformFile);
  for(int i = 0; i < movingFileList.size(); i++)
  {

    ///////////////MOVING IMAGE////////////////////

    DataContainerArray::Pointer dca = DataContainerArray::New();

    itkImageReader->setDataContainerArray(dca);

    QVariant dcName;
    QString pathname = "_INTERNAL_USE_ONLY_MovingImageDataContainerName";
    dcName.setValue(pathname);
    itkImageReader->setProperty("DataContainerName", dcName);

    QVariant amName;
    QString name = "_INTERNAL_USE_ONLY_attributeMatrixName";
    amName.setValue(name);
    itkImageReader->setProperty("CellAttributeMatrixName", amName);

    QVariant imDAName;
    name = "_INTERNAL_USE_ONLY_imageDataArrayName";
    imDAName.setValue(name);
    itkImageReader->setProperty("ImageDataArrayName", imDAName);

    QVariant var;
    var.setValue(movingFileList[i]);
    itkImageReader->setProperty("FileName", var);

    // itkImageReader->preflight();
    itkImageReader->execute();

    if(itkImageReader->getErrorCondition() < 0)
    {
      setErrorCondition(itkImageReader->getErrorCondition());
      notifyErrorMessage(getHumanLabel(), "Error Reading Input Image.", getErrorCondition());

      return;
    }

    DataArrayPath movingpath("_INTERNAL_USE_ONLY_MovingImageDataContainerName", "_INTERNAL_USE_ONLY_attributeMatrixName", "_INTERNAL_USE_ONLY_imageDataArrayName");
    setMovingImageArrayPath(movingpath);

    ///////////////FIXED IMAGE////////////////////

    pathname = "_INTERNAL_USE_ONLY_FixedImageDataContainerName";
    dcName.setValue(pathname);
    itkImageReader->setProperty("DataContainerName", dcName);

    var.setValue(fixedFileList[i]);
    itkImageReader->setProperty("FileName", var);

    // itkImageReader->preflight();
    itkImageReader->execute();

    if(itkImageReader->getErrorCondition() < 0)
    {
      setErrorCondition(itkImageReader->getErrorCondition());
      notifyErrorMessage(getHumanLabel(), "Error Reading Input Image.", getErrorCondition());

      return;
    }

    DataArrayPath fixedpath("_INTERNAL_USE_ONLY_FixedImageDataContainerName", "_INTERNAL_USE_ONLY_attributeMatrixName", "_INTERNAL_USE_ONLY_imageDataArrayName");
    setFixedImageArrayPath(fixedpath);

    m_FixedImagePtr = dca->getPrereqIDataArrayFromPath<IDataArray, AbstractFilter>(this, getFixedImageArrayPath());
    m_MovingImagePtr = dca->getPrereqIDataArrayFromPath<IDataArray, AbstractFilter>(this, getMovingImageArrayPath()); /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */

    if(getErrorCondition() < 0)
    {
      return;
    }

    EXECUTE_FUNCTION_TEMPLATE_NO_BOOL(this, convertToDouble, m_FixedImagePtr.lock(), m_FixedImagePtr, m_FixedImageArrayPath, "_INTERNAL_USE_ONLY_FixedImageDouble", dca);
    EXECUTE_FUNCTION_TEMPLATE_NO_BOOL(this, convertToDouble, m_MovingImagePtr.lock(), m_MovingImagePtr, m_MovingImageArrayPath, "_INTERNAL_USE_ONLY_MovingImageDouble", dca);

    registerImagePair2D(dca, i, fileId);
  }

  H5Utilities::closeFile(fileId);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ITKPairwiseImageRegistration::execute()
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

  if(m_OperationMode == 0)
  {
    SinglePairRegistration();
  }

  if(m_OperationMode == 1)
  {
    SeriesPairRegistration();
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AbstractFilter::Pointer ITKPairwiseImageRegistration::newFilterInstance(bool copyFilterParameters) const
{
  ITKPairwiseImageRegistration::Pointer filter = ITKPairwiseImageRegistration::New();
  if(copyFilterParameters)
  {
    copyFilterParameterInstanceVariables(filter.get());
  }
  return filter;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString ITKPairwiseImageRegistration::getCompiledLibraryName() const
{
  return DREAM3DReviewConstants::DREAM3DReviewBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString ITKPairwiseImageRegistration::getBrandingString() const
{
  return "DREAM3DReview";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString ITKPairwiseImageRegistration::getFilterVersion() const
{
  QString version;
  QTextStream vStream(&version);
  vStream << DREAM3DReview::Version::Major() << "." << DREAM3DReview::Version::Minor() << "." << DREAM3DReview::Version::Patch();
  return version;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString ITKPairwiseImageRegistration::getGroupName() const
{
  return SIMPL::FilterGroups::Unsupported;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString ITKPairwiseImageRegistration::getSubGroupName() const
{
  return DREAM3DReviewConstants::FilterSubGroups::RegistrationFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString ITKPairwiseImageRegistration::getHumanLabel() const
{
  return "ITKPairwiseImageRegistration";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QUuid ITKPairwiseImageRegistration::getUuid()
{
  return QUuid("{0a283663-e1c8-596d-bcd0-fd27b747551d}");
}
