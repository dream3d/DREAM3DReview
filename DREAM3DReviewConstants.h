#pragma once

#include <QtCore/QString>

/** @brief This macro is used to shorten the code needed to go from std::string to QString. Helpful in other codes that use
 * QString instead of std::string
 */
#define S2Q(var) QString::fromStdString((var))

/**
 * @brief This namespace is used to define some Constants for the plugin itself.
 */
namespace DREAM3DReviewConstants
{

const QString DREAM3DReviewPluginFile("DREAM3DReviewPlugin");
const QString DREAM3DReviewPluginDisplayName("DREAM3D Review Plugin");
const QString DREAM3DReviewBaseName("DREAM3DReview");

namespace FilterGroups
{
const QString DREAM3DReviewFilters("DREAM3D Review");
}

namespace FilterSubGroups
{
const QString RotationTransformationFilters("Rotation/Transforming");
const QString StatisticsFilters("Statistics");
const QString RegistrationFilters("Registration");
const QString MemoryManagementFilters("Memory/Management");
const QString GenerationFilters("Generation");
const QString CropCutFilters("Cropping/Cutting");
const QString ClusteringFilters("Clustering");
const QString GeometryFilters("Geometry");
const QString DimensionalityReductionFilters("Dimensionality Reduction");
const QString ThresholdFilters("Threshold");
const QString Coarsening("Coarsening");
const QString InterpolationFilters("InterpolationFilters");
const QString PointCloudFilters("PointCloudFilters");

} // namespace FilterSubGroups
} // namespace DREAM3DReviewConstants
/**
 * @brief This namespace is used to define some Constants for the plugin itself.
 */
namespace MASSIFUtilitiesConstants
{
const QString MASSIFUtilitiesPluginFile("MASSIFUtilitiesPlugin");
const QString MASSIFUtilitiesPluginDisplayName("MASSIF Utilities Plugin");
const QString MASSIFUtilitiesBaseName("MASSIFUtilities");

namespace FilterGroups
{
const QString MASSIFUtilitiesFilters("MASSIFUtilities");
}

namespace ImportMassifData
{
const QString MassifDC = "MassifDataContainer";
const QString MassifAM = "MassifAttributeMatrix";
const QString DField = "D11";
const QString EField = "El11";
const QString SField = "S11";
const QString DimGrpName = "Dimension";
const QString OriginGrpName = "Origin";
const QString SpacingGrpName = "Spacing";
const QString DCGrpName = "3Ddatacontainer";
const QString GeometryGrpName = "Geometry";
const QString EVM = "EVM";
const QString SVM = "SVM";
const QString GrainID = "Grainid";
const QString Phase = "Phase";
const QString Datapoint = "Datapoint";
const QString DFieldsGrpName = "Dfields";
const QString EFieldsGrpName = "Elastic strain";
const QString SFieldsGrpName = "Sfields";
const QString EulerAngleGrpName = "Eulerangle";
const QString Phi1 = "Phi1";
const QString Phi = "Phi";
const QString Phi2 = "phi2";

const int MaxStepNumber = 999999;
} // namespace ImportMassifData
} // namespace MASSIFUtilitiesConstants

namespace AnisotropyConstants
{
const QString AnisotropyPluginFile("AnisotropyPlugin");
const QString AnisotropyPluginDisplayName("Anisotropy");
const QString AnisotropyBaseName("Anisotropy");

const QString Version("1.0");
const QString CompatibilityVersion("1.0");

const QString VendorName("Czech Academy of Sciences, Institute of Physics, Group of Bulk Nanomaterials and Interfaces");
const QString URL("http://ams.fzu.cz");
const QString Copyright("(C) 2016 Czech Academy of Sciences, v.v.i.");

namespace FilterSubGroups
{
const QString AnisotropicAlignment("Anisotropic Alignment");
}

} // namespace AnisotropyConstants

/**
 * @brief This namespace is used to define some Constants for the plugin itself.
 */
namespace HEDMAnalysisConstants
{
const QString HEDMAnalysisPluginFile("HEDMAnalysisPlugin");
const QString HEDMAnalysisPluginDisplayName("HEDMAnalysis");
const QString HEDMAnalysisBaseName("HEDMAnalysis");

namespace FilterGroups
{
const QString HEDMAnalysisFilters("HEDM Analysis");
}
} // namespace HEDMAnalysisConstants

/**
 * @brief This namespace is used to define some Constants for the plugin itself.
 */
namespace DDDAnalysisToolboxConstants
{
const QString DDDAnalysisToolboxPluginFile("DDDAnalysisToolboxPlugin");
const QString DDDAnalysisToolboxPluginDisplayName("DDDAnalysisToolbox");
const QString DDDAnalysisToolboxBaseName("DDDAnalysisToolbox");

namespace FilterGroups
{
const QString DDDAnalyticsToolboxFilters("DDD Analytics");
}
} // namespace DDDAnalysisToolboxConstants

namespace TransformationPhaseConstants
{
const QString TransformationPhasePluginFile("TransformationPhasePlugin");
const QString TransformationPhasePluginDisplayName("Transformation Phase");
// const QString TransformationPhaseBaseName("TransformationPhase");

namespace FilterGroups
{
const QString TransformationPhaseFilters("Transformation Phase");
}

const QString Initiators("Initiators");
const QString HardFeatures("HardFeatures");
const QString SoftFeatures("SoftFeatures");
const QString HardSoftGroups("HardSoftGroups");
const QString SelectedFeatures("SelectedFeatures");
const QString SubsurfaceFeatures("SubsurfaceFeatures");

const QString SurfaceMeshCSLBoundary("SurfaceMeshCSLBoundary");
const QString SurfaceMeshCSLBoundaryIncoherence("SurfaceMeshCSLBoundaryIncoherence");

namespace CSL
{

const float Sigma3 = 3.0f;
const float Sigma5 = 5.0f;
const float Sigma7 = 7.0f;
const float Sigma9 = 9.0f;
const float Sigma11a = 11.0f;
const float Sigma11b = 11.5f;
const float Sigma13a = 13.0f;
const float Sigma13b = 13.5f;
const float Sigma15 = 15.0f;
const float Sigma17a = 17.0f;
const float Sigma17b = 17.5f;
const float Sigma19a = 19.0f;
const float Sigma19b = 19.5f;
const float Sigma21a = 21.0f;
const float Sigma21b = 21.5f;
const float Sigma23 = 23.0f;
const float Sigma25a = 25.0f;
const float Sigma25b = 25.5f;
const float Sigma27a = 27.0f;
const float Sigma27b = 27.5f;
const float Sigma29a = 29.0f;
const float Sigma29b = 29.5f;

} // namespace CSL

static const float CSLAxisAngle[21][5] = {{3.0f, 60.0f, 1.0f, 1.0f, 1.0f},   {5.0f, 36.87f, 1.0f, 0.0f, 0.0f},  {7.0f, 38.21f, 1.0f, 1.0f, 0.0f},  {9.0f, 38.94f, 1.0f, 1.0f, 0.0f},
                                          {11.0f, 50.48f, 1.0f, 1.0f, 0.0f}, {13.0f, 22.62f, 1.0f, 0.0f, 0.0f}, {13.5f, 27.8f, 1.0f, 1.0f, 1.0f},  {15.0f, 48.19f, 2.0f, 1.0f, 0.0f},
                                          {17.0f, 28.07f, 1.0f, 0.0f, 0.0f}, {17.5f, 61.93f, 2.0f, 2.0f, 1.0f}, {19.0f, 26.53f, 1.0f, 1.0f, 0.0f}, {19.5f, 46.83f, 1.0f, 1.0f, 1.0f},
                                          {21.0f, 21.79f, 1.0f, 1.0f, 1.0f}, {21.5f, 44.4f, 2.0f, 1.0f, 1.0f},  {23.0f, 40.45f, 3.0f, 1.0f, 1.0f}, {25.0f, 16.25f, 1.0f, 0.0f, 0.0f},
                                          {25.5f, 51.68f, 3.0f, 3.0f, 1.0f}, {27.0f, 31.58f, 1.0f, 1.0f, 0.0f}, {27.5f, 35.42f, 2.0f, 1.0f, 0.0f}, {29.0f, 43.61f, 1.0f, 0.0f, 0.0f},
                                          {29.5f, 46.39f, 2.0f, 2.0f, 1.0f}};
} // namespace TransformationPhaseConstants

namespace UUtahDMREFConstants
{
const QString UUtahDMREFPluginFile("UUtahDMREFPlugin");
const QString UUtahDMREFPluginDisplayName("UUtahDMREF");
const QString UUtahDMREFBaseName("UUtahDMREF");

namespace FilterGroups
{
const QString UUtahDMREFFilters("UUtah DMREF");
}

namespace FilterSubGroups
{
const QString PackingFilters("Packing");
}
} // namespace UUtahDMREFConstants

/**
 * @brief Use this namespace to define any custom GUI widgets that collect FilterParameters
 * for a filter. Do NOT define general reusable widgets here.
 */
namespace FilterParameterWidgetType
{
/* const QString SomeCustomWidget("SomeCustomWidget"); */
}
