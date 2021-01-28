/* ============================================================================
 * Software developed by US federal government employees (including military personnel)
 * as part of their official duties is not subject to copyright protection and is
 * considered "public domain" (see 17 USC Section 105). Public domain software can be used
 * by anyone for any purpose, and cannot be released under a copyright license
 * (including typical open source software licenses).
 *
 * This source code file was originally written by United States DoD employees. The
 * original source code files are released into the Public Domain.
 *
 * Subsequent changes to the codes by others may elect to add a copyright and license
 * for those changes.
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
#ifndef _interpolatepointcloudtoregulargrid_h_
#define _interpolatepointcloudtoregulargrid_h_

#include <memory>

#include "SIMPLib/SIMPLib.h"
#include "SIMPLib/DataArrays/DataArray.hpp"
#include "SIMPLib/DataArrays/NeighborList.hpp"
#include "SIMPLib/FilterParameters/FloatVec3FilterParameter.h"
#include "SIMPLib/Filtering/AbstractFilter.h"
#include "SIMPLib/Geometry/IGeometry.h"
#include "SIMPLib/Geometry/ImageGeom.h"

#include "DREAM3DReview/DREAM3DReviewDLLExport.h"

/**
 * @brief The InterpolatePointCloudToRegularGrid class. See [Filter documentation](@ref interpolatepointcloudtoregulargrid) for details.
 */
class DREAM3DReview_EXPORT InterpolatePointCloudToRegularGrid : public AbstractFilter
{
  Q_OBJECT
  PYB11_BEGIN_BINDINGS(InterpolatePointCloudToRegularGrid SUPERCLASS AbstractFilter)
  PYB11_FILTER()
  PYB11_SHARED_POINTERS(InterpolatePointCloudToRegularGrid)
  PYB11_FILTER_NEW_MACRO(InterpolatePointCloudToRegularGrid)

  PYB11_END_BINDINGS()

public:
  using Self = InterpolatePointCloudToRegularGrid;
  using Pointer = std::shared_ptr<Self>;
  using ConstPointer = std::shared_ptr<const Self>;
  using WeakPointer = std::weak_ptr<Self>;
  using ConstWeakPointer = std::weak_ptr<const Self>;
  static Pointer NullPointer();

  static std::shared_ptr<InterpolatePointCloudToRegularGrid> New();

  /**
   * @brief Returns the name of the class for InterpolatePointCloudToRegularGrid
   */
  QString getNameOfClass() const override;
  /**
   * @brief Returns the name of the class for InterpolatePointCloudToRegularGrid
   */
  static QString ClassName();

  ~InterpolatePointCloudToRegularGrid() override;

  /**
   * @brief Setter property for DataContainerName
   */
  void setDataContainerName(const DataArrayPath& value);
  /**
   * @brief Getter property for DataContainerName
   * @return Value of DataContainerName
   */
  DataArrayPath getDataContainerName() const;
  Q_PROPERTY(DataArrayPath DataContainerName READ getDataContainerName WRITE setDataContainerName)

  /**
   * @brief Setter property for ArraysToInterpolate
   */
  void setArraysToInterpolate(const QVector<DataArrayPath>& value);
  /**
   * @brief Getter property for ArraysToInterpolate
   * @return Value of ArraysToInterpolate
   */
  QVector<DataArrayPath> getArraysToInterpolate() const;
  Q_PROPERTY(QVector<DataArrayPath> ArraysToInterpolate READ getArraysToInterpolate WRITE setArraysToInterpolate)

  /**
   * @brief Setter property for ArraysToCopy
   */
  void setArraysToCopy(const QVector<DataArrayPath>& value);
  /**
   * @brief Getter property for ArraysToCopy
   * @return Value of ArraysToCopy
   */
  QVector<DataArrayPath> getArraysToCopy() const;
  Q_PROPERTY(QVector<DataArrayPath> ArraysToCopy READ getArraysToCopy WRITE setArraysToCopy)

  /**
   * @brief Setter property for VoxelIndicesArrayPath
   */
  void setVoxelIndicesArrayPath(const DataArrayPath& value);
  /**
   * @brief Getter property for VoxelIndicesArrayPath
   * @return Value of VoxelIndicesArrayPath
   */
  DataArrayPath getVoxelIndicesArrayPath() const;
  Q_PROPERTY(DataArrayPath VoxelIndicesArrayPath READ getVoxelIndicesArrayPath WRITE setVoxelIndicesArrayPath)

  /**
   * @brief Setter property for InterpolatedDataContainerName
   */
  void setInterpolatedDataContainerName(const DataArrayPath& value);
  /**
   * @brief Getter property for InterpolatedDataContainerName
   * @return Value of InterpolatedDataContainerName
   */
  DataArrayPath getInterpolatedDataContainerName() const;
  Q_PROPERTY(DataArrayPath InterpolatedDataContainerName READ getInterpolatedDataContainerName WRITE setInterpolatedDataContainerName)

  /**
   * @brief Setter property for InterpolatedAttributeMatrixName
   */
  void setInterpolatedAttributeMatrixName(const QString& value);
  /**
   * @brief Getter property for InterpolatedAttributeMatrixName
   * @return Value of InterpolatedAttributeMatrixName
   */
  QString getInterpolatedAttributeMatrixName() const;
  Q_PROPERTY(QString InterpolatedAttributeMatrixName READ getInterpolatedAttributeMatrixName WRITE setInterpolatedAttributeMatrixName)

  /**
   * @brief Setter property for InterpolationTechnique
   */
  void setInterpolationTechnique(int value);
  /**
   * @brief Getter property for InterpolationTechnique
   * @return Value of InterpolationTechnique
   */
  int getInterpolationTechnique() const;
  Q_PROPERTY(int InterpolationTechnique READ getInterpolationTechnique WRITE setInterpolationTechnique)

  /**
   * @brief Setter property for KernelSize
   */
  void setKernelSize(const FloatVec3Type& value);
  /**
   * @brief Getter property for KernelSize
   * @return Value of KernelSize
   */
  FloatVec3Type getKernelSize() const;
  Q_PROPERTY(FloatVec3Type KernelSize READ getKernelSize WRITE setKernelSize)

  /**
   * @brief Setter property for Sigmas
   */
  void setSigmas(const FloatVec3Type& value);
  /**
   * @brief Getter property for Sigmas
   * @return Value of Sigmas
   */
  FloatVec3Type getSigmas() const;
  Q_PROPERTY(FloatVec3Type Sigmas READ getSigmas WRITE setSigmas)

  /**
   * @brief Setter property for UseMask
   */
  void setUseMask(bool value);
  /**
   * @brief Getter property for UseMask
   * @return Value of UseMask
   */
  bool getUseMask() const;
  Q_PROPERTY(bool UseMask READ getUseMask WRITE setUseMask)

  /**
   * @brief Setter property for MaskArrayPath
   */
  void setMaskArrayPath(const DataArrayPath& value);
  /**
   * @brief Getter property for MaskArrayPath
   * @return Value of MaskArrayPath
   */
  DataArrayPath getMaskArrayPath() const;
  Q_PROPERTY(DataArrayPath MaskArrayPath READ getMaskArrayPath WRITE setMaskArrayPath)

  /**
   * @brief Setter property for StoreKernelDistances
   */
  void setStoreKernelDistances(bool value);
  /**
   * @brief Getter property for StoreKernelDistances
   * @return Value of StoreKernelDistances
   */
  bool getStoreKernelDistances() const;
  Q_PROPERTY(bool StoreKernelDistances READ getStoreKernelDistances WRITE setStoreKernelDistances)

  /**
   * @brief Setter property for KernelDistancesArrayName
   */
  void setKernelDistancesArrayName(const QString& value);
  /**
   * @brief Getter property for KernelDistancesArrayName
   * @return Value of KernelDistancesArrayName
   */
  QString getKernelDistancesArrayName() const;
  Q_PROPERTY(QString KernelDistancesArrayName READ getKernelDistancesArrayName WRITE setKernelDistancesArrayName)

  /**
   * @brief Setter property for KernelDistancesArrayName
   */
  void setInterpolatedSuffix(const QString& value);
  /**
   * @brief Getter property for KernelDistancesArrayName
   * @return Value of KernelDistancesArrayName
   */
  QString getInterpolatedSuffix() const;
  Q_PROPERTY(QString InterpolatedSuffix READ getInterpolatedSuffix WRITE setInterpolatedSuffix)

  /**
   * @brief Setter property for KernelDistancesArrayName
   */
  void setCopySuffix(const QString& value);
  /**
   * @brief Getter property for KernelDistancesArrayName
   * @return Value of KernelDistancesArrayName
   */
  QString getCopySuffix() const;
  Q_PROPERTY(QString CopySuffix READ getCopySuffix WRITE setCopySuffix)

  /**
   * @brief getCompiledLibraryName Reimplemented from @see AbstractFilter class
   */
  QString getCompiledLibraryName() const override;

  /**
   * @brief getBrandingString Returns the branding string for the filter, which is a tag
   * used to denote the filter's association with specific plugins
   * @return Branding string
   */
  QString getBrandingString() const override;

  /**
   * @brief getFilterVersion Returns a version string for this filter. Default
   * value is an empty string.
   * @return
   */
  QString getFilterVersion() const override;

  /**
   * @brief newFilterInstance Reimplemented from @see AbstractFilter class
   */
  AbstractFilter::Pointer newFilterInstance(bool copyFilterParameters) const override;

  /**
   * @brief getGroupName Reimplemented from @see AbstractFilter class
   */
  QString getGroupName() const override;

  /**
   * @brief getSubGroupName Reimplemented from @see AbstractFilter class
   */
  QString getSubGroupName() const override;

  /**
   * @brief getHumanLabel Reimplemented from @see AbstractFilter class
   */
  QString getHumanLabel() const override;

  /**
   * @brief setupFilterParameters Reimplemented from @see AbstractFilter class
   */
  void setupFilterParameters() override;

  /**
   * @brief execute Reimplemented from @see AbstractFilter class
   */
  void execute() override;

  /**
   * @brief getUuid Return the unique identifier for this filter.
   * @return A QUuid object.
   */
  QUuid getUuid() const override;

protected:
  InterpolatePointCloudToRegularGrid();

  /**
   * @brief determineKernel Determines the kernel to be used for interpolation
   * @param kernelNumVoxels Voxel extents of the kernel
   * @return Float vector for the kernel values
   */
  void determineKernel(int64_t kernelNumVoxels[3]);

  /**
   * @brief determineKernelDistances Determines the Euclidean distance between each
   * voxel in the kernel and the kernel center
   * @param kernelNumVoxels Voxel extents of the kernel
   * @return Float vector for the kernel distances
   */
  void determineKernelDistances(int64_t kernelNumVoxels[3], FloatVec3Type res);

  /**
   * @brief mapKernelDistances Maps kernel point distances during interpolation
   * @param kernel Voxel dimensions of the kernel
   * @param dims Total dimensions of the interpolation grid
   * @param curX Current x position
   * @param curY Current y position
   * @param curZ Current z position
   */
  void mapKernelDistances(int64_t kernel[3], size_t dims[3], size_t curX, size_t curY, size_t curZ);
  /**
   * @brief dataCheck Checks for the appropriate parameter values and availability of arrays
   */
  void dataCheck() override;

  /**
   * @brief Initializes all the private instance variables.
   */
  void initialize();

private:
  std::weak_ptr<DataArray<MeshIndexType>> m_VoxelIndicesPtr;
  MeshIndexType* m_VoxelIndices = nullptr;
  std::weak_ptr<DataArray<bool>> m_MaskPtr;
  bool* m_Mask = nullptr;

  DataArrayPath m_DataContainerName = {};
  QVector<DataArrayPath> m_ArraysToInterpolate = {};
  QVector<DataArrayPath> m_ArraysToCopy = {};
  DataArrayPath m_VoxelIndicesArrayPath = {};
  DataArrayPath m_InterpolatedDataContainerName = {};
  QString m_InterpolatedAttributeMatrixName = {};
  int m_InterpolationTechnique = 0;
  FloatVec3Type m_KernelSize = {1.0f, 1.0f, 1.0f};
  FloatVec3Type m_Sigmas = {1.0f, 1.0f, 1.0f};
  bool m_UseMask = false;
  DataArrayPath m_MaskArrayPath = {};
  bool m_StoreKernelDistances = false;
  QString m_KernelDistancesArrayName = {};

  QString m_InterpolatedSuffix = " [Interpolated]";
  QString m_CopySuffix = " [Copied]";

  NeighborList<float>::WeakPointer m_KernelDistances = NeighborList<float>::NullPointer();

  QList<QString> m_AttrMatList;
  std::vector<IDataArray::WeakPointer> m_SourceArraysToInterpolate;
  std::vector<IDataArray::WeakPointer> m_SourceArraysToCopy;
  std::vector<IDataArray::WeakPointer> m_DynamicArraysToInterpolate;
  std::vector<IDataArray::WeakPointer> m_DynamicArraysToCopy;
  std::vector<float> m_Kernel;
  std::vector<float> m_KernelValDistances;

public:
  InterpolatePointCloudToRegularGrid(const InterpolatePointCloudToRegularGrid&) = delete;            // Copy Constructor Not Implemented
  InterpolatePointCloudToRegularGrid(InterpolatePointCloudToRegularGrid&&) = delete;                 // Move Constructor Not Implemented
  InterpolatePointCloudToRegularGrid& operator=(const InterpolatePointCloudToRegularGrid&) = delete; // Copy Assignment Not Implemented
  InterpolatePointCloudToRegularGrid& operator=(InterpolatePointCloudToRegularGrid&&) = delete;      // Move Assignment Not Implemented
};

#endif /* _InterpolatePointCloudToRegularGrid_H_ */
